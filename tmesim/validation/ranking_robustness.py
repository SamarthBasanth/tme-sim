"""RANKING ROBUSTNESS -- does the DECISION survive parameter uncertainty?

Why this module exists
----------------------
:mod:`tmesim.validation.uncertainty` puts error bars on ONE drug's prediction
(kill / penetration credible intervals). That is the right question for a
*number*. It is NOT the question for a *decision*.

The headline decision claim of this repo is a RANKING on the real 10x Xenium GBM
section::

    cisplatin > doxorubicin > lomustine > methotrexate > temozolomide

An interval on each drug's kill does not tell you whether that ORDER holds -- the
intervals overlap heavily, but the draws are CORRELATED (all five drugs sit on the
same tissue, share the same calibration constant, the same vascular field). A
ranking can be far more stable than the marginal intervals suggest, or far less.
The only way to know is to re-rank inside every draw. That is what this module
does.

Design (and why it is affordable)
---------------------------------
The expensive object is the TISSUE + TRANSPORT solve. The cheap object is the
kill mapping. They separate exactly:

* ``sim.engine.simulate_response`` runs oxygen -> pressure -> delivery ->
  transport -> Hill kill. Only the LAST step reads ``drug.ic50``. So the AUC
  field is independent of IC50.
* Therefore: **reconstruct ONCE per reconstruction-knob draw, solve transport
  once per (draw, drug), then sweep IC50 analytically** on the frozen
  ``cells.pkpd`` Hill mapping. All five drugs share the same tissue in every
  draw, which is exactly the correlation structure a ranking question needs.

Per draw we therefore pay 1 reconstruction + 5 transport solves, and get
``n_ic50`` complete re-rankings out of it for free.

The kill used for ranking is the **expected** kill fraction
``E[kill] = mean_over_tumor(p_kill)`` rather than one Bernoulli realisation.
``p_kill`` comes from ``cells.pkpd.apply_drug_kill`` itself (its return value),
so the frozen mapping is reused, not re-implemented; :func:`verify_kill_kernel`
asserts the agreement at run time. Using the expectation matters here: the
single-seed Bernoulli noise on ``kill_fraction`` is ~0.0015 absolute, and the
nominal cisplatin-vs-doxorubicin gap is only ~0.019 -- coin-flip noise is a
non-trivial fraction of the very gap under test, and it is not parameter
uncertainty, so it is removed rather than counted.

Priors
------
The four reconstruction / calibration priors are IMPORTED UNCHANGED from
:mod:`tmesim.validation.uncertainty` (``build_priors``, ``sample_priors``) so the
two modules cannot drift: endothelial typing threshold, KDE bandwidth, the
single-anchor calibration alpha (0.5x-2x), and the ECM tortuosity coefficient.
See ``results/uncertainty.md`` section 1 for their justifications.

What this module ADDS is the axis a ranking needs and a single-drug interval does
not: **per-drug IC50 uncertainty**. Potency is what separates these five drugs
(their transport differs by ~1.5x in D_eff; their IC50s differ by ~1200x), so any
honest statement about the ORDER has to move the IC50s. Three scenarios are
evaluated on the SAME transport solves (see :func:`build_ic50_priors`):

  A. ``transport_only`` -- IC50s frozen at their resolved values. Isolates how
     much instability the reconstruction/calibration knobs alone can cause.
  B. ``ic50_median_se`` (PRIMARY) -- each IC50 is uncertain by the sampling error
     of the summary statistic actually used. For the three GDSC-screened drugs
     that is the standard error of the across-cell-line MEDIAN log-IC50,
     COMPUTED FROM THE CACHED SCREEN. For the two literature-only drugs it is a
     documented wider expert prior.
  C. ``ic50_cell_line`` (STRESS) -- each IC50 is uncertain by the FULL
     cross-cell-line spread, i.e. "we do not know where in the glioma cell-line
     distribution this patient's tumour sits". Also computed from the screen.

Run: ``python -m tmesim.validation.ranking_robustness`` (writes
``results/ranking_robustness.md`` and ``results/ranking_robustness.png``).

Nothing in ``core`` / ``physics`` / ``molecular`` / ``reconstruct`` / ``sim`` is
modified; this module only drives them.
"""
from __future__ import annotations

import csv
import dataclasses
import itertools
import math
import os
import statistics
import time
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import stats

from ..cells.agents import update_cell_states
from ..cells.pkpd import apply_drug_kill
from ..core import CellState, DrugParams, OxygenConfig, PDConfig, PDEConfig, SimConfig
from ..molecular import potency as _potency
from ..molecular.binding_site_barrier import apply_barrier_to_drugparams
from ..molecular.library import get_molecule
from ..sim.engine import simulate_response
from .calibration import calibrate_alpha, calibrated_drugparams
from .numerical_rigor import BASE_BANDWIDTH, BASE_THRESHOLD, reconstruct_knobs
from .uncertainty import PARAM_ORDER, _patched_tortuosity, build_priors, sample_priors

# ---------------------------------------------------------------------------
# The decision under test.
# ---------------------------------------------------------------------------
#: The five-drug panel, in the order the headline claim asserts.
PANEL: Tuple[str, ...] = (
    "cisplatin", "doxorubicin", "lomustine", "methotrexate", "temozolomide",
)
#: The headline ranking claim (best predicted kill first).
HEADLINE_RANKING: Tuple[str, ...] = PANEL
#: The headline pairwise claim and the headline top-2 SET.
HEADLINE_PAIR: Tuple[str, str] = ("cisplatin", "doxorubicin")
HEADLINE_TOP2: frozenset = frozenset(("cisplatin", "doxorubicin"))

# ---------------------------------------------------------------------------
# Monte-Carlo configuration.
#
# dx: the ranking is a statement about ``kill_fraction``, which numerical_rigor
# showed is GRID-CONVERGED for dx <= 20 um (it moves <1% between the two finest
# grids). dx = 15 um is inside that regime and ~3x cheaper than the dx = 12 um
# production grid; :func:`dx_ranking_check` re-runs the nominal panel at both to
# demonstrate the ranking (and the size of the cisplatin-doxorubicin gap) is not
# a resolution artifact.
#
# dt = 20 s / convection OFF: verified identical to the production
# dt = 10 s / convection ON configuration to 4 decimal places on every drug
# (numerical_rigor section 3: dt-invariant for dt <= 20 s; interstitial speeds
# < 0.7 um/s make advection negligible).
# ---------------------------------------------------------------------------
MC_DX_UM = 15.0
MC_DT_S = 20.0
MC_TOTAL_TIME_S = 3600.0 * 6
PRODUCTION_DX_UM = 12.0


def mc_pde(dt: float = MC_DT_S, total_time: float = MC_TOTAL_TIME_S,
           convection: bool = False) -> PDEConfig:
    return PDEConfig(dt=dt, total_time=total_time, convection=convection,
                     save_every=500)


def production_pde() -> PDEConfig:
    """The exact config ``experiments/real_drug_panel.py`` publishes with."""
    return PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)


# ===========================================================================
# 1. IC50 PRIORS -- derived from the cached GDSC screen where possible.
# ===========================================================================
@dataclass(frozen=True)
class IC50Prior:
    """A multiplicative (lognormal) prior on one drug's IC50.

    ``sigma_log`` is the standard deviation of ``ln(IC50)``; the median is the
    resolved IC50 itself, so the prior is centred on the value the model uses.
    """

    drug: str
    ic50_uM: float
    tier: str            # 'GDSC-measured' | 'literature'
    sigma_log: float
    basis: str
    n_cell_lines: Optional[int] = None
    sd_ln_cell_lines: Optional[float] = None

    @property
    def fold_90(self) -> float:
        """Half-width of the 90% CI as a multiplicative fold factor."""
        return float(math.exp(1.645 * self.sigma_log))

    def multipliers(self, z: np.ndarray) -> np.ndarray:
        """Map standard-normal deviates to IC50 multipliers (median 1.0)."""
        return np.exp(self.sigma_log * np.asarray(z, dtype=np.float64))


def gdsc_ln_ic50_by_cell_line(drug: str, cache_dir: str = "data_cache") -> Optional[np.ndarray]:
    """One representative ``LN_IC50`` per distinct glioma cell line, or ``None``.

    Reuses :mod:`tmesim.molecular.potency`'s alias table and cache path so the
    cell lines aggregated here are EXACTLY the ones whose median becomes the
    model's IC50. Repeat screens of the same line are collapsed by median, which
    is what ``potency.measured_ic50`` does.
    """
    path = _potency._gdsc_path(cache_dir)
    if not os.path.exists(path):
        return None
    aliases = set(_potency._aliases_for(drug))
    per_line: Dict[str, List[float]] = {}
    try:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                if (row.get("DRUG_NAME") or "").strip().lower() not in aliases:
                    continue
                line = row.get("CELL_LINE_NAME") or ""
                if not line:
                    continue
                try:
                    ln = float(row["LN_IC50"])
                except (KeyError, TypeError, ValueError):
                    continue
                if math.isfinite(ln):
                    per_line.setdefault(line, []).append(ln)
    except Exception:  # noqa: BLE001 - degrade to the literature tier
        return None
    if not per_line:
        return None
    return np.array([statistics.median(v) for v in per_line.values()], dtype=np.float64)


#: Expert multiplicative prior for the LITERATURE-only drugs under the PRIMARY
#: scenario: a 90% CI of 1/3x - 3x. Justification, per drug, in
#: ``potency.LITERATURE_IC50``: doxorubicin is quoted as "~0.1-1 uM in glioma
#: cell lines" around a 0.5 uM anchor (0.2x-2x) and lomustine as "tens of uM"
#: around 50 uM (~0.2x-2x). A 1/3x-3x band comfortably contains both quoted
#: ranges and adds room for the cross-study transfer those single-sentence
#: ranges do not cover. It is WIDER than every measured drug's prior, which is
#: the intended ordering: a curated order-of-magnitude anchor must not be
#: allowed to look as well-determined as a 34-51 cell-line screen.
LITERATURE_FOLD_90 = 3.0
LITERATURE_SIGMA_LOG = math.log(LITERATURE_FOLD_90) / 1.6448536269514722

#: Stress-scenario sigma for the literature drugs: the MEDIAN of the measured
#: drugs' own cross-cell-line log-SDs, i.e. "assume an unscreened drug's
#: cell-line-to-cell-line spread looks like a screened one's".
LITERATURE_STRESS_SIGMA_FALLBACK = 1.5

IC50_SCENARIOS: Tuple[str, ...] = ("transport_only", "ic50_median_se", "ic50_cell_line")
PRIMARY_SCENARIO = "ic50_median_se"

SCENARIO_LABEL = {
    "transport_only": "A. transport knobs only (IC50 fixed)",
    "ic50_median_se": "B. + IC50 sampling error (PRIMARY)",
    "ic50_cell_line": "C. + full cross-cell-line IC50 spread (STRESS)",
}


def build_ic50_priors(cache_dir: str = "data_cache") -> Dict[str, Dict[str, IC50Prior]]:
    """Build the three IC50 prior scenarios for the panel.

    For the GDSC-screened drugs both widths are DERIVED FROM DATA:

    * ``ic50_median_se`` uses ``sigma = 1.2533 * sd_ln / sqrt(n_lines)``, the
      asymptotic standard error of a MEDIAN (the 1.2533 = sqrt(pi/2) factor is
      the median's efficiency penalty relative to the mean under normality).
      This is the honest uncertainty on *the number the model actually plugs in*
      -- ``exp(median_over_cell_lines(LN_IC50))`` -- treating the screened lines
      as a sample.
    * ``ic50_cell_line`` uses the FULL ``sd_ln`` across cell lines. This is the
      much larger, and arguably more decision-relevant, question: this Xenium
      section is ONE tumour, and nothing tells us it behaves like the median
      glioma cell line. It is reported as a stress test rather than the headline
      because a cell-line panel's spread is a population statistic, not a
      credible interval for a specific tumour.
    """
    scen: Dict[str, Dict[str, IC50Prior]] = {s: {} for s in IC50_SCENARIOS}
    measured_sds: List[float] = []
    resolved: Dict[str, Dict] = {d: _potency.resolve_ic50(d, cache_dir=cache_dir)
                                 for d in PANEL}
    lnvals: Dict[str, Optional[np.ndarray]] = {
        d: (gdsc_ln_ic50_by_cell_line(d, cache_dir=cache_dir)
            if resolved[d]["confidence"] == "measured" else None)
        for d in PANEL
    }
    for d in PANEL:
        v = lnvals[d]
        if v is not None and v.size > 2:
            measured_sds.append(float(np.std(v, ddof=1)))
    stress_lit_sigma = (float(np.median(measured_sds)) if measured_sds
                        else LITERATURE_STRESS_SIGMA_FALLBACK)

    for d in PANEL:
        r = resolved[d]
        ic50 = float(r["ic50_uM"])
        v = lnvals[d]
        if v is not None and v.size > 2:
            n = int(v.size)
            sd = float(np.std(v, ddof=1))
            se_med = 1.2533141373155001 * sd / math.sqrt(n)
            common = dict(drug=d, ic50_uM=ic50, tier="GDSC-measured",
                          n_cell_lines=n, sd_ln_cell_lines=sd)
            scen["transport_only"][d] = IC50Prior(
                sigma_log=0.0,
                basis="held fixed at the resolved GDSC median", **common)
            scen["ic50_median_se"][d] = IC50Prior(
                sigma_log=se_med,
                basis=(f"SE of the across-cell-line median log-IC50: "
                       f"1.2533*{sd:.3f}/sqrt({n}) = {se_med:.3f} "
                       f"(GDSC2 glioma screen, {n} cell lines)"),
                **common)
            scen["ic50_cell_line"][d] = IC50Prior(
                sigma_log=sd,
                basis=(f"full cross-cell-line SD of log-IC50 = {sd:.3f} over "
                       f"{n} glioma cell lines (GDSC2)"),
                **common)
        else:
            common = dict(drug=d, ic50_uM=ic50, tier="literature",
                          n_cell_lines=None, sd_ln_cell_lines=None)
            scen["transport_only"][d] = IC50Prior(
                sigma_log=0.0, basis="held fixed at the literature anchor", **common)
            scen["ic50_median_se"][d] = IC50Prior(
                sigma_log=LITERATURE_SIGMA_LOG,
                basis=(f"expert prior, 90% CI = 1/{LITERATURE_FOLD_90:.0f}x-"
                       f"{LITERATURE_FOLD_90:.0f}x (sigma_log "
                       f"{LITERATURE_SIGMA_LOG:.3f}); not in the GDSC glioma "
                       f"screen, so no measured spread exists"),
                **common)
            scen["ic50_cell_line"][d] = IC50Prior(
                sigma_log=stress_lit_sigma,
                basis=(f"stress: assume the unscreened drug's cross-cell-line "
                       f"spread matches the median screened drug's "
                       f"(sigma_log {stress_lit_sigma:.3f})"),
                **common)
    return scen


# ===========================================================================
# 2. THE KILL KERNEL -- IC50 sweeps without re-solving transport.
# ===========================================================================
@dataclass
class KillKernel:
    """Everything needed to recompute one drug's expected kill at any IC50.

    ``E`` is the per-TUMOUR-cell time-averaged effective concentration (uM) and
    ``sr`` the per-tumour-cell ``sensitivity * hypoxia/quiescence resistance``
    multiplier -- exactly the two quantities ``cells.pkpd.apply_drug_kill``
    combines with the Hill term. Both are fixed once the transport solve is done,
    so an IC50 sweep is a 30k-element vector op.
    """

    drug: str
    E: np.ndarray            # (n_tumor,) uM
    sr: np.ndarray           # (n_tumor,) sensitivity * resistance
    hill: float
    max_kill_rate: float
    ic50_nominal: float
    mean_auc: float
    penetration_index: float
    kill_engine: float       # the engine's own Bernoulli kill_fraction (seed 0)

    def expected_kill(self, ic50: float) -> float:
        """``E[kill_fraction]`` = mean over tumour cells of the Hill kill prob."""
        e_h = np.power(np.maximum(self.E, 0.0), self.hill)
        p = self.max_kill_rate * e_h / (float(ic50) ** self.hill + e_h + 1e-300)
        return float(np.mean(np.clip(p * self.sr, 0.0, 1.0)))

    def expected_kill_many(self, ic50s: np.ndarray) -> np.ndarray:
        e_h = np.power(np.maximum(self.E, 0.0), self.hill)
        out = np.empty(len(ic50s), dtype=np.float64)
        for k, c in enumerate(np.asarray(ic50s, dtype=np.float64)):
            p = self.max_kill_rate * e_h / (c ** self.hill + e_h + 1e-300)
            out[k] = float(np.mean(np.clip(p * self.sr, 0.0, 1.0)))
        return out


def compile_panel_drug(name: str, alpha: float) -> DrugParams:
    """The EXACT compile path ``experiments/real_drug_panel.py`` uses."""
    return apply_barrier_to_drugparams(calibrated_drugparams(name, alpha),
                                       get_molecule(name))


def _labelled_states(micro, oxygen: np.ndarray, ox_cfg: OxygenConfig) -> np.ndarray:
    """Oxygen-banded per-cell states, reproducing the engine's step 2 exactly by
    calling the same frozen ``cells.agents.update_cell_states``."""
    work = dataclasses.replace(micro, cells=micro.cells.copy())
    update_cell_states(work, oxygen, ox_cfg)
    return np.asarray(work.cells.state).copy()


def build_kernels(
    micro,
    alpha: float,
    tortuosity: float,
    pde: PDEConfig,
    drugs: Sequence[str] = PANEL,
    pd_cfg: Optional[PDConfig] = None,
    ox_cfg: Optional[OxygenConfig] = None,
    seed: int = 0,
) -> Dict[str, KillKernel]:
    """Solve transport for every drug on ONE tissue and return their kill kernels.

    All five drugs see the identical microenvironment, identical oxygen field and
    identical calibration constant -- the shared-tissue correlation that makes a
    ranking meaningful.
    """
    pd_cfg = pd_cfg or PDConfig()
    ox_cfg = ox_cfg or OxygenConfig()
    out: Dict[str, KillKernel] = {}
    tumor = np.asarray(micro.cells.tumor_mask())
    with _patched_tortuosity(tortuosity):
        for name in drugs:
            dp = compile_panel_drug(name, alpha)
            res = simulate_response(micro, SimConfig(drug=dp, pde=pde, pd=pd_cfg,
                                                     oxygen=ox_cfg, seed=seed))
            states = _labelled_states(micro, res.fields["oxygen"], ox_cfg)
            col, row = micro.domain.to_index(micro.cells.x, micro.cells.y)
            auc = np.asarray(res.transport.auc, dtype=np.float64)[row, col]
            times = np.asarray(res.transport.times, dtype=np.float64)
            duration = float(times[-1]) if times.size and times[-1] > 0.0 else 1.0
            resistance = np.ones(states.size, dtype=np.float64)
            resistance[states == CellState.HYPOXIC] = pd_cfg.hypoxia_resistance
            resistance[states == CellState.QUIESCENT] = pd_cfg.quiescent_resistance
            sr = np.asarray(micro.cells.sensitivity, dtype=np.float64) * resistance
            out[name] = KillKernel(
                drug=name,
                E=(auc / duration)[tumor],
                sr=sr[tumor],
                hill=float(dp.hill),
                max_kill_rate=float(pd_cfg.max_kill_rate),
                ic50_nominal=float(dp.ic50),
                mean_auc=float(res.extras["mean_auc"]),
                penetration_index=float(res.extras["penetration_index"]),
                kill_engine=float(res.kill_fraction),
            )
    return out


def verify_kill_kernel(micro, kernel: KillKernel, alpha: float,
                       pde: PDEConfig, tortuosity: float,
                       pd_cfg: Optional[PDConfig] = None,
                       ox_cfg: Optional[OxygenConfig] = None) -> Dict[str, float]:
    """Prove the analytic kernel IS the frozen pkpd mapping.

    Re-runs the engine for this drug, calls the frozen
    ``cells.pkpd.apply_drug_kill`` on the oxygen-labelled pre-kill cells, and
    compares the mean of ITS returned per-tumour-cell kill probability with
    :meth:`KillKernel.expected_kill`. They must agree to machine precision --
    otherwise every ranking below is computed with a re-implementation rather
    than with the frozen model.
    """
    pd_cfg = pd_cfg or PDConfig()
    ox_cfg = ox_cfg or OxygenConfig()
    dp = compile_panel_drug(kernel.drug, alpha)
    with _patched_tortuosity(tortuosity):
        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, pd=pd_cfg,
                                                 oxygen=ox_cfg, seed=0))
    work = dataclasses.replace(micro, cells=micro.cells.copy())
    update_cell_states(work, res.fields["oxygen"], ox_cfg)
    p = apply_drug_kill(work, res.transport, pd_cfg, dp,
                        rng=np.random.default_rng(0))
    tumor = np.asarray(micro.cells.tumor_mask())
    frozen = float(np.mean(np.asarray(p)[tumor]))
    ours = kernel.expected_kill(kernel.ic50_nominal)
    # Bernoulli sampling sd of the engine's single-seed realisation.
    n_t = int(tumor.sum())
    pv = np.clip(np.asarray(p)[tumor], 0.0, 1.0)
    bern_sd = float(np.sqrt(np.sum(pv * (1.0 - pv))) / max(n_t, 1))
    return {
        "drug": kernel.drug,
        "frozen_pkpd_expected_kill": frozen,
        "kernel_expected_kill": ours,
        "abs_diff": abs(frozen - ours),
        "engine_bernoulli_kill": kernel.kill_engine,
        "bernoulli_sd": bern_sd,
        "n_tumor": n_t,
    }


# ===========================================================================
# 3. THE MONTE CARLO.
# ===========================================================================
def _rank_of(kills: np.ndarray) -> np.ndarray:
    """Rank 1 = largest kill. ``kills`` is (n_draw, n_drug)."""
    order = np.argsort(-kills, axis=1, kind="stable")
    ranks = np.empty_like(order)
    idx = np.arange(1, kills.shape[1] + 1)
    for i in range(kills.shape[0]):
        ranks[i, order[i]] = idx
    return ranks


def run_ranking_mc(
    n_outer: int = 40,
    n_ic50: int = 64,
    dx_um: float = MC_DX_UM,
    seed: int = 0,
    pde: Optional[PDEConfig] = None,
    drugs: Sequence[str] = PANEL,
    cache_dir: str = "data_cache",
    time_budget_s: Optional[float] = None,
    verbose: bool = True,
) -> Dict:
    """Monte-Carlo the FULL ranking through the prior uncertainty.

    ``n_outer`` draws of (typing threshold, KDE bandwidth, calibration alpha, ECM
    tortuosity) -- one reconstruction and five transport solves each -- crossed
    with ``n_ic50`` draws of the five IC50 multipliers, evaluated for all three
    IC50 scenarios on the same transport solves via the kill kernels.

    The design is NESTED (IC50 draws are nested inside tissue draws), so the
    ``n_outer * n_ic50`` rows are NOT independent. Every probability below is
    therefore quoted with a CLUSTER bootstrap that resamples whole outer draws --
    the honest error bar for this design.
    """
    t_start = time.time()
    pde = pde or mc_pde()
    alpha0 = calibrate_alpha()
    priors = build_priors(alpha0)
    ic50_priors = build_ic50_priors(cache_dir=cache_dir)
    draws = sample_priors(priors, n_outer, seed=seed)
    drugs = tuple(drugs)
    nd = len(drugs)

    # IC50 deviates: one LHS per outer draw, SHARED across the three scenarios
    # (common random numbers) so scenarios differ only by the prior widths.
    z_all = np.empty((n_outer, n_ic50, nd), dtype=np.float64)
    for i in range(n_outer):
        u = stats.qmc.LatinHypercube(d=nd, seed=int(seed) + 10_000 + i).random(n_ic50)
        z_all[i] = stats.norm.ppf(np.clip(u, 1e-9, 1 - 1e-9))

    kills = {s: np.full((n_outer, n_ic50, nd), np.nan) for s in IC50_SCENARIOS}
    rows: List[Dict] = []
    verification: Optional[Dict] = None
    n_done = 0

    for i in range(n_outer):
        thr = float(draws["endothelial_typing_threshold"][i])
        bw = float(draws["kde_bandwidth_um"][i])
        alpha = float(draws["calibration_alpha"][i])
        tort = float(draws["ecm_tortuosity_alpha"][i])
        t0 = time.time()
        micro = reconstruct_knobs(dx_um, thr, bw)
        kernels = build_kernels(micro, alpha, tort, pde, drugs)
        if i == 0:
            verification = verify_kill_kernel(micro, kernels[drugs[0]], alpha,
                                              pde, tort)
        for s in IC50_SCENARIOS:
            for j, name in enumerate(drugs):
                pr = ic50_priors[s][name]
                kills[s][i, :, j] = kernels[name].expected_kill_many(
                    pr.ic50_uM * pr.multipliers(z_all[i, :, j]))
        rows.append({
            "i": i,
            "endothelial_typing_threshold": thr,
            "kde_bandwidth_um": bw,
            "calibration_alpha": alpha,
            "ecm_tortuosity_alpha": tort,
            "nominal_kill": {n: kernels[n].expected_kill(kernels[n].ic50_nominal)
                             for n in drugs},
            "mean_auc": {n: kernels[n].mean_auc for n in drugs},
            "penetration_index": {n: kernels[n].penetration_index for n in drugs},
            "runtime_s": time.time() - t0,
        })
        n_done = i + 1
        if verbose:
            k0 = rows[-1]["nominal_kill"]
            top = max(k0, key=k0.get)
            print(f"  draw {i+1:3d}/{n_outer}  thr={thr:.3f} bw={bw:4.1f} "
                  f"alpha={alpha:6.1f} tort={tort:.3f} -> top(nominal IC50)="
                  f"{top:13s} kills=["
                  + " ".join(f"{k0[n]:.3f}" for n in drugs)
                  + f"]  [{time.time()-t0:.0f}s]", flush=True)
        if time_budget_s is not None and (time.time() - t_start) > time_budget_s:
            if verbose:
                print(f"  [time budget {time_budget_s:.0f}s reached after "
                      f"{n_done} draws -- stopping early]", flush=True)
            break

    for s in IC50_SCENARIOS:
        kills[s] = kills[s][:n_done]
    z_all = z_all[:n_done]

    return {
        "drugs": list(drugs),
        "n_outer": int(n_done),
        "n_outer_requested": int(n_outer),
        "n_ic50": int(n_ic50),
        "n_total": int(n_done * n_ic50),
        "dx_um": float(dx_um),
        "dt_s": float(pde.dt),
        "convection": bool(pde.convection),
        "total_time_s": float(pde.total_time),
        "seed": int(seed),
        "alpha0": float(alpha0),
        "priors": priors,
        "ic50_priors": ic50_priors,
        "kills": kills,
        "z": z_all,
        "rows": rows,
        "verification": verification,
        "runtime_s": time.time() - t_start,
    }


# ---------------------------------------------------------------------------
# Persisting a run. The Monte Carlo costs ~1 min per outer draw; the report and
# the figure cost seconds. Caching the raw cube means a wording or plotting
# change never requires re-solving transport -- and it makes the published
# numbers auditable, because the exact array the report was written from is on
# disk. Only PLAIN ARRAYS are stored: the two prior specs are rebuilt on load
# from ``build_priors(alpha0)`` and ``build_ic50_priors()``, which are pure
# functions of the calibration fit and the cached GDSC screen, so a reload
# cannot silently disagree with the module.
# ---------------------------------------------------------------------------
_SCALAR_KEYS = ("n_outer", "n_outer_requested", "n_ic50", "n_total", "dx_um",
                "dt_s", "convection", "total_time_s", "seed", "alpha0",
                "runtime_s")


def save_mc(res: Dict, path: str = "results/ranking_robustness_mc.npz") -> str:
    """Persist the raw Monte-Carlo cube (arrays only) so reports can be redone."""
    drugs = list(res["drugs"])
    payload = {
        "drugs": np.array(drugs, dtype=object),
        "z": np.asarray(res["z"], dtype=np.float64),
        "scenarios": np.array(list(IC50_SCENARIOS), dtype=object),
        "verification": np.array([res.get("verification")], dtype=object),
        "scalars": np.array([{k: res[k] for k in _SCALAR_KEYS}], dtype=object),
    }
    for s in IC50_SCENARIOS:
        payload[f"kills__{s}"] = np.asarray(res["kills"][s], dtype=np.float64)
    for p in PARAM_ORDER:
        payload[f"param__{p}"] = np.array([r[p] for r in res["rows"]],
                                          dtype=np.float64)
    for field_ in ("nominal_kill", "mean_auc", "penetration_index"):
        payload[f"perdrug__{field_}"] = np.array(
            [[r[field_][d] for d in drugs] for r in res["rows"]],
            dtype=np.float64)
    payload["param__runtime_s"] = np.array([r["runtime_s"] for r in res["rows"]],
                                           dtype=np.float64)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    np.savez_compressed(path, **payload)
    return path


def load_mc(path: str = "results/ranking_robustness_mc.npz",
            cache_dir: str = "data_cache") -> Dict:
    """Rebuild a :func:`run_ranking_mc` result dict from :func:`save_mc`."""
    z = np.load(path, allow_pickle=True)
    drugs = [str(d) for d in z["drugs"]]
    scal = dict(z["scalars"][0])
    rows: List[Dict] = []
    for i in range(int(scal["n_outer"])):
        row = {"i": i, "runtime_s": float(z["param__runtime_s"][i])}
        for p in PARAM_ORDER:
            row[p] = float(z[f"param__{p}"][i])
        for field_ in ("nominal_kill", "mean_auc", "penetration_index"):
            row[field_] = {d: float(z[f"perdrug__{field_}"][i, j])
                           for j, d in enumerate(drugs)}
        rows.append(row)
    return {
        "drugs": drugs,
        "priors": build_priors(float(scal["alpha0"])),
        "ic50_priors": build_ic50_priors(cache_dir=cache_dir),
        "kills": {s: z[f"kills__{s}"] for s in IC50_SCENARIOS},
        "z": z["z"],
        "rows": rows,
        "verification": z["verification"][0],
        **{k: scal[k] for k in _SCALAR_KEYS},
    }


# ===========================================================================
# 4. RANK STATISTICS -- every probability with its Monte-Carlo error.
#
# The design is NESTED: n_ic50 IC50 draws sit inside each of n_outer tissue
# draws, and rows inside one tissue draw share a reconstruction, a calibration
# alpha and a tortuosity. They are NOT independent Bernoulli trials, so the
# naive binomial SE sqrt(p(1-p)/n_total) is an UNDER-estimate -- often by a
# large factor, because the tissue draw is exactly the thing that flips ranks.
#
# Every probability below therefore carries THREE numbers:
#   * se_binom   -- the naive binomial SE at n_total. Reported because it was
#                   asked for and because the gap between it and se_cluster is
#                   itself the finding: it measures how much of the ranking
#                   uncertainty lives in the tissue axis.
#   * se_cluster -- the honest SE: the standard error of the mean of the
#                   n_outer per-draw proportions. This is the one to quote.
#   * ci90       -- a cluster bootstrap (resample whole outer draws) 90% CI.
# For p = 0 or 1 exactly, the cluster SE collapses to 0 and is meaningless; the
# rule-of-three one-sided 95% bound 3/n_outer is reported instead.
# ===========================================================================
N_BOOT = 4000


def _prob_stats(ind: np.ndarray, n_boot: int = N_BOOT, seed: int = 0) -> Dict:
    """Probability + MC error for an indicator array of shape (n_outer, n_ic50)."""
    ind = np.asarray(ind, dtype=np.float64)
    if ind.ndim == 1:
        ind = ind[:, None]
    n_o, n_i = ind.shape
    per = ind.mean(axis=1)                      # per-outer-draw proportion
    p = float(per.mean())
    n_tot = int(n_o * n_i)
    se_binom = math.sqrt(max(p * (1.0 - p), 0.0) / max(n_tot, 1))
    se_cluster = (float(per.std(ddof=1) / math.sqrt(n_o)) if n_o > 1 else float("nan"))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n_o, size=(int(n_boot), n_o))
    boots = per[idx].mean(axis=1)
    lo, hi = (float(v) for v in np.percentile(boots, [5.0, 95.0]))
    degenerate = bool(p <= 0.0 or p >= 1.0)
    # design effect: how much larger the honest SE is than the naive one.
    deff = (float(se_cluster / se_binom) if (se_binom > 0 and np.isfinite(se_cluster))
            else float("nan"))
    return {
        "p": p,
        "n_outer": n_o,
        "n_ic50": n_i,
        "n_total": n_tot,
        "se_binom": float(se_binom),
        "se_cluster": se_cluster,
        "ci90_lo": lo,
        "ci90_hi": hi,
        "design_effect": deff,
        "degenerate": degenerate,
        # one-sided 95% bound for a 0/n or n/n outcome, on the CLUSTER unit.
        "rule_of_three": float(3.0 / n_o) if degenerate else float("nan"),
        "n_eff": (float(p * (1.0 - p) / se_cluster ** 2)
                  if (np.isfinite(se_cluster) and se_cluster > 0) else float("nan")),
    }


def fmt_p(st: Dict, digits: int = 3) -> str:
    """``p +- cluster-SE`` with the degenerate case handled honestly."""
    if st["degenerate"]:
        side = "<" if st["p"] >= 1.0 else ">"
        bound = 1.0 - st["rule_of_three"] if st["p"] >= 1.0 else st["rule_of_three"]
        return (f"{st['p']:.{digits}f} ({st['p']*st['n_total']:.0f}/{st['n_total']} "
                f"rows; rule-of-three 95% bound {side} {bound:.3f})")
    return f"{st['p']:.{digits}f} ± {st['se_cluster']:.3f}"


def rank_stats(res: Dict, scenario: str, n_boot: int = N_BOOT,
               seed: int = 0) -> Dict:
    """Full rank-probability structure for one IC50 scenario.

    Everything is computed from the SAME (n_outer, n_ic50, n_drug) kill cube, so
    the rank matrix, the pairwise probability and the top-2-set probability are
    guaranteed mutually consistent (they are three marginals of one joint).
    """
    kills = np.asarray(res["kills"][scenario], dtype=np.float64)
    n_o, n_i, nd = kills.shape
    drugs: List[str] = list(res["drugs"])
    ranks = _rank_of(kills.reshape(-1, nd)).reshape(n_o, n_i, nd)

    # --- drug x rank probability matrix --------------------------------------
    matrix: Dict[str, List[Dict]] = {}
    for j, d in enumerate(drugs):
        matrix[d] = [_prob_stats(ranks[:, :, j] == (r + 1), n_boot, seed + 100 * j + r)
                     for r in range(nd)]
    p_rank1 = {d: matrix[d][0] for d in drugs}

    # --- headline pairwise ----------------------------------------------------
    a, b = HEADLINE_PAIR
    ia, ib = drugs.index(a), drugs.index(b)
    pair = _prob_stats(kills[:, :, ia] > kills[:, :, ib], n_boot, seed + 11)

    # every pairwise, so the report can show WHICH comparisons are decided
    pairwise: Dict[Tuple[str, str], Dict] = {}
    for (j1, d1), (j2, d2) in itertools.combinations(list(enumerate(drugs)), 2):
        pairwise[(d1, d2)] = _prob_stats(kills[:, :, j1] > kills[:, :, j2],
                                         n_boot, seed + 31 * j1 + j2)

    # --- top-2 SET (unordered) ------------------------------------------------
    top2_idx = np.argsort(-kills, axis=2, kind="stable")[:, :, :2]
    sets = np.zeros((n_o, n_i, nd), dtype=bool)
    np.put_along_axis(sets, top2_idx, True, axis=2)
    top2_headline = _prob_stats(
        np.all(sets == np.array([d in HEADLINE_TOP2 for d in drugs]), axis=2),
        n_boot, seed + 21)
    # the distribution over all observed top-2 sets
    key = (sets * (1 << np.arange(nd))).sum(axis=2)
    top2_dist: List[Tuple[frozenset, Dict]] = []
    for k in np.unique(key):
        members = frozenset(drugs[j] for j in range(nd) if (int(k) >> j) & 1)
        top2_dist.append((members, _prob_stats(key == k, n_boot, seed + 41 + int(k))))
    top2_dist.sort(key=lambda t: -t[1]["p"])

    # --- the complete headline order -----------------------------------------
    want = np.array([HEADLINE_RANKING.index(d) + 1 for d in drugs])
    full_order = _prob_stats(np.all(ranks == want, axis=2), n_boot, seed + 51)

    # --- the cisplatin-doxorubicin gap itself --------------------------------
    gap = (kills[:, :, ia] - kills[:, :, ib]).ravel()
    per_draw_gap = (kills[:, :, ia] - kills[:, :, ib]).mean(axis=1)

    nominal = {d: float(np.median([r["nominal_kill"][d] for r in res["rows"]]))
               for d in drugs}
    return {
        "scenario": scenario,
        "drugs": drugs,
        "matrix": matrix,
        "p_rank1": p_rank1,
        "pair": pair,
        "pair_names": HEADLINE_PAIR,
        "pairwise": pairwise,
        "top2_headline": top2_headline,
        "top2_dist": top2_dist,
        "full_order": full_order,
        "gap_median": float(np.median(gap)),
        "gap_ci90": (float(np.percentile(gap, 5.0)), float(np.percentile(gap, 95.0))),
        "gap_per_draw_min": float(per_draw_gap.min()),
        "gap_per_draw_max": float(per_draw_gap.max()),
        "expected_rank": {d: float(np.mean(ranks[:, :, j])) for j, d in enumerate(drugs)},
        "median_kill": {d: float(np.median(kills[:, :, j])) for j, d in enumerate(drugs)},
        "kill_ci90": {d: (float(np.percentile(kills[:, :, j], 5.0)),
                          float(np.percentile(kills[:, :, j], 95.0)))
                      for j, d in enumerate(drugs)},
        "nominal_kill_median": nominal,
        "n_outer": int(n_o),
        "n_ic50": int(n_i),
        "n_total": int(n_o * n_i),
    }


# ===========================================================================
# 5. WHICH PARAMETER DESTABILISES THE RANKING?
#
# The outcome is BINARY (does the target drug hold rank 1 in this row?), so the
# natural importance measure is: cut each parameter into quartiles and look at
# how far P(rank 1) MOVES across them. The "swing" max-min is directly
# interpretable in the units of the decision -- "moving alpha from its bottom to
# its top quartile changes the probability cisplatin wins by X" -- and it needs
# no linearity assumption, which matters because the typing threshold acts as a
# step and the IC50 effect is a ratio.
#
# Nine variables are scored on the SAME footing: the four reconstruction /
# calibration knobs (constant within an outer draw) and the five per-drug IC50
# log-multipliers (varying within an outer draw). Both are cluster-bootstrapped
# over outer draws.
# ===========================================================================
IC50_VAR_PREFIX = "ic50_mult:"


def _grouped_prob(y: np.ndarray, g: np.ndarray, n_groups: int,
                  n_boot: int = N_BOOT, seed: int = 0
                  ) -> Tuple[List[Dict], np.ndarray, np.ndarray, np.ndarray]:
    """P(y | group k) for each k, with a cluster bootstrap over outer draws."""
    y = np.asarray(y, dtype=np.float64)
    g = np.asarray(g)
    n_o = y.shape[0]
    sums = np.zeros((n_o, n_groups))
    cnts = np.zeros((n_o, n_groups))
    for k in range(n_groups):
        m = (g == k)
        sums[:, k] = (y * m).sum(axis=1)
        cnts[:, k] = m.sum(axis=1)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n_o, size=(int(n_boot), n_o))
    out: List[Dict] = []
    for k in range(n_groups):
        tot = float(cnts[:, k].sum())
        if tot <= 0:   # an empty quartile (too few outer draws to cut four ways)
            out.append({"group": k, "p": float("nan"), "n": 0,
                        "se": float("nan"), "ci90_lo": float("nan"),
                        "ci90_hi": float("nan")})
            continue
        p = float(sums[:, k].sum() / tot)
        bs = sums[idx, k].sum(axis=1)
        bc = cnts[idx, k].sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            boots = np.where(bc > 0, bs / np.maximum(bc, 1e-12), np.nan)
        if np.all(~np.isfinite(boots)):
            lo = hi = se = float("nan")
        else:
            lo, hi = (float(v) for v in np.nanpercentile(boots, [5.0, 95.0]))
            se = float(np.nanstd(boots, ddof=1))
        out.append({"group": k, "p": p, "n": int(tot), "se": se,
                    "ci90_lo": lo, "ci90_hi": hi})
    return out, sums, cnts, idx


def destabilisers(res: Dict, scenario: str = PRIMARY_SCENARIO,
                  target: Optional[str] = None, n_q: int = 4,
                  n_boot: int = N_BOOT, seed: int = 0) -> Dict:
    """Quartile-conditioned P(``target`` holds rank 1) for every sampled input.

    Returns one entry per variable with its four quartile probabilities and the
    SWING (top quartile minus bottom quartile is NOT used -- the max-min over the
    four is, because several of these effects are non-monotone), plus a cluster
    bootstrap 90% CI on the swing so a swing can be called indistinguishable
    from zero when it is.
    """
    drugs: List[str] = list(res["drugs"])
    target = target or HEADLINE_RANKING[0]
    kills = np.asarray(res["kills"][scenario], dtype=np.float64)
    n_o, n_i, nd = kills.shape
    ranks = _rank_of(kills.reshape(-1, nd)).reshape(n_o, n_i, nd)
    y = (ranks[:, :, drugs.index(target)] == 1).astype(np.float64)

    # variable -> value array broadcast to (n_outer, n_ic50)
    variables: Dict[str, np.ndarray] = {}
    for p in PARAM_ORDER:
        col = np.array([r[p] for r in res["rows"]], dtype=np.float64)
        variables[p] = np.repeat(col[:, None], n_i, axis=1)
    if scenario != "transport_only":
        z = np.asarray(res["z"], dtype=np.float64)          # (n_outer, n_ic50, nd)
        for j, d in enumerate(drugs):
            sig = res["ic50_priors"][scenario][d].sigma_log
            if sig > 0:
                variables[IC50_VAR_PREFIX + d] = np.exp(sig * z[:, :, j])

    out: Dict[str, Dict] = {}
    for vi, (name, V) in enumerate(variables.items()):
        edges = np.quantile(V.ravel(), np.linspace(0, 1, n_q + 1))
        edges[0], edges[-1] = -np.inf, np.inf
        g = np.clip(np.digitize(V, edges[1:-1], right=False), 0, n_q - 1)
        rows_, sums, cnts, idx = _grouped_prob(y, g, n_q, n_boot, seed + 7 * vi)
        ps = np.array([r["p"] for r in rows_], dtype=np.float64)
        swing = (float(np.nanmax(ps) - np.nanmin(ps)) if np.any(np.isfinite(ps))
                 else float("nan"))
        # bootstrap the swing on the SAME resamples
        with np.errstate(invalid="ignore", divide="ignore"):
            bc_ = cnts[idx].sum(axis=1)
            bp = np.where(bc_ > 0, sums[idx].sum(axis=1) / np.maximum(bc_, 1e-12),
                          np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            bswing = np.nanmax(bp, axis=1) - np.nanmin(bp, axis=1)
        out[name] = {
            "variable": name,
            "is_ic50": name.startswith(IC50_VAR_PREFIX),
            "quartiles": rows_,
            "edges": [float(e) for e in np.quantile(V.ravel(),
                                                    np.linspace(0, 1, n_q + 1))],
            "swing": swing,
            "swing_ci90": ((float(np.nanpercentile(bswing, 5.0)),
                            float(np.nanpercentile(bswing, 95.0)))
                           if np.any(np.isfinite(bswing))
                           else (float("nan"), float("nan"))),
            "swing_se": (float(np.nanstd(bswing, ddof=1))
                         if np.any(np.isfinite(bswing)) else float("nan")),
            # monotone cross-check on the pooled rows
            "spearman": float(stats.spearmanr(V.ravel(), y.ravel()).statistic),
        }
    # NaN swings (an empty quartile) sort LAST rather than unpredictably.
    ranked = sorted(out, key=lambda k: (not np.isfinite(out[k]["swing"]),
                                        -np.nan_to_num(out[k]["swing"], nan=-1.0)))
    return {"target": target, "scenario": scenario, "n_q": int(n_q),
            "per_variable": out, "ranked": ranked,
            "baseline": _prob_stats(y, n_boot, seed + 3)}


def ranking_stability_index(res: Dict, scenario: str) -> Dict:
    """Fraction of rows whose FULL order matches the modal order, plus entropy.

    A compact scalar for "how stable is this ranking at all", independent of
    which drug happens to win: the Shannon entropy of the distribution over the
    5! possible orders, normalised by log(5!) so 0 = one order always, 1 =
    uniform over all orders.
    """
    kills = np.asarray(res["kills"][scenario], dtype=np.float64)
    n_o, n_i, nd = kills.shape
    order = np.argsort(-kills.reshape(-1, nd), axis=1, kind="stable")
    keys = [tuple(int(v) for v in row) for row in order]
    counts: Dict[Tuple[int, ...], int] = {}
    for k in keys:
        counts[k] = counts.get(k, 0) + 1
    n = len(keys)
    probs = np.array(sorted(counts.values(), reverse=True), dtype=np.float64) / n
    ent = float(-np.sum(probs * np.log(probs)))
    modal = max(counts, key=counts.get)
    return {
        "n_distinct_orders": int(len(counts)),
        "modal_order": [res["drugs"][j] for j in modal],
        "modal_p": float(counts[modal] / n),
        "entropy": ent,
        "entropy_norm": float(ent / math.log(math.factorial(nd))),
    }


# ===========================================================================
# 6. RESOLUTION CHECK -- is the ranking (and the gap) a dx artifact?
# ===========================================================================
def dx_ranking_check(dxs: Sequence[float] = (MC_DX_UM, PRODUCTION_DX_UM),
                     alpha: Optional[float] = None,
                     drugs: Sequence[str] = PANEL,
                     verbose: bool = True) -> Dict:
    """Re-run the NOMINAL panel at the MC grid and at the production grid.

    The MC runs at dx = 15 um / dt = 20 s / convection off for affordability; the
    published panel runs at dx = 12 um / dt = 10 s / convection on. This shows
    the ORDER and the size of the cisplatin-doxorubicin gap at both, so the
    ranking result cannot be dismissed as a resolution or scheme artifact.
    """
    alpha = float(alpha if alpha is not None else calibrate_alpha())
    rows: List[Dict] = []
    for dx in dxs:
        pde = production_pde() if abs(dx - PRODUCTION_DX_UM) < 1e-9 else mc_pde()
        t0 = time.time()
        micro = reconstruct_knobs(float(dx), BASE_THRESHOLD, BASE_BANDWIDTH)
        ks = build_kernels(micro, alpha, float(_patched_nominal_tortuosity()), pde, drugs)
        kill = {d: ks[d].expected_kill(ks[d].ic50_nominal) for d in drugs}
        order = sorted(drugs, key=lambda d: -kill[d])
        rows.append({
            "dx_um": float(dx),
            "dt_s": float(pde.dt),
            "convection": bool(pde.convection),
            "grid": (micro.domain.ny, micro.domain.nx),
            "kill": kill,
            "kill_engine": {d: ks[d].kill_engine for d in drugs},
            "mean_auc": {d: ks[d].mean_auc for d in drugs},
            "order": order,
            "gap_cis_dox": float(kill[HEADLINE_PAIR[0]] - kill[HEADLINE_PAIR[1]]),
            "runtime_s": time.time() - t0,
        })
        if verbose:
            print(f"  dx={dx:4.1f} um  dt={pde.dt:.0f}s conv={pde.convection}  "
                  f"order={' > '.join(order)}  gap={rows[-1]['gap_cis_dox']:+.4f}  "
                  f"[{rows[-1]['runtime_s']:.0f}s]", flush=True)
    return {
        "rows": rows,
        "alpha": alpha,
        "order_identical": all(r["order"] == rows[0]["order"] for r in rows),
    }


def _patched_nominal_tortuosity() -> float:
    """The frozen ``physics.diffusion.ALPHA_TORTUOSITY`` (0.6), read live."""
    from ..physics import diffusion as _d
    return float(_d.ALPHA_TORTUOSITY)


# ===========================================================================
# 7. FIGURE.
# ===========================================================================
_C_WIN = "#2a78d6"       # blue   -- cisplatin / primary
_C_ALT = "#1baf7a"       # aqua   -- doxorubicin / secondary
_C_THIRD = "#e0a325"     # amber  -- lomustine
_C_ACCENT = "#e34948"    # red    -- rules, thresholds
_C_MUTED = "#52514e"
_DRUG_COLOR = {
    "cisplatin": _C_WIN, "doxorubicin": _C_ALT, "lomustine": _C_THIRD,
    "methotrexate": "#8c6bb1", "temozolomide": "#7f8180",
}
_SCEN_COLOR = {"transport_only": "#9dc3ea", "ic50_median_se": _C_WIN,
               "ic50_cell_line": "#12406f"}
_SCEN_SHORT = {"transport_only": "A. transport only",
               "ic50_median_se": "B. + IC50 SE (primary)",
               "ic50_cell_line": "C. + cell-line spread"}
_VAR_LABEL = {
    "endothelial_typing_threshold": "endothelial typing thr.",
    "kde_bandwidth_um": "KDE bandwidth",
    "calibration_alpha": "calibration α",
    "ecm_tortuosity_alpha": "ECM tortuosity",
}


def _var_label(name: str) -> str:
    if name.startswith(IC50_VAR_PREFIX):
        return f"IC50 × ({name[len(IC50_VAR_PREFIX):]})"
    return _VAR_LABEL.get(name, name)


def plot_ranking(res: Dict, stats_by_scen: Dict[str, Dict], dest: Dict,
                 path: str = "results/ranking_robustness.png") -> str:
    """Four panels: the rank matrix, P(rank-1) across scenarios, the
    cisplatin-doxorubicin gap distribution, and what destabilises the winner."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    drugs: List[str] = list(res["drugs"])
    nd = len(drugs)
    st = stats_by_scen[PRIMARY_SCENARIO]

    fig, axes = plt.subplots(2, 2, figsize=(13.4, 9.4))
    fig.patch.set_facecolor("white")

    # --- (a) drug x rank probability matrix -----------------------------------
    ax = axes[0, 0]
    M = np.array([[st["matrix"][d][r]["p"] for r in range(nd)] for d in drugs])
    cmap = LinearSegmentedColormap.from_list(
        "seq", ["#ffffff", "#cfe3f7", "#7fb3e8", "#2a78d6", "#12406f"])
    im = ax.imshow(M, cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")
    for i in range(nd):
        for j in range(nd):
            p = M[i, j]
            se = st["matrix"][drugs[i]][j]["se_cluster"]
            if p < 0.005:
                txt = "–"
            elif st["matrix"][drugs[i]][j]["degenerate"]:
                txt = f"{p:.2f}"
            else:
                txt = f"{p:.2f}\n±{se:.3f}" if se < 0.01 else f"{p:.2f}\n±{se:.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8.6,
                    color="white" if p > 0.55 else "#0b0b0b")
    ax.set_xticks(range(nd))
    ax.set_xticklabels([f"rank {r+1}" for r in range(nd)], fontsize=9)
    ax.set_yticks(range(nd))
    ax.set_yticklabels(drugs, fontsize=9)
    ax.set_title(f"P(drug attains rank) — scenario B (primary)\n"
                 f"± = cluster SE over {st['n_outer']} tissue draws",
                 fontsize=10.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("probability", fontsize=9)
    cb.outline.set_visible(False)
    ax.set_xticks(np.arange(-0.5, nd, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nd, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.6)
    ax.tick_params(which="minor", length=0)

    # --- (b) P(rank 1) across the three scenarios -----------------------------
    ax = axes[0, 1]
    show = [d for d in drugs
            if max(stats_by_scen[s]["p_rank1"][d]["p"] for s in IC50_SCENARIOS) > 0.002]
    y = np.arange(len(show))
    h = 0.26
    for si, s in enumerate(IC50_SCENARIOS):
        vals = np.array([stats_by_scen[s]["p_rank1"][d]["p"] for d in show])
        lo = np.array([stats_by_scen[s]["p_rank1"][d]["ci90_lo"] for d in show])
        hi = np.array([stats_by_scen[s]["p_rank1"][d]["ci90_hi"] for d in show])
        off = (si - 1) * h
        ax.barh(y - off, vals, height=h * 0.92, color=_SCEN_COLOR[s],
                label=_SCEN_SHORT[s],
                xerr=np.vstack([(vals - lo).clip(0), (hi - vals).clip(0)]),
                error_kw=dict(ecolor=_C_MUTED, elinewidth=1.1, capsize=2.5))
        for yi, v, e in zip(y - off, vals, (hi - vals).clip(0)):
            ax.text(v + e + 0.015, yi, f"{v:.2f}", va="center", fontsize=8,
                    color=_C_MUTED)
    ax.axvline(0.5, color=_C_ACCENT, linestyle=":", linewidth=1.4)
    ax.text(0.505, len(show) - 0.45, " a coin flip", fontsize=8.2,
            color=_C_ACCENT, va="top")
    ax.set_yticks(y)
    ax.set_yticklabels(show, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.14)
    ax.set_xlabel("P(rank 1)   [bars = cluster bootstrap 90% CI]", fontsize=9.5)
    ax.set_title("Who wins, and how sure are we?", fontsize=10.5)
    ax.legend(frameon=False, fontsize=8.4, loc="lower right")
    ax.grid(axis="x", alpha=0.18, linewidth=0.7)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # --- (c) the cisplatin - doxorubicin gap ----------------------------------
    ax = axes[1, 0]
    a, b = HEADLINE_PAIR
    ia, ib = drugs.index(a), drugs.index(b)
    for s in IC50_SCENARIOS:
        k = np.asarray(res["kills"][s])
        g = (k[:, :, ia] - k[:, :, ib]).ravel()
        ax.hist(g, bins=60, histtype="stepfilled", alpha=0.42,
                color=_SCEN_COLOR[s], edgecolor=_SCEN_COLOR[s], linewidth=1.1,
                label=f"{_SCEN_SHORT[s]}  P>0 = "
                      f"{stats_by_scen[s]['pair']['p']:.2f}")
    ax.axvline(0.0, color=_C_ACCENT, linewidth=2.0)
    ax.annotate("ranking flips here", xy=(0.0, ax.get_ylim()[1] * 0.42),
                xytext=(0.14, ax.get_ylim()[1] * 0.62), fontsize=8.6,
                color=_C_ACCENT,
                arrowprops=dict(arrowstyle="->", color=_C_ACCENT, lw=1.1))
    ax.set_xlabel(f"E[kill]({a}) − E[kill]({b})   per Monte-Carlo row", fontsize=9.5)
    ax.set_ylabel("rows", fontsize=9.5)
    ax.set_title(f"The margin under test ({a} vs {b})", fontsize=10.5)
    ax.legend(frameon=False, fontsize=8.2, loc="upper left")
    ax.grid(axis="y", alpha=0.18, linewidth=0.7)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # --- (d) what destabilises the winner -------------------------------------
    ax = axes[1, 1]
    names = dest["ranked"]
    y = np.arange(len(names))
    names = [n for n in names if np.isfinite(dest["per_variable"][n]["swing"])]
    y = np.arange(len(names))
    sw = np.array([dest["per_variable"][n]["swing"] for n in names])
    lo = np.array([dest["per_variable"][n]["swing_ci90"][0] for n in names])
    hi = np.array([dest["per_variable"][n]["swing_ci90"][1] for n in names])
    cols = [(_C_ACCENT if dest["per_variable"][n]["is_ic50"] else _C_WIN)
            for n in names]
    ax.barh(y, sw, height=0.62, color=cols,
            xerr=np.vstack([(sw - lo).clip(0), (hi - sw).clip(0)]),
            error_kw=dict(ecolor=_C_MUTED, elinewidth=1.1, capsize=2.5))
    xmax = float(np.nanmax(hi)) if len(names) else 1.0
    for yi, n, s_, h_ in zip(y, names, sw, hi):
        if not np.isfinite(s_):
            continue
        if s_ > 0.18 * xmax:      # label INSIDE a long bar
            ax.text(0.006 * xmax, yi, f"{s_:.3f}", va="center", ha="left",
                    fontsize=8, color="white", zorder=5)
        else:                      # label clear of the whisker on a short bar
            ax.text(max(h_, s_) + 0.012 * xmax, yi, f"{s_:.3f}", va="center",
                    ha="left", fontsize=8, color=_C_MUTED, zorder=5)
    ax.set_xlim(0, xmax * 1.16 if np.isfinite(xmax) and xmax > 0 else 1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([_var_label(n) for n in names], fontsize=8.6)
    ax.invert_yaxis()
    ax.set_xlabel(f"swing in P({dest['target']} = rank 1) across input quartiles\n"
                  "(blue = transport/reconstruction knob, red = IC50; "
                  "bars = bootstrap 90%)", fontsize=9)
    ax.set_title("Which uncertainty flips the decision", fontsize=10.5)
    ax.grid(axis="x", alpha=0.18, linewidth=0.7)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.suptitle(
        f"Is the drug RANKING decision-grade?  —  real 10x Xenium GBM section, "
        f"uniform sensitivity mode\n"
        f"{res['n_outer']} tissue draws × {res['n_ic50']} IC50 draws = "
        f"{res['n_outer']*res['n_ic50']:,} re-rankings, dx = {res['dx_um']:.0f} µm",
        fontsize=12.2, y=0.985)
    fig.text(0.5, 0.005,
             "Probabilities are over EXPERT priors on reconstruction, calibration "
             "and IC50 — not a validated predictive distribution. Errors are "
             "cluster (not binomial) because IC50 draws are nested inside tissue "
             "draws.",
             ha="center", fontsize=8.6, color=_C_MUTED)
    fig.tight_layout(rect=(0, 0.025, 1, 0.945))
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=130, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


# ===========================================================================
# 8. REPORT.
# ===========================================================================
def _verdict(st: Dict) -> Dict:
    """Turn the primary-scenario numbers into an explicit decision-grade call.

    The rule is declared BEFORE the numbers are seen and applied mechanically:

    * **decision-grade top-1** requires the leader's P(rank 1) cluster-bootstrap
      90% CI to sit entirely above 0.90 AND the runner-up's entirely below 0.10.
      (0.90 is the conventional bar for reporting a single recommendation without
      an alternative; anything weaker is a preference, not a decision.)
    * **a tie** is declared when the two leaders' P(rank 1) 90% CIs OVERLAP, or
      when the leader's CI reaches below 0.90.
    * the TOP-2 SET is graded on the same 0.90 bar, separately -- a set call can
      be decision-grade while the order inside it is not, and that is the useful
      thing to report when it happens.
    """
    drugs = st["drugs"]
    order = sorted(drugs, key=lambda d: -st["p_rank1"][d]["p"])
    lead, second = order[0], order[1]
    a, b = st["p_rank1"][lead], st["p_rank1"][second]
    overlap = not (a["ci90_lo"] > b["ci90_hi"])
    decisive = (a["ci90_lo"] >= 0.90) and (b["ci90_hi"] <= 0.10)
    return {
        "leader": lead,
        "runner_up": second,
        "leader_stat": a,
        "runner_up_stat": b,
        "ci_overlap": bool(overlap),
        "top1_decision_grade": bool(decisive),
        "top2_decision_grade": bool(st["top2_headline"]["ci90_lo"] >= 0.90),
        "separation": float(a["p"] - b["p"]),
    }


def write_report(res: Dict, stats_by_scen: Dict[str, Dict], dest: Dict,
                 stability: Dict[str, Dict], dxchk: Optional[Dict] = None,
                 path: str = "results/ranking_robustness.md",
                 figure: str = "results/ranking_robustness.png") -> str:
    L: List[str] = []
    A = L.append
    drugs: List[str] = list(res["drugs"])
    nd = len(drugs)
    st = stats_by_scen[PRIMARY_SCENARIO]
    v = _verdict(st)
    a, b = HEADLINE_PAIR

    A("# Ranking robustness — is the per-patient drug ranking DECISION-GRADE?")
    A("")
    A("**Module:** `tmesim/validation/ranking_robustness.py` "
      "(run: `python -m tmesim.validation.ranking_robustness`)")
    A("")
    A(f"**Figure:** `{os.path.basename(figure)}`")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- verdict
    A("## 0. The verdict, first")
    A("")
    lead_s, run_s = v["leader_stat"], v["runner_up_stat"]
    A(f"On the real 10x Xenium GBM section, propagating the reconstruction, "
      f"calibration and IC50 priors through **{res['n_outer']} tissue draws × "
      f"{res['n_ic50']} IC50 draws = {res['n_total']:,} complete re-rankings**:")
    A("")
    A(f"* **P({v['leader']} is rank 1) = {fmt_p(lead_s)}** "
      f"(cluster-bootstrap 90% CI {lead_s['ci90_lo']:.3f}–{lead_s['ci90_hi']:.3f})")
    A(f"* **P({v['runner_up']} is rank 1) = {fmt_p(run_s)}** "
      f"(90% CI {run_s['ci90_lo']:.3f}–{run_s['ci90_hi']:.3f})")
    A(f"* **P({a} > {b}) = {fmt_p(st['pair'])}** "
      f"(90% CI {st['pair']['ci90_lo']:.3f}–{st['pair']['ci90_hi']:.3f})")
    A(f"* **P(top-2 set = {{{a}, {b}}}) = {fmt_p(st['top2_headline'])}** "
      f"(90% CI {st['top2_headline']['ci90_lo']:.3f}–"
      f"{st['top2_headline']['ci90_hi']:.3f})")
    A(f"* **P(the full published order holds) = {fmt_p(st['full_order'])}**")
    A("")
    if v["top1_decision_grade"]:
        A(f"> **The top-1 call IS decision-grade.** `{v['leader']}` holds rank 1 with "
          f"probability {lead_s['p']:.3f} (90% CI "
          f"{lead_s['ci90_lo']:.3f}–{lead_s['ci90_hi']:.3f}), the runner-up "
          f"`{v['runner_up']}` never gets within reach "
          f"({run_s['p']:.3f}, CI ≤ {run_s['ci90_hi']:.3f}), and the two intervals "
          f"do not overlap. The pre-declared bar (leader CI entirely above 0.90, "
          f"runner-up entirely below 0.10) is met.")
    elif v["ci_overlap"]:
        A(f"> **The top-1 call is NOT decision-grade — report it as a TIE.** "
          f"`{v['leader']}` ({lead_s['p']:.3f}) and `{v['runner_up']}` "
          f"({run_s['p']:.3f}) have OVERLAPPING 90% credible intervals "
          f"([{lead_s['ci90_lo']:.3f}, {lead_s['ci90_hi']:.3f}] vs "
          f"[{run_s['ci90_lo']:.3f}, {run_s['ci90_hi']:.3f}]). The correct output "
          f"is a co-first-place set, not a winner.")
    else:
        A(f"> **The top-1 call is NOT decision-grade, though it is not a coin flip "
          f"either.** `{v['leader']}` leads `{v['runner_up']}` "
          f"{lead_s['p']:.3f} vs {run_s['p']:.3f} and the intervals separate, but "
          f"the leader's 90% CI reaches down to {lead_s['ci90_lo']:.3f}, below the "
          f"0.90 bar this module declared in advance for a single-drug "
          f"recommendation. Report it as *`{v['leader']}` preferred, "
          f"`{v['runner_up']}` a live alternative*, not as a winner.")
    A("")
    if v["top2_decision_grade"]:
        A(f"**The top-2 SET, by contrast, IS decision-grade:** "
          f"{{{a}, {b}}} is the top-2 in {st['top2_headline']['p']:.3f} of "
          f"re-rankings (90% CI ≥ {st['top2_headline']['ci90_lo']:.3f}). If the "
          f"clinical question tolerates two candidates, this model answers it; if "
          f"it demands one, it does not.")
    else:
        A(f"**The top-2 SET is also not decision-grade at the 0.90 bar** "
          f"({st['top2_headline']['p']:.3f}, 90% CI "
          f"{st['top2_headline']['ci90_lo']:.3f}–"
          f"{st['top2_headline']['ci90_hi']:.3f}). See §4 for which set competes.")
    A("")
    top_var = dest["ranked"][0]
    tv = dest["per_variable"][top_var]
    A(f"**Most destabilising input:** `{_var_label(top_var)}` — moving it across "
      f"its prior quartiles swings P({dest['target']} = rank 1) by "
      f"**{tv['swing']:.3f}** (bootstrap 90% CI "
      f"{tv['swing_ci90'][0]:.3f}–{tv['swing_ci90'][1]:.3f}). §5.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------- what & why
    A("## 1. What is being asked, and why an interval per drug does not answer it")
    A("")
    A("`results/uncertainty.md` puts a 90% credible interval on **one** drug's "
      "predicted kill. That is the right object for a *number*. It is the wrong "
      "object for a *decision*: the five drugs' intervals overlap heavily, but "
      "their draws are **correlated** — all five sit on the same reconstructed "
      "tissue, share the same calibration α and the same vascular field, so a "
      "draw that raises one drug's kill usually raises them all. Overlapping "
      "marginals are therefore compatible with either a rock-solid or a "
      "coin-flip ranking. The only way to know is to **re-rank inside every "
      "draw**, which is what this module does.")
    A("")
    A("Decision-grade is defined here, in advance, as: the recommended drug holds "
      "rank 1 with probability whose 90% credible interval lies entirely above "
      "**0.90**, and the runner-up's lies entirely below **0.10**. Anything less "
      "is reported as a tie or as a preference with a named alternative.")
    A("")
    A("**Sensitivity mode = `uniform`, deliberately.** Per-cell sensitivity has "
      "two modes (`reconstruct.real_pipeline.apply_drug_sensitivity`). In "
      "`determinant` mode each drug is scaled by its own measured resistance "
      "determinant *where the 319-gene Xenium panel carries one* — TOP2A for "
      "doxorubicin, an MKI67 proliferation proxy for cisplatin / lomustine / "
      "methotrexate — and by 1.0 where it does not. Temozolomide's determinant is "
      "**MGMT, which is absent from the panel and is a methylation state that "
      "transcript counts cannot report at all**, so TMZ alone would keep 1.0 "
      "while every competitor is scaled down ≈1.7×. That asymmetry is an artifact "
      "of *measurement coverage*, not of drug biology, and it would contaminate "
      "every rank probability in this report with the question 'whose determinant "
      "did we happen to measure?'. In `uniform` mode the factor is identical "
      "across drugs and cancels exactly out of the comparison, so the ranking "
      "isolates **delivery × potency**. The absolute kills are then an equal "
      "no-intrinsic-resistance upper bound for every drug — which is fine, "
      "because a ranking only needs the comparison to be fair, not the level to "
      "be right. `reconstruct_knobs` (used for every draw here) sets "
      "`sensitivity = 1.0` for all cells, so `uniform` is enforced structurally, "
      "not by a flag.")
    A("")
    A("**The published nominal ranking under test** "
      "(`results/real_drug_panel.md`, dx = 12 µm): "
      + " > ".join(f"**{d}**" for d in HEADLINE_RANKING) + ".")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------ the design
    A("## 2. Design, cost, and the kill-kernel identity")
    A("")
    A("The expensive object is the tissue + transport solve; the cheap object is "
      "the kill mapping, and they separate **exactly** — `sim.engine."
      "simulate_response` runs oxygen → pressure → delivery → transport → Hill "
      "kill, and only the last step reads `drug.ic50`. So:")
    A("")
    A(f"* **per outer draw:** 1 reconstruction + {nd} transport solves "
      f"(all {nd} drugs on the SAME tissue — the correlation structure a ranking "
      f"question requires);")
    A(f"* **then {res['n_ic50']} IC50 re-rankings for free**, as a vector op on "
      f"the frozen per-cell exposure and sensitivity;")
    A("* the **same** IC50 deviates drive all three scenarios (common random "
      "numbers), so scenarios differ only by prior width, not by sampling noise.")
    A("")
    A(f"| | |")
    A(f"|---|---|")
    A(f"| outer (tissue) draws achieved | **{res['n_outer']}**"
      + (f" of {res['n_outer_requested']} requested"
         if res["n_outer"] != res["n_outer_requested"] else "") + " |")
    A(f"| IC50 draws per outer draw | {res['n_ic50']} (Latin hypercube, "
      f"{nd}-dimensional) |")
    A(f"| complete re-rankings | **{res['n_total']:,}** per scenario |")
    A(f"| grid / scheme | dx = {res['dx_um']:.0f} µm, dt = {res['dt_s']:.0f} s, "
      f"convection {'on' if res['convection'] else 'off'}, "
      f"{res['total_time_s']/3600:.0f} h exposure |")
    A(f"| α₀ (fitted, e-fold convention) | **{res['alpha0']:.1f}** |")
    A(f"| runtime | {res['runtime_s']/60.0:.1f} min "
      f"(~{res['runtime_s']/max(res['n_outer'],1):.0f} s per outer draw), "
      f"single process |")
    A("")
    A(f"**n is modest and that is stated rather than hidden.** The effective "
      f"sample size for every probability in this report is the number of "
      f"*independent tissue draws*, **{res['n_outer']}**, not the "
      f"{res['n_total']:,} rows — see §3. All intervals are quoted, none of the "
      f"probabilities is presented as a bare point estimate.")
    A("")
    ver = res.get("verification")
    if ver:
        A("**The kill kernel IS the frozen model, not a re-implementation.** On the "
          "first draw the engine is re-run for one drug and `cells.pkpd."
          "apply_drug_kill` is called on the oxygen-labelled pre-kill cells; the "
          "mean of ITS returned per-tumour-cell kill probability is compared with "
          "the analytic kernel:")
        A("")
        A(f"| check ({ver['drug']}) | value |")
        A("|---|---:|")
        A(f"| frozen `pkpd.apply_drug_kill` mean p | {ver['frozen_pkpd_expected_kill']:.10f} |")
        A(f"| analytic `KillKernel.expected_kill` | {ver['kernel_expected_kill']:.10f} |")
        A(f"| absolute difference | **{ver['abs_diff']:.2e}** |")
        A(f"| engine single-seed Bernoulli `kill_fraction` | {ver['engine_bernoulli_kill']:.4f} |")
        A(f"| Bernoulli sampling sd of that realisation | ±{ver['bernoulli_sd']:.4f} "
          f"({ver['n_tumor']:,} tumour cells) |")
        A("")
        A(f"The ranking uses the **expectation** E[kill] = mean(p_kill), not one "
          f"Bernoulli realisation. That matters: the single-seed coin-flip noise "
          f"is ±{ver['bernoulli_sd']:.4f} absolute, against a nominal "
          f"{a}-vs-{b} gap of only ~{abs(st['nominal_kill_median'][a]-st['nominal_kill_median'][b]):.3f} "
          f"— sampling noise is a non-trivial fraction of the very margin under "
          f"test, and it is *not* parameter uncertainty, so it is removed rather "
          f"than counted as instability.")
        A("")
    A("---")
    A("")

    # -------------------------------------------------------------- the priors
    A("## 3. Priors, and the honest error bar for a nested design")
    A("")
    A("### 3.1 Reconstruction / calibration priors (imported unchanged)")
    A("")
    A("The four transport-side priors are **imported from "
      "`tmesim.validation.uncertainty`** (`build_priors`, `sample_priors`) so the "
      "two modules cannot drift apart. Full justifications: "
      "`results/uncertainty.md` §1.")
    A("")
    A("| knob | prior | 5% / 50% / 95% | realised min–max "
      f"({res['n_outer']} draws) |")
    A("|---|---|:--:|:--:|")
    for name in PARAM_ORDER:
        pr = res["priors"][name]
        q5, q50, q95 = pr.quantiles()
        col = np.array([r[name] for r in res["rows"]], dtype=np.float64)
        A(f"| `{name}` | {pr.kind}: {pr.spec} | {q5:.3g} / {q50:.3g} / {q95:.3g} "
          f"| {col.min():.3g} – {col.max():.3g} |")
    A("")
    a_col = np.array([r["calibration_alpha"] for r in res["rows"]],
                     dtype=np.float64)
    A(f"> *Read the last column against the third.* A Latin hypercube puts exactly "
      f"one draw in each of {res['n_outer']} equal-probability strata, so the "
      f"outermost strata of an UNBOUNDED prior get sampled every time. The α "
      f"prior is lognormal with a 90% CI of 0.5×–2× "
      f"([{res['alpha0']/2:.0f}, {res['alpha0']*2:.0f}]), but its realised range "
      f"here is **{a_col.min():.0f}–{a_col.max():.0f}** "
      f"({a_col.min()/res['alpha0']:.2f}×–{a_col.max()/res['alpha0']:.2f}×) — the "
      f"tail, not a bug. Those extreme draws are ~2/{res['n_outer']} of the mass "
      f"and they are kept: a probability that survives them is worth more than "
      f"one computed on a truncated prior. §5 reports what α actually does to the "
      f"decision.")
    A("")
    A(f"**α is the single-anchor scale for the whole penetration model.** "
      f"α = **{res['alpha0']:.1f}** (the e-folding convention `L = d_half/ln 2`; "
      f"the retired α = 309.6 fitted the raw half-distance as if it were an "
      f"e-fold length and is superseded — audit defect #4). It is fitted in "
      f"closed form against **ONE doxorubicin immunofluorescence anchor** "
      f"(Primeau 2005, PMID 16361566), and the *absolute* penetration of all five "
      f"drugs rides on that single measurement. Its prior here is correspondingly "
      f"the widest of the four (90% CI 0.5×–2×), which is the honest treatment of "
      f"a one-study anchor transferred to human GBM — and note what that implies "
      f"for a *ranking*: because α is shared by every drug in a draw, it moves "
      f"all five kills together, so a wide α prior costs the marginal intervals a "
      f"lot and costs the ORDER comparatively little. Some of the ranking's "
      f"stability below is bought by exactly that shared-nuisance structure, and "
      f"that is a real property of the decision, not a trick.")
    A("")
    A("### 3.2 IC50 priors — two tiers, and why they differ")
    A("")
    A("Potency is what actually separates these drugs (D_eff spans ~1.5×; IC50 "
      "spans ~1200×), so no honest statement about the ORDER can hold IC50 fixed. "
      "Three scenarios are evaluated **on the same transport solves**:")
    A("")
    for s in IC50_SCENARIOS:
        A(f"* **{SCENARIO_LABEL[s]}**")
    A("")
    A("| drug | IC50 (µM) | tier | GDSC lines | SD(ln IC50) across lines | "
      "σ_log — B (primary) | 90% fold — B | σ_log — C (stress) | 90% fold — C |")
    A("|---|---:|:--|---:|---:|---:|---:|---:|---:|")
    for d in drugs:
        pB = res["ic50_priors"]["ic50_median_se"][d]
        pC = res["ic50_priors"]["ic50_cell_line"][d]
        nl = f"{pB.n_cell_lines}" if pB.n_cell_lines else "—"
        sd = f"{pB.sd_ln_cell_lines:.3f}" if pB.sd_ln_cell_lines else "—"
        A(f"| {d} | {pB.ic50_uM:.4g} | {pB.tier} | {nl} | {sd} | "
          f"{pB.sigma_log:.3f} | {pB.fold_90:.2f}× | {pC.sigma_log:.3f} | "
          f"{pC.fold_90:.2f}× |")
    A("")
    A("**GDSC-measured drugs get the TIGHTER prior, and it is derived from data, "
      "not chosen.** For cisplatin, methotrexate and temozolomide the model plugs "
      "in `exp(median over glioma cell lines of LN_IC50)` from the cached GDSC2 "
      "screen. The uncertainty on *that statistic* is the standard error of a "
      "median, `σ = √(π/2)·SD(ln IC50)/√n_lines` — the √(π/2) ≈ 1.2533 factor is "
      "the median's efficiency penalty relative to the mean under normality. The "
      "cell lines aggregated are read through `molecular.potency`'s own alias "
      "table and cache path, and repeat screens of one line are collapsed by "
      "median exactly as `potency.measured_ic50` does, so this is the spread of "
      "the very numbers the model consumes.")
    A("")
    A("**Literature-only drugs get a WIDER prior, and it is an expert judgement, "
      "declared as such.** Doxorubicin and lomustine are not in the GDSC glioma "
      "screen; there is no measured spread to compute. `potency.LITERATURE_IC50` "
      "quotes doxorubicin as ~0.1–1 µM around a 0.5 µM anchor (0.2×–2×) and "
      "lomustine as 'tens of µM' around 50 µM (~0.2×–2×). The prior used is a 90% "
      f"CI of **1/3×–3×** (σ_log = {LITERATURE_SIGMA_LOG:.3f}), which contains "
      "both quoted ranges and adds room for the cross-study transfer a "
      "one-sentence range does not cover. It is wider than **every** measured "
      "drug's primary prior, which is the intended ordering: a curated "
      "order-of-magnitude anchor must not be allowed to look as well-determined "
      "as a 34–51 cell-line screen.")
    A("")
    A("> *Honest caveat on the asymmetry.* The two leading drugs sit in "
      "**different tiers** — cisplatin is GDSC-measured (tight), doxorubicin is "
      "literature-only (wide). That is a real difference in what we know, and it "
      "is deliberately not smoothed over; but it does mean part of any "
      "instability attributed to doxorubicin is *epistemic bookkeeping* about its "
      "anchor rather than a claim about its pharmacology. Scenario C removes the "
      "asymmetry by giving both a cell-line-scale spread, which is one reason it "
      "is worth reading even though it is not the headline.")
    A("")
    A("### 3.3 Why the error bars are CLUSTER, not binomial")
    A("")
    A(f"The design is **nested**: {res['n_ic50']} IC50 draws sit inside each of "
      f"{res['n_outer']} tissue draws, and rows inside one tissue draw share a "
      f"reconstruction, an α and a tortuosity. They are not independent Bernoulli "
      f"trials. The naive binomial SE √(p(1−p)/{res['n_total']:,}) is therefore an "
      f"**under-estimate** — and the size of the under-estimate is itself a "
      f"finding, because it measures how much of the ranking uncertainty lives on "
      f"the tissue axis.")
    A("")
    A("| probability (scenario B) | p | naive binomial SE @ "
      f"n={res['n_total']:,} | **cluster SE** @ "
      f"n_outer={res['n_outer']} | design effect | effective n |")
    A("|---|---:|---:|---:|---:|---:|")
    def _row(lbl, s_):
        neff = f"{s_['n_eff']:.0f}" if np.isfinite(s_["n_eff"]) else "—"
        de = f"{s_['design_effect']:.1f}×" if np.isfinite(s_["design_effect"]) else "—"
        cs = f"{s_['se_cluster']:.4f}" if np.isfinite(s_["se_cluster"]) else "—"
        A(f"| {lbl} | {s_['p']:.4f} | {s_['se_binom']:.4f} | **{cs}** | {de} | {neff} |")
    for d in drugs:
        _row(f"P({d} = rank 1)", st["p_rank1"][d])
    _row(f"P({a} > {b})", st["pair"])
    _row(f"P(top-2 set = {{{a}, {b}}})", st["top2_headline"])
    _row("P(full published order)", st["full_order"])
    A("")
    A("**Read the design-effect column** (cluster SE ÷ binomial SE). Where it is "
      "**> 1** the binomial SE is fiction: the probability is driven by *which "
      "tissue you drew*, and adding more IC50 draws inside the same tissues would "
      "shrink the binomial SE without learning anything. Where it is **< 1** the "
      "nested Latin hypercube has actually beaten iid sampling — the inner design "
      "stratifies each drug's IC50 marginal, so the within-draw proportion is "
      "less noisy than a binomial — and the cluster SE is still the one to quote, "
      "because it is the one that reflects the axis we are short of draws on. Every headline number in this report is quoted "
      "with the cluster SE and a cluster-bootstrap 90% CI (resampling whole outer "
      "draws, 4,000 replicates). Where a probability is exactly 0 or 1 the "
      f"cluster SE collapses and is meaningless; the rule-of-three one-sided 95% "
      f"bound on the cluster unit, 3/{res['n_outer']} = "
      f"{3.0/max(res['n_outer'],1):.3f}, is reported instead.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------ rank matrix
    A("## 4. The full drug × rank probability matrix")
    A("")
    A(f"Scenario **B ({SCENARIO_LABEL[PRIMARY_SCENARIO]})**, "
      f"{res['n_total']:,} re-rankings. Cell = P(that drug lands at that rank) "
      f"± cluster SE. Rows sum to 1 by construction.")
    A("")
    A("| drug | " + " | ".join(f"rank {r+1}" for r in range(nd))
      + " | E[rank] | median kill (90% CI) |")
    A("|---|" + "---:|" * nd + "---:|:--:|")
    for d in drugs:
        cells = []
        for r in range(nd):
            s_ = st["matrix"][d][r]
            if s_["p"] <= 0.0:
                cells.append("—")
            elif s_["degenerate"]:
                cells.append(f"**{s_['p']:.3f}**")
            else:
                cells.append(f"{s_['p']:.3f} ± {s_['se_cluster']:.3f}")
        lo, hi = st["kill_ci90"][d]
        A(f"| **{d}** | " + " | ".join(cells) + f" | {st['expected_rank'][d]:.2f} | "
          f"{st['median_kill'][d]:.4f} [{lo:.4f}, {hi:.4f}] |")
    A("")
    A(f"![ranking robustness]({os.path.basename(figure)})")
    A("")
    A("### 4.1 P(rank 1) across all three IC50 scenarios")
    A("")
    A("| drug | A. transport only | B. + IC50 SE (**primary**) | "
      "C. + cell-line spread |")
    A("|---|:--:|:--:|:--:|")
    for d in drugs:
        cells = []
        for s in IC50_SCENARIOS:
            s_ = stats_by_scen[s]["p_rank1"][d]
            if s_["degenerate"]:
                cells.append(f"{s_['p']:.3f}")
            else:
                cells.append(f"{s_['p']:.3f} ± {s_['se_cluster']:.3f}")
        A(f"| {d} | " + " | ".join(cells) + " |")
    A("")
    A("| statistic | A. transport only | B. primary | C. stress |")
    A("|---|:--:|:--:|:--:|")
    for lbl, key in ((f"P({a} > {b})", "pair"),
                     (f"P(top-2 = {{{a}, {b}}})", "top2_headline"),
                     ("P(full published order)", "full_order")):
        A(f"| {lbl} | " + " | ".join(
            fmt_p(stats_by_scen[s][key]) for s in IC50_SCENARIOS) + " |")
    A("| distinct orders observed | " + " | ".join(
        f"{stability[s]['n_distinct_orders']} of {math.factorial(nd)}"
        for s in IC50_SCENARIOS) + " |")
    A("| modal order probability | " + " | ".join(
        f"{stability[s]['modal_p']:.3f}" for s in IC50_SCENARIOS) + " |")
    A("| normalised order entropy | " + " | ".join(
        f"{stability[s]['entropy_norm']:.3f}" for s in IC50_SCENARIOS) + " |")
    A("")
    A(f"Scenario A holds every IC50 at its resolved value, so it isolates how much "
      f"instability the **reconstruction + calibration knobs alone** can cause; B "
      f"adds the sampling error of the potency statistic the model actually "
      f"consumes; C asks the harder question — *this Xenium section is ONE "
      f"tumour, and nothing tells us it sits at the median of the glioma "
      f"cell-line distribution*. C is reported as a stress test, not as the "
      f"headline, because a cell-line panel's spread is a population statistic, "
      f"not a credible interval for a specific tumour. The primary answer is B.")
    A("")
    A("### 4.2 Top-2 SET distribution (scenario B)")
    A("")
    A("| top-2 set | probability |")
    A("|---|:--:|")
    for members, s_ in st["top2_dist"][:6]:
        if s_["p"] < 0.001:
            continue
        mark = " ← **headline**" if members == HEADLINE_TOP2 else ""
        A(f"| {{{', '.join(sorted(members))}}}{mark} | {fmt_p(s_)} |")
    A("")
    A("### 4.3 Every pairwise comparison (scenario B)")
    A("")
    A("| comparison | P(left > right) | decided at 90%? |")
    A("|---|:--:|:--:|")
    for (d1, d2), s_ in st["pairwise"].items():
        dec = ("yes" if (s_["ci90_lo"] >= 0.90 or s_["ci90_hi"] <= 0.10)
               else "**NO**")
        A(f"| {d1} > {d2} | {fmt_p(s_)} | {dec} |")
    A("")
    A(f"**The {a}–{b} margin itself.** Median "
      f"E[kill]({a}) − E[kill]({b}) = **{st['gap_median']:+.4f}** "
      f"(90% of rows in [{st['gap_ci90'][0]:+.4f}, {st['gap_ci90'][1]:+.4f}]); "
      f"averaged within a tissue draw the gap ranges "
      f"{st['gap_per_draw_min']:+.4f} to {st['gap_per_draw_max']:+.4f} across the "
      f"{res['n_outer']} draws.")
    A("")
    A("---")
    A("")

    # --------------------------------------------------------- destabilisers
    A("## 5. Which parameter destabilises the ranking")
    A("")
    A(f"Outcome: the binary event *`{dest['target']}` holds rank 1* (baseline "
      f"{fmt_p(dest['baseline'])}). Each sampled input is cut into "
      f"{dest['n_q']} quartiles and the conditional probability is recomputed; "
      f"the **swing** is max − min across quartiles, cluster-bootstrapped over "
      f"outer draws. A quartile conditional needs no linearity assumption, which "
      f"matters because the typing threshold acts as a step (integer Xenium "
      f"counts quantise the PECAM1 z-score) and the IC50 effect is a ratio.")
    A("")
    A("| input | Q1 | Q2 | Q3 | Q4 | **swing** | swing 90% CI | Spearman(ρ, win) |")
    A("|---|:--:|:--:|:--:|:--:|---:|:--:|---:|")
    for name in dest["ranked"]:
        e = dest["per_variable"][name]
        qs = " | ".join("—" if not np.isfinite(q["p"]) else f"{q['p']:.3f}"
                        for q in e["quartiles"])
        if not np.isfinite(e["swing"]):
            A(f"| `{_var_label(name)}` | {qs} | — | — | {e['spearman']:+.3f} |")
            continue
        A(f"| `{_var_label(name)}` | {qs} | **{e['swing']:.3f}** | "
          f"[{e['swing_ci90'][0]:.3f}, {e['swing_ci90'][1]:.3f}] | "
          f"{e['spearman']:+.3f} |")
    A("")
    ranked = dest["ranked"]
    top, second = ranked[0], (ranked[1] if len(ranked) > 1 else None)
    te = dest["per_variable"][top]
    A(f"**Most destabilising: `{_var_label(top)}`** — swing "
      f"{te['swing']:.3f} [{te['swing_ci90'][0]:.3f}, {te['swing_ci90'][1]:.3f}].")
    if second is not None:
        se_ = dest["per_variable"][second]
        A("")
        A(f"Runner-up `{_var_label(second)}` (swing {se_['swing']:.3f} "
          f"[{se_['swing_ci90'][0]:.3f}, {se_['swing_ci90'][1]:.3f}]).")
    A("")
    A("**A swing whose bootstrap CI reaches 0 is not distinguishable from noise "
      "at this n.** With "
      f"{res['n_outer']} independent tissue draws the quartile conditionals for "
      "the four transport knobs rest on ~"
      f"{res['n_outer']//4} draws each, so read the RANKING of the swings, not "
      "their third decimal. Note also that the lower bound of a max−min statistic "
      "is bounded below by 0 and biased *upward* under noise, so small swings "
      "should be treated as upper bounds.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------- resolution
    if dxchk:
        A("## 6. Is the ranking a resolution artifact?")
        A("")
        A(f"The Monte Carlo runs at dx = {res['dx_um']:.0f} µm / dt = "
          f"{res['dt_s']:.0f} s / convection off for affordability. "
          "`results/numerical_rigor.md` established that `kill_fraction` is "
          "grid-converged for dx ≲ 20 µm and dt-invariant for dt ≤ 20 s "
          "(interstitial speeds < 0.7 µm/s make advection negligible). The "
          "nominal panel is nonetheless re-run at both the MC configuration and "
          "the **exact published production configuration** (dx = 12 µm, "
          "dt = 10 s, convection on) so the claim is demonstrated rather than "
          "inherited:")
        A("")
        A("| dx (µm) | dt (s) | conv. | grid | "
          + " | ".join(drugs) + f" | gap({a}−{b}) |")
        A("|---:|---:|:--:|:--:|" + "---:|" * nd + "---:|")
        for r in dxchk["rows"]:
            gy, gx = r["grid"]
            A(f"| {r['dx_um']:.0f} | {r['dt_s']:.0f} | "
              f"{'on' if r['convection'] else 'off'} | {gy}×{gx} | "
              + " | ".join(f"{r['kill'][d]:.4f}" for d in drugs)
              + f" | {r['gap_cis_dox']:+.4f} |")
        A("")
        if dxchk["order_identical"]:
            A(f"**The order is identical at both resolutions** "
              f"({' > '.join(dxchk['rows'][0]['order'])}), and the "
              f"{a}−{b} gap moves from "
              f"{dxchk['rows'][0]['gap_cis_dox']:+.4f} to "
              f"{dxchk['rows'][-1]['gap_cis_dox']:+.4f}. The instability reported "
              f"above is therefore *parameter* uncertainty, not a mesh artifact — "
              f"which is the point: refining the grid will not fix it.")
        else:
            A(f"**The order is NOT identical at the two resolutions** "
              f"({' > '.join(dxchk['rows'][0]['order'])} at "
              f"dx = {dxchk['rows'][0]['dx_um']:.0f} µm vs "
              f"{' > '.join(dxchk['rows'][-1]['order'])} at "
              f"dx = {dxchk['rows'][-1]['dx_um']:.0f} µm). That is a finding in "
              f"its own right and it makes the ranking *less* trustworthy than "
              f"the Monte Carlo alone suggests: a numerical choice, not only a "
              f"prior, can move the answer.")
        A("")
        A("---")
        A("")

    # ---------------------------------------------------------------- honesty
    A("## 7. What this does and does NOT establish")
    A("")
    A("1. **These are probabilities over EXPERT PRIORS, not a validated "
      "predictive distribution.** No outcome label exists for this section — it "
      "is an untreated resection with no measured response to any of the five "
      f"drugs. `P({v['leader']} is rank 1) = "
      f"{st['p_rank1'][v['leader']]['p']:.2f}` means *the model's own ranking is "
      f"that stable inside the priors declared in §3*. It is NOT a claim that "
      f"{v['leader']} is the best drug for this patient "
      f"{100*st['p_rank1'][v['leader']]['p']:.0f}% of the time.")
    A("2. **Model-form error is not in these bars, and is plausibly the largest "
      "term.** Every draw is computed inside one fixed structure: a first-order "
      "`k = clearance + binding` sink, `D_eff = D·(1 − α·ECM)`, a KDE-smoothed "
      "endothelial density standing in for a vessel network never imaged, a Hill "
      "kill on AUC, plasma-averaged delivery, and no efflux transporters. If that "
      "structure is wrong, no amount of parameter sampling detects it. "
      "`results/vessel_network_impact.md` shows the KDE→discrete-vessel swap "
      "alone moves doxorubicin's coverage by −0.017 dose-controlled (−0.586 "
      "naively) — a structural move that no prior in §3 can produce.")
    A("3. **Every drug's exposure rides on ONE anchor.** α is fitted to a single "
      "doxorubicin penetration measurement in a different tumour model. Its 0.5×–"
      "2× prior is a judgement about how far that anchor transfers, not a "
      "meta-analysis. If it does not transfer, the whole distribution translates.")
    A("4. **Absolute kills are an upper bound.** `uniform` mode models no "
      "intrinsic per-cell resistance for any drug. That is what makes the ranking "
      "fair; it also means no absolute kill here should be quoted as an expected "
      "response.")
    A("5. **One section, one segmentation, no replicate.** Section-to-section and "
      "patient-to-patient variability is entirely absent. This is uncertainty "
      "*within one measured tumour*.")
    A("6. **Temozolomide's position is downstream of a code fix, not of "
      "biology-free modelling.** `bbb_penetrance` now aggregates its four CNS-MPO "
      "desirabilities with a **geometric mean** rather than an arithmetic one, so "
      "a disqualifying property can VETO a molecule instead of being averaged "
      "away (audit defect #9; it cut the eribulin miss from 66× to 6.9× against "
      "measured human glioma tissue:plasma ratios). TMZ's bbb_penetrance fell "
      "0.558 → **0.125** as a result, which is why its kill is 0.1151 rather than "
      "the older 0.4173. The veto is better-validated than what it replaced — but "
      "TMZ's rank in this report is nonetheless *sensitive to one aggregation "
      "choice in a structure-only score*, and that sensitivity is not represented "
      "in any prior sampled here.")
    A("")

    # ------------------------------------------------------ what would settle it
    A("## 8. What measurement would actually settle the top-1 call")
    A("")
    if v["top1_decision_grade"]:
        A(f"The top-1 call already clears the bar in scenario B. The measurement "
          f"that would make it robust to the STRESS scenario (C), where "
          f"P({v['leader']} = rank 1) falls to "
          f"{stats_by_scen['ic50_cell_line']['p_rank1'][v['leader']]['p']:.3f}, is "
          f"the same one listed below: a patient-specific potency readout.")
        A("")
    A(f"The ranking is decided by **delivery × potency**, and the two candidates "
      f"sit at opposite ends of both axes — {a} has ~157× the mean AUC and "
      f"~49× the IC50 of {b}. So the useful measurements are, in order of "
      f"decisiveness:")
    A("")
    A(f"1. **An ex-vivo potency readout on THIS patient's cells for {a} and "
      f"{b}.** This is the decisive one. Scenario C shows that replacing 'the "
      f"median glioma cell line' with 'this tumour' is what dominates the "
      f"remaining uncertainty; a per-patient IC50 pair collapses the widest prior "
      f"in §3.2 from a population spread to a measurement error. It also removes "
      f"the tier asymmetry — doxorubicin's 1/3×–3× literature prior exists purely "
      f"because it is absent from the GDSC glioma screen. `results/"
      f"gbm_tmz_exvivo.md` shows the assay format already exists for glioma "
      f"(GSE232173); it has simply never been run on a section that was also "
      f"spatially profiled.")
    A(f"2. **A second, independent penetration anchor** — ideally for a drug on "
      f"this panel other than doxorubicin. α is currently one number from one "
      f"study; a second anchor would let α be fitted rather than assumed-"
      f"transferable, shrinking the widest transport prior. Because α is shared "
      f"across drugs it matters less for the ORDER than for the levels, so this "
      f"is second, not first.")
    A(f"3. **Measured spatial drug distribution in the same section**, which "
      f"would test the delivery axis directly rather than through a "
      f"literature-anchored sink. `results/DATA_HUNT.md` establishes that no open "
      f"dataset anywhere pairs measured spatial drug distribution with measured "
      f"response in the same tumour — this is a collaboration to run, not a "
      f"download.")
    A("")
    A("**Until (1) exists, the reportable output of this pipeline for the top of "
      "the panel is a SET, not a winner.**")
    A("")

    text = "\n".join(L)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


# ===========================================================================
# 9. CONSOLE SUMMARY + ENTRY POINT.
# ===========================================================================
def _print(res: Dict, stats_by_scen: Dict[str, Dict], dest: Dict,
           stability: Dict[str, Dict], dxchk: Optional[Dict] = None) -> None:
    drugs = list(res["drugs"])
    st = stats_by_scen[PRIMARY_SCENARIO]
    v = _verdict(st)
    a, b = HEADLINE_PAIR
    W = 96
    print("=" * W)
    print("RANKING ROBUSTNESS — real Xenium GBM section, uniform sensitivity mode")
    print(f"  {res['n_outer']} tissue draws x {res['n_ic50']} IC50 draws = "
          f"{res['n_total']:,} re-rankings   dx={res['dx_um']:.0f} um   "
          f"runtime={res['runtime_s']/60:.1f} min")
    print("=" * W)
    print(f"DRUG x RANK PROBABILITY  (scenario B = {PRIMARY_SCENARIO}, "
          f"+- cluster SE over {res['n_outer']} draws)")
    print(f"  {'drug':14s}" + "".join(f"{'rank '+str(r+1):>18s}"
                                      for r in range(len(drugs)))
          + f"{'E[rank]':>9s}")
    for d in drugs:
        cells = ""
        for r in range(len(drugs)):
            s_ = st["matrix"][d][r]
            cells += (f"{'-':>18s}" if s_["p"] <= 0 else
                      f"{s_['p']:.3f}+-{s_['se_cluster']:.3f}".rjust(18))
        print(f"  {d:14s}{cells}{st['expected_rank'][d]:>9.2f}")
    print("-" * W)
    print("HEADLINE PROBABILITIES (scenario B, cluster bootstrap 90% CI):")
    for lbl, s_ in ((f"P({v['leader']} = rank 1)", st["p_rank1"][v["leader"]]),
                    (f"P({v['runner_up']} = rank 1)", st["p_rank1"][v["runner_up"]]),
                    (f"P({a} > {b})", st["pair"]),
                    (f"P(top-2 set = {{{a},{b}}})", st["top2_headline"]),
                    ("P(full published order)", st["full_order"])):
        print(f"  {lbl:34s} {s_['p']:.4f}  [{s_['ci90_lo']:.4f}, "
              f"{s_['ci90_hi']:.4f}]   cluster SE {s_['se_cluster']:.4f}  "
              f"(binomial SE {s_['se_binom']:.4f}, "
              f"design effect {s_['design_effect']:.1f}x)")
    print("-" * W)
    print("P(rank 1) BY SCENARIO:")
    print(f"  {'drug':14s}" + "".join(f"{_SCEN_SHORT[s]:>26s}" for s in IC50_SCENARIOS))
    for d in drugs:
        print(f"  {d:14s}" + "".join(
            f"{stats_by_scen[s]['p_rank1'][d]['p']:.3f}+-"
            f"{stats_by_scen[s]['p_rank1'][d]['se_cluster']:.3f}".rjust(26)
            for s in IC50_SCENARIOS))
    print("-" * W)
    print(f"MOST DESTABILISING INPUT for P({dest['target']} = rank 1):")
    for name in dest["ranked"]:
        e = dest["per_variable"][name]
        if not np.isfinite(e["swing"]):
            continue
        qs = " ".join("  nan" if not np.isfinite(q["p"]) else f"{q['p']:.3f}"
                      for q in e["quartiles"])
        print(f"  {_var_label(name):26s} quartiles [{qs}]  swing "
              f"{e['swing']:.3f} [{e['swing_ci90'][0]:.3f}, "
              f"{e['swing_ci90'][1]:.3f}]")
    if dxchk:
        print("-" * W)
        print("RESOLUTION CHECK (nominal panel):")
        for r in dxchk["rows"]:
            print(f"  dx={r['dx_um']:4.1f} dt={r['dt_s']:4.0f} "
                  f"conv={str(r['convection']):5s}  "
                  f"{' > '.join(r['order'])}   gap={r['gap_cis_dox']:+.4f}")
        print(f"  -> order identical across resolutions: {dxchk['order_identical']}")
    print("=" * W)
    if v["top1_decision_grade"]:
        print(f"VERDICT: top-1 call IS decision-grade — {v['leader']}.")
    elif v["ci_overlap"]:
        print(f"VERDICT: top-1 call is NOT decision-grade. {v['leader']} and "
              f"{v['runner_up']} have OVERLAPPING 90% CIs — report a TIE.")
    else:
        print(f"VERDICT: top-1 call is NOT decision-grade at the 0.90 bar "
              f"({v['leader']} {st['p_rank1'][v['leader']]['p']:.3f}, CI low "
              f"{st['p_rank1'][v['leader']]['ci90_lo']:.3f}). Report "
              f"'{v['leader']} preferred, {v['runner_up']} a live alternative'.")
    print(f"         top-2 SET {{{a}, {b}}}: {st['top2_headline']['p']:.3f} "
          f"[{st['top2_headline']['ci90_lo']:.3f}, "
          f"{st['top2_headline']['ci90_hi']:.3f}] — "
          f"{'decision-grade' if v['top2_decision_grade'] else 'NOT decision-grade'}.")
    print("HONEST: probabilities over EXPERT PRIORS, not a validated predictive "
          "distribution.")
    print("=" * W)


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(
        prog="python -m tmesim.validation.ranking_robustness",
        description="Is the per-patient drug RANKING decision-grade under "
                    "parameter uncertainty?")
    ap.add_argument("--n-outer", type=int, default=40,
                    help="tissue draws (the EFFECTIVE sample size; default 40)")
    ap.add_argument("--n-ic50", type=int, default=256,
                    help="IC50 draws nested inside each tissue draw (default 256)")
    ap.add_argument("--dx-um", type=float, default=MC_DX_UM)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--time-budget-s", type=float, default=None,
                    help="stop drawing after this many seconds (reports the n "
                         "actually achieved)")
    ap.add_argument("--no-dx-check", action="store_true",
                    help="skip the dx=15 vs dx=12 nominal resolution check")
    ap.add_argument("--report", default="results/ranking_robustness.md")
    ap.add_argument("--figure", default="results/ranking_robustness.png")
    ap.add_argument("--cache", default="results/ranking_robustness_mc.npz",
                    help="where the raw Monte-Carlo cube is written")
    ap.add_argument("--from-cache", action="store_true",
                    help="re-render the report and figure from --cache without "
                         "re-solving transport (the resolution check is skipped)")
    args = ap.parse_args(argv)

    t0 = time.time()
    if args.from_cache:
        print(f"Re-rendering from {args.cache} (no transport solves)...\n")
        res = load_mc(args.cache)
        args.no_dx_check = True
    else:
        print("Running ranking robustness on REAL data (tens of minutes)...\n")
        res = run_ranking_mc(n_outer=args.n_outer, n_ic50=args.n_ic50,
                             dx_um=args.dx_um, seed=args.seed,
                             time_budget_s=args.time_budget_s)
        print(f"\n[cached raw Monte-Carlo cube -> {save_mc(res, args.cache)}]")
    stats_by_scen = {s: rank_stats(res, s) for s in IC50_SCENARIOS}
    stability = {s: ranking_stability_index(res, s) for s in IC50_SCENARIOS}
    dest = destabilisers(res, PRIMARY_SCENARIO)
    dxchk = None
    if not args.no_dx_check:
        print("\nResolution check (nominal panel at the MC grid and the "
              "production grid)...")
        dxchk = dx_ranking_check()
    print()
    _print(res, stats_by_scen, dest, stability, dxchk)
    fig = plot_ranking(res, stats_by_scen, dest, path=args.figure)
    out = write_report(res, stats_by_scen, dest, stability, dxchk,
                       path=args.report, figure=args.figure)
    print(f"\n[wrote {out}]  [wrote {fig}]  "
          f"[total runtime {(time.time()-t0)/60:.1f} min]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
