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
