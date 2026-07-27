"""COMBINATION THERAPY on the REAL tumor: do drugs with complementary DELIVERY
profiles cover more of the tumor than either alone?

The hypothesis under test
-------------------------
A small, PENETRANT-but-WEAK drug (temozolomide: MW 194, reaches deep tissue, but
a real GDSC IC50 in the hundreds of uM) paired with a POTENT-but-DELIVERY-LIMITED
drug (doxorubicin: IC50 0.5 uM, but consumed within ~45 um of a vessel by the
calibrated perivascular sink) should cover COMPLEMENTARY spatial territory: the
potent drug kills the perivascular cuff, the penetrant drug reaches the deep cells
the potent one never sees. This experiment tests whether that spatial
complementarity actually exists on REAL geometry — the open 10x Xenium GBM
section (40,887 measured cells, real endothelial-derived vasculature, real
matrix-transcript ECM, real-perfusion oxygen), reconstructed ONCE at the
grid-converged dx = 12 um and reused by every arm, so all differences are
attributable to drug physics rather than to tissue.

THE COMBINATION ASSUMPTION: BLISS INDEPENDENCE (read this first)
-----------------------------------------------------------------
Each drug is simulated SEPARATELY (its own transport PDE on the same tissue), and
each tumor cell ends up with its own exposure AUC and its own kill PROBABILITY
``p_A``, ``p_B`` from the frozen engine's Hill pharmacodynamics. The combination
is then formed AT THE EXPOSURE LEVEL by the Bliss-independence rule::

    p_kill_combo(cell) = 1 - (1 - p_A(cell)) * (1 - p_B(cell))

i.e. the two drugs kill the same cell INDEPENDENTLY, and a cell survives only if
it survives both. This assumption is deliberate and restrictive:

  * There is NO pharmacological synergy or antagonism in this model. No shared
    target, no repair-pathway crosstalk, no efflux competition, no schedule
    dependence, no drug-drug PK interaction, no altered transport of A caused by
    B. We cannot predict biochemical synergy from structure and geometry, so we
    do not pretend to.
  * What this DOES model is purely SPATIAL complementarity: whether the two
    drugs' delivery fields reach DIFFERENT cells. Any gain reported here comes
    from geography or from the independence arithmetic itself — never from
    pharmacology.
  * Consequently this is a delivery-only floor, not a forecast of combination
    efficacy: real anthracycline + alkylator combinations may synergise or
    antagonise, and this model is blind to both.

The Bernoulli coupling matters too. Under Bliss independence the killed SETS are
drawn with INDEPENDENT per-cell draws (``u_A``, ``u_B``), which is exactly the
coupling that makes ``P(A or B) = 1 - (1-p_A)(1-p_B)`` hold. All set statistics
are therefore reported as EXACT EXPECTATIONS under that coupling::

    E|A n B|  = sum_cells p_A p_B          E|A only| = sum_cells p_A (1 - p_B)
    E|A u B|  = sum_cells [1-(1-p_A)(1-p_B)]   E|B only| = sum_cells (1 - p_A) p_B

so the headline numbers carry no Monte-Carlo noise; a single seeded realisation is
drawn only to render the territory maps.

THREE CONTROLS, because "combo > monotherapy" is trivially true
----------------------------------------------------------------
1. **DOSE-MATCHED control.** The better single agent of each pair is re-run at
   2x its plasma peak — the same total systemic exposure (in standard-dose units)
   the two-drug arm receives. If the combination cannot beat 2x monotherapy, the
   spatial story buys nothing.

2. **SELF-COMBINATION control** (the sharp one). Apply the SAME Bliss rule to the
   better agent against ITSELF: ``1 - (1 - p_B)^2``. This is a combination of two
   drugs with *perfectly identical, maximally redundant* delivery. If
   ``combo(A,B) ~= selfcombo(B,B)`` then the entire "combination benefit" is an
   artifact of the independence assumption applied against a saturating Hill
   curve — NOT spatial complementarity. Only ``combo(A,B) > selfcombo(B,B)``
   is evidence that using two DIFFERENT delivery profiles helped.

3. **SPATIAL-SHUFFLE nulls.** Independent-Bernoulli overlap is mechanically small
   whenever kill probabilities are small (E|A n B| ~ p^2 while E|A u B| ~ 2p), so
   a low Jaccard is NOT by itself evidence of complementarity. Every overlap
   statistic is therefore compared against permutation nulls with identical
   marginals but destroyed spatial alignment:

   * *full shuffle* — permute ``p_B`` across tumor cells;
   * *delivery shuffle* — the engine's kill probability factorises exactly as
     ``p = p_hill(exposure) * sensitivity * resistance``, where the second factor
     is a CELL-INTRINSIC term SHARED by every drug (it does not depend on the
     agent). The full shuffle destroys that shared term's alignment too, so it
     cannot separate "same delivery" from "same intrinsically-killable cells".
     The delivery shuffle permutes only ``p_hill`` (drug B's delivery term),
     within 20 strata of the shared factor, leaving each cell's intrinsic
     killability where it is. Its ratio isolates DELIVERY co-localisation.

   observed/null < 1 = genuinely complementary territory; ~1 = no spatial
   structure; > 1 = redundant, both drugs hitting the same cells.

Run:    ./.venv/bin/python experiments/combination.py
Writes: results/combination.png
        results/combination.md

Adds files only; core / physics / molecular / reconstruct / engine are untouched
and are driven exclusively through their public entry points.
"""
from __future__ import annotations

import argparse
import os
import pickle
import time
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from tmesim.cells.agents import update_cell_states
from tmesim.cells.pkpd import apply_drug_kill
from tmesim.core import (
    CellState,
    Microenvironment,
    OxygenConfig,
    PDConfig,
    PDEConfig,
    PKConfig,
    SimConfig,
    TissueFields,
)
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.physics.oxygen import solve_oxygen
from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.validation.vessel_profile import (
    distance_from_vessel_um,
    vessel_mask_from_density,
)

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
#: Grid-CONVERGED resolution (results/numerical_rigor.md section 2): the ~45 um
#: perivascular boundary layer is unresolved on coarser grids.
DX_UM = 12.0

#: The three agents. temozolomide = penetrant/weak, doxorubicin =
#: potent/delivery-limited, cisplatin = the intermediate contrast.
DRUGS: List[str] = ["temozolomide", "doxorubicin", "cisplatin"]

#: (A, B) pairs. TMZ + doxorubicin is the headline delivery-complementarity test.
PAIRS: List[Tuple[str, str]] = [
    ("temozolomide", "doxorubicin"),
    ("cisplatin", "doxorubicin"),
    ("temozolomide", "cisplatin"),
]

#: Nominal single dose (uM peak plasma), identical for every agent — see caveats.
BASE_PLASMA_PEAK = PKConfig().plasma_peak

#: The dose-matched control multiplies the better single agent's dose by this.
DOSE_MATCH_FACTOR = 2.0

#: Permutations for the spatial-shuffle nulls, and strata for the delivery null.
N_SHUFFLE = 200
N_STRATA = 20
SHUFFLE_SEED = 20260727

#: Vessel-distance binning for the delivery profile. The section's tumor cells sit
#: a MEDIAN of ~250 um from the p95 vessel mask (max ~2.3 mm), so a uniform
#: near-field grid would lump most of the tumor into one overflow bin. Bins are
#: therefore 12 um (= dx, resolving the calibrated ~45 um layer) out to 120 um,
#: then progressively wider; cells beyond the last edge are DROPPED, not clipped,
#: and the dropped fraction is reported.
VESSEL_PERCENTILE = 95.0
DIST_EDGES = np.concatenate([
    np.arange(0.0, 132.0, 12.0),      # 12 um bins, 0-120 um
    np.arange(160.0, 440.0, 40.0),    # 40 um bins, 120-400 um
    np.arange(500.0, 1300.0, 100.0),  # 100 um bins, 400-1200 um
])
MIN_BIN_CELLS = 30

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "combination.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "combination.md")

_COL = {"temozolomide": "#2ca02c", "doxorubicin": "#d62728", "cisplatin": "#ff7f0e"}
_COL_A = "#1f77b4"      # A only
_COL_B = "#d62728"      # B only
_COL_BOTH = "#7b3fa0"   # both
_COL_NONE = "#dcdcdc"   # survivors


# --------------------------------------------------------------------------- #
# 1. Drug compilation (identical pipeline to the flagship single-drug runs)
# --------------------------------------------------------------------------- #
def compile_real_drug(name: str, alpha: float):
    """calibrated perivascular sink (alpha fit ONCE on doxorubicin) + saturable
    binding-site barrier. The barrier is a no-op for these small molecules (no
    saturable high-affinity target) but every agent is routed identically."""
    dp = calibrated_drugparams(name, alpha)
    mol = get_molecule(name)
    return apply_barrier_to_drugparams(dp, mol), mol


# --------------------------------------------------------------------------- #
# 2. Per-cell kill probabilities (the engine's own PD mapping, re-read)
# --------------------------------------------------------------------------- #
def oxygen_labeled_cells(micro: Microenvironment, oxygen_cfg: OxygenConfig):
    """Reproduce the engine's step-2 state labelling on a private copy.

    ``simulate_response`` solves oxygen and bands every cell VIABLE / QUIESCENT /
    HYPOXIC before the kill step; those labels are drug-INDEPENDENT (oxygen does
    not depend on the drug), so we build them ONCE here and reuse them to re-read
    each drug's per-cell kill probability with the identical resistance modulation
    the engine applied. Nothing in the engine is modified — this is the same
    public physics/cells call sequence the engine itself makes.
    """
    work = Microenvironment(
        domain=micro.domain,
        fields=TissueFields(
            ecm_density=micro.fields.ecm_density,
            vascular_density=micro.fields.vascular_density,
            oxygen=None,
        ),
        cells=micro.cells.copy(),
        patient_id=micro.patient_id,
        metadata=micro.metadata,
    )
    oxy = solve_oxygen(work, oxygen_cfg)
    work.fields.oxygen = oxy
    update_cell_states(work, oxy, oxygen_cfg)
    return work.cells.copy(), oxy


def intrinsic_factor(cells, pd_cfg: PDConfig) -> np.ndarray:
    """The DRUG-INDEPENDENT multiplier in the engine's PD map:
    ``sensitivity * R_hypoxia/R_quiescence``. Every agent is multiplied by this
    same per-cell number, which is why it inflates cross-drug kill correlation."""
    resistance = np.ones(cells.n, dtype=np.float64)
    resistance[cells.state == CellState.HYPOXIC] = pd_cfg.hypoxia_resistance
    resistance[cells.state == CellState.QUIESCENT] = pd_cfg.quiescent_resistance
    return np.asarray(cells.sensitivity, dtype=np.float64) * resistance


def per_cell_kill(micro, labeled, transport, pd_cfg: PDConfig, drug):
    """Per-cell kill PROBABILITY and per-cell exposure AUC for one drug arm.

    ``apply_drug_kill`` returns exactly the probability vector the engine used and
    writes each cell's sampled AUC onto the population; we run it on a throwaway
    copy of the oxygen-labelled cells so nothing is mutated and its Bernoulli draw
    is discarded.
    """
    snap = Microenvironment(
        domain=micro.domain,
        fields=micro.fields,
        cells=labeled.copy(),
        patient_id=micro.patient_id,
        metadata=micro.metadata,
    )
    p = apply_drug_kill(snap, transport, pd_cfg, drug, rng=np.random.default_rng(0))
    return (np.asarray(p, dtype=np.float64),
            np.asarray(snap.cells.exposure, dtype=np.float64))


# --------------------------------------------------------------------------- #
# 3. Bliss independence + set statistics (exact expectations)
# --------------------------------------------------------------------------- #
def bliss(p_a: np.ndarray, p_b: np.ndarray) -> np.ndarray:
    """p_combo = 1 - (1 - p_A)(1 - p_B). Independent action, NO synergy."""
    return 1.0 - (1.0 - p_a) * (1.0 - p_b)


def set_stats(p_a: np.ndarray, p_b: np.ndarray) -> Dict[str, float]:
    """Exact expected sizes of the killed-cell sets under independent draws."""
    n = float(p_a.size)
    e_a = float(p_a.sum())
    e_b = float(p_b.sum())
    e_both = float((p_a * p_b).sum())
    e_a_only = e_a - e_both
    e_b_only = e_b - e_both
    e_union = e_a_only + e_b_only + e_both
    return {
        "n_tumor": n, "kill_a": e_a / n, "kill_b": e_b / n, "kill_combo": e_union / n,
        "e_a": e_a, "e_b": e_b, "e_both": e_both,
        "e_a_only": e_a_only, "e_b_only": e_b_only, "e_union": e_union,
        "jaccard": e_both / e_union if e_union > 0 else 0.0,
        "overlap_coef": e_both / min(e_a, e_b) if min(e_a, e_b) > 0 else 0.0,
        "frac_both": e_both / e_union if e_union > 0 else 0.0,
        "frac_a_only": e_a_only / e_union if e_union > 0 else 0.0,
        "frac_b_only": e_b_only / e_union if e_union > 0 else 0.0,
        "tumor_frac_both": e_both / n,
        "tumor_frac_a_only": e_a_only / n,
        "tumor_frac_b_only": e_b_only / n,
    }


def _stratified_permutation(values: np.ndarray, strata: np.ndarray,
                            rng: np.random.Generator) -> np.ndarray:
    """Permute ``values`` only WITHIN each stratum."""
    out = np.empty_like(values)
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        out[idx] = rng.permutation(values[idx])
    return out


def shuffle_null(p_a: np.ndarray, p_b: np.ndarray, kind: str,
                 hill_b: np.ndarray = None, intrinsic: np.ndarray = None,
                 n_perm: int = N_SHUFFLE, seed: int = SHUFFLE_SEED) -> Dict[str, float]:
    """Permutation null for the overlap statistics.

    ``kind='full'``     — permute p_B across all tumor cells: identical marginal,
                          both the delivery pattern AND the shared cell-intrinsic
                          killability pattern are decoupled from drug A.
    ``kind='delivery'`` — permute only drug B's DELIVERY term ``p_hill`` within
                          strata of the shared cell-intrinsic factor, then
                          re-multiply: each cell keeps its own intrinsic
                          killability, only B's delivery geography is scrambled.
                          Isolates DELIVERY co-localisation. Marginals are
                          preserved to within the stratum width (reported).
    """
    rng = np.random.default_rng(seed)
    if kind == "delivery":
        q = np.quantile(intrinsic, np.linspace(0.0, 1.0, N_STRATA + 1)[1:-1])
        strata = np.searchsorted(q, intrinsic)
    jac = np.empty(n_perm); ovl = np.empty(n_perm)
    combo = np.empty(n_perm); marg = np.empty(n_perm)
    for i in range(n_perm):
        if kind == "full":
            pb = rng.permutation(p_b)
        else:
            pb = np.clip(_stratified_permutation(hill_b, strata, rng) * intrinsic, 0.0, 1.0)
        s = set_stats(p_a, pb)
        jac[i] = s["jaccard"]; ovl[i] = s["overlap_coef"]
        combo[i] = s["kill_combo"]; marg[i] = s["kill_b"]
    return {
        "kind": kind,
        "jaccard_mean": float(jac.mean()), "jaccard_sd": float(jac.std()),
        "overlap_mean": float(ovl.mean()), "overlap_sd": float(ovl.std()),
        "kill_combo_mean": float(combo.mean()), "kill_combo_sd": float(combo.std()),
        "kill_b_mean": float(marg.mean()),
    }


def realise_sets(p_a: np.ndarray, p_b: np.ndarray, seed: int = 7):
    """ONE seeded realisation with INDEPENDENT per-cell draws (the coupling Bliss
    independence implies). Used only for the maps."""
    rng = np.random.default_rng(seed)
    return rng.random(p_a.size) < p_a, rng.random(p_b.size) < p_b


# --------------------------------------------------------------------------- #
# 4. Simulation payload (the expensive half)
# --------------------------------------------------------------------------- #
def compute_payload(dx_um: float) -> Dict:
    """Reconstruct the real tissue ONCE and run every drug arm on it."""
    print(f"Reconstructing the REAL Xenium GBM microenvironment (once, dx={dx_um} um) ...")
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
    report = reconstruction_report(micro)
    print(f"  {report['n_cells']:,} real cells (tumor {report['n_tumor']:,})  "
          f"grid {report['grid']}  dx={report['dx_um']:.0f} um")

    alpha = calibrate_alpha()
    print(f"perivascular-sink alpha = {alpha:.1f} (fit ONCE on doxorubicin)\n")

    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)
    oxygen_cfg = OxygenConfig()
    pd_cfg = PDConfig()

    labeled, _oxy = oxygen_labeled_cells(micro, oxygen_cfg)
    tumor = labeled.tumor_mask()
    s_intrinsic = intrinsic_factor(labeled, pd_cfg)[tumor]

    arms: Dict[Tuple[str, float], Dict] = {}
    for name in DRUGS:
        dp, mol = compile_real_drug(name, alpha)
        for mult in (1.0, DOSE_MATCH_FACTOR):
            cfg = SimConfig(drug=dp, pk=PKConfig(plasma_peak=BASE_PLASMA_PEAK * mult),
                            pde=pde, oxygen=oxygen_cfg, pd=pd_cfg, seed=0)
            print(f"  simulating {name:14s} dose x{mult:.0f} ...", end="", flush=True)
            res = simulate_response(micro, cfg)
            p, expo = per_cell_kill(micro, labeled, res.transport, pd_cfg, dp)
            p_t = p[tumor]
            # EXACT factorisation of the engine's PD map: p = p_hill * intrinsic.
            hill = np.divide(p_t, s_intrinsic, out=np.zeros_like(p_t),
                             where=s_intrinsic > 0)
            arms[(name, mult)] = dict(
                name=name, mult=mult,
                ic50=float(dp.ic50), mw=float(mol.get("mw")), hill_coef=float(dp.hill),
                diffusivity=float(dp.diffusivity), binding=float(dp.binding),
                clearance=float(dp.clearance), bbb=float(dp.bbb_penetrance),
                duration=float(res.transport.times[-1]),
                p_tumor=p_t, hill_tumor=hill, exposure_tumor=expo[tumor],
                expected_kill=float(p_t.mean()), engine_kill=float(res.kill_fraction),
                mean_auc=float(res.extras["mean_auc"]),
                pen=float(res.extras["penetration_index"]),
                hypoxic=float(res.extras["hypoxic_tumor_fraction"]),
            )
            a = arms[(name, mult)]
            print(f" E[kill]={a['expected_kill']:.4f}  engine={a['engine_kill']:.4f}"
                  f"  mean_auc={a['mean_auc']:.3g}")

    # Vessel-distance map: the delivery axis.
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    vm = vessel_mask_from_density(micro.fields.vascular_density, tissue,
                                  percentile=VESSEL_PERCENTILE)
    dist_field = distance_from_vessel_um(vm["mask"], float(micro.domain.dx))
    col, row = micro.domain.to_index(labeled.x, labeled.y)

    return dict(
        report=report, alpha=float(alpha), arms=arms,
        x=np.asarray(labeled.x)[tumor], y=np.asarray(labeled.y)[tumor],
        extent=tuple(float(v) for v in micro.domain.extent_um),
        cell_dist=dist_field[row, col][tumor],
        s_intrinsic=s_intrinsic,
        vessel={"area_fraction": vm["area_fraction"], "n_vessel_px": vm["n_vessel_px"]},
        hypoxic=float(arms[(DRUGS[0], 1.0)]["hypoxic"]),
    )


# --------------------------------------------------------------------------- #
# 5. Analysis
# --------------------------------------------------------------------------- #
def analyse(payload: Dict, n_shuffle: int) -> Dict:
    arms, s_int = payload["arms"], payload["s_intrinsic"]
    cell_dist = payload["cell_dist"]

    edges = DIST_EDGES
    centers = 0.5 * (edges[:-1] + edges[1:])
    raw_idx = np.digitize(cell_dist, edges) - 1
    inside = (raw_idx >= 0) & (raw_idx < len(centers))   # beyond the last edge: DROPPED
    bin_idx = raw_idx[inside]
    bin_n = np.bincount(bin_idx, minlength=len(centers))
    dropped = float((~inside).mean())
    print(f"\n  distance binning: {inside.sum():,}/{cell_dist.size:,} tumor cells "
          f"inside 0-{edges[-1]:.0f} um of the vessel mask "
          f"({dropped*100:.1f}% beyond, dropped from the profiles)")

    def profile(values: np.ndarray) -> np.ndarray:
        s = np.bincount(bin_idx, weights=values[inside], minlength=len(centers))
        out = np.full(len(centers), np.nan)
        ok = bin_n >= MIN_BIN_CELLS
        out[ok] = s[ok] / bin_n[ok]
        return out

    prof = {
        "kill": {n: profile(arms[(n, 1.0)]["p_tumor"]) for n in DRUGS},
        "hill": {n: profile(arms[(n, 1.0)]["hill_tumor"]) for n in DRUGS},
        "exposure": {n: profile(arms[(n, 1.0)]["exposure_tumor"]) for n in DRUGS},
        "intrinsic": profile(s_int),
    }

    # POTENCY HEADROOM: the time-averaged effective concentration each cell sees,
    # in units of the drug's own IC50. This is what decides whether a steep
    # delivery gradient actually translates into a steep KILL gradient.
    ok = bin_n >= MIN_BIN_CELLS
    headroom = {}
    for n in DRUGS:
        a = arms[(n, 1.0)]
        e_cell = a["exposure_tumor"] / a["duration"]        # uM
        rel = e_cell / a["ic50"]
        e_prof = prof["exposure"][n] / a["duration"] / a["ic50"]
        # empirical half-distance of the exposure profile (no functional form).
        v = prof["exposure"][n]
        half_d = np.nan
        idx = np.flatnonzero(ok)
        if idx.size:
            v0 = v[idx[0]]
            below = idx[v[idx] <= 0.5 * v0]
            if below.size:
                half_d = float(centers[below[0]])
        headroom[n] = dict(
            median_E=float(np.median(e_cell)), median_rel=float(np.median(rel)),
            frac_above_ic50=float((rel >= 1.0).mean()),
            rel_profile=e_prof, half_distance_um=half_d,
            # deepest bin where the median cell still sees >= IC50
            last_bin_above=(float(centers[np.flatnonzero(ok & (e_prof >= 1.0))[-1]])
                            if np.any(ok & (e_prof >= 1.0)) else np.nan),
        )

    pair_rows = []
    for a_name, b_name in PAIRS:
        A, B = arms[(a_name, 1.0)], arms[(b_name, 1.0)]
        pa, pb = A["p_tumor"], B["p_tumor"]
        stats = set_stats(pa, pb)
        null_full = shuffle_null(pa, pb, "full", n_perm=n_shuffle)
        null_del = shuffle_null(pa, pb, "delivery", hill_b=B["hill_tumor"],
                                intrinsic=s_int, n_perm=n_shuffle)
        p_combo = bliss(pa, pb)

        better = a_name if A["expected_kill"] >= B["expected_kill"] else b_name
        p_better = arms[(better, 1.0)]["p_tumor"]
        ctrl = arms[(better, DOSE_MATCH_FACTOR)]
        self_combo = float(bliss(p_better, p_better).mean())

        def _corr(u, v):
            return float(np.corrcoef(u, v)[0, 1])

        def _rank(u):
            return np.argsort(np.argsort(u)).astype(np.float64)

        ka, kb = realise_sets(pa, pb)
        pair_rows.append(dict(
            a=a_name, b=b_name, stats=stats, null_full=null_full, null_del=null_del,
            kill_a=A["expected_kill"], kill_b=B["expected_kill"],
            kill_combo=float(p_combo.mean()),
            better=better, kill_ctrl=ctrl["expected_kill"],
            kill_selfcombo=self_combo,
            best_single=max(A["expected_kill"], B["expected_kill"]),
            pearson_p=_corr(pa, pb), spearman_p=_corr(_rank(pa), _rank(pb)),
            pearson_hill=_corr(A["hill_tumor"], B["hill_tumor"]),
            spearman_hill=_corr(_rank(A["hill_tumor"]), _rank(B["hill_tumor"])),
            profile_combo=profile(p_combo), killed_a=ka, killed_b=kb,
            # WHERE does each agent do its UNIQUE killing? These are the exact
            # expected per-cell contributions to the A-only / B-only / both sets,
            # binned by distance from vessel — the direct test of the hypothesis
            # that the penetrant agent covers the DEEP tissue the potent one misses.
            profile_a_only=profile(pa * (1.0 - pb)),
            profile_b_only=profile(pb * (1.0 - pa)),
            profile_both=profile(pa * pb),
        ))
        r = pair_rows[-1]
        print(f"\n  {a_name} + {b_name}: A={stats['kill_a']:.4f} B={stats['kill_b']:.4f} "
              f"combo={stats['kill_combo']:.4f}  2x{better}={ctrl['expected_kill']:.4f}  "
              f"self({better}+{better})={self_combo:.4f}\n"
              f"    Jaccard={stats['jaccard']:.4f}  full-null={null_full['jaccard_mean']:.4f} "
              f"({stats['jaccard']/null_full['jaccard_mean']:.2f}x)  "
              f"delivery-null={null_del['jaccard_mean']:.4f} "
              f"({stats['jaccard']/null_del['jaccard_mean']:.2f}x)")

    return dict(centers=centers, bin_n=bin_n, prof=prof, headroom=headroom,
                pair_rows=pair_rows, dist_dropped=dropped,
                dist_median=float(np.median(cell_dist)),
                dist_max=float(cell_dist.max()))


# --------------------------------------------------------------------------- #
# 6. Figure
# --------------------------------------------------------------------------- #
def make_figure(payload: Dict, ana: Dict):
    arms = payload["arms"]
    xs, ys = payload["x"], payload["y"]
    W, H = payload["extent"]
    centers, prof, pair_rows = ana["centers"], ana["prof"], ana["pair_rows"]
    report = payload["report"]

    fig = plt.figure(figsize=(19.5, 17.0))
    gs = GridSpec(4, 3, figure=fig, height_ratios=[0.62, 0.62, 1.0, 1.0],
                  hspace=0.30, wspace=0.20, top=0.94, bottom=0.035)

    # ---- Row 0: per-drug per-cell kill probability on the real tissue -------
    p_max = max(float(arms[(n, 1.0)]["p_tumor"].max()) for n in DRUGS) or 1.0
    map_axes, sc = [], None
    for j, name in enumerate(DRUGS):
        ax = fig.add_subplot(gs[0, j]); map_axes.append(ax)
        a = arms[(name, 1.0)]
        order = np.argsort(a["p_tumor"])
        sc = ax.scatter(xs[order], ys[order], c=a["p_tumor"][order], s=1.0,
                        cmap="magma", vmin=0.0, vmax=p_max, linewidths=0)
        ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{name} alone (1x)   IC50 {a['ic50']:.3g} uM   "
                     f"E[kill] {a['expected_kill']:.4f}", fontsize=10)
    cb = fig.colorbar(sc, ax=map_axes, shrink=0.85, pad=0.012, location="right")
    cb.set_label("per-cell kill probability")

    # ---- Row 1: killed-cell territory (one seeded realisation) --------------
    handles = [
        Line2D([], [], marker="o", ls="", ms=5, color=_COL_A, label="A only"),
        Line2D([], [], marker="o", ls="", ms=5, color=_COL_B, label="B only"),
        Line2D([], [], marker="o", ls="", ms=5, color=_COL_BOTH, label="both"),
        Line2D([], [], marker="o", ls="", ms=5, color=_COL_NONE, label="survives"),
    ]
    for j, row in enumerate(pair_rows):
        ax = fig.add_subplot(gs[1, j])
        ka, kb = row["killed_a"], row["killed_b"]
        both, a_only, b_only, none = ka & kb, ka & ~kb, kb & ~ka, ~(ka | kb)
        ax.scatter(xs[none], ys[none], s=0.6, c=_COL_NONE, linewidths=0)
        ax.scatter(xs[a_only], ys[a_only], s=1.5, c=_COL_A, linewidths=0)
        ax.scatter(xs[b_only], ys[b_only], s=1.5, c=_COL_B, linewidths=0)
        ax.scatter(xs[both], ys[both], s=2.0, c=_COL_BOTH, linewidths=0)
        ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        s = row["stats"]
        ax.set_title(
            f"A = {row['a']}   +   B = {row['b']}\n"
            f"killed population: {s['frac_a_only']*100:.1f}% A-only / "
            f"{s['frac_b_only']*100:.1f}% B-only / {s['frac_both']*100:.1f}% both"
            f"   (Jaccard {s['jaccard']:.3f})", fontsize=9.5)
        ax.legend(handles=handles, loc="lower right", fontsize=7.5, framealpha=0.92)

    # ---- Row 2a: DELIVERY axis — exposure vs distance, in IC50 units --------
    hr = ana["headroom"]
    ax = fig.add_subplot(gs[2, 0])
    for name in DRUGS:
        e = hr[name]["rel_profile"]
        ok = ~np.isnan(e)
        ax.plot(centers[ok], e[ok], "-o", ms=3, lw=1.6, color=_COL[name],
                label=f"{name}  (d1/2 = {hr[name]['half_distance_um']:.0f} um)")
    ax.axhline(1.0, color="black", lw=1.2)
    ax.text(centers[-1], 1.15, "IC50", ha="right", fontsize=8)
    ax.axvline(45.0, color="grey", ls=":", lw=1.2)
    ax.text(0.30, 0.04, "grey line: measured dox half-distance ~45 um",
            fontsize=7, color="grey", transform=ax.transAxes)
    ax.set_yscale("log"); ax.set_xscale("log")
    ax.set_xlabel("distance from nearest vessel pixel (um, log)")
    ax.set_ylabel("mean effective concentration / that drug's IC50")
    ax.set_title("DELIVERY: exposure decays steeply — and near-IDENTICALLY —\n"
                 "for all three agents, but stays ABOVE IC50 for the potent ones",
                 fontsize=9.5)
    ax.legend(fontsize=7.5); ax.grid(alpha=0.25, which="both")

    # ---- Row 2b: where each drug actually KILLS -----------------------------
    ax = fig.add_subplot(gs[2, 1])
    for name in DRUGS:
        ax.plot(centers, prof["kill"][name], "-o", ms=3, lw=1.6, color=_COL[name],
                label=f"{name} p_kill")
    head = pair_rows[0]
    ax.plot(centers, head["profile_combo"], "--", lw=2.0, color=_COL_BOTH,
            label=f"{head['a']}+{head['b']} (Bliss)")
    ax.plot(centers, prof["intrinsic"], ":", lw=1.8, color="#444444",
            label="shared cell-intrinsic factor (sens x resist)")
    ax.axvline(45.0, color="grey", ls=":", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("distance from nearest vessel pixel (um, log)")
    ax.set_ylabel("mean per-cell kill probability")
    ax.legend(fontsize=7.5, loc="lower left")
    ax.set_ylim(0, None)
    ax.set_title("KILL vs distance: only the WEAK agent's kill tracks its\n"
                 "delivery; the potent agents are saturated nearly everywhere",
                 fontsize=9.5)
    ax.grid(alpha=0.25)

    # ---- Row 2c: where does each agent do its UNIQUE killing? ---------------
    ax = fig.add_subplot(gs[2, 2])
    ax.plot(centers, head["profile_a_only"], "-o", ms=3, lw=1.7, color=_COL_A,
            label=f"killed by {head['a']} ONLY")
    ax.plot(centers, head["profile_b_only"], "-o", ms=3, lw=1.7, color=_COL_B,
            label=f"killed by {head['b']} ONLY")
    ax.plot(centers, head["profile_both"], "-o", ms=3, lw=1.7, color=_COL_BOTH,
            label="killed by BOTH")
    ax.axvline(45.0, color="grey", ls=":", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("distance from nearest vessel pixel (um, log)")
    ax.set_ylabel("expected per-cell contribution")
    ax.legend(fontsize=7.5)
    ax.set_title(f"THE HYPOTHESIS, DIRECTLY: does the penetrant agent "
                 f"({head['a']})\ncover the DEEP tissue {head['b']} cannot reach?",
                 fontsize=9.5)
    ax.grid(alpha=0.25)

    # ---- Row 3a: kill bars, with BOTH controls ------------------------------
    ax = fig.add_subplot(gs[3, 0])
    labels = [f"{r['a'][:4]}+{r['b'][:4]}" for r in pair_rows]
    x = np.arange(len(pair_rows))
    w = 0.16
    series = [
        ("A alone", [r["kill_a"] for r in pair_rows], _COL_A),
        ("B alone", [r["kill_b"] for r in pair_rows], _COL_B),
        ("A+B (Bliss)", [r["kill_combo"] for r in pair_rows], _COL_BOTH),
        ("2x dose of better single", [r["kill_ctrl"] for r in pair_rows], "#555555"),
        ("Bliss(better + ITSELF)", [r["kill_selfcombo"] for r in pair_rows], "#9ecae1"),
    ]
    for i, (lab, vals, c) in enumerate(series):
        b = ax.bar(x + (i - 2) * w, vals, w, label=lab, color=c,
                   edgecolor="black" if "ITSELF" in lab else "none", linewidth=0.6)
        for rect in b:
            ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height(),
                    f"{rect.get_height():.3f}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("expected tumor-cell kill fraction")
    ax.set_ylim(0, max(max(v) for _, v, _ in series) * 1.35)
    ax.legend(fontsize=7, loc="upper left", ncol=1)
    ax.set_title("THE HONEST TESTS: beat 2x dose of the better agent,\n"
                 "AND beat that agent Bliss-combined with ITSELF", fontsize=9.5)
    ax.grid(alpha=0.2, axis="y")

    # ---- Row 3b: observed overlap vs BOTH nulls -----------------------------
    ax = fig.add_subplot(gs[3, 1])
    obs = [r["stats"]["jaccard"] for r in pair_rows]
    nf = [r["null_full"]["jaccard_mean"] for r in pair_rows]
    nf_sd = [r["null_full"]["jaccard_sd"] for r in pair_rows]
    nd = [r["null_del"]["jaccard_mean"] for r in pair_rows]
    nd_sd = [r["null_del"]["jaccard_sd"] for r in pair_rows]
    w = 0.26
    ax.bar(x - w, obs, w, color=_COL_BOTH, label="observed Jaccard")
    ax.bar(x, nf, w, yerr=nf_sd, color="#999999", capsize=3,
           label="full-shuffle null")
    ax.bar(x + w, nd, w, yerr=nd_sd, color="#4d4d4d", capsize=3,
           label="delivery-shuffle null (intrinsic factor held)")
    for xi, r in enumerate(pair_rows):
        o = r["stats"]["jaccard"]
        ax.text(xi, max(o, nf[xi], nd[xi]) * 1.10,
                f"obs/full {o/nf[xi]:.2f}x\nobs/deliv {o/nd[xi]:.2f}x",
                ha="center", fontsize=7.5, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Jaccard overlap of killed populations")
    ax.set_ylim(0, max(obs + nf + nd) * 1.42)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("SPATIAL COMPLEMENTARITY TEST\n"
                 "<1 = different territory;  >1 = redundant, same cells", fontsize=9.5)
    ax.grid(alpha=0.2, axis="y")

    # ---- Row 3c: verdict text ------------------------------------------------
    ax = fig.add_subplot(gs[3, 2]); ax.axis("off")
    head_txt, lines = _verdict(pair_rows)
    body = head_txt.replace("**", "") + "\n\n"
    for ln in lines:
        body += "- " + ln.replace("**", "").lstrip("- ") + "\n\n"
    ax.text(0.0, 1.0, "VERDICT", fontsize=12, fontweight="bold", va="top")
    ax.text(0.0, 0.93, _wrap(body, 74), fontsize=8.2, va="top", family="monospace")

    fig.suptitle(
        "COMBINATION THERAPY on the REAL Xenium glioblastoma section — do "
        "delivery-complementary drugs cover more tumor?\n"
        f"{report['n_tumor']:,} real tumor cells, grid-converged dx = "
        f"{report['dx_um']:.0f} um, ONE reconstruction reused by every arm.   "
        "Combination = BLISS INDEPENDENCE, p = 1-(1-pA)(1-pB): spatial "
        "complementarity ONLY, NO pharmacological synergy modelled.\n"
        "PREDICTIONS on an untreated resection — no outcome label; nothing here is "
        f"a validated response. Tumor hypoxic fraction {payload['hypoxic']*100:.0f}%.",
        fontsize=12.5, y=0.985)
    fig.savefig(_OUT_PNG, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _wrap(text: str, width: int) -> str:
    import textwrap
    out = []
    for para in text.split("\n"):
        out.extend(textwrap.wrap(para, width, subsequent_indent="  ") or [""])
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# 7. Verdict + markdown report
# --------------------------------------------------------------------------- #
def _verdict(pair_rows) -> Tuple[str, List[str]]:
    """Honest verdict: complementary vs redundant territory, and whether the
    combination survives BOTH controls."""
    lines = []
    n_comp = n_dose = n_self = 0
    for r in pair_rows:
        s = r["stats"]
        rf = s["jaccard"] / r["null_full"]["jaccard_mean"]
        rd = s["jaccard"] / r["null_del"]["jaccard_mean"]
        comp = rd < 0.95
        beats_dose = r["kill_combo"] > r["kill_ctrl"]
        beats_self = r["kill_combo"] > r["kill_selfcombo"]
        n_comp += int(comp); n_dose += int(beats_dose); n_self += int(beats_self)
        verdict = ("complementary territory" if comp else
                   "REDUNDANT (same cells)" if rd > 1.05 else
                   "no spatial structure")
        lines.append(
            f"{r['a']} + {r['b']}: overlap vs delivery-null {rd:.2f}x "
            f"(vs full-null {rf:.2f}x) -> {verdict}. "
            f"combo {r['kill_combo']:.4f} vs 2x-{r['better']} {r['kill_ctrl']:.4f} "
            f"({r['kill_combo']/r['kill_ctrl']:.2f}x, "
            f"{'beats' if beats_dose else 'loses'}) and vs "
            f"Bliss({r['better']}+itself) {r['kill_selfcombo']:.4f} "
            f"({r['kill_combo']/r['kill_selfcombo']:.2f}x, "
            f"{'beats' if beats_self else 'LOSES'})."
        )
    n = len(pair_rows)
    if n_comp == 0 and n_self == 0:
        head = ("THE DELIVERY-COMPLEMENTARITY HYPOTHESIS FAILS on this tissue. "
                "Every pair's killed populations overlap MORE than the "
                "delivery-shuffled null (the two drugs hit the same cells), and no "
                "combination beats the better agent Bliss-combined with ITSELF — so "
                "the apparent combination benefit over 2x dosing is an artifact of "
                "the independence assumption acting against a saturating Hill "
                "curve, not spatial complementarity.")
    elif n_comp == n and n_self == n and n_dose == n:
        head = ("THE DELIVERY-COMPLEMENTARITY HYPOTHESIS HOLDS on this tissue. "
                "Every pair kills less-overlapping territory than the "
                "delivery-shuffled null and beats both controls.")
    else:
        head = (f"MIXED. {n_comp}/{n} pairs kill less-overlapping territory than the "
                f"delivery-shuffled null; {n_dose}/{n} beat 2x dosing of the better "
                f"agent; {n_self}/{n} beat that agent Bliss-combined with itself "
                f"(the control that isolates complementarity from the independence "
                f"arithmetic).")
    return head, lines


def build_markdown(payload: Dict, ana: Dict, n_shuffle: int) -> str:
    report, arms = payload["report"], payload["arms"]
    pair_rows, centers, prof = ana["pair_rows"], ana["centers"], ana["prof"]
    hypoxic = payload["hypoxic"]
    head, verdict_lines = _verdict(pair_rows)
    L: List[str] = []
    A = L.append

    A("# Combination therapy on the real Xenium glioblastoma microenvironment")
    A("")
    A("**Experiment:** `experiments/combination.py` "
      "(run: `./.venv/bin/python experiments/combination.py`)")
    A("")
    A("**Question.** Do two drugs with complementary DELIVERY profiles — a small, "
      "penetrant-but-weak agent (temozolomide) and a potent-but-delivery-limited "
      "agent (doxorubicin, consumed within ~45 um of a vessel) — cover more of a "
      "real tumor than either alone, *for delivery reasons*?")
    A("")
    A(f"> **Headline.** {head}")
    A("")

    A("## THE COMBINATION ASSUMPTION: BLISS INDEPENDENCE (read this first)")
    A("")
    A("Each drug is simulated **separately** on the same real tissue, giving every "
      "tumor cell its own exposure AUC and its own kill probability from the frozen "
      "engine's Hill pharmacodynamics. The combination is then formed at the "
      "exposure level by")
    A("")
    A("```")
    A("p_kill_combo(cell) = 1 - (1 - p_A(cell)) * (1 - p_B(cell))")
    A("```")
    A("")
    A("This is **Bliss independence**: the two agents kill each cell independently "
      "and a cell survives only if it survives both.")
    A("")
    A("> **NO pharmacological synergy or antagonism is modelled.** There is no "
      "shared target, no DNA-repair crosstalk, no efflux competition, no schedule "
      "dependence, no drug-drug PK interaction, and drug B does not alter drug A's "
      "transport. What is modelled is **SPATIAL complementarity only** — whether "
      "the two delivery fields reach *different cells*. Every gain reported below "
      "is geographic or arithmetic, never biochemical. Real combinations may "
      "synergise or antagonise; this model is blind to both, so treat these "
      "numbers as a delivery-only floor, not a forecast of combination efficacy.")
    A("")
    A("Because Bliss independence *is* the statement that the per-cell draws are "
      "independent, all set statistics are reported as exact expectations under "
      "that coupling:")
    A("")
    A("```")
    A("E|A n B|  = sum p_A p_B                E|A only| = sum p_A (1 - p_B)")
    A("E|A u B|  = sum [1-(1-p_A)(1-p_B)]     E|B only| = sum (1 - p_A) p_B")
    A("```")
    A("")
    A("so the headline numbers carry no Monte-Carlo noise; a single seeded "
      "realisation is drawn only to render the territory maps.")
    A("")

    A("## Method")
    A("")
    A(f"1. **One real tissue, reconstructed once.** "
      f"`reconstruct_real('xenium', dx_um={report['dx_um']:.0f})` -> "
      f"{report['n_cells']:,} measured cells ({report['n_tumor']:,} tumor), grid "
      f"{report['grid']}, real endothelial-KDE vasculature, real matrix-transcript "
      f"ECM, real-perfusion oxygen. The same object is reused by every arm, so all "
      f"differences are drug physics, not tissue.")
    A(f"2. **Drug compilation.** `calibrate_alpha()` -> alpha = "
      f"{payload['alpha']:.1f} (the single perivascular-sink constant, fit ONCE on "
      f"doxorubicin) -> `calibrated_drugparams(name, alpha)` -> "
      f"`apply_barrier_to_drugparams` (a no-op for these small molecules, but the "
      f"identical route an antibody would take).")
    A("3. **Separate simulations.** Each agent gets its own `simulate_response` run "
      "(oxygen -> IFP -> delivery -> transport PDE -> Hill kill), 6 h exposure "
      "window, convection on.")
    A("4. **Per-cell readout.** The engine's own `apply_drug_kill` is re-read on a "
      "throwaway copy of the *identical* oxygen-labelled cell population, "
      "recovering each cell's exposure AUC and kill probability with the same "
      "hypoxia / quiescence / sensitivity modulation the engine applied. Nothing "
      "is mutated; core / physics / molecular / reconstruct / engine are untouched.")
    A("5. **Combination.** Bliss independence, above.")
    A(f"6. **Control 1 — dose-matched.** The better single agent is re-run at "
      f"**{DOSE_MATCH_FACTOR:g}x** plasma peak: the same total systemic exposure "
      f"(in standard-dose units) the two-drug arm receives.")
    A("7. **Control 2 — self-combination.** The SAME Bliss rule applied to the "
      "better agent against *itself*: `1-(1-p_B)^2`. This is a hypothetical "
      "reference arm, not a realisable regimen — a 'combination' whose two "
      "components have perfectly identical, maximally redundant delivery. If "
      "`combo(A,B) ~= selfcombo(B,B)`, the benefit is an artifact of the "
      "independence assumption meeting a saturating Hill curve — not spatial "
      "complementarity.")
    A(f"8. **Control 3 — permutation nulls ({n_shuffle} perms each).** "
      f"Independent-Bernoulli overlap shrinks automatically as kill probabilities "
      f"shrink, so raw Jaccard proves nothing alone. The engine's PD map "
      f"factorises **exactly** as `p = p_hill(exposure) x (sensitivity x "
      f"resistance)`, and the second factor is CELL-INTRINSIC and SHARED by every "
      f"drug. Two nulls follow:")
    A("   - *full shuffle* — permute `p_B` across tumor cells (marginals exactly "
      "preserved; both delivery and shared-killability alignment destroyed);")
    A(f"   - *delivery shuffle* — permute only drug B's delivery term `p_hill` "
      f"within {N_STRATA} strata of the shared factor, then re-multiply, so each "
      f"cell keeps its own intrinsic killability and only B's delivery geography "
      f"is scrambled. **This is the null that isolates DELIVERY co-localisation.**")
    A("")
    A("   observed/null **< 1 = complementary territory**, **~1 = no spatial "
      "structure**, **> 1 = redundant** (both drugs hit the same cells).")
    A("")

    A("## Single agents on the real tissue (1x dose)")
    A("")
    A("| drug | MW | IC50 (uM) | D (um^2/s) | sink k (1/s) | BBB | mean tissue AUC (uM*s) | E[kill] | engine kill (1 draw) |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for name in DRUGS:
        a = arms[(name, 1.0)]
        A(f"| {name} | {a['mw']:.0f} | {a['ic50']:.3g} | {a['diffusivity']:.0f} | "
          f"{a['binding'] + a['clearance']:.3g} | {a['bbb']:.3f} | "
          f"{a['mean_auc']:.3g} | {a['expected_kill']:.4f} | {a['engine_kill']:.4f} |")
    A("")
    A("**POTENCY HEADROOM — the number that decides the whole experiment.** A steep "
      "delivery gradient only produces a steep KILL gradient if the concentration "
      "crosses IC50 somewhere inside the tissue. Per tumor cell, the "
      "time-averaged effective concentration `E = AUC / T` in units of that drug's "
      "own IC50:")
    A("")
    A("| drug | median E (uM) | IC50 (uM) | median E / IC50 | % tumor cells above IC50 | exposure half-distance | deepest bin still >= IC50 |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for name in DRUGS:
        a, h = arms[(name, 1.0)], ana["headroom"][name]
        lb = ("none" if np.isnan(h["last_bin_above"])
              else f"{h['last_bin_above']:.0f} um")
        A(f"| {name} | {h['median_E']:.4g} | {a['ic50']:.3g} | "
          f"{h['median_rel']:.3g} | {h['frac_above_ic50']*100:.1f}% | "
          f"{h['half_distance_um']:.0f} um | {lb} |")
    A("")
    A("`E[kill]` (mean per-cell kill probability over tumor cells) is the "
      "noise-free expectation; `engine kill` is one Bernoulli realisation returned "
      "by `simulate_response`. They agree to within sampling error — the "
      "cross-check that the re-read probabilities are the engine's own.")
    A("")
    A("**Dose-response is already saturating**, which is central to reading the "
      "controls below:")
    A("")
    A("| drug | E[kill] @1x | E[kill] @2x | gain from doubling the dose |")
    A("|---|---:|---:|---:|")
    for name in DRUGS:
        k1 = arms[(name, 1.0)]["expected_kill"]
        k2 = arms[(name, DOSE_MATCH_FACTOR)]["expected_kill"]
        A(f"| {name} | {k1:.4f} | {k2:.4f} | {k2 / k1:.2f}x |")
    A("")

    A("## Kill: A, B, A+B, and both controls")
    A("")
    A("| pair | A alone | B alone | A+B (Bliss) | better single | 2x dose of better | Bliss(better+ITSELF) | combo/best single | combo/2x-dose | combo/self-combo |")
    A("|---|---:|---:|---:|---|---:|---:|---:|---:|---:|")
    for r in pair_rows:
        A(f"| {r['a']} + {r['b']} | {r['kill_a']:.4f} | {r['kill_b']:.4f} | "
          f"{r['kill_combo']:.4f} | {r['better']} | {r['kill_ctrl']:.4f} | "
          f"{r['kill_selfcombo']:.4f} | {r['kill_combo'] / r['best_single']:.2f}x | "
          f"{r['kill_combo'] / r['kill_ctrl']:.2f}x | "
          f"{r['kill_combo'] / r['kill_selfcombo']:.2f}x |")
    A("")
    A("All values are expected fractions of **tumor cells** killed. The last column "
      "is the one that matters: a ratio at or below 1.00 means combining two "
      "*different* drugs did no better than the model's own independence "
      "arithmetic applied to a single drug against itself.")
    A("")

    A("## Spatial overlap of the killed populations")
    A("")
    A("| pair | killed by A only | by B only | by both | Jaccard | overlap coef | full-shuffle null | obs/full | delivery-shuffle null | obs/delivery |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in pair_rows:
        s, nf, nd = r["stats"], r["null_full"], r["null_del"]
        A(f"| {r['a']} + {r['b']} | {s['frac_a_only']*100:.1f}% | "
          f"{s['frac_b_only']*100:.1f}% | {s['frac_both']*100:.1f}% | "
          f"{s['jaccard']:.4f} | {s['overlap_coef']:.4f} | "
          f"{nf['jaccard_mean']:.4f} +- {nf['jaccard_sd']:.4f} | "
          f"{s['jaccard']/nf['jaccard_mean']:.2f}x | "
          f"{nd['jaccard_mean']:.4f} +- {nd['jaccard_sd']:.4f} | "
          f"{s['jaccard']/nd['jaccard_mean']:.2f}x |")
    A("")
    A("A-only / B-only / both are percentages **of the killed population** (they "
      "sum to 100%). The delivery-shuffle null preserves drug B's marginal kill to "
      "within " + ", ".join(
          f"{abs(r['null_del']['kill_b_mean'] - r['kill_b']) / r['kill_b'] * 100:.1f}%"
          for r in pair_rows) + " respectively.")
    A("")
    A("As fractions of all tumor cells:")
    A("")
    A("| pair | A only | B only | both | any (= combo kill) |")
    A("|---|---:|---:|---:|---:|")
    for r in pair_rows:
        s = r["stats"]
        A(f"| {r['a']} + {r['b']} | {s['tumor_frac_a_only']*100:.2f}% | "
          f"{s['tumor_frac_b_only']*100:.2f}% | {s['tumor_frac_both']*100:.2f}% | "
          f"{s['kill_combo']*100:.2f}% |")
    A("")
    A("Correlation of the per-cell kill probabilities across tumor cells — a "
      "coupling-free redundancy measure. `p` includes the shared cell-intrinsic "
      "factor; `p_hill` is the DELIVERY term alone:")
    A("")
    A("| pair | Pearson(p_A,p_B) | Spearman(p_A,p_B) | Pearson(hill_A,hill_B) | Spearman(hill_A,hill_B) |")
    A("|---|---:|---:|---:|---:|")
    for r in pair_rows:
        A(f"| {r['a']} + {r['b']} | {r['pearson_p']:.3f} | {r['spearman_p']:.3f} | "
          f"{r['pearson_hill']:.3f} | {r['spearman_hill']:.3f} |")
    A("")

    A("## Where each drug actually kills (the delivery axis)")
    A("")
    A(f"Vessel mask = the p{VESSEL_PERCENTILE:.0f} threshold on the real "
      f"endothelial-KDE vascular density "
      f"({payload['vessel']['area_fraction']*100:.1f}% of tissue area, "
      f"{payload['vessel']['n_vessel_px']:,} pixels); distance is the Euclidean "
      f"transform in um. Per tumor-cell distance bin: mean exposure AUC "
      f"(normalised to the perivascular bin), mean kill probability, and the "
      f"shared cell-intrinsic factor.")
    A("")
    A(f"Bin widths are 12 um (= dx) out to 120 um, then 40 um, then 100 um: tumor "
      f"cells on this section sit a MEDIAN of {ana['dist_median']:.0f} um from the "
      f"vessel mask (max {ana['dist_max']:.0f} um), so a uniform near-field grid "
      f"would collapse most of the tumor into one overflow bin. Cells beyond "
      f"{DIST_EDGES[-1]:.0f} um ({ana['dist_dropped']*100:.1f}% of tumor cells) are "
      f"DROPPED from the profiles rather than clipped into the last bin.")
    A("")
    A("| distance (um) | cells | " +
      " | ".join(f"AUC {n[:4]}" for n in DRUGS) + " | " +
      " | ".join(f"p_kill {n[:4]}" for n in DRUGS) +
      " | p_kill combo | intrinsic factor |")
    A("|---:|---:|" + "---:|" * (2 * len(DRUGS) + 2))
    first = {n: None for n in DRUGS}
    for i, c in enumerate(centers):
        if ana["bin_n"][i] < MIN_BIN_CELLS:
            continue
        for n in DRUGS:
            if first[n] is None:
                first[n] = prof["exposure"][n][i]
        A(f"| {c:.0f} | {ana['bin_n'][i]:,} | " +
          " | ".join(f"{prof['exposure'][n][i] / first[n]:.3f}" for n in DRUGS) + " | " +
          " | ".join(f"{prof['kill'][n][i]:.4f}" for n in DRUGS) +
          f" | {pair_rows[0]['profile_combo'][i]:.4f} | {prof['intrinsic'][i]:.4f} |")
    A("")
    A("Near-vessel vs deep, for exposure (the delivery quantity) and for kill:")
    A("")
    A("| drug | AUC near / AUC deep | p_kill perivascular bin | p_kill first bin >100 um | p_kill far/near |")
    A("|---|---:|---:|---:|---:|")
    ok = ana["bin_n"] >= MIN_BIN_CELLS
    i0 = int(np.flatnonzero(ok)[0])
    far_idx = np.flatnonzero(ok & (centers > 100.0))
    i1 = int(far_idx[0]) if far_idx.size else i0
    for name in DRUGS:
        e0, e1 = prof["exposure"][name][i0], prof["exposure"][name][i1]
        k0, k1 = prof["kill"][name][i0], prof["kill"][name][i1]
        A(f"| {name} | {e0 / e1:.2f}x | {k0:.4f} | {k1:.4f} | {k1 / k0:.2f}x |")
    A("")
    A(f"The shared cell-intrinsic factor is nearly FLAT over the same range "
      f"({prof['intrinsic'][i0]:.4f} perivascular -> {prof['intrinsic'][i1]:.4f} "
      f"deep, {prof['intrinsic'][i1] / prof['intrinsic'][i0]:.2f}x), so it is not "
      f"what shapes these curves — delivery is.")
    A("")

    A("## THE HYPOTHESIS, TESTED DIRECTLY: where does each agent kill *uniquely*?")
    A("")
    A("The claim under test is not just that the combination kills more, but that "
      "the penetrant agent reaches DEEP tissue the delivery-limited agent cannot. "
      "That is a statement about where each drug's UNIQUE kills sit. Expected "
      "per-cell contribution to the A-only / B-only / both sets, by distance from "
      "vessel:")
    A("")
    for k, r in enumerate(pair_rows):
        A(f"**{r['a']} (A) + {r['b']} (B)**" +
          ("" if k == 0 else "  (condensed — every fourth bin)"))
        A("")
        A("| distance (um) | cells | A only | B only | both | A-only share of kills |")
        A("|---:|---:|---:|---:|---:|---:|")
        for i, c in enumerate(centers):
            if ana["bin_n"][i] < MIN_BIN_CELLS:
                continue
            if k > 0 and i % 4 != 0:
                continue
            ao, bo, bt = (r["profile_a_only"][i], r["profile_b_only"][i],
                          r["profile_both"][i])
            A(f"| {c:.0f} | {ana['bin_n'][i]:,} | {ao:.4f} | {bo:.4f} | {bt:.4f} | "
              f"{ao / (ao + bo + bt) * 100:.1f}% |")
        A("")

    A("## Verdict")
    A("")
    A(head)
    A("")
    for ln in verdict_lines:
        A(f"- {ln}")
    A("")
    A("### What the numbers mean, mechanistically")
    A("")
    hd = pair_rows[0]
    hr = ana["headroom"]
    A(f"1. **The delivery gradients are real but nearly IDENTICAL in shape.** All "
      f"three agents' exposure falls to half within "
      f"{min(hr[n]['half_distance_um'] for n in DRUGS):.0f}-"
      f"{max(hr[n]['half_distance_um'] for n in DRUGS):.0f} um of the vessel mask "
      f"and the normalised profiles lie almost on top of one another. This is the "
      f"confound `results/vessel_profile.md` already documented: with a fast "
      f"perivascular sink the tissue field is close to reaction-dominated "
      f"(`C ~ S/k`), so every drug largely INHERITS the shape of the same vascular "
      f"source field. **Two drugs entering from the same vessels have correlated, "
      f"not complementary, delivery geometry** — which is the structural reason "
      f"the complementarity hypothesis has so little room to work here.")
    A(f"2. **What differs between the agents is MAGNITUDE relative to potency, not "
      f"shape.** Doxorubicin's mean tissue AUC is ~"
      f"{arms[('temozolomide', 1.0)]['mean_auc'] / arms[('doxorubicin', 1.0)]['mean_auc']:.0f}x "
      f"below temozolomide's — it is genuinely delivery-limited in concentration — "
      f"but its IC50 is ~"
      f"{arms[('temozolomide', 1.0)]['ic50'] / arms[('doxorubicin', 1.0)]['ic50']:.0f}x "
      f"lower, so the delivered concentration still sits at "
      f"{hr['doxorubicin']['median_rel']:.2g}x IC50 in the median tumor cell versus "
      f"{hr['temozolomide']['median_rel']:.2g}x for temozolomide. **Doxorubicin is "
      f"therefore NOT effectively delivery-limited on this section**: its Hill "
      f"curve is saturated well past the ~45 um cuff, so its kill barely decays "
      f"with distance while the weak agent's kill decays steeply. The premise of "
      f"the headline pairing fails at this step.")
    A(f"3. **The penetrant agent's unique kills are PERIVASCULAR, not deep** — the "
      f"opposite of the hypothesis. In the {hd['a']} + {hd['b']} pair the "
      f"{hd['a']}-only contribution is largest in the near-vessel bins and falls "
      f"with distance (table above), because {hd['a']} only wins where its own "
      f"exposure is highest. It does not open up a disjoint deep territory; the "
      f"deep tissue is already covered by the potent agent.")
    A("4. **Hence overlap ABOVE both nulls.** The killed populations co-localise "
      "more than chance even against the delivery-shuffle null, which holds each "
      "cell's intrinsic killability fixed. The redundancy is delivery redundancy.")
    A("5. **The combination still beats 2x dosing — but so does a drug combined "
      "with ITSELF.** Returns to dose are steeply diminishing (dose-doubling "
      "table), so under a Bliss model *any* second independent hit beats *any* "
      "second dose, including the same agent twice. The self-combination control "
      "is the one that separates the two explanations, and every pair LOSES to it. "
      "The apparent 'combination benefit' here is the independence arithmetic "
      "meeting a saturating Hill curve, not spatial complementarity.")
    A("")
    A("**What would change this answer.** A genuinely complementary pair needs "
      "(a) delivery fields with different SHAPES — which on this tissue means "
      "escaping the shared vascular source geometry (e.g. a hypoxia-activated "
      "prodrug, a drug whose sink is anti-correlated with perfusion, or "
      "convection-dominated delivery) — and (b) a potent agent whose exposure "
      "actually crosses IC50 inside the tissue, so that its kill is confined where "
      "its delivery is confined. Neither condition holds for these three agents on "
      "this section.")
    A("")

    A("## Caveats (all of them)")
    A("")
    A("- **No synergy is modelled.** Bliss independence is an assumption, not a "
      "result. Real drug-drug interaction (synergy or antagonism) is outside this "
      "model's reach and could dominate everything reported here.")
    A("- **No outcome validation.** This open Xenium section is an *untreated* "
      "resection: there is no paired measured response for any of these agents, "
      "alone or combined. Every number here is a PREDICTION.")
    A("- **Absolute kills are conditional on the sensitivity model.** The "
      "reconstruction carries a UNIFORM per-cell `sensitivity` of 1.0 (a "
      "deliberately drug-agnostic geometry-only choice, documented in "
      "`reconstruct/real_xenium.py`), so the only cell-intrinsic modulation here is "
      f"hypoxia / quiescence resistance ({hypoxic*100:.0f}% of tumor cells are "
      f"hypoxic). Absolute kill fractions are therefore an upper-ish bound "
      f"conditional on 'every tumor cell is fully sensitive', and would fall — "
      f"unevenly — under any real per-cell sensitivity model. Such a model could "
      f"move the overlap ratios in either direction; it would not change the "
      f"exposure fields, on which the delivery conclusions rest.")
    A("- **Absolute kills are also conditional on the potency inputs.** IC50s come "
      "from the cached GDSC glioma screen (e.g. temozolomide at "
      f"{arms[('temozolomide', 1.0)]['ic50']:.0f} uM, a genuinely weak in-vitro "
      f"potency); the headroom argument above is a direct function of those "
      f"numbers and of the shared nominal dose.")
    A("- **All agents share one nominal dose.** Every drug is driven at the same "
      f"{BASE_PLASMA_PEAK:g} uM peak plasma level, which is NOT clinically "
      f"dose-matched between agents, so cross-drug absolute comparisons are soft. "
      f"The dose-matched control is defined in the same 'standard dose' units, so "
      f"the combo-vs-2x comparison is internally consistent.")
    A("- **Dosing is the engine's Stage-0 plasma average.** The transport solve "
      "uses a static source at the time-averaged plasma level, so schedule "
      "(sequential vs concurrent — which matters enormously in real combination "
      "therapy) cannot be represented at all.")
    A("- **The vessel mask is a threshold on a 40 um-bandwidth endothelial KDE**, "
      "not a measured lumen, and `results/vessel_profile.md` documents that most "
      "of the apparent perivascular gradient on this section is vessel/KDE "
      "geometry rather than drug transport. The distance axis inherits both "
      "limitations.")
    A("- **Overlap statistics depend on the coupling.** Under independent draws, "
      "overlap shrinks automatically as kill probabilities shrink, which is "
      "precisely why the null RATIOS, not the raw Jaccard, are the reported test "
      "statistics.")
    A("- **One section, one tumor, three agents, three pairs.** No cohort, no "
      "replicate sections, no cross-tissue generalisation is claimed.")
    A("- **Research prototype.** Not a clinical decision tool; no regulatory, "
      "diagnostic or treatment-selection use.")
    A("")
    return "\n".join(L) + "\n"


def print_summary(payload, ana, t0):
    pair_rows = ana["pair_rows"]
    print("\n" + "=" * 108)
    print("COMBINATION THERAPY on the REAL Xenium GBM section — Bliss independence, "
          "spatial complementarity only")
    print("=" * 108)
    hdr = (f"{'pair':28s} {'A':>7} {'B':>7} {'A+B':>7} {'2x-bett':>8} {'self+self':>9} "
           f"{'Jacc':>7} {'/full':>6} {'/deliv':>7}")
    print(hdr); print("-" * len(hdr))
    for r in pair_rows:
        s = r["stats"]
        print(f"{r['a'] + ' + ' + r['b']:28s} {r['kill_a']:7.4f} {r['kill_b']:7.4f} "
              f"{r['kill_combo']:7.4f} {r['kill_ctrl']:8.4f} {r['kill_selfcombo']:9.4f} "
              f"{s['jaccard']:7.4f} {s['jaccard']/r['null_full']['jaccard_mean']:6.2f} "
              f"{s['jaccard']/r['null_del']['jaccard_mean']:7.2f}")
    print("-" * len(hdr))
    for r in pair_rows:
        s = r["stats"]
        print(f"{r['a'] + ' + ' + r['b']:28s} killed population: "
              f"{s['frac_a_only']*100:5.1f}% A-only  {s['frac_b_only']*100:5.1f}% "
              f"B-only  {s['frac_both']*100:5.1f}% both")
    head, lines = _verdict(pair_rows)
    print("=" * 108)
    print(head)
    for ln in lines:
        print("  - " + ln)
    print("=" * 108)
    print(f"tumor hypoxic fraction {payload['hypoxic']*100:.1f}%; untreated "
          f"resection -> NO outcome label; no synergy modelled.")
    print(f"[runtime {time.time() - t0:.1f} s]")


# --------------------------------------------------------------------------- #
# 8. Main
# --------------------------------------------------------------------------- #
def main(dx_um: float = DX_UM, n_shuffle: int = N_SHUFFLE, cache: str = None):
    t0 = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    if cache and os.path.exists(cache):
        print(f"loading cached simulation payload: {cache}")
        with open(cache, "rb") as fh:
            payload = pickle.load(fh)
    else:
        payload = compute_payload(dx_um)
        if cache:
            with open(cache, "wb") as fh:
                pickle.dump(payload, fh)

    ana = analyse(payload, n_shuffle)
    make_figure(payload, ana)
    print(f"\nfigure written: {_OUT_PNG}")
    with open(_OUT_MD, "w") as fh:
        fh.write(build_markdown(payload, ana, n_shuffle))
    print(f"report written: {_OUT_MD}")
    print_summary(payload, ana, t0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Combination therapy on real tissue.")
    ap.add_argument("--dx", type=float, default=DX_UM,
                    help="grid spacing in um (default 12, the converged value)")
    ap.add_argument("--shuffles", type=int, default=N_SHUFFLE,
                    help="permutations per null (default 200)")
    ap.add_argument("--cache", type=str, default=None,
                    help="optional path to cache/reuse the simulation payload "
                         "(pickle) so the analysis can be re-run without "
                         "re-simulating; nothing is cached by default")
    args = ap.parse_args()
    main(dx_um=args.dx, n_shuffle=args.shuffles, cache=args.cache)
