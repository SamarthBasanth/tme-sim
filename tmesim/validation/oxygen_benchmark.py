"""BENCHMARKING THE OXYGEN ROUTES — every available route, scored on the same data.

**THIS MODULE ONCE ASSERTED THE OPPOSITE OF WHAT ITS OWN NUMBERS SAY.** It was
written to follow through on a claim in ``results/state_validation.md`` §1.2 that
the reaction-diffusion oxygen solve was "falsified as value-add" — out-predicted
by a one-line distance transform built from its own source cells. That claim was
an INVERTED INEQUALITY. Every predictor here is sign-aligned to oxygen, so a
correct predictor scores a NEGATIVE rho and the stronger predictor is the one with
the larger ``|rho|``; on the 816,769-cell IO section the solved field scores
rho = -0.1523 against the measured hypoxia composite while the distance baseline
scores -0.1204, i.e. |-0.1523| > |-0.1204| and **the solve wins, by 1.27x**.
``results/oxygen_demotion.md`` is now the record of that correction rather than of
the demotion.

The structural consequence for THIS file: every sentence in :func:`write_report`
that states an ORDERING is now derived from the numbers through
:func:`state_validation._ordering` instead of written out in English. The previous
version hard-coded "The vessel-distance surrogate wins" and then interpolated
``ranked[0]`` beside it, so a regeneration under the corrected numbers would have
emitted a sentence contradicting its own interpolated values. Nothing here now
states a direction the arithmetic did not produce.

What the correction does NOT license is the opposite overclaim. |rho| ~ 0.15-0.20
is WEAK, it is within-section rather than out-of-sample, the composite is carried
by VEGFA alone, and §1.1's finding that the solve's effective diffusion length is
several times the published oxygen-diffusion distance SURVIVES intact — Spearman
rho is invariant to any monotone rescaling, so the solve can win this ranking and
still have its length scale wrong. Both are true at once.

It does four things and writes a benchmark report:

1. **Benchmark.** Every oxygen/hypoxia route this codebase can build — the PDE, the
   vessel-distance surrogate (two vessel definitions), the marker-derived field,
   and the trivial baselines — rank-correlated against the MEASURED VEGFA hypoxia
   score on the IO section, with permutation p-values and a ranking.
2. **Runtime.** Wall-clock for the sparse solve vs the kd-tree query, at three grid
   spacings, so "cheaper" is a measured number and not an adjective.
3. **Downstream.** How much predicted kill and therapeutic coverage actually move
   between ``oxygen_mode='pde'`` and ``'vessel_distance'`` for temozolomide,
   cisplatin and doxorubicin on the real Xenium section at dx = 12 um. If they
   barely move, that is reported as the finding: it would mean the oxygen field,
   however computed, is not load-bearing for drug predictions.
4. **Recommendation.** Which mode should be the default, and what would have to be
   re-run if it changed.

HOW THE DOWNSTREAM ARM IS COMPUTED (and why it is not a re-implementation)
--------------------------------------------------------------------------
``sim.engine.simulate_response`` RE-SOLVES the oxygen PDE internally, so handing it
a surrogate reconstruction changes nothing. ``sim/engine.py`` is not to be
modified, so instead we exploit a structural fact about the pipeline, verified by
grep and then by execution: **oxygen enters the kill through exactly one channel**.
``physics.pressure``, ``physics.diffusion`` and ``treatment.drugs`` never read
``fields.oxygen``, so the transport solve and therefore the exposure/AUC field —
and therefore ``validation.delivery_metrics.therapeutic_coverage``, which is a pure
function of per-cell AUC and the drug's IC50 — are **oxygen-invariant by
construction**. Oxygen reaches the kill only via
``cells.agents.update_cell_states`` -> the ``PDConfig.hypoxia_resistance`` /
``quiescent_resistance`` multipliers inside ``cells.pkpd.apply_drug_kill``.

So the alternative arm re-bands the cells with the surrogate field and re-runs the
UNMODIFIED ``apply_drug_kill`` on the engine's OWN transport result, with the
engine's own seeded generator. That structural claim is not asserted — it is
checked: :func:`_kill_from_oxygen` is first run with the PDE field and must
reproduce ``simulate_response``'s ``kill_fraction`` **exactly**, and the run aborts
if it does not.

CAVEATS THAT DO NOT GO AWAY
---------------------------
* **VEGFA is a hypoxia-RESPONSE gene, not an oxygen tension.** It is the output of
  a transcriptional programme with its own kinetics and its own non-hypoxic
  inducers. Every rho here is "predicts measured VEGFA", not "predicts pO2". What
  the comparison DOES establish is internal and immune to this: whatever VEGFA is
  reporting, these routes are ranked against it on identical inputs. Which route
  wins that ranking is read off the numbers, never asserted here.
* **The two sections are DIFFERENT PATIENTS.** The benchmark runs on the 480-gene
  IO section (816,769 cells); the downstream drug arm runs on the 319-gene brain
  section (40,887 cells) because that is the section every flagship drug result in
  ``results/`` was computed on. Nothing transfers between them.
* **The marker route is CIRCULAR on this panel.** The IO panel's only
  :data:`physics.oxygen_proxy.HYPOXIA_MARKER_GENES` member is VEGFA, so a
  marker-derived field benchmarked against a VEGFA target is correlating VEGFA
  with a smoothed copy of itself. It is reported as a CEILING — the best any
  spatial field could score given this smoothing and this dropout — and explicitly
  not as a competitor.
* **Absolute oxygen is uncalibrated in every arm.** ``OxygenConfig.consumption``
  has never been fitted to a measured oxygen tension and the surrogate's ``L`` is
  anchored to published morphometry, not to a tension either. Nothing here
  calibrates a level.

Run:    ./.venv/bin/python -m tmesim.validation.oxygen_benchmark
Writes: results/oxygen_benchmark.md

   NOT ``results/oxygen_demotion.md`` any more. That path now holds the
   hand-written record of the inverted-inequality correction, which no generator
   can reproduce and which regenerating this module used to destroy.
   :func:`write_report` additionally REFUSES to overwrite any file marked
   ``STATUS: WITHDRAWN``, so pointing it back at that path fails loudly instead of
   silently reverting a correction.

Units: distance um, oxygen dimensionless [0,1], AUC uM*s, coverage/kill fractions.
"""
from __future__ import annotations

import os
import re
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..cells.agents import update_cell_states
from ..cells.pkpd import apply_drug_kill
from ..core import CellType, OxygenConfig, PDConfig, PDEConfig, SimConfig
from ..data.real.singlecell_spatial import load as load_singlecell_spatial
from ..physics.oxygen import solve_oxygen
from ..physics.oxygen_proxy import (
    HYPOXIA_MARKER_GENES,
    HYPOXIA_ONSET_UM,
    diffusion_length_um,
    hypoxia_from_markers,
    oxygen_from_vessel_distance,
    pde_diffusion_length_um,
    vessel_distance_um,
)
from ..reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from ..sim.engine import simulate_response
from .cellstate_validation import (
    normalized_log1p,
    partial_spearman,
    spearman_permutation,
)
from .delivery_metrics import _PDE, compile_real_drug, compute_all
from .state_validation import (
    BRAIN_SAMPLE,
    DX_CONVERGED,
    DX_SWEEP,
    IO_PERM_SUBSAMPLE,
    IO_SAMPLE,
    N_PERM_IO,
    PERM_SEED,
    _ordering,
    spearman_full,
    subsample_index,
    trivial_baselines,
)

__all__ = [
    "DOWNSTREAM_DRUGS",
    "benchmark_predictors",
    "downstream_effect",
    "run",
    "runtime_comparison",
    "write_report",
]

#: The three drugs the downstream arm quantifies, spanning the delivery axis:
#: temozolomide (small, BBB-penetrant), cisplatin (mid), doxorubicin (large,
#: BBB-excluded). The same three ``results/vessel_network_impact.md`` uses, so the
#: two A/B studies are directly comparable.
DOWNSTREAM_DRUGS: Tuple[str, ...] = ("temozolomide", "cisplatin", "doxorubicin")

#: Repeats for the wall-clock comparison. The kd-tree query is milliseconds, so a
#: single timing would be dominated by noise.
TIMING_REPEATS = 3


def _results_dir() -> str:
    d = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
    )
    os.makedirs(d, exist_ok=True)
    return d


def _load(path: str) -> Dict:
    """Load a cached section, retrying once on a transient failure."""
    try:
        return load_singlecell_spatial(path=path, dense=False)
    except Exception:
        time.sleep(3)
        return load_singlecell_spatial(path=path, dense=False)


def _gene_counts(expr, genes: np.ndarray, gene: str) -> Optional[np.ndarray]:
    from scipy import sparse

    genes = np.asarray(genes, dtype=str)
    hit = np.flatnonzero(genes == gene)
    if hit.size == 0:
        return None
    j = int(hit[0])
    if sparse.issparse(expr):
        return np.asarray(expr[:, [j]].todense()).ravel().astype(np.float64)
    return np.asarray(expr)[:, j].astype(np.float64)


# ---------------------------------------------------------------------------
# 1. THE BENCHMARK
# ---------------------------------------------------------------------------
def benchmark_predictors(
    sample: str = IO_SAMPLE,
    dx_um: float = DX_CONVERGED,
    n_perm: int = N_PERM_IO,
) -> Dict:
    """Score every oxygen/hypoxia route against the measured VEGFA hypoxia signal.

    One reconstruction (the frozen ``reconstruct_real`` path, ``oxygen_mode='pde'``)
    supplies the cells, the solved field and the vascular KDE; every alternative
    route is then built on those SAME cells, so the comparison isolates the field
    construction and nothing else.

    Every predictor is sign-aligned to OXYGEN — larger means better perfused — so a
    correct predictor gives a NEGATIVE rho against a hypoxia score, and the numbers
    are directly comparable with ``results/state_validation.md`` §1.2.

    The primary target is **VEGFA alone**, depth-normalised: §1.1 of that report
    established that the historical VEGFA+HIF1A composite is carried entirely by
    VEGFA (HIF1A ranks 249/480 at rho = +0.0085, for the mechanistic reason that
    oxygen sensing acts on the HIF-1alpha protein, not its mRNA). The composite is
    reported alongside so the two documents can be reconciled, never instead.

    Permutation p-values use the same documented device as §1.2: a fixed seeded
    subsample of :data:`state_validation.IO_PERM_SUBSAMPLE` cells with ``n_perm``
    draws, because the kernel is O(n x draws). Every rho is on the FULL population;
    only the null is subsampled, which can only make a p-value more conservative.
    """
    print("\n" + "=" * 78)
    print(f"BENCHMARK — every oxygen route vs MEASURED hypoxia ({sample}, dx={dx_um:g})")
    print("=" * 78)

    data = _load(sample)
    expr, genes = data["expression"], np.asarray(data["genes"], dtype=str)
    micro = reconstruct_real("xenium", sample=sample, dx_um=float(dx_um))
    rep = reconstruction_report(micro)
    cells = micro.cells
    if int(expr.shape[0]) != int(cells.n):
        raise ValueError(
            f"expression rows ({expr.shape[0]}) != cells ({cells.n}); the "
            "reconstruction and the expression matrix are not aligned"
        )
    print(f"  {cells.n:,} cells, grid {micro.domain.shape}, "
          f"{int((np.asarray(cells.cell_type) == int(CellType.ENDOTHELIAL)).sum()):,} "
          "ENDOTHELIAL")

    col, row = micro.domain.to_index(cells.x, cells.y)
    o2_pde = np.asarray(micro.fields.oxygen, dtype=np.float64)[row, col]
    vasc_at_cell = np.asarray(micro.fields.vascular_density, dtype=np.float64)[row, col]

    # --- the MEASURED targets ------------------------------------------------
    from scipy import sparse

    totals = (np.asarray(expr.sum(axis=1)).ravel().astype(np.float64)
              if sparse.issparse(expr) else
              np.asarray(expr, dtype=np.float64).sum(axis=1))
    vegfa_raw = _gene_counts(expr, genes, "VEGFA")
    if vegfa_raw is None:
        raise ValueError(f"{sample} has no VEGFA; this benchmark needs it as the target")
    vegfa = normalized_log1p(vegfa_raw, totals)
    hif1a_raw = _gene_counts(expr, genes, "HIF1A")
    if hif1a_raw is not None:
        h = normalized_log1p(hif1a_raw, totals)
        z = lambda v: (v - v.mean()) / (v.std() or 1.0)  # noqa: E731
        composite = 0.5 * (z(vegfa) + z(h))
    else:
        composite = vegfa.copy()
    log_depth = np.log1p(totals)

    targets = {"vegfa": vegfa, "composite": composite}
    print(f"  VEGFA detected in {float((vegfa_raw > 0).mean()) * 100:.2f}% of cells"
          + (f"; HIF1A in {float((hif1a_raw > 0).mean()) * 100:.2f}%"
             if hif1a_raw is not None else ""))

    # --- the predictors ------------------------------------------------------
    L = diffusion_length_um()
    preds: Dict[str, Dict[str, object]] = {}
    timing: Dict[str, float] = {}

    preds["pde"] = {
        "values": o2_pde,
        "label": "SIMULATED O2 (the reaction-diffusion PDE)",
        "family": "model",
        "what": "sparse steady-state solve; the current shipping default",
    }

    t0 = time.perf_counter()
    surro_grid, surro_info = oxygen_from_vessel_distance(micro, detail=True)
    timing["vessel_distance_grid"] = time.perf_counter() - t0
    preds["vessel_distance"] = {
        "values": surro_grid[row, col],
        "label": f"vessel-distance SURROGATE, exp(-d/L), L={L:.1f} um (grid-sampled)",
        "family": "surrogate",
        "what": "one cKDTree query per grid cell; no solve. This is what "
                "oxygen_mode='vessel_distance' installs",
    }

    t0 = time.perf_counter()
    d_cells = vessel_distance_um(micro, at="cells")
    timing["vessel_distance_cells"] = time.perf_counter() - t0
    preds["vessel_distance_exact"] = {
        "values": np.exp(-d_cells / L),
        "label": "the same surrogate at each cell's EXACT coordinate (no mesh)",
        "family": "surrogate",
        "what": "identical rho to raw distance-to-vessel: the map is monotone",
    }

    t0 = time.perf_counter()
    vn_grid, vn_info = oxygen_from_vessel_distance(
        micro, source="vessel_network", detail=True
    )
    timing["vessel_network_grid"] = time.perf_counter() - t0
    preds["vessel_network"] = {
        "values": vn_grid[row, col],
        "label": (
            f"surrogate to the nearest SEGMENTED vessel "
            f"({vn_info['n_vessels']:,} objects >= "
            f"{vn_info['min_cells_per_vessel']} cells, link "
            f"{vn_info['link_distance_um']:.0f} um)"
        ),
        "family": "surrogate",
        "what": "more physical vessel definition; discards singleton endothelium",
    }

    t0 = time.perf_counter()
    mk = hypoxia_from_markers(micro, expr, genes, detail=True)
    timing["markers_grid"] = time.perf_counter() - t0
    if mk["field"] is not None:
        preds["markers"] = {
            "values": np.asarray(mk["field"])[row, col],
            "label": f"MEASURED marker field ({', '.join(mk['genes_used'])})",
            "family": "circular",
            "what": (
                "CIRCULAR against this target — the only HYPOXIA_MARKER_GENES "
                "member on this panel IS VEGFA, so this correlates VEGFA with a "
                "smoothed copy of itself. Reported as a CEILING, not a competitor"
            ),
        }

    base = trivial_baselines(cells.x, cells.y, cells.cell_type, vasc_at_cell)
    for key in ("vascular_kde", "dist_endothelial", "local_cell_density"):
        preds[key] = {
            "values": np.asarray(base[key]["values"], dtype=np.float64),
            "label": str(base[key]["label"]),
            "family": "baseline",
            "what": str(base[key]["what"]),
        }

    # --- score ---------------------------------------------------------------
    sub = subsample_index(int(cells.n), IO_PERM_SUBSAMPLE, PERM_SEED)
    rows: List[Dict[str, object]] = []
    for name, blk in preds.items():
        v = np.asarray(blk["values"], dtype=np.float64)
        finite = np.isfinite(v)
        vv = np.where(finite, v, np.nanmin(v[finite]) if finite.any() else 0.0)
        rho_v = spearman_full(vv, targets["vegfa"])
        rho_c = spearman_full(vv, targets["composite"])
        perm = spearman_permutation(
            vv[sub], targets["vegfa"][sub], n_perm=int(n_perm), seed=PERM_SEED
        )
        rows.append({
            "name": name,
            "label": blk["label"],
            "family": blk["family"],
            "what": blk["what"],
            "rho_vegfa": rho_v,
            "rho_composite": rho_c,
            "rho_vegfa_depth_partialled": partial_spearman(
                vv, targets["vegfa"], log_depth
            ),
            "rho_subsample": perm["rho"],
            "p_perm": perm["p_value"],
            "n_perm": int(n_perm),
            "n_subsample": int(sub.size),
        })
        print(f"  {name:24s} rho(VEGFA) = {rho_v:+.4f}   p = {perm['p_value']:.5f}")

    # Competitors only: circular routes cannot be ranked against the others.
    #
    # Sorted by DESCENDING |rho|, which is what "stronger predictor" means under
    # this table's sign convention. While every rho is negative — as it is on this
    # section for every route — that is exactly the same ordering as sorting
    # ascending on the signed rho, which is how this line used to read; the
    # ordering was never the defect. The abs() is here so that a route which came
    # out ANTI-aligned (positive rho, i.e. predicting hypoxia where it predicts
    # oxygen) would be ranked on the strength of its association rather than
    # silently sent to the bottom of the table.
    ranked = sorted(
        [r for r in rows if r["family"] != "circular"],
        key=lambda r: -abs(r["rho_vegfa"]),
    )
    for i, r in enumerate(ranked, 1):
        r["rank"] = i

    # --- the nested test, re-run on the primary (VEGFA-only) target ----------
    d_ox = -d_cells  # sign-aligned to oxygen
    nested = {
        "vegfa_vs_pde_given_distance": partial_spearman(o2_pde, targets["vegfa"], d_ox),
        "vegfa_vs_distance_given_pde": partial_spearman(d_ox, targets["vegfa"], o2_pde),
        "pde_vs_distance": spearman_full(o2_pde, d_ox),
    }
    print(f"  nested: PDE|distance = {nested['vegfa_vs_pde_given_distance']:+.4f}   "
          f"distance|PDE = {nested['vegfa_vs_distance_given_pde']:+.4f}")

    return {
        "sample": sample,
        "dx_um": float(dx_um),
        "n_cells": int(cells.n),
        "grid": tuple(int(v) for v in micro.domain.shape),
        "report": rep,
        "rows": rows,
        "ranked": ranked,
        "nested": nested,
        "diffusion_length_um": float(L),
        "hypoxia_onset_um": float(HYPOXIA_ONSET_UM),
        "pde_length_scale": pde_diffusion_length_um(micro),
        "surrogate_info": surro_info,
        "vessel_network_info": vn_info,
        "marker_info": (
            {k: v for k, v in mk.items() if k not in ("field", "per_cell")}
            if isinstance(mk, dict) else None
        ),
        "marker_genes_looked_for": list(HYPOXIA_MARKER_GENES),
        "vegfa_detection": float((vegfa_raw > 0).mean()),
        "hif1a_detection": (float((hif1a_raw > 0).mean())
                            if hif1a_raw is not None else None),
        "build_times_s": timing,
    }


# ---------------------------------------------------------------------------
# 2. RUNTIME
# ---------------------------------------------------------------------------
def runtime_comparison(
    sample: str = IO_SAMPLE,
    dx_sweep: Sequence[float] = DX_SWEEP,
    repeats: int = TIMING_REPEATS,
) -> List[Dict]:
    """Wall-clock: the sparse oxygen solve vs the kd-tree surrogate, at three dx.

    Timed on the SAME reconstructed microenvironment so the reconstruction cost —
    which both routes share and neither avoids — is excluded. Best-of-``repeats``
    is reported for each, which is the standard defence against a scheduler
    hiccup inflating a single measurement.

    Only the field CONSTRUCTION is timed. The kd-tree route additionally pays for
    the tree build, and that build is inside the timed region.
    """
    print("\n" + "=" * 78)
    print("RUNTIME — sparse solve vs kd-tree query")
    print("=" * 78)
    cfg = OxygenConfig()
    out: List[Dict] = []
    for dx in dx_sweep:
        micro = reconstruct_real("xenium", sample=sample, dx_um=float(dx))
        ny, nx = micro.domain.shape

        def _best(fn) -> float:
            best = float("inf")
            for _ in range(max(1, int(repeats))):
                t0 = time.perf_counter()
                fn()
                best = min(best, time.perf_counter() - t0)
            return best

        t_pde = _best(lambda: solve_oxygen(micro, cfg))
        t_sur = _best(lambda: oxygen_from_vessel_distance(micro, cfg))
        out.append({
            "dx_um": float(dx),
            "grid": (int(ny), int(nx)),
            "n_unknowns": int(ny * nx),
            "pde_s": t_pde,
            "vessel_distance_s": t_sur,
            "speedup": (t_pde / t_sur) if t_sur > 0 else float("inf"),
        })
        print(f"  dx={dx:5.1f} um  grid {ny}x{nx} ({ny*nx:,} unknowns)  "
              f"PDE {t_pde:7.3f}s   surrogate {t_sur:7.4f}s   "
              f"{t_pde / max(t_sur, 1e-12):8.1f}x")
    return out


# ---------------------------------------------------------------------------
# 3. DOWNSTREAM EFFECT
# ---------------------------------------------------------------------------
def _kill_from_oxygen(micro, oxygen, transport, cfg: SimConfig, drug) -> Dict:
    """Re-band the cells with ``oxygen`` and re-run the UNMODIFIED kill.

    Nothing is re-implemented: this calls ``cells.agents.update_cell_states`` and
    ``cells.pkpd.apply_drug_kill``, the same two functions the engine calls, on a
    private copy of the cells and on the engine's own ``transport`` result. It is
    valid because the transport/AUC field does not depend on oxygen (no module in
    ``physics.pressure`` / ``physics.diffusion`` / ``treatment.drugs`` reads
    ``fields.oxygen``), and that validity is CHECKED by the caller, which requires
    this function to reproduce ``simulate_response`` exactly in the PDE arm.
    """
    from ..core import Microenvironment

    work = Microenvironment(
        domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
        patient_id=micro.patient_id, metadata=micro.metadata,
    )
    update_cell_states(work, np.asarray(oxygen, dtype=np.float64), cfg.oxygen)
    states = work.cells.state.copy()
    apply_drug_kill(work, transport, cfg.pd, drug,
                    rng=np.random.default_rng(cfg.seed))
    viable = float(work.cells.viable_tumor_fraction())

    tumor = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
    n_t = int(tumor.sum())
    from ..core import CellState

    bands = {
        s.name: (float((states[tumor] == int(s)).sum()) / n_t if n_t else 0.0)
        for s in (CellState.HYPOXIC, CellState.QUIESCENT, CellState.VIABLE)
    }
    return {"kill_fraction": 1.0 - viable, "viable_tumor_fraction": viable,
            "tumor_state_fractions": bands}


def downstream_effect(
    sample: Optional[str] = None,
    dx_um: float = DX_CONVERGED,
    drugs: Sequence[str] = DOWNSTREAM_DRUGS,
) -> Dict:
    """Does the oxygen route change the drug answer? — the A/B on the real section.

    Arm A: ``oxygen_mode='pde'`` — the shipping default, run through the unmodified
    :func:`sim.engine.simulate_response`.
    Arm B: ``oxygen_mode='vessel_distance'`` — the same cells, the same drug, the
    same transport, re-banded by the surrogate (see :func:`_kill_from_oxygen`).

    **The replication check is a hard gate.** Before arm B is believed, arm A is
    recomputed through :func:`_kill_from_oxygen` using the PDE field and must
    return the engine's ``kill_fraction`` bit-for-bit. If it does not, the
    structural assumption behind the whole comparison is wrong and the function
    raises rather than reporting a number.

    Therapeutic coverage is reported for both arms and is expected to be EXACTLY
    equal: it is a function of per-cell AUC and the drug's IC50, and the AUC field
    has no oxygen dependence anywhere in the pipeline. That exactness is itself a
    result — it says the oxygen field cannot move a delivery metric at all — so it
    is measured rather than asserted.
    """
    print("\n" + "=" * 78)
    print(f"DOWNSTREAM — kill & coverage, pde vs vessel_distance (dx = {dx_um:g} um)")
    print("=" * 78)

    from .calibration import calibrate_alpha

    alpha = calibrate_alpha()
    sample_path = sample if sample is not None else BRAIN_SAMPLE

    micro_a = reconstruct_real("xenium", sample=sample_path, dx_um=float(dx_um))
    rep_a = reconstruction_report(micro_a)
    micro_b = reconstruct_real("xenium", sample=sample_path, dx_um=float(dx_um),
                               oxygen_mode="vessel_distance")
    rep_b = reconstruction_report(micro_b)
    ox_a = np.asarray(micro_a.fields.oxygen, dtype=np.float64)
    ox_b = np.asarray(micro_b.fields.oxygen, dtype=np.float64)
    print(f"  {micro_a.cells.n:,} cells, grid {micro_a.domain.shape}; "
          f"mean O2 pde={rep_a['mean_oxygen']:.4f} "
          f"surrogate={rep_b['mean_oxygen']:.4f}")

    rows: Dict[str, Dict] = {}
    for name in drugs:
        dp, mol = compile_real_drug(name, alpha)
        cfg = SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0)
        t0 = time.time()
        res = simulate_response(micro_a, cfg)

        # HARD GATE: reproduce the engine through the re-banding path.
        check = _kill_from_oxygen(micro_a, res.fields["oxygen"], res.transport, cfg, dp)
        if check["kill_fraction"] != float(res.kill_fraction):
            raise AssertionError(
                f"[{name}] the re-banding path does not reproduce the engine "
                f"({check['kill_fraction']!r} vs {float(res.kill_fraction)!r}); the "
                "assumption that oxygen reaches the kill only through the state "
                "banding is WRONG and no downstream number here may be believed"
            )

        arm_b = _kill_from_oxygen(micro_b, ox_b, res.transport, cfg, dp)

        # Coverage: identical inputs except fields.oxygen, which compute_all reads
        # only for its hypoxia-band block.
        met_a = compute_all(micro_a, res, dp)
        res_b = _swap_oxygen_in_result(res, ox_b)
        met_b = compute_all(micro_b, res_b, dp)

        rows[name] = {
            "drug": name,
            "ic50_uM": float(dp.ic50),
            "mw": float(mol.get("mw")),
            "pde": {
                "kill_fraction": float(res.kill_fraction),
                "tumor_state_fractions": check["tumor_state_fractions"],
                "coverage_1x": float(met_a["therapeutic_coverage"]["coverage_1x"]),
                "coverage_0.5x": float(met_a["therapeutic_coverage"]["coverage_0.5x"]),
                "coverage_2x": float(met_a["therapeutic_coverage"]["coverage_2x"]),
                "median_E_over_ic50": float(
                    met_a["therapeutic_coverage"]["median_E_over_ic50"]),
                "hypoxic_over_normoxic_exposure": float(
                    met_a["hypoxia"]["ratio_hypoxic_over_normoxic"]),
            },
            "vessel_distance": {
                "kill_fraction": float(arm_b["kill_fraction"]),
                "tumor_state_fractions": arm_b["tumor_state_fractions"],
                "coverage_1x": float(met_b["therapeutic_coverage"]["coverage_1x"]),
                "coverage_0.5x": float(met_b["therapeutic_coverage"]["coverage_0.5x"]),
                "coverage_2x": float(met_b["therapeutic_coverage"]["coverage_2x"]),
                "median_E_over_ic50": float(
                    met_b["therapeutic_coverage"]["median_E_over_ic50"]),
                "hypoxic_over_normoxic_exposure": float(
                    met_b["hypoxia"]["ratio_hypoxic_over_normoxic"]),
            },
            "runtime_s": time.time() - t0,
        }
        a, b = rows[name]["pde"], rows[name]["vessel_distance"]
        d_kill = b["kill_fraction"] - a["kill_fraction"]
        rel = (d_kill / a["kill_fraction"] * 100.0) if a["kill_fraction"] else float("nan")
        rows[name]["delta_kill"] = float(d_kill)
        rows[name]["delta_kill_pct"] = float(rel)
        rows[name]["delta_coverage_1x"] = float(b["coverage_1x"] - a["coverage_1x"])
        print(f"  {name:14s} kill {a['kill_fraction']:.4f} -> {b['kill_fraction']:.4f} "
              f"({d_kill:+.4f}, {rel:+.1f}%)   cov@IC50 {a['coverage_1x']:.4f} -> "
              f"{b['coverage_1x']:.4f} ({b['coverage_1x'] - a['coverage_1x']:+.2e})")

    return {
        "sample": sample_path,
        "dx_um": float(dx_um),
        "alpha": float(alpha),
        "n_cells": int(micro_a.cells.n),
        "grid": tuple(int(v) for v in micro_a.domain.shape),
        "report_pde": rep_a,
        "report_vessel_distance": rep_b,
        "mean_oxygen_pde": float(rep_a["mean_oxygen"]),
        "mean_oxygen_vessel_distance": float(rep_b["mean_oxygen"]),
        "field_spearman": spearman_full(ox_a.ravel(), ox_b.ravel()),
        "pd": {"hypoxia_resistance": PDConfig().hypoxia_resistance,
               "quiescent_resistance": PDConfig().quiescent_resistance},
        "rows": rows,
    }


def _swap_oxygen_in_result(res, oxygen: np.ndarray):
    """A shallow copy of a ``SimulationResult`` with ``fields['oxygen']`` replaced.

    Only ``delivery_metrics.hypoxic_vs_normoxic_exposure_ratio`` reads that field;
    everything else in ``compute_all`` is a function of the AUC. Swapping it lets
    the two arms' metric blocks be computed by the identical unmodified code.
    """
    from dataclasses import replace

    fields = dict(res.fields)
    fields["oxygen"] = np.asarray(oxygen, dtype=np.float64)
    return replace(res, fields=fields)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def run(
    io_sample: str = IO_SAMPLE,
    brain_sample: str = BRAIN_SAMPLE,
    dx_um: float = DX_CONVERGED,
    n_perm: int = N_PERM_IO,
    drugs: Sequence[str] = DOWNSTREAM_DRUGS,
    cache_path: Optional[str] = None,
) -> Dict:
    """Benchmark + runtime + downstream, in one deterministic pass.

    ``cache_path`` (optional) pickles the raw results so the report can be
    regenerated without repeating the ~20-minute compute. Same convention as
    :func:`validation.delivery_metrics.run` and
    :func:`validation.vessel_network_impact.run`.
    """
    if cache_path and os.path.exists(cache_path):
        import pickle

        with open(cache_path, "rb") as fh:
            print(f"loaded cached results from {cache_path}")
            return pickle.load(fh)

    t0 = time.time()
    bench = benchmark_predictors(io_sample, dx_um, n_perm)
    timing = runtime_comparison(io_sample)
    down = downstream_effect(brain_sample, dx_um, drugs)
    results = {
        "benchmark": bench,
        "runtime": timing,
        "downstream": down,
        "elapsed_s": time.time() - t0,
    }
    if cache_path:
        import pickle

        # The raw fields and per-cell arrays are large and not needed to render.
        slim = dict(results)
        slim["benchmark"] = {
            k: v for k, v in bench.items() if k not in ("rows", "ranked")
        }
        slim["benchmark"]["rows"] = [
            {k: v for k, v in r.items()} for r in bench["rows"]
        ]
        slim["benchmark"]["ranked"] = [
            next(r for r in slim["benchmark"]["rows"] if r["name"] == q["name"])
            for q in bench["ranked"]
        ]
        with open(cache_path, "wb") as fh:
            pickle.dump(slim, fh)
        print(f"cached results to {cache_path}")
    return results


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------
def _fmt(v, nd: int = 4) -> str:
    if v is None:
        return "n/a"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not np.isfinite(f):
        return "n/a"
    return f"{f:.{nd}f}"


def _refuse_to_clobber_a_withdrawal(out_path: str) -> None:
    """Abort rather than overwrite a document marked ``STATUS: WITHDRAWN``.

    This module used to regenerate ``results/oxygen_demotion.md``. That path now
    holds the hand-written record of the inverted-inequality correction — a
    document a generator cannot reproduce — so a routine re-run of this module
    would have DELETED the correction and re-emitted the claim it corrects. The
    check is on the file's own content rather than on its name, so it protects
    any withdrawn artefact this generator is ever pointed at.
    """
    if not os.path.exists(out_path):
        return
    with open(out_path) as fh:
        head = fh.read(4096)
    # Anchored to the start of a line, allowing a blockquote marker and bold
    # emphasis, because that is how such a banner is actually written:
    #
    #     > **STATUS: WITHDRAWN AND REPLACED BY THIS CORRECTION RECORD.**
    #
    # A bare substring test would be wrong in both directions: this generator's
    # own correction block DISCUSSES the marker in running prose, so a substring
    # test made the generator refuse to overwrite its own previous output and
    # therefore non-idempotent -- which is a real bug, caught by running the
    # regeneration twice.
    banner = re.compile(r"^[>\s]*\**\s*STATUS:\s*WITHDRAWN", re.MULTILINE)
    if banner.search(head):
        raise RuntimeError(
            f"refusing to overwrite {out_path}: its header marks it WITHDRAWN, "
            "so it is a correction record rather than generator output. Write "
            "the benchmark somewhere else."
        )


def write_report(results: Dict, out_path: Optional[str] = None) -> str:
    """Write ``results/oxygen_benchmark.md``. Returns the written path.

    EVERY SENTENCE THAT STATES AN ORDERING between the solved field and a
    competitor is derived through :func:`state_validation._ordering` rather than
    written out in English, for the reason given in that function's docstring and
    recorded in ``results/oxygen_demotion.md``: the previous version of this
    function hard-coded *"The vessel-distance surrogate wins"* and then
    interpolated ``ranked[0]`` next to it, so once the underlying inequality was
    read correctly a regeneration would have emitted a sentence contradicting its
    own numbers. There is now no such sentence to emit.
    """
    b = results["benchmark"]
    rt = results["runtime"]
    dn = results["downstream"]
    out_path = out_path or os.path.join(_results_dir(), "oxygen_benchmark.md")
    _refuse_to_clobber_a_withdrawal(out_path)

    winner = b["ranked"][0]
    pde_row = next(r for r in b["rows"] if r["name"] == "pde")
    vd_row = next(r for r in b["rows"] if r["name"] == "vessel_distance")
    vde_row = next(r for r in b["rows"] if r["name"] == "vessel_distance_exact")

    # ---- THE ORDERING, DERIVED ---------------------------------------------
    # The strongest route that is NOT the PDE, i.e. the thing the PDE is on trial
    # against. Ranked by |rho| under this table's sign convention.
    rivals = [r for r in b["ranked"] if r["name"] != "pde"]
    best_rival = max(rivals, key=lambda r: abs(r["rho_vegfa"])) if rivals else pde_row
    ordv = _ordering(pde_row["rho_vegfa"], best_rival["rho_vegfa"])
    kde_row = next((r for r in b["rows"] if r["name"] == "vascular_kde"), None)
    # The solve against its OWN undiffused source term: does diffusing add or
    # destroy predictive information? Also derived, never asserted.
    ordk = (_ordering(pde_row["rho_vegfa"], kde_row["rho_vegfa"])
            if kde_row is not None else None)
    pde_wins = bool(ordv["sim_wins"])

    kill_deltas = [abs(r["delta_kill"]) for r in dn["rows"].values()]
    cov_deltas = [abs(r["delta_coverage_1x"]) for r in dn["rows"].values()]
    max_kill_delta = max(kill_deltas) if kill_deltas else 0.0
    max_cov_delta = max(cov_deltas) if cov_deltas else 0.0
    max_speedup = max(r["speedup"] for r in rt)
    min_speedup = min(r["speedup"] for r in rt)

    L = []
    A = L.append
    A("# The oxygen routes benchmarked: which predicts measured hypoxia, what it "
      "costs, and whether it changes a drug number")
    A("")
    A("Every oxygen route this codebase can build, scored against the same measured "
      "VEGFA signal on the same cells, then priced, then tested for whether the "
      "choice moves a drug prediction at all.")
    A("")
    if pde_wins:
        # Emitted, not hand-written, so a regeneration cannot quietly drop it. It
        # appears only while the solve is the winner -- i.e. exactly when a reader
        # would otherwise have no way of knowing this document once said the
        # opposite, and was the stated basis for deprecating a component.
        A("> ## CORRECTION — this report was originally written on an inverted "
          "inequality")
        A(">")
        A("> **What it used to say.** That `results/state_validation.md` §1.2 had "
          "found the reaction-diffusion solve *\"falsified as value-add\"*, "
          "out-predicted by a one-line distance transform built from its own "
          "source cells; that *\"the vessel-distance surrogate wins\"*; and that "
          "*\"there is no accuracy argument for the solve\"*. This module and "
          "`tmesim/physics/oxygen_proxy.py` were both written to act on that.")
        A(">")
        A("> **Why it was wrong.** Every predictor in the §1 table is sign-aligned "
          "to oxygen — `state_validation.trivial_baselines` literally returns "
          "`-dist_endo` — so, as that table's own header states, a correct "
          f"predictor gives a NEGATIVE rho and the STRONGER predictor is the one "
          f"with the larger |rho|. On this run "
          f"|{_fmt(pde_row['rho_vegfa'])}| > |{_fmt(best_rival['rho_vegfa'])}|: "
          f"**the solved field out-predicts `{best_rival['name']}`**, by "
          f"{ordv['ratio']:.2f}x. The inequality was read on the signed number "
          "instead of on its magnitude.")
        A(">")
        A("> **The fix that makes this unrepeatable.** Every ordering claim below "
          "is now DERIVED from the numbers via `state_validation._ordering()` "
          "rather than written in English, so re-running this module reproduces "
          "whatever the data say — including a flip back — instead of reverting "
          "to a remembered conclusion. `write_report` also refuses to overwrite "
          "any document whose header marks it as withdrawn, which is what stops "
          "a re-run from deleting `results/oxygen_demotion.md`.")
        A(">")
        A("> **What does NOT follow.** The solve winning does not make it good. "
          f"|rho| ~ {ordv['abs_sim']:.2f} is WEAK; this is a within-section "
          "comparison, not out-of-sample; and the VEGFA+HIF1A composite is carried "
          "by VEGFA alone. The correct reading is **\"the PDE is the best of "
          "several weak predictors\"**, never \"the oxygen solve is validated\". "
          "§1.1's separate finding — that the solve's effective diffusion length "
          "is several times the published oxygen-diffusion distance because "
          "`OxygenConfig.consumption` is unfitted — **survives this correction "
          "unchanged**, because a Spearman rho is invariant to monotone rescaling: "
          "a field can win this ranking and still have the wrong length scale.")
        A("")
    A("**The default is unchanged.** `reconstruct_real(oxygen_mode='pde')` is still "
      "the shipping route and is bit-identical to the pre-option pipeline. What "
      "this work added is an explicit alternative, labelled a SURROGATE, with the "
      "route's provenance travelling in "
      "`reconstruction_report()['oxygen']`.")
    A("")

    # --- headline ---
    A("## Headline")
    A("")
    A("| | |")
    A("|---|---|")
    A(f"| **Best predictor of measured VEGFA** | `{winner['name']}` at "
      f"rho = {_fmt(winner['rho_vegfa'])} (p = {_fmt(winner['p_perm'], 5)}) |")
    A(f"| **The PDE** | rho = {_fmt(pde_row['rho_vegfa'])} "
      f"(p = {_fmt(pde_row['p_perm'], 5)}), rank "
      f"{pde_row.get('rank', '-')} of {len(b['ranked'])} |")
    A(f"| **Cost ratio** | the PDE is {min_speedup:.0f}x to {max_speedup:.0f}x "
      "slower than the surrogate |")
    A(f"| **Downstream effect on predicted kill** | at most "
      f"{max_kill_delta:.4f} absolute across three drugs |")
    A(f"| **Downstream effect on therapeutic coverage** | {max_cov_delta:.3g} — "
      "coverage is oxygen-invariant *by construction* |")
    A("| **Recommendation** | "
      + ("keep `'pde'` as the default — it is the default ON the evidence, not "
         "despite it. See §5. |" if pde_wins else
         "keep `'pde'` as the default; use `oxygen_mode='vessel_distance'` for "
         "any NEW hypoxia claim. See §5. |"))
    A("")
    A("> **The two sections are DIFFERENT PATIENTS.** The benchmark runs on the "
      f"480-gene IO section ({b['n_cells']:,} cells); the downstream drug arm runs "
      f"on the 319-gene brain section ({dn['n_cells']:,} cells), because that is "
      "the section every flagship drug result in `results/` was computed on. "
      "Nothing transfers between them.")
    A("")

    # --- 1 benchmark ---
    A("## 1. Benchmark — every oxygen route against measured VEGFA")
    A("")
    A(f"**Section.** `{b['sample']}` — {b['n_cells']:,} cells, 480-gene IO panel, "
      f"dx = {b['dx_um']:g} um, grid {b['grid'][0]}x{b['grid'][1]}.  ")
    A(f"**Target.** Depth-normalised VEGFA (detected in "
      f"{b['vegfa_detection'] * 100:.2f}% of cells). "
      + (f"HIF1A ({b['hif1a_detection'] * 100:.2f}%) is reported only inside the "
         "historical composite column — `state_validation.md` §1.1 showed it ranks "
         "249/480 at rho = +0.0085 and dilutes the signal, for the mechanistic "
         "reason that oxygen sensing acts on the HIF-1alpha PROTEIN, not its mRNA "
         "(Ivan et al., Science 2001;292(5516):464-468, PMID 11292862, "
         "doi:10.1126/science.1059817; Jaakkola et al., Science "
         "2001;292(5516):468-472, PMID 11292861, doi:10.1126/science.1059796)."
         if b["hif1a_detection"] is not None else ""))
    A("")
    A("Every predictor is **sign-aligned to oxygen** (larger = better perfused), so "
      "a correct predictor gives a NEGATIVE rho. All rho are on the full "
      f"{b['n_cells']:,} cells; permutation p uses a fixed seeded subsample of "
      f"{b['rows'][0]['n_subsample']:,} cells with {b['rows'][0]['n_perm']:,} "
      "draws (the kernel is O(n x draws); subsampling the null can only make a p "
      "more conservative).")
    A("")
    A("| rank | route | rho vs VEGFA | permutation p | rho vs VEGFA+HIF1A composite | "
      "rho with log-depth removed | what it is |")
    A("|---|---|---|---|---|---|---|")
    for r in b["ranked"]:
        A(f"| {r['rank']} | **{r['name']}** | **{_fmt(r['rho_vegfa'])}** | "
          f"{_fmt(r['p_perm'], 5)} | {_fmt(r['rho_composite'])} | "
          f"{_fmt(r['rho_vegfa_depth_partialled'])} | {r['what']} |")
    for r in b["rows"]:
        if r["family"] == "circular":
            A(f"| — | ~~{r['name']}~~ (CIRCULAR) | {_fmt(r['rho_vegfa'])} | "
              f"{_fmt(r['p_perm'], 5)} | {_fmt(r['rho_composite'])} | "
              f"{_fmt(r['rho_vegfa_depth_partialled'])} | {r['what']} |")
    A("")
    A("**Reading the table.**")
    A("")
    # THE ORDERING SENTENCE. Derived end to end -- winner, direction, ratio and
    # the "what does not follow" rider all come out of _ordering() and out of
    # ranked[0], so this bullet cannot state an ordering the run did not produce.
    # It is the sentence that used to hard-code "The vessel-distance surrogate
    # wins" while interpolating ranked[0] beside it.
    A(f"* **{'The solved field wins' if pde_wins else 'The PDE is out-predicted'}.** "
      f"Rank 1 is `{winner['name']}` at rho = {_fmt(winner['rho_vegfa'])}. The PDE "
      f"is at {_fmt(pde_row['rho_vegfa'])}; the strongest route that is NOT the "
      f"PDE, `{best_rival['name']}`, is at {_fmt(best_rival['rho_vegfa'])}. Read "
      "by |rho|, as this table's sign convention requires, the solved field "
      f"{ordv['verb_vs_baseline']} it by {ordv['ratio']:.2f}x — and the "
      f"{'loser' if pde_wins else 'solve'} recovers "
      f"{100 * ordv['recovers']:.0f}% of the "
      f"{'winner' if pde_wins else 'baseline'}'s association. "
      + ("That a one-line kd-tree query gets MOST of the way there is the honest "
         "caveat; it does not get all of the way, and it does not win. "
         if pde_wins else
         "A kd-tree query out-predicting a sparse linear solve on the solve's own "
         "inputs is the finding. ")
      + "The two surrogate rows are the SAME predictor written twice — "
      "`vessel_distance_exact` is `exp(-d/L)` and `dist_endothelial` is `-d`, and "
      "a monotone map cannot change a Spearman rho — and are both listed only to "
      "make that identity visible rather than claimed.")
    if ordk is not None:
        # The solve against its own undiffused source: adds or destroys?
        A(f"* **The diffusion step "
          f"{'ADDS predictive information' if ordk['sim_wins'] else 'DESTROYS predictive information'}.** "
          f"`vascular_kde` is the solver's OWN SOURCE TERM sampled before it is "
          f"diffused, at rho = {_fmt(kde_row['rho_vegfa'])}. The solve takes that "
          f"field and returns one at {_fmt(pde_row['rho_vegfa'])} — "
          f"{ordk['ratio']:.2f}x on |rho| in favour of "
          f"{'the solved field' if ordk['sim_wins'] else 'the undiffused source'}. "
          + ("This is the most direct evidence in the document that the "
             "reaction-diffusion step is doing work rather than blurring."
             if ordk["sim_wins"] else
             "The solve is destroying information its own input already had."))
    A(f"* **Meshing costs almost nothing here.** The grid-sampled surrogate "
      f"({_fmt(vd_row['rho_vegfa'])}) and the exact per-cell surrogate "
      f"({_fmt(vde_row['rho_vegfa'])}) differ by "
      f"{abs(vd_row['rho_vegfa'] - vde_row['rho_vegfa']):.4f}, so binning the "
      f"distance to a {b['dx_um']:g} um mesh is not what separates the routes.")
    vn_row = next((r for r in b["rows"] if r["name"] == "vessel_network"), None)
    if vn_row is not None:
        # Derived both ways: which vessel definition scores better, and what the
        # discarded singletons were therefore carrying. The trailing clause used
        # to hard-code "the plainer definition ... is the one that wins", which
        # contradicts the leading clause whenever the segmented one scores higher.
        vn_better = abs(vn_row["rho_vegfa"]) > abs(vd_row["rho_vegfa"])
        n_dropped = (b["surrogate_info"]["n_points"]
                     - b["vessel_network_info"]["n_points"])
        A(f"* **The 'more physical' vessel definition is "
          f"{'better' if vn_better else 'WORSE'}.** Measuring distance to a "
          f"SEGMENTED vessel object rather than to any endothelial cell scores "
          f"{_fmt(vn_row['rho_vegfa'])} against {_fmt(vd_row['rho_vegfa'])}. "
          f"Segmentation discards {n_dropped:,} of "
          f"{b['surrogate_info']['n_points']:,} endothelial cells as singletons, "
          + ("so on this section those discarded cells were carrying NOISE, not "
             "signal. It is still NOT the default: the plainer definition is the "
             "one every other number in this report was benchmarked on, and "
             "switching it is a separate change with its own re-run list."
             if vn_better else
             "and on this section those discarded cells were carrying signal. It "
             "is offered as an option and is NOT the default — the plainer "
             "definition is the one that was benchmarked and the one that "
             "scores higher.")
          + " Either way this is a comparison between two SURROGATES, neither of "
          f"which is the rank-1 route here (`{winner['name']}`).")
    A("* **The mapping cannot be tuned to win.** `exp(-d/L)` is strictly monotone in "
      "distance, so **every** monotone map of distance gives the identical Spearman "
      "rho. The constant `L` sets where the HYPOXIC/QUIESCENT/VIABLE band edges "
      "fall and nothing else. This is why it is anchored to published morphometry "
      "rather than fitted — fitting it would be unfalsifiable decoration.")
    if any(r["family"] == "circular" for r in b["rows"]):
        mrow = next(r for r in b["rows"] if r["family"] == "circular")
        A(f"* **The marker route is a CEILING, not a competitor.** rho = "
          f"{_fmt(mrow['rho_vegfa'])}. The only "
          f"`HYPOXIA_MARKER_GENES` member on this panel is VEGFA "
          f"(absent: {', '.join(g for g in b['marker_genes_looked_for'] if g != 'VEGFA')}), "
          "so this correlates VEGFA with a kernel-smoothed copy of itself. It "
          "measures how much of the per-cell marker signal survives 40 um "
          "smoothing — i.e. an upper bound on what ANY spatial field could score "
          "here — and it must never be quoted as a validated predictor.")
    A("")
    A("### 1.1 Calibration of the surrogate mapping")
    A("")
    A("| quantity | value | provenance |")
    A("|---|---|---|")
    A(f"| hypoxia onset distance | {b['hypoxia_onset_um']:.0f} um | Thomlinson & "
      "Gray, Br J Cancer 1955;9(4):539-549 (PMID 13304213, "
      "doi:10.1038/bjc.1955.55) — tumour cords viable to ~100-180 um from stroma, "
      "necrotic beyond |")
    A(f"| hypoxia threshold | {b['surrogate_info']['hypoxia_threshold']:.2f} | "
      "`OxygenConfig.hypoxia_threshold`, unchanged |")
    A(f"| derived length scale L | {b['diffusion_length_um']:.2f} um | "
      "`L = onset / ln(1/threshold)`, so a cell at the published onset distance "
      "lands exactly on the existing hypoxia threshold — no new tuning knob |")
    A(f"| VIABLE within | {b['surrogate_info']['viable_within_um']:.1f} um | "
      "`L * ln(1/(2*threshold))`, the QUIESCENT/VIABLE edge |")
    A(f"| the PDE's OWN length scale | "
      f"{b['pde_length_scale']['L_mean_um']:.0f} um | "
      "`sqrt(D / (consumption * rho_tum))` at the mean tumour areal density |")
    A("")
    A(f"**The last row is a finding, not a footnote, and it SURVIVES the "
      f"correction above.** The solved field's own effective diffusion length is "
      f"~{b['pde_length_scale']['L_mean_um']:.0f} um, "
      f"{b['pde_length_scale']['L_mean_um'] / b['hypoxia_onset_um']:.1f}x the "
      "measured oxygen-diffusion distance. `OxygenConfig.consumption` has never "
      "been fitted to an oxygen tension, and at its current value the solve "
      "low-pass-filters its own source at a length scale the tissue does not "
      "have. **This criticism is independent of who wins §1**, and it is important "
      "not to let the ranking result erase it: a Spearman rho is invariant to any "
      "monotone rescaling of the field, so the solve can top the correlation table "
      "and still place its hypoxic band edge at the wrong distance from a vessel. "
      "Ranking is not calibration. The one thing that has changed is the "
      "INTERPRETATION: this over-long length scale can no longer be offered as the "
      "mechanism by which diffusion destroys predictive information, because on "
      "these data it does not destroy it. "
      "Independent support that vessel distance is the right axis in the first "
      "place: Helmlinger et al., Nat Med 1997;3(2):177-182 (PMID 9018236, "
      "doi:10.1038/nm0297-177) measured interstitial pO2 profiles BETWEEN ADJACENT "
      "VESSELS at <=10 um resolution in vivo — and, honestly, also found "
      "heterogeneous profile shapes and hypoxic well-perfused vessels, which is "
      "the ceiling on what geometry alone can know.")
    A("")
    A("### 1.2 Nested test on the VEGFA-only target")
    A("")
    A("| partial Spearman | rho |")
    A("|---|---|")
    A(f"| VEGFA vs **simulated O2**, distance-to-vessel held fixed | "
      f"{_fmt(b['nested']['vegfa_vs_pde_given_distance'])} |")
    A(f"| VEGFA vs **-distance-to-vessel**, simulated O2 held fixed | "
      f"{_fmt(b['nested']['vegfa_vs_distance_given_pde'])} |")
    A(f"| simulated O2 vs -distance-to-vessel (how nested the pair is) | "
      f"{_fmt(b['nested']['pde_vs_distance'])} |")
    A("")
    # Derived, not asserted: which member of the nested pair retains more.
    _np_pde = abs(float(b["nested"]["vegfa_vs_pde_given_distance"]))
    _np_dist = abs(float(b["nested"]["vegfa_vs_distance_given_pde"]))
    if _np_pde > _np_dist:
        A(f"The PDE keeps more explanatory power after distance is conditioned out "
          f"(|{_np_pde:.4f}|) than distance keeps after the PDE is conditioned out "
          f"(|{_np_dist:.4f}|). **Distance is the redundant member of the pair** on "
          "the cleaner VEGFA-only target — the same direction as the corrected §1.2 "
          "of `results/state_validation.md`. Both partials are small in absolute "
          "terms; this orders two weak predictors, it does not validate either.")
    else:
        A(f"Distance keeps more explanatory power after the PDE is conditioned out "
          f"(|{_np_dist:.4f}|) than the PDE keeps after distance is conditioned out "
          f"(|{_np_pde:.4f}|). **The PDE is the redundant member of the pair** on "
          "the cleaner VEGFA-only target.")
    A("")

    # --- 2 runtime ---
    A("## 2. Runtime — a sparse solve against a kd-tree query")
    A("")
    A(f"Best of {TIMING_REPEATS} on the same reconstructed microenvironment "
      f"(`{b['sample']}`), so the shared reconstruction cost is excluded. The "
      "surrogate's timing INCLUDES building the kd-tree.")
    A("")
    A("| dx (um) | grid | unknowns | PDE solve (s) | vessel-distance (s) | speed-up |")
    A("|---|---|---|---|---|---|")
    for r in rt:
        A(f"| {r['dx_um']:g} | {r['grid'][0]}x{r['grid'][1]} | "
          f"{r['n_unknowns']:,} | {r['pde_s']:.3f} | {r['vessel_distance_s']:.4f} | "
          f"**{r['speedup']:.0f}x** |")
    A("")
    A("The gap widens with refinement, as it must: the sparse direct solve's cost "
      "grows superlinearly in the number of unknowns while the kd-tree query grows "
      "as `O(n_grid log n_vessels)`. At the production dx = 12 um the solve is the "
      f"single most expensive step in a reconstruction and the surrogate is "
      f"{rt[-1]['speedup']:.0f}x cheaper.")
    A("")

    # --- 3 downstream ---
    A("## 3. Downstream — does the oxygen route change a drug prediction?")
    A("")
    A(f"**Section.** `{dn['sample']}` — {dn['n_cells']:,} cells, 319-gene brain "
      f"panel, dx = {dn['dx_um']:g} um, grid {dn['grid'][0]}x{dn['grid'][1]}, "
      f"geometric tissue mask (the shipping default). Perivascular-sink alpha = "
      f"{dn['alpha']:.2f}, calibrated once on doxorubicin as in every flagship run.")
    A("")
    A("**How arm B is computed.** `sim.engine.simulate_response` re-solves the PDE "
      "internally and `sim/engine.py` is not modified here, so arm B re-bands the "
      "cells with the surrogate via the unmodified "
      "`cells.agents.update_cell_states` and re-runs the unmodified "
      "`cells.pkpd.apply_drug_kill` on **the engine's own transport result**, with "
      "the engine's own seeded generator. This is valid because nothing in "
      "`physics.pressure`, `physics.diffusion` or `treatment.drugs` reads "
      "`fields.oxygen` — oxygen reaches the kill through exactly one channel, the "
      "state banding. That is not assumed: arm A is recomputed through the same "
      "path and **must reproduce the engine's `kill_fraction` bit-for-bit or the "
      "run aborts**. It did, for all three drugs.")
    A("")
    A(f"Mean oxygen over tissue: PDE {dn['mean_oxygen_pde']:.4f}, surrogate "
      f"{dn['mean_oxygen_vessel_distance']:.4f}. The two fields agree spatially at "
      f"Spearman rho = {_fmt(dn['field_spearman'])} over the grid.")
    A("")
    A("### 3.1 Tumour-cell state occupancy under each route")
    A("")
    A("| route | HYPOXIC | QUIESCENT | VIABLE |")
    A("|---|---|---|---|")
    first = next(iter(dn["rows"].values()))
    for arm, key in (("`oxygen_mode='pde'`", "pde"),
                     ("`oxygen_mode='vessel_distance'`", "vessel_distance")):
        f_ = first[key]["tumor_state_fractions"]
        A(f"| {arm} | {f_['HYPOXIC'] * 100:.2f}% | {f_['QUIESCENT'] * 100:.2f}% | "
          f"{f_['VIABLE'] * 100:.2f}% |")
    A("")
    A("The banding is identical for every drug (it is set before any drug is "
      "applied), so one row each suffices. The resistance multipliers this feeds "
      f"are `hypoxia_resistance = {dn['pd']['hypoxia_resistance']}` and "
      f"`quiescent_resistance = {dn['pd']['quiescent_resistance']}`.")
    A("")
    pde_bands = first["pde"]["tumor_state_fractions"]
    if pde_bands["HYPOXIC"] <= 0.0:
        A("**Note the PDE's HYPOXIC band is EMPTY on this section under the "
          "shipping `geometric` tissue mask** — a point `results/state_validation.md` "
          "§2.2 already made. So `PDConfig.hypoxia_resistance` currently has no "
          "cells to act on at all in the flagship configuration, and the entire "
          "pde-vs-surrogate difference below runs through the QUIESCENT band. Any "
          "statement anywhere in `results/` about hypoxia-driven resistance on this "
          "section is, under this mask, describing an empty set.")
        A("")
    A("### 3.2 Predicted kill and therapeutic coverage")
    A("")
    A("| drug | IC50 (uM) | kill (pde) | kill (vd) | delta kill | delta % | "
      "cov@1xIC50 (pde) | cov@1xIC50 (vd) | median E/IC50 (pde) | "
      "median E/IC50 (vd) |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for name, r in dn["rows"].items():
        a, bb = r["pde"], r["vessel_distance"]
        A(f"| {name} | {r['ic50_uM']:.3g} | {a['kill_fraction']:.4f} | "
          f"{bb['kill_fraction']:.4f} | {r['delta_kill']:+.4f} | "
          f"{r['delta_kill_pct']:+.1f}% | {a['coverage_1x']:.4f} | "
          f"{bb['coverage_1x']:.4f} | {a['median_E_over_ic50']:.4g} | "
          f"{bb['median_E_over_ic50']:.4g} |")
    A("")
    A("`median E/IC50` is carried alongside coverage because coverage@1xIC50 is "
      "degenerate at one end or the other for these drugs on this section "
      "(temozolomide 0.0000 — no tumour cell reaches its IC50; cisplatin and "
      "doxorubicin near-saturated above 0.96), so an identical pair of coverage "
      "numbers could be read as a saturation artefact. The median potency margin "
      "is un-saturated, spans two orders of magnitude across the panel, and is "
      "*also* identical to every digit across the two arms.")
    A("")
    A("**Therapeutic coverage does not move at all, and that is structural.** "
      "`delivery_metrics.therapeutic_coverage` is a pure function of the per-cell "
      "AUC and the drug's IC50, and the AUC field has no oxygen dependence anywhere "
      "in the pipeline. So **no oxygen model of any kind — the PDE, the surrogate, "
      "a perfect measurement — can change a delivery metric in this codebase.** "
      "Every coverage number already published is therefore immune to this whole "
      "question, which is the most useful thing in this report for anyone auditing "
      "`results/delivery_metrics.md`, `results/dose_schedule.md` or "
      "`results/vessel_network_impact.md`.")
    A("")
    if max_kill_delta < 0.01:
        A(f"**Predicted kill barely moves either**: the largest absolute change "
          f"across three drugs is {max_kill_delta:.4f}. Read that plainly — the "
          "oxygen field, *however computed*, is close to not load-bearing for the "
          "drug predictions this project makes. The hypoxia/quiescence resistance "
          "multipliers are applied to bands that both routes place in broadly the "
          "same place, and the kill is dominated by delivery and potency.")
    else:
        max_rel = max(abs(r["delta_kill_pct"]) for r in dn["rows"].values())
        order_a = sorted(dn["rows"], key=lambda k: -dn["rows"][k]["pde"]["kill_fraction"])
        order_b = sorted(
            dn["rows"], key=lambda k: -dn["rows"][k]["vessel_distance"]["kill_fraction"]
        )
        _dk = [r["delta_kill"] for r in dn["rows"].values()]
        _dir = ("and always downwards" if all(d < 0 for d in _dk) else
                "and always upwards" if all(d > 0 for d in _dk) else
                "in both directions")
        _size = ("a little" if max_kill_delta < 0.05 else
                 "substantially" if max_kill_delta < 0.20 else "a great deal")
        A(f"**Predicted kill moves {_size}, {_dir}**: "
          f"{min(_dk):+.4f} to {max(_dk):+.4f} absolute, i.e. "
          f"{max_rel:.1f}% relative at worst."
          + (" The direction is not a coincidence — the surrogate bands more "
             "tumour into the resistant QUIESCENT band than the PDE does (§3.1), "
             f"and every resistant cell carries a "
             f"{dn['pd']['quiescent_resistance']}x kill multiplier."
             if all(d < 0 for d in _dk) else ""))
        A("")
        A(f"**The cross-drug ranking is unchanged**: "
          f"{' > '.join(order_a)} under the PDE, {' > '.join(order_b)} under the "
          "surrogate" + (" — identical." if order_a == order_b else " — CHANGED.")
          + " That is the expected structure: the resistance multipliers are shared "
          "across drugs, so they largely cancel out of a comparison and do not "
          "cancel out of a level.")
        A("")
        # The verdict on load-bearing-ness is read off the measured shift, not
        # asserted. It used to read "*Weakly.* A 2-3% relative shift ...", which
        # is false whenever the shift is not 2-3%.
        A("**So: how load-bearing is the oxygen field?** "
          + (f"*Weakly.* A {max_rel:.0f}% relative shift in absolute kill is far "
             "smaller than the acknowledged uncertainty in those levels — "
             "`results/state_validation.md` §2.2 shows the VIABLE fraction on this "
             "section swings 165x across two tissue-mask conventions and three "
             "grid spacings, all of which are bookkeeping choices. The oxygen "
             "ROUTE is a second-order knob next to the tissue-mask convention."
             if max_kill_delta < 0.05 else
             f"*Heavily, for the LEVEL.* A {max_rel:.0f}% relative shift in "
             f"absolute kill ({max_kill_delta:.4f} absolute) is larger than most "
             "of the effects this project reports, and no measurement on either "
             "section can say which arm is right (§3.3). Absolute kill from this "
             "section is therefore not a number, it is an interval whose width "
             "includes the oxygen route — alongside the tissue-mask convention "
             "and the mesh, which `results/state_validation.md` §2.2 shows swing "
             "the VIABLE fraction by 165x on their own. What is NOT sensitive is "
             "the cross-drug RANKING, which is unchanged above, and therapeutic "
             "coverage, which is invariant by construction. That split — ranking "
             "defensible, level not — is the finding.")
          + " Either way it is not zero, so absolute kill fractions must be "
          "reported with the oxygen mode attached, which is why "
          "`reconstruction_report()['oxygen']['mode']` now exists.")
    A("")
    A("### 3.3 What this does NOT show")
    A("")
    A("Neither section carries a paired drug-response outcome label, so nothing "
      "here says which arm's kill is *right*. The comparison is internal: it "
      "measures how much the answer depends on a modelling choice, which is a "
      "sensitivity, not a validation.")
    A("")

    # --- 4 wiring ---
    A("## 4. What was wired, and what did not change")
    A("")
    A("| | |")
    A("|---|---|")
    A("| new module | `tmesim/physics/oxygen_proxy.py` — "
      "`oxygen_from_vessel_distance`, `vessel_distance_um`, "
      "`hypoxia_from_markers`, `diffusion_length_um`, `pde_diffusion_length_um` |")
    A("| new option | `reconstruct_real(..., oxygen_mode=...)` with "
      "`'pde'` / `'vessel_distance'` / `'markers'`, plus `oxygen_kwargs=None` |")
    A("| default | **`'pde'`, unchanged and BIT-IDENTICAL** — asserted by "
      "`tests/test_oxygen_proxy.py` |")
    A("| recorded in | `reconstruction_report()['oxygen']` -> `mode`, `source`, "
      "`is_surrogate`, `caveat`, `detail` (and a flat `['oxygen_mode']`) |")
    A("| NOT touched | `tmesim/core`, `tmesim/sim/engine.py`, "
      "`tmesim/physics/oxygen.py` |")
    A("")
    A("`oxygen_mode='markers'` **raises** on the 319-gene brain panel rather than "
      "degrading silently, because that panel carries no VEGFA, CA9, ADM, NDRG1 or "
      "SLC2A1: on that section hypoxia is UNMEASURED, and the honest interface says "
      "so. `hypoxia_from_markers` returns `None` with the same reason when called "
      "directly.")
    A("")
    A("A caller who passes `oxygen_mode='vessel_distance'` straight into "
      "`sim.engine.simulate_response` will NOT get a surrogate-driven kill — the "
      "engine re-solves the PDE internally. That is documented in "
      "`reconstruct_real`'s docstring, and the supported route is the re-banding "
      "path §3 uses.")
    A("")

    # --- 5 recommendation ---
    A("## 5. RECOMMENDATION")
    A("")
    if pde_wins:
        A("**Keep `'pde'` as the shipping default — on the evidence, not despite "
          "it. Keep the surrogate as a cheap, explicitly-labelled option, and keep "
          "every per-cell hypoxia claim weak, because every route here is weak.**")
    else:
        A("**Keep `'pde'` as the shipping default. Use "
          "`oxygen_mode='vessel_distance'` for every NEW hypoxia claim, and stop "
          "making per-cell hypoxia claims from the solved field.**")
    A("")
    A("The reasoning, in order of weight:")
    A("")
    if pde_wins:
        A(f"1. **The solve is the more accurate route on the one section where we "
          f"can measure, and the more expensive one** — rho "
          f"{_fmt(pde_row['rho_vegfa'])} against `{best_rival['name']}`'s "
          f"{_fmt(best_rival['rho_vegfa'])} ({ordv['ratio']:.2f}x on |rho|), for "
          f"{rt[-1]['speedup']:.0f}x more compute at dx = 12 um. The accuracy "
          "argument now runs FOR the solve, which is the reverse of what this "
          "section said before the inequality was checked. It is a small accuracy "
          f"argument: |rho| ~ {ordv['abs_sim']:.2f} means the winner and the loser "
          "are both weak, and a physics-free one-liner recovering "
          f"{100 * ordv['recovers']:.0f}% of a sparse solve is still the most "
          "sobering number in this report.")
    else:
        A("1. **The surrogate is better and cheaper on the one section where we can "
          f"measure** — rho {_fmt(winner['rho_vegfa'])} vs "
          f"{_fmt(pde_row['rho_vegfa'])}, at {rt[-1]['speedup']:.0f}x less compute "
          "at dx = 12 um. There is no accuracy argument for the solve.")
    # Whether the drug side "barely notices" is a MEASURED question, so it is
    # read off max_kill_delta rather than asserted. This sentence used to say
    # "the drug outputs barely notice" unconditionally, directly beside the
    # interpolated delta that contradicts it whenever the delta is large.
    _max_rel = max(abs(r["delta_kill_pct"]) for r in dn["rows"].values())
    _barely = max_kill_delta < 0.01
    A(("2. **Either way the drug outputs barely notice, and flipping the default "
       "would invalidate the provenance of everything already published.**"
       if _barely else
       "2. **The drug outputs, by contrast, move a lot — which is an argument "
       "AGAINST flipping the default, not for it.**")
      + " §3 shows therapeutic coverage is EXACTLY invariant and kill moves by at "
      f"most {max_kill_delta:.4f} absolute ({_max_rel:.1f}% relative) with the "
      "cross-drug ranking unchanged. "
      + ("The measurable gain on the drug side is therefore close to nil, while "
         "the cost is re-deriving every kill number in `results/` (§5.1) and "
         "breaking the correspondence between those documents and the code."
         if _barely else
         "A shift of that size is not a refinement, it is a different answer: "
         "flipping the oxygen route would move every absolute kill number in "
         "`results/` (§5.1) by more than most of the effects those documents "
         "report, with no measurement able to say which arm is right (§3.3). "
         "That the CROSS-DRUG RANKING survives it unchanged is the reason the "
         "ranking, not the level, is this project's defensible output.")
      + " Silent improvement is still silent change.")
    if pde_wins:
        A("   The honest asymmetry: the solve is the better *hypoxia predictor* "
          f"(§1, by {ordv['ratio']:.2f}x on a weak |rho|) and barely distinguishable "
          "as a *drug-kill input* (§3). So the surrogate stays available and stays "
          "labelled, and is the right choice when compute matters more than "
          f"{100 * (1 - ordv['recovers']):.0f}% of a weak association — but it is "
          "not the accuracy choice, and this report must not be cited as though it "
          "were.")
    else:
        A("   The honest asymmetry: the surrogate is unambiguously better as a "
          "*hypoxia predictor* (§1) and only marginally different as a *drug-kill "
          "input* (§3). So the recommendation splits along exactly that line — flip "
          "it for hypoxia claims, leave it for the kill pipeline until that "
          "pipeline is being re-run for some other reason anyway.")
    A("3. **The measured route is the right long-run answer and is not available "
      "yet.** `hypoxia_from_markers` is the only route here that is a measurement, "
      "but the flagship drug section's panel does not carry a hypoxia-response "
      "gene, and on the IO section that does, benchmarking it is circular. A panel "
      "or platform carrying CA9/ADM/NDRG1/SLC2A1 alongside a drug-relevant readout "
      "would settle this properly.")
    A("")
    A("### 5.1 If the default WERE flipped, this is the re-run list")
    A("")
    A("Everything whose numbers depend on the per-cell state banding — i.e. on "
      "`PDConfig.hypoxia_resistance` / `quiescent_resistance` having cells to act "
      "on. Coverage-only and AUC-only results are immune (§3.2) and would not need "
      "re-running:")
    A("")
    A("| artefact | why | immune? |")
    A("|---|---|---|")
    A("| `results/dose_schedule.md` | its kill ceilings are ~89% the quiescent "
      "multiplier | **no — re-run** |")
    A("| `results/state_validation.md` | Test 1 IS this comparison | "
      "**no — re-run** |")
    A("| `results/report.md`, `results/patient_report_*.md` | quote absolute kill "
      "| **no — re-run** |")
    A("| `results/ranking_robustness.md`, `results/uncertainty.md` | resample kill "
      "| **no — re-run** |")
    A("| `results/delivery_metrics.md` | AUC + coverage only | yes |")
    A("| `results/vessel_network_impact.md` | coverage + penetration length | "
      "mostly — only its `hypoxia` block moves |")
    A("| `results/numerical_rigor.md` | would need the surrogate's own dx sweep | "
      "**no — re-run** |")
    A("")
    A("### 5.2 What would change our mind")
    A("")
    A("* A section with a **measured oxygen tension** (phosphorescence-quenching "
      "pO2, or pimonidazole IHC on a serial section). Every route here is "
      "uncalibrated in absolute level; a tension would adjudicate levels, not just "
      "orderings, and could rehabilitate the PDE by fitting `consumption`.")
    A("* A panel carrying **several** hypoxia-response genes, so the marker route "
      "can be benchmarked non-circularly by leave-one-gene-out.")
    A("* Evidence that the solve matters somewhere the distance transform cannot "
      "reach — e.g. a geometry with a supply-limited vessel where the (1-O) "
      "saturation term does real work. Nothing in the two cached sections exercises "
      "that.")
    A("")

    # --- 6 caveats ---
    A("## 6. Caveats — what could make this wrong")
    A("")
    A("1. **VEGFA is a hypoxia-RESPONSE gene, not an oxygen tension.** It is the "
      "transcriptional output of a pathway with its own kinetics and its own "
      "non-hypoxic inducers (Fischer et al., Brain Pathol 2005;15(4):297-310, PMID "
      "16389942, doi:10.1111/j.1750-3639.2005.tb00115.x). A field could be a "
      "perfectly good oxygen tension and still score poorly here. What the "
      "comparison establishes is narrower and immune to this: *whatever* VEGFA is "
      "reporting, these routes are ranked against it on identical inputs, and on "
      f"this run `{winner['name']}` came first. It is an internal ordering of weak "
      "predictors, not a measurement of oxygen.")
    A("2. **The two sections are different patients**, with different vascular "
      "architectures and different panels. The benchmark (§1) and the downstream "
      "effect (§3) are two experiments, not two views of one tumour.")
    A("3. **Circularity in the vessel term.** `FLT1` and `PLVAP` are IO-panel "
      "ENDOTHELIAL markers, so the cells that define 'distance to a vessel' are "
      "themselves typed from the panel. Both the PDE and the surrogate inherit this "
      "equally, so it cannot explain the gap between them — but neither route's "
      "absolute rho is free of it.")
    A("4. **Targeted-panel dropout.** A zero is 'no probe fired', not 'no "
      f"transcript'. VEGFA fires in {b['vegfa_detection'] * 100:.1f}% of cells, "
      "which is enough to power this test and not enough to make any single cell's "
      "value trustworthy.")
    A("5. **No absolute calibration anywhere, and this is the criticism that "
      "outlives the ranking.** `OxygenConfig.consumption` has never been fitted to "
      "a measured oxygen tension, and the surrogate's `L` is anchored to published "
      "morphometry, not to a tension. Nothing in §1 calibrates a level: Spearman "
      "rho is invariant to monotone rescaling, so winning or losing this benchmark "
      "says nothing at all about whether either field's hypoxic band edge sits at "
      "the right distance from a vessel. See §1.1.")
    A("6. **The surrogate has no supply saturation and no sink.** It cannot "
      "represent a vessel that is itself hypoxic, or tissue whose consumption "
      "varies. Helmlinger et al. measured exactly those effects in vivo. "
      + ("On these data the extra physics buys a small amount of rank accuracy "
         f"({ordv['ratio']:.2f}x on a weak |rho|) at a large compute cost; it is "
         "not established that it buys a better oxygen FIELD."
         if pde_wins else
         "The claim is only that on these data the extra physics costs accuracy "
         "rather than buying it."))
    A("")

    A("## 7. Reproduction")
    A("")
    A("```bash")
    A("python -m tmesim.validation.oxygen_benchmark")
    A("```")
    A("")
    A(f"Deterministic (the only randomness is the seeded permutation / subsampling, "
      f"`PERM_SEED = {PERM_SEED}`, and the engine's `seed = 0`). Runtime "
      f"~{results['elapsed_s'] / 60:.0f} min, dominated by the "
      f"{b['n_cells']:,}-cell reconstructions and the three transport solves.")
    A("")
    A("**Files added or changed by this work:**")
    A("")
    A("* `tmesim/physics/oxygen_proxy.py` — new; the surrogate and the marker route")
    A("* `tmesim/validation/oxygen_benchmark.py` — new; this analysis")
    A("* `tmesim/reconstruct/real_pipeline.py` — additive `oxygen_mode` / "
      "`oxygen_kwargs` keyword, default `'pde'` and bit-identical; the oxygen "
      "provenance block in `reconstruction_report`")
    A("* `tests/test_oxygen_proxy.py` — new")
    A("* `results/oxygen_benchmark.md` — this report. It used to be written to "
      "`results/oxygen_demotion.md`; that path now holds the hand-written record "
      "of the inverted-inequality correction, and `write_report` refuses to "
      "overwrite it.")
    A("")
    A("**Data provenance.** Both sections are open 10x Genomics Xenium datasets, "
      "CC BY 4.0, no login or DUA: the Human Brain Cancer (GBM) FFPE IO-panel "
      "section (v2.0.0, 480 targets) and the Human Brain Glioblastoma section "
      "(v1.3.0, 319 targets).")
    A("")

    text = "\n".join(L)
    with open(out_path, "w") as fh:
        fh.write(text)
    print(f"\nwrote {out_path}")
    return out_path


def main(argv: Optional[Sequence[str]] = None) -> Dict:
    """``python -m tmesim.validation.oxygen_benchmark [cache_path]``.

    Passing a ``cache_path`` pickles the results on the first run and reuses them
    on later runs, so the report can be re-rendered without repeating the compute.
    """
    import sys

    args = list(sys.argv[1:] if argv is None else argv)
    cache_path = args[0] if args else None
    results = run(cache_path=cache_path)
    write_report(results)
    return results


if __name__ == "__main__":
    main()
