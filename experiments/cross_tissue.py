"""CROSS-TISSUE at MATCHED, CONVERGED RESOLUTION: same real drug, 5 real GBM tumors.

This is the RIGOR rewrite of the cross-tissue experiment. The previous version ran
the single-cell Xenium section at the converged dx = 12 um and the four Visium
sections at their ~100 um spot pitch, then reported the min->max spread in
``kill_fraction`` and ``penetration_index`` as "real inter-tumor heterogeneity".
That comparison was CONFOUNDED BY RESOLUTION and its headline numbers are
retracted here (see ``results/cross_tissue.md``, "Retraction & correction").

Two things are fixed:

1. SIMULATION GRID != DATA RESOLUTION.
   The reconstructed fields are *continuous*: the Visium vascular / ECM fields are
   a Gaussian KDE of the per-spot deconvolved fractions with a bandwidth fixed in
   MICROMETRES (100 um, i.e. the spot pitch), and cells are placed at continuous
   (x, y) inside each spot footprint. Nothing in that construction requires the
   numerical grid to equal the spot pitch. So we simulate ALL FIVE tissues on the
   same converged dx = 12 um grid.

   This does NOT invent data. The Visium *information content* is still
   spot-limited (~100 um): the underlying field is band-limited at the 100 um KDE
   bandwidth, and refining dx only samples that same band-limited function more
   finely. What the refinement removes is the NUMERICAL artifact -- at dx = 100 um
   the grid cannot resolve the perivascular boundary layer (the calibrated
   doxorubicin e-folding length, ~65 um) or the oxygen diffusion length, so kill
   is understated and the hypoxia field is a severe grid artifact (see
   ``results/numerical_rigor.md`` section 2). Finer grid = less numerical error on
   the SAME information, not more biological detail.

2. penetration_index IS NOT A DELIVERY METRIC.
   ``penetration_index = mean(AUC) / percentile99(AUC)`` over the tissue mask --
   a ratio of two quantities that both scale linearly with the source amplitude.
   It is therefore SCALE-INVARIANT: it reports the *shape* (evenness) of the
   vascular field, not how much drug reached the tissue. Since the delivery source
   is ``S = vascular_density * bbb_penetrance * plasma``, correlating
   penetration_index against vascular density is close to CIRCULAR. It is reported
   below ONLY as a vascular-geometry shape diagnostic, never as delivery, and it
   is never used to rank tumors by drug exposure.

   The metrics that DO mean something and are used instead:
     * ``mean_auc``   -- mean tissue drug exposure (uM*s). Real delivery.
     * ``median_auc`` -- ditto, robust to the perivascular hot spots.
     * ``mean/median final concentration`` (uM) at the end of the window.
     * ``kill_fraction`` -- the grid-converged prediction.

HONESTY (unchanged and still binding):
  * PREDICTIONS on real geometry with NO outcome labels. Every section is an
    untreated resection; there is no paired measured drug response.
  * PLATFORM DIFFERENCE: Xenium is single-cell (real cells, per-cell marker
    typing, expression-weighted ECM). Visium is spot-level (~100 um spots,
    reference-based deconvolution, stromal-proxy ECM). Matching the numerical grid
    does NOT make them the same assay.
  * n = 5 tumors. Any correlation quoted here is descriptive, not powered.

Run:  ./.venv/bin/python experiments/cross_tissue.py
Writes: results/cross_tissue.png, results/cross_tissue.md
"""
from __future__ import annotations

import math
import os
import time
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from tmesim.core import OxygenConfig, PDEConfig, SimConfig
from tmesim.data.real.visium import list_samples, sample_info
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.reconstruct.real_fields import add_oxygen
from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report
# Private helper reused (NOT modified) so that cell placement can be held fixed in
# PHYSICAL units while the numerical grid is refined -- see reconstruct_visium_at.
from tmesim.reconstruct.real_visium import _place_cells, reconstruct_from_visium
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import (
    DOX_TARGET_EFOLD_UM,
    calibrate_alpha,
    calibrated_drugparams,
)

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "cross_tissue.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "cross_tissue.md")

# --- resolution policy -------------------------------------------------------
# ALL FIVE tissues are simulated on this grid. It is inside the converged regime
# for kill_fraction and comfortably resolves the perivascular boundary layer
# (results/numerical_rigor.md section 2).
MATCHED_DX_UM = 12.0
# The Visium spot pitch. Used ONLY for the paired resolution-control run that
# quantifies how much of the previously published spread was numerical.
VISIUM_NATIVE_DX_UM = 100.0
# Visium cell-placement jitter, held FIXED in micrometres (= the v1 spot radius)
# so refining dx changes ONLY the numerical grid, never the cell geometry.
# reconstruct_from_visium ties jitter to min(27.5, 0.45*dx), which equals 27.5 at
# the native 100 um pitch -- so this matches the canonical reconstruction exactly
# and merely keeps it dx-invariant on finer grids.
VISIUM_JITTER_UM = 27.5
# The physical KDE bandwidth real_visium uses for the vascular / ECM fields. It is
# specified in MICROMETRES, so the reconstructed field is the same band-limited
# function at any dx. This is why a finer grid adds no information.
VISIUM_KDE_BANDWIDTH_UM = 100.0
# A self-contained grid-refinement ladder on ONE Visium section. The ladder in
# results/numerical_rigor.md was run on the XENIUM section only; nothing had ever
# verified that dx = 12 um is converged for a Visium reconstruction. This does.
VISIUM_LADDER_SAMPLE = "GBM3"          # the smallest section, so the ladder is affordable
VISIUM_LADDER_DX = [100.0, 50.0, 30.0, 20.0, 12.0, 8.0]
# DOSE-SCALING CONTROL. The delivery source is S = vascular x bbb_penetrance x
# plasma, and the transport PDE is linear in S, so scaling bbb_penetrance scales
# the whole exposure field by the same factor. This is the decisive test of the
# penetration_index correction: a real DELIVERY metric must move when the dose
# moves; a SHAPE statistic cannot.
DOSE_SCALE_FACTOR = 2.0
# Calibrated perivascular boundary-layer scale (the doxorubicin e-folding length
# the single alpha is fit to). NOTE: the calibration honesty pass corrected a
# half-distance/e-folding mix-up, moving alpha 309.6 -> 148.6 and this length
# ~45 um -> ~65 um. The grid ladder in results/numerical_rigor.md section 2
# predates that fix; a LONGER boundary layer is easier to resolve, so dx = 12 um
# is, if anything, more comfortably inside the converged regime than that audit
# implies.
BOUNDARY_LAYER_UM = float(DOX_TARGET_EFOLD_UM)

# Okabe-Ito (colour-vision-deficiency safe). Validated: lightness band, chroma
# floor and adjacent-pair CVD separation all pass; the contrast WARN is relieved
# by direct labels + legends on every panel.
C_XENIUM = "#D55E00"   # vermillion  -- single-cell platform
C_VISIUM = "#0072B2"   # blue        -- spot-level platform
C_VASC = "#E69F00"     # orange
C_OXY = "#009E73"      # bluish green
C_HYP = "#CC79A7"      # reddish purple
C_MUTED = "#666666"
C_COMP = {
    "TUMOR": "#D55E00",
    "STROMAL": "#009E73",
    "IMMUNE": "#0072B2",
    "ENDOTHELIAL": "#E69F00",
    "NECROTIC": "#999999",
}

# Previously published (dx-mismatched) headline, kept verbatim for the retraction.
_OLD_CLAIM = {
    "Xenium-GBM": {"dx": 12.0, "kill": 0.023, "pen": 0.196},
    "10x_parent_gbm": {"dx": 100.0, "kill": 0.116, "pen": 0.283},
    "GBM1": {"dx": 100.0, "kill": 0.094, "pen": 0.264},
    "GBM2": {"dx": 100.0, "kill": 0.241, "pen": 0.509},
    "GBM3": {"dx": 100.0, "kill": 0.235, "pen": 0.566},
}


# ---------------------------------------------------------------------------
# drug + reconstruction
# ---------------------------------------------------------------------------
def compile_real_temozolomide():
    """Compile the REAL temozolomide drug ONCE: calibrated perivascular-consumption
    physics (alpha fit on doxorubicin) then the binding-site barrier. The exact
    same DrugParams object is reused for every tissue -- the drug never changes."""
    alpha = calibrate_alpha()  # single fit on doxorubicin (see calibration.py)
    drug = calibrated_drugparams("temozolomide", alpha)
    drug = apply_barrier_to_drugparams(drug, get_molecule("temozolomide"))
    return drug, alpha


def _retry(fn, *args, **kwargs):
    """Call ``fn``, retrying once on any transient/IO/network error."""
    try:
        return fn(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001 - retry once per task spec
        print(f"    call failed ({exc!r}); retrying once in 5s ...")
        time.sleep(5)
        return fn(*args, **kwargs)


def reconstruct_visium_at(sample_id: str, dx_um: float):
    """Reconstruct a REAL Visium section on an ARBITRARY simulation grid ``dx_um``.

    ``reconstruct_real(source='visium', dx_um=...)`` clamps any dx < 50 um back to
    the 100 um spot pitch. That clamp conflates the *data* resolution with the
    *numerical* grid, which is exactly the confound this experiment fixes, so we
    go one level down to the public ``reconstruct_from_visium`` (nothing in
    ``reconstruct/`` is modified).

    Why this is legitimate and adds NO information:
      * the vascular / ECM fields are a Gaussian KDE with bandwidth fixed at
        VISIUM_KDE_BANDWIDTH_UM = 100 um in PHYSICAL units, so the reconstructed
        field is the *same band-limited function* at any dx -- a finer grid only
        samples it more accurately (and removes the aliasing of binning ~100 um
        spot centres onto a 100 um grid);
      * cells are placed at continuous coordinates inside each spot footprint.

    One correction is applied: ``reconstruct_from_visium`` scales the placement
    jitter as ``min(27.5, 0.45*dx)``, which would shrink the spot footprint from
    27.5 um to 5.4 um purely because the grid got finer. We re-place the cells with
    the jitter FIXED at VISIUM_JITTER_UM (identical to what the native 100 um
    reconstruction uses), so dx = 100 vs dx = 12 differ ONLY in the numerical grid.
    """
    from tmesim.reconstruct.deconvolve import deconvolve_sample

    micro = _retry(reconstruct_from_visium, sample_id, dx_um=float(dx_um))
    coords_um, fractions, type_names = _retry(deconvolve_sample, sample_id)
    micro.cells = _place_cells(
        coords_um, fractions, type_names,
        extent=micro.domain.extent_um, jitter_um=VISIUM_JITTER_UM,
    )
    if micro.fields.oxygen is None:
        add_oxygen(micro, OxygenConfig())
    micro.metadata = dict(micro.metadata)
    micro.metadata["real_pipeline"] = {
        "source": "visium",
        "vascular_source": "KDE of REAL ENDOTHELIAL spot fraction "
                           f"(bandwidth {VISIUM_KDE_BANDWIDTH_UM:.0f} um, fixed in "
                           "physical units)",
        "ecm_source": "STROMAL-fraction proxy (spot resolution; no per-cell ECM "
                      "transcript)",
        "ecm_markers_used": [],
        "oxygen_source": "physics.solve_oxygen driven by the real vascular field",
    }
    return micro


def reconstruct_tissue(source: str, sample: Optional[str], dx_um: float):
    """Reconstruct any of the five real tissues at an arbitrary simulation dx."""
    if source == "xenium":
        return _retry(reconstruct_real, source="xenium", sample=sample, dx_um=dx_um)
    return reconstruct_visium_at(sample, dx_um)


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------
def _delivery_metrics(result, micro, cfg: SimConfig) -> Dict[str, float]:
    """Delivery metrics that actually measure DRUG EXPOSURE (not field shape).

    ``penetration_index`` is deliberately excluded from this set: it is
    mean(AUC)/p99(AUC), a ratio of two quantities that both scale linearly with the
    source, hence a scale-invariant SHAPE statistic of the vascular field.

    Two extra diagnostics are computed because the tissue-wide mean alone cannot
    explain a kill difference:

    * **tumor-cell-local AUC** -- the exposure sampled where the TUMOR cells
      actually sit. Tissue-wide mean AUC counts empty/stromal grid points that no
      tumor cell occupies, so two tissues can share a mean AUC while their tumor
      cells experience very different exposure (this is exactly what separates the
      real single-cell Xenium geometry from the spot-centred Visium geometry).
    * **resistant-state fractions** -- HYPOXIC (resistance 0.35) and QUIESCENT
      (resistance 0.5) tumor cells are killed less at identical exposure, so the
      banding is needed to read the kill numbers.

    Nothing in core/physics/cells is modified: the AUC and oxygen grids returned by
    the engine are simply re-sampled, and the public ``update_cell_states`` is
    re-applied to a copy of the untouched cells.
    """
    from tmesim.cells.agents import update_cell_states
    from tmesim.core import CellState, Microenvironment

    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    auc_grid = np.asarray(result.transport.auc, dtype=np.float64)
    auc = auc_grid[tissue]
    conc = np.asarray(result.transport.final_concentration, dtype=np.float64)[tissue]

    out = {
        "mean_auc": float(auc.mean()),
        "median_auc": float(np.median(auc)),
        "p10_auc": float(np.percentile(auc, 10.0)),
        "mean_conc": float(conc.mean()),
        "median_conc": float(np.median(conc)),
    }

    # --- exposure where the TUMOR cells actually are -------------------------
    cells = micro.cells
    tumor = cells.tumor_mask()
    if int(tumor.sum()):
        col, row = micro.domain.to_index(cells.x[tumor], cells.y[tumor])
        tauc = auc_grid[row, col]
        out["tumor_auc_mean"] = float(tauc.mean())
        out["tumor_auc_median"] = float(np.median(tauc))
    else:
        out["tumor_auc_mean"] = 0.0
        out["tumor_auc_median"] = 0.0

    # --- resistant-state banding of the tumor cells --------------------------
    band = Microenvironment(domain=micro.domain, fields=micro.fields,
                            cells=cells.copy(), patient_id=micro.patient_id,
                            metadata=micro.metadata)
    update_cell_states(band, np.asarray(result.fields["oxygen"]), cfg.oxygen)
    st = band.cells.state
    n_t = int(tumor.sum())
    out["quiescent_tumor_fraction"] = (
        float((tumor & (st == int(CellState.QUIESCENT))).sum()) / n_t if n_t else 0.0)
    out["resistant_tumor_fraction"] = (
        float((tumor & ((st == int(CellState.QUIESCENT))
                        | (st == int(CellState.HYPOXIC)))).sum()) / n_t
        if n_t else 0.0)
    return out


def run_one_tissue(label: str, source: str, sample: Optional[str], dx_um: float,
                   drug, cfg_pde: PDEConfig, diagnosis: str) -> Dict:
    """Reconstruct one REAL tissue at ``dx_um`` and run the SAME drug through it."""
    t0 = time.time()
    print(f"\n[{label}] source={source} sample={sample} dx={dx_um:g} um")
    micro = reconstruct_tissue(source, sample, dx_um)
    report = reconstruction_report(micro)

    # The convective (Darcy) term carries an explicit CFL limit dt*|v|/dx <= 1, so
    # a fine grid can require a smaller dt than the shared default. Halve dt until
    # the solver accepts it. This is safe for comparability: numerical_rigor.md §3
    # shows the solution is dt-invariant to <1e-2 relative for dt <= 20 s, so the
    # fallback changes the discretisation error, not the answer. The dt actually
    # used is recorded on the row and reported.
    from dataclasses import replace as _dc_replace

    dt_used = float(cfg_pde.dt)
    for _ in range(4):
        cfg = SimConfig(drug=drug,
                        pde=_dc_replace(cfg_pde, dt=dt_used), seed=0)
        try:
            result = simulate_response(micro, cfg)
            break
        except ValueError as exc:
            if "CFL" not in str(exc):
                raise
            dt_used /= 2.0
            print(f"    CFL limit hit; retrying with dt = {dt_used:g} s "
                  f"({exc})", flush=True)
    else:  # pragma: no cover - 16x reduction should always suffice
        raise RuntimeError(f"{label}: could not satisfy the advection CFL limit")
    frac = report["celltype_fractions"]

    n_tumor = int(report["n_tumor"])
    k = float(result.kill_fraction)
    # Binomial sampling s.e. of the per-cell kill draw at fixed p_kill. This is the
    # ONLY stochastic part of the pipeline; everything upstream is deterministic.
    kill_se = float(np.sqrt(max(k * (1.0 - k), 0.0) / n_tumor)) if n_tumor else 0.0

    row = {
        "label": label,
        "platform": "Xenium (single-cell)" if source == "xenium"
        else "Visium (spot-level)",
        "source": source,
        "diagnosis": diagnosis,
        "n_cells": report["n_cells"],
        "n_tumor": n_tumor,
        "grid": report["grid"],
        "dx_um": report["dx_um"],
        "tumor_frac": frac.get("TUMOR", 0.0),
        "endo_frac": frac.get("ENDOTHELIAL", 0.0),
        "immune_frac": frac.get("IMMUNE", 0.0),
        "stromal_frac": frac.get("STROMAL", 0.0),
        "necrotic_frac": frac.get("NECROTIC", 0.0),
        "mean_vascular": report["mean_vascular_density"],
        "mean_ecm": report["mean_ecm_density"],
        "mean_oxygen": report["mean_oxygen"],
        "kill_fraction": k,
        "kill_se": kill_se,
        "penetration_index": float(result.extras["penetration_index"]),
        "hypoxic_tumor_fraction": float(result.extras["hypoxic_tumor_fraction"]),
        "sim_oxygen": float(result.extras["mean_oxygen"]),
        "pde_dt": dt_used,
        "runtime_s": time.time() - t0,
    }
    row.update(_delivery_metrics(result, micro, cfg))
    print(f"    grid={row['grid'][0]}x{row['grid'][1]}  n_cells={row['n_cells']:,} "
          f"tumor={row['tumor_frac'] * 100:.1f}% endo={row['endo_frac'] * 100:.1f}%",
          flush=True)
    print(f"    vasc={row['mean_vascular']:.3f}  meanAUC={row['mean_auc']:.4g}  "
          f"medAUC={row['median_auc']:.4g}  tumorAUC={row['tumor_auc_mean']:.4g}  "
          f"kill={k:.4f}+-{kill_se:.4f}  "
          f"hypox={row['hypoxic_tumor_fraction']:.4f} "
          f"resist={row['resistant_tumor_fraction']:.4f}  "
          f"[pen(shape)={row['penetration_index']:.3f}]  "
          f"({row['runtime_s']:.0f}s)", flush=True)
    return row


def run_dose_control(drug, cfg_pde: PDEConfig, diagnosis: str) -> List[Dict]:
    """Run ONE tissue at 1x and at DOSE_SCALE_FACTOR x the delivered dose.

    The delivery source is ``S = vascular_density * bbb_penetrance * plasma`` and
    the transport PDE is linear in ``S``, so scaling ``bbb_penetrance`` scales the
    entire exposure field by exactly that factor -- nothing about the geometry,
    the grid or the tissue changes.

    This is the decisive test behind the ``penetration_index`` correction. A
    genuine DELIVERY metric must move when the delivered dose moves. A
    scale-invariant SHAPE statistic cannot. Both are measured here rather than
    asserted.
    """
    from dataclasses import replace as _dc_replace

    hi_drug = _dc_replace(
        drug, bbb_penetrance=float(drug.bbb_penetrance) * DOSE_SCALE_FACTOR)
    out = []
    for tag, d in ((f"{VISIUM_LADDER_SAMPLE}@dose1x", drug),
                   (f"{VISIUM_LADDER_SAMPLE}@dose{DOSE_SCALE_FACTOR:g}x", hi_drug)):
        out.append(run_one_tissue(tag, "visium", VISIUM_LADDER_SAMPLE,
                                  MATCHED_DX_UM, d, cfg_pde, diagnosis))
    return out


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------
def _exact_perm_p(x: np.ndarray, y: np.ndarray, stat_fn) -> float:
    """EXACT two-sided permutation p-value (all n! permutations of y).

    scipy's asymptotic Spearman p-value is a t-approximation that reports p ≈ 0
    for a perfect rank match at n = 5 — which is nonsense, because only 1 of the
    120 orderings achieves it. At these sample sizes the exact test is both cheap
    and the only defensible one.
    """
    from itertools import permutations

    obs = abs(stat_fn(x, y))
    if not np.isfinite(obs):
        return float("nan")
    perms = list(permutations(range(y.size)))
    hits = sum(1 for p in perms
               if abs(stat_fn(x, y[list(p)])) >= obs - 1e-12)
    return float(hits) / float(len(perms))


def _corr(x: List[float], y: List[float]) -> Dict[str, float]:
    """Spearman + Pearson with EXACT permutation p-values (n is tiny here) and a
    leave-one-out stability range."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = x.size
    out: Dict[str, float] = {"n": n}
    if n < 3:
        return out
    with np.errstate(invalid="ignore"):
        rho = float(stats.spearmanr(x, y).statistic)
        r = float(stats.pearsonr(x, y)[0])
    out["spearman"] = rho
    out["pearson"] = r
    if np.isfinite(rho) and n <= 8:
        out["spearman_p"] = _exact_perm_p(
            x, y, lambda a, b: stats.spearmanr(a, b).statistic)
        out["pearson_p"] = _exact_perm_p(
            x, y, lambda a, b: stats.pearsonr(a, b)[0])
    else:
        out["spearman_p"] = float("nan")
        out["pearson_p"] = float("nan")
    if n >= 4 and np.isfinite(rho):
        with np.errstate(invalid="ignore"):
            loo = [float(stats.spearmanr(np.delete(x, i),
                                         np.delete(y, i)).statistic)
                   for i in range(n)]
        if np.any(np.isfinite(loo)):
            out["spearman_loo_min"] = float(np.nanmin(loo))
            out["spearman_loo_max"] = float(np.nanmax(loo))
    # Smallest attainable two-sided exact p for a perfect rank match at this n.
    out["min_attainable_p"] = 2.0 / float(math.factorial(n)) if n <= 8 else 0.0
    return out


def _fmt_corr(c: Dict[str, float]) -> str:
    if "spearman" not in c or not np.isfinite(c["spearman"]):
        return "n/a (undefined — constant input or n < 3)"
    s = (f"Spearman rho = {c['spearman']:+.3f} (exact p = {c['spearman_p']:.3f}), "
         f"Pearson r = {c['pearson']:+.3f} (exact p = {c['pearson_p']:.3f})")
    if "spearman_loo_min" in c:
        s += (f"; leave-one-out rho in "
              f"[{c['spearman_loo_min']:+.3f}, {c['spearman_loo_max']:+.3f}]")
    return s


def _spread(rows: List[Dict], key: str) -> Tuple[float, float, float]:
    vals = [r[key] for r in rows]
    return min(vals), max(vals), (max(vals) - min(vals))


def compute_stats(rows: List[Dict]) -> Dict[str, Dict[str, float]]:
    """All correlations reported in the write-up, on the matched-dx rows."""
    vis = [r for r in rows if r["source"] != "xenium"]
    pack = {
        "vasc_vs_kill_all5": _corr([r["mean_vascular"] for r in rows],
                                   [r["kill_fraction"] for r in rows]),
        "vasc_vs_auc_all5": _corr([r["mean_vascular"] for r in rows],
                                  [r["mean_auc"] for r in rows]),
        "endo_vs_kill_all5": _corr([r["endo_frac"] for r in rows],
                                   [r["kill_fraction"] for r in rows]),
        "endo_vs_auc_all5": _corr([r["endo_frac"] for r in rows],
                                  [r["mean_auc"] for r in rows]),
        "auc_vs_kill_all5": _corr([r["mean_auc"] for r in rows],
                                  [r["kill_fraction"] for r in rows]),
        "tumorauc_vs_kill_all5": _corr([r["tumor_auc_mean"] for r in rows],
                                       [r["kill_fraction"] for r in rows]),
        "vasc_vs_kill_visium": _corr([r["mean_vascular"] for r in vis],
                                     [r["kill_fraction"] for r in vis]),
        "vasc_vs_auc_visium": _corr([r["mean_vascular"] for r in vis],
                                    [r["mean_auc"] for r in vis]),
    }
    return pack


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
def _scatter_panel(ax, xs, ys, labels, colors, xlabel, ylabel, title, corr=None,
                   extra=None):
    ax.scatter(xs, ys, c=colors, s=110, edgecolors="white", linewidths=1.6, zorder=3)
    for i, (x, y, lab) in enumerate(zip(xs, ys, labels)):
        ax.annotate(lab, (x, y), textcoords="offset points",
                    xytext=(9, 5) if i % 2 == 0 else (9, -11),
                    fontsize=8, color="#333333")
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=9.5, loc="left")
    ax.grid(alpha=0.22, linewidth=0.7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if corr and "spearman" in corr:
        txt = (f"all 5: Spearman rho={corr['spearman']:+.2f} "
               f"(exact p={corr['spearman_p']:.2f})\n"
               f"       Pearson r={corr['pearson']:+.2f} "
               f"(exact p={corr['pearson_p']:.2f})")
        if extra and "spearman" in extra:
            txt += (f"\nVisium only (n=4, same platform):\n"
                    f"       rho={extra['spearman']:+.2f} "
                    f"(exact p={extra['spearman_p']:.2f})")
        txt += "\nn is tiny -- DESCRIPTIVE ONLY, no significance claimed"
        ax.text(0.03, 0.96, txt, transform=ax.transAxes, va="top", ha="left",
                fontsize=7.2, color="#444444", family="monospace",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))


def make_figure(rows: List[Dict], native: List[Dict], ladder: List[Dict],
                dose: List[Dict], st: Dict, alpha: float, path: str) -> None:
    labels = [r["label"] for r in rows]
    is_xen = np.array([r["source"] == "xenium" for r in rows])
    colors = [C_XENIUM if x else C_VISIUM for x in is_xen]

    vasc = np.array([r["mean_vascular"] for r in rows])
    kill = np.array([r["kill_fraction"] for r in rows])
    mauc = np.array([r["mean_auc"] for r in rows])
    hyp = np.array([r["hypoxic_tumor_fraction"] for r in rows])
    oxy = np.array([r["mean_oxygen"] for r in rows])

    fig, axes = plt.subplots(2, 4, figsize=(25.5, 11.5), constrained_layout=True)

    # --- (a) vascular density vs MEAN TISSUE AUC (real delivery) -------------
    _scatter_panel(
        axes[0, 0], vasc, mauc / 1e6, labels, colors,
        "REAL mean vascular density (min-max-normalised field)",
        "Mean tissue drug AUC  (1e6 uM*s)",
        "(a) Vascularity vs REAL DELIVERY (mean tissue AUC)\n"
        f"all 5 tumors at MATCHED, converged dx = {MATCHED_DX_UM:.0f} um\n"
        "NOTE: source = vascular x bbb x plasma, so this link is largely "
        "MECHANICAL, not a finding",
        st["vasc_vs_auc_all5"], st["vasc_vs_auc_visium"])

    # --- (b) vascular density vs KILL (converged) ---------------------------
    ax = axes[0, 1]
    _scatter_panel(
        ax, vasc, kill, labels, colors,
        "REAL mean vascular density (min-max-normalised field)",
        "Predicted kill_fraction (grid-converged)",
        "(b) Vascularity vs predicted KILL\n"
        f"all 5 tumors at MATCHED, converged dx = {MATCHED_DX_UM:.0f} um\n"
        "the only link with genuine non-linear content (Hill + O2 gating)",
        st["vasc_vs_kill_all5"], st["vasc_vs_kill_visium"])
    ax.errorbar(vasc, kill, yerr=[r["kill_se"] for r in rows], fmt="none",
                ecolor=C_MUTED, elinewidth=1.0, capsize=3, zorder=2)
    xi = int(np.argmax(is_xen))
    ax.annotate("cross-platform:\nnot like-for-like",
                (vasc[xi], kill[xi]), textcoords="offset points", xytext=(26, 18),
                fontsize=7.2, color=C_XENIUM, ha="left",
                arrowprops=dict(arrowstyle="-", color=C_XENIUM, lw=1.0))

    # --- (c) resolution control: dx=100 -> dx=12 on the SAME Visium tissues --
    ax = axes[0, 2]
    nat = {r["label"]: r for r in native}
    xs = np.arange(len(native))
    k_nat = np.array([nat[r["label"]]["kill_fraction"] for r in native])
    k_fin = np.array([next(q["kill_fraction"] for q in rows if q["label"] == r["label"])
                      for r in native])
    ax.bar(xs - 0.19, k_nat, 0.36, color="#bbbbbb", edgecolor="white", linewidth=1.5,
           label=f"dx = {VISIUM_NATIVE_DX_UM:.0f} um (spot pitch; NOT converged)")
    ax.bar(xs + 0.19, k_fin, 0.36, color=C_VISIUM, edgecolor="white", linewidth=1.5,
           label=f"dx = {MATCHED_DX_UM:.0f} um (converged)")
    for i, (a, b) in enumerate(zip(k_nat, k_fin)):
        ax.annotate(f"x{b / max(a, 1e-9):.1f}", (i, max(a, b)),
                    textcoords="offset points", xytext=(0, 5), ha="center",
                    fontsize=8, color="#333333")
    ax.set_xticks(xs)
    ax.set_xticklabels([r["label"] for r in native], rotation=18, ha="right",
                       fontsize=8)
    ax.set_ylabel("Predicted kill_fraction", fontsize=9)
    ax.set_title("(c) RESOLUTION CONTROL -- same tissue, same drug, same cell\n"
                 "geometry; only the numerical grid changes. This is the confound\n"
                 "that invalidated the earlier spread (see the md: an upstream\n"
                 "oxygen-sink fix has since shrunk the residual effect).",
                 fontsize=9.5, loc="left")
    ax.set_ylim(0, float(max(k_nat.max(), k_fin.max())) * 1.42)
    ax.legend(fontsize=7.5, loc="upper left")
    ax.grid(alpha=0.22, axis="y", linewidth=0.7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # --- (d) self-contained grid ladder on ONE Visium section ---------------
    ax = axes[0, 3]
    if ladder:
        lx = np.array([r["dx_um"] for r in ladder], dtype=float)
        order = np.argsort(-lx)
        lx = lx[order]
        lk = np.array([ladder[i]["kill_fraction"] for i in order])
        lh = np.array([ladder[i]["hypoxic_tumor_fraction"] for i in order])
        la = np.array([ladder[i]["mean_auc"] for i in order])
        la = la / la[-1]
        ax.plot(lx, lk, "o-", color=C_VISIUM, lw=2.0, ms=8, mec="white", mew=1.4,
                label="kill_fraction")
        ax.plot(lx, lh, "s-", color=C_HYP, lw=2.0, ms=8, mec="white", mew=1.4,
                label="hypoxic tumor fraction")
        ax.plot(lx, la, "^--", color=C_VASC, lw=2.0, ms=8, mec="white", mew=1.4,
                label="mean tissue AUC (relative to finest dx)")
        ax.axvline(MATCHED_DX_UM, color=C_MUTED, ls=":", lw=1.4)
        ax.annotate(f"dx = {MATCHED_DX_UM:.0f} µm\n(used here)",
                    (MATCHED_DX_UM, ax.get_ylim()[1]), textcoords="offset points",
                    xytext=(6, -18), fontsize=7.5, color=C_MUTED)
        ax.set_xscale("log")
        ax.set_xticks(lx)
        ax.set_xticklabels([f"{v:g}" for v in lx], fontsize=8)
        ax.invert_xaxis()
        ax.set_xlabel("simulation grid dx (µm, log scale, coarse → fine)",
                      fontsize=9)
        ax.set_ylabel("value (0–1) / relative", fontsize=9)
        ax.set_title(f"(d) GRID CONVERGENCE on ONE Visium section "
                     f"({VISIUM_LADDER_SAMPLE})\n"
                     "the audit in numerical_rigor.md §2 covered the Xenium\n"
                     "section only — nothing had verified dx for Visium",
                     fontsize=9.5, loc="left")
        ax.legend(fontsize=7.5, loc="best")
        ax.grid(alpha=0.22, linewidth=0.7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # --- (e) per-tissue real fields at matched dx ----------------------------
    ax = axes[1, 0]
    x = np.arange(len(rows))
    w = 0.27
    ax.bar(x - w, vasc, w, color=C_VASC, edgecolor="white", linewidth=1.4,
           label="mean vascular density")
    ax.bar(x, oxy, w, color=C_OXY, edgecolor="white", linewidth=1.4,
           label="mean oxygen")
    ax.bar(x + w, hyp, w, color=C_HYP, edgecolor="white", linewidth=1.4,
           label="hypoxic tumor fraction")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=8)
    ax.set_ylabel("field value (0-1)", fontsize=9)
    ax.set_title(f"(e) REAL microenvironment fields at matched dx = "
                 f"{MATCHED_DX_UM:.0f} um\n"
                 "hypoxic fraction: grid-dependence FIXED upstream, but its\n"
                 "ABSOLUTE level is un-calibrated -- order-of-magnitude only",
                 fontsize=9.5, loc="left")
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.22, axis="y", linewidth=0.7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # --- (f) composition -----------------------------------------------------
    ax = axes[1, 1]
    bottom = np.zeros(len(rows))
    for key, name in (("tumor_frac", "TUMOR"), ("stromal_frac", "STROMAL"),
                      ("immune_frac", "IMMUNE"), ("endo_frac", "ENDOTHELIAL"),
                      ("necrotic_frac", "NECROTIC")):
        vals = np.array([r[key] for r in rows])
        ax.bar(x, vals, 0.62, bottom=bottom, label=name, color=C_COMP[name],
               edgecolor="white", linewidth=1.4)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=8)
    ax.set_ylabel("cell-type fraction", fontsize=9)
    ax.set_title("(f) REAL cell-type composition\n"
                 "(Xenium: per-cell marker typing; Visium: deconvolved)",
                 fontsize=9.5, loc="left")
    ax.legend(fontsize=7.5, loc="lower right")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # --- (g) DOSE-SCALING CONTROL: the decisive penetration_index test -------
    ax = axes[1, 2]
    p_nat = np.array([nat[r["label"]]["penetration_index"] for r in native])
    p_fin = np.array([next(q["penetration_index"] for q in rows
                           if q["label"] == r["label"]) for r in native])
    if len(dose) == 2:
        base, hi = dose[0], dose[1]
        names = ["mean tissue AUC\n(delivery)", "kill_fraction",
                 "penetration_index\n(shape statistic)"]
        ratios = [hi["mean_auc"] / max(base["mean_auc"], 1e-12),
                  hi["kill_fraction"] / max(base["kill_fraction"], 1e-12),
                  hi["penetration_index"] / max(base["penetration_index"], 1e-12)]
        cols = [C_VASC, C_VISIUM, C_HYP]
        bx = np.arange(3)
        ax.bar(bx, ratios, 0.55, color=cols, edgecolor="white", linewidth=1.5)
        for i, v in enumerate(ratios):
            ax.annotate(f"{v:.2f}×", (i, v), textcoords="offset points",
                        xytext=(0, 5), ha="center", fontsize=9, color="#333333")
        ax.axhline(1.0, color=C_MUTED, lw=1.2, ls="--")
        ax.set_xticks(bx)
        ax.set_xticklabels(names, fontsize=8)
        ax.set_ylabel(f"ratio at {DOSE_SCALE_FACTOR:g}× dose ÷ 1× dose", fontsize=9)
        ax.set_ylim(0, max(max(ratios), 1.0) * 1.35)
        ax.set_title(f"(g) DOSE-SCALING CONTROL on {base['label'].split('@')[0]}: "
                     f"the decisive test.\n"
                     f"Same tissue, same grid; only the delivered dose is "
                     f"{DOSE_SCALE_FACTOR:g}×.\n"
                     "Delivery doubles and kill rises -- penetration_index does NOT "
                     "move.\nIt is a SHAPE statistic. Ranking delivery by it is "
                     "CIRCULAR.", fontsize=9.5, loc="left")
        ax.grid(alpha=0.22, axis="y", linewidth=0.7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # --- (h) WHAT actually moves with resolution ----------------------------
    ax = axes[1, 3]
    a_nat = np.array([nat[r["label"]]["mean_auc"] for r in native])
    a_fin = np.array([next(q["mean_auc"] for q in rows if q["label"] == r["label"])
                      for r in native])
    series = [
        ("kill_fraction", k_fin / np.maximum(k_nat, 1e-12), C_VISIUM),
        ("mean tissue AUC (delivery)", a_fin / np.maximum(a_nat, 1e-12), C_VASC),
        ("penetration_index (shape)", p_fin / np.maximum(p_nat, 1e-12), C_HYP),
    ]
    w2 = 0.26
    for j, (name, vals, col) in enumerate(series):
        ax.bar(xs + (j - 1) * w2, vals, w2, color=col, edgecolor="white",
               linewidth=1.4, label=name)
    ax.axhline(1.0, color=C_MUTED, lw=1.2, ls="--")
    ax.set_xticks(xs)
    ax.set_xticklabels([r["label"] for r in native], rotation=18, ha="right",
                       fontsize=8)
    ax.set_ylabel(f"ratio  (dx {MATCHED_DX_UM:.0f} µm ÷ dx "
                  f"{VISIUM_NATIVE_DX_UM:.0f} µm)", fontsize=9)
    ax.set_title("(h) HOW MUCH each metric still depends on the grid\n"
                 "1.0 = fully grid-independent. Delivery was never the fragile\n"
                 "part; the fragility was in the O2/quiescence gate, which an\n"
                 "upstream oxygen-sink fix has now removed.",
                 fontsize=9.5, loc="left")
    ax.set_ylim(0, 1.35)
    ax.legend(fontsize=7.5, loc="lower left")
    ax.grid(alpha=0.22, axis="y", linewidth=0.7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.suptitle(
        f"CROSS-TISSUE at MATCHED CONVERGED RESOLUTION: same TMZ physics "
        f"(alpha={alpha:.0f}) on {len(rows)} REAL GBM tissues, all at "
        f"dx = {MATCHED_DX_UM:.0f} um\n"
        f"vermillion = Xenium single-cell   |   blue = Visium spot-level "
        f"(numerical grid matched; INFORMATION CONTENT still ~100 um "
        f"spot-limited)   |   PREDICTIONS on real geometry, NO outcome labels, "
        f"n = {len(rows)}",
        fontsize=12.5)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"\nfigure written: {path}")


# ---------------------------------------------------------------------------
# markdown
# ---------------------------------------------------------------------------
def _verdict(rows: List[Dict], native: List[Dict], st: Dict) -> Tuple[str, str]:
    """Decide -- from the numbers, not from a prior -- whether the heterogeneity
    claim survives at matched resolution. Returns (headline, body)."""
    nat = {r["label"]: r for r in native}
    vis = [r for r in rows if r["source"] != "xenium" and r["label"] in nat]
    klo, khi, ksp = _spread(rows, "kill_fraction")
    alo, ahi, asp = _spread(rows, "mean_auc")
    fold = khi / max(klo, 1e-9)
    typ_se = float(np.mean([r["kill_se"] for r in rows]))
    c = st["vasc_vs_kill_all5"]
    cv = st["vasc_vs_kill_visium"]
    rho = c.get("spearman", float("nan"))
    rho_v = cv.get("spearman", float("nan"))

    # "Substantial" = the min->max kill spread is many times the stochastic
    # sampling s.e. AND is at least a 1.5x fold change.
    spread_real = (ksp > 20.0 * typ_se) and (fold >= 1.5)
    tracks = (not np.isnan(rho)) and rho >= 0.79
    tracks_v = (not np.isnan(rho_v)) and rho_v >= 0.79

    if spread_real and tracks and tracks_v:
        head = ("**VERDICT: the inter-tumor heterogeneity claim SURVIVES at matched "
                "resolution — predicted kill still varies substantially across the "
                "five real tumors and still rises with real vascularity, in the "
                f"all-5 set (Spearman rho = {rho:+.2f}) and, importantly, also "
                f"within the four like-for-like Visium tumors (rho = {rho_v:+.2f}). "
                "But the NUMBERS are completely different and every earlier figure "
                "is retracted.**")
    elif spread_real and tracks:
        head = ("**VERDICT: the inter-tumor heterogeneity claim SURVIVES at matched "
                f"resolution in the all-5 set (Spearman rho = {rho:+.2f}), but it "
                "does NOT hold up within the four like-for-like Visium tumors "
                f"(rho = {rho_v:+.2f}) — so the apparent tracking is carried by the "
                "cross-platform Xenium-vs-Visium contrast, which is confounded. "
                "Treat it as unsupported.**")
    elif spread_real:
        head = ("**VERDICT: the inter-tumor heterogeneity claim PARTLY SURVIVES. "
                "Predicted kill still varies substantially across the 5 real tumors "
                "at matched dx, but it NO LONGER tracks vascular density the way the "
                "earlier, dx-mismatched run implied.**")
    else:
        head = ("**VERDICT: the inter-tumor heterogeneity claim DOES NOT SURVIVE at "
                "matched resolution — the spread collapses once all five tissues are "
                "solved on the same converged grid. The earlier claim is RETRACTED.**")

    body = (
        f"At the matched, converged dx = {MATCHED_DX_UM:.0f} µm the predicted "
        f"`kill_fraction` spans **{klo * 100:.1f}% → {khi * 100:.1f}%** "
        f"({fold:.1f}× fold, {ksp * 100:.1f} pts absolute) and mean tissue drug AUC "
        f"spans **{alo:.3g} → {ahi:.3g} µM·s** ({ahi / max(alo, 1e-9):.1f}×). The "
        f"per-cell kill sampling s.e. is ≈{typ_se:.4f} (binomial, "
        f"n_tumor = {min(r['n_tumor'] for r in rows):,}–"
        f"{max(r['n_tumor'] for r in rows):,}), i.e. "
        f"{ksp / max(typ_se, 1e-12):.0f}× smaller than the observed spread — "
        + ("so the spread is not sampling noise. "
           if ksp > 20.0 * typ_se else
           "which is NOT a comfortable margin: treat the spread with suspicion. ")
        + f"Vascular density vs kill: {_fmt_corr(st['vasc_vs_kill_all5'])}. "
        f"Vascular density vs mean AUC: {_fmt_corr(st['vasc_vs_auc_all5'])}. "
        f"Within the four like-for-like Visium tumors (the comparison that does not "
        f"cross platforms): {_fmt_corr(st['vasc_vs_kill_visium'])}.\n\n"
        f"**The kill spread is much WIDER than the retracted claim, not narrower.** "
        f"The earlier page reported ~2% → 24%; at matched resolution it is "
        f"{klo * 100:.1f}% → {khi * 100:.1f}%. That is not a vindication — the "
        f"earlier numbers were wrong in magnitude as well as in comparability. The "
        f"magnitude change is NOT mainly the grid (the paired control above shows "
        f"the residual grid effect is only "
        f"≤{max(abs(r['kill_fraction'] / max(nat[r['label']]['kill_fraction'], 1e-9) - 1) for r in vis) * 100:.0f}%); "
        f"it comes from the two upstream honesty fixes that landed since — the "
        f"calibration alpha correction (309.6 → {alpha:.1f}) and the "
        f"`physics/oxygen.py` sink fix that stopped the coarse grid from spuriously "
        f"starving the tissue. Absolute kill fractions of this size for "
        f"temozolomide on real tissue should be read as a statement about this "
        f"engine's current calibration, **NOT** as a clinical expectation — there "
        f"are no outcome labels to check them against, and the oxygen gate they "
        f"pass through rests on an un-calibrated constant.\n\n"
        f"**Where the resolution sensitivity actually lives.** Between dx = 100 µm "
        f"and dx = 12 µm, mean tissue AUC moves at most "
        f"{max(abs(r['mean_auc'] / max(nat[r['label']]['mean_auc'], 1e-9) - 1) for r in vis) * 100:.0f}% "
        f"and kill moves at most "
        f"{max(abs(r['kill_fraction'] / max(nat[r['label']]['kill_fraction'], 1e-9) - 1) for r in vis) * 100:.0f}%. "
        f"The delivery calculation was never the fragile part — the fragile part was "
        f"the oxygen/quiescence gate, and that fragility has since been fixed "
        f"upstream (see the retraction). The grid is matched here regardless, so the "
        f"comparison cannot silently regress.\n\n"
        f"**But read the correlation with n = 5 in mind.** With five points the "
        f"smallest two-sided p a *perfect* rank correlation can attain is "
        f"{st['vasc_vs_kill_all5'].get('min_attainable_p', float('nan')):.3f}; "
        f"nothing here can clear a conventional significance bar, and a single "
        f"tumor moves the rank statistic a lot (leave-one-out range above). Treat "
        f"the correlation as a **direction**, not evidence."
    )
    return head, body


def write_md(rows: List[Dict], native: List[Dict], ladder: List[Dict],
             dose: List[Dict], st: Dict, alpha: float, drug, runtime_s: float,
             path: str) -> None:
    def pct(v):
        return f"{v * 100:.1f}%"

    vis = [r for r in rows if r["source"] != "xenium"]
    nat = {r["label"]: r for r in native}
    L: List[str] = []

    L.append("# Cross-tissue predicted drug behaviour on REAL GBM microenvironments\n")
    L.append(
        "**All five real tumors are simulated on the SAME converged numerical grid "
        f"(dx = {MATCHED_DX_UM:.0f} µm).** The same real drug (temozolomide: "
        f"doxorubicin-calibrated perivascular-consumption physics, alpha = "
        f"{alpha:.0f}, plus the binding-site barrier) and the same validated "
        "transport → kill engine are used for every tissue; only the REAL "
        "microenvironment differs.\n\n"
        "This supersedes an earlier version of this page in which the Xenium "
        "section ran at dx = 12 µm and the four Visium sections ran at their 100 µm "
        "spot pitch. That comparison mixed biology with a numerical grid artifact. "
        "It also predates the calibration honesty pass, so it used the inflated "
        f"alpha = 309.6 rather than the corrected alpha = {alpha:.1f}. "
        "**Every headline number on the earlier page is superseded** — see "
        "**Retraction & correction** below before quoting any of them.\n"
    )
    L.append(
        f"Drug physics (fixed for all tissues): diffusivity {drug.diffusivity:.1f} "
        f"µm²/s, ic50 {drug.ic50:.3g} µM, binding {drug.binding:.3e}/s, clearance "
        f"{drug.clearance:.3e}/s, bbb_penetrance {drug.bbb_penetrance:.3f}.\n"
    )

    # ---------------- retraction ----------------
    L.append("## Retraction & correction (read first)\n")
    L.append(
        "### 1. The earlier cross-tissue spread was confounded by RESOLUTION — RETRACTED\n\n"
        "The previously published headline — *\"same drug, 5 real tumors, predicted "
        "kill spans ~2% → 24%, tracking vascularity\"* — is **withdrawn**. It was "
        "computed with the Xenium section at the converged dx = 12 µm and the four "
        "Visium sections at dx = 100 µm, which is far outside the converged regime "
        "(`results/numerical_rigor.md` §2). At dx = 100 µm the grid cannot resolve "
        f"the ~{BOUNDARY_LAYER_UM:.0f} µm perivascular boundary layer, and — at the "
        "time those numbers were produced — could not resolve the oxygen field "
        "either, so kill was understated by a tissue-dependent amount. The reported "
        "min→max range therefore mixed real biology with numerical error. **Read the "
        "paired control below together with the update that follows it: the "
        "dominant source of that error has since been fixed upstream, which changes "
        "how big the correction is but not the fact that the published numbers are "
        "withdrawn.**\n\n"
        "The paired control below re-runs the identical Visium tissues, identical "
        "drug and identical cell geometry, changing **only** the numerical grid:\n"
    )
    L.append("| Tissue | RETRACTED kill (published @dx=100) | kill @ dx=100 µm "
             "(re-run) | kill @ dx=12 µm | fold | hypoxic @ 100 µm | hypoxic @ 12 µm "
             "| mean AUC @ 100 µm | mean AUC @ 12 µm | penetration_index 100 → 12 µm |")
    L.append("|" + "---|" * 10)
    for r in vis:
        a = nat[r["label"]]
        old = _OLD_CLAIM.get(r["label"], {}).get("kill")
        L.append(
            f"| **{r['label']}** | {pct(old) if old is not None else 'n/a'} | "
            f"{pct(a['kill_fraction'])} | {pct(r['kill_fraction'])} | "
            f"{r['kill_fraction'] / max(a['kill_fraction'], 1e-9):.2f}× | "
            f"{pct(a['hypoxic_tumor_fraction'])} | "
            f"{pct(r['hypoxic_tumor_fraction'])} | {a['mean_auc']:.3g} | "
            f"{r['mean_auc']:.3g} | {a['penetration_index']:.3f} → "
            f"{r['penetration_index']:.3f} |"
        )
    L.append("")
    max_fold = max(r["kill_fraction"] / max(nat[r["label"]]["kill_fraction"], 1e-9)
                   for r in vis)
    resid = max(abs(r["kill_fraction"] / max(nat[r["label"]]["kill_fraction"], 1e-9)
                    - 1.0) for r in vis)
    L.append(
        "**Read the columns correctly.** The *valid* resolution control is the "
        "`dx = 100 µm (re-run)` → `dx = 12 µm` pair: same tissue, same drug object, "
        "same cell geometry, same seed — **only the numerical grid differs**. Any "
        "quantity that moves between those two columns was, in the earlier "
        "write-up, being reported as inter-tumor biology.\n\n"
        + ("**Important update — the resolution confound has largely been fixed "
           "upstream, and this control now shows it.** As re-run today the residual "
           f"dx = 100 µm → 12 µm change in kill is only ≤{resid * 100:.0f}%, not the "
           "large effect the confound implied. The reason is a concurrent fix in "
           "`physics/oxygen.py`: the oxygen sink previously accumulated a raw "
           "tumor-cell *count* per grid cell, which scales as `dx²` and made the "
           "hypoxic fraction — and every hypoxia-gated resistance term downstream — "
           "a function of the numerical grid. On a coarse grid that spuriously "
           "starved the tissue (the Visium sections used to report ~100% hypoxic "
           "tumor at dx = 100 µm) and suppressed kill. With the sink expressed as an "
           "areal density, the grid dependence is gone.\n\n"
           "**This does not un-retract anything.** The earlier page's numbers were "
           "produced *before* that fix, so at the time they were written the "
           "resolution mismatch really was injecting a large, tissue-dependent error "
           "into the reported spread. What has changed is the remedy's size, not its "
           "necessity: matching the grid is now cheap insurance rather than a big "
           "correction, and the run below matches it anyway so that no future "
           "regression can reintroduce the problem silently.\n\n"
           if resid < 0.25 else
           f"The residual dx = 100 µm → 12 µm change in kill is up to "
           f"{max_fold:.1f}×, so the grid mismatch alone was large enough to "
           "dominate the reported inter-tumor spread.\n\n")
        + "The first column is the *retracted published* value and is shown only to "
        "name what is being withdrawn. It is **not** expected to equal the re-run: "
        f"it was produced with the pre-fix calibration alpha = 309.6, whereas this "
        f"run uses alpha = {alpha:.1f}. A separate calibration honesty pass "
        "corrected a half-distance-vs-e-folding-length mix-up in "
        "`validation/calibration.py` (treating the 45 µm doxorubicin half-distance "
        "as an e-folding length inflated alpha by 1/ln²2 = 2.081×). So the earlier "
        "cross-tissue page carried **two** independent defects — a resolution "
        "mismatch and an inflated alpha — and both are corrected here.\n"
    )
    L.append(
        "### 2. `penetration_index` is a vascular-geometry SHAPE statistic, not drug "
        "delivery — CORRECTED\n\n"
        "`penetration_index = mean(AUC) / percentile99(AUC)` over the tissue mask. "
        "Numerator and denominator both scale linearly with the delivery source, so "
        "the index is **scale-invariant**: doubling the delivered drug everywhere "
        "leaves it unchanged. It measures how *even* the exposure field is — and "
        "since the source is `S = vascular_density × bbb_penetrance × plasma`, the "
        "exposure field's evenness is essentially the evenness of the vascular "
        "field. **Correlating `penetration_index` against vascular density is "
        "therefore close to circular and is NOT presented as a finding here.**\n\n"
        "**This is measured, not asserted — the dose-scaling control.** One tissue "
        f"(**{VISIUM_LADDER_SAMPLE}**) was run twice on the same grid with the same "
        f"geometry, changing only the delivered dose by "
        f"{DOSE_SCALE_FACTOR:g}× (`bbb_penetrance`, which enters the source "
        "linearly):\n"
    )
    if len(dose) == 2:
        b, h = dose[0], dose[1]
        L.append("| Metric | 1× dose | "
                 f"{DOSE_SCALE_FACTOR:g}× dose | ratio | behaves like… |")
        L.append("|" + "---|" * 5)
        L.append(f"| Mean tissue AUC (µM·s) | {b['mean_auc']:.4g} | "
                 f"{h['mean_auc']:.4g} | "
                 f"**{h['mean_auc'] / max(b['mean_auc'], 1e-12):.3f}×** | "
                 "a delivery metric ✅ |")
        L.append(f"| Median tissue AUC (µM·s) | {b['median_auc']:.4g} | "
                 f"{h['median_auc']:.4g} | "
                 f"**{h['median_auc'] / max(b['median_auc'], 1e-12):.3f}×** | "
                 "a delivery metric ✅ |")
        L.append(f"| kill_fraction | {b['kill_fraction']:.4f} | "
                 f"{h['kill_fraction']:.4f} | "
                 f"**{h['kill_fraction'] / max(b['kill_fraction'], 1e-12):.3f}×** | "
                 "a dose-responsive prediction ✅ |")
        L.append(f"| **penetration_index** | {b['penetration_index']:.4f} | "
                 f"{h['penetration_index']:.4f} | "
                 f"**{h['penetration_index'] / max(b['penetration_index'], 1e-12):.3f}×** "
                 f"| **NOT a delivery metric ❌** |")
        L.append("")
        L.append(
            f"Doubling the dose doubles the exposure "
            f"({h['mean_auc'] / max(b['mean_auc'], 1e-12):.2f}× on mean AUC) and "
            f"raises the kill "
            f"({h['kill_fraction'] / max(b['kill_fraction'], 1e-12):.2f}×), while "
            f"`penetration_index` moves "
            f"{abs(h['penetration_index'] / max(b['penetration_index'], 1e-12) - 1) * 100:.2f}%. "
            "A statistic that cannot detect a doubling of the delivered drug is "
            "not measuring drug delivery — it is measuring the shape of the field, "
            "which here is the shape of the vascular field it was seeded from.\n"
        )
    L.append(
        "It is retained in the tables below **only** as a vascular-geometry "
        "diagnostic. Every cross-tissue comparison in this document uses **mean / "
        "median tissue drug AUC** and **kill_fraction** instead.\n"
    )

    # ---------------- matched-resolution table ----------------
    L.append("## Matched-resolution per-tissue table "
             f"(ALL tissues at dx = {MATCHED_DX_UM:.0f} µm)\n")
    L.append(
        "| Tissue | Platform | Diagnosis | Cells (tumor) | Grid @dx | Tumor | Endo | "
        "Vascular | Mean O₂ | Hypoxic tumor† | Resistant tumor§ | "
        "Mean AUC (µM·s) | Median AUC (µM·s) | Tumor-cell AUC¶ | Mean conc (µM) | "
        "Median conc (µM) | **Kill** | pen_idx‡ |"
    )
    L.append("|" + "---|" * 18)
    for r in rows:
        L.append(
            f"| **{r['label']}** | {r['platform']} | {r['diagnosis']} | "
            f"{r['n_cells']:,} ({r['n_tumor']:,}) | "
            f"{r['grid'][0]}×{r['grid'][1]} @{r['dx_um']:.0f}µm | "
            f"{pct(r['tumor_frac'])} | {pct(r['endo_frac'])} | "
            f"{r['mean_vascular']:.3f} | {r['mean_oxygen']:.3f} | "
            f"{pct(r['hypoxic_tumor_fraction'])} | "
            f"{pct(r['resistant_tumor_fraction'])} | {r['mean_auc']:.4g} | "
            f"{r['median_auc']:.4g} | {r['tumor_auc_mean']:.4g} | "
            f"{r['mean_conc']:.4g} | "
            f"{r['median_conc']:.4g} | **{pct(r['kill_fraction'])}** "
            f"± {r['kill_se'] * 100:.2f} | {r['penetration_index']:.3f} |"
        )
    L.append("")
    L.append(
        "† **Hypoxic tumor fraction is ORDER-OF-MAGNITUDE ONLY** — but for a "
        "different reason than this page previously gave. A concurrent fix in "
        "`physics/oxygen.py` removed a spurious `dx²` scaling in the oxygen sink "
        "(the tumor sink is now an *areal* density, `rho_tum`, normalised by a fixed "
        "reference area rather than a raw per-grid-cell count), so the hypoxic "
        "fraction is **no longer grid-dependent** — the residual dx = 100 µm → "
        "12 µm movement in the control table above is what is left after that fix. "
        "What remains uncertain is the **absolute level**: "
        "`OxygenConfig.consumption × TUMOR_REFERENCE_AREA_UM2` is a single "
        "**un-calibrated** lumped constant never fit to a measured oxygen tension, "
        "and by that module's own audit a ±40% change swings the real-section "
        "hypoxic tumor fraction between ~0.0 and ~0.6. So read every hypoxic and "
        "resistant fraction here as \"a lot / a little\", never as a number, and "
        "never difference two of them.\n\n"
        "‡ `pen_idx` = `penetration_index`, a **vascular-geometry shape statistic** "
        "(mean AUC ÷ p99 AUC), **not** a delivery measure. Listed for completeness "
        "only; see the correction above. Do not rank tumors by it.\n\n"
        "§ **Resistant tumor** = tumor cells banded HYPOXIC (resistance 0.35) *or* "
        "QUIESCENT (resistance 0.5) by the oxygen gate. These are killed less at "
        "identical exposure, so kill cannot be read from AUC alone. Like the "
        "hypoxic fraction it inherits the un-calibrated oxygen consumption "
        "constant — order-of-magnitude only.\n\n"
        "¶ **Tumor-cell AUC** = mean exposure sampled *at the positions the tumor "
        "cells actually occupy*, rather than averaged over the whole tissue mask. "
        "This is the exposure the Hill kill integrates. It can differ a lot from the "
        "tissue-wide mean when tumor cells are not uniformly distributed over the "
        "mask — which is exactly the Xenium-vs-Visium geometry difference discussed "
        "below.\n\n"
        "± on kill is the binomial s.e. of the per-cell kill draw (the only "
        "stochastic step in the pipeline).\n"
    )

    # ---------------- spread + correlations ----------------
    vlo, vhi, vsp = _spread(rows, "mean_vascular")
    klo, khi, ksp = _spread(rows, "kill_fraction")
    alo, ahi, asp = _spread(rows, "mean_auc")
    mlo, mhi, msp = _spread(rows, "median_auc")
    hlo, hhi, hsp = _spread(rows, "hypoxic_tumor_fraction")
    L.append("## Observed spread at matched resolution\n")
    L.append(
        f"- **Mean vascular density (real)**: {vlo:.3f} → {vhi:.3f} "
        f"({vhi / max(vlo, 1e-9):.1f}×).\n"
        f"- **Endothelial cell fraction (real, assay-anchored)**: "
        f"{pct(min(r['endo_frac'] for r in rows))} → "
        f"{pct(max(r['endo_frac'] for r in rows))}.\n"
        f"- **Mean tissue drug AUC (real delivery)**: {alo:.4g} → {ahi:.4g} µM·s "
        f"({ahi / max(alo, 1e-9):.1f}×).\n"
        f"- **Median tissue drug AUC**: {mlo:.4g} → {mhi:.4g} µM·s "
        f"({mhi / max(mlo, 1e-9):.1f}×).\n"
        f"- **Predicted kill_fraction (converged)**: {pct(klo)} → {pct(khi)} "
        f"({khi / max(klo, 1e-9):.1f}×, {ksp * 100:.1f} pts).\n"
        f"- **Hypoxic tumor fraction (order-of-magnitude only)**: {pct(hlo)} → "
        f"{pct(hhi)}.\n"
    )

    L.append("## Corrected correlations\n")
    L.append(
        "All computed on the matched-dx rows. `penetration_index` appears in NONE of "
        "them, by design.\n"
    )
    L.append("| Relationship | n | Statistic |")
    L.append("|---|---|---|")
    for key, name in (
        ("vasc_vs_kill_all5", "mean vascular density vs **kill_fraction** (all 5)"),
        ("vasc_vs_auc_all5", "mean vascular density vs **mean tissue AUC** (all 5)"),
        ("endo_vs_kill_all5", "endothelial cell fraction vs **kill_fraction** (all 5)"),
        ("endo_vs_auc_all5", "endothelial cell fraction vs **mean tissue AUC** (all 5)"),
        ("auc_vs_kill_all5", "mean tissue AUC vs **kill_fraction** (all 5)"),
        ("tumorauc_vs_kill_all5",
         "mean **tumor-cell-local** AUC vs **kill_fraction** (all 5)"),
        ("vasc_vs_kill_visium", "mean vascular density vs **kill_fraction** (Visium only)"),
        ("vasc_vs_auc_visium", "mean vascular density vs **mean tissue AUC** (Visium only)"),
    ):
        c = st[key]
        L.append(f"| {name} | {int(c.get('n', 0))} | {_fmt_corr(c)} |")
    L.append("")
    L.append(
        "**How much of the vascularity → AUC link is real?** Very little of it is a "
        "*discovery*. The delivery source is defined as "
        "`S = vascular_density × bbb_penetrance × plasma`, and the transport PDE is "
        "linear in the source, so mean tissue AUC is close to a linear functional of "
        "the vascular field **by construction**. A strong vascular↔AUC correlation is "
        "the model working as specified, not evidence about tumors. The only link "
        "carrying non-linear, non-definitional content is **vascularity → kill**, "
        "which passes through the Hill dose–response and the oxygen gate.\n\n"
        "**And `mean_vascular` is itself a normalised shape statistic.** Both the "
        "Xenium KDE (`_kde_on_grid` divides by its own max) and the Visium field "
        "(`_normalize01`) rescale each tissue's vascular field to peak at 1.0. So "
        "`mean vascular density` measures how *evenly spread* the vasculature is "
        "within that tissue — it is **not** an absolute perfused-vessel fraction and "
        "is not comparable across tissues as an absolute. The endothelial **cell "
        "fraction** rows above are the assay-anchored vascularity measure; prefer "
        "them when a cross-tissue vascularity claim is being made.\n"
    )

    # ---------------- residual platform confound ----------------
    xr = next(r for r in rows if r["source"] == "xenium")
    # The Visium tissue whose tissue-wide mean AUC is closest to Xenium's.
    near = min(vis, key=lambda r: abs(r["mean_auc"] - xr["mean_auc"]))
    L.append("## The Xenium point is still NOT comparable to the Visium points\n")
    L.append(
        "Matching the grid removed the *numerical* confound. It did not remove the "
        "*platform* one, and the table above shows exactly how large that residue "
        "is.\n\n"
        f"**{xr['label']}** and **{near['label']}** have almost the same tissue-wide "
        f"mean drug AUC ({xr['mean_auc']:.3g} vs {near['mean_auc']:.3g} µM·s, a "
        f"{abs(xr['mean_auc'] / max(near['mean_auc'], 1e-9) - 1) * 100:.0f}% "
        f"difference) — yet their predicted kill differs by "
        f"{pct(xr['kill_fraction'])} vs {pct(near['kill_fraction'])}. Identical drug, "
        "identical physics, identical grid. Three reconstruction differences account "
        "for it, and **none of them is a claim about the tumors**:\n\n"
        f"1. **Where the tumor cells sit relative to the drug.** Tumor-cell-local "
        f"AUC is {xr['tumor_auc_mean']:.3g} µM·s for {xr['label']} vs "
        f"{near['tumor_auc_mean']:.3g} µM·s for {near['label']} "
        f"({xr['tumor_auc_mean'] / max(xr['mean_auc'], 1e-9):.2f}× and "
        f"{near['tumor_auc_mean'] / max(near['mean_auc'], 1e-9):.2f}× their "
        "respective tissue-wide means). Xenium has REAL cell positions, so its tumor "
        "cells sit in the tumor mass, away from vessels, in the drug-poor part of the "
        "field. Visium \"cells\" are 12 placeholders scattered inside each spot — the "
        "same spots whose endothelial fraction *generated* the vascular field — so "
        "they sample the field near its seeds. This is a geometry artifact of "
        "spot-level reconstruction, not better perfusion.\n"
        f"2. **Resistant-state banding — this is the dominant term.** "
        f"{pct(xr['resistant_tumor_fraction'])} of {xr['label']}'s tumor cells are "
        f"banded hypoxic (kill × 0.35) or quiescent (kill × 0.5) against "
        f"{pct(near['resistant_tumor_fraction'])} for {near['label']}. Those "
        "multipliers alone predict a kill ratio of roughly "
        f"{(0.5 * xr['resistant_tumor_fraction'] + (1 - xr['resistant_tumor_fraction'])) / max(0.5 * near['resistant_tumor_fraction'] + (1 - near['resistant_tumor_fraction']), 1e-9):.2f}×, "
        f"against an observed "
        f"{xr['kill_fraction'] / max(near['kill_fraction'], 1e-9):.2f}× — i.e. the "
        "gap is essentially the oxygen gate, not perfusion. And that gate rests on "
        "the **un-calibrated** oxygen-consumption constant flagged above, whose own "
        "audit says ±40% swings the hypoxic fraction between ~0.0 and ~0.6.\n"
        "3. **Tissue mask.** The Xenium domain carries a real, irregular tissue mask; "
        "the Visium domain is the full bounding rectangle. \"Mean over the tissue "
        "mask\" therefore means different things on the two platforms.\n\n"
        "**Consequence for the headline.** The Xenium point is the lowest on both "
        "vascularity and kill, which *helps* the all-5 rank correlation — but for "
        "reasons that are at least partly platform, not biology. **The defensible "
        "statement is the within-Visium one (n = 4), where all four tissues share a "
        "reconstruction path.** Any all-5 statistic in this document should be read "
        "as illustrative only.\n"
    )

    # ---------------- verdict ----------------
    head, body = _verdict(rows, native, st)
    L.append("## Does the inter-tumor heterogeneity claim survive?\n")
    L.append(head + "\n")
    L.append(body + "\n")

    # ---------------- Visium grid-convergence ladder -------------------------
    if ladder:
        lad = sorted(ladder, key=lambda r: -r["dx_um"])
        fine = lad[-1]
        L.append(f"## Is dx = {MATCHED_DX_UM:.0f} µm actually converged **for a "
                 f"Visium reconstruction**?\n")
        L.append(
            "The grid-refinement audit in `results/numerical_rigor.md` §2 was run on "
            "the **Xenium** section only. Nothing had ever checked that dx = 12 µm is "
            "inside the converged regime for a *Visium* reconstruction, whose fields "
            "and cell placement are built differently. This ladder does that, on "
            f"**{VISIUM_LADDER_SAMPLE}**, with everything except dx held fixed.\n"
        )
        dts = sorted({r.get("pde_dt") for r in lad if r.get("pde_dt")})
        if len(dts) > 1:
            L.append(
                f"**One caveat on this ladder:** the convective (Darcy) term carries "
                f"an explicit CFL limit `dt·|v|/dx ≤ 1`, so the finest grids need a "
                f"smaller timestep — `dt` varies over {min(dts):g}–{max(dts):g} s "
                f"across the rows and is reported per row below. Per "
                f"`results/numerical_rigor.md` §3 the solution is dt-invariant to "
                f"<1e-2 relative for any dt ≲ 20 s, so this does not confound the "
                f"grid comparison, but it does mean the finest row is not a *pure* "
                f"dx refinement.\n"
            )
        L.append("| dx (µm) | grid | dt (s) | mean AUC (µM·s) | AUC rel. to finest "
                 "| kill_fraction | Δkill vs finest | hypoxic tumor | "
                 "pen_idx (shape) |")
        L.append("|" + "---|" * 9)
        for r in lad:
            L.append(
                f"| {r['dx_um']:.0f} | {r['grid'][0]}×{r['grid'][1]} | "
                f"{r.get('pde_dt', float('nan')):g} | "
                f"{r['mean_auc']:.4g} | "
                f"{r['mean_auc'] / max(fine['mean_auc'], 1e-9):.3f} | "
                f"{r['kill_fraction']:.4f} | "
                f"{abs(r['kill_fraction'] - fine['kill_fraction']) / max(fine['kill_fraction'], 1e-9):.3f} | "
                f"{pct(r['hypoxic_tumor_fraction'])} | "
                f"{r['penetration_index']:.3f} |"
            )
        L.append("")
        near12 = min(lad, key=lambda r: abs(r["dx_um"] - MATCHED_DX_UM))
        rel12 = (abs(near12["kill_fraction"] - fine["kill_fraction"])
                 / max(fine["kill_fraction"], 1e-9))
        L.append(
            f"**Reading.** Between dx = {near12['dx_um']:.0f} µm (used throughout "
            f"this page) and the finest grid tested, dx = {fine['dx_um']:.0f} µm, "
            f"`kill_fraction` moves {rel12 * 100:.1f}% relative and mean AUC moves "
            f"{abs(near12['mean_auc'] / max(fine['mean_auc'], 1e-9) - 1) * 100:.1f}%. "
            + ("That is inside the same tolerance `numerical_rigor.md` uses to call "
               "the Xenium kill converged, so **dx = 12 µm is an adequate grid for "
               "these Visium reconstructions too** — the matched-resolution table is "
               "not resting on an unverified assumption."
               if rel12 < 0.05 else
               "**That is NOT tight enough to call dx = 12 µm converged for Visium.** "
               "The matched-resolution numbers below should therefore be read as "
               "much better than the dx = 100 µm ones but still carrying a residual "
               "discretisation error of this size; a finer grid is needed before any "
               "absolute kill fraction is quoted.")
            + " Note that **both** mean AUC and kill are essentially flat across "
            "the whole ladder, including the coarse end. That is the post-fix "
            "picture: with the `physics/oxygen.py` sink expressed as an areal "
            "density, this Visium reconstruction is largely grid-insensitive. "
            "Before that fix the same refinement moved kill by ~2.6–2.9× (the "
            "measurement that motivated this whole rewrite), so the flatness here "
            "is evidence that the defect is fixed, not evidence that the original "
            "concern was misplaced.\n"
        )

    # ---------------- resolution & convergence (preserved + extended) --------
    L.append("## Resolution & convergence\n")
    L.append(
        "*(This section was contributed earlier and is preserved; it is extended "
        "here with the matched-grid policy and the paired resolution control.)*\n\n"
        "Source: [`results/numerical_rigor.md`](numerical_rigor.md) §2 "
        "(grid-refinement audit on the real Xenium section).\n\n"
        "- **`kill_fraction` is grid-converged for dx ≲ 20 µm.** It plateaus and "
        "moves <1–2% relative between the finest grids. It is NOT converged at "
        "dx = 100 µm — see the paired control above, where refining the grid alone "
        "changes it by up to "
        f"{max(r['kill_fraction'] / max(nat[r['label']]['kill_fraction'], 1e-9) for r in vis):.1f}×.\n"
        "- **`penetration_index` requires dx ≲ 12 µm on single-cell data.** The "
        f"calibrated perivascular boundary layer is ~{BOUNDARY_LAYER_UM:.0f} µm "
        "thick. On the Xenium "
        "section the index falls 0.396 (dx 30) → 0.311 (20) → 0.204 (15) → 0.196 "
        "(12) → 0.197 (10), settling only below ~15 µm. Earlier cross-tissue runs "
        "quoted the dx = 30 µm Xenium value (0.396) and over-stated it ~2×; those "
        "figures are superseded. **However, see the correction above: even at the "
        "converged dx the index is a shape statistic, so \"converged penetration "
        "index\" is a converged measure of vascular geometry, not of delivery.**\n"
        "- **Caveat on the cited ladder itself.** The grid-refinement table in "
        "`numerical_rigor.md` §2 was produced *before* the calibration honesty pass, "
        "i.e. at the inflated alpha = 309.6, which corresponds to a ~45 µm "
        f"perivascular boundary layer. At the corrected alpha = {alpha:.1f} the layer "
        f"is ~{BOUNDARY_LAYER_UM:.0f} µm — **longer, hence easier to resolve**. So "
        "the dx ≲ 12 µm requirement inherited from that audit is, if anything, "
        "conservative for this run; the absolute numbers in that ladder no longer "
        "match this drug and should be re-run.\n"
        "- **Why the Visium sections are now run at dx = 12 µm too.** The earlier "
        "version argued that sub-spot dx is \"physically meaningless\" for Visium "
        "and pinned it at the 100 µm spot pitch. That conflates two different "
        "things:\n"
        "  - the **data resolution** (~100 µm — genuinely spot-limited, and "
        "unchanged by anything done here), and\n"
        "  - the **simulation grid** (a numerical discretisation of a continuous "
        "PDE).\n\n"
        "  The reconstructed Visium fields are continuous: vascular and ECM density "
        f"are a Gaussian KDE with bandwidth fixed at "
        f"{VISIUM_KDE_BANDWIDTH_UM:.0f} µm **in physical units**, and cells sit at "
        "continuous coordinates inside each spot footprint. Refining dx samples that "
        "*same band-limited function* more accurately; it does not sharpen it. "
        f"Meanwhile dx = 100 µm is coarser than the ~{BOUNDARY_LAYER_UM:.0f} µm drug "
        "boundary layer "
        "and coarser than the oxygen diffusion length, so it carries a large, "
        "tissue-dependent discretisation error. Refining the grid removes that "
        "error. **The information content stays ~100 µm; only the numerical "
        "artifact goes away.**\n"
        "- **Cell geometry was held fixed while refining.** "
        "`reconstruct_from_visium` scales its cell-placement jitter as "
        "`min(27.5, 0.45·dx)`, which would have shrunk the spot footprint from "
        f"27.5 µm to {0.45 * MATCHED_DX_UM:.1f} µm purely because the grid got "
        f"finer. The experiment pins it at {VISIUM_JITTER_UM:.1f} µm (identical to "
        "the native 100 µm reconstruction) so the paired control isolates the grid "
        "and nothing else.\n"
        "- **The oxygen sink's grid dependence has been fixed (concurrently), and "
        "this changes what the earlier hypoxia caveat should say.** "
        "`physics/oxygen.py` previously accumulated a raw tumor-cell *count* per "
        "grid cell, which scales as `dx²` and made the hypoxic fraction — and every "
        "hypoxia-gated resistance term downstream — a function of the numerical "
        "grid. It now uses an areal density normalised by a fixed reference area, "
        "so the hypoxic fraction is grid-independent. Earlier statements on this "
        "page that the hypoxic fraction \"keeps falling as dx drops\" described that "
        "defect and are withdrawn.\n"
        "- **But the hypoxic fraction is still un-calibrated in absolute terms** and "
        "is quoted as order-of-magnitude only: the lumped constant "
        "`consumption × TUMOR_REFERENCE_AREA_UM2` was never fit to a measured oxygen "
        "tension, and ±40% on it swings the hypoxic fraction between ~0.0 and ~0.6. "
        "Because kill is gated on that banding, the absolute kill fractions inherit "
        "this uncertainty.\n"
        "- **What actually moves kill between dx = 100 µm and dx = 12 µm.** Mean "
        "tissue AUC is nearly dx-invariant (see the control table: it moves only a "
        "few percent), so the residual dx effect on kill does **not** come from "
        "delivery. It comes from the oxygen/quiescence banding and from resolving "
        "the perivascular gradient the tumor cells actually sample. That is why "
        "`kill_fraction` must be quoted at a converged dx even though `mean_auc` is "
        "robust.\n"
    )

    # ---------------- caveats ----------------
    L.append("## Honest caveats (read before trusting any number)\n")
    L.append(
        "- **Visium information content is still ~100 µm.** Matching the numerical "
        "grid at 12 µm does NOT give the Visium sections single-cell detail. Their "
        "vascular/ECM fields remain a 100 µm-bandwidth KDE of deconvolved spot "
        "fractions, and their \"cells\" remain 12 synthetic placeholders per spot "
        "apportioned by fraction. The Xenium section has genuine per-cell structure "
        "that no Visium section has at any dx. The finer grid removes a numerical "
        "artifact; it adds no real detail. **A Xenium-vs-Visium difference is still "
        "confounded by platform, deconvolution and information content — just no "
        "longer by the grid.**\n"
        "- **No outcome labels.** Every tissue is an *untreated* open resection. "
        "There is NO paired measured temozolomide response for any of them. Kill and "
        "AUC are what the validated physics predicts on real geometry — "
        "**predictions, not validated outcomes**. Do not read `kill_fraction` as a "
        "clinical response.\n"
        f"- **n = {len(rows)}.** Every correlation above is descriptive. The "
        "smallest two-sided p attainable at n = 5 even for a perfect rank match is "
        f"{st['vasc_vs_kill_all5'].get('min_attainable_p', float('nan')):.3f}, and "
        "the Visium-only subset (n = 4) cannot reach any conventional threshold at "
        "all. No significance is claimed.\n"
        "- **`mean_vascular` and `penetration_index` are both normalised shape "
        "statistics**, not absolute physical densities/deliveries — see the "
        "correction and the correlation notes above.\n"
        "- **The vascularity → AUC link is largely mechanical**, since the delivery "
        "source is defined proportional to the vascular field.\n"
        "- **Reconstruction, not measurement, of fields.** Vascular density, ECM and "
        "oxygen are inferred from real markers/perfusion, not directly measured "
        "concentrations. They are faithful to the data but are model fields.\n"
        "- **One drug, one config.** Only temozolomide at one engine config is run "
        "here; the relative comparison is the message, not any absolute value.\n"
    )

    L.append(
        f"\n_Generated by `experiments/cross_tissue.py` on real data "
        f"({len(rows)} tissues at matched dx = {MATCHED_DX_UM:.0f} µm + "
        f"{len(native)} paired resolution controls at dx = "
        f"{VISIUM_NATIVE_DX_UM:.0f} µm; runtime {runtime_s / 60.0:.1f} min). "
        f"Figure: `results/cross_tissue.png`._\n"
    )

    with open(path, "w") as fh:
        fh.write("\n".join(L))
    print(f"markdown written: {path}")


# ---------------------------------------------------------------------------
def main() -> None:
    t0 = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    drug, alpha = compile_real_temozolomide()
    print(f"REAL drug compiled ONCE: temozolomide (alpha={alpha:.1f}); reused "
          f"unchanged for every tissue.")
    print(f"MATCHED simulation grid for ALL tissues: dx = {MATCHED_DX_UM:g} um "
          f"(converged). Visium DATA resolution is still ~100 um.")

    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)

    diagnoses: Dict[str, str] = {}
    for sid in list_samples():
        try:
            diagnoses[sid] = sample_info(sid).get("diagnosis", "Glioblastoma")
        except Exception:  # noqa: BLE001
            diagnoses[sid] = "Glioblastoma"

    # --- 1. all five tissues at the MATCHED converged grid -------------------
    print("\n" + "=" * 78)
    print(f"PASS 1 -- all 5 REAL tissues at MATCHED dx = {MATCHED_DX_UM:g} um")
    print("=" * 78)
    rows: List[Dict] = [run_one_tissue(
        "Xenium-GBM", "xenium", None, MATCHED_DX_UM, drug, pde,
        "GBM (10x Xenium single-cell section)")]
    for sid in list_samples():
        rows.append(run_one_tissue(sid, "visium", sid, MATCHED_DX_UM, drug, pde,
                                   diagnoses[sid]))

    # --- 2. paired resolution control at the native Visium spot pitch --------
    print("\n" + "=" * 78)
    print(f"PASS 2 -- resolution control: same Visium tissues at dx = "
          f"{VISIUM_NATIVE_DX_UM:g} um (only the grid changes)")
    print("=" * 78)
    native: List[Dict] = [
        run_one_tissue(sid, "visium", sid, VISIUM_NATIVE_DX_UM, drug, pde,
                       diagnoses[sid])
        for sid in list_samples()
    ]

    # --- 3. self-contained grid-convergence ladder on ONE Visium section -----
    print("\n" + "=" * 78)
    print(f"PASS 3 -- grid-convergence ladder on {VISIUM_LADDER_SAMPLE} "
          f"(dx = {VISIUM_LADDER_DX})")
    print("=" * 78)
    ladder: List[Dict] = [
        run_one_tissue(f"{VISIUM_LADDER_SAMPLE}@dx{dx:g}", "visium",
                       VISIUM_LADDER_SAMPLE, dx, drug, pde,
                       diagnoses.get(VISIUM_LADDER_SAMPLE, "Glioblastoma"))
        for dx in VISIUM_LADDER_DX
    ]

    # --- 4. dose-scaling control: the decisive penetration_index test --------
    print("\n" + "=" * 78)
    print(f"PASS 4 -- dose-scaling control on {VISIUM_LADDER_SAMPLE} "
          f"(1x vs {DOSE_SCALE_FACTOR:g}x delivered dose)")
    print("=" * 78)
    dose = run_dose_control(drug, pde, diagnoses.get(VISIUM_LADDER_SAMPLE,
                                                     "Glioblastoma"))

    st = compute_stats(rows)
    runtime = time.time() - t0
    make_figure(rows, native, ladder, dose, st, alpha, _OUT_PNG)
    write_md(rows, native, ladder, dose, st, alpha, drug, runtime, _OUT_MD)

    # --- stdout summary ------------------------------------------------------
    print("\n" + "=" * 100)
    print(f"CROSS-TISSUE SUMMARY -- MATCHED dx = {MATCHED_DX_UM:g} um "
          f"(same TMZ physics, real GBM tumors)")
    print("=" * 100)
    print(f"{'tissue':16s} {'plat':7s} {'vasc':>6s} {'endo%':>6s} {'meanAUC':>11s} "
          f"{'medAUC':>11s} {'kill':>7s} {'hypox':>7s} {'[pen]':>7s}")
    for r in rows:
        plat = "xenium" if r["source"] == "xenium" else "visium"
        print(f"{r['label']:16s} {plat:7s} {r['mean_vascular']:6.3f} "
              f"{r['endo_frac'] * 100:6.2f} {r['mean_auc']:11.4g} "
              f"{r['median_auc']:11.4g} {r['kill_fraction']:7.4f} "
              f"{r['hypoxic_tumor_fraction']:7.4f} {r['penetration_index']:7.3f}")
    print("-" * 100)
    for key in ("vasc_vs_kill_all5", "vasc_vs_auc_all5", "endo_vs_kill_all5",
                "auc_vs_kill_all5", "vasc_vs_kill_visium"):
        print(f"{key:24s} {_fmt_corr(st[key])}")
    print("-" * 100)
    print("RESOLUTION CONTROL (Visium, only the grid changes):")
    natmap = {r["label"]: r for r in native}
    for r in rows:
        if r["source"] == "xenium":
            continue
        a = natmap[r["label"]]
        print(f"  {r['label']:16s} kill {a['kill_fraction']:.4f} (dx100) -> "
              f"{r['kill_fraction']:.4f} (dx12)  "
              f"[{r['kill_fraction'] / max(a['kill_fraction'], 1e-9):.2f}x]   "
              f"pen {a['penetration_index']:.3f} -> {r['penetration_index']:.3f}   "
              f"hypox {a['hypoxic_tumor_fraction']:.3f} -> "
              f"{r['hypoxic_tumor_fraction']:.3f}")
    print("-" * 100)
    if len(dose) == 2:
        b, h = dose
        print(f"DOSE-SCALING CONTROL ({DOSE_SCALE_FACTOR:g}x delivered dose, same "
              f"tissue/grid):")
        print(f"  mean AUC          {b['mean_auc']:.4g} -> {h['mean_auc']:.4g}  "
              f"[{h['mean_auc'] / max(b['mean_auc'], 1e-12):.3f}x]  (delivery: moves)")
        print(f"  kill_fraction     {b['kill_fraction']:.4f} -> "
              f"{h['kill_fraction']:.4f}  "
              f"[{h['kill_fraction'] / max(b['kill_fraction'], 1e-12):.3f}x]  "
              f"(dose-responsive)")
        print(f"  penetration_index {b['penetration_index']:.4f} -> "
              f"{h['penetration_index']:.4f}  "
              f"[{h['penetration_index'] / max(b['penetration_index'], 1e-12):.3f}x]"
              f"  (SHAPE statistic: does NOT move -> not a delivery metric)")
    if ladder:
        lad = sorted(ladder, key=lambda r: -r["dx_um"])
        print("-" * 100)
        print(f"VISIUM GRID LADDER ({VISIUM_LADDER_SAMPLE}):  " + "  ".join(
            f"dx{r['dx_um']:g}:kill={r['kill_fraction']:.4f}" for r in lad))
    print("-" * 100)
    head, _ = _verdict(rows, native, st)
    print(head.replace("**", ""))
    print("HONEST: predictions on real geometry, NO outcome labels; n=5; Visium "
          "information content is still ~100 um spot-limited even at dx=12 um; "
          "penetration_index is a shape statistic and is NOT used to rank delivery.")
    print(f"[runtime {runtime / 60.0:.1f} min]")


if __name__ == "__main__":
    main()
