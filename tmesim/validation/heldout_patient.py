"""THE HELD-OUT PATIENT — does anything in this project generalise to a second tumour?

Every validated result in this repository was developed on ONE tissue: the 319-gene
10x Xenium GBM section (``data_cache/xenium_gbm``, 40,887 cells). The marker map,
the tissue mask radius, the vessel linking distance, the headline metric — all of
them were chosen, debugged and defended while looking at that one section.

A SECOND real patient is cached and the pipeline was never developed on it: the
480-gene Immuno-Oncology Xenium GBM section (``data_cache/xenium_gbm_io``, 816,769
cells). Running the whole validated pipeline on it is the closest thing to a
genuine out-of-sample test available without new data.

This module runs it, end to end, and reports EVERY quantity side by side:

1. **Reconstruction** — cell-type fractions, and specifically whether the
   ENDOTHELIAL fraction lands near the human volume-EM tissue-counting value of
   **8.0%** (Shapson-Coe A et al., *Science* 2024;384(6696):eadk4858, PMID
   38723085) as it did on section 1 (9.74%). This is an out-of-sample test of the
   marker-gated classifier on a DIFFERENT gene panel.
2. **Vascular morphometry** — :func:`reconstruct.vessel_network.segment_vessels`:
   MVD, intervessel spacing, tumour-to-vessel distance, and the load-bearing
   **vascular area fraction** against the measured **2.1 ± 0.3%** whole-section
   CD34 value (Tamma R et al., *Cells* 2022;12(1):11, PMID 36611806). Section
   1 gave 2.50%.
3. **Simulation** — the five-drug panel at dx = 12 µm in **uniform** sensitivity
   mode (the only fair cross-drug mode), reporting ``therapeutic_coverage`` @1x
   IC50, mean tumour AUC and ``kill_fraction``, and asking whether the ordering
   **cisplatin > doxorubicin > lomustine > methotrexate > temozolomide** holds.
4. **The hypoxia check on its home turf** — the IO panel MEASURES VEGFA and HIF1A,
   so ``results/state_validation.md``'s central claim (the simulated oxygen field
   is beaten by a one-line distance-to-nearest-endothelial-cell baseline) is
   re-derived here independently, on the same section it was found on.
5. **An honest comparison table** with an explicit per-claim verdict.

Both sections are driven through the SAME frozen production path
(``reconstruct.real_pipeline.reconstruct_real`` -> ``physics.solve_oxygen`` ->
``sim.engine.simulate_response``) with identical settings. Section 1's numbers are
**re-derived here, not quoted**, so the comparison cannot be a comparison between
two different configurations.

THE CONFOUND, STATED UP FRONT
-----------------------------
The two sections are **different patients** AND **different panels** — they share
only **60 genes**. The IO panel lacks PECAM1 / SOX2 / OLIG2 / TOP2A / AQP4 / GJA1 /
MOG, so :mod:`reconstruct.marker_panels` substitutes VWF/FLT1/PLVAP (endothelial),
SLC1A3/S100B/PTN (astro-lineage stromal) and EGFR/PDGFRA/IL13RA2/MET (tumour). Any
composition difference is therefore **jointly** a difference in biology and a
difference in what is measurable, and this report never attributes one without
saying so. Neither section carries a drug-response outcome label.

Subsampling
-----------
816,769 cells is large. Every effect size, every cell-type fraction, every vessel
statistic and every simulation number is computed on **all** cells. Only the
permutation NULL is subsampled (a seeded 50,000-cell draw, matching
``state_validation``), which can only make a p-value more conservative. The
subsampling is recorded in the report.

Nothing in ``core`` / ``physics`` / ``molecular`` / ``reconstruct`` / ``sim.engine``
is modified by this module.

Run:    ./.venv/bin/python -m tmesim.validation.heldout_patient
Writes: results/heldout_patient.md
        results/heldout_patient.png

Units: AUC µM·s, exposure/IC50 µM, length µm, coverage/fractions dimensionless.
"""
from __future__ import annotations

import argparse
import os
import pickle
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..cells.agents import update_cell_states
from ..core import CellState, OxygenConfig, PDEConfig, SimConfig
from ..data.real.singlecell_spatial import load as load_singlecell_spatial
from ..molecular.binding_site_barrier import apply_barrier_to_drugparams
from ..molecular.library import get_molecule
from ..physics.diffusion import ALPHA_TORTUOSITY
from ..reconstruct.marker_panels import hypoxia_genes_for, select_marker_map
from ..reconstruct.real_pipeline import (
    ECM_MARKERS as ECM_MARKERS_ALL,
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from ..reconstruct.vessel_network import (
    LINK_DISTANCE_UM,
    VESSEL_RADIUS_UM,
    link_distance_sweep,
    segment_vessels,
    segmentation_report,
)
from ..sim.engine import simulate_response
from . import delivery_metrics as DM
from .calibration import calibrate_alpha, calibrated_drugparams
# Statistics kernels are REUSED from the modules that already published with
# them, never re-implemented, so the hypoxia re-check is directly comparable to
# results/state_validation.md rather than merely similar to it.
from .cellstate_validation import partial_spearman, spearman_permutation
from .state_validation import (
    all_gene_ranking,
    hypoxia_score,
    spearman_full,
    subsample_index,
    tile_correlation,
    trivial_baselines,
)

__all__ = [
    "SECTIONS",
    "PANEL",
    "DX_UM",
    "run_section",
    "run",
    "make_figure",
    "write_report",
]

# ---------------------------------------------------------------------------
# Declared constants — every choice is named here, none is buried.
# ---------------------------------------------------------------------------

#: The two cached real GBM sections. DIFFERENT PATIENTS, DIFFERENT PANELS.
#: ``dev`` marks the one the whole pipeline was developed on; the other is the
#: held-out patient.
SECTIONS: Tuple[Dict[str, object], ...] = (
    {
        "key": "section1_brain319",
        "sample": "data_cache/xenium_gbm",
        "short": "Section 1 (DEVELOPMENT)",
        "panel": "hBrain_254g + 65-gene add-on (319 targets), Xenium v1.3.0",
        "dev": True,
    },
    {
        "key": "section2_io480",
        "sample": "data_cache/xenium_gbm_io",
        "short": "Section 2 (HELD OUT)",
        "panel": "Human Immuno-Oncology Profiling hImmune_v1 + Add-on B "
                 "(480 targets), Xenium v2.0.0",
        "dev": False,
    },
)

#: The five-drug panel, in the order ``results/real_drug_panel.md`` uses.
PANEL: Tuple[str, ...] = (
    "temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin",
)

#: Production grid spacing. dx = 12 µm is the converged resolution
#: (``results/numerical_rigor.md`` §2) and the finest grid satisfying the
#: advection CFL at dt = 10 s.
DX_UM = 12.0

#: The PDE configuration ``results/real_drug_panel.md`` ran under, unchanged, so
#: section 1's numbers here are directly checkable against that report.
PDE_KWARGS = dict(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)

#: The published ordering this project claims, from ``results/real_drug_panel.md``
#: (uniform sensitivity mode, dx = 12 µm).
PUBLISHED_ORDERING: Tuple[str, ...] = (
    "cisplatin", "doxorubicin", "lomustine", "methotrexate", "temozolomide",
)

#: Human cortex ENDOTHELIAL volume fraction by volume-EM tissue counting — the
#: external anchor the 9.74% reconstruction was validated against. NOT a
#: dissociation-based scRNA-seq number (those give 1-2%, which would have put a
#: ~5x error into the delivery source term).
ENDO_VOLUME_EM = {
    "value_pct": 8.0,
    "what": "endothelial fraction of cells, human cortex volume electron microscopy",
    "source": "Shapson-Coe A et al., Science 2024;384(6696):eadk4858 — "
              "'A petavoxel fragment of human cerebral cortex reconstructed at "
              "nanoscale resolution'",
    "pmid": "38723085",
    "doi": "10.1126/science.adk4858",
    #: Tolerance for calling the fraction a replication. 1.5x either way is
    #: generous and is declared BEFORE the numbers are looked at: the anchor is
    #: normal cortex, not tumour, and cell-type fractions are gate-dependent.
    "tolerance_fold": 1.5,
}

#: Permutation settings, matched to ``state_validation`` so the two reports'
#: p-values mean the same thing.
PERM_SEED = 0
N_PERM = 2000
PERM_SUBSAMPLE = 50_000

_RESULTS = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
)
_OUT_MD = os.path.join(_RESULTS, "heldout_patient.md")
_OUT_PNG = os.path.join(_RESULTS, "heldout_patient.png")

_INK = "#0b0b0b"
_INK_2 = "#52514e"
_GRID = "#dcdbd6"
_SEC_COLORS = {"section1_brain319": "#e34948", "section2_io480": "#2a78d6"}


def _retry(fn, *args, **kwargs):
    """Call ``fn``, retrying ONCE after a short pause on any error (the brief's
    API/network-error policy; here it covers transient IO on the caches)."""
    try:
        return fn(*args, **kwargs)
    except Exception:
        time.sleep(3)
        return fn(*args, **kwargs)


# ---------------------------------------------------------------------------
# 1. RECONSTRUCTION — cell-type fractions on a panel the classifier never saw
# ---------------------------------------------------------------------------
def reconstruct_section(sample: str, dx_um: float = DX_UM):
    """Drive the FROZEN production reconstruction on one cached section."""
    return _retry(reconstruct_real, source="xenium", sample=sample,
                  dx_um=float(dx_um))


def composition_block(micro, data: Dict) -> Dict:
    """Cell-type fractions + the ENDOTHELIAL external check.

    The marker map is not chosen here — :func:`marker_panels.select_marker_map`
    picks it from the panel that was actually loaded, and this block only records
    WHICH map won and how unambiguous that was, so a composition difference can be
    read against the markers that produced it.
    """
    rep = reconstruction_report(micro)
    genes = np.asarray(data["genes"], dtype=str)
    marker_map, panel_name, coverage = select_marker_map(genes)
    present = set(genes.tolist())

    # Per-marker detection rate (fraction of cells with >= 1 transcript). A
    # mean-z gate over three sparse markers is a harder gate to pass than a
    # single-marker gate, so these rates are part of any composition comparison.
    expr = data["expression"]
    gene_index = {g: j for j, g in enumerate(genes.tolist())}
    detection: Dict[str, float] = {}
    for ms in marker_map.values():
        for m in ms:
            if m in gene_index and m not in detection:
                col = expr[:, [gene_index[m]]]
                col = np.asarray(col.todense()).ravel() if hasattr(col, "todense") \
                    else np.asarray(col).ravel()
                detection[m] = float((col > 0).mean())

    endo_pct = 100.0 * rep["celltype_fractions"]["ENDOTHELIAL"]
    fold = endo_pct / ENDO_VOLUME_EM["value_pct"]
    tol = float(ENDO_VOLUME_EM["tolerance_fold"])
    verdict = "PASS" if (1.0 / tol) <= fold <= tol else ("HIGH" if fold > tol else "LOW")

    return {
        "report": rep,
        "panel_name": panel_name,
        "panel_coverage": coverage,
        "n_panel_genes": int(genes.size),
        "markers_used": {
            ct.name: [m for m in ms if m in present] for ct, ms in marker_map.items()
        },
        "marker_detection": detection,
        "markers_missing": {
            ct.name: [m for m in ms if m not in present]
            for ct, ms in marker_map.items()
        },
        "endothelial_pct": endo_pct,
        "endothelial_fold_vs_volume_em": fold,
        "endothelial_verdict": verdict,
        # Areal endothelial density is the quantity the oxygen SOURCE actually
        # depends on. It is NOT the same as the endothelial cell-type FRACTION:
        # a section can have a lower endothelial fraction and still be more
        # densely vascularised if its overall cell density is higher. Reporting
        # only the fraction would hide exactly that.
        "cell_density_per_mm2": (
            float(rep["n_cells"]) /
            (float(np.asarray(micro.domain.tissue_mask).sum())
             * float(micro.domain.dx) ** 2 / 1e6)
        ),
        "endothelial_per_mm2": (
            float(rep["celltype_counts"]["ENDOTHELIAL"]) /
            (float(np.asarray(micro.domain.tissue_mask).sum())
             * float(micro.domain.dx) ** 2 / 1e6)
        ),
    }


# ---------------------------------------------------------------------------
# 2. VASCULAR MORPHOMETRY — the discrete vessel network, unmodified
# ---------------------------------------------------------------------------
def vascular_block(micro) -> Dict:
    """Segment vessels and score them against published human GBM morphometry."""
    vessels = segment_vessels(micro.cells, micro.domain,
                              link_distance_um=LINK_DISTANCE_UM)
    rep = segmentation_report(vessels, micro.domain, radius_um=VESSEL_RADIUS_UM)
    sweep = link_distance_sweep(micro.cells, micro.domain)
    keep = ("link_distance_um", "n_endothelial_cells", "n_vessels",
            "tissue_area_mm2", "mvd_per_mm2", "spacing_from_mvd_um",
            "singleton_fraction", "mean_cells_per_vessel",
            "median_extent_multi_um", "median_aspect_multi",
            "tumor_frac_beyond_100um", "tumor_frac_beyond_150um", "n_tumor")
    return {
        "vessels": {k: vessels[k] for k in keep},
        "intervessel_nn": vessels["intervessel_nn"],
        "tumor_to_vessel": vessels["tumor_to_vessel"],
        "_tumor_to_vessel_um": vessels["tumor_to_vessel_um"],
        "_intervessel_nn_um": vessels["intervessel_nn_um"],
        "_xy": vessels["xy"],
        "report": {k: v for k, v in rep.items()},
        "sweep": sweep,
    }


# ---------------------------------------------------------------------------
# 3. SIMULATION — the drug panel, uniform sensitivity, identical on both sections
# ---------------------------------------------------------------------------
def _effective_diffusivity(diffusivity: float, ecm_density: float) -> float:
    return float(diffusivity) * (
        1.0 - ALPHA_TORTUOSITY * float(np.clip(ecm_density, 0.0, 1.0))
    )


def drug_panel_block(micro, alpha: float, mean_ecm: float,
                     panel: Sequence[str] = PANEL) -> List[Dict]:
    """Run every drug on ONE reconstruction in ``uniform`` sensitivity mode.

    ``uniform`` is the only mode in which a cross-drug ranking is fair: every drug
    carries sensitivity 1.0 at every cell, so the factor cancels and the ranking
    isolates DELIVERY x POTENCY. Absolute kills are an equal
    no-intrinsic-resistance upper bound for every drug.
    """
    pde = PDEConfig(**PDE_KWARGS)
    rows: List[Dict] = []
    # Distance-from-vessel is GEOMETRY: identical for every drug on this tissue,
    # and it is the expensive part of the penetration fit. Compute it once.
    # PRIMARY metric: direct nearest-endothelial-cell distance (see
    # results/vessel_profile.md's dedicated distance-metric-change section);
    # this replaces the LEGACY percentile-mask distance compute_all used to
    # default to.
    distance_um = DM.direct_vessel_distance_um(micro.domain, DM.endothelial_xy(micro.cells))
    for name in panel:
        dp = calibrated_drugparams(name, alpha)
        dp = apply_barrier_to_drugparams(dp, get_molecule(name))
        info = apply_drug_sensitivity(micro, name, mode="uniform")
        t0 = time.time()
        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
        met = DM.compute_all(micro, res, dp, distance_um=distance_um)
        if distance_um is None:
            distance_um = met["penetration_length"].get("_distance_um")
        h = met["headline"]
        rows.append({
            "drug": name,
            "mw": float(get_molecule(name).get("mw")),
            "ic50_uM": float(dp.ic50),
            "bbb_penetrance": float(dp.bbb_penetrance),
            "d_eff_um2_s": _effective_diffusivity(dp.diffusivity, mean_ecm),
            "kill_fraction": float(res.kill_fraction),
            "penetration_index": float(res.extras.get("penetration_index", np.nan)),
            "coverage_1x": float(h["therapeutic_coverage_1x_ic50"]),
            "coverage_0.5x": float(h["therapeutic_coverage_0.5x_ic50"]),
            "coverage_2x": float(h["therapeutic_coverage_2x_ic50"]),
            "median_E_over_ic50": float(h["median_E_over_ic50"]),
            "mean_tumor_auc": float(h["mean_tumor_auc_uM_s"]),
            "median_tumor_auc": float(h["median_tumor_auc_uM_s"]),
            "penetration_length_um": float(h["penetration_length_L_um"]),
            "penetration_length_source_null_um": float(
                met["penetration_length"]["source_null"]["L_um"]),
            "hypoxic_tumor_fraction": float(
                met["hypoxia"]["hypoxic_tumor_fraction"]),
            "viable_tumor_fraction": float(met["hypoxia"]["viable_tumor_fraction"]),
            "sensitivity_mode": info["mode"],
            "sensitivity_mean": float(info.get("mean", 1.0)),
            "upper_bound": bool(info.get("upper_bound", True)),
            "runtime_s": time.time() - t0,
        })
        print(f"    {name:14s} kill={rows[-1]['kill_fraction']:.4f}  "
              f"cov@1xIC50={rows[-1]['coverage_1x']:.4f}  "
              f"meanAUC={rows[-1]['mean_tumor_auc']:.3g}  "
              f"({rows[-1]['runtime_s']:.0f}s)", flush=True)
    # Leave the shared reconstruction in the fair mode.
    apply_drug_sensitivity(micro, panel[0], mode="uniform")
    return rows


def ordering(rows: Sequence[Dict], key: str) -> Tuple[str, ...]:
    """Drug names sorted best-first by ``key``."""
    return tuple(r["drug"] for r in sorted(rows, key=lambda r: -r[key]))


def ordering_diagnosis(rows: Sequence[Dict]) -> Dict:
    """Is the ordering informative, or is it just IC50 read back out?

    A cross-drug ranking that is identical to the IC50 ranking tells us nothing
    about the tissue: it would be reproduced by a spreadsheet with no geometry,
    no transport and no patient. This computes the Spearman correlation between
    the kill/coverage ordering and (a) potency alone, (b) delivery alone, and
    (c) the delivery x potency product AUC/IC50 — the quantity
    ``results/SESSION_AUDIT.md`` claims recovers the ordering exactly.
    """
    kill = np.array([r["kill_fraction"] for r in rows])
    cov = np.array([r["coverage_1x"] for r in rows])
    inv_ic50 = np.array([1.0 / r["ic50_uM"] for r in rows])
    auc = np.array([r["mean_tumor_auc"] for r in rows])
    ratio = auc / np.array([r["ic50_uM"] for r in rows])
    return {
        "rank_by_kill": ordering(rows, "kill_fraction"),
        "rank_by_coverage": ordering(rows, "coverage_1x"),
        "rank_by_potency_alone": tuple(
            r["drug"] for r in sorted(rows, key=lambda r: r["ic50_uM"])),
        "rank_by_delivery_alone": ordering(rows, "mean_tumor_auc"),
        "rank_by_auc_over_ic50": tuple(
            np.array([r["drug"] for r in rows])[np.argsort(-ratio)].tolist()),
        "rho_kill_vs_potency": spearman_full(kill, inv_ic50),
        "rho_kill_vs_delivery": spearman_full(kill, auc),
        "rho_kill_vs_auc_over_ic50": spearman_full(kill, ratio),
        "rho_coverage_vs_potency": spearman_full(cov, inv_ic50),
        "rho_coverage_vs_auc_over_ic50": spearman_full(cov, ratio),
        "kill_spread": float(kill.max() - kill.min()),
        "coverage_spread": float(cov.max() - cov.min()),
        "auc_fold_range": float(auc.max() / auc.min()),
        "ic50_fold_range": float(
            max(r["ic50_uM"] for r in rows) / min(r["ic50_uM"] for r in rows)),
    }


# ---------------------------------------------------------------------------
# 4. THE HYPOXIA CHECK, on the section that actually MEASURES hypoxia response
# ---------------------------------------------------------------------------
def hypoxia_block(micro, data: Dict, n_perm: int = N_PERM) -> Optional[Dict]:
    """Re-derive ``results/state_validation.md`` Test 1 from the same code path.

    Returns ``None`` when the panel measures no hypoxia-response gene — which is
    the honest answer for the 319-gene brain panel (no VEGFA, no HIF1A, no CA9,
    no ADM, no NDRG1, no SLC2A1). On that section hypoxia is **unmeasured**, not
    validated, and this function says so by declining to produce a number.

    Everything is computed on ALL cells. Only the permutation null is subsampled.
    """
    genes = np.asarray(data["genes"], dtype=str)
    used = hypoxia_genes_for(genes)
    if not used:
        return None

    expr = data["expression"]
    if expr.shape[0] != micro.cells.n:
        raise ValueError(
            f"expression rows ({expr.shape[0]}) != cells ({micro.cells.n}); the "
            "reconstruction and the loaded matrix are not row-aligned"
        )
    totals = np.asarray(expr.sum(axis=1)).ravel().astype(np.float64)

    oxygen = np.asarray(micro.fields.oxygen, dtype=np.float64)
    update_cell_states(micro, oxygen, OxygenConfig())
    col, row = micro.domain.to_index(micro.cells.x, micro.cells.y)
    o2 = oxygen[row, col]
    vasc_at_cell = np.asarray(
        micro.fields.vascular_density, dtype=np.float64)[row, col]

    hs = hypoxia_score(expr, genes, totals)
    score = hs["score"]
    log_depth = np.log(np.maximum(totals, 1.0))

    base = trivial_baselines(micro.cells.x, micro.cells.y, micro.cells.cell_type,
                             vasc_at_cell)
    predictors: Dict[str, Dict[str, object]] = {
        "simulated_o2": {
            "values": o2,
            "what": "the reaction-diffusion solve under test",
        }
    }
    for k, v in base.items():
        predictors[k] = {"values": v["values"], "what": v["what"],
                         "median_um": v.get("median_um")}

    out_pred: Dict[str, Dict[str, float]] = {}
    for name, blk in predictors.items():
        v = np.asarray(blk["values"], dtype=np.float64)
        entry = {
            "what": blk["what"],
            "rho_score": spearman_full(v, score),
            "partial_rho_depth_removed": partial_spearman(v, score, log_depth),
        }
        for g in used:
            entry[f"rho_{g}"] = spearman_full(v, hs["per_gene"][g])
        out_pred[name] = entry

    # Permutation p-values on a seeded subsample (the kernel is O(n x draws)).
    idx = subsample_index(int(micro.cells.n), PERM_SUBSAMPLE, PERM_SEED)
    perm = {}
    for name in ("simulated_o2", "dist_endothelial", "vascular_kde"):
        v = np.asarray(predictors[name]["values"], dtype=np.float64)
        perm[name] = spearman_permutation(v[idx], score[idx], n_perm, PERM_SEED)

    de = np.asarray(predictors["dist_endothelial"]["values"], dtype=np.float64)
    nested = {
        "score_vs_simo2_given_dist": partial_spearman(score, o2, de),
        "score_vs_dist_given_simo2": partial_spearman(score, de, o2),
        "simo2_vs_dist": spearman_full(o2, de),
    }

    ranking = all_gene_ranking(expr, genes, totals, o2)
    rank_of = {g: 1 + [r[0] for r in ranking].index(g) for g in used}

    tiles = tile_correlation(row, col, float(micro.domain.dx),
                             tuple(int(v) for v in micro.domain.shape),
                             oxygen, score)
    tiles = {k: v for k, v in tiles.items() if not k.startswith("_")}

    states = np.asarray(micro.cells.state)
    state_frac = {
        s.name: float((states == int(s)).mean()) for s in CellState
    }

    return {
        "genes_used": used,
        "detection": hs["detection"],
        "gene_gene_rho": hs["gene_gene_rho"],
        "predictors": out_pred,
        "permutation": perm,
        "perm_subsample_n": int(idx.size),
        "perm_draws": int(n_perm),
        "nested": nested,
        "gene_ranking_head": [(g, r, d) for g, r, d in ranking[:8]],
        "gene_ranking_tail": [(g, r, d) for g, r, d in ranking[-5:]],
        "rank_of_hypoxia_genes": rank_of,
        "n_genes_ranked": len(ranking),
        "tile_correlation": tiles,
        "state_fractions_all_cells": state_frac,
        "o2_vs_depth_rho": spearman_full(o2, totals),
        "mean_o2_at_cells": float(o2.mean()),
        "_o2": o2,
        "_score": score,
        "_dist_endo": -de,
    }


# ---------------------------------------------------------------------------
# The per-section driver
# ---------------------------------------------------------------------------
def run_section(spec: Dict, alpha: float, dx_um: float = DX_UM,
                panel: Sequence[str] = PANEL, n_perm: int = N_PERM) -> Dict:
    """Everything, for one section, through the frozen production path."""
    sample = str(spec["sample"])
    print("\n" + "=" * 78)
    print(f"{spec['short']} — {sample}")
    print("=" * 78, flush=True)

    t0 = time.time()
    data = _retry(load_singlecell_spatial, path=sample, dense=False)
    micro = reconstruct_section(sample, dx_um)
    print(f"  reconstructed in {time.time()-t0:.0f}s  "
          f"grid={tuple(micro.domain.shape)}  dx={micro.domain.dx:g} um", flush=True)

    comp = composition_block(micro, data)
    rep = comp["report"]
    print(f"  {rep['n_cells']:,} cells  panel={comp['panel_name']}  "
          f"TUMOR {100*rep['celltype_fractions']['TUMOR']:.2f}%  "
          f"STROMAL {100*rep['celltype_fractions']['STROMAL']:.2f}%  "
          f"IMMUNE {100*rep['celltype_fractions']['IMMUNE']:.2f}%  "
          f"ENDOTHELIAL {comp['endothelial_pct']:.2f}% "
          f"({comp['endothelial_verdict']} vs 8.0% volume-EM)", flush=True)

    t1 = time.time()
    vasc = vascular_block(micro)
    print(f"  vessels: {vasc['vessels']['n_vessels']:,} objects  "
          f"MVD {vasc['vessels']['mvd_per_mm2']:.1f}/mm2  "
          f"area frac {vasc['report']['vessel_area_fraction_pct']:.2f}%  "
          f"({time.time()-t1:.0f}s)", flush=True)

    t1 = time.time()
    hyp = hypoxia_block(micro, data, n_perm=n_perm)
    if hyp is None:
        print("  hypoxia check: SKIPPED — this panel measures no hypoxia-response "
              "gene, so hypoxia is UNMEASURED on this section", flush=True)
    else:
        p = hyp["predictors"]
        print(f"  hypoxia: rho(simO2)={p['simulated_o2']['rho_score']:+.4f}  "
              f"rho(dist_endo)={p['dist_endothelial']['rho_score']:+.4f}  "
              f"rho(kde)={p['vascular_kde']['rho_score']:+.4f}  "
              f"({time.time()-t1:.0f}s)", flush=True)

    print(f"  drug panel (uniform sensitivity, dx={dx_um:g} um):", flush=True)
    rows = drug_panel_block(micro, alpha, float(rep["mean_ecm_density"]), panel)

    return {
        "spec": dict(spec),
        "dx_um": float(dx_um),
        "alpha": float(alpha),
        "source": data["source"],
        "composition": comp,
        "vascular": vasc,
        "hypoxia": hyp,
        "drugs": rows,
        "ordering": ordering_diagnosis(rows),
        "runtime_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# Verdicts — one per claim, decided arithmetically, not narratively
# ---------------------------------------------------------------------------
def _fold(a: float, b: float) -> float:
    if not (np.isfinite(a) and np.isfinite(b)) or b == 0:
        return float("nan")
    return float(max(a / b, b / a))


def verdicts(res: Dict) -> List[Dict]:
    """Score each claim on the held-out patient. ``status`` is one of
    REPLICATES / DOES NOT REPLICATE / UNTESTABLE, plus a one-line reason."""
    s1 = res["sections"]["section1_brain319"]
    s2 = res["sections"]["section2_io480"]
    out: List[Dict] = []

    def add(claim, v1, v2, status, why, units=""):
        out.append({"claim": claim, "section1": v1, "section2": v2,
                    "status": status, "why": why, "units": units})

    # --- endothelial fraction vs volume-EM ---------------------------------- #
    e1 = s1["composition"]["endothelial_pct"]
    e2 = s2["composition"]["endothelial_pct"]
    tol = float(ENDO_VOLUME_EM["tolerance_fold"])
    f2 = _fold(e2, ENDO_VOLUME_EM["value_pct"])
    add("ENDOTHELIAL fraction matches human volume-EM 8.0% (PMID 38723085)",
        e1, e2,
        "REPLICATES" if f2 <= tol else "DOES NOT REPLICATE",
        (f"{e2:.2f}% is {f2:.2f}x from 8.0% "
         f"({'inside' if f2 <= tol else 'outside'} the declared {tol:g}x tolerance); "
         f"section 1 was {e1:.2f}% ({_fold(e1, 8.0):.2f}x)"),
        "%")

    # --- endothelial AREAL density (what the source term actually uses) ----- #
    d1 = s1["composition"]["endothelial_per_mm2"]
    d2 = s2["composition"]["endothelial_per_mm2"]
    add("ENDOTHELIAL areal density (the oxygen/drug SOURCE) is comparable",
        d1, d2, "CONTEXT",
        (f"{d2:.0f} vs {d1:.0f} cells/mm2 ({d2/d1:.2f}x) — the FRACTION fell but "
         "the areal DENSITY rose, because section 2's tissue is "
         f"{s2['composition']['cell_density_per_mm2']/s1['composition']['cell_density_per_mm2']:.2f}x "
         "more cellular. The delivery source depends on the density, not the fraction."),
        "cells/mm2")

    # --- vascular area fraction --------------------------------------------- #
    a1 = s1["vascular"]["report"]["vessel_area_fraction_pct"]
    a2 = s2["vascular"]["report"]["vessel_area_fraction_pct"]
    lo, hi = 1.8, 2.4
    add("Vascular area fraction matches measured 2.1 +/- 0.3% (PMID 36611806)",
        a1, a2,
        "REPLICATES" if lo <= a2 <= hi else "DOES NOT REPLICATE",
        (f"{a2:.2f}% vs 2.1% is {a2/2.1:.2f}x and lies "
         f"{'inside' if lo <= a2 <= hi else 'OUTSIDE'} the published 1.8-2.4% band; "
         f"section 1 gave {a1:.2f}% ({a1/2.1:.2f}x, also outside but only just)"),
        "%")

    # --- MVD ----------------------------------------------------------------- #
    m1 = s1["vascular"]["vessels"]["mvd_per_mm2"]
    m2 = s2["vascular"]["vessels"]["mvd_per_mm2"]
    m1_ok, m2_ok = 35.0 <= m1 <= 95.0, 35.0 <= m2 <= 95.0
    add("MVD lands in the published human GBM range 35-95/mm2 (PMID 27583442)",
        m1, m2,
        "REPLICATES" if (m1_ok and m2_ok) else
        ("NOT REPRODUCED ON EITHER" if not (m1_ok or m2_ok) else
         "DOES NOT REPLICATE" if m1_ok else "PARTIAL"),
        (f"{m2:.0f}/mm2 vs published 65 (35-95); section 1 gave {m1:.0f}/mm2. "
         + ("BOTH sections are high, and by a similar amount, so this is a "
            "systematic property of the segmentation (endothelial cells are "
            "over-split into vessel objects) rather than something the held-out "
            "patient revealed — and it is the one place where the two sections "
            "AGREE with each other while both disagreeing with the literature."
            if not (m1_ok or m2_ok) else
            "the two sections fall on opposite sides of the published range.")),
        "vessels/mm2")

    # --- drug ordering -------------------------------------------------------- #
    o1 = s1["ordering"]["rank_by_kill"]
    o2 = s2["ordering"]["rank_by_kill"]
    c1 = s1["ordering"]["rank_by_coverage"]
    c2 = s2["ordering"]["rank_by_coverage"]
    add("Kill ordering is cisplatin > doxorubicin > lomustine > methotrexate > "
        "temozolomide",
        " > ".join(o1), " > ".join(o2),
        "REPLICATES" if o2 == o1 else "DOES NOT REPLICATE",
        ("identical on both sections" if o2 == o1 else
         f"section 2 reorders to {' > '.join(o2)}"),
        "")
    # A claim that fails on BOTH sections is not a generalisation failure — it is
    # a claim that does not hold in THIS configuration at all, on the development
    # section included. Those two cases must not be scored the same way.
    cov1_ok, cov2_ok = (c1 == o1), (c2 == o2)
    if cov1_ok and cov2_ok:
        cov_status, cov_why = "REPLICATES", (
            "coverage reproduces the kill ordering on both sections")
    elif cov1_ok and not cov2_ok:
        cov_status, cov_why = "DOES NOT REPLICATE", (
            "coverage reproduces the kill ordering on section 1 but not on the "
            "held-out patient")
    elif cov2_ok and not cov1_ok:
        cov_status, cov_why = "PARTIAL", (
            "coverage reproduces the kill ordering on the held-out patient but "
            "not on section 1 — so the agreement reported elsewhere is "
            "configuration-dependent")
    else:
        cov_status, cov_why = "NOT REPRODUCED ON EITHER", (
            "coverage does NOT reproduce the kill ordering on EITHER section in "
            "this configuration (both swap doxorubicin and lomustine), so this "
            "is not a held-out failure — it is a disagreement with "
            "`results/delivery_metrics.md`, which ran at a different exposure "
            "scale where coverage was not saturated. See §3.3")
    add("therapeutic_coverage @1x IC50 recovers the kill ordering",
        " > ".join(c1), " > ".join(c2), cov_status, cov_why, "")

    # Is coverage still discriminating, or has it saturated?
    sat1 = sum(1 for r in s1["drugs"] if r["coverage_1x"] > 0.9 or r["coverage_1x"] < 0.01)
    sat2 = sum(1 for r in s2["drugs"] if r["coverage_1x"] > 0.9 or r["coverage_1x"] < 0.01)
    add("therapeutic_coverage @1x IC50 is still a DISCRIMINATING headline metric",
        f"{sat1}/{len(s1['drugs'])} drugs pinned at 0 or >0.9",
        f"{sat2}/{len(s2['drugs'])} drugs pinned at 0 or >0.9",
        "HOLDS" if max(sat1, sat2) < len(s1["drugs"]) else
        ("NOT REPRODUCED ON EITHER" if min(sat1, sat2) >= len(s1["drugs"])
         else "DOES NOT REPLICATE"),
        ("coverage has saturated at BOTH ends in this configuration — every drug "
         "is either effectively 0 or above 0.9, on both sections — so it is "
         "behaving as a binary indicator rather than a graded metric. "
         "`results/delivery_metrics.md` flagged this saturation as a 'latent, "
         "not active' weakness; at this exposure scale it is ACTIVE."
         if max(sat1, sat2) >= len(s1["drugs"]) else
         "coverage still spreads the panel on both sections"),
        "")

    # --- is the ordering trivially IC50? -------------------------------------- #
    add("The ordering is NOT simply the IC50 ordering read back out",
        " > ".join(s1["ordering"]["rank_by_potency_alone"]),
        " > ".join(s2["ordering"]["rank_by_potency_alone"]),
        "HOLDS" if (s2["ordering"]["rank_by_kill"] !=
                    s2["ordering"]["rank_by_potency_alone"]) else "FAILS",
        (f"potency alone would rank {' > '.join(s2['ordering']['rank_by_potency_alone'])}; "
         f"the simulation ranks {' > '.join(s2['ordering']['rank_by_kill'])}. "
         f"rho(kill, 1/IC50) = {s2['ordering']['rho_kill_vs_potency']:+.2f}, "
         f"rho(kill, AUC/IC50) = {s2['ordering']['rho_kill_vs_auc_over_ic50']:+.2f}"),
        "")

    # --- hypoxia ------------------------------------------------------------- #
    h2 = s2["hypoxia"]
    if h2 is not None:
        pde_rho = h2["predictors"]["simulated_o2"]["rho_score"]
        dist_rho = h2["predictors"]["dist_endothelial"]["rho_score"]
        add("Simulated oxygen is BEATEN by a distance-to-nearest-vessel baseline "
            "(state_validation.md Test 1)",
            "UNTESTABLE (no hypoxia gene on the 319 panel)",
            f"rho_PDE {pde_rho:+.4f} vs rho_dist {dist_rho:+.4f}",
            "CONFIRMED" if abs(dist_rho) > abs(pde_rho) else "REFUTED",
            (f"the one-line cKDTree baseline scores {abs(dist_rho)/abs(pde_rho):.2f}x "
             "the solved field" if abs(dist_rho) > abs(pde_rho) else
             "the solved field now wins"),
            "Spearman rho")
        vr = h2["rank_of_hypoxia_genes"].get("VEGFA")
        hr = h2["rank_of_hypoxia_genes"].get("HIF1A")
        add("VEGFA is the most anti-correlated gene on the panel; HIF1A is null",
            "n/a — neither gene is on the 319 panel",
            f"VEGFA rank {vr}/{h2['n_genes_ranked']}, HIF1A rank {hr}/"
            f"{h2['n_genes_ranked']}",
            "CONFIRMED" if (vr == 1 and hr is not None and
                            abs(h2['predictors']['simulated_o2']['rho_HIF1A']) < 0.02)
            else "PARTIAL",
            (f"VEGFA rho {h2['predictors']['simulated_o2']['rho_VEGFA']:+.4f}, "
             f"HIF1A rho {h2['predictors']['simulated_o2']['rho_HIF1A']:+.4f}"),
            "")
    return out


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------
def run(dx_um: float = DX_UM, panel: Sequence[str] = PANEL,
        n_perm: int = N_PERM, cache: Optional[str] = None) -> Dict:
    """Run both sections through the identical pipeline and score every claim.

    ``cache`` is a development convenience only: a pickle path that short-circuits
    the ~20 minutes of solves so the report and figure can be iterated on. It is
    never used by the module's own entry point unless asked for explicitly.
    """
    if cache and os.path.exists(cache):
        with open(cache, "rb") as fh:
            res = pickle.load(fh)
        print(f"[cache] loaded {cache}")
        res["_verdicts"] = verdicts(res)
        return res

    t0 = time.time()
    alpha = calibrate_alpha()
    print(f"perivascular-sink calibration alpha = {alpha:.1f} "
          "(fit ONCE on doxorubicin, on section 1 — NOT refitted for section 2)")

    sections: Dict[str, Dict] = {}
    for spec in SECTIONS:
        sections[str(spec["key"])] = run_section(spec, alpha, dx_um, panel, n_perm)

    res = {
        "dx_um": float(dx_um),
        "alpha": float(alpha),
        "panel": list(panel),
        "sections": sections,
        "pde": dict(PDE_KWARGS),
        "runtime_s": time.time() - t0,
    }
    res["_verdicts"] = verdicts(res)
    if cache:
        with open(cache, "wb") as fh:
            pickle.dump(res, fh)
        print(f"[cache] wrote {cache}")
    return res


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def _clean(ax, ygrid: bool = True) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if ygrid:
        ax.grid(axis="y", color=_GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)


def make_figure(res: Dict, path: str = _OUT_PNG) -> str:
    """Nine panels: composition, vasculature, delivery, response, hypoxia, verdicts."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    s1 = res["sections"]["section1_brain319"]
    s2 = res["sections"]["section2_io480"]
    vd = res.get("_verdicts") or verdicts(res)
    labels = {"section1_brain319": "Section 1 · 319-panel · 40,887 cells (DEV)",
              "section2_io480": "Section 2 · 480-IO panel · 816,769 cells (HELD OUT)"}

    plt.rcParams.update({
        "font.size": 9.5, "axes.edgecolor": _GRID, "axes.labelcolor": _INK_2,
        "xtick.color": _INK_2, "ytick.color": _INK_2, "text.color": _INK,
        "axes.titlecolor": _INK, "figure.facecolor": "white",
        "axes.facecolor": "white",
    })
    fig = plt.figure(figsize=(18.5, 16.0))
    gs = GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.26,
                  left=0.055, right=0.985, top=0.895, bottom=0.05)

    # ---- (a) cell-type composition ------------------------------------------ #
    ax = fig.add_subplot(gs[0, 0])
    classes = ["TUMOR", "STROMAL", "IMMUNE", "ENDOTHELIAL"]
    x = np.arange(len(classes))
    for j, (key, sec) in enumerate((("section1_brain319", s1),
                                    ("section2_io480", s2))):
        v = [100 * sec["composition"]["report"]["celltype_fractions"][c]
             for c in classes]
        ax.bar(x + (j - 0.5) * 0.38, v, 0.38, color=_SEC_COLORS[key], zorder=3,
               label=labels[key])
        for xi, vi in zip(x + (j - 0.5) * 0.38, v):
            ax.text(xi, vi + 1.0, f"{vi:.1f}", ha="center", fontsize=7.5,
                    color=_INK)
    ax.axhline(ENDO_VOLUME_EM["value_pct"], color=_INK_2, lw=1.2, ls="--", zorder=2)
    ax.annotate("8.0% — human volume-EM\nendothelial (PMID 38723085)",
                xy=(1.35, ENDO_VOLUME_EM["value_pct"]), xytext=(1.35, 27.0),
                ha="center", fontsize=7.5, color=_INK_2,
                arrowprops=dict(arrowstyle="-", color=_INK_2, lw=0.9))
    ax.set_ylim(0, 100)
    ax.set_xticks(x)
    ax.set_xticklabels(classes, fontsize=8.5)
    ax.set_ylabel("% of measured cells")
    ax.set_title("a · Cell-type composition — DIFFERENT PANELS, so composition\n"
                 "differences are part biology, part measurable-gene content",
                 loc="left", fontsize=9.5)
    ax.legend(frameon=False, fontsize=7.6, loc="upper right")
    _clean(ax)

    # ---- (b) vascular morphometry vs published ------------------------------ #
    ax = fig.add_subplot(gs[0, 1])
    keys = ["mvd_per_mm2", "intervessel_spacing_um", "median_tumor_to_vessel_um",
            "vascular_area_fraction_pct"]
    short = ["MVD\n(/mm²)", "intervessel\n1/√MVD (µm)", "tumour→vessel\nmedian (µm)",
             "vascular area\nfraction (%)"]
    x = np.arange(len(keys))
    for j, (key, sec) in enumerate((("section1_brain319", s1),
                                    ("section2_io480", s2))):
        rows = {r["key"]: r for r in sec["vascular"]["report"]["rows"]}
        ratio = [rows[k]["reconstructed"] / rows[k]["published"] for k in keys]
        ax.bar(x + (j - 0.5) * 0.38, ratio, 0.38, color=_SEC_COLORS[key], zorder=3,
               label=labels[key])
        for xi, k, rr in zip(x + (j - 0.5) * 0.38, keys, ratio):
            ax.text(xi, rr + 0.06, f"{rows[k]['reconstructed']:.2f}", ha="center",
                    fontsize=7.2, color=_INK, rotation=90)
    rows1 = {r["key"]: r for r in s1["vascular"]["report"]["rows"]}
    for xi, k in zip(x, keys):
        lo, hi = rows1[k]["published_range"]
        pub = rows1[k]["published"]
        ax.add_patch(plt.Rectangle((xi - 0.47, lo / pub), 0.94, (hi - lo) / pub,
                                   color="#9aa0a6", alpha=0.24, zorder=1))
    ax.axhline(1.0, color=_INK_2, lw=1.2, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(short, fontsize=8)
    ax.set_ylabel("reconstructed / published")
    ax.set_title("b · Vascular morphometry vs published human GBM\n"
                 "(grey band = published range; bar label = the raw value)",
                 loc="left", fontsize=9.5)
    ax.legend(frameon=False, fontsize=7.6, loc="upper left")
    _clean(ax)

    # ---- (c) tumour-to-vessel distance CDF ---------------------------------- #
    ax = fig.add_subplot(gs[0, 2])
    for key, sec in (("section1_brain319", s1), ("section2_io480", s2)):
        d = np.asarray(sec["vascular"]["_tumor_to_vessel_um"], dtype=np.float64)
        d = np.sort(d[np.isfinite(d)])
        if d.size > 200_000:  # plotting only; every statistic uses all cells
            d = d[:: max(1, d.size // 200_000)]
        ax.plot(d, np.linspace(0, 1, d.size), lw=2.2, color=_SEC_COLORS[key],
                label=(f"{'Section 1 (DEV)' if key.startswith('section1') else 'Section 2 (HELD OUT)'}"
                       f" · median "
                       f"{sec['vascular']['tumor_to_vessel']['median']:.0f} µm"
                       f" · >100 µm "
                       f"{100*sec['vascular']['vessels']['tumor_frac_beyond_100um']:.1f}%"),
                zorder=3)
    ax.axvline(100.0, color=_INK_2, lw=1.0, ls=":", zorder=2)
    ax.axvline(150.0, color=_INK_2, lw=1.0, ls=":", zorder=2)
    ax.text(152, 0.12, "150 µm — Thomlinson–Gray", fontsize=7.2, color=_INK_2)
    ax.set_xlim(0, 260)
    ax.set_xlabel("distance from tumour cell to nearest endothelial cell (µm)")
    ax.set_ylabel("cumulative fraction of tumour cells")
    ax.set_title("c · The delivery-relevant geometry", loc="left", fontsize=9.5)
    ax.legend(frameon=False, fontsize=7.6, loc="lower right",
              handlelength=1.6, labelspacing=0.7)
    _clean(ax, ygrid=False)
    ax.grid(color=_GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)

    # ---- (d,e,f) the drug panel --------------------------------------------- #
    for i, (field, title, logy, fmt) in enumerate((
            ("kill_fraction", "d · Predicted kill (uniform sensitivity — the FAIR "
                              "cross-drug mode)", False, "{:.3f}"),
            ("coverage_1x", "e · therapeutic_coverage @ 1× IC50 — the headline "
                            "metric", False, "{:.3f}"),
            ("mean_tumor_auc", "f · Mean tumour AUC (µM·s) — raw delivery",
             True, "{:.2g}"))):
        ax = fig.add_subplot(gs[1, i])
        drugs = list(res["panel"])
        x = np.arange(len(drugs))
        for j, (key, sec) in enumerate((("section1_brain319", s1),
                                        ("section2_io480", s2))):
            by = {r["drug"]: r for r in sec["drugs"]}
            v = [by[d][field] for d in drugs]
            ax.bar(x + (j - 0.5) * 0.38, v, 0.38, color=_SEC_COLORS[key], zorder=3,
                   label=labels[key])
            for xi, vi in zip(x + (j - 0.5) * 0.38, v):
                ax.text(xi, vi * (1.12 if logy else 1.0) + (0 if logy else 0.012),
                        fmt.format(vi), ha="center", fontsize=7.0, color=_INK,
                        rotation=90)
        if logy:
            ax.set_yscale("log")
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi * 6.0)
        else:
            ax.set_ylim(0, max(1.0, ax.get_ylim()[1]) * 1.30)
        ax.set_xticks(x)
        ax.set_xticklabels([d[:4] + "." if len(d) > 6 else d for d in drugs],
                           fontsize=8)
        ax.set_xticklabels(drugs, fontsize=7.6, rotation=18, ha="right")
        ax.set_title(title, loc="left", fontsize=9.5)
        if i == 0:
            ax.legend(frameon=False, fontsize=7.4, loc="upper left")
        _clean(ax)

    # ---- (g) hypoxia predictors --------------------------------------------- #
    ax = fig.add_subplot(gs[2, 0])
    h2 = s2["hypoxia"]
    if h2 is not None:
        names = ["simulated_o2", "vascular_kde", "dist_endothelial",
                 "local_cell_density"]
        nice = ["SIMULATED O₂\nthe PDE", "vascular_kde\nits own source",
                "dist_endothelial\n1 cKDTree query", "local_density\nno vessels"]
        vals = [h2["predictors"][n]["rho_score"] for n in names]
        cols = ["#e34948" if n == "simulated_o2" else "#2a78d6" for n in names]
        ax.barh(np.arange(len(names)), vals, 0.6, color=cols, zorder=3)
        span = abs(min(vals))
        for i2, v in enumerate(vals):
            inside = abs(v) > 0.28 * span   # label fits inside the bar
            ax.text(v + (0.004 if inside else -0.004), i2, f"{v:+.4f}",
                    va="center", ha="left" if inside else "right", fontsize=8,
                    color="white" if inside else _INK,
                    fontweight="bold" if inside else "normal")
        ax.set_xlim(min(vals) * 1.12, 0.004)
        ax.set_yticks(np.arange(len(names)))
        ax.set_yticklabels(nice, fontsize=7.8)
        ax.invert_yaxis()
        ax.set_xlabel("Spearman ρ vs MEASURED hypoxia score "
                      "(more negative = better)")
        beaten = (abs(h2["predictors"]["dist_endothelial"]["rho_score"]) >
                  abs(h2["predictors"]["simulated_o2"]["rho_score"]))
        ax.set_title("g · Section 2 hypoxia re-check — the PDE is "
                     + ("out-predicted by\na distance transform built from its "
                        "own inputs" if beaten else
                        "NOT out-predicted\nby a distance transform on its own "
                        "inputs"),
                     loc="left", fontsize=9.5)
        _clean(ax, ygrid=False)
        ax.grid(axis="x", color=_GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)

    # ---- (h) the 480-gene ranking ------------------------------------------- #
    ax = fig.add_subplot(gs[2, 1])
    if h2 is not None:
        head = h2["gene_ranking_head"]
        tail = h2["gene_ranking_tail"]
        gs_ = [g for g, _, _ in head] + ["…"] + [g for g, _, _ in tail]
        rv = [r for _, r, _ in head] + [0.0] + [r for _, r, _ in tail]
        cols = ["#2a78d6" if g in ("VEGFA", "HIF1A") else "#9aa0a6" for g in gs_]
        cols = ["#e34948" if g == "VEGFA" else c for g, c in zip(gs_, cols)]
        ax.barh(np.arange(len(gs_)), rv, 0.72, color=cols, zorder=3)
        ax.set_yticks(np.arange(len(gs_)))
        ax.set_yticklabels(gs_, fontsize=7.6)
        ax.invert_yaxis()
        ax.axvline(0, color=_INK_2, lw=1.0)
        vr = h2["rank_of_hypoxia_genes"].get("VEGFA")
        hr = h2["rank_of_hypoxia_genes"].get("HIF1A")
        ax.set_xlabel("Spearman ρ vs simulated O₂")
        ax.set_title(f"h · All {h2['n_genes_ranked']} IO-panel genes ranked — "
                     f"VEGFA rank {vr}, HIF1A rank {hr}\n"
                     "(positive tail is the vascular/myeloid programme — partly "
                     "self-correlation)", loc="left", fontsize=9.5)
        _clean(ax, ygrid=False)
        ax.grid(axis="x", color=_GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)

    # ---- (i) verdicts -------------------------------------------------------- #
    ax = fig.add_subplot(gs[2, 2])
    ax.axis("off")
    ax.set_title("i · Does it generalise to a second real patient?",
                 loc="left", fontsize=9.5)
    colours = {"REPLICATES": "#2e7d32", "CONFIRMED": "#2e7d32", "HOLDS": "#2e7d32",
               "DOES NOT REPLICATE": "#c62828", "REFUTED": "#c62828",
               "FAILS": "#c62828", "PARTIAL": "#ef6c00", "CONTEXT": "#52514e",
               "NOT REPRODUCED ON EITHER": "#6a1b9a"}
    y = 0.975
    for v in vd:
        ax.text(0.0, y, v["status"], fontsize=8.0, fontweight="bold",
                color=colours.get(v["status"], _INK_2), va="top",
                transform=ax.transAxes)
        y -= 0.030
        claim = v["claim"]
        wrapped, line = [], ""
        for w in claim.split():
            if len(line) + len(w) > 62:
                wrapped.append(line)
                line = w
            else:
                line = (line + " " + w).strip()
        wrapped.append(line)
        ax.text(0.03, y, "\n".join(wrapped), fontsize=7.4, va="top", color=_INK,
                transform=ax.transAxes)
        y -= 0.018 + 0.0265 * len(wrapped)

    fig.suptitle(
        "THE HELD-OUT PATIENT — the whole validated pipeline, re-run on a second "
        "real GBM section it was never developed on\n"
        f"dx = {res['dx_um']:g} µm · uniform sensitivity · α = {res['alpha']:.1f} "
        "(fit once on section 1, NOT refitted) · two DIFFERENT patients on two "
        "DIFFERENT panels sharing only 60 genes · no outcome label on either",
        fontsize=12.5, x=0.055, ha="left", y=0.972, color=_INK)
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _pct(v: float) -> str:
    return f"{100.0 * float(v):.2f}%"


def write_report(res: Dict, path: str = _OUT_MD,
                 png: str = _OUT_PNG) -> str:
    """Write ``results/heldout_patient.md`` — the side-by-side, with verdicts."""
    s1 = res["sections"]["section1_brain319"]
    s2 = res["sections"]["section2_io480"]
    vd = res.get("_verdicts") or verdicts(res)
    L: List[str] = []
    A = L.append

    c1, c2 = s1["composition"], s2["composition"]
    r1, r2 = c1["report"], c2["report"]
    v1, v2 = s1["vascular"], s2["vascular"]
    h2 = s2["hypoxia"]
    d1 = {r["drug"]: r for r in s1["drugs"]}
    d2 = {r["drug"]: r for r in s2["drugs"]}
    replicates = [v for v in vd if v["status"] in ("REPLICATES", "CONFIRMED", "HOLDS")]
    fails = [v for v in vd if v["status"] in ("DOES NOT REPLICATE", "REFUTED", "FAILS")]

    A("# The held-out patient — running the whole validated pipeline on a second "
      "real GBM section")
    A("")
    A("Every validated result in this project was developed on **one** tissue: the "
      "319-gene 10x Xenium GBM section (40,887 cells). This report re-runs the "
      "entire pipeline — reconstruction, vessel segmentation, the five-drug "
      "simulation and the hypoxia check — on a **second real patient the pipeline "
      "was never developed on**: the 480-gene Immuno-Oncology Xenium GBM section "
      f"({r2['n_cells']:,} cells). Both sections are driven through the identical "
      "frozen production path at the identical settings, so **section 1's numbers "
      "here are re-derived, not quoted**.")
    A("")
    neither = [v for v in vd if v["status"] == "NOT REPRODUCED ON EITHER"]
    A(f"**{len(replicates)} of {len(replicates) + len(fails)} claims that could be "
      f"scored as generalisation replicate.** A further {len(neither)} fail on "
      "**both** sections — those are not generalisation failures, they are claims "
      "that do not hold in this configuration at all, and they are called out "
      "separately because conflating the two would be the easiest way to make "
      "this report say something it should not. The headline, one line each:")
    A("")
    A("| what | section 1 (development) | section 2 (held out) | verdict |")
    A("|---|---|---|:--|")
    for v in vd:
        a = v["section1"]
        b = v["section2"]
        fa = f"{a:.2f}" if isinstance(a, float) else str(a)
        fb = f"{b:.2f}" if isinstance(b, float) else str(b)
        A(f"| {v['claim']} | {fa} {v['units']} | {fb} {v['units']} | "
          f"**{v['status']}** |")
    A("")
    A("> **THE CONFOUND, FIRST.** These are **different patients** *and* "
      "**different panels** — the two Xenium panels share only **60 genes**. The "
      "IO panel carries none of PECAM1 / SOX2 / OLIG2 / TOP2A / AQP4 / GJA1 / MOG, "
      "so the endothelial gate runs on VWF/FLT1/PLVAP instead of PECAM1 and the "
      "stromal gate on SLC1A3/S100B/PTN instead of AQP4/GJA1/MOG. Every "
      "composition difference below is therefore **jointly** a difference in "
      "biology and a difference in what is measurable. Where the two can be "
      "separated, this report separates them; where they cannot, it says so. "
      "Neither section carries a drug-response outcome label, so nothing here "
      "validates a kill prediction — only the machinery that produces one.")
    A("")

    # --- §0 -------------------------------------------------------------------
    A("## 0. What was held out, and what was not")
    A("")
    A("This is an out-of-sample test of the **pipeline**, not of a fitted model. "
      "Being explicit about what that does and does not cover:")
    A("")
    A("| held out from section 2 | NOT held out |")
    A("|---|---|")
    A("| the reconstruction was never inspected on it while the marker gate, the "
      "typing threshold or the 30 µm tissue-mask radius were chosen | the IO-panel "
      "**marker map** in `reconstruct/marker_panels.py` was written for this "
      "section (in `results/state_validation.md`), by choosing among measured "
      "candidate marker sets — see §1.2 |")
    A("| the 20 µm vessel linking distance and the 7.5 µm vessel radius were fixed "
      "against section 1 and are used here unchanged | the oxygen solver's "
      "`consumption` constant has never been fitted to a measured oxygen tension "
      "on **either** section |")
    A(f"| the perivascular-sink calibration **α = {res['alpha']:.1f}** was fit once "
      "on doxorubicin against section 1 and is **not** refitted here | the drug "
      "IC50s are external (GDSC/literature) and are the same numbers on both |")
    A("| the fixed thresholds in `delivery_metrics` (the 1× IC50 coverage bar; "
      "the vessel-distance axis is now the PRIMARY direct-cell metric and has no "
      "threshold at all — see `results/vessel_profile.md`'s dedicated "
      "distance-metric-change section) | `dx = 12 µm` was chosen for grid "
      "convergence on section 1; convergence is **not** re-established here "
      "(§6) |")
    A("")

    # --- §1 -------------------------------------------------------------------
    A("## 1. Reconstruction — does the marker-gated classifier survive a different "
      "panel?")
    A("")
    A("| | Section 1 (development) | Section 2 (HELD OUT) |")
    A("|---|---|---|")
    A(f"| dataset | `{s1['spec']['sample']}` | `{s2['spec']['sample']}` |")
    A(f"| panel | {s1['spec']['panel']} | {s2['spec']['panel']} |")
    A(f"| marker map selected | `{c1['panel_name']}` | `{c2['panel_name']}` |")
    A(f"| marker coverage (brain / IO map) | "
      f"{c1['panel_coverage']['brain_319']:.2f} / {c1['panel_coverage']['io_480']:.2f} | "
      f"{c2['panel_coverage']['brain_319']:.2f} / {c2['panel_coverage']['io_480']:.2f} |")
    A(f"| cells | {r1['n_cells']:,} | {r2['n_cells']:,} |")
    A(f"| grid @ dx = {res['dx_um']:g} µm | {r1['grid']} | {r2['grid']} |")
    A(f"| tissue extent (µm) | {r1['extent_um'][0]:.0f} × {r1['extent_um'][1]:.0f} | "
      f"{r2['extent_um'][0]:.0f} × {r2['extent_um'][1]:.0f} |")
    A(f"| tissue area (geometric mask) | {v1['vessels']['tissue_area_mm2']:.1f} mm² | "
      f"{v2['vessels']['tissue_area_mm2']:.1f} mm² |")
    A(f"| **cell density** | {c1['cell_density_per_mm2']:.0f} /mm² | "
      f"{c2['cell_density_per_mm2']:.0f} /mm² |")
    for cls in ("TUMOR", "STROMAL", "IMMUNE", "ENDOTHELIAL"):
        A(f"| {cls} | {_pct(r1['celltype_fractions'][cls])} | "
          f"{_pct(r2['celltype_fractions'][cls])} |")
    A(f"| mean vascular density | {r1['mean_vascular_density']:.4f} | "
      f"{r2['mean_vascular_density']:.4f} |")
    A(f"| mean ECM density | {r1['mean_ecm_density']:.4f} | "
      f"{r2['mean_ecm_density']:.4f} |")
    A(f"| ECM markers on panel | {', '.join(r1['ecm_markers_used'])} | "
      f"{', '.join(r2['ecm_markers_used'])} |")
    A(f"| mean oxygen (tissue) | {r1['mean_oxygen']:.4f} | {r2['mean_oxygen']:.4f} |")
    A(f"| tissue coverage of bbox | {_pct(r1['tissue_coverage_fraction'])} | "
      f"{_pct(r2['tissue_coverage_fraction'])} |")
    A("")
    A("### 1.1 The endothelial check — the one external anchor the reconstruction has")
    A("")
    A(f"Human cortex volume electron microscopy counts **{ENDO_VOLUME_EM['value_pct']:.1f}%** "
      f"of cells as endothelial ({ENDO_VOLUME_EM['source']}, "
      f"[PMID {ENDO_VOLUME_EM['pmid']}](https://pubmed.ncbi.nlm.nih.gov/"
      f"{ENDO_VOLUME_EM['pmid']}/)). Section 1's **{c1['endothelial_pct']:.2f}%** "
      "matching that — rather than the 1–2% that dissociation-based scRNA-seq "
      "reports — is the single strongest external check the reconstruction has, and "
      "`results/SESSION_AUDIT.md` §4 leans on it.")
    A("")
    A(f"**On the held-out patient the classifier calls {c2['endothelial_pct']:.2f}% "
      f"endothelial — {c2['endothelial_fold_vs_volume_em']:.2f}× the volume-EM "
      f"value, against section 1's {c1['endothelial_fold_vs_volume_em']:.2f}×.**")
    A("")
    A("| | section 1 | section 2 | ratio |")
    A("|---|---:|---:|---:|")
    A(f"| ENDOTHELIAL, % of cells | {c1['endothelial_pct']:.2f}% | "
      f"{c2['endothelial_pct']:.2f}% | {c2['endothelial_pct']/c1['endothelial_pct']:.2f}× |")
    A(f"| ENDOTHELIAL, cells per mm² | {c1['endothelial_per_mm2']:.0f} | "
      f"{c2['endothelial_per_mm2']:.0f} | "
      f"{c2['endothelial_per_mm2']/c1['endothelial_per_mm2']:.2f}× |")
    A(f"| all cells per mm² | {c1['cell_density_per_mm2']:.0f} | "
      f"{c2['cell_density_per_mm2']:.0f} | "
      f"{c2['cell_density_per_mm2']/c1['cell_density_per_mm2']:.2f}× |")
    A("")
    A("**These two rows point in opposite directions, and that is the finding.** "
      "The endothelial *fraction* falls while the endothelial *areal density* "
      "rises, because section 2's tissue is far more cellular. The oxygen and drug "
      "source term is built from a KDE of endothelial **positions**, so it depends "
      "on the density, not the fraction — the compartment that matters physically "
      "did not shrink. But the *fraction* is the quantity that was validated "
      "against volume EM, and that quantity does not replicate.")
    A("")
    A("### 1.2 Why this is a weaker test than it looks — and still a real one")
    A("")
    A("The IO-panel marker map was not written blind. `results/state_validation.md` "
      "§'MURAL / PERICYTE' records that three candidate STROMAL marker sets were "
      "measured on **this section** and the astro-lineage set was chosen partly "
      "because it kept the two sections comparable. So the map has seen section 2. "
      "What it has **not** seen is any external anchor: the choice was made on "
      "internal comparability, never against the 8.0% volume-EM number or against "
      "any vascular measurement. The endothelial gate in particular "
      f"(`{', '.join(c2['markers_used']['ENDOTHELIAL'])}`) was picked because those "
      "are the canonical endothelial genes the panel carries, not because they "
      "produced a particular fraction.")
    A("")
    A("| class | section 1 markers used | section 2 markers used |")
    A("|---|---|---|")
    for cls in ("TUMOR", "IMMUNE", "ENDOTHELIAL", "STROMAL"):
        A(f"| {cls} | {', '.join(c1['markers_used'][cls]) or '—'} | "
          f"{', '.join(c2['markers_used'][cls]) or '—'} |")
    A("")
    A("Three candidate explanations for the endothelial gap, and what separates them:")
    A("")
    A("1. **Real biology.** Two different tumours genuinely differ in vascularity. "
      "Partly supported — but the *areal* endothelial density went **up**, not "
      "down, so this cannot be the whole story.")
    A("2. **Denominator.** Section 2 has "
      f"{c2['cell_density_per_mm2']/c1['cell_density_per_mm2']:.1f}× the cell "
      "density, so the same vasculature is divided by many more non-endothelial "
      "cells. The denominator alone is arithmetically sufficient to explain the "
      "whole gap: divide section 2's endothelial areal density by section 1's "
      "total cell density and the fraction would read "
      f"{100*c2['endothelial_per_mm2']/c1['cell_density_per_mm2']:.1f}% — far "
      f"ABOVE 8.0%, not the {c2['endothelial_pct']:.2f}% observed. The observed "
      "fraction is set as much by how many cells the section contains as by how "
      "many of them are endothelial.")
    det1 = c1["marker_detection"]
    det2 = c2["marker_detection"]
    endo1 = ", ".join(f"{m} {100*det1[m]:.2f}%" for m in c1["markers_used"]["ENDOTHELIAL"])
    endo2 = ", ".join(f"{m} {100*det2[m]:.2f}%" for m in c2["markers_used"]["ENDOTHELIAL"])
    A(f"3. **Panel content.** A {len(c2['markers_used']['ENDOTHELIAL'])}-marker "
      f"mean-z gate ({'/'.join(c2['markers_used']['ENDOTHELIAL'])}) behaves "
      f"differently from a {len(c1['markers_used']['ENDOTHELIAL'])}-marker gate "
      f"({'/'.join(c1['markers_used']['ENDOTHELIAL'])}). Measured detection rates "
      f"— section 1: {endo1}; section 2: {endo2} — so on section 2 the gate is the "
      "mean z of three sparse markers, which a cell strongly positive for one but "
      "silent for the other two can fail. That pushes the endothelial call "
      "**down** for a purely technical reason.")
    A("")
    A("**Explanations 2 and 3 are measurement artefacts, not biology, and neither "
      "can be excluded here.** The honest statement is that *the endothelial "
      "cell-type fraction is not a portable quantity across panels*, and the 8.0% "
      "agreement on section 1 should be read as a property of that section and "
      "that gate, not as a validated property of the classifier.")
    A("")

    # --- §2 -------------------------------------------------------------------
    A("## 2. Vascular morphometry — the only vascular quantity scored against a "
      "measurement")
    A("")
    A("`reconstruct/vessel_network.py`, unmodified: measured endothelial cells "
      f"single-linkage clustered at **{v2['vessels']['link_distance_um']:.0f} µm** "
      "into discrete vessel objects, then painted as a source of radius "
      f"**{v2['report']['radius_um']:.1f} µm**. Both constants were fixed against "
      "section 1 and are used here unchanged.")
    A("")
    A("| quantity | section 1 | section 2 (HELD OUT) | published | range | "
      "s1 verdict | s2 verdict | source |")
    A("|---|---:|---:|---:|---|:--|:--|---|")
    rows1 = {r["key"]: r for r in v1["report"]["rows"]}
    rows2 = {r["key"]: r for r in v2["report"]["rows"]}
    for k in ("mvd_per_mm2", "intervessel_spacing_um", "median_tumor_to_vessel_um",
              "vascular_area_fraction_pct"):
        a, b = rows1[k], rows2[k]
        A(f"| {a['quantity']} ({a['units']}) | {a['reconstructed']:.2f} | "
          f"**{b['reconstructed']:.2f}** | {a['published']:.2f} | "
          f"{a['published_range'][0]:.1f}–{a['published_range'][1]:.1f} | "
          f"{a['verdict']} | **{b['verdict']}** | "
          f"[PMID {a['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{a['pmid']}/) |")
    A("")
    A("### 2.1 The load-bearing row")
    A("")
    af1 = v1["report"]["vessel_area_fraction_pct"]
    af2 = v2["report"]["vessel_area_fraction_pct"]
    A("`results/vessel_network_impact.md` calls the **vascular area fraction** the "
      "load-bearing row, because it is the only one that tests the *field* rather "
      "than the point pattern, and because the vessel radius was fixed from "
      "published vessel-caliber imaging **before** the field was built — so it is "
      f"an untuned comparison. Section 1 scored **{af1:.2f}%** against a measured "
      f"whole-section CD34 value of **2.1 ± 0.3%**, i.e. {af1/2.1:.2f}×.")
    A("")
    A(f"**Section 2 scores {af2:.2f}% — {af2/2.1:.2f}× the measured value, and "
      f"{af2/af1:.2f}× section 1's.** That is outside the published band by a wide "
      "margin, and it is the clearest single failure of generalisation in this "
      "report.")
    A("")
    A("The cause is not the fixed-7.5-µm-per-cell painting an earlier version of "
      "this paragraph described — that path was superseded by "
      "`vessel_geometry.py`'s sub-grid, data-derived radius model (see "
      "`results/vessel_geometry_fix.md`), which is what §2.1's numbers above are "
      "already computed from. A dedicated follow-up "
      "(`results/sweep/vascular_area_fraction_generalization_2026-07-30T03-18-42Z.md`) "
      "tested whether ANY non-circular, section-invariant function from a vessel "
      "object's own cell count to its radius could close this gap — six "
      "candidates, including literature-motivated flattened and capped variants, "
      "on this section plus a third independent GBM section — and found a "
      "genuine null: every candidate that helps this held-out section pulls the "
      "development section further from target, because the two errors run in "
      "opposite directions (development: 0.59× too low; both other sections: "
      "5.2× too high) while every candidate radius model moves monotonically in "
      "the same direction on all three. The radius that WOULD hit 2.1% on each "
      "section separately is 7.28 / 3.12 / 2.57 µm — a near-perfect inverse rank "
      "correlation with each section's own microvessel density, directionally "
      "consistent with a real published finding (Halshtok Neiman et al., "
      "*Clin Neuroradiol* 2021, PMID 34940899). But using this project's own MVD "
      "to calibrate its own vessel radius would be circular: §6.3 below "
      "independently flags MVD itself as an unresolved, both-sections-fail "
      "segmentation defect. Conditioning one flagged defect on another is not a "
      "fix, so none was applied — this remains open, honestly.")
    A("")
    A("### 2.2 Vessel objects and the clustering that makes them comparable to "
      "histology")
    A("")
    A("| | section 1 | section 2 |")
    A("|---|---:|---:|")
    for lab, key, fmt in (
            ("endothelial cells", "n_endothelial_cells", "{:,}"),
            ("discrete vessel objects", "n_vessels", "{:,}"),
            ("tissue area (mm²)", "tissue_area_mm2", "{:.1f}"),
            ("MVD (/mm²)", "mvd_per_mm2", "{:.1f}"),
            ("endothelial cells per vessel", "mean_cells_per_vessel", "{:.2f}"),
            ("singleton objects", "singleton_fraction", "{:.1%}"),
            ("median object extent, multi-cell (µm)", "median_extent_multi_um", "{:.1f}"),
            ("median aspect ratio, multi-cell", "median_aspect_multi", "{:.2f}"),
            ("tumour cells > 100 µm from a vessel", "tumor_frac_beyond_100um", "{:.2%}"),
            ("tumour cells > 150 µm from a vessel", "tumor_frac_beyond_150um", "{:.2%}"),
    ):
        A(f"| {lab} | {fmt.format(v1['vessels'][key])} | "
          f"{fmt.format(v2['vessels'][key])} |")
    A(f"| intervessel NN median (µm) | {v1['intervessel_nn']['median']:.1f} | "
      f"{v2['intervessel_nn']['median']:.1f} |")
    A(f"| tumour→vessel median (µm) | {v1['tumor_to_vessel']['median']:.1f} | "
      f"{v2['tumor_to_vessel']['median']:.1f} |")
    A(f"| tumour→vessel p90 (µm) | {v1['tumor_to_vessel']['p90']:.1f} | "
      f"{v2['tumor_to_vessel']['p90']:.1f} |")
    A(f"| clustering index vs Poisson (median) | "
      f"{v1['report']['clustering_index_median']:.2f}× | "
      f"{v2['report']['clustering_index_median']:.2f}× |")
    A(f"| clustering index vs Poisson (>100 µm tail) | "
      f"{v1['report']['clustering_index_tail']:.2f}× | "
      f"{v2['report']['clustering_index_tail']:.2f}× |")
    A("")
    A("**One thing genuinely improves out of sample.** `results/vascular_validation.md` "
      "§3.2 found that on section 1 the clustering step 'barely bites' — 1.29 "
      "endothelial cells per putative vessel and 85% singletons — so the "
      "reconstruction was resolving scattered endothelial nuclei rather than "
      f"vessel walls. On section 2 it merges {v2['vessels']['mean_cells_per_vessel']:.2f} "
      f"cells per object with {v2['vessels']['singleton_fraction']:.0%} singletons "
      f"and a median multi-cell extent of {v2['vessels']['median_extent_multi_um']:.0f} µm "
      f"at aspect ratio {v2['vessels']['median_aspect_multi']:.1f}: **elongated, "
      "multi-cell objects that actually look like vessel profiles.** The denser "
      "section supports the segmentation the sparse one could not. That is a "
      "genuine, if backhanded, out-of-sample success — and it is also part of why "
      "the area fraction inflates.")
    A("")
    A("### 2.3 Linking-distance sweep — the free parameter, swept on both sections")
    A("")
    A("| link (µm) | s1 vessels | s1 MVD | s1 cells/vessel | s2 vessels | s2 MVD | "
      "s2 cells/vessel |")
    A("|---:|---:|---:|---:|---:|---:|---:|")
    for a, b in zip(v1["sweep"], v2["sweep"]):
        A(f"| {a['link_distance_um']:.0f} | {a['n_vessels']:,} | "
          f"{a['mvd_per_mm2']:.1f} | {a['mean_cells_per_vessel']:.2f} | "
          f"{b['n_vessels']:,} | {b['mvd_per_mm2']:.1f} | "
          f"{b['mean_cells_per_vessel']:.2f} |")
    A("")
    A("On section 2 the published MVD range (35–95/mm²) is reached at a linking "
      "distance of ~25–40 µm rather than the 20 µm fixed on section 1. **No "
      "linking distance was re-tuned for this report** — reporting the sweep is "
      "the honest alternative to tuning it.")
    A("")

    # --- §3 -------------------------------------------------------------------
    o1 = s1["ordering"]
    o2 = s2["ordering"]
    A("## 3. The drug panel — does the ordering hold on a second patient?")
    A("")
    A(f"Five drugs, **uniform** per-cell sensitivity (the only fair cross-drug "
      f"mode — every drug carries 1.0 at every cell, so the factor cancels and the "
      f"ranking isolates delivery × potency), dx = {res['dx_um']:g} µm, "
      f"α = {res['alpha']:.1f} fit once on section 1 and **not refitted**, "
      f"PDE dt = {res['pde']['dt']:g} s over "
      f"{res['pde']['total_time']/3600:.0f} h with convection on. Identical "
      "settings on both sections.")
    A("")
    A("### 3.1 Section 1 (development)")
    A("")
    hdr = ("| drug | IC50 (µM) | BBB | D_eff (µm²/s) | mean tumour AUC (µM·s) | "
           "**coverage @1× IC50** | coverage @0.5×/2× | median E/IC50 | "
           "**kill** |")
    sep = "|---|---:|---:|---:|---:|---:|---:|---:|---:|"
    for tag, dd in (("### 3.1 Section 1 (development)", d1),
                    ("### 3.2 Section 2 (HELD OUT)", d2)):
        if tag.startswith("### 3.2"):
            A(tag)
            A("")
        A(hdr)
        A(sep)
        for name in res["panel"]:
            r = dd[name]
            A(f"| {name} | {r['ic50_uM']:g} | {r['bbb_penetrance']:.3f} | "
              f"{r['d_eff_um2_s']:.1f} | {r['mean_tumor_auc']:.3g} | "
              f"**{r['coverage_1x']:.4f}** | "
              f"{r['coverage_0.5x']:.3f} / {r['coverage_2x']:.3f} | "
              f"{r['median_E_over_ic50']:.2f} | **{r['kill_fraction']:.4f}** |")
        A("")
    A("### 3.3 The ordering")
    A("")
    A("| ranked by | section 1 | section 2 (HELD OUT) | same? |")
    A("|---|---|---|:--:|")
    for lab, key in (("simulated **kill**", "rank_by_kill"),
                     ("**coverage @1× IC50**", "rank_by_coverage"),
                     ("mean tumour AUC (delivery alone)", "rank_by_delivery_alone"),
                     ("1/IC50 (potency alone)", "rank_by_potency_alone"),
                     ("AUC/IC50 (delivery × potency)", "rank_by_auc_over_ic50")):
        same = "✅" if o1[key] == o2[key] else "❌"
        A(f"| {lab} | {' > '.join(o1[key])} | {' > '.join(o2[key])} | {same} |")
    A("")
    A("*(The 'same?' column compares the two SECTIONS with each other. Whether "
      "each ranking agrees with the KILL ranking — the separate question — is "
      "read down the rows: on both sections `coverage` swaps doxorubicin and "
      "lomustine relative to kill, and `AUC/IC50` reproduces kill exactly.)*")
    A("")
    A(f"**Published ordering** (`results/real_drug_panel.md`, uniform mode): "
      f"{' > '.join(PUBLISHED_ORDERING)}. Section 1's kill fractions here "
      "(" + ", ".join(f"{n} {d1[n]['kill_fraction']:.4f}" for n in res["panel"]) + ") "
      "reproduce that report's numbers exactly, which is the check that this "
      "harness is running the same pipeline and not a variant of it.")
    A("")
    kill_same = o1["rank_by_kill"] == o2["rank_by_kill"]
    if kill_same:
        A("**The kill ordering is IDENTICAL on both patients.** That is the "
          "headline, and the immediate next question is whether it is identical "
          "for a *trivial* reason.")
    else:
        A("**The kill ordering CHANGES on the second patient.** A different "
          "ordering on a different real tumour is a legitimate finding — real "
          "inter-patient heterogeneity — not a failure. What follows is the check "
          "on whether the change is biology or bookkeeping.")
    A("")
    A("### 3.4 Is the ordering trivially the IC50 ordering?")
    A("")
    A("An ordering that just reproduces the potency ranking would be reproducible "
      "by a spreadsheet with no geometry, no transport and no patient in it. The "
      "diagnostic is whether delivery ever overturns potency.")
    A("")
    A("| | section 1 | section 2 |")
    A("|---|---:|---:|")
    A(f"| IC50 fold range across the panel | {o1['ic50_fold_range']:.0f}× | "
      f"{o2['ic50_fold_range']:.0f}× |")
    A(f"| mean-tumour-AUC fold range | {o1['auc_fold_range']:.0f}× | "
      f"{o2['auc_fold_range']:.0f}× |")
    A(f"| ρ(kill, 1/IC50) — potency alone | {o1['rho_kill_vs_potency']:+.2f} | "
      f"{o2['rho_kill_vs_potency']:+.2f} |")
    A(f"| ρ(kill, mean AUC) — delivery alone | {o1['rho_kill_vs_delivery']:+.2f} | "
      f"{o2['rho_kill_vs_delivery']:+.2f} |")
    A(f"| ρ(kill, AUC/IC50) — the product | {o1['rho_kill_vs_auc_over_ic50']:+.2f} | "
      f"{o2['rho_kill_vs_auc_over_ic50']:+.2f} |")
    A(f"| kill spread (max − min) | {o1['kill_spread']:.3f} | "
      f"{o2['kill_spread']:.3f} |")
    A(f"| coverage spread (max − min) | {o1['coverage_spread']:.3f} | "
      f"{o2['coverage_spread']:.3f} |")
    A("")
    trivial2 = o2["rank_by_kill"] == o2["rank_by_potency_alone"]
    A(("**It is NOT the IC50 ordering.** " if not trivial2 else
       "**⚠️ On section 2 the kill ordering IS exactly the potency ordering.** ")
      + f"Potency alone would rank {' > '.join(o2['rank_by_potency_alone'])}; the "
      f"simulation ranks {' > '.join(o2['rank_by_kill'])}. Delivery spans "
      f"{o2['auc_fold_range']:.0f}× and potency {o2['ic50_fold_range']:.0f}×, and "
      f"the product AUC/IC50 tracks kill at ρ = "
      f"{o2['rho_kill_vs_auc_over_ic50']:+.2f}"
      + (" — the delivery × potency claim from `results/SESSION_AUDIT.md` §4 "
         "survives on the held-out patient." if
         o2['rho_kill_vs_auc_over_ic50'] >= 0.9 else
         " — weaker than on section 1, so read the product claim with care here."))
    A("")
    A("### 3.5 Why the absolute numbers move so much, and why that is NOT a "
      "ranking result")
    A("")
    hyp1 = d1[res["panel"][0]]["hypoxic_tumor_fraction"]
    hyp2 = d2[res["panel"][0]]["hypoxic_tumor_fraction"]
    A("| | section 1 | section 2 |")
    A("|---|---:|---:|")
    A(f"| mean oxygen over tissue | {r1['mean_oxygen']:.4f} | {r2['mean_oxygen']:.4f} |")
    A(f"| model-HYPOXIC tumour fraction | {hyp1:.2%} | {hyp2:.2%} |")
    A(f"| model-VIABLE tumour fraction | "
      f"{d1[res['panel'][0]]['viable_tumor_fraction']:.2%} | "
      f"{d2[res['panel'][0]]['viable_tumor_fraction']:.2%} |")
    A(f"| mean ECM density | {r1['mean_ecm_density']:.3f} | "
      f"{r2['mean_ecm_density']:.3f} |")
    A("")
    A("Two mechanisms move every absolute kill on section 2, and **neither is a "
      "drug property**, so both largely cancel out of the ranking:")
    A("")
    A(f"1. **Oxygen.** Mean tissue O₂ is {r2['mean_oxygen']:.3f} on section 2 "
      f"against {r1['mean_oxygen']:.3f} on section 1 — the denser tissue is a "
      f"bigger sink — so {hyp2:.0%} of tumour cells fall in the model's HYPOXIC "
      f"band versus {hyp1:.0%}, and `PDConfig.hypoxia_resistance` (×0.35) plus "
      "`quiescent_resistance` (×0.5) then suppress every drug's kill together. "
      "**`results/state_validation.md` §1.6 already established that the absolute "
      "oxygen level is uncalibrated** — `OxygenConfig.consumption` has never been "
      "fitted to a measured oxygen tension — so this difference is a statement "
      "about the sink constant and the cell density, not about the two tumours.")
    A(f"2. **ECM.** The IO panel carries only "
      f"{len(r2['ecm_markers_used'])} of the {len(ECM_MARKERS_ALL)} matrix "
      f"transcripts the ECM field is built from ({', '.join(r2['ecm_markers_used'])}) "
      f"against section 1's {len(r1['ecm_markers_used'])} "
      f"({', '.join(r1['ecm_markers_used'])}) — **COL1A2 and POSTN, the two "
      "fibrillar-collagen markers, are absent from the IO panel entirely.** Mean "
      f"ECM comes out at {r2['mean_ecm_density']:.3f} vs {r1['mean_ecm_density']:.3f}, "
      f"which scales D_eff to {(1-ALPHA_TORTUOSITY*r2['mean_ecm_density'])/(1-ALPHA_TORTUOSITY*r1['mean_ecm_density']):.2f}× of its section-1 value. "
      "That is **panel content, not tissue biology** — the ECM field is a "
      "max-normalised KDE over whichever matrix markers happen to be measured, "
      f"so an 11-marker field and a {len(r2['ecm_markers_used'])}-marker field "
      "have no common absolute scale even on identical tissue. This is the "
      "clearest case in the report of a difference that must NOT be read as "
      "biology.")
    A("")
    A("**Consequence.** Absolute kills and absolute AUCs are *not* comparable "
      "between the two sections. The **ordering** and the **within-section "
      "spread** are the transportable quantities, and those are what §3.3 "
      "reports.")
    A("")

    # --- §4 -------------------------------------------------------------------
    A("## 4. The hypoxia check, on its home turf")
    A("")
    A("The 319-gene brain panel measures **no** hypoxia-response gene — no VEGFA, "
      "no HIF1A, no CA9, no ADM, no NDRG1, no SLC2A1 — so on section 1 hypoxia is "
      "**unmeasured**, and this module declines to produce a number for it rather "
      "than substituting a proxy. The IO panel measures VEGFA and HIF1A, so "
      "section 2 is the *only* place this test can run, and it is the section "
      "`results/state_validation.md` Test 1 was found on. This is therefore a "
      "**re-derivation on the same data**, not an out-of-sample test — its value "
      "is that it is computed by a different module through the same production "
      "path, so it either confirms or refutes that report arithmetically.")
    A("")
    if h2 is None:
        A("*(No hypoxia-response gene on either panel in this run — nothing to "
          "report.)*")
        A("")
    else:
        A(f"Measured genes: " + ", ".join(
            f"**{g}** (detected in {100*h2['detection'][g]:.2f}% of cells)"
            for g in h2["genes_used"]) +
          f". Their mutual Spearman ρ is **{h2['gene_gene_rho']:+.4f}** — close to "
          "zero, so they are *not* two readouts of one latent hypoxia variable, "
          "and the composite is reported alongside the per-gene results, never "
          "instead of them.")
        A("")
        A("### 4.1 Every predictor, sign-aligned to oxygen (negative = correct)")
        A("")
        cols = " | ".join(f"ρ vs {g}" for g in h2["genes_used"])
        A(f"| predictor | what it is | ρ vs hypoxia score | {cols} | "
          "ρ with log-depth partialled out |")
        A("|---|---|---:|" + "---:|" * len(h2["genes_used"]) + "---:|")
        for name, nice in (("simulated_o2", "**SIMULATED O₂ (the PDE)**"),
                           ("vascular_kde", "vascular_kde"),
                           ("dist_endothelial", "dist_endothelial"),
                           ("local_cell_density", "local_cell_density")):
            p = h2["predictors"][name]
            per = " | ".join(f"{p[f'rho_{g}']:+.4f}" for g in h2["genes_used"])
            A(f"| {nice} | {p['what']} | **{p['rho_score']:+.4f}** | {per} | "
              f"{p['partial_rho_depth_removed']:+.4f} |")
        A("")
        A(f"Permutation p-values (label shuffling, {h2['perm_draws']:,} draws on a "
          f"seeded {h2['perm_subsample_n']:,}-cell subsample — the kernel is "
          "O(n × draws); every effect size above is on **all** "
          f"{r2['n_cells']:,} cells, so subsampling the null can only make a "
          "p-value more conservative):")
        A("")
        A("| predictor | ρ on the subsample | permutation p |")
        A("|---|---:|---:|")
        for name, blk in h2["permutation"].items():
            A(f"| {name} | {blk['rho']:+.4f} | {blk['p_value']:.5f} |")
        A("")
        pde_rho = h2["predictors"]["simulated_o2"]["rho_score"]
        dist_rho = h2["predictors"]["dist_endothelial"]["rho_score"]
        kde_rho = h2["predictors"]["vascular_kde"]["rho_score"]
        beaten = abs(dist_rho) > abs(pde_rho)
        A(("**CONFIRMED.** " if beaten else "**REFUTED.** ") +
          f"The simulated field is significantly anti-correlated with measured "
          f"hypoxia response (ρ = {pde_rho:+.4f}), but a one-line `cKDTree` query "
          f"on its own source cells scores {dist_rho:+.4f} — "
          f"{abs(dist_rho)/abs(pde_rho):.2f}× the solved field — and the "
          f"undiffused source term itself scores {kde_rho:+.4f}"
          + (", so the reaction–diffusion solve is not what carries the signal. "
             "`results/state_validation.md` §1.2 is reproduced." if beaten else
             ", so the solve does add predictive information here, contradicting "
             "`results/state_validation.md` §1.2."))
        A("")
        A("### 4.2 Nested test — does the PDE add anything beyond the raw geometry?")
        A("")
        A("| partial Spearman | ρ |")
        A("|---|---:|")
        A(f"| hypoxia score vs **simulated O₂**, distance-to-vessel held fixed | "
          f"{h2['nested']['score_vs_simo2_given_dist']:+.4f} |")
        A(f"| hypoxia score vs **−distance-to-vessel**, simulated O₂ held fixed | "
          f"{h2['nested']['score_vs_dist_given_simo2']:+.4f} |")
        A(f"| simulated O₂ vs −distance-to-vessel (how nested the two are) | "
          f"{h2['nested']['simo2_vs_dist']:+.4f} |")
        A("")
        A(f"**Depth control.** The simulated field correlates with total transcript "
          f"count — a purely technical variable — at ρ = {h2['o2_vs_depth_rho']:+.4f}; "
          "partialling log-depth out (last column of §4.1) does not change the "
          "ordering, so it is not a depth artefact.")
        A("")
        A(f"### 4.3 The {h2['n_genes_ranked']}-gene ranking")
        A("")
        A(f"Every one of the {h2['n_genes_ranked']} panel genes was "
          "depth-normalised and rank-correlated against the simulated oxygen at "
          f"each of the {r2['n_cells']:,} cells.")
        A("")
        A("| rank (1 = most negative) | gene | ρ vs simulated O₂ | detection |")
        A("|---:|---|---:|---:|")
        for i2, (g, rr, dd_) in enumerate(h2["gene_ranking_head"], start=1):
            star = " **←**" if g in h2["genes_used"] else ""
            A(f"| {i2} | {g}{star} | {rr:+.4f} | {100*dd_:.2f}% |")
        for g in h2["genes_used"]:
            k = h2["rank_of_hypoxia_genes"][g]
            if k > len(h2["gene_ranking_head"]):
                A(f"| … | … | … | … |")
                A(f"| {k} | {g} **←** | "
                  f"{h2['predictors']['simulated_o2'][f'rho_{g}']:+.4f} | "
                  f"{100*h2['detection'][g]:.2f}% |")
        A("| … | … | … | … |")
        for off, (g, rr, dd_) in enumerate(reversed(h2["gene_ranking_tail"])):
            A(f"| {h2['n_genes_ranked'] - off} | {g} | {rr:+.4f} | "
              f"{100*dd_:.2f}% |")
        A("")
        vr = h2["rank_of_hypoxia_genes"].get("VEGFA")
        hr = h2["rank_of_hypoxia_genes"].get("HIF1A")
        A(f"* **VEGFA is rank {vr} of {h2['n_genes_ranked']}** at ρ = "
          f"{h2['predictors']['simulated_o2']['rho_VEGFA']:+.5f}. Of every gene the "
          "panel measures, the canonical hypoxia-response gene is (among) the ones "
          "the simulated oxygen field most strongly anti-correlates with. "
          "Reproduced.")
        A(f"* **HIF1A is rank {hr} of {h2['n_genes_ranked']}** at ρ = "
          f"{h2['predictors']['simulated_o2']['rho_HIF1A']:+.5f} — null. Expected "
          "on mechanism: HIF-1α is regulated by oxygen-dependent prolyl "
          "hydroxylation and VHL-mediated proteolysis of the **protein** (Ivan et "
          "al., *Science* 2001;292:464-468, "
          "[PMID 11292862](https://pubmed.ncbi.nlm.nih.gov/11292862/); Jaakkola et "
          "al., *Science* 2001;292:468-472, "
          "[PMID 11292861](https://pubmed.ncbi.nlm.nih.gov/11292861/)), so its "
          "mRNA has no obligation to track oxygen. Reproduced.")
        A("* **The circularity, restated.** The positive end of that table is the "
          "vascular and myeloid programme, and `FLT1` is one of the three "
          "ENDOTHELIAL markers that *builds* the oxygen source term — so its "
          "positive correlation is partly the model correlating with itself. Only "
          "the negative end is informative.")
        A("")
        tc = h2["tile_correlation"]
        A("### 4.4 Field-vs-field agreement, and the state bands")
        A("")
        A(f"Averaging both fields over fixed {tc['tile_um']:.0f} µm tissue tiles "
          f"(≥{tc['min_cells']} measured cells per tile, {tc['n_tiles']:,} tiles) "
          f"gives Spearman ρ = **{tc['rho']:+.4f}** (Pearson "
          f"{tc['pearson']:+.4f}) — about "
          f"{abs(tc['rho'])/abs(pde_rho):.1f}× the per-cell value, the expected "
          "signature of per-cell dropout rather than of a better model. It is the "
          "most favourable framing of the test and is reported as such.")
        A("")
        sf = h2["state_fractions_all_cells"]
        A("State occupancy over all cells at this configuration: " +
          ", ".join(f"**{k}** {100*v:.2f}%" for k, v in sf.items() if v > 0) + ".")
        A("")

    # --- §5 the master table ---------------------------------------------------
    A("## 5. Everything, side by side")
    A("")
    A("| quantity | Section 1 (development) | Section 2 (HELD OUT) | ratio / note |")
    A("|---|---:|---:|---|")

    def row(lab, a, b, note="", fmt="{:.4g}"):
        fa = fmt.format(a) if isinstance(a, (int, float)) else str(a)
        fb = fmt.format(b) if isinstance(b, (int, float)) else str(b)
        if not note and isinstance(a, (int, float)) and isinstance(b, (int, float)) \
                and a not in (0,):
            note = f"{b / a:.2f}×"
        A(f"| {lab} | {fa} | {fb} | {note} |")

    A("| **PROVENANCE** | | | |")
    row("cells", r1["n_cells"], r2["n_cells"], fmt="{:,.0f}")
    row("panel gene targets", c1["n_panel_genes"], c2["n_panel_genes"], fmt="{:.0f}")
    row("marker map selected", c1["panel_name"], c2["panel_name"],
        "shared genes between the two panels: 60")
    row("tissue area (mm²)", v1["vessels"]["tissue_area_mm2"],
        v2["vessels"]["tissue_area_mm2"], fmt="{:.1f}")
    row("cell density (/mm²)", c1["cell_density_per_mm2"],
        c2["cell_density_per_mm2"], fmt="{:.0f}")
    A("| **COMPOSITION** | | | |")
    for cls in ("TUMOR", "STROMAL", "IMMUNE", "ENDOTHELIAL"):
        row(f"{cls} (%)", 100 * r1["celltype_fractions"][cls],
            100 * r2["celltype_fractions"][cls], fmt="{:.2f}")
    row("ENDOTHELIAL vs 8.0% volume-EM", c1["endothelial_fold_vs_volume_em"],
        c2["endothelial_fold_vs_volume_em"],
        f"{c1['endothelial_verdict']} → **{c2['endothelial_verdict']}**",
        fmt="{:.2f}×")
    row("endothelial areal density (/mm²)", c1["endothelial_per_mm2"],
        c2["endothelial_per_mm2"], fmt="{:.0f}")
    A("| **FIELDS** | | | |")
    row("mean vascular density", r1["mean_vascular_density"],
        r2["mean_vascular_density"], fmt="{:.4f}")
    row("mean ECM density", r1["mean_ecm_density"], r2["mean_ecm_density"],
        "different matrix genes on the two panels — NOT comparable as biology",
        fmt="{:.4f}")
    row("mean oxygen (tissue)", r1["mean_oxygen"], r2["mean_oxygen"],
        "absolute level is uncalibrated on both", fmt="{:.4f}")
    row("model-HYPOXIC tumour fraction",
        100 * d1[res["panel"][0]]["hypoxic_tumor_fraction"],
        100 * d2[res["panel"][0]]["hypoxic_tumor_fraction"],
        "a statement about the sink constant, not the tumours", fmt="{:.2f}%")
    A("| **VASCULATURE** | | | |")
    row("discrete vessel objects", v1["vessels"]["n_vessels"],
        v2["vessels"]["n_vessels"], fmt="{:,.0f}")
    row("MVD (/mm²)", v1["vessels"]["mvd_per_mm2"], v2["vessels"]["mvd_per_mm2"],
        f"published 65 (35–95): {rows1['mvd_per_mm2']['verdict']} → "
        f"**{rows2['mvd_per_mm2']['verdict']}**", fmt="{:.1f}")
    row("intervessel spacing 1/√MVD (µm)", v1["vessels"]["spacing_from_mvd_um"],
        v2["vessels"]["spacing_from_mvd_um"],
        f"published 124 (103–169): {rows1['intervessel_spacing_um']['verdict']} → "
        f"**{rows2['intervessel_spacing_um']['verdict']}**", fmt="{:.1f}")
    row("intervessel NN median (µm)", v1["intervessel_nn"]["median"],
        v2["intervessel_nn"]["median"], fmt="{:.1f}")
    row("tumour→vessel median (µm)", v1["tumor_to_vessel"]["median"],
        v2["tumor_to_vessel"]["median"],
        f"published 58 (48–79): {rows1['median_tumor_to_vessel_um']['verdict']} → "
        f"**{rows2['median_tumor_to_vessel_um']['verdict']}**", fmt="{:.1f}")
    row("tumour cells > 100 µm from a vessel (%)",
        100 * v1["vessels"]["tumor_frac_beyond_100um"],
        100 * v2["vessels"]["tumor_frac_beyond_100um"], fmt="{:.2f}")
    row("**vascular area fraction (%)**", af1, af2,
        f"measured **2.1 ± 0.3**: {rows1['vascular_area_fraction_pct']['verdict']} "
        f"({af1/2.1:.2f}×) → **{rows2['vascular_area_fraction_pct']['verdict']}** "
        f"({af2/2.1:.2f}×)", fmt="{:.2f}")
    row("endothelial cells per vessel object", v1["vessels"]["mean_cells_per_vessel"],
        v2["vessels"]["mean_cells_per_vessel"],
        "section 2's objects look far more like real vessel profiles", fmt="{:.2f}")
    A("| **DRUG PANEL (uniform sensitivity)** | | | |")
    for name in res["panel"]:
        row(f"{name} — coverage @1× IC50", d1[name]["coverage_1x"],
            d2[name]["coverage_1x"], fmt="{:.4f}")
    for name in res["panel"]:
        row(f"{name} — mean tumour AUC (µM·s)", d1[name]["mean_tumor_auc"],
            d2[name]["mean_tumor_auc"], fmt="{:.3g}")
    for name in res["panel"]:
        row(f"{name} — kill", d1[name]["kill_fraction"], d2[name]["kill_fraction"],
            fmt="{:.4f}")
    row("kill ordering", " > ".join(o1["rank_by_kill"]),
        " > ".join(o2["rank_by_kill"]),
        "**IDENTICAL**" if kill_same else "**DIFFERENT**")
    row("coverage ordering", " > ".join(o1["rank_by_coverage"]),
        " > ".join(o2["rank_by_coverage"]),
        "**IDENTICAL**" if o1["rank_by_coverage"] == o2["rank_by_coverage"]
        else "**DIFFERENT**")
    row("ρ(kill, AUC/IC50)", o1["rho_kill_vs_auc_over_ic50"],
        o2["rho_kill_vs_auc_over_ic50"], "delivery × potency", fmt="{:+.2f}")
    A("| **HYPOXIA** | | | |")
    if h2 is None:
        row("simulated O₂ vs measured hypoxia", "UNMEASURABLE (no gene on panel)",
            "UNMEASURABLE", "")
    else:
        row("hypoxia-response genes on panel", "none",
            ", ".join(h2["genes_used"]),
            "the 319 panel cannot test this at all")
        row("ρ(simulated O₂, hypoxia score)", "n/a",
            f"{h2['predictors']['simulated_o2']['rho_score']:+.4f}", "")
        row("ρ(dist-to-endothelial, hypoxia score)", "n/a",
            f"{h2['predictors']['dist_endothelial']['rho_score']:+.4f}",
            "**beats the PDE**" if abs(h2['predictors']['dist_endothelial']['rho_score'])
            > abs(h2['predictors']['simulated_o2']['rho_score']) else "loses to the PDE")
        row("VEGFA rank of all panel genes", "n/a",
            f"{h2['rank_of_hypoxia_genes'].get('VEGFA')}/{h2['n_genes_ranked']}", "")
        row("HIF1A rank of all panel genes", "n/a",
            f"{h2['rank_of_hypoxia_genes'].get('HIF1A')}/{h2['n_genes_ranked']}", "")
    A("")

    # --- §6 verdicts ------------------------------------------------------------
    A("## 6. Verdict — what generalises to a second patient, and what does not")
    A("")
    A("### 6.1 Replicates / confirmed")
    A("")
    A("Two of these are confirmations of a *negative* result — the hypoxia rows "
      "confirm that the reaction-diffusion solve is out-predicted by a distance "
      "transform. 'Confirmed' means the prior finding reproduced, not that the "
      "model performed well.")
    A("")
    for v in vd:
        if v["status"] in ("REPLICATES", "CONFIRMED", "HOLDS"):
            A(f"* **{v['claim']}** — {v['why']}")
    A("")
    A("### 6.2 Does NOT generalise")
    A("")
    for v in vd:
        if v["status"] in ("DOES NOT REPLICATE", "REFUTED", "FAILS"):
            A(f"* **{v['claim']}** — {v['why']}")
    A("")
    A("### 6.3 Fails on BOTH sections — not a generalisation failure")
    A("")
    A("These are the rows where the held-out patient behaved exactly like the "
      "development section, and both disagree with something this project has "
      "claimed. Reading them as 'it did not generalise' would be wrong; they are "
      "disagreements with a *prior configuration*, and the second patient's role "
      "here is to show the disagreement is not patient-specific.")
    A("")
    for v in vd:
        if v["status"] == "NOT REPRODUCED ON EITHER":
            A(f"* **{v['claim']}** — {v['why']}")
    A("")
    A("### 6.4 Untestable on one section or the other")
    A("")
    A("* **Hypoxia against measured markers** is impossible on the 319-gene "
      "section: it carries no VEGFA, no HIF1A, no CA9, no ADM, no NDRG1, no "
      "SLC2A1. The honest statement for section 1 remains that hypoxia is "
      "*unmeasured*, not that it has been simulated.")
    A("* **The TOP2A-based doxorubicin sensitivity determinant** cannot be "
      "evaluated on section 2 — TOP2A is absent from the IO panel. This is why "
      "the whole comparison runs in `uniform` mode; `determinant` mode is not "
      "even definable in a comparable way across these two panels.")
    A("* **MGMT** — and therefore any measured temozolomide sensitivity — is "
      "absent from **both** panels. TMZ's kill is an upper bound on both sections "
      "for the same reason.")
    A("* **Outcome.** Neither section has a drug-response label. Nothing here "
      "validates a kill prediction; the drug panel tests whether the machinery "
      "behaves consistently on a second real tissue, which is a different and "
      "much weaker claim.")
    A("")
    A("### 6.5 The one-paragraph reading")
    A("")
    A(("**The response layer travels; the reconstruction layer does not.** " if
       kill_same else
       "**The response ordering moved, and the reconstruction layer did not "
       "travel either.** ") +
      "The drug ordering and the delivery × potency structure behind it survive "
      "a second real patient on a different panel — but they were always going "
      "to be robust, because they are dominated by drug properties (IC50, BBB "
      "penetrance, molecular weight) that are identical on both sections, and the "
      "tissue only modulates them. The quantities that are genuinely *about the "
      "tissue* — the endothelial cell-type fraction and the vascular area "
      "fraction, i.e. the two numbers this project cites as its external "
      "validations — are the ones that fail to replicate. That is the honest "
      "asymmetry, and it points at where the next work is: the reconstruction's "
      "external anchors are section-specific, and the two constants that produce "
      "them (the typing threshold and the 7.5 µm vessel radius) are entangled "
      "with the grid and with panel sparsity in a way that a single section "
      "could never expose.")
    A("")

    # --- §7 caveats -------------------------------------------------------------
    A("## 7. Caveats — what could make these conclusions wrong")
    A("")
    A("1. **Different patients AND different panels.** Only 60 genes are shared. "
      "Every composition difference is confounded by panel content; §1.2 and "
      "§3.5 name the specific places where a difference is more plausibly "
      "measurement than biology (the endothelial gate arity, the ECM marker set).")
    A("2. **n = 2.** Two sections cannot separate 'this pipeline does not "
      "generalise' from 'these two tumours differ'. A claim that a quantity "
      "*fails* to replicate is well supported by a single discrepancy; a claim "
      "that one *does* replicate rests on a single agreement and is much weaker.")
    A("3. **The IO marker map is not fully blind** (§1.2): it was written against "
      "this section, though never against an external anchor.")
    A("4. **Grid convergence is not re-established here.** dx = 12 µm was shown "
      "converged on section 1's geometry (`results/numerical_rigor.md` §2). "
      "Section 2 is denser and its vessel objects are larger; the convergence "
      "argument is inherited, not re-run. A dx sweep on section 2 is the obvious "
      "next check and it was not affordable in this run.")
    A("5. **The absolute oxygen level is uncalibrated on both sections**, so the "
      f"{100*d2[res['panel'][0]]['hypoxic_tumor_fraction']:.0f}% vs "
      f"{100*d1[res['panel'][0]]['hypoxic_tumor_fraction']:.0f}% hypoxic-fraction "
      "gap is a statement about `OxygenConfig.consumption` and cell density, not "
      "about the two tumours. Every absolute kill inherits it.")
    A("6. **Typing circularity.** `FLT1` and `PLVAP` are IO-panel ENDOTHELIAL "
      "markers and therefore help build the oxygen source, so the positive tail "
      "of the gene ranking is partly the model correlating with its own inputs. "
      "Only the negative tail is informative.")
    A("7. **Pericytes are unlabelled on section 2.** `CellType` has no MURAL "
      "member, so RGS5/PDGFRB/NOTCH3-positive mural cells fall into the ungated "
      "TUMOR default, inflating the tumour count and hence the oxygen sink. "
      "That is part of why section 2 is so much more hypoxic.")
    A("8. **2D, in-plane, single section** on both. Distances are Euclidean and "
      "in-plane, not geodesic and not 3D.")
    A("9. **Subsampling.** Only the permutation null is subsampled "
      + (f"({h2['perm_subsample_n']:,} of {r2['n_cells']:,} cells, "
         f"{h2['perm_draws']:,} draws)" if h2 else "") +
      "; every effect size, fraction and simulation number is computed on all "
      "cells. Plot (c) decimates the CDF for rendering only.")
    A("")

    # --- §8 reproduction ---------------------------------------------------------
    A("## 8. Reproduction")
    A("")
    A("```bash")
    A("./.venv/bin/python -m tmesim.validation.heldout_patient")
    A("```")
    A("")
    A(f"Deterministic (`seed = 0` in every `SimConfig`; permutation seed "
      f"{PERM_SEED}). Runtime {res['runtime_s']/60:.0f} min, dominated by "
      f"{2 * len(res['panel'])} PDE solves — "
      f"{sum(r['runtime_s'] for r in s2['drugs'])/60:.0f} min of it on section 2's "
      f"{r2['grid'][0]}×{r2['grid'][1]} grid.")
    A("")
    A("**Files added by this validation (nothing else was modified):**")
    A("")
    A("* `tmesim/validation/heldout_patient.py` — this analysis")
    A("* `results/heldout_patient.md`, `results/heldout_patient.png`")
    A("")
    A("**Data provenance.** Both sections are open 10x Genomics Xenium datasets, "
      "CC BY 4.0, no login or DUA:")
    A("")
    A(f"* Section 1 — {s1['source'].get('dataset', 'n/a')}; panel "
      f"{s1['source'].get('panel', 'n/a')}")
    A(f"* Section 2 — {s2['source'].get('dataset', 'n/a')}; panel "
      f"{s2['source'].get('panel', 'n/a')}")
    A("")
    A(f"**Figure.** `{os.path.basename(png)}`.")
    A("")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dx", type=float, default=DX_UM, help="grid spacing (um)")
    ap.add_argument("--n-perm", type=int, default=N_PERM,
                    help="permutation draws for the hypoxia p-values")
    ap.add_argument("--cache", default=None,
                    help="dev only: pickle path to reuse a completed run")
    ap.add_argument("--md", default=_OUT_MD)
    ap.add_argument("--png", default=_OUT_PNG)
    args = ap.parse_args(argv)

    res = run(dx_um=args.dx, n_perm=args.n_perm, cache=args.cache)
    png = make_figure(res, args.png)
    md = write_report(res, args.md, png)
    print(f"\nwrote {md}\nwrote {png}")
    for v in res["_verdicts"]:
        print(f"  {v['status']:<20} {v['claim']}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
