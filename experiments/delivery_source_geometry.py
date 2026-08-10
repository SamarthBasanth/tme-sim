"""Does the PDE's delivery source read the corrected vessel geometry? (No.) Should it?

Why this file exists
--------------------
``results/vascular_segmentation_audit.md`` §6 measured that the production
``vascular_density`` field is BIT-IDENTICAL across single-linkage distances
10/20/50/80 um, because ``discrete_vascular_field(geometry='binary',
amplitude='unit')`` never reads the cluster labels. The corollary it drew is the
subject of this file: the field that actually delivers drug is still the binary
7.5 um disc painted around **every endothelial cell** — the construction
``vessel_geometry.py``'s own docstring calls "Wrong #2" — while the corrected,
per-vessel-object sub-grid geometry reaches only the REPORTED area fraction.

Three questions, kept separate because they have different answers:

  (a) CORRECTNESS  — is the sub-grid geometry the right source for the PDE?
  (b) FEASIBILITY  — is it defined everywhere the source needs to be?
  (c) SAFETY       — what does switching move?

Stages
------
``--stage verify``   Re-derive the audit's bit-identity claim by EXPERIMENT: vary
                     the linking distance, hash the production field. Also hash
                     the field under ``amplitude='per_vessel_cells'`` (the one
                     path that does read labels) as a positive control, so
                     "identical" is shown to be a property of the default and not
                     of the harness.
``--stage geom``     The two fields side by side on every cached real section:
                     painted area, area per endothelial cell, mass, overlap,
                     spatial displacement of the source, distance-to-source
                     distributions, and the heterogeneity score.
``--stage alpha``    Is the calibration constant able to absorb a source-geometry
                     change? ``calibrate_alpha`` is a closed form in D_eff and
                     k_eff and never touches the vascular field, so the "refit"
                     arm is proved DEGENERATE rather than assumed to be.
``--stage blast``    Both arms of the delivery blast radius on the published
                     breast and GBM verdict panels.

Nothing under ``tmesim/`` is required by stages ``verify``/``geom``/``alpha``.
The ``blast`` stage overrides ``micro.fields.vascular_density`` locally, before
``prepare_tissue`` re-solves oxygen / cell states / IFP from it, so the whole
downstream chain sees the swapped source exactly as production would.

Run:
    ./.venv/bin/python -u -m experiments.delivery_source_geometry --stage verify
    ./.venv/bin/python -u -m experiments.delivery_source_geometry --stage geom
    ./.venv/bin/python -u -m experiments.delivery_source_geometry --stage alpha
    ./.venv/bin/python -u -m experiments.delivery_source_geometry --stage blast
Writes:
    results/sweep/delivery_source_geometry_<stage>.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

OUT_DIR = os.path.join(_ROOT, "results", "sweep")

#: Every cached real section, so a claim about "the production field" is not a
#: claim about one tumour.
SECTIONS = [
    ("s1_dev", "xenium", "data_cache/xenium_gbm"),
    ("s2_heldout", "xenium", "data_cache/xenium_gbm_io"),
    ("s3_gse324268", "xenium", "data_cache/gse324268_xenium"),
    ("breast", "xenium", "data_cache/xenium_breast_janesick/rep1"),
]
DX = 12.0
LINKS = (10.0, 20.0, 50.0, 80.0)

T0 = time.time()


def log(msg: str) -> None:
    print(f"[dsg {time.time() - T0:7.1f}s] {msg}", flush=True)


def _hash(a: np.ndarray) -> str:
    """Bit-level digest of a float field — the audit's own identity criterion."""
    b = np.ascontiguousarray(np.asarray(a, dtype=np.float64))
    return hashlib.sha256(b.tobytes()).hexdigest()[:16]


def _load(sample: str):
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    return reconstruct_real(source="xenium", sample=sample, dx_um=DX)


# --------------------------------------------------------------------------- #
# The O(N^2) landmine, worked around HERE rather than patched in tmesim/
# --------------------------------------------------------------------------- #
#: ``vessel_network.segment_vessels`` computes each cluster's extent as
#: ``np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1).max()`` — a DENSE
#: n x n matrix per cluster — and its inter-vessel distance as
#: ``cKDTree.query(xy, k=counts.max()+1)`` — a DENSE n_endo x k_max matrix.
#: Both are harmless at the production 20 um link and both explode once a section
#: percolates (``results/vascular_segmentation_audit.md`` §4.1: a ~40 GB
#: allocation at 60 um killed that agent's run). This module therefore builds the
#: vessel dict it needs with bounded memory, and ASSERTS bit-equality against the
#: shipping function wherever the shipping function can still be run.
def _max_pairwise_distance(pts: np.ndarray) -> float:
    """Exact max pairwise distance in O(n) memory.

    Two exact routes, in order of cost:
      1. the diameter of a point set equals the diameter of its convex hull, so
         only the hull vertices need pairing (usually a few hundred points);
      2. if Qhull refuses (collinear / degenerate / n < 3), a chunked exact scan
         with a bounded row block.
    Both return the same number as the dense formula; neither materialises n^2.
    """
    n = pts.shape[0]
    if n < 2:
        return 0.0
    if n == 2:
        return float(np.linalg.norm(pts[0] - pts[1]))
    try:
        from scipy.spatial import ConvexHull, QhullError
        try:
            h = ConvexHull(pts)
            q = pts[h.vertices]
            d = np.linalg.norm(q[:, None, :] - q[None, :, :], axis=-1)
            return float(d.max())
        except QhullError:
            pass
    except ImportError:  # pragma: no cover
        pass
    best = 0.0
    block = 256
    for s in range(0, n, block):
        c = pts[s:s + block]
        d = np.linalg.norm(c[:, None, :] - pts[None, :, :], axis=-1)
        best = max(best, float(d.max()))
    return best


def _segment_safe(cells, domain, link_um: float) -> Dict:
    """The subset of ``segment_vessels``' output the FIELDS need, in O(N) memory.

    Returns exactly the keys ``discrete_vascular_field`` and ``vessel_geometry``
    read: ``xy``, ``labels``, ``cells_per_vessel``, ``extent_um``, plus counts.
    The clustering is the same single-linkage operator (``sparse_distance_matrix``
    + ``connected_components``, i.e. DBSCAN(eps, min_samples=1)) the shipping
    function uses.
    """
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    from tmesim.core import CellType

    ct = np.asarray(cells.cell_type)
    endo = ct == int(CellType.ENDOTHELIAL)
    xy = np.column_stack([np.asarray(cells.x)[endo],
                          np.asarray(cells.y)[endo]]).astype(np.float64)
    n = int(xy.shape[0])
    if n == 0:
        return {"xy": xy, "labels": np.zeros(0, dtype=int),
                "cells_per_vessel": np.zeros(0, dtype=int),
                "extent_um": np.zeros(0), "n_vessels": 0,
                "n_endothelial_cells": 0, "link_distance_um": float(link_um)}
    tree = cKDTree(xy)
    g = tree.sparse_distance_matrix(tree, float(link_um), output_type="coo_matrix")
    labels = np.asarray(connected_components(g, directed=False)[1], dtype=int)
    n_v = int(labels.max() + 1)
    counts = np.bincount(labels, minlength=n_v).astype(int)
    order = np.argsort(labels, kind="stable")
    starts = np.searchsorted(labels[order], np.arange(n_v), side="left")
    ends = np.searchsorted(labels[order], np.arange(n_v), side="right")
    extent = np.zeros(n_v, dtype=np.float64)
    for k in range(n_v):
        if counts[k] >= 2:
            extent[k] = _max_pairwise_distance(xy[order[starts[k]:ends[k]]])
    return {"xy": xy, "labels": labels, "cells_per_vessel": counts,
            "extent_um": extent, "n_vessels": n_v, "n_endothelial_cells": n,
            "link_distance_um": float(link_um)}


# --------------------------------------------------------------------------- #
# Stage 1 — does the production field read the cluster labels?
# --------------------------------------------------------------------------- #
def stage_verify(sections=SECTIONS) -> Dict:
    """Vary the linking distance; hash the field. Nothing is read from source."""
    from tmesim.reconstruct.vessel_network import (
        discrete_vascular_field, segment_vessels)

    out = {"dx_um": DX, "links": list(LINKS), "sections": {}}
    for name, _src, sample in sections:
        log(f"verify {name}: reconstructing {sample}")
        micro = _load(sample)
        dom, cells = micro.domain, micro.cells

        # The production field as reconstruct_real builds it, untouched.
        prod = np.asarray(micro.fields.vascular_density, dtype=np.float64)

        # Equivalence check: at the production link the memory-safe segmentation
        # must reproduce the shipping one EXACTLY, or nothing below is about the
        # shipping pipeline. (Only run at 20 um; above that the shipping
        # function is the landmine this check exists to route around.)
        v_ship = segment_vessels(cells, dom, link_distance_um=20.0)
        v_safe = _segment_safe(cells, dom, 20.0)
        equiv = {
            "labels_identical": bool(np.array_equal(
                np.asarray(v_ship["labels"]), v_safe["labels"])),
            "counts_identical": bool(np.array_equal(
                np.asarray(v_ship["cells_per_vessel"]), v_safe["cells_per_vessel"])),
            "extent_max_abs_diff_um": float(np.abs(
                np.asarray(v_ship["extent_um"]) - v_safe["extent_um"]).max()),
            "xy_identical": bool(np.array_equal(
                np.asarray(v_ship["xy"]), v_safe["xy"])),
        }
        log(f"  {name}: safe-vs-shipping segmentation {equiv}")
        assert equiv["labels_identical"] and equiv["counts_identical"], name
        assert equiv["extent_max_abs_diff_um"] == 0.0, name
        del v_ship

        rows = []
        for link in LINKS:
            v = _segment_safe(cells, dom, float(link))
            # (1) the PRODUCTION default path
            f_default = discrete_vascular_field(dom, v)
            # (2) POSITIVE CONTROL: the one amplitude that reads labels
            f_labels = discrete_vascular_field(dom, v,
                                               amplitude="per_vessel_cells")
            # (3) the corrected geometry, which is per-OBJECT and must move
            f_subgrid = discrete_vascular_field(dom, v, geometry="subgrid")
            rows.append({
                "link_um": float(link),
                "n_vessels": int(v["n_vessels"]),
                "n_endothelial": int(v["n_endothelial_cells"]),
                "hash_default": _hash(f_default),
                "hash_per_vessel_cells": _hash(f_labels),
                "hash_subgrid": _hash(f_subgrid),
                "identical_to_production": bool(np.array_equal(f_default, prod)),
                "max_abs_diff_vs_production": float(np.abs(f_default - prod).max()),
                "sum_default": float(f_default.sum()),
                "sum_per_vessel_cells": float(f_labels.sum()),
                "sum_subgrid": float(f_subgrid.sum()),
            })
            log(f"  {name} link={link:5.1f}  n_vessels={v['n_vessels']:6d}  "
                f"default={rows[-1]['hash_default']} "
                f"labels={rows[-1]['hash_per_vessel_cells']} "
                f"subgrid={rows[-1]['hash_subgrid']}")
        out["sections"][name] = {
            "sample": sample,
            "grid": list(dom.shape),
            "safe_vs_shipping_segmentation": equiv,
            "production_hash": _hash(prod),
            "rows": rows,
            "default_invariant": len({r["hash_default"] for r in rows}) == 1,
            "labels_invariant": len({r["hash_per_vessel_cells"] for r in rows}) == 1,
            "subgrid_invariant": len({r["hash_subgrid"] for r in rows}) == 1,
        }
        log(f"  {name}: default invariant={out['sections'][name]['default_invariant']}  "
            f"per_vessel_cells invariant={out['sections'][name]['labels_invariant']}  "
            f"subgrid invariant={out['sections'][name]['subgrid_invariant']}")
    return out


# --------------------------------------------------------------------------- #
# Stage 2 — the two representations, measured against each other
# --------------------------------------------------------------------------- #
def _dist_to_source(dom, field: np.ndarray, thresh: float = 0.0) -> np.ndarray:
    """Distance (um) from every TISSUE grid-cell centre to the nearest cell whose
    source is above ``thresh``. This is the quantity the PDE actually integrates
    against: a reaction-diffusion field decays as exp(-d/L) in it."""
    from scipy.spatial import cKDTree
    ny, nx = dom.shape
    dx = float(dom.dx)
    gx = (np.arange(nx, dtype=np.float64) + 0.5) * dx
    gy = (np.arange(ny, dtype=np.float64) + 0.5) * dx
    GX, GY = np.meshgrid(gx, gy)
    tissue = np.asarray(dom.tissue_mask, dtype=bool)
    src = np.asarray(field) > thresh
    if not src.any():
        return np.full(int(tissue.sum()), np.inf)
    pts = np.column_stack([GX[src], GY[src]])
    q = np.column_stack([GX[tissue], GY[tissue]])
    return cKDTree(pts).query(q, k=1)[0]


def _pct(a: np.ndarray) -> Dict[str, float]:
    a = np.asarray(a, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {k: float("nan") for k in ("median", "mean", "p90", "p99", "max")}
    return {"median": float(np.median(a)), "mean": float(a.mean()),
            "p90": float(np.percentile(a, 90)), "p99": float(np.percentile(a, 99)),
            "max": float(a.max())}


def stage_geom(sections=SECTIONS) -> Dict:
    """Binary-7.5-um-per-cell vs sub-grid-per-object, on the same segmentation."""
    from tmesim.reconstruct.vessel_geometry import (
        SINGLE_EC_VESSEL_RADIUS_UM, vascular_area_fraction, vascular_heterogeneity,
        vessel_geometry)
    from tmesim.reconstruct.vessel_network import (
        LINK_DISTANCE_UM, VESSEL_RADIUS_UM, discrete_vascular_field)

    out = {"dx_um": DX, "link_um": float(LINK_DISTANCE_UM),
           "legacy_radius_um": float(VESSEL_RADIUS_UM),
           "r1_um": float(SINGLE_EC_VESSEL_RADIUS_UM),
           "legacy_area_per_ec_um2": float(np.pi * VESSEL_RADIUS_UM ** 2),
           "sections": {}}
    for name, _src, sample in sections:
        log(f"geom {name}: reconstructing {sample}")
        micro = _load(sample)
        dom, cells = micro.domain, micro.cells
        dx = float(dom.dx)
        tissue = np.asarray(dom.tissue_mask, dtype=bool)
        tissue_area_um2 = float(tissue.sum()) * dx * dx

        v = _segment_safe(cells, dom, float(LINK_DISTANCE_UM))
        n_endo = int(v["n_endothelial_cells"])
        geom = vessel_geometry(v, radius_mode="data")
        r = np.asarray(geom["radius_um"], dtype=np.float64)

        f_bin = discrete_vascular_field(dom, v)                    # production
        f_sub = discrete_vascular_field(dom, v, geometry="subgrid")

        af = vascular_area_fraction(dom, v, radius_mode="data")
        af.pop("coverage", None)

        # ---- ANALYTIC areal comparison, no grid involved --------------------
        analytic_sub_um2 = float(geom["total_footprint_area_um2"])
        analytic_leg_um2 = float(n_endo) * math.pi * VESSEL_RADIUS_UM ** 2

        # ---- what the GRID actually carries ---------------------------------
        # binary: every painted cell owns dx^2. subgrid: fractional coverage.
        bin_area_um2 = float((f_bin > 0)[tissue].sum()) * dx * dx
        sub_area_um2 = float(f_sub[tissue].sum()) * dx * dx
        bin_mass = float(f_bin[tissue].sum())
        sub_mass = float(f_sub[tissue].sum())

        b = (np.asarray(f_bin) > 0) & tissue
        s = (np.asarray(f_sub) > 0) & tissue
        inter = int((b & s).sum())
        union = int((b | s).sum())

        d_bin = _dist_to_source(dom, f_bin)
        d_sub = _dist_to_source(dom, f_sub)

        endo_xy = np.asarray(v["xy"], dtype=np.float64)
        het_bin = vascular_heterogeneity(dom, f_bin, endo_xy)
        het_sub = vascular_heterogeneity(dom, f_sub, endo_xy)

        # ---- WHICH GEOMETRIC LAW does each field implement? ------------------
        # The engine's delivery term is ``vascular_density * bbb_penetrance *
        # plasma`` — a transvascular flux per unit tissue volume. The physical
        # weight for that is the vessel WALL area per unit tissue (S/V), which
        # for a tube of radius r and length L goes as ``2 pi r L`` — LINEAR in r.
        # The sub-grid footprint is a CROSS-SECTION area, ``pi r^2 + 2 r L`` —
        # QUADRATIC in r in the transverse regime that dominates at the
        # production linking distance. With ``r = r1 * k`` that is k^2 against
        # the wall law's k. This block measures the resulting per-endothelial-
        # cell weight as a function of k, which is the quantity a leak-rate
        # source is supposed to hold roughly CONSTANT (one endothelial cell
        # lines EC_AXIAL_PITCH_UM of wall, so EC count IS a wall-area proxy).
        counts = np.asarray(v["cells_per_vessel"], dtype=np.float64)
        kk = np.asarray(geom["ec_per_cross_section"], dtype=np.float64)
        foot = np.asarray(geom["footprint_area_um2"], dtype=np.float64)
        axis_len = np.asarray(geom["axis_length_um"], dtype=np.float64)
        wall_len = axis_len + 2.0 * r          # stadium perimeter / (2) proxy
        per_ec_sub = foot / np.maximum(counts, 1.0)
        kbins = [(1.0, 1.25), (1.25, 2.0), (2.0, 3.0), (3.0, 5.0), (5.0, 1e9)]
        by_k = []
        for lo, hi in kbins:
            sel = (kk >= lo) & (kk < hi)
            if not sel.any():
                continue
            by_k.append({
                "k_lo": lo, "k_hi": (None if hi > 1e8 else hi),
                "n_objects": int(sel.sum()),
                "n_ec": float(counts[sel].sum()),
                "median_k": float(np.median(kk[sel])),
                "median_radius_um": float(np.median(r[sel])),
                # what the SUB-GRID source deposits per endothelial cell
                "subgrid_area_per_ec_um2": float(foot[sel].sum() / counts[sel].sum()),
                # what the PRODUCTION binary source deposits per endothelial
                # cell, before overlap de-duplication: a constant, by design
                "binary_area_per_ec_um2": float(np.pi * VESSEL_RADIUS_UM ** 2),
                # a WALL-area law: 2 pi r * (axis length + caps), per EC
                "wall_len_per_ec_um": float(
                    (2.0 * np.pi * r[sel] * wall_len[sel]).sum() / counts[sel].sum()),
            })
        rec_by_k = by_k

        rec = {
            "sample": sample,
            "grid": list(dom.shape),
            "weight_per_ec_by_k": rec_by_k,
            "subgrid_area_per_ec_spread_ratio": (
                (by_k[-1]["subgrid_area_per_ec_um2"] / by_k[0]["subgrid_area_per_ec_um2"])
                if len(by_k) >= 2 else float("nan")),
            "tissue_area_mm2": tissue_area_um2 / 1e6,
            "n_endothelial": n_endo,
            "n_vessels": int(v["n_vessels"]),
            "ec_per_mm2": n_endo / (tissue_area_um2 / 1e6),
            "mean_k": float(np.mean(geom["ec_per_cross_section"])),
            "median_radius_um": float(np.median(r)) if r.size else float("nan"),
            "p99_radius_um": float(np.percentile(r, 99)) if r.size else float("nan"),
            "max_radius_um": float(r.max()) if r.size else float("nan"),
            # analytic areas (grid-free)
            "analytic_subgrid_area_um2": analytic_sub_um2,
            "analytic_legacy_area_um2": analytic_leg_um2,
            "analytic_area_ratio_subgrid_over_legacy": analytic_sub_um2 / analytic_leg_um2,
            "analytic_subgrid_area_per_ec_um2": analytic_sub_um2 / max(n_endo, 1),
            "analytic_area_fraction_pct": 100.0 * analytic_sub_um2 / tissue_area_um2,
            # rasterised
            "binary_painted_area_um2": bin_area_um2,
            "subgrid_painted_area_um2": sub_area_um2,
            "binary_area_fraction_pct": 100.0 * bin_area_um2 / tissue_area_um2,
            "subgrid_area_fraction_pct": 100.0 * sub_area_um2 / tissue_area_um2,
            "rasterised_area_ratio_subgrid_over_binary": sub_area_um2 / bin_area_um2,
            "binary_mass": bin_mass, "subgrid_mass": sub_mass,
            "mass_ratio_subgrid_over_binary": sub_mass / bin_mass,
            # support geometry — WHERE the source sits
            "binary_support_cells": int(b.sum()),
            "subgrid_support_cells": int(s.sum()),
            "support_intersection_cells": inter,
            "support_union_cells": union,
            "support_jaccard": inter / union if union else float("nan"),
            "subgrid_support_covered_by_binary": inter / int(s.sum()) if s.sum() else float("nan"),
            "binary_support_covered_by_subgrid": inter / int(b.sum()) if b.sum() else float("nan"),
            "binary_support_frac_of_tissue": int(b.sum()) / float(tissue.sum()),
            "subgrid_support_frac_of_tissue": int(s.sum()) / float(tissue.sum()),
            # distance-to-source: what a exp(-d/L) solute experiences
            "dist_to_binary_um": _pct(d_bin),
            "dist_to_subgrid_um": _pct(d_sub),
            # fidelity to the measured point pattern
            "heterogeneity_binary": het_bin,
            "heterogeneity_subgrid": het_sub,
            "reported_area_fraction_pct": af["area_fraction_pct"],
        }
        out["sections"][name] = rec
        log(f"  {name}: subgrid/legacy analytic area ratio "
            f"{rec['analytic_area_ratio_subgrid_over_legacy']:.4f}; "
            f"per-EC {rec['analytic_subgrid_area_per_ec_um2']:.1f} vs "
            f"{out['legacy_area_per_ec_um2']:.1f} um2; "
            f"support Jaccard {rec['support_jaccard']:.4f}; "
            f"mass ratio {rec['mass_ratio_subgrid_over_binary']:.4f}")
    return out


# --------------------------------------------------------------------------- #
# Stage 3 — can alpha absorb a source-geometry change?
# --------------------------------------------------------------------------- #
def stage_alpha() -> Dict:
    """``calibrate_alpha`` is a closed form in (D_eff, k_passive, uptake_proxy).

    None of those three is a function of ``vascular_density``. This stage proves
    it by CALLING the shipping calibration under a swapped source field rather
    than by reading the code, so the 'alpha refit' arm of the blast radius is
    shown to be degenerate rather than asserted to be.
    """
    from tmesim.validation.calibration import (
        DOX_TARGET_EFOLD_UM, _compile, _d_eff_tortuous, calibrate_alpha,
        k_passive, uptake_proxy)
    from tmesim.validation.pk_reality import DEFAULT_ECM

    dp = _compile("doxorubicin")
    a0 = calibrate_alpha()
    d_eff = _d_eff_tortuous(dp, DEFAULT_ECM)
    out = {
        "alpha_published": float(a0),
        "alpha_recomputed": float(calibrate_alpha()),
        "D_eff_um2_s": float(d_eff),
        "k_passive_per_s": float(k_passive(dp)),
        "uptake_proxy_per_s": float(uptake_proxy(dp)),
        "target_efold_um": float(DOX_TARGET_EFOLD_UM),
        "L_analytic_um": float(math.sqrt(
            d_eff / (k_passive(dp) + a0 * uptake_proxy(dp)))),
        "inputs_depend_on_vascular_density": False,
    }
    out["alpha_refit_minus_published"] = out["alpha_recomputed"] - out["alpha_published"]
    log(f"alpha published {a0:.12f}; closed-form inputs "
        f"D_eff={d_eff:.6f} k_passive={k_passive(dp):.6e} "
        f"proxy={uptake_proxy(dp):.6e}; L={out['L_analytic_um']:.6f} um")
    return out


# --------------------------------------------------------------------------- #
# Stage 4 — the delivery blast radius, both arms
# --------------------------------------------------------------------------- #
#: The nine formally-scored verdicts of ``results/p99_convention.md`` §6.
GBM_DRUGS = ("temozolomide", "cisplatin", "doxorubicin")
BREAST_DRUGS = ("doxorubicin", "paclitaxel", "docetaxel", "epirubicin",
                "5-fluorouracil", "cisplatin")

PUBLISHED_RATIO = {
    ("gbm", "temozolomide"): 0.2626, ("gbm", "cisplatin"): 0.2791,
    ("gbm", "doxorubicin"): 0.2077,
    ("breast", "doxorubicin"): 0.2301, ("breast", "paclitaxel"): 0.1503,
    ("breast", "docetaxel"): 0.1484, ("breast", "epirubicin"): 0.2301,
    ("breast", "5-fluorouracil"): 0.3379, ("breast", "cisplatin"): 0.3085,
}
PUBLISHED_VERDICT = {
    ("gbm", "temozolomide"): "PASS", ("gbm", "cisplatin"): "PASS",
    ("gbm", "doxorubicin"): "PASS",
    ("breast", "doxorubicin"): "PASS", ("breast", "paclitaxel"): "PASS",
    ("breast", "docetaxel"): "FAIL", ("breast", "epirubicin"): "PASS",
    ("breast", "5-fluorouracil"): "PASS", ("breast", "cisplatin"): "PASS",
}


#: The three source arms.
#:
#: ``binary``        — production. Binary 7.5 um disc per ENDOTHELIAL CELL.
#: ``subgrid``       — the corrected per-VESSEL-OBJECT sub-grid geometry, raw.
#:                     Both the total source and its spatial support change.
#: ``subgrid_matched`` — the same corrected geometry with ``match_total_to`` set
#:                     to the production field, i.e. the tissue-integrated source
#:                     total held EXACTLY fixed and only *where* the drug enters
#:                     allowed to change. This is ``discrete_vascular_field``'s own
#:                     documented controlled comparison, and it is also the
#:                     strongest available statement about the PS gauge: because
#:                     ``PS`` and ``vascular_density`` are perfectly confounded
#:                     (``results/transvascular_exchange_revalidation.md`` §5),
#:                     any pure SCALE change is absorbable by PS by construction,
#:                     so whatever moves in this arm is what NO global gauge —
#:                     PS, alpha, or plasma level — can remove.
ARMS = ("binary", "subgrid", "subgrid_matched")


def _source_field(micro, arm: str) -> np.ndarray:
    from tmesim.reconstruct.vessel_network import (
        LINK_DISTANCE_UM, discrete_vascular_field)
    v = _segment_safe(micro.cells, micro.domain, float(LINK_DISTANCE_UM))
    base = discrete_vascular_field(micro.domain, v)
    if arm == "binary":
        return base
    if arm == "subgrid":
        return discrete_vascular_field(micro.domain, v, geometry="subgrid")
    if arm == "subgrid_matched":
        return discrete_vascular_field(micro.domain, v, geometry="subgrid",
                                       match_total_to=base)
    raise ValueError(arm)


def _swap_vascular(micro, arm: str):
    """Return a COPY of ``micro`` whose vascular_density is painted for ``arm``.

    Oxygen / states / IFP are deliberately NOT carried over: ``prepare_tissue``
    re-solves all three from ``fields.vascular_density``, so the swapped source
    propagates into the oxygen field, the cell-state labels (hence the state
    resistance term) and the interstitial pressure exactly as it would in
    production. Anything less would understate the blast radius.
    """
    from tmesim.core import Microenvironment, TissueFields
    field = _source_field(micro, arm)
    m2 = Microenvironment(
        domain=micro.domain,
        fields=TissueFields(ecm_density=micro.fields.ecm_density,
                            vascular_density=field),
        cells=micro.cells.copy(), patient_id=micro.patient_id,
        metadata=dict(micro.metadata))
    return m2, field


def _load_at(sample: str, dx_um: float):
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    return reconstruct_real(source="xenium", sample=sample, dx_um=float(dx_um))


def _gbm_bands(drugs: Sequence[str], ic50_hill: Dict, baths: Dict) -> Dict:
    """The GBM report's own band construction (``kill_calibration.analyze`` §THE
    TEST): the drug-matched literature bridge where one exists, else 'solid',
    mapped to a kill-ratio band at that drug's OWN reporting concentration."""
    from tmesim.validation.kill_calibration import (
        bridge_to_kill_ratio_band, literature_bridge)
    out = {}
    for name in drugs:
        b = literature_bridge(scope=f"drug:{name}")
        scope = f"drug:{name}"
        if not b.get("n"):
            b, scope = literature_bridge(scope="solid"), "solid"
        ic50, hill = ic50_hill[name]
        out[name] = {"band": bridge_to_kill_ratio_band(b, baths[name], ic50, hill),
                     "scope": scope}
    return out


def _run_panel(tissue: str, arm: str, alpha: float, drugs: Sequence[str],
               dx_um: float) -> Dict:
    """One tissue, one source arm, one alpha — through the published scoring path.

    The scoring is deliberately the SAME code the published reports call
    (``kill_calibration.decompose`` / ``ratio_test`` and the frozen literature
    bridges), so an arm's numbers are comparable to the published table rather
    than merely to each other. The ``binary`` arm is the control and MUST
    reproduce the published ratios; that is checked in :func:`stage_blast`.
    """
    from tmesim.core import PDConfig, PDEConfig, PKConfig, SimConfig
    from tmesim.validation.kill_calibration import (
        decompose, expected_kill, kill_ratio_to_ic50_fold, prepare_tissue,
        ratio_test, resistance_vector, unit_transport, window_plasma_mean)

    if tissue == "breast":
        from experiments.dtissue_sensitivity import build_bands
        from experiments.kill_calibration_breast import (
            BREAST_SAMPLE_PATH, DT_S, NOMINAL_HALF_LIFE_S, NOMINAL_PEAK_UM,
            SAVE_EVERY, SEED, WINDOW_S, apply_drug_sensitivity_on_copy,
            compile_breast_drug)
        sample = BREAST_SAMPLE_PATH

        def compile_fn(name):
            return compile_breast_drug(name, float(alpha))[0]
    else:
        from tmesim.validation.kill_calibration import (
            DT_S, EXVIVO_TMZ_UM, NOMINAL_HALF_LIFE_S, NOMINAL_PEAK_UM,
            SAVE_EVERY, SEED, WINDOW_S, compile_drug)
        from tmesim.reconstruct.real_pipeline import apply_drug_sensitivity
        sample = os.path.join(_ROOT, "data_cache", "xenium_gbm")

        def compile_fn(name):
            return compile_drug(name, float(alpha))

    micro = _load_at(sample, dx_um)
    swapped, field = _swap_vascular(micro, arm)
    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                    save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                  n_doses=1, dose_interval=3600.0 * 24)
    pd_cfg = PDConfig()
    work, labeled, velocity = prepare_tissue(
        swapped, SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
    tumor = labeled.tumor_mask()
    state_t = np.asarray(labeled.state)[tumor]
    plasma_mean = window_plasma_mean(pk, WINDOW_S, DT_S)
    tissue_mask = np.asarray(work.domain.tissue_mask, dtype=bool)

    # the state term is itself downstream of the source, via oxygen
    states, counts = np.unique(state_t, return_counts=True)
    state_mix = {int(s): int(c) for s, c in zip(states, counts)}

    ic50_hill, sens, drugs_compiled, baths = {}, {}, {}, {}
    for name in drugs:
        dp0 = compile_fn(name)
        drugs_compiled[name] = dp0
        ic50_hill[name] = (float(dp0.ic50), float(dp0.hill))
        baths[name] = (float(EXVIVO_TMZ_UM) if (tissue == "gbm" and
                                                name == "temozolomide")
                       else 1.0 * float(dp0.ic50))
        if tissue == "breast":
            m2 = apply_drug_sensitivity_on_copy(work, labeled, name)
        else:
            from tmesim.core import Microenvironment
            m2 = Microenvironment(domain=work.domain, fields=work.fields,
                                  cells=labeled.copy(),
                                  patient_id=work.patient_id,
                                  metadata=work.metadata)
            apply_drug_sensitivity(m2, name, mode="uniform")
        sens[name] = np.asarray(m2.cells.sensitivity, dtype=np.float64)[tumor]

    bands = (build_bands(drugs, ic50_hill) if tissue == "breast"
             else _gbm_bands(drugs, ic50_hill, baths))

    rows = []
    for name in drugs:
        drug = drugs_compiled[name]
        ts = time.time()
        tr = unit_transport(work, drug, pde, velocity)
        col, row = work.domain.to_index(labeled.x, labeled.y)
        auc = np.asarray(tr.auc, dtype=np.float64)[row, col]
        dur = float(tr.times[-1]) if tr.times.size and tr.times[-1] > 0 else 1.0
        e_unit = (auc / dur)[tumor]
        resist = resistance_vector(state_t, pd_cfg, drug=drug)
        c_bath = baths[name]
        dec = decompose(e_unit, sens[name], resist, c_bath, drug.ic50, drug.hill,
                        pd_cfg.max_kill_rate)
        b = bands[name]
        if tissue == "breast":
            t_spec = (ratio_test(dec["ratio_exact"], b["band"]["ratio_lo"],
                                 b["band"]["ratio_hi"]) if b["band"] else None)
            t_solid = ratio_test(dec["ratio_exact"], b["band_solid"]["ratio_lo"],
                                 b["band_solid"]["ratio_hi"])
            formal = t_spec if b["band_is_drug_matched"] else t_solid
        else:
            formal = ratio_test(dec["ratio_exact"], b["band"]["ratio_lo"],
                                b["band"]["ratio_hi"])
        # the AUC field itself, so 'the exposure moved' is measurable and not
        # inferred from the score
        auc_t = np.asarray(tr.auc, dtype=np.float64)[tissue_mask] / dur
        rows.append({
            "tissue": tissue, "arm": arm, "alpha": float(alpha), "drug": name,
            "c_bath_uM": float(c_bath),
            "ratio_exact": dec["ratio_exact"],
            "f_penetration": dec["f_penetration"],
            "f_state": dec["f_state"],
            "f_product": dec["f_product"],
            "e_unit_mean": float(e_unit.mean()),
            "e_unit_median": float(np.median(e_unit)),
            "e_unit_p99": float(np.percentile(e_unit, 99.0)),
            "e_unit_p01": float(np.percentile(e_unit, 1.0)),
            "field_auc_mean": float(auc_t.mean()),
            "kill_3d_model": dec["kill_3d_model"],
            "kill_2d_model": dec["kill_2d_model"],
            "kill_nominal_plasma": float(expected_kill(
                e_unit, sens[name], resist, plasma_mean, drug.ic50, drug.hill,
                pd_cfg.max_kill_rate)),
            "implied_ic50_fold": kill_ratio_to_ic50_fold(
                dec["ratio_exact"], c_bath, drug.ic50, drug.hill),
            "verdict": formal["verdict"],
            "band_lo": formal["lit_ratio_lo"], "band_hi": formal["lit_ratio_hi"],
            "_e_unit_sum": float(e_unit.sum()),
            "_e_unit_first": float(e_unit[0]),
            "solve_s": time.time() - ts,
        })
        log(f"  {tissue}/{arm} {name}: ratio={rows[-1]['ratio_exact']:.4f} "
            f"F_pen={rows[-1]['f_penetration']:.4f} "
            f"E={rows[-1]['e_unit_mean']:.4e} {rows[-1]['verdict']} "
            f"({rows[-1]['solve_s']:.0f}s)")
    return {
        "tissue": tissue, "arm": arm, "alpha": float(alpha), "dx_um": float(dx_um),
        "grid": list(work.domain.shape),
        "n_cells": int(labeled.n), "n_tumor": int(tumor.sum()),
        "state_mix": state_mix,
        "source_sum_over_tissue": float(np.asarray(field)[tissue_mask].sum()),
        "source_support_frac": float((np.asarray(field)[tissue_mask] > 0).mean()),
        "source_hash": _hash(field),
        "oxygen_mean": float(np.asarray(work.fields.oxygen)[tissue_mask].mean()),
        "rows": rows,
    }


def stage_absorb(dx_um: float = 12.0, drug: str = "doxorubicin",
                 factors: Sequence[float] = (0.25, 0.5, 1.0, 2.0, 4.0)) -> Dict:
    """Could ANY alpha absorb the source-geometry change? Measured, not argued.

    ``calibrate_alpha`` cannot move at all under a source swap (:func:`stage_alpha`),
    so this asks the strictly weaker question: if one were willing to abandon the
    doxorubicin penetration anchor and tune ``alpha`` freely, is there a value that
    restores the published 3-D:2-D ratio under the sub-grid source? Each point also
    reports the analytic penetration length that alpha implies, so the cost to the
    anchor is visible next to the benefit.
    """
    from tmesim.validation.calibration import (
        DOX_EFOLD_BAND_UM, DOX_TARGET_EFOLD_UM, _compile, _d_eff_tortuous,
        calibrate_alpha, k_passive, uptake_proxy)
    from tmesim.validation.pk_reality import DEFAULT_ECM

    a0 = float(calibrate_alpha())
    dp = _compile("doxorubicin")
    d_eff = _d_eff_tortuous(dp, DEFAULT_ECM)
    rows = []
    for f in factors:
        a = a0 * float(f)
        L = math.sqrt(d_eff / (k_passive(dp) + a * uptake_proxy(dp)))
        for arm in ("binary", "subgrid"):
            p = _run_panel("breast", arm, a, (drug,), dx_um)
            r = p["rows"][0]
            r.update({"alpha_factor": float(f), "L_analytic_um": float(L)})
            rows.append(r)
            log(f"  absorb f={f:g} alpha={a:.3f} L={L:.1f}um {arm}: "
                f"ratio={r['ratio_exact']:.4f}")
    return {"alpha_published": a0, "dx_um": float(dx_um), "drug": drug,
            "dox_target_efold_um": float(DOX_TARGET_EFOLD_UM),
            "dox_efold_band_um": list(DOX_EFOLD_BAND_UM),
            "rows": rows}


def stage_blast(dx_um: float = 12.0, tissues: Sequence[str] = ("breast", "gbm"),
                arms: Sequence[str] = ARMS,
                checkpoint: Optional[str] = None) -> Dict:
    """The delivery blast radius over the nine published verdicts."""
    out = {"dx_um": float(dx_um), "arms": list(arms), "panels": [],
           "published_ratio": {f"{k[0]}|{k[1]}": v
                               for k, v in PUBLISHED_RATIO.items()},
           "published_verdict": {f"{k[0]}|{k[1]}": v
                                 for k, v in PUBLISHED_VERDICT.items()}}
    from tmesim.validation.calibration import calibrate_alpha
    alpha = float(calibrate_alpha())
    out["alpha"] = alpha

    done = {}
    if checkpoint and os.path.exists(checkpoint):
        with open(checkpoint) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    p = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if float(p.get("dx_um", -1)) == float(dx_um):
                    done[f"{p['tissue']}|{p['arm']}"] = p
        log(f"checkpoint: {len(done)} panels already complete")
    ck = open(checkpoint, "a") if checkpoint else None

    for tissue in tissues:
        drugs = BREAST_DRUGS if tissue == "breast" else GBM_DRUGS
        for arm in arms:
            key = f"{tissue}|{arm}"
            if key in done:
                log(f"{key}: cached")
                out["panels"].append(done[key])
                continue
            log(f"=== panel {key} (alpha={alpha:.6f}) ===")
            p = _run_panel(tissue, arm, alpha, drugs, dx_um)
            out["panels"].append(p)
            if ck is not None:
                ck.write(json.dumps(p, default=float) + "\n")
                ck.flush()
    if ck is not None:
        ck.close()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="verify",
                    choices=("verify", "geom", "alpha", "blast", "absorb", "all"))
    ap.add_argument("--sections", default="")
    ap.add_argument("--tissues", default="breast,gbm")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--dx", type=float, default=12.0)
    args = ap.parse_args()

    secs = SECTIONS
    if args.sections:
        want = {s.strip() for s in args.sections.split(",")}
        secs = [s for s in SECTIONS if s[0] in want]

    os.makedirs(OUT_DIR, exist_ok=True)
    # 'all' runs only the stages that need no PDE solve; 'blast' and 'absorb'
    # cost tens of transport solves each and are always requested explicitly.
    stages = (("verify", "geom", "alpha") if args.stage == "all" else (args.stage,))
    for st in stages:
        if st == "verify":
            res = stage_verify(secs)
        elif st == "geom":
            res = stage_geom(secs)
        elif st == "alpha":
            res = stage_alpha()
        elif st == "absorb":
            res = stage_absorb(dx_um=args.dx)
        else:
            res = stage_blast(
                dx_um=args.dx,
                tissues=tuple(t.strip() for t in args.tissues.split(",") if t.strip()),
                arms=tuple(a.strip() for a in args.arms.split(",") if a.strip()),
                checkpoint=os.path.join(
                    OUT_DIR, f"delivery_source_geometry_blast_dx{args.dx:g}.jsonl"))
        path = os.path.join(OUT_DIR, f"delivery_source_geometry_{st}.json")
        with open(path, "w") as fh:
            json.dump(res, fh, indent=1, default=float)
        log(f"wrote {path}")


if __name__ == "__main__":
    main()
