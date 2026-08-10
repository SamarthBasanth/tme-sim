"""Audit of the vascular area-fraction / MVD chain — diagnosis, not tuning.

Measures, on every cached real GBM section, the quantities needed to decide WHY
``vessel_network.segmentation_report`` reports a vascular area fraction of
1.25 % (development section) and 10.98 % (held-out section) against a measured
2.1 +/- 0.3 %, and an MVD of 111 / 233 /mm^2 against a measured 65 (35-95).

Nothing here tunes anything and no constant is changed. Every number printed is
a measurement of the existing frozen pipeline.

Three blocks:

1. ``production_stats``  — the chain at the frozen 20 um linking distance:
   endothelial density, MVD, per-object ``k`` and radius distributions, the
   painted area per endothelial cell, and an INTERNAL consistency check (does
   the painted footprint fit inside the endothelial nuclei that generated it?).
   None of these reference a published target.

2. ``link_sweep`` — MVD, mean ``k``, giant-component fraction and area fraction
   versus the single-linkage distance. The linking distance is the parameter the
   measurement's own counting convention constrains (PMID 27583442: "if
   interconnected, stained objects were considered a single vessel"), and this
   package already carries the measured endothelial axial pitch that fixes it
   (``vessel_geometry.EC_AXIAL_PITCH_UM`` = 50 um, PMID 38723085).

3. ``blast_radius`` — is the PRODUCTION delivery source field affected by any of
   this? (It is painted by ``discrete_vascular_field(geometry='binary',
   amplitude='unit')``, which uses endothelial positions only.)

Run:  ./.venv/bin/python -m experiments.vascular_segmentation_audit
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from tmesim.reconstruct.real_fields import vascular_density_from_cells
from tmesim.reconstruct.vessel_geometry import (
    EC_AXIAL_PITCH_UM,
    MEAN_VESSEL_RADIUS_UM,
    SINGLE_EC_VESSEL_RADIUS_UM,
    _subgrid_coverage,
    vascular_area_fraction,
    vessel_geometry,
)
from tmesim.reconstruct.vessel_network import (
    LINK_DISTANCE_UM,
    discrete_vascular_field,
    segment_vessels,
)

SECTIONS = [
    ("s1_dev", "data_cache/xenium_gbm"),
    ("s2_heldout", "data_cache/xenium_gbm_io"),
    ("s3_gse324268", "data_cache/gse324268_xenium"),
]
DX = 12.0
# Bounded: above ~50 um the denser sections percolate and segment_vessels'
# per-cluster O(n^2) extent computation becomes intractable, which is itself
# reported as a finding rather than worked around.
LINKS = (10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0)
T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _labels(xy: np.ndarray, eps: float) -> np.ndarray:
    """Same single-linkage as segment_vessels, without the O(n^2) extent step."""
    tree = cKDTree(xy)
    g = tree.sparse_distance_matrix(tree, float(eps), output_type="coo_matrix")
    return np.asarray(connected_components(g, directed=False)[1], dtype=int)


def audit(name: str, sample: str) -> dict:
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    log(f"{name}: reconstructing {sample} at dx={DX}")
    micro = reconstruct_real(source="xenium", sample=sample, dx_um=DX)
    dom, cells = micro.domain, micro.cells
    tissue = np.asarray(dom.tissue_mask, dtype=bool)
    tissue_area_mm2 = float(tissue.sum()) * DX ** 2 / 1e6
    n_cells = int(np.asarray(cells.cell_type).size)

    log(f"{name}: segmenting at the production {LINK_DISTANCE_UM:g} um")
    v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
    xy = np.asarray(v["xy"], dtype=np.float64)
    n_endo = int(v["n_endothelial_cells"])
    counts = np.asarray(v["cells_per_vessel"], dtype=np.float64)
    extent = np.asarray(v["extent_um"], dtype=np.float64)

    # nearest OTHER endothelial cell: the length scale the 20 um link must span
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]

    geom = vessel_geometry(v, radius_mode="data")
    r = np.asarray(geom["radius_um"], dtype=np.float64)
    k = np.asarray(geom["ec_per_cross_section"], dtype=np.float64)
    foot = np.asarray(geom["footprint_area_um2"], dtype=np.float64)

    af = vascular_area_fraction(dom, v, radius_mode="data")
    af.pop("coverage", None)

    # ---- INTERNAL consistency, no published number involved -----------------
    # vessel_geometry's OWN documented model states that a transverse cut has
    # "E ~= 2r (the cells sit across one lumen)". So the object's own extent E is
    # an independent estimate of its radius, which the code computes and then
    # discards in favour of r = r1*k. Where r1*k > E/2 the painted disc is wider
    # than the nuclei that generated it -- geometrically impossible for a single
    # vessel profile. For k nuclei tiling a circle of radius r the largest
    # pairwise chord is 2 r sin(pi*floor(k/2)/k), so the extent-implied radius is
    # r_E = E / (2 sin(pi*floor(k/2)/k)); this is exact for k=2 and tightest for
    # even k. Objects with a single nucleus carry no extent information.
    kk = np.maximum(np.rint(counts), 1.0)
    chord = 2.0 * np.sin(np.pi * np.floor(kk / 2.0) / np.maximum(kk, 1.0))
    r_from_extent = np.where(chord > 0, extent / np.maximum(chord, 1e-9), np.inf)
    multi = counts >= 2
    ratio_r = np.where(np.isfinite(r_from_extent) & (r_from_extent > 0),
                       r / np.maximum(r_from_extent, 1e-9), np.nan)
    rm = ratio_r[multi & np.isfinite(ratio_r)]

    painted_um2 = af["area_fraction_pct"] / 100.0 * tissue_area_mm2 * 1e6
    per_ec_um2 = painted_um2 / max(n_endo, 1)

    out = {
        "name": name,
        "sample": sample,
        "production_stats": {
            "n_cells": n_cells,
            "n_endothelial": n_endo,
            "endothelial_pct_of_cells": 100.0 * n_endo / max(n_cells, 1),
            "tissue_area_mm2": tissue_area_mm2,
            "ec_density_per_mm2": n_endo / tissue_area_mm2,
            "n_vessels": int(v["n_vessels"]),
            "mvd_per_mm2": float(v["mvd_per_mm2"]),
            "singleton_fraction": float(v["singleton_fraction"]),
            "mean_cells_per_vessel": float(v["mean_cells_per_vessel"]),
            "median_nn_endo_um": float(np.median(nn)),
            "frac_endo_with_neighbour_within_link": float(
                (nn <= LINK_DISTANCE_UM).mean()),
            "frac_endo_with_neighbour_within_axial_pitch": float(
                (nn <= EC_AXIAL_PITCH_UM).mean()),
            "mean_k": float(k.mean()),
            "median_k": float(np.median(k)),
            "p95_k": float(np.percentile(k, 95)),
            "max_k": float(k.max()),
            "frac_k_gt2": float((k > 2).mean()),
            "frac_k_gt6": float((k > 6).mean()),
            # the k^2 amplification: painted area per EC relative to the k=1 floor
            "k2_amplification_mean_k2_over_mean_k": float(
                (k ** 2 * 1.0).sum() / max(counts.sum(), 1.0)),
            "median_radius_um": float(np.median(r)),
            "mean_radius_um": float(r.mean()),
            "p99_radius_um": float(np.percentile(r, 99)),
            "max_radius_um": float(r.max()),
            "area_fraction_pct": float(af["area_fraction_pct"]),
            "analytic_area_fraction_pct": float(af["analytic_area_fraction_pct"]),
            "painted_um2_per_ec": per_ec_um2,
            "equiv_radius_per_ec_um": float(np.sqrt(per_ec_um2 / np.pi)),
            # internal consistency
            "frac_multi_objects_overpainted": (
                float((rm > 1.0).mean()) if rm.size else float("nan")),
            "median_r_over_r_from_extent_multi": (
                float(np.median(rm)) if rm.size else float("nan")),
            "p90_r_over_r_from_extent_multi": (
                float(np.percentile(rm, 90)) if rm.size else float("nan")),
            "n_multi_objects": int(multi.sum()),
        },
    }

    log(f"{name}: link sweep")
    sweep = []
    for ld in LINKS:
        lab = _labels(xy, ld)
        nv = int(lab.max() + 1)
        c = np.bincount(lab, minlength=nv).astype(np.float64)
        vv = segment_vessels(cells, dom, link_distance_um=ld)
        gg = vessel_geometry(vv, radius_mode="data")
        cov, cell_area = _subgrid_coverage(dom, gg)
        ta = float(cell_area[tissue].sum())
        pct = 100.0 * float((cov * cell_area)[tissue].sum()) / ta
        kkk = np.asarray(gg["ec_per_cross_section"], dtype=np.float64)
        row = {
            "link_um": ld,
            "n_vessels": nv,
            "mvd_per_mm2": nv / tissue_area_mm2,
            "singleton_fraction": float((c == 1).mean()),
            "giant_component_frac_of_ec": float(c.max() / n_endo),
            "mean_k": float(kkk.mean()),
            "area_fraction_pct": pct,
        }
        sweep.append(row)
        log(f"    link={ld:4.0f}um  MVD={row['mvd_per_mm2']:7.1f}  "
            f"giant={row['giant_component_frac_of_ec']*100:5.1f}%  "
            f"mean_k={row['mean_k']:5.2f}  AF={pct:7.3f}%")
    out["link_sweep"] = sweep

    # ---- blast radius: is the PRODUCTION delivery source touched by any of it?
    log(f"{name}: blast-radius check")
    prod = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    same = []
    for ld in (10.0, 20.0, 50.0):
        vv = segment_vessels(cells, dom, link_distance_um=ld)
        f = discrete_vascular_field(dom, vv)
        same.append({"link_um": ld,
                     "bit_identical_to_production": bool(np.array_equal(f, prod)),
                     "max_abs_diff": float(np.abs(f - prod).max())})
    sub = discrete_vascular_field(dom, v, geometry="subgrid")
    b = prod[tissue]
    s = sub[tissue]
    out["blast_radius"] = {
        "production_field_is_binary": bool(set(np.unique(prod).tolist()) <= {0.0, 1.0}),
        "production_matches_vascular_density_from_cells": bool(np.array_equal(
            prod, vascular_density_from_cells(dom, cells))),
        "link_distance_invariance": same,
        # what WOULD change if the corrected sub-grid geometry were promoted
        "binary_total_over_tissue": float(b.sum()),
        "subgrid_total_over_tissue": float(s.sum()),
        "subgrid_over_binary_total": float(s.sum() / max(b.sum(), 1e-12)),
        "binary_nonzero_frac_of_tissue": float((b > 0).mean()),
        "subgrid_nonzero_frac_of_tissue": float((s > 0).mean()),
        "pearson_r_binary_vs_subgrid": float(np.corrcoef(b, s)[0, 1]),
    }
    log(f"{name}: done")
    return out


def main() -> int:
    rows = []
    for name, sample in SECTIONS:
        try:
            rows.append(audit(name, sample))
        except Exception as exc:  # pragma: no cover
            log(f"{name}: FAILED {exc!r}")
    print("===JSON===")
    print(json.dumps(rows, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
