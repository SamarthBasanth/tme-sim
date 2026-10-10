"""Three post-hoc DIAGNOSTICS. None changes a pre-registered number.

D1  Why is only part of the 542-gene panel found in the scRNA reference?
D2  What is in the own-reference ENDOTHELIAL class, given cos(IMMUNE,ENDO)=0.978?
D3  Is the region-mismatch band narrowed by boundary-clipped windows?
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import io_loaders as L  # noqa: E402
from tmesim.core.microenvironment import CellType  # noqa: E402
from run_test import (  # noqa: E402
    GRID_CELL_AREA_UM2, GRID_UM, TYPE_ORDER, WINDOW_AREA_TOL, WINDOW_MIN_CELLS,
    dedup_first, l1, occupancy, pca_align, rect_sum, summed_area,
)

OUT = "/tmp/vmp/out"
A = json.load(open(os.path.join(OUT, "stage_a.json")))

print("=" * 78)
print("D1 -- panel genes not found in the B189 scRNA symbol index")
panel = L.xenium_panel_genes(
    "/tmp/vmp/dl_xenium/GSM9809297_WU1221_Edge_cell_feature_matrix.zarr.zip")
sc, sc_genes, _ = L.read_10x("/tmp/vmp/ex_scrna")
sc_dedup, sc_pos = dedup_first(sc_genes)
pos_of = {g: int(p) for g, p in zip(sc_dedup, sc_pos)}
present = [g for g in panel if g in pos_of]
absent = [g for g in panel if g not in pos_of]
print(f"  panel {len(panel)}   present {len(present)}   absent {len(absent)}")
print(f"  absent (first 60): {absent[:60]}")
print(f"  stage_a recorded panel_genes_in_scrna = {A['panel_genes_in_scrna']}")
print("  brain-panel markers present: "
      + str({m: (m in pos_of) for m in
             ["EGFR", "SOX2", "OLIG2", "MKI67", "PTPRC", "CD68", "AIF1", "CD163",
              "PECAM1", "MOG", "AQP4", "GJA1"]}))

print("=" * 78)
print("D2 -- composition of the own-reference ENDOTHELIAL class")
labels = list(np.load(os.path.join(OUT, "scrna_labels.npy"), allow_pickle=True))
lab = np.asarray(labels, dtype=object)
groups = {"Neoplastic": "TUMOR", "Astrocyte": "STROMAL",
          "Immune cell": "IMMUNE", "Vascular": "ENDOTHELIAL"}
markers = {"ENDO": ["PECAM1"], "IMMUNE": ["PTPRC", "CD68", "AIF1", "CD163"],
           "TUMOR": ["EGFR", "SOX2", "OLIG2", "MKI67"],
           "STROMAL": ["MOG", "AQP4", "GJA1"],
           "MURAL": ["RGS5", "PDGFRB", "NOTCH3", "ACTA2"]}
print("  detection rate (% of cells with count > 0) by assigned class")
hdr = "  {:<12}".format("class") + "".join(f"{g:>10s}" for grp in markers.values()
                                           for g in grp)
print(hdr)
for dk, name in groups.items():
    sel = lab == dk
    if not sel.any():
        continue
    row = f"  {name:<12}"
    for grp in markers.values():
        for g in grp:
            if g not in pos_of:
                row += f"{'n/a':>10s}"
                continue
            col = np.asarray(sc[sel][:, pos_of[g]].todense()).ravel()
            row += f"{100.0 * np.mean(col > 0):9.1f}%"
    print(row + f"   (n={int(sel.sum()):,})")

print("=" * 78)
print("D3 -- region band with boundary-clipped windows EXCLUDED")
panel_list = list(panel)
for pid, xg, xs in (("P1", "GSM9809293", "ST073021_2"),
                    ("P2", "GSM9809294", "ST073021_4"),
                    ("P3", "GSM9809296", "WU1220_3"),
                    ("P4", "GSM9809298", "WU1221_Core2"),
                    ("P5", "GSM9809297", "WU1221_Edge")):
    from tmesim.reconstruct import marker_panels
    from tmesim.reconstruct.real_xenium import (
        MATCHED_NULL_POOL_SIZE, TYPING_AGGREGATION, TYPING_THRESHOLD, phenotype_cells,
    )
    z = f"/tmp/vmp/dl_xenium/{xg}_{xs}_cell_feature_matrix.zarr.zip"
    c = f"/tmp/vmp/dl_xenium/{xg}_{xs}_cells.parquet.gz"
    xc, xgenes = L.read_xenium_matrix(z)
    xy, _ = L.read_xenium_cells(c)
    mmap, _pn, _cv = marker_panels.select_marker_map(xgenes)
    ct, _u, _m = phenotype_cells(xc, np.asarray(xgenes, dtype=str),
                                 threshold=TYPING_THRESHOLD, marker_map=mmap,
                                 aggregation=TYPING_AGGREGATION,
                                 pool_size=MATCHED_NULL_POOL_SIZE)
    xy_al = pca_align(xy)
    sec = A["sections"][pid]
    a_vis = sec["visium_area_mm2"] * 1e6
    aspect = sec["visium_aspect_pca"]

    ij, (ny, nx), _o = occupancy(xy_al)
    occ = np.zeros((ny, nx)); occ[ij[:, 1], ij[:, 0]] = 1.0
    s_occ = summed_area(occ)
    per = {}
    for n in TYPE_ORDER:
        g = np.zeros((ny, nx)); sel = ct == int(getattr(CellType, n))
        np.add.at(g, (ij[sel, 1], ij[sel, 0]), 1.0); per[n] = summed_area(g)
    n_all = np.zeros((ny, nx)); np.add.at(n_all, (ij[:, 1], ij[:, 0]), 1.0)
    s_all = summed_area(n_all)
    whole = np.array([rect_sum(per[n], 0, ny, 0, nx) for n in TYPE_ORDER])
    whole = whole / whole.sum()

    target = a_vis / GRID_CELL_AREA_UM2
    bh = np.sqrt(max(target, 1.0) / max(aspect, 1e-9)); bw = bh * max(aspect, 1e-9)
    rng = np.random.default_rng(20261010)
    scales = np.linspace(0.6, 4.0, 60)
    d_all, d_unclipped = [], []
    att = 0
    while att < 200_000 and len(d_all) < 2000:
        att += 1
        cx, cy = rng.uniform(0, nx), rng.uniform(0, ny)
        got = None
        for s in scales:
            j0r, j1r = cx - .5 * bw * s, cx + .5 * bw * s
            i0r, i1r = cy - .5 * bh * s, cy + .5 * bh * s
            j0, j1 = int(max(0, np.floor(j0r))), int(min(nx, np.ceil(j1r)))
            i0, i1 = int(max(0, np.floor(i0r))), int(min(ny, np.ceil(i1r)))
            if j1 <= j0 or i1 <= i0:
                continue
            oc = rect_sum(s_occ, i0, i1, j0, j1)
            if oc >= target * (1 - WINDOW_AREA_TOL):
                if oc <= target * (1 + WINDOW_AREA_TOL):
                    clipped = (j0r < 0) or (i0r < 0) or (j1r > nx) or (i1r > ny)
                    got = (i0, i1, j0, j1, clipped)
                break
        if got is None:
            continue
        i0, i1, j0, j1, clipped = got
        if rect_sum(s_all, i0, i1, j0, j1) < WINDOW_MIN_CELLS:
            continue
        v = np.array([rect_sum(per[n], i0, i1, j0, j1) for n in TYPE_ORDER])
        if v.sum() <= 0:
            continue
        d = l1(v / v.sum(), whole)
        d_all.append(d)
        if not clipped:
            d_unclipped.append(d)
    d_all = np.asarray(d_all); d_unc = np.asarray(d_unclipped)
    msg = (f"  {pid:>3s}  all n={len(d_all):>4d} p50={np.percentile(d_all,50):.4f} "
           f"p95={np.percentile(d_all,95):.4f} max={d_all.max():.4f}")
    if d_unc.size >= 30:
        msg += (f"   | unclipped n={len(d_unc):>4d} "
                f"p50={np.percentile(d_unc,50):.4f} "
                f"p95={np.percentile(d_unc,95):.4f} max={d_unc.max():.4f}")
    else:
        msg += f"   | unclipped n={len(d_unc)} (too few to quote)"
    print(msg)
