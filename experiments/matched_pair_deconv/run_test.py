"""Matched-pair, own-reference deconvolution test. Protocol: PREREGISTRATION.md.

Two stages so the TRUTH half is read after the prediction half is frozen:

    stage_a   prediction side only -- B189 scRNA reference, Visium NNLS.
    stage_b   truth side -- Xenium labels, the region-mismatch band, pseudo-spots.

Every scoring component is the repo's shipped default, imported and called
unmodified. ``tmesim/**`` is not edited by anything here.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
if REPO not in sys.path:
    sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import io_loaders as L  # noqa: E402
from tmesim.core.microenvironment import CellType  # noqa: E402
from tmesim.reconstruct import marker_panels  # noqa: E402
from tmesim.reconstruct.deconvolve import (  # noqa: E402
    build_signatures,
    deconvolve_spot,
    iter_spot_vectors,
    map_label,
)
from tmesim.reconstruct.real_xenium import (  # noqa: E402
    MATCHED_NULL_POOL_SIZE,
    TYPING_AGGREGATION,
    TYPING_THRESHOLD,
    phenotype_cells,
)

# --------------------------------------------------------------------------- #
# Constants fixed by PREREGISTRATION.md (sections 2.3, 3.1, 4). None is tuned.
# --------------------------------------------------------------------------- #
HEX_CELL_AREA_UM2 = 100.0 * 100.0 * np.sqrt(3.0) / 2.0   # 8660.254..., = §1.4
GRID_UM = 100.0                                          # occupancy grid pitch
GRID_CELL_AREA_UM2 = GRID_UM * GRID_UM
N_WINDOWS_TARGET = 2000
N_ATTEMPTS_CAP = 200_000
WINDOW_AREA_TOL = 0.10
WINDOW_MIN_CELLS = 200
RNG_SEED = 20261010
SPOT_RADIUS_UM = 27.5                                    # Visium 55 um diameter
PSEUDOSPOT_MIN_CELLS = 5
TAU_RETRACTED = 0.1778        # L1 of the 4.44 pp published in revalidation.md
TAU_CURRENT = 0.4884          # L1 of the 12.21 pp published in revalidation.md

#: Darmanis vocabulary, so the repo's own ``map_label`` performs the collapse.
CELLTYPE_TO_DARMANIS = {
    int(CellType.TUMOR): "Neoplastic",
    int(CellType.IMMUNE): "Immune cell",
    int(CellType.ENDOTHELIAL): "Vascular",
    int(CellType.STROMAL): "Astrocyte",
}
TYPE_ORDER = ["TUMOR", "STROMAL", "IMMUNE", "ENDOTHELIAL"]

PAIRS = [
    # (pair id, visium dir stem, xenium gsm, xenium stem, donor, arm)
    ("P1", "GSM9497624_TWBK-ST073021-2", "GSM9809293", "ST073021_2", "B189", "own"),
    ("P2", "GSM9497625_TWBK-ST073021-4", "GSM9809294", "ST073021_4", "B189", "own"),
    ("P3", "GSM9497628_TWKI-WU1220-3", "GSM9809296", "WU1220_3", "WU1220", "foreign"),
    ("P4", "GSM9497629_TWKI-WU1221-Core2", "GSM9809298", "WU1221_Core2", "WU1221", "foreign"),
    ("P5", "GSM9497630_TWKI-WU1221-Edge1", "GSM9809297", "WU1221_Edge", "WU1221", "foreign"),
]


# --------------------------------------------------------------------------- #
# Signature construction
# --------------------------------------------------------------------------- #
def dedup_first(symbols: List[str]) -> Tuple[List[str], np.ndarray]:
    """First-occurrence dedup, as the repo does. Returns (symbols, positions)."""
    seen, out, pos = set(), [], []
    for i, s in enumerate(symbols):
        if s not in seen:
            seen.add(s)
            out.append(s)
            pos.append(i)
    return out, np.asarray(pos, dtype=np.int64)


def signatures_sparse(
    counts_cells_x_genes: sparse.csr_matrix,
    ref_genes: List[str],
    labels: List[str],
    shared_genes: List[str],
) -> Tuple[np.ndarray, List[str]]:
    """``build_signatures`` arithmetic, done sparsely.

    Identical formula::

        sig[g, k] = mean_{cell c in class k} ( X[g, c] / sum_{g' in found} X[g', c] )

    where ``found`` is the subset of ``shared_genes`` present in the
    (first-occurrence-deduped) reference index. Rows for shared genes absent
    from the reference stay zero, exactly as the repo's version leaves them.
    Verified against the repo's own function in :func:`check_signature_equivalence`.
    """
    ref_dedup, ref_pos = dedup_first(list(ref_genes))
    pos_of = {g: int(p) for g, p in zip(ref_dedup, ref_pos)}
    found = [g for g in shared_genes if g in pos_of]
    row_of = {g: i for i, g in enumerate(shared_genes)}
    found_rows = np.asarray([row_of[g] for g in found], dtype=np.int64)
    found_cols = np.asarray([pos_of[g] for g in found], dtype=np.int64)

    X = counts_cells_x_genes[:, found_cols].tocsr().astype(np.float64)  # cells x found
    tot = np.asarray(X.sum(axis=1)).ravel()
    inv = np.zeros_like(tot)
    nz = tot > 0
    inv[nz] = 1.0 / tot[nz]
    Xn = sparse.diags(inv) @ X

    mapped = [map_label(s) for s in labels]
    out_names, cols = [], []
    for ct in (CellType.TUMOR, CellType.STROMAL, CellType.IMMUNE, CellType.ENDOTHELIAL):
        sel = np.asarray([m == ct for m in mapped], dtype=bool)
        if not sel.any():
            continue
        out_names.append(ct.name)
        cols.append(np.asarray(Xn[sel].mean(axis=0)).ravel())
    sig = np.zeros((len(shared_genes), len(cols)), dtype=np.float64)
    for k, c in enumerate(cols):
        sig[found_rows, k] = c
    return sig, out_names


def check_signature_equivalence(
    counts_cells_x_genes, ref_genes, labels, shared_genes, n_cells=1200, seed=7,
) -> float:
    """Max abs difference between the sparse build and the repo's own.

    Run on a subsample so the dense genes x cells frame the repo's function
    needs actually fits; the arithmetic being checked does not depend on which
    cells are in it.
    """
    rng = np.random.default_rng(seed)
    n = counts_cells_x_genes.shape[0]
    idx = np.sort(rng.choice(n, size=min(n_cells, n), replace=False))
    sub = counts_cells_x_genes[idx]
    sub_labels = [labels[i] for i in idx]

    mine, mine_names = signatures_sparse(sub, ref_genes, sub_labels, shared_genes)

    ref_dedup, ref_pos = dedup_first(list(ref_genes))
    keep = np.asarray([g in set(shared_genes) for g in ref_dedup], dtype=bool)
    rows = ref_pos[keep]
    syms = [g for g, k in zip(ref_dedup, keep) if k]
    dense = np.asarray(sub[:, rows].todense()).T          # genes x cells
    bcs = [f"c{i}" for i in range(dense.shape[1])]
    expr_df = pd.DataFrame(dense, index=syms, columns=bcs)
    labels_df = pd.DataFrame({"cell_type": sub_labels}, index=bcs)
    theirs, their_names = build_signatures(
        shared_genes, reference=(expr_df, labels_df, syms)
    )
    if list(mine_names) != list(their_names):
        raise AssertionError(f"type order differs: {mine_names} vs {their_names}")
    return float(np.max(np.abs(mine - theirs)))


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def pca_align(xy: np.ndarray) -> np.ndarray:
    """Rotate so PC1 is the x axis. Same device as the 2026-10-09 footprints."""
    c = xy - xy.mean(axis=0, keepdims=True)
    _u, _s, vt = np.linalg.svd(c, full_matrices=False)
    return c @ vt.T


def occupancy(xy: np.ndarray, grid: float = GRID_UM):
    """Grid-cell indices per point, plus the integer grid shape and origin."""
    o = np.floor(xy.min(axis=0) / grid) * grid
    ij = np.floor((xy - o) / grid).astype(np.int64)
    ny = int(ij[:, 1].max()) + 1
    nx = int(ij[:, 0].max()) + 1
    return ij, (ny, nx), o


def summed_area(a: np.ndarray) -> np.ndarray:
    """Integral image with a zero first row/column (O(1) rectangle sums)."""
    s = np.zeros((a.shape[0] + 1, a.shape[1] + 1), dtype=np.float64)
    s[1:, 1:] = np.cumsum(np.cumsum(np.asarray(a, dtype=np.float64), axis=0), axis=1)
    return s


def rect_sum(s: np.ndarray, i0: int, i1: int, j0: int, j1: int) -> float:
    return float(s[i1, j1] - s[i0, j1] - s[i1, j0] + s[i0, j0])


def rect_sum_vec(s, i0, i1, j0, j1):
    return s[i1, j1] - s[i0, j1] - s[i1, j0] + s[i0, j0]


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def l1(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.abs(np.asarray(a) - np.asarray(b)).sum())


def nontumour(v: np.ndarray, names: List[str]) -> np.ndarray:
    """Renormalised (STROMAL, IMMUNE, ENDOTHELIAL) subvector."""
    idx = [names.index(n) for n in ("STROMAL", "IMMUNE", "ENDOTHELIAL")]
    sub = np.asarray(v, dtype=np.float64)[idx]
    t = sub.sum()
    return sub / t if t > 0 else sub


def align_to_order(vec: np.ndarray, names: List[str]) -> np.ndarray:
    """Reorder a composition vector onto TYPE_ORDER, zero-filling absent types."""
    out = np.zeros(len(TYPE_ORDER))
    for n, v in zip(names, np.asarray(vec, dtype=np.float64)):
        out[TYPE_ORDER.index(n)] = v
    return out


# --------------------------------------------------------------------------- #
# STAGE A -- prediction side
# --------------------------------------------------------------------------- #
def stage_a(args) -> None:
    panel = L.xenium_panel_genes(args.panel_zarr)
    print(f"[panel] {len(panel)} features from {os.path.basename(args.panel_zarr)}")

    print("[scrna] reading B189 ...")
    sc, sc_genes, _sc_bcs = L.read_10x(args.scrna)
    print(f"[scrna] {sc.shape[0]} cells x {sc.shape[1]} genes")

    # --- labels: same classifier, same markers, same 542-gene space as Xenium
    sc_dedup, sc_pos = dedup_first(sc_genes)
    pos_of = {g: int(p) for g, p in zip(sc_dedup, sc_pos)}
    panel_present = [g for g in panel if g in pos_of]
    panel_cols = np.asarray([pos_of[g] for g in panel_present], dtype=np.int64)
    sc_panel = sc[:, panel_cols].tocsr()
    mmap, pname, cov = marker_panels.select_marker_map(panel_present)
    print(f"[labels] marker map {pname!r} coverage {cov}")
    ct, used_map, missing = phenotype_cells(
        sc_panel, np.asarray(panel_present, dtype=str),
        threshold=TYPING_THRESHOLD, marker_map=mmap,
        aggregation=TYPING_AGGREGATION, pool_size=MATCHED_NULL_POOL_SIZE,
    )
    labels = [CELLTYPE_TO_DARMANIS[int(v)] for v in ct]
    counts = {n: int(np.sum(ct == int(getattr(CellType, n)))) for n in TYPE_ORDER}
    print(f"[labels] scRNA class counts {counts}  missing markers {missing}")

    out: Dict = {
        "panel_n_features": len(panel),
        "panel_genes_in_scrna": len(panel_present),
        "marker_map_panel": pname,
        "marker_map_coverage": {k: float(v) for k, v in cov.items()},
        "marker_map_used": {k.name: v for k, v in used_map.items()},
        "markers_missing_from_panel": list(missing),
        "scrna_n_cells": int(sc.shape[0]),
        "scrna_n_genes": int(sc.shape[1]),
        "scrna_class_counts": counts,
        "typing": {"threshold": TYPING_THRESHOLD, "aggregation": TYPING_AGGREGATION,
                   "pool_size": MATCHED_NULL_POOL_SIZE},
        "sections": {},
    }

    # --- Darmanis reference symbols, for the C6 contrast
    from tmesim.data.real.scref_gbm import reference_gene_symbols
    darmanis_syms = set(reference_gene_symbols())

    equiv_done = False
    for pid, vstem, _xg, _xs, donor, arm in PAIRS:
        vdir_m = os.path.join(args.visium, f"{vstem}_filtered_feature_bc_matrix")
        vdir_s = os.path.join(args.visium, f"{vstem}_spatial")
        print(f"[{pid}] {vstem}")
        vx, vgenes, vbcs = L.read_10x(vdir_m)
        coords_map, um_px = L.read_visium_spatial(vdir_s)
        keep = np.asarray([b in coords_map for b in vbcs], dtype=bool)
        vx = vx[keep].tocsr()
        bcs = [b for b, k in zip(vbcs, keep) if k]
        coords = np.array([coords_map[b] for b in bcs], dtype=np.float64)
        print(f"[{pid}] {vx.shape[0]} in-tissue spots x {vx.shape[1]} genes, "
              f"{um_px:.5f} um/px")

        rec: Dict = {
            "pair": pid, "donor": donor, "arm": arm, "visium_stem": vstem,
            "n_spots_matrix": int(len(vbcs)),
            "n_spots_in_tissue_used": int(vx.shape[0]),
            "n_genes_visium": int(vx.shape[1]),
            "um_per_fullres_px": um_px,
            "visium_area_mm2": float(vx.shape[0] * HEX_CELL_AREA_UM2 / 1e6),
        }
        al = pca_align(coords)
        ext = al.max(axis=0) - al.min(axis=0)
        rec["visium_extent_um_pca"] = [float(ext[0]), float(ext[1])]
        rec["visium_aspect_pca"] = float(max(ext) / max(min(ext), 1e-9))

        for refname, refsyms in (("own_b189", set(sc_dedup)),
                                 ("darmanis", darmanis_syms)):
            seen, shared_genes, shared_cols = set(), [], []
            for j, g in enumerate(vgenes):
                if g in refsyms and g not in seen:
                    seen.add(g)
                    shared_genes.append(g)
                    shared_cols.append(j)
            if refname == "own_b189":
                sig, names = signatures_sparse(sc, sc_genes, labels, shared_genes)
                if not equiv_done:
                    d = check_signature_equivalence(sc, sc_genes, labels, shared_genes)
                    out["signature_equivalence_max_abs_diff"] = d
                    print(f"[check] sparse-vs-repo signature max|diff| = {d:.3e}")
                    equiv_done = True
            else:
                sig, names = build_signatures(shared_genes)

            sub = vx[:, np.asarray(shared_cols, dtype=int)].tocsr()
            fr = np.zeros((sub.shape[0], sig.shape[1]), dtype=np.float64)
            for i, spot in enumerate(iter_spot_vectors(sub, len(shared_genes))):
                fr[i] = deconvolve_spot(spot, sig)
            comp = align_to_order(fr.mean(axis=0), names)
            rec[refname] = {
                "n_shared_genes": len(shared_genes),
                "type_names": names,
                "composition_pct": {n: float(100.0 * c) for n, c in zip(TYPE_ORDER, comp)},
                "composition": [float(c) for c in comp],
            }
            print(f"[{pid}/{refname}] genes={len(shared_genes):>6}  "
                  + "  ".join(f"{n}={100*c:5.2f}%" for n, c in zip(TYPE_ORDER, comp)))

        np.save(os.path.join(args.out, f"{pid}_visium_coords_um.npy"), coords)
        out["sections"][pid] = rec

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "stage_a.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    np.save(os.path.join(args.out, "scrna_labels.npy"), np.asarray(labels, dtype=object),
            allow_pickle=True)
    print("[stage_a] written")


# --------------------------------------------------------------------------- #
# STAGE B -- truth side, region band, pseudo-spots
# --------------------------------------------------------------------------- #
def region_band(xy_al, ct, a_vis_um2, aspect, rng):
    """Sub-window ensemble per PREREGISTRATION.md §2.3.

    A window is an axis-aligned rectangle in the PCA-aligned Xenium frame with
    the Visium extent's aspect ratio, grown about its centre until the OCCUPIED
    tissue area inside it matches ``a_vis_um2``; accepted when the match is
    within +/-10 % and the window holds >= 200 cells.
    """
    ij, (ny, nx), _o = occupancy(xy_al)
    occ = np.zeros((ny, nx), dtype=np.float64)
    occ[ij[:, 1], ij[:, 0]] = 1.0
    per_type = {}
    for n in TYPE_ORDER:
        g = np.zeros((ny, nx), dtype=np.float64)
        sel = ct == int(getattr(CellType, n))
        np.add.at(g, (ij[sel, 1], ij[sel, 0]), 1.0)
        per_type[n] = summed_area(g)
    s_occ = summed_area(occ)
    n_all = np.zeros((ny, nx), dtype=np.float64)
    np.add.at(n_all, (ij[:, 1], ij[:, 0]), 1.0)
    s_all = summed_area(n_all)

    a_xen = float(occ.sum()) * GRID_CELL_AREA_UM2
    target_cells = a_vis_um2 / GRID_CELL_AREA_UM2          # occupied grid cells
    # base half-extents in grid units for area a_vis with the Visium aspect
    base_h = np.sqrt(max(target_cells, 1.0) / max(aspect, 1e-9))
    base_w = base_h * max(aspect, 1e-9)

    comps, l1_wf, attempts, accepted = [], [], 0, 0
    whole = np.array([rect_sum(per_type[n], 0, ny, 0, nx) for n in TYPE_ORDER])
    whole = whole / max(whole.sum(), 1e-12)
    scales = np.linspace(0.6, 4.0, 60)
    while accepted < N_WINDOWS_TARGET and attempts < N_ATTEMPTS_CAP:
        attempts += 1
        cx = rng.uniform(0, nx)
        cy = rng.uniform(0, ny)
        got = None
        for s in scales:
            j0 = int(max(0, np.floor(cx - 0.5 * base_w * s)))
            j1 = int(min(nx, np.ceil(cx + 0.5 * base_w * s)))
            i0 = int(max(0, np.floor(cy - 0.5 * base_h * s)))
            i1 = int(min(ny, np.ceil(cy + 0.5 * base_h * s)))
            if j1 <= j0 or i1 <= i0:
                continue
            occ_cells = rect_sum(s_occ, i0, i1, j0, j1)
            if occ_cells >= target_cells * (1.0 - WINDOW_AREA_TOL):
                if occ_cells <= target_cells * (1.0 + WINDOW_AREA_TOL):
                    got = (i0, i1, j0, j1)
                break
        if got is None:
            continue
        i0, i1, j0, j1 = got
        if rect_sum(s_all, i0, i1, j0, j1) < WINDOW_MIN_CELLS:
            continue
        v = np.array([rect_sum(per_type[n], i0, i1, j0, j1) for n in TYPE_ORDER])
        t = v.sum()
        if t <= 0:
            continue
        v = v / t
        comps.append(v)
        l1_wf.append(l1(v, whole))
        accepted += 1
    return {
        "whole_field": whole, "windows": np.asarray(comps),
        "l1_window_vs_whole": np.asarray(l1_wf),
        "attempts": attempts, "accepted": accepted,
        "xenium_area_mm2": a_xen / 1e6,
        "grid_shape": (ny, nx),
    }


def pseudospots(xy_al, ct, xcounts, panel_genes, sig, names):
    """C7: Visium-geometry pseudo-spots on the Xenium cells, same signature."""
    pitch = GRID_UM
    x0, y0 = xy_al.min(axis=0)
    x1, y1 = xy_al.max(axis=0)
    rows = int(np.ceil((y1 - y0) / (pitch * np.sqrt(3.0) / 2.0))) + 1
    cols = int(np.ceil((x1 - x0) / pitch)) + 2
    centres = []
    for r in range(rows):
        yy = y0 + r * pitch * np.sqrt(3.0) / 2.0
        off = 0.0 if r % 2 == 0 else pitch / 2.0
        for c in range(cols):
            centres.append((x0 + off + c * pitch, yy))
    centres = np.asarray(centres)

    from scipy.spatial import cKDTree
    tree = cKDTree(xy_al)
    members = tree.query_ball_point(centres, r=SPOT_RADIUS_UM)

    truth, pred = [], []
    for mem in members:
        if len(mem) < PSEUDOSPOT_MIN_CELLS:
            continue
        mem = np.asarray(mem, dtype=np.int64)
        tv = np.array([np.sum(ct[mem] == int(getattr(CellType, n))) for n in TYPE_ORDER],
                      dtype=np.float64)
        if tv.sum() <= 0:
            continue
        agg = np.asarray(xcounts[mem].sum(axis=0), dtype=np.float64).ravel()
        pv = align_to_order(deconvolve_spot(agg, sig), names)
        truth.append(tv / tv.sum())
        pred.append(pv)
    truth = np.asarray(truth)
    pred = np.asarray(pred)
    return truth, pred


def stage_b(args) -> None:
    with open(os.path.join(args.out, "stage_a.json")) as fh:
        A = json.load(fh)
    panel = L.xenium_panel_genes(args.panel_zarr)
    sc, sc_genes, _ = L.read_10x(args.scrna)
    labels = list(np.load(os.path.join(args.out, "scrna_labels.npy"),
                          allow_pickle=True))
    rng = np.random.default_rng(RNG_SEED)
    results: Dict = {"pairs": {}, "prereg_sha256": args.prereg_sha}

    for pid, vstem, xg, xs, donor, arm in PAIRS:
        zpath = os.path.join(args.xenium, f"{xg}_{xs}_cell_feature_matrix.zarr.zip")
        cpath = os.path.join(args.xenium, f"{xg}_{xs}_cells.parquet.gz")
        if not (os.path.exists(zpath) and os.path.exists(cpath)):
            print(f"[{pid}] SKIP -- xenium files absent")
            continue
        print(f"[{pid}] xenium {xs}")
        xc, xgenes = L.read_xenium_matrix(zpath)
        xy, _ids = L.read_xenium_cells(cpath)
        if xy.shape[0] != xc.shape[0]:
            raise ValueError(f"{pid}: {xy.shape[0]} centroids vs {xc.shape[0]} rows")
        if list(xgenes) != list(panel):
            raise ValueError(f"{pid}: panel differs from {args.panel_zarr}")
        mmap, pname, _cov = marker_panels.select_marker_map(xgenes)
        ct, _um, _ms = phenotype_cells(
            xc, np.asarray(xgenes, dtype=str), threshold=TYPING_THRESHOLD,
            marker_map=mmap, aggregation=TYPING_AGGREGATION,
            pool_size=MATCHED_NULL_POOL_SIZE,
        )
        xy_al = pca_align(xy)
        sec = A["sections"][pid]
        band = region_band(xy_al, ct, sec["visium_area_mm2"] * 1e6,
                           sec["visium_aspect_pca"], rng)

        rec: Dict = {
            "donor": donor, "arm": arm, "xenium_stem": xs, "xenium_gsm": xg,
            "n_cells": int(xc.shape[0]), "marker_map_panel": pname,
            "xenium_area_mm2": band["xenium_area_mm2"],
            "visium_area_mm2": sec["visium_area_mm2"],
            "area_ratio_xen_over_vis": band["xenium_area_mm2"] / sec["visium_area_mm2"],
            "xenium_cell_density_per_mm2": xc.shape[0] / band["xenium_area_mm2"],
            "windows_accepted": band["accepted"], "windows_attempted": band["attempts"],
            "xenium_whole_field_pct": {n: float(100.0 * v)
                                       for n, v in zip(TYPE_ORDER, band["whole_field"])},
        }
        W = band["windows"]
        if W.size:
            rec["l1_wf_50"] = float(np.percentile(band["l1_window_vs_whole"], 50))
            rec["l1_wf_95"] = float(np.percentile(band["l1_window_vs_whole"], 95))
            Wnt = np.asarray([nontumour(w, TYPE_ORDER) for w in W])
            wnt = nontumour(band["whole_field"], TYPE_ORDER)
            l1nt = np.abs(Wnt - wnt[None, :]).sum(axis=1)
            rec["l1nt_wf_50"] = float(np.percentile(l1nt, 50))
            rec["l1nt_wf_95"] = float(np.percentile(l1nt, 95))

        for refname in ("own_b189", "darmanis"):
            if refname not in sec:
                continue
            p = np.asarray(sec[refname]["composition"], dtype=np.float64)
            r: Dict = {
                "visium_composition_pct": sec[refname]["composition_pct"],
                "l1_full": l1(p, band["whole_field"]),
                "l1nt_full": l1(nontumour(p, TYPE_ORDER),
                                nontumour(band["whole_field"], TYPE_ORDER)),
            }
            if W.size:
                d = np.abs(W - p[None, :]).sum(axis=1)
                k = int(np.argmin(d))
                r["l1_best"] = float(d[k])
                r["best_window_pct"] = {n: float(100.0 * v)
                                        for n, v in zip(TYPE_ORDER, W[k])}
                Wnt = np.asarray([nontumour(w, TYPE_ORDER) for w in W])
                dnt = np.abs(Wnt - nontumour(p, TYPE_ORDER)[None, :]).sum(axis=1)
                r["l1nt_best"] = float(np.min(dnt))
                r["l1_full_vs_band95"] = ("inside" if r["l1_full"] <= rec["l1_wf_95"]
                                          else "outside")
            r["mean_abs_diff_pp_full"] = 100.0 * r["l1_full"] / 4.0
            if "l1_best" in r:
                r["mean_abs_diff_pp_best"] = 100.0 * r["l1_best"] / 4.0
                r["C1_fail_worse_than_current"] = bool(r["l1_best"] > TAU_CURRENT)
                r["C2_meets_retracted_era"] = bool(r["l1_best"] <= TAU_RETRACTED)
            rho = spearmanr(p, band["whole_field"]).statistic
            r["spearman_rho_full"] = float(rho) if rho == rho else None
            rec[refname] = r

        # C7, own reference only and only where the donor's own scRNA exists
        if arm == "own":
            sc_dedup, sc_pos = dedup_first(sc_genes)
            pos_of = {g: int(p) for g, p in zip(sc_dedup, sc_pos)}
            shared = [g for g in xgenes if g in pos_of]
            sig7, names7 = signatures_sparse(sc, sc_genes, labels, shared)
            cols7 = np.asarray([list(xgenes).index(g) for g in shared], dtype=np.int64)
            truth, pred = pseudospots(xy_al, ct, xc[:, cols7].tocsr(), shared,
                                      sig7, names7)
            if truth.size:
                per = np.abs(truth - pred)
                rec["C7_pseudospots"] = {
                    "n_pseudospots": int(truth.shape[0]),
                    "n_genes": len(shared),
                    "mean_l1_per_spot": float(per.sum(axis=1).mean()),
                    "median_l1_per_spot": float(np.median(per.sum(axis=1))),
                    "pooled_truth_pct": {n: float(100.0 * v) for n, v
                                         in zip(TYPE_ORDER, truth.mean(axis=0))},
                    "pooled_pred_pct": {n: float(100.0 * v) for n, v
                                        in zip(TYPE_ORDER, pred.mean(axis=0))},
                    "pooled_l1": l1(truth.mean(axis=0), pred.mean(axis=0)),
                    "pooled_mean_abs_diff_pp": 100.0 * l1(truth.mean(axis=0),
                                                          pred.mean(axis=0)) / 4.0,
                }
        results["pairs"][pid] = rec
        print(json.dumps(rec, indent=2, sort_keys=True)[:2000])

    with open(os.path.join(args.out, "stage_b.json"), "w") as fh:
        json.dump(results, fh, indent=2, sort_keys=True)
    print("[stage_b] written")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["stage_a", "stage_b"])
    ap.add_argument("--scrna", default="/tmp/vmp/ex_scrna")
    ap.add_argument("--visium", default="/tmp/vmp/ex_visium")
    ap.add_argument("--xenium", default="/tmp/vmp/dl_xenium")
    ap.add_argument("--panel-zarr", dest="panel_zarr",
                    default="/tmp/vmp/dl_xenium/"
                            "GSM9809297_WU1221_Edge_cell_feature_matrix.zarr.zip")
    ap.add_argument("--out", default="/tmp/vmp/out")
    ap.add_argument("--prereg-sha", dest="prereg_sha", default="")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    (stage_a if args.stage == "stage_a" else stage_b)(args)


if __name__ == "__main__":
    main()
