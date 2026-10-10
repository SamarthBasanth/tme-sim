"""Minimal readers for the GEO files this test uses. No repo file is modified.

Everything here is a reader. No scoring, no thresholds, no choices that could
move a result; those all live in PREREGISTRATION.md and in run_test.py.
"""
from __future__ import annotations

import gzip
import io
import json
import os
import zipfile
from typing import Dict, List, Tuple

import numpy as np
from scipy import sparse
from scipy.io import mmread


# --------------------------------------------------------------------------- #
# 10x Market-Matrix triplets (scRNA and Visium both ship this layout)
# --------------------------------------------------------------------------- #
def _find_triplet(root: str) -> str:
    """Directory inside ``root`` holding matrix.mtx.gz / features / barcodes."""
    for dirpath, _dirnames, filenames in os.walk(root):
        if "matrix.mtx.gz" in filenames:
            return dirpath
    raise FileNotFoundError(f"no matrix.mtx.gz under {root}")


def read_10x(root: str) -> Tuple[sparse.csr_matrix, List[str], List[str]]:
    """Return ``(counts_cells_x_genes, gene_symbols, barcodes)``.

    The .mtx on disk is genes x cells; it is transposed once here so every
    caller sees the (observations, features) orientation that
    ``phenotype_cells`` and the deconvolution both want.
    """
    d = _find_triplet(root)
    with gzip.open(os.path.join(d, "matrix.mtx.gz"), "rb") as fh:
        m = mmread(io.BytesIO(fh.read()))
    m = sparse.csr_matrix(m)                      # genes x cells
    syms: List[str] = []
    with gzip.open(os.path.join(d, "features.tsv.gz"), "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            syms.append(parts[1] if len(parts) > 1 else parts[0])
    with gzip.open(os.path.join(d, "barcodes.tsv.gz"), "rt") as fh:
        bcs = [ln.strip() for ln in fh if ln.strip()]
    if m.shape != (len(syms), len(bcs)):
        raise ValueError(f"{d}: mtx {m.shape} vs {len(syms)} genes {len(bcs)} cells")
    return sparse.csr_matrix(m.T), syms, bcs


# --------------------------------------------------------------------------- #
# Visium spatial
# --------------------------------------------------------------------------- #
def read_visium_spatial(root: str) -> Tuple[Dict[str, np.ndarray], float]:
    """``({barcode: [x_um, y_um]}, um_per_fullres_px)`` for IN-TISSUE spots.

    Spot pitch is 100 um and spot diameter 55 um, so the scale factor is
    ``55 / spot_diameter_fullres`` -- the same conversion the 2026-10-09
    footprint measurement used.
    """
    pos_path = None
    scale_path = None
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn in ("tissue_positions_list.csv", "tissue_positions.csv"):
                pos_path = os.path.join(dirpath, fn)
            if fn == "scalefactors_json.json":
                scale_path = os.path.join(dirpath, fn)
    if pos_path is None or scale_path is None:
        raise FileNotFoundError(f"spatial files missing under {root}")
    with open(scale_path) as fh:
        sf = json.load(fh)
    um_per_px = 55.0 / float(sf["spot_diameter_fullres"])

    coords: Dict[str, np.ndarray] = {}
    with open(pos_path) as fh:
        first = fh.readline()
        rows = []
        if not first.lower().startswith("barcode"):
            rows.append(first)
        rows.extend(fh.readlines())
    for line in rows:
        p = line.strip().split(",")
        if len(p) < 6:
            continue
        if int(p[1]) != 1:                        # in_tissue == 1 only
            continue
        row_px, col_px = float(p[4]), float(p[5])
        coords[p[0]] = np.array([col_px * um_per_px, row_px * um_per_px])
    return coords, um_per_px


# --------------------------------------------------------------------------- #
# Xenium: cell_feature_matrix.zarr.zip and cells.parquet.gz
# --------------------------------------------------------------------------- #
#: Feature types the repo's own loader keeps when ``genes_only=True``
#: (``tmesim/data/real/singlecell_spatial.py``: "restrict to real Gene
#: Expression targets, dropping negative-control and unassigned codeword
#: features"). The GSE335275 panel declares 362 ``gene`` features inside 542
#: feature columns; the other 180 are 111 unassigned codewords, 41 negative
#: control codewords, 20 negative control probes, 7 deprecated codewords and
#: 1 ``aggregate_gene`` pseudo-feature (``Total transcripts``).
GENE_FEATURE_TYPES = ("gene",)


def read_xenium_matrix(
    zip_path: str, genes_only: bool = True,
) -> Tuple[sparse.csr_matrix, List[str]]:
    """``(counts_cells_x_genes, feature_keys)`` from a Xenium zarr.zip.

    The archive stores the matrix twice -- a CSR over features and a ``csc/``
    group over cells. ``csc/indptr`` has ``n_cells + 1`` entries, so that group
    is the cells-major view and is the one read here.

    ``genes_only`` matches the default of the repo's own Xenium loader: real
    gene targets only, control and codeword features dropped.
    """
    import zarr

    store = zarr.ZipStore(zip_path, mode="r")
    g = zarr.open(store, mode="r")
    att = dict(g["cell_features"].attrs)
    keys = [str(k) for k in att["feature_keys"]]
    n_cells = int(att["number_cells"])
    n_feat = int(att["number_features"])
    cf = g["cell_features"]["csc"]
    data = np.asarray(cf["data"][:], dtype=np.float64)
    indices = np.asarray(cf["indices"][:], dtype=np.int64)
    indptr = np.asarray(cf["indptr"][:], dtype=np.int64)
    if indptr.size != n_cells + 1:
        raise ValueError(f"{zip_path}: csc/indptr {indptr.size} != n_cells+1 {n_cells+1}")
    m = sparse.csr_matrix((data, indices, indptr), shape=(n_cells, n_feat))
    store.close()
    if len(keys) != n_feat:
        raise ValueError(f"{zip_path}: {len(keys)} keys vs {n_feat} features")
    if genes_only:
        ftypes = [str(t) for t in att["feature_types"]]
        keep = np.asarray([t in GENE_FEATURE_TYPES for t in ftypes], dtype=bool)
        m = m[:, keep].tocsr()
        keys = [k for k, kp in zip(keys, keep) if kp]
    return m, keys


def read_xenium_cells(parquet_gz: str) -> Tuple[np.ndarray, List[str]]:
    """``(centroids_um[n,2], cell_ids)`` from ``cells.parquet.gz``."""
    import pyarrow.parquet as pq

    with gzip.open(parquet_gz, "rb") as fh:
        raw = fh.read()
    tbl = pq.read_table(io.BytesIO(raw))
    cols = set(tbl.column_names)
    if not {"x_centroid", "y_centroid"} <= cols:
        raise ValueError(f"{parquet_gz}: columns {sorted(cols)}")
    x = np.asarray(tbl.column("x_centroid").to_numpy(), dtype=np.float64)
    y = np.asarray(tbl.column("y_centroid").to_numpy(), dtype=np.float64)
    ids = [str(v) for v in tbl.column("cell_id").to_pylist()] if "cell_id" in cols else []
    return np.column_stack([x, y]), ids


def xenium_panel_genes(zip_path: str, genes_only: bool = True) -> List[str]:
    """``feature_keys`` only -- reads the .zattrs and nothing else."""
    z = zipfile.ZipFile(zip_path)
    att = json.loads(z.read("cell_features/.zattrs").decode())
    keys = [str(k) for k in att["feature_keys"]]
    if genes_only:
        ftypes = [str(t) for t in att["feature_types"]]
        keys = [k for k, t in zip(keys, ftypes) if t in GENE_FEATURE_TYPES]
    return keys
