"""Loader for a REAL single-cell-resolution spatial glioblastoma dataset.

Unlike the Visium ingest (``tmesim/data/spatial_ingest.py``), which pools
~1-10 cells into each ~55 um spot, this dataset resolves **individual cells**:
every row is one segmented cell with a measured ``(x, y)`` centroid in
micrometers and a targeted-panel expression vector. It is the ideal input for a
cell-resolved microenvironment (true single-cell positions + counts).

Dataset (OPEN, CC BY 4.0, no login, no DUA)
-------------------------------------------
10x Genomics **Xenium Human Brain — Glioblastoma section** (Xenium In Situ,
Xenium Analyzer v1.3.0). One adult human GBM FFPE section (Avaden Biosciences),
"Xenium Human Brain Preview Data".
  Landing: https://www.10xgenomics.com/datasets/xenium-human-brain-preview-data-1-standard
  Files under cf.10xgenomics.com/samples/xenium/1.3.0/
      Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon/
Panel: Xenium Human Brain base panel (hBrain_254g, 254 genes) + custom 65-gene
add-on = 319 gene-expression targets (plus negative-control / codeword features).

What is real here
-----------------
REAL (measured bytes on disk):
  * per-cell (x, y) centroids in micrometers from ``cells.csv.gz``;
  * per-cell targeted-panel counts from the MEX sparse matrix
    (``cell_feature_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}``).
Column i of the matrix corresponds to row i of ``cells.csv.gz`` (verified: the
barcodes file and the cells file share identical cell_id order).

Honest limits: this open section carries NO paired drug-response outcome label
(same gap documented in ``docs/real_data/singlecell_spatial.md``). Negative-
control features are dropped by default so ``genes`` is the real 319-gene panel.

Parsed with stdlib + scipy/pandas only (no scanpy/anndata/h5py needed).
"""
from __future__ import annotations

import gzip
import os
import tarfile
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

# ---------------------------------------------------------------------------
# Provenance (verified reachable 2026-07-26; individual files all HTTP 200)
# ---------------------------------------------------------------------------
XENIUM_GBM_BASE = (
    "https://cf.10xgenomics.com/samples/xenium/1.3.0/"
    "Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon/"
    "Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon"
)
XENIUM_GBM_FILES = {
    # local basename -> (url suffix, verified content-length bytes)
    "Xenium_GBM_cells.csv.gz": ("_cells.csv.gz", 2_167_253),
    "Xenium_GBM_cell_feature_matrix.tar.gz": ("_cell_feature_matrix.tar.gz", 14_279_458),
    "Xenium_GBM_gene_panel.json": ("_gene_panel.json", 152_134),
}
XENIUM_GBM_LANDING = (
    "https://www.10xgenomics.com/datasets/xenium-human-brain-preview-data-1-standard"
)

_DEFAULT_CACHE = os.path.join("data_cache", "xenium_gbm")
_CELLS_FILE = "Xenium_GBM_cells.csv.gz"
_MATRIX_TAR = "Xenium_GBM_cell_feature_matrix.tar.gz"
_MATRIX_DIR = "cell_feature_matrix"

# Feature-type strings in the Xenium features.tsv third column.
_GENE_TYPE = "Gene Expression"


def _read_mex(matrix_dir: str):
    """Read a 10x MEX triple into (csr features x cells, features_df, barcodes)."""
    mtx = os.path.join(matrix_dir, "matrix.mtx.gz")
    feats = os.path.join(matrix_dir, "features.tsv.gz")
    bcs = os.path.join(matrix_dir, "barcodes.tsv.gz")
    for p in (mtx, feats, bcs):
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing MEX component: {p}")
    with gzip.open(mtx, "rb") as fh:
        mat = sparse.csr_matrix(mmread(fh))  # features x cells
    features = pd.read_csv(
        feats, sep="\t", header=None, names=["ensembl", "symbol", "type"]
    )
    barcodes = pd.read_csv(bcs, sep="\t", header=None, names=["cell_id"])
    return mat, features, barcodes


def load(
    path: str = _DEFAULT_CACHE,
    genes_only: bool = True,
    dense: bool = False,
) -> Dict[str, object]:
    """Load the real Xenium single-cell-resolution GBM section.

    Parameters
    ----------
    path : str
        Directory holding the cached Xenium files (``cells.csv.gz`` and the
        extracted / to-be-extracted ``cell_feature_matrix`` MEX triple). Default
        ``data_cache/xenium_gbm``.
    genes_only : bool
        If True (default) restrict ``expression``/``genes`` to real
        ``Gene Expression`` targets, dropping negative-control and unassigned
        codeword features. If False, keep all matrix features.
    dense : bool
        If True return ``expression`` as a dense ``float32`` ndarray
        (n_cells x n_genes). If False (default) return a ``scipy`` CSR matrix in
        the same orientation to stay memory-light.

    Returns
    -------
    dict with keys
        cell_xy_um : ndarray (n_cells, 2) float64
            Per-cell (x, y) centroid in **micrometers** (measured; Xenium global
            coordinate frame). Row order matches ``expression`` rows.
        expression : csr_matrix or ndarray (n_cells, n_genes)
            Per-cell targeted-panel raw counts (cells as rows, genes as columns).
        genes : ndarray (n_genes,) of str
            Gene symbols (column order of ``expression``).
        cell_ids : ndarray (n_cells,) of str
            Xenium cell ids, aligned to the rows.
        cell_meta : DataFrame (n_cells, ...)
            Full ``cells.csv`` per-cell metadata (transcript_counts, cell_area,
            nucleus_area, control counts, ...), aligned to the rows.
        coordinate_units : "micrometer"
        source : provenance dict.

    Raises
    ------
    FileNotFoundError
        If the cache is missing. Download the three files listed in
        ``XENIUM_GBM_FILES`` into ``path`` first
        (see ``docs/real_data/singlecell_spatial.md``).
    """
    cells_path = os.path.join(path, _CELLS_FILE)
    if not os.path.exists(cells_path):
        raise FileNotFoundError(
            f"Missing cached cells file: {cells_path}. "
            "Download the Xenium GBM files into data_cache/xenium_gbm/ first "
            "(see docs/real_data/singlecell_spatial.md)."
        )

    # Ensure the MEX triple is extracted from the tarball.
    matrix_dir = os.path.join(path, _MATRIX_DIR)
    if not os.path.exists(os.path.join(matrix_dir, "matrix.mtx.gz")):
        tar_path = os.path.join(path, _MATRIX_TAR)
        if not os.path.exists(tar_path):
            raise FileNotFoundError(
                f"Missing matrix (neither extracted dir {matrix_dir} nor tar "
                f"{tar_path}). See docs/real_data/singlecell_spatial.md."
            )
        with tarfile.open(tar_path, "r:gz") as tf:
            tf.extractall(path)

    cells = pd.read_csv(cells_path)
    mat, features, barcodes = _read_mex(matrix_dir)  # mat: features x cells

    # Align cells.csv rows to matrix columns by cell_id (they share order, but
    # reindex defensively so the loader is robust to any reordering).
    if not (len(cells) == mat.shape[1] == len(barcodes)):
        raise ValueError(
            f"cell-count mismatch: cells.csv={len(cells)}, "
            f"matrix cols={mat.shape[1]}, barcodes={len(barcodes)}"
        )
    if not (cells["cell_id"].values == barcodes["cell_id"].values).all():
        order = {cid: j for j, cid in enumerate(barcodes["cell_id"].values)}
        col_idx = cells["cell_id"].map(order).to_numpy()
        mat = mat[:, col_idx]

    # Select feature rows.
    if genes_only:
        keep = (features["type"] == _GENE_TYPE).to_numpy()
    else:
        keep = np.ones(len(features), dtype=bool)
    mat = mat[keep, :]
    genes = features.loc[keep, "symbol"].to_numpy().astype(str)

    # Orient as cells x genes.
    expr = mat.T.tocsr()  # cells x genes
    if dense:
        expr = expr.toarray().astype(np.float32)

    cell_xy_um = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=np.float64)

    return {
        "cell_xy_um": cell_xy_um,
        "expression": expr,
        "genes": genes,
        "cell_ids": cells["cell_id"].to_numpy().astype(str),
        "cell_meta": cells,
        "coordinate_units": "micrometer",
        "source": {
            "dataset": "10x Xenium Human Brain — Glioblastoma section (v1.3.0)",
            "accession": "10x Genomics Xenium Human Brain Preview Data",
            "sample": "Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon",
            "assay": "Xenium In Situ (single-cell-resolution spatial)",
            "panel": "hBrain_254g base + 65-gene add-on (319 gene targets)",
            "license": "CC BY 4.0 (open, no login/DUA)",
            "landing_page": XENIUM_GBM_LANDING,
            "has_outcome_label": False,
        },
    }


if __name__ == "__main__":  # pragma: no cover - smoke test
    d = load()
    xy = d["cell_xy_um"]
    print(f"n_cells       : {xy.shape[0]:,}")
    print(f"coord units   : {d['coordinate_units']}")
    print(f"x range um    : {xy[:, 0].min():.1f} .. {xy[:, 0].max():.1f}")
    print(f"y range um    : {xy[:, 1].min():.1f} .. {xy[:, 1].max():.1f}")
    print(f"panel genes   : {len(d['genes'])}")
    print(f"expression    : {type(d['expression']).__name__} {d['expression'].shape}")
    print(f"sample genes  : {list(d['genes'][:8])}")
