"""Load REAL 10x Visium (spot-level spatial transcriptomics) glioblastoma
sections from the local cache under ``data_cache/visium_gbm/``.

This generalises the single-sample reader in
``tmesim/data/spatial_ingest.py`` to the several REAL Visium GBM sections we
have downloaded (no synthetic stand-ins). Every sample here is genuine
Space Ranger output: a MEX sparse expression matrix
(``matrix.mtx.gz`` + ``features.tsv.gz`` + ``barcodes.tsv.gz``) plus the
``spatial/`` folder carrying per-spot pixel coordinates
(``tissue_positions*.csv``) and the ``scalefactors_json.json`` needed to convert
pixels to micrometers.

Samples (see ``docs/real_data/visium_gbm.md`` for full provenance):

  * ``10x_parent_gbm`` - 10x Genomics public "Human Glioblastoma, Whole
    Transcriptome Analysis" (Visium v1, Space Ranger 1.2.0). Open, no login.
  * ``GBM1`` / ``GBM2`` / ``GBM3`` - three independent primary human GBM
    sections from GEO series **GSE194329** (Ren et al.; Space Ranger 1.1.0),
    supplementary Space Ranger outputs, openly downloadable from NCBI GEO.

Public API
----------
``list_samples()``                 -> list of sample-id strings that are present.
``load_sample(sample_id)``         -> dict with (at least) the keys
                                      ``coords_um``, ``expression``, ``genes``,
                                      ``barcodes``.

Only scipy/numpy/stdlib are used; no scanpy/anndata required.
"""
from __future__ import annotations

import gzip
import json
import os
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Cache location + sample registry
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
_CACHE = os.path.join(_REPO_ROOT, "data_cache", "visium_gbm")

# Visium v1 geometry: 55 um spot diameter. um-per-pixel = 55 / spot_diameter_px.
_SPOT_DIAMETER_UM = 55.0

# GBM marker panel used only for the presence/QC report (not for phenotyping here).
GBM_MARKERS: Dict[str, List[str]] = {
    "tumor": ["SOX2", "OLIG2", "EGFR", "MKI67", "GFAP"],
    "immune": ["CD68", "AIF1", "CD3E", "PTPRC"],
    "endothelial": ["PECAM1", "VEGFA"],
    "oligodendrocyte": ["MBP", "PLP1"],
}

# sample_id -> provenance + on-disk location (relative to _CACHE). ``subdir=""``
# means the spatial/ and matrix dirs live at the cache root (the 10x Parent sample).
_SAMPLES: Dict[str, Dict] = {
    "10x_parent_gbm": {
        "subdir": "",
        "source": "10x Genomics public Visium (Space Ranger 1.2.0)",
        "accession": "Parent_Visium_Human_Glioblastoma",
        "landing": (
            "https://www.10xgenomics.com/datasets/"
            "human-glioblastoma-whole-transcriptome-analysis-1-standard-1-2-0"
        ),
        "open": True,
        "diagnosis": "Glioblastoma (primary resection, fresh-frozen)",
    },
    "GBM1": {
        "subdir": "GBM1",
        "source": "GEO GSE194329 (Space Ranger 1.1.0) supplementary Space Ranger output",
        "accession": "GSM5833533 (GBM1)",
        "landing": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194329",
        "open": True,
        "diagnosis": "Glioblastoma, IDH-wildtype (primary)",
    },
    "GBM2": {
        "subdir": "GBM2",
        "source": "GEO GSE194329 (Space Ranger 1.1.0) supplementary Space Ranger output",
        "accession": "GSM5833534 (GBM2)",
        "landing": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194329",
        "open": True,
        "diagnosis": "Glioblastoma, IDH-wildtype (recurrent)",
    },
    "GBM3": {
        "subdir": "GBM3",
        "source": "GEO GSE194329 (Space Ranger 1.1.0) supplementary Space Ranger output",
        "accession": "GSM5833535 (GBM3)",
        "landing": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194329",
        "open": True,
        "diagnosis": "Glioblastoma, IDH-wildtype (primary)",
    },
}


# ---------------------------------------------------------------------------
# Directory discovery (robust to Space Ranger version differences)
# ---------------------------------------------------------------------------
def _sample_root(sample_id: str) -> str:
    if sample_id not in _SAMPLES:
        raise KeyError(
            f"unknown sample_id {sample_id!r}; known: {sorted(_SAMPLES)}"
        )
    return os.path.join(_CACHE, _SAMPLES[sample_id]["subdir"])


def _find_dir(
    root: str,
    required_file_options: Tuple[Tuple[str, ...], ...],
    exclude_names: Optional[set] = None,
) -> Optional[str]:
    """Return the first directory at/under ``root`` (depth <= 2) that contains
    at least one file from each option-group in ``required_file_options``.

    ``exclude_names`` are directory basenames NOT descended into (used to keep a
    root-level sample from picking up a sibling sample's subdirectory)."""
    exclude_names = exclude_names or set()
    candidates = [root]
    if os.path.isdir(root):
        for name in sorted(os.listdir(root)):
            if name in exclude_names:
                continue
            p = os.path.join(root, name)
            if os.path.isdir(p):
                candidates.append(p)
                for sub in sorted(os.listdir(p)):
                    sp = os.path.join(p, sub)
                    if os.path.isdir(sp):
                        candidates.append(sp)
    for d in candidates:
        if not os.path.isdir(d):
            continue
        listing = set(os.listdir(d))
        if all(any(f in listing for f in group) for group in required_file_options):
            return d
    return None


def _locate(sample_id: str) -> Tuple[str, str]:
    """Return ``(spatial_dir, matrix_dir)`` for a sample. Raises FileNotFoundError."""
    root = _sample_root(sample_id)
    # A root-level sample (subdir="") shares the cache dir with the per-sample
    # subdirs; never descend into a sibling sample's folder or the raw downloads.
    exclude = {s["subdir"] for sid, s in _SAMPLES.items() if sid != sample_id and s["subdir"]}
    exclude.add("_geo_downloads")
    spatial_dir = _find_dir(
        root,
        (
            ("scalefactors_json.json",),
            ("tissue_positions_list.csv", "tissue_positions.csv"),
        ),
        exclude_names=exclude,
    )
    matrix_dir = _find_dir(
        root,
        (("matrix.mtx.gz", "matrix.mtx"),
         ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"),
         ("barcodes.tsv.gz", "barcodes.tsv")),
        exclude_names=exclude,
    )
    if spatial_dir is None or matrix_dir is None:
        raise FileNotFoundError(
            f"sample {sample_id!r}: could not locate "
            f"{'spatial/ ' if spatial_dir is None else ''}"
            f"{'matrix/ ' if matrix_dir is None else ''}under {root!r}. "
            "Expected Space Ranger output (spatial/ + filtered_feature_bc_matrix/)."
        )
    return spatial_dir, matrix_dir


def _has_local(sample_id: str) -> bool:
    try:
        _locate(sample_id)
        return True
    except (FileNotFoundError, KeyError):
        return False


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------
def _open_maybe_gz(path_no_ext: str, *names: str):
    """Open the first existing file among ``names`` (relative to a dir), handling
    optional .gz transparently in text mode. Returns a file handle."""
    for name in names:
        full = os.path.join(path_no_ext, name)
        if os.path.exists(full):
            if full.endswith(".gz"):
                return gzip.open(full, "rt")
            return open(full, "rt")
    raise FileNotFoundError(f"none of {names} in {path_no_ext}")


def _read_positions(spatial_dir: str) -> "Dict[str, Tuple[float, float]]":
    """Read tissue positions -> {barcode: (x_px, y_px)} for in-tissue spots only.

    Handles both the headerless ``tissue_positions_list.csv`` (Space Ranger < 1.2)
    and the headered ``tissue_positions.csv`` (>= 1.2). Column order is identical:
    barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres.
    """
    coords: Dict[str, Tuple[float, float]] = {}
    fh = _open_maybe_gz(
        spatial_dir,
        "tissue_positions_list.csv",
        "tissue_positions.csv",
        "tissue_positions_list.csv.gz",
        "tissue_positions.csv.gz",
    )
    with fh:
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 6:
                continue
            if parts[0].strip().lower() in ("barcode", "barcodes"):
                continue  # header row
            barcode, in_tissue = parts[0], parts[1].strip()
            if in_tissue not in ("1", "1.0"):
                continue
            try:
                pxl_row = float(parts[4])  # y
                pxl_col = float(parts[5])  # x
            except ValueError:
                continue
            coords[barcode] = (pxl_col, pxl_row)
    return coords


def _read_umpp(spatial_dir: str) -> Tuple[float, float]:
    """Return ``(umpp, spot_diameter_fullres)`` from ``scalefactors_json.json``."""
    with open(os.path.join(spatial_dir, "scalefactors_json.json")) as fh:
        sf = json.load(fh)
    diam = float(sf["spot_diameter_fullres"])
    return _SPOT_DIAMETER_UM / diam, diam


def _read_features(matrix_dir: str) -> Tuple[List[str], List[str]]:
    """Return ``(ensembl_ids, symbols)`` in matrix-row order."""
    ensembl: List[str] = []
    symbols: List[str] = []
    fh = _open_maybe_gz(
        matrix_dir, "features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"
    )
    with fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if not cols or not cols[0]:
                continue
            ensembl.append(cols[0])
            symbols.append(cols[1] if len(cols) >= 2 else cols[0])
    return ensembl, symbols


def _read_barcodes(matrix_dir: str) -> List[str]:
    fh = _open_maybe_gz(matrix_dir, "barcodes.tsv.gz", "barcodes.tsv")
    with fh:
        return [ln.strip() for ln in fh if ln.strip()]


def _read_matrix(matrix_dir: str):
    """Read the MEX matrix as a scipy CSR of shape (n_genes, n_spots)."""
    from scipy.io import mmread

    for name in ("matrix.mtx.gz", "matrix.mtx"):
        full = os.path.join(matrix_dir, name)
        if os.path.exists(full):
            if full.endswith(".gz"):
                with gzip.open(full, "rb") as fh:
                    return mmread(fh).tocsr()
            return mmread(full).tocsr()
    raise FileNotFoundError(f"no matrix.mtx(.gz) in {matrix_dir}")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def list_samples(present_only: bool = True) -> List[str]:
    """List available Visium GBM sample ids.

    Parameters
    ----------
    present_only : if True (default) only return samples whose bytes are locally
        cached and locatable; if False return every registered sample id.
    """
    ids = list(_SAMPLES.keys())
    if not present_only:
        return ids
    return [s for s in ids if _has_local(s)]


def sample_info(sample_id: str) -> Dict:
    """Return the static provenance record for ``sample_id`` (no file access)."""
    if sample_id not in _SAMPLES:
        raise KeyError(f"unknown sample_id {sample_id!r}; known: {sorted(_SAMPLES)}")
    rec = dict(_SAMPLES[sample_id])
    rec["sample_id"] = sample_id
    rec["present"] = _has_local(sample_id)
    return rec


def load_sample(sample_id: str, marker_report: bool = True) -> Dict:
    """Load a REAL Visium GBM section from the local cache.

    Returns a dict with the keys required by the build plan plus useful extras:

    ``coords_um``   : (n_spots, 2) float64 array of in-tissue spot centroids in
                      micrometers, origin shifted to (0, 0), aligned row-for-row
                      with ``barcodes`` and the rows of ``expression``.
    ``expression``  : scipy.sparse CSR of shape (n_spots, n_genes), integer UMI
                      counts (spots as rows, genes as columns).
    ``genes``       : list[str] of gene SYMBOLS, length n_genes (column order).
    ``barcodes``    : list[str] of spot barcodes, length n_spots (row order).

    extras: ``ensembl_ids`` (list, aligned to ``genes``), ``umpp``
    (um-per-pixel), ``spot_diameter_fullres`` (px), ``coords_px`` (n_spots, 2),
    ``n_spots``, ``n_genes``, ``total_counts`` (per-spot UMI), ``source``,
    ``accession``, ``diagnosis``, and (if ``marker_report``) ``markers_present``.

    Raises
    ------
    FileNotFoundError if the sample's bytes are not cached (see
    ``docs/real_data/visium_gbm.md`` for how to fetch).
    """
    from scipy.sparse import csr_matrix

    spatial_dir, matrix_dir = _locate(sample_id)

    coords_px = _read_positions(spatial_dir)
    umpp, spot_diam = _read_umpp(spatial_dir)
    ensembl_ids, symbols = _read_features(matrix_dir)
    barcodes_all = _read_barcodes(matrix_dir)
    mat = _read_matrix(matrix_dir)  # (n_genes, n_spots) CSR

    n_genes_file, n_spots_file = mat.shape
    if n_genes_file != len(symbols) or n_spots_file != len(barcodes_all):
        # mmread gives (rows, cols); MEX is genes x spots. Guard against a
        # transposed/edge encoding by aligning on the counts that match.
        if mat.shape == (len(barcodes_all), len(symbols)):
            mat = mat.T.tocsr()
            n_genes_file, n_spots_file = mat.shape

    # keep only spots present in BOTH the matrix and the in-tissue position list
    keep_idx = [i for i, bc in enumerate(barcodes_all) if bc in coords_px]
    if not keep_idx:
        raise ValueError(
            f"sample {sample_id!r}: no barcodes shared between matrix and "
            "in-tissue tissue_positions"
        )
    keep_idx = np.asarray(keep_idx, dtype=int)
    barcodes = [barcodes_all[i] for i in keep_idx]

    xs_px = np.array([coords_px[b][0] for b in barcodes], dtype=np.float64)
    ys_px = np.array([coords_px[b][1] for b in barcodes], dtype=np.float64)
    coords_px_arr = np.column_stack([xs_px, ys_px])

    x_um = (xs_px - xs_px.min()) * umpp
    y_um = (ys_px - ys_px.min()) * umpp
    coords_um = np.column_stack([x_um, y_um])

    # expression -> (n_spots, n_genes)
    expr = mat[:, keep_idx].T.tocsr()  # (n_spots, n_genes)
    total_counts = np.asarray(expr.sum(axis=1)).ravel().astype(np.float64)

    out: Dict = {
        "sample_id": sample_id,
        "coords_um": coords_um,
        "coords_px": coords_px_arr,
        "expression": expr,
        "genes": symbols,
        "ensembl_ids": ensembl_ids,
        "barcodes": barcodes,
        "umpp": umpp,
        "spot_diameter_fullres": spot_diam,
        "n_spots": int(expr.shape[0]),
        "n_genes": int(expr.shape[1]),
        "total_counts": total_counts,
        "source": _SAMPLES[sample_id]["source"],
        "accession": _SAMPLES[sample_id]["accession"],
        "diagnosis": _SAMPLES[sample_id]["diagnosis"],
        "spatial_dir": spatial_dir,
        "matrix_dir": matrix_dir,
    }

    if marker_report:
        sym_to_col = {}
        for j, s in enumerate(symbols):
            sym_to_col.setdefault(s, j)
        present: Dict[str, List[str]] = {}
        for cls, genes in GBM_MARKERS.items():
            present[cls] = [g for g in genes if g in sym_to_col]
        out["markers_present"] = present

    return out


if __name__ == "__main__":
    print("Available Visium GBM samples (present locally):", list_samples())
    for sid in list_samples():
        d = load_sample(sid)
        res_um = 100.0  # Visium v1 center-to-center spot pitch is 100 um
        print(
            f"\n[{sid}] {d['accession']}\n"
            f"  source            : {d['source']}\n"
            f"  n_spots           : {d['n_spots']}\n"
            f"  n_genes           : {d['n_genes']}\n"
            f"  umpp              : {d['umpp']:.4f} um/px "
            f"(spot_diameter={d['spot_diameter_fullres']:.1f} px = 55 um)\n"
            f"  spot pitch        : ~{res_um:.0f} um center-to-center (Visium v1)\n"
            f"  tissue extent     : {d['coords_um'][:,0].max():.0f} x "
            f"{d['coords_um'][:,1].max():.0f} um\n"
            f"  median UMI/spot   : {np.median(d['total_counts']):.0f}\n"
            f"  markers present   : "
            + ", ".join(
                f"{k}={len(v)}/{len(GBM_MARKERS[k])}"
                for k, v in d["markers_present"].items()
            )
        )
