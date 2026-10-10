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


#: Entries compared per chunk by :func:`_is_csr_canonical_order`. 2**20 int32
#: entries is a 4 MiB slice and a 1 MiB boolean temporary, so the ordering
#: check costs megabytes on a matrix whose index arrays are hundreds of them.
_ORDER_CHECK_CHUNK = 1 << 20


def _read_matrix_coo(matrix_dir: str):
    """Read the MEX matrix as scipy COORDINATE format, (n_genes, n_spots).

    This is ``mmread``'s own output, NOT converted. The conversion is what
    costs the memory: ``mmread`` on GBM1 produces a 232.4 MiB COO at a
    MEASURED peak of 286.4 MiB in isolation and 294.5 MiB inside the loader
    -- near-optimal either way -- and every byte above that in the old reader
    was spent turning that COO into the (n_spots, n_genes) CSR
    the pipeline wants. :func:`_spots_by_genes_csr` does that turn without
    copying the two big arrays. See results/visium_load_memory_2026-10-10.md.

    A Matrix Market file written in ``array`` (dense) rather than
    ``coordinate`` form makes ``mmread`` return a plain ndarray; it is
    converted here so the one caller only ever sees coordinate format. No
    Space Ranger output is written that way (the old reader would have raised
    ``AttributeError`` on it).
    """
    from scipy.io import mmread
    from scipy.sparse import coo_matrix, issparse

    for name in ("matrix.mtx.gz", "matrix.mtx"):
        full = os.path.join(matrix_dir, name)
        if os.path.exists(full):
            if full.endswith(".gz"):
                with gzip.open(full, "rb") as fh:
                    m = mmread(fh)
            else:
                m = mmread(full)
            if not issparse(m):
                return coo_matrix(m)
            return m.tocoo() if m.format != "coo" else m
    raise FileNotFoundError(f"no matrix.mtx(.gz) in {matrix_dir}")


def _is_csr_canonical_order(rows, cols, chunk: int = _ORDER_CHECK_CHUNK) -> bool:
    """Are ``(rows, cols)`` already in CSR canonical order?

    Canonical means ``rows`` non-decreasing AND ``cols`` strictly increasing
    within each row -- i.e. sorted with no duplicate ``(row, col)`` pair. That
    is exactly the condition under which ``coo_matrix(...).tocsr()`` would
    neither reorder nor sum anything, so a CSR assembled straight from these
    arrays is BIT-IDENTICAL to the one that conversion produces. When it is
    False the caller must do the real conversion instead.

    Compared in chunks so the check never allocates an array proportional to
    ``nnz``. ``np.diff(rows)`` would cost 58 MiB on GBM1 and ``np.bincount``
    silently casts its input to ``intp``, which is 116 MiB; both were measured
    and both are why this is a loop.
    """
    n = int(rows.size)
    if n < 2:
        return True
    for lo in range(1, n, chunk):
        hi = min(lo + chunk, n)
        r, r_prev = rows[lo:hi], rows[lo - 1:hi - 1]
        if np.any(r < r_prev):
            return False
        same_row = r == r_prev
        c, c_prev = cols[lo:hi], cols[lo - 1:hi - 1]
        if np.any(c[same_row] <= c_prev[same_row]):
            return False
    return True


def _spots_by_genes_csr(coo, spots_axis: int, keep_idx, n_genes: int):
    """Turn ``mmread``'s COO into the ``(n_kept_spots, n_genes)`` CSR.

    The old reader did this in three steps, each a full copy of the matrix:

        mat = mmread(...).tocsr()        # genes x spots  (+116 MiB)
        mat[:, keep_idx]                 # a CSR column gather
        .T.tocsr()                       # CSC -> CSR     (+174 MiB)

    MEASURED on GBM1: +630.4 MiB over imports for a 174.3 MiB result, against
    +296.4 MiB for this function. Two of those steps are avoidable outright
    and one is avoidable on every section that exists here:

    * **The column gather selects everything.** ``barcodes.tsv`` in a
      *filtered* Space Ranger matrix already lists only in-tissue spots, so
      ``keep_idx`` is ``arange(n_spots)`` for all four built-in sections
      (MEASURED: 3,468 / 4,764 / 2,891 / 1,775 of 3,468 / 4,764 / 2,891 /
      1,775). It is still computed and still applied when it is NOT the
      identity -- a *raw* matrix carries ~1.5x more barcodes than there are
      in-tissue spots -- but then it is applied to the COO's flat arrays,
      where it costs the size of the output rather than a sparse gather.

    * **The transpose needs no data movement.** A MEX file lists its entries
      column-major (spot-major), rows ascending within a column, so in the
      TRANSPOSED orientation they are already in CSR canonical order. The
      gene index array is therefore already the CSR ``indices``, the value
      array is already the CSR ``data``, and ``indptr`` is a ``searchsorted``
      over the spot index array. Nothing is copied.

    That ordering is a property of the file, not a guarantee, so it is
    CHECKED (:func:`_is_csr_canonical_order`) rather than assumed, and a file
    that fails the check falls through to the ordinary conversion, which
    sorts and sums duplicates exactly as the old reader did.
    """
    from scipy.sparse import coo_matrix, csr_matrix

    spot_of = coo.col if spots_axis == 1 else coo.row
    gene_of = coo.row if spots_axis == 1 else coo.col
    data = coo.data
    n_spots_file = int(coo.shape[spots_axis])
    n_keep = int(keep_idx.size)

    if n_keep == n_spots_file and np.array_equal(
        keep_idx, np.arange(n_spots_file)
    ):
        rows = spot_of           # the identity: no gather, no copy
    else:
        # New position of each file column, -1 for a column being dropped.
        # int32 keeps the (nnz,) gather the same width as mmread's own index
        # arrays; int64 would double it.
        pos = np.full(n_spots_file, -1, dtype=np.int32)
        pos[keep_idx] = np.arange(n_keep, dtype=np.int32)
        rows = pos[spot_of]
        kept = rows >= 0
        if not bool(kept.all()):
            rows = rows[kept]
            gene_of = gene_of[kept]
            data = data[kept]
        del pos, kept

    if _is_csr_canonical_order(rows, gene_of):
        # indptr[i] = first entry whose row is >= i. Valid because `rows` is
        # non-decreasing, which is what the check above establishes. The
        # needle shares `rows`' dtype so that searchsorted cannot promote --
        # and copy -- the big array.
        indptr = np.searchsorted(
            rows, np.arange(n_keep + 1, dtype=rows.dtype), side="left"
        )
        return csr_matrix(
            (data, gene_of, indptr), shape=(n_keep, n_genes), copy=False
        )
    return coo_matrix((data, (rows, gene_of)), shape=(n_keep, n_genes)).tocsr()


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


def register_sample(
    sample_id: str,
    path: str,
    source: str = "ingested by path (unregistered provenance)",
    accession: str = "unregistered",
    diagnosis: str = "unknown",
    open_access: bool = False,
) -> Dict:
    """Add a Space Ranger directory to the sample registry AT RUNTIME.

    This is what lets a NEW patient specimen reach the reconstruction pipeline
    without editing code. ``reconstruct.real_visium.reconstruct_from_visium``
    takes a ``sample_id``, not a path, and resolves it through this module's
    registry; registering the directory under an id therefore makes the whole
    downstream chain work unchanged.

    ``path`` is stored as the sample's ``subdir``. ``os.path.join`` returns an
    absolute right-hand operand unchanged, so an absolute ``path`` bypasses
    ``_CACHE`` entirely and a relative one is still resolved under it - which is
    exactly the behaviour the four built-in samples rely on.

    Raises ``KeyError`` if ``sample_id`` is already registered with a DIFFERENT
    location, so an ingested sample can never silently shadow a built-in one and
    change what a published number was measured from. Re-registering the same
    location is a no-op.
    """
    existing = _SAMPLES.get(sample_id)
    if existing is not None:
        if os.path.abspath(str(existing["subdir"]) or _CACHE) != os.path.abspath(path):
            raise KeyError(
                f"sample_id {sample_id!r} is already registered at "
                f"{existing['subdir']!r}; refusing to repoint it to {path!r}. "
                "Use a different id."
            )
        return dict(existing, sample_id=sample_id)
    _SAMPLES[sample_id] = {
        "subdir": path,
        "source": source,
        "accession": accession,
        "landing": None,
        "open": bool(open_access),
        "diagnosis": diagnosis,
        "ingested": True,
    }
    return dict(_SAMPLES[sample_id], sample_id=sample_id)


def unregister_sample(sample_id: str) -> bool:
    """Remove a RUNTIME-registered sample. Built-in samples cannot be removed."""
    rec = _SAMPLES.get(sample_id)
    if rec is None or not rec.get("ingested"):
        return False
    del _SAMPLES[sample_id]
    return True


def sample_info(sample_id: str) -> Dict:
    """Return the static provenance record for ``sample_id`` (no file access)."""
    if sample_id not in _SAMPLES:
        raise KeyError(f"unknown sample_id {sample_id!r}; known: {sorted(_SAMPLES)}")
    rec = dict(_SAMPLES[sample_id])
    rec["sample_id"] = sample_id
    rec["present"] = _has_local(sample_id)
    return rec


def load_from_dirs(
    spatial_dir: str,
    matrix_dir: str,
    sample_id: str = "unregistered",
    marker_report: bool = True,
) -> Dict:
    """Read a Space Ranger output from an EXPLICIT ``(spatial_dir, matrix_dir)``.

    This is the whole reader; :func:`load_sample` is this function plus the
    registry lookup that finds the two directories and stamps provenance onto
    the result. Split out so the generic ingester
    (:mod:`tmesim.data.ingest`) can read a Space Ranger bundle from ANY
    directory through exactly the same code path — byte-identical output, no
    second implementation to drift.

    The returned dict carries no ``source`` / ``accession`` / ``diagnosis``:
    those are properties of the sample registry, not of the bytes, and a
    directory handed in by a user has none. :func:`load_sample` adds them.
    """
    from scipy.sparse import csr_matrix  # noqa: F401  (parity with load_sample)

    coords_px = _read_positions(spatial_dir)
    umpp, spot_diam = _read_umpp(spatial_dir)
    ensembl_ids, symbols = _read_features(matrix_dir)
    barcodes_all = _read_barcodes(matrix_dir)
    mat = _read_matrix_coo(matrix_dir)  # coordinate format, (n_genes, n_spots)

    # Which AXIS of the file is the spots. MEX is genes x spots, so normally
    # axis 1; the guard below catches a transposed/edge encoding by aligning on
    # the counts that match. The old reader physically transposed the matrix
    # here (`mat.T.tocsr()`, a full copy); an axis number does the same job,
    # because nothing between here and `_spots_by_genes_csr` reads `mat`.
    spots_axis = 1
    n_genes_file, n_spots_file = mat.shape
    if n_genes_file != len(symbols) or n_spots_file != len(barcodes_all):
        if mat.shape == (len(barcodes_all), len(symbols)):
            spots_axis = 0
            n_spots_file, n_genes_file = mat.shape

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
    expr = _spots_by_genes_csr(mat, spots_axis, keep_idx, n_genes_file)
    # Drop the COO. On the identity selection `expr` holds two of its three
    # arrays, so this frees exactly the third (the spot index: 58.1 MiB on
    # GBM1) -- and on a non-identity selection it frees all three.
    del mat
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
    spatial_dir, matrix_dir = _locate(sample_id)
    out = load_from_dirs(
        spatial_dir, matrix_dir, sample_id=sample_id, marker_report=marker_report
    )
    out["source"] = _SAMPLES[sample_id]["source"]
    out["accession"] = _SAMPLES[sample_id]["accession"]
    out["diagnosis"] = _SAMPLES[sample_id]["diagnosis"]
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
