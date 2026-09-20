"""Load the REAL GBM single-cell RNA-seq reference atlas.

Source: Darmanis et al. 2017, *Cell Reports* 21(5):1399-1410,
"Single-Cell RNAseq analysis of diffuse neoplastic infiltrating cells at the
migrating front of human glioblastoma."  GEO accession **GSE84465**
(public / open access; NCBI GEO data are freely redistributable).

This is a full-length (SmartSeq2, plate-based) scRNA-seq atlas of 4 primary
human GBM tumours, sampling both the tumour core and the peritumoral
("periphery") brain, with per-cell cell-type annotations published by the
authors in the GEO series-matrix metadata.  It serves as the **deconvolution
reference**: a genes x cells expression matrix plus a real cell-type label per
cell, so spatial (Visium) spots can be resolved into real cell-type fractions
and malignant vs. non-malignant content.

Cached files (see ``docs/real_data/scref_gbm.md``):
  * ``GSE84465_GBM_All_data.csv.gz``     - genes x cells raw HTSeq counts.
  * ``GSE84465_series_matrix.txt.gz``    - per-cell annotations (cell type,
    tissue [Tumor/Periphery], patient id, neoplastic flag, tSNE cluster).

The two are joined on the cell barcode ``<plate_id>.<well>`` which is the column
name in the counts matrix (e.g. ``1001000173.G8``).

Only scipy/pandas/stdlib are used; no scanpy/anndata required.

MEMORY
------
This loader is the largest single allocation in the product, and the grid
ceiling (``TMESIM_MAX_GRID_CELLS``) structurally cannot catch it: the cost is a
function of the REFERENCE PANEL, not of the requested grid, so shrinking
``dx_um`` does not shrink it by one byte. Every Visium section goes through
here. See :func:`load_reference` and :func:`check_memory_budget`.
"""
from __future__ import annotations

import gzip
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Default cache location (repo-root/data_cache/scref_gbm).
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
_DEFAULT_CACHE = os.path.join(_REPO_ROOT, "data_cache", "scref_gbm")

_COUNTS_FILE = "GSE84465_GBM_All_data.csv.gz"
_META_FILE = "GSE84465_series_matrix.txt.gz"

# HTSeq summary rows that trail the real gene rows in the counts CSV; not genes.
_NON_GENE_ROWS = {
    "no_feature",
    "ambiguous",
    "too_low_aQual",
    "not_aligned",
    "alignment_not_unique",
}

# --------------------------------------------------------------------------- #
# Memory: the measured cost of this loader, and the refusal built on it
# --------------------------------------------------------------------------- #
# EVERY number in this block is a measurement, taken on 2026-09-20 on darwin
# 24.6.0 / CPython 3.9 / pandas 2.3.3 / numpy 2.0.2, on the real cached
# GSE84465 counts file (20.5 MB gzip -> 23,465 x 3,589). They are recorded in
# results/scref_memory_2026-09-20.md together with the scripts that produced
# them. Peak is ``resource.getrusage(RUSAGE_SELF).ru_maxrss``, which on darwin
# is in BYTES (on Linux it is KiB); :func:`_maxrss_bytes` normalises that.
#
# NOTHING here is fitted or extrapolated from a different quantity. The one
# derived figure, ``MEASURED_PEAK_BYTES_PER_ELEMENT``, is the measured peak
# divided by the measured element count of the same load -- it exists so that a
# DIFFERENT reference file is sized from its own shape rather than from a
# constant that silently stops being true.

#: Measured peak RSS rise of ``load_reference()`` AS IT NOW STANDS, in MiB.
#: Five cold runs: 609.8 / 576.1 / 593.2 / 592.5 / 590.4. The MAXIMUM is taken
#: rather than the mean, because this number is a refusal threshold and a
#: threshold set at the average is wrong half the time in the direction that
#: kills the container.
MEASURED_LOAD_PEAK_MIB = 610.0

#: Measured peak RSS rise of the code this replaced, in MiB. Five cold runs:
#: 1586.3 / 1612.8 / 1850.2 / 1960.6 / 1599.5 -- the spread is the allocator
#: not returning freed arenas, and the maximum is again what is recorded.
#: Kept because a refusal's arithmetic is only meaningful next to what it saved.
MEASURED_LEGACY_PEAK_MIB = 1960.6

#: Measured resident size of the returned int32 frame, in MiB (was 643.8 as
#: int64). This is what a process PAYS TO KEEP the reference, and it is why
#: this module does NOT cache the frame -- see :func:`load_reference`.
MEASURED_RESIDENT_MIB = 322.6

#: The shape the three figures above were measured on.
MEASURED_SHAPE = (23465, 3589)

#: Measured peak bytes per matrix element: MEASURED_LOAD_PEAK_MIB over the
#: element count of the same load. Used to size a reference whose shape differs
#: from ``MEASURED_SHAPE``; on the shipped reference it reproduces the measured
#: 610 MiB exactly, because that is where it came from.
MEASURED_PEAK_BYTES_PER_ELEMENT = (
    MEASURED_LOAD_PEAK_MIB * 1024 * 1024 / float(MEASURED_SHAPE[0] * MEASURED_SHAPE[1])
)

#: Row-chunk for the bounded parse. MEASURED at 500 / 2,000 / 8,000 rows:
#: peak +479.5 / +451.8 / +585.9 MiB in 5.51 / 3.42 / 3.11 s. 2,000 is the
#: measured memory optimum and within 10 % of the best time.
_PARSE_CHUNK_ROWS = 2000

#: Operator escape hatch. Set to 1/true/yes to run the load anyway on a host
#: the refusal below judges too small. Named in the refusal message, because a
#: refusal with no documented way past it is a bug report waiting to happen.
_OVERRIDE_ENV = "TMESIM_SCREF_ALLOW_OVERSIZE"


class ReferenceTooLargeForHost(MemoryError):
    """The scRNA reference will not fit in this host's declared memory.

    A ``MemoryError`` subclass on purpose: it IS the memory error, raised
    before the allocation instead of by the OOM killer afterwards, and any
    caller already handling ``MemoryError`` keeps handling it.
    """


def _maxrss_bytes() -> int:
    """High-water RSS of this process in BYTES.

    ``ru_maxrss`` units are platform-dependent -- BYTES on darwin/macOS, KiB on
    Linux -- and reporting one as the other is a 1,024x lie in a message whose
    only job is to be believed.
    """
    import resource

    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(v) if sys.platform == "darwin" else int(v) * 1024


def _current_rss_bytes() -> Optional[int]:
    """CURRENT resident size in bytes, or ``None`` if it cannot be read.

    ``/proc/self/statm`` first (Linux, i.e. the container that matters) because
    it is a read of a pseudo-file; ``ps`` second (darwin) because there is no
    ``/proc`` there. ``None`` rather than a guess when neither works: the
    headroom calculation below says so instead of inventing a number.
    """
    try:
        with open("/proc/self/statm") as fh:
            pages = int(fh.read().split()[1])
        return pages * os.sysconf("SC_PAGE_SIZE")
    except (OSError, ValueError, IndexError):
        pass
    try:
        import subprocess

        out = subprocess.check_output(
            ["ps", "-o", "rss=", "-p", str(os.getpid())], stderr=subprocess.DEVNULL
        )
        return int(out.strip()) * 1024  # ps reports KiB
    except Exception:
        return None


def _host_memory_limit() -> Dict:
    """What the host grants this process, via the server's own detector.

    Imported lazily and defensively: ``app.server.account`` imports FastAPI at
    module scope, and ``tmesim`` must stay importable on a core-only install
    with no ``app`` extra (the same reason ``tmesim/protocol/executor.py``
    imports ``app.server.cache`` inside a function). There is deliberately no
    second copy of the cgroup-reading logic here -- a duplicated limit detector
    that drifts from the one the rest of the server uses is precisely the
    defect this repository keeps finding.

    On an install where the import fails there is no server, therefore no
    visitor-driven OOM, and the answer is an honest "not declared".
    """
    try:
        from app.server.account import memory_limit
    except Exception:
        return {"bytes": None, "source": "not declared (app.server.account unavailable)"}
    try:
        return memory_limit()
    except Exception:
        return {"bytes": None, "source": "not declared (detector raised)"}


def memory_headroom() -> Dict:
    """How much room this process has left, and how confident that number is.

    Returns ``{"limit_bytes", "rss_bytes", "headroom_bytes", "source",
    "override"}``. ``limit_bytes`` is ``None`` when nothing declares a limit,
    and then ``headroom_bytes`` is ``None`` too -- an undeclared limit is
    reported as unknown, never as generous.

    Shared by this module's refusal and by the Visium path's
    (``tmesim/reconstruct/deconvolve.py``) so that both quote the same detected
    limit from the same detector.
    """
    override = str(os.environ.get(_OVERRIDE_ENV, "")).strip().lower() in (
        "1", "true", "yes", "on"
    )
    mem = _host_memory_limit()
    limit = mem.get("bytes")
    rss = _current_rss_bytes()
    if not limit:
        headroom = None
    elif rss is None:
        headroom = limit
    else:
        headroom = limit - rss
    return {
        "limit_bytes": limit or None,
        "rss_bytes": rss,
        "headroom_bytes": headroom,
        "source": mem.get("source"),
        "override": override,
    }


def estimate_load_bytes(n_genes: int, n_cells: int) -> int:
    """Measured peak bytes the bounded parse needs for an ``n_genes x n_cells``
    counts matrix. See the constants above for provenance."""
    return int(round(MEASURED_PEAK_BYTES_PER_ELEMENT * float(n_genes) * float(n_cells)))


def check_memory_budget(n_genes: int, n_cells: int) -> None:
    """Refuse the load if this host's declared memory cannot hold its peak.

    This exists because ``TMESIM_MAX_GRID_CELLS`` **cannot** catch this
    allocation. That ceiling counts VOXELS; this cost is a function of the
    reference panel and is identical at every ``dx_um``, so a visitor who picks
    a Visium section defeats the grid ceiling without trying.

    Refuses only when a limit is actually DECLARED (cgroup v2 -> v1 ->
    ``TMESIM_HOST_MEMORY_MB``, read by ``app/server/account.py::memory_limit``)
    and the measured peak does not fit in the headroom. An undeclared limit is
    never treated as small: this raises nothing rather than refusing a run on a
    machine it knows nothing about.

    Raises
    ------
    ReferenceTooLargeForHost
        With the measured requirement, the detected limit and its source, the
        headroom, and the two things an operator can do about it.
    """
    if str(os.environ.get(_OVERRIDE_ENV, "")).strip().lower() in ("1", "true", "yes", "on"):
        return

    mem = _host_memory_limit()
    limit = mem.get("bytes")
    if not limit:
        return  # nothing declared -> nothing to compare against, so no refusal

    need = estimate_load_bytes(n_genes, n_cells)
    rss = _current_rss_bytes()
    headroom = (limit - rss) if rss is not None else limit
    if need <= headroom:
        return

    mib = 1024.0 * 1024.0
    if rss is None:
        where = (f"the whole {limit / mib:,.0f} MiB limit (this process's current "
                 f"RSS could not be read, so none of it is assumed to be in use)")
    else:
        where = (f"{headroom / mib:,.0f} MiB of headroom ({limit / mib:,.0f} MiB limit "
                 f"minus {rss / mib:,.0f} MiB already resident)")

    raise ReferenceTooLargeForHost(
        f"Refusing to load the GBM scRNA-seq deconvolution reference "
        f"(counts matrix {n_genes:,} rows x {n_cells:,} cells): parsing it peaks "
        f"at a MEASURED "
        f"{need / mib:,.0f} MiB and this host declares {where}. "
        f"Limit source: {mem.get('source')}. "
        f"\n\nThis is NOT a grid-size problem and TMESIM_MAX_GRID_CELLS cannot "
        f"prevent it: every Visium section loads this same reference panel at "
        f"every dx_um, so choosing a coarser grid does not reduce it. "
        f"\n\nOptions: (1) run the Visium sections on an instance with more "
        f"memory -- the reference load alone needs ~{MEASURED_LOAD_PEAK_MIB:,.0f} MiB "
        f"and the full Visium reconstruction was measured far above that "
        f"(see results/scref_memory_2026-09-20.md); (2) use a Xenium section, "
        f"which does not touch this reference at all; or (3) set "
        f"{_OVERRIDE_ENV}=1 to attempt the load anyway and accept that the "
        f"container may be OOM-killed."
    )


#: One-entry memo for :func:`_counts_shape`, keyed by (path, size, mtime_ns).
#: This one IS safe to cache: it is ~27k short strings (a couple of MB), not
#: the 322.6 MiB counts matrix, and the key changes if the file does. Without
#: it a single Visium run reads the symbol column three times at 0.72 s each.
_SHAPE_MEMO: Dict[Tuple, Tuple] = {}


def _counts_shape(path: str) -> Tuple[int, int, List[str], List[str]]:
    """``(n_rows, n_cols, gene_symbols, cell_barcodes)`` read CHEAPLY.

    MEASURED: +30.0 MiB and 0.72 s, against +1,444 MiB and 1.94 s to parse the
    whole matrix. The symbols returned here were verified ``array_equal`` to
    ``pd.read_csv(...).index`` on the real file, so a caller that needs only
    the gene vocabulary never has to pay for the counts.
    """
    try:
        st = os.stat(path)
        key = (os.path.abspath(path), st.st_size, st.st_mtime_ns)
    except OSError:
        key = None
    if key is not None and key in _SHAPE_MEMO:
        return _SHAPE_MEMO[key]
    with gzip.open(path, "rt", errors="replace") as fh:
        cells = [c.strip('"') for c in fh.readline().split()]
    idx = pd.read_csv(path, sep=r"\s+", index_col=0, usecols=[0]).index
    genes = [str(g) for g in idx]
    got = (len(genes), len(cells), genes, cells)
    if key is not None:
        _SHAPE_MEMO.clear()  # one reference at a time; never grows unbounded
        _SHAPE_MEMO[key] = got
    return got


def _read_counts_bounded(
    path: str,
    shape: Optional[Tuple] = None,
    keep_rows: Optional[np.ndarray] = None,
    keep_cols: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Parse the counts matrix with a BOUNDED peak, bit-exactly.

    The replaced code was::

        expr = pd.read_csv(path, sep=r"\\s+", index_col=0)   # dense int64
        expr = expr.loc[~expr.index.isin(_NON_GENE_ROWS)]   # full-matrix copy
        expr = expr[common]                                 # full-matrix copy

    -- a MEASURED +1,444 MiB for the parse alone (2.2x the 643.8 MiB frame it
    returned, because the C parser builds 3,589 per-column arrays and then
    consolidates them into one block) and then two more whole-matrix copies to
    drop 5 rows and reorder the columns.

    Three changes, all value-preserving:

    * **int32 instead of int64.** The largest count in the file is 77,122, so
      every value is exactly representable; the frame's only consumers convert
      to float64 before doing arithmetic, and numpy/pandas integer reductions
      accumulate in int64 regardless of the input width, so no sum narrows.
    * **chunked into a preallocated array** instead of parsed whole and then
      consolidated, which removes the duplicate-the-whole-matrix step.
    * **the row and column selections applied DURING the fill**, so the two
      post-hoc copies never happen. ``keep_cols`` is positional, which also
      makes a duplicated barcode select one column rather than silently
      fanning out the way label-based ``expr[common]`` would.

    MEASURED for the parse in isolation: +451.8 MiB peak, against +1,444 to
    +1,506 MiB for the bare ``read_csv`` it replaces.

    Verified against the original on the real cached file:
    ``np.array_equal(new.astype(int64), old)`` is **True**, and so is every
    downstream deconvolution output (see the results file).
    """
    n_genes, n_cells, _, cells = shape if shape is not None else _counts_shape(path)
    if keep_rows is None:
        keep_rows = np.ones(n_genes, dtype=bool)
    if keep_cols is None:
        keep_cols = np.arange(n_cells, dtype=np.intp)
    keep_cols = np.asarray(keep_cols, dtype=np.intp)

    out = np.empty((int(keep_rows.sum()), keep_cols.size), dtype=np.int32)
    dtypes = {c: np.int32 for c in cells}
    symbols: List[str] = []
    src = 0  # row offset in the FILE
    dst = 0  # row offset in ``out``
    for chunk in pd.read_csv(
        path, sep=r"\s+", index_col=0, dtype=dtypes, chunksize=_PARSE_CHUNK_ROWS
    ):
        n = chunk.shape[0]
        mask = keep_rows[src:src + n]
        if mask.any():
            sel = chunk.to_numpy()[mask][:, keep_cols]
            out[dst:dst + sel.shape[0]] = sel
            symbols.extend(
                str(s) for s, k in zip(chunk.index, mask) if k
            )
            dst += sel.shape[0]
        src += n
        del chunk
    if src != n_genes:  # pragma: no cover - a malformed cache, not a code path
        raise ValueError(
            f"{path}: parsed {src} rows but the gene-symbol column has {n_genes}"
        )
    return pd.DataFrame(
        out,
        index=pd.Index(symbols),
        columns=[cells[i] for i in keep_cols],
        copy=False,
    )


def reference_fingerprint(cache_dir: str = None) -> Tuple:
    """A cheap identity for the cached reference: (path, size, mtime_ns) x2.

    Used as a cache key by the Visium path so that REPLACING the reference
    files invalidates a memoised result instead of serving a stale one.
    """
    cache_dir = _DEFAULT_CACHE if cache_dir is None else cache_dir
    out = []
    for name in (_COUNTS_FILE, _META_FILE):
        p = os.path.join(cache_dir, name)
        try:
            st = os.stat(p)
            out.append((p, st.st_size, st.st_mtime_ns))
        except OSError:
            out.append((p, None, None))
    return tuple(out)


def reference_gene_symbols(cache_dir: str = _DEFAULT_CACHE) -> np.ndarray:
    """The reference's gene vocabulary, WITHOUT parsing the counts matrix.

    Returns exactly ``load_reference(cache_dir)[0].index`` as an object array
    (HTSeq summary rows dropped, order preserved, duplicates kept) for a
    MEASURED +30 MiB instead of +452 MiB.

    This exists because the Visium path used to load the entire 23,460 x 3,589
    matrix purely to build ``set(expr_df.index)`` -- ~1.5 GiB of parsing for a
    list of strings.
    """
    counts_path = os.path.join(cache_dir, _COUNTS_FILE)
    if not os.path.exists(counts_path):
        raise FileNotFoundError(
            f"Missing cached reference file: {counts_path}. "
            "Download GSE84465 into data_cache/scref_gbm/ first "
            "(see docs/real_data/scref_gbm.md)."
        )
    symbols = _counts_shape(counts_path)[2]
    return np.asarray(
        [g for g in symbols if g not in _NON_GENE_ROWS], dtype=object
    )


def _parse_series_matrix(path: str) -> pd.DataFrame:
    """Parse per-cell annotations from a GEO series-matrix .txt.gz.

    Returns a DataFrame indexed by the ``<plate_id>.<well>`` cell barcode with
    columns: cell_type, tissue, patient_id, neoplastic, tsne_cluster, plate_id,
    well, gsm.
    """
    titles: List[str] = []
    gsms: List[str] = []
    # characteristic key -> list of per-sample values
    chars: Dict[str, List[str]] = {}

    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [c.strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                gsms = [c.strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [c.strip('"') for c in line.rstrip("\n").split("\t")[1:]]
                # Each characteristics line is "key: value"; key is stable per line.
                key = None
                for v in vals:
                    if ": " in v:
                        key = v.split(": ", 1)[0].strip()
                        break
                if key is None:
                    continue
                parsed = [v.split(": ", 1)[1] if ": " in v else "" for v in vals]
                # A GSE can repeat the same key across lines; keep them distinct.
                base = key
                n = 1
                while key in chars:
                    n += 1
                    key = f"{base}_{n}"
                chars[key] = parsed

    n_samples = len(titles)
    plate = chars.get("plate id", [""] * n_samples)
    well = chars.get("well", [""] * n_samples)
    barcodes = [f"{p}.{w}" for p, w in zip(plate, well)]

    df = pd.DataFrame(
        {
            "cell_type": chars.get("cell type", [""] * n_samples),
            "tissue": chars.get("tissue", [""] * n_samples),
            "patient_id": chars.get("patient id", [""] * n_samples),
            "neoplastic": chars.get("neoplastic", [""] * n_samples),
            "tsne_cluster": chars.get("tsne cluster", [""] * n_samples),
            "plate_id": plate,
            "well": well,
            "gsm": gsms if len(gsms) == n_samples else [""] * n_samples,
        },
        index=barcodes,
    )
    df.index.name = "cell_id"
    return df


def load_reference(
    cache_dir: str = _DEFAULT_CACHE,
    max_cells: Optional[int] = None,
    random_state: int = 0,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Load the Darmanis GBM scRNA-seq reference from the cache.

    Parameters
    ----------
    cache_dir : str
        Directory holding the two cached GSE84465 files.
    max_cells : int, optional
        If given, return a random subsample of this many cells (stratified is
        not applied; a documented uniform subsample). ``None`` returns all.
    random_state : int
        Seed for the optional subsample.

    Returns
    -------
    expression : pandas.DataFrame
        genes x cells raw counts (index = gene symbol, columns = cell barcode).
    cell_labels : pandas.DataFrame
        per-cell annotations aligned to ``expression.columns`` (cell_type,
        tissue, patient_id, neoplastic, tsne_cluster, ...).
    gene_names : numpy.ndarray
        1-D array of gene symbols (== ``expression.index``).

    Raises
    ------
    FileNotFoundError
        If the cache files are missing (run the acquisition step first).
    ReferenceTooLargeForHost
        If this host DECLARES a memory limit that the measured peak of this
        load does not fit in. Raised BEFORE the allocation, with the measured
        requirement and the detected limit in the message.

    Notes
    -----
    **This function is deliberately NOT cached.** Memoising it would remove the
    repeat parse (4.5 s) but pin a MEASURED 322.7 MiB resident for the life of
    the process, which on the 512 MB instance this product deploys to is worse
    than the problem it solves. The Visium path instead caches the small thing
    it actually wants -- see
    ``tmesim/reconstruct/deconvolve.py::build_signatures``, whose memoised
    result is a few hundred kB.
    """
    counts_path = os.path.join(cache_dir, _COUNTS_FILE)
    meta_path = os.path.join(cache_dir, _META_FILE)
    for p in (counts_path, meta_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"Missing cached reference file: {p}. "
                "Download GSE84465 into data_cache/scref_gbm/ first "
                "(see docs/real_data/scref_gbm.md)."
            )

    # REFUSE BEFORE ALLOCATING. The shape is read cheaply (+30 MiB) so the
    # refusal can quote the real matrix rather than a hardcoded guess, and so
    # that a host too small for the parse never starts it.
    shape = _counts_shape(counts_path)
    n_genes, n_cells = shape[0], shape[1]
    check_memory_budget(n_genes, n_cells)

    # The annotations are parsed FIRST (a 145 kB file) so that both selections
    # -- drop the HTSeq summary rows, keep only barcodes the metadata knows --
    # are known before the counts matrix is allocated, and can therefore be
    # applied while it is filled. Doing them afterwards, as this used to, cost
    # two extra whole-matrix copies: MEASURED +852 MiB peak for the function
    # against +452 MiB for the parse it wraps.
    meta = _parse_series_matrix(meta_path)
    symbols, cells = shape[2], shape[3]
    keep_rows = np.array([s not in _NON_GENE_ROWS for s in symbols], dtype=bool)
    meta_ids = set(meta.index)
    keep_cols = np.array(
        [j for j, c in enumerate(cells) if c in meta_ids], dtype=np.intp
    )

    # Counts: first column is the gene symbol (unnamed in the header). Despite
    # the ``.csv`` name the file is whitespace-delimited.
    expr = _read_counts_bounded(
        counts_path, shape=shape, keep_rows=keep_rows, keep_cols=keep_cols
    )
    expr.index = expr.index.astype(str)
    labels = meta.loc[list(expr.columns)].copy()

    if max_cells is not None and max_cells < expr.shape[1]:
        rng = np.random.default_rng(random_state)
        keep = rng.choice(expr.shape[1], size=max_cells, replace=False)
        keep.sort()
        cols = expr.columns[keep]
        expr = expr[cols]
        labels = labels.loc[cols]

    gene_names = np.asarray(expr.index, dtype=object)
    return expr, labels, gene_names


if __name__ == "__main__":
    X, y, genes = load_reference()
    print(f"expression: {X.shape[0]} genes x {X.shape[1]} cells")
    print(f"cell-type vocabulary: {sorted(y['cell_type'].unique())}")
    print(y["cell_type"].value_counts().to_string())
    print(f"tissue split: {dict(y['tissue'].value_counts())}")
    print(f"patients: {sorted(y['patient_id'].unique())}")
