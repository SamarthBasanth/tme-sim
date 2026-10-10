"""Spot deconvolution: real Visium spots -> real cell-type fractions.

This turns each REAL 10x Visium spot's expression profile into a vector of real
cell-type fractions, using the REAL Darmanis GBM single-cell atlas
(GSE84465, ``tmesim.data.real.scref_gbm``) as the deconvolution reference.
There is NO synthetic data anywhere in this module.

Method (deliberately simple + honest):

1. Load the Darmanis genes x cells reference and its per-cell ``cell_type``
   labels. Collapse the 7 published labels onto the coarse ``CellType`` classes
   of the core model:

       Neoplastic                               -> TUMOR
       Immune cell                              -> IMMUNE
       Vascular                                 -> ENDOTHELIAL
       Astrocyte / Oligodendrocyte / OPC / Neuron -> STROMAL

2. Restrict both reference and spot to the SHARED gene set (symbols present in
   both the Visium panel and the reference), library-size-normalize each single
   cell over those shared genes, then average per collapsed class to get one
   mean-expression *signature* per cell type.
3. For each spot, solve ``spot_expr ~= signatures @ w`` with non-negative least
   squares (:func:`scipy.optimize.nnls`) and normalize ``w`` to sum to 1 -> the
   spot's cell-type fractions.

Public API
----------
``build_signatures(shared_genes)``  -> (signatures[G x K], type_names)
``deconvolve_spot(spot_expr, signatures)`` -> fractions[K] (sum 1)
``deconvolve_sample(sample_id)``    -> (coords_um[S,2], fractions[S,K], type_names)
"""
from __future__ import annotations

from collections import OrderedDict
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.optimize import nnls

from ..core.microenvironment import CellType
from ..data.real.scref_gbm import (_OVERRIDE_ENV, load_reference,
                                   memory_headroom, reference_fingerprint,
                                   reference_gene_symbols)
from ..data.real.visium import load_sample

# --------------------------------------------------------------------------- #
# What a Visium deconvolution actually costs, and the refusal built on it
# --------------------------------------------------------------------------- #
# MEASURED 2026-09-20, darwin 24.6.0 / CPython 3.9 / pandas 2.3.3 / numpy
# 2.0.2, peak RSS rise of one cold `deconvolve_sample(sample_id)` per process
# (``resource.getrusage(RUSAGE_SELF).ru_maxrss``, which is BYTES on darwin and
# KiB on Linux). Three cold runs per section before, five or six after; the
# MAXIMUM observed is recorded. Full table and method in
# results/scref_signatures_2026-09-20.md.
#
# THESE NUMBERS WERE RE-MEASURED AND THEY MOVED, so the reasoning attached to
# them moved too. The code they described densified a 17,828 x 3,589 reference
# submatrix to float64 twice in `build_signatures`, and a 4,764 x 17,828 spot
# matrix twice more here. All four are gone -- streamed in chunks,
# bit-identically -- and the peak fell by 2.2x to 3.3x. Before and after are
# both measured HERE, same harness, same host, same session, three cold runs
# before and five or six after, maximum of each:
#
#     section         | before  | after   | before/after
#     GBM1            | 2716.8  | 1227.9  | 2.2x
#     GBM2            | 3099.8  |  939.0  | 3.3x
#     GBM3            | 2567.3  |  922.2  | 2.8x
#     10x_parent_gbm  | 2923.5  | 1215.7  | 2.4x
#
# The "before" column is NOT the table in results/scref_memory_2026-09-20.md
# (GBM1 3,454.8 / GBM2 3,151.8 / GBM3 3,156.1 / 10x_parent_gbm 3,146.1). That
# file measured the same code and got 13-21 % more, and the run-to-run spread
# within each of my own triples is up to 1,400 MiB, so the two are consistent
# without being the same measurement. Do not average them and do not compute a
# ratio across the two tables; the ratios above are within-harness.
#
# WHAT THE DOMINANT TERM IS NOW, because the refusal message has to name it
# correctly and the old one named `build_signatures`. Staged, cold, GBM1:
#
#     imports                     109.7 MiB peak
#     load_sample(visium)         743.4     (+633.7)   <-- now the largest
#     reference_gene_symbols      772.7     ( +29.3)
#     build_signatures           1280.8     (+507.8 transient, 184 resident)
#     CSR column subset + NNLS   1280.8     (  +0.0)
#
# `load_sample` is now the biggest single term and it lives in
# `tmesim/data/real/visium.py`, which this module does not own. A Visium
# section still needs ~0.9 to 1.2 GB on this machine, so a 512 MB container
# still cannot serve one at any `dx_um` -- the gap is 1.7x to 2.4x rather than
# 6x, but it is a gap, and the refusal below must not imply it has closed.
#
# RE-MEASURED AGAIN 2026-10-10, same method, after `load_sample` was rebuilt in
# `tmesim/data/real/visium.py` to assemble its CSR out of mmread's COO instead
# of copying the matrix three times. That term fell 2.0x to 2.2x (GBM1 644.2 ->
# 299.3 MiB) and the output is bit-identical -- `deconvolve_sample`'s `coords`
# and `fractions` are `np.array_equal` to the pre-change values on all four
# sections, max|d| exactly 0. The whole-call peaks below fell with it, 8 to 12
# cold runs per section, maximum recorded:
#
#     section         | 2026-09-20 | 2026-10-10 | runs
#     GBM1            |    1227.9  |     895.6  |  8
#     GBM2            |     939.0  |     782.3  |  8
#     GBM3            |     922.2  |     787.4  | 12
#     10x_parent_gbm  |    1215.7  |     896.7  | 12
#
# The table had to move: 1227.9 against a real 895.6 is 332 MiB of refusal
# aimed at hosts that can now serve the section. The gap to 512 MB is now
# 1.53x to 1.75x, NOT the "1.7x to 2.4x" the paragraph above states -- that
# figure is the 2026-09-20 one and is superseded here, not contradicted.
#
# THE DOMINANT TERM MOVED AGAIN, and it is no longer in the Visium loader.
# Staged, cold, GBM1, after the change:
#
#     imports                     109.3 MiB peak
#     load_sample(visium)         398.3     (+289.0)
#     reference_gene_symbols      425.3     ( +27.0)
#     build_signatures            952.8     (+527.5)  <-- now the largest
#     CSR column subset + NNLS    952.8     (  +0.0)
#
# `build_signatures`' +527.5 MiB is `load_reference()` assembling the whole
# 23,460 x 3,589 int32 reference frame (321 MiB resident) plus its parse
# transient. It is a high-water rise here and was not one before only because
# the old loader left ~450 MiB of freed-but-retained allocator arenas behind
# for it to land in.
#
# AND 512 MB IS NOT REACHABLE. That is now measured rather than projected.
# results/visium_load_memory_2026-10-10.md runs the pipeline with the
# reference STREAMED -- never assembled, two chunked passes straight off the
# gzip, which is the cheapest shape this computation has -- and with the full
# spot x gene matrix dropped the moment the shared-gene submatrix exists.
# Absolute peak RSS of that floor: GBM1 686.6, GBM2 542.1, GBM3 531.5,
# 10x_parent_gbm 707.2 MiB. Every one is over 512, so no rewrite of
# `build_signatures` gets there either, and the refusal stays.
MEASURED_VISIUM_PEAK_MIB: Dict[str, float] = {
    "GBM1": 895.6,
    "GBM2": 782.3,
    "GBM3": 787.4,
    "10x_parent_gbm": 896.7,
}

#: The CHEAPEST Visium deconvolution ever measured here, across every section
#: and every run (GBM3, 565.8 MiB -- a single low run against a mode of 750 to
#: 790 for that same section; the run-to-run spread is the allocator, not the
#: workload). Used as the threshold for a sample with no measurement of its own
#: -- an ingested section, say. It is
#: deliberately the minimum and not the mean: this is the figure below which
#: NOTHING measured here has ever fitted, so refusing under it never refuses a
#: run that was observed to succeed. It is the one constant here where the
#: minimum is right and the maximum would be wrong, and it moved DOWN (869.9 ->
#: 565.8) for the same reason the table above did.
MEASURED_VISIUM_MIN_MIB = 565.8


class DeconvolutionTooLargeForHost(MemoryError):
    """A Visium deconvolution will not fit in this host's declared memory.

    ``MemoryError`` subclass for the same reason as
    ``scref_gbm.ReferenceTooLargeForHost``: it IS the memory error, raised
    before the allocation rather than by the OOM killer afterwards.

    NOT YET WIRED TO THE CLIENT. ``app/server/jobs.py`` renders a failed run's
    error verbatim only for classes named in ``RUN_REFUSALS``; this one is not
    in that set, so a visitor currently gets the redacted sentence plus a
    reference id while the text below goes to stderr and
    ``GET /api/runs/{id}/traceback``. The container is saved either way -- the
    refusal happens before the allocation -- but the reader is not served.
    Adding this name and ``ReferenceTooLargeForHost`` to ``RUN_REFUSALS``
    closes it; see results/scref_memory_2026-09-20.md §4.5.
    """


def check_visium_memory_budget(sample_id: str) -> None:
    """Refuse a Visium deconvolution this host cannot hold. See the table above.

    This is the refusal that ``TMESIM_MAX_GRID_CELLS`` cannot make. That
    ceiling counts VOXELS, and every megabyte measured above is spent before a
    grid exists -- on the Visium section itself and on the scRNA reference.
    Shrinking ``dx_um`` does not move any of it, so a visitor who picks a
    Visium section walks straight past the grid ceiling.

    Refuses only against a DECLARED limit (cgroup v2 -> v1 ->
    ``TMESIM_HOST_MEMORY_MB``, read by ``app/server/account.py::memory_limit``).
    Undeclared memory raises nothing.
    """
    h = memory_headroom()
    if h["override"] or h["headroom_bytes"] is None:
        return
    mib = 1024.0 * 1024.0
    need_mib = MEASURED_VISIUM_PEAK_MIB.get(sample_id, MEASURED_VISIUM_MIN_MIB)
    have_mib = h["headroom_bytes"] / mib
    if need_mib <= have_mib:
        return

    measured = (f"MEASURED at {need_mib:,.0f} MiB for this section"
                if sample_id in MEASURED_VISIUM_PEAK_MIB else
                f"not measured for this section; the cheapest Visium section "
                f"measured here needed {need_mib:,.0f} MiB")
    raise DeconvolutionTooLargeForHost(
        f"Refusing to deconvolve Visium sample {sample_id!r}: peak memory is "
        f"{measured}, and this host has {have_mib:,.0f} MiB of headroom "
        f"({h['limit_bytes'] / mib:,.0f} MiB limit, source {h['source']}, minus "
        f"{h['rss_bytes'] / mib:,.0f} MiB already resident)."
        f"\n\nThis is NOT a grid-size problem. TMESIM_MAX_GRID_CELLS counts "
        f"voxels and every byte above is spent before a grid exists, so no "
        f"choice of dx_um avoids it. The measured cost is dominated by reading "
        f"the single-cell reference panel the section is deconvolved against "
        f"-- a 23,460 x 3,589 counts matrix -- not by the requested resolution "
        f"and no longer by the Visium section itself."
        f"\n\nOptions: (1) choose a Xenium section -- the Xenium path does not "
        f"load this reference at all, so none of the cost above applies to it "
        f"(its own cost was not measured in this pass); "
        f"(2) run Visium on an instance with more memory; or (3) set "
        f"{_OVERRIDE_ENV}=1 to attempt it anyway and accept that the container "
        f"may be OOM-killed. Measurements: "
        f"results/visium_load_memory_2026-10-10.md."
    )

# ---------------------------------------------------------------------------
# Darmanis label -> coarse CellType mapping
# ---------------------------------------------------------------------------
# Keys are lower-cased label substrings so the (mis-spelled) "Astocyte" label in
# the GEO metadata is still caught.  Order of the canonical output classes:
_TYPE_ORDER: List[CellType] = [
    CellType.TUMOR,
    CellType.STROMAL,
    CellType.IMMUNE,
    CellType.ENDOTHELIAL,
]


def map_label(label: str) -> "CellType | None":
    """Map a Darmanis ``cell_type`` string onto a coarse core ``CellType``.

    Returns ``None`` for labels that carry no assignment (none currently, but
    keeps the collapse robust to unexpected vocabulary).
    """
    s = str(label).strip().lower()
    if "neoplastic" in s:
        return CellType.TUMOR
    if "immune" in s:
        return CellType.IMMUNE
    if "vascular" in s or "endothel" in s:
        return CellType.ENDOTHELIAL
    if any(k in s for k in ("astocyte", "astrocyte", "oligodendrocyte", "opc", "neuron")):
        return CellType.STROMAL
    return None


def _normalize_columns_to_unit_sum(mat: np.ndarray) -> np.ndarray:
    """Library-size normalize each COLUMN (cell) to sum 1 over the given genes.

    Cells with zero counts over the shared genes are left as all-zero.

    The division is written with ``out=``/``where=`` rather than as
    ``out[:, nz] = mat[:, nz] / colsums[nz]``. That expression built TWO extra
    whole-matrix temporaries -- the boolean-fancy-index copy ``mat[:, nz]`` and
    the quotient -- each the size of ``mat`` (488 MiB on the shipped
    reference), for a form that computes the identical numbers. Every element
    is still the single IEEE-754 double division ``mat[i, j] / colsums[j]``,
    which is correctly rounded and therefore bit-identical to what the slice
    form produced; ``where=nz`` leaves the zero-sum columns at the ``0.0`` that
    ``np.zeros_like`` already put there, which is what the slice form did too.
    """
    colsums = mat.sum(axis=0)
    out = np.zeros_like(mat, dtype=np.float64)
    nz = colsums > 0
    np.divide(mat, colsums, out=out, where=nz)
    return out


# --------------------------------------------------------------------------- #
# Signature memoisation
# --------------------------------------------------------------------------- #
# WHAT IS CACHED AND WHY IT IS THIS AND NOT THE REFERENCE.
#
# ``load_reference()`` peaks at a MEASURED 452 MiB and stays 322.7 MiB resident
# (results/scref_memory_2026-09-20.md). Memoising IT would pin that 322.7 MiB
# for the life of a process that deploys to a 512 MB container -- the fix that
# looks obvious and is wrong. What the Visium path actually wants out of that
# 452 MiB is this: a ``(G_shared, K)`` float64 matrix with K = 4, i.e. a few
# hundred kB. That is what is cached.
#
# KEY. ``(reference_fingerprint(), tuple(shared_genes))``. The fingerprint is
# (path, size, mtime_ns) of both cached files, so replacing the reference
# invalidates rather than serves a stale signature. ``shared_genes`` is part of
# the key because it HAS to be: ``_normalize_columns_to_unit_sum`` divides each
# cell by its total over the shared genes, so the signature for gene set A is
# not a row-slice of the signature for gene set B. A per-section cache is the
# largest correct unit, not a conservative choice.
#
# COPY-ON-READ. Callers get a fresh array every time. Handing out the cached
# object would make one caller's in-place edit everyone else's result, which is
# a correctness bug rather than an optimisation; ~640 kB a copy is not worth
# arguing about.
#
# A caller that supplies its own ``reference=`` (the held-out split in
# ``tmesim/validation/reconstruction_fidelity.py``) is NEVER cached: that frame
# is a subset chosen by the caller and the fingerprint cannot describe it.
_SIGNATURE_CACHE: "OrderedDict[Tuple, Tuple[np.ndarray, List[str]]]" = OrderedDict()

#: Enough for every built-in Visium section plus room; each entry is ~640 kB.
_SIGNATURE_CACHE_MAX = 8


def signature_cache_clear() -> None:
    """Drop every memoised signature. For tests and for an operator REPL."""
    _SIGNATURE_CACHE.clear()


def signature_cache_info() -> Dict:
    """``{'entries': n, 'max': m}`` -- what is currently memoised."""
    return {"entries": len(_SIGNATURE_CACHE), "max": _SIGNATURE_CACHE_MAX}


def build_signatures(
    shared_genes: Sequence[str],
    reference: "Tuple | None" = None,
) -> Tuple[np.ndarray, List[str]]:
    """Build per-cell-type mean-expression signatures over ``shared_genes``.

    Memoised when ``reference`` is omitted; see the block above for the key,
    why the reference frame itself is deliberately not the thing cached, and
    why the result is copied on the way out.

    Parameters
    ----------
    shared_genes : ordered gene symbols. The returned signature rows are aligned
        to this exact order. Symbols not found in the reference contribute an
        all-zero row (they carry no discriminative information).
    reference : optional pre-loaded ``(expr_df, labels_df, gene_names)`` tuple as
        returned by :func:`tmesim.data.real.scref_gbm.load_reference`; loaded
        on demand if omitted.

    Returns
    -------
    signatures : (n_shared_genes, K) float array, one column per cell type.
        Each column is the mean of the library-size-normalized single-cell
        profiles of that class over ``shared_genes`` (columns are on a
        comparable scale, so NNLS weights read as fractions).
    type_names : list[str] of length K — the ``CellType`` names, in the order of
        the signature columns. Only classes actually present in the reference
        are included.
    """
    if reference is not None:
        return _build_signatures_uncached(shared_genes, reference)

    key = (reference_fingerprint(), tuple(str(g) for g in shared_genes))
    hit = _SIGNATURE_CACHE.get(key)
    if hit is None:
        # The reference is loaded, consumed and dropped inside this call, so
        # its 452 MiB peak is not held alongside anything else.
        sig, names = _build_signatures_uncached(shared_genes, load_reference())
        _SIGNATURE_CACHE[key] = (sig, names)
        while len(_SIGNATURE_CACHE) > _SIGNATURE_CACHE_MAX:
            _SIGNATURE_CACHE.popitem(last=False)
        hit = (sig, names)
    else:
        _SIGNATURE_CACHE.move_to_end(key)
    sig, names = hit
    return sig.copy(), list(names)


def _build_signatures_uncached(
    shared_genes: Sequence[str],
    reference: Tuple,
) -> Tuple[np.ndarray, List[str]]:
    """:func:`build_signatures` with the memo removed. Unchanged arithmetic.

    "Unchanged" is meant literally and is TESTED, not asserted: the streamed
    path below is byte-for-byte identical to the dense one it replaced
    (``np.array_equal`` on the signature matrix, and on ``deconvolve_sample``'s
    output end to end). See :func:`_signatures_streamed` for why that is not
    luck, and the dense fallback for what happens when the argument it rests on
    does not hold.
    """
    expr_df, labels_df, _ = reference

    # Collapse labels -> coarse CellType (drop unmapped cells).
    mapped = labels_df["cell_type"].map(lambda v: map_label(v))
    keep = mapped.notna()
    mapped = mapped[keep]
    cells_kept = labels_df.index[keep]

    # Restrict reference to shared genes, preserving requested row order. Genes
    # absent from the reference get an all-zero row.
    shared_genes = list(shared_genes)
    ref_index = expr_df.index
    # First occurrence wins for any duplicate symbol in the reference index.
    dup = ref_index.duplicated()
    present_set = set(ref_index[~dup])

    G = len(shared_genes)
    row_of = {g: i for i, g in enumerate(shared_genes)}
    found = [g for g in shared_genes if g in present_set]
    found_rows = np.array([row_of[g] for g in found], dtype=int)

    type_names: List[str] = []
    groups: List[np.ndarray] = []
    mapped_vals = mapped.to_numpy()
    for ct in _TYPE_ORDER:
        sel = mapped_vals == ct
        if not sel.any():
            continue
        # Positions into ``cells_kept``, i.e. into the selected column block.
        groups.append(np.flatnonzero(sel))
        type_names.append(ct.name)
    if not groups:
        return np.zeros((G, 0)), []

    gene_pos, cell_pos = _positional_selection(expr_df, found, cells_kept)
    signatures = None
    if gene_pos is not None:
        signatures = _signatures_streamed(
            expr_df, gene_pos, cell_pos, groups, found_rows, G
        )
    if signatures is None:
        signatures = _signatures_dense(
            expr_df, dup, cells_kept, found, groups, found_rows, G
        )
    return signatures, type_names


def _positional_selection(
    expr_df, found: List[str], cells_kept,
) -> "Tuple[np.ndarray | None, np.ndarray | None]":
    """Resolve the label selections to POSITIONS, or refuse to.

    The dense path expressed the same selection with pandas labels::

        ref_sub = expr_df.loc[~expr_df.index.duplicated()]   # whole-matrix copy
        ref_sub = ref_sub[cells_kept]                        # whole-matrix copy
        ref_sub.loc[found]                                   # and a third

    On the shipped reference the first two select EVERYTHING -- the index has
    zero duplicated symbols and every one of the 3,589 cells carries a mappable
    label -- so they were two 321 MiB copies of a 23,460 x 3,589 int32 frame
    that changed nothing (MEASURED: 999.6 then 1,320.9 MiB resident, up from
    678.3). Positions let the one selection that does narrow the matrix be
    applied while the float64 block is filled, instead of three times before.

    Returns ``(None, None)`` -- meaning "use the dense path" -- when the label
    selection is not unambiguously positional: duplicate column labels, where
    ``ref_sub[cells_kept]`` would fan a barcode out into several columns, or a
    requested cell that is not a column at all. Row labels cannot fail here
    (``found`` was built from this index), but the lookup is guarded anyway.
    """
    cols = expr_df.columns
    if not cols.is_unique:
        return None, None
    cell_pos = np.asarray(cols.get_indexer(cells_kept), dtype=np.intp)
    if cell_pos.size and int(cell_pos.min()) < 0:
        return None, None
    first: Dict[str, int] = {}
    for i, g in enumerate(expr_df.index):
        if g not in first:
            first[g] = i          # first occurrence wins, as ~duplicated() does
    try:
        gene_pos = np.fromiter(
            (first[g] for g in found), dtype=np.intp, count=len(found)
        )
    except KeyError:  # pragma: no cover - `found` is built from this index
        return None, None
    return gene_pos, cell_pos


#: Gene rows per streamed chunk. The chunk is the only thing the streamed build
#: holds beyond the (G, 4) result, so this constant IS the peak.
#:
#: MEASURED, not chosen: cold peak of ``build_signatures(shared_genes)`` in its
#: own process, four runs per value (the figure below is the MAXIMUM, and it
#: includes ``load_reference()``'s own ~610 MiB, which is the floor underneath
#: all of them):
#:
#:     chunk rows |   250 |   500 |  1000 |  2000 |  4000 |  8000 | dense
#:     peak (MiB) | 579.0 | 620.1 | 601.9 | 684.3 | 786.1 | 1064.9| 1612.2
#:     wall   (s) |  4.05 |  4.13 |  4.23 |  4.50 |  3.98 |  4.26 |   4.20
#:
#: 1,000 is the measured optimum that is not also the slowest. Every value in
#: that sweep, including 8,000 and the dense fallback, produced a
#: BIT-IDENTICAL signature matrix -- which is the point: the result does not
#: depend on the chunking, so this constant is a memory knob and nothing else.
#:
#: That is true because :func:`_chunk_bounds` keeps it true. The floor is TWO
#: rows, not one, and a lone trailing row is absorbed; see that function and
#: premise 2 of :func:`_signatures_streamed` for the measured reason.
_SIGNATURE_CHUNK_ROWS = 1000

#: Above this, a float64 cannot hold every integer, and the exactness argument
#: in :func:`_signatures_streamed` stops being true. Checked, not assumed.
_EXACT_INT_LIMIT = 2.0 ** 53


def _chunk(expr_df, gene_pos, cell_pos, a: int, b: int, identity_cells: bool):
    """Gene rows ``gene_pos[a:b]`` x columns ``cell_pos``, as a numpy array.

    ``DataFrame.take`` rather than ``expr_df.to_numpy()[gene_pos[a:b]]``: on
    the shipped single-block frame ``to_numpy()`` happens to be a free view,
    but on a frame with more than one block it silently materialises the WHOLE
    matrix, which is the allocation this function exists to avoid. ``take``
    costs one chunk whatever the block layout is.
    """
    ch = expr_df.take(gene_pos[a:b], axis=0).to_numpy()
    if not identity_cells:
        ch = ch[:, cell_pos]
    return ch


def _signatures_streamed(
    expr_df, gene_pos: np.ndarray, cell_pos: np.ndarray,
    groups: List[np.ndarray], found_rows: np.ndarray, G: int,
) -> "np.ndarray | None":
    """Per-type mean signatures without ever holding the dense float64 matrix.

    Returns ``None`` if the exactness argument below does not hold, in which
    case the caller runs the dense path instead.

    WHY THIS IS BIT-IDENTICAL, AND WHERE THAT STOPS BEING TRUE
    ----------------------------------------------------------
    Floating-point addition is not associative, so a reduction rewritten into a
    different ORDER is not entitled to the same bits. Two reorderings happen
    here and each needs its own argument; neither is "it was close enough".

    1. **The column sums are accumulated over gene chunks in int64, not summed
       once in float64.** This is exact rather than approximately equal. The
       reference holds non-negative integer counts, so every partial sum -- in
       any order -- is an integer, and an integer below 2**53 is represented
       EXACTLY by a float64. A sum of exactly-representable integers whose
       every partial result is also exactly representable is computed with zero
       rounding error, and a computation with no rounding error cannot depend
       on the order it was done in. So the chunked int64 total converts to the
       same float64 that ``mat.sum(axis=0)`` produced. MEASURED on the shipped
       reference: largest column sum 1,445,849, i.e. 2**53 / 1,445,849 = 6.2e9
       times of headroom, and ``np.array_equal`` against the dense float64
       column sums is True.

       The two premises are CHECKED each call, not assumed: integer dtype and
       no negative entry (without which partial sums are not bounded by the
       total), plus the 2**53 ceiling on the result. Any failure returns
       ``None`` and the dense path runs.

    2. **Each per-type mean is computed over a chunk of GENE ROWS at a time.**
       The reduction runs along axis 1 (cells) and row i of the result depends
       on row i alone, so splitting the rows leaves every individual reduction
       with the same values in the same order. What is NOT split is the
       reduction itself -- a chunked sum along axis 1 would be a genuine
       reordering and is not done.

       That argument is sound and it is STILL not sufficient, which is the
       sharpest thing measured in this pass. ``block[:, idx]`` is an
       advanced-index copy, and numpy 2.0.2 picks its memory ORDER from the
       shape: with two or more rows the copy comes back non-C-contiguous, and
       ``.mean(axis=1)`` then reduces along a strided axis sequentially; with
       exactly ONE row the copy is contiguous and the same call reduces
       pairwise instead. Same values, same order of operands, different
       summation tree. MEASURED on the full 17,828-gene panel, a one-row chunk
       against the shipped one: 45.7 % of signature entries differ, max
       absolute 2.4e-17, max relative 3.6e-15 (7 ulp) -- and it does NOT stop
       there, because 69.9 % of GBM3's 7,100 spot fractions then differ too,
       max absolute 2.3e-15, max relative 6.0e-13. Invisible to any tolerance,
       real in the bits, and it propagates. So a chunk of fewer than 2 rows is
       refused below rather than clamped, and the standing guard against the
       property silently changing under a numpy upgrade is the test that
       compares this function to :func:`_signatures_dense` on the real
       reference with ``np.array_equal`` -- not this paragraph.

    The element-wise division is not a reduction and needs no argument:
    ``value / colsum`` is one correctly-rounded IEEE-754 double division either
    way.
    """
    if not all(np.issubdtype(dt, np.integer) for dt in expr_df.dtypes):
        return None          # premise 1 needs integers
    n = int(gene_pos.size)
    identity_cells = bool(np.array_equal(cell_pos, np.arange(cell_pos.size)))
    bounds = _chunk_bounds(n, _SIGNATURE_CHUNK_ROWS)

    # Pass 1: column sums, exactly, in int64.
    colsums = np.zeros(int(cell_pos.size), dtype=np.int64)
    for a, b in bounds:
        ch = _chunk(expr_df, gene_pos, cell_pos, a, b, identity_cells)
        if ch.size and int(ch.min()) < 0:
            return None      # premise 1 needs non-negative counts
        colsums += ch.sum(axis=0, dtype=np.int64)
        del ch
    cs = colsums.astype(np.float64)
    if cs.size and float(cs.max()) >= _EXACT_INT_LIMIT:
        return None          # premise 1 needs every total below 2**53
    nz = cs > 0

    # Pass 2: normalize each chunk in place and reduce it into the signatures.
    signatures = np.zeros((G, len(groups)), dtype=np.float64)
    for a, b in bounds:
        block = _chunk(expr_df, gene_pos, cell_pos, a, b, identity_cells)
        block = block.astype(np.float64)
        np.divide(block, cs, out=block, where=nz)
        if not nz.all():
            block[:, ~nz] = 0.0   # zero-count cells stay all-zero, as before
        rows = found_rows[a:b]
        for k, idx in enumerate(groups):
            signatures[rows, k] = block[:, idx].mean(axis=1)
        del block
    return signatures


def _chunk_bounds(n: int, step: int) -> List[Tuple[int, int]]:
    """Row-chunk ``(start, stop)`` pairs. NO CHUNK MAY HAVE EXACTLY ONE ROW.

    That invariant is the whole function. Premise 2 of
    :func:`_signatures_streamed` holds only while every chunk provokes the same
    memory-order choice out of ``block[:, idx]`` that the full matrix does, and
    numpy 2.0.2 makes a different choice for a one-row array than for any
    other. Two ways to end up with a lone row, both closed here:

    * ``step`` set below 2. Raised to 2 rather than refused: this is a MEMORY
      knob, and answering "you asked for the smallest possible chunk" with the
      dense 1.6 GB path would be the opposite of what was asked for.
    * ``n % step == 1``, which leaves a single trailing row -- the quiet one.
      With the shipped 1,000 and the shipped 17,828-gene panel the remainder is
      828 and nothing goes wrong, so a panel of 17,001 genes would have been
      the first anyone heard of it. The final chunk absorbs the orphan instead.

    ``n == 1`` is the one case where a one-row chunk is correct, because the
    dense path is then a one-row matrix too and makes the same choice.
    """
    if n <= 0:
        return []
    if n == 1:
        return [(0, 1)]
    step = max(int(step), 2)
    bounds: List[Tuple[int, int]] = []
    a = 0
    while a < n:
        b = min(a + step, n)
        if n - b == 1:
            b = n          # absorb the orphan rather than leave it alone
        bounds.append((a, b))
        a = b
    return bounds


def _signatures_dense(
    expr_df, dup: np.ndarray, cells_kept, found: List[str],
    groups: List[np.ndarray], found_rows: np.ndarray, G: int,
) -> np.ndarray:
    """The original whole-matrix build. Kept as the fallback, arithmetic verbatim.

    This is what ran before the streamed path existed and it is what still runs
    when :func:`_signatures_streamed` declines (see its docstring for the three
    conditions). It peaks at a MEASURED 2,503 MiB on the shipped reference
    against the streamed path's 748 MiB, which is the whole reason it is no
    longer the default -- but it is the definition of the answer, so it stays.
    """
    ref_sub = expr_df.loc[~dup]     # dedupe rows
    ref_sub = ref_sub[cells_kept]   # keep mapped cells
    ref_found = ref_sub.loc[found].to_numpy(dtype=np.float64)  # (n_found, n_cells)
    del ref_sub

    # Library-size normalize each cell over the shared genes (comparable columns).
    ref_norm = _normalize_columns_to_unit_sum(ref_found)
    del ref_found

    signatures = np.zeros((G, len(groups)), dtype=np.float64)
    for k, idx in enumerate(groups):
        signatures[found_rows, k] = ref_norm[:, idx].mean(axis=1)
    return signatures


def deconvolve_spot(spot_expr: np.ndarray, signatures: np.ndarray) -> np.ndarray:
    """Deconvolve one spot's expression into cell-type fractions via NNLS.

    Parameters
    ----------
    spot_expr : (G,) expression vector over the shared gene set (raw counts are
        fine — absolute scale does not affect the normalized fractions).
    signatures : (G, K) signature matrix from :func:`build_signatures`.

    Returns
    -------
    fractions : (K,) non-negative array summing to 1. If NNLS returns an
        all-zero solution (e.g. an empty spot), a uniform vector is returned.
    """
    b = np.asarray(spot_expr, dtype=np.float64).ravel()
    # Library-size normalize the spot to the signature scale (unit sum). This
    # does not change the normalized fractions but keeps NNLS well-conditioned.
    bsum = b.sum()
    if bsum > 0:
        b = b / bsum
    # scipy 1.13's nnls emits spurious FP warnings from its internal residual-norm
    # matmul even when the returned solution is correct; silence just those.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        w, _ = nnls(signatures, b)
    total = w.sum()
    if total <= 0:
        K = signatures.shape[1]
        return np.full(K, 1.0 / K) if K else w
    return w / total


def deconvolve_sample(
    sample_id: str,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Deconvolve every spot of a REAL Visium sample into cell-type fractions.

    Returns
    -------
    coords_um : (S, 2) spot centroids in micrometers (from the loader).
    fractions : (S, K) per-spot cell-type fractions (each row sums to 1).
    type_names : list[str] of length K, aligned to the fraction columns.

    Raises
    ------
    DeconvolutionTooLargeForHost
        Before loading anything, if this host declares a memory limit smaller
        than the measured peak of this section. See
        :func:`check_visium_memory_budget`.
    """
    # REFUSE FIRST. `load_sample` alone was measured at +400 to +620 MiB, so
    # the check has to precede it or it is checking after the fact.
    check_visium_memory_budget(sample_id)

    sample = load_sample(sample_id)
    coords_um = np.asarray(sample["coords_um"], dtype=np.float64)
    genes = list(sample["genes"])
    expr = sample["expression"].tocsr()  # (S, G_visium)

    # Shared gene set = Visium symbols also present in the reference.
    #
    # This used to be `load_reference()` -- a MEASURED +1,444 MiB parse of a
    # 23,460 x 3,589 counts matrix -- to build a set of gene symbols and
    # nothing else, and it then held that frame resident while
    # build_signatures made three more copies of it. `reference_gene_symbols()`
    # reads only the symbol column: +30 MiB, and verified array_equal to the
    # index the full parse returns.
    ref_syms = set(reference_gene_symbols())
    # Dedupe Visium symbols (first occurrence) to a well-defined column mapping.
    seen = set()
    shared_genes: List[str] = []
    shared_cols: List[int] = []
    for j, g in enumerate(genes):
        if g in ref_syms and g not in seen:
            seen.add(g)
            shared_genes.append(g)
            shared_cols.append(j)
    if not shared_genes:
        raise ValueError(f"sample {sample_id!r}: no genes shared with reference")

    # No `reference=` argument: build_signatures loads, consumes and frees the
    # reference inside the call (and memoises the small result), instead of
    # this frame staying alive across it.
    signatures, type_names = build_signatures(shared_genes)

    # Per-spot NNLS, one spot's dense vector at a time.
    #
    # THE SECOND DOUBLE DENSIFICATION, and it was the larger one. This read
    #
    #     expr[:, cols].toarray().astype(np.float64)
    #
    # which builds the whole (S, G_shared) matrix as int64 and then AGAIN as
    # float64: MEASURED 648 MiB each, 1,296 MiB together, for GBM1's 4,764
    # spots x 17,828 shared genes -- and the loop below then reads one row of
    # it at a time. The CSR it came from holds 15.2 M nonzeros in a 174 M-entry
    # matrix, so ~91 % of those bytes were zeros with nothing to say.
    #
    # Scattering one row into a reused buffer is bit-identical rather than
    # merely close: `toarray()` writes each stored value into a zero matrix and
    # so does this, and the int64 -> float64 conversion is the same
    # element-wise cast whether it happens on the dense matrix or on the value
    # being written. `sum_duplicates()` is what makes that true -- an
    # un-canonicalised CSR row can list an index twice, which `toarray()` SUMS
    # and a scatter assignment would silently resolve to the last write. It is
    # a no-op on a canonical matrix.
    expr_sub = expr[:, np.asarray(shared_cols, dtype=int)].tocsr()
    S = expr_sub.shape[0]
    K = signatures.shape[1]
    fractions = np.zeros((S, K), dtype=np.float64)
    for i, spot in enumerate(iter_spot_vectors(expr_sub, len(shared_genes))):
        fractions[i] = deconvolve_spot(spot, signatures)

    return coords_um, fractions, type_names


def iter_spot_vectors(csr, n_genes: int):
    """Yield each CSR row as a dense float64 vector, reusing one buffer.

    Stands in for ``csr.toarray().astype(np.float64)`` and yields exactly its
    rows. Factored out of :func:`deconvolve_sample` so that the test which
    pins it against ``toarray`` exercises THIS code rather than a copy of it
    that could drift -- the failure mode where a check passes against both the
    fixed and the broken version.

    THE BUFFER IS REUSED AND IS NOT YOURS. Each yielded array is the same
    object, cleared and refilled; a caller that wants to keep one must copy it.
    ``deconvolve_spot`` consumes it immediately and does not retain it.
    """
    csr.sum_duplicates()   # toarray() SUMS repeated entries; a scatter would
                           # keep only the last. No-op when already canonical.
    indptr, indices, data = csr.indptr, csr.indices, csr.data
    buf = np.zeros(int(n_genes), dtype=np.float64)
    for i in range(csr.shape[0]):
        lo, hi = int(indptr[i]), int(indptr[i + 1])
        cols_i = indices[lo:hi]
        buf[cols_i] = data[lo:hi]
        yield buf
        buf[cols_i] = 0.0   # reset only what was touched: O(nnz), not O(G)


if __name__ == "__main__":
    sid = "10x_parent_gbm"
    coords, fr, names = deconvolve_sample(sid)
    print(f"[{sid}] spots={fr.shape[0]}  types={names}")
    print("mean cell-type fractions across spots:")
    for name, m in zip(names, fr.mean(axis=0)):
        print(f"  {name:12s} {m:.3f}")
