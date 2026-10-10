"""The Visium CSR must be assembled, not copied — and not one bit may move.

WHY THIS FILE EXISTS
--------------------
``visium.load_from_dirs`` used to turn ``mmread``'s output into the
``(n_spots, n_genes)`` CSR the pipeline wants in three full-matrix steps:
``.tocsr()``, a column gather that MEASURABLY selected every column, and a
``.T.tocsr()``. MEASURED on GBM1: +630.4 MiB over imports for a 174.5 MiB
result. ``_spots_by_genes_csr`` does the same turn by reading the COO's own
arrays as the CSR's ``data`` and ``indices`` and computing ``indptr`` with a
``searchsorted`` — +296.4 MiB, nothing copied.

That is only correct because a MEX file lists its entries column-major with
rows ascending inside each column, which in the TRANSPOSED orientation is
exactly CSR canonical order. It is a property of the files, not a guarantee of
the format, so the code CHECKS it and falls back to the real conversion when it
does not hold. These tests pin both halves: that the assembled matrix is
bit-identical to the conversion it replaces on the real sections, and that the
check actually catches the orderings it claims to.

The one-row trap recorded in results/scref_signatures_2026-09-20.md §3.3 does
not apply here — nothing in this path is a floating-point reduction; the values
are integers that are moved, not summed. What CAN move here is an index, and an
index that moves silently relabels a spot or a gene, so the comparisons below
are ``np.array_equal`` on the raw ``data`` / ``indices`` / ``indptr`` arrays
rather than on anything reduced.

EVERY TEST HERE WAS CONFIRMED RED by planting the specific defect it targets,
confirming the plant parses and is reached, then reverting. The plant and its
OBSERVED failure are named in each docstring.

The tests marked ``needs_cache`` need the real Visium sections and skip without
them. The rest are pure scipy/numpy and always run.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from tmesim.data.real import visium


_HAVE_SECTIONS = bool(visium.list_samples())
needs_cache = pytest.mark.skipif(
    not _HAVE_SECTIONS, reason="no Visium GBM section is cached locally"
)
# The smallest cached section, so the real-data checks stay cheap.
_SMALLEST = "GBM3" if "GBM3" in visium.list_samples() else (
    visium.list_samples()[0] if _HAVE_SECTIONS else None
)


def _old_way(coo, keep_idx):
    """The three-copy computation ``_spots_by_genes_csr`` replaced, verbatim.

    This is the DEFINITION of the answer for every comparison below. It is
    written out here rather than imported because the code it came from no
    longer exists; a test that compared the new path to itself would pass
    against any rewrite at all.
    """
    return coo.tocsr()[:, keep_idx].T.tocsr()


# --------------------------------------------------------------------------- #
# 1. The real sections: assembled == converted, bit for bit
# --------------------------------------------------------------------------- #
@needs_cache
def test_assembled_csr_is_bit_identical_to_the_conversion_it_replaced():
    """On a real section, the assembled CSR equals the old three-copy result.

    Not ``allclose``: ``array_equal`` on ``data``, ``indices`` and ``indptr``,
    plus the dtypes, because the dtype is part of what downstream reads.

    PLANT that made this RED: ``side="left"`` -> ``side="right"`` in the
    ``searchsorted`` that builds ``indptr``, which shifts every row boundary by
    the width of the first row. Observed: ``ValueError: index pointer should
    start with 0``, raised out of scipy's own ``check_format`` before the
    comparison was reached -- the constructor rejects that particular damage
    on its own, which is worth knowing and is not a reason to trust it with
    the rest.
    """
    spatial_dir, matrix_dir = visium._locate(_SMALLEST)
    coords_px = visium._read_positions(spatial_dir)
    barcodes_all = visium._read_barcodes(matrix_dir)
    _ens, symbols = visium._read_features(matrix_dir)
    keep_idx = np.asarray(
        [i for i, bc in enumerate(barcodes_all) if bc in coords_px], dtype=int
    )

    coo = visium._read_matrix_coo(matrix_dir)
    got = visium._spots_by_genes_csr(coo, 1, keep_idx, len(symbols))
    want = _old_way(coo, keep_idx)

    assert got.shape == want.shape
    assert np.array_equal(got.indptr, want.indptr), "indptr"
    assert np.array_equal(got.indices, want.indices), "indices"
    assert np.array_equal(got.data, want.data), "data"
    assert got.data.dtype == want.data.dtype
    assert got.indices.dtype == want.indices.dtype
    assert got.indptr.dtype == want.indptr.dtype


@needs_cache
def test_the_real_sections_take_the_assembled_path_and_select_every_spot():
    """Both claims the memory figure rests on, measured rather than assumed.

    If a real file were NOT in canonical order the fallback would run and the
    saving would be zero while every other test here stayed green; if the
    column selection were not the identity the gather would still copy. Both
    are properties of the cached files, so they are asserted against the files.

    PLANT that made this RED: ``_is_csr_canonical_order`` returning ``False``
    unconditionally. Observed:
    ``AssertionError: 10x_parent_gbm is not in canonical order``.
    """
    # SCOPED TO THE SECTIONS THE MEMORY FIGURE WAS MEASURED ON, not to
    # whatever happens to be registered when this runs.
    #
    # It used to iterate `visium.list_samples()`, and that failed in the full
    # suite on `gbmspace_at10_bra_5_fo_1_1 is not in canonical order` — for two
    # separate reasons, both worth recording:
    #
    # 1. `list_samples()` reads the module-level `_SAMPLES`, which
    #    `register_sample` MUTATES. Some other test in the suite registers a
    #    GBM-Space export and does not restore it, so this test's subject
    #    depended on what ran before it. An isolation defect, not this file's.
    #
    # 2. More usefully: that export is genuinely NOT in canonical order. So the
    #    fast path's precondition does not hold for every registrable sample,
    #    and `_is_csr_canonical_order`'s fallback is load-bearing IN PRODUCTION
    #    rather than only in principle. The memory saving reported in
    #    results/visium_load_memory_2026-10-10.md is therefore a property of
    #    these four sections and must not be generalised to any Visium input.
    #    `test_a_non_canonical_matrix_takes_the_fallback_and_is_still_correct`
    #    is what covers the other branch.
    measured = ("GBM1", "GBM2", "GBM3", "10x_parent_gbm")
    for sid in [s for s in measured if s in visium.list_samples()]:
        spatial_dir, matrix_dir = visium._locate(sid)
        coords_px = visium._read_positions(spatial_dir)
        barcodes_all = visium._read_barcodes(matrix_dir)
        keep_idx = np.asarray(
            [i for i, bc in enumerate(barcodes_all) if bc in coords_px], dtype=int
        )
        assert np.array_equal(keep_idx, np.arange(len(barcodes_all))), (
            f"{sid}: the column selection is no longer the identity"
        )
        coo = visium._read_matrix_coo(matrix_dir)
        assert visium._is_csr_canonical_order(coo.col, coo.row), (
            f"{sid} is not in canonical order"
        )
        del coo


@needs_cache
def test_the_assembly_does_not_copy_the_two_large_arrays():
    """The memory claim itself: ``data`` and ``indices`` are the COO's own.

    This is the whole point of the rewrite, and it is the one property that no
    value comparison can see — a version that copied both would stay
    bit-identical and cost the memory the change was made to stop paying.

    PLANT that made this RED: ``data = coo.data`` -> ``data = coo.data.copy()``
    in ``_spots_by_genes_csr``. Observed:
    ``AssertionError: data was copied``.
    """
    spatial_dir, matrix_dir = visium._locate(_SMALLEST)
    coords_px = visium._read_positions(spatial_dir)
    barcodes_all = visium._read_barcodes(matrix_dir)
    _ens, symbols = visium._read_features(matrix_dir)
    keep_idx = np.asarray(
        [i for i, bc in enumerate(barcodes_all) if bc in coords_px], dtype=int
    )
    coo = visium._read_matrix_coo(matrix_dir)
    data_in, rows_in = coo.data, coo.row
    expr = visium._spots_by_genes_csr(coo, 1, keep_idx, len(symbols))
    assert np.shares_memory(expr.data, data_in), "data was copied"
    assert np.shares_memory(expr.indices, rows_in), "indices were copied"


# --------------------------------------------------------------------------- #
# 2. The ordering check decides which path runs, so it has to be right
# --------------------------------------------------------------------------- #
def test_canonical_order_check_accepts_only_sorted_unique_pairs():
    """Sorted-and-unique is True; every way of not being so is False.

    The three failures are the three things the fallback exists for: a row that
    goes backwards, a column that goes backwards inside a row, and a repeated
    ``(row, col)`` pair — which ``tocsr`` would SUM and a direct assembly would
    leave as two entries.

    PLANT that made this RED: dropping the ``c[same_row] <= c_prev[same_row]``
    comparison (keeping only the row check). Observed:
    ``AssertionError: a column that goes backwards inside a row is not
    canonical`` -- the first of the two column cases; the duplicate-pair case
    on the next line goes with it.
    """
    r = np.array([0, 0, 1, 2, 2], dtype=np.int32)
    c = np.array([1, 4, 0, 2, 9], dtype=np.int32)
    assert visium._is_csr_canonical_order(r, c)

    assert not visium._is_csr_canonical_order(
        np.array([0, 1, 0], dtype=np.int32), np.array([0, 0, 1], dtype=np.int32)
    ), "a row that goes backwards is not canonical"
    assert not visium._is_csr_canonical_order(
        np.array([0, 0], dtype=np.int32), np.array([4, 1], dtype=np.int32)
    ), "a column that goes backwards inside a row is not canonical"
    assert not visium._is_csr_canonical_order(
        np.array([0, 0], dtype=np.int32), np.array([3, 3], dtype=np.int32)
    ), "a repeated (row, col) pair is not canonical"

    # Degenerate inputs: nothing to compare is trivially in order.
    assert visium._is_csr_canonical_order(
        np.zeros(0, np.int32), np.zeros(0, np.int32)
    )
    assert visium._is_csr_canonical_order(
        np.zeros(1, np.int32), np.zeros(1, np.int32)
    )


def test_canonical_order_check_sees_a_violation_that_straddles_a_chunk():
    """The check is chunked, so the pair it compares across a seam must count.

    This is the quiet failure mode of any windowed scan: the violation sits
    between the last element of one chunk and the first of the next, and a loop
    that compares only WITHIN each window never looks at it.

    PLANT that made this RED: windowing that compares only WITHIN each window
    -- ``for lo in range(0, n, chunk)`` with ``r, r_prev = rows[lo + 1:hi],
    rows[lo:hi - 1]`` (and the same for the columns) instead of reaching back
    to ``lo - 1``. Observed with chunk=2:
    ``AssertionError: a backwards row across a chunk seam must be caught``.

    THIS TEST WAS VACUOUS ON ITS FIRST ATTEMPT and the plant is what caught it.
    It used ``rows = [0, 1, 2, 1]``, where the violating pair is (index 2,
    index 3) -- which falls INSIDE a window under both the correct and the
    planted indexing, so the plant stayed green. The violating pair has to sit
    on the planted windowing's seam, i.e. at (index 1, index 2).
    """
    # chunk=2: pairs (0,1) and (2,3) are inside windows; (1,2) is the seam.
    rows = np.array([0, 2, 1, 3], dtype=np.int32)
    cols = np.array([0, 0, 0, 0], dtype=np.int32)
    assert not visium._is_csr_canonical_order(rows, cols, chunk=2), (
        "a backwards row across a chunk seam must be caught"
    )
    # Same for a repeated pair across the seam.
    rows = np.array([0, 1, 1, 2], dtype=np.int32)
    cols = np.array([0, 5, 5, 0], dtype=np.int32)
    assert not visium._is_csr_canonical_order(rows, cols, chunk=2), (
        "a repeated pair across a chunk seam must be caught"
    )
    # And a legitimate matrix must survive every chunk size, seams included.
    rows = np.array([0, 0, 1, 1, 2, 3, 3], dtype=np.int32)
    cols = np.array([0, 3, 1, 2, 5, 0, 4], dtype=np.int32)
    for chunk in (1, 2, 3, 7, 4096):
        assert visium._is_csr_canonical_order(rows, cols, chunk=chunk), chunk


# --------------------------------------------------------------------------- #
# 3. The fallback: a file that is NOT in canonical order
# --------------------------------------------------------------------------- #
def _toy_coo(shape, entries, dtype=np.int64):
    """A COO built from raw arrays, with NO canonicalisation on the way in.

    ``coo_matrix((data, (row, col)))`` stores what it is given; the duplicate
    and the out-of-order entry below both survive to the call, which is the
    whole point. (The analogous test in ``test_signature_streaming.py`` was
    vacuous on its first attempt because ``.tocsr()`` canonicalised the fixture
    before the code under test ever saw it.)
    """
    rows = np.asarray([e[0] for e in entries], dtype=np.int32)
    cols = np.asarray([e[1] for e in entries], dtype=np.int32)
    data = np.asarray([e[2] for e in entries], dtype=dtype)
    return sp.coo_matrix((data, (rows, cols)), shape=shape)


def test_a_non_canonical_matrix_falls_back_and_sums_its_duplicates():
    """Out of order and with a repeated entry: the fallback must do the work.

    ``tocsr()`` SUMS a repeated ``(row, col)``; a direct assembly would keep
    both entries and every later reader would see whichever it happened to
    scatter last. The assembled path must therefore not be taken here, and the
    answer must equal the conversion's.

    PLANT that made this RED: ``if _is_csr_canonical_order(...)`` -> ``if
    True``, forcing the assembled path. Observed: ``[[5, 11, 9], [0, 0, 0]]``
    against the conversion's ``[[5, 7, 0], [0, 4, 9]]`` -- every value in the
    wrong spot, 5 stored entries instead of 4, the duplicate left unsummed.
    """
    # genes x spots, deliberately NOT in column-major order, with (gene 1,
    # spot 0) listed twice.
    coo = _toy_coo(
        (3, 2),
        [(1, 1, 4), (0, 0, 5), (1, 0, 3), (2, 1, 9), (1, 0, 4)],
    )
    keep = np.arange(2)
    got = visium._spots_by_genes_csr(coo, 1, keep, 3)
    want = _old_way(coo, keep)
    assert np.array_equal(got.toarray(), want.toarray())
    assert np.array_equal(got.indptr, want.indptr), "indptr"
    assert np.array_equal(got.indices, want.indices), "indices"
    assert np.array_equal(got.data, want.data), "data"
    # The duplicate really was summed: spot 0, gene 1 = 3 + 4.
    assert got.toarray()[0, 1] == 7


# --------------------------------------------------------------------------- #
# 4. The spot selection, which NO cached section exercises
# --------------------------------------------------------------------------- #
def test_a_non_identity_spot_selection_matches_the_gather_it_replaced():
    """Dropping and reordering spots must still equal ``mat[:, keep].T``.

    Every cached section selects all of its spots in order, so this branch is
    dead on this machine's data and alive on anyone's RAW Space Ranger matrix,
    where ``barcodes.tsv`` carries ~1.5x more barcodes than there are in-tissue
    spots. Untested, it would fail for the first person to ingest one.

    PLANT that made this RED: ``pos[keep_idx] = np.arange(n_keep)`` ->
    ``pos[keep_idx] = keep_idx``, i.e. keeping the FILE's spot numbering
    instead of renumbering into the kept block. Observed on ``keep = [0, 2,
    4]``: ``[[1, 0, 2, 0], [0, 0, 0, 0], [4, 0, 0, 5]]`` against the gather's
    ``[[1, 0, 2, 0], [4, 0, 0, 5], [0, 0, 7, 8]]`` -- spot 4's entries dropped
    off the end of the renumbered block and spot 2's landed one row late. A
    silent relabelling of which spot is which, which is the worst shape this
    defect could take.
    """
    coo = _toy_coo(
        (4, 5),  # 4 genes x 5 spots
        [(0, 0, 1), (2, 0, 2), (1, 1, 3), (0, 2, 4), (3, 2, 5),
         (1, 3, 6), (2, 4, 7), (3, 4, 8)],
    )
    for keep in (
        np.array([0, 2, 4]),        # a proper subset
        np.array([1, 3]),           # a different subset
        np.array([4, 0]),           # out of order
        np.array([2]),              # a single spot
    ):
        got = visium._spots_by_genes_csr(coo, 1, keep, 4)
        want = _old_way(coo, keep)
        assert np.array_equal(got.toarray(), want.toarray()), keep
        assert np.array_equal(got.indptr, want.indptr), (keep, "indptr")
        assert np.array_equal(got.indices, want.indices), (keep, "indices")
        assert np.array_equal(got.data, want.data), (keep, "data")


def test_a_transposed_file_is_read_through_the_other_axis():
    """``spots_axis=0`` must give the same matrix as transposing the file did.

    The old reader handled a ``(spots, genes)`` encoding with a physical
    ``mat.T.tocsr()``; this one passes an axis number instead, and an axis
    number that is ignored would silently transpose the whole section.

    PLANT that made this RED: ``spot_of = coo.col if spots_axis == 1 else
    coo.row`` -> ``spot_of = coo.col`` (axis argument ignored). Observed:
    ``ValueError: axis 0 index 2 exceeds matrix dimension 2`` -- the gene index
    used as a spot index, out of bounds for the 2-spot result.
    """
    entries = [(0, 0, 1), (0, 2, 2), (1, 1, 3)]
    spots_by_genes = _toy_coo((2, 3), entries)          # spots x genes
    genes_by_spots = _toy_coo((3, 2), [(c, r, v) for r, c, v in entries])
    keep = np.arange(2)

    got = visium._spots_by_genes_csr(spots_by_genes, 0, keep, 3)
    want = _old_way(genes_by_spots, keep)
    assert np.array_equal(got.toarray(), want.toarray())


# --------------------------------------------------------------------------- #
# 5. End to end, on a real section
# --------------------------------------------------------------------------- #
@needs_cache
def test_load_sample_still_returns_a_canonical_int_csr():
    """The contract every downstream consumer reads, pinned.

    ``deconvolve.iter_spot_vectors`` scatters CSR rows into a dense buffer and
    is only equivalent to ``toarray()`` on a canonical matrix; ~30 call sites
    index ``sample["expression"]`` by gene column. Shape, orientation, dtype
    and canonical form are therefore part of the output, not an implementation
    detail.

    PLANT that made this RED: ``shape=(n_keep, n_genes)`` ->
    ``shape=(n_genes, n_keep)`` in the assembled branch. Observed:
    ``ValueError: index pointer size (1776) should be (36602)`` out of scipy's
    ``check_format`` -- a 1,775-spot section claiming 36,601 rows.
    """
    d = visium.load_sample(_SMALLEST)
    e = d["expression"]
    assert e.format == "csr"
    assert e.shape == (d["n_spots"], d["n_genes"])
    assert e.shape == (len(d["barcodes"]), len(d["genes"]))
    assert np.issubdtype(e.data.dtype, np.integer)
    assert e.has_canonical_format
    assert int(e.data.min()) >= 0
    # total_counts is the row sum of exactly this matrix.
    assert np.array_equal(
        d["total_counts"], np.asarray(e.sum(axis=1)).ravel().astype(np.float64)
    )
