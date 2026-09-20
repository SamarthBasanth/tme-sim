"""The streamed signature build must not change one bit of the answer.

WHY THIS FILE EXISTS
--------------------
``build_signatures`` used to densify a 17,828 x 3,589 submatrix to float64
twice and ``deconvolve_sample`` densified a 4,764 x 17,828 spot matrix twice
more. Both are now streamed in chunks, which took the MEASURED cold peak of a
Visium deconvolution from 2,567-3,100 MiB to 884-1,225 MiB.

That rewrite reorders arithmetic, and floating-point addition is not
associative, so a reordered reduction is NOT entitled to the same bits. One of
the two reorderings here is exact by a proof (integers below 2**53); the other
is exact only CONDITIONALLY, and this pass found the condition by tripping over
it -- a chunk of exactly one row silently changes the answer, and the change
reaches the published fractions. See ``_signatures_streamed`` and
``_chunk_bounds``. A proof nobody can re-run is worth nothing, so these tests
re-run it against the real GSE84465 reference: streamed output vs the dense
definition, ``np.array_equal``, no tolerance.

EVERY TEST HERE WAS CONFIRMED RED by planting the specific defect it targets,
before being confirmed green. The plant and the OBSERVED failure are named in
each docstring, because a check that passes against both the fixed and the
broken code is worse than no check -- this project has already shipped eight of
those, and writing this file produced two more before the plants caught them
(both are recorded, in ``test_spot_scatter_...`` and
``test_the_measured_table_...``).

The measurement-backed tests need the real 21 MB GSE84465 cache and skip
without it. The pure-numpy ones do not.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from tmesim.data.real import scref_gbm
from tmesim.reconstruct import deconvolve

_CACHE = scref_gbm._DEFAULT_CACHE
_COUNTS = os.path.join(_CACHE, scref_gbm._COUNTS_FILE)
_META = os.path.join(_CACHE, scref_gbm._META_FILE)
_HAVE_CACHE = os.path.exists(_COUNTS) and os.path.exists(_META)
needs_cache = pytest.mark.skipif(
    not _HAVE_CACHE, reason=f"real GSE84465 cache not present at {_CACHE}"
)


@pytest.fixture(autouse=True)
def _isolate(monkeypatch):
    monkeypatch.delenv("TMESIM_HOST_MEMORY_MB", raising=False)
    monkeypatch.setenv(scref_gbm._OVERRIDE_ENV, "1")
    deconvolve.signature_cache_clear()
    yield
    deconvolve.signature_cache_clear()


@pytest.fixture(scope="module")
def reference():
    if not _HAVE_CACHE:
        pytest.skip("real GSE84465 cache not present")
    os.environ[scref_gbm._OVERRIDE_ENV] = "1"
    expr, labels, _ = scref_gbm.load_reference()
    return expr, labels


def _force_dense(monkeypatch):
    """Make ``_build_signatures_uncached`` take the dense fallback."""
    monkeypatch.setattr(
        deconvolve, "_positional_selection", lambda *a, **k: (None, None)
    )


def _forbid_dense(monkeypatch):
    """Fail loudly if the dense fallback is taken when it should not be."""
    def boom(*a, **k):
        raise AssertionError("dense fallback taken; the streamed path was expected")
    monkeypatch.setattr(deconvolve, "_signatures_dense", boom)


# --------------------------------------------------------------------------- #
# 1. The headline: streamed == dense, on the real panel, bit for bit
# --------------------------------------------------------------------------- #
@needs_cache
def test_streamed_signatures_are_bit_identical_to_the_dense_definition(
        reference, monkeypatch):
    """The whole justification for the rewrite, on the real reference.

    PLANT that made this RED: replace the element-wise division

        np.divide(block, cs, out=block, where=nz)

    with the reciprocal-multiply form ``block *= 1.0 / cs``, which any reviewer
    would wave through as the same thing and which is faster. It is not the
    same thing -- ``x * (1/y)`` rounds twice where ``x / y`` rounds once -- and
    ``np.array_equal`` went False. Observed:
    ``assert False + where False = np.array_equal(...)`` at line 111.

    NEGATIVE CONTROL, recorded because it says what the exactness rests on:
    accumulating the column sums in chunked FLOAT64 instead of int64 leaves
    this test GREEN. That plant was tried first, expecting red. The reason it
    passes is the same reason the int64 version is allowed -- the values are
    integers below 2**53, so float64 addition of them carries no rounding error
    and is order-independent either way. The int64 accumulator is chosen for
    the guard it makes checkable, not because float64 would have been wrong.
    """
    expr, labels = reference
    genes = list(expr.index[~expr.index.duplicated()])

    _forbid_dense(monkeypatch)
    streamed, names_s = deconvolve._build_signatures_uncached(
        genes, (expr, labels, None)
    )
    monkeypatch.undo()
    _force_dense(monkeypatch)
    dense, names_d = deconvolve._build_signatures_uncached(
        genes, (expr, labels, None)
    )

    assert names_s == names_d == ["TUMOR", "STROMAL", "IMMUNE", "ENDOTHELIAL"]
    assert streamed.shape == (len(genes), 4)
    # Not allclose. Bits.
    assert np.array_equal(streamed, dense)
    # And the thing being compared is not trivially zero.
    assert np.count_nonzero(streamed) > 10_000


@needs_cache
def test_chunking_is_a_memory_knob_and_not_a_numerics_knob(reference, monkeypatch):
    """The same signatures at every chunk size, awkward remainders included.

    The sizes below are chosen against the specific hazard this pass found, not
    for roundness. 1 exercises the floor; 1499 over 1500 genes leaves a single
    trailing row; 7 and 512 are ordinary. All must agree, bit for bit.

    TWO PLANTS, both RED:

    * lower the ``_chunk_bounds`` floor from 2 rows to 1. ``chunk=1`` then
      diverged from every other size -- ``assert np.array_equal(sig, out[0])``
      False at line 134.
    * delete the orphan-row absorption (``if n - b == 1: b = n``). The
      three-gene block at the end went red: ``assert np.array_equal(two[0],
      three[0])`` False.

    The second plant is the reason that three-gene block exists. With only the
    1,500-gene panel the SAME plant left this test GREEN, because the orphaned
    row there happens to be one of the ~18 % whose value does not move. A check
    that depends on which row lands last is not a check.
    """
    expr, labels = reference
    genes = list(expr.index[:1500])
    _forbid_dense(monkeypatch)

    out = []
    for chunk in (1, 2, 7, 512, 749, 1499, 4000):
        monkeypatch.setattr(deconvolve, "_SIGNATURE_CHUNK_ROWS", chunk)
        sig, _ = deconvolve._build_signatures_uncached(genes, (expr, labels, None))
        out.append(sig)
    for sig in out[1:]:
        assert np.array_equal(sig, out[0])
    assert np.count_nonzero(out[0]) > 1000

    # A remainder of exactly one row, on three genes where that row's value is
    # KNOWN to move if it is left alone (the 1,500-gene panel above happens to
    # end on a row that does not, which is why it cannot carry this check).
    small = list(expr.index[:3])
    two, three = [], []
    for chunk, sink in ((2, two), (3, three)):
        monkeypatch.setattr(deconvolve, "_SIGNATURE_CHUNK_ROWS", chunk)
        sig, _ = deconvolve._build_signatures_uncached(small, (expr, labels, None))
        sink.append(sig)
    assert np.array_equal(two[0], three[0])
    assert np.count_nonzero(two[0]) > 0


def test_no_chunk_may_ever_contain_exactly_one_row():
    """The invariant ``_chunk_bounds`` exists to hold, over every remainder.

    A one-row chunk is not a rounding curiosity. MEASURED on the full
    17,828-gene panel, forcing one: 45.7 % of signature entries change (max
    relative 3.6e-15) and it propagates -- 69.9 % of GBM3's 7,100 spot
    fractions change, max relative 6.0e-13. No tolerance would ever see it.

    PLANT that made this RED: delete the orphan-row absorption. Observed:
    ``AssertionError: lone row at n=3 step=0: [(0, 2), (2, 3)]``.
    """
    for n in range(0, 60):
        for step in (0, 1, 2, 3, 7, 16, 59, 60, 61):
            bounds = deconvolve._chunk_bounds(n, step)
            assert bool(bounds) == (n > 0)
            # exact cover, in order, no gaps or overlaps
            if bounds:
                assert [a for a, _ in bounds] == [0] + [b for _, b in bounds[:-1]]
                assert bounds[-1][1] == n
            for a, b in bounds:
                assert b > a
                if n > 1:
                    assert b - a >= 2, f"lone row at n={n} step={step}: {bounds}"
    # n == 1 is the one legitimate one-row chunk: the dense path is one row too.
    assert deconvolve._chunk_bounds(1, 1000) == [(0, 1)]
    assert deconvolve._chunk_bounds(0, 1000) == []


@needs_cache
def test_a_lone_trailing_row_would_change_the_answer_if_it_were_allowed(
        reference, monkeypatch):
    """Proof that the invariant above is load-bearing and not decoration.

    Reaches past ``_chunk_bounds`` to produce the forbidden shape directly, and
    shows the signatures move. If numpy ever stopped distinguishing a one-row
    advanced-index copy from a multi-row one, this test would go red and
    ``_chunk_bounds`` could be simplified -- which is the point of pinning it.

    No plant is needed in the product: this test IS the plant, reaching past
    ``_chunk_bounds`` to build the forbidden shape and asserting the damage.
    Its counterpart is the three-gene block in ``test_chunking_...``, which
    asserts the shipped path does NOT produce that shape.
    """
    expr, labels = reference
    genes = list(expr.index[:301])
    _forbid_dense(monkeypatch)

    monkeypatch.setattr(deconvolve, "_SIGNATURE_CHUNK_ROWS", 100)
    good, _ = deconvolve._build_signatures_uncached(genes, (expr, labels, None))

    # 300 + 1: the orphan the real code absorbs.
    monkeypatch.setattr(
        deconvolve, "_chunk_bounds",
        lambda n, step: [(0, 100), (100, 200), (200, 300), (300, 301)],
    )
    orphaned, _ = deconvolve._build_signatures_uncached(genes, (expr, labels, None))

    assert np.array_equal(good[:300], orphaned[:300]), "only the lone row may move"
    assert not np.array_equal(good, orphaned), (
        "a one-row chunk no longer changes the result; _chunk_bounds may be "
        "simplifiable -- re-measure before deleting it"
    )


# --------------------------------------------------------------------------- #
# 2. The premise the exactness argument rests on, measured on the real file
# --------------------------------------------------------------------------- #
@needs_cache
def test_int64_column_sums_are_the_same_float64_as_the_dense_pairwise_sum(reference):
    """Premise 1 of ``_signatures_streamed``, checked rather than believed.

    Every column total is an integer far below 2**53, so the float64 sum of the
    dense matrix carries no rounding error and cannot depend on the order it
    was computed in. That is what lets the chunked int64 accumulation stand in
    for it.

    PLANT that made this RED: none was needed for the equality itself -- it is
    exact, and no reordering of it is not. What this test guards is the
    HEADROOM claim, so the plant was to lower ``_EXACT_INT_LIMIT`` to 2**20.
    Observed: ``assert (1048576.0 / 487894.0) > 1000000.0`` -- i.e. the
    assertion is reached and it is about the real column sums of the real file
    (487,894 over this 6,000-gene slice; 1,445,849 over the full 17,828-gene
    panel the product actually uses).
    """
    genes = list(reference[0].index[:6000])
    sub = reference[0].loc[genes]
    dense_f8 = sub.to_numpy(dtype=np.float64)
    colsums_dense = dense_f8.sum(axis=0)

    acc = np.zeros(sub.shape[1], dtype=np.int64)
    for a in range(0, len(genes), 997):   # a deliberately ragged chunk
        acc += sub.iloc[a:a + 997].to_numpy().sum(axis=0, dtype=np.int64)

    assert np.array_equal(acc.astype(np.float64), colsums_dense)
    assert float(acc.max()) < deconvolve._EXACT_INT_LIMIT
    # The headroom is not marginal, and the number is real, not a placeholder.
    assert acc.max() > 0
    assert deconvolve._EXACT_INT_LIMIT / float(acc.max()) > 1e6


# --------------------------------------------------------------------------- #
# 3. When the premises fail, the streamed path must decline, not approximate
# --------------------------------------------------------------------------- #
def _toy():
    """A tiny integer reference with two of each mappable class."""
    rng = np.random.default_rng(3)
    genes = [f"G{i}" for i in range(40)]
    cells = [f"c{j}" for j in range(8)]
    counts = rng.integers(0, 50, size=(40, 8)).astype(np.int32)
    expr = pd.DataFrame(counts, index=pd.Index(genes), columns=cells)
    labels = pd.DataFrame(
        {"cell_type": ["Neoplastic", "Neoplastic", "Immune cell", "Immune cell",
                       "Vascular", "Vascular", "Astrocyte", "Neuron"]},
        index=cells,
    )
    return expr, labels, genes


@pytest.mark.parametrize("break_it", ["float_dtype", "negative_value", "dup_columns"])
def test_streamed_path_declines_when_its_premises_do_not_hold(break_it, monkeypatch):
    """Three ways to invalidate the exactness argument; all must fall back.

    A silent approximation here would be the worst outcome of the whole change,
    so the fallback is asserted to be TAKEN -- not merely that the answer looks
    right.

    PLANT that made this RED: delete the guard for the case under test. Each
    deletion failed its own parametrisation and no other:

    * ``np.issubdtype(..., np.integer)`` removed ->
      ``AssertionError: float_dtype: streamed path ran on an invalid premise``
    * ``ch.min() < 0`` removed ->
      ``AssertionError: negative_value: streamed path ran on an invalid premise``
    * ``cols.is_unique`` removed in ``_positional_selection`` ->
      ``pandas.errors.InvalidIndexError: Reindexing only valid with uniquely
      valued Index objects``, raised out of ``get_indexer``. That one does not
      even reach the assertion: without the guard the function cannot compute
      an answer at all, which is the strongest form of "the guard is load-
      bearing" available.
    """
    expr, labels, genes = _toy()
    if break_it == "float_dtype":
        expr = expr.astype(np.float64)
    elif break_it == "negative_value":
        expr = expr.copy()
        expr.iloc[0, 0] = -5
    elif break_it == "dup_columns":
        expr = pd.concat([expr, expr.iloc[:, :2]], axis=1)
        labels = pd.concat([labels, labels.iloc[:2]], axis=0)

    took_dense = {"yes": False}
    real_dense = deconvolve._signatures_dense

    def spy(*a, **k):
        took_dense["yes"] = True
        return real_dense(*a, **k)
    monkeypatch.setattr(deconvolve, "_signatures_dense", spy)

    sig, names = deconvolve._build_signatures_uncached(genes, (expr, labels, None))
    assert took_dense["yes"], f"{break_it}: streamed path ran on an invalid premise"

    monkeypatch.undo()
    _force_dense(monkeypatch)
    ref_sig, ref_names = deconvolve._build_signatures_uncached(
        genes, (expr, labels, None)
    )
    assert names == ref_names
    assert np.array_equal(sig, ref_sig)
    assert np.count_nonzero(sig) > 0


def test_a_valid_toy_reference_takes_the_streamed_path(monkeypatch):
    """The counterpart: with the premises intact, the fallback must NOT run.

    Without this, the parametrised test above would still pass if the streamed
    path had been deleted outright and everything always fell back -- which is
    precisely the vacuous-green shape this project keeps finding.

    PLANT that made this RED: make ``_positional_selection`` return
    ``(None, None)`` unconditionally (i.e. pretend the streamed path can never
    apply). The ``_forbid_dense`` guard fired: AssertionError "dense fallback
    taken; the streamed path was expected".
    """
    expr, labels, genes = _toy()
    _forbid_dense(monkeypatch)
    sig, names = deconvolve._build_signatures_uncached(genes, (expr, labels, None))
    assert names == ["TUMOR", "STROMAL", "IMMUNE", "ENDOTHELIAL"]
    assert np.count_nonzero(sig) > 0


def test_first_occurrence_still_wins_for_a_duplicated_gene_symbol(monkeypatch):
    """Positional row lookup must reproduce ``~index.duplicated()``, not ``.loc``.

    ``ref_sub.loc[found]`` on a de-duplicated frame takes the FIRST row for each
    symbol. The positional mapping has to make the same choice; taking the last
    would be an easy and completely silent mistake, because the two agree on
    every frame that has no duplicates -- including the shipped reference.

    PLANT that made this RED: change ``if g not in first`` to an unconditional
    ``first[g] = i`` in ``_positional_selection`` (last occurrence wins). The
    doubled rows then read from the ``+7`` copy: ``assert np.array_equal(
    sig_dup, sig_ref)`` False, with all 40 gene rows moved, max absolute
    6.8e-3 and max relative 2.6 -- a 262 % error, not a rounding one.
    """
    expr, labels, genes = _toy()
    tail = expr.iloc[:6].copy()
    tail.iloc[:, :] = tail.to_numpy() + 7      # a visibly different second copy
    doubled = pd.concat([expr, tail], axis=0)  # 6 symbols now appear twice

    _forbid_dense(monkeypatch)
    sig_dup, _ = deconvolve._build_signatures_uncached(
        genes, (doubled, labels, None)
    )
    sig_ref, _ = deconvolve._build_signatures_uncached(genes, (expr, labels, None))
    assert np.array_equal(sig_dup, sig_ref)

    # ...and the +7 copy really would have changed the answer if it had won.
    monkeypatch.undo()
    _forbid_dense(monkeypatch)
    last_wins = pd.concat([tail, expr.iloc[6:]], axis=0)
    sig_other, _ = deconvolve._build_signatures_uncached(
        genes, (last_wins.reindex(genes), labels, None)
    )
    assert not np.array_equal(sig_other, sig_ref)


# --------------------------------------------------------------------------- #
# 4. The two smaller rewrites, each against the expression it replaced
# --------------------------------------------------------------------------- #
def test_normalize_columns_matches_the_slice_form_it_replaced():
    """``np.divide(out=, where=)`` vs ``out[:, nz] = mat[:, nz] / colsums[nz]``.

    PLANT that made this RED: drop ``where=nz`` from the ``np.divide`` call.
    Observed: ``RuntimeWarning: invalid value encountered in divide`` at the
    planted line, and ``array_equal`` False -- the zero-count columns became
    ``nan`` (0/0) instead of staying 0.0.
    """
    rng = np.random.default_rng(11)
    mat = rng.integers(0, 900, size=(300, 40)).astype(np.float64)
    mat[:, 5] = 0.0     # a cell with no counts over these genes
    mat[:, 19] = 0.0
    before = mat.copy()

    colsums = mat.sum(axis=0)
    nz = colsums > 0
    expected = np.zeros_like(mat, dtype=np.float64)
    expected[:, nz] = mat[:, nz] / colsums[nz]

    got = deconvolve._normalize_columns_to_unit_sum(mat)
    assert np.array_equal(got, expected)
    assert (got[:, 5] == 0).all() and (got[:, 19] == 0).all()
    assert np.array_equal(mat, before), "the input must not be mutated"


def test_spot_scatter_matches_toarray_including_duplicate_entries():
    """The per-spot buffer must reproduce ``toarray().astype(float64)`` exactly.

    The dangerous case is a non-canonical CSR: ``toarray`` SUMS repeated
    (row, col) entries and a scatter assignment would keep only the last.
    ``sum_duplicates()`` is what closes that, so the matrix here is built with
    a deliberate duplicate.

    PLANT that made this RED: remove the ``sum_duplicates()`` call from
    ``iter_spot_vectors``. Observed: row 1 came back ``[9., 0., 4., 0.]`` --
    the last write -- against ``toarray``'s ``[9., 0., 11., 0.]``.

    That plant passed against the FIRST version of this test, which built the
    matrix with ``coo_matrix(...).tocsr()``. The conversion canonicalises, so
    the duplicate never survived to the call and the test proved nothing. It is
    built from raw CSR arrays now, and asserts the duplicate is still there.
    """
    # Built from raw CSR arrays on purpose. Going through `coo.tocsr()` would
    # canonicalise on the way and the duplicate would never reach the code --
    # which is exactly how the first version of this test passed against a
    # build with `sum_duplicates()` deleted.
    indptr = np.array([0, 1, 4, 5])
    indices = np.array([3, 2, 2, 0, 1])   # row 1 lists column 2 twice: 7 + 4
    data = np.array([5, 7, 4, 9, 1], dtype=np.int64)
    m = sp.csr_matrix((data, indices, indptr), shape=(3, 4))
    assert not m.has_canonical_format, "the duplicate must survive to the call"
    dense = m.toarray().astype(np.float64)

    got = np.array([v.copy() for v in deconvolve.iter_spot_vectors(m, 4)])

    assert np.array_equal(got, dense)
    assert got.dtype == np.float64
    assert dense[1, 2] == 11.0, "the duplicate must be summed, not overwritten"


def test_the_spot_buffer_is_cleared_between_rows():
    """A stale entry would leak one spot's counts into the next one's NNLS.

    PLANT that made this RED: delete the ``buf[cols_i] = 0.0`` reset at the end
    of ``iter_spot_vectors``. Observed: spot 1 came back ``[0., 5., 7., 0.]``
    -- spot 0's count of 5 still sitting in the buffer.
    """
    m = sp.csr_matrix(np.array([[0, 5, 0, 0],
                                [0, 0, 7, 0],
                                [0, 0, 0, 0]], dtype=np.int64))
    got = [v.copy() for v in deconvolve.iter_spot_vectors(m, 4)]
    assert np.array_equal(got[0], [0., 5., 0., 0.])
    assert np.array_equal(got[1], [0., 0., 7., 0.])
    assert np.array_equal(got[2], [0., 0., 0., 0.])


# --------------------------------------------------------------------------- #
# 5. The refusal thresholds are measurements, and they still refuse 512 MB
# --------------------------------------------------------------------------- #
def test_the_measured_table_still_refuses_every_section_on_512_mb(monkeypatch):
    """The peak fell 2.2x-3.4x. It did NOT fall below 512 MB, and must not claim to.

    TWO PLANTS, both RED with ``Failed: DID NOT RAISE``: the whole
    ``MEASURED_VISIUM_PEAK_MIB`` table set to 100.0, and separately
    ``MEASURED_VISIUM_MIN_MIB`` set to 100.0 -- i.e. pretending the rewrite got
    a section under the live container's limit.

    A first attempt at 400.0 left this test GREEN and is recorded because it is
    instructive: 400 MiB still exceeds the headroom a 512 MB host has left once
    the test process is resident (~156 MiB), so the refusal fired anyway and
    the plant proved nothing. A threshold plant has to clear the headroom, not
    just the limit.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit",
        lambda: {"bytes": 512 * 1024 * 1024, "source": "test"},
    )
    monkeypatch.delenv(scref_gbm._OVERRIDE_ENV, raising=False)
    for sid in list(deconvolve.MEASURED_VISIUM_PEAK_MIB) + ["some_ingested_section"]:
        with pytest.raises(deconvolve.DeconvolutionTooLargeForHost) as e:
            deconvolve.check_visium_memory_budget(sid)
        assert "512" in str(e.value)
        # The message must not still blame build_signatures: that was the
        # dominant term before this rewrite and it is not any more.
        assert "build_signatures" not in str(e.value)

    # Every threshold is a real measurement, not a round placeholder.
    vals = list(deconvolve.MEASURED_VISIUM_PEAK_MIB.values())
    assert all(v != round(v, -2) for v in vals)
    assert deconvolve.MEASURED_VISIUM_MIN_MIB <= min(vals)


@needs_cache
def test_a_host_with_room_is_not_refused(monkeypatch):
    """The other half: the lowered thresholds must let a big host through.

    Without this, the test above would still pass if the refusal had been
    hard-wired to always raise.

    PLANT that made this RED: replace the ``if need_mib <= have_mib: return``
    early exit in ``check_visium_memory_budget`` with ``pass``. Observed:
    ``DeconvolutionTooLargeForHost: ... peak memory is MEASURED at 1,228 MiB
    for this section, and this host has 8,036 MiB of headroom``.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit",
        lambda: {"bytes": 8192 * 1024 * 1024, "source": "test"},
    )
    monkeypatch.delenv(scref_gbm._OVERRIDE_ENV, raising=False)
    for sid in deconvolve.MEASURED_VISIUM_PEAK_MIB:
        deconvolve.check_visium_memory_budget(sid)   # must not raise
