"""The scRNA reference must not be able to OOM-kill the container.

WHY THIS FILE EXISTS
--------------------
``TMESIM_MAX_GRID_CELLS`` refuses a run whose GRID is too big for the host.
It cannot refuse this one. ``tmesim/data/real/scref_gbm.py::load_reference``
allocates as a function of the REFERENCE PANEL, not of the grid, so its cost is
identical at every ``dx_um``; a visitor who picks a Visium section walks past
the voxel ceiling without trying. MEASURED 2026-09-20 on this machine:

    load_reference()        peak +1,961 MiB before, +610 MiB after
    deconvolve_sample()     peak ~2,900-3,455 MiB, essentially unchanged

so on a 512 MB instance the Visium path cannot run at all and the only honest
answer is a refusal that says so. Every test below was confirmed RED against
the pre-fix code by planting the specific defect it targets; the plants are
named in each docstring so the next person can repeat the exercise rather than
trust it.

The measurement-backed tests need the real 21 MB GSE84465 cache and skip
without it. The refusal tests do NOT -- they drive the detector directly, so
the guard that protects the live site is checked even in a bare checkout.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from tmesim.data.real import scref_gbm
from tmesim.reconstruct import deconvolve

_CACHE = scref_gbm._DEFAULT_CACHE
_COUNTS = os.path.join(_CACHE, scref_gbm._COUNTS_FILE)
_META = os.path.join(_CACHE, scref_gbm._META_FILE)

_HAVE_CACHE = os.path.exists(_COUNTS) and os.path.exists(_META)
needs_cache = pytest.mark.skipif(
    not _HAVE_CACHE,
    reason=f"real GSE84465 cache not present at {_CACHE}",
)


@pytest.fixture(autouse=True)
def _isolate(monkeypatch):
    """No ambient host-size env, and no signatures memoised across tests."""
    monkeypatch.delenv("TMESIM_HOST_MEMORY_MB", raising=False)
    monkeypatch.delenv(scref_gbm._OVERRIDE_ENV, raising=False)
    deconvolve.signature_cache_clear()
    yield
    deconvolve.signature_cache_clear()


def _limit(monkeypatch, mb, source="test"):
    """Force a declared host limit without needing a cgroup."""
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit",
        lambda: {"bytes": int(mb * 1024 * 1024), "source": source},
    )


# --------------------------------------------------------------------------- #
# 1. The refusal exists, fires, and is not a warning
# --------------------------------------------------------------------------- #
@needs_cache
def test_load_reference_refuses_on_a_512_mb_host(monkeypatch):
    """The exact deployment that is live: a 512 MB container.

    PLANTED to confirm this is not vacuous: deleting the
    ``check_memory_budget(...)`` call from ``load_reference`` makes this FAIL
    with "DID NOT RAISE ReferenceTooLargeForHost" -- i.e. the pre-fix code,
    which really does go on to allocate, is red here.
    """
    _limit(monkeypatch, 512, source="cgroup v2 memory.max")
    with pytest.raises(scref_gbm.ReferenceTooLargeForHost):
        scref_gbm.load_reference()


@needs_cache
def test_refusal_names_the_measured_requirement_and_the_detected_limit(monkeypatch):
    """A refusal that does not say WHAT and HOW MUCH is a 500 with extra steps.

    PLANTED: replacing the message body with a bare "out of memory" makes every
    assertion below fail. Both numbers are checked against the objects they
    came from, not against a copy of the string, so an edit that changes the
    measured constant without changing the message is also red.
    """
    _limit(monkeypatch, 512, source="cgroup v2 memory.max")
    with pytest.raises(scref_gbm.ReferenceTooLargeForHost) as e:
        scref_gbm.load_reference()
    msg = str(e.value)

    # The measured requirement, as the module records it.
    assert f"{scref_gbm.MEASURED_LOAD_PEAK_MIB:,.0f} MiB" in msg
    # The detected limit AND where it was detected.
    assert "512 MiB" in msg
    assert "cgroup v2 memory.max" in msg
    # The thing a reader will otherwise get wrong.
    assert "TMESIM_MAX_GRID_CELLS" in msg
    # A way out.
    assert scref_gbm._OVERRIDE_ENV in msg


def test_no_refusal_when_the_host_declares_nothing(monkeypatch):
    """An undetected limit must never be treated as a small one.

    PLANTED: making ``check_memory_budget`` fall back to a default limit when
    none is declared turns this red immediately.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit",
        lambda: {"bytes": None, "source": "not declared"},
    )
    scref_gbm.check_memory_budget(23465, 3589)  # must not raise


def test_override_env_lets_an_operator_through(monkeypatch):
    """The escape hatch the refusal advertises has to actually exist."""
    _limit(monkeypatch, 64)
    with pytest.raises(scref_gbm.ReferenceTooLargeForHost):
        scref_gbm.check_memory_budget(23465, 3589)
    monkeypatch.setenv(scref_gbm._OVERRIDE_ENV, "1")
    scref_gbm.check_memory_budget(23465, 3589)  # must not raise


def test_requirement_scales_with_the_reference_and_is_not_a_constant(monkeypatch):
    """A bigger reference must ask for more, or the estimate is decoration.

    PLANTED: hardcoding ``estimate_load_bytes`` to return the measured 610 MiB
    regardless of shape makes the strict inequality fail.
    """
    small = scref_gbm.estimate_load_bytes(1_000, 100)
    shipped = scref_gbm.estimate_load_bytes(*scref_gbm.MEASURED_SHAPE)
    bigger = scref_gbm.estimate_load_bytes(
        scref_gbm.MEASURED_SHAPE[0] * 2, scref_gbm.MEASURED_SHAPE[1]
    )
    assert small < shipped < bigger
    # On the shipped reference it must reproduce the recorded measurement.
    assert shipped == pytest.approx(
        scref_gbm.MEASURED_LOAD_PEAK_MIB * 1024 * 1024, rel=1e-9
    )


# --------------------------------------------------------------------------- #
# 2. The Visium path refuses BEFORE it allocates
# --------------------------------------------------------------------------- #
def test_visium_refusal_happens_before_any_data_is_loaded(monkeypatch):
    """"Refuse before allocating" is the whole claim, so it is the assertion.

    ``load_sample`` alone was MEASURED at +400 to +620 MiB. If the budget check
    ran after it, the refusal would arrive having already spent half the
    container. Here ``load_sample`` is booby-trapped: reaching it fails the
    test rather than merely being slower.

    PLANTED: moving ``check_visium_memory_budget(sample_id)`` to AFTER the
    ``load_sample`` call makes this red with "load_sample was reached".
    """
    _limit(monkeypatch, 512, source="cgroup v2 memory.max")

    def _trap(*a, **k):
        raise AssertionError("load_sample was reached before the budget check")

    monkeypatch.setattr(deconvolve, "load_sample", _trap)
    with pytest.raises(deconvolve.DeconvolutionTooLargeForHost) as e:
        deconvolve.deconvolve_sample("GBM3")
    msg = str(e.value)
    assert f"{deconvolve.MEASURED_VISIUM_PEAK_MIB['GBM3']:,.0f} MiB" in msg
    assert "512 MiB" in msg
    assert "cgroup v2 memory.max" in msg
    assert "TMESIM_MAX_GRID_CELLS" in msg


def test_visium_refusal_covers_a_section_it_has_no_measurement_for(monkeypatch):
    """An ingested section is not in the table and must still be refused.

    PLANTED: making the lookup ``MEASURED_VISIUM_PEAK_MIB[sample_id]`` instead
    of ``.get(sample_id, MEASURED_VISIUM_MIN_MIB)`` turns this into a KeyError
    -- which is exactly the hole this asserts is closed.
    """
    _limit(monkeypatch, 512)
    with pytest.raises(deconvolve.DeconvolutionTooLargeForHost) as e:
        deconvolve.check_visium_memory_budget("an_ingested_section")
    assert "not measured for this section" in str(e.value)
    assert f"{deconvolve.MEASURED_VISIUM_MIN_MIB:,.0f} MiB" in str(e.value)


def test_visium_runs_when_the_host_is_big_enough(monkeypatch):
    """The refusal must not be a blanket ban.

    PLANTED: dropping the ``need_mib <= have_mib`` early return makes this red.
    """
    _limit(monkeypatch, 16_384)
    deconvolve.check_visium_memory_budget("GBM3")  # must not raise


# --------------------------------------------------------------------------- #
# 3. The cheap accessor really is the same data as the expensive one
# --------------------------------------------------------------------------- #
@needs_cache
def test_gene_symbols_are_identical_to_the_full_load(monkeypatch):
    """The Visium path swapped a +1,444 MiB parse for a +30 MiB one. Same answer?

    PLANTED two ways, both red: (a) dropping the ``_NON_GENE_ROWS`` filter from
    ``reference_gene_symbols`` leaves the 5 HTSeq summary rows in and the
    lengths differ by 5; (b) returning ``sorted(symbols)`` passes a set
    comparison and fails this one, which is why the assertion is
    ``array_equal`` on the ORDER and not ``set(...) == set(...)``.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit", lambda: {"bytes": None, "source": "test"}
    )
    cheap = scref_gbm.reference_gene_symbols()
    expr, _labels, gene_names = scref_gbm.load_reference()
    full = np.asarray(expr.index, dtype=object)
    assert np.array_equal(cheap, full)
    assert np.array_equal(cheap, gene_names)
    assert not set(cheap) & scref_gbm._NON_GENE_ROWS


@needs_cache
def test_counts_are_bit_exact_against_the_original_read(monkeypatch):
    """int32 + chunked must equal int64 + monolithic, value for value.

    The comparison is made against the ORIGINAL expression --
    ``pd.read_csv(path, sep=r"\\s+", index_col=0)`` -- on the first rows of the
    real file, not against a stored golden array, so it is still a real
    comparison if the cached reference is ever replaced.

    PLANTED: changing the parse dtype to ``np.int16`` wraps the counts above
    32,767 (the file's maximum is 77,122) and this goes red on
    ``array_equal``. Changing ``_PARSE_CHUNK_ROWS`` does NOT make it red, and
    should not: that is a performance knob, not a numerics one.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit", lambda: {"bytes": None, "source": "test"}
    )
    n = 800
    original = pd.read_csv(_COUNTS, sep=r"\s+", index_col=0, nrows=n)
    original.index = original.index.astype(str)

    expr, _labels, _g = scref_gbm.load_reference()
    assert expr.dtypes.iloc[0] == np.int32

    shared_rows = [g for g in original.index if g in set(expr.index)][:n]
    got = expr.loc[shared_rows, list(expr.columns)].to_numpy()
    want = original.loc[shared_rows, list(expr.columns)].to_numpy()
    assert got.shape == want.shape and got.size > 0
    assert np.array_equal(got.astype(np.int64), want)


# --------------------------------------------------------------------------- #
# 4. The cache is on the small object, and handing it out is safe
# --------------------------------------------------------------------------- #
@needs_cache
def test_signature_cache_hands_out_a_copy(monkeypatch):
    """A cached array a caller can mutate is a correctness bug, not a cache.

    PLANTED: returning ``hit`` directly instead of ``sig.copy(), list(names)``
    makes the second call see the poisoned value and this goes red.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit", lambda: {"bytes": None, "source": "test"}
    )
    genes = list(scref_gbm.reference_gene_symbols()[:300])
    first, names = deconvolve.build_signatures(genes)
    first[:] = -12345.0
    names.append("NOT_A_CELL_TYPE")

    second, names2 = deconvolve.build_signatures(genes)
    assert not np.any(second == -12345.0)
    assert "NOT_A_CELL_TYPE" not in names2


@needs_cache
def test_signature_cache_is_keyed_on_the_reference_files(monkeypatch):
    """Replace the reference and a stale signature must NOT be served.

    PLANTED: dropping ``reference_fingerprint()`` from the cache key makes the
    second call a cache hit, ``calls`` stays 1, and this goes red.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit", lambda: {"bytes": None, "source": "test"}
    )
    genes = list(scref_gbm.reference_gene_symbols()[:300])

    calls = {"n": 0}
    real_load = scref_gbm.load_reference

    def _counted(*a, **k):
        calls["n"] += 1
        return real_load(*a, **k)

    monkeypatch.setattr(deconvolve, "load_reference", _counted)

    deconvolve.build_signatures(genes)
    assert calls["n"] == 1
    deconvolve.build_signatures(genes)
    assert calls["n"] == 1, "a second identical request must be served from cache"

    # Now pretend the cached files changed on disk.
    monkeypatch.setattr(
        deconvolve, "reference_fingerprint", lambda *a, **k: ("CHANGED",)
    )
    deconvolve.build_signatures(genes)
    assert calls["n"] == 2, "a changed reference must invalidate the memo"


@needs_cache
def test_signature_cache_does_not_pin_the_reference_frame(monkeypatch):
    """What is cached must be the small thing, not the 322.6 MiB frame.

    The reference frame is 23,460 x 3,589; a signature is (G, 4). If a future
    edit "simplifies" this by memoising ``load_reference`` instead, the cached
    payload grows by three orders of magnitude and this goes red.

    PLANTED: storing ``(expr_df, labels_df)`` in ``_SIGNATURE_CACHE`` alongside
    the signature takes the measured payload past the bound and fails.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit", lambda: {"bytes": None, "source": "test"}
    )
    genes = list(scref_gbm.reference_gene_symbols()[:300])
    deconvolve.build_signatures(genes)
    assert deconvolve.signature_cache_info()["entries"] == 1

    payload = 0
    for key, (sig, names) in deconvolve._SIGNATURE_CACHE.items():
        assert isinstance(sig, np.ndarray), "the cached value must be the signature"
        payload += sig.nbytes
    # Two orders of magnitude below the measured resident frame.
    assert payload < scref_gbm.MEASURED_RESIDENT_MIB * 1024 * 1024 / 100


@needs_cache
def test_signature_cache_is_bounded(monkeypatch):
    """An unbounded per-gene-set cache is a slow memory leak on a long-lived
    server. The bound is asserted here so removing it is visible.

    PLANTED: deleting the ``while len(...) > _SIGNATURE_CACHE_MAX`` eviction
    loop makes ``entries`` grow past the maximum and this goes red.
    """
    monkeypatch.setattr(
        scref_gbm, "_host_memory_limit", lambda: {"bytes": None, "source": "test"}
    )
    all_genes = list(scref_gbm.reference_gene_symbols())
    for i in range(deconvolve._SIGNATURE_CACHE_MAX + 3):
        deconvolve.build_signatures(all_genes[i:i + 50])
    info = deconvolve.signature_cache_info()
    assert info["entries"] <= deconvolve._SIGNATURE_CACHE_MAX


# --------------------------------------------------------------------------- #
# 5. The platform trap the brief called out by name
# --------------------------------------------------------------------------- #
def test_maxrss_is_reported_in_bytes_on_darwin_and_kib_elsewhere(monkeypatch):
    """``ru_maxrss`` is BYTES on macOS and KiB on Linux. Getting this wrong is a
    1,024x error in a number whose only job is to be believed.

    PLANTED: removing the ``sys.platform`` branch from ``_maxrss_bytes`` makes
    one of the two assertions below fail whichever way it is removed.
    """
    import resource

    class _RU:
        ru_maxrss = 1_048_576  # 1 MiB in bytes, or 1 GiB if read as KiB

    monkeypatch.setattr(resource, "getrusage", lambda who: _RU())

    monkeypatch.setattr(scref_gbm.sys, "platform", "darwin")
    assert scref_gbm._maxrss_bytes() == 1_048_576

    monkeypatch.setattr(scref_gbm.sys, "platform", "linux")
    assert scref_gbm._maxrss_bytes() == 1_048_576 * 1024
