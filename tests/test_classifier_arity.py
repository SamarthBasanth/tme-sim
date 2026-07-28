"""Pins the cell-type classifier's ARITY fix on a REAL section (no synthetic data).

Background. ``reconstruct.real_xenium.phenotype_cells`` used to gate the RAW mean
marker z-score at a fixed 1.0. A mean of ``n`` unit-variance z-scores has SD
``sqrt((1 + (n-1)*rho)/n) <= 1``, so a fixed gate was mechanically easier to clear
the FEWER markers a type had — which inflated the one-marker ``ENDOTHELIAL`` class
on the Xenium brain panel and deflated it on the three-marker IO panel, swinging
the reconstructed endothelial fraction 9.74 / 4.18 / 12.19 % over three real GBM
sections. See ``results/classifier_fix.md``.

These tests pin, on the real 40,887-cell Xenium brain section:

1. the BUG still reproduces exactly under ``aggregation='mean_z'`` (so the legacy
   numbers stay regenerable and the defect stays documented);
2. the arity relationship itself — measured ``sd(mean-z)`` equals the closed form
   and the one-marker class is strictly the widest — so the bug cannot silently
   return through a marker-map edit;
3. the FIXED default no longer orders score spread by marker count, and produces
   the composition recorded in ``results/classifier_fix.md``;
4. the whole thing is deterministic and ``marker_scores`` is consistent with the
   labels ``phenotype_cells`` derives from it.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core import CellType
from tmesim.reconstruct.marker_panels import BRAIN_PANEL_MARKER_MAP
from tmesim.reconstruct.real_xenium import (
    TYPING_AGGREGATION,
    TYPING_AGGREGATIONS,
    TYPING_THRESHOLD,
    marker_scores,
    phenotype_cells,
)

#: Section 1 — the 10x Xenium Human Brain 319-gene GBM section (real bytes).
_SECTION = "data_cache/xenium_gbm"

#: Composition under the LEGACY, arity-biased rule. Unchanged from the numbers
#: recorded throughout the repo before the fix; kept so `aggregation='mean_z'`
#: provably still regenerates history.
LEGACY_PCT = {"TUMOR": 73.80, "STROMAL": 7.92, "IMMUNE": 8.54, "ENDOTHELIAL": 9.74}

#: Composition under the FIXED default (`matched_null_z`), from
#: `results/classifier_fix.md`. Deterministic: no RNG anywhere in the path.
FIXED_PCT = {"TUMOR": 58.24, "STROMAL": 11.78, "IMMUNE": 19.48, "ENDOTHELIAL": 10.50}

#: Percentage-point tolerance. The computation is deterministic on fixed bytes;
#: this only absorbs BLAS summation-order jitter.
TOL_PP = 0.05


@pytest.fixture(scope="module")
def xenium():
    """Raw real Xenium section-1 counts + gene symbols (retries once)."""
    from tmesim.data.real.singlecell_spatial import load

    try:
        d = load(path=_SECTION)
    except Exception:
        import time

        time.sleep(3)
        d = load(path=_SECTION)
    return d["expression"], np.asarray(d["genes"], dtype=str)


def _pct(labels: np.ndarray) -> dict:
    return {c.name: 100.0 * float((labels == int(c)).mean()) for c in CellType}


# --------------------------------------------------------------------------- #
# 1. the bug is still reproducible under the legacy rule
# --------------------------------------------------------------------------- #
def test_legacy_mean_z_reproduces_the_historical_composition(xenium):
    expression, genes = xenium
    labels, _used, _missing = phenotype_cells(
        expression, genes, aggregation="mean_z"
    )
    got = _pct(labels)
    for name, want in LEGACY_PCT.items():
        assert got[name] == pytest.approx(want, abs=TOL_PP), got


# --------------------------------------------------------------------------- #
# 2. the arity relationship itself
# --------------------------------------------------------------------------- #
def test_mean_z_score_spread_is_a_function_of_marker_count(xenium):
    """sd(mean-z) == sqrt((1+(n-1)rho)/n): the bug, as arithmetic."""
    expression, genes = xenium
    scores, types, used_map, _ = marker_scores(
        expression, genes, aggregation="mean_z"
    )

    sds = {}
    for i, ct in enumerate(types):
        markers = used_map[ct]
        n = len(markers)
        assert n > 0, ct
        col = scores[:, i]
        sds[ct] = float(col.std())

        # Closed form, using the measured mean pairwise marker correlation.
        per_marker = np.column_stack([
            marker_scores(expression, genes, aggregation="mean_z",
                          marker_map={ct: [m]})[0][:, 0]
            for m in markers
        ])
        rho = (float(np.corrcoef(per_marker.T)[np.triu_indices(n, 1)].mean())
               if n > 1 else 0.0)
        assert sds[ct] == pytest.approx(np.sqrt((1 + (n - 1) * rho) / n), abs=1e-6)

    # A single-marker type's score IS one unit-variance z-score, and it is
    # strictly the widest-spread class -- which is why a fixed gate favoured it.
    assert len(used_map[CellType.ENDOTHELIAL]) == 1
    assert sds[CellType.ENDOTHELIAL] == pytest.approx(1.0, abs=1e-9)
    assert sds[CellType.ENDOTHELIAL] == max(sds.values())
    for ct, sd in sds.items():
        if len(used_map[ct]) > 1:
            assert sd < 0.9, (ct, sd, len(used_map[ct]))

    # ... and the gate-clearing rate follows the spread, not the biology.
    clear = {ct: float((scores[:, i] >= TYPING_THRESHOLD).mean())
             for i, ct in enumerate(types)}
    assert clear[CellType.ENDOTHELIAL] == max(clear.values())


def test_matched_null_z_does_not_order_score_spread_by_marker_count(xenium):
    """The fix: the one-marker class is no longer the widest-spread class."""
    expression, genes = xenium
    scores, types, used_map, _ = marker_scores(expression, genes,
                                               aggregation="matched_null_z")
    sds = {ct: float(scores[:, i].std()) for i, ct in enumerate(types)}
    n_markers = {ct: len(used_map[ct]) for ct in types}

    assert n_markers[CellType.ENDOTHELIAL] == 1
    assert sds[CellType.ENDOTHELIAL] < max(sds.values()), sds

    # No monotone relationship between marker count and spread survives: at least
    # one multi-marker type is WIDER and at least one is NARROWER than the
    # single-marker one (under 'mean_z' every multi-marker type is narrower).
    multi = [sds[ct] for ct in types if n_markers[ct] > 1]
    assert max(multi) > sds[CellType.ENDOTHELIAL]
    assert min(multi) < sds[CellType.ENDOTHELIAL]


def test_max_z_is_arity_biased_in_the_opposite_direction(xenium):
    """Why 'max_z' was not chosen: E[max of n] grows with n."""
    expression, genes = xenium
    scores, types, used_map, _ = marker_scores(expression, genes, aggregation="max_z")
    clear = {ct: float((scores[:, i] >= TYPING_THRESHOLD).mean())
             for i, ct in enumerate(types)}
    single = clear[CellType.ENDOTHELIAL]
    multi = [clear[ct] for ct in types if len(used_map[ct]) > 1]
    # The single-marker class now clears LEAST, by a wide margin.
    assert single == min(clear.values())
    assert min(multi) > 2.0 * single, clear


# --------------------------------------------------------------------------- #
# 3. the fixed default, pinned
# --------------------------------------------------------------------------- #
def test_default_aggregation_is_the_matched_null():
    assert TYPING_AGGREGATION == "matched_null_z"
    assert TYPING_AGGREGATION in TYPING_AGGREGATIONS


def test_fixed_classifier_composition_on_section_one(xenium):
    expression, genes = xenium
    labels, used_map, missing = phenotype_cells(expression, genes)
    got = _pct(labels)

    for name, want in FIXED_PCT.items():
        assert got[name] == pytest.approx(want, abs=TOL_PP), got

    # Still a GBM section: malignant cells dominate, every niche compartment is
    # populated, and no cell is phenotyped NECROTIC (no necrosis marker on panel).
    assert got["TUMOR"] > 50.0
    assert all(got[c] > 1.0 for c in ("STROMAL", "IMMUNE", "ENDOTHELIAL"))
    assert got["NECROTIC"] == 0.0
    assert sum(got.values()) == pytest.approx(100.0, abs=1e-6)

    # The marker map is untouched by the fix -- only the scoring changed.
    assert used_map[CellType.ENDOTHELIAL] == ["PECAM1"]
    assert used_map == {ct: list(ms) for ct, ms in BRAIN_PANEL_MARKER_MAP.items()}
    assert missing == []


def test_endothelial_fraction_moved_toward_the_independent_estimates(xenium):
    """Both external references sit near 8-10%; the fix must not overshoot.

    NOT a tuning target -- the gate was frozen at 1.0 for every candidate. This
    only pins that the corrected value stays in the single-digit-to-low-teens
    band the reconstruction is documented to produce.
    """
    expression, genes = xenium
    fixed = _pct(phenotype_cells(expression, genes)[0])["ENDOTHELIAL"]
    assert 5.0 < fixed < 15.0, fixed


# --------------------------------------------------------------------------- #
# 4. contract: determinism, validation, score/label consistency
# --------------------------------------------------------------------------- #
def test_phenotyping_is_deterministic(xenium):
    expression, genes = xenium
    a, _, _ = phenotype_cells(expression, genes)
    b, _, _ = phenotype_cells(expression, genes)
    assert np.array_equal(a, b)


def test_unknown_aggregation_raises(xenium):
    expression, genes = xenium
    with pytest.raises(ValueError, match="unknown aggregation"):
        phenotype_cells(expression, genes, aggregation="mean")


@pytest.mark.parametrize("rule", TYPING_AGGREGATIONS)
def test_labels_follow_the_gated_argmax_of_marker_scores(xenium, rule):
    expression, genes = xenium
    labels, _u, _m = phenotype_cells(expression, genes, aggregation=rule)
    scores, types, _u2, _m2 = marker_scores(expression, genes, aggregation=rule)

    gated = np.where(scores >= TYPING_THRESHOLD, scores, -np.inf)
    expected = np.array([int(types[i]) for i in gated.argmax(axis=1)], dtype=np.int8)
    expected[~np.isfinite(gated.max(axis=1))] = int(CellType.TUMOR)
    assert np.array_equal(labels, expected)

    # Cells clearing no gate fall back to the malignant background.
    ungated = ~np.isfinite(gated.max(axis=1))
    assert (labels[ungated] == int(CellType.TUMOR)).all()


def test_every_aggregation_produces_a_valid_label_vector(xenium):
    expression, genes = xenium
    n = expression.shape[0]
    valid = {int(c) for c in CellType}
    for rule in TYPING_AGGREGATIONS:
        labels, _u, _m = phenotype_cells(expression, genes, aggregation=rule)
        assert labels.shape == (n,)
        assert labels.dtype == np.int8
        assert set(np.unique(labels)).issubset(valid), rule
        assert (labels == int(CellType.TUMOR)).any(), rule
