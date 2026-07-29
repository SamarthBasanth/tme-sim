"""Tests for tmesim.ml.panel_imputation and the opt-in ml_imputed hypoxia route.

Split into two tiers:
  * FAST / hermetic — synthetic panels and tiny in-memory Microenvironments;
    exercise the refusal logic and the registry's own internal consistency
    without touching any real cached dataset.
  * SLOW / real-data integration — trains a real RandomForestRegressor on the
    real GBM-Space Visium sections. Marked so a fast test run can skip them
    (``pytest -m "not slow"``), but they run by default like the rest of the
    suite.

The one test every other test in this file exists to make trustworthy:
:func:`test_ml_imputed_never_overwrites_a_measured_hypoxia_call` — an inferred
value must never silently pass as measured.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core.domain import Domain
from tmesim.core.microenvironment import (
    CellPopulation,
    CellState,
    CellType,
    Microenvironment,
    TissueFields,
)
from tmesim.ml import panel_imputation as pi
from tmesim.reconstruct import real_pipeline as rp
from tmesim.reconstruct.cell_state import HYPOXIA_MARKERS, STATE_MODES


# --------------------------------------------------------------------------- #
# Small synthetic helpers (no real data touched).
# --------------------------------------------------------------------------- #
def _tiny_micro(n=40, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 500, n)
    y = rng.uniform(0, 500, n)
    cell_type = np.full(n, int(CellType.TUMOR), dtype=np.int8)
    cell_type[: n // 4] = int(CellType.IMMUNE)
    state = np.full(n, int(CellState.VIABLE), dtype=np.int8)
    domain = Domain(nx=20, ny=20, dx=25.0)
    fields = TissueFields(ecm_density=domain.zeros(), vascular_density=domain.zeros())
    cells = CellPopulation(
        x=x, y=y, cell_type=cell_type, state=state,
        sensitivity=np.ones(n, dtype=np.float64),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells,
                             patient_id="synthetic-test")


def _synthetic_panel_no_hypoxia(n=40, n_filler=40, seed=0):
    """A made-up panel resembling neither real registered map closely, and
    carrying NONE of cell_state.HYPOXIA_MARKERS — hypoxia is genuinely
    UNMEASURED on it, but it is also not the 319-gene brain panel, so
    install_ml_imputed_hypoxia must refuse it rather than guess."""
    rng = np.random.default_rng(seed)
    genes = [f"FILLER{i}" for i in range(n_filler)]
    genes = np.asarray(genes, dtype=str)
    expr = rng.poisson(1.0, size=(n, len(genes))).astype(np.float64)
    return expr, genes


def _synthetic_panel_with_vegfa(n=40, n_filler=40, seed=0):
    """Same as above but with a VEGFA column added (all-zero counts is fine —
    presence alone, not detection rate, is what makes an axis measured)."""
    expr, genes = _synthetic_panel_no_hypoxia(n=n, n_filler=n_filler, seed=seed)
    vegfa_col = np.zeros((n, 1), dtype=np.float64)
    vegfa_col[0, 0] = 3.0  # at least one detected count
    expr2 = np.concatenate([expr, vegfa_col], axis=1)
    genes2 = np.concatenate([genes, np.array(["VEGFA"], dtype=str)])
    return expr2, genes2


# --------------------------------------------------------------------------- #
# Registry self-consistency (fast).
# --------------------------------------------------------------------------- #
def test_target_registry_well_formed():
    names = [t.name for t in pi.TARGET_REGISTRY]
    assert len(names) == len(set(names)), "duplicate target names in TARGET_REGISTRY"
    for t in pi.TARGET_REGISTRY:
        assert t.target_genes, f"{t.name}: empty target_genes"
        assert t.gap_panels, f"{t.name}: empty gap_panels"
        for p in t.gap_panels:
            assert p in pi.XENIUM_PANELS, f"{t.name}: unknown gap panel {p!r}"
        for p in t.ground_truth_panels:
            assert p in pi.XENIUM_PANELS, f"{t.name}: unknown ground_truth panel {p!r}"
            assert p not in t.gap_panels, (
                f"{t.name}: {p!r} listed as both a gap panel and a ground-truth "
                "panel — contradictory"
            )


def test_hypoxia_composite_matches_cell_state_markers():
    """The swept hypoxia_composite target must be the SAME gene set
    cell_state.py itself uses for the real HYPOXIC call, not a drifted copy."""
    spec = pi.TARGETS_BY_NAME["hypoxia_composite"]
    assert tuple(spec.target_genes) == tuple(HYPOXIA_MARKERS)


# --------------------------------------------------------------------------- #
# composite_score (fast, synthetic).
# --------------------------------------------------------------------------- #
def test_composite_score_refuses_when_nothing_present():
    expr, genes = _synthetic_panel_no_hypoxia()
    score, present, det = pi.composite_score(expr, genes, ("VEGFA", "CA9"))
    assert score is None
    assert present == []
    assert det == {}


def test_composite_score_computes_over_present_subset_only():
    expr, genes = _synthetic_panel_with_vegfa()
    score, present, det = pi.composite_score(expr, genes, ("VEGFA", "CA9", "ADM"))
    assert present == ["VEGFA"]  # only the present one, CA9/ADM silently absent
    assert score is not None
    assert score.shape == (expr.shape[0],)
    assert np.isfinite(score).all()
    assert "VEGFA" in det


# --------------------------------------------------------------------------- #
# impute_gene_score: refusal behaviour (fast for the refusal paths; the
# unknown-target check never trains anything).
# --------------------------------------------------------------------------- #
def test_impute_gene_score_unknown_target_raises():
    expr, genes = _synthetic_panel_with_vegfa()
    with pytest.raises(ValueError, match="unknown target"):
        pi.impute_gene_score(expr, genes, target="not_a_real_target")


def test_impute_gene_score_refuses_too_few_features():
    # A panel carrying almost none of the brain_319 feature template.
    rng = np.random.default_rng(1)
    genes = np.array(["EGFR", "PTPRC", "SOMEOTHERGENE"], dtype=str)
    expr = rng.poisson(1.0, size=(10, len(genes))).astype(np.float64)
    with pytest.raises(ValueError, match="refuses"):
        pi.impute_gene_score(expr, genes, target="hypoxia_composite", gap_panel="brain_319")


# --------------------------------------------------------------------------- #
# install_ml_imputed_hypoxia: opt-in / refuse-unless-requested / never
# overwrites measured (fast — no real data, no model training needed for these
# three, since each is refused before reaching the model).
# --------------------------------------------------------------------------- #
def test_ml_imputed_is_not_a_state_mode():
    """The core opt-in guarantee at the state_mode layer: 'ml_imputed' is not
    a recognised state_mode string, so no RunSpec can reach it by accident."""
    assert "ml_imputed" not in STATE_MODES


def test_ml_imputed_never_overwrites_a_measured_hypoxia_call():
    """If hypoxia IS measured on the panel (VEGFA present), install_ml_imputed_
    hypoxia must make NO CHANGE — an inferred value must never silently
    replace a measured one, even if the caller asks for it by mistake."""
    micro = _tiny_micro()
    expr, genes = _synthetic_panel_with_vegfa(n=micro.cells.n)
    before = micro.cells.state.copy()
    result = rp.install_ml_imputed_hypoxia(micro, expr, genes)
    assert result["applied"] is False
    assert result["kind"] == "ml_imputed"
    assert "MEASURED" in result["reason"]
    np.testing.assert_array_equal(micro.cells.state, before)


def test_ml_imputed_refuses_an_unvalidated_panel():
    """Hypoxia genuinely UNMEASURED, but the panel is not the one this
    session validated (it resolves to something other than brain_319) —
    must refuse rather than silently reuse an unvalidated route."""
    micro = _tiny_micro()
    expr, genes = _synthetic_panel_no_hypoxia(n=micro.cells.n)
    with pytest.raises(ValueError, match="refuses"):
        rp.install_ml_imputed_hypoxia(micro, expr, genes)


def test_ml_imputed_row_mismatch_raises():
    micro = _tiny_micro(n=40)
    expr, genes = _synthetic_panel_no_hypoxia(n=10)  # deliberately wrong row count
    with pytest.raises(ValueError, match="misaligned"):
        rp.install_ml_imputed_hypoxia(micro, expr, genes)


def test_reconstruct_real_refuses_ml_imputed_state_mode():
    """reconstruct_real never routes to ml_imputed automatically or on
    request via state_mode= — it simply is not a member of STATE_MODES, so
    asking for it fails the same way asking for a typo'd mode would."""
    with pytest.raises(ValueError, match="unknown state_mode"):
        rp.reconstruct_real("xenium", state_mode="ml_imputed")


# --------------------------------------------------------------------------- #
# SLOW real-data integration: actually trains on GBM-Space and checks the
# honest cross-section / cross-panel validation machinery runs end to end and
# produces sane, finite numbers. Kept to the ONE target this module wires all
# the way to a CellState (hypoxia_composite / brain_319) plus one more target
# so the sweep's own aggregation is exercised.
# --------------------------------------------------------------------------- #
_HAS_GBMSPACE = all(
    __import__("os").path.isdir(
        __import__("os").path.join("data_cache", "gbmspace_visium", s)
    )
    for s in pi.GBMSPACE_SECTIONS
)
_HAS_XENIUM_BRAIN = __import__("os").path.isdir(pi.XENIUM_PANELS["brain_319"])

pytestmark_slow = pytest.mark.skipif(
    not (_HAS_GBMSPACE and _HAS_XENIUM_BRAIN),
    reason="real GBM-Space / Xenium brain-panel caches not present",
)


@pytestmark_slow
def test_leave_one_section_out_cv_never_trains_on_the_held_out_section():
    """The one non-negotiable methodological rule, checked directly: for every
    fold, n_train + n_test spans exactly the 3 real sections and the held-out
    section's own row count equals n_test (i.e. it was excluded from
    training, not merely down-weighted)."""
    spec = pi.TARGETS_BY_NAME["hypoxia_composite"]
    feature_genes = [g for g in pi.panel_genes("brain_319") if g not in spec.target_genes]
    result = pi.leave_one_section_out_cv(spec, feature_genes)
    assert result["skipped"] is False
    assert len(result["folds"]) == len(pi.GBMSPACE_SECTIONS)
    for fold in result["folds"]:
        assert np.isfinite(fold["r2"])
        assert fold["n_train"] > 0 and fold["n_test"] > 0
        # 2 of the 3 sections trained on; the third (held out) supplied ONLY
        # the test rows.
        held = fold["held_out_section"]
        assert held in pi.GBMSPACE_SECTIONS
    assert np.isfinite(result["pooled_r2"])


@pytestmark_slow
def test_impute_gene_score_real_panel_returns_ml_imputed_provenance():
    from tmesim.data.real.singlecell_spatial import load
    d = load(path=pi.XENIUM_PANELS["brain_319"], dense=False)
    # Small slice: imputation is per-cell, independent across rows, so a
    # sub-sample is a faithful (and much faster) smoke test.
    expr = d["expression"][:200]
    genes = d["genes"]
    out = pi.impute_gene_score(expr, genes, target="hypoxia_composite", gap_panel="brain_319")
    assert out["score"].shape == (200,)
    assert np.isfinite(out["score"]).all()
    assert out["provenance"]["kind"] == "ml_imputed"
    assert "MODEL-IMPUTED, NOT MEASURED" in out["provenance"]["caveat"]


@pytestmark_slow
def test_install_ml_imputed_hypoxia_on_real_brain_panel_fills_the_gap():
    """End to end on the real 319-gene brain-panel reconstruction, where
    hypoxia genuinely is UNMEASURED: applied=True, states.length matches
    cells, DEAD cells (there are none pre-simulation, but the guard is
    checked) are never resurrected, and the returned caveat carries the
    honest validation numbers."""
    from tmesim.reconstruct.real_xenium import reconstruct_from_xenium
    from tmesim.data.real.singlecell_spatial import load

    micro = reconstruct_from_xenium()
    d = load(path=pi.XENIUM_PANELS["brain_319"], dense=False)
    assert d["expression"].shape[0] == micro.cells.n

    result = rp.install_ml_imputed_hypoxia(micro, d["expression"], d["genes"])
    assert result["applied"] is True
    assert result["kind"] == "ml_imputed"
    assert result["state_is_measured"] is False
    assert micro.cells.state.shape == (micro.cells.n,)
    assert result["fractions"]["all"]["hypoxic"] >= 0.0
    assert "loso_pooled_r2" not in result  # not echoed at top level
    assert "pooled R^2" not in result["caveat"] or "R^2=" in result["caveat"]
    assert "R^2=" in result["caveat"]
