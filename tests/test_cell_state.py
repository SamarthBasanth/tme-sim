"""Tests for the MEASURED marker-derived cell state, and for the default staying put.

Two jobs, and the second is the more important one:

1. :mod:`tmesim.reconstruct.cell_state` behaves — it REFUSES on a panel that
   cannot support a call (the 319-gene brain panel carries no hypoxia-response
   gene, so HYPOXIC must be declared UNMEASURED and assigned to nobody), it falls
   back to a neighbourhood aggregation rather than fabricating a per-cell call out
   of a 0.35% detection rate, it is deterministic, it is completely independent of
   the simulation mesh, and a row misalignment is a hard error.
2. **``reconstruct_real``'s default is BIT-IDENTICAL.** The whole point of adding
   ``state_mode`` was to make a measured alternative available WITHOUT moving any
   published number, so the default path is asserted not to touch
   ``cells.state`` at all (``np.array_equal`` against a pre-option-equivalent
   reconstruction), and the historical ``metadata['real_pipeline']`` /
   ``reconstruction_report`` keys are asserted still present with their historical
   values.

Everything runs on the real cached 10x Xenium GBM section at a coarse grid, via
the session-scoped ``real_micro`` fixture where possible, so the suite stays fast.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core import CellState, CellType
from tmesim.data.real.singlecell_spatial import load as load_singlecell_spatial
from tmesim.reconstruct.cell_state import (
    CELL_CYCLE_SPECIFIC_MARKERS,
    HYPOXIA_MARKERS,
    HYPOXIA_SPECIFIC_MARKERS,
    LOW_SPECIFICITY_MARKERS,
    PER_CELL_DETECTION_FLOOR,
    PROLIFERATION_MARKERS,
    STATE_MODES,
    marker_detection,
    neighbourhood_positive_fraction,
    panel_coverage,
    state_fractions_from_states,
    state_from_markers,
    state_from_markers_micro,
)
from tmesim.reconstruct.real_pipeline import (
    STATE_MODE_DEFAULT,
    reconstruct_real,
    reconstruction_report,
)

from .conftest import REAL_DX_UM, REAL_SOURCE

BRAIN_SAMPLE = "data_cache/xenium_gbm"


@pytest.fixture(scope="module")
def brain_panel():
    """Raw expression + panel genes for the 319-gene flagship section."""
    try:
        return load_singlecell_spatial(path=BRAIN_SAMPLE, dense=False)
    except Exception:
        import time

        time.sleep(3)
        return load_singlecell_spatial(path=BRAIN_SAMPLE, dense=False)


# --------------------------------------------------------------------------- #
# The marker sets are declared consistently
# --------------------------------------------------------------------------- #
def test_marker_sets_are_consistent():
    """The specific and low-specificity proliferation sets partition the full set,
    and HIF1A is absent from the hypoxia set on mechanism."""
    assert set(CELL_CYCLE_SPECIFIC_MARKERS) | set(LOW_SPECIFICITY_MARKERS) == set(
        PROLIFERATION_MARKERS
    )
    assert not set(CELL_CYCLE_SPECIFIC_MARKERS) & set(LOW_SPECIFICITY_MARKERS)
    assert set(HYPOXIA_SPECIFIC_MARKERS) <= set(HYPOXIA_MARKERS)
    assert "HIF1A" not in HYPOXIA_MARKERS
    assert 0.0 < PER_CELL_DETECTION_FLOOR < 1.0
    assert STATE_MODE_DEFAULT in STATE_MODES


# --------------------------------------------------------------------------- #
# The refusals — the whole point of the module
# --------------------------------------------------------------------------- #
def test_hypoxia_is_declared_unmeasured_on_the_brain_panel(brain_panel):
    """The 319-gene panel carries no VEGFA/CA9/ADM/NDRG1/SLC2A1, so HYPOXIC must be
    UNMEASURED and assigned to NO cell — not silently computed as zero."""
    cov = panel_coverage(brain_panel["expression"], brain_panel["genes"])
    h = cov["hypoxia"]
    assert h["present"] == []
    assert set(h["missing"]) == set(HYPOXIA_MARKERS)
    assert h["call_level"] == "UNMEASURED"
    assert "UNMEASURED" in h["reason"]
    assert "HYPOXIC" in cov["unmeasured_states"]

    xy = np.asarray(brain_panel["cell_xy_um"], dtype=np.float64)
    out = state_from_markers(
        brain_panel["expression"], brain_panel["genes"],
        x=xy[:, 0], y=xy[:, 1], scoring_population="all",
    )
    assert not (out["states"] == int(CellState.HYPOXIC)).any(), (
        "a cell was labelled HYPOXIC on a panel that measures no hypoxia gene"
    )
    assert out["fractions"]["all"]["hypoxic"] == 0.0
    assert "UNMEASURED" in out["provenance"]


def test_per_cell_proliferation_call_is_refused_at_brain_panel_dropout(brain_panel):
    """MKI67 ~0.35% and TOP2A ~1.08% on this panel, against a published 19-40%
    Ki-67 index. A per-cell call there is probe sensitivity, so the module must
    fall back to the neighbourhood level and SAY so."""
    cov = panel_coverage(brain_panel["expression"], brain_panel["genes"])
    p = cov["proliferation"]
    assert set(p["present"]) == set(CELL_CYCLE_SPECIFIC_MARKERS)
    assert p["union_detection"] < PER_CELL_DETECTION_FLOOR
    assert p["call_level"] == "neighbourhood"
    assert "BELOW" in p["reason"]
    assert "CONVENTION" in p["binarisation"]


def test_neighbourhood_axis_is_unmeasured_without_coordinates(brain_panel):
    """No coordinates -> the neighbourhood aggregation is impossible -> the axis is
    declared UNMEASURED rather than guessed."""
    out = state_from_markers(brain_panel["expression"], brain_panel["genes"])
    p = out["coverage"]["proliferation"]
    assert p["call_level"] == "UNMEASURED"
    assert "coordinates" in p["reason"]
    # Nothing is labelled QUIESCENT when proliferation is unmeasured.
    assert not (out["states"] == int(CellState.QUIESCENT)).any()


def test_low_specificity_marker_is_excluded_by_default(brain_panel):
    """PCNA is on this panel at ~46% and would buy a per-cell call by measuring
    something other than proliferation. It must be off by default and must change
    the call level when switched on."""
    genes = np.asarray(brain_panel["genes"], dtype=str)
    assert all(g in genes for g in LOW_SPECIFICITY_MARKERS)

    default = panel_coverage(brain_panel["expression"], genes)
    assert "PCNA" not in default["proliferation"]["present"]
    assert default["proliferation"]["call_level"] == "neighbourhood"

    with_pcna = panel_coverage(
        brain_panel["expression"], genes,
        proliferation_markers=tuple(CELL_CYCLE_SPECIFIC_MARKERS)
        + tuple(LOW_SPECIFICITY_MARKERS),
    )
    assert "PCNA" in with_pcna["proliferation"]["present"]
    assert with_pcna["proliferation"]["union_detection"] > PER_CELL_DETECTION_FLOOR
    assert with_pcna["proliferation"]["call_level"] == "per_cell"


# --------------------------------------------------------------------------- #
# Mechanics
# --------------------------------------------------------------------------- #
def test_marker_detection_reports_presence_and_rates(brain_panel):
    md = marker_detection(brain_panel["expression"], brain_panel["genes"],
                          ("MKI67", "NOT_A_REAL_GENE"))
    assert md["present"] == ["MKI67"]
    assert md["missing"] == ["NOT_A_REAL_GENE"]
    assert 0.0 < md["detection"]["MKI67"] < 1.0
    assert md["union_detection"] == pytest.approx(md["detection"]["MKI67"])
    assert md["positive"].dtype == bool

    empty = marker_detection(brain_panel["expression"], brain_panel["genes"],
                             ("NOT_A_REAL_GENE",))
    assert empty["present"] == [] and empty["positive"] is None
    assert empty["union_detection"] == 0.0


def test_neighbourhood_fraction_on_a_hand_checkable_layout():
    """Three cells at 0, 10 and 1000 um. With R = 50 the first two see each other
    and the third is alone; below the min-cells floor every window falls back to
    the section rate."""
    x = np.array([0.0, 10.0, 1000.0])
    y = np.zeros(3)
    pos = np.array([True, False, False])

    agg = neighbourhood_positive_fraction(x, y, pos, radius_um=50.0, min_cells=1)
    assert agg["section_rate"] == pytest.approx(1.0 / 3.0)
    assert agg["fraction"][0] == pytest.approx(0.5)   # {self, neighbour}, 1 positive
    assert agg["fraction"][1] == pytest.approx(0.5)
    assert agg["fraction"][2] == pytest.approx(0.0)   # alone, not positive
    assert agg["n_unsupported"] == 0

    # Raise the floor above every window size -> all windows fall back.
    fallback = neighbourhood_positive_fraction(x, y, pos, radius_um=50.0,
                                               min_cells=99)
    assert np.allclose(fallback["fraction"], fallback["section_rate"])
    assert fallback["n_unsupported"] == 3

    with pytest.raises(ValueError):
        neighbourhood_positive_fraction(x, y, pos, radius_um=0.0)


def test_states_are_deterministic_and_reproducible(real_micro, brain_panel):
    a = state_from_markers_micro(real_micro, brain_panel["expression"],
                                 brain_panel["genes"])
    b = state_from_markers_micro(real_micro, brain_panel["expression"],
                                 brain_panel["genes"])
    assert np.array_equal(a["states"], b["states"])
    assert a["states"].dtype == np.int8
    assert set(np.unique(a["states"])) <= {int(s) for s in CellState}


def test_state_fractions_sum_to_one(real_micro, brain_panel):
    out = state_from_markers_micro(real_micro, brain_panel["expression"],
                                   brain_panel["genes"])
    for key in ("all", "tumor", "scored"):
        fr = out["fractions"][key]
        assert sum(fr.values()) == pytest.approx(1.0)
    # And it agrees with the standalone helper on the same mask.
    tumor = np.asarray(real_micro.cells.cell_type) == int(CellType.TUMOR)
    assert state_fractions_from_states(out["states"], tumor) == out["fractions"]["tumor"]


def test_dead_cells_are_never_resurrected(real_micro, brain_panel):
    cells = real_micro.cells.copy()
    cells.state[:100] = int(CellState.DEAD)
    xy_x, xy_y = cells.x, cells.y
    out = state_from_markers(
        brain_panel["expression"], brain_panel["genes"], x=xy_x, y=xy_y,
        cell_type=cells.cell_type, state=cells.state,
    )
    assert np.all(out["states"][:100] == int(CellState.DEAD))


def test_misalignment_is_a_hard_error(real_micro, brain_panel):
    expr = brain_panel["expression"]
    with pytest.raises(ValueError, match="not aligned"):
        state_from_markers(expr[:-1], brain_panel["genes"],
                           x=real_micro.cells.x, y=real_micro.cells.y)
    with pytest.raises(ValueError, match="not aligned"):
        state_from_markers_micro(real_micro, expr[:-1], brain_panel["genes"])


def test_marker_states_are_independent_of_the_simulation_mesh(brain_panel):
    """The property the whole recommendation rests on: the marker route never
    touches a grid, so refining the mesh cannot move a single label. Asserted with
    ``array_equal``, not a tolerance."""
    fine = reconstruct_real(source=REAL_SOURCE, dx_um=40.0)
    coarse = reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM)
    a = state_from_markers_micro(fine, brain_panel["expression"],
                                 brain_panel["genes"])
    b = state_from_markers_micro(coarse, brain_panel["expression"],
                                 brain_panel["genes"])
    assert fine.domain.shape != coarse.domain.shape
    assert np.array_equal(a["states"], b["states"])


def test_conservative_hypoxia_set_is_selectable(brain_panel):
    """The CA9-only variant is reachable through the public keyword; on this panel
    it is also absent, so the axis stays honestly UNMEASURED."""
    out = state_from_markers(
        brain_panel["expression"], brain_panel["genes"],
        hypoxia_markers=HYPOXIA_SPECIFIC_MARKERS,
    )
    assert out["coverage"]["hypoxia"]["call_level"] == "UNMEASURED"
    assert out["coverage"]["hypoxia"]["requested"] == list(HYPOXIA_SPECIFIC_MARKERS)


# --------------------------------------------------------------------------- #
# THE DEFAULT DOES NOT MOVE
# --------------------------------------------------------------------------- #
def test_default_state_mode_is_bit_identical(real_micro):
    """``state_mode='oxygen'`` (the default) must not touch ``cells.state`` at all.

    The pre-option pipeline never set states in ``reconstruct_real``; the engine
    banded them later. So the default is bit-identical iff an explicit
    ``state_mode='oxygen'`` reconstruction and the fixture (built with no
    ``state_mode`` argument at all) carry identical state arrays, and iff those
    states are whatever the cell reconstructor produced.
    """
    explicit = reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                               state_mode="oxygen")
    assert np.array_equal(explicit.cells.state, real_micro.cells.state)
    assert explicit.cells.state.dtype == real_micro.cells.state.dtype
    # Every other field is untouched too.
    assert np.array_equal(explicit.fields.oxygen, real_micro.fields.oxygen)
    assert np.array_equal(explicit.fields.vascular_density,
                          real_micro.fields.vascular_density)
    assert np.array_equal(explicit.fields.ecm_density, real_micro.fields.ecm_density)


def test_default_reconstruction_report_keys_are_unchanged(real_micro):
    """The historical report keys keep their historical values, and the new
    ``cell_state`` block is additive and labelled NOT-a-measurement."""
    rep = reconstruction_report(real_micro)
    for key in ("patient_id", "source", "n_cells", "n_tumor", "grid", "dx_um",
                "celltype_fractions", "mean_vascular_density", "mean_ecm_density",
                "mean_oxygen", "oxygen", "oxygen_mode", "sensitivity"):
        assert key in rep, f"historical report key {key!r} disappeared"
    assert rep["oxygen_mode"] == "pde"

    cs = rep["cell_state"]
    assert rep["state_mode"] == cs["mode"] == STATE_MODE_DEFAULT == "oxygen"
    assert cs["is_measured"] is False
    assert "NOT A MEASUREMENT" in cs["caveat"]
    assert cs["detail"] is None
    assert cs["unmeasured_states"] is None
    # The historical metadata block still carries the historical keys.
    rp = real_micro.metadata["real_pipeline"]
    assert rp["source"] == "xenium"
    assert rp["oxygen_mode"] == "pde"
    assert rp["state_mode"] == "oxygen"


def test_marker_state_mode_is_selectable_and_reported():
    micro = reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                             state_mode="markers")
    rep = reconstruction_report(micro)
    cs = rep["cell_state"]
    assert rep["state_mode"] == "markers"
    assert cs["is_measured"] is True
    assert "HYPOXIC" in cs["unmeasured_states"]
    assert cs["fractions_tumor"]["hypoxic"] == 0.0
    assert sum(cs["fractions_tumor"].values()) == pytest.approx(1.0)
    # The states really were installed, and they differ from the default's.
    default = reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM)
    assert not np.array_equal(micro.cells.state, default.cells.state)
    # The oxygen field is NOT affected by the state route.
    assert np.array_equal(micro.fields.oxygen, default.fields.oxygen)


def test_unknown_state_mode_and_stray_kwargs_are_errors():
    with pytest.raises(ValueError, match="unknown state_mode"):
        reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM, state_mode="nope")
    with pytest.raises(TypeError, match="state_kwargs"):
        reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                         state_mode="oxygen", state_kwargs={"radius_um": 60.0})


def test_marker_mode_refuses_without_per_cell_expression(real_micro):
    """The Visium path has no per-cell transcripts, so ``state_mode='markers'``
    must raise with a reason rather than silently degrading to the oxygen route."""
    from tmesim.reconstruct.real_pipeline import _install_states

    before = real_micro.cells.state.copy()
    with pytest.raises(ValueError, match="needs per-cell expression"):
        _install_states(real_micro, "markers")
    # And it failed before touching anything.
    assert np.array_equal(real_micro.cells.state, before)


def test_state_kwargs_reach_the_marker_route():
    """A non-default window radius must change the labels, which is the measured
    proof that the convention is real and that the keyword is wired through."""
    tight = reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                             state_mode="markers", state_kwargs={"radius_um": 50.0})
    wide = reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                            state_mode="markers", state_kwargs={"radius_um": 150.0})
    assert not np.array_equal(tight.cells.state, wide.cells.state)
    t = reconstruction_report(tight)["cell_state"]["fractions_tumor"]["viable"]
    w = reconstruction_report(wide)["cell_state"]["fractions_tumor"]["viable"]
    assert w > t, "a wider window must call more cells proliferative at this sparsity"


# --------------------------------------------------------------------------- #
# The KILL MULTIPLIERS themselves — the largest unvalidated priors in the engine
#
# `results/cell_state_multipliers_2026-09-20.md` is the claim surface these four
# tests defend. Until they existed nothing in the suite pinned the multipliers,
# their count, their per-drug resolution, or the boundedness statement that
# `results/kill_calibration.md` section 5 derives FROM them ("relabelling cells
# can move absolute kill only 2.00x"). That sentence is arithmetic on the
# multiplier VALUES, so a value change silently invalidates a published number
# with no test going red. These make it go red.
# --------------------------------------------------------------------------- #
def _uniform_exposure_micro(states, exposure_ratio=1.0, ic50=5.0, sensitivity=1.0):
    """4 cells, one per state, at IDENTICAL exposure — isolates the multiplier."""
    from tmesim.core import (
        CellPopulation, Domain, Microenvironment, TissueFields, TransportResult,
    )

    n = len(states)
    domain = Domain(nx=n, ny=1, dx=10.0)
    duration = 1000.0
    auc = np.full((1, n), exposure_ratio * ic50 * duration)
    transport = TransportResult(
        auc=auc, final_concentration=np.zeros((1, n)),
        times=np.array([0.0, duration]),
        mean_concentration_history=np.array([0.0, 0.0]),
    )
    cells = CellPopulation(
        x=(np.arange(n) + 0.5) * 10.0,
        y=np.full(n, 5.0),
        cell_type=np.full(n, int(CellType.TUMOR), dtype=np.int8),
        state=np.asarray(states, dtype=np.int8),
        sensitivity=np.full(n, float(sensitivity)),
    )
    micro = Microenvironment(
        domain=domain,
        fields=TissueFields(ecm_density=domain.zeros(),
                            vascular_density=domain.zeros()),
        cells=cells,
    )
    return micro, transport


def test_the_kill_path_applies_exactly_three_state_multipliers_and_no_fourth():
    """Cell state reaches the kill through ONE channel with exactly three values.

    Every downstream bound in `results/` (and every report that swaps states on a
    frozen transport result) rests on this. At identical exposure and identical
    sensitivity the per-cell kill probability must be the Hill term times exactly
    `1.0 / effective_hypoxia_resistance / quiescent_resistance / 0.0`. A fourth
    distinct value, or any state-dependent term the multiplier does not explain,
    breaks the decomposition rather than merely shifting it.
    """
    from tmesim.cells.pkpd import apply_drug_kill, effective_hypoxia_resistance
    from tmesim.core import DrugParams, PDConfig

    order = [CellState.VIABLE, CellState.HYPOXIC, CellState.QUIESCENT,
             CellState.DEAD]
    micro, transport = _uniform_exposure_micro([int(s) for s in order])
    drug = DrugParams(name="probe", ic50=5.0, hill=1.5)
    pd_cfg = PDConfig()

    diag = {}
    p = apply_drug_kill(micro, transport, pd_cfg, drug,
                        rng=np.random.default_rng(0), diagnostics=diag)
    p_hill = float(diag["p_hill"][0])
    assert p_hill > 0.0

    expected = {
        CellState.VIABLE: 1.0,
        CellState.HYPOXIC: effective_hypoxia_resistance(drug, pd_cfg),
        CellState.QUIESCENT: float(pd_cfg.quiescent_resistance),
        CellState.DEAD: 0.0,
    }
    for i, s in enumerate(order):
        assert p[i] == pytest.approx(p_hill * expected[s], rel=1e-12), (
            f"{s.name}: kill probability is not p_hill x {expected[s]}"
        )
    # Exactly three non-zero multipliers, and they are the documented ones.
    live = sorted({round(float(p[i]) / p_hill, 12) for i, s in enumerate(order)
                   if s is not CellState.DEAD})
    assert live == sorted({1.0, expected[CellState.HYPOXIC],
                           expected[CellState.QUIESCENT]}), live


def test_the_shipped_multiplier_values_are_pinned_to_the_published_claim():
    """The VALUES `results/` quotes. Change one and this goes red by design.

    `results/kill_calibration.md` section 5 and
    `results/cell_state_multipliers_2026-09-20.md` both quote x0.80 / x0.50 and
    derive a 2.00x bound from them; `results/cell_state_impact.md` section 3.1
    quotes the global x0.35. None of those documents has a test behind it, so a
    revision to `PDConfig` would leave three published numbers wrong and every
    check green. This is that check.
    """
    from tmesim.core import PDConfig
    from tmesim.core.config import DRUG_LIBRARY, TMZ_HYPOXIA_RESISTANCE

    pd_cfg = PDConfig()
    assert pd_cfg.quiescent_resistance == 0.5, (
        "quiescent_resistance moved; results/cell_state_multipliers_2026-09-20.md "
        "and the 2.00x bound in results/kill_calibration.md section 5 are now wrong"
    )
    assert pd_cfg.hypoxia_resistance == 0.35, (
        "the GLOBAL hypoxia_resistance moved; results/cell_state_impact.md "
        "section 3.1 quotes x0.35 for cisplatin and doxorubicin"
    )
    assert pd_cfg.max_kill_rate == 1.0
    assert TMZ_HYPOXIA_RESISTANCE == 0.80, (
        "temozolomide's per-drug override moved; the x0.8 in "
        "results/kill_calibration.md section 5 is now wrong"
    )
    assert DRUG_LIBRARY["temozolomide-like"].hypoxia_resistance == 0.80
    # And it really is the ONLY drug that overrides.
    overriding = {n for n, d in DRUG_LIBRARY.items()
                  if d.hypoxia_resistance is not None}
    assert overriding == {"temozolomide-like"}, overriding


def test_hypoxia_multiplier_resolution_is_per_drug_not_global():
    """A per-drug override must win, and must be what the kill path actually used.

    `cells.pkpd.effective_hypoxia_resistance` is the only sanctioned resolver;
    reading `pd.hypoxia_resistance` directly silently ignores temozolomide's
    override. The diagnostics the kill emits must report the multiplier it
    applied, not the global default, or every attribution built on them is
    describing a run that did not happen.
    """
    from tmesim.cells.pkpd import apply_drug_kill, effective_hypoxia_resistance
    from tmesim.core import PDConfig
    from tmesim.core.config import DRUG_LIBRARY

    pd_cfg = PDConfig()
    tmz = DRUG_LIBRARY["temozolomide-like"]
    dox = DRUG_LIBRARY["doxorubicin-like"]
    assert effective_hypoxia_resistance(tmz, pd_cfg) == 0.80
    assert effective_hypoxia_resistance(dox, pd_cfg) == pd_cfg.hypoxia_resistance

    for drug in (tmz, dox):
        micro, transport = _uniform_exposure_micro(
            [int(CellState.HYPOXIC), int(CellState.VIABLE)], ic50=drug.ic50)
        diag = {}
        p = apply_drug_kill(micro, transport, pd_cfg, drug,
                            rng=np.random.default_rng(0), diagnostics=diag)
        eff = effective_hypoxia_resistance(drug, pd_cfg)
        assert diag["hypoxia_resistance"] == eff
        assert p[0] == pytest.approx(p[1] * eff, rel=1e-12)


def test_relabelling_bound_is_one_over_the_smallest_reachable_multiplier():
    """The BOUNDEDNESS claim, as arithmetic rather than as prose.

    `results/kill_calibration.md` section 5 states that because the multipliers
    are bounded, no revision of WHICH cells carry which label can move absolute
    kill by more than 2.00x. The general statement is
    `1 / min(reachable multipliers)`, it is TIGHT (both endpoints are attainable
    by labelling every cell), and it is DRUG-SPECIFIC: 2.00x holds for
    temozolomide (min(0.80, 0.50)) and is 2.857x for a drug on the global 0.35.
    The published 2.00x is therefore not a property of the engine but of one
    drug's override, which is exactly the thing an unvalidated prior can change.
    """
    from tmesim.cells.pkpd import apply_drug_kill, effective_hypoxia_resistance
    from tmesim.core import PDConfig
    from tmesim.core.config import DRUG_LIBRARY

    pd_cfg = PDConfig()
    n = 64
    rng = np.random.default_rng(3)
    for key, expected_bound in (("temozolomide-like", 2.0),
                                ("doxorubicin-like", 1.0 / 0.35)):
        drug = DRUG_LIBRARY[key]
        r_hyp = effective_hypoxia_resistance(drug, pd_cfg)
        floor = min(r_hyp, float(pd_cfg.quiescent_resistance))
        assert 1.0 / floor == pytest.approx(expected_bound, rel=1e-12)

        def mean_p(states):
            micro, transport = _uniform_exposure_micro(
                states, exposure_ratio=0.05, ic50=drug.ic50)
            return float(apply_drug_kill(
                micro, transport, pd_cfg, drug,
                rng=np.random.default_rng(0)).mean())

        all_viable = mean_p([int(CellState.VIABLE)] * n)
        assert all_viable > 0.0
        # Both endpoints attainable -> the bound is TIGHT, not merely an upper
        # limit, so a report may quote it as the movement and not only as a cap.
        worst = int(CellState.QUIESCENT if pd_cfg.quiescent_resistance <= r_hyp
                    else CellState.HYPOXIC)
        all_worst = mean_p([worst] * n)
        assert all_viable / all_worst == pytest.approx(1.0 / floor, rel=1e-12)
        # And no mixed labelling can beat it.
        for _ in range(20):
            mixed = rng.choice(
                [int(CellState.VIABLE), int(CellState.HYPOXIC),
                 int(CellState.QUIESCENT)], size=n).tolist()
            assert all_viable / mean_p(mixed) <= 1.0 / floor + 1e-12
