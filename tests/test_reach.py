"""THE SCALE-FREE REACH MAP (``app/server/reach.py``).

The point of ``reach`` is that its headline is a RATIO INSIDE ONE MAP, so the
per-compound plasma factor — the axis that broke the cross-compound ordering
(``results/pk_ranking_robustness.md``: real Cmax across the retained GBM
compounds spans 113.1x) — divides out of it. That claim is not an argument here,
it is :meth:`TestPlasmaScaleInvariance`: the AUC field is multiplied by 5,000
and every point of the ``F(x)`` ladder and every normalised falloff bin must be
unchanged.

Everything runs on the REAL Xenium GBM section from ``conftest.real_micro``
(dx = 80 um) with a short transport window, so the whole module is a few
seconds. No synthetic tissue is used for the physics; the degenerate-input tests
edit a COPY of that real result to manufacture the specific undefined case they
name.
"""
from __future__ import annotations

import copy
import os
import sys

import numpy as np
import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app.server.reach import (  # noqa: E402
    CONVENTION_NOTE,
    HEADLINE_X,
    REACH_LADDER,
    compute_reach,
)
from tmesim.core import (  # noqa: E402
    CellPopulation,
    CellType,
    Microenvironment,
    OxygenConfig,
    PDEConfig,
    SimConfig,
)
from tmesim.sim.engine import simulate_response  # noqa: E402
from tmesim.validation.delivery_metrics import penetration_length_fit  # noqa: E402
from tmesim.validation.vessel_profile import (  # noqa: E402
    direct_vessel_distance_um,
    endothelial_xy,
)

#: Plasma-scale invariance tolerance, ABSOLUTE, on quantities of order 1.
#:
#: The reach ratio is AUC_i / E_ref and both are linear in the plasma level, so
#: the factor cancels ANALYTICALLY exactly; what is asserted here is that it also
#: cancels NUMERICALLY — that nothing in the implementation (a threshold, a bin
#: edge, an accumulation order) has quietly made the answer depend on the level.
#: MEASURED on this fixture at factors 5000, 1/5000 and 113.1 (the real-Cmax
#: span of results/pk_ranking_robustness.md): every F(x) is BIT-IDENTICAL
#: (max |dF| = 0.0) and the worst normalised falloff bin moves by 3.3e-16, about
#: 1.5 ulp of 1.0 — pure float64 round-off in the rescale itself. 1e-14 is
#: therefore ~30x the largest deviation actually observed, tight enough that any
#: real level-dependence fails this test immediately.
INVARIANCE_ATOL = 1e-14

#: Short transport window; the reach map is a property of whatever exposure
#: field it is handed, so a 10-minute window exercises every code path a 6 h one
#: does at a fraction of the cost.
WINDOW_S = 600.0


@pytest.fixture(scope="module")
def real_result(real_micro, real_drug):
    """One real simulation on the real section, shared by every test here."""
    cfg = SimConfig(
        drug=real_drug,
        pde=PDEConfig(dt=10.0, total_time=WINDOW_S, convection=True,
                      save_every=500),
        seed=0,
    )
    return simulate_response(real_micro, cfg)


@pytest.fixture(scope="module")
def distance_um(real_micro):
    """The PRIMARY distance-from-nearest-endothelial-cell map, built once."""
    return direct_vessel_distance_um(real_micro.domain,
                                     endothelial_xy(real_micro.cells))


@pytest.fixture(scope="module")
def reach(real_micro, real_result, distance_um):
    out = compute_reach(real_micro, real_result, distance_um=distance_um)
    assert out["null_reason"] is None, out["null_reason"]
    return out


def _scaled_result(res, factor: float):
    """A copy of ``res`` whose AUC field is multiplied by ``factor``.

    This is EXACTLY what a different per-compound plasma level does to the
    field and nothing else: the delivery source is one multiplication by the
    plasma level and the transport PDE is linear in C, measured to a worst
    relative deviation of 1.41e-14 over a 5,200x span of plasma factors
    (``results/pk_ranking_robustness.md`` sec. 1.2). The oxygen field is
    deliberately NOT scaled — it does not depend on the drug dose, and the
    viable population must therefore be the same set of cells in both arms.
    """
    out = copy.copy(res)
    out.fields = dict(res.fields)
    out.fields["auc"] = np.asarray(res.fields["auc"], dtype=np.float64) * factor
    return out


def _retyped(micro, cell_type: CellType):
    """The same real tissue with every cell relabelled — used only to build the
    "no tumour cells" case out of a real reconstruction."""
    cells = micro.cells
    return Microenvironment(
        domain=micro.domain,
        fields=micro.fields,
        cells=CellPopulation(
            x=cells.x, y=cells.y,
            cell_type=np.full(cells.n, int(cell_type), dtype=np.int8),
            state=cells.state, sensitivity=cells.sensitivity,
        ),
        patient_id=micro.patient_id,
    )


# --------------------------------------------------------------------------- #
# The property the whole design exists for
# --------------------------------------------------------------------------- #
class TestPlasmaScaleInvariance:
    """Multiply the exposure field by an arbitrary constant; nothing moves."""

    FACTOR = 5000.0

    def test_ladder_is_unchanged_under_a_5000x_plasma_rescale(
            self, real_micro, real_result, distance_um, reach):
        scaled = compute_reach(real_micro, _scaled_result(real_result, self.FACTOR),
                               distance_um=distance_um)
        assert scaled["null_reason"] is None
        base_F = np.asarray(reach["curve"]["under_reached_fraction"])
        scaled_F = np.asarray(scaled["curve"]["under_reached_fraction"])
        assert scaled["curve"]["x"] == reach["curve"]["x"]
        # Tolerance INVARIANCE_ATOL = 1e-12, absolute, on a quantity in [0, 1].
        np.testing.assert_allclose(scaled_F, base_F, rtol=0.0,
                                   atol=INVARIANCE_ATOL)
        assert scaled["headline_under_reached_fraction"] == pytest.approx(
            reach["headline_under_reached_fraction"], abs=INVARIANCE_ATOL)
        # The same cells are being counted in both arms.
        assert scaled["n_viable"] == reach["n_viable"]
        assert scaled["n_tumor"] == reach["n_tumor"]

    def test_falloff_is_unchanged_under_a_5000x_plasma_rescale(
            self, real_micro, real_result, distance_um, reach):
        scaled = compute_reach(real_micro, _scaled_result(real_result, self.FACTOR),
                               distance_um=distance_um)
        base = reach["falloff"]
        got = scaled["falloff"]
        assert got["distance_um"] == base["distance_um"]
        assert got["count"] == base["count"]
        np.testing.assert_allclose(got["normalised_mean_auc"],
                                   base["normalised_mean_auc"],
                                   rtol=0.0, atol=INVARIANCE_ATOL)

    def test_e_ref_itself_does_scale(self, real_micro, real_result, distance_um,
                                     reach):
        """The control: E_ref is an ABSOLUTE exposure and must move with the
        dose. If it did not, the invariance above would be measuring nothing."""
        scaled = compute_reach(real_micro, _scaled_result(real_result, self.FACTOR),
                               distance_um=distance_um)
        assert scaled["E_ref_uM_s"] == pytest.approx(
            reach["E_ref_uM_s"] * self.FACTOR, rel=1e-12)

    def test_a_downward_rescale_moves_nothing_either(
            self, real_micro, real_result, distance_um, reach):
        """1/5000 as well as 5000: the cancellation is not a large-number
        artefact."""
        scaled = compute_reach(real_micro,
                               _scaled_result(real_result, 1.0 / self.FACTOR),
                               distance_um=distance_um)
        np.testing.assert_allclose(
            scaled["curve"]["under_reached_fraction"],
            reach["curve"]["under_reached_fraction"],
            rtol=0.0, atol=INVARIANCE_ATOL)
        np.testing.assert_allclose(
            scaled["falloff"]["normalised_mean_auc"],
            reach["falloff"]["normalised_mean_auc"],
            rtol=0.0, atol=INVARIANCE_ATOL)


# --------------------------------------------------------------------------- #
# The curve is a distribution function
# --------------------------------------------------------------------------- #
class TestCurveShape:

    def test_ladder_is_the_documented_one(self, reach):
        assert reach["curve"]["x"] == [float(v) for v in REACH_LADDER]

    def test_F_is_monotone_non_decreasing_in_x(self, reach):
        F = np.asarray(reach["curve"]["under_reached_fraction"])
        assert np.all(np.diff(F) >= 0.0), F

    def test_F_is_a_fraction(self, reach):
        F = np.asarray(reach["curve"]["under_reached_fraction"])
        assert np.all(F >= 0.0) and np.all(F <= 1.0), F

    def test_headline_is_picked_off_the_ladder_not_recomputed(self, reach):
        i = reach["curve"]["x"].index(HEADLINE_X)
        assert reach["headline_x"] == HEADLINE_X
        # Exact equality: the headline IS the ladder entry, not a second
        # evaluation of the same definition.
        assert (reach["headline_under_reached_fraction"]
                == reach["curve"]["under_reached_fraction"][i])

    def test_counts_are_consistent(self, reach):
        assert 0 < reach["n_viable"] <= reach["n_tumor"]


# --------------------------------------------------------------------------- #
# E_ref is the SAME number the penetration length was fitted from
# --------------------------------------------------------------------------- #
class TestPerivascularReference:

    def test_e_ref_is_penetration_length_fits_perivascular_mean(
            self, real_micro, real_result, distance_um, reach):
        pen = penetration_length_fit(real_result.fields["auc"], real_micro,
                                     distance_um=distance_um)
        # Equality, not approximation: E_ref is REUSED from this dict, never
        # recomputed, so the reach number and the published penetration length
        # cannot drift apart.
        assert reach["E_ref_uM_s"] == pen["drug"]["perivascular_mean"]

    def test_passing_the_precomputed_fit_changes_nothing(
            self, real_micro, real_result, distance_um, reach):
        pen = penetration_length_fit(real_result.fields["auc"], real_micro,
                                     distance_um=distance_um)
        passed = compute_reach(real_micro, real_result, distance_um=distance_um,
                               penetration=pen)
        assert passed["E_ref_uM_s"] == reach["E_ref_uM_s"]
        assert (passed["curve"]["under_reached_fraction"]
                == reach["curve"]["under_reached_fraction"])

    def test_first_falloff_bin_is_the_reference_itself(self, reach):
        assert reach["falloff"]["normalised_mean_auc"][0] == 1.0

    def test_falloff_bins_all_meet_the_minimum_count(self, reach):
        counts = reach["falloff"]["count"]
        assert counts and min(counts) >= reach["falloff"]["min_bin_count"]

    def test_falloff_arrays_are_parallel(self, reach):
        f = reach["falloff"]
        assert (len(f["distance_um"]) == len(f["normalised_mean_auc"])
                == len(f["count"]))


# --------------------------------------------------------------------------- #
# The threshold is a CHOICE, and the payload says so
# --------------------------------------------------------------------------- #
class TestReportingConvention:

    def test_the_convention_sentence_travels_with_the_number(self, reach):
        assert reach["convention"] == CONVENTION_NOTE
        low = CONVENTION_NOTE.lower()
        assert "reporting convention" in low
        assert "not a physical constant" in low
        assert "0.1" in CONVENTION_NOTE

    def test_the_whole_ladder_is_returned_beside_the_headline(self, reach):
        # The curve is returned precisely so the threshold is visibly one point
        # on it; a payload with only the headline would rebuild 1x IC50.
        assert len(reach["curve"]["x"]) == len(REACH_LADDER) > 1

    def test_the_scale_invariance_disclosure_states_both_strengths(self, reach):
        s = reach["scale_invariance"]
        assert "1.41e-14" in s and "pk_ranking_robustness.md" in s
        assert "1e-9" in s and "capacity_ranking_robustness.md" in s
        assert "descriptor_falsification.md" in s
        # CONVERTED 2026-08-21, and made stricter. This line was
        # `assert "NOT exact over a finite window" in s`, which pinned the
        # correct SIGN of the capacity claim but let it stay qualitative. The
        # capacity residual has since been MEASURED on the flagship section, so
        # the disclosure now has to carry the numbers rather than the adjective:
        # a qualitative "not exact" is what let this note ship claiming a regime
        # the product is not in.
        assert "BOUNDED RESIDUAL" in s          # not an invariance, at this window
        assert "REGIME CLAIM WAS NOT" in s      # correct physics, wrong regime
        # The measured residual, at both readings, because they differ by an
        # order of magnitude and the panel must not present them as equal.
        assert "0.0011" in s                    # the ladder, over R = 1..39
        assert "0.3772" in s                    # the falloff SHAPE, at R = 39
        assert "0.80" in s                      # the ladder past the bound
        # The window is NOT converged, with the decay law that says so.
        assert "0.01302" in s and "281/T" in s

    def test_the_disclosure_names_the_permeability_blindness(self, reach):
        """The map is EXACTLY invariant to bbb_penetrance, so it cannot see a
        barrier at all. That is the price of the scale-freedom that makes it
        trustworthy, and a reader asking "does my compound get into the brain"
        must not read this figure as the answer."""
        s = reach["scale_invariance"]
        assert "bbb_penetrance" in s
        assert "sqrt(D_eff/k)" in s
        low = s.lower()
        assert "permeability" in low or "barrier" in low

    def test_units_are_declared(self, reach):
        assert reach["units"]["E_ref"] == "uM*s"
        assert "dimensionless" in reach["units"]["reach_ratio"]


# --------------------------------------------------------------------------- #
# REFUSAL: an undefined quantity is never given a number
# --------------------------------------------------------------------------- #
def _assert_refused(out):
    assert out["null_reason"], "an undefined case returned no reason"
    assert isinstance(out["null_reason"], str) and len(out["null_reason"]) > 20
    assert out["headline_under_reached_fraction"] is None
    assert out["curve"] is None
    assert out["falloff"] is None
    # The convention still travels, so a consumer rendering the refusal shows
    # the same sentence it would have shown beside the number.
    assert out["convention"] == CONVENTION_NOTE


class TestRefusals:

    def test_no_tumour_cells(self, real_micro, real_result, distance_um):
        out = compute_reach(_retyped(real_micro, CellType.ENDOTHELIAL),
                            real_result, distance_um=distance_um)
        _assert_refused(out)
        assert "TUMOR cells" in out["null_reason"]
        assert out["n_tumor"] == 0
        assert out["E_ref_uM_s"] is None

    def test_no_viable_tumour_cells(self, real_micro, real_result, distance_um):
        """Every tumour cell below the viable oxygen band: the population the
        fraction is taken OVER is empty, and it does not fall back to all
        tumour cells."""
        starved = copy.copy(real_result)
        starved.fields = dict(real_result.fields)
        starved.fields["oxygen"] = np.zeros_like(
            np.asarray(real_result.fields["oxygen"], dtype=np.float64))
        out = compute_reach(real_micro, starved, distance_um=distance_um)
        _assert_refused(out)
        assert "VIABLE" in out["null_reason"]
        assert out["n_viable"] == 0
        assert out["n_tumor"] > 0

    def test_e_ref_not_positive(self, real_micro, real_result, distance_um):
        """A zero exposure field has a zero perivascular reference, and every
        reach ratio would be a division by it."""
        zeroed = _scaled_result(real_result, 0.0)
        out = compute_reach(real_micro, zeroed, distance_um=distance_um)
        _assert_refused(out)
        assert "E_ref" in out["null_reason"]

    def test_missing_oxygen_field(self, real_micro, real_result, distance_um):
        blind = copy.copy(real_result)
        blind.fields = {k: v for k, v in real_result.fields.items()
                        if k != "oxygen"}
        out = compute_reach(real_micro, blind, distance_um=distance_um)
        _assert_refused(out)
        assert "oxygen" in out["null_reason"]

    def test_headline_off_the_ladder_is_refused_not_interpolated(
            self, real_micro, real_result, distance_um):
        """The headline is read OFF the ladder. Hand a ladder without it and the
        answer is a refusal, never a value computed some other way."""
        out = compute_reach(real_micro, real_result, distance_um=distance_um,
                            ladder=(0.2, 0.5))
        _assert_refused(out)
        assert "ladder" in out["null_reason"]


# --------------------------------------------------------------------------- #
# The viable band is metric 5's band, not a second opinion
# --------------------------------------------------------------------------- #
def test_viable_band_matches_delivery_metrics(real_micro, real_result, reach):
    from tmesim.validation.delivery_metrics import (
        hypoxic_vs_normoxic_exposure_ratio,
    )

    m5 = hypoxic_vs_normoxic_exposure_ratio(real_micro, real_result)
    assert reach["n_viable"] == m5["n_viable"]
    assert reach["n_tumor"] == int(m5["n_hypoxic"] + m5["n_normoxic"])
    assert reach["hypoxia_threshold"] == OxygenConfig().hypoxia_threshold


# --------------------------------------------------------------------------- #
# The wire contract the Studio consumes
# --------------------------------------------------------------------------- #
def test_result_is_json_serialisable(reach):
    import json

    from app.server.jobs import jsonable

    doc = json.loads(json.dumps(jsonable(reach)))
    assert doc["headline_x"] == HEADLINE_X
    assert doc["curve"]["under_reached_fraction"][
        doc["curve"]["x"].index(HEADLINE_X)] == \
        doc["headline_under_reached_fraction"]


# --------------------------------------------------------------------------- #
# The run pipeline: reach travels with a run exactly as extra_metrics does
# --------------------------------------------------------------------------- #
SECTION = "xenium_319"
RUN_DX_UM = 30.0


def _section_available() -> bool:
    from app.server import cache as cache_mod

    return any(s["id"] == SECTION and s["available"]
               for s in cache_mod.list_sections())


def _run_spec(**overrides):
    spec = {
        "section": SECTION,
        "dx_um": RUN_DX_UM,
        "drug": "temozolomide",
        "sensitivity_mode": "uniform",
        "vascular_mode": "discrete",
        "oxygen_mode": "pde",
        "mgmt": None,
        "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
               "n_doses": 1, "dose_interval_s": 86400.0},
        "total_time_s": WINDOW_S,
    }
    spec.update(overrides)
    return spec


def _drive(manager, spec, timeout=300.0):
    import time as _time

    job = manager.submit(spec)
    deadline = _time.time() + timeout
    while job.status not in ("done", "failed", "cancelled"):
        if _time.time() > deadline:
            raise AssertionError(f"run {job.run_id} did not finish")
        _time.sleep(0.1)
    assert job.status == "done", job.error
    return job


@pytest.mark.skipif(not _section_available(),
                    reason=f"section {SECTION!r} is not cached in data_cache/")
class TestRunWiring:
    """``reach`` is threaded exactly the way ``extra_metrics`` is: onto the run
    body, out through ``result_view()``, into ``result.json``, and back in
    through ``_restore``."""

    def test_reach_is_on_the_run_and_survives_a_restart(self, tmp_path):
        from app.server.jobs import JobManager

        root = str(tmp_path / "runs")
        manager = JobManager(root=root)
        try:
            job = _drive(manager, _run_spec())
            view = job.result_view()
            reach = view["reach"]
            assert reach["null_reason"] is None, reach["null_reason"]
            assert 0.0 <= reach["headline_under_reached_fraction"] <= 1.0
            assert reach["curve"]["x"] == [float(v) for v in REACH_LADDER]
            assert reach["convention"] == CONVENTION_NOTE
            run_id = job.run_id
        finally:
            manager.shutdown()

        # A brand-new manager over the same directory = a server restart.
        restarted = JobManager(root=root, start_worker=False)
        try:
            back = restarted.get(run_id)
            assert back is not None
            assert back.result_view()["reach"] == reach
        finally:
            restarted.shutdown()

    def test_a_reconstruction_only_run_refuses_rather_than_omitting(
            self, tmp_path):
        """No drug means no exposure field. The key is still present, carrying
        the reason — a consumer never has to tell 'absent' from 'undefined'."""
        from app.server.jobs import JobManager

        manager = JobManager(root=str(tmp_path / "runs"))
        try:
            job = _drive(manager, _run_spec(drug=None))
            reach = job.result_view()["reach"]
            _assert_refused(reach)
            assert "drug=null" in reach["null_reason"]
        finally:
            manager.shutdown()
