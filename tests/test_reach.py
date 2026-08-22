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
    FINITE_WINDOW_NOTE,
    HEADLINE_X,
    REACH_LADDER,
    STEADY_SCALE_INVARIANCE_NOTE,
    STEADY_STATE_NOTE,
    WINDOW_MODES,
    _advection_matrix,
    compute_reach,
    steady_state_concentration,
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
#:
#: THIS IS A FINITE-WINDOW CONSTANT and it now says so. Every test that uses it
#: is testing ``window_mode = "finite_window"``; the steady-state tests below
#: take no window at all, which is the point of them.
WINDOW_S = 600.0

#: How far a march has to run before ``C(T)`` may be compared with ``C_ss``.
#:
#: The field converges at the reaction rate ``k = clearance + binding``, so the
#: error falls like ``exp(-k*T)`` and this is expressed in TIME CONSTANTS rather
#: than in seconds — a fixed number of seconds would be a different amount of
#: convergence for every compound, which is exactly the confusion the whole
#: change is about. 25 time constants predicts ``exp(-25) = 1.4e-11``; the
#: assertion is loosened to 1e-6 because the march also carries its own
#: backward-Euler time-discretisation error at ``dt = 10 s``, which does NOT
#: vanish with T and is the larger of the two terms. MEASURED on this fixture:
#: the agreement reached is reported by the test itself when it fails.
STEADY_MARCH_TAUS = 25.0
STEADY_MARCH_ATOL_REL = 1e-6

#: Capacity-invariance tolerance for the STEADY map, relative to the falloff's
#: own range. The steady problem does not contain R at all, so the only thing
#: that moves is float64 round-off in forming ``D/R``, ``k/R`` and ``S/R`` and
#: in the sparse factorisation of a matrix whose entries have been scaled by R.
#: That is why this is 1e-9 rather than 0: the algebra is exact, the arithmetic
#: is not. 1e-9 is the same strength ``results/capacity_ranking_robustness.md``
#: sec. 2b states its own steady-state verification at.
CAPACITY_ATOL = 1e-9


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
# THE STEADY-STATE OBJECT
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def velocity(real_micro):
    """The interstitial velocity, from the microenvironment alone."""
    from tmesim.physics.pressure import solve_ifp

    return solve_ifp(real_micro)[1]


@pytest.fixture(scope="module")
def css(real_micro, real_drug, velocity):
    """``C_ss`` for the real section and the real drug, convection ON."""
    from tmesim.treatment.drugs import delivery_source

    return steady_state_concentration(
        real_micro, real_drug, delivery_source(real_micro, real_drug, 1.0),
        velocity=velocity)


@pytest.fixture(scope="module")
def steady_reach(real_micro, real_result, distance_um, real_drug, velocity):
    out = compute_reach(real_micro, real_result, distance_um=distance_um,
                        drug=real_drug, velocity=velocity, convection=True)
    assert out["null_reason"] is None, out["null_reason"]
    return out


class TestSteadyStateOperator:
    """The one piece of the engine's discretisation re-expressed here.

    ``steady_state_concentration`` IMPORTS ``_build_diffusion_operator`` and
    ``_effective_diffusivity``, so those cannot drift. The upwind advection
    term is the exception: the engine applies it as a function of ``c`` and the
    steady solve needs it as a matrix. That re-expression is the only place a
    silent disagreement could live, so it is measured against the engine's own
    function rather than reasoned about.
    """

    def test_the_matrix_applies_the_engines_own_advection(self):
        from tmesim.physics.diffusion import _advective_divergence

        rng = np.random.default_rng(0)
        worst = 0.0
        for _ in range(8):
            ny, nx = int(rng.integers(4, 22)), int(rng.integers(4, 22))
            dx = float(rng.uniform(5.0, 40.0))
            tissue = rng.random((ny, nx)) > 0.25
            vx, vy = rng.normal(size=(ny, nx)), rng.normal(size=(ny, nx))
            c = rng.random((ny, nx))
            c[~tissue] = 0.0
            direct = _advective_divergence(c, vx, vy, tissue, dx)
            viamat = (_advection_matrix(vx, vy, tissue, dx)
                      @ c.ravel()).reshape(ny, nx)
            scale = max(float(np.abs(direct).max()), 1e-300)
            worst = max(worst, float(np.abs(direct - viamat).max()) / scale)
        # MEASURED at 2.140e-16 over 6 cases in
        # results/steady_state_reach.md sec. 0; 1e-13 is ~470x that, tight
        # enough that any genuine disagreement in signs, faces or the upwind
        # switch fails immediately.
        assert worst < 1e-13, worst

    def test_a_zero_velocity_field_contributes_nothing(self):
        z = np.zeros((5, 6))
        m = _advection_matrix(z, z, np.ones((5, 6), dtype=bool), 10.0)
        assert m.nnz == 0 or float(np.abs(m.toarray()).max()) == 0.0


class TestSteadyStateField:
    """``C_ss`` is the fixed point of the shipped stepper, measured not argued."""

    def test_it_is_where_a_long_march_goes(self, real_micro, real_drug, css,
                                           velocity):
        """March far past equilibrium and land on the direct solve.

        This is the check that matters: the algebra says the direct solve is
        the stepper's fixed point, and a fixed point of the WRONG operator
        would still look like a perfectly good linear solve. Only running the
        actual stepper can tell the difference.
        """
        from tmesim.core import PDEConfig
        from tmesim.physics.diffusion import simulate_drug_transport
        from tmesim.treatment.drugs import delivery_source

        k = float(real_drug.clearance) + float(real_drug.binding)
        total = STEADY_MARCH_TAUS / k
        tr = simulate_drug_transport(
            real_micro, real_drug, delivery_source(real_micro, real_drug, 1.0),
            PDEConfig(dt=10.0, total_time=total, convection=True,
                      save_every=10 ** 9),
            velocity=velocity)
        tissue = np.asarray(real_micro.domain.tissue_mask, dtype=bool)
        scale = float(css[tissue].max())
        err = float(np.abs(tr.final_concentration[tissue]
                           - css[tissue]).max()) / scale
        assert err < STEADY_MARCH_ATOL_REL, (
            f"the direct steady solve and a {STEADY_MARCH_TAUS:g}-time-constant "
            f"march disagree by {err:.3e} of max C_ss")

    def test_the_unclamped_solve_is_still_non_negative(self, css):
        """The stepper clamps C at zero every step; the direct solve clamps
        nothing. If the discrete steady state were not non-negative the two
        routes could not agree, so this is not cosmetic."""
        assert float(css.min()) >= 0.0, float(css.min())

    def test_it_is_exactly_linear_in_the_source(self, real_micro, real_drug,
                                                velocity, css):
        """Every scalar in the source — the plasma level, ``bbb_penetrance`` —
        multiplies the whole field, which is why the map needs no dose. Solved
        with a unit source in production, so this is the property that licenses
        that."""
        from tmesim.treatment.drugs import delivery_source

        f = 3721.0
        scaled = steady_state_concentration(
            real_micro, real_drug,
            delivery_source(real_micro, real_drug, 1.0) * f, velocity=velocity)
        ref = css * f
        denom = max(float(np.abs(ref).max()), 1e-300)
        assert float(np.abs(scaled - ref).max()) / denom < 1e-12

    def test_solve_ifp_from_micro_reproduces_the_engines_velocity(
            self, real_micro, velocity):
        """``compute_reach`` rebuilds the velocity from ``micro`` when it is not
        given one, and that has to be the SAME array the engine's transport call
        received. ``solve_ifp`` reads the grid, ``vascular_density``,
        ``ecm_density`` and the cells' ``cell_type``; the engine sets the oxygen
        field and updates cell STATE before calling it, and neither is an input.
        Bit-identical, not close."""
        from tmesim.core import OxygenConfig
        from tmesim.cells.agents import update_cell_states
        from tmesim.physics.oxygen import solve_oxygen
        from tmesim.physics.pressure import solve_ifp
        from tmesim.sim.engine import _working_copy

        work = _working_copy(real_micro)
        ox = solve_oxygen(work, OxygenConfig())
        work.fields.oxygen = ox
        update_cell_states(work, ox, OxygenConfig())
        assert np.array_equal(solve_ifp(work)[1], velocity)


class TestSteadyCapacityInvariance:
    """The measurement the whole change rests on.

    A capacity term ``R`` multiplies the TIME DERIVATIVE and nothing else, so
    it is absent from the steady problem. That is a prediction of exactly zero
    movement, and predictions of zero are the ones worth testing — the
    finite-window falloff moves 0.3772 at ``R = 39``
    (``results/capacity_ranking_robustness.md``), so there is a large number
    here to fail to reproduce.

    ``R`` is emulated exactly as ``experiments/capacity_ranking_robustness.py``
    emulates it — ``D/R``, ``binding/R``, ``clearance/R``, ``bbb_penetrance/R``
    — and **convection is OFF**, which is not a convenience. Multiply the steady
    equation of the emulated system through by ``R``:

    ``div(D grad C) - R*div(v C) + S - k*C = 0``

    The advection term does NOT lose its ``R``, because the emulation cannot
    divide the velocity field (``capacity_ranking_robustness.md`` sec. 2(d): it
    is built inside ``simulate_response`` and is unreachable from outside). So
    with convection on, the published emulation is not a capacity change at
    all — it is a capacity change AND a factor-``R`` boost to the Péclet number,
    and the map moves accordingly. MEASURED: 1.044e-01 at ``R = 39`` on this
    fixture with convection on, against the ``< 1e-9`` below with it off.

    ``test_dividing_the_velocity_too_restores_the_invariance`` closes the loop:
    scale ``v`` by ``1/R`` as well and the invariance returns with convection
    on, which is what identifies the residual as the emulation's artefact
    rather than the map's. In the real physics ``R`` multiplies ``dC/dt`` and
    the steady problem has no ``dC/dt``, convection or not.
    """

    RS = (39.0, 4184.0)

    #: What the published emulation does to the steady map when convection is
    #: left ON. Not a tolerance — a measured artefact of the emulation, pinned
    #: here so that it cannot be mistaken for a property of the map.
    #:
    #: This is THIS FIXTURE's value (dx = 80 um). On the flagship at dx = 12 the
    #: same sweep moves the normalised falloff by 1.226 at R = 39 and 4.380e+03
    #: at R = 4184 (results/steady_state_reach.md sec. 5). The artefact is
    #: strongly grid-dependent because it is a Peclet-number boost, which is
    #: part of why it is easy to walk into and hard to recognise.
    EMULATION_WITH_CONVECTION_MOVE = 1.044e-01

    @staticmethod
    def _retard(dp, R):
        import dataclasses

        return dataclasses.replace(
            dp, diffusivity=dp.diffusivity / R, binding=dp.binding / R,
            clearance=dp.clearance / R, bbb_penetrance=dp.bbb_penetrance / R)

    def _reach_at(self, real_micro, real_result, distance_um, dp,
                  convection=False, velocity=None):
        out = compute_reach(real_micro, real_result, distance_um=distance_um,
                            drug=dp, velocity=velocity, convection=convection)
        assert out["null_reason"] is None, out["null_reason"]
        return out

    def test_the_published_emulation_is_not_capacity_when_convection_is_on(
            self, real_micro, real_result, distance_um, real_drug, velocity):
        """The counter-case, kept as an assertion rather than as a remark.

        If someone later runs the capacity sweep against the steady map with
        convection left on, they will get a large residual and conclude the
        invariance claim is false. It is not — the emulation is. This pins the
        size of that trap.
        """
        base = self._reach_at(real_micro, real_result, distance_um, real_drug,
                              convection=True, velocity=velocity)
        got = self._reach_at(real_micro, real_result, distance_um,
                             self._retard(real_drug, 39.0),
                             convection=True, velocity=velocity)
        move = float(np.abs(
            np.asarray(got["curve"]["under_reached_fraction"])
            - np.asarray(base["curve"]["under_reached_fraction"])).max())
        assert move == pytest.approx(self.EMULATION_WITH_CONVECTION_MOVE,
                                     abs=5e-4), move

    def test_dividing_the_velocity_too_restores_the_invariance(
            self, real_micro, real_result, distance_um, real_drug, velocity):
        """A FAITHFUL capacity emulation with convection on: ``v/R`` as well.

        ``R*dC/dt = div(D grad C) - div(v C) + S - k*C`` divided through by
        ``R`` scales the velocity along with everything else. Do that and the
        residual above collapses to round-off, which is what identifies it as
        the emulation's artefact.
        """
        R = 39.0
        base = self._reach_at(real_micro, real_result, distance_um, real_drug,
                              convection=True, velocity=velocity)
        got = self._reach_at(real_micro, real_result, distance_um,
                             self._retard(real_drug, R),
                             convection=True, velocity=velocity / R)
        for key, arr in (("curve", "under_reached_fraction"),
                         ("falloff", "normalised_mean_auc")):
            move = float(np.abs(np.asarray(got[key][arr])
                                - np.asarray(base[key][arr])).max())
            assert move < CAPACITY_ATOL, f"{arr} moved {move:.3e} at R = {R:g}"

    def test_the_falloff_shape_does_not_move(self, real_micro, real_result,
                                             distance_um, real_drug):
        base = self._reach_at(real_micro, real_result, distance_um, real_drug)
        b = np.asarray(base["falloff"]["normalised_mean_auc"])
        for R in self.RS:
            got = self._reach_at(real_micro, real_result, distance_um,
                                 self._retard(real_drug, R))
            g = np.asarray(got["falloff"]["normalised_mean_auc"])
            assert g.shape == b.shape
            move = float(np.abs(g - b).max())
            assert move < CAPACITY_ATOL, (
                f"R = {R:g} moved the STEADY falloff by {move:.3e}; the steady "
                f"problem does not contain R, so this is not a tolerance "
                f"question")

    def test_the_under_reached_ladder_does_not_move(self, real_micro,
                                                    real_result, distance_um,
                                                    real_drug):
        base = self._reach_at(real_micro, real_result, distance_um, real_drug)
        b = np.asarray(base["curve"]["under_reached_fraction"])
        for R in self.RS:
            got = self._reach_at(real_micro, real_result, distance_um,
                                 self._retard(real_drug, R))
            move = float(np.abs(
                np.asarray(got["curve"]["under_reached_fraction"]) - b).max())
            assert move < CAPACITY_ATOL, f"R = {R:g} moved the ladder {move:.3e}"

    def test_the_finite_window_map_DOES_move(self, real_micro, real_result,
                                             distance_um, real_drug):
        """The control. If the finite-window arm did not move either, the test
        above would be measuring nothing — it would be consistent with the
        capacity emulation simply not being applied."""
        from tmesim.core import PDEConfig, SimConfig
        from tmesim.sim.engine import simulate_response

        base = compute_reach(real_micro, real_result, distance_um=distance_um)
        b = np.asarray(base["falloff"]["normalised_mean_auc"])
        res = simulate_response(
            real_micro,
            SimConfig(drug=self._retard(real_drug, 39.0),
                      pde=PDEConfig(dt=10.0, total_time=WINDOW_S,
                                    convection=True, save_every=500),
                      seed=0))
        got = compute_reach(real_micro, res, distance_um=distance_um)
        move = float(np.abs(
            np.asarray(got["falloff"]["normalised_mean_auc"]) - b).max())
        assert move > 1e3 * CAPACITY_ATOL, (
            "the finite-window falloff did not move under R = 39, so the "
            f"capacity emulation is not reaching the solve (moved {move:.3e})")


class TestSteadyDegenerateTissue:
    """Where a one-shot direct solve returns a number it should not.

    A time march on a zero source produces zeros and the refusal rule catches
    it. A direct solve on a degenerate system produces whatever the linear
    solver produces, and the failure mode is a PLAUSIBLE FINITE NUMBER rather
    than an error. These are the cases where that could happen.
    """

    def test_absent_vasculature_refuses_and_returns_no_number(
            self, real_micro, real_result, distance_um, real_drug, velocity):
        base = np.array(real_micro.fields.vascular_density, copy=True)
        try:
            real_micro.fields.vascular_density = base * 0.0
            out = compute_reach(real_micro, real_result,
                                distance_um=distance_um, drug=real_drug,
                                velocity=velocity, convection=True)
        finally:
            real_micro.fields.vascular_density = base
        _assert_refused(out)
        assert "E_ref" in out["null_reason"]

    def test_sparse_vasculature_does_NOT_refuse_and_that_is_correct(
            self, real_micro, real_result, distance_um, real_drug, velocity):
        """Kept as an assertion because it is counter-intuitive.

        Scaling the vessel density down by 1e-3 leaves a well-defined map: the
        field is exactly linear in the source, so the RATIOS this module
        publishes are unchanged. Refusing here would be discarding a number
        that is defined, which the module docstring calls the opposite of the
        refusal rule. Measured on the flagship at dx = 30 the same holds
        geometrically — 4 endothelial cells still gives a usable first bin and
        a real perivascular mean (results/steady_state_reach.md sec. 9c).
        """
        base = np.array(real_micro.fields.vascular_density, copy=True)
        ref = compute_reach(real_micro, real_result, distance_um=distance_um,
                            drug=real_drug, velocity=velocity, convection=True)
        try:
            real_micro.fields.vascular_density = base * 1e-3
            out = compute_reach(real_micro, real_result,
                                distance_um=distance_um, drug=real_drug,
                                velocity=velocity, convection=True)
        finally:
            real_micro.fields.vascular_density = base
        assert out["null_reason"] is None
        np.testing.assert_allclose(out["falloff"]["normalised_mean_auc"],
                                   ref["falloff"]["normalised_mean_auc"],
                                   rtol=0.0, atol=1e-12)
        assert (out["curve"]["under_reached_fraction"]
                == ref["curve"]["under_reached_fraction"])

    def test_where_the_denormal_cliff_IS(self, real_micro, real_result,
                                         distance_um, real_drug, velocity):
        """The one silent-wrong regime this map has, PINNED rather than guarded.

        Scaling the source until ``E_ref`` goes denormal eventually returns a
        plausible finite falloff that is simply wrong, with no refusal. It is
        deliberately not guarded: a positivity floor above ``E_ref > 0`` would
        be a threshold chosen to defend a case that cannot occur, and the
        realised ``E_ref`` on real tissue is O(1-10) because the plasma scalar
        — the one axis that could be tiny — divides out by construction.

        So the cliff's LOCATION is asserted instead. This survives someone
        later changing what feeds ``E_ref``; a magic floor would be inherited
        and never re-derived. MEASURED on the flagship: exact to 2.220e-16 at
        1e-300, 3.804e-13 at 1e-310, 4.745e-08 at 1e-315, 5.050e-02 at 1e-320
        (results/steady_state_reach.md sec. 9d).
        """
        base = np.array(real_micro.fields.vascular_density, copy=True)
        ref = compute_reach(real_micro, real_result, distance_um=distance_um,
                            drug=real_drug, velocity=velocity, convection=True)
        rf = np.asarray(ref["falloff"]["normalised_mean_auc"])
        try:
            # Still exact 250 orders down: the ratio is scale-free until the
            # floats stop being able to represent the scale.
            real_micro.fields.vascular_density = base * 1e-250
            far = compute_reach(real_micro, real_result,
                                distance_um=distance_um, drug=real_drug,
                                velocity=velocity, convection=True)
            assert far["null_reason"] is None
            assert float(np.abs(
                np.asarray(far["falloff"]["normalised_mean_auc"]) - rf).max()
            ) < 1e-15
            # And broken by 1e-320, silently. If this ever STOPS being broken,
            # the cliff has moved and the sentence in the results file is stale.
            real_micro.fields.vascular_density = base * 1e-320
            gone = compute_reach(real_micro, real_result,
                                 distance_um=distance_um, drug=real_drug,
                                 velocity=velocity, convection=True)
        finally:
            real_micro.fields.vascular_density = base
        if gone["null_reason"] is None:
            moved = float(np.abs(
                np.asarray(gone["falloff"]["normalised_mean_auc"]) - rf).max())
            assert moved > 1e-3, (
                "the denormal cliff has MOVED: scaling the source to 1e-320 no "
                f"longer breaks the falloff (moved {moved:.3e}). "
                "results/steady_state_reach.md sec. 9d locates it here and is "
                "now stale.")


class TestWindowMode:
    """A payload always says WHICH OBJECT it is. That is the defect this repo
    spent the week fixing: records that were true about a configuration nobody
    wrote down."""

    def test_the_two_modes_are_labelled(self, reach, steady_reach):
        assert reach["window_mode"] == "finite_window"
        assert steady_reach["window_mode"] == "steady_state"
        assert set(WINDOW_MODES) == {"steady_state", "finite_window"}

    def test_the_finite_window_payload_carries_its_own_window(self, reach):
        assert reach["window_s"] == pytest.approx(WINDOW_S)

    def test_the_steady_payload_has_no_window_at_all(self, steady_reach):
        # Not zero, not a default: there is no window in the object.
        assert steady_reach["window_s"] is None

    def test_the_units_change_because_the_quantity_changes(self, reach,
                                                           steady_reach):
        """``E_ref`` is an EXPOSURE in one mode and a CONCENTRATION in the
        other. The wire key is still ``E_ref_uM_s`` for compatibility with
        ``app/web`` and ``app/server/report.py``, so the payload has to carry
        the true unit somewhere a consumer can read it."""
        assert reach["E_ref_units"] == "uM*s"
        assert reach["units"]["E_ref"] == "uM*s"
        assert steady_reach["E_ref_units"] == "uM"
        assert steady_reach["units"]["E_ref"] == "uM"

    def test_the_mode_is_stated_in_words_not_only_in_a_field(self, reach,
                                                            steady_reach):
        assert reach["window_note"] == FINITE_WINDOW_NOTE
        assert steady_reach["window_note"] == STEADY_STATE_NOTE
        assert "FINITE WINDOW" in FINITE_WINDOW_NOTE
        assert "STEADY STATE" in STEADY_STATE_NOTE
        # The reframing itself, in the sentence a reader meets.
        assert "GEOMETRY PERMITS" in STEADY_STATE_NOTE
        # And what it costs, in the same breath.
        assert "DIFFERENT QUESTION" in STEADY_STATE_NOTE

    def test_a_caller_supplied_fit_is_discarded_under_steady_state(
            self, real_micro, real_result, distance_um, real_drug, velocity):
        """The production caller passes a ``penetration_length_fit`` of the
        finite-window AUC profile. Normalising the STEADY field by a reference
        taken off a different field would be silently wrong, so it is dropped
        and the payload says it was."""
        pen = penetration_length_fit(real_result.fields["auc"], real_micro,
                                     distance_um=distance_um)
        out = compute_reach(real_micro, real_result, distance_um=distance_um,
                            penetration=pen, drug=real_drug,
                            velocity=velocity, convection=True)
        assert out["null_reason"] is None
        assert out["penetration_reused"] is False
        assert out["E_ref_uM_s"] != pen["drug"]["perivascular_mean"]

    def test_the_finite_window_mode_still_reuses_it(self, real_micro,
                                                    real_result, distance_um):
        pen = penetration_length_fit(real_result.fields["auc"], real_micro,
                                     distance_um=distance_um)
        out = compute_reach(real_micro, real_result, distance_um=distance_um,
                            penetration=pen)
        assert out["penetration_reused"] is True
        assert out["E_ref_uM_s"] == pen["drug"]["perivascular_mean"]

    def test_omitting_convection_refuses_rather_than_defaulting(
            self, real_micro, real_result, distance_um, real_drug):
        """``convection`` must match the run's ``PDEConfig.convection`` and is
        not carried on a ``SimulationResult``, so there is nothing to infer it
        from. A default of ``True`` would agree with ``jobs.PDE_CONVECTION``
        today and stop agreeing silently the moment either moved — which is the
        same shape as the emulation trap two tests up. So it REFUSES."""
        out = compute_reach(real_micro, real_result, distance_um=distance_um,
                            drug=real_drug)
        _assert_refused(out)
        assert "convection" in out["null_reason"]
        assert out["window_mode"] == "steady_state"

    def test_the_finite_window_mode_does_not_require_it(self, reach):
        """It builds no field, so there is nothing for convection to mean."""
        assert reach["null_reason"] is None

    def test_a_refusal_still_says_which_object_failed(self, real_micro,
                                                      real_result, distance_um,
                                                      real_drug):
        out = compute_reach(_retyped(real_micro, CellType.ENDOTHELIAL),
                            real_result, distance_um=distance_um,
                            drug=real_drug, convection=True)
        _assert_refused(out)
        assert out["window_mode"] == "steady_state"
        assert out["window_note"] == STEADY_STATE_NOTE
        assert out["scale_invariance"] == STEADY_SCALE_INVARIANCE_NOTE


class TestSteadyDisclosure:
    """The third version of this note in two days. It has to be RIGHT rather
    than merely current, so every clause that quotes a measurement is pinned
    at the magnitude, here, at the site that quotes it."""

    def test_it_states_the_capacity_claim_as_exact_and_says_why(self,
                                                               steady_reach):
        s = steady_reach["scale_invariance"]
        assert "THE MAP IS EXACTLY CAPACITY-INVARIANT" in s
        # The REASON, which is what makes it checkable rather than assertable:
        # R multiplies the time derivative, and there is no time derivative.
        assert "TIME DERIVATIVE AND NOTHING ELSE" in s
        assert "1e-9" in s and "capacity_ranking_robustness.md" in s
        # The magnitudes, pinned at the site that quotes them.
        assert "0.000e+00" in s and "4.452e-14" in s and "0.3731" in s

    def test_the_two_capacity_sentences_name_their_own_subjects(self,
                                                               steady_reach):
        """(2) is about the MODEL and (2b) is about an INSTRUMENT. Both are
        true at once, and a reader who did not follow the argument must not be
        able to read the second as a walk-back of the first — which is what
        happens if either sentence leaves its subject implicit."""
        s = steady_reach["scale_invariance"]
        assert "A STATEMENT ABOUT THE MAP" in s
        assert "A STATEMENT ABOUT AN INSTRUMENT, NOT A WALK-BACK" in s
        assert "THE PUBLISHED CAPACITY SWEEP IS NOT A FAITHFUL CAPACITY CHANGE" in s
        # The trap's size, and the fact the shipping config walks into it.
        # The flagship magnitudes; 1.044e-01 below is this COARSE fixture's
        # own value for the same sweep, and the two differ by an order of
        # magnitude because the Peclet boost is grid-dependent — which is
        # itself a reason the trap is easy to walk into.
        assert "1.226" in s and "4.380e+03" in s
        assert "PDE_CONVECTION" in s
        assert "different subjects" in s

    def test_the_supersession_pointer_is_present_and_not_hollow(
            self, steady_reach):
        """A note that silently drops its own superseded claim is how version 2
        happened, so the claim is NAMED with its numbers. But full archaeology
        in the payload buys non-silence at the cost of length, and a disclosure
        nobody finishes protects nobody — so the derivation lives in results/
        and this is a pointer. A pointer is only honest if it resolves, so the
        target is read here rather than trusted."""
        s = steady_reach["scale_invariance"]
        assert "SUPERSESSION" in s
        assert "0.01302" in s and "281/T" in s
        assert "ESTIMATOR, not the field" in s
        assert "results/steady_state_reach.md sec. 1" in s

        target = os.path.join(_REPO, "results", "steady_state_reach.md")
        assert os.path.exists(target), (
            "the disclosure points a reader at results/steady_state_reach.md "
            "and it does not exist")
        doc = open(target, encoding="utf-8").read()
        # The derivation the pointer promises, with the magnitudes the payload
        # no longer carries.
        for token in ("AUC(x,T) = C_ss(x)*T + J(x,T)", "281/T", "0.01302",
                      "1.112e-05", "1997.0", "10.82"):
            assert token in doc, f"{token!r} missing from {target}"

    def test_the_untouched_disclosures_survive_verbatim(self, steady_reach):
        """Permeability blindness, the descriptor-path limit and the plasma
        bound are unaffected by the steady-state change and must not have been
        weakened while the capacity clause was rewritten."""
        s = steady_reach["scale_invariance"]
        assert "1.41e-14" in s and "pk_ranking_robustness.md" in s
        assert "bbb_penetrance" in s and "sqrt(D_eff/k)" in s
        assert "LIMIT AS MUCH AS A PROPERTY" in s
        assert "descriptor_falsification.md" in s
        assert "an axis not named here was not tested" in s

    def test_the_limits_come_first_and_name_the_absolute_problem(
            self, steady_reach):
        """The note exists because two previous versions were confidently
        wrong, so it opens with what it does NOT buy. A reader who stops
        halfway must have read the limits, not three exactness claims."""
        s = steady_reach["scale_invariance"]
        assert s.startswith("WHAT THIS MAP DOES NOT BUY")
        # The limits precede the first exactness claim, not merely exist.
        assert s.index("WHAT THIS MAP DOES NOT BUY") < s.index("PLASMA/DOSE")
        # Invariance is not accuracy.
        assert "exactly as wrong" in s
        # The absolute claim is bad in a direction that FLATTERS the compound,
        # and the magnitude is pinned here rather than gestured at.
        assert "14x too shallow" in s and "3.7x too compressed" in s
        assert "kyle2004_penetration.md" in s
        assert "may not" in s
        # And robustness is not validation.
        assert "robustness is not accuracy" in s

    def test_the_finite_window_note_is_corrected_too(self, reach):
        """The old note is still shipped in finite-window mode, so its false
        sentence had to be fixed there as well — not only in the new one."""
        s = reach["scale_invariance"]
        assert "THIS WINDOW IS NOT CONVERGED" not in s
        assert "ESTIMATOR" in s
        assert "manufactures one" in s
        # The falloff residual REPRODUCED on re-measurement (0.37723 against a
        # shipped 0.3772, results/reach_invariances_measure.md), so it is still
        # a live claim and is pinned as one.
        assert "0.3772" in s

    def test_the_refuted_ladder_numbers_appear_ONLY_AS_RETRACTED(self, reach):
        """0.0011 and 0.80 were REFUTED on 2026-08-22 and must never again read
        as claims.

        This test exists because the obvious version of it is vacuous. An
        earlier revision here asserted `"0.0011" in s` under the comment "the
        ladder, over R = 1..39". That assertion still PASSES against the
        corrected note — because the string survives inside the sentence that
        RETRACTS it — so it would have gone on being green while asserting the
        opposite of the truth, and would equally have stayed green if someone
        later restored 0.0011 as a live bound. A substring check cannot tell a
        claim from its retraction, so the retraction itself is what gets pinned.
        """
        s = reach["scale_invariance"]
        # The numbers are still present: a refuted claim is retracted, never
        # deleted, or the next reader re-derives it from the old document.
        assert "0.0011" in s and "0.80" in s
        # But ONLY inside the retraction, and the retraction has to be findable
        # from the numbers themselves rather than sitting elsewhere in a long
        # disclosure.
        i = s.index("0.0011")
        j = s.index("0.80")
        span = s[min(i, j) - 200:max(i, j) + 400]
        assert "a previous version of this note said" in span
        assert "no arm and no reading reproduces any of those three" in span
        assert "reach_invariances_measure.md" in span
        # And the replacement measurement, which is what a reader should carry
        # away instead: the HEADLINE rung is bit-identical while the WHOLE
        # ladder moves 0.4021, i.e. the retracted single number understated the
        # ladder by 365x by reporting the best rung as if it were all of them.
        assert "0.4021" in s
        assert "365x" in s


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
        # The falloff SHAPE residual at R = 39. Still a live claim: re-measured
        # 2026-08-22 at 0.37723 against the shipped 0.3772
        # (results/reach_invariances_measure.md).
        assert "0.3772" in s
        # 0.0011 and 0.80 were REFUTED on the same day and are asserted by
        # test_the_refuted_ladder_numbers_appear_ONLY_AS_RETRACTED instead.
        # They are deliberately NOT asserted here as bare substrings: this test
        # would pass either way, because the strings survive inside the
        # retraction, and a check that cannot fail when its subject is inverted
        # is not a check. Their replacement IS asserted, because that is the
        # number a reader should leave with.
        assert "0.4021" in s                    # the WHOLE ladder, R = 1 to 39
        # The superseded convergence claim, kept so it can be found and
        # recognised as superseded. CORRECTED 2026-08-22: this pair used to be
        # asserted under the comment "the window is NOT converged, with the
        # decay law that says so", which is the claim the note now RETRACTS —
        # 281/T was measuring the AUC estimator, and the field at that same
        # window is converged to 1.112e-05. The strings stay; what they mean
        # here does not.
        assert "0.01302" in s and "281/T" in s
        assert "ESTIMATOR" in s

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
            # THE PRODUCTION PATH PRODUCES THE STEADY-STATE OBJECT. Asserted
            # here rather than inferred from the fact that jobs.py passes a
            # drug: `steady = drug is not None` means one argument is the whole
            # cutover, and a change to that call would flip which object the
            # product emits without touching this module at all.
            assert reach["window_mode"] == "steady_state"
            assert reach["window_s"] is None
            assert reach["E_ref_units"] == "uM"
            assert reach["window_note"] == STEADY_STATE_NOTE
            assert reach["scale_invariance"] == STEADY_SCALE_INVARIANCE_NOTE
            run_id = job.run_id
        finally:
            manager.shutdown()

        # A brand-new manager over the same directory = a server restart.
        restarted = JobManager(root=root, start_worker=False)
        try:
            back = restarted.get(run_id)
            assert back is not None
            assert back.result_view()["reach"] == reach
            # WHICH OBJECT IT IS reaches the STORED record, not only the live
            # response. A run's window_mode is a fact about that run in the same
            # category as its claim-set stamp, and a record that does not carry
            # it is a record that was true about a configuration nobody wrote
            # down — which is the defect this field exists to close.
            restored = back.result_view()["reach"]
            assert restored["window_mode"] == "steady_state"
            assert restored["window_note"] == STEADY_STATE_NOTE
        finally:
            restarted.shutdown()

    def test_a_record_written_before_window_mode_existed_is_not_assumed(self):
        """Runs stored before 2026-08-22 have no ``window_mode``. They PREDATE
        the field and are silent about it, which is not the same as saying
        ``finite_window`` — so a reader must distinguish ABSENT from either
        value rather than defaulting a missing key to the historical
        behaviour. This pins the distinction as a contract; ``WINDOW_MODES``
        deliberately does not contain a value meaning 'unrecorded'."""
        legacy = {"null_reason": None, "headline_x": HEADLINE_X}
        assert "window_mode" not in legacy
        assert legacy.get("window_mode") is None
        assert None not in WINDOW_MODES
        assert "unrecorded" not in WINDOW_MODES

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
