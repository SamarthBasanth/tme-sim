"""PROPERTY tests for the physical invariants of the transport solve.

Every test in this file is of the form *"this must hold for ANY input in this
class"*. None of them asserts that a quantity equals a stored number for one
run — ``tests/test_physics.py`` already does that, and a fixture cannot tell you
whether the scheme is conservative or only whether it still is at that one point.

The five invariants, the derivations that fix their thresholds, and the full
measured ladders live in ``experiments/conservation_audit.py`` and
``results/conservation_audit.md``. This file imports that driver's primitives so
the test and the record cannot drift apart: ``march_mirror`` (a step-by-step
mirror of ``diffusion.simulate_drug_transport`` with the zero-clamp made optional
and per-step accounting added), ``divergent_cross_velocity`` (the minimal
reproducing input for the one violation found), and the frozen ``PREREG``
thresholds.

TWO OF THESE TESTS ARE RED ON PURPOSE, AND SAY SO.
``test_cfl_guard_is_not_sufficient_for_positivity`` and
``test_apply_drug_kill_default_rng_is_unseeded`` assert that a measured DEFECT
is still present. They are not aspirational: they pin the size and the shape of
something the engine does today, so that if it is ever fixed the test fails
loudly and the record gets updated rather than quietly going stale. Each one
names the fix that would make it fail.

Runtime: the whole file is a few seconds on controlled geometry plus one coarse
real section from the shared ``real_micro`` fixture.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO = str(Path(__file__).resolve().parents[1])
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from experiments.conservation_audit import (  # noqa: E402
    ALPHA_PINNED,
    PREREG,
    assert_mirror_is_production,
    cross_undershoot_law,
    divergent_cross_velocity,
    flat_micro,
    march_mirror,
    shipped_cfl,
    strict_positivity_cfl,
)
from tmesim.core import DrugParams, PDEConfig, PKConfig, SimConfig  # noqa: E402
from tmesim.physics.diffusion import (  # noqa: E402
    _effective_diffusivity,
    simulate_drug_transport,
)


# --------------------------------------------------------------------------- #
# The constant this repo pins everywhere, asserted before and after.
# --------------------------------------------------------------------------- #
def test_calibrate_alpha_is_pinned():
    from tmesim.validation.calibration import calibrate_alpha

    assert calibrate_alpha() == ALPHA_PINNED


def test_zz_calibrate_alpha_is_still_pinned():
    """Same assertion, named to sort last: nothing in this file moved it."""
    from tmesim.validation.calibration import calibrate_alpha

    assert calibrate_alpha() == ALPHA_PINNED


# --------------------------------------------------------------------------- #
# The mirror is the load-bearing instrument. If it is not the solver, every
# measurement made through it is about something else.
# --------------------------------------------------------------------------- #
def test_mirror_reproduces_the_solver_bit_for_bit():
    """Static source / time-resolved profile x convection off / on: all four."""
    assert_mirror_is_production(verbose=False)


def test_mirror_without_the_clamp_differs_only_where_the_clamp_fires():
    """The clamp is the ONLY difference between the two mirror modes.

    On an input where the scheme never undershoots the two must be bit-identical;
    on the divergent-cross input they must differ. This is what makes
    ``clamp=False`` a measurement of the scheme rather than of a second scheme.
    """
    micro = flat_micro(24, 24, 15.0, ecm=0.3)
    drug = DrugParams(name="probe", diffusivity=250.0, clearance=2e-4, binding=4e-4)
    src = micro.domain.zeros()
    src[::6, ::6] = 1.0
    pde = PDEConfig(dt=8.0, total_time=480.0, convection=False, save_every=10 ** 9)
    a = march_mirror(micro, drug, src, pde, clamp=True)
    b = march_mirror(micro, drug, src, pde, clamp=False)
    assert np.array_equal(a["c"], b["c"])
    assert np.array_equal(a["auc"], b["auc"])
    assert a["clamp_mass"] == 0.0


# --------------------------------------------------------------------------- #
# (1) MASS BALANCE -- as a budget derived from the discretisation.
# --------------------------------------------------------------------------- #
#   sum(c^{n+1}) + dt*sum(k*c^{n+1}) = sum(c^n) + dt*P_n*sum(S)
# because L has zero row sums and the upwind face fluxes are antisymmetric and
# zeroed on every boundary face. Cumulatively, from c^0 = 0:
#   final mass == INJECTED - dt*sum_n sum(k*c^{n+1})
# The loss quadrature is a RIGHT-ENDPOINT sum because the sink is implicit; the
# returned `auc` is trapezoidal, so k*sum(auc) is the wrong rule and misses by
# O(dt). Both halves are tested.
_BUDGET_CASES = [
    # (label, D, clearance, binding, dt, T, convection, time_resolved, ecm)
    ("inert tracer", 300.0, 0.0, 0.0, 10.0, 600.0, False, False, 0.30),
    ("reactive", 300.0, 2e-4, 5e-4, 10.0, 600.0, False, False, 0.30),
    ("fine dt", 300.0, 2e-4, 5e-4, 0.5, 120.0, False, False, 0.30),
    ("coarse dt", 300.0, 2e-4, 5e-4, 60.0, 1800.0, False, False, 0.30),
    ("fast sink", 50.0, 1e-2, 1e-2, 10.0, 600.0, False, False, 0.55),
    ("slow diffusion", 3.0, 3e-2, 1e-3, 5.0, 300.0, False, False, 0.80),
    ("convection", 300.0, 2e-4, 5e-4, 10.0, 600.0, True, False, 0.30),
    ("time-resolved", 300.0, 2e-4, 5e-4, 10.0, 600.0, True, True, 0.30),
    ("zero ecm", 300.0, 2e-4, 5e-4, 10.0, 600.0, False, False, 0.0),
    ("dense ecm", 300.0, 2e-4, 5e-4, 10.0, 600.0, False, False, 0.95),
]


def _budget(label, D, cl, bd, dt, T, conv, resolved, ecm, hole=False):
    rng = np.random.default_rng(abs(hash(label)) % (2 ** 32))
    ny = nx = 30
    dx = 15.0
    mask = None
    if hole:
        mask = np.ones((ny, nx), dtype=bool)
        mask[7:12, 8:15] = False
    micro = flat_micro(ny, nx, dx, ecm=ecm, tissue_mask=mask)
    drug = DrugParams(name=label, diffusivity=D, clearance=cl, binding=bd)
    src = micro.domain.zeros()
    src[::6, ::7] = 1.2
    vel = None
    if conv:
        vel = np.zeros((2, ny, nx))
        vel[0] = 0.04 * rng.standard_normal((ny, nx))
        vel[1] = 0.04 * rng.standard_normal((ny, nx))
    n_steps = int(round(T / dt))
    profile = (np.abs(rng.standard_normal(n_steps)) + 0.02) if resolved else None
    pde = PDEConfig(dt=dt, total_time=T, convection=conv, save_every=10 ** 9)
    return march_mirror(micro, drug, src, pde, velocity=vel,
                        source_profile=profile, clamp=True), (cl + bd)


@pytest.mark.parametrize("case", _BUDGET_CASES, ids=[c[0] for c in _BUDGET_CASES])
def test_mass_budget_closes_for_any_input(case):
    """PROPERTY: final mass == injected - removed, to solver tolerance, for ANY
    drug / grid / dt / window / source profile / convection setting."""
    m, _ = _budget(*case)
    residual = abs(m["c"].sum() - (m["injected"] - m["removed"]))
    rel = residual / m["injected"]
    tol = PREREG["A1_mass_budget_backward_euler"]["threshold"]
    assert rel < tol, (
        f"{case[0]}: mass budget leaks {rel:.3e} of the injected drug "
        f"(threshold {tol:.0e}). final={m['c'].sum():.6e} "
        f"injected={m['injected']:.6e} removed={m['removed']:.6e}")


def test_mass_budget_closes_across_an_internal_no_flux_hole():
    """PROPERTY: an interior non-tissue region is a no-flux boundary, not a
    leak. Real sections are full of them."""
    m, _ = _budget("internal hole", 300.0, 2e-4, 5e-4, 10.0, 600.0,
                   True, True, 0.30, hole=True)
    rel = abs(m["c"].sum() - (m["injected"] - m["removed"])) / m["injected"]
    assert rel < PREREG["A1_mass_budget_backward_euler"]["threshold"]


def test_mass_budget_absorbs_the_transvascular_exchange_term():
    """PROPERTY: ``exchange_k`` is a first-order sink like any other, so the same
    budget closes with it on. (It ships OFF; this says the budget is not
    special-cased to the shipped default.)"""
    ny = nx = 30
    micro = flat_micro(ny, nx, 15.0, ecm=0.3)
    drug = DrugParams(name="ex", diffusivity=300.0, clearance=2e-4, binding=5e-4)
    src = micro.domain.zeros()
    src[::6, ::7] = 1.2
    ex = np.full((ny, nx), 4e-3)
    pde = PDEConfig(dt=10.0, total_time=600.0, convection=False, save_every=10 ** 9)
    m = march_mirror(micro, drug, src, pde, exchange_k=ex, clamp=True)
    rel = abs(m["c"].sum() - (m["injected"] - m["removed"])) / m["injected"]
    assert rel < PREREG["A1_mass_budget_backward_euler"]["threshold"]


def test_the_trapezoidal_budget_does_not_close():
    """PROPERTY (the negative control): closing the budget with ``k * sum(auc)``
    instead of the right-endpoint sum is the WRONG quadrature and must miss by
    O(dt). If this ever passes, the derivation behind every other budget test in
    this file is wrong and they are all measuring the wrong thing.

    The O(dt) is checked directly: halving dt must roughly halve the error.
    """
    errs = {}
    for dt in (2.5, 5.0, 10.0, 20.0):
        m, k = _budget("trap", 300.0, 2e-4, 5e-4, dt, 600.0, False, False, 0.30)
        wrong = abs(m["c"].sum() - (m["injected"] - k * m["auc"].sum()))
        errs[dt] = wrong / m["injected"]
    floor = PREREG["A2_trapezoid_budget_does_not_close"]["threshold"]
    assert min(errs.values()) > floor, (
        f"the trapezoidal budget closed to {min(errs.values()):.3e}, which "
        f"falsifies the right-endpoint derivation: {errs}")
    # first order in dt: the ratio across a 2x step change sits near 2
    for coarse, fine in ((5.0, 2.5), (10.0, 5.0), (20.0, 10.0)):
        ratio = errs[coarse] / errs[fine]
        assert 1.7 < ratio < 2.3, f"not O(dt): errs={errs}, ratio={ratio:.3f}"


def test_the_zero_clamp_is_not_reached_on_any_physical_input(real_micro, real_drug):
    """PROPERTY: the solver's ``np.maximum(c, 0)`` creates mass whenever it
    fires. On every physical input — ladder geometry and a real section — it
    must fire exactly zero times, so it is never a silent source."""
    for case in _BUDGET_CASES:
        m, _ = _budget(*case)
        assert m["clamp_mass"] == 0.0, f"{case[0]}: clamp created mass"
        assert m["min_raw"] == 0.0, f"{case[0]}: scheme undershot"

    from tmesim.physics.pressure import solve_ifp
    from tmesim.treatment.drugs import delivery_source

    micro = dataclasses.replace(real_micro)
    ifp, velocity = solve_ifp(micro)
    src = delivery_source(micro, real_drug, 1.0)
    pde = PDEConfig(dt=10.0, total_time=600.0, convection=True, save_every=10 ** 9)
    m = march_mirror(micro, real_drug, src, pde, velocity=velocity, clamp=False)
    assert m["min_raw"] == 0.0
    rel = abs(m["c"].sum() - (m["injected"] - m["removed"])) / m["injected"]
    assert rel < PREREG["A1_mass_budget_backward_euler"]["threshold"]


# --------------------------------------------------------------------------- #
# (2) NON-NEGATIVITY -- across a ladder of Peclet numbers.
# --------------------------------------------------------------------------- #
# A = I - dt*L + dt*K is an M-matrix (positive diagonal, non-positive
# off-diagonals, strictly diagonally dominant), so A^-1 >= 0 and the implicit
# half cannot undershoot at ANY dt or Peclet number. Only the explicit upwind
# advection can push the right-hand side negative.
_PECLET_LADDER = (0.0, 0.1, 1.0, 10.0, 100.0, 1000.0)


@pytest.mark.parametrize("peclet", _PECLET_LADDER)
def test_no_undershoot_at_any_peclet_number_with_smooth_flow(peclet):
    """PROPERTY: the RAW (pre-clamp) field never goes negative under a smooth
    flow field, at any Peclet number, with CFL held at 0.9.

    This is measured on the unclamped mirror: the shipped result is clamped and
    therefore cannot answer the question it is being asked."""
    ny = nx = 32
    dx = 15.0
    micro = flat_micro(ny, nx, dx, ecm=0.35)
    drug = DrugParams(name="ladder", diffusivity=300.0, clearance=2e-4, binding=5e-4)
    d_eff = float(_effective_diffusivity(drug, micro.fields.ecm_density).max())
    src = micro.domain.zeros()
    src[ny // 2, nx // 2] = 2.0
    v = peclet * d_eff / dx
    vel = None
    dt = 1.0
    if v > 0:
        vel = np.zeros((2, ny, nx))
        vel[0] = v
        dt = 0.9 * dx / v
    pde = PDEConfig(dt=dt, total_time=120 * dt, convection=v > 0,
                    save_every=10 ** 9)
    m = march_mirror(micro, drug, src, pde, velocity=vel, clamp=False)
    rel = m["min_raw"] / m["max_c"] if m["max_c"] > 0 else 0.0
    assert rel >= PREREG["B1_nonneg_smooth_flow"]["threshold"], (
        f"Pe={peclet}: relative undershoot {rel:.3e} "
        f"(min_raw={m['min_raw']:.3e}, max_c={m['max_c']:.3e})")


@pytest.mark.parametrize("peclet", (1.0, 10.0, 100.0, 1000.0))
def test_no_undershoot_under_smooth_diverging_flow(peclet):
    """PROPERTY: same, for a smoothly DIVERGING radial field at CFL 0.99 — the
    configuration in which a conservative upwind scheme is most exposed."""
    ny = nx = 32
    dx = 15.0
    micro = flat_micro(ny, nx, dx, ecm=0.35)
    drug = DrugParams(name="radial", diffusivity=300.0, clearance=2e-4, binding=5e-4)
    d_eff = float(_effective_diffusivity(drug, micro.fields.ecm_density).max())
    X, Y = micro.domain.coords_um()
    cx, cy = X.mean(), Y.mean()
    radius = max(float(np.abs(X - cx).max()), float(np.abs(Y - cy).max()))
    vmax = peclet * d_eff / dx
    vel = np.zeros((2, ny, nx))
    vel[0] = vmax * (X - cx) / radius
    vel[1] = vmax * (Y - cy) / radius
    speed = float(np.max(np.abs(vel[0]) + np.abs(vel[1])))
    dt = 0.99 * dx / speed
    src = micro.domain.zeros()
    src[ny // 2, nx // 2] = 2.0
    pde = PDEConfig(dt=dt, total_time=120 * dt, convection=True, save_every=10 ** 9)
    m = march_mirror(micro, drug, src, pde, velocity=vel, clamp=False)
    rel = m["min_raw"] / m["max_c"] if m["max_c"] > 0 else 0.0
    assert rel >= PREREG["B1_nonneg_smooth_flow"]["threshold"], (
        f"Pe={peclet}: relative undershoot {rel:.3e}")


# ---- RED ON PURPOSE ------------------------------------------------------- #
@pytest.mark.parametrize("cfl", (0.6, 0.75, 0.9, 1.0))
def test_cfl_guard_is_not_sufficient_for_positivity(cfl):
    """DEFECT PIN (currently red-by-assertion, i.e. it asserts the bug is here).

    ``diffusion.simulate_drug_transport`` accepts any ``dt`` with
    ``dt*max(|vx| + |vy|)/dx <= 1``. That is a max over CELLS of the POINTWISE
    SUM. The upwind stencil's positivity condition is
    ``dt*(max|vx| + max|vy|)/dx <= 1/2`` with the two maxima taken
    INDEPENDENTLY, because cell P can lose through all four of its faces at once
    while no single cell carries both maxima. The guard is therefore up to 2x
    too loose, and on the minimal field that realises the factor exactly the
    scheme undershoots by ``-(2*CFL - 1)`` of the peak — **-100% at the guard's
    own accepted limit** — after which the zero-clamp silently puts that mass
    back.

    THE FIX THAT WOULD BREAK THIS TEST: tighten the guard in
    ``simulate_drug_transport`` to reject ``dt*(max|vx| + max|vy|)/dx > 1/2``
    (or switch to a positivity-preserving advection step). When that lands,
    ``simulate_drug_transport`` will raise here, this test will fail, and
    ``results/conservation_audit.md`` sec. 2 must be superseded by a new dated
    record rather than edited.
    """
    ny = nx = 9
    dx = 15.0
    v = 30.0
    i = j = 4
    micro = flat_micro(ny, nx, dx, ecm=0.0)
    drug = DrugParams(name="tracer", diffusivity=1e-6, clearance=0.0, binding=0.0)
    vel = divergent_cross_velocity(ny, nx, v, i, j)
    src = micro.domain.zeros()
    src[i, j] = 1.0
    speed = float(np.max(np.abs(vel[0]) + np.abs(vel[1])))
    dt = cfl * dx / speed
    n_steps = 8
    pde = PDEConfig(dt=dt, total_time=n_steps * dt, convection=True,
                    save_every=10 ** 9)
    profile = np.zeros(n_steps)
    profile[0] = 1.0

    # the guard accepts it ...
    assert shipped_cfl(vel, dt, dx) <= 1.0 + 1e-12
    simulate_drug_transport(micro, drug, src, pde, velocity=vel,
                            source_profile=profile)
    # ... and the strict bound says it should not have
    assert strict_positivity_cfl(vel, dt, dx) > 0.5

    raw = march_mirror(micro, drug, src, pde, velocity=vel,
                       source_profile=profile, clamp=False)
    rel = raw["min_raw"] / raw["max_c"]
    assert rel < PREREG["B2_cfl_guard_sufficient_for_positivity"]["threshold"], (
        f"CFL={cfl}: expected a large undershoot, measured {rel:.3e}. If the "
        f"guard was tightened, supersede results/conservation_audit.md sec. 2.")
    assert rel == pytest.approx(cross_undershoot_law(cfl), abs=1e-5), (
        f"CFL={cfl}: undershoot {rel:.6f} no longer matches the derived law "
        f"{cross_undershoot_law(cfl):.6f}; the mechanism has changed.")

    clamped = march_mirror(micro, drug, src, pde, velocity=vel,
                           source_profile=profile, clamp=True)
    created = clamped["clamp_mass"] / clamped["injected"]
    assert created == pytest.approx(-cross_undershoot_law(cfl), abs=1e-5), (
        f"CFL={cfl}: the clamp created {created:.4f} of the injected mass")


def test_positivity_holds_at_or_below_the_strict_half_cfl_bound():
    """PROPERTY (the other side of the same coin): on the SAME adversarial field,
    the derived bound ``CFL <= 1/2`` is sufficient — no undershoot at all."""
    ny = nx = 9
    dx = 15.0
    micro = flat_micro(ny, nx, dx, ecm=0.0)
    drug = DrugParams(name="tracer", diffusivity=1e-6, clearance=0.0, binding=0.0)
    vel = divergent_cross_velocity(ny, nx, 30.0, 4, 4)
    src = micro.domain.zeros()
    src[4, 4] = 1.0
    speed = float(np.max(np.abs(vel[0]) + np.abs(vel[1])))
    for cfl in (0.1, 0.25, 0.4, 0.5):
        dt = cfl * dx / speed
        n_steps = 8
        pde = PDEConfig(dt=dt, total_time=n_steps * dt, convection=True,
                        save_every=10 ** 9)
        profile = np.zeros(n_steps)
        profile[0] = 1.0
        raw = march_mirror(micro, drug, src, pde, velocity=vel,
                           source_profile=profile, clamp=False)
        assert raw["min_raw"] == 0.0, f"CFL={cfl} undershot: {raw['min_raw']:.3e}"


def test_production_velocity_stays_inside_the_strict_bound_at_the_default_dt(real_micro):
    """PROPERTY: the shipped ``PDEConfig()`` default ``dt = 0.5 s`` keeps a real
    section's Darcy velocity far inside the derived ``CFL <= 1/2`` bound.

    This is the check that says whether the B2 defect is REACHED in production.
    The app runs ``PDE_DT_S = 10.0`` instead; that arm is measured per grid in
    ``results/conservation_audit.md`` sec. 2.3 and is NOT asserted here, because
    it depends on ``dx`` and on the section."""
    from tmesim.physics.pressure import solve_ifp

    _, velocity = solve_ifp(real_micro)
    dx = float(real_micro.domain.dx)
    assert strict_positivity_cfl(velocity, PDEConfig().dt, dx) < 0.5


# --------------------------------------------------------------------------- #
# (3) MONOTONICITY -- in the obvious parameters, for every compound.
# --------------------------------------------------------------------------- #
_COMPOUNDS = ("temozolomide", "doxorubicin", "lomustine", "methotrexate",
              "vincristine", "cisplatin")


def _real_setup(real_micro):
    from tmesim.physics.pressure import solve_ifp

    _, velocity = solve_ifp(real_micro)
    return velocity


@pytest.mark.parametrize("name", _COMPOUNDS)
def test_more_dose_never_less_exposure_anywhere(real_micro, name):
    """PROPERTY: for every compound, raising the dose raises the exposure at
    EVERY grid point, and does so by exactly the dose ratio (the PDE is linear
    in the source, so this is an equality, not an inequality)."""
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.treatment.drugs import delivery_source
    from tmesim.validation.calibration import calibrated_drugparams

    velocity = _real_setup(real_micro)
    drug = apply_barrier_to_drugparams(calibrated_drugparams(name),
                                       library.get_molecule(name))
    src = delivery_source(real_micro, drug, 1.0)
    if float(np.abs(src).max()) == 0.0:
        pytest.skip(f"{name}: bbb_penetrance is 0, the source is identically zero")
    pde = PDEConfig(dt=10.0, total_time=600.0, convection=True, save_every=10 ** 9)
    base = simulate_drug_transport(real_micro, drug, src, pde,
                                   velocity=velocity).auc
    prev = None
    tol = PREREG["C1_dose_monotone_and_linear"]["threshold"]
    for lam in (1e-3, 0.1, 1.0, 10.0, 1e3):
        a = simulate_drug_transport(real_micro, drug, lam * src, pde,
                                    velocity=velocity).auc
        scaled = lam * base
        mask = np.abs(scaled) > 0
        dev = float((np.abs(a - scaled)[mask] / np.abs(scaled)[mask]).max())
        assert dev < tol, f"{name} lambda={lam}: {dev:.3e} off an exact rescale"
        if prev is not None:
            assert np.all(a >= prev), (
                f"{name}: raising the dose to {lam} LOWERED exposure at "
                f"{int((a < prev).sum())} grid points")
        prev = a


@pytest.mark.parametrize("name", _COMPOUNDS)
def test_longer_window_never_less_auc_anywhere(real_micro, name):
    """PROPERTY: AUC is a cumulative integral of a non-negative field whose
    trajectory does not know where the window ends, so extending the window can
    only add. Checked on the tissue as well as grid-wide, because grid-wide the
    inequality is satisfied trivially by the non-tissue zeros."""
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.treatment.drugs import delivery_source
    from tmesim.validation.calibration import calibrated_drugparams

    velocity = _real_setup(real_micro)
    drug = apply_barrier_to_drugparams(calibrated_drugparams(name),
                                       library.get_molecule(name))
    src = delivery_source(real_micro, drug, 1.0)
    if float(np.abs(src).max()) == 0.0:
        pytest.skip(f"{name}: bbb_penetrance is 0")
    tissue = np.asarray(real_micro.domain.tissue_mask, dtype=bool)
    prev = None
    for T in (300.0, 600.0, 1200.0):
        pde = PDEConfig(dt=10.0, total_time=T, convection=True, save_every=10 ** 9)
        a = simulate_drug_transport(real_micro, drug, src, pde,
                                    velocity=velocity).auc
        if prev is not None:
            assert np.all(a >= prev), f"{name}: AUC fell when the window grew"
            assert float((a - prev)[tissue].min()) >= 0.0
        prev = a


@pytest.mark.parametrize("name", _COMPOUNDS)
def test_total_exposure_is_exactly_invariant_to_diffusivity(real_micro, name):
    """PROPERTY, and the reason the naive diffusivity claim is false.

    At the fixed point ``1^T (K - L + Adv) c* = 1^T S`` reduces to
    ``sum(k*c*) = sum(S)`` (zero row sums again), so with uniform ``k`` the TOTAL
    exposure is ``sum(S)/k`` — independent of the diffusivity. Raising D
    therefore CANNOT raise exposure everywhere: it redistributes. The property
    that survives is measured per distance band in
    ``results/conservation_audit.md`` sec. 3.3.
    """
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.treatment.drugs import delivery_source
    from tmesim.validation.calibration import calibrated_drugparams

    velocity = _real_setup(real_micro)
    drug = apply_barrier_to_drugparams(calibrated_drugparams(name),
                                       library.get_molecule(name))
    src = delivery_source(real_micro, drug, 1.0)
    if float(np.abs(src).max()) == 0.0:
        pytest.skip(f"{name}: bbb_penetrance is 0")
    tissue = np.asarray(real_micro.domain.tissue_mask, dtype=bool)
    # A long window so the field is near its fixed point and the identity bites.
    pde = PDEConfig(dt=20.0, total_time=20000.0, convection=True,
                    save_every=10 ** 9)
    totals = []
    fields = []
    for factor in (0.25, 1.0, 4.0):
        dr = dataclasses.replace(drug, diffusivity=drug.diffusivity * factor)
        a = simulate_drug_transport(real_micro, dr, src, pde,
                                    velocity=velocity).auc
        fields.append(a)
        totals.append(float(a[tissue].sum()))
    spread = (max(totals) - min(totals)) / max(totals)
    assert spread < 1e-6, (
        f"{name}: total exposure moved by {spread:.3e} across a 16x change in "
        f"diffusivity; the conservation identity sum(k*c*) = sum(S) is broken")
    # and therefore: the naive "higher D never gives less anywhere" is FALSE
    assert np.any(fields[-1][tissue] < fields[0][tissue]), (
        f"{name}: higher D raised exposure everywhere while leaving the total "
        f"fixed, which is arithmetically impossible — one of the two "
        f"measurements is wrong")


def test_the_far_field_repair_does_not_transfer_to_real_tissue(real_micro):
    """MEASURED LIMIT, pinned so it is not quietly claimed away.

    In single-vessel geometry "higher D never gives less beyond r*" is a real
    property (test above). On a real section it is NOT: even the outermost
    decile of distance-to-nearest-vessel still contains cells that LOSE exposure
    when D rises, because with many vessels at many spacings the distance to the
    NEAREST vessel does not order cells by their position relative to the
    aggregate crossover. Measured at two grid spacings in
    ``results/conservation_audit.md`` sec. 3.3 and stable across them.

    If this test starts failing, the far-field repair has become claimable on
    real tissue — which would be good news, and must be recorded as a new dated
    measurement rather than assumed.
    """
    from app.server.reach import steady_state_concentration
    from experiments.conservation_audit import vessel_distance_field
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.physics.pressure import solve_ifp
    from tmesim.treatment.drugs import delivery_source
    from tmesim.validation.calibration import calibrated_drugparams

    _, velocity = solve_ifp(real_micro)
    drug = apply_barrier_to_drugparams(calibrated_drugparams("temozolomide"),
                                       library.get_molecule("temozolomide"))
    src = delivery_source(real_micro, drug, 1.0)
    tissue = np.asarray(real_micro.domain.tissue_mask, dtype=bool)
    dist = vessel_distance_field(real_micro, src)
    lo = steady_state_concentration(real_micro, drug, src, velocity=velocity)
    hi = steady_state_concentration(
        real_micro, dataclasses.replace(drug, diffusivity=4 * drug.diffusivity),
        src, velocity=velocity)
    down = (hi < lo) & tissue
    far = tissue & (dist >= np.quantile(dist[tissue], 0.9))
    frac = float((down & far).sum()) / float(far.sum())
    assert frac > 0.0, (
        "the outermost distance decile no longer loses exposure when D rises; "
        "supersede results/conservation_audit.md sec. 3.3 with a new record")


@pytest.mark.parametrize("k", (7e-4, 5e-3))
def test_single_vessel_diffusivity_crossover_is_a_single_clean_radius(k):
    """PROPERTY: in single-vessel geometry the D-crossover exists, is SINGLE, and
    the region beyond it is most of the domain — so "higher D never gives less
    beyond r*" is a real statement there, not one that is true by how r* was
    defined. Measured at the steady state, which is the exact fixed point of the
    same stepper (see ``test_the_march_converges_to_the_steady_solve``), so no
    finite window is entangled in the answer."""
    from app.server.reach import steady_state_concentration

    ny = nx = 81
    dx = 10.0
    micro = flat_micro(ny, nx, dx, ecm=0.0)
    X, Y = micro.domain.coords_um()
    ctr = ny // 2
    r = np.hypot(X - X[ctr, ctr], Y - Y[ctr, ctr])
    src = micro.domain.zeros()
    src[ctr, ctr] = 1.0
    fields = {}
    for D in (50.0, 100.0, 200.0):
        drug = DrugParams(name="one-vessel", diffusivity=D,
                          clearance=0.2 * k, binding=0.8 * k)
        fields[D] = steady_state_concentration(micro, drug, src)
        # the identity that forces a crossing to exist: sum(k*C_ss) == sum(S)
        assert float(fields[D].sum()) * k == pytest.approx(1.0, rel=1e-9)
    ds = sorted(fields)
    for i in range(len(ds) - 1):
        lo, hi = fields[ds[i]], fields[ds[i + 1]]
        ray = hi[ctr, ctr:] - lo[ctr, ctr:]
        sign = np.sign(ray[np.abs(ray) > 0])
        crossings = int((np.diff(sign) != 0).sum())
        assert crossings == 1, (
            f"k={k} D {ds[i]}->{ds[i + 1]}: {crossings} sign changes along the "
            f"ray; the crossover is not a single radius")
        down = hi < lo
        r_star = float(r[down].max())
        beyond = (r > r_star)
        assert beyond.mean() > 0.5, (
            f"k={k}: only {beyond.mean():.4f} of the domain lies beyond r*, so "
            f"'monotone beyond r*' would be a near-vacuous claim")
        assert np.all(hi[beyond] >= lo[beyond])


@pytest.mark.parametrize("name", ("temozolomide", "doxorubicin", "vincristine"))
def test_kill_fraction_is_monotone_in_dose_at_a_fixed_seed(real_micro, name):
    """PROPERTY: with the seed fixed, the Bernoulli draws are the same vector, and
    the Hill map is monotone in exposure, so the killed SET can only grow."""
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.calibration import calibrated_drugparams

    drug = apply_barrier_to_drugparams(calibrated_drugparams(name),
                                       library.get_molecule(name))
    prev = -1.0
    for peak in (1.0, 25.0, 250.0):
        cfg = SimConfig(drug=drug, pk=PKConfig(plasma_peak=peak),
                        pde=PDEConfig(dt=10.0, total_time=600.0,
                                      convection=True, save_every=60),
                        seed=3)
        kill = float(simulate_response(real_micro, cfg).kill_fraction)
        assert kill >= prev, f"{name}: kill fell from {prev} to {kill} at {peak} uM"
        prev = kill


# --------------------------------------------------------------------------- #
# (4) DETERMINISM.
# --------------------------------------------------------------------------- #
def test_same_spec_same_seed_is_bit_identical(real_micro):
    """PROPERTY: two calls with the same spec produce identical bytes."""
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.calibration import calibrated_drugparams

    drug = apply_barrier_to_drugparams(calibrated_drugparams("temozolomide"),
                                       library.get_molecule("temozolomide"))
    cfg = SimConfig(drug=drug,
                    pde=PDEConfig(dt=10.0, total_time=600.0, convection=True,
                                  save_every=60), seed=7)
    a = simulate_response(real_micro, cfg)
    b = simulate_response(real_micro, cfg)
    assert np.array_equal(a.transport.auc, b.transport.auc)
    assert np.array_equal(a.transport.final_concentration,
                          b.transport.final_concentration)
    assert np.array_equal(a.viable_history, b.viable_history)
    assert a.kill_fraction == b.kill_fraction


def test_determinism_is_independent_of_dict_and_set_iteration_order(tmp_path):
    """PROPERTY: the result does not depend on ``PYTHONHASHSEED``.

    Run in FRESH subprocesses, because ``PYTHONHASHSEED`` is read once at
    interpreter start — setting it inside the test process would prove nothing.
    Two seeds here; the full four-seed sweep including the reconstruction,
    oxygen, ECM and cell arrays is in the audit driver."""
    import json
    import os
    import subprocess

    worker = tmp_path / "worker.py"
    worker.write_text(
        "import hashlib, json, sys\n"
        f"sys.path.insert(0, {_REPO!r})\n"
        "import numpy as np\n"
        "from tmesim.reconstruct.real_pipeline import reconstruct_real\n"
        "from tmesim.sim.engine import simulate_response\n"
        "from tmesim.core import SimConfig, PDEConfig\n"
        "from tmesim.molecular import library\n"
        "from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams\n"
        "from tmesim.validation.calibration import calibrated_drugparams\n"
        "h = lambda a: hashlib.sha256(np.ascontiguousarray(\n"
        "    np.asarray(a, dtype=np.float64)).tobytes()).hexdigest()\n"
        "m = reconstruct_real(source='xenium', dx_um=80.0)\n"
        "d = apply_barrier_to_drugparams(calibrated_drugparams('temozolomide'),\n"
        "                                library.get_molecule('temozolomide'))\n"
        "cfg = SimConfig(drug=d, pde=PDEConfig(dt=10.0, total_time=600.0,\n"
        "                convection=True, save_every=60), seed=7)\n"
        "r = simulate_response(m, cfg)\n"
        "print('RESULT ' + json.dumps({'auc': h(r.transport.auc),\n"
        "    'oxygen': h(m.fields.oxygen), 'cells_x': h(m.cells.x),\n"
        "    'kill': repr(r.kill_fraction)}))\n"
    )
    seen = []
    for seed in ("0", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        out = subprocess.run([sys.executable, str(worker)], env=env,
                             capture_output=True, text=True, timeout=900)
        lines = [ln for ln in out.stdout.splitlines() if ln.startswith("RESULT ")]
        assert lines, (f"worker failed under PYTHONHASHSEED={seed}:\n"
                       + "\n".join(out.stderr.splitlines()[-10:]))
        seen.append(json.loads(lines[-1][len("RESULT "):]))
    assert seen[0] == seen[1], f"hash-seed dependent result: {seen}"


# ---- RED ON PURPOSE ------------------------------------------------------- #
def test_apply_drug_kill_default_rng_is_unseeded(real_micro):
    """DEFECT PIN (asserts the bug is still here).

    ``cells.pkpd.apply_drug_kill`` documents ``rng=None`` as "a fresh
    ``np.random.default_rng()``" — i.e. seeded from OS entropy. A caller that
    omits ``rng`` therefore gets a different kill every time. This test asserts
    that is STILL true, so the day it is fixed the record gets updated instead
    of going stale.

    THE FIX THAT WOULD BREAK THIS TEST: give ``apply_drug_kill`` a default seed,
    or make ``rng`` required. The engine path is safe today and is asserted
    separately below; this is about the direct callers.
    """
    from tmesim.cells.agents import update_cell_states
    from tmesim.cells.pkpd import apply_drug_kill
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.calibration import calibrated_drugparams

    drug = apply_barrier_to_drugparams(calibrated_drugparams("temozolomide"),
                                       library.get_molecule("temozolomide"))
    cfg = SimConfig(drug=drug,
                    pde=PDEConfig(dt=10.0, total_time=600.0, convection=True,
                                  save_every=60), seed=7)
    res = simulate_response(real_micro, cfg)

    kills = []
    for _ in range(6):
        work = dataclasses.replace(real_micro, cells=real_micro.cells.copy())
        update_cell_states(work, res.fields["oxygen"], cfg.oxygen)
        apply_drug_kill(work, res.transport, cfg.pd, drug)      # no rng
        kills.append(1.0 - float(work.cells.viable_tumor_fraction()))
    assert len(set(kills)) > 1, (
        "apply_drug_kill(rng=None) is now deterministic — the defect this test "
        "pins has been fixed. Supersede results/conservation_audit.md sec. 4 "
        "with a new dated record and delete this test.")


def test_the_unseeded_default_is_not_reachable_through_the_engine(real_micro):
    """PROPERTY: ``simulate_response`` always passes a seeded generator, so the
    defect above cannot reach a shipped run."""
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.calibration import calibrated_drugparams

    drug = apply_barrier_to_drugparams(calibrated_drugparams("temozolomide"),
                                       library.get_molecule("temozolomide"))
    cfg = SimConfig(drug=drug,
                    pde=PDEConfig(dt=10.0, total_time=600.0, convection=True,
                                  save_every=60), seed=7)
    kills = {float(simulate_response(real_micro, cfg).kill_fraction)
             for _ in range(3)}
    assert len(kills) == 1, f"simulate_response is non-deterministic: {kills}"


# --------------------------------------------------------------------------- #
# (5) SCALE INVARIANCE -- re-derived, not read off the stored records.
# --------------------------------------------------------------------------- #
def _scale_fixture():
    ny = nx = 33
    dx = 12.0
    micro = flat_micro(ny, nx, dx, ecm=0.25)
    drug = DrugParams(name="scale", diffusivity=120.0, clearance=1e-4, binding=6e-4)
    src = micro.domain.zeros()
    src[::6, ::6] = 0.8
    rng = np.random.default_rng(11)
    vel = np.zeros((2, ny, nx))
    vel[0] = 0.02 * rng.standard_normal((ny, nx))
    vel[1] = 0.02 * rng.standard_normal((ny, nx))
    return micro, drug, src, vel


def _relmax(a, b):
    den = np.abs(b)
    mask = den > 0
    return float((np.abs(a - b)[mask] / den[mask]).max()) if mask.any() else 0.0


@pytest.mark.parametrize("lam", (1e-3, 0.1, 10.0, 1e3))
def test_auc_is_exactly_proportional_to_any_scalar_source_factor(lam):
    """PROPERTY: this is the invariance the reach map rests on. The plasma level
    and ``bbb_penetrance`` both enter as a single multiplication on the source,
    so a ratio between two points of the map cannot depend on either."""
    micro, drug, src, vel = _scale_fixture()
    pde = PDEConfig(dt=4.0, total_time=1200.0, convection=True, save_every=10 ** 9)
    base = simulate_drug_transport(micro, drug, src, pde, velocity=vel).auc
    a = simulate_drug_transport(micro, drug, lam * src, pde, velocity=vel).auc
    assert _relmax(a, lam * base) < PREREG["E1_source_scale_exact"]["threshold"]


def test_superposition():
    """PROPERTY: the transport operator is linear, so two vessels are the sum of
    each vessel alone. A clamp that fired would break this; it does not."""
    micro, drug, _, vel = _scale_fixture()
    s1 = micro.domain.zeros()
    s1[::6, ::6] = 0.8
    s2 = micro.domain.zeros()
    s2[3::8, 2::8] = 1.7
    pde = PDEConfig(dt=4.0, total_time=1200.0, convection=True, save_every=10 ** 9)
    a1 = simulate_drug_transport(micro, drug, s1, pde, velocity=vel).auc
    a2 = simulate_drug_transport(micro, drug, s2, pde, velocity=vel).auc
    a12 = simulate_drug_transport(micro, drug, s1 + s2, pde, velocity=vel).auc
    assert _relmax(a1 + a2, a12) < PREREG["E2_superposition"]["threshold"]


@pytest.mark.parametrize("s", (0.25, 0.5, 2.0, 4.0))
def test_capacity_identity_holds_without_a_capacity_term(s):
    """PROPERTY: a tissue-capacity term ``R`` multiplies the time derivative and
    nothing else, so the capacity-``R`` problem IS the ordinary problem with
    ``(D, k, v, S)`` and ``(T, dt)`` rescaled. Re-derived here rather than taken
    from ``results/reach_invariances_measure.md`` sec. 3:

        auc(D, k, v, S; T, dt) == s * auc(sD, sk, sv, sS; T/s, dt/s)

    with ``s = 1/R``. No capacity term is implemented and none is needed to test
    the identity it would satisfy."""
    micro, drug, src, vel = _scale_fixture()
    T, dt = 1200.0, 4.0
    base = simulate_drug_transport(
        micro, drug, src,
        PDEConfig(dt=dt, total_time=T, convection=True, save_every=10 ** 9),
        velocity=vel).auc
    scaled_drug = dataclasses.replace(drug, diffusivity=drug.diffusivity * s,
                                      clearance=drug.clearance * s,
                                      binding=drug.binding * s)
    scaled = simulate_drug_transport(
        micro, scaled_drug, s * src,
        PDEConfig(dt=dt / s, total_time=T / s, convection=True,
                  save_every=10 ** 9),
        velocity=s * vel).auc
    assert _relmax(s * scaled, base) < PREREG["E3_capacity_time_rescale"]["threshold"]


@pytest.mark.parametrize("s", (1e-3, 0.5, 2.0, 1e3))
def test_steady_state_is_exactly_capacity_invariant(s):
    """PROPERTY: ``(K - L + Adv) c* = S`` is homogeneous of degree one in a
    common factor on ``(D, k, v, S)``, so ``c*`` is EXACTLY invariant to it. That
    is what "capacity is absent from the steady problem" means, stated as
    something that can fail."""
    from app.server.reach import steady_state_concentration

    micro, drug, src, vel = _scale_fixture()
    c0 = steady_state_concentration(micro, drug, src, velocity=vel)
    scaled_drug = dataclasses.replace(drug, diffusivity=drug.diffusivity * s,
                                      clearance=drug.clearance * s,
                                      binding=drug.binding * s)
    c1 = steady_state_concentration(micro, scaled_drug, s * src, velocity=s * vel)
    tol = PREREG["E4_steady_state_capacity_invariant"]["threshold"]
    assert _relmax(c1, c0) < tol


def test_the_march_converges_to_the_steady_solve():
    """PROPERTY: the steady solve is the march's own fixed point, so a long
    enough march must reach it. If it does not, one of the two is not the
    operator it claims to be."""
    from app.server.reach import steady_state_concentration

    micro, drug, src, vel = _scale_fixture()
    c0 = steady_state_concentration(micro, drug, src, velocity=vel)
    res = simulate_drug_transport(
        micro, drug, src,
        PDEConfig(dt=10.0, total_time=86400.0, convection=True,
                  save_every=10 ** 9),
        velocity=vel)
    tol = PREREG["E5_march_reaches_steady_state"]["threshold"]
    assert _relmax(res.final_concentration, c0) < tol
