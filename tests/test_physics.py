"""Physics numerics tests: drug transport, oxygen, and interstitial pressure.

These pin the load-bearing qualitative + quantitative properties the Stage 0
kill-test relies on:

1. transport mass conservation under no-flux with no source/sink,
2. transport non-negativity,
3. oxygen: hypoxic tumor core vs oxygenated periphery,
4. pressure: elevated core IFP vs the margin (the Jain effect).

Everything runs on small grids / short windows to stay fast.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core import DrugParams, OxygenConfig, PDEConfig
from tmesim.physics.diffusion import simulate_drug_transport
from tmesim.physics.oxygen import label_hypoxia, solve_oxygen
from tmesim.physics.pressure import solve_ifp
from tmesim.core import CellState


# --------------------------------------------------------------------------- #
# (1) Mass conservation
# --------------------------------------------------------------------------- #
def test_transport_mass_conservation(inert_drug, flat_micro_factory):
    """With clearance=binding=0 and no-flux boundaries, nothing escapes: the
    final tissue mass must equal exactly the injected mass (n_steps*dt*sum S)."""
    micro = flat_micro_factory(ny=40, nx=40, dx=15.0, ecm_value=0.3)
    pde = PDEConfig(dt=10.0, total_time=1200.0, convection=False, save_every=20)

    source = micro.domain.zeros()
    source[20, 20] = 2.0  # a single interior vessel, uM/s

    res = simulate_drug_transport(micro, inert_drug, source, pde)

    n_steps = int(round(pde.total_time / pde.dt))
    injected = n_steps * pde.dt * float(source.sum())
    final_mass = float(res.final_concentration.sum())
    rel_err = abs(final_mass - injected) / injected

    assert rel_err < 1e-9, f"mass not conserved: rel_err={rel_err:.2e}"


def test_transport_zero_source_stays_empty(inert_drug, flat_micro_factory):
    """No source, zero initial condition -> the field stays identically zero."""
    micro = flat_micro_factory(ny=30, nx=30, dx=15.0, ecm_value=0.5)
    pde = PDEConfig(dt=10.0, total_time=600.0, convection=False, save_every=20)
    res = simulate_drug_transport(micro, inert_drug, micro.domain.zeros(), pde)
    assert float(np.abs(res.final_concentration).max()) == 0.0
    assert float(np.abs(res.auc).max()) == 0.0


# --------------------------------------------------------------------------- #
# (2) Non-negativity
# --------------------------------------------------------------------------- #
def test_transport_nonnegativity(flat_micro_factory):
    """Concentration and accumulated AUC must never go negative, even with a
    reactive drug and multiple scattered sources."""
    micro = flat_micro_factory(ny=48, nx=48, dx=15.0, ecm_value=0.4)
    drug = DrugParams(name="tmz", diffusivity=300.0, clearance=3e-4, binding=1e-4)
    pde = PDEConfig(dt=10.0, total_time=1800.0, convection=False, save_every=30)

    source = micro.domain.zeros()
    source[::12, ::12] = 1.0  # scattered vessels

    res = simulate_drug_transport(micro, drug, source, pde)

    assert res.final_concentration.min() >= -1e-9
    assert res.auc.min() >= -1e-9
    # And the drug actually reached somewhere (non-trivial run).
    assert res.auc.max() > 0.0
    assert np.all(np.isfinite(res.final_concentration))


# --------------------------------------------------------------------------- #
# (3) Oxygen: hypoxic core
# --------------------------------------------------------------------------- #
def test_oxygen_core_more_hypoxic_than_periphery(radial_micro_factory):
    micro, r, core_radius = radial_micro_factory(ny=60, nx=60, n_tumor=5000,
                                                 seed=0)
    cfg = OxygenConfig()
    oxygen = solve_oxygen(micro, cfg)

    assert oxygen.shape == micro.domain.shape
    assert oxygen.min() >= 0.0 and oxygen.max() <= 1.0

    core_mask = r < 0.5 * core_radius
    periph_mask = r > 1.3 * core_radius
    core_o2 = float(oxygen[core_mask].mean())
    periph_o2 = float(oxygen[periph_mask].mean())
    assert core_o2 < periph_o2, (
        f"tumor core O2 {core_o2:.3f} is not below periphery {periph_o2:.3f}")


def test_oxygen_is_deterministic(radial_micro_factory):
    micro, _, _ = radial_micro_factory(ny=48, nx=48, n_tumor=3000, seed=1)
    cfg = OxygenConfig()
    o1 = solve_oxygen(micro, cfg)
    o2 = solve_oxygen(micro, cfg)
    assert np.array_equal(o1, o2)


def test_oxygen_sink_is_grid_independent(radial_micro_factory):
    """The O2 sink must be an areal DENSITY, not a per-grid-cell COUNT.

    Same physical tissue (identical 900x900 um extent, identical 12,000 tumor cells
    at identical positions) solved on four grids spanning a 4x range of dx. A raw
    count sink scales as dx^2 and used to swing the hypoxic fraction from 0.998 at
    dx = 30 um to 0.000 at dx = 10 um; the areal density is invariant.
    """
    cfg = OxygenConfig()
    hypoxic, mean_o2 = [], []
    for dx, n in ((30.0, 30), (15.0, 60), (10.0, 90), (7.5, 120)):
        micro, _, _ = radial_micro_factory(ny=n, nx=n, dx=dx, n_tumor=12000, seed=0)
        assert micro.domain.extent_um == (900.0, 900.0)  # same physical tissue
        oxygen = solve_oxygen(micro, cfg)
        # Sample at the (grid-independent) real cell positions.
        col, row = micro.domain.to_index(micro.cells.x, micro.cells.y)
        local = oxygen[row, col]
        hypoxic.append(float((local < cfg.hypoxia_threshold).mean()))
        mean_o2.append(float(local.mean()))

    assert max(mean_o2) - min(mean_o2) < 0.01, f"mean O2 varies with dx: {mean_o2}"
    assert max(hypoxic) - min(hypoxic) < 0.02, f"hypoxic fraction varies with dx: {hypoxic}"


def test_label_hypoxia_marks_core_cells(radial_micro_factory):
    # Dense tumor core so oxygen genuinely drops below the hypoxia threshold.
    micro, _, _ = radial_micro_factory(ny=60, nx=60, n_tumor=12000, seed=0)
    cfg = OxygenConfig()
    oxygen = solve_oxygen(micro, cfg)
    label_hypoxia(micro, oxygen, cfg)
    n_hypoxic = int((micro.cells.state == int(CellState.HYPOXIC)).sum())
    assert n_hypoxic > 0, "no cells were labeled hypoxic despite a hypoxic core"


# --------------------------------------------------------------------------- #
# (3b) Oxygen: opt-in Michaelis-Menten sink form (default stays "linear")
# --------------------------------------------------------------------------- #
def test_default_sink_form_is_linear_and_unchanged(radial_micro_factory):
    """``OxygenConfig()`` must keep solving the ORIGINAL linear-sink system --
    the new ``sink_form`` option must not perturb the default in any way."""
    micro, _, _ = radial_micro_factory(ny=48, nx=48, n_tumor=6000, seed=2)
    cfg = OxygenConfig()
    assert cfg.sink_form == "linear"
    oxygen_default = solve_oxygen(micro, cfg)
    oxygen_explicit = solve_oxygen(micro, OxygenConfig(sink_form="linear"))
    assert np.array_equal(oxygen_default, oxygen_explicit)


def test_michaelis_menten_sink_is_an_available_opt_in_option(radial_micro_factory):
    """The MM sink form must be selectable, produce a valid, deterministic
    field, and (holding consumption fixed) draw MORE oxygen down than the
    linear sink -- consistent with a non-self-throttling sink term."""
    micro, _, _ = radial_micro_factory(ny=48, nx=48, n_tumor=6000, seed=2)
    cfg_linear = OxygenConfig()
    cfg_mm = OxygenConfig(sink_form="michaelis_menten", mm_half_saturation=0.05)

    o_linear = solve_oxygen(micro, cfg_linear)
    o_mm = solve_oxygen(micro, cfg_mm)
    o_mm_again = solve_oxygen(micro, cfg_mm)

    assert o_mm.shape == micro.domain.shape
    assert o_mm.min() >= 0.0 and o_mm.max() <= 1.0
    assert np.array_equal(o_mm, o_mm_again), "MM sink solve is not deterministic"
    # Not bit-identical to the linear sink (it is a materially different PDE).
    assert not np.array_equal(o_linear, o_mm)
    # Where there is tumor, the non-self-throttling MM sink pulls O2 lower.
    tum_mask = micro.cells.tumor_mask()
    col, row = micro.domain.to_index(micro.cells.x[tum_mask], micro.cells.y[tum_mask])
    assert float(o_mm[row, col].mean()) < float(o_linear[row, col].mean())


def test_sink_form_rejects_unknown_values(radial_micro_factory):
    micro, _, _ = radial_micro_factory(ny=30, nx=30, n_tumor=2000, seed=3)
    cfg = OxygenConfig(sink_form="quadratic")
    with pytest.raises(ValueError):
        solve_oxygen(micro, cfg)


def test_michaelis_menten_raises_rather_than_silently_returning_unconverged(
    radial_micro_factory,
):
    """A production nonlinear solver must fail loudly, not return a field that
    never reached tolerance."""
    micro, _, _ = radial_micro_factory(ny=30, nx=30, n_tumor=2000, seed=4)
    cfg = OxygenConfig(sink_form="michaelis_menten", mm_half_saturation=0.02,
                       mm_max_iter=1, mm_tol=1e-14)
    with pytest.raises(RuntimeError):
        solve_oxygen(micro, cfg)


# --------------------------------------------------------------------------- #
# (4) Pressure: elevated core IFP (Jain effect)
# --------------------------------------------------------------------------- #
def test_pressure_core_ifp_exceeds_margin(radial_micro_factory):
    micro, r, core_radius = radial_micro_factory(ny=60, nx=60, n_tumor=3000,
                                                 seed=0)
    ifp, velocity = solve_ifp(micro)

    assert ifp.shape == micro.domain.shape
    assert velocity.shape == (2,) + micro.domain.shape
    assert ifp.min() >= 0.0  # pressure is non-negative

    core_mask = r < 0.35 * core_radius
    margin_mask = (r > 0.8 * core_radius) & (r < 1.1 * core_radius)
    core_ifp = float(ifp[core_mask].mean())
    margin_ifp = float(ifp[margin_mask].mean())
    assert core_ifp > margin_ifp, (
        f"core IFP {core_ifp:.3f} should exceed margin IFP {margin_ifp:.3f}")


# --------------------------------------------------------------------------- #
# (5) Transvascular exchange (PDEConfig.transvascular_exchange / exchange_k)
# --------------------------------------------------------------------------- #
def test_exchange_k_none_is_bit_identical(inert_drug, flat_micro_factory):
    """``exchange_k=None`` must reproduce the one-way path EXACTLY.

    The whole backward-compatibility argument for adding the term rests on
    this: every stored reference result in the repo was produced on the
    one-way path, so the default must not move a single bit.
    """
    micro = flat_micro_factory()
    ny, nx = micro.domain.ny, micro.domain.nx
    src = np.zeros((ny, nx))
    src[ny // 2, nx // 2] = 1e-3
    pde = PDEConfig(dt=1.0, total_time=60.0, convection=False, save_every=10)

    a = simulate_drug_transport(micro, inert_drug, src, pde)
    b = simulate_drug_transport(micro, inert_drug, src, pde, exchange_k=None)
    assert np.array_equal(a.auc, b.auc)


def test_exchange_k_zero_is_bit_identical(inert_drug, flat_micro_factory):
    """A zero exchange field is also a no-op: ``k + 0.0 == k`` in IEEE-754."""
    micro = flat_micro_factory()
    ny, nx = micro.domain.ny, micro.domain.nx
    src = np.zeros((ny, nx))
    src[ny // 2, nx // 2] = 1e-3
    pde = PDEConfig(dt=1.0, total_time=60.0, convection=False, save_every=10)

    a = simulate_drug_transport(micro, inert_drug, src, pde)
    b = simulate_drug_transport(micro, inert_drug, src, pde,
                                exchange_k=np.zeros((ny, nx)))
    assert np.array_equal(a.auc, b.auc)


def test_exchange_bounds_tissue_below_plasma(flat_micro_factory):
    """THE point of the term: with exchange on, tissue concentration can never
    exceed the plasma level, for any dose, at any time.

    The one-way form has no such bound -- its steady state is
    ``PS*Cp/(clearance+binding)``, which for a slowly-cleared drug sits far
    ABOVE plasma (measured: 5-FU at 7.8x plasma in results/breast_kp_validation.md,
    against 0.61x measured in humans). This test pins the fix, not the symptom:
    it asserts the bound holds by construction, on a drug with deliberately
    tiny loss terms where the one-way form blows past plasma.
    """
    micro = flat_micro_factory()
    ny, nx = micro.domain.ny, micro.domain.nx
    micro.fields.vascular_density = np.full((ny, nx), 0.02)

    # Almost no clearance/binding -> the one-way steady state is enormous.
    drug = DrugParams(name="slow", diffusivity=200.0, molecular_weight=300.0,
                      ic50=1.0, hill=1.0, bbb_penetrance=0.5,
                      binding=1e-7, clearance=1e-7)
    plasma = 10.0
    ps = micro.fields.vascular_density * drug.bbb_penetrance   # 1/s
    src = ps * plasma
    pde = PDEConfig(dt=1.0, total_time=3600.0, convection=False, save_every=100)

    one_way = simulate_drug_transport(micro, drug, src, pde)
    exchange = simulate_drug_transport(micro, drug, src, pde, exchange_k=ps)

    tis = micro.domain.tissue_mask
    c_one = one_way.auc[tis] / pde.total_time
    c_exc = exchange.auc[tis] / pde.total_time

    # The one-way form must actually breach plasma here, or the test is not
    # exercising the failure it claims to.
    assert c_one.max() > plasma, (
        f"one-way path did not exceed plasma ({c_one.max():.4g} <= {plasma}); "
        "test no longer exercises the defect it was written for")
    # The exchange form must not, ever.
    assert c_exc.max() <= plasma + 1e-9, (
        f"exchange path breached plasma: {c_exc.max():.6g} > {plasma}")


def test_exchange_k_rejects_negative(inert_drug, flat_micro_factory):
    """A negative exchange rate would pump drug in against its own gradient."""
    micro = flat_micro_factory()
    ny, nx = micro.domain.ny, micro.domain.nx
    pde = PDEConfig(dt=1.0, total_time=10.0, convection=False, save_every=10)
    bad = np.full((ny, nx), -1e-4)
    with pytest.raises(ValueError, match="non-negative"):
        simulate_drug_transport(micro, inert_drug, np.zeros((ny, nx)), pde,
                                exchange_k=bad)


# --------------------------------------------------------------------------- #
# (6) The transvascular rate scale k_0 (treatment.drugs.PS_SCALE_S_INV)
# --------------------------------------------------------------------------- #
def test_ps_scale_default_is_bit_identical(flat_micro_factory):
    """The default 1.0 must not move a single bit.

    PS_SCALE_S_INV exists to make a dimensional error visible, not to change
    any number. Every stored reference in the repo was produced with the
    implicit factor of 1, so the explicit default has to reproduce it exactly.
    """
    from tmesim.core import DrugParams
    from tmesim.treatment.drugs import PS_SCALE_S_INV, delivery_source

    assert PS_SCALE_S_INV == 1.0, (
        "the default must stay 1.0 -- changing it rescales every absolute "
        "delivery number in results/ and is a re-validation exercise")

    micro = flat_micro_factory()
    ny, nx = micro.domain.ny, micro.domain.nx
    micro.fields.vascular_density = np.full((ny, nx), 0.03)
    drug = DrugParams(name="x", diffusivity=200.0, molecular_weight=300.0,
                      ic50=1.0, hill=1.0, bbb_penetrance=0.4)

    a = delivery_source(micro, drug, 10.0)
    b = delivery_source(micro, drug, 10.0, ps_scale=None)
    c = delivery_source(micro, drug, 10.0, ps_scale=1.0)
    assert np.array_equal(a, b)
    assert np.array_equal(a, c)


def test_ps_scale_is_linear_and_reaches_exchange(flat_micro_factory):
    """The scale multiplies the source, and the SAME scale reaches the
    back-flux field -- otherwise source and exchange would silently model
    different permeabilities."""
    from tmesim.core import DrugParams
    from tmesim.treatment.drugs import delivery_source, exchange_coefficient

    micro = flat_micro_factory()
    ny, nx = micro.domain.ny, micro.domain.nx
    micro.fields.vascular_density = np.full((ny, nx), 0.03)
    drug = DrugParams(name="x", diffusivity=200.0, molecular_weight=300.0,
                      ic50=1.0, hill=1.0, bbb_penetrance=0.4)

    base = delivery_source(micro, drug, 10.0)
    assert np.allclose(delivery_source(micro, drug, 10.0, ps_scale=0.25),
                       base * 0.25, rtol=0, atol=0)
    # exchange_coefficient IS delivery_source at Cp=1, scale included.
    assert np.array_equal(
        exchange_coefficient(micro, drug, ps_scale=0.25),
        delivery_source(micro, drug, 1.0, ps_scale=0.25))

    with pytest.raises(ValueError, match="non-negative"):
        delivery_source(micro, drug, 10.0, ps_scale=-1.0)


def test_ps_scale_calibration_residual_equals_half_spread():
    """The calibration's key identity, which is the whole argument for why a
    better k_0 is NOT the remaining fix.

    A single scale can only remove the systematic bias. What is left is
    exactly the square root of the disagreement between the per-drug implied
    scales -- i.e. the engine's cross-drug permeability ORDERING error, which
    no rescale can touch.
    """
    from tmesim.validation.ps_scale import calibrate_ps_scale

    c = calibrate_ps_scale()
    assert c["n_anchors"] >= 2
    # The fit must improve on the implicit factor of 1.
    assert c["gm_fold_after"] < c["gm_fold_before"]
    # And the residual must be exactly the half-spread, to floating tolerance.
    assert c["gm_fold_after"] == pytest.approx(c["irreducible_residual"],
                                               rel=1e-9)
    assert c["irreducible_residual"] == pytest.approx(
        c["disagreement_fold"] ** 0.5, rel=1e-12)
