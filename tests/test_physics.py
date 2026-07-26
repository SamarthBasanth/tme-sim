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


def test_label_hypoxia_marks_core_cells(radial_micro_factory):
    # Dense tumor core so oxygen genuinely drops below the hypoxia threshold.
    micro, _, _ = radial_micro_factory(ny=60, nx=60, n_tumor=12000, seed=0)
    cfg = OxygenConfig()
    oxygen = solve_oxygen(micro, cfg)
    label_hypoxia(micro, oxygen, cfg)
    n_hypoxic = int((micro.cells.state == int(CellState.HYPOXIC)).sum())
    assert n_hypoxic > 0, "no cells were labeled hypoxic despite a hypoxic core"


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
