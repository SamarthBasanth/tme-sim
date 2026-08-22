"""The entry/spread split: is it OFF, is it exact when off, and does it BITE when on?

WHAT THIS FILE GUARDS, AND WHY EACH CHECK CAN FAIL
==================================================
``results/permeability_split.md`` rests on one exact property and one measured
consequence, and both are asserted here against the live solver rather than
quoted:

* **the one-way law is exactly homogeneous of degree 1 in the entry term** --
  scaling ``PS`` scales the whole field and leaves the SHAPE bit-invariant, so
  the entry term cannot be part of any gradient error;
* **the split breaks that homogeneity** -- which is the only reason the split
  can change an absolute per-cell verdict at all.

The second one is the check that matters, because it is the one with a live
failure mode. ``tmesim.validation.propagation`` documents the exact defect:
a caller that sets ``PDEConfig.transvascular_exchange`` but never delivers
``exchange_k`` to the solver produces a run that is BIT-IDENTICAL to the flag
being off, and reports "this physics change has no effect" with total
confidence. A test that only asserted "OFF and ON differ in the counts" would
be green on a broken build for exactly the wrong reason; this file instead
asserts the SHAPE-invariance property directly, which an inert flag cannot
satisfy.

DEFECT PLANTED AND WATCHED TO BITE
-----------------------------------
Each of the three load-bearing checks below was run against a deliberately
broken build before being committed, and the failures are recorded in
``results/permeability_split.md`` section 6 with their verbatim messages:

1. ``delivery_pair`` made to ignore ``split=True`` and return
   ``(source, None)`` -- the INERT defect -- fails
   ``test_split_on_carries_the_same_ps`` and
   ``test_split_breaks_the_homogeneity``;
2. ``SPLIT_DEFAULT`` flipped to ``True`` -- the shipped-default defect -- fails
   ``test_the_shipped_default_is_off``;
3. ``exchange_coefficient`` given a different ``ps_scale`` from the source --
   the source/back-flux drift ``drugs.exchange_coefficient`` warns about --
   fails ``test_split_on_carries_the_same_ps``.

Configuration is stated in the same place as every number: the PDE checks use a
1-D slab of 25 columns at ``dx = 8 um``, ``ecm_density = 0.3`` uniform,
``dt = 2 s``, ``total_time = 1200 s``, convection OFF, and route through
``sim.engine.simulate_response`` so the ``transvascular_exchange`` flag reaches
the solver by the ONE production wiring point rather than by hand.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from tmesim.core import (
    CellPopulation,
    CellState,
    CellType,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    PKConfig,
    SimConfig,
    TissueFields,
)
from tmesim.sim.engine import simulate_response
from tmesim.treatment.drugs import (
    delivery_source,
    exchange_coefficient,
    plasma_profile,
)
from tmesim.treatment.transvascular import (
    SPLIT_DEFAULT,
    boundary_layer_um,
    delivery_pair,
    entry_elasticity,
    entry_transport_regime,
    implied_permeability_um_s,
    one_way_plasma_multiple,
    ps_from_permeability,
    saturation_fraction,
    surface_to_volume_per_um,
)

# --- the slab, stated once and used by every PDE check --------------------- #
SLAB_NX = 25
SLAB_NY = 3
SLAB_DX_UM = 8.0
SLAB_ECM = 0.3
SLAB_DT_S = 2.0
SLAB_TOTAL_S = 1200.0


def _slab_micro():
    """One-sided planar source at column 0, sealed elsewhere (Neumann)."""
    dom = Domain(nx=SLAB_NX, ny=SLAB_NY, dx=SLAB_DX_UM)
    vasc = np.zeros(dom.shape, dtype=np.float64)
    vasc[:, 0] = 1.0
    fields = TissueFields(
        ecm_density=np.full(dom.shape, SLAB_ECM), vascular_density=vasc)
    cells = CellPopulation(
        x=np.array([0.0]), y=np.array([0.0]),
        cell_type=np.array([CellType.TUMOR]), state=np.array([CellState.VIABLE]),
        sensitivity=np.array([1.0]))
    return Microenvironment(domain=dom, fields=fields, cells=cells,
                            patient_id="split_slab")


def _patch_micro():
    """A small 2-D patch with a graded vessel field, for the field builders."""
    dom = Domain(nx=11, ny=7, dx=15.0)
    vasc = np.linspace(0.0, 0.6, 11)[None, :].repeat(7, axis=0)
    fields = TissueFields(ecm_density=dom.ones() * 0.25, vascular_density=vasc)
    cells = CellPopulation(
        x=np.array([15.0]), y=np.array([15.0]),
        cell_type=np.array([CellType.TUMOR]), state=np.array([CellState.VIABLE]),
        sensitivity=np.array([0.5]))
    return Microenvironment(domain=dom, fields=fields, cells=cells)


DRUG = DrugParams(name="split-probe", bbb_penetrance=0.25, diffusivity=40.0,
                  clearance=4e-4, binding=1.6e-3, ic50=10.0)


def _profile(drug, exchange: bool) -> np.ndarray:
    """Depth profile of the AUC across the slab, averaged along y."""
    micro = _slab_micro()
    pde = PDEConfig(dt=SLAB_DT_S, total_time=SLAB_TOTAL_S, convection=False,
                    save_every=100000, transvascular_exchange=exchange)
    res = simulate_response(micro, SimConfig(drug=drug, pde=pde, seed=0))
    return np.asarray(res.fields["auc"], dtype=np.float64).mean(axis=0)


# =========================================================================== #
# 1. The shipped position: OFF, and exactly OFF.
# =========================================================================== #
def test_the_shipped_default_is_off():
    """Both switches that select the delivery law are False as shipped.

    Fails if anyone flips ``SPLIT_DEFAULT`` or
    ``PDEConfig.transvascular_exchange`` without re-validating every absolute
    number in ``results/``.
    """
    assert SPLIT_DEFAULT is False
    assert PDEConfig().transvascular_exchange is False
    # and the module-level default really is what the keyword takes
    micro = _patch_micro()
    _, ex = delivery_pair(micro, DRUG, 5.0)
    assert ex is None


def test_split_off_is_bit_identical_to_delivery_source():
    micro = _patch_micro()
    for cp in (0.0, 1.0, 7.25):
        src, ex = delivery_pair(micro, DRUG, cp)
        assert ex is None
        ref = delivery_source(micro, DRUG, cp)
        assert np.array_equal(src, ref), f"split=False drifted at Cp={cp}"


def test_split_on_does_not_touch_the_source():
    """The split adds a back-flux; it must not silently restate the source."""
    micro = _patch_micro()
    ref = delivery_source(micro, DRUG, 3.5)
    src, ex = delivery_pair(micro, DRUG, 3.5, split=True)
    assert np.array_equal(src, ref)
    assert ex is not None


def test_split_on_carries_the_same_ps():
    """``exchange_k * Cp == source``, EXACTLY, for every argument combination.

    This is the check that catches both the INERT defect (``exchange_k`` never
    built) and the drift defect (source and back-flux built from different
    permeabilities, which ``drugs.exchange_coefficient`` warns about in prose).
    """
    micro = _patch_micro()
    for cp in (1.0, 4.0, 100.0):
        for disruption in (None, 0.0, 0.35):
            for efflux in (False, True):
                src, ex = delivery_pair(
                    micro, DRUG, cp, split=True,
                    bbb_disruption=disruption, efflux=efflux)
                assert ex is not None, "split=True must build an exchange field"
                assert np.array_equal(ex * cp, src), (cp, disruption, efflux)
                assert np.array_equal(
                    ex,
                    exchange_coefficient(micro, DRUG,
                                         bbb_disruption=disruption,
                                         efflux=efflux))
    # Cp = 0 is an ordinary point of a dosing profile and must stay finite.
    _, ex0 = delivery_pair(micro, DRUG, 0.0, split=True)
    assert np.all(np.isfinite(ex0)) and np.all(ex0 >= 0.0)


def test_exchange_field_is_non_negative_and_finite():
    micro = _patch_micro()
    _, ex = delivery_pair(micro, DRUG, 2.0, split=True)
    assert np.all(np.isfinite(ex))
    assert np.all(ex >= 0.0), "a negative exchange rate would pump drug inward"


# =========================================================================== #
# 2. The exact property the whole result rests on, asserted on the solver.
# =========================================================================== #
def test_one_way_law_is_exactly_homogeneous_in_the_entry_term():
    """Scaling PS scales the field and leaves the SHAPE invariant.

    Slab: 25 x 3 at dx = 8 um, ecm 0.3, dt = 2 s, 1200 s window, convection
    OFF, exchange OFF, through ``simulate_response``. ``bbb_penetrance`` is the
    only per-compound factor in ``PS``, so scaling it IS scaling the entry term.
    """
    base = _profile(DRUG, exchange=False)
    assert base[0] > 0.0
    shape0 = base / base[0]
    for c in (1e-4, 1e2, 1e5):
        scaled = _profile(
            dataclasses.replace(DRUG, bbb_penetrance=DRUG.bbb_penetrance * c),
            exchange=False)
        # the LEVEL moves by exactly c ...
        assert np.allclose(scaled, base * c, rtol=1e-11, atol=0.0), c
        # ... and the SHAPE does not move at all.
        assert np.allclose(scaled / scaled[0], shape0, rtol=1e-11, atol=0.0), c


def test_split_breaks_the_homogeneity():
    """With the split ON the shape MUST depend on the entry term.

    This is the check with a live failure mode: if the flag ever stops reaching
    the solver (``tmesim.validation.propagation``'s INERT category), the ON arm
    silently becomes the OFF arm, the shape becomes invariant again, and this
    assertion fires. A test that merely compared ON against OFF would not --
    it would compare a run to itself and pass.

    Same slab and same configuration as the test above, exchange ON.
    """
    lo = _profile(dataclasses.replace(DRUG, bbb_penetrance=DRUG.bbb_penetrance * 1e-4),
                  exchange=True)
    hi = _profile(dataclasses.replace(DRUG, bbb_penetrance=DRUG.bbb_penetrance * 1e4),
                  exchange=True)
    shape_lo = lo / lo[0]
    shape_hi = hi / hi[0]
    rel = float(np.max(np.abs(shape_hi - shape_lo) / np.maximum(shape_lo, 1e-300)))
    assert rel > 1e-3, (
        "the split did not change the gradient SHAPE across 8 decades of PS "
        f"(max relative shape change {rel:.3g}). Either exchange_k never "
        "reached the solver -- the INERT defect in "
        "tmesim.validation.propagation -- or the split is not wired at all."
    )
    # and the level must SATURATE rather than scale: 8 decades of PS must not
    # buy 8 decades of exposure.
    assert hi[0] / lo[0] < 1e7, (hi[0], lo[0])


def test_split_bounds_the_perivascular_level_and_the_one_way_law_does_not():
    """The plasma bound, on the solver rather than in closed form.

    The bound is ``C <= Cp``, so over the window the exposure of a tissue that
    exactly tracked plasma is ``mean_t(Cp) * total_time``. That is computed from
    the run's OWN ``PKConfig`` on the solver's OWN step grid rather than
    hard-coded: an earlier revision of this test assumed ``Cp = 1`` uM and
    failed at 58254 against 1200, which was the test being wrong and not the
    physics. Same slab and configuration as above; the default
    ``PKConfig(plasma_peak=50, half_life=14400, n_doses=1)`` gives a window mean
    of 48.586 uM over 1200 s.
    """
    open_wall = dataclasses.replace(DRUG, bbb_penetrance=DRUG.bbb_penetrance * 1e5)
    off = _profile(open_wall, exchange=False)
    on = _profile(open_wall, exchange=True)
    n_steps = int(round(SLAB_TOTAL_S / SLAB_DT_S))
    plasma_mean = float(plasma_profile(
        PKConfig(), np.arange(n_steps + 1, dtype=np.float64) * SLAB_DT_S).mean())
    plasma_auc = SLAB_TOTAL_S * plasma_mean
    assert off[0] > 10.0 * plasma_auc, (
        "the one-way law was expected to overshoot the plasma bound at a wide "
        f"open wall; it reached {off[0]:.4g} against a plasma AUC of "
        f"{plasma_auc:.4g}")
    assert on[0] <= plasma_auc * (1.0 + 1e-9), (
        f"the split law exceeded the plasma bound: {on[0]:.6g} > {plasma_auc:.6g}")


# =========================================================================== #
# 3. The closed forms that state what the split buys.
# =========================================================================== #
def test_entry_elasticity_is_exactly_one_under_the_one_way_law():
    for ps in (0.0, 1e-8, 1e-3, 1.0, 1e9):
        assert entry_elasticity(ps, 2e-3, split=False) == 1.0


def test_entry_elasticity_under_the_split_is_k_over_ps_plus_k():
    k = 2e-3
    for ps in (1e-6, 1e-4, 1e-3, k, 1e-1, 1.0, 1e3):
        assert entry_elasticity(ps, k, split=True) == pytest.approx(
            k / (ps + k), rel=1e-14)
    # monotone decreasing in PS, and it reaches zero
    vals = [entry_elasticity(p, k, split=True)
            for p in (1e-6, 1e-3, 1e-1, 1.0, 1e3, 1e9)]
    assert all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)), vals
    assert vals[0] == pytest.approx(1.0, rel=1e-2)
    assert vals[-1] < 1e-9


def test_saturation_is_bounded_and_the_one_way_multiple_is_not():
    k = 2e-3
    for ps in (0.0, 1e-6, 1e-3, 1.0, 1e6, 1e12):
        assert 0.0 <= saturation_fraction(ps, k) <= 1.0
    assert saturation_fraction(0.0, k) == 0.0
    assert one_way_plasma_multiple(1.0, k) == pytest.approx(500.0, rel=1e-12)
    assert one_way_plasma_multiple(1e6, k) > 1e8
    assert one_way_plasma_multiple(1.0, 0.0) == float("inf")


def test_boundary_layer_shortens_only_under_the_split():
    d, k = 40.0, 2e-3
    for ps in (1e-4, 1e-2, 1.0):
        assert boundary_layer_um(d, ps, k, split=False) == pytest.approx(
            np.sqrt(d / k), rel=1e-14)
    lens = [boundary_layer_um(d, p, k, split=True) for p in (0.0, 1e-3, 1e-1, 10.0)]
    assert lens[0] == pytest.approx(np.sqrt(d / k), rel=1e-14)
    assert all(lens[i] > lens[i + 1] for i in range(len(lens) - 1)), lens


def test_regime_diagnostic_moves_with_the_spacing_it_is_told():
    d, k = 40.0, 2e-3
    near = entry_transport_regime(d, 1e-3, k, spacing_um=10.0)
    far = entry_transport_regime(d, 1e-3, k, spacing_um=400.0)
    assert near["regime"] == "entry_limited"
    assert far["regime"] == "transport_limited"
    # quadratic in the spacing, and the spacing is reported so it cannot be lost
    assert far["damkohler"] / near["damkohler"] == pytest.approx(1600.0, rel=1e-12)
    assert near["spacing_um"] == 10.0 and far["spacing_um"] == 400.0


# =========================================================================== #
# 4. The decomposition PS = P * (S/V) and the parameter it exposes.
# =========================================================================== #
def test_ps_decomposition_round_trips_to_the_engine_field():
    """The implied permeability reproduces the engine's own PS field exactly.

    It must, and the fact that it does for EVERY radius is the point: the
    radius cancels, so the shipped law pins ``P`` only up to a vessel geometry
    nobody in this repo has measured.
    """
    micro = _patch_micro()
    engine_ps = exchange_coefficient(micro, DRUG)
    for radius in (1.5, 4.0, 12.0):
        p = implied_permeability_um_s(DRUG.bbb_penetrance, radius)
        assert p == pytest.approx(DRUG.bbb_penetrance * radius / 2.0, rel=1e-14)
        rebuilt = ps_from_permeability(p, micro.fields.vascular_density, radius)
        assert np.allclose(rebuilt, engine_ps, rtol=1e-12, atol=0.0)


def test_surface_to_volume_is_the_cylinder_geometry():
    phi = np.array([0.0, 0.05, 0.5])
    assert np.allclose(surface_to_volume_per_um(phi, 4.0), 2.0 * phi / 4.0)
    with pytest.raises(ValueError):
        surface_to_volume_per_um(phi, 0.0)
    with pytest.raises(ValueError):
        surface_to_volume_per_um(np.array([-0.1]), 4.0)


def test_rate_arguments_reject_negatives():
    for fn in (saturation_fraction, one_way_plasma_multiple):
        with pytest.raises(ValueError):
            fn(-1.0, 1e-3)
        with pytest.raises(ValueError):
            fn(1e-3, -1.0)
    with pytest.raises(ValueError):
        entry_elasticity(-1.0, 1e-3)
    with pytest.raises(ValueError):
        boundary_layer_um(0.0, 1e-3, 1e-3)
    with pytest.raises(ValueError):
        implied_permeability_um_s(-0.1, 4.0)
