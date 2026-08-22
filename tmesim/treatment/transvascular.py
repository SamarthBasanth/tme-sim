"""ENTRY versus SPREAD: separating "how much gets in" from "where it goes once in".

WHAT THIS MODULE IS FOR, AND THE RESULT THAT DEFINES ITS SCOPE
==============================================================
The shipped delivery law is one-way (``treatment.drugs.delivery_source``)::

    dC/dt = div(D_eff grad C) - div(vC) + PS*Cp - (clearance + binding)*C

``PS = vascular_density * bbb_penetrance * PS_SCALE_S_INV`` is the ENTRY term
and ``D_eff`` is the SPREAD term. The question this module was built to answer
is whether conflating them is part of the engine's ~14x perivascular-gradient
error against Kyle 2004. **It is not, and the reason is exact rather than
empirical.**

THE HOMOGENEITY RESULT — measured, and it settles the question
--------------------------------------------------------------
Under the one-way law the solved field is EXACTLY positively homogeneous of
degree 1 in the entry term::

    C(x, t; c*PS)  ==  c * C(x, t; PS)     for every c > 0

because ``PS`` enters the linear PDE only as a multiplicative factor on the
source and nowhere else. Measured on the frozen solver, arm-A slab (140 um,
dx = 2 um, ecm 0.0875, 1 h window, convection OFF, exchange OFF, routed through
``sim.engine.simulate_response``): sweeping ``c`` over twelve decades
(1e-6 to 1e+6) moves the perivascular AUC exactly proportionally
(1.850775e-02 -> 1.850775e+10) while the near:far ratio stays
2.997687165949104 -> 2.9976871659491144, i.e. invariant to ~1e-14.

**The consequence, stated in the direction that matters.** The entry term
contributes NOTHING to gradient shape. It cannot be part of a shape error of
any size. What sets shape is ``lambda = sqrt(D_eff / (clearance + binding))``,
and that same pair also sets the perivascular LEVEL once ``PS`` has fixed the
scale. So the parameter doing two jobs is ``D_eff``/``k_loss``, not ``PS``.

**And the structural reason the split cannot rescue the gradient.** ``PS`` is
supported only where ``vascular_density > 0``. Kyle's gradient is measured
where it is not. Turning the exchange term on adds a sink at the vessels and
changes nothing in the tissue interior: same slab, same dx/window/convection,
sweeping PS_wall from 2.06e-6 to 2.06e+3 1/s with
``PDEConfig.transvascular_exchange=True`` moves near:far from 2.9977 to 2.9931
-- 0.15 %, in the SHALLOWER direction, against a gap of 14x. The split is
structurally incapable of reaching the tissue interior. See
``results/permeability_split.md`` sections 2 and 3.

WHAT THE SPLIT DOES BUY, WHICH IS A DIFFERENT THING ENTIRELY
-------------------------------------------------------------
Homogeneity of degree 1 is not a harmless property. It says that every
ABSOLUTE quantity the engine reports is proportional to ``PS``, and ``PS``
carries an implicit, underived factor of exactly 1 s^-1
(``treatment.drugs.PS_SCALE_S_INV``, and
``results/transvascular_exchange_revalidation.md`` section 6.2 for the
dimensional audit). So:

* every RATIO WITHIN ONE MAP is exactly invariant to that unmeasured scale --
  these are the claims that survive;
* every ABSOLUTE THRESHOLD verdict -- per-cell DOSE / DELIVERY / RESISTANCE,
  which compares local exposure against ``IC50`` -- scales with it, without
  bound.

The exchange (split) law ``S = PS*(Cp - C)`` replaces the unbounded scaling
with a SATURATING one: ``C_ss = PS*Cp/(PS + k_loss) <= Cp``. The elasticity of
exposure with respect to the entry term drops from exactly 1 to
``k_loss/(PS + k_loss)``, which goes to 0 as the wall opens. That is the whole
case for the split, and it is a claim about FALSIFIABILITY -- an unmeasured
parameter stops propagating without bound into the one class of output a user
reads as an answer. :func:`entry_elasticity` and :func:`saturation_fraction`
are those two statements as code.

HONEST SCOPE — the parameter this trades for
---------------------------------------------
A genuine transvascular permeability is ``PS [1/s] = P [um/s] * (S/V) [1/um]``:
a wall permeability times a vessel surface-to-volume density.
:func:`ps_from_permeability` and :func:`implied_permeability_um_s` write that
decomposition down. **The repository has neither factor.** It has
``bbb_penetrance``, a dimensionless brain:plasma PARTITION coefficient that
three separate modules state is not a permeability, multiplied by an
undocumented 1 s^-1. And ``S/V`` needs a vessel radius that this repo does not
carry as a measured quantity. So the split as implemented does NOT introduce a
measured entry parameter -- it makes the unmeasured one saturate instead of
diverge. That is a real gain and it is a smaller one than "we now model entry
properly". Read :func:`implied_permeability_um_s`'s docstring before quoting
any ``P``.

WHAT THIS MODULE DOES NOT DO
-----------------------------
It never calls the transport solver. It builds fields and returns closed-form
scalars; the caller wires them. That is deliberate: ``tmesim/sim/engine.py`` is
and must remain the ONLY production place that turns
``PDEConfig.transvascular_exchange`` into an ``exchange_k`` argument
(``tests/test_propagation.py::test_there_is_exactly_one_production_wiring_point``).
Adding a second wiring point here would be the exact defect
``tmesim.validation.propagation`` exists to prevent.

DEFAULT: OFF. :data:`SPLIT_DEFAULT` is ``False`` and :func:`delivery_pair` with
``split=False`` returns ``delivery_source(...)`` bit-for-bit and ``None`` for
the exchange field, so this module changes no shipped number until a caller
opts in. Whether it SHOULD be on by default is a recommendation with evidence
in ``results/permeability_split.md`` section 7; it is not switched here.

Units: concentration uM, time s, length um, permeability um/s, PS 1/s.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from tmesim.core import DrugParams, Microenvironment

from .drugs import delivery_source, exchange_coefficient

__all__ = [
    "SPLIT_DEFAULT",
    "delivery_pair",
    "entry_elasticity",
    "saturation_fraction",
    "one_way_plasma_multiple",
    "ps_from_permeability",
    "implied_permeability_um_s",
    "surface_to_volume_per_um",
    "boundary_layer_um",
    "entry_transport_regime",
]

#: The shipped default for the entry/spread split. ``False`` -- the one-way law.
#: This mirrors ``tmesim.core.PDEConfig.transvascular_exchange``, which is also
#: ``False``, and it is NOT independently switchable in production: the engine
#: reads the PDEConfig flag, not this one. This constant exists so a caller of
#: :func:`delivery_pair` that is outside the engine (an experiment, a probe)
#: has a single named place to read the shipped position from rather than
#: hard-coding ``False`` and drifting.
SPLIT_DEFAULT: bool = False


# --------------------------------------------------------------------------- #
# 1. The field builder. OFF is bit-identical to the shipped law.
# --------------------------------------------------------------------------- #
def delivery_pair(
    micro: Microenvironment,
    drug: DrugParams,
    plasma_conc: float,
    *,
    split: bool = SPLIT_DEFAULT,
    bbb_disruption: "np.ndarray | float | None" = None,
    efflux: bool = False,
    ps_scale: "float | None" = None,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """The ``(source, exchange_k)`` pair for one delivery law, in one call.

    The two fields MUST be built from the same ``PS``. Building them separately
    -- ``delivery_source`` here, ``exchange_coefficient`` three functions away --
    is how a source and a back-flux end up using different permeabilities, which
    ``drugs.exchange_coefficient``'s own docstring already warns about. This
    function makes them one object so they cannot drift.

    Parameters
    ----------
    split : bool
        ``False`` (the shipped default, :data:`SPLIT_DEFAULT`) returns
        ``(delivery_source(...), None)``. The source is the SAME OBJECT the
        shipped path builds -- same call, same arguments -- and ``None`` is what
        ``simulate_drug_transport(exchange_k=...)`` short-circuits on, so the
        default arm is bit-for-bit the one-way law.
        ``True`` returns ``(PS*Cp, PS)``, the gradient-driven exchange law
        ``S = PS*(Cp - C)`` written as source-minus-first-order-loss.
    bbb_disruption, efflux, ps_scale
        Passed through unchanged to both fields, which is the point.

    Returns
    -------
    (source, exchange_k)
        ``source`` is ``(ny, nx)`` in uM/s. ``exchange_k`` is ``(ny, nx)`` in
        1/s when ``split`` is True, and ``None`` otherwise.

    Notes
    -----
    Because ``delivery_source`` is exactly linear in ``plasma_conc``,
    ``exchange_k`` here is numerically ``source / plasma_conc`` wherever
    ``plasma_conc != 0`` -- but it is COMPUTED as
    ``exchange_coefficient(...)``, i.e. the same expression at
    ``plasma_conc = 1``, rather than by dividing. Dividing would produce
    ``nan``/``inf`` at ``plasma_conc = 0``, which is a perfectly ordinary point
    of a dosing profile (before the first dose, and between doses of a drug that
    has washed out), and a ``nan`` in ``exchange_k`` is rejected by the solver.
    """
    source = delivery_source(
        micro, drug, plasma_conc,
        bbb_disruption=bbb_disruption, efflux=efflux, ps_scale=ps_scale,
    )
    if not split:
        return source, None
    ps = exchange_coefficient(
        micro, drug,
        bbb_disruption=bbb_disruption, efflux=efflux, ps_scale=ps_scale,
    )
    return source, ps


# --------------------------------------------------------------------------- #
# 2. What the split changes, in closed form. These are the claims.
# --------------------------------------------------------------------------- #
def entry_elasticity(ps_s_inv: float, k_loss_s_inv: float,
                     split: bool = SPLIT_DEFAULT) -> float:
    """``d ln(exposure) / d ln(PS)`` -- how far an entry-term error propagates.

    This is the number that says how much an unmeasured entry parameter matters
    to an ABSOLUTE output, and the two laws give categorically different answers.

    Well-mixed (zero-dimensional) steady state, which is the limit in which the
    two laws can be compared with no geometry in the way:

    * one-way   ``C_ss = PS*Cp / k``            -> elasticity **exactly 1**
    * split     ``C_ss = PS*Cp / (PS + k)``     -> elasticity ``k/(PS + k) < 1``

    Elasticity 1 means a factor-``f`` error in ``PS`` is a factor-``f`` error in
    every absolute exposure, forever, with no bound. That is the shipped
    position, and ``PS`` carries an underived 1 s^-1
    (``drugs.PS_SCALE_S_INV``). Under the split the same error is damped by
    ``k/(PS+k)`` and vanishes entirely once the wall is open compared with the
    tissue's own loss rate.

    **This is a statement about the well-mixed limit, not about the solved
    field.** In a real geometry the interstitium far from a vessel sees
    ``PS = 0`` and is therefore still elasticity-1 in the ENTRY of its feeding
    vessel; the damping is local to where ``PS`` is supported. Do not quote this
    as a whole-tissue number -- quote it as the perivascular one it is. The
    solved-field version is measured in ``results/permeability_split.md``
    section 4, not asserted here.

    Returns 1.0 for ``split=False`` regardless of ``ps_s_inv``, which is the
    exact homogeneity result and not an approximation.
    """
    k = float(k_loss_s_inv)
    if k < 0.0:
        raise ValueError(f"k_loss_s_inv must be >= 0, got {k_loss_s_inv!r}")
    ps = float(ps_s_inv)
    if ps < 0.0:
        raise ValueError(f"ps_s_inv must be >= 0, got {ps_s_inv!r}")
    if not split:
        return 1.0
    total = ps + k
    if total <= 0.0:
        # Both rates zero: nothing enters and nothing leaves. The derivative is
        # not defined; 1.0 is the one-way answer and is the conservative choice
        # (it claims the MOST propagation, never less).
        return 1.0
    return k / total


def saturation_fraction(ps_s_inv: float, k_loss_s_inv: float) -> float:
    """``C_ss / Cp`` under the SPLIT law: ``PS/(PS + k)``, always in [0, 1].

    The tissue:plasma ratio the split law permits at steady state in the
    well-mixed limit. Bounded by 1 by construction, for every drug at every
    dose -- which is the property the one-way law does not have and the
    property that
    ``results/breast_kp_validation.md`` measured it failing (5-fluorouracil
    predicted at 7.80x plasma against 0.61x measured).
    """
    ps = float(ps_s_inv)
    k = float(k_loss_s_inv)
    if ps < 0.0 or k < 0.0:
        raise ValueError(f"rates must be >= 0, got ps={ps!r}, k={k!r}")
    total = ps + k
    return 0.0 if total <= 0.0 else ps / total


def one_way_plasma_multiple(ps_s_inv: float, k_loss_s_inv: float) -> float:
    """``C_ss / Cp`` under the ONE-WAY law: ``PS/k``, unbounded above.

    Provided next to :func:`saturation_fraction` so the two can be read side by
    side. Returns ``inf`` at zero loss, which is the correct statement that a
    one-way source into a lossless tissue has no steady state at all.
    """
    ps = float(ps_s_inv)
    k = float(k_loss_s_inv)
    if ps < 0.0 or k < 0.0:
        raise ValueError(f"rates must be >= 0, got ps={ps!r}, k={k!r}")
    if k <= 0.0:
        return float("inf") if ps > 0.0 else 0.0
    return ps / k


def boundary_layer_um(d_eff_um2_s: float, ps_s_inv: float,
                      k_loss_s_inv: float,
                      split: bool = SPLIT_DEFAULT) -> float:
    """Perivascular decay length ``sqrt(D_eff / (PS + k))``, um.

    Under the one-way law (``split=False``) ``PS`` is absent from the sink and
    this is ``sqrt(D_eff/k)`` -- the length the GBM penetration anchor is fit
    to, and the reason a scale error in ``PS`` of ANY magnitude was invisible to
    that anchor (``tmesim/validation/ps_scale.py``, "WHY IT WAS INVISIBLE").
    Under the split the wall's own conductance shortens it. The discretisation
    consequence -- this length falling below ``dx`` -- is measured in
    ``experiments/boundary_layer_resolution.py``, which owns that question.
    """
    d = float(d_eff_um2_s)
    if d <= 0.0:
        raise ValueError(f"d_eff_um2_s must be positive, got {d_eff_um2_s!r}")
    k = float(k_loss_s_inv)
    ps = float(ps_s_inv)
    if k < 0.0 or ps < 0.0:
        raise ValueError(f"rates must be >= 0, got ps={ps!r}, k={k!r}")
    total = k + (ps if split else 0.0)
    if total <= 0.0:
        return float("inf")
    return float(np.sqrt(d / total))


def entry_transport_regime(d_eff_um2_s: float, ps_s_inv: float,
                           k_loss_s_inv: float,
                           spacing_um: float) -> Dict[str, float]:
    """Is delivery ENTRY-limited or TRANSPORT-limited, at this vessel spacing?

    The dimensionless group is the ratio of the two resistances a molecule
    meets: crossing the wall, then crossing the tissue to the next cell.

        R_entry     ~ 1 / PS                      [s]
        R_transport ~ spacing^2 / D_eff           [s]
        Da          = R_transport / R_entry = PS * spacing^2 / D_eff

    ``Da >> 1`` -- the wall is fast compared with the tissue, so the tissue is
    the bottleneck: TRANSPORT-limited, and an error in ``PS`` barely matters.
    ``Da << 1`` -- the wall is the bottleneck: ENTRY-limited, and an error in
    ``PS`` is the whole answer.

    **This is a diagnostic, not a claim.** ``spacing_um`` is the caller's choice
    of what "the tissue" means (a Krogh radius, a mean intercapillary distance,
    a slab thickness) and the verdict moves with it, quadratically. It is
    reported with the number it was computed from so that cannot be lost.
    """
    d = float(d_eff_um2_s)
    if d <= 0.0:
        raise ValueError(f"d_eff_um2_s must be positive, got {d_eff_um2_s!r}")
    s = float(spacing_um)
    if s <= 0.0:
        raise ValueError(f"spacing_um must be positive, got {spacing_um!r}")
    ps = float(ps_s_inv)
    if ps < 0.0:
        raise ValueError(f"ps_s_inv must be >= 0, got {ps_s_inv!r}")
    da = ps * s * s / d
    return {
        "damkohler": da,
        "spacing_um": s,
        "ps_s_inv": ps,
        "d_eff_um2_s": d,
        "k_loss_s_inv": float(k_loss_s_inv),
        "entry_resistance_s": (float("inf") if ps <= 0.0 else 1.0 / ps),
        "transport_resistance_s": s * s / d,
        "regime": ("entry_limited" if da < 1.0 else "transport_limited"),
    }


# --------------------------------------------------------------------------- #
# 3. The decomposition PS = P * (S/V), and the parameter it exposes.
# --------------------------------------------------------------------------- #
def surface_to_volume_per_um(vascular_density, vessel_radius_um: float):
    """Vessel surface area per unit tissue volume, 1/um: ``2*phi / R``.

    For cylinders of radius ``R`` occupying volume fraction ``phi``, the surface
    per unit tissue volume is ``(2*pi*R*L)/(pi*R^2*L) * phi = 2*phi/R``. Pure
    geometry, no fitting.

    ``vessel_radius_um`` IS NOT MEASURED IN THIS REPOSITORY. It is the caller's
    parameter and there is no default, deliberately: a default would be quoted
    as though it came from somewhere. See :func:`implied_permeability_um_s`.
    """
    r = float(vessel_radius_um)
    if r <= 0.0:
        raise ValueError(f"vessel_radius_um must be positive, got {r!r}")
    phi = np.asarray(vascular_density, dtype=np.float64)
    if np.any(phi < 0.0):
        raise ValueError("vascular_density must be non-negative")
    return 2.0 * phi / r


def ps_from_permeability(permeability_um_s: float, vascular_density,
                         vessel_radius_um: float):
    """``PS [1/s] = P [um/s] * (S/V) [1/um]`` -- the dimensionally honest entry term.

    This is what the engine's ``PS`` field WOULD be if the repo had a measured
    wall permeability. It does not; see :func:`implied_permeability_um_s`. The
    function exists so that the moment a measured ``P`` appears, the wiring is
    already written and the units are already right.
    """
    p = float(permeability_um_s)
    if p < 0.0:
        raise ValueError(f"permeability_um_s must be >= 0, got {p!r}")
    return p * surface_to_volume_per_um(vascular_density, vessel_radius_um)


def implied_permeability_um_s(bbb_penetrance: float, vessel_radius_um: float,
                              ps_scale_s_inv: float = 1.0) -> float:
    """The wall permeability the SHIPPED law implicitly asserts, um/s.

    Inverting ``PS = phi * bbb_penetrance * ps_scale`` against
    ``PS = P * 2*phi/R`` gives ``P = bbb_penetrance * ps_scale * R / 2``, and
    the vascular density ``phi`` cancels exactly -- so the implied permeability
    is a per-compound constant, independent of where in the tissue you evaluate
    it. That cancellation is the only reason this number is well defined at all.

    **WHAT THIS NUMBER IS, SAID PLAINLY.** It is not a measurement and it is not
    a prediction. It is the answer to "if the shipped delivery law were a
    permeability law, what permeability would it be claiming?", and it is
    reported so that the claim is at least VISIBLE. Every factor in it is
    contested: ``bbb_penetrance`` is a dimensionless brain:plasma PARTITION
    coefficient that ``core/config.py``, ``reconstruct/barrier.py`` and
    ``molecular/permeability.py`` all state is not a permeability;
    ``ps_scale`` is an undocumented 1 s^-1 by default
    (``drugs.PS_SCALE_S_INV``); and ``vessel_radius_um`` is supplied by the
    caller because this repository does not carry a measured one.

    **Do not quote the output of this function as the engine's permeability.**
    Quote it as what the engine would have to mean, and quote the radius you
    passed in the same sentence.
    """
    r = float(vessel_radius_um)
    if r <= 0.0:
        raise ValueError(f"vessel_radius_um must be positive, got {r!r}")
    b = float(bbb_penetrance)
    if b < 0.0:
        raise ValueError(f"bbb_penetrance must be >= 0, got {b!r}")
    s = float(ps_scale_s_inv)
    if s < 0.0:
        raise ValueError(f"ps_scale_s_inv must be >= 0, got {s!r}")
    return b * s * r / 2.0


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    from tmesim.core import (
        CellPopulation, CellState, CellType, Domain, TissueFields,
    )

    print("=" * 78)
    print("transvascular self-test: entry/spread split")
    print("=" * 78)

    dom = Domain(nx=12, ny=9, dx=15.0)
    vasc = np.zeros(dom.shape)
    vasc[:, :6] = 0.05
    vasc[:, 6:] = 0.40
    fields = TissueFields(ecm_density=dom.ones() * 0.3, vascular_density=vasc)
    cells = CellPopulation(
        x=np.array([10.0]), y=np.array([10.0]),
        cell_type=np.array([CellType.TUMOR]), state=np.array([CellState.VIABLE]),
        sensitivity=np.array([0.5]))
    micro = Microenvironment(domain=dom, fields=fields, cells=cells)
    drug = DrugParams(name="probe", bbb_penetrance=0.2, diffusivity=30.0,
                      clearance=3e-4, binding=2e-3)

    # [1] OFF is bit-identical to the shipped one-way law.
    src_off, ex_off = delivery_pair(micro, drug, 12.0)
    ref = delivery_source(micro, drug, 12.0)
    assert ex_off is None
    assert np.array_equal(src_off, ref), "split=False must be bit-identical"
    assert SPLIT_DEFAULT is False
    print(f"\n[1] split=False -> exchange_k is None and source is bit-identical "
          f"to delivery_source (max|diff| = {np.abs(src_off - ref).max():.17g})")

    # [2] ON returns the SAME PS in both slots.
    src_on, ex_on = delivery_pair(micro, drug, 12.0, split=True)
    assert np.array_equal(src_on, ref), "the source must not change with the flag"
    assert np.allclose(ex_on * 12.0, src_on, rtol=0, atol=0), "PS*Cp != source"
    print(f"[2] split=True  -> source unchanged; exchange_k*Cp == source exactly; "
          f"PS range {ex_on.min():.4g}-{ex_on.max():.4g} 1/s")

    # [3] the elasticity claim, which is the whole case for the split
    k = drug.clearance + drug.binding
    print(f"\n[3] entry elasticity d ln(C_ss)/d ln(PS), k_loss = {k:.4g} 1/s")
    print(f"    {'PS (1/s)':>12}{'one-way':>10}{'split':>10}"
          f"{'C_ss/Cp one-way':>18}{'C_ss/Cp split':>16}")
    for ps in (1e-6, 1e-4, 1e-3, k, 1e-1, 1.0, 1e2):
        e0 = entry_elasticity(ps, k, split=False)
        e1 = entry_elasticity(ps, k, split=True)
        assert e0 == 1.0
        assert 0.0 <= e1 <= 1.0
        assert saturation_fraction(ps, k) <= 1.0
        print(f"    {ps:>12.4g}{e0:>10.4f}{e1:>10.4f}"
              f"{one_way_plasma_multiple(ps, k):>18.4g}"
              f"{saturation_fraction(ps, k):>16.4f}")
    # monotone: a more open wall damps the entry error harder
    es = [entry_elasticity(p, k, split=True) for p in (1e-6, 1e-3, 1e-1, 1.0, 1e3)]
    assert all(es[i] >= es[i + 1] for i in range(len(es) - 1)), es
    assert entry_elasticity(1e12, k, split=True) < 1e-6
    print("    -> one-way is elasticity 1 at EVERY PS; split -> 0 as the wall opens")

    # [4] the one-way law has no plasma bound and the split does
    assert one_way_plasma_multiple(1.0, k) > 100.0
    assert saturation_fraction(1e12, k) <= 1.0
    print(f"[4] at PS = 1 1/s the one-way law asserts "
          f"{one_way_plasma_multiple(1.0, k):.1f}x plasma; the split asserts "
          f"{saturation_fraction(1.0, k):.4f}x, bounded by 1 by construction")

    # [5] boundary layer and regime
    d_eff = 30.0
    l_off = boundary_layer_um(d_eff, 0.02, k, split=False)
    l_on = boundary_layer_um(d_eff, 0.02, k, split=True)
    assert l_on < l_off
    print(f"[5] decay length at D_eff = {d_eff} um2/s, PS = 0.02 1/s: "
          f"one-way {l_off:.2f} um -> split {l_on:.2f} um")
    reg = entry_transport_regime(d_eff, 0.02, k, spacing_um=50.0)
    print(f"    Damkohler at 50 um spacing = {reg['damkohler']:.3g} "
          f"-> {reg['regime']}")

    # [6] the decomposition and the parameter it exposes
    print("\n[6] PS = P * (S/V). The implied P is independent of phi:")
    for radius in (2.0, 4.0, 8.0):
        p_imp = implied_permeability_um_s(drug.bbb_penetrance, radius)
        back = ps_from_permeability(p_imp, vasc, radius)
        direct = exchange_coefficient(micro, drug)
        assert np.allclose(back, direct, rtol=1e-12, atol=0), (back, direct)
        print(f"    R = {radius:>4.1f} um -> implied P = {p_imp:.5g} um/s "
              f"(round-trip to PS: max rel err "
              f"{np.max(np.abs(back - direct) / np.maximum(direct, 1e-300)):.3g})")
    print("    the radius is the CALLER's number; this repo does not measure one")

    print("\nAll asserts passed.")
