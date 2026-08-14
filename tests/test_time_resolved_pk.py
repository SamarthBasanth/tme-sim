"""Time-resolved pharmacokinetics: the drug layer's schedule capability.

Before this, ``sim.engine.simulate_response`` collapsed the multi-dose plasma
profile to a single time-averaged level and drove the transport solver with a
STATIC source. Two schedules with equal cumulative plasma exposure were then
bit-identical inputs, so fractionation and infusion duration produced a kill
spread of exactly 0.0 — a structural inability to answer a scheduling question,
not a finding about biology (see ``results/dose_schedule.md``).

``PDEConfig.time_resolved=True`` (and, at the solver level, the ``source_profile``
argument of :func:`tmesim.physics.diffusion.simulate_drug_transport`) drives the
source with the real ``plasma(t)`` instead. This module pins:

**Backward compatibility (must be BIT-identical, not merely close)**
  1. the default static path reproduces a stored pre-change reference exactly,
  2. an all-ones ``source_profile`` is bit-identical to ``source_profile=None``,
  3. ``PDEConfig.time_resolved`` defaults to False.

**The physics (checks that fail if the implementation is wrong)**
  4. *Convergence.* As the plasma level is held constant, the time-resolved path
     converges to the static-source result.
  5. *Mass balance.* Conservation still holds exactly with a time-varying source,
     and matched-cumulative-exposure schedules inject exactly the same total drug.
  6. *Depth / direction.* At matched total exposure over a fixed window the
     schedule genuinely changes the exposure field, and it does so in a
     DERIVABLE direction — see :func:`test_earlier_delivery_reaches_deeper`.
  7. *Capability.* Fractionation at matched cumulative exposure no longer gives a
     spread of exactly zero.

Why check 6 goes the opposite way to the naive expectation
----------------------------------------------------------
The usual transport-limited intuition says a long slow infusion should reach
DEEPER than a short bolus. For this (linear) transport model that is provably
false, and the proof also says exactly how big the schedule effect can be.
Integrate ``dC/dt = div(D grad C) - div(vC) + S - kC`` over the window ``[0, T]``
and write ``A = AUC = int_0^T C dt``::

    C(T) - C(0) = div(D grad A) - div(v A) - k A + int_0^T S dt

With ``C(0) = 0`` and ``M := k - div(D grad .) + div(v .)``::

    A = M^{-1} ( S_total - C(T) ),      S_total = int_0^T S dt

``S_total`` is IDENTICAL for any two schedules at matched cumulative plasma
exposure. So the entire schedule effect flows through ``C(T)``, the drug still
in flight at the end of the window. A bolus given early has washed out by ``T``
(``C(T) ~ 0``) and therefore realises the full steady-state exposure field; a
sustained infusion still running at ``T`` leaves ``C(T) > 0`` un-counted, so it
loses exposure. Two consequences, both tested here:

* Earlier delivery yields MORE exposure and, because the un-counted term
  ``M^{-1}C(T)`` is itself a doubly-smoothed (flatter) field, also a relatively
  FLATTER profile — i.e. the bolus reaches deeper, not the infusion.
* The effect is bounded by the terminal residual and therefore VANISHES as the
  window grows past washout. Schedule matters exactly insofar as the window
  truncates the response.

The 2026-08-12 re-pin of ``tests/data/static_source_reference.json``
-------------------------------------------------------------------
The stored reference exists so that nobody can change static-source behaviour
without it showing up. That only works if a re-pin is loud. This one was, and is
recorded here as well as in the reference file's own ``_README._repin_history``.

**Cause.** ``tmesim/physics/pressure.py`` smoothed the tumour-occupancy field
``phi`` with ``_TUMOR_SIGMA = 5.0``, commented "grid cells" but handed straight
to ``scipy.ndimage.gaussian_filter``, whose ``sigma`` is in PIXELS. ``phi`` sets
the vessel-leak gain and the lymphatic suppression, so it is a MODEL parameter:
the physical smoothing length was therefore ``5 * dx`` — 60 um at the production
grid but 75 um at ``dx = 15`` and 400 um at ``dx = 80``. Refining the mesh
silently redefined the model. It is now ``_TUMOR_SIGMA_UM = 60.0``, converted to
pixels at the point of use.

**What moved, and why that is the predicted result rather than a regression.**
``phi`` reaches a transport result by exactly one route — the Darcy velocity
``v = -K grad p`` — so it can only bite when convection is enabled. That is
precisely the observed split: the three ``dx = 15`` convection-ON cases and the
``dx = 80`` real-tissue case moved (by 1-3 ULP in the transport sums; more in
``penetration_index``, which is a shape statistic); ``weak_10dose_noconv``
(convection OFF) and ``solver_flat`` (raw solver, no IFP at all) are byte-for-byte
unchanged and were not regenerated. ``plasma_mean``, a pure-PK quantity with no
spatial dependence, is unchanged in every case. The superseded values are kept
verbatim in the reference file under
``_superseded_2026-08-12_grid_cell_tumor_sigma``.

**The guard was strengthened, not weakened.** ``60.0 / 12.0 == 5.0`` exactly in
IEEE-754, so the production grid is the fixed point of the change and every
shipped ``dx = 12`` result is bit-identical — but until now ``dx = 12`` was the
one grid with *no* bitwise pin, while the grids that did have pins were the ones
whose values were wrong. :func:`test_production_dx12_static_path_matches_stored_reference_bitwise`
closes that hole. Its values were generated BEFORE the fix and verified
identical AFTER it, so it pins the invariance claim itself. Full write-up:
``results/tumor_sigma_units.md``.

The 2026-08-13 re-pin — the two-route transport model became the default
-----------------------------------------------------------------------
**Cause.** ``compile_drug.build_profile`` now applies a two-route cell term
(paracellular between cells plus transcellular through membranes), so
temozolomide's compiled diffusivity went ``277.985 -> 42.4245 um^2/s`` and the
single perivascular-sink constant ``alpha``, fitted in closed form from
``L = sqrt(D_eff/k_eff)`` against the same measured 64.92 um anchor, went
``148.60815765445545 -> 22.223477711172563``. Both cases below use a COMPILED,
alpha-calibrated temozolomide, so both had to move.

**Why the new values are right rather than merely different.** The superseded
values encode a diffusivity law whose only solute input was molecular weight,
which was scored blind against 71 measured compounds and REFUTED at R^2 = 0.022
(``results/diffusivity_falsification.md``). The law now in force scores
R^2 = 0.718 on the same 71 compounds, blind and unfitted, and brings the engine's
dimensionless multicellular-layer hindrance factor inside its measured band for
the first time (1.77-fold -> 10.33-fold against a measured 12-fold, PMID
9806526). These cases now run a transport coefficient that has been tested
externally and largely agrees; they did not before.

**What moved, and why that is the predicted result.** Only a case whose drug
comes from the molecular compile can move, and only those did:
``real_xenium_dx80_tmz_3dose`` and ``production_xenium_dx12_tmz_3dose`` were
regenerated; ``default_1dose``, ``weak_1dose``, ``weak_3dose``,
``weak_10dose_noconv`` and ``solver_flat`` all use hand-built ``DrugParams``
literals with no compile in the path and are byte-for-byte unchanged; and
:func:`test_production_dx12_ifp_solve_matches_stored_reference_bitwise`, which
has no drug in it at all, passed untouched. ``plasma_mean`` — pure PK, no spatial
dependence — is BIT-IDENTICAL in the cases that moved, which is the positive
control.

**The magnitudes are a finding, not a footnote.** Unlike the 2026-08-12 ULP-level
re-pin, ``kill_fraction`` moves 2.76x at ``dx = 12`` and 5.29x at ``dx = 80``.
Diffusivity FELL 6.6x and exposure ROSE, because the Primeau anchor pins the
ratio ``D/k`` and the fitted loss rate fell with ``D``. The p99-normalised ratio
the breast panel scores is immune to that common-mode change — 0 of 12 of those
verdicts moved — but threshold statistics such as kill fraction are not. Full
write-up: ``results/two_route_promotion.md``; superseded values under
``_superseded_2026-08-13_mw_only_diffusivity``.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pytest

from tmesim.core import (
    CellPopulation,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    PKConfig,
    SimConfig,
    TissueFields,
)
from tmesim.physics.diffusion import simulate_drug_transport
from tmesim.sim.engine import _build_selftest_micro, simulate_response
from tmesim.treatment.drugs import delivery_source, plasma_profile

_REF_PATH = os.path.join(os.path.dirname(__file__), "data",
                         "static_source_reference.json")

# A deliberately weak drug so the kill fraction sits mid-curve and is genuinely
# sensitive to the exposure field (the library default saturates on this micro).
WEAK_DRUG = DrugParams(name="weak", diffusivity=200.0, ic50=4000.0, hill=1.2,
                       bbb_penetrance=0.15, clearance=3e-3, binding=0.0)

#: The exact cases captured in the stored reference, BEFORE the time-resolved
#: source existed. Any change to these definitions invalidates the reference.
_REFERENCE_CASES = [
    ("default_1dose", None,
     PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=1,
              dose_interval=3600.0 * 12),
     PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=200)),
    ("weak_1dose", WEAK_DRUG,
     PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=1,
              dose_interval=3600.0 * 12),
     PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=200)),
    ("weak_3dose", WEAK_DRUG,
     PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=3,
              dose_interval=3600.0 * 2),
     PDEConfig(dt=20.0, total_time=3600.0 * 6, convection=True, save_every=50)),
    ("weak_10dose_noconv", WEAK_DRUG,
     PKConfig(plasma_peak=20.0, half_life=1800.0, n_doses=10,
              dose_interval=1200.0),
     PDEConfig(dt=10.0, total_time=3600.0 * 4, convection=False, save_every=100)),
]


@pytest.fixture(scope="module")
def reference():
    with open(_REF_PATH) as f:
        return json.load(f)


def _hex(x: str) -> float:
    """Reference values are stored as IEEE-754 hex so the comparison is exact."""
    return float.fromhex(x)


def _slab_micro(nx: int = 160, ny: int = 3, dx: float = 10.0):
    """A 1-D transport slab: a single vessel column on the left edge.

    Deliberately simple geometry so the perivascular-to-deep profile is a clean
    1-D decay and a "depth" statistic means what it says.
    """
    domain = Domain(nx=nx, ny=ny, dx=dx)
    vasc = np.zeros((ny, nx))
    vasc[:, 0] = 1.0
    fields = TissueFields(ecm_density=np.zeros((ny, nx)), vascular_density=vasc)
    cells = CellPopulation(
        x=np.array([nx * dx / 2]), y=np.array([ny * dx / 2]),
        cell_type=np.array([0]), state=np.array([0]),
        sensitivity=np.array([0.5]),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells)


def _relative_depth_um(auc_row: np.ndarray, dx: float, frac: float = 0.1) -> float:
    """Distance at which AUC has fallen to ``frac`` of its perivascular value.

    Normalised by the perivascular value, so this measures the SHAPE of the
    profile, not its amplitude — raising the dose cannot move it.
    """
    target = auc_row[0] * frac
    below = np.nonzero(auc_row < target)[0]
    assert below.size, "profile never falls below the depth threshold"
    i = int(below[0])
    a, b = float(auc_row[i - 1]), float(auc_row[i])
    return (i - 1 + (a - target) / (a - b)) * dx


def _matched_peak(target_mean: float, half_life: float, total_time: float,
                  dt: float) -> float:
    """Per-dose peak that makes the window-mean plasma equal ``target_mean``."""
    n = int(round(total_time / dt))
    unit = plasma_profile(
        PKConfig(plasma_peak=1.0, half_life=half_life, n_doses=1,
                 dose_interval=total_time),
        np.arange(n, dtype=float) * dt,
    ).mean()
    return float(target_mean) / float(unit)


# --------------------------------------------------------------------------- #
# (1-3) Backward compatibility — BIT-identical, not approximately equal
# --------------------------------------------------------------------------- #
def test_time_resolved_defaults_to_off():
    """The new capability must be opt-in, or every stored result in the repo
    silently changes meaning."""
    assert PDEConfig().time_resolved is False
    assert SimConfig().pde.time_resolved is False


@pytest.mark.parametrize("label,drug,pk,pde", _REFERENCE_CASES,
                         ids=[c[0] for c in _REFERENCE_CASES])
def test_static_path_matches_stored_reference_bitwise(reference, label, drug,
                                                      pk, pde):
    """The default (``time_resolved=False``) engine path reproduces values
    captured BEFORE the time-resolved source was added, bit for bit."""
    micro = _build_selftest_micro()
    kw = dict(pk=pk, pde=pde, seed=0)
    if drug is not None:
        kw["drug"] = drug
    res = simulate_response(micro, SimConfig(**kw))
    ref = reference[label]

    assert float(res.kill_fraction) == _hex(ref["kill_fraction"])
    assert float(res.extras["mean_auc"]) == _hex(ref["mean_auc"])
    assert float(res.extras["plasma_mean"]) == _hex(ref["plasma_mean"])
    assert float(res.extras["penetration_index"]) == _hex(ref["penetration_index"])
    assert float(np.asarray(res.transport.auc).sum()) == _hex(ref["auc_sum"])
    assert (float(np.asarray(res.transport.final_concentration).sum())
            == _hex(ref["final_c_sum"]))
    assert [float(v) for v in res.viable_history] == [
        _hex(v) for v in ref["viable_history"]
    ]


def test_static_solver_matches_stored_reference_bitwise(reference,
                                                        flat_micro_factory):
    """Same pin one level down, on the raw solver."""
    micro = flat_micro_factory(ny=40, nx=40, dx=15.0, ecm_value=0.3)
    drug = DrugParams(name="tmz", diffusivity=300.0, clearance=3e-4, binding=1e-4)
    pde = PDEConfig(dt=10.0, total_time=1800.0, convection=False, save_every=30)
    source = micro.domain.zeros()
    source[::12, ::12] = 1.0

    res = simulate_drug_transport(micro, drug, source, pde)
    ref = reference["solver_flat"]
    assert float(res.auc.sum()) == _hex(ref["auc_sum"])
    assert float(res.auc.max()) == _hex(ref["auc_max"])
    assert float(res.final_concentration.sum()) == _hex(ref["final_c_sum"])
    assert float(res.penetration_index) == _hex(ref["penetration_index"])
    assert [float(v) for v in res.mean_concentration_history] == [
        _hex(v) for v in ref["mean_hist"]
    ]


def test_real_tissue_static_path_matches_stored_reference_bitwise(reference,
                                                                  real_micro,
                                                                  real_drug):
    """And on the REAL Xenium section, through the whole engine."""
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=3,
                  dose_interval=3600.0 * 4)
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 12, convection=True,
                    save_every=100)
    res = simulate_response(real_micro, SimConfig(drug=real_drug, pk=pk, pde=pde,
                                                  seed=0))
    ref = reference["real_xenium_dx80_tmz_3dose"]
    assert float(res.kill_fraction) == _hex(ref["kill_fraction"])
    assert float(res.extras["mean_auc"]) == _hex(ref["mean_auc"])
    assert float(np.asarray(res.transport.auc).sum()) == _hex(ref["auc_sum"])
    assert (float(np.asarray(res.transport.final_concentration).sum())
            == _hex(ref["final_c_sum"]))


# --------------------------------------------------------------------------- #
# (3b) The PRODUCTION grid — added 2026-08-12, see the module docstring
# --------------------------------------------------------------------------- #
#: The exact production configuration: the same reconstruction call the landing
#: asset and every shipped result use. ``dt`` must stay at 10 s — at dx = 12 um
#: the advection CFL limit on this section is ~18.9 s, so dt = 20 raises.
_DX12_PDE = dict(dt=10.0, total_time=1800.0, convection=True, save_every=30)


@pytest.fixture(scope="module")
def production_micro_dx12():
    """The real Xenium GBM section on the PRODUCTION grid (dx = 12 um).

    Identical arguments to ``experiments/export_landing_asset.py``, so this pins
    the configuration the shipped IFP figures are actually computed on. Costs a
    couple of seconds to assemble; module-scoped so both dx = 12 pins share it.
    """
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    return reconstruct_real(source="xenium", sample="data_cache/xenium_gbm",
                            dx_um=12.0, oxygen_mode="pde", state_mode="oxygen")


def test_production_dx12_ifp_solve_matches_stored_reference_bitwise(
        reference, production_micro_dx12):
    """The IFP solve itself, on the production grid, pinned bit for bit.

    This is the level the tumour-sigma units fix acts at, so it is the level the
    production guard belongs at: if ``_TUMOR_SIGMA_UM`` is ever changed, or the
    um -> pixel conversion is dropped and the constant goes back to meaning grid
    cells, these three numbers move immediately — no transport solve required to
    see it.
    """
    from tmesim.physics.pressure import solve_ifp

    p, velocity = solve_ifp(production_micro_dx12)
    mask = np.asarray(production_micro_dx12.domain.tissue_mask, dtype=bool)
    ref = reference["production_xenium_dx12_tmz_3dose"]

    assert float(p[mask].max()) == _hex(ref["ifp_max_in_tissue_mmhg"])
    assert float(p[mask].mean()) == _hex(ref["ifp_mean_in_tissue_mmhg"])
    assert (float(np.hypot(velocity[0], velocity[1])[mask].max())
            == _hex(ref["ifp_speed_max_um_s"]))


def test_production_dx12_static_path_matches_stored_reference_bitwise(
        reference, production_micro_dx12):
    """End-to-end on the PRODUCTION grid, convection ON.

    The reference file used to pin only dx = 15 and dx = 80 — neither of which
    anything ships on. This closes that hole: the grid every published result in
    ``results/`` is computed on now has a bitwise pin too, running the same
    Darcy-convection path that the tumour-sigma fix touches.
    """
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.validation.calibration import calibrated_drugparams

    drug = apply_barrier_to_drugparams(calibrated_drugparams("temozolomide"),
                                       library.get_molecule("temozolomide"))
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=3,
                  dose_interval=3600.0 * 4)
    res = simulate_response(
        production_micro_dx12,
        SimConfig(drug=drug, pk=pk, pde=PDEConfig(**_DX12_PDE), seed=0),
    )
    ref = reference["production_xenium_dx12_tmz_3dose"]

    assert float(res.kill_fraction) == _hex(ref["kill_fraction"])
    assert float(res.extras["mean_auc"]) == _hex(ref["mean_auc"])
    assert float(res.extras["plasma_mean"]) == _hex(ref["plasma_mean"])
    assert float(res.extras["penetration_index"]) == _hex(ref["penetration_index"])
    assert float(np.asarray(res.transport.auc).sum()) == _hex(ref["auc_sum"])
    assert (float(np.asarray(res.transport.final_concentration).sum())
            == _hex(ref["final_c_sum"]))
    assert [float(v) for v in res.viable_history] == [
        _hex(v) for v in ref["viable_history"]
    ]


def test_tumor_occupancy_smoothing_is_a_physical_length():
    """The units fix, asserted directly so it cannot silently regress.

    ``_TUMOR_SIGMA_UM`` is a MODEL parameter (it sets the vessel-leak gain and
    the lymphatic suppression), so its value must be a physical length and the
    pixel sigma handed to ``gaussian_filter`` must scale as ``1/dx``. Asserting
    the conversion rather than the constant means a deliberate recalibration of
    the 60 um is still allowed; reverting to grid-cell units is not.

    Also pins the fixed point that makes the whole change safe: at the
    production grid the conversion is EXACT in IEEE-754, so ``gaussian_filter``
    receives a bit-identical argument and no shipped result moves.
    """
    from tmesim.physics import pressure as _pressure

    assert not hasattr(_pressure, "_TUMOR_SIGMA"), (
        "the grid-cell tumour-sigma constant is back; the smoothing length is a "
        "physical length, not a pixel count"
    )
    sigma_um = _pressure._TUMOR_SIGMA_UM
    assert sigma_um > 0.0
    # dx = 12 um is the exact fixed point -> production is bit-identical.
    assert sigma_um / 12.0 == 5.0
    assert float.hex(sigma_um / 12.0) == float.hex(5.0)
    # The physical smoothing length is the SAME on every mesh.
    for dx in (8.0, 12.0, 15.0, 18.0, 24.0, 80.0):
        assert (sigma_um / dx) * dx == pytest.approx(sigma_um, rel=1e-12)


def test_unit_source_profile_is_bit_identical_to_none(flat_micro_factory):
    """``source_profile = ones`` must reproduce the static path EXACTLY (not just
    to a tolerance): ``dt * 1.0 == dt`` in IEEE-754, so the arithmetic is the
    same arithmetic. This is what makes the new branch provably inert."""
    micro = flat_micro_factory(ny=40, nx=40, dx=15.0, ecm_value=0.3)
    drug = DrugParams(name="tmz", diffusivity=300.0, clearance=3e-4, binding=1e-4)
    pde = PDEConfig(dt=10.0, total_time=1800.0, convection=False, save_every=30)
    source = micro.domain.zeros()
    source[::12, ::12] = 1.0
    n_steps = int(round(pde.total_time / pde.dt))

    a = simulate_drug_transport(micro, drug, source, pde)
    b = simulate_drug_transport(micro, drug, source, pde,
                                source_profile=np.ones(n_steps))
    assert np.array_equal(a.auc, b.auc)
    assert np.array_equal(a.final_concentration, b.final_concentration)
    assert np.array_equal(a.mean_concentration_history,
                          b.mean_concentration_history)
    assert a.penetration_index == b.penetration_index


def test_source_profile_validation(flat_micro_factory):
    micro = flat_micro_factory(ny=16, nx=16, dx=15.0, ecm_value=0.3)
    drug = DrugParams(name="tmz", diffusivity=300.0, clearance=3e-4)
    pde = PDEConfig(dt=10.0, total_time=200.0, convection=False, save_every=5)
    source = micro.domain.ones()
    n_steps = int(round(pde.total_time / pde.dt))

    with pytest.raises(ValueError, match="at least one scale factor per step"):
        simulate_drug_transport(micro, drug, source, pde,
                                source_profile=np.ones(n_steps - 1))
    with pytest.raises(ValueError, match="1-D"):
        simulate_drug_transport(micro, drug, source, pde,
                                source_profile=np.ones((n_steps, 2)))
    with pytest.raises(ValueError, match="non-finite"):
        bad = np.ones(n_steps)
        bad[3] = np.nan
        simulate_drug_transport(micro, drug, source, pde, source_profile=bad)
    with pytest.raises(ValueError, match="non-negative"):
        neg = np.ones(n_steps)
        neg[3] = -0.5
        simulate_drug_transport(micro, drug, source, pde, source_profile=neg)
    # A LONGER profile is fine — only the first n_steps entries are consumed.
    long_p = np.concatenate([np.ones(n_steps), np.full(50, 1e6)])
    ok = simulate_drug_transport(micro, drug, source, pde, source_profile=long_p)
    base = simulate_drug_transport(micro, drug, source, pde)
    assert np.array_equal(ok.auc, base.auc)


# --------------------------------------------------------------------------- #
# (4) Convergence: constant plasma -> the static-source answer
# --------------------------------------------------------------------------- #
def test_time_resolved_converges_to_static_for_constant_plasma():
    """PHYSICS CHECK (b). Hold the plasma level constant and the time-resolved
    path must reproduce the static-source result. Verified two ways: with a
    literally constant profile (exact to round-off) and with a plasma half-life
    made long compared with the window, which is the physical limit the schedule
    sweep actually approaches."""
    micro = _slab_micro()
    drug = DrugParams(name="tl", diffusivity=100.0, clearance=1e-3, binding=0.0)
    dt, total_time = 5.0, 3600.0 * 6
    pde = PDEConfig(dt=dt, total_time=total_time, convection=False,
                    save_every=200)
    n_steps = int(round(total_time / dt))
    unit = delivery_source(micro, drug, 1.0)
    level = 10.0

    static = simulate_drug_transport(micro, drug, unit * level, pde)

    # (i) literally constant profile
    flat = simulate_drug_transport(micro, drug, unit, pde,
                                   source_profile=np.full(n_steps, level))
    assert np.allclose(flat.auc, static.auc, rtol=1e-12, atol=0.0)

    # (ii) the physical limit: half-life >> window, peak rescaled to matched
    #      cumulative exposure. Convergence must be MONOTONE in half-life.
    errs = []
    for half_life_h in (4.0, 48.0, 480.0, 4800.0):
        hl = half_life_h * 3600.0
        pk = PKConfig(plasma_peak=_matched_peak(level, hl, total_time, dt),
                      half_life=hl, n_doses=1, dose_interval=total_time)
        prof = plasma_profile(pk, np.arange(n_steps, dtype=float) * dt)
        assert abs(float(prof.mean()) - level) < 1e-9 * level  # matched exposure
        res = simulate_drug_transport(micro, drug, unit, pde, source_profile=prof)
        errs.append(abs(float(res.auc.sum()) / float(static.auc.sum()) - 1.0))

    assert errs == sorted(errs, reverse=True), (
        f"convergence to the static limit is not monotone in half-life: {errs}"
    )
    assert errs[-1] < 1e-4, (
        f"a 4800 h half-life over a 6 h window is effectively a constant plasma "
        f"level but still differs from the static source by {errs[-1]:.2e}"
    )


# --------------------------------------------------------------------------- #
# (5) Mass balance with a time-varying source
# --------------------------------------------------------------------------- #
def test_mass_balance_with_time_varying_source(inert_drug, flat_micro_factory):
    """PHYSICS CHECK (c). With clearance = binding = 0 and no-flux boundaries,
    nothing escapes, so the final tissue mass must equal exactly the mass
    injected by the (now time-varying) source: ``dt * sum_n p[n] * sum(S)``."""
    micro = flat_micro_factory(ny=40, nx=40, dx=15.0, ecm_value=0.3)
    pde = PDEConfig(dt=10.0, total_time=3600.0, convection=False, save_every=20)
    source = micro.domain.zeros()
    source[20, 20] = 2.0
    n_steps = int(round(pde.total_time / pde.dt))

    pk = PKConfig(plasma_peak=50.0, half_life=600.0, n_doses=4,
                  dose_interval=600.0)
    prof = plasma_profile(pk, np.arange(n_steps, dtype=float) * pde.dt)
    assert prof.max() > 1.5 * prof.min(), "profile must actually vary"

    res = simulate_drug_transport(micro, inert_drug, source, pde,
                                  source_profile=prof)
    injected = pde.dt * float(prof.sum()) * float(source.sum())
    final_mass = float(res.final_concentration.sum())

    assert abs(final_mass - injected) / injected < 1e-9
    assert res.final_concentration.min() >= -1e-9
    assert res.auc.min() >= -1e-9


def test_matched_cumulative_exposure_injects_identical_total(inert_drug,
                                                             flat_micro_factory):
    """The identity the whole "matched cumulative exposure" comparison rests on:
    a time-resolved schedule whose window-mean plasma equals ``plasma_mean``
    delivers exactly as much drug as the static source at ``plasma_mean``. If
    this ever drifted, every schedule comparison would be confounded by dose."""
    micro = flat_micro_factory(ny=30, nx=30, dx=15.0, ecm_value=0.2)
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 2, convection=False,
                    save_every=50)
    source = micro.domain.zeros()
    source[15, 15] = 1.0
    n_steps = int(round(pde.total_time / pde.dt))
    times = np.arange(n_steps, dtype=float) * pde.dt

    level = 7.5
    static = simulate_drug_transport(micro, inert_drug, source * level, pde)
    static_mass = float(static.final_concentration.sum())

    for n_doses, interval, half_life in ((1, pde.total_time, 3600.0),
                                         (3, 2400.0, 1800.0),
                                         (10, 720.0, 600.0)):
        peak = level / plasma_profile(
            PKConfig(plasma_peak=1.0, half_life=half_life, n_doses=n_doses,
                     dose_interval=interval), times).mean()
        prof = plasma_profile(
            PKConfig(plasma_peak=peak, half_life=half_life, n_doses=n_doses,
                     dose_interval=interval), times)
        res = simulate_drug_transport(micro, inert_drug, source, pde,
                                      source_profile=prof)
        mass = float(res.final_concentration.sum())
        assert abs(mass - static_mass) / static_mass < 1e-9, (
            f"n_doses={n_doses} delivered {mass:.6e} vs static {static_mass:.6e}"
        )


# --------------------------------------------------------------------------- #
# (6) Depth: the schedule changes the exposure field, in a derivable direction
# --------------------------------------------------------------------------- #
def test_earlier_delivery_reaches_deeper():
    """PHYSICS CHECK (a), with the direction ASSERTED.

    Naive expectation: at matched total exposure a long slow infusion reaches
    DEEPER than a short bolus. For this linear transport model that is provably
    the wrong way round (see the module docstring): the entire schedule effect is
    the terminal residual ``C(T)``, and the earlier the drug is given the smaller
    that residual is. So the BOLUS wins — in absolute exposure everywhere, and
    in relative depth too, because the un-counted term is a doubly-smoothed
    (flatter) field whose removal steepens what is left.

    This test asserts that direction on a transport-limited slab, and asserts the
    accompanying bound: the effect must SHRINK as the observation window grows
    past washout, which is what identifies it as a residual effect rather than a
    numerical artefact.
    """
    micro = _slab_micro()
    dx = 10.0
    # Transport-limited: L = sqrt(D/k) = 316 um against a 1600 um slab, and a
    # loss time-scale 1/k = 1000 s comparable to the short window below.
    drug = DrugParams(name="tl", diffusivity=100.0, clearance=1e-3, binding=0.0)
    dt = 5.0
    unit = delivery_source(micro, drug, 1.0)
    level = 10.0
    mid = micro.domain.shape[0] // 2

    def profile_result(half_life_h, total_time):
        n = int(round(total_time / dt))
        hl = half_life_h * 3600.0
        pk = PKConfig(plasma_peak=_matched_peak(level, hl, total_time, dt),
                      half_life=hl, n_doses=1, dose_interval=total_time)
        prof = plasma_profile(pk, np.arange(n, dtype=float) * dt)
        assert abs(float(prof.mean()) - level) < 1e-9 * level
        pde = PDEConfig(dt=dt, total_time=total_time, convection=False,
                        save_every=200)
        return simulate_drug_transport(micro, drug, unit, pde,
                                       source_profile=prof)

    # --- short window (2 h): comparable to washout, so schedule bites ---------
    short_T = 3600.0 * 2
    bolus = profile_result(0.05, short_T)      # 3 min half-life  -> all up front
    sustained = profile_result(48.0, short_T)  # 48 h half-life   -> ~flat infusion

    d_bolus = _relative_depth_um(bolus.auc[mid], dx)
    d_sustained = _relative_depth_um(sustained.auc[mid], dx)

    # 1. The schedule changes the answer AT ALL (this is the whole point).
    assert not np.allclose(bolus.auc, sustained.auc), (
        "matched-exposure schedules still produce an identical exposure field — "
        "the time-resolved source is not actually driving the solve"
    )
    # 2. Direction, on total exposure.
    assert bolus.auc.sum() > sustained.auc.sum() * 1.01
    # 3. Direction, on SHAPE — the depth at which exposure falls to 10% of its
    #    perivascular value, which is invariant to amplitude.
    assert d_bolus > d_sustained + 1.0, (
        f"bolus depth {d_bolus:.1f} um is not deeper than sustained "
        f"{d_sustained:.1f} um"
    )
    # 4. Same conclusion from the solver's own shape statistic.
    assert bolus.penetration_index > sustained.penetration_index
    # 5. The mechanism: the sustained infusion leaves drug in flight at T.
    assert (sustained.final_concentration.sum()
            > 10.0 * bolus.final_concentration.sum())

    # --- long window (96 h): past washout, so the effect must collapse --------
    long_T = 3600.0 * 96
    bolus_l = profile_result(0.05, long_T)
    sustained_l = profile_result(48.0, long_T)
    gap_short = bolus.auc.sum() / sustained.auc.sum() - 1.0
    gap_long = bolus_l.auc.sum() / sustained_l.auc.sum() - 1.0
    assert gap_long > 0.0
    assert gap_long < 0.1 * gap_short, (
        f"the schedule effect ({gap_long:.4f}) did not collapse relative to the "
        f"short window ({gap_short:.4f}) — it is not a terminal-residual effect"
    )


# --------------------------------------------------------------------------- #
# (7) Capability: the previously-zero spread is no longer zero
# --------------------------------------------------------------------------- #
def test_fractionation_at_matched_exposure_is_no_longer_degenerate(real_micro,
                                                                   real_drug):
    """The regression this whole change exists to fix, on the REAL section.

    With the static source, 1 / 3 / 10 fractions at matched cumulative plasma
    exposure were bit-identical inputs and the kill spread was exactly
    ``0.00e+00``. Assert (i) that the static path STILL behaves that way — the
    old behaviour is preserved bit-for-bit when the flag is off — and (ii) that
    the time-resolved path does not, and moves in the derivable direction: more
    fractions means later average delivery means a bigger terminal residual means
    less exposure realised inside the window.

    Run on ``real_micro`` rather than the 300-cell hand-built micro because kill
    fraction is quantised at 1/n_tumor and the schedule effect here is ~1%.
    """
    total_time, dt = 3600.0 * 12, 20.0
    half_life = 3600.0 * 4
    times = np.arange(int(round(total_time / dt)), dtype=float) * dt
    target_mean = plasma_profile(
        PKConfig(plasma_peak=300.0, half_life=half_life, n_doses=1,
                 dose_interval=total_time), times).mean()

    schedules = [(1, total_time), (3, 3600.0 * 4), (10, 3600.0 * 1.2)]
    for flag in (False, True):
        pde = PDEConfig(dt=dt, total_time=total_time, convection=True,
                        save_every=200, time_resolved=flag)
        kills, aucs = [], []
        for n, interval in schedules:
            unit = plasma_profile(
                PKConfig(plasma_peak=1.0, half_life=half_life, n_doses=n,
                         dose_interval=interval), times).mean()
            pk = PKConfig(plasma_peak=float(target_mean / unit),
                          half_life=half_life, n_doses=n, dose_interval=interval)
            res = simulate_response(real_micro,
                                    SimConfig(drug=real_drug, pk=pk, pde=pde,
                                              seed=0))
            kills.append(float(res.kill_fraction))
            aucs.append(float(np.asarray(res.transport.auc).sum()))

        if not flag:
            # Static source: matched cumulative exposure -> identical inputs.
            assert max(aucs) - min(aucs) == 0.0
            assert max(kills) - min(kills) == 0.0
        else:
            # Time-resolved: fractionation is now visible, and monotone in n.
            assert aucs == sorted(aucs, reverse=True), aucs
            assert max(aucs) - min(aucs) > 0.0
            assert kills == sorted(kills, reverse=True), kills
            assert max(kills) - min(kills) > 0.0, (
                "fractionation still produces a zero kill spread with "
                "time_resolved=True"
            )


def test_time_resolved_engine_reports_matched_total_exposure():
    """The engine must deliver the same total drug on both paths, so that a
    static-vs-resolved comparison is a comparison of SCHEDULE and not of dose."""
    micro = _build_selftest_micro()
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0 * 2, n_doses=4,
                  dose_interval=3600.0 * 2)
    kw = dict(drug=WEAK_DRUG, pk=pk, seed=0)
    base = PDEConfig(dt=20.0, total_time=3600.0 * 8, convection=True,
                     save_every=100)
    res_static = simulate_response(micro, SimConfig(pde=base, **kw))
    res_tr = simulate_response(
        micro,
        SimConfig(pde=PDEConfig(dt=base.dt, total_time=base.total_time,
                                convection=True, save_every=base.save_every,
                                time_resolved=True), **kw),
    )
    assert res_static.extras["time_resolved"] == 0.0
    assert res_tr.extras["time_resolved"] == 1.0
    # Same window-mean plasma reported on both paths -> same total delivered.
    assert res_tr.extras["plasma_delivered_mean"] == pytest.approx(
        res_static.extras["plasma_mean"], rel=1e-12
    )
    # ... but a genuinely different exposure field.
    assert res_tr.extras["mean_auc"] != res_static.extras["mean_auc"]


def test_time_resolved_is_deterministic():
    micro = _build_selftest_micro()
    pde = PDEConfig(dt=20.0, total_time=3600.0 * 4, convection=True,
                    save_every=100, time_resolved=True)
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0, n_doses=3,
                  dose_interval=3600.0)
    cfg = SimConfig(drug=WEAK_DRUG, pk=pk, pde=pde, seed=0)
    a = simulate_response(micro, cfg)
    b = simulate_response(micro, cfg)
    assert a.kill_fraction == b.kill_fraction
    assert np.array_equal(a.transport.auc, b.transport.auc)
