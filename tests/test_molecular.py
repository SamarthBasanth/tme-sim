"""Fast focused tests for the molecular layer (tmesim.molecular).

Covers the whole pipeline: descriptors -> Stokes-Einstein free diffusivity ->
Brownian-dynamics effective diffusivity -> permeability / BBB -> binding ->
the compile bridge into a frozen ``core.DrugParams``, INCLUDING an end-to-end
run of a compiled drug through the tissue engine.

The suite is deliberately fast: the Brownian sims use small particle/step
counts (they are the only expensive piece), fixed seeds, and short PDE windows.
Everything here is deterministic, so a value that passes once passes forever.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from tmesim.molecular import (
    binding,
    compile_drug,
    descriptors,
    library,
    permeability,
    stokes_einstein,
)
from tmesim.molecular.brownian import simulate_effective_diffusivity
from tmesim.molecular.types import Molecule
from tmesim.core import DrugParams, PDEConfig, SimConfig, SimulationResult
from tmesim.sim.engine import simulate_response

# Small Brownian settings: keep every simulation to a few hundredths of a second.
# Deterministic given the seed, so these values are stable across runs.
_FAST_BROWNIAN = {"n_particles": 150, "steps": 600, "seed": 0}


# --------------------------------------------------------------------------- #
# descriptors
# --------------------------------------------------------------------------- #
def test_temozolomide_mw_from_library():
    tmz = library.get_molecule("temozolomide")
    assert tmz.get("mw") == pytest.approx(194.0, abs=2.0)


def test_compute_descriptors_bad_smiles_raises():
    with pytest.raises(ValueError):
        descriptors.compute_descriptors("Xy!!zzz-not-a-smiles")


def test_from_smiles_never_crashes_and_fills_radius():
    # Good SMILES: radius must be filled.
    good = descriptors.from_smiles("temozolomide", "CN1C(=O)N2C=NC(=C2N=N1)C(=O)N")
    assert good.get("radius_hydrodynamic_nm") is not None
    assert good.get("radius_hydrodynamic_nm") > 0.0

    # Bad SMILES + unknown name: must degrade gracefully (no crash) and still
    # populate a radius rather than raising.
    junk = descriptors.from_smiles("totally_unknown_drug", "%%%not-a-smiles%%%")
    assert isinstance(junk, Molecule)
    assert "radius_hydrodynamic_nm" in junk.descriptors


def test_antibody_entry_flagged():
    bev = library.get_molecule("bevacizumab")
    assert bool(bev.get("is_antibody")) is True
    assert bev.get("mw") > 100_000.0


# --------------------------------------------------------------------------- #
# stokes_einstein
# --------------------------------------------------------------------------- #
def test_free_diffusivity_scales_inverse_radius():
    radii = [0.4, 0.8, 1.6, 3.2]
    ds = [stokes_einstein.free_diffusivity(r) for r in radii]
    # strictly decreasing in radius
    assert all(ds[i] > ds[i + 1] for i in range(len(ds) - 1))
    # D * r constant  =>  D ~ 1/r
    products = [d * r for d, r in zip(ds, radii)]
    assert max(products) - min(products) < 1e-6 * max(products)


def test_free_diffusivity_scales_inverse_viscosity():
    etas = [1.0e-3, 2.0e-3, 4.0e-3]
    ds = [stokes_einstein.free_diffusivity(1.0, viscosity_pa_s=e) for e in etas]
    assert all(ds[i] > ds[i + 1] for i in range(len(ds) - 1))
    products = [d * e for d, e in zip(ds, etas)]
    assert max(products) - min(products) < 1e-6 * max(products)


def test_small_molecule_D_much_greater_than_antibody():
    r_small = stokes_einstein.hydrodynamic_radius_nm(300.0)
    r_ab = stokes_einstein.hydrodynamic_radius_nm(150000.0, is_antibody=True)
    d_small = stokes_einstein.free_diffusivity(r_small)
    d_ab = stokes_einstein.free_diffusivity(r_ab)
    assert d_small > 8.0 * d_ab


def test_small_molecule_D_in_plausible_hundreds():
    r = stokes_einstein.hydrodynamic_radius_nm(300.0)
    d = stokes_einstein.free_diffusivity(r)
    # ~300 g/mol solute in interstitial fluid: hundreds of um^2/s.
    assert 100.0 <= d <= 600.0


# --------------------------------------------------------------------------- #
# brownian
# --------------------------------------------------------------------------- #
def test_zero_obstacle_recovers_free_D():
    free_D = 300.0
    res = simulate_effective_diffusivity(
        free_D, obstacle_fraction=0.0, n_particles=200, steps=800, seed=0)
    # With no obstacles the effective D must recover free_D within ~10%.
    assert abs(res["effective_D"] - free_D) / free_D < 0.12


def test_hindrance_decreases_with_obstacle_fraction():
    free_D = 300.0
    lo = simulate_effective_diffusivity(
        free_D, obstacle_fraction=0.0, n_particles=200, steps=800, seed=0)
    hi = simulate_effective_diffusivity(
        free_D, obstacle_fraction=0.4, n_particles=200, steps=800, seed=0)
    assert hi["hindrance"] < lo["hindrance"]


def test_msd_derived_effective_D_positive():
    res = simulate_effective_diffusivity(
        300.0, obstacle_fraction=0.3, n_particles=150, steps=600, seed=0)
    assert res["effective_D"] > 0.0
    msd = np.asarray(res["msd"])
    assert msd.size > 0
    assert np.all(np.isfinite(msd))
    assert msd[-1] > 0.0


# --------------------------------------------------------------------------- #
# permeability
# --------------------------------------------------------------------------- #
def test_bbb_antibody_near_zero_and_below_small_lipophilic():
    bev = library.get_molecule("bevacizumab")
    lomustine = library.get_molecule("lomustine")  # small, lipophilic, CNS-like
    bbb_ab = permeability.bbb_penetrance(bev.descriptors)
    bbb_sm = permeability.bbb_penetrance(lomustine.descriptors)
    assert bbb_ab == pytest.approx(0.0, abs=1e-9)
    assert bbb_sm > bbb_ab
    assert bbb_sm > 0.1  # a genuinely lipophilic small molecule crosses well


def test_membrane_permeability_positive_and_finite():
    lomustine = library.get_molecule("lomustine")
    p = permeability.membrane_permeability(lomustine.descriptors)
    assert math.isfinite(p)
    assert p > 0.0


# --------------------------------------------------------------------------- #
# binding
# --------------------------------------------------------------------------- #
def test_binding_outputs_in_documented_ranges():
    for name in ("temozolomide", "doxorubicin", "lomustine"):
        d = library.get_molecule(name).descriptors
        ff = binding.plasma_free_fraction(d)
        bind = binding.tissue_binding_rate(d)
        clr = binding.clearance_rate(d)
        assert 0.01 <= ff <= 1.0
        assert 1e-5 <= bind <= 5e-3
        assert 1e-5 <= clr <= 1e-3


def test_lipophilic_lower_free_fraction_than_hydrophilic():
    lipophilic = {"mw": 540.0, "logp": 4.5, "hbd": 4}
    hydrophilic = {"mw": 150.0, "logp": -1.0, "hbd": 1}
    ff_lip = binding.plasma_free_fraction(lipophilic)
    ff_hyd = binding.plasma_free_fraction(hydrophilic)
    assert ff_lip < ff_hyd


# --------------------------------------------------------------------------- #
# compile_drug (the bridge)
# --------------------------------------------------------------------------- #
def _finite(x) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(x)


def test_molecule_to_drugparams_valid_and_sane():
    tmz = library.get_molecule("temozolomide")
    dp = compile_drug.molecule_to_drugparams(tmz, brownian_kw=_FAST_BROWNIAN)
    assert isinstance(dp, DrugParams)
    for fld in ("diffusivity", "molecular_weight", "ic50", "hill",
                "binding", "clearance", "bbb_penetrance"):
        assert _finite(getattr(dp, fld)), fld
    assert 0.1 <= dp.diffusivity <= 1000.0
    assert 50.0 <= dp.molecular_weight <= 2000.0
    assert dp.ic50 > 0.0
    assert 0.0 <= dp.bbb_penetrance <= 1.0
    assert 1e-5 <= dp.binding <= 5e-3
    assert 1e-5 <= dp.clearance <= 1e-3


def test_doxorubicin_lower_diffusivity_and_bbb_than_temozolomide():
    tmz = compile_drug.molecule_to_drugparams(
        library.get_molecule("temozolomide"), brownian_kw=_FAST_BROWNIAN)
    dox = compile_drug.molecule_to_drugparams(
        library.get_molecule("doxorubicin"), brownian_kw=_FAST_BROWNIAN)
    # Doxorubicin is larger => hindered more => lower effective diffusivity.
    assert dox.diffusivity < tmz.diffusivity
    # Doxorubicin (large, polar, many donors) crosses no better than TMZ.
    assert dox.bbb_penetrance <= tmz.bbb_penetrance


def test_compiled_drug_runs_through_tissue_engine(real_micro):
    # Compile a real molecule into DrugParams and drive the FROZEN engine with it
    # on the REAL reconstructed microenvironment.
    dp = compile_drug.molecule_to_drugparams(
        library.get_molecule("temozolomide"), brownian_kw=_FAST_BROWNIAN)

    cfg = SimConfig(seed=0)
    cfg.drug = dp
    cfg.pde = PDEConfig(dt=20.0, total_time=800.0, convection=True, save_every=10)

    result = simulate_response(real_micro, cfg)

    assert isinstance(result, SimulationResult)
    assert 0.0 <= float(result.predicted_response) <= 1.0
    assert 0.0 <= float(result.viable_tumor_fraction) <= 1.0
    # The engine actually saw our compiled BBB value.
    assert float(result.extras["bbb_penetrance"]) == pytest.approx(
        dp.bbb_penetrance, abs=1e-9)


# --------------------------------------------------------------------------- #
# determinism
# --------------------------------------------------------------------------- #
def test_compile_is_deterministic_with_fixed_seed():
    kw = {"n_particles": 150, "steps": 600, "seed": 42}
    dp_a = compile_drug.molecule_to_drugparams(
        library.get_molecule("temozolomide"), brownian_kw=kw)
    dp_b = compile_drug.molecule_to_drugparams(
        library.get_molecule("temozolomide"), brownian_kw=kw)
    assert dp_a.diffusivity == dp_b.diffusivity


#: Sentinel for "this descriptor was never computed", as distinct from a
#: descriptor whose value happens to be zero. Conflating those two is the
#: exact defect the class below guards.
class _Absent:
    def __repr__(self):  # pragma: no cover - debugging aid only
        return "<absent>"


_ABSENT = _Absent()


class _FakeMolecule:
    """The minimum surface ``build_profile`` touches: a ``descriptors`` dict.

    A real ``Molecule`` needs RDKit and a SMILES parse; this test is about one
    branch reached from ``descriptors``, and building a real molecule would
    make the test depend on a chemistry toolkit to check a dictionary lookup.
    """

    def __init__(self, descriptors):
        self.descriptors = dict(descriptors)
        self.name = "probe"
        self.smiles = "CCO"

    def get(self, key, default=None):
        """``build_profile`` reads mw / is_antibody / radius through this as
        well as through ``.descriptors``; a stand-in that answers only one of
        the two fails on an AttributeError that says nothing about the branch
        under test."""
        return self.descriptors.get(key, default)


class TestAnAbsentDescriptorIsRefusedNotDefaulted:
    """`compile_drug.build_profile` defaulted a missing logP/HBD to zero.

    The two-route law's membrane term is ``log K = logP - 1.34*HBD``, so
    ``(0, 0)`` is ``log K = 0.00`` — the **85th percentile** of membrane
    openness across the scored benchmark. Defaulting silently asserted that an
    unknown molecule crosses membranes more easily than four out of five known
    ones.

    MEASURED 2026-09-20 at MW 263 Da: 45.68 → 189.34 um^2/s, a **4.14x
    inflation** — carried to the user with the law's full "scored against the
    benchmark, r^2 = 0.718" provenance stapled to it, because the provenance
    string does not know the descriptors were absent.

    The benchmark experiment already excludes no-logP compounds as
    UNSCOREABLE. The shipping path compiled them confidently. Both positions
    cannot be right, and the experiment's is the correct one.
    """

    def _mol(self, **over):
        """A descriptor dict shaped like the real one, minus what we drop."""
        d = {"mw": 263.0, "logp": 1.2, "hbd": 1}
        for k, v in over.items():
            if v is _ABSENT:
                d.pop(k, None)
            else:
                d[k] = v
        return d

    @pytest.mark.parametrize("missing", ["logp", "hbd"])
    def test_a_missing_descriptor_raises_rather_than_compiling(self, missing):
        from tmesim.molecular import compile_drug as cd

        mol = _FakeMolecule(self._mol(**{missing: _ABSENT}))
        with pytest.raises(ValueError) as exc:
            cd.build_profile(mol, two_route=True)
        msg = str(exc.value)
        assert missing in msg
        # The refusal must say what to DO, not merely that it refused.
        assert "two_route=False" in msg, \
            "a refusal that names no alternative is a dead end"

    def test_a_genuine_zero_is_not_treated_as_absent(self):
        """`or 0.0` swallowed a real 0.0 into the same branch as a missing
        key, so a compound whose logP genuinely IS zero was indistinguishable
        from one whose logP was never computed. Those are different facts."""
        from tmesim.molecular import compile_drug as cd

        prof = cd.build_profile(_FakeMolecule(self._mol(logp=0.0, hbd=0)),
                                two_route=True)
        assert float(prof.effective_diffusivity) > 0.0

    def test_the_paracellular_path_is_still_open_without_descriptors(self):
        """The refusal is scoped to the law that NEEDS the descriptors. A
        caller who wants the paracellular-only answer may still have it — and
        gets it labelled as that, rather than as the two-route result."""
        from tmesim.molecular import compile_drug as cd

        prof = cd.build_profile(_FakeMolecule(self._mol(logp=_ABSENT)),
                                two_route=False)
        assert float(prof.effective_diffusivity) > 0.0
