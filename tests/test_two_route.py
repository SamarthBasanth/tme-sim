"""Tests for the two-route (paracellular + transcellular) diffusivity model.

Three things are pinned here, in descending order of importance:

1. **It is default-off and additive.** `DrugParams.diffusivity` must be
   bit-identical with the two-route module present, and `calibrate_alpha()` must
   still return the published 148.60815765445545. A stale agent broke exactly
   this by changing a molecular input that flows into doxorubicin's binding;
   alpha absorbed it and drifted, breaking two bitwise PK pins.
2. **The limit reduces correctly.** At zero membrane permeability the two-route
   law must reproduce `ecm_mesh.cell_exclusion` to machine precision -- that is
   the real check on the algebra, and it is what makes the new law a strict
   generalisation of the shipped one rather than a competitor to it.
3. **It produces spread at fixed molecular size**, which is the property the
   whole exercise exists to add and which the impermeable-cell law cannot have.
"""
from __future__ import annotations

import math

import pytest

from tmesim.molecular import ecm_mesh, two_route
from tmesim.molecular.two_route import (
    DEFAULT_CELL_RADIUS_UM, DEFAULT_CELL_VOLUME_FRACTION,
    cytoplasm_diffusivity_um2s, hasselman_johnson_hindrance,
    membrane_permeability_um_s, two_route_ceiling, two_route_diffusivity,
    two_route_factor,
)

PHI = DEFAULT_CELL_VOLUME_FRACTION
SIGMA_M = 265.0   # um^2/s, an ECM-hindered interstitial D for a small molecule
D_CYTO = 200.0    # um^2/s


# --------------------------------------------------------------------------- #
# 1. DEFAULT-OFF: nothing the engine ships may move.
# --------------------------------------------------------------------------- #
class TestDefaultOff:
    def test_alpha_is_unchanged(self):
        """The one number that must never drift."""
        from tmesim.validation.calibration import calibrate_alpha
        assert calibrate_alpha() == 148.60815765445545

    def test_compile_drug_does_not_import_two_route(self):
        """The two-route model must not be wired in. If someone wires it, this
        test fails and they are forced to re-run the blast radius."""
        import inspect

        from tmesim.molecular import compile_drug
        src = inspect.getsource(compile_drug)
        assert "two_route" not in src, (
            "compile_drug now references two_route. The two-route model is "
            "DEFAULT-OFF by decision (results/two_route_transport.md 7); wiring "
            "it in changes every DrugParams.diffusivity and moves alpha.")

    def test_drugparams_diffusivity_is_bit_identical(self):
        """BITWISE pin on the quantity the two-route model would change.

        184.75 um^2/s is the value results/cellular_exclusion.md 1 and
        results/diffusivity_falsification.md were computed against, recorded to
        full precision here. This is the number that flows into doxorubicin's
        `L = sqrt(D_eff/k)` and therefore into alpha; if it moves, alpha moves.
        """
        from tmesim.molecular import compile_drug
        from tmesim.molecular.library import get_molecule

        dp = compile_drug.molecule_to_drugparams(get_molecule("doxorubicin"))
        assert dp.diffusivity == 184.74861430191334

    def test_existing_ecm_mesh_functions_unchanged(self):
        """The additive wiring must not perturb the shipped ECM legs."""
        a = ecm_mesh.ogston_renkin_hindrance(0.4)
        assert a["effective_hindrance"] == pytest.approx(
            ecm_mesh.ogston_partition(0.4) * a["hindrance_diffusion"], rel=0, abs=0)
        assert ecm_mesh.tissue_hindrance(0.4, 0.0)["total_hindrance"] == \
            pytest.approx(a["effective_hindrance"], rel=0, abs=0)


# --------------------------------------------------------------------------- #
# 2. THE LIMITS. This is the check on the algebra.
# --------------------------------------------------------------------------- #
class TestLimits:
    def test_zero_permeability_reduces_to_cell_exclusion_exactly(self):
        shipped = ecm_mesh.cell_exclusion(PHI)["hindrance"]
        got = two_route_factor(PHI, SIGMA_M, D_CYTO, 0.0)["hindrance"]
        assert abs(got - shipped) < 1e-15

    @pytest.mark.parametrize("phi", [0.0, 0.1, 0.4, 0.6, 0.7892, 0.9])
    def test_zero_permeability_reduces_at_every_phi(self, phi):
        shipped = ecm_mesh.cell_exclusion(phi)["hindrance"]
        got = two_route_factor(phi, SIGMA_M, D_CYTO, 0.0)["hindrance"]
        assert abs(got - shipped) < 1e-15

    def test_tiny_permeability_approaches_the_impermeable_limit(self):
        shipped = ecm_mesh.cell_exclusion(PHI)["hindrance"]
        assert two_route_factor(PHI, SIGMA_M, D_CYTO, 1e-12)["hindrance"] == \
            pytest.approx(shipped, abs=1e-9)

    def test_no_cells_is_a_no_op(self):
        assert two_route_factor(0.0, SIGMA_M, D_CYTO, 1.0)["hindrance"] == \
            pytest.approx(1.0, abs=1e-15)

    def test_zero_contrast_tissue_is_transparent(self):
        got = two_route_factor(PHI, SIGMA_M, SIGMA_M, float("inf"))["hindrance"]
        assert got == pytest.approx(1.0, abs=1e-15)

    @pytest.mark.parametrize("p", [1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e4])
    def test_identity_with_hasselman_johnson_1987(self, p):
        """Maxwell-with-series-resistance == Hasselman & Johnson (1987) Eqn 8.

        doi:10.1177/002199838702100602. Two independent algebraic forms of the
        same published result must agree to machine precision.
        """
        mine = two_route_factor(PHI, SIGMA_M, D_CYTO, p)["hindrance"]
        theirs = hasselman_johnson_hindrance(
            PHI, SIGMA_M, D_CYTO, DEFAULT_CELL_RADIUS_UM, p)
        assert abs(mine - theirs) < 1e-14

    def test_dilute_maxwell_slope_survives_the_generalisation(self):
        """At phi -> 0 with impermeable cells the slope is the classic -3/2."""
        h = two_route_factor(1e-4, SIGMA_M, D_CYTO, 0.0)["hindrance"]
        assert abs((1.0 - h) / 1e-4 - 1.5) < 1e-3


# --------------------------------------------------------------------------- #
# 3. MONOTONICITY AND BOUNDS -- the model must be physical everywhere.
# --------------------------------------------------------------------------- #
class TestPhysicality:
    def test_monotone_increasing_in_membrane_permeability(self):
        hs = [two_route_factor(PHI, SIGMA_M, D_CYTO, p)["hindrance"]
              for p in (0.0, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e3, 1e6)]
        assert all(hs[i] <= hs[i + 1] + 1e-15 for i in range(len(hs) - 1))

    def test_bounded_by_the_infinite_permeability_ceiling(self):
        ceil_ = two_route_ceiling(PHI, SIGMA_M, D_CYTO)["hindrance"]
        for p in (0.0, 1.0, 1e3, 1e9):
            assert two_route_factor(PHI, SIGMA_M, D_CYTO, p)["hindrance"] <= \
                ceil_ + 1e-12

    def test_never_exceeds_one_when_cells_are_slower_than_interstitium(self):
        for p in (0.0, 1.0, 1e6):
            assert 0.0 < two_route_factor(PHI, SIGMA_M, D_CYTO, p)["hindrance"] <= 1.0

    def test_monotone_decreasing_in_cell_volume_fraction(self):
        hs = [two_route_factor(f, SIGMA_M, D_CYTO, 10.0)["hindrance"]
              for f in (0.0, 0.2, 0.4, 0.6, 0.8)]
        assert all(hs[i] > hs[i + 1] for i in range(len(hs) - 1))

    @pytest.mark.parametrize("bad", [
        dict(cell_volume_fraction=1.0), dict(cell_volume_fraction=-0.1),
        dict(interstitial_D_um2s=0.0), dict(membrane_permeability_um_s=-1.0),
        dict(cell_radius_um=0.0),
    ])
    def test_rejects_unphysical_inputs(self, bad):
        kw = dict(cell_volume_fraction=PHI, interstitial_D_um2s=SIGMA_M,
                  cytoplasm_D_um2s=D_CYTO, membrane_permeability_um_s=1.0)
        kw.update(bad)
        with pytest.raises(ValueError):
            two_route_factor(**kw)


# --------------------------------------------------------------------------- #
# 4. THE POINT OF THE EXERCISE: spread at FIXED molecular size.
# --------------------------------------------------------------------------- #
class TestSpreadAtFixedSize:
    def test_produces_spread_where_the_impermeable_law_is_flat(self):
        """Same radius, different chemistry. `cell_exclusion` is flat by
        construction; the two-route law must not be."""
        r = 0.45
        probes = [(-2.6, 2), (-0.4, 2), (0.0, 1), (1.5, 1), (2.9, 1)]
        vals = [two_route_diffusivity(SIGMA_M, r, lp, hd)["hindrance"]
                for lp, hd in probes]
        assert max(vals) / min(vals) > 2.0
        # the shipped law is exactly flat over the same probes
        flat = ecm_mesh.cell_exclusion(PHI)["hindrance"]
        assert all(v >= flat - 1e-15 for v in vals)

    def test_more_donors_never_helps(self):
        a = two_route_diffusivity(SIGMA_M, 0.45, 0.0, 1)["hindrance"]
        b = two_route_diffusivity(SIGMA_M, 0.45, 0.0, 3)["hindrance"]
        assert b < a

    def test_higher_logp_never_hurts(self):
        a = two_route_diffusivity(SIGMA_M, 0.45, 0.0, 1)["hindrance"]
        b = two_route_diffusivity(SIGMA_M, 0.45, 2.0, 1)["hindrance"]
        assert b > a

    def test_doxorubicin_is_effectively_paracellular_only(self):
        """6 H-bond donors at logP 0 -> the transcellular route is shut, so the
        two-route model reduces to the impermeable-cell answer for the drug alpha
        is calibrated on. This is why the blast radius is small."""
        res = two_route_diffusivity(184.75, 0.5393, 0.0, 6)
        assert res["transcellular_gain"] == pytest.approx(1.0, abs=5e-4)


# --------------------------------------------------------------------------- #
# 5. THE MEMBRANE LEG.
# --------------------------------------------------------------------------- #
class TestMembraneLeg:
    def test_constants_are_permeability_pys_published_ones(self):
        """No new tunable membrane constant was introduced."""
        from tmesim.molecular import permeability as perm
        assert two_route.ETA_MEMBRANE_PA_S == perm._ETA_MEMBRANE_PA_S
        assert two_route.BILAYER_THICKNESS_NM == perm._BILAYER_THICKNESS_NM
        assert two_route.LOGK_MEMBRANE_MAX == perm._LOGP_KP_MAX
        assert two_route.ETA_WATER_PA_S == perm._ETA_WATER_PA_S

    def test_upper_partition_clamp_binds(self):
        hi = membrane_permeability_um_s(9.0, 0, 0.45)
        at = membrane_permeability_um_s(3.0, 0, 0.45)
        assert hi["clamped"] and hi["permeability_um_s"] == \
            pytest.approx(at["permeability_um_s"], rel=1e-12)

    def test_no_lower_clamp(self):
        """A floor at K = 1e-3 would force the transcellular route open for
        exactly the hydrophilic compounds whose exclusion from it is the
        mechanism under test."""
        assert membrane_permeability_um_s(-2.0, 4, 0.45)["log_k_membrane"] == \
            pytest.approx(-6.0)

    def test_cytoplasm_diffusivity_uses_water_not_interstitial_viscosity(self):
        """Kao 1993's 0.27 is against WATER; using the engine's interstitial
        1.2e-3 Pa.s would apply an interstitial claim to the cytoplasm."""
        from tmesim.molecular import stokes_einstein
        d = cytoplasm_diffusivity_um2s(0.45)
        d_water = stokes_einstein.free_diffusivity(
            0.45, viscosity_pa_s=two_route.ETA_WATER_PA_S)
        assert d == pytest.approx(0.27 * d_water, rel=1e-12)

    def test_antibody_sized_solute_has_a_shut_membrane(self):
        res = two_route_diffusivity(5.0, 5.5, -1.0, 12)
        assert res["transcellular_gain"] == pytest.approx(1.0, abs=1e-9)


# --------------------------------------------------------------------------- #
# 6. THE ecm_mesh WIRING.
# --------------------------------------------------------------------------- #
class TestEcmMeshWiring:
    def test_two_route_wiring_reduces_to_tissue_hindrance_when_membrane_is_shut(self):
        """A compound with no lipophilicity and many donors must reproduce the
        impermeable-cell `tissue_hindrance` to within a part in 1e3."""
        a = ecm_mesh.tissue_hindrance_two_route(0.45, 0.7892, logp=-5.0, hbd=8)
        b = ecm_mesh.tissue_hindrance(0.45, 0.7892)
        assert a["total_hindrance"] == pytest.approx(b["total_hindrance"], rel=1e-6)

    def test_two_route_wiring_beats_paracellular_for_a_lipophilic_compound(self):
        a = ecm_mesh.tissue_hindrance_two_route(0.45, 0.7892, logp=2.5, hbd=0)
        assert a["total_hindrance"] > a["paracellular_only_hindrance"]
        assert a["transcellular_gain"] > 1.0

    def test_zero_cells_reproduces_the_pure_ecm_law(self):
        a = ecm_mesh.tissue_hindrance_two_route(0.45, 0.0, logp=1.0, hbd=1)
        b = ecm_mesh.ogston_renkin_hindrance(0.45)
        assert a["total_hindrance"] == pytest.approx(b["effective_hindrance"],
                                                     rel=1e-12)


# --------------------------------------------------------------------------- #
# 7. THE PRE-REGISTRATION must still describe the live code.
# --------------------------------------------------------------------------- #
def test_parameters_match_the_pre_registration():
    """If a constant is later tuned to the Pruijn measurement, this fails."""
    import json
    import os

    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "results", "sweep", "tworoute_prereg.json")
    if not os.path.exists(path):
        pytest.skip("pre-registration file not present")
    p = json.load(open(path))["parameters"]
    assert two_route.DEFAULT_CELL_VOLUME_FRACTION == p["cell_volume_fraction"]["value"]
    assert two_route.DEFAULT_CELL_RADIUS_UM == p["cell_radius_um"]["value"]
    assert two_route.CYTOPLASM_D_OVER_WATER == p["cytoplasm_D_over_water"]["value"]
    assert two_route.ETA_MEMBRANE_PA_S == p["membrane_viscosity_pa_s"]["value"]
    assert two_route.BILAYER_THICKNESS_NM == p["bilayer_thickness_nm"]["value"]
    assert two_route.LOGK_MEMBRANE_MAX == p["logk_membrane_max"]["value"]
    assert two_route.HBD_LOGP_DECREMENT == p["hbd_logp_decrement"]["value"]
