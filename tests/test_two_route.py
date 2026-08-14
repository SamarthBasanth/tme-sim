"""Tests for the two-route (paracellular + transcellular) diffusivity model.

Three things are pinned here, in descending order of importance:

1. **It is the SHIPPING transport model, and its two numbers are pinned.**
   `DrugParams.diffusivity` for doxorubicin and `calibrate_alpha()` are both
   bitwise pins. THEY WERE RE-PINNED when the model was promoted to the default
   (`results/two_route_promotion.md`); the superseded values are kept beside the
   new ones in this module, in `tests/data/two_route_promotion_reference.json`
   and in `results/two_route_promotion.md`, and the assertions were changed
   rather than removed. `test_two_route_is_wired_into_compile_drug` replaces the
   old `test_compile_drug_does_not_import_two_route`, which existed to force
   exactly this re-run and did its job.

   HISTORY, so a future reader can tell a deliberate re-pin from a drift:

   | | superseded | current |
   |---|---|---|
   | `calibrate_alpha()` | 148.60815765445545 | 22.223477711172563 |
   | doxorubicin `DrugParams.diffusivity` | 184.74861430191334 | 27.92558262333494 |

   Cause: `compile_drug.build_profile` now multiplies the interstitial `D` by the
   two-route cell factor, which for doxorubicin is 0.15115449027238376 -- almost
   entirely the IMPERMEABLE-cell Maxwell factor, because doxorubicin's six
   hydrogen-bond donors shut its membrane and its transcellular gain is 1.0000002.
   `alpha` follows because it is fitted in closed form from `L = sqrt(D_eff/k)`.
   A stale agent once broke these pins by changing a molecular input that flows
   into doxorubicin's binding; the pins exist for that, and a change here that is
   NOT accompanied by a results/ document is that failure, not this one.
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
# 1. THE SHIPPING DEFAULT, and the two bitwise pins that move with it.
#
# These were `TestDefaultOff` until results/two_route_promotion.md. The class was
# renamed and its assertions re-pinned rather than deleted; the superseded values
# are asserted to be recoverable, which is a stronger guard than the old one.
# --------------------------------------------------------------------------- #
#: Re-pinned at the promotion. See this module's docstring for the cause.
ALPHA_SUPERSEDED = 148.60815765445545
ALPHA_CURRENT = 22.223477711172563
DOX_D_SUPERSEDED = 184.74861430191334
DOX_D_CURRENT = 27.92558262333494
#: doxorubicin's two-route factor. Overwhelmingly the impermeable-cell Maxwell
#: term: its transcellular gain is 1.0000002 because six donors shut its membrane.
DOX_TWO_ROUTE_SCALE = 0.15115449027238376


class TestShippingDefault:
    def test_alpha_is_pinned_at_its_repinned_value(self):
        """BITWISE pin. Re-pinned once, deliberately, at the promotion."""
        from tmesim.validation.calibration import calibrate_alpha
        assert calibrate_alpha() == ALPHA_CURRENT

    def test_two_route_is_wired_into_compile_drug(self):
        """The contract this test asserts was INVERTED at the promotion.

        It used to assert that `compile_drug` contained no reference to
        `two_route` at all, so that anyone wiring it in would be forced to re-run
        the blast radius before the suite went green. That is exactly what
        happened (results/two_route_promotion.md), so the test now pins the new
        contract: the model IS wired in, it IS on by default, and turning it off
        is still possible and still recovers the old number exactly.
        """
        import inspect

        from tmesim.molecular import compile_drug
        src = inspect.getsource(compile_drug)
        assert "two_route" in src, (
            "compile_drug no longer references two_route. The two-route model is "
            "the SHIPPING default (results/two_route_promotion.md); removing it "
            "silently returns every DrugParams.diffusivity to the MW-only model "
            "that results/diffusivity_falsification.md falsified.")
        assert compile_drug._TWO_ROUTE_DEFAULT is True
        assert (inspect.signature(compile_drug.build_profile)
                .parameters["two_route"].default
                is compile_drug._TWO_ROUTE_DEFAULT)

    def test_drugparams_diffusivity_is_bit_identical(self):
        """BITWISE pin on the quantity the two-route model changed.

        27.9256 um^2/s is the shipping value since the promotion. 184.7486 is what
        results/cellular_exclusion.md 1 and results/diffusivity_falsification.md
        were computed against, and is asserted below to be exactly recoverable.
        This is the number that flows into doxorubicin's `L = sqrt(D_eff/k)` and
        therefore into alpha; if it moves, alpha moves.
        """
        from tmesim.molecular import compile_drug
        from tmesim.molecular.library import get_molecule

        dp = compile_drug.molecule_to_drugparams(get_molecule("doxorubicin"))
        assert dp.diffusivity == DOX_D_CURRENT

    def test_the_superseded_value_is_still_exactly_recoverable(self):
        """Every published number predating the promotion must stay reproducible.

        This is what makes the re-pin auditable rather than a reset: the old
        model is not gone, it is one keyword away, and it is asserted BITWISE.
        """
        from tmesim.molecular import compile_drug
        from tmesim.molecular.library import get_molecule

        prof = compile_drug.build_profile(get_molecule("doxorubicin"),
                                          two_route=False)
        assert prof.effective_diffusivity == DOX_D_SUPERSEDED
        # and it is carried on every compile, so a caller never has to re-run one
        on = compile_drug.build_profile(get_molecule("doxorubicin"))
        assert on.interstitial_diffusivity == DOX_D_SUPERSEDED

    def test_the_repin_is_exactly_the_two_route_factor(self):
        """The new pins must be the OLD pins times the law, not new numbers.

        If this fails, something other than the two-route term moved D, and the
        re-pin in this module is not explained by the change it claims.
        """
        from tmesim.molecular import compile_drug
        from tmesim.molecular.library import get_molecule

        prof = compile_drug.build_profile(get_molecule("doxorubicin"))
        assert prof.two_route["hindrance"] == DOX_TWO_ROUTE_SCALE
        assert prof.effective_diffusivity == DOX_D_SUPERSEDED * DOX_TWO_ROUTE_SCALE
        # doxorubicin's membrane is SHUT: the move is the impermeable-cell term,
        # not the transcellular route. This is why the blast radius is small.
        assert prof.two_route["transcellular_gain"] == pytest.approx(1.0, abs=5e-4)
        assert prof.two_route["hindrance"] == pytest.approx(
            ecm_mesh.cell_exclusion(two_route.DEFAULT_CELL_VOLUME_FRACTION)
            ["hindrance"], rel=1e-4)

    def test_alpha_follows_D_through_its_own_closed_form(self):
        """alpha is a GAUGE COORDINATE (results/alpha_identifiability.md), so its
        value moving is expected. What must hold is that it moved for the stated
        reason and by the stated arithmetic, and NOT in proportion to D -- the
        un-amplified passive clearance term is subtracted first, which is why the
        alpha ratio (0.14954) is not the D ratio (0.15115)."""
        from tmesim.validation import calibration as C

        dp = C._compile("doxorubicin")
        d_eff = C._d_eff_tortuous(dp, C.DEFAULT_ECM)
        k_target = d_eff / (C.DOX_TARGET_EFOLD_UM ** 2)
        assert (k_target - C.k_passive(dp)) / C.uptake_proxy(dp) == ALPHA_CURRENT
        assert ALPHA_CURRENT / ALPHA_SUPERSEDED != pytest.approx(
            DOX_D_CURRENT / DOX_D_SUPERSEDED, rel=1e-6)

    def test_existing_ecm_mesh_functions_unchanged(self):
        """The wiring must not perturb the shipped ECM legs."""
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
        mechanism under test.

        Asserted against the LAW rather than against a literal, so that promoting
        `HBD_LOGP_DECREMENT` (1.0 -> 1.34, results/two_route_promotion.md) cannot
        make this test pass or fail for the wrong reason. What it pins is that
        nothing arrests the descent at permeability.py's -3.
        """
        got = membrane_permeability_um_s(-2.0, 4, 0.45)["log_k_membrane"]
        assert got == pytest.approx(-2.0 - 4.0 * two_route.HBD_LOGP_DECREMENT)
        assert got < -3.0 and not membrane_permeability_um_s(
            -2.0, 4, 0.45)["clamped"]

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
# 6b. EFFLUX composed with the two-route law.
#
# The point of these tests is that the composition is an ALGEBRAIC IDENTITY with
# no new free parameter, and that the identity has a falsifiable consequence:
# efflux cannot move a drug whose membrane is shut. See
# results/two_route_efflux.md.
# --------------------------------------------------------------------------- #
class TestEffluxComposition:
    def test_efflux_is_opt_in_and_bit_identical_when_off(self):
        """Omitting `efflux_ratio` and passing None must be BITWISE identical."""
        a = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
        b = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=None)
        assert a["hindrance"] == b["hindrance"]
        assert a["effective_D_um2s"] == b["effective_D_um2s"]

    def test_efflux_ratio_of_one_is_exactly_a_no_op(self):
        a = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
        b = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=1.0)
        assert a["hindrance"] == b["hindrance"]

    def test_efflux_never_boosts(self):
        """ER < 1 means net ABSORPTIVE transport, which this model does not
        represent. It must clamp, exactly as `efflux._apply_form` does, and never
        make a drug faster."""
        a = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
        for er in (0.01, 0.5, 0.999):
            b = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=er)
            assert b["hindrance"] == a["hindrance"]

    def test_the_slab_identity_that_makes_this_not_a_modelling_choice(self):
        """(C1-C2)/J = 2/P + (L/sigma_c)*ER for a cell with efflux at both faces.

        This is the whole justification for `K -> K/ER`: efflux multiplies the
        INTERIOR resistance by the efflux ratio and leaves the membrane
        resistance alone. If this fails, the composition is unjustified.
        """
        worst = 0.0
        for p_mem in (1e-3, 0.1, 1.0, 10.0, 1e3):
            for sig_c in (0.5, 20.0, 500.0):
                for length in (0.5, 5.0, 20.0):
                    for e_pump in (0.0, 0.5, 10.0, 500.0):
                        j = p_mem / (2.0 + (p_mem + e_pump) * length / sig_c)
                        ca = (p_mem - j) / (p_mem + e_pump)
                        cb = j / (p_mem + e_pump)
                        # all three steady-state flux statements hold at this J
                        assert p_mem * (1.0 - ca) - e_pump * ca == pytest.approx(j)
                        assert (sig_c / length) * (ca - cb) == pytest.approx(j)
                        assert p_mem * cb + e_pump * cb == pytest.approx(j)
                        er = 1.0 + e_pump / p_mem
                        want = 2.0 / p_mem + (length / sig_c) * er
                        worst = max(worst, abs(1.0 / j - want) / want)
        assert worst < 1e-12, worst

    def test_efflux_equals_reducing_the_partition(self):
        """`efflux_ratio=ER` must be identical to passing `cell_partition=K/ER`."""
        for er in (2.0, 7.5, 39.0, 100.0):
            a = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=er)
            b = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4,
                                      cell_partition=1.0 / er)
            assert a["hindrance"] == b["hindrance"]
        assert two_route.efflux_reduced_partition(1.0, 4.0) == 0.25
        assert two_route.efflux_reduced_partition(2.0, 4.0) == 0.5
        assert two_route.efflux_reduced_partition(1.0, 0.3) == 1.0   # clamped

    def test_a_shut_membrane_is_immune_to_efflux(self):
        """THE FALSIFIABLE PREDICTION. Doxorubicin (logP 0, 6 donors) has a
        membrane resistance ratio ~2e5, so the cell interior contributes nothing
        to the series and no efflux ratio, however large, can move it.

        This is why the previous transport model -- which had no membrane and no
        cell inclusions at all -- could not have been helped by efflux: the term
        was identically zero for EVERY drug, not just for this one.
        """
        a = two_route_diffusivity(SIGMA_M, 0.539, 0.0, 6)
        assert a["membrane_resistance_ratio"] > 1e4
        for er in (10.0, 1e3, 1e6):
            b = two_route_diffusivity(SIGMA_M, 0.539, 0.0, 6, efflux_ratio=er)
            assert b["hindrance"] == pytest.approx(a["hindrance"], rel=1e-4)
            assert b["transcellular_gain"] == pytest.approx(
                a["transcellular_gain"], rel=1e-4)
        # and the gain was never more than a part in 1e4 above 1 to begin with
        assert a["transcellular_gain"] < 1.0001

    def test_an_open_membrane_is_strongly_moved_by_efflux(self):
        """A compound whose membrane is genuinely open: its transcellular gain
        must collapse toward 1 as ER rises, monotonically.

        The probe is a bare lipophile (logP 3.74, ZERO donors) rather than
        paclitaxel-like (3.74, four donors), because at the promoted decrement
        paclitaxel's membrane is no longer open -- see the test below, which pins
        that as the finding it is rather than letting it silently break this one.
        """
        a = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 0)
        assert a["membrane_resistance_ratio"] < 1.0
        assert a["transcellular_gain"] > 3.0
        gains = [two_route_diffusivity(SIGMA_M, 0.626, 3.74, 0,
                                       efflux_ratio=er)["transcellular_gain"]
                 for er in (1.0, 2.0, 5.0, 20.0, 100.0, 1e4)]
        assert all(gains[i] >= gains[i + 1] for i in range(len(gains) - 1)), gains
        assert gains[-1] == pytest.approx(1.0, abs=1e-3)

    def test_the_sourced_decrement_closes_paclitaxels_membrane(self):
        """THE FINDING THAT COUPLED THE TWO BLOCKERS, pinned so it cannot rot.

        results/two_route_efflux.md 4.1: at the pre-registered convention (1.0)
        paclitaxel's membrane was nearly open -- resistance ratio 0.366, a 4.94x
        transcellular gain, and the largest single number in the blast radius,
        which is what made efflux look necessary. At the SOURCED decrement (1.34)
        4 donors cost 5.36 log units instead of 4.00, the resistance ratio becomes
        8.31, and the outlier largely disappears WITHOUT any efflux term. The 4.94x
        was substantially an artefact of an under-sourced donor penalty.
        """
        at_convention = two_route_diffusivity(
            SIGMA_M, 0.626, 3.74, 4,
            hbd_decrement=two_route.HBD_LOGP_DECREMENT_PREREGISTERED)
        shipping = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
        assert at_convention["membrane_resistance_ratio"] < 1.0
        assert shipping["membrane_resistance_ratio"] > 8.0
        assert at_convention["transcellular_gain"] > 3.0
        assert shipping["transcellular_gain"] < 2.0

    def test_efflux_rejects_nonsense(self):
        with pytest.raises(ValueError):
            two_route.efflux_reduced_partition(1.0, 0.0)
        with pytest.raises(ValueError):
            two_route.efflux_reduced_partition(1.0, -3.0)
        with pytest.raises(ValueError):
            two_route.efflux_reduced_partition(-1.0, 2.0)

    def test_none_of_the_pruijn_compounds_carries_an_efflux_ratio(self):
        """The commensurability fact, pinned so it cannot rot silently.

        If someone later adds an SN-series compound to EFFLUX_TABLE, the claim in
        results/two_route_efflux.md that efflux is unscoreable against Pruijn
        stops being true and this test forces a re-read.
        """
        from tmesim.molecular.efflux import EFFLUX_TABLE
        from experiments.diffusivity_falsification import PRUIJN_2008_TABLE_S1
        names = {row[0].strip().lower() for row in PRUIJN_2008_TABLE_S1}
        assert not (names & set(EFFLUX_TABLE)), (
            "a Pruijn compound now has an efflux assignment; "
            "results/two_route_efflux.md 2 must be re-run")


# --------------------------------------------------------------------------- #
# 7. THE PRE-REGISTRATION must still describe the live code.
# --------------------------------------------------------------------------- #
def test_parameters_match_the_pre_registration():
    """If a constant is later tuned to the Pruijn measurement, this fails.

    The ORIGINAL pre-registration is checked in full and is NEVER edited. The one
    parameter that moved at the promotion is checked here against the historical
    value that live code still carries by name, and against its own registration
    in `test_the_promoted_parameter_matches_its_own_registration`. The
    pre-registration was EXTENDED, not bypassed, and neither half is skipped.
    """
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
    assert (two_route.HBD_LOGP_DECREMENT_PREREGISTERED
            == p["hbd_logp_decrement"]["value"] == 1.0)


def test_the_promoted_parameter_matches_its_own_registration():
    """The sourced decrement must agree with the registration that sourced it.

    The value 1.34 is NOT free: it is 3.621 * 0.37, from two open-access Abraham
    LSER coefficient sets. If someone edits the constant without editing the
    registration -- or edits the registration to chase a score -- this fails.
    """
    import json
    import os

    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "results", "sweep", "tworoute_promotion_prereg.json")
    if not os.path.exists(path):
        pytest.skip("promotion registration file not present")
    ch = json.load(open(path))["parameters_changed"]["hbd_logp_decrement"]
    assert two_route.HBD_LOGP_DECREMENT == ch["value"] == 1.34
    assert two_route.HBD_LOGP_DECREMENT_PREREGISTERED == ch["superseded_value"]
    # the derivation is arithmetic, so check the arithmetic rather than the prose
    assert round(3.621 * 0.37, 3) == 1.34
    # and it must NOT be the value that maximises the Pruijn score, which is the
    # single most attackable property the parent file had.
    assert two_route.HBD_LOGP_DECREMENT != 1.0


# --------------------------------------------------------------------------- #
# 8. INTRACELLULAR RETARDATION -- the CAPACITY term.
#
# The point of these tests is that retardation is a statement about CAPACITY and
# not about conductance, that it is therefore arithmetically incapable of moving
# `transcellular_gain`, and that it is bitwise off by default. See
# results/intracellular_retardation.md.
# --------------------------------------------------------------------------- #
class TestIntracellularRetardation:
    def test_retardation_is_opt_in_and_bit_identical_when_off(self):
        """Omitting `bound_to_free_ratio` and passing None must be BITWISE identical."""
        a = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
        b = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, bound_to_free_ratio=None)
        assert a["hindrance"] == b["hindrance"]
        assert a["effective_D_um2s"] == b["effective_D_um2s"]
        assert a["retardation_factor"] == 1.0
        assert a["apparent_D_um2s"] == a["effective_D_um2s"]
        assert a["apparent_hindrance"] == a["hindrance"]

    def test_zero_binding_is_exactly_a_no_op(self):
        a = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
        b = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, bound_to_free_ratio=0.0)
        assert b["retardation_factor"] == 1.0
        assert b["apparent_D_um2s"] == a["effective_D_um2s"]

    def test_retardation_cannot_move_a_conductance(self):
        """THE test. `transcellular_gain` is a ratio of two CONDUCTANCES, and
        reversible binding does not touch a conductance. No value of K_b may
        move it. This is what makes "retardation reduces paclitaxel's gain"
        false as stated -- see results/intracellular_retardation.md."""
        base = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
        for kb in (0.0, 1.0, 1e2, 1e4, 1e8):
            r = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4,
                                      bound_to_free_ratio=kb)
            assert r["transcellular_gain"] == base["transcellular_gain"]
            assert r["hindrance"] == base["hindrance"]
            assert r["effective_D_um2s"] == base["effective_D_um2s"]

    def test_retardation_identity_against_the_capacity_ratio(self):
        """R must be exactly s_eff(K_b)/s_eff(0), the ratio of two volume averages.

        If this fails the term is not a capacity ratio and its derivation in
        two_route.retardation_factor does not apply.
        """
        worst = 0.0
        for phi in (0.0, 0.1, 0.4, PHI, 0.95):
            for k in (0.5, 1.0, 3.0):
                for kb in (0.0, 0.3, 3.0, 30.0, 300.0, 3000.0):
                    s0 = (1.0 - phi) + phi * k
                    s1 = (1.0 - phi) + phi * k * (1.0 + kb)
                    direct = s1 / s0
                    mine = two_route.retardation_factor(
                        kb, cell_volume_fraction=phi, cell_partition=k)
                    worst = max(worst, abs(mine - direct) / direct)
        assert worst < 1e-14, worst

    def test_retardation_is_monotone_and_never_speeds_anything_up(self):
        prev = 0.0
        for kb in (0.0, 0.1, 1.0, 10.0, 1000.0):
            r = two_route.retardation_factor(kb)
            assert r >= 1.0
            assert r >= prev
            prev = r

    def test_no_cells_means_no_retardation_however_hard_the_drug_binds(self):
        """phi -> 0 leaves nothing to bind to, so R must be exactly 1."""
        for kb in (0.0, 10.0, 1e6):
            assert two_route.retardation_factor(
                kb, cell_volume_fraction=0.0) == 1.0

    def test_rejects_nonsense_inputs(self):
        with pytest.raises(ValueError):
            two_route.retardation_factor(-1.0)
        with pytest.raises(ValueError):
            two_route.retardation_factor(1.0, cell_volume_fraction=1.0)
        with pytest.raises(ValueError):
            two_route.retardation_factor(1.0, cell_partition=-1.0)

    def test_structural_rule_reproduces_its_own_derivation(self):
        """K_b = 125 * 10**(0.072 logP^2 + 0.067 logP - 1.126), and the value at
        logP = 0 is the one the pre-registration calls out as out-of-domain."""
        for lp in (-2.62, -0.8, 0.0, 2.0, 3.74):
            expect = 125.0 * 10.0 ** (0.072 * lp * lp + 0.067 * lp - 1.126)
            assert abs(two_route.binding_capacity_from_logp(lp) - expect) < 1e-12
        assert abs(two_route.binding_capacity_from_logp(0.0) - 9.3521) < 1e-3

    def test_equilibration_time_is_infinite_for_a_sealed_cell(self):
        """A cell that admits nothing never loads its bound pool. This is the
        guard on the local-equilibrium assumption that the retardation factor
        rests on."""
        assert two_route.equilibration_time_s(10.0, 0.0) == float("inf")
        fast = two_route.equilibration_time_s(10.0, 1e4)
        slow = two_route.equilibration_time_s(10.0, 1e-4)
        assert slow > fast
        # more binding takes proportionally longer to load
        assert abs(two_route.equilibration_time_s(1.0, 1.0)
                   / two_route.equilibration_time_s(0.0, 1.0) - 2.0) < 1e-12

    def test_retardation_composes_with_efflux_without_interfering(self):
        """Efflux acts on the partition K (a conductance leg) and retardation on
        the capacity. Both may be on at once and neither may corrupt the other."""
        a = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=39.0)
        b = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=39.0,
                                  bound_to_free_ratio=30.0)
        assert b["hindrance"] == a["hindrance"]
        assert b["transcellular_gain"] == a["transcellular_gain"]
        assert b["apparent_D_um2s"] < a["apparent_D_um2s"]


def test_retardation_parameters_match_the_pre_registration():
    """The retardation pre-registration was written BEFORE the term was scored.
    If a constant is later tuned to Pruijn, this fails."""
    import json
    import os

    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "results", "sweep", "retard_prereg.json")
    if not os.path.exists(path):
        pytest.skip("retardation pre-registration file not present")
    reg = json.load(open(path))
    rule = reg["quantities_added"]["structural_rule_for_K_b"]["value"]
    assert "125" in rule and "0.072" in rule and "0.067" in rule and "1.126" in rule
    assert two_route.KILFORD_MG_EQUIV_PER_ML_CELL == 125.0
    assert two_route.HALLIFAX_C2 == 0.072
    assert two_route.HALLIFAX_C1 == 0.067
    assert two_route.HALLIFAX_C0 == -1.126
    # the term must still be OFF by default
    assert (reg["quantities_added"]["structural_rule_for_K_b"]
            ["adopted_in_code_as_default"]) is False
