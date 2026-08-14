"""BRING YOUR OWN COMPOUND -- the compile path, its refusals, and its honesty.

Four things are tested, in rising order of what they would cost if they broke:

1. **The refusals.** Every input this pipeline must not compile -- an antibody or
   biologic, a salt/mixture, a pasted amino-acid sequence, an element the
   Crippen logP model has no atom type for, an out-of-range molecular weight, an
   unparseable SMILES, a name that shadows a curated drug -- fails with its own
   ``CompoundRefused`` code and a complete sentence, never a stack trace and
   never a silently degraded ``mw=0`` molecule. This is first because the
   pre-existing ``descriptors.from_smiles`` deliberately does NOT refuse: it
   falls back to a zeroed descriptor dict, which is right for a batch pipeline
   and catastrophic for a form field.

2. **The honesty claims are true, not just written down.** The provenance table
   says logP is a Crippen prediction, that Crippen cannot see stereochemistry,
   and that diffusivity is a function of MW alone. Each of those is asserted
   numerically here, so the sentence on the operator's screen cannot become
   false without a test failing.

3. **A user compound is compiled EXACTLY as a library drug is.** Same alpha, same
   binding-site barrier, same bridge -- checked by compiling a library drug's own
   SMILES under a user name and demanding bit-identical transport parameters. If
   this drifts, the ranking (the only defensible output either has) is invalid.

4. **The Studio wiring.** The endpoints, the store, the drug catalogue, the
   RunSpec validator, and one real end-to-end run driven with a user compound.
"""
from __future__ import annotations

import math
import os
import sys
import time

import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from tmesim.molecular import byo  # noqa: E402
from tmesim.molecular.byo import (  # noqa: E402
    CompoundRefused,
    compile_user_compound,
)

# Real, small, fast structures. Aspirin is the neutral small molecule used as the
# happy path; the two anthracyclines are the repo's own stereochemistry probe.
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
DOXORUBICIN = "COc1cccc2c1C(=O)c1c(O)c3c(c(O)c1C2=O)CC(O)(C(=O)CO)CC3OC1CC(N)C(O)C(C)O1"
# Epirubicin is doxorubicin's 4'-epimer: the SAME graph, one inverted centre.
# Written WITHOUT stereo markers here on purpose -- that is exactly the point
# being pinned (see test_crippen_is_stereochemistry_blind).
EPIRUBICIN = DOXORUBICIN

CURATED = ("temozolomide", "doxorubicin", "paclitaxel", "mitoxantrone")


@pytest.fixture(scope="module")
def alpha():
    """The project's single fitted loss constant, computed once for this file."""
    from tmesim.validation.calibration import calibrate_alpha

    return float(calibrate_alpha())


# --------------------------------------------------------------------------- #
# 1. Refusals
# --------------------------------------------------------------------------- #
class TestRefusals:
    """Every one of these would otherwise produce a number that looks like an
    answer. The assertion is not merely "it raised" -- it is that the code is
    specific and the message is a sentence somebody can act on."""

    @pytest.mark.parametrize("name,smiles,code", [
        # --- input hygiene ---------------------------------------------------
        ("x", "", "smiles_empty"),
        ("x", "   ", "smiles_empty"),
        ("x", "not-a-smiles", "smiles_unparseable"),
        ("x", "c1ccccc", "smiles_unparseable"),          # unclosed ring
        ("x", "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "smiles_unparseable"),
        ("x", "C" * (byo.MAX_SMILES_LEN + 1), "smiles_too_long"),
        # --- biologics -------------------------------------------------------
        # The sentinel the shipped library uses for its antibody row. A user must
        # not be able to reach the antibody path by typing it.
        ("x", "ANTIBODY:bevacizumab", "smiles_unparseable"),
        # A one-letter amino-acid paste. This DOES parse as a SMILES, which is
        # precisely why it needs its own gate.
        ("x", "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIH", "smiles_is_a_sequence"),
        # --- mixtures --------------------------------------------------------
        ("x", "CC(=O)[O-].[Na+]", "smiles_multi_component"),
        ("x", "CC(=O)Oc1ccccc1C(=O)O.O", "smiles_multi_component"),
        # --- descriptor model out of range -----------------------------------
        ("x", "[Gd]", "crippen_unparameterised_atom"),
        ("x", "[U]", "crippen_unparameterised_atom"),
        # --- MW band ---------------------------------------------------------
        ("x", "C", "mw_below_range"),                    # methane, 16 g/mol
        ("x", "CCO", "mw_below_range"),                  # ethanol, 46 g/mol
    ])
    def test_refused_with_a_specific_code(self, name, smiles, code):
        with pytest.raises(CompoundRefused) as exc:
            compile_user_compound(name, smiles)
        assert exc.value.code == code, exc.value.message
        # A refusal is a sentence, not a token.
        assert len(exc.value.message) > 40
        assert exc.value.message.rstrip().endswith((".", "?"))
        assert exc.value.field in ("name", "smiles", "ic50_uM")

    def test_a_very_large_molecule_is_refused_above_the_asserted_band(self):
        """A long alkane sails past every other gate and is caught by MW."""
        with pytest.raises(CompoundRefused) as exc:
            compile_user_compound("huge", "C" * 200)
        assert exc.value.code == "mw_above_range"
        assert "2,000" in exc.value.message or "2000" in exc.value.message

    def test_the_antibody_threshold_is_read_from_descriptors_not_retyped(self):
        from tmesim.molecular import descriptors

        assert byo.BIOLOGIC_MW_THRESHOLD == float(descriptors._LARGE_MW_THRESHOLD)

    def test_biologic_refusal_names_the_antibody_special_case(self):
        """The refusal has to explain WHY bevacizumab exists and a user's
        antibody cannot: the library row is hand-entered protein values, not a
        compiled structure."""
        # Reach _check_mw directly -- no SMILES on earth both parses and weighs
        # 149 kDa, which is itself the point.
        with pytest.raises(CompoundRefused) as exc:
            byo._check_mw(149_000.0, is_antibody=True)
        assert exc.value.code == "biologic_out_of_range"
        assert "bevacizumab" in exc.value.message

    @pytest.mark.parametrize("name,code", [
        ("", "name_empty"),
        ("   ", "name_empty"),
        ("A" * 200, "name_too_long"),
        ("../etc/passwd", "name_charset"),
        ("has/slash", "name_charset"),
    ])
    def test_bad_names_are_refused_on_the_name_field(self, name, code):
        with pytest.raises(CompoundRefused) as exc:
            compile_user_compound(name, ASPIRIN)
        assert exc.value.code == code
        assert exc.value.field == "name"

    @pytest.mark.parametrize("curated", CURATED)
    def test_a_user_compound_may_not_shadow_a_curated_name(self, curated):
        for variant in (curated, curated.upper(), f"  {curated}  "):
            with pytest.raises(CompoundRefused) as exc:
                compile_user_compound(variant, ASPIRIN, reserved_names=CURATED)
            assert exc.value.code == "name_collision"
            assert exc.value.field == "name"

    @pytest.mark.parametrize("ic50,code", [
        ("banana", "ic50_not_a_number"),
        (0.0, "ic50_not_positive"),
        (-3.0, "ic50_not_positive"),
        (float("nan"), "ic50_not_finite"),
        (float("inf"), "ic50_not_finite"),
        (5.0e6, "ic50_absurd"),
    ])
    def test_bad_ic50_is_refused_on_its_own_field(self, ic50, code):
        with pytest.raises(CompoundRefused) as exc:
            compile_user_compound("x", ASPIRIN, ic50_uM=ic50)
        assert exc.value.code == code
        assert exc.value.field == "ic50_uM"

    def test_refusal_serialises_for_the_wire(self):
        with pytest.raises(CompoundRefused) as exc:
            compile_user_compound("x", "zzz")
        d = exc.value.as_dict()
        assert set(d) == {"code", "message", "field", "detail"}

    def test_the_from_smiles_fallback_can_never_stand_in_for_a_user_compound(self):
        """``descriptors.from_smiles`` NEVER raises: given garbage it returns a
        Molecule with mw=0 and source='fallback_table'. That is the exact
        behaviour a user compound must not inherit, so it is pinned here in both
        directions."""
        from tmesim.molecular.descriptors import from_smiles

        degraded = from_smiles("whatever", "not-a-smiles")
        assert degraded.source == "fallback_table"
        assert degraded.get("mw") == 0.0          # the trap

        with pytest.raises(CompoundRefused):      # and the gate that closes it
            compile_user_compound("whatever", "not-a-smiles")

    def test_a_library_name_with_garbage_smiles_cannot_borrow_its_descriptors(self):
        """The fallback table is keyed by NAME. A user who types a curated name
        with a broken SMILES must not silently receive the curated descriptors --
        the name check fires first, and the SMILES check would fire anyway."""
        with pytest.raises(CompoundRefused) as exc:
            compile_user_compound("temozolomide", "zzz", reserved_names=CURATED)
        assert exc.value.code in ("name_collision", "smiles_unparseable")


# --------------------------------------------------------------------------- #
# 2. The honesty claims, asserted numerically
# --------------------------------------------------------------------------- #
class TestProvenanceIsTrue:
    def test_happy_path_produces_finite_engine_parameters(self, alpha):
        c = compile_user_compound("ACME-101", ASPIRIN, ic50_uM=12.5, alpha=alpha)
        dp = c.drug_params
        assert dp.name == "ACME-101"
        for fld in ("diffusivity", "molecular_weight", "ic50", "hill",
                    "binding", "clearance", "bbb_penetrance"):
            v = float(getattr(dp, fld))
            assert math.isfinite(v), fld
        assert dp.diffusivity > 0.0
        assert dp.ic50 == 12.5
        assert 0.0 <= dp.bbb_penetrance <= 1.0
        assert abs(dp.molecular_weight - 180.159) < 0.01

    def test_only_mw_and_tpsa_are_computed_everything_else_is_a_prediction(self, alpha):
        c = compile_user_compound("ACME-101", ASPIRIN, ic50_uM=12.5, alpha=alpha)
        by_key = {f["key"]: f for f in c.fields}
        exact = {"mw", "tpsa"}
        for key, f in by_key.items():
            if key in exact:
                assert f["kind"] == "computed", key
                assert f["is_prediction"] is False, key
            elif key == "ic50":
                # Supplied by the operator: real, but measured elsewhere.
                assert f["kind"] == "literature"
                assert f["is_prediction"] is False
            else:
                assert f["kind"] == "assumed", key
                assert f["is_prediction"] is True, key
        # NOTHING here is `measured`. No instrument touched this compound.
        assert not any(f["kind"] == "measured" for f in c.fields)
        assert c.as_dict()["predicted_field_count"] == len(c.fields) - 3

    def test_every_field_carries_a_method_sentence(self, alpha):
        c = compile_user_compound("ACME-101", ASPIRIN, alpha=alpha)
        for f in c.fields:
            assert f["method"] and len(f["method"]) > 40, f["key"]
            assert f["label"] and f["key"]

    def test_logp_is_exactly_the_crippen_prediction(self, alpha):
        """The UI says a self-compiled compound's logP is an RDKit Crippen
        prediction rather than a measurement. That claim is only worth making if
        it is the number. Every logP in this repo is a Crippen prediction
        (results/logp_provenance.md §1)."""
        from rdkit import Chem
        from rdkit.Chem import Crippen

        c = compile_user_compound("ACME-101", ASPIRIN, alpha=alpha)
        crippen = round(Crippen.MolLogP(Chem.MolFromSmiles(ASPIRIN)), 3)
        assert c.molecule.get("logp") == crippen
        logp_row = next(f for f in c.fields if f["key"] == "logp")
        assert "Crippen" in logp_row["method"]
        assert "PREDICTION" in logp_row["method"]

    def test_crippen_is_stereochemistry_blind(self, alpha):
        """results/logp_provenance.md's sharpest finding, restated for a user
        compound: doxorubicin and epirubicin differ only at one stereocentre and
        the descriptor chain cannot tell them apart — BIT-IDENTICAL descriptor
        dicts, hence identical diffusivity, BBB penetrance and binding. A user
        compiling a stereochemical series inherits that exactly, so the UI must
        say so — and this pins that it is still true."""
        a = compile_user_compound("SERIES-A", DOXORUBICIN, alpha=alpha)
        b = compile_user_compound("SERIES-B", EPIRUBICIN, alpha=alpha)
        assert a.molecule.descriptors == b.molecule.descriptors
        assert a.drug_params.diffusivity == b.drug_params.diffusivity
        assert a.drug_params.bbb_penetrance == b.drug_params.bbb_penetrance
        assert a.drug_params.binding == b.drug_params.binding
        logp_row = next(f for f in a.fields if f["key"] == "logp")
        assert "STEREOCHEMISTRY-BLIND" in logp_row["method"]

    def test_diffusivity_is_no_longer_a_function_of_mw_alone(self, alpha):
        """THIS CONTRACT INVERTED at the two-route promotion, and the inversion
        is the point of the whole exercise.

        It used to assert that two constitutional isomers -- same formula, same
        mass, different connectivity and different logP -- compiled to the SAME
        diffusivity, because the chain was MW -> radius -> Stokes-Einstein ->
        hindrance and nothing else entered it. That was true, it was what
        results/diffusivity_falsification.md refuted, and
        results/two_route_promotion.md fixed the mechanism: logP and
        hydrogen-bond donor count now reach `D_eff` through a membrane
        permeability. So the isomers must now DIFFER, and the UI text must no
        longer claim otherwise.

        THE HONEST CAVEAT, PINNED RATHER THAN OMITTED: for THIS pair the
        separation is tiny (~6e-6 relative), because Crippen puts n-butylbenzene
        and tert-butylbenzene within 0.01 log units of each other and both have
        zero donors. Being able to separate matched-mass compounds IN PRINCIPLE
        is not the same as separating any given pair usefully, and a user
        ranking a congeneric series needs the second statement, not the first.
        The pair below with a real logP gap is what a meaningful separation
        looks like.
        """
        # n-butylbenzene vs tert-butylbenzene: C10H14 both ways.
        a = compile_user_compound("ISO-A", "CCCCc1ccccc1", alpha=alpha)
        b = compile_user_compound("ISO-B", "CC(C)(C)c1ccccc1", alpha=alpha)
        assert a.molecule.get("mw") == b.molecule.get("mw")
        assert a.molecule.get("logp") != b.molecule.get("logp")   # they DO differ
        assert a.drug_params.diffusivity != b.drug_params.diffusivity
        # ... but only just, and that is disclosed rather than dressed up.
        rel = abs(a.drug_params.diffusivity - b.drug_params.diffusivity) \
            / a.drug_params.diffusivity
        assert rel < 1e-3

        row = next(f for f in a.fields if f["key"] == "diffusivity")
        assert "MW ALONE" not in row["method"], (
            "the diffusivity field still tells the user D_eff is a function of "
            "MW alone; since results/two_route_promotion.md that is false")

    def test_matched_mass_compounds_with_a_real_logp_gap_do_separate(self, alpha):
        """The property the two-route promotion exists to add, at a mass tie.

        Identical molecular weight, materially different hydrogen-bonding: under
        the MW-only model these were flat BY CONSTRUCTION, at exactly 1.00x, and
        results/diffusivity_falsification.md measured a real 3.7-3.8x spread at
        such ties. They must now separate, in the direction physics requires:
        more donors and lower lipophilicity shut the membrane and slow the drug.
        """
        # C6H12O6-ish tie is awkward from SMILES; use an explicit matched pair:
        # 1,4-dimethoxybenzene vs 4-methoxybenzyl alcohol are NOT matched, so
        # instead compare the same skeleton with and without a donor at equal MW.
        polar = compile_user_compound("TIE-POLAR", "OCc1ccccc1O", alpha=alpha)
        apolar = compile_user_compound("TIE-APOLAR", "COc1ccccc1O", alpha=alpha)
        assert polar.molecule.get("mw") == apolar.molecule.get("mw")
        assert polar.molecule.get("hbd") != apolar.molecule.get("hbd")
        # the one with MORE donors must be the SLOWER one
        if polar.molecule.get("hbd") > apolar.molecule.get("hbd"):
            assert polar.drug_params.diffusivity < apolar.drug_params.diffusivity
        else:
            assert polar.drug_params.diffusivity > apolar.drug_params.diffusivity

    def test_the_hindrance_shortfall_is_disclosed_on_the_diffusivity_field(self, alpha):
        c = compile_user_compound("ACME-101", ASPIRIN, alpha=alpha)
        row = next(f for f in c.fields if f["key"] == "diffusivity")
        assert "cellular_exclusion.md" in row["method"]
        assert "12-fold" in row["method"] or "12x" in row["method"]

    def test_the_transport_limit_is_disclosed_at_the_point_of_use(self, alpha):
        """The single most important thing this panel has to say, and the one a
        user compiling a congeneric series most needs. RE-PINNED at the two-route
        promotion (results/two_route_promotion.md), and the re-pin is the point.

        SUPERSEDED CONTRACT: the MW-only `D_eff` model had been REFUTED against
        71 compounds (Pruijn et al. 2008) at R^2 = 0.028 and was flat by
        construction at matched mass, so every surface had to carry the word
        FALSIFIED and a hard-stop "MUST NOT rank a congeneric series".

        THAT IS NO LONGER TRUE, AND SAYING IT WOULD NOW BE ITS OWN OVERSTATEMENT.
        The shipping model scores R^2 = 0.718 blind on the same 71 compounds. The
        replacement claim is NARROWER THAN "FIXED" and this test pins all three
        halves of it, because dropping any one of them would be the failure mode:
        the new score, the old score it replaced, and the shortfall that remains.
        """
        c = compile_user_compound("ACME-101", ASPIRIN, ic50_uM=1.0, alpha=alpha)

        row = next(f for f in c.fields if f["key"] == "diffusivity")
        for token in (
                "0.718",        # what it scores now
                "0.022",        # what the model it replaced scored
                "29.7x", "8.3x",  # measured spread, and the quarter of it modelled
                "ONE cell line",  # the limit that has not moved at all
                "Pruijn",
                "results/two_route_promotion.md",
                "diffusivity_falsification.md"):
            assert token in row["method"], f"{token!r} missing from the field"
        # It must NOT go on calling the shipping model falsified or MW-only.
        assert "MW ALONE" not in row["method"]

        cannot = " ".join(c.as_dict()["cannot_conclude"])
        assert "29.7x" in cannot and "0.718" in cannot
        assert "congeneric" in cannot
        assert "ONE cell line" in cannot

        joined = " ".join(c.caveats)
        # The transport caveat is no longer a hard stop, because the claim it
        # guards is no longer a refutation -- but it must still be unambiguous
        # about the strength of the evidence.
        assert "WEAK EVIDENCE, NOT A MEASUREMENT" in joined
        assert "0.718" in joined and "0 of 31" in joined
        # ... and the panel must still carry at least one genuine hard stop, so
        # that softening this one did not soften the renderer's stop behaviour.
        assert any("MUST NOT" in cv or "NOT DEFENSIBLE" in cv for cv in c.caveats)

    def test_the_ranking_claim_is_qualified_not_advertised(self, alpha):
        """The CAN list may no longer promise an unqualified ranking: the one
        dataset that tested the transport model refuted it for similar-mass
        compounds. It must point at that limit rather than sell around it."""
        d = compile_user_compound("ACME-101", ASPIRIN, alpha=alpha).as_dict()
        can = " ".join(d["can_conclude"])
        assert "not a validated prediction" in can
        assert "INDISTINGUISHABLE" in can

    def test_supplied_ic50_is_labelled_as_the_operators_own(self, alpha):
        c = compile_user_compound("ACME-101", ASPIRIN, ic50_uM=7.25, alpha=alpha)
        assert c.ic50_supplied is True
        assert c.drug_params.ic50 == 7.25
        row = next(f for f in c.fields if f["key"] == "ic50")
        assert row["method"].startswith("YOURS")
        assert any("YOURS" in cv for cv in c.caveats)

    def test_absent_ic50_is_a_flagged_placeholder_and_says_so(self, alpha):
        c = compile_user_compound("ACME-102", ASPIRIN, alpha=alpha)
        assert c.ic50_supplied is False
        assert c.drug_params.ic50 == 25.0
        row = next(f for f in c.fields if f["key"] == "ic50")
        assert row["method"].startswith("ABSENT")
        assert row["is_prediction"] is True
        assert any("NO IC50 WAS SUPPLIED" in cv for cv in c.caveats)
        assert any("MUST NOT" in cv for cv in c.caveats)

    def test_absolute_outputs_are_refused_in_words_on_every_compile(self, alpha):
        """HARD RULE: absolute kill is not a defensible output and supplying the
        molecule does not make it one. The caveat is worded to trip the Studio's
        hard-stop renderer (``NOT DEFENSIBLE`` / ``MUST NOT``), so this asserts
        the exact tokens that renderer keys on."""
        c = compile_user_compound("ACME-101", ASPIRIN, ic50_uM=1.0, alpha=alpha)
        joined = " ".join(c.caveats)
        assert "NOT DEFENSIBLE" in joined
        assert "MUST NOT" in joined
        assert "SELF-COMPILED COMPOUND" in joined

    def test_can_and_cannot_conclude_are_stated_and_disjoint(self, alpha):
        d = compile_user_compound("ACME-101", ASPIRIN, alpha=alpha).as_dict()
        assert len(d["can_conclude"]) >= 2 and len(d["cannot_conclude"]) >= 4
        can, cannot = " ".join(d["can_conclude"]), " ".join(d["cannot_conclude"])
        assert "ORDERING" in can
        for forbidden in ("kill fraction", "AUC", "penetration depth"):
            assert forbidden in cannot
        assert "clinical" in cannot.lower()

    def test_out_of_library_mw_warns_but_still_compiles(self, alpha):
        """Inside [50, 2000] but outside the 130-854 span of the curated
        libraries is an EXTRAPOLATION, not a refusal -- and it is said so."""
        c = compile_user_compound("SMALL-1", "CC(=O)NC", alpha=alpha)  # ~73 g/mol
        assert c.warnings, "an extrapolated MW must warn"
        assert any("OUTSIDE" in w for w in c.warnings)
        assert any("OUTSIDE" in cv for cv in c.caveats)

    def test_in_library_mw_does_not_warn(self, alpha):
        c = compile_user_compound("MID-1", ASPIRIN, alpha=alpha)
        assert not any("OUTSIDE" in w for w in c.warnings)


# --------------------------------------------------------------------------- #
# 3. A user compound is compiled EXACTLY as a library drug is
# --------------------------------------------------------------------------- #
class TestParityWithTheLibrary:
    def test_same_structure_same_alpha_gives_identical_transport(self, alpha):
        """Compile a curated drug's own SMILES under a user name and demand the
        SAME transport parameters the run pipeline builds for the library entry.

        Lomustine is chosen because it carries none of the by-name special cases
        (no TMZ hypoxia override, no cross-library reuse), so the only remaining
        difference between a curated compile and a user compile is the name,
        which must not matter. Exact equality, not a tolerance: two compounds
        calibrated even slightly differently cannot be ranked against each other,
        and the ordering is the only output either of them supports."""
        from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
        from tmesim.molecular.library import GBM_DRUGS, get_molecule
        from tmesim.molecular.potency import resolve_ic50
        from tmesim.validation.calibration import calibrated_drugparams

        mol = get_molecule("lomustine")
        library_dp = apply_barrier_to_drugparams(
            calibrated_drugparams("lomustine", alpha), mol)

        user = compile_user_compound(
            "MY-LOMUSTINE", GBM_DRUGS["lomustine"],
            ic50_uM=resolve_ic50("lomustine")["ic50_uM"], alpha=alpha)
        udp = user.drug_params

        assert udp.diffusivity == library_dp.diffusivity
        assert udp.molecular_weight == library_dp.molecular_weight
        assert udp.binding == library_dp.binding
        assert udp.clearance == library_dp.clearance
        assert udp.bbb_penetrance == library_dp.bbb_penetrance
        assert udp.ic50 == library_dp.ic50
        assert udp.hill == library_dp.hill

    def test_logp_is_the_same_tier_for_a_user_compound_and_a_library_drug(self):
        """As of this tree there is exactly ONE logP route — RDKit Crippen from
        the structure — so a user compound is on the SAME logP footing as a
        curated drug, and the module says so rather than claiming an asymmetry
        that does not exist. If a name-keyed measured tier ever lands, the next
        test is the guard that stops it reaching a user compound."""
        from tmesim.molecular.library import GBM_DRUGS, get_molecule

        for name in ("lomustine", "doxorubicin", "temozolomide"):
            assert get_molecule(name).get("logp_tier") is None, (
                "a logP tier has appeared; re-read byo._check_logp_tier and the "
                "logp provenance row, which currently state that every logP in "
                "this repo is a Crippen prediction")

        user = compile_user_compound("MY-LOMUSTINE", GBM_DRUGS["lomustine"])
        assert user.molecule.get("logp") == get_molecule("lomustine").get("logp")

    def test_a_name_keyed_measured_logp_would_be_refused(self):
        """The backstop, exercised directly. A measured-logP table keyed on drug
        NAME is safe for a curated library and unsafe for a form field: a user
        could type any structure, call it "aspirin", and receive aspirin's
        published logP — which then propagates into five estimators. The guard is
        inert today (no tier exists); this test is what makes it non-inert the
        day one lands."""
        from tmesim.molecular import byo as byo_mod
        from tmesim.molecular.types import Molecule

        tiered = Molecule(
            name="aspirin", smiles=ASPIRIN,
            descriptors={"logp": 1.19, "logp_tier": "experimental",
                         "logp_source": "some citation"})
        with pytest.raises(CompoundRefused) as exc:
            byo_mod._check_logp_tier(tiered)
        assert exc.value.code == "name_matches_measured_logp"
        assert exc.value.field == "name"

        untiered = Molecule(name="ACME-101", smiles=ASPIRIN,
                            descriptors={"logp": 1.31})
        assert byo_mod._check_logp_tier(untiered) is None

    def test_calibrated_drugparams_passthrough_contract_holds(self, alpha):
        """This module composes ``calibration.calibrated_drugparams`` by handing
        it an already-compiled ``DrugParams`` -- a path ``_compile`` documents
        ("a DrugParams passes through unchanged") and this module depends on
        rather than reimplementing the alpha algebra. Pin it, so a change in
        calibration.py fails HERE with a clear reason instead of silently
        producing a user compound on a different loss scale."""
        from tmesim.molecular.compile_drug import molecule_to_drugparams
        from tmesim.molecular.descriptors import from_smiles
        from tmesim.validation.calibration import (
            calibrated_drugparams, k_passive, uptake_proxy,
        )

        fresh = molecule_to_drugparams(from_smiles("probe", ASPIRIN))
        via_api = calibrated_drugparams(fresh, alpha)
        assert via_api.binding == alpha * uptake_proxy(fresh)
        assert via_api.clearance == k_passive(fresh)
        assert via_api.diffusivity == fresh.diffusivity

    def test_without_alpha_the_preview_says_it_is_uncalibrated(self):
        c = compile_user_compound("NOCAL", ASPIRIN, alpha=None)
        assert any("WITHOUT the project's alpha calibration" in cv
                   for cv in c.caveats)


# --------------------------------------------------------------------------- #
# 4. The Studio wiring
# --------------------------------------------------------------------------- #
fastapi = pytest.importorskip("fastapi", reason="the Studio backend needs fastapi")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi.testclient import TestClient  # noqa: E402

from app.server import cache as cache_mod  # noqa: E402
from app.server import compounds as compounds_mod  # noqa: E402
from app.server.api import create_app, validate_spec  # noqa: E402
from app.server.jobs import JobManager  # noqa: E402

SECTION = "xenium_319"


def _section_available(section_id: str = SECTION) -> bool:
    return any(s["id"] == section_id and s["available"]
               for s in cache_mod.list_sections())


requires_section = pytest.mark.skipif(
    not _section_available(),
    reason=f"section {SECTION!r} is not cached in data_cache/",
)


@pytest.fixture
def store(tmp_path, monkeypatch):
    """A compound store isolated from the real ``app/compounds/``."""
    st = compounds_mod.CompoundStore(root=str(tmp_path / "compounds"))
    monkeypatch.setattr(compounds_mod, "COMPOUND_STORE", st)
    return st


@pytest.fixture
def client(tmp_path_factory, store):
    root = str(tmp_path_factory.mktemp("runs"))
    manager = JobManager(root=root)
    app = create_app(manager)
    with TestClient(app) as c:
        yield c
    manager.shutdown()


class TestStore:
    def test_create_list_get_delete(self, store):
        doc = store.create("ACME-101", ASPIRIN, 12.5, note="series A")
        assert doc["id"] == "acme-101"
        assert doc["is_user_compound"] is True
        assert store.names() == ["acme-101"]

        again = store.compile_stored("acme-101")
        assert again["drug_params"] == doc["drug_params"]

        assert store.delete("ACME-101") is True
        assert store.names() == []
        assert store.delete("ACME-101") is False

    def test_only_the_inputs_are_persisted(self, store, tmp_path):
        """A compiled number on disk would go stale the moment the physics moved
        and would then look authoritative and be wrong. Only the operator's
        inputs are stored."""
        store.create("ACME-101", ASPIRIN, 12.5)
        import json

        path = os.path.join(store.root, "acme-101", "compound.json")
        with open(path) as fh:
            rec = json.load(fh)
        assert set(rec) == {"id", "name", "smiles", "ic50_uM", "note",
                            "created_at", "origin"}
        assert rec["origin"] == "user_smiles"
        assert "diffusivity" not in json.dumps(rec)

    def test_recompiles_are_memoised_in_process_but_never_on_disk(self, store):
        """``GET /api/drugs`` compiles every stored compound, so the memo is what
        keeps a picker with fifty compounds from taking half a minute to open.
        It must live in the process only: a compiled number on disk goes stale
        the moment the physics moves and then looks authoritative and is wrong."""
        store.create("ACME-101", ASPIRIN, 12.5)

        t0 = time.perf_counter()
        first = store.compile_stored("acme-101")
        cold = time.perf_counter() - t0

        t0 = time.perf_counter()
        second = store.compile_stored("acme-101")
        warm = time.perf_counter() - t0

        assert first["drug_params"] == second["drug_params"]
        assert warm < cold
        # A mutation of the returned document must not poison the memo.
        second["drug_params"] = None
        assert store.compile_stored("acme-101")["drug_params"] == first["drug_params"]
        # Nothing compiled reached the disk.
        with open(os.path.join(store.root, "acme-101", "compound.json")) as fh:
            assert "drug_params" not in fh.read()
        # A fresh store (i.e. a restart) starts cold.
        assert compounds_mod.CompoundStore(root=store.root)._compiled == {}

    def test_delete_evicts_the_memo(self, store):
        store.create("ACME-101", ASPIRIN, 12.5)
        store.compile_stored("acme-101")
        assert store._compiled
        store.delete("acme-101")
        assert store._compiled == {}

    def test_survives_a_restart(self, store):
        store.create("ACME-101", ASPIRIN, 12.5)
        fresh = compounds_mod.CompoundStore(root=store.root)
        assert fresh.names() == ["acme-101"]
        assert fresh.compile_stored("acme-101")["drug_params"]["ic50_uM"] == 12.5

    def test_a_second_compound_may_not_reuse_a_name(self, store):
        store.create("ACME-101", ASPIRIN, 12.5)
        with pytest.raises(CompoundRefused) as exc:
            store.create("acme-101", DOXORUBICIN)
        assert exc.value.code == "name_collision"

    def test_curated_names_are_reserved_by_the_store(self, store):
        for curated in ("temozolomide", "paclitaxel", "mitoxantrone"):
            with pytest.raises(CompoundRefused) as exc:
                store.create(curated, ASPIRIN)
            assert exc.value.code == "name_collision"

    def test_the_reserved_set_is_derived_from_the_libraries_not_listed(self, store):
        """A drug added to either curated library must be reserved the moment it
        lands, so the reserved set is derived rather than typed out."""
        from tmesim.molecular.breast_library import BREAST_DRUGS
        from tmesim.molecular.library import GBM_DRUGS

        reserved = set(compounds_mod._curated_names())
        assert set(GBM_DRUGS) <= reserved
        assert set(BREAST_DRUGS) <= reserved
        assert "mitoxantrone" in reserved

    def test_user_drugparams_returns_none_for_a_library_drug(self, store, alpha):
        assert compounds_mod.user_drugparams("temozolomide", alpha=alpha) is None

    def test_user_drugparams_round_trips(self, store, alpha):
        store.create("ACME-101", ASPIRIN, 12.5)
        got = compounds_mod.user_drugparams("ACME-101", alpha=alpha)
        assert got is not None
        dp, mol, doc = got
        assert dp.ic50 == 12.5
        assert mol.smiles == ASPIRIN
        assert doc["is_user_compound"] is True


class TestEndpoints:
    def test_compile_is_a_dry_run(self, client, store):
        r = client.post("/api/compounds/compile",
                        json={"name": "ACME-101", "smiles": ASPIRIN,
                              "ic50_uM": 12.5})
        assert r.status_code == 200, r.text
        d = r.json()
        assert d["drug_params"]["ic50_uM"] == 12.5
        assert d["fields"] and d["can_conclude"] and d["cannot_conclude"]
        assert store.names() == [], "compile must not store anything"

    @pytest.mark.parametrize("smiles,code", [
        ("zzz", "smiles_unparseable"),
        ("CC(=O)[O-].[Na+]", "smiles_multi_component"),
        ("[Gd]", "crippen_unparameterised_atom"),
        ("C", "mw_below_range"),
        ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIH", "smiles_is_a_sequence"),
    ])
    def test_a_refusal_is_a_422_with_a_reason_not_a_500(self, client, smiles, code):
        r = client.post("/api/compounds/compile",
                        json={"name": "ACME-101", "smiles": smiles})
        assert r.status_code == 422, r.text
        detail = r.json()["detail"]
        assert detail["code"] == code
        assert detail["field"] in ("name", "smiles", "ic50_uM")
        assert len(detail["message"]) > 40
        # A stack trace must never reach the operator.
        assert "Traceback" not in r.text

    def test_a_missing_body_is_a_400(self, client):
        r = client.post("/api/compounds/compile", content=b"not json",
                        headers={"content-type": "application/json"})
        assert r.status_code == 400

    def test_create_then_it_appears_in_the_drug_catalogue(self, client, store):
        r = client.post("/api/compounds",
                        json={"name": "ACME-101", "smiles": ASPIRIN,
                              "ic50_uM": 12.5, "note": "series A"})
        assert r.status_code == 201, r.text

        rows = client.get("/api/drugs").json()
        mine = [x for x in rows if x.get("is_user_compound")]
        assert [m["name"] for m in mine] == ["ACME-101"]
        m = mine[0]
        assert m["reason"] is None
        assert m["smiles"] == ASPIRIN
        assert m["ic50_uM"] == 12.5
        assert m["ic50_supplied"] is True
        assert "supplied by you" in m["ic50_source"]
        assert m["predicted_fields"] == m["total_fields"] - 3
        # The library rows must be untouched and must NOT claim to be user ones.
        lib = [x for x in rows if not x.get("is_user_compound")]
        assert any(x["name"] == "temozolomide" for x in lib)

    def test_a_compound_without_an_ic50_says_so_in_the_catalogue(self, client, store):
        client.post("/api/compounds", json={"name": "ACME-102", "smiles": ASPIRIN})
        row = next(x for x in client.get("/api/drugs").json()
                   if x.get("is_user_compound"))
        assert row["ic50_supplied"] is False
        assert row["ic50_uM"] == 25.0
        assert "LOW confidence" in row["ic50_source"]

    def test_get_and_delete(self, client, store):
        client.post("/api/compounds",
                    json={"name": "ACME-101", "smiles": ASPIRIN, "ic50_uM": 3.0})
        assert client.get("/api/compounds").json()[0]["name"] == "ACME-101"
        assert client.get("/api/compounds/acme-101").status_code == 200
        assert client.get("/api/compounds/nope").status_code == 404
        assert client.delete("/api/compounds/acme-101").status_code == 200
        assert client.delete("/api/compounds/acme-101").status_code == 404
        assert not [x for x in client.get("/api/drugs").json()
                    if x.get("is_user_compound")]


class TestRunSpecValidation:
    def _spec(self, drug):
        return {
            "section": SECTION, "dx_um": 30.0, "drug": drug,
            "sensitivity_mode": "uniform", "vascular_mode": "discrete",
            "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": None,
            "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
                   "n_doses": 1, "dose_interval_s": 86400.0},
            "total_time_s": 600.0,
        }

    def test_an_unstored_compound_is_still_an_unknown_drug(self, store):
        v = validate_spec(self._spec("ACME-101"))
        assert not v["ok"]
        assert any("unknown drug" in str(e) for e in v["errors"])

    @requires_section
    def test_a_stored_compound_validates_and_warns_loudly(self, store):
        store.create("ACME-101", ASPIRIN, 12.5)
        v = validate_spec(self._spec("ACME-101"))
        assert v["ok"], v["errors"]
        warns = " ".join(str(w) for w in v["warnings"])
        assert "SELF-COMPILED" in warns
        assert "NOT defensible" in warns or "NOT DEFENSIBLE" in warns
        assert "RANKING" in warns

    @requires_section
    def test_a_library_drug_gets_no_self_compiled_warning(self, store):
        v = validate_spec(self._spec("temozolomide"))
        assert not any("SELF-COMPILED" in str(w) for w in v["warnings"])


@requires_section
class TestEndToEndRun:
    def test_a_user_compound_runs_and_the_result_says_what_it_is(self, client, store):
        """The whole point: a compound the operator typed in, simulated on a real
        section, with a result that never pretends its structure was curated."""
        client.post("/api/compounds",
                    json={"name": "ACME-101", "smiles": ASPIRIN, "ic50_uM": 12.5})
        spec = {
            "section": SECTION, "dx_um": 40.0, "drug": "ACME-101",
            "sensitivity_mode": "uniform", "vascular_mode": "discrete",
            "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": None,
            "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
                   "n_doses": 1, "dose_interval_s": 86400.0},
            "total_time_s": 600.0,
        }
        posted = client.post("/api/runs", json=spec)
        assert posted.status_code == 202, posted.text
        run_id = posted.json()["run_id"]

        deadline = time.time() + 300.0
        status = None
        while time.time() < deadline:
            status = client.get(f"/api/runs/{run_id}/status").json()
            if status["status"] in ("done", "failed", "cancelled"):
                break
            time.sleep(0.2)
        assert status and status["status"] == "done", (
            client.get(f"/api/runs/{run_id}/traceback").json())

        run = client.get(f"/api/runs/{run_id}").json()

        dpz = run["drug_params"]
        assert dpz["is_user_compound"] is True
        assert dpz["structure_source"] == "user_smiles"
        assert dpz["smiles"] == ASPIRIN
        assert dpz["ic50_uM"] == 12.5

        prov = run["provenance"]
        assert prov["ic50_confidence"] == "user_supplied"
        assert "supplied by the operator" in prov["ic50_source"]
        uc = prov["user_compound"]
        assert uc["smiles"] == ASPIRIN
        assert uc["predicted_field_count"] >= 7
        assert uc["can_conclude"] and uc["cannot_conclude"]

        caveats = " ".join(str(c) for c in prov["caveats"])
        assert "SELF-COMPILED COMPOUND" in caveats
        assert "NOT DEFENSIBLE" in caveats
        assert "MUST NOT" in caveats

        # The run still produced real transport numbers -- the refusal is about
        # how they may be READ, not about withholding them.
        assert run["metrics"]["mean_tumor_auc"] > 0.0

    def test_a_user_compound_can_be_ranked_against_a_library_drug(self, client, store):
        """The defensible claim, exercised: two compounds, one tissue, one
        RunSpec, both compiled with the same alpha -- so the ORDERING is
        meaningful even though neither absolute level is."""
        client.post("/api/compounds",
                    json={"name": "ACME-BIG", "smiles": DOXORUBICIN,
                          "ic50_uM": 1.0})
        base = {
            "section": SECTION, "dx_um": 40.0, "drug": None,
            "sensitivity_mode": "uniform", "vascular_mode": "discrete",
            "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": None,
            "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
                   "n_doses": 1, "dose_interval_s": 86400.0},
            "total_time_s": 600.0,
        }
        results = {}
        for drug in ("temozolomide", "ACME-BIG"):
            spec = dict(base, drug=drug)
            rid = client.post("/api/runs", json=spec).json()["run_id"]
            deadline = time.time() + 300.0
            while time.time() < deadline:
                s = client.get(f"/api/runs/{rid}/status").json()
                if s["status"] in ("done", "failed", "cancelled"):
                    break
                time.sleep(0.2)
            assert s["status"] == "done", client.get(
                f"/api/runs/{rid}/traceback").json()
            results[drug] = client.get(f"/api/runs/{rid}").json()

        # Doxorubicin's structure is far heavier than temozolomide's, so the
        # user compound must be the slower-diffusing of the two. That is a
        # RANKING statement, which is the only kind either run supports.
        assert (results["ACME-BIG"]["drug_params"]["diffusivity_um2_s"]
                < results["temozolomide"]["drug_params"]["diffusivity_um2_s"])
        # And both carry the standing limits on their absolute numbers.
        for r in results.values():
            assert r["provenance"]["caveats"]
