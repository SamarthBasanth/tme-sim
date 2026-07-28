"""Tests for TPL — ``tmesim/protocol/`` and the ``/api/protocol/*`` routes.

Grouped by the thing being pinned:

* :class:`TestLexer` — tokens, positions, comments, strings.
* :class:`TestParserErrors` — every refusal carries a LINE and a COLUMN.
* :class:`TestUnits` — the frozen unit table and explicit conversion.
* :class:`TestDoseToPk` — the conversion and every assumption it states.
* :class:`TestSweepExpansion` — cross-product, counting, and the budget refusal.
* :class:`TestCapability` — a mode the panel cannot support is an ERROR.
* :class:`TestDrugSources` — smiles / library / file, and per-field provenance.
* :class:`TestRoutes` — the five HTTP routes, without running a solve.

Only :class:`TestExecution` runs a PDE solve, and it uses a coarse grid and
a short observation window so the suite stays fast.
"""
from __future__ import annotations

import json
import math
import os
import tempfile

import pytest

from tmesim.protocol import (
    METRICS,
    UNITS,
    Quantity,
    convert,
    convert_dose,
    parse,
    plan,
    schedule_times,
    tokenize,
    validate_source,
)
from tmesim.protocol.dosing import BODY_WEIGHT_KG, BSA_M2, PK_LITERATURE
from tmesim.protocol.errors import TplError
from tmesim.protocol.executor import DEFAULT_MAX_RUNS

# --------------------------------------------------------------------------- #
# Shared scripts
# --------------------------------------------------------------------------- #
TMZ_SMILES = "CN1C(=O)N2C=NC(=C2N=N1)C(=O)N"

STUPP = f'''drug TMZ from smiles "{TMZ_SMILES}"
  ic50 from gdsc
end

protocol stupp_adjuvant
  route oral
  dose 150 mg/m2
  every 24 h
  repeat 5
  then rest 23 d
  cycles 6
end

run TMZ with stupp_adjuvant on xenium_319
  dx 30 um
  state_mode oxygen
  sensitivity uniform
  report coverage, kill, mean_auc, penetration_length
end
'''


def errors(source: str):
    return [d for d in validate_source(source)["diagnostics"]
            if d["severity"] == "error"]


def first_error(source: str):
    errs = errors(source)
    assert errs, f"expected an error, got none for:\n{source}"
    return errs[0]


# --------------------------------------------------------------------------- #
class TestLexer:
    def test_tokens_carry_one_based_line_and_column(self):
        toks = tokenize("drug TMZ from library \"cisplatin\"\n  ic50 12 uM\nend\n")
        by_text = {t.text: t for t in toks if t.kind != "NEWLINE"}
        assert by_text["drug"].line == 1 and by_text["drug"].col == 1
        assert by_text["TMZ"].col == 6
        assert by_text["cisplatin"].kind == "STRING"
        assert by_text["ic50"].line == 2 and by_text["ic50"].col == 3
        assert by_text["uM"].line == 2
        # end_col is exclusive
        assert by_text["drug"].end_col == 5

    def test_compound_units_lex_as_one_word(self):
        kinds = {t.text for t in tokenize("dose 150 mg/m2\n")}
        assert "mg/m2" in kinds
        assert "um2/s" in {t.text for t in tokenize("diffusivity 12 um2/s\n")}

    def test_comments_and_blank_lines_are_skipped(self):
        toks = tokenize("# a comment\n\ndose 1 mg   # trailing\n")
        texts = [t.text for t in toks if t.kind not in ("NEWLINE", "EOF")]
        assert texts == ["dose", "1", "mg"]

    def test_percent_comma_and_comparators(self):
        assert [t.kind for t in tokenize("expect kill >= 5 %\n")][:5] == [
            "WORD", "WORD", "COMPARATOR", "NUMBER", "PERCENT"]
        assert "COMMA" in [t.kind for t in tokenize("report a, b\n")]

    def test_scientific_and_negative_numbers(self):
        assert [t.text for t in tokenize("clearance 3e-4\n")][1] == "3e-4"
        assert [t.text for t in tokenize("expect kill > -1\n")][3] == "-1"

    def test_unterminated_string_is_a_positioned_error(self):
        with pytest.raises(TplError) as exc:
            tokenize('drug X from smiles "CN1C\n')
        d = exc.value.diagnostics[0]
        assert d.code == "unterminated-string"
        assert d.line == 1 and d.col == 20

    def test_unexpected_character_is_a_positioned_error(self):
        with pytest.raises(TplError) as exc:
            tokenize("dose 1 mg\ndose 2 @ mg\n")
        d = exc.value.diagnostics[0]
        assert d.line == 2 and d.col == 8 and d.code == "unexpected-char"


# --------------------------------------------------------------------------- #
class TestParserErrors:
    def test_the_spec_example_parses_clean(self):
        prog, diags = parse(STUPP)
        assert diags.ok, diags.to_dicts()
        assert [d.name for d in prog.drugs] == ["TMZ"]
        assert prog.protocols[0].name == "stupp_adjuvant"
        assert prog.protocols[0].repeat.count == 5
        assert prog.protocols[0].cycles.count == 6
        assert prog.runs[0].section == "xenium_319"
        assert prog.runs[0].options.report == [
            'coverage', 'kill', 'mean_auc', 'penetration_length']
        assert set(prog.runs[0].options.report) == set(METRICS)

    def test_every_ast_node_is_positioned(self):
        prog, _ = parse(STUPP)
        p = prog.protocols[0]
        for node in [prog.drugs[0], *prog.drugs[0].fields, p, *p.statements(),
                     prog.runs[0]]:
            assert node.line > 0, node
            assert node.col > 0, node
            assert node.end_col > node.col, node

    def test_dose_line_and_column_are_the_dose_statement(self):
        prog, _ = parse(STUPP)
        dose = prog.protocols[0].dose
        assert dose.line == 7 and dose.col == 3
        assert dose.amount.unit == "mg/m2" and dose.amount.value == 150.0

    @pytest.mark.parametrize("bad,code,line", [
        (STUPP.replace("dose 150 mg/m2", "dose 150"), "missing-unit", 7),
        (STUPP.replace("dose 150 mg/m2", "dose 150 mg/m3"), "unknown-unit", 7),
        (STUPP.replace("every 24 h", "every 24 um"), "wrong-dimension", 8),
        (STUPP.replace("repeat 5", "repeat 5 h"), "unexpected-unit", 9),
        (STUPP.replace("route oral", "route inhaled"), "unknown-route", 6),
        (STUPP.replace("  dx 30 um", "  dx 30 um\n  bogus 1"),
         "unknown-option", 17),
        (STUPP.replace("report coverage", "report bogus"), "unknown-metric", 19),
        (STUPP.replace("state_mode oxygen", "state_mode aura"),
         "unknown-state-mode", 17),
        (STUPP.replace("from smiles", "from inchi"), "unknown-drug-source", 1),
        (STUPP.replace("  ic50 from gdsc", "  potency 1 uM"), "unknown-field", 2),
        (STUPP.replace("\nend\n\nrun", "\n\nrun"), "unclosed-block", 15),
    ])
    def test_bad_syntax_is_one_positioned_error(self, bad, code, line):
        d = first_error(bad)
        assert d["code"] == code, d
        assert d["line"] == line, d
        assert d["col"] > 0
        assert d["message"] and not d["message"].startswith("Traceback")

    def test_a_second_statement_still_parses_after_a_bad_one(self):
        src = ('protocol broken\n  dose 5 furlongs\nend\n\n'
               'protocol fine\n  dose 5 mg\nend\n')
        prog, diags = parse(src)
        assert not diags.ok
        assert [p.name for p in prog.protocols] == ["fine"]

    def test_unknown_statement_names_the_valid_ones(self):
        d = first_error("import numpy\n")
        assert d["code"] == "unknown-statement"
        assert "protocol" in d["message"] and "sweep" in d["message"]

    def test_duplicate_field_and_duplicate_statement_are_refused(self):
        assert first_error(STUPP.replace(
            "  ic50 from gdsc", "  ic50 1 uM\n  ic50 2 uM"))[
                "code"] == "duplicate-field"
        assert first_error(STUPP.replace(
            "  repeat 5", "  repeat 5\n  repeat 6"))[
                "code"] == "duplicate-statement"

    def test_expect_is_parsed_and_never_blocks(self):
        prog, diags = parse(STUPP + "expect kill > 0.05\nexpect coverage >= 10 %\n")
        assert diags.ok
        assert [e.metric for e in prog.expects] == ["kill", "coverage"]
        assert prog.expects[0].comparator == ">"
        assert prog.expects[1].value.unit == "%"

    def test_compare_needs_two_drugs(self):
        src = STUPP + "compare TMZ with stupp_adjuvant on xenium_319\n  end\n"
        assert first_error(src)["code"] == "incomplete-compare"


# --------------------------------------------------------------------------- #
class TestUnits:
    def test_the_frozen_unit_table_is_exactly_the_spec_list(self):
        assert set(UNITS) == {
            "mg/m2", "mg/kg", "mg", "uM", "nM", "mM", "um", "mm", "s", "min",
            "h", "d", "um2/s", "%"}

    @pytest.mark.parametrize("value,src,dst,expected", [
        (24, "h", "s", 86400.0),
        (1, "d", "s", 86400.0),
        (90, "min", "s", 5400.0),
        (2.5, "mm", "um", 2500.0),
        (500, "nM", "uM", 0.5),
        (1.5, "mM", "uM", 1500.0),
    ])
    def test_conversions_are_exact(self, value, src, dst, expected):
        out, record = convert(value, src, dst)
        assert out == pytest.approx(expected)
        assert src in record and dst in record  # the arithmetic is reported

    def test_percent_canonicalises_to_a_fraction(self):
        out, record = Quantity(31.0, "%").canonical()
        assert out == pytest.approx(0.31)
        assert "fraction" in record

    def test_case_matters_because_mM_and_mm_are_different_dimensions(self):
        assert UNITS["mM"].dimension == "concentration"
        assert UNITS["mm"].dimension == "length"
        with pytest.raises(TplError):
            convert(1.0, "mm", "uM")

    def test_cross_dimension_conversion_is_refused_not_coerced(self):
        with pytest.raises(TplError) as exc:
            convert(1.0, "h", "um")
        assert exc.value.diagnostics[0].code == "wrong-dimension"

    def test_dimensionless_quantity_cannot_be_converted(self):
        with pytest.raises(TplError) as exc:
            Quantity(1.0, None).to("s")
        assert exc.value.diagnostics[0].code == "missing-unit"


# --------------------------------------------------------------------------- #
class TestDoseToPk:
    def _tmz(self, **kw):
        args = dict(dose=Quantity(150.0, "mg/m2"), route="oral",
                    every_s=86400.0, repeat=5, rest_s=23 * 86400.0, cycles=6,
                    molecular_weight=194.154, drug_key="temozolomide",
                    max_window_s=86400.0)
        args.update(kw)
        return convert_dose(**args)

    def test_schedule_expansion_matches_a_28_day_stupp_cycle(self):
        times, duration = schedule_times(86400.0, 5, 23 * 86400.0, 6)
        assert len(times) == 30
        assert times[:5] == [0.0, 86400.0, 172800.0, 259200.0, 345600.0]
        assert times[5] == 28 * 86400.0        # cycle length = 5 d + 23 d
        assert duration == 6 * 28 * 86400.0

    def test_mg_per_m2_uses_the_stated_body_surface_area(self):
        c = self._tmz()
        assert c.amount_mg == pytest.approx(150.0 * BSA_M2)
        assert any(a.quantity == "body_surface_area" and a.value == BSA_M2
                   and a.provenance == "assumed" for a in c.assumptions)

    def test_mg_per_kg_uses_the_stated_body_weight(self):
        c = self._tmz(dose=Quantity(2.0, "mg/kg"))
        assert c.amount_mg == pytest.approx(2.0 * BODY_WEIGHT_KG)
        assert any(a.quantity == "body_weight" for a in c.assumptions)

    def test_the_peak_is_the_one_compartment_bolus_relation(self):
        c = self._tmz()
        expected = (1000.0 * c.amount_mg * c.bioavailability
                    / (c.volume_of_distribution_L * 194.154))
        assert c.plasma_peak_uM == pytest.approx(expected)
        assert "1000" in c.formula and "Vd_L" in c.formula
        assert len(c.conversions) >= 2

    def test_every_assumption_is_named_sourced_and_tagged(self):
        c = self._tmz()
        names = {a.quantity for a in c.assumptions}
        assert names == {"body_surface_area", "bioavailability",
                         "volume_of_distribution", "plasma_half_life"}
        for a in c.assumptions:
            assert a.provenance in ("measured", "literature", "computed",
                                    "assumed", "flagged-default")
            assert len(a.source) > 20, a          # a real sentence, not a word

    def test_temozolomide_pk_is_literature_with_a_citation(self):
        c = self._tmz()
        assert c.half_life_s == pytest.approx(1.8 * 3600.0)
        hl = next(a for a in c.assumptions if a.quantity == "plasma_half_life")
        assert hl.provenance == "literature"
        assert "PMID" in hl.source
        assert c.volume_of_distribution_L == PK_LITERATURE[
            "temozolomide"]["vd_L"]

    def test_an_uncurated_drug_falls_to_flagged_and_assumed_defaults(self):
        c = self._tmz(drug_key="novel_compound_47")
        hl = next(a for a in c.assumptions if a.quantity == "plasma_half_life")
        vd = next(a for a in c.assumptions
                  if a.quantity == "volume_of_distribution")
        assert hl.provenance == "flagged-default"
        assert vd.provenance == "assumed"
        assert "no curated" in vd.source

    def test_an_iv_route_is_bioavailability_one_by_definition(self):
        c = self._tmz(route="iv")
        f = next(a for a in c.assumptions if a.quantity == "bioavailability")
        assert f.value == 1.0 and f.provenance == "computed"
        assert "by definition" in f.source

    def test_oral_is_refused_for_a_drug_that_is_not_given_orally(self):
        with pytest.raises(ValueError, match="not given orally"):
            convert_dose(dose=Quantity(75.0, "mg/m2"), route="oral",
                         every_s=86400.0, repeat=1, rest_s=0.0, cycles=1,
                         molecular_weight=298.0, drug_key="cisplatin")

    def test_a_missing_molecular_weight_is_refused_not_defaulted(self):
        with pytest.raises(ValueError, match="molecular weight"):
            self._tmz(molecular_weight=0.0)

    def test_time_resolved_is_always_set_and_explained(self):
        c = self._tmz()
        assert c.time_resolved is True
        assert "time_resolved" in c.time_resolved_reason
        assert "dose_schedule" in c.time_resolved_reason

    def test_the_window_truncation_is_reported_with_its_dose_count(self):
        c = self._tmz()
        assert c.protocol_duration_s == 6 * 28 * 86400.0
        assert c.window_s == 86400.0
        assert c.truncated is True
        assert c.n_administrations_total == 30
        assert c.n_administrations_outside_window == 29
        assert any("outside" in n.lower() for n in c.notes)

    def test_a_bolus_is_one_dose_at_the_stated_interval(self):
        c = self._tmz(rest_s=0.0, cycles=1, max_window_s=5 * 86400.0)
        assert c.infusion is None
        assert c.n_doses == 5
        assert c.dose_interval_s == 86400.0

    def test_over_produces_a_sustained_profile_not_a_bolus(self):
        c = self._tmz(repeat=1, rest_s=0.0, cycles=1,
                      over=Quantity(24.0, "h"), max_window_s=86400.0)
        assert c.infusion is not None
        assert c.infusion["mode"].startswith("sub-dose expansion")
        k = c.infusion["n_subdoses_per_administration"]
        assert k > 1
        assert c.n_doses == k                       # drug ARRIVES over T
        assert c.dose_interval_s == pytest.approx(86400.0 / k)
        bolus = self._tmz(repeat=1, rest_s=0.0, cycles=1)
        assert c.plasma_peak_uM == pytest.approx(bolus.plasma_peak_uM / k)
        # Same total drug, spread out: peak x doses is conserved.
        assert (c.plasma_peak_uM * c.n_doses
                == pytest.approx(bolus.plasma_peak_uM * bolus.n_doses))

    def test_a_repeating_infusion_falls_back_to_a_labelled_approximation(self):
        c = self._tmz(repeat=3, rest_s=0.0, cycles=1,
                      over=Quantity(1.0, "h"), max_window_s=3 * 86400.0)
        assert c.infusion["mode"].startswith("peak attenuation")
        assert 0.0 < c.infusion["attenuation_factor"] < 1.0
        assert "APPROXIMATE" in c.infusion["mode"]
        assert "OVER-stated" in c.infusion["explanation"]
        # It reproduces the exact one-compartment infusion Cmax.
        lam = math.log(2.0) / c.half_life_s
        expected = (1.0 - math.exp(-lam * 3600.0)) / (lam * 3600.0)
        assert c.infusion["attenuation_factor"] == pytest.approx(expected)


# --------------------------------------------------------------------------- #
class TestSweepExpansion:
    def _n(self, sweep_line, extra=""):
        src = STUPP.replace("run TMZ", sweep_line + "\n" + extra + "run TMZ")
        v = validate_source(src)
        return v

    def test_step_range_is_inclusive_of_both_ends(self):
        v = self._n("sweep dose from 50 to 400 mg/m2 step 50")
        assert v["ok"], v["diagnostics"]
        assert v["n_runs"] == 8

    def test_in_n_steps_is_linear(self):
        p = plan(STUPP.replace("run TMZ",
                               "sweep every from 6 to 48 h in 4 steps\nrun TMZ"),
                 compile_drugs=False)
        axis = p.sweep_axes[0]
        assert axis["n_points"] == 4
        assert axis["values"] == ["6 h", "20 h", "34 h", "48 h"]

    def test_a_fixed_sweep_is_one_point_and_an_override(self):
        p = plan(STUPP.replace("run TMZ", "sweep dx 12 um\nrun TMZ"),
                 compile_drugs=False)
        assert p.sweep_axes[0]["kind"] == "fixed"
        assert p.n_runs == 1
        assert p.runs[0].dx_um == 12.0        # the sweep overrides `dx 30 um`

    def test_two_axes_are_a_cross_product(self):
        v = self._n("sweep dose from 50 to 200 mg/m2 step 50",
                    "sweep every from 12 to 24 h in 2 steps\n")
        assert v["n_runs"] == 8               # 4 x 2

    def test_a_compare_multiplies_by_the_drug_count(self):
        src = (f'drug TMZ from smiles "{TMZ_SMILES}"\nend\n'
               'drug CIS from library "cisplatin"\nend\n'
               'protocol p\n  route iv\n  dose 75 mg/m2\nend\n'
               'sweep dose from 50 to 100 mg/m2 step 25\n'
               'compare TMZ, CIS with p on xenium_319\n  dx 30 um\n'
               '  rank by coverage\nend\n')
        p = plan(src, compile_drugs=False)
        assert p.ok, p.diagnostics.to_dicts()
        assert p.n_runs == 6                  # 3 dose points x 2 drugs

    def test_the_cross_product_shares_one_reconstruction(self):
        p = plan(STUPP.replace(
            "run TMZ", "sweep dose from 50 to 400 mg/m2 step 50\nrun TMZ"),
            compile_drugs=False)
        assert p.n_runs == 8
        assert p.estimated_cost["n_reconstructions"] == 1
        assert len({r.recon_key for r in p.runs}) == 1

    def test_a_dx_sweep_needs_one_reconstruction_per_dx(self):
        p = plan(STUPP.replace(
            "run TMZ", "sweep dx from 30 to 60 um in 2 steps\nrun TMZ"),
            compile_drugs=False)
        assert p.n_runs == 2
        assert p.estimated_cost["n_reconstructions"] == 2

    def test_over_budget_is_refused_with_a_line_number_and_the_fix(self):
        src = STUPP.replace(
            "run TMZ", "sweep dose from 10 to 1000 mg/m2 step 10\nrun TMZ")
        d = first_error(src)
        assert d["code"] == "over-budget"
        assert d["line"] == 15                 # the `sweep` line
        assert "100 runs" in d["message"]
        assert str(DEFAULT_MAX_RUNS) in d["message"]
        assert d["hint"] == "allow runs 100"

    def test_allow_runs_lifts_the_budget(self):
        src = STUPP.replace(
            "run TMZ",
            "sweep dose from 10 to 1000 mg/m2 step 10\nallow runs 100\nrun TMZ")
        v = validate_source(src)
        assert v["ok"], v["diagnostics"]
        assert v["n_runs"] == 100

    def test_allow_runs_below_the_count_is_still_refused(self):
        src = STUPP.replace(
            "run TMZ",
            "sweep dose from 10 to 1000 mg/m2 step 10\nallow runs 20\nrun TMZ")
        d = first_error(src)
        assert d["code"] == "over-allowance"
        assert "you allowed" in d["message"]

    def test_a_non_sweepable_parameter_is_refused(self):
        d = first_error(STUPP.replace("run TMZ", "sweep colour 1 um\nrun TMZ"))
        assert d["code"] == "unknown-sweep-param"
        assert "dose" in d["message"]

    def test_an_absurd_step_is_refused_before_it_is_materialised(self):
        d = first_error(STUPP.replace(
            "run TMZ", "sweep dose from 1 to 1000000 mg/m2 step 0.01\nrun TMZ"))
        assert d["code"] == "absurd-sweep"

    def test_a_backwards_step_is_refused(self):
        d = first_error(STUPP.replace(
            "run TMZ", "sweep dose from 400 to 50 mg/m2 step 50\nrun TMZ"))
        assert d["code"] == "empty-sweep"


# --------------------------------------------------------------------------- #
class TestCapability:
    def test_marker_state_mode_on_the_319_panel_is_an_error(self):
        d = first_error(STUPP.replace("state_mode oxygen", "state_mode markers"))
        assert d["code"] == "capability"
        assert d["line"] == 17
        assert "hypoxia" in d["message"]
        assert "cell_state_impact" in d["message"]

    def test_the_default_oxygen_route_is_allowed_with_a_caveat(self):
        v = validate_source(STUPP)
        assert v["ok"]
        warns = " ".join(d["message"] for d in v["diagnostics"]
                         if d["severity"] == "warning")
        assert "ABSOLUTE KILL MUST NOT BE QUOTED" in warns
        assert "UPPER BOUND" in warns
        assert "NOT FOR CLINICAL USE" in warns

    def test_an_unknown_specimen_is_refused_with_the_known_list(self):
        d = first_error(STUPP.replace("xenium_319", "not_a_section"))
        assert d["code"] == "unknown-section"
        assert "xenium_319" in d["message"]

    def test_determinant_sensitivity_is_refused_on_a_visium_section(self):
        src = STUPP.replace("xenium_319", "visium_gbm1").replace(
            "sensitivity uniform", "sensitivity determinant")
        codes = [d["code"] for d in errors(src)]
        # Either the section is not cached in this checkout (refused with a
        # reason) or the mode is refused; both are the required behaviour.
        assert codes and set(codes) <= {"capability", "section-unavailable"}

    def test_the_mode_tuples_match_the_server_registry(self):
        from app.server.cache import SENSITIVITY_MODES, STATE_MODES

        from tmesim.protocol.ast import (
            SENSITIVITY_MODES as TPL_SENS,
            STATE_MODES as TPL_STATE,
        )
        assert tuple(TPL_STATE) == tuple(STATE_MODES)
        assert tuple(TPL_SENS) == tuple(SENSITIVITY_MODES)


# --------------------------------------------------------------------------- #
class TestDrugSources:
    def test_ic50_from_gdsc_on_a_smiles_drug_is_measured(self):
        p = plan(STUPP)
        assert p.ok, p.diagnostics.to_dicts()
        d = p.drugs["TMZ"]
        assert d.identity == "temozolomide"
        assert "canonical-SMILES" in d.identity_route
        assert d.provenance["ic50"]["provenance"] == "measured"
        assert "GDSC" in d.provenance["ic50"]["source"]
        assert d.ic50_resolution["n_cell_lines"] > 0

    def test_ic50_from_gdsc_errors_when_the_drug_is_absent_from_the_screen(self):
        src = STUPP.replace(f'smiles "{TMZ_SMILES}"', 'library "doxorubicin"')
        d = first_error(src)
        assert d["code"] == "gdsc-absent"
        assert d["line"] == 2
        assert "literature" in d["message"]     # names what it refused to accept

    def test_ic50_from_gdsc_errors_when_no_identity_can_be_established(self):
        src = STUPP.replace(TMZ_SMILES, "CCCCCCCCCCCCCCCCCCBr")
        d = first_error(src)
        assert d["code"] == "gdsc-no-identity"

    def test_an_unknown_library_drug_is_refused_with_the_available_list(self):
        src = STUPP.replace(f'smiles "{TMZ_SMILES}"', 'library "aspirin"')
        d = first_error(src)
        assert d["code"] == "unknown-drug"
        assert "temozolomide" in d["message"]

    def test_an_unparseable_smiles_is_refused(self):
        d = first_error(STUPP.replace(TMZ_SMILES, "Xy!"))
        assert d["code"] in ("bad-smiles", "unexpected-char")

    def test_per_field_provenance_covers_every_drugparams_field(self):
        p = plan(STUPP)
        prov = p.drugs["TMZ"].provenance
        for fld in ("ic50", "hill", "diffusivity", "binding", "clearance",
                    "bbb_penetrance", "molecular_weight"):
            assert fld in prov, fld
            assert prov[fld]["provenance"] in (
                "measured", "computed", "user-supplied", "flagged-default")
            assert len(prov[fld]["source"]) > 15
        assert prov["hill"]["provenance"] == "flagged-default"
        assert prov["diffusivity"]["provenance"] == "computed"
        assert prov["molecular_weight"]["value"] == pytest.approx(194.154, abs=0.1)

    def test_an_explicit_override_wins_and_is_labelled_user_supplied(self):
        src = STUPP.replace("  ic50 from gdsc",
                            "  ic50 12 uM\n  bbb_penetrance 0.31")
        p = plan(src)
        assert p.ok, p.diagnostics.to_dicts()
        d = p.drugs["TMZ"]
        assert d.drug_params.ic50 == pytest.approx(12.0)
        assert d.drug_params.bbb_penetrance == pytest.approx(0.31)
        assert d.provenance["ic50"]["provenance"] == "user-supplied"
        assert d.provenance["bbb_penetrance"]["provenance"] == "user-supplied"

    def test_nM_ic50_is_converted_and_the_conversion_is_reported(self):
        p = plan(STUPP.replace("  ic50 from gdsc", "  ic50 500 nM"))
        d = p.drugs["TMZ"]
        assert d.drug_params.ic50 == pytest.approx(0.5)
        assert "nM" in d.provenance["ic50"]["source"]
        assert "uM" in d.provenance["ic50"]["source"]

    def test_gdsc_and_an_explicit_ic50_together_are_refused(self):
        src = STUPP.replace("  ic50 from gdsc", "  ic50 from gdsc\n  ic50 12 uM")
        assert first_error(src)["code"] == "conflicting-ic50"

    def test_bbb_penetrance_out_of_range_is_refused(self):
        src = STUPP.replace("  ic50 from gdsc", "  bbb_penetrance 1.4")
        assert first_error(src)["code"] == "out-of-range"

    def test_a_missing_upload_names_where_it_looked(self):
        src = STUPP.replace(f'smiles "{TMZ_SMILES}"',
                            'file "uploads/does_not_exist.sdf"')
        d = first_error(src)
        assert d["code"] == "file-not-found"
        assert "uploads/does_not_exist.sdf" in d["message"]

    def test_a_path_outside_the_repository_is_refused(self):
        src = STUPP.replace(f'smiles "{TMZ_SMILES}"', 'file "/etc/passwd"')
        assert first_error(src)["code"] == "file-not-found"

    @pytest.mark.parametrize("name,content", [
        ("tpl_test.smi", f"{TMZ_SMILES} temozolomide\n"),
        ("tpl_test.json", json.dumps({"smiles": TMZ_SMILES})),
    ])
    def test_a_structure_file_compiles_like_a_smiles_drug(self, name, content):
        from tmesim.protocol.resolve import uploads_dir

        path = os.path.join(uploads_dir(), name)
        with open(path, "w") as fh:
            fh.write(content)
        try:
            src = STUPP.replace(f'smiles "{TMZ_SMILES}"',
                                f'file "uploads/{name}"')
            p = plan(src)
            assert p.ok, p.diagnostics.to_dicts()
            assert p.drugs["TMZ"].identity == "temozolomide"
            assert p.drugs["TMZ"].provenance["ic50"]["provenance"] == "measured"
        finally:
            os.remove(path)

    def test_a_parameter_file_without_a_structure_must_supply_what_it_needs(self):
        from tmesim.protocol.resolve import uploads_dir

        path = os.path.join(uploads_dir(), "tpl_params_only.json")
        with open(path, "w") as fh:
            json.dump({"ic50": 3.0}, fh)
        try:
            src = STUPP.replace(f'smiles "{TMZ_SMILES}"',
                                'file "uploads/tpl_params_only.json"').replace(
                                    "  ic50 from gdsc", "")
            d = first_error(src)
            assert d["code"] == "incomplete-parameters"
            assert "diffusivity" in d["message"]
        finally:
            os.remove(path)


# --------------------------------------------------------------------------- #
class TestPlanAndProvenance:
    def test_the_plan_reports_the_cost_without_executing(self):
        p = plan(STUPP.replace(
            "run TMZ", "sweep dose from 100 to 250 mg/m2 step 50\nrun TMZ"))
        assert p.ok and p.n_runs == 4
        cost = p.estimated_cost
        assert cost["n_reconstructions"] == 1
        assert cost["total_pde_timesteps"] > 0
        assert cost["wall_time_s"] is None       # never invented
        assert "not predicted" in cost["wall_time_reason"]

    def test_each_run_carries_its_own_dose_conversion(self):
        p = plan(STUPP.replace(
            "run TMZ", "sweep dose from 100 to 200 mg/m2 step 50\nrun TMZ"))
        peaks = [r.dose_conversion.plasma_peak_uM for r in p.runs]
        assert len(peaks) == 3
        assert peaks == sorted(peaks)            # monotone in dose
        assert peaks[2] / peaks[0] == pytest.approx(2.0)

    def test_a_program_with_no_run_is_refused(self):
        d = first_error(STUPP.split("run TMZ")[0])
        assert d["code"] == "no-run"

    def test_a_protocol_without_a_dose_is_refused(self):
        src = STUPP.replace("  dose 150 mg/m2\n", "")
        assert first_error(src)["code"] == "incomplete-protocol"

    def test_repeat_without_every_is_refused(self):
        src = STUPP.replace("  every 24 h\n", "")
        d = first_error(src)
        assert d["code"] == "incomplete-protocol"
        assert "t = 0" in d["message"]


# --------------------------------------------------------------------------- #
class TestRoutes:
    @pytest.fixture(scope="class")
    def client(self):
        from fastapi.testclient import TestClient

        from app.server.api import create_app
        from app.server.jobs import JobManager

        manager = JobManager(root=tempfile.mkdtemp(prefix="tpl-runs-"),
                             start_worker=False)
        with TestClient(create_app(manager)) as c:
            yield c

    def test_validate_returns_positioned_diagnostics(self, client):
        r = client.post("/api/protocol/validate", json={"source": STUPP})
        assert r.status_code == 200 and r.json()["ok"] is True
        r = client.post("/api/protocol/validate",
                        json={"source": STUPP.replace("dose 150 mg/m2",
                                                      "dose 150")})
        body = r.json()
        assert r.status_code == 200 and body["ok"] is False
        bad = [d for d in body["diagnostics"] if d["severity"] == "error"][0]
        assert bad["line"] == 7 and bad["col"] > 0 and bad["code"] == "missing-unit"

    def test_validate_accepts_a_plain_text_body(self, client):
        r = client.post("/api/protocol/validate", content=STUPP.encode(),
                        headers={"Content-Type": "text/plain"})
        assert r.status_code == 200 and r.json()["ok"] is True

    def test_validate_never_500s_on_garbage(self, client):
        for body in ("", "{{{", "drug\n", "\x00"):
            r = client.post("/api/protocol/validate", json={"source": body})
            assert r.status_code == 200
            assert "diagnostics" in r.json()

    def test_plan_expands_without_executing(self, client):
        src = STUPP.replace("run TMZ",
                            "sweep dose from 100 to 250 mg/m2 step 50\nrun TMZ")
        r = client.post("/api/protocol/plan", json={"source": src})
        body = r.json()
        assert r.status_code == 200 and body["ok"] is True
        assert body["n_runs"] == 4 and len(body["runs"]) == 4
        assert body["estimated_cost"]["n_reconstructions"] == 1
        assert body["drugs"]["TMZ"]["provenance"]["ic50"]["provenance"] == "measured"
        assert body["runs"][0]["dose_conversion"]["pk"]["plasma_peak"] > 0
        # nothing ran
        assert client.get("/api/runs").json() == []

    def test_run_refuses_a_bad_script_with_a_400_and_positions(self, client):
        r = client.post("/api/protocol/run",
                        json={"source": STUPP.replace("dose 150 mg/m2",
                                                      "dose 150")})
        assert r.status_code == 400
        diags = r.json()["detail"]["diagnostics"]
        assert diags[0]["line"] == 7 and diags[0]["code"] == "missing-unit"

    def test_run_refuses_an_over_budget_sweep(self, client):
        src = STUPP.replace(
            "run TMZ", "sweep dose from 10 to 1000 mg/m2 step 10\nrun TMZ")
        r = client.post("/api/protocol/run", json={"source": src})
        assert r.status_code == 400
        codes = [d["code"] for d in r.json()["detail"]["diagnostics"]]
        assert "over-budget" in codes

    def test_run_enqueues_on_the_shared_queue(self, client):
        r = client.post("/api/protocol/run",
                        json={"source": STUPP, "max_window_s": 600})
        assert r.status_code == 202
        body = r.json()
        assert body["run_id"].startswith("r_")
        assert body["n_runs"] == 1
        assert body["warnings"]                   # the standing caveats travel
        # It is a job on the SAME manager the /api/runs path uses.
        assert any(j["run_id"] == body["run_id"]
                   for j in client.get("/api/runs").json())
        got = client.get(f"/api/protocol/runs/{body['run_id']}").json()
        assert got["status"] in ("queued", "running", "done")
        assert got["source"] == STUPP

    def test_a_missing_protocol_run_is_a_404(self, client):
        assert client.get("/api/protocol/runs/r_99999").status_code == 404

    def test_a_runspec_run_is_not_served_as_a_protocol_run(self, client):
        spec = client.get("/api/defaults").json()
        rid = client.post("/api/runs", json=spec).json()["run_id"]
        r = client.get(f"/api/protocol/runs/{rid}")
        assert r.status_code == 409
        assert "not a protocol run" in r.json()["detail"]["message"]

    def test_upload_stores_a_structure_and_previews_it(self, client):
        r = client.post("/api/protocol/upload",
                        json={"filename": "tpl_route_test.smi",
                              "content": f"{TMZ_SMILES} temozolomide\n"})
        assert r.status_code == 200
        body = r.json()
        assert body["path"] == "uploads/tpl_route_test.smi"
        assert body["preview"]["kind"] == "structure"
        assert body["preview"]["descriptors"]["mw"] == pytest.approx(194.15,
                                                                     abs=0.1)
        assert body["use_in_script"].startswith("drug ")
        # the returned path is usable verbatim
        src = STUPP.replace(f'smiles "{TMZ_SMILES}"', f'file "{body["path"]}"')
        assert client.post("/api/protocol/validate",
                           json={"source": src}).json()["ok"] is True

    def test_upload_accepts_base64(self, client):
        import base64

        r = client.post("/api/protocol/upload", json={
            "filename": "tpl_b64.json",
            "content_base64": base64.b64encode(
                json.dumps({"smiles": TMZ_SMILES}).encode()).decode()})
        assert r.status_code == 200
        assert r.json()["preview"]["kind"] == "parameters"

    def test_upload_refuses_an_unsupported_extension(self, client):
        r = client.post("/api/protocol/upload",
                        json={"filename": "evil.py", "content": "import os"})
        assert r.status_code == 400
        assert ".smi" in str(r.json()["detail"]["allowed"])

    def test_upload_refuses_an_unparseable_file_and_does_not_keep_it(self, client):
        from tmesim.protocol.resolve import uploads_dir

        r = client.post("/api/protocol/upload",
                        json={"filename": "tpl_bad.smi", "content": "Xy!!zz\n"})
        assert r.status_code == 400
        assert not os.path.exists(os.path.join(uploads_dir(), "tpl_bad.smi"))

    def test_upload_path_traversal_is_neutralised(self, client):
        r = client.post("/api/protocol/upload",
                        json={"filename": "../../../../etc/passwd.smi",
                              "content": f"{TMZ_SMILES}\n"})
        assert r.status_code == 200
        assert r.json()["path"].startswith("uploads/")
        assert ".." not in r.json()["path"]


# --------------------------------------------------------------------------- #
class TestExecution:
    """One real end-to-end sweep. Short window, coarse grid, 2 points.

    Deliberately cheap (dx = 60 um, a 5-10 minute window) so the suite stays
    fast; the full-resolution end-to-end verification is a script, not a test.
    """

    def test_a_dose_sweep_reconstructs_once_and_returns_a_row_per_point(self):
        from tmesim.protocol import execute_plan

        src = STUPP.replace(
            "run TMZ", "sweep dose from 100 to 200 mg/m2 step 100\nrun TMZ"
        ).replace("dx 30 um", "dx 60 um")
        p = plan(src, max_window_s=600.0)
        assert p.ok, p.diagnostics.to_dicts()
        assert p.n_runs == 2

        result = execute_plan(p)
        assert result["reconstructions"]["n_built"] == 1
        assert result["reconstructions"]["n_cache_hits"] == 1
        assert [c["status"] for c in result["grid"]] == ["done", "done"]

        rows = result["rows"]
        assert len(rows) == 2
        for row in rows:
            assert row["metrics"]["kill_fraction"] is not None
            assert row["provenance"]["dose_to_pk"]["time_resolved"] is True
            assert row["provenance"]["solver"]["pde_time_resolved"] is True
            assert row["provenance"]["drug"]["provenance"]["ic50"][
                "provenance"] == "measured"
            assert any("NOT FOR CLINICAL USE" in c
                       for c in row["provenance"]["caveats"])
        # more drug -> more exposure
        assert (rows[1]["metrics"]["mean_tumor_auc"]
                > rows[0]["metrics"]["mean_tumor_auc"])

    def test_expectations_are_scored_and_never_block(self):
        from tmesim.protocol import execute_plan

        src = (STUPP.replace("dx 30 um", "dx 60 um")
               + "expect kill > 0.99\nexpect coverage >= 0\n")
        p = plan(src, max_window_s=300.0)
        result = execute_plan(p)
        assert result["grid"][0]["status"] == "done"
        exps = result["rows"][0]["expectations"]
        assert [e["passed"] for e in exps] == [False, True]
