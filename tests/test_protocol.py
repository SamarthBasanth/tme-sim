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
from tmesim.protocol.executor import DEFAULT_MAX_RUNS as _DEFAULT_MAX_RUNS  # noqa: F401
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

    @pytest.mark.parametrize("bad,code,marker", [
        (STUPP.replace("dose 150 mg/m2", "dose 150"), "missing-unit",
         "dose 150"),
        (STUPP.replace("dose 150 mg/m2", "dose 150 mg/m3"), "unknown-unit",
         "mg/m3"),
        (STUPP.replace("every 24 h", "every 24 um"), "wrong-dimension",
         "every 24 um"),
        (STUPP.replace("repeat 5", "repeat 5 h"), "unexpected-unit",
         "repeat 5 h"),
        (STUPP.replace("route oral", "route inhaled"), "unknown-route",
         "inhaled"),
        (STUPP.replace("  dx 30 um", "  dx 30 um\n  bogus 1"),
         "unknown-option", "bogus 1"),
        (STUPP.replace("report coverage", "report bogus"), "unknown-metric",
         "report bogus"),
        (STUPP.replace("state_mode oxygen", "state_mode aura"),
         "unknown-state-mode", "state_mode aura"),
        (STUPP.replace("from smiles", "from inchi"), "unknown-drug-source",
         "from inchi"),
        (STUPP.replace("  ic50 from gdsc", "  potency 1 uM"), "unknown-field",
         "potency 1 uM"),
    ])
    def test_bad_syntax_is_one_positioned_error(self, bad, code, marker):
        expected_line = next(i for i, ln in enumerate(bad.split("\n"), 1)
                             if marker in ln)
        d = first_error(bad)
        assert d["code"] == code, d
        assert d["line"] == expected_line, d
        assert d["col"] > 0
        assert d["end_col"] > d["col"]
        assert d["message"] and not d["message"].startswith("Traceback")

    def test_an_unclosed_block_is_reported(self):
        d = first_error(STUPP.replace("  cycles 6\nend\n", "  cycles 6\n"))
        assert d["code"] == "unclosed-block"
        assert d["line"] > 0

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

    def test_compare_needs_two_runs(self):
        # WAS `test_compare_needs_two_drugs`. The refusal survives with its code
        # intact, but it now counts RUNS, not drugs: `compare` is two or more
        # runs of ONE compound, so a compare with no sweep plans one run and one
        # run is not a comparison.
        src = STUPP + "compare TMZ with stupp_adjuvant on xenium_319\n  end\n"
        d = first_error(src)
        assert d["code"] == "incomplete-compare"
        assert "only one run" in d["message"]
        assert "sweep" in d["message"]

    def test_compare_refuses_two_distinct_compounds_and_cites_the_evidence(self):
        src = (STUPP + 'drug CIS from library "cisplatin"\nend\n'
               "compare TMZ, CIS with stupp_adjuvant on xenium_319\n  end\n")
        d = first_error(src)
        assert d["code"] == "withdrawn-ranking"
        assert "WITHDRAWN" in d["message"]
        # The refusal must carry its evidence, not just its verdict.
        assert "results/pk_ranking_robustness.md" in d["message"]
        assert "results/capacity_ranking_robustness.md" in d["message"]

    def test_compare_refuses_one_compound_named_twice(self):
        src = (STUPP + "sweep dose from 50 to 100 mg/m2 step 25\n"
               "compare TMZ, TMZ with stupp_adjuvant on xenium_319\n  end\n")
        d = first_error(src)
        assert d["code"] == "incomplete-compare"
        assert "named 2 times" in d["message"]

    def test_rank_by_is_refused_at_validation_and_names_the_withdrawal(self):
        src = STUPP + ("sweep dose from 50 to 100 mg/m2 step 25\n"
                       "compare TMZ with stupp_adjuvant on xenium_319\n"
                       "  dx 30 um\n  rank by coverage\nend\n")
        errs = errors(src)
        d = next(e for e in errs if e["code"] == "withdrawn-ranking")
        assert "rank by" in d["message"] and "WITHDRAWN" in d["message"]
        assert "results/pk_ranking_robustness.md" in d["message"]
        assert "results/capacity_ranking_robustness.md" in d["message"]
        # Positioned, like every other refusal.
        assert d["line"] == next(i for i, ln in enumerate(src.split("\n"), 1)
                                 if ln.strip().startswith("rank by"))

    def test_rank_by_is_not_a_silent_no_op(self):
        src = STUPP + ("sweep dose from 50 to 100 mg/m2 step 25\n"
                       "compare TMZ with stupp_adjuvant on xenium_319\n"
                       "  dx 30 um\n  rank by coverage\nend\n")
        assert validate_source(src)["ok"] is False

    def test_run_options_carry_no_ranking_field(self):
        # The AST field is gone, so nothing downstream can resurrect an order.
        prog, diags = parse(STUPP)
        assert diags.ok
        opts = prog.runs[0].options.to_dict()
        assert "rank_by" not in opts
        assert set(opts) >= {"dx", "state_mode", "sensitivity", "report"}


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

    def test_the_stated_patient_constants_are_the_ones_in_force(self):
        """The operator-facing sentence used to HARDCODE "1.73 m2" while the
        arithmetic used ``BSA_M2``. Planting ``BSA_M2 = 1.99`` left all 132
        tests green -- the only test that pinned "1.73 m2" was matching the
        prose, not the value -- so the run could have been charged at one body
        size and described at another.
        """
        import tmesim.protocol.dosing as D

        # The LITERAL, not `== BSA_M2`. 1.73 m2 / 70 kg are documented values,
        # not free parameters: API_CONTRACT.md:287 publishes
        # "150 mg/m2 x 1.73 m2 = 259.5 mg" as the conversion string this
        # backend returns. Every assertion phrased against the imported
        # constant moves with it and cannot fail -- that is exactly how
        # BSA_M2 = 1.99 stayed green across all 132 tests.
        assert BSA_M2 == 1.73
        assert BODY_WEIGHT_KG == 70.0
        # The module docstring is a THIRD copy, and no f-string can interpolate
        # a docstring, so it is pinned rather than generated.
        assert f"BSA = {BSA_M2:g} m2" in D.__doc__
        assert f"weight = {BODY_WEIGHT_KG:g} kg" in D.__doc__

        c = self._tmz()
        bsa = next(a for a in c.assumptions
                   if a.quantity == "body_surface_area")
        assert bsa.value == BSA_M2
        assert f"{BSA_M2:g} m2" in bsa.source
        assert f"{BSA_M2}" in " ".join(c.conversions)
        w = self._tmz(dose=Quantity(2.0, "mg/kg"))
        bw = next(a for a in w.assumptions if a.quantity == "body_weight")
        assert bw.value == BODY_WEIGHT_KG
        assert f"{BODY_WEIGHT_KG:g} kg" in bw.source

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

    def test_a_compare_over_one_compound_plans_the_whole_sweep(self):
        # WAS `test_a_compare_multiplies_by_the_drug_count`, which planned
        # 3 dose points x 2 drugs = 6. The drug factor is gone: a compare names
        # one compound and the runs come from the sweep alone.
        src = (f'drug TMZ from smiles "{TMZ_SMILES}"\nend\n'
               'protocol p\n  route iv\n  dose 75 mg/m2\nend\n'
               'sweep dose from 50 to 100 mg/m2 step 25\n'
               'compare TMZ with p on xenium_319\n  dx 30 um\nend\n')
        p = plan(src, compile_drugs=False)
        assert p.ok, p.diagnostics.to_dicts()
        assert p.n_runs == 3                  # 3 dose points x 1 compound
        assert {r.drug_alias for r in p.runs} == {"TMZ"}

    def test_a_compare_over_two_compounds_plans_nothing(self):
        # The other half of the old multiply-by-drug-count test: the two-drug
        # program that used to plan 6 runs now plans 0 and refuses.
        src = (f'drug TMZ from smiles "{TMZ_SMILES}"\nend\n'
               'drug CIS from library "cisplatin"\nend\n'
               'protocol p\n  route iv\n  dose 75 mg/m2\nend\n'
               'sweep dose from 50 to 100 mg/m2 step 25\n'
               'compare TMZ, CIS with p on xenium_319\n  dx 30 um\nend\n')
        p = plan(src, compile_drugs=False)
        assert not p.ok
        assert p.n_runs == 0
        codes = [d["code"] for d in p.diagnostics.to_dicts()]
        assert "withdrawn-ranking" in codes

    def test_two_compounds_in_one_program_plan_nothing(self):
        """The `compare A, B` guard is scoped to ONE STATEMENT. On 2026-09-03
        four rewrites were measured walking around it, every one planning
        cleanly with ok = True and putting two compounds' metrics in ONE grid:

          two `run` blocks (2 runs) | `compare` + `run` (6) |
          two `compare` blocks (6)  | `sweep` x two `run` blocks (8)

        A plan is one job with one `grid`, one `rows` list and one sweep
        snapshot, and the editor paints that snapshot on a single colour scale
        spanning every cell (app/web/js/ui/protocol.js:793-801) — so two
        compounds in one plan IS the withdrawn ordering, however it was
        written. The rule is therefore on the plan.
        """
        head = (f'drug TMZ from smiles "{TMZ_SMILES}"\nend\n'
                'drug CIS from library "cisplatin"\nend\n'
                'protocol p\n  route iv\n  dose 75 mg/m2\nend\n')
        tails = {
            "two run blocks":
                'run TMZ with p on xenium_319\n  dx 30 um\nend\n'
                'run CIS with p on xenium_319\n  dx 30 um\nend\n',
            "compare + run":
                'sweep dose from 50 to 100 mg/m2 step 25\n'
                'compare TMZ with p on xenium_319\n  dx 30 um\nend\n'
                'run CIS with p on xenium_319\n  dx 30 um\nend\n',
            "two compare blocks":
                'sweep dose from 50 to 100 mg/m2 step 25\n'
                'compare TMZ with p on xenium_319\n  dx 30 um\nend\n'
                'compare CIS with p on xenium_319\n  dx 30 um\nend\n',
            "sweep x two run blocks":
                'sweep dose from 50 to 200 mg/m2 step 50\n'
                'run TMZ with p on xenium_319\n  dx 30 um\nend\n'
                'run CIS with p on xenium_319\n  dx 30 um\nend\n',
        }
        for label, tail in tails.items():
            p = plan(head + tail, compile_drugs=False)
            assert not p.ok, label
            assert p.n_runs == 0, label
            assert {r.drug_alias for r in p.runs} == set(), label
            d = next((x for x in p.diagnostics.to_dicts()
                      if x["code"] == "withdrawn-ranking"), None)
            assert d is not None, (label, p.diagnostics.to_dicts())
            assert "TMZ" in d["message"] and "CIS" in d["message"], label
            assert "results/pk_ranking_robustness.md" in d["message"], label
            assert "results/capacity_ranking_robustness.md" in d["message"], label
            assert d["line"] > 0 and d["col"] > 0, label

    def test_one_compound_across_several_run_blocks_is_still_allowed(self):
        # The rule is one COMPOUND, not one run statement. Two runs of the same
        # compound at different dx are a parameter comparison, which is the
        # product; this is what must NOT be caught by the guard above.
        src = (f'drug TMZ from smiles "{TMZ_SMILES}"\nend\n'
               'drug CIS from library "cisplatin"\nend\n'
               'protocol p\n  route iv\n  dose 75 mg/m2\nend\n'
               'run TMZ with p on xenium_319\n  dx 30 um\nend\n'
               'run TMZ with p on xenium_319\n  dx 60 um\nend\n')
        p = plan(src, compile_drugs=False)
        assert p.ok, p.diagnostics.to_dicts()
        assert p.n_runs == 2
        assert {r.drug_alias for r in p.runs} == {"TMZ"}
        # CIS is defined but never run: defining a second compound is not the
        # violation, putting it in the grid is.
        assert not [d for d in p.diagnostics.errors
                    if d.code == "withdrawn-ranking"]

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
        assert d["line"] == next(i for i, ln in enumerate(src.split("\n"), 1)
                                 if ln.startswith("sweep dose"))
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


class TestWindowTruncation:
    """``tmesim/protocol/dosing.py`` says the window truncation "is reported
    (with how many doses fell outside it) rather than hidden". These trace that
    claim from the arithmetic to every route that can be asked about the run.

    It used to stop short of two of them: the sentence lived only on a
    ``DoseConversion``, and ``POST /api/protocol/run`` plans with
    ``compile_drugs=False`` and so never builds one. An operator who went
    validate -> run was told nothing about 29 of 30 doses being dropped.
    """

    def test_the_planner_states_it_without_compiling_a_drug(self):
        p = plan(STUPP, compile_drugs=False)
        assert p.ok, p.diagnostics.to_dicts()
        assert p.drugs["TMZ"].drug_params is None     # nothing was compiled
        # The record rides on the RUN, not on the DoseConversion: the window
        # cuts the schedule, and the schedule does not depend on the compound.
        assert p.runs[0].window["truncated"] is True
        assert p.runs[0].window["n_administrations_outside_window"] == 29
        t = p.truncation
        assert t["truncated"] is True
        assert t["n_runs"] == 1 and t["n_runs_truncated"] == 1
        assert t["runs"][0]["n_administrations_total"] == 30
        assert t["runs"][0]["n_administrations_outside_window"] == 29
        assert "29 of 30 administrations fall OUTSIDE" in t["message"]

    def test_it_is_the_same_sentence_the_conversion_puts_on_the_result(self):
        # One arithmetic, one sentence. Two copies could drift, and then the
        # operator is shown one window and the solver runs another.
        planner = plan(STUPP, compile_drugs=False).truncation["message"]
        compiled = plan(STUPP, compile_drugs=True).runs[0].dose_conversion
        assert planner in compiled.notes
        assert compiled.n_administrations_outside_window == 29

    def test_a_schedule_that_fits_the_window_is_not_flagged(self):
        p = plan(STUPP, compile_drugs=False, max_window_s=200 * 86400.0)
        assert p.truncation["truncated"] is False
        assert p.truncation["n_runs_truncated"] == 0
        assert p.truncation["message"] == ""

    def test_a_sweep_that_changes_the_schedule_counts_the_truncated_points(self):
        # `sweep cycles` varies the protocol duration, so some points are cut
        # and some are not: the aggregate has to count, not just flag.
        src = STUPP.replace("run TMZ",
                            "sweep cycles from 1 to 2 in 2 steps\nrun TMZ")
        p = plan(src, compile_drugs=False, max_window_s=40 * 86400.0)
        assert p.n_runs == 2
        flags = [r.window["truncated"] for r in p.runs]
        assert flags == [False, True], flags
        assert p.truncation["n_runs_truncated"] == 1
        assert p.truncation["worst_run_index"] == 1


# --------------------------------------------------------------------------- #
class TestCapability:
    def test_marker_state_mode_on_the_319_panel_is_an_error(self):
        src = STUPP.replace("state_mode oxygen", "state_mode markers")
        d = first_error(src)
        assert d["code"] == "capability"
        assert d["line"] == next(i for i, ln in enumerate(src.split("\n"), 1)
                                 if "state_mode markers" in ln)
        assert "hypoxia" in d["message"]
        assert "cell_state_impact" in d["message"]

    def test_a_program_naming_no_metric_does_not_open_on_absolute_kill(self):
        """`report` is optional, and what you get for saying nothing becomes
        the sweep grid's opening column: the editor takes ``metrics[0]``
        (app/web/js/ui/protocol.js:758), the editor's metric list is built from
        ``plan.runs[0].report``, and that defaults to ``list(METRICS)``.

        With `kill` first, a program that named no metric opened on a
        colour-scaled map of ABSOLUTE KILL — the one number the SAME plan's
        capability caveat says must never be quoted from this panel. The metric
        is not removed; it is no longer what you get for saying nothing.
        """
        from app.server.protocol_wire import _metric_rows

        src = STUPP.replace(
            "  report coverage, kill, mean_auc, penetration_length\n", "")
        p = plan(src, compile_drugs=False)
        assert p.ok, p.diagnostics.to_dicts()
        assert p.runs[0].report == list(METRICS)     # this IS the default
        assert p.runs[0].report[0] != "kill"
        assert "kill" in p.runs[0].report            # still available
        assert _metric_rows(p)[0]["key"] != "kill"
        # An explicit request still wins: the user named it.
        explicit = plan(src.replace("  dx 30 um", "  dx 30 um\n  report kill"),
                        compile_drugs=False)
        assert _metric_rows(explicit)[0]["key"] == "kill"
        # The caveat that makes the default matter is on the same plan.
        assert any("ABSOLUTE KILL MUST NOT BE QUOTED" in d.message
                   for d in p.diagnostics.warnings)

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
        assert codes
        # The old assertion was `set(codes) <= {"capability",
        # "section-unavailable"}`, which passes when the section is simply not
        # cached -- i.e. it went green without ever reaching the determinant
        # guard it is named for. A checkout that cannot reach the subject must
        # SKIP, not pass.
        if set(codes) == {"section-unavailable"}:
            pytest.skip("visium_gbm1 is not cached in this checkout, so the "
                        "determinant refusal cannot be reached here")
        assert "capability" in codes
        d = next(x for x in errors(src) if x["code"] == "capability")
        assert "per-cell expression" in d["message"]

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
        assert d["line"] == 2               # the `ic50 from gdsc` line
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
            assert "molecular_weight" in d["message"]
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

    def test_a_two_compound_program_is_refused_at_every_route(self, client):
        """No route may hand back a payload with two compounds' numbers in it.

        `compare TMZ, CIS` was already refused; this is the same program
        written as two `run` blocks crossed with a `sweep`, which planned
        cleanly as an 8-cell two-compound dose ladder until 2026-09-03.
        """
        src = (f'drug TMZ from smiles "{TMZ_SMILES}"\nend\n'
               'drug CIS from library "cisplatin"\nend\n'
               'protocol p\n  route iv\n  dose 75 mg/m2\nend\n'
               'sweep dose from 50 to 200 mg/m2 step 50\n'
               'run TMZ with p on xenium_319\n  dx 30 um\nend\n'
               'run CIS with p on xenium_319\n  dx 30 um\nend\n')
        n_before = len(client.get("/api/runs").json())

        v = client.post("/api/protocol/validate", json={"source": src}).json()
        assert v["ok"] is False
        assert any(e["code"] == "withdrawn-ranking" for e in v["errors"])

        pl = client.post("/api/protocol/plan", json={"source": src}).json()
        assert pl["ok"] is False
        assert pl["runs"] == [] and pl["run_count"] == 0 and pl["n_runs"] == 0
        assert pl["refusal"]["code"] == "withdrawn-ranking"

        r = client.post("/api/protocol/run", json={"source": src})
        assert r.status_code == 400
        detail = r.json()["detail"]
        assert "withdrawn-ranking" in [d["code"] for d in detail["diagnostics"]]
        # nothing was enqueued, so no grid exists to compare them in
        assert len(client.get("/api/runs").json()) == n_before

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
class TestEditorWire:
    """The second projection: `app/web/js/protocol/api.js`'s frozen wire shape.

    Both projections ride in one payload under disjoint keys (see
    ``app/server/protocol_wire.py``); these pin the mapping so a change to one
    cannot silently break the other client.
    """

    @pytest.fixture(scope="class")
    def client(self):
        from fastapi.testclient import TestClient

        from app.server.api import create_app
        from app.server.jobs import JobManager

        manager = JobManager(root=tempfile.mkdtemp(prefix="tpl-wire-"),
                             start_worker=False)
        with TestClient(create_app(manager)) as c:
            yield c

    def test_validate_carries_both_projections(self, client):
        body = client.post("/api/protocol/validate",
                           json={"source": STUPP}).json()
        assert body["ok"] is True
        assert body["n_runs"] == body["run_count"] == 1
        assert body["symbols"] == {"drugs": ["TMZ"],
                                   "protocols": ["stupp_adjuvant"],
                                   "sweeps": []}
        assert body["errors"] == []
        assert body["warnings"] and all("line" in w for w in body["warnings"])

    def test_the_two_column_conventions_differ_by_exactly_one(self, client):
        bad = STUPP.replace("dose 150 mg/m2", "dose 150")
        body = client.post("/api/protocol/validate", json={"source": bad}).json()
        contract = [d for d in body["diagnostics"] if d["severity"] == "error"][0]
        editor = body["errors"][0]
        assert contract["line"] == editor["line"]
        assert editor["col"] == contract["col"] - 1        # 1-based -> 0-based
        assert editor["end_col"] == contract["end_col"] - 1
        assert editor["remedy"] == contract["hint"]

    def test_plan_axes_and_coords_index_the_grid(self, client):
        src = STUPP.replace("run TMZ",
                            "sweep dose from 100 to 200 mg/m2 step 50\nrun TMZ")
        body = client.post("/api/protocol/plan", json={"source": src}).json()
        assert body["run_count"] == 3
        assert body["axes"] == [{"name": "dose", "unit": "mg/m2",
                                 "values": [100.0, 150.0, 200.0]}]
        assert [r["coords"] for r in body["runs"]] == [[0], [1], [2]]
        assert [m["key"] for m in body["metrics"]] == [
            "coverage", "kill", "mean_auc", "penetration_length"]
        assert body["refusal"] is None
        assert body["estimated_cost_s"] > 0
        assert "ESTIMATE, not a measurement" in body["cost_note"]
        # the contract's own field still refuses to guess
        assert body["estimated_cost"]["wall_time_s"] is None

    def test_a_fixed_sweep_is_not_a_grid_dimension(self, client):
        src = STUPP.replace("run TMZ", "sweep dx 12 um\nrun TMZ")
        body = client.post("/api/protocol/plan", json={"source": src}).json()
        assert body["axes"] == []          # one point is not an axis
        assert body["run_count"] == 1

    def test_the_dose_to_pk_assumptions_are_renderable_sentences(self, client):
        body = client.post("/api/protocol/plan", json={"source": STUPP}).json()
        assumptions = body["assumptions"]
        assert assumptions and all(isinstance(a, str) for a in assumptions)
        joined = " ".join(assumptions)
        assert f"{BSA_M2:g} m2" in joined              # the BSA IN FORCE
        assert "PMID" in joined                        # a citation travels
        assert "time_resolved" in joined
        assert body["caveats"] and all(isinstance(c, str) for c in body["caveats"])

    def test_over_budget_becomes_a_refusal_with_insertable_tpl(self, client):
        src = STUPP.replace(
            "run TMZ", "sweep dose from 10 to 1000 mg/m2 step 10\nrun TMZ")
        body = client.post("/api/protocol/plan", json={"source": src}).json()
        ref = body["refusal"]
        assert body["ok"] is False
        assert ref["code"] == "max_runs_exceeded"
        assert ref["run_count"] == 100 and ref["max_runs"] == DEFAULT_MAX_RUNS
        assert ref["remedy_text"] == "allow runs 100"
        assert ref["remedy_line"] == ref["line"] > 0

    def test_the_truncation_reaches_plan_run_and_the_snapshot(self, client):
        """``tmesim/protocol/dosing.py`` says the truncation "is reported (with
        how many doses fell outside it) rather than hidden".

        It used to reach ``/plan``'s ``assumptions`` list and nothing else:
        ``/validate`` and ``/run`` plan with ``compile_drugs=False`` and so never
        built the ``DoseConversion`` that carried the sentence, and on the
        results panel the editor renders ``assumptions`` inside a COLLAPSED
        <details> (app/web/js/ui/protocol.js:809) while ``caveats`` sit beside
        the numbers (:790).
        """
        pl = client.post("/api/protocol/plan", json={"source": STUPP}).json()
        t = pl["truncation"]
        assert t["truncated"] is True
        assert t["n_runs"] == 1 and t["n_runs_truncated"] == 1
        assert "29 of 30 administrations fall OUTSIDE" in t["message"]
        assert any("fall OUTSIDE" in c for c in pl["caveats"])
        assert any("fall OUTSIDE" in a for a in pl["assumptions"])

        # /run is where the operator commits the compute.
        r = client.post("/api/protocol/run",
                        json={"source": STUPP, "max_window_s": 300})
        assert r.status_code == 202
        body = r.json()
        assert body["truncation"]["truncated"] is True
        assert body["truncation"]["n_runs_truncated"] == 1
        rid = body["run_id"]

        # ...and a poll of the STILL-QUEUED sweep says it too, before any
        # result document exists to carry it.
        got = client.get(f"/api/protocol/runs/{rid}").json()
        assert got["status"] in ("queued", "running", "done")
        assert got["truncation"]["truncated"] is True
        sweep = client.get(f"/api/protocol/sweeps/{rid}").json()
        assert sweep["truncation"]["truncated"] is True

    def test_a_rejected_plan_still_has_the_shape_of_a_plan(self, client):
        body = client.post("/api/protocol/plan",
                           json={"source": "drug\n"}).json()
        assert body["ok"] is False
        for key in ("runs", "axes", "metrics", "assumptions", "caveats",
                    "errors", "warnings", "diagnostics"):
            assert isinstance(body[key], list), key
        assert body["run_count"] == 0

    def test_the_rejected_plan_shape_holds_on_the_EXCEPTION_path_too(self):
        """`_diag_only(plan_shape=True)` is the fallback body for a plan that
        RAISED. The sibling test above never reaches it: `plan()` does not raise
        for malformed input, it returns a not-ok Plan, so every garbage body the
        suite posts ("", "{{{", "\\x00", an unterminated string) comes back
        through the NORMAL projection. Breaking `_diag_only`'s shape left the
        whole file green. Call it directly, since no route can.
        """
        from app.server.protocol_api import _diag_only

        diag = [{"severity": "error", "line": 1, "col": 1, "end_col": 2,
                 "code": "internal", "message": "boom", "hint": None}]
        out = _diag_only(diag, plan_shape=True)
        assert out["ok"] is False and out["run_count"] == 0
        for key in ("runs", "axes", "metrics", "assumptions", "caveats",
                    "errors", "warnings", "diagnostics", "expects"):
            assert isinstance(out[key], list), key
        assert isinstance(out["truncation"], dict)
        assert out["errors"] and out["errors"][0]["code"] == "internal"

    def test_confirm_run_count_mismatch_is_a_409(self, client):
        r = client.post("/api/protocol/run",
                        json={"source": STUPP, "confirm_run_count": 7})
        assert r.status_code == 409
        assert "confirmed 7 runs" in r.json()["detail"]["message"]
        assert r.json()["detail"]["run_count"] == 1

    def test_confirm_run_count_that_agrees_is_accepted(self, client):
        r = client.post("/api/protocol/run",
                        json={"source": STUPP, "confirm_run_count": 1,
                              "max_window_s": 300})
        assert r.status_code == 202
        body = r.json()
        assert body["sweep_id"] == body["run_id"]
        assert body["run_count"] == 1
        assert body["runs"][0]["run_id"] is None
        assert "one job" in body["runs"][0]["run_id_reason"].lower()

    def test_the_sweeps_alias_serves_the_same_job(self, client):
        rid = client.post("/api/protocol/run",
                          json={"source": STUPP, "max_window_s": 300}
                          ).json()["run_id"]
        sweep = client.get(f"/api/protocol/sweeps/{rid}").json()
        assert sweep["sweep_id"] == rid
        assert sweep["status"] in ("queued", "running", "done")
        assert set(sweep["counts"]) >= {"queued", "running", "done", "failed",
                                        "cancelled"}
        # /runs/{id} carries BOTH projections
        both = client.get(f"/api/protocol/runs/{rid}").json()
        assert both["run_id"] == both["sweep_id"] == rid
        assert "grid" in both and "runs" in both

    def test_cancelling_a_sweep_is_cooperative_and_says_so(self, client):
        rid = client.post("/api/protocol/run",
                          json={"source": STUPP, "max_window_s": 300}
                          ).json()["run_id"]
        out = client.delete(f"/api/protocol/sweeps/{rid}").json()
        assert out["cancelled"] is True
        assert "next run boundary" in out["reason"]
        assert client.delete("/api/protocol/sweeps/r_99999").status_code == 404

    def test_multipart_upload_is_parsed_without_python_multipart(self, client):
        import importlib

        with pytest.raises(ImportError):
            importlib.import_module("multipart")     # the dependency is absent
        body = (b'--BND\r\nContent-Disposition: form-data; name="file"; '
                b'filename="tpl_mp.smi"\r\n\r\n'
                + TMZ_SMILES.encode() + b" temozolomide\n\r\n--BND--\r\n")
        r = client.post("/api/protocol/upload", content=body,
                        headers={"content-type":
                                 "multipart/form-data; boundary=BND"})
        assert r.status_code == 200
        out = r.json()
        assert out["file"] == {"name": "tpl_mp.smi",
                               "size": len(TMZ_SMILES) + len(" temozolomide\n"),
                               "format": "smi"}
        fields = out["drug"]["fields"]
        assert set(fields) == {"ic50", "hill", "diffusivity", "binding",
                               "clearance", "bbb_penetrance",
                               "molecular_weight"}
        assert fields["ic50"]["provenance"] == "measured"
        assert fields["hill"]["provenance"] == "flagged-default"
        assert out["drug"]["identity"] == "temozolomide"
        assert out["caveats"]

    def test_a_multipart_body_with_no_file_part_is_a_400(self, client):
        r = client.post("/api/protocol/upload", content=b"--BND--\r\n",
                        headers={"content-type":
                                 "multipart/form-data; boundary=BND"})
        assert r.status_code == 400


class TestUploadBodyLimits:
    """H3. The parser buffers, splits, and copies; measured at ~5.7x the
    request size. The 2 MB upload cap used to be checked only AFTER all of
    that, so a ~500 MB POST cost ~2.8 GB before anything refused it, on an
    unauthenticated route with no rate limit. The cap has to bind before the
    allocation to be a cap."""

    @pytest.fixture(scope="class")
    def client(self):
        from fastapi.testclient import TestClient

        from app.server.api import create_app
        from app.server.jobs import JobManager

        manager = JobManager(root=tempfile.mkdtemp(prefix="tpl-cap-"),
                             start_worker=False)
        with TestClient(create_app(manager)) as c:
            yield c

    def test_an_oversized_body_is_refused_before_it_is_buffered(self, client):
        from app.server.protocol_api import _MAX_BODY_BYTES

        big = b"x" * (_MAX_BODY_BYTES + 1)
        r = client.post("/api/protocol/upload", content=big,
                        headers={"content-type": "text/plain"})
        assert r.status_code == 413
        assert "limit" in json.dumps(r.json())

    def test_the_cap_is_the_documented_size_and_not_merely_self_consistent(self):
        """Both oversized-body tests build their payload from
        ``_MAX_BODY_BYTES + 1``, so the cap's VALUE is untested: raising the
        constant to 8 MB leaves them green (measured), and so would 4 GB. This
        is the only body limit in the product (there is no reverse proxy in
        deployment mode 1a, protocol_api.py:68-74), so the number is pinned."""
        from app.server.protocol_api import _MAX_BODY_BYTES, _MAX_UPLOAD_BYTES

        assert _MAX_BODY_BYTES == 4 * 1024 * 1024
        assert _MAX_UPLOAD_BYTES == 2 * 1024 * 1024
        # the multipart envelope has to fit around a full-size payload
        assert _MAX_BODY_BYTES > _MAX_UPLOAD_BYTES

    def test_the_cap_binds_with_no_content_length_to_trust(self, client):
        """The Content-Length branch is the cheap path and it is CLIENT-SUPPLIED.
        A chunked body carries no length to refuse on, so only the counter
        inside the stream loop can stop it -- and that counter is the whole
        basis for the claim that the cap binds BEFORE the allocation. The
        sibling test above is satisfied by the header branch alone."""
        from app.server.protocol_api import _MAX_BODY_BYTES

        def chunks():
            for _ in range((_MAX_BODY_BYTES // 65536) + 2):
                yield b"x" * 65536

        r = client.post("/api/protocol/upload", content=chunks(),
                        headers={"content-type": "text/plain"})
        assert r.status_code == 413
        assert "limit" in json.dumps(r.json())

    def test_an_oversized_multipart_upload_is_refused(self, client):
        from app.server.protocol_api import _MAX_BODY_BYTES

        body = (b'--BND\r\nContent-Disposition: form-data; name="file"; '
                b'filename="huge.smi"\r\n\r\n'
                + b"C" * (_MAX_BODY_BYTES + 1) + b"\r\n--BND--\r\n")
        r = client.post("/api/protocol/upload", content=body,
                        headers={"content-type":
                                 "multipart/form-data; boundary=BND"})
        assert r.status_code == 413

    def test_an_empty_boundary_is_refused_rather_than_split_on(self, client):
        """`boundary=` made `sep` == b"--", splitting the body on every "--"
        in it: 4.2 million fragment objects out of 8 MB. RFC 2046 requires
        1-70 characters from a fixed charset."""
        body = b"-" * 4096
        r = client.post("/api/protocol/upload", content=body,
                        headers={"content-type":
                                 "multipart/form-data; boundary="})
        assert r.status_code == 400
        assert "boundary" in json.dumps(r.json())

    @pytest.mark.parametrize("boundary", ["a" * 71, "bad\\boundary", "b<ad>"])
    def test_a_malformed_boundary_is_refused(self, client, boundary):
        """The body must be WELL-FORMED under the bad boundary, or the 400 comes
        from "no file part" and _BOUNDARY_RE is never reached. It used to be
        b"--x--\r\n", which has no part at all: deleting the regex left this
        test green, which is the whole defect this file keeps growing."""
        sep = f"--{boundary}".encode()
        body = (sep + b'\r\nContent-Disposition: form-data; name="file"; '
                b'filename="tpl_bnd.smi"\r\n\r\n'
                + TMZ_SMILES.encode() + b"\r\n" + sep + b"--\r\n")
        r = client.post("/api/protocol/upload", content=body,
                        headers={"content-type":
                                 f"multipart/form-data; boundary={boundary}"})
        assert r.status_code == 400, r.text
        assert "boundary" in json.dumps(r.json())

    def test_a_legitimate_upload_is_unaffected(self, client):
        body = (b'--BND\r\nContent-Disposition: form-data; name="file"; '
                b'filename="tpl_cap.smi"\r\n\r\n'
                + TMZ_SMILES.encode() + b" temozolomide\n\r\n--BND--\r\n")
        r = client.post("/api/protocol/upload", content=body,
                        headers={"content-type":
                                 "multipart/form-data; boundary=BND"})
        assert r.status_code == 200
        assert r.json()["file"]["name"] == "tpl_cap.smi"


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
