"""Tests for ``vivum validate`` — the external-anchor scorecard.

Deliberately **fast and data-free**: not one test here touches ``data_cache/``
or runs a PDE. The real anchors are exercised by running the command; what is
tested here is the machinery that decides what a result MEANS:

1. **Argparse plumbing** — the subcommand is registered alongside the existing
   ones and every flag parses to the expected default.
2. **The PASS/FAIL rule in isolation** (:func:`scorecard.classify`), which is
   the only place a band comparison is spelled out.
3. **The PASS / FAIL / NOT TESTABLE / DEFERRED classification**, driven through
   :func:`scorecard.run_suite` with a synthetic registry. This is the part that
   matters most: the whole point of the command is that "unreachable" and
   "deferred" must never be silently converted into "passing", so there are
   explicit tests that a missing data file, a runner that declines, a runner
   that *crashes*, and a ``--quick`` deferral all produce a non-PASS status
   carrying a stated reason.
4. **Registry hygiene** — every registered anchor names its source module, its
   report, and only conditionals that actually exist.

The synthetic registry is injected by monkeypatching
``scorecard.build_registry``, so these tests stay valid as anchors are added or
removed.
"""
from __future__ import annotations

import json
import os
import sys

import pytest

import tmesim.cli as cli
from tmesim.validation import scorecard as SC


# --------------------------------------------------------------------------- #
# Synthetic registry helpers (no tissue, no PDE, no data_cache)
# --------------------------------------------------------------------------- #
def _row(group: str, status: str, **kw) -> SC.Row:
    kw.setdefault("scope", SC.SCOPE_ANALYTIC)
    kw.setdefault("quantity", "a model number vs a measured number")
    kw.setdefault("predicted", "1.0")
    kw.setdefault("measured", "1.0")
    kw.setdefault("pmids", ("12345678",))
    return SC.Row(group=group, label="probe", status=status, **kw)


def _group(key: str, runner, tier: str = SC.TIER_QUICK, needs=()) -> SC.AnchorGroup:
    return SC.AnchorGroup(
        key=key, title=f"synthetic {key}", tissue="test", tier=tier,
        scope=SC.SCOPE_ANALYTIC, source=f"tests::{key}", report="tests",
        runner=runner, needs=tuple(needs), est="0 s",
    )


def _install(monkeypatch, groups) -> None:
    monkeypatch.setattr(SC, "build_registry", lambda exchange=False: list(groups))


# --------------------------------------------------------------------------- #
# 1. Argparse plumbing
# --------------------------------------------------------------------------- #
class TestArgparsePlumbing:
    def test_validate_is_registered_alongside_the_others(self):
        parser = cli.build_parser()
        sub = next(a for a in parser._subparsers._group_actions
                   if hasattr(a, "choices"))
        assert "validate" in sub.choices

    def test_help_exits_cleanly(self):
        with pytest.raises(SystemExit) as exc:
            cli.build_parser().parse_args(["validate", "--help"])
        assert exc.value.code == 0

    def test_defaults(self):
        args = cli.build_parser().parse_args(["validate"])
        assert args.quick is False
        assert args.json is False
        assert args.only is None
        assert args.list_anchors is False
        assert args.exchange is False
        assert args.func is cli.cmd_validate

    @pytest.mark.parametrize("flag,attr", [
        ("--quick", "quick"),
        ("--json", "json"),
        ("--list", "list_anchors"),
        ("--exchange", "exchange"),
    ])
    def test_boolean_flags_parse(self, flag, attr):
        args = cli.build_parser().parse_args(["validate", flag])
        assert getattr(args, attr) is True

    def test_only_takes_a_comma_separated_list(self):
        args = cli.build_parser().parse_args(
            ["validate", "--only", "gbm-penetration,ps-scale"])
        assert args.only == "gbm-penetration,ps-scale"

    def test_list_runs_nothing(self, monkeypatch, capsys):
        # --list must be pure metadata: if it ever started the suite, a user
        # asking "what does this check?" would pay for a full validation run.
        def _boom(*a, **kw):
            raise AssertionError("run_suite must not be called for --list")

        monkeypatch.setattr(SC, "run_suite", _boom)
        assert cli.main(["validate", "--list"]) == 0
        out = capsys.readouterr().out
        assert "Registered external anchors" in out
        assert "DELIBERATELY NOT registered" in out


# --------------------------------------------------------------------------- #
# 2. The band rule, in isolation
# --------------------------------------------------------------------------- #
class TestClassify:
    @pytest.mark.parametrize("value,expected", [
        (0.5, SC.FAIL), (0.99, SC.FAIL),
        (1.0, SC.PASS), (1.5, SC.PASS), (2.0, SC.PASS),
        (2.01, SC.FAIL), (100.0, SC.FAIL),
    ])
    def test_bounds_are_inclusive(self, value, expected):
        assert SC.classify(value, 1.0, 2.0) == expected

    def test_reversed_bounds_are_tolerated(self):
        # A caller passing (hi, lo) should not silently fail everything.
        assert SC.classify(1.5, 2.0, 1.0) == SC.PASS
        assert SC.classify(3.0, 2.0, 1.0) == SC.FAIL

    def test_nan_prediction_never_passes(self):
        assert SC.classify(float("nan"), 0.0, 1e9) == SC.FAIL


# --------------------------------------------------------------------------- #
# 3. PASS / FAIL / NOT TESTABLE / DEFERRED classification
# --------------------------------------------------------------------------- #
class TestClassification:
    def test_pass_only_exits_zero(self, monkeypatch):
        _install(monkeypatch, [_group("ok", lambda: [_row("ok", SC.PASS)])])
        res = SC.run_suite()
        assert res["counts"][SC.PASS] == 1
        assert res["exit_code"] == 0

    def test_any_fail_exits_nonzero(self, monkeypatch):
        _install(monkeypatch, [
            _group("ok", lambda: [_row("ok", SC.PASS)]),
            _group("bad", lambda: [_row("bad", SC.FAIL)]),
        ])
        res = SC.run_suite()
        assert res["counts"][SC.FAIL] == 1
        assert res["exit_code"] == 1

    def test_missing_data_is_not_testable_not_a_pass(self, monkeypatch):
        def _never_called():
            raise AssertionError("runner must not run when its data is absent")

        _install(monkeypatch, [
            _group("needsdata", _never_called,
                   needs=("data_cache/definitely-not-here-9f2a",)),
        ])
        res = SC.run_suite()
        assert res["counts"][SC.NOT_TESTABLE] == 1
        assert res["counts"][SC.PASS] == 0
        assert res["exit_code"] == 0
        reason = res["rows"][0]["reason"]
        assert "REQUIRED DATA ABSENT" in reason
        assert "definitely-not-here-9f2a" in reason

    def test_runner_declining_is_not_testable_with_its_reason(self, monkeypatch):
        def _decline():
            raise SC.NotTestable("no measured target exists for this quantity")

        _install(monkeypatch, [_group("decline", _decline)])
        res = SC.run_suite()
        assert res["rows"][0]["status"] == SC.NOT_TESTABLE
        assert "no measured target exists" in res["rows"][0]["reason"]

    def test_a_crash_is_not_testable_never_a_pass(self, monkeypatch):
        # An anchor that blows up produced NO verdict. Swallowing the exception
        # into a pass (or dropping the row) is the exact failure mode this
        # command exists to prevent.
        def _crash():
            raise ValueError("upstream module changed shape")

        _install(monkeypatch, [_group("crash", _crash)])
        res = SC.run_suite()
        row = res["rows"][0]
        assert row["status"] == SC.NOT_TESTABLE
        assert "CRASHED" in row["reason"]
        assert "ValueError" in row["reason"]
        assert res["counts"][SC.PASS] == 0

    def test_quick_defers_full_tier_and_never_calls_it(self, monkeypatch):
        def _never_called():
            raise AssertionError("a full-tier runner must not run under --quick")

        _install(monkeypatch, [
            _group("cheap", lambda: [_row("cheap", SC.PASS)]),
            _group("dear", _never_called, tier=SC.TIER_FULL),
        ])
        res = SC.run_suite(quick=True)
        assert res["counts"][SC.DEFERRED] == 1
        assert res["counts"][SC.PASS] == 1
        deferred = [r for r in res["rows"] if r["status"] == SC.DEFERRED][0]
        assert "NOT a pass" in deferred["reason"]

    def test_quick_does_not_turn_a_full_tier_failure_into_a_pass(self, monkeypatch):
        _install(monkeypatch, [
            _group("dear", lambda: [_row("dear", SC.FAIL)], tier=SC.TIER_FULL),
        ])
        quick = SC.run_suite(quick=True)
        full = SC.run_suite(quick=False)
        assert quick["counts"][SC.PASS] == 0      # deferred, not passed
        assert quick["counts"][SC.DEFERRED] == 1
        assert full["counts"][SC.FAIL] == 1
        assert full["exit_code"] == 1

    def test_every_verdictless_row_states_a_reason(self, monkeypatch):
        def _decline():
            raise SC.NotTestable("because")

        _install(monkeypatch, [
            _group("a", _decline),
            _group("b", lambda: [_row("b", SC.PASS)]),
            _group("c", lambda: [_row("c", SC.FAIL)], tier=SC.TIER_FULL),
            _group("d", lambda: [], needs=("data_cache/nope-4b1c",)),
        ])
        res = SC.run_suite(quick=True)
        for r in res["rows"]:
            if r["status"] in (SC.NOT_TESTABLE, SC.DEFERRED):
                assert r["reason"], f"{r['group']} has no stated reason"

    def test_not_testable_alone_does_not_fail_the_suite(self, monkeypatch):
        def _decline():
            raise SC.NotTestable("structurally unmeasurable")

        _install(monkeypatch, [_group("nt", _decline)])
        assert SC.run_suite()["exit_code"] == 0

    def test_analytic_rows_are_counted_separately_from_pde_rows(self, monkeypatch):
        # The PDE row must declare its delivery law, or _stamp_physics
        # downgrades it and it never reaches the scored count at all — see
        # TestPhysicsProvenance.
        law = SC.physics_label(SC.repo_default_exchange())
        _install(monkeypatch, [
            _group("an", lambda: [_row("an", SC.PASS, scope=SC.SCOPE_ANALYTIC)]),
            _group("pd", lambda: [_row("pd", SC.PASS, scope=SC.SCOPE_PDE,
                                       physics=law)]),
        ])
        res = SC.run_suite()
        assert res["n_scored"] == 2
        assert res["n_scored_not_pde"] == 1


# --------------------------------------------------------------------------- #
# 4. --only, --json and exit codes through the CLI
# --------------------------------------------------------------------------- #
class TestCliBehaviour:
    def test_only_selects_a_subset(self, monkeypatch, capsys):
        _install(monkeypatch, [
            _group("keep", lambda: [_row("keep", SC.PASS)]),
            _group("drop", lambda: [_row("drop", SC.FAIL)]),
        ])
        rc = cli.main(["validate", "--only", "keep"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "keep" in out
        assert "drop" not in out

    def test_unknown_only_key_is_a_clean_exit_2(self, monkeypatch, capsys):
        _install(monkeypatch, [_group("keep", lambda: [_row("keep", SC.PASS)])])
        rc = cli.main(["validate", "--only", "not-an-anchor"])
        assert rc == 2
        assert "not-an-anchor" in capsys.readouterr().out

    def test_failing_anchor_exits_one_through_the_cli(self, monkeypatch):
        _install(monkeypatch, [_group("bad", lambda: [_row("bad", SC.FAIL)])])
        assert cli.main(["validate"]) == 1

    def test_json_is_parseable_and_carries_the_scorecard(self, monkeypatch, capsys):
        _install(monkeypatch, [
            _group("ok", lambda: [_row("ok", SC.PASS)]),
            _group("bad", lambda: [_row("bad", SC.FAIL)]),
        ])
        rc = cli.main(["validate", "--json"])
        assert rc == 1
        payload = json.loads(capsys.readouterr().out)
        assert payload["exit_code"] == 1
        assert payload["counts"][SC.FAIL] == 1
        assert {r["group"] for r in payload["rows"]} == {"ok", "bad"}
        for r in payload["rows"]:
            assert r["status"] in SC.STATUS_ORDER
            assert r["quantity"]
        assert {g["key"] for g in payload["groups"]} == {"ok", "bad"}

    def test_json_emits_no_progress_chatter(self, monkeypatch, capsys):
        _install(monkeypatch, [_group("ok", lambda: [_row("ok", SC.PASS)])])
        cli.main(["validate", "--json"])
        # Must be a single parseable document, not JSON with log lines around it.
        json.loads(capsys.readouterr().out)

    def test_failures_are_repeated_in_a_dedicated_block(self, monkeypatch, capsys):
        _install(monkeypatch, [_group(
            "bad",
            lambda: [_row("bad", SC.FAIL, predicted="7.80", measured="0.61")],
        )])
        cli.main(["validate"])
        out = capsys.readouterr().out
        assert "FAILURES (1)" in out
        # once in the row block, once in the failures block
        assert out.count("7.80") >= 2

    def test_quick_mode_prints_the_not_a_green_suite_banner(self, monkeypatch,
                                                            capsys):
        _install(monkeypatch, [
            _group("dear", lambda: [_row("dear", SC.FAIL)], tier=SC.TIER_FULL),
        ])
        rc = cli.main(["validate", "--quick"])
        assert rc == 0            # nothing was scored, so nothing failed
        out = capsys.readouterr().out
        assert "QUICK MODE" in out
        assert "not a green suite" in out
        assert "DEFERRED IS NOT A PASS" in out

    def test_conditionals_cited_by_a_row_are_printed_in_full(self, monkeypatch,
                                                            capsys):
        _install(monkeypatch, [_group(
            "c", lambda: [_row("c", SC.PASS, conditional=("ps-scale",))])])
        cli.main(["validate"])
        out = capsys.readouterr().out
        assert "CONDITIONALS" in out
        assert "PS_SCALE_S_INV" in out


# --------------------------------------------------------------------------- #
# 5. Registry hygiene — cheap invariants over the REAL registry (no runs)
# --------------------------------------------------------------------------- #
class TestRegistryHygiene:
    def test_keys_are_unique(self):
        keys = [g.key for g in SC.build_registry()]
        assert len(keys) == len(set(keys))

    def test_every_anchor_names_its_source_report_and_tier(self):
        for g in SC.build_registry():
            assert g.source, f"{g.key} does not say which module computes it"
            assert g.report, f"{g.key} does not say where it is published"
            assert g.tier in (SC.TIER_QUICK, SC.TIER_FULL), g.key
            assert g.scope in (SC.SCOPE_PDE, SC.SCOPE_ANALYTIC, SC.SCOPE_GEOMETRY)
            assert callable(g.runner), g.key

    def test_quick_tier_needs_no_real_tissue_section(self):
        # The point of --quick is that it runs without the big Xenium/Visium
        # caches. Literature tables under data_cache/hunt/ are fine.
        for g in SC.build_registry():
            if g.tier != SC.TIER_QUICK:
                continue
            for need in g.needs:
                assert "xenium" not in need.lower(), (
                    f"{g.key} is quick-tier but needs a tissue section: {need}")
                assert "visium" not in need.lower(), g.key

    def test_conditional_keys_all_resolve(self):
        # A row citing "[foo]" that CONDITIONALS does not define would print a
        # caveat marker with no caveat behind it.
        import inspect

        src = inspect.getsource(SC)
        for key in SC.CONDITIONALS:
            assert SC.CONDITIONALS[key].strip(), key
        # every literal conditional tuple entry in the module must be a known key
        for key in ("alpha-fit", "ps-scale", "leak-fit", "analytic-only",
                    "gdsc-potency", "untreated-section"):
            assert key in SC.CONDITIONALS
            assert f'"{key}"' in src

    def test_out_of_scope_entries_carry_a_reason(self):
        assert SC.OUT_OF_SCOPE
        for what, why in SC.OUT_OF_SCOPE:
            assert what.strip() and why.strip()

    def test_run_suite_rejects_an_unknown_key(self):
        with pytest.raises(KeyError):
            SC.run_suite(only=["definitely-not-an-anchor"])

    def test_geometry_anchors_are_registered_and_marked_as_such(self):
        # Rows that score the reconstructed TISSUE rather than the transport
        # must be distinguishable from rows that score the solver, or the
        # summary's "n of m scored rows never invoke the PDE" line is wrong.
        scopes = {g.scope for g in SC.build_registry()}
        assert SC.SCOPE_GEOMETRY in scopes
        assert SC.SCOPE_PDE in scopes
        assert SC.SCOPE_ANALYTIC in scopes

    def test_both_tissues_are_represented(self):
        tissues = {g.tissue for g in SC.build_registry()}
        assert {"GBM", "breast"} <= tissues

    def test_mcl_coefficient_anchors_are_analytic_and_never_enter_the_solver(
            self):
        """The Hicks MCL rows score the solver's INPUT COEFFICIENTS, not any
        solved field, and must never be tagged as PDE evidence.

        This is the one mis-tag that would be actively harmful rather than
        merely wrong: three more PDE-scoped rows would move the summary's
        "N of N scored rows never invoke the transport PDE" line without a
        single extra solve, which is the precise disclosure that line exists to
        make. So the scope tag is asserted here, and so is the fact behind it --
        that running the anchor enters neither transport entry point.
        """
        group = next(g for g in SC.build_registry()
                     if g.key == "mcl-transport-coefficients")
        assert group.scope == SC.SCOPE_ANALYTIC
        assert group.tier == SC.TIER_QUICK
        assert not group.needs, "an analytic coefficient row needs no tissue"

        from tmesim.physics import diffusion as _diff
        from tmesim.sim import engine as _eng

        seen = []
        solver_files = {os.path.abspath(_diff.__file__),
                        os.path.abspath(_eng.__file__)}

        def _profile(frame, event, arg):
            if event == "call" and frame.f_code.co_filename in solver_files:
                seen.append(frame.f_code.co_name)
            return None

        sys.setprofile(_profile)
        try:
            rows = SC._stamp_physics(list(group.runner()))
        finally:
            sys.setprofile(None)

        assert not seen, f"the transport solver was entered: {sorted(set(seen))}"
        assert len(rows) == 3, "diffusivity, hindrance and k are separate rows"
        for r in rows:
            assert r.scope == SC.SCOPE_ANALYTIC, r.label
            assert r.physics == SC.PHYSICS_NOT_APPLICABLE, r.label
            assert "coefficient-not-field" in r.conditional, r.label
            assert "mcl-no-vasculature" in r.conditional, r.label

    def test_mcl_verdicts_are_scored_on_the_engine_as_it_ships(self):
        """The cellular-exclusion term is merged but DEFAULT-OFF, so it must not
        supply a verdict -- and the hindrance row must record that the engine
        lands SHORT rather than being marked a pass because it improved."""
        rows = SC._run_mcl_transport_coefficients()
        hind = next(r for r in rows if r.label.startswith("hindrance"))
        v = hind.values
        # Scored on the shipping compiled value, not on the cell-term-on arm.
        assert hind.status == SC.FAIL
        assert v["shipped_fold"] < v["cell_term_on_armB_fold"] < v["measured_fold"]
        assert v["armB_shortfall_x"] > 1.0, (
            "arm B must be recorded as SHORT of the measurement")
        assert "cell-exclusion-default-off" in hind.conditional

        diff = next(r for r in rows
                    if r.label.startswith("effective diffusivity"))
        lo, hi = diff.values["measured_band_um2s"]
        # The cell-term-on arm lands INSIDE the measured band; the shipping
        # value does not, and the shipping value is the one that is scored.
        assert lo <= diff.values["cell_term_on_armB_D_eff_um2s"] <= hi
        assert not (lo <= diff.values["shipped_D_eff_um2s"] <= hi)
        assert diff.status == SC.FAIL

    def test_mcl_consumption_rate_declines_to_score_and_says_which_way(self):
        """NOT TESTABLE here removes a would-be PASS, not a FAIL, and the row
        has to say so -- otherwise it reads as a softened failure."""
        row = next(r for r in SC._run_mcl_transport_coefficients()
                   if r.label.startswith("first-order"))
        assert row.status == SC.NOT_TESTABLE
        assert row.reason and row.bar == "n/a"
        assert "REMOVES a pass" in row.reason
        v = row.values
        # The two readings of the engine's own k straddle the measurement by an
        # order of magnitude, which is WHY no single verdict is available: one
        # would read as a clean failure and the other as a comfortable pass.
        assert max(v["raw_ratio_to_measured"]) < 0.1
        assert min(v["calibrated_ratio_to_measured"]) > 0.3
        assert v["k_calibrated_per_s"] > 10 * v["k_raw_per_s"]


class TestRowHonestyInvariants:
    """Invariants that keep a row from quietly overstating itself.

    Driven through the synthetic registry so they cost nothing, but they
    encode the rules the real rows are written to.
    """

    def test_a_scored_row_must_name_its_bar_and_cite_a_pmid(self, monkeypatch):
        # A PASS or FAIL with no bar is unfalsifiable, and one with no PMID is
        # not an EXTERNAL anchor at all.
        _install(monkeypatch, [
            _group("a", lambda: [_row("a", SC.PASS, bar="inside 1-2")]),
            _group("b", lambda: [_row("b", SC.FAIL, bar="inside 1-2")]),
        ])
        for r in SC.run_suite()["rows"]:
            assert r["bar"], f"{r['group']} is scored but names no bar"
            assert r["pmids"], f"{r['group']} is scored but cites no PMID"

    def test_the_real_registry_obeys_the_same_invariants_where_it_can_run(self):
        # Only the anchors that need no data at all, so this stays data-free.
        res = SC.run_suite(only=["ps-scale", "immune-sink"])
        assert res["rows"]
        for r in res["rows"]:
            if r["status"] in (SC.PASS, SC.FAIL):
                assert r["bar"] and r["pmids"]
            else:
                assert r["reason"]
            assert r["quantity"] and "MEASURED" in r["quantity"].upper()

    def test_deferred_and_not_testable_never_reach_the_pass_count(self,
                                                                  monkeypatch):
        def _decline():
            raise SC.NotTestable("nope")

        _install(monkeypatch, [
            _group("nt", _decline),
            _group("df", lambda: [_row("df", SC.PASS)], tier=SC.TIER_FULL),
            _group("md", lambda: [], needs=("data_cache/absent-7c3e",)),
        ])
        res = SC.run_suite(quick=True)
        assert res["counts"][SC.PASS] == 0
        assert res["n_scored"] == 0
        assert res["exit_code"] == 0

    def test_bridge_pmids_are_recovered_from_report_local_row_tuples(self):
        # The breast kill report builds several bands from row tuples that
        # carry no ``scope``, so the scope-based lookup finds nothing for them.
        # A scored row with no PMID is not an EXTERNAL anchor, so this path
        # must not silently degrade to "n/a".
        bridge = {"n": 2, "rows": [
            {"pmid": "29221130", "fold": 3.0},
            {"pmid": "33431915", "fold": 5.0},
            {"pmid": "29221130", "fold": 4.0},   # duplicate, must collapse
            {"fold": 6.0},                        # no pmid, must be skipped
        ]}
        assert SC._rows_pmids(bridge) == ("29221130", "33431915")
        assert SC._rows_pmids(None) == ()
        assert SC._rows_pmids({"n": 0, "rows": []}) == ()

    def test_conditional_keys_on_real_rows_are_all_defined(self):
        res = SC.run_suite(only=["ps-scale", "immune-sink"])
        for key in res["conditionals"]:
            assert key in SC.CONDITIONALS
        for r in res["rows"]:
            for key in r["conditional"]:
                assert key in SC.CONDITIONALS, f"{r['group']} cites unknown {key}"


class TestPhysicsProvenance:
    """The delivery law must be visible on every row that could depend on it.

    The breast interstitial:plasma anchor moves by ~10x between the one-way
    leak and gradient-driven exchange, so a scorecard that prints a bare number
    is ambiguous by an order of magnitude. These tests pin the rule that a
    solver-derived verdict without a declared law is downgraded rather than
    allowed to pass.
    """

    def test_repo_default_is_read_from_pdeconfig_not_hardcoded(self):
        from tmesim.core import PDEConfig

        assert SC.repo_default_exchange() is bool(
            PDEConfig().transvascular_exchange)

    def test_physics_label_names_the_law_and_flags_non_default(self):
        default = SC.repo_default_exchange()
        assert "transvascular_exchange" in SC.physics_label(default)
        assert "NOT THE SHIPPING DEFAULT" not in SC.physics_label(default)
        assert "NOT THE SHIPPING DEFAULT" in SC.physics_label(not default)
        assert SC.physics_label(None) == SC.PHYSICS_NOT_APPLICABLE

    def test_non_pde_rows_are_stamped_not_applicable(self, monkeypatch):
        _install(monkeypatch, [
            _group("an", lambda: [_row("an", SC.PASS, scope=SC.SCOPE_ANALYTIC)]),
            _group("ge", lambda: [_row("ge", SC.PASS, scope=SC.SCOPE_GEOMETRY)]),
        ])
        for r in SC.run_suite()["rows"]:
            assert r["physics"] == SC.PHYSICS_NOT_APPLICABLE

    def test_a_pde_verdict_without_declared_physics_is_not_a_pass(self,
                                                                  monkeypatch):
        # The load-bearing one: a solver number with no declared delivery law
        # must NOT be reported as passing.
        _install(monkeypatch, [
            _group("pde", lambda: [_row("pde", SC.PASS, scope=SC.SCOPE_PDE)]),
        ])
        res = SC.run_suite()
        row = res["rows"][0]
        assert row["status"] == SC.NOT_TESTABLE
        assert res["counts"][SC.PASS] == 0
        assert "PHYSICS PROVENANCE MISSING" in row["reason"]
        assert row["bar"] == ""

    def test_a_pde_failure_without_declared_physics_is_also_downgraded(
            self, monkeypatch):
        _install(monkeypatch, [
            _group("pde", lambda: [_row("pde", SC.FAIL, scope=SC.SCOPE_PDE)]),
        ])
        res = SC.run_suite()
        assert res["rows"][0]["status"] == SC.NOT_TESTABLE
        assert res["exit_code"] == 0     # nothing was scored, so nothing failed

    def test_a_pde_row_that_declares_its_physics_keeps_its_verdict(self,
                                                                   monkeypatch):
        law = SC.physics_label(SC.repo_default_exchange())
        _install(monkeypatch, [_group(
            "pde",
            lambda: [_row("pde", SC.PASS, scope=SC.SCOPE_PDE, physics=law)])])
        res = SC.run_suite()
        assert res["rows"][0]["status"] == SC.PASS
        assert res["rows"][0]["physics"] == law

    def test_a_verdictless_pde_row_is_not_downgraded_further(self, monkeypatch):
        # immune-sink is PDE-scoped but produces no number; it must stay
        # NOT TESTABLE with its own reason, not acquire the provenance one.
        _install(monkeypatch, [_group(
            "pde",
            lambda: [_row("pde", SC.NOT_TESTABLE, scope=SC.SCOPE_PDE,
                          reason="no measured target exists")])])
        row = SC.run_suite()["rows"][0]
        assert row["status"] == SC.NOT_TESTABLE
        assert row["reason"] == "no measured target exists"
        assert row["physics"] == SC.PHYSICS_NOT_APPLICABLE

    def test_suite_result_carries_the_run_physics(self, monkeypatch):
        _install(monkeypatch, [_group("ok", lambda: [_row("ok", SC.PASS)])])
        default = SC.run_suite()["physics"]
        assert default["transvascular_exchange_requested"] is False
        assert default["is_shipping_default"] is (
            SC.repo_default_exchange() is False)
        assert default["law"] == SC.PHYSICS_LAWS[False]
        other = SC.run_suite(exchange=True)["physics"]
        assert other["transvascular_exchange_requested"] is True
        assert other["law"] == SC.PHYSICS_LAWS[True]

    def test_cli_prints_the_physics_banner(self, monkeypatch, capsys):
        _install(monkeypatch, [_group("ok", lambda: [_row("ok", SC.PASS)])])
        cli.main(["validate"])
        out = capsys.readouterr().out
        assert "PHYSICS" in out
        assert "transvascular_exchange" in out

    def test_cli_shouts_when_the_run_is_not_the_shipping_default(
            self, monkeypatch, capsys):
        _install(monkeypatch, [_group("ok", lambda: [_row("ok", SC.PASS)])])
        cli.main(["validate", "--exchange"])
        out = capsys.readouterr().out
        assert "NOT THE SHIPPING DEFAULT" in out

    def test_json_carries_physics_on_the_run_and_on_every_row(self, monkeypatch,
                                                              capsys):
        _install(monkeypatch, [_group("ok", lambda: [_row("ok", SC.PASS)])])
        cli.main(["validate", "--json"])
        payload = json.loads(capsys.readouterr().out)
        assert "physics" in payload
        assert payload["physics"]["law"]
        for r in payload["rows"]:
            assert r["physics"], f"{r['group']} carries no physics provenance"

    def test_real_data_free_anchors_all_declare_physics(self):
        # ps-scale and immune-sink need no data_cache, so this exercises the
        # stamping on REAL registry rows rather than synthetic ones.
        for r in SC.run_suite(only=["ps-scale", "immune-sink"])["rows"]:
            assert r["physics"], f"{r['group']} carries no physics provenance"
