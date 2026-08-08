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
        _install(monkeypatch, [
            _group("an", lambda: [_row("an", SC.PASS, scope=SC.SCOPE_ANALYTIC)]),
            _group("pd", lambda: [_row("pd", SC.PASS, scope=SC.SCOPE_PDE)]),
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
