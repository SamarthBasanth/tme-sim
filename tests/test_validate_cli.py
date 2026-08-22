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

One documented exception to "data-free": ``TestBrainSlicePreregisteredRow``
reads ``data_cache/brain_slice/vu_brain.csv``, a 30 kB committed literature
table — not a tissue section, no PDE, ~3 s. It is there because that row is the
enforcement of a PRE-REGISTERED failure, and a pre-registered failure that
nothing tests is the exact defect the row exists to close.

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
        """RE-PINNED at the two-route promotion (results/two_route_promotion.md).

        SUPERSEDED CONTRACT: the cell term was merged but DEFAULT-OFF, so the
        hindrance row had to be FAIL (shipping 1.77-fold against a measured
        12-fold) and had to record that even the cell-term-on arm landed SHORT.
        The promotion put a PERMEABLE-cell term into the shipping compile, so:

        | | superseded | current |
        |---|---|---|
        | hindrance, shipping fold | 1.7732 | 10.3274 |
        | hindrance verdict | FAIL | **PASS** |
        | D_eff, shipping um^2/s | 287.32 | 49.33 |
        | D_eff verdict | FAIL | FAIL (direction flipped) |

        The verdict move is REAL and is the point of the change. What this test
        now pins is that it was not bought by moving the bar: the band is still
        the repo's own 3x ORDER_OF_MAGNITUDE_FACTOR either side of the measured
        12-fold, chosen while the engine was failing it, and the ABSOLUTE row is
        still scored and still FAILS.
        """
        rows = SC._run_mcl_transport_coefficients()
        hind = next(r for r in rows if r.label.startswith("hindrance"))
        v = hind.values
        assert hind.status == SC.PASS
        # ... and the pass must be inside the UNCHANGED band, not a widened one.
        blo, bhi = v["band_fold"]
        assert (blo, bhi) == (4.0, 36.0)
        assert blo <= v["shipped_fold"] <= bhi
        # the promotion is what moved it, and the old value is still recorded
        assert v["superseded_mw_only_fold"] < blo < v["shipped_fold"]
        # the pass is still SHORT of the central measured value; say so
        assert v["shipped_fold"] < v["measured_fold"]

        diff = next(r for r in rows
                    if r.label.startswith("effective diffusivity"))
        lo, hi = diff.values["measured_band_um2s"]
        # The ABSOLUTE row is still scored, still fails, and the miss flipped
        # sides -- which is why "the transport model is now right" is not a
        # conclusion this group supports.
        assert not (lo <= diff.values["shipped_D_eff_um2s"] <= hi)
        assert diff.status == SC.FAIL
        assert diff.values["superseded_mw_only_D_eff_um2s"] > hi
        assert diff.values["shipped_D_eff_um2s"] < lo

    def test_mcl_consumption_rate_declines_to_score_and_says_which_way(self):
        """NOT TESTABLE has to say which way the omission cuts -- and AT THE
        PROMOTION IT REVERSED, which is the whole reason this test exists.

        SUPERSEDED: at alpha = 148.6082 the calibrated reading was 0.49-0.82x
        measured, i.e. inside the repo's 3x band, so declining to score removed a
        would-be PASS. At alpha = 22.2235 it is 0.080-0.134x measured, outside
        that band, so declining to score now removes a would-be FAILURE. Leaving
        the old sentence in place would have been a disclosure that flattered us
        while claiming to do the opposite.
        """
        row = next(r for r in SC._run_mcl_transport_coefficients()
                   if r.label.startswith("first-order"))
        assert row.status == SC.NOT_TESTABLE
        assert row.reason and row.bar == "n/a"
        assert "removes a FAILURE" in row.reason, (
            "the row must state which way declining to score cuts, and since "
            "results/two_route_promotion.md it cuts the other way")
        v = row.values
        # Both readings are now BELOW the measurement -- they no longer straddle
        # it -- so the reason can no longer rest on 'straddling'. It rests on
        # commensurability, which is reason (1) and is untouched.
        assert max(v["raw_ratio_to_measured"]) < 0.1
        assert max(v["calibrated_ratio_to_measured"]) < 0.3
        assert v["k_calibrated_per_s"] > 5 * v["k_raw_per_s"]


class TestBrainSlicePreregisteredRow:
    """The pre-registered rat brain-slice tissue-capacity anchor (TEST B).

    Why this class exists, stated plainly
    -------------------------------------
    ``results/transport_anchor_prereg.md`` fixed the dataset cascade, the
    statistic, the 3.0 bar, the bootstrap and the withdrawal conditions BEFORE
    the prediction was computed; ``results/brain_slice_anchor.md`` executed it
    and it FAILED (ratio 12.05 against a bar of 3.0). Prereg sec. 8.5 requires
    that failure to appear in ``vivum validate`` and to make it exit 1. Until
    it did, the failure was recorded in a markdown file while the scorecard
    reported as though nothing had happened — a check that had stopped
    checking. These tests pin the three properties that keep it checking:

    1. the scope tag is ``analytic`` and the solver really is never entered,
       so the row can never be quoted as PDE evidence (prereg sec. 8.5 says
       ``SCOPE_ANALYTIC``, **never** ``SCOPE_PDE``);
    2. ``GMFE_floor`` — the ENGINE-INDEPENDENT headline, the best any constant
       predictor could do — is not droppable from the row;
    3. the verdict is recomputed live from the shipping engine and still
       matches the published execution record, so the row cannot go stale
       against the document, and a softened verdict cannot pass unnoticed.

    Not data-free (see the module docstring): it reads the 30 kB committed
    literature table the anchor is scored against. No tissue, no PDE.
    """

    #: Published in ``results/brain_slice_anchor.md`` sec. 0.2 / sec. 3. The row
    #: recomputes these live; if it ever stops reproducing them, either the
    #: engine moved (which must be re-registered under prereg sec. 9.5) or the
    #: row stopped scoring what the document says it scores.
    PUBLISHED = {
        "engine_vu_brain_ml_per_g": 1.000,
        "gmfe_engine_x": 57.50,
        "gmfe_floor_x": 4.77,
        "ratio": 12.045,
        "gmfe_floor_optimal_constant_ml_per_g": 39.0,
        "n": 33,
    }

    @pytest.fixture(scope="class")
    def row(self):
        need = os.path.join(SC.REPO_ROOT, "data_cache", "brain_slice",
                            "vu_brain.csv")
        if not os.path.exists(need):
            pytest.skip(f"frozen literature table absent: {need}")
        rows = SC._stamp_physics(list(SC._run_brain_slice_capacity()))
        assert len(rows) == 1, "TEST B is one comparison, not several"
        return rows[0]

    def test_it_is_registered_at_all(self):
        group = next((g for g in SC.build_registry()
                      if g.key == "brain-slice-tissue-capacity"), None)
        assert group is not None, (
            "prereg sec. 8.5 requires this anchor to be registered in "
            "vivum validate; without it the FAIL is recorded but not enforced")
        assert group.tier == SC.TIER_QUICK
        assert group.needs, "the row must declare the table it scores against"

    def test_it_cites_the_prereg_and_the_execution_record(self, row):
        """A reader must be able to see the criterion was fixed in advance."""
        group = next(g for g in SC.build_registry()
                     if g.key == "brain-slice-tissue-capacity")
        assert "transport_anchor_prereg.md" in group.report
        assert "brain_slice_anchor.md" in group.report
        # and in the row itself, not only in the group header
        assert "transport_anchor_prereg.md" in row.bar
        assert row.values["prereg"].endswith("transport_anchor_prereg.md")
        assert row.values["execution_record"].endswith("brain_slice_anchor.md")
        assert row.pmids, "an EXTERNAL anchor must cite its sources"

    def test_scope_is_analytic_and_the_solver_is_never_entered(self, row):
        """It scores an INPUT COEFFICIENT. Tagging it PDE would move the
        summary's 'N of N scored rows never invoke the transport PDE' line
        without a single solve — the mis-tag prereg sec. 8.5 forbids by name."""
        group = next(g for g in SC.build_registry()
                     if g.key == "brain-slice-tissue-capacity")
        assert group.scope == SC.SCOPE_ANALYTIC
        assert row.scope == SC.SCOPE_ANALYTIC
        assert row.physics == SC.PHYSICS_NOT_APPLICABLE
        assert "coefficient-not-field" in row.conditional
        # verified, not asserted: the runner's own live call counters
        assert row.values["solver_calls"] == {
            "simulate_drug_transport": 0, "simulate_response": 0}

    def test_the_solver_is_never_entered_under_a_profiler_too(self, row):
        """Call counters catch the two named entry points; the profiler catches
        indirect and renamed calls that a wrapper would miss."""
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
            SC._run_brain_slice_capacity()
        finally:
            sys.setprofile(None)
        assert not seen, f"the transport solver was entered: {sorted(set(seen))}"

    def test_the_verdict_is_FAIL_and_it_makes_the_suite_exit_one(self, row):
        assert row.status == SC.FAIL
        lo, hi = row.values["ratio_ci90"]
        bar = row.values["the_test_bar"]
        assert bar == 3.0, "the pre-registered bar may not be widened"
        assert lo > bar, (
            "prereg sec. 3.6: the verdict stands only because the 90 % BCa "
            "interval lies ENTIRELY above the bar; a straddle is B-C "
            "UNDERPOWERED, which is NOT a soft pass")
        assert row.values["verdict_bucket"].startswith("B-B")

    def test_the_row_carries_verdict_ratio_AND_floor(self, row):
        """GMFE_floor is the engine-independent headline (prereg sec. 3.5a) and
        must not be droppable: it is what says the best ANY constant could do
        is still ~4.8x, i.e. what quantifies the missing per-compound term."""
        v = row.values
        for key in ("gmfe_floor_x", "gmfe_floor_ci90",
                    "gmfe_floor_optimal_constant_ml_per_g", "ratio",
                    "ratio_ci90", "gmfe_engine_x", "verdict_bucket"):
            assert key in v, f"{key} dropped from the row"
        # ...and visibly, in the text a reader actually sees, not only in the
        # machine-readable values.
        assert "GMFE_floor" in row.measured and "GMFE_floor" in row.note
        assert f"{v['gmfe_floor_x']:.2f}x" in row.note
        assert f"{v['ratio']:.3f}" in row.predicted
        assert v["gmfe_floor_x"] == pytest.approx(v["gmfe_engine_x"] / v["ratio"])
        # the floor is a MEASURED property: no engine quantity may enter it
        assert v["gmfe_floor_x"] < v["gmfe_engine_x"]
        assert v["gmfe_floor_x"] > v["reproducibility_floor_x"], (
            "if the floor ever fell to the assay's own noise floor, a single "
            "constant WOULD be an adequate model and this row's reading "
            "changes; that is a finding, not a silent pass")

    def test_recentring_the_constant_is_not_reported_as_a_fix(self, row):
        # The trap prereg sec. 3.5a exists to close: a re-centred constant
        # passes THE TEST by construction and is still wrong by GMFE_floor.
        assert "RE-CENTRING IS THEREFORE NOT A FIX" in row.note
        assert "BY CONSTRUCTION" in row.note

    def test_it_reproduces_the_published_execution_record(self, row):
        """Live recomputation, pinned to results/brain_slice_anchor.md."""
        for key, want in self.PUBLISHED.items():
            got = row.values[key]
            assert got == pytest.approx(want, rel=2e-3), (
                f"{key}: live {got} vs published {want}. Either the engine "
                f"moved (prereg sec. 9.5 requires re-registration, not a "
                f"quiet edit) or this row stopped scoring what "
                f"results/brain_slice_anchor.md says it scores.")
        assert row.values["frozen_constants_drift"] == "", (
            "a prereg sec. 2.3 frozen constant has moved since execution")

    def test_the_engine_prediction_is_the_same_for_every_compound(self, row):
        # The whole content of the failure: one capacity for 33 compounds whose
        # measured capacities span three orders of magnitude.
        assert row.values["s_eff"] == 1.0
        assert row.values["n_below_engine"] == row.values["n"]
        assert row.values["within_3x_pct"] == 0.0
        lo, hi = row.values["measured_range_ml_per_g"]
        assert hi / lo > 1000
        assert "NOT COMPUTABLE" in row.values["spearman_rho"], (
            "a zero-variance prediction has no ranks; reporting rho = 0.00 "
            "would be inventing a number (prereg sec. 3.5)")

    def test_BW1_is_recorded_as_triggered(self, row):
        lo, hi = row.values["measured_range_ml_per_g"]
        pred = row.values["engine_vu_brain_ml_per_g"]
        assert row.values["bw1_reachability_triggered"] is (
            not (lo <= pred <= hi))
        assert row.values["bw1_reachability_triggered"] is True
        assert "BW1" in row.note and "WITHDRAWN, not narrowed" in row.note

    def test_both_judgement_calls_are_recorded_not_papered_over(self, row):
        """Neither is resolved by this registration; both must stay visible.

        (1) GATE 1's engine quantity is not on prereg sec. 3.3's list, and a
        reader could hold that a capacity fixed at 1 by NORMALISATION is not a
        quantity the engine expresses — which would make this NOT TESTABLE
        (B-D) rather than FAIL. The row takes the less flattering reading and
        says the other was available.
        (2) Prereg sec. 4.2 is internally inconsistent about BW1's range. The
        operative sentence was applied; the inconsistency needs amendment under
        sec. 9.5, which an implementer may not do.
        """
        assert "ALTERNATIVE READING" in row.note
        assert "B-D" in row.note and "NOT TESTABLE" in row.note
        assert "LESS FLATTERING READING" in row.note
        assert "INTERNALLY INCONSISTENT" in row.note
        assert "sec. 9.5" in row.note
        # the alternative reading must not be laundered as a pass
        assert "Both readings are non-passes" in row.note
        # and the row must say the test fails under either reading of sec. 4.2
        assert "NOTHING TURNS ON IT" in row.note

    def test_it_does_not_license_moving_cell_partition(self, row):
        """prereg sec. 3.2: Vu,brain is the LUMP, cell_partition is one
        component of it. This FAIL constrains neither K nor the solver."""
        assert "capacity-normalised-to-one" in row.conditional
        assert "brain-slice-not-tumour" in row.conditional
        text = SC.CONDITIONALS["capacity-normalised-to-one"]
        assert "DOES NOT LICENSE MOVING cell_partition" in text
        assert "K remains unconstrained by any measurement" in text
        bound = SC.CONDITIONALS["brain-slice-not-tumour"]
        assert "does NOT license" in bound and "human tumour" in bound

    def test_the_registered_row_reaches_the_suite_and_fails_it(self):
        res = SC.run_suite(only=["brain-slice-tissue-capacity"])
        assert res["counts"][SC.FAIL] == 1
        assert res["counts"][SC.PASS] == 0
        assert res["exit_code"] == 1, (
            "prereg sec. 8.5: an anchor that fails must make vivum validate "
            "exit 1")
        # it is scored, and it is NOT solver evidence
        assert res["n_scored"] == 1
        assert res["n_scored_not_pde"] == 1
        assert res["rows"][0]["bar"] and res["rows"][0]["pmids"]
        assert "MEASURED" in res["rows"][0]["quantity"].upper()
        for key in res["rows"][0]["conditional"]:
            assert key in SC.CONDITIONALS


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


# --------------------------------------------------------------------------- #
# 5. Zero-width bands, and the two paths that must not disagree about them
# --------------------------------------------------------------------------- #
#
# The defect these guard against, recorded in results/zero_width_bands.md:
# ``vivum validate`` scored the breast kill panel 3 PASS / 3 FAIL while
# results/kill_calibration_breast.md published 5 PASS / 1 FAIL, because the
# scorecard and the report generator each decided SEPARATELY which band a drug
# is scored against. Two of the scorecard's FAILs were rendered against
# zero-width bands (0.1636-0.1636 and 0.0412-0.0412) built by taking a min-max
# over a SINGLE literature row.
class TestZeroWidthBands:
    """A band with no width is not a bar. It must never produce a verdict."""

    @pytest.mark.parametrize("predicted", [0.0, 0.0412, 0.1636, 0.5, 1.0, 1e9])
    def test_a_degenerate_band_is_never_scored(self, predicted):
        from tmesim.validation import kill_calibration as KC

        t = KC.ratio_test(predicted, 0.0412, 0.0412)
        assert t["verdict"] == SC.NOT_TESTABLE
        assert t["verdict"] != SC.FAIL
        assert t["degenerate"] is True
        assert t["reason"], "a verdictless row must say why"

    def test_a_real_interval_is_still_scored_both_ways(self):
        from tmesim.validation import kill_calibration as KC

        assert KC.ratio_test(0.5, 0.2, 0.8)["verdict"] == SC.PASS
        assert KC.ratio_test(0.9, 0.2, 0.8)["verdict"] == SC.FAIL
        assert KC.ratio_test(0.1, 0.2, 0.8)["verdict"] == SC.FAIL

    def test_band_degeneracy_covers_crossed_and_non_finite_edges(self):
        from tmesim.validation import kill_calibration as KC

        assert KC.band_is_degenerate(0.3, 0.3) is True
        assert KC.band_is_degenerate(float("nan"), 0.3) is True
        assert KC.band_is_degenerate(float("inf"), 0.3) is True
        assert KC.band_is_degenerate(0.2, 0.8) is False
        # A crossed pair still encloses an interval and is tolerated, exactly
        # as SC.classify tolerates reversed bounds.
        assert KC.band_is_degenerate(0.8, 0.2) is False

    def test_one_study_is_a_point_estimate_not_an_interval(self):
        from tmesim.validation import kill_calibration as KC

        bridge = {"n": 1, "rows": [{"pmid": "41422245", "fold": 25.0}],
                  "fold_median": 25.0, "fold_min": 25.0, "fold_max": 25.0}
        band = KC.bridge_to_kill_ratio_band(bridge, 1.0, 1.0, 1.0)
        assert band["degenerate"] is True
        assert band["n_rows"] == 1
        s = KC.score_ratio_row(0.1503, specific_bridge=bridge, specific_band=band)
        assert s["status"] == SC.NOT_TESTABLE
        assert s["scored"] is False
        assert "n=1" in s["reason"]

    def test_the_broad_band_is_never_substituted_to_rescue_a_point_estimate(self):
        # The honest fix must not become "widen the bar until it passes": the
        # broad cross-tumour band is CONTEXT on a single-study row, never the
        # scored band, even when the prediction would land inside it.
        from tmesim.validation import kill_calibration as KC

        one = {"n": 1, "fold_median": 25.0, "fold_min": 25.0, "fold_max": 25.0}
        broad = {"n": 20, "fold_median": 4.0, "fold_min": 0.4, "fold_max": 25.0}
        band = KC.bridge_to_kill_ratio_band(one, 1.0, 1.0, 1.0)
        ctx = KC.bridge_to_kill_ratio_band(broad, 1.0, 1.0, 1.0)
        s = KC.score_ratio_row(0.5, specific_bridge=one, specific_band=band,
                               context_bridge=broad, context_band=ctx)
        assert s["context_test"]["verdict"] == SC.PASS   # it WOULD have passed
        assert s["status"] == SC.NOT_TESTABLE            # and it still is not scored
        assert s["band"] is band and s["test"] is None

    def test_two_studies_do_define_an_interval(self):
        from tmesim.validation import kill_calibration as KC

        bridge = {"n": 2, "fold_median": 33.0, "fold_min": 8.3, "fold_max": 58.0}
        band = KC.bridge_to_kill_ratio_band(bridge, 1.0, 1.0, 1.0)
        assert band["degenerate"] is False
        s = KC.score_ratio_row(band["ratio_central"], specific_bridge=bridge,
                               specific_band=band)
        assert s["scored"] is True
        assert s["status"] in (SC.PASS, SC.FAIL)

    def test_the_kill_row_renderer_refuses_a_degenerate_band(self):
        # Belt and braces: even if a caller hands _kill_rows a pre-computed
        # verdict against a zero-width band, the renderer downgrades it.
        t = {"predicted_ratio": 0.34, "lit_ratio_lo": 0.1636,
             "lit_ratio_hi": 0.1636, "verdict": SC.FAIL, "direction": "outside",
             "implied_fold": 3.77, "fold_lo": 7.5, "fold_hi": 7.5,
             "scope": "drug:5-fluorouracil"}
        rows = SC._kill_rows("g", "tissue", {"5-fluorouracil": t}, {}, ())
        assert [r.status for r in rows] == [SC.NOT_TESTABLE]
        assert rows[0].reason


class TestBreastKillPanelCannotDiverge:
    """The scorecard and the published report must render ONE decision.

    Data-free: a synthetic ``res`` is fed to the real
    ``build_ratio_section``, and the verdict the DOCUMENT would print is
    compared cell-for-cell with the status ``vivum validate`` assigns. No
    tissue, no PDE, no data_cache.
    """

    @staticmethod
    def _fake_res():
        kcb = SC._import_experiment("kill_calibration_breast")
        # Ratios chosen to straddle the bands: whatever they score, the two
        # paths have to score them the SAME.
        per_drug = {}
        for i, name in enumerate(kcb.DRUGS):
            ratio = 0.05 + 0.06 * i
            per_drug[name] = {
                "ic50_uM": 1.0, "hill": 1.0,
                "headline": {"ratio_exact": ratio, "implied_ic50_fold": 1.0 / ratio,
                             "label": name, "c_bath_uM": 1.0},
            }
        return {"per_drug": per_drug}

    def _section(self):
        rep = SC._import_experiment("kill_calibration_breast_report")
        return rep, rep.build_ratio_section(self._fake_res())

    def test_document_verdict_matches_scorecard_status_for_every_drug(
            self, monkeypatch):
        rep, section = self._section()
        kcb = SC._import_experiment("kill_calibration_breast")
        res = self._fake_res()

        monkeypatch.setattr(kcb, "run", lambda **kw: res)
        monkeypatch.setattr(SC.os.path, "isdir", lambda p: True)
        sc_rows = {r.label: r.status for r in SC._run_breast_kill_3d2d()}

        assert set(sc_rows) == set(section)
        for name, row in section.items():
            line = rep.ratio_table_row(name, row)
            cell = line.rsplit("|", 2)[1].replace("*", "")
            doc_verdict = ("NOT TESTABLE" if "NOT TESTABLE" in cell
                           else cell.strip())
            assert doc_verdict == sc_rows[name], (
                f"{name}: report says {doc_verdict!r}, scorecard says "
                f"{sc_rows[name]!r} — the two paths have diverged again")

    def test_no_scorecard_row_is_scored_against_a_zero_width_band(
            self, monkeypatch):
        from tmesim.validation import kill_calibration as KC

        kcb = SC._import_experiment("kill_calibration_breast")
        res = self._fake_res()
        monkeypatch.setattr(kcb, "run", lambda **kw: res)
        monkeypatch.setattr(SC.os.path, "isdir", lambda p: True)
        for r in SC._run_breast_kill_3d2d():
            if r.status in (SC.PASS, SC.FAIL):
                lo = r.values["lit_ratio_lo"]
                hi = r.values["lit_ratio_hi"]
                assert not KC.band_is_degenerate(lo, hi), (
                    f"{r.label} carries a verdict against a zero-width band")

    def test_single_study_drugs_are_not_testable_not_failed(self, monkeypatch):
        # Pins the specific rows the defect was found on. If new literature
        # gives either drug a second study, this test SHOULD be updated -- and
        # the panel tally with it.
        _rep, section = self._section()
        for name in ("paclitaxel", "5-fluorouracil"):
            s = section[name]["scored"]
            assert s["n_rows"] == 1
            assert s["status"] == SC.NOT_TESTABLE
            assert s["scored"] is False
