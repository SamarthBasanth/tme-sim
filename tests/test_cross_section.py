"""ONE COMPOUND, MANY SECTIONS — the acceptance gate for the cross-section view.

What this file is guarding
--------------------------
``app/server/cross_section.py`` and ``app/web/js/ui/crosssection.js`` build the
one comparison this product can defend: the SAME compound on two or more pieces
of real tissue. The comparison it must never build — two DIFFERENT compounds on
one axis — was withdrawn on 2026-08-21 and is enforced at three layers already
(``projects._kill_comparability`` raises, ``protocol.executor._check_one_compound``
plans zero runs). This view is the third surface, and the tests below are
ordered by how much damage the defect they catch would do:

1. the axis carries exactly ONE compound, audited on the OUTPUT;
2. a spot-resolution section contributes NO per-cell statistic, with the
   refusal imported from ``app/server/api.py`` rather than paraphrased;
3. comparability is the BACKEND's verdict — literally
   ``projects._kill_comparability`` — and a refused row renders no numbers;
4. the DISAGREEMENT is computed, because a flip a reader has to spot for
   themselves is a flip that does not exist on screen.

``calibrate_alpha()`` is asserted first and last, per this repo's standing
rule: nothing in this feature touches the calibration chain, and a test file
that cannot demonstrate that is not evidence of it.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

pytest.importorskip("fastapi", reason="the Studio backend is an optional extra")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi.testclient import TestClient  # noqa: E402

from app.server import cross_section as xs  # noqa: E402
from app.server.projects import (  # noqa: E402
    _HEADLINE_METRICS,
    _drug_of,
    _kill_comparability,
)

#: The value this repo's calibration chain must produce. Asserted at the top and
#: the bottom of this file.
ALPHA = 22.223477711172563

JS = os.path.join(_ROOT, "app", "web", "js", "ui", "crosssection.js")


def test_calibrate_alpha_is_unchanged_at_the_START_of_this_file():
    from tmesim.validation.calibration import calibrate_alpha

    assert calibrate_alpha() == ALPHA


# --------------------------------------------------------------------------- #
# Fakes — a server made of runs, with no disk and no real sections
# --------------------------------------------------------------------------- #
LADDER = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]


class FakeJob:
    def __init__(self, view):
        self._view = view
        self.status = view.get("status", "done")

    def result_view(self):
        return dict(self._view)


class FakeJobs:
    """The shape ``cross_section`` uses of a ``JobManager``: ``list()``."""

    def __init__(self, views):
        self._jobs = [FakeJob(v) for v in views]

    def list(self, limit=None):
        return self._jobs[:limit] if limit else list(self._jobs)


def make_run(run_id, section, drug="temozolomide", *, fractions=None,
             ladder=None, n_viable=1000, n_tumor=1200, state_measured=True,
             window_mode="finite_window", reach="ok", null_reason=None,
             metrics=None, created_at=100.0, dx_um=12.0):
    """One finished run, shaped exactly like ``Job.result_view()``."""
    if reach is None:
        reach_block = None
    elif null_reason is not None:
        reach_block = {"null_reason": null_reason, "curve": None,
                       "window_mode": window_mode, "n_viable": None,
                       "n_tumor": None}
    else:
        reach_block = {
            "null_reason": None,
            "curve": {"x": list(ladder or LADDER),
                      "under_reached_fraction": list(
                          fractions or [0.01] * len(ladder or LADDER))},
            "window_mode": window_mode,
            "n_viable": n_viable,
            "n_tumor": n_tumor,
            "E_ref_uM_s": 1.234e5,
            "E_ref_units": "uM*s",
            "scale_invariance": "SCALE INVARIANCE NOTE (verbatim from the run)",
            "window_note": "WINDOW NOTE (verbatim from the run)",
            "convention": "CONVENTION NOTE (verbatim from the run)",
        }
    return {
        "run_id": run_id,
        "status": "done",
        "created_at": created_at,
        "spec": {"section": section, "drug": drug, "dx_um": dx_um,
                 "state_mode": "markers", "sensitivity_mode": "measured"},
        "metrics": metrics or {"kill_fraction": 0.12, "mean_tumor_auc": 3.4e6,
                               "therapeutic_coverage_1x": 0.3,
                               "penetration_length_um": 44.0},
        "provenance": {"state_is_measured": state_measured},
        "report": {},
        "reach": reach_block,
    }


def capability_for(section, *, panel="hBrain_254g", absolute_ok=True,
                   source="xenium"):
    return {
        "section": section,
        "label": f"{section} label",
        "available": True,
        "source": source,
        "panel": panel,
        "state_routes": {"markers": {"supported": True}},
        "supports_cell_state": True,
        "absolute_kill_supported": absolute_ok,
        "absolute_kill_reason": ("panel carries a hypoxia transcript"
                                 if absolute_ok else
                                 "no hypoxia-response transcript on this panel"),
    }


def build(views, *, compound="temozolomide", sources=None, caps=None, **kw):
    """Build a payload from fake runs, fake capabilities and fake platforms."""
    sources = sources or {}
    caps = caps or {}

    def capability(section):
        if section in caps:
            return caps[section]
        return capability_for(section, source=sources.get(section, "xenium"))

    def record(section):
        return {"source": sources.get(section, "xenium"), "label": section}

    return xs.cross_section_payload(
        FakeJobs(views), compound=compound, capability=capability,
        section_record=record, **kw)


# --------------------------------------------------------------------------- #
# 1. THE AXIS CARRIES ONE COMPOUND
# --------------------------------------------------------------------------- #
class TestOneCompoundOnOneAxis:
    """The withdrawn output, and the guard that cannot be walked past."""

    def test_the_axis_declares_what_is_fixed_and_what_varies(self):
        p = build([make_run("r1", "sec_a"), make_run("r2", "sec_b")])
        assert p["axis"]["fixed"] == "compound"
        assert p["axis"]["varies"] == "section"
        assert p["compound"] == "temozolomide"
        assert p["axis"]["n_sections"] == 2

    def test_a_second_compounds_runs_are_not_selected(self):
        p = build([make_run("r1", "sec_a"),
                   make_run("r2", "sec_b"),
                   make_run("r3", "sec_c", drug="doxorubicin")])
        assert [s["section"] for s in p["sections"]] == ["sec_a", "sec_b"]
        assert xs.assert_one_compound_axis(p) == "temozolomide"

    def test_every_section_row_carries_the_RUNS_OWN_compound(self):
        """Not the requested name.

        A row stamped with what the CALLER asked for agrees with itself no
        matter which run produced it, which is a guard that cannot fail. The
        run below names the compound with different case and spacing, and the
        row must still be read back off the run.
        """
        p = build([make_run("r1", "sec_a", drug="  Temozolomide "),
                   make_run("r2", "sec_b")])
        assert [s["compound"] for s in p["sections"]] == ["temozolomide"] * 2
        assert [s["spec"]["drug"] for s in p["sections"]][0] == "  Temozolomide "

    def test_the_guard_refuses_a_payload_that_names_two_compounds(self):
        p = build([make_run("r1", "sec_a"), make_run("r2", "sec_b")])
        p["sections"][1]["compound"] = "doxorubicin"
        with pytest.raises(xs.CrossSectionError) as exc:
            xs.assert_one_compound_axis(p)
        assert exc.value.status == 500
        assert "WITHDRAWN" in exc.value.message
        assert "pk_ranking_robustness.md" in exc.value.message
        assert "capacity_ranking_robustness.md" in exc.value.message
        assert exc.value.detail["compounds"] == ["doxorubicin", "temozolomide"]

    def test_the_guard_is_mechanical_it_walks_to_any_depth(self):
        """It audits the SHAPE, not the fields this module happens to write."""
        with pytest.raises(xs.CrossSectionError):
            xs.assert_one_compound_axis(
                {"compound": "a",
                 "something": {"new": [{"deeper": {"drugs": ["b"]}}]}})

    def test_a_two_compound_request_is_refused_with_the_withdrawal(self):
        for raw in ("temozolomide,doxorubicin", "temozolomide doxorubicin"):
            with pytest.raises(xs.CrossSectionError) as exc:
                build([make_run("r1", "sec_a")], compound=raw)
            assert exc.value.status == 400
            assert "WITHDRAWN" in exc.value.message

    def test_the_payload_builder_runs_the_guard_on_its_own_output(self):
        """The builder's last act. Asserted by source, because a guard that is
        defined and never called is the defect this test exists to catch."""
        src = open(os.path.join(_ROOT, "app", "server", "cross_section.py")).read()
        body = src.split("def cross_section_payload")[1]
        assert "assert_one_compound_axis(payload)" in body
        assert body.index("assert_one_compound_axis(payload)") < body.index("return payload")

    def test_the_names_this_module_borrows_from_projects_still_exist(self):
        """A private import is a coupling, so it is pinned rather than hoped
        for: a rename in projects.py fails HERE, not at request time."""
        assert callable(_kill_comparability) and callable(_drug_of)
        assert "kill_fraction" in _HEADLINE_METRICS


# --------------------------------------------------------------------------- #
# 2. SPOT SECTIONS CONTRIBUTE NO PER-CELL STATISTIC
# --------------------------------------------------------------------------- #
class TestSpotResolutionSectionsAreRefusedNotRendered:

    def test_which_platforms_are_spot_is_imported_not_assumed(self):
        from app.server.api import SPOT_SPATIAL_REFUSAL, _SPOT_SOURCES

        sources, refusal = xs._spot_rules()
        assert sources is _SPOT_SOURCES
        assert refusal is SPOT_SPATIAL_REFUSAL
        assert "visium" in _SPOT_SOURCES

    def test_a_visium_section_puts_no_number_on_the_ladder(self):
        from app.server.api import SPOT_SPATIAL_REFUSAL

        p = build([make_run("r1", "sec_x"), make_run("r2", "sec_v")],
                  sources={"sec_v": "visium"})
        refused = {r["section"]: r for r in p["reach"]["sections_refused"]}
        assert "sec_v" in refused
        assert refused["sec_v"]["kind"] == "platform"
        # VERBATIM. Not summarised, not truncated, not re-worded.
        assert refused["sec_v"]["reason"] == SPOT_SPATIAL_REFUSAL
        assert "sec_v" not in p["reach"]["sections_on_axis"]
        for row in p["reach"]["rows"]:
            assert "sec_v" not in row["under_reached_fraction"]

    def test_the_refusal_names_every_quantity_it_withholds(self):
        p = build([make_run("r1", "sec_v")], sources={"sec_v": "visium"})
        section = p["sections"][0]
        assert section["spot"]["per_cell_statistics_refused"] is True
        names = [q["quantity"] for q in section["spot"]["refused_quantities"]]
        assert "reach ladder F(x)" in names
        assert "penetration_length_um" in names
        assert section["E_ref"] is None

    def test_penetration_length_is_withheld_for_a_spot_section(self):
        """It is FIT on the vessel-distance profile the refusal names, so it
        cannot travel while the ladder it shares a fit with does not."""
        p = build([make_run("r1", "sec_x"), make_run("r2", "sec_v")],
                  sources={"sec_v": "visium"}, metric="penetration_length_um")
        assert p["absolute"]["comparable"] is False
        assert p["absolute"]["values"] == {}
        assert any("sec_v" in b for b in p["absolute"]["blockers"])

    def test_it_agrees_with_the_API_about_which_sections_those_are(self):
        """A MIRROR IS A SECOND COPY AND A SECOND COPY DRIFTS.

        ``app/server/api.py::_with_caveats`` is what the section picker uses.
        This view must refuse exactly the sections that picker refuses, so both
        are run over the same records and compared.
        """
        from app.server.api import _with_caveats

        for source in ("visium", "xenium", "ingested", ""):
            record = {"id": "s", "source": source}
            theirs = _with_caveats(dict(record))
            mine = xs._spot_verdict(record, {})
            assert mine["spatial_spread_valid"] == theirs["spatial_spread_valid"], source
            blocked = any(c["capability"] == "spatial_spread"
                          and c["status"] == "NO" for c in theirs["caveats"])
            assert mine["per_cell_statistics_refused"] == blocked, source

    def test_an_unnameable_platform_is_not_established_rather_than_valid(self):
        p = xs.cross_section_payload(
            FakeJobs([make_run("r1", "sec_q")]), compound="temozolomide",
            capability=lambda s: {"label": s, "available": True},
            section_record=lambda s: {})
        assert p["sections"][0]["spot"]["spatial_spread_valid"] is None
        assert p["sections"][0]["spot"]["per_cell_statistics_refused"] is False


# --------------------------------------------------------------------------- #
# 3. COMPARABILITY IS A BACKEND VERDICT
# --------------------------------------------------------------------------- #
class TestComparabilityIsNotADisplayChoice:

    def test_the_absolute_block_uses_projects_own_verdict_function(self):
        runs = [make_run("r1", "sec_a"), make_run("r2", "sec_b")]
        p = build(runs)
        entries = [
            {"specimen": {"id": s, "label": f"{s} label", "section": s},
             "capability": capability_for(s), "run": r}
            for s, r in zip(("sec_a", "sec_b"), runs)
        ]
        theirs = _kill_comparability(entries)
        assert p["absolute"]["comparable"] == theirs["comparable"]
        assert p["absolute"]["caveats"] == theirs["caveats"]
        assert p["absolute"]["state_measured"] == theirs["state_measured"]

    def test_a_refused_row_renders_no_numbers_at_all(self):
        """The refusal takes the place of the number. A value shown under a
        'NOT COMPARABLE' chip is still a value shown."""
        p = build([make_run("r1", "sec_a"), make_run("r2", "sec_b")],
                  caps={"sec_b": capability_for("sec_b", absolute_ok=False)})
        assert p["absolute"]["comparable"] is False
        assert p["absolute"]["withheld"] is True
        assert p["absolute"]["values"] == {}
        assert p["absolute"]["blockers"]

    def test_a_mixed_measured_and_assumed_state_blocks_the_row(self):
        p = build([make_run("r1", "sec_a", state_measured=True),
                   make_run("r2", "sec_b", state_measured=False)])
        assert p["absolute"]["comparable"] is False
        assert any("MEASURED" in b and "ASSUMED" in b
                   for b in p["absolute"]["blockers"])

    def test_absolute_kill_is_never_flagged_as_the_headline(self):
        p = build([make_run("r1", "sec_a"), make_run("r2", "sec_b")])
        assert p["absolute"]["is_headline"] is False
        assert "under-reached" in p["reading_order"][1] or "F(x)" in p["reading_order"][1]

    def test_an_unknown_metric_is_refused_rather_than_defaulted(self):
        with pytest.raises(xs.CrossSectionError) as exc:
            build([make_run("r1", "sec_a")], metric="auc_over_vibes")
        assert exc.value.status == 400


# --------------------------------------------------------------------------- #
# 4. THE DISAGREEMENT IS THE OUTPUT
# --------------------------------------------------------------------------- #
class TestWhereTheAnswerIsNotStableAcrossTissue:

    def test_curves_that_cross_produce_a_flip_verdict(self):
        """A compound more under-reached in tissue A at one rung and in tissue
        B at another cannot order those two tissues at all."""
        a = [0.30, 0.28, 0.20, 0.10, 0.05, 0.02, 0.01]
        b = [0.10, 0.12, 0.18, 0.22, 0.30, 0.40, 0.50]
        p = build([make_run("r1", "sec_a", fractions=a),
                   make_run("r2", "sec_b", fractions=b)])
        assert p["stability"]["any_order_flips"] is True
        assert p["stability"]["level"] == "crit"
        pair = p["stability"]["pairs"][0]
        assert pair["order_flips"] is True
        assert pair["higher_throughout"] is None
        assert pair["flip_rungs"]

    def test_curves_that_do_not_cross_hold_and_name_the_higher_section(self):
        a = [0.30, 0.28, 0.20, 0.10, 0.05, 0.02, 0.01]
        b = [0.20, 0.18, 0.10, 0.05, 0.02, 0.01, 0.005]
        p = build([make_run("r1", "sec_a", fractions=a),
                   make_run("r2", "sec_b", fractions=b)])
        assert p["stability"]["any_order_flips"] is False
        assert p["stability"]["pairs"][0]["higher_throughout"] == "sec_a"
        assert p["stability"]["level"] == "ok"

    def test_a_fraction_and_its_count_pointing_opposite_ways_is_reported(self):
        """results/heldout_patient.md §1.1's shape: the FRACTION falls while
        the count it is a fraction OF rises, because the denominator moved."""
        a = [0.20] * 7           # 20% of 1,000 viable  =   200 cells
        b = [0.10] * 7           # 10% of 9,000 viable  =   900 cells
        p = build([make_run("r1", "sec_a", fractions=a, n_viable=1000),
                   make_run("r2", "sec_b", fractions=b, n_viable=9000)])
        pair = p["stability"]["pairs"][0]
        assert p["stability"]["any_fraction_count_divergence"] is True
        assert len(pair["fraction_count_divergences"]) == 7
        one = pair["fraction_count_divergences"][0]
        assert one["fractions"]["sec_a"] > one["fractions"]["sec_b"]
        assert one["counts"]["sec_a"] < one["counts"]["sec_b"]

    def test_every_fraction_carries_its_own_denominator(self):
        p = build([make_run("r1", "sec_a", n_viable=1000, n_tumor=1200),
                   make_run("r2", "sec_b", n_viable=9000, n_tumor=11000)])
        row = p["reach"]["rows"][0]
        assert row["of_viable_tumor_cells"] == {"sec_a": 1000, "sec_b": 9000}
        assert row["under_reached_cells"]["sec_a"] == round(
            row["under_reached_fraction"]["sec_a"] * 1000)
        assert "round(F(x)" in row["cells_are_derived"]

    def test_the_denominator_caveat_travels_with_the_ladder(self):
        p = build([make_run("r1", "sec_a"), make_run("r2", "sec_b")])
        joined = " ".join(p["reach"]["caveats"])
        assert "heldout_patient.md" in joined
        assert "0.67x" in joined and "4.19x" in joined
        assert "cross_tissue.md" in joined

    def test_no_cross_section_noise_floor_is_claimed(self):
        p = build([make_run("r1", "sec_a"), make_run("r2", "sec_b")])
        note = p["stability"]["noise_floor"]
        assert "NO CROSS-SECTION NOISE FLOOR HAS BEEN MEASURED" in note
        assert "nothing here is a significance test" in note

    def test_the_convention_point_is_marked_but_not_headlined(self):
        p = build([make_run("r1", "sec_a"), make_run("r2", "sec_b")])
        marked = [r for r in p["reach"]["rows"] if r["is_convention_point"]]
        assert len(marked) == 1 and marked[0]["x"] == pytest.approx(0.1)
        assert "headline" not in p["stability"]["verdict"].lower()


# --------------------------------------------------------------------------- #
# 5. REFUSALS, NEVER DEFAULTS
# --------------------------------------------------------------------------- #
class TestRefusalsRatherThanDefaults:

    def test_two_window_modes_are_not_placed_on_one_axis(self):
        p = build([make_run("r1", "sec_a", window_mode="steady_state"),
                   make_run("r2", "sec_b", window_mode="finite_window")])
        assert p["reach"]["comparable"] is False
        assert any("window_mode" in b for b in p["reach"]["blockers"])
        assert p["stability"]["verdict"] == "REFUSED"

    def test_two_ladders_are_refused_rather_than_interpolated(self):
        p = build([make_run("r1", "sec_a"),
                   make_run("r2", "sec_b", ladder=[0.01, 0.1, 1.0],
                            fractions=[0.3, 0.2, 0.1])])
        assert p["reach"]["comparable"] is False
        assert any("interpolat" in b for b in p["reach"]["blockers"])

    def test_a_null_reach_is_carried_as_its_own_reason_not_as_zero(self):
        p = build([make_run("r1", "sec_a"),
                   make_run("r2", "sec_b", null_reason="no viable tumour cells")])
        refused = {r["section"]: r for r in p["reach"]["sections_refused"]}
        assert refused["sec_b"]["reason"] == "no viable tumour cells"
        assert refused["sec_b"]["kind"] == "null"

    def test_a_run_with_no_reach_map_says_so(self):
        p = build([make_run("r1", "sec_a"), make_run("r2", "sec_b", reach=None)])
        refused = {r["section"]: r for r in p["reach"]["sections_refused"]}
        assert refused["sec_b"]["kind"] == "absent"
        assert "predates" in refused["sec_b"]["reason"]

    def test_only_one_runspec_signature_is_placed_on_the_axis(self):
        """results/cross_tissue.md retracted its own cross-tissue spread
        because the sections had been run at different grids."""
        p = build([make_run("r1", "sec_a", dx_um=12.0),
                   make_run("r2", "sec_b", dx_um=12.0),
                   make_run("r3", "sec_c", dx_um=20.0)])
        assert [s["section"] for s in p["sections"]] == ["sec_a", "sec_b"]
        assert len(p["matched_sets"]) == 2
        assert sum(1 for s in p["matched_sets"] if s["is_primary"]) == 1
        assert "cross_tissue.md" in p["signature_note"]

    @pytest.mark.parametrize("order", ["oldest_first", "newest_first"])
    def test_a_second_run_of_one_section_is_reported_not_dropped_silently(self, order):
        """Both arrival orders. ``JobManager.list()`` is NEWEST FIRST, so the
        displaced run is usually the one arriving — the case a check written
        only in the replacing branch keeps silent."""
        runs = [make_run("r1", "sec_a", created_at=1.0),
                make_run("r2", "sec_a", created_at=2.0)]
        if order == "newest_first":
            runs.reverse()
        p = build(runs + [make_run("r3", "sec_b")])
        assert p["sections"][0]["run_id"] == "r2"
        assert any("r1" in problem for problem in p["problems"]), p["problems"]

    def test_one_section_is_not_a_comparison(self):
        p = build([make_run("r1", "sec_a")])
        assert p["reach"]["comparable"] is False
        assert p["absolute"]["withheld"] is True
        assert p["stability"]["verdict"] in ("REFUSED", "NOT A COMPARISON")


# --------------------------------------------------------------------------- #
# 6. THE ROUTES
# --------------------------------------------------------------------------- #
def _mount(app, router):
    """Apply the production mount line to an app that is already built.

    In ``create_app`` the line sits beside the other ``include_router`` calls,
    BEFORE ``_mount_web(app)`` — starlette matches in registration order and the
    UI is mounted at ``/``, so a router appended after it would be shadowed and
    every route below would 404. These tests therefore prepend, which is the
    same position the production line occupies.
    """
    if any(getattr(r, "path", "") == "/api/cross-section" for r in app.routes):
        return app
    app.router.routes = list(router.routes) + list(app.router.routes)
    return app


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    """A real app with this router mounted the way ``create_app`` will mount it.

    The production mount is ONE line inside ``app/server/api.py::create_app``,
    beside the other ``include_router`` calls::

        app.include_router(build_cross_section_router(jobs=jobs))

    That file belongs to another agent, so this fixture applies the same line to
    its own app instance: the routes below are exercised for real, and the mount
    is the only thing left to land.
    """
    from app.server.api import create_app
    from app.server.cross_section import build_cross_section_router
    from app.server.jobs import JobManager

    manager = JobManager(root=str(tmp_path_factory.mktemp("runs")))
    app = create_app(manager)
    _mount(app, build_cross_section_router(jobs=manager))
    with TestClient(app) as c:
        yield c
    manager.shutdown()


class TestTheRoute:
    """The routes, against a real app with the one-line mount applied."""

    def _get(self, client, path):
        return client.get(path)

    def test_no_compound_is_a_refusal_that_says_what_to_do_next(self, client):
        r = self._get(client, "/api/cross-section")
        assert r.status_code == 400
        assert "compounds" in str(r.json())

    def test_two_compounds_are_refused_with_the_withdrawal(self, client):
        r = self._get(client, "/api/cross-section?compound=a,b")
        assert r.status_code == 400
        assert "WITHDRAWN" in str(r.json())

    def test_an_unknown_metric_is_a_400(self, client):
        r = self._get(client, "/api/cross-section?compound=x&metric=vibes")
        assert r.status_code == 400

    def test_the_compound_list_is_a_set_never_an_ordering(self, client):
        r = self._get(client, "/api/cross-section/compounds")
        assert r.status_code == 200
        body = r.json()
        assert [c["compound"] for c in body["compounds"]] == sorted(
            c["compound"] for c in body["compounds"])
        assert "ALPHABETICAL" in body["ordering"]
        for row in body["compounds"]:
            assert set(row) >= {"compound", "n_sections"}
            assert not any(k in row for k in ("rank", "best", "score"))

    def test_it_is_gated_like_every_other_api_route(self, tmp_path, monkeypatch):
        from app.server.api import create_app
        from app.server.jobs import JobManager

        from app.server.cross_section import build_cross_section_router

        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "xs-key")
        manager = JobManager(root=str(tmp_path))
        app = create_app(manager)
        _mount(app, build_cross_section_router(jobs=manager))
        c = TestClient(app)
        # The gate is middleware on every /api/ path, so a router mounted after
        # create_app is authenticated by construction rather than by remembering.
        assert c.get("/api/cross-section/compounds").status_code == 401
        assert c.get("/api/cross-section?compound=x").status_code == 401
        ok = c.get("/api/cross-section/compounds",
                   headers={"x-api-key": "xs-key"})
        assert ok.status_code == 200


# --------------------------------------------------------------------------- #
# 7. THE BROWSER SIDE
# --------------------------------------------------------------------------- #
needs_node = pytest.mark.skipif(shutil.which("node") is None,
                                reason="node is not installed")


class TestTheViewCannotRenderTwoCompounds:

    def test_the_module_exists(self):
        assert os.path.exists(JS)

    @needs_node
    def test_it_parses_as_an_es_module(self):
        """`node --check` alone parses as CommonJS and has passed here on a
        real ESM error that broke the whole app."""
        with open(JS) as fh:
            src = fh.read()
        r = subprocess.run(["node", "--input-type=module", "--check"],
                           input=src, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr

    @needs_node
    def test_the_browser_guard_refuses_a_two_compound_payload(self, tmp_path):
        """The guard is executed, not grepped for."""
        harness = tmp_path / "guard.mjs"
        harness.write_text(
            f"const m = await import({_js_url()!r});\n"
            "const one = m.assertOneCompoundAxis("
            "{compound: 'tmz', sections: [{compound: 'tmz'}]});\n"
            "if (one !== 'tmz') { console.log('BAD:' + one); process.exit(1); }\n"
            "try {\n"
            "  m.assertOneCompoundAxis({compound: 'tmz', reach: {rows: "
            "[{deep: {drugs: ['doxorubicin']}}]}});\n"
            "  console.log('GUARD DID NOT FIRE'); process.exit(1);\n"
            "} catch (e) {\n"
            "  if (e.name !== 'CrossCompoundAxisRefused') { console.log('WRONG:'"
            " + e.name); process.exit(1); }\n"
            "  if (!/WITHDRAWN/.test(e.message)) { console.log('NO WITHDRAWAL');"
            " process.exit(1); }\n"
            "  console.log('GUARD FIRED ' + e.compounds.join(','));\n"
            "}\n")
        r = subprocess.run(["node", str(harness)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        assert "GUARD FIRED doxorubicin,tmz" in r.stdout

    def test_render_runs_the_guard_before_it_builds_a_row(self):
        src = open(JS).read()
        body = src.split("function render()")[1]
        assert "assertOneCompoundAxis(payload)" in body, (
            "render() does not call the axis guard. A guard that is defined, "
            "exported and tested but never CALLED on the render path is the "
            "exact defect this assertion exists to catch.")
        assert (body.index("assertOneCompoundAxis(payload)")
                < body.index("1 · Sections")), (
            "the guard runs after rows are built; it must refuse before the "
            "axis exists, not after part of it is on screen")

    def test_both_sides_of_the_wire_audit_the_same_keys(self):
        """A guard that runs only on the server is one a stale build walks
        past; two guards that disagree about what a compound key IS are one
        guard."""
        src = open(JS).read()
        found = re.search(r"const COMPOUND_KEYS = \[([^\]]+)\]", src)
        assert found, "the browser guard's key list could not be located"
        keys = tuple(k.strip().strip("'\"") for k in found.group(1).split(","))
        assert keys == tuple(xs._COMPOUND_KEYS)

    def test_it_defines_no_colour_of_its_own(self):
        """Workspace contract §7: every colour is a token, so the theme toggle
        works on this dialog too."""
        src = open(JS).read()
        block = src.split("const css = `")[1].split("`;")[0]
        assert not re.search(r"#[0-9a-fA-F]{3,8}\b", block)
        assert not re.search(r"\b(rgb|hsl)a?\(", block)
        tokens = set(re.findall(r"var\((--[a-z0-9-]+)\)", block))
        with open(os.path.join(_ROOT, "app", "web", "css", "tokens.css")) as fh:
            declared = fh.read()
        missing = [t for t in sorted(tokens) if f"{t}:" not in declared]
        assert not missing, f"not tokens: {missing}"

    def test_every_class_it_defines_is_prefixed(self):
        src = open(JS).read()
        block = src.split("const css = `")[1].split("`;")[0]
        classes = set(re.findall(r"\.([a-zA-Z][\w-]*)", block))
        assert all(c.startswith("xs-") for c in classes), sorted(classes)


# --------------------------------------------------------------------------- #
# `n_runs` is not the workspace over-count wearing the same clothes
# --------------------------------------------------------------------------- #
class TestNRunsIsNotTheWorkspaceOvercount:
    """``app/server/workspace.py`` summed each project's ``n_runs`` and reported
    **960 runs on a store holding 32**, because ``docs/WORKSPACE_CONTRACT.md``
    §2 lets one run id belong to any number of projects. ``available_compounds``
    has a ``len(slot["runs"])`` that LOOKS identical.

    It is not, and this class is the reason the next reader must not "fix" it
    into a bug. Two independent facts, each asserted separately so a failure
    says which one stopped being true.
    """

    def test_no_run_id_is_ever_in_two_slots(self):
        """Fact 1. A slot is keyed by ``drug`` and ``drug`` is a pure function
        of ONE run's spec, so the slots PARTITION the store — unlike projects,
        which own ids they did not create and which another project may own.

        Asserted through its one OBSERVABLE consequence rather than by
        re-deriving the slots here: if any run reached two slots, the rows
        would sum above the store's distinct run count. That over-count is
        precisely the 960-on-32 defect, so the equality below is the thing
        workspace.py could not have claimed.
        """
        views = [make_run("r_00001", "sec_a", drug="temozolomide"),
                 make_run("r_00002", "sec_b", drug="temozolomide"),
                 make_run("r_00003", "sec_a", drug="irinotecan"),
                 make_run("r_00004", "sec_b", drug="irinotecan"),
                 make_run("r_00005", "sec_c", drug="lomustine")]
        rows = xs.available_compounds(FakeJobs(views))["compounds"]

        # The premise, checked rather than assumed: more than one slot, or
        # "nothing was double-counted" is true for an uninteresting reason.
        assert len(rows) == 3, [r["compound"] for r in rows]
        assert {_drug_of(v) for v in views} == {r["compound"] for r in rows}

        distinct = len({v["run_id"] for v in views})
        assert sum(r["n_runs"] for r in rows) == distinct == 5, (
            f"rows sum to {sum(r['n_runs'] for r in rows)} on a store holding "
            f"{distinct} distinct run(s) -- a run reached more than one slot")

    def test_there_is_no_cross_slot_sum_to_be_wrong(self):
        """Fact 2. workspace.py's defect was in the SUM, not in the per-row
        ``n_runs`` — those rows were right there too. This function emits rows
        and no totals block, so the defective construct does not exist here."""
        out = xs.available_compounds(FakeJobs(
            [make_run("r_00001", "sec_a", drug="temozolomide")]))
        assert "totals" not in out
        assert set(out) == {"compounds", "ordering", "axis", "generated_at"}
        row_keys = set(out["compounds"][0])
        assert not (row_keys & {"total_runs", "runs_total"}), row_keys

    def test_a_run_seen_twice_is_counted_once(self):
        """The ``set`` is load-bearing and must survive a "simplification".

        ``JobManager.list()`` is ``[self._jobs[i] for i in
        reversed(self._order)]`` and ``_restore`` appends to ``_order`` using
        the ``run_id`` written INSIDE each ``result.json``. Two run directories
        declaring the same id therefore hand the SAME ``Job`` back twice. Built
        here on a REAL ``JobManager`` rather than a fake, because the claim is
        about that class's restore path and a fake would only restate it.
        """
        import json
        import tempfile

        from app.server.jobs import JobManager

        tmp = tempfile.mkdtemp()
        try:
            runs = os.path.join(tmp, "runs")
            for d in ("r_00001", "r_00002"):
                os.makedirs(os.path.join(runs, d))
                with open(os.path.join(runs, d, "result.json"), "w") as fh:
                    json.dump({"run_id": "r_00001", "status": "done",
                               "spec": {"section": "sec_a",
                                        "drug": "temozolomide"},
                               "metrics": {"kill_fraction": 0.1}}, fh)
            manager = JobManager(root=runs, start_worker=False)
            try:
                listed = manager.list()
                # The premise, checked rather than assumed: if list() did not
                # actually yield the run twice, this test proves nothing.
                assert [j.run_id for j in listed] == ["r_00001", "r_00001"], \
                    f"the duplicate was not constructed: {[j.run_id for j in listed]}"

                (row,) = xs.available_compounds(manager)["compounds"]
                assert row["n_runs"] == 1, (
                    "n_runs counted one run twice -- the accumulator stopped "
                    "being a set of ids")
                assert row["n_sections"] == 1
            finally:
                manager.shutdown()
        finally:
            shutil.rmtree(tmp)


def test_calibrate_alpha_is_unchanged_at_the_END_of_this_file():
    from tmesim.validation.calibration import calibrate_alpha

    assert calibrate_alpha() == ALPHA


def _js_url():
    return "file://" + JS
