"""Tests for ``GET /api/workspace`` (``app/server/workspace.py``).

The home page's single boot call, frozen by ``docs/WORKSPACE_CONTRACT.md`` §1.
Four things are checked, in rising order of what they would cost if they broke:

1. **It renders at all.** A fresh install has zero projects, zero runs and
   possibly zero cached sections. That is the FIRST screen anyone sees, and it
   must be a page that says "nothing here yet", not a 500. The same goes for a
   project with no runs, a run whose project was deleted, a run id that no
   longer resolves, and a project document that cannot be read.

2. **It says what is true.** ``last_opened_at`` is null, never now, for every
   project that predates the field (contract §3). A specimen label that does
   not resolve renders as a refusal that says what to do (contract §6), never
   as ``None`` and never as an exception. A thumbnail is a measured
   composition or an explained absence, never a fabricated picture (§4).

3. **It is FAST, and that is the point of the endpoint.** ``GET /api/sections``
   was measured at 12.7 s cold on a 0.1-CPU instance because it reads every
   section's metadata. This route is now the product's first paint, so three
   tests below are performance guards written as correctness guards: they
   break the expensive call and assert the payload is UNCHANGED. If anyone
   reintroduces ``list_sections()``, ``drug_rows()`` or a capability refresh
   into this path, they go red.

4. **It is gated.** Like every other ``/api/`` route, by being under ``/api/``.

No test here simulates anything. Runs are written as ``result.json`` documents
and loaded through ``JobManager``'s own restore path, so the whole file runs in
well under a second.
"""
from __future__ import annotations

import json
import os
import sys

import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

fastapi = pytest.importorskip("fastapi", reason="the app/server backend needs fastapi")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi.testclient import TestClient  # noqa: E402

from app.server import api as api_mod  # noqa: E402
from app.server import cache as cache_mod  # noqa: E402
from app.server import projects as projects_mod  # noqa: E402
from app.server import workspace as ws  # noqa: E402
from app.server.api import create_app  # noqa: E402
from app.server.jobs import JobManager  # noqa: E402
from app.server.projects import ProjectStore  # noqa: E402

#: Two ids the section REGISTRY always knows. Nothing here needs their data to
#: be cached: ``add_specimen`` only requires registry membership, and
#: ``capability_report`` degrades to "nothing is known" without a data_cache.
SECTION_A = "xenium_319"
SECTION_B = "xenium_io"

#: Every environment variable ``create_app`` reads. Cleared around any test
#: that sets one, because leaking one of these into a later module has broken
#: this suite twice (see tests/test_app_server.py's note).
_APP_ENV = (
    "TMESIM_STUDIO_API_KEY",
    "TMESIM_STUDIO_ALLOWED_HOSTS",
    "TMESIM_STUDIO_ALLOWED_ORIGINS",
    "TMESIM_STUDIO_DATA_ROOTS",
)

#: The contract's RunCard, exactly. Not a superset: a RunCard is deliberately
#: NOT a RunResult, and a field that creeps in here is a field the home page
#: will start rendering without anyone deciding that it should.
RUN_CARD_KEYS = {
    "run_id", "project_id", "project_name", "specimen_id", "specimen_label",
    "section", "drug", "status", "created_at", "finished_at", "runtime_s",
}

PAYLOAD_KEYS = {"projects", "recent_runs", "sections_available",
                "drugs_available", "totals", "server"}


# --------------------------------------------------------------------------- #
# Fixtures / builders
# --------------------------------------------------------------------------- #
@pytest.fixture
def runs_root(tmp_path):
    root = tmp_path / "runs"
    root.mkdir()
    return root


@pytest.fixture
def store(tmp_path):
    return ProjectStore(root=str(tmp_path / "projects"))


def write_run(runs_root, run_id, section=SECTION_A, drug="temozolomide", *,
              status="done", created_at=1000.0, started_at=1001.0,
              finished_at=1042.0, report=None):
    """A completed run as it exists on disk, so ``JobManager`` restores it.

    Deliberately goes through the real restore path rather than poking the
    manager's internals: a RunCard is built from a restored ``Job``, and a
    restored ``Job`` is the only kind this endpoint will ever see in
    production (the server is restarted far more often than a run is watched
    live).
    """
    d = runs_root / run_id
    d.mkdir(parents=True, exist_ok=True)
    doc = {
        "run_id": run_id,
        "status": status,
        "stage": status,
        "progress": 1.0 if status == "done" else 0.0,
        "created_at": created_at,
        "started_at": started_at,
        "finished_at": finished_at,
        "spec": {"section": section, "drug": drug, "dx_um": 30.0,
                 "total_time_s": 600.0},
        "metrics": {"mean_tumor_auc": 1.0e6},
    }
    if report is not None:
        doc["report"] = report
    (d / "result.json").write_text(json.dumps(doc))
    return run_id


def a_report(n_cells=40887):
    """A reconstruction report carrying a MEASURED composition — the only
    thing contract §4 allows a thumbnail to be built from."""
    return {"n_cells": n_cells,
            "celltype_fractions": {"TUMOR": 0.738, "ENDOTHELIAL": 0.0974,
                                   "IMMUNE": 0.1646}}


def manager(runs_root):
    return JobManager(root=str(runs_root), start_worker=False)


def attach(store, mgr, project_id, specimen_id, run_id):
    return store.attach_run(project_id, specimen_id, run_id,
                            mgr.get(run_id).result_view())


def payload(mgr, store):
    return ws.workspace_payload(mgr, store)


# --------------------------------------------------------------------------- #
# 1. It renders at all
# --------------------------------------------------------------------------- #
def test_empty_workspace_renders_rather_than_500(runs_root, store):
    """A fresh install: no projects, no runs, maybe no cached sections.

    This is the first screen of the whole product. It has to be a page.
    """
    out = payload(manager(runs_root), store)

    assert set(out) == PAYLOAD_KEYS
    assert out["projects"] == []
    assert out["recent_runs"] == []
    assert out["totals"] == {"projects": 0, "runs": 0, "specimens": 0}
    # The contract's two keys, pinned exactly as sec. 5 declares them.
    assert out["server"]["contract"] == "v1"
    assert out["server"]["stub"] is False
    # Plus the additive disclosure (contract amendment 2026-09-20). Asserted
    # rather than tolerated: a field the page uses to tell a user their work
    # will be destroyed must not be able to go missing silently.
    storage = out["server"]["storage"]
    assert set(storage) == {"persistent", "detail"}
    assert isinstance(storage["persistent"], bool)
    assert storage["detail"] and isinstance(storage["detail"], str)
    assert set(out["server"]) == {"contract", "stub", "storage"}
    # Counts are ints even when they are zero -- `null` here would make the
    # home page render "sections available: null".
    assert isinstance(out["sections_available"], int)
    assert isinstance(out["drugs_available"], int)
    assert out["sections_available"] >= 0
    assert out["drugs_available"] >= 0


def test_empty_workspace_over_http_is_200(tmp_path, monkeypatch):
    for var in _APP_ENV:
        monkeypatch.delenv(var, raising=False)
    mgr = JobManager(root=str(tmp_path / "runs"), start_worker=False)
    app = create_app(mgr, ProjectStore(root=str(tmp_path / "projects")))
    with TestClient(app) as c:
        r = c.get("/api/workspace")
    assert r.status_code == 200, r.text
    assert set(r.json()) == PAYLOAD_KEYS


def test_project_with_no_finished_runs_still_renders(runs_root, store):
    """A project can exist for days before anything finishes in it.

    Three states are checked at once: a specimen with no runs at all, a
    specimen whose only run is still queued, and what the card says about a
    composition it cannot draw.
    """
    mgr = manager(runs_root)
    p = store.create("has nothing yet")
    store.add_specimen(p["id"], SECTION_A, label="patient A")

    out = payload(mgr, store)
    assert len(out["projects"]) == 1
    assert out["projects"][0]["n_runs"] == 0
    assert out["recent_runs"] == []
    assert out["totals"] == {"projects": 1, "runs": 0, "specimens": 1}
    # Contract §4, built by ProjectStore: not a grey rectangle that reads as a
    # failed image, and never an empty chart.
    thumb = out["projects"][0]["thumbnail"]
    assert thumb["kind"] == "none"
    assert "no run" in thumb["reason"]
    # Computed, and the answer is none. Distinct from the key being absent.
    assert out["projects"][0]["compounds"] == []

    # Now a run that exists but has not finished.
    write_run(runs_root, "r_00001", status="queued", started_at=None,
              finished_at=None)
    mgr2 = manager(runs_root)
    spec_id = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]
    attach(store, mgr2, p["id"], spec_id, "r_00001")

    out = payload(mgr2, store)
    assert len(out["recent_runs"]) == 1
    card = out["recent_runs"][0]
    assert card["status"] == "queued"
    # A job that never started has no runtime. 0.0 would be a claim that it ran
    # instantly, which is the kind of number this repo refuses to invent.
    assert card["runtime_s"] is None
    assert card["finished_at"] is None
    thumb = out["projects"][0]["thumbnail"]
    assert thumb["kind"] == "none"
    # The store's refusal states the real reason rather than a bare "none".
    assert "queued" in thumb["reason"]
    assert out["projects"][0]["compounds"] == ["temozolomide"]


def test_a_deleted_projects_runs_never_appear(runs_root, store):
    """THE join guard. Deleting a project does NOT delete its runs (contract
    §2) -- they stay in app/runs/ and stay in the job manager. They must not
    keep appearing in a cross-project strip of a workspace they are no longer
    part of.
    """
    mgr = manager(runs_root)
    keep = store.create("kept")
    doomed = store.create("doomed")
    store.add_specimen(keep["id"], SECTION_A, label="kept patient")
    store.add_specimen(doomed["id"], SECTION_B, label="doomed patient")
    s_keep = store.get(keep["id"], refresh_capability=False)["specimens"][0]["id"]
    s_doom = store.get(doomed["id"], refresh_capability=False)["specimens"][0]["id"]

    write_run(runs_root, "r_00001", section=SECTION_A, created_at=1000.0,
              finished_at=1100.0, report=a_report())
    write_run(runs_root, "r_00002", section=SECTION_B, created_at=2000.0,
              finished_at=2100.0, report=a_report())
    mgr = manager(runs_root)
    attach(store, mgr, keep["id"], s_keep, "r_00001")
    attach(store, mgr, doomed["id"], s_doom, "r_00002")

    before = payload(mgr, store)
    assert {c["run_id"] for c in before["recent_runs"]} == {"r_00001", "r_00002"}

    deleted = store.delete(doomed["id"])
    # The run is retained -- that is the contract, and it is why this guard is
    # needed at all.
    assert deleted["runs_retained"] == ["r_00002"]
    assert mgr.get("r_00002") is not None

    after = payload(mgr, store)
    assert {c["run_id"] for c in after["recent_runs"]} == {"r_00001"}
    assert all(c["project_id"] != doomed["id"] for c in after["recent_runs"])
    assert after["totals"] == {"projects": 1, "runs": 1, "specimens": 1}


def test_a_dangling_run_id_is_skipped_not_fatal(runs_root, store):
    """A project document outlives the run directory it points at.

    The association is a plain id; nothing stops app/runs/ being pruned. A
    dangling id has no timestamp at all, so it cannot be placed in a strip
    ordered by time -- it is skipped, and it still counts in `totals.runs`,
    which reports what the project cards report.
    """
    mgr = manager(runs_root)
    p = store.create("pruned runs")
    store.add_specimen(p["id"], SECTION_A, label="patient A")
    sid = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]

    write_run(runs_root, "r_00001", report=a_report())
    mgr = manager(runs_root)
    attach(store, mgr, p["id"], sid, "r_00001")

    # The run directory goes away; the association does not.
    for f in (runs_root / "r_00001").iterdir():
        f.unlink()
    (runs_root / "r_00001").rmdir()
    mgr = manager(runs_root)
    assert mgr.get("r_00001") is None

    out = payload(mgr, store)
    assert out["recent_runs"] == []
    assert out["totals"]["runs"] == 1          # what the card says
    assert out["projects"][0]["n_runs"] == 1   # the card itself


def test_totals_runs_counts_a_shared_run_once(runs_root, store):
    """``totals.runs`` is how many runs EXIST, not how many references do.

    Contract §2: "a project owns nothing but ids". A run lives in ``app/runs/``
    and any number of projects may hold its id -- contract §1 makes that
    ordinary on purpose ("a run belongs to the store, and two specimens may
    reference the same run id"). Summing each project's ``n_runs`` therefore
    counts a shared run once per project holding it, and the header reports a
    workspace with more runs in it than the server has.

    Built WITHOUT ``duplicate()``, so what is covered is cross-project sharing
    itself and not one route's behaviour.
    """
    write_run(runs_root, "r_00001", report=a_report())
    mgr = manager(runs_root)

    for name in ("arm A", "arm B"):
        p = store.create(name)
        store.add_specimen(p["id"], SECTION_A, label="patient A")
        sid = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]
        attach(store, mgr, p["id"], sid, "r_00001")

    out = payload(mgr, store)

    # ONE run exists on this server, and both projects can reach it.
    assert sorted(os.listdir(runs_root)) == ["r_00001"]
    assert len(mgr.list()) == 1
    assert [r["n_runs"] for r in out["projects"]] == [1, 1]
    assert len({c["run_id"] for c in out["recent_runs"]}) == 1

    assert out["totals"]["runs"] == 1, (
        "the header counted one run once per project that holds its id")
    # Projects are not shared, and a specimen is a per-project record.
    assert out["totals"]["projects"] == 2
    assert out["totals"]["specimens"] == 2


def test_totals_runs_does_not_move_when_a_project_is_duplicated(runs_root,
                                                                store):
    """``duplicate`` is the route that makes sharing routine, and it is where
    this was found: duplicating a project holding two runs reported FOUR runs
    on a server holding two.

    "Runs are NOT re-executed and NOT copied" (contract §1). A number that
    doubles when nothing was created is not a count of runs.
    """
    for rid in ("r_00001", "r_00002"):
        write_run(runs_root, rid, report=a_report())
    mgr = manager(runs_root)

    p = store.create("original")
    store.add_specimen(p["id"], SECTION_A, label="patient A")
    sid = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]
    for rid in ("r_00001", "r_00002"):
        attach(store, mgr, p["id"], sid, rid)

    before = payload(mgr, store)
    assert before["totals"]["runs"] == 2

    store.duplicate(p["id"])
    after = payload(mgr, store)

    # Nothing was created on disk or in the manager.
    assert sorted(os.listdir(runs_root)) == ["r_00001", "r_00002"]
    assert len(mgr.list()) == 2
    assert after["totals"]["runs"] == 2, (
        "duplicating a project invented runs in the workspace header")
    # The CARDS are unchanged and are not wrong: each project really does hold
    # two run ids. The header is simply not their sum, and must not be.
    assert [r["n_runs"] for r in after["projects"]] == [2, 2]
    assert after["totals"]["specimens"] == 2


def test_totals_specimens_is_a_sum_and_has_to_stay_one(runs_root, store):
    """The mirror-image error, guarded so the fix above cannot travel to a
    number where distinct-counting would be WRONG.

    A specimen is not a shared object: it lives inside the project document,
    ``ProjectStore`` has no specimen store, and its id is allocated per project
    -- ``add_specimen`` takes max+1 of THAT project's ids. So ``s_00001``
    exists in very nearly every project on the server, and two projects'
    ``s_00001`` are two different patients. Contract §1 draws the line in the
    duplicate row itself: specimens are COPIED, runs are not.

    Counting distinct specimen ids would report one specimen for a workspace
    holding two.
    """
    a = store.create("unrelated A")
    b = store.create("unrelated B")
    store.add_specimen(a["id"], SECTION_A, label="patient A")
    store.add_specimen(b["id"], SECTION_B, label="patient B")

    sa = store.get(a["id"], refresh_capability=False)["specimens"][0]
    sb = store.get(b["id"], refresh_capability=False)["specimens"][0]
    # Two different patients on two different sections -- and the SAME id.
    assert sa["id"] == sb["id"] == "s_00001"
    assert sa["section"] != sb["section"]

    out = payload(manager(runs_root), store)
    assert out["totals"]["specimens"] == 2, (
        "two unrelated specimens were collapsed into one by their id")


def test_an_unreadable_project_document_does_not_500(runs_root, store):
    """``ProjectStore.list()`` reports a corrupt document as a row with
    ``unreadable`` set rather than dropping it. The workspace must carry that
    row through, count it, and say why it has no picture."""
    mgr = manager(runs_root)
    good = store.create("fine")
    bad = store.create("corrupt")
    path = os.path.join(store.root, bad["id"], "project.json")
    with open(path, "w") as fh:
        fh.write("{ this is not json")

    out = payload(mgr, store)
    assert out["totals"]["projects"] == 2
    rows = {p["id"]: p for p in out["projects"]}
    assert rows[good["id"]]["name"] == "fine"
    assert rows[bad["id"]]["unreadable"]
    # The store cannot build card art for a document it could not parse, and
    # this endpoint does not invent one.
    assert rows[bad["id"]]["thumbnail"] is None
    assert rows[bad["id"]]["last_opened_at"] is None
    # ABSENT, not []: this server cannot say what it has run.
    assert "compounds" not in rows[bad["id"]]


# --------------------------------------------------------------------------- #
# 2. It says what is true
# --------------------------------------------------------------------------- #
def test_a_project_that_vanishes_between_the_list_and_the_read_is_skipped(
        runs_root):
    """The workspace is a SNAPSHOT of a store two other routes are mutating.

    ``store.list()`` and ``store.get()`` are separate reads under separate
    acquisitions of the store lock, so a DELETE landing between them is
    ordinary, not exceptional -- and a home page that 500s because a project
    was deleted while it was loading is a home page that 500s.

    Reached deliberately: a corrupt document does NOT reach this branch
    (``store.list()`` reports it as an ``unreadable`` row and the loop stops
    there), so a test built on a corrupt file leaves this code untested while
    reporting green.
    """
    mgr = manager(runs_root)

    class _Vanishing:
        root = "/nowhere"

        def __init__(self, exc):
            self.exc = exc

        def list(self, thumbnails=True, resolve_run=None):
            return [{"id": "p_00001", "name": "here a moment ago",
                     "notes": "", "created_at": 10.0, "updated_at": 20.0,
                     "last_opened_at": None, "n_specimens": 2, "n_runs": 3,
                     "sections": [SECTION_A]},
                    {"id": "p_00002", "name": "still here", "notes": "",
                     "created_at": 5.0, "updated_at": 5.0,
                     "last_opened_at": None, "n_specimens": 0, "n_runs": 0,
                     "sections": []}]

        def get(self, pid, refresh_capability=True):
            if pid == "p_00001":
                raise self.exc
            return {"id": pid, "specimens": []}

    logged = []
    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(ws, "log_internal",
                        lambda exc, ctx: logged.append(ctx) or "ref")
    try:
        for exc, is_fault in (
                (projects_mod.ProjectError("no such project", status=404), False),
                (RuntimeError("the disk went away mid-read"), True)):
            logged.clear()
            out = payload(mgr, _Vanishing(exc))
            names = [p["name"] for p in out["projects"]]
            assert "still here" in names
            assert out["recent_runs"] == []
            # The row survives -- it is what `store.list()` said -- but this
            # server could not read its specimens, so it cannot claim a
            # compound set.
            gone = [p for p in out["projects"] if p["id"] == "p_00001"][0]
            assert "compounds" not in gone
            assert names.count("here a moment ago") == 1
            # Its row survives with its own `n_runs: 3`, but `totals.runs`
            # counts DISTINCT ids read off the documents -- and this project's
            # document could not be read, so it contributes none. The header
            # under-reports rather than guessing, and it does so for exactly
            # the project whose cards and compound set it has already dropped.
            assert gone["n_runs"] == 3
            assert out["totals"]["runs"] == 0
            # THE DIFFERENCE BETWEEN THE TWO BRANCHES. A deleted project is
            # traffic and is not logged; a fault removes a project's runs from
            # the page and MUST leave a record, or the page is wrong and
            # nothing anywhere says so.
            if is_fault:
                assert logged and "p_00001" in logged[0], (
                    "an unexpected read failure was swallowed without a "
                    "server-side record")
            else:
                assert logged == [], (
                    f"a project deleted mid-load logged {logged!r}; deleting a "
                    f"project is ordinary traffic, not a fault")
    finally:
        monkeypatch.undo()


def test_run_card_is_the_contracts_shape_and_no_more(runs_root, store):
    mgr = manager(runs_root)
    p = store.create("shape")
    store.add_specimen(p["id"], SECTION_A, label="patient A")
    sid = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]
    write_run(runs_root, "r_00001", drug="doxorubicin", report=a_report())
    mgr = manager(runs_root)
    attach(store, mgr, p["id"], sid, "r_00001")

    card = payload(mgr, store)["recent_runs"][0]
    assert set(card) == RUN_CARD_KEYS
    assert card["run_id"] == "r_00001"
    assert card["project_id"] == p["id"]
    assert card["project_name"] == "shape"
    assert card["specimen_id"] == sid
    assert card["specimen_label"] == "patient A"
    assert card["section"] == SECTION_A
    assert card["drug"] == "doxorubicin"
    assert card["status"] == "done"
    assert card["created_at"] == 1000.0
    assert card["finished_at"] == 1042.0
    assert card["runtime_s"] == pytest.approx(41.0)


def test_an_unresolvable_specimen_label_is_a_refusal_not_a_crash(
        runs_root, store):
    """The store will not let you blank a label. A document written by an
    older layout, restored from a backup or edited by hand can still arrive
    without one, and the home page must render a sentence that says what to
    do -- not ``None``, not a bare dash, and not a traceback."""
    mgr = manager(runs_root)
    p = store.create("hand edited")
    store.add_specimen(p["id"], SECTION_A, label="patient A")
    write_run(runs_root, "r_00001", report=a_report())
    mgr = manager(runs_root)
    sid = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]
    attach(store, mgr, p["id"], sid, "r_00001")

    path = os.path.join(store.root, p["id"], "project.json")
    with open(path) as fh:
        doc = json.load(fh)
    doc["specimens"][0].pop("label")
    with open(path, "w") as fh:
        json.dump(doc, fh)

    card = payload(mgr, store)["recent_runs"][0]
    label = card["specimen_label"]
    assert isinstance(label, str) and label.strip()
    assert sid in label                  # WHICH specimen
    assert SECTION_A in label            # what it is bound to
    assert "rename" in label             # what to do about it

    # And with the section gone too, which is the other half of "no longer
    # resolves" -- a record that cannot be run at all.
    doc["specimens"][0].pop("section")
    with open(path, "w") as fh:
        json.dump(doc, fh)
    label = payload(mgr, store)["recent_runs"][0]["specimen_label"]
    assert "incomplete" in label and sid in label


def test_last_opened_at_is_null_never_now(runs_root, store):
    """Contract §3. Every project that predates the field has never been
    opened by a human, and the only honest value is null -- which the card
    renders as "not opened yet". Filling it in with ``updated_at`` would make
    a background job that attached a run look like a person.

    TWO cases, and the second is the one this module is responsible for.
    ``ProjectStore._summary`` publishes the field today, so the first case
    tests the store's null and not this file's handling of it -- a guard that
    passes while testing someone else's code. The second gives this module a
    row with the key genuinely ABSENT, which is what every project document
    written before the field looks like, and is the only case
    ``workspace_payload`` itself decides.
    """
    mgr = manager(runs_root)
    store.create("never opened")
    row = payload(mgr, store)["projects"][0]
    assert "last_opened_at" in row
    assert row["last_opened_at"] is None
    # ... and specifically not smuggled in from either timestamp it sits next to
    assert row["updated_at"] is not None
    assert row["created_at"] is not None

    class _StoreWithoutTheField:
        """A summary row exactly as it looked before contract §3."""

        root = "/nowhere"

        def list(self, thumbnails=True, resolve_run=None):
            return [{"id": "p_00001", "name": "predates the field",
                     "notes": "", "created_at": 10.0, "updated_at": 20.0,
                     "n_specimens": 0, "n_runs": 0, "sections": []}]

        def get(self, pid, refresh_capability=True):
            return {"id": pid, "specimens": []}

    row = payload(mgr, _StoreWithoutTheField())["projects"][0]
    assert "last_opened_at" in row, (
        "the key must be present even when the store does not publish it -- "
        "the card distinguishes null from absent")
    assert row["last_opened_at"] is None, (
        f"a project that has never been opened reported "
        f"{row['last_opened_at']!r}; null is the only honest value and now is "
        f"the worst possible wrong one")


def test_a_store_that_publishes_last_opened_at_wins(runs_root, store):
    """``app/server/projects.py`` owns this field. When it starts publishing
    one, that value travels untouched and is NOT overwritten with null."""
    mgr = manager(runs_root)
    p = store.create("opened")
    path = os.path.join(store.root, p["id"], "project.json")
    with open(path) as fh:
        doc = json.load(fh)
    doc["last_opened_at"] = 5000.0
    with open(path, "w") as fh:
        json.dump(doc, fh)

    row = payload(mgr, store)["projects"][0]
    assert row["last_opened_at"] == 5000.0


def test_projects_sort_by_last_opened_then_updated(runs_root, store):
    """Contract §3: last_opened_at, falling back to updated_at when null.
    ``ProjectStore.list()`` sorts by ``created_at``, so this endpoint MUST
    re-sort -- the two orders are different and §3 is the load-bearing one."""
    mgr = manager(runs_root)
    rows = [
        {"id": "p_00001", "name": "oldest touch",  "updated_at": 100.0,
         "created_at": 900.0, "n_specimens": 0, "n_runs": 0},
        {"id": "p_00002", "name": "newest touch",  "updated_at": 300.0,
         "created_at": 100.0, "n_specimens": 0, "n_runs": 0},
        {"id": "p_00003", "name": "opened by hand", "updated_at": 200.0,
         "created_at": 200.0, "n_specimens": 0, "n_runs": 0,
         "last_opened_at": 400.0},
    ]

    class _Fake:
        root = "/nowhere"

        def list(self, thumbnails=True, resolve_run=None):
            return [dict(r) for r in rows]

        def get(self, pid, refresh_capability=True):
            return {"id": pid, "specimens": []}

    out = payload(mgr, _Fake())
    assert [p["name"] for p in out["projects"]] == [
        "opened by hand",   # last_opened_at 400 beats every updated_at
        "newest touch",     # updated_at 300
        "oldest touch",     # updated_at 100
    ]


def test_recent_runs_are_newest_first_and_capped_at_twelve(runs_root, store):
    """Contract §1: across ALL projects, newest first, <= 12."""
    mgr = manager(runs_root)
    p = store.create("many runs")
    store.add_specimen(p["id"], SECTION_A, label="patient A")
    sid = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]
    for i in range(1, 21):
        write_run(runs_root, f"r_{i:05d}", created_at=1000.0 + i,
                  started_at=1000.0 + i, finished_at=1000.0 + i + 0.5,
                  report=a_report())
    mgr = manager(runs_root)
    for i in range(1, 21):
        attach(store, mgr, p["id"], sid, f"r_{i:05d}")

    out = payload(mgr, store)
    assert len(out["recent_runs"]) == ws.MAX_RECENT_RUNS == 12
    finished = [c["finished_at"] for c in out["recent_runs"]]
    assert finished == sorted(finished, reverse=True)
    # The newest twelve, not the first twelve the store happened to list.
    assert out["recent_runs"][0]["run_id"] == "r_00020"
    assert out["recent_runs"][-1]["run_id"] == "r_00009"
    # The cap is on the STRIP, not on the count.
    assert out["totals"]["runs"] == 20


def test_recent_runs_span_projects(runs_root, store):
    """It is a FLATTENED cross-project view, not a per-project list."""
    mgr = manager(runs_root)
    a = store.create("project A")
    b = store.create("project B")
    store.add_specimen(a["id"], SECTION_A, label="A1")
    store.add_specimen(b["id"], SECTION_B, label="B1")
    sa = store.get(a["id"], refresh_capability=False)["specimens"][0]["id"]
    sb = store.get(b["id"], refresh_capability=False)["specimens"][0]["id"]
    write_run(runs_root, "r_00001", section=SECTION_A, finished_at=10.0)
    write_run(runs_root, "r_00002", section=SECTION_B, finished_at=20.0)
    mgr = manager(runs_root)
    attach(store, mgr, a["id"], sa, "r_00001")
    attach(store, mgr, b["id"], sb, "r_00002")

    cards = payload(mgr, store)["recent_runs"]
    assert [c["project_name"] for c in cards] == ["project B", "project A"]


def test_thumbnail_is_resolved_against_the_live_job_manager(runs_root, store):
    """Contract §4 is built by ``ProjectStore``, NOT by a second copy here.

    What this module is responsible for is the ``resolve_run`` hook. Passed
    one, the store reads the live JobManager; omitted, it reads
    ``app/runs/<id>/result.json`` from the PROCESS-WIDE runs directory -- the
    real one, not this test's. So the n_cells assertion below is the guard:
    drop the resolver and the thumbnail comes back from whatever run happens to
    be sitting in the developer's own app/runs/.

    The store picks the FIRST run that can supply a composition, not the
    newest, deliberately: a card that repaints every time a background job
    finishes is a card that flickers for no information."""
    mgr = manager(runs_root)
    p = store.create("has a composition")
    store.add_specimen(p["id"], SECTION_A, label="patient A")
    sid = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]
    write_run(runs_root, "r_00001", created_at=10.0, finished_at=20.0,
              report=a_report(n_cells=111))
    write_run(runs_root, "r_00002", created_at=30.0, finished_at=40.0,
              report=a_report(n_cells=222))
    # No report at all: cannot support a thumbnail, must not be chosen even
    # though it is the newest.
    write_run(runs_root, "r_00003", created_at=50.0, finished_at=60.0)
    mgr = manager(runs_root)
    for rid in ("r_00001", "r_00002", "r_00003"):
        attach(store, mgr, p["id"], sid, rid)

    thumb = payload(mgr, store)["projects"][0]["thumbnail"]
    assert thumb["kind"] == "composition"
    assert thumb["section"] == SECTION_A
    # 111 is THIS test's first composition-capable run. 40887 would mean the
    # store fell back to the real app/runs/ because the resolver went missing.
    assert thumb["n_cells"] == 111
    assert thumb["type_fractions"]["TUMOR"] == pytest.approx(0.738)
    # No ordering of anything, and no kill fraction anywhere near a card.
    assert "kill_fraction" not in json.dumps(thumb)


def test_compounds_is_a_deduplicated_alphabetical_set(runs_root, store):
    """The compounds a project has RUN, as a SET.

    Cross-compound ordering by predicted effect was withdrawn 2026-08-21 and
    the prohibition is absolute -- no metric, no "top", no badge. So the only
    order this field may carry is the alphabet, and the card says so on screen.
    This test pins the alphabet, which is what makes "not ordered by anything
    else" checkable at all.
    """
    mgr = manager(runs_root)
    p = store.create("three compounds, four runs")
    store.add_specimen(p["id"], SECTION_A, label="patient A")
    sid = store.get(p["id"], refresh_capability=False)["specimens"][0]["id"]
    for i, drug in enumerate(
            ("temozolomide", "doxorubicin", "temozolomide", "Carboplatin"), 1):
        write_run(runs_root, f"r_{i:05d}", drug=drug, created_at=1000.0 + i,
                  finished_at=1000.0 + i)
    mgr = manager(runs_root)
    for i in range(1, 5):
        attach(store, mgr, p["id"], sid, f"r_{i:05d}")

    compounds = payload(mgr, store)["projects"][0]["compounds"]
    # Deduplicated (temozolomide was run twice) and alphabetical regardless of
    # case -- "Carboplatin" must not sort ahead of everything because it is
    # capitalised.
    assert compounds == ["Carboplatin", "doxorubicin", "temozolomide"]
    # The runs themselves were attached in a different order, so this is not
    # attachment order wearing a hat.
    assert [c["drug"] for c in payload(mgr, store)["recent_runs"]] != compounds


def test_compounds_absent_and_empty_are_different_states(runs_root, store):
    """[] means "computed, and nothing has been run". An ABSENT key means
    "this server cannot tell you". The card renders them differently and it is
    right to: collapsing them would turn "I do not know" into "there is none".

    A partial answer is treated as not knowing. A project with one pruned run
    genuinely does not know its own compound set, and a list of the four it
    CAN see looks exactly like a complete list of four.
    """
    mgr = manager(runs_root)
    empty = store.create("nothing run yet")
    store.add_specimen(empty["id"], SECTION_A, label="A")
    partial = store.create("one run pruned")
    store.add_specimen(partial["id"], SECTION_B, label="B")
    sid = store.get(partial["id"],
                    refresh_capability=False)["specimens"][0]["id"]

    write_run(runs_root, "r_00001", section=SECTION_B, drug="temozolomide")
    write_run(runs_root, "r_00002", section=SECTION_B, drug="doxorubicin")
    mgr = manager(runs_root)
    attach(store, mgr, partial["id"], sid, "r_00001")
    attach(store, mgr, partial["id"], sid, "r_00002")

    rows = {p["id"]: p for p in payload(mgr, store)["projects"]}
    assert rows[empty["id"]]["compounds"] == []
    assert rows[partial["id"]]["compounds"] == ["doxorubicin", "temozolomide"]

    # Now prune one run document. The association survives; the run does not.
    for f in (runs_root / "r_00002").iterdir():
        f.unlink()
    (runs_root / "r_00002").rmdir()
    mgr = manager(runs_root)

    rows = {p["id"]: p for p in payload(mgr, store)["projects"]}
    assert rows[empty["id"]]["compounds"] == []
    assert "compounds" not in rows[partial["id"]], (
        "a project with an unresolvable run does not know its own compound "
        "set and must not publish a partial one as if it were complete")


# --------------------------------------------------------------------------- #
# 3. Performance guards, written as correctness guards
# --------------------------------------------------------------------------- #
def test_cheap_section_count_matches_full_enumeration():
    """THE DRIFT GUARD for the mirrored dispatch in ``workspace._available``.

    ``sections_available()`` deliberately runs only the availability half of
    ``describe_section`` -- 18-33x cheaper, because it stats filenames instead
    of reading every panel's gene list. That makes it a SECOND COPY of the
    dispatch rule, and a second copy is correct on the day it is written.
    """
    full = {s["id"]: bool(s["available"]) for s in cache_mod.list_sections()}
    cheap = {}
    for sid, rec in cache_mod.all_sections().items():
        try:
            cheap[sid] = bool(ws._available(rec)[0])
        except Exception:  # noqa: BLE001
            cheap[sid] = False

    assert set(cheap) == set(full)
    disagree = {k for k in full if full[k] != cheap[k]}
    assert not disagree, (
        f"the cheap availability dispatch disagrees with describe_section on "
        f"{sorted(disagree)}")
    assert ws.sections_available() == sum(1 for v in full.values() if v)


def test_a_section_whose_data_is_absent_is_not_counted(monkeypatch, tmp_path):
    """The half of the drift guard above that the real registry cannot test.

    Every section on a developer's machine with a full ``data_cache/`` is
    available, so "cheap == full" is satisfied by an implementation that
    simply returns True -- a check that passes while testing nothing, which is
    this repo's most-repeated defect. So the registry is given sections whose
    data is NOT there, one per branch of the dispatch, and the count must not
    move.
    """
    real = cache_mod.all_sections()
    missing = os.path.join(str(tmp_path), "not-a-cache-dir")
    assert not os.path.exists(missing)

    fake = dict(real)
    fake["probe_xenium_absent"] = {
        "id": "probe_xenium_absent", "label": "absent xenium", "source": "xenium",
        "sample": "probe", "cache_dir": missing, "panel": None,
    }
    fake["probe_visium_absent"] = {
        "id": "probe_visium_absent", "label": "absent visium", "source": "visium",
        "sample": "no-such-visium-sample", "cache_dir": missing, "panel": None,
    }
    fake["probe_ingested_absent"] = {
        "id": "probe_ingested_absent", "label": "absent ingested",
        "source": "xenium", "sample": "probe", "cache_dir": missing,
        "panel": None, "ingested": True, "platform": "xenium", "runnable": True,
    }
    # Both namespaces: `cache.list_sections` reads the module global, and
    # `workspace.sections_available` reads the name bound at import.
    monkeypatch.setattr(cache_mod, "all_sections", lambda: fake)
    monkeypatch.setattr(ws, "all_sections", lambda: fake)

    before = sum(1 for r in real.values()
                 if ws._available(r)[0])
    assert ws.sections_available() == before, (
        "three sections whose data is not on this disk were counted as "
        "runnable")

    full = {s["id"]: bool(s["available"]) for s in cache_mod.list_sections()}
    for probe in ("probe_xenium_absent", "probe_visium_absent",
                  "probe_ingested_absent"):
        assert probe in full
        assert full[probe] is False, f"{probe} should be unavailable"
        assert ws._available(fake[probe])[0] is False, (
            f"the cheap dispatch calls {probe} runnable and describe_section "
            f"does not")


def test_workspace_never_enumerates_section_metadata(runs_root, store,
                                                     monkeypatch):
    """PERFORMANCE GUARD. ``_section_metadata`` is the expensive half -- it is
    what made ``GET /api/sections`` 12.7 s cold on a 0.1-CPU instance.

    Break it and the count must not move. If anyone rewrites
    ``sections_available()`` in terms of ``list_sections()``, every section
    reports unavailable (``describe_section`` catches a metadata failure and
    flips ``available`` to False) and this assertion fails.
    """
    expected = ws.sections_available()

    def boom(*a, **k):
        raise AssertionError("_section_metadata must not be on the boot path")

    monkeypatch.setattr(cache_mod, "_section_metadata", boom)
    monkeypatch.setattr(cache_mod, "_meta_cache", None)

    out = payload(manager(runs_root), store)
    assert out["sections_available"] == expected


def test_workspace_never_compiles_the_drug_library(runs_root, store,
                                                   monkeypatch):
    """PERFORMANCE GUARD. ``drug_rows()`` compiles ten molecules through the
    calibration chain: MEASURED 6.5 s cold here, ~15 s on the deploy target,
    which is why ``create_app`` warms it on a background thread. The home page
    needs a COUNT, and ``_known_drugs()`` produces the same names in 0.114 s
    cold / 0.000 s warm.
    """
    expected = len(api_mod._known_drugs())
    assert expected > 0

    def boom(*a, **k):
        raise AssertionError("drug_rows() must not be on the boot path")

    monkeypatch.setattr(api_mod, "drug_rows", boom)
    out = payload(manager(runs_root), store)
    assert out["drugs_available"] == expected


def test_workspace_never_refreshes_specimen_capability(runs_root, store,
                                                       monkeypatch):
    """PERFORMANCE GUARD. ``ProjectStore.get()`` re-derives every specimen's
    capability report by default, which reads that section's panel gene list
    per specimen. Nothing on a project card renders one, so every read on this
    path passes ``refresh_capability=False``.
    """
    p = store.create("capability")
    store.add_specimen(p["id"], SECTION_A, label="patient A")
    store.add_specimen(p["id"], SECTION_B, label="patient B")

    calls = []
    real = projects_mod.capability_report
    monkeypatch.setattr(projects_mod, "capability_report",
                        lambda sid: (calls.append(sid), real(sid))[1])

    out = payload(manager(runs_root), store)
    assert out["totals"]["specimens"] == 2
    assert calls == [], f"capability_report was called for {calls}"


# --------------------------------------------------------------------------- #
# 4. It is gated
# --------------------------------------------------------------------------- #
def test_every_workspace_route_is_gated(tmp_path, monkeypatch):
    """The gate keys off ``get_route_path(...).startswith("/api/")``, so this
    router is authenticated by construction. Asserted rather than assumed --
    the gate is the reason a remote deployment is not an open store, and
    "it should be covered by the prefix" is exactly the kind of reasoning this
    repo has been wrong about before.

    EVERY route the router mounts, not just the one it is named after: the two
    project mutations MUTATE, and an unauthenticated POST that duplicates or
    restamps someone's project is a worse hole than an unauthenticated read.
    """
    for var in _APP_ENV:
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "workspace-gate-key")
    mgr = JobManager(root=str(tmp_path / "runs"), start_worker=False)
    store = ProjectStore(root=str(tmp_path / "projects"))
    p = store.create("gated")
    app = create_app(mgr, store)
    key = {"x-api-key": "workspace-gate-key"}
    routes = [("GET", "/api/workspace"),
              ("POST", f"/api/projects/{p['id']}/open"),
              ("POST", f"/api/projects/{p['id']}/duplicate")]
    with TestClient(app) as c:
        for method, path in routes:
            assert c.request(method, path).status_code == 401, path
            assert c.request(method, path,
                             headers={"x-api-key": "wrong"}).status_code == 401
        ok = c.get("/api/workspace", headers=key)
        assert ok.status_code == 200
        assert set(ok.json()) == PAYLOAD_KEYS
        assert c.post(f"/api/projects/{p['id']}/open", headers=key,
                      json={}).status_code == 200
        assert c.post(f"/api/projects/{p['id']}/duplicate", headers=key,
                      json={}).status_code == 201


# --------------------------------------------------------------------------- #
# 5. The home page's two project mutations
# --------------------------------------------------------------------------- #
@pytest.fixture
def http(tmp_path, monkeypatch):
    """A TestClient with an isolated store AND runs dir, and no API key."""
    for var in _APP_ENV:
        monkeypatch.delenv(var, raising=False)
    runs = tmp_path / "runs"
    mgr = JobManager(root=str(runs), start_worker=False)
    store = ProjectStore(root=str(tmp_path / "projects"))
    app = create_app(mgr, store)
    with TestClient(app) as c:
        c._store = store      # type: ignore[attr-defined]
        c._manager = mgr      # type: ignore[attr-defined]
        c._runs = runs        # type: ignore[attr-defined]
        yield c


def test_open_stamps_recency_without_touching_updated_at(http):
    """Contract §3, and the whole reason ``last_opened_at`` is a second field.

    ``updated_at`` answers "did the project DOCUMENT change". Opening is not a
    change. Bump both and a workspace sorted by recency puts a project on top
    because a background job attached a run to it, which is precisely the
    behaviour the field was added to stop.
    """
    created = http.post("/api/projects", json={"name": "recency"}).json()
    pid = created["id"]
    assert created["last_opened_at"] is None

    before = http.get(f"/api/projects/{pid}").json()
    r = http.post(f"/api/projects/{pid}/open", json={})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["ok"] is True
    assert isinstance(body["last_opened_at"], float)
    assert body["updated_at"] == before["updated_at"], (
        "opening a project moved updated_at; opening is not a change to the "
        "document and contract §3 is the distinction")

    after = http.get(f"/api/projects/{pid}").json()
    assert after["last_opened_at"] == body["last_opened_at"]
    assert after["updated_at"] == before["updated_at"]


def test_open_is_idempotent_and_reaches_the_workspace(http):
    pid = http.post("/api/projects", json={"name": "twice"}).json()["id"]
    first = http.post(f"/api/projects/{pid}/open", json={}).json()
    second = http.post(f"/api/projects/{pid}/open", json={}).json()
    assert second["ok"] is True
    assert second["last_opened_at"] >= first["last_opened_at"]

    row = [p for p in http.get("/api/workspace").json()["projects"]
           if p["id"] == pid][0]
    assert row["last_opened_at"] == second["last_opened_at"]


def test_open_sorts_a_project_to_the_top_of_the_workspace(http):
    """The end-to-end point of the whole mechanism: opening a project is what
    puts it at the top of the home page, and nothing else does."""
    old = http.post("/api/projects", json={"name": "opened by a human"}).json()
    new = http.post("/api/projects", json={"name": "created later"}).json()

    # Newest-created is on top while nothing has been opened.
    assert [p["name"] for p in http.get("/api/workspace").json()["projects"]] \
        == ["created later", "opened by a human"]

    http.post(f"/api/projects/{old['id']}/open", json={})
    assert [p["name"] for p in http.get("/api/workspace").json()["projects"]] \
        == ["opened by a human", "created later"]


def test_open_on_a_project_that_is_gone_is_a_404_with_a_sentence(http):
    """Not a 500, not a traceback, and not a 200.

    The Studio fires this on boot and swallows the rejection, because recency
    is an annotation. That swallow is a CLIENT decision; this route must not
    make it for every client by answering 200 to an id that does not exist --
    "recorded" and "there was nothing to record" would become the same answer.
    """
    r = http.post("/api/projects/p_99999/open", json={})
    assert r.status_code == 404, r.text
    assert "p_99999" in json.dumps(r.json())

    # A malformed id must not reach the filesystem as a path segment either.
    for bad in ("..", "%2e%2e", "p_1", "not-an-id"):
        assert http.post(f"/api/projects/{bad}/open",
                         json={}).status_code in (404, 405), bad


def test_open_and_duplicate_accept_no_body_at_all(http):
    """Both take an optional body, so an ABSENT one is a well-formed request
    -- but a body that is present and is not a JSON object still 400s exactly
    as it does on every other route."""
    pid = http.post("/api/projects", json={"name": "bodyless"}).json()["id"]
    assert http.post(f"/api/projects/{pid}/open").status_code == 200
    assert http.post(f"/api/projects/{pid}/duplicate").status_code == 201
    for bad in ("[]", '"a string"', "17", "{not json"):
        r = http.post(f"/api/projects/{pid}/open", content=bad,
                      headers={"content-type": "application/json"})
        assert r.status_code == 400, (bad, r.status_code)


def test_duplicate_copies_specimens_and_run_ids_but_runs_are_not_copied(
        tmp_path, monkeypatch):
    """Contract §1 and §2. A run belongs to ``app/runs/``; a project owns
    nothing but ids. After a duplicate, two specimens in two projects point at
    the SAME run id -- that is what makes the copy instant, and it is why
    deleting either project still deletes no run.

    The project is given REAL attached runs first. Duplicating a project with
    no runs proves nothing about whether runs are copied, and a test that
    proves nothing while passing is the defect this file is most on guard
    against.
    """
    for var in _APP_ENV:
        monkeypatch.delenv(var, raising=False)
    runs = tmp_path / "runs"
    runs.mkdir()
    write_run(runs, "r_00001", section=SECTION_A, drug="temozolomide",
              report=a_report())
    write_run(runs, "r_00002", section=SECTION_A, drug="doxorubicin",
              report=a_report())
    mgr = JobManager(root=str(runs), start_worker=False)
    store = ProjectStore(root=str(tmp_path / "projects"))
    app = create_app(mgr, store)

    with TestClient(app) as http:
        src = http.post("/api/projects", json={"name": "original"}).json()
        http.post(f"/api/projects/{src['id']}/specimens",
                  json={"section": SECTION_A, "label": "patient A"})
        doc = http.get(f"/api/projects/{src['id']}").json()
        sid = doc["specimens"][0]["id"]
        for rid in ("r_00001", "r_00002"):
            att = http.post(
                f"/api/projects/{src['id']}/specimens/{sid}/runs",
                json={"run_id": rid})
            assert att.status_code in (200, 201), att.text
        doc = http.get(f"/api/projects/{src['id']}").json()
        assert doc["specimens"][0]["run_ids"] == ["r_00001", "r_00002"]

        projects_before = set(os.listdir(store.root))
        runs_before = sorted(os.listdir(runs))
        n_jobs_before = len(mgr.list())
        assert runs_before == ["r_00001", "r_00002"]
        assert n_jobs_before == 2

        r = http.post(f"/api/projects/{src['id']}/duplicate", json={})
        assert r.status_code == 201, r.text
        copy = r.json()
        assert copy["id"] != src["id"]
        assert copy["name"] == "original (copy)"
        assert copy["duplicated_from"] == src["id"]
        # Nobody has opened the copy.
        assert copy["last_opened_at"] is None
        # Specimens copied verbatim, run REFERENCES included -- the same ids,
        # now reachable from two projects.
        assert [s["id"] for s in copy["specimens"]] == [sid]
        assert copy["specimens"][0]["label"] == "patient A"
        assert copy["specimens"][0]["section"] == SECTION_A
        assert copy["specimens"][0]["run_ids"] == ["r_00001", "r_00002"]

        # NO RUN WAS CREATED, COPIED OR RE-EXECUTED. Checked three ways,
        # because this is the clause of contract §2 that the word "duplicate"
        # most invites a reader to assume the opposite of.
        assert sorted(os.listdir(runs)) == runs_before
        assert len(mgr.list()) == n_jobs_before
        assert set(os.listdir(store.root)) - projects_before == {copy["id"]}

        # The original is untouched. Compared on the fields a duplicate could
        # plausibly disturb -- not on the whole specimen, whose `capability`
        # block is re-derived on every read and is not what this test is about.
        after = http.get(f"/api/projects/{src['id']}").json()["specimens"]
        keys = ("id", "label", "section", "run_ids", "created_at")
        assert [tuple(x[k] for k in keys) for x in after] == \
               [tuple(x[k] for k in keys) for x in doc["specimens"]]

        # Both projects in the workspace, both counted, and the SAME two runs
        # now appear once per project in the cross-project strip -- which is
        # the honest rendering of one run reachable from two places.
        ws_out = http.get("/api/workspace").json()
        assert ws_out["totals"]["projects"] == 2
        assert ws_out["totals"]["specimens"] == 2
        # TWO, not four. This assertion read `== 4` until 2026-09-20 and was
        # green over the defect it should have caught: the header summed each
        # project's `n_runs`, so duplicating a project doubled the workspace's
        # run count while the three assertions above proved no run had been
        # created. `totals.runs` counts DISTINCT ids -- see
        # test_totals_runs_does_not_move_when_a_project_is_duplicated.
        assert ws_out["totals"]["runs"] == len(runs_before) == 2
        # The per-project rows are unchanged and are NOT wrong: each project
        # really does hold both ids. The header is simply not their sum.
        assert [row["n_runs"] for row in ws_out["projects"]] == [2, 2]
        by_project = {}
        for card in ws_out["recent_runs"]:
            by_project.setdefault(card["project_id"], set()).add(card["run_id"])
        assert by_project == {src["id"]: {"r_00001", "r_00002"},
                              copy["id"]: {"r_00001", "r_00002"}}
        for row in ws_out["projects"]:
            assert row["compounds"] == ["doxorubicin", "temozolomide"]


def test_duplicate_takes_a_name_and_refuses_an_empty_one(http):
    src = http.post("/api/projects", json={"name": "original"}).json()
    named = http.post(f"/api/projects/{src['id']}/duplicate",
                      json={"name": "arm B"})
    assert named.status_code == 201
    assert named.json()["name"] == "arm B"
    blank = http.post(f"/api/projects/{src['id']}/duplicate",
                      json={"name": "   "})
    assert blank.status_code == 400, blank.text
    assert http.post("/api/projects/p_99999/duplicate",
                     json={}).status_code == 404
