"""``GET /api/workspace`` — the home page's single boot call.

``docs/WORKSPACE_CONTRACT.md`` §1 freezes the payload. This module builds it
and nothing else.

Why one call
------------

The home page needs four unrelated things: the projects, a cross-project strip
of recent runs, and the two counts that tell a new arrival whether this server
has any data at all. Four round-trips on a cold 0.1-CPU instance is four
chances to look broken, and this page is now the product's first paint.

Why it is not ``GET /api/sections`` with extras
-----------------------------------------------

**MEASURED, and this is the whole design of the file.** ``list_sections()``
calls ``describe_section`` for every registered section, and that reads
``features.tsv.gz`` / the per-cell centroid table for each one to fill in
``n_cells``, ``extent_um``, ``n_genes`` and the panel name. On a warm
``app/runs/_section_meta.json`` memo that is ~0.03 s; with the memo cold it is
0.84 s here for 38 sections, and ``GET /api/sections`` was measured at 12.7 s
cold on the 0.1-CPU instance this deploys to.

The home page does not render any of those numbers. It renders ONE integer:
how many sections this server can run. That integer needs the availability
half of ``describe_section`` and none of the metadata half, and the
availability half is 18-33x cheaper because it only stats filenames:

    availability-only  0.0266 s   (38 sections)
    list_sections()    0.8782 s   (38 sections, memo cold)   -> 33.0x

So :func:`sections_available` mirrors ``describe_section``'s dispatch and calls
only the ``_*_available`` half.

A MIRROR IS A SECOND COPY, AND A SECOND COPY DRIFTS. That is a defect this repo
has shipped before, so the mirror is not trusted — it is CHECKED.
``tests/test_workspace.py::test_cheap_section_count_matches_full_enumeration``
runs both paths over the real registry and fails if they disagree on a single
id. If ``describe_section`` grows a fourth source, or reorders its dispatch,
the suite goes red rather than this endpoint going quietly wrong.

What is NOT read here
---------------------

* ``drug_rows()`` — compiles ten molecules through the calibration chain.
  MEASURED 6.5 s cold here (~15 s on the deploy target, which is why
  ``create_app``'s lifespan warms it on a background thread). ``_known_drugs()``
  returns the same names — the exact set a ``RunSpec``'s ``drug`` field may
  take — in 0.114 s cold and 0.000 s warm. The home page shows a count, not a
  library.
* ``ProjectStore.get(refresh_capability=True)`` — re-derives every specimen's
  capability report, which reads each section's panel gene list. Every read
  here passes ``refresh_capability=False``. A capability block is not on a
  project card.
* ``Job.result_view()`` — copies the whole run body (report, metrics,
  timelapse, provenance), backfills the delivery law and recomputes the
  claim-set verdict. A ``RunCard`` needs seven scalars off the ``Job`` object
  itself, so nothing here calls it, and neither does the thumbnail resolver.

What is NOT built here
----------------------

The ``thumbnail`` (contract §4) is built by ``ProjectStore`` and published on
every row of ``store.list()``. This module supplies the ``resolve_run`` hook
that points that machinery at the LIVE job manager instead of at
``app/runs/<id>/result.json``, and otherwise passes the block through
untouched. Two implementations of "what may be drawn on a card" is how a rule
that says *never fabricate a picture of tissue* gets weakened in the copy
nobody remembered to update.

``compounds`` IS built here, because ``ProjectSummary`` has no such field and
the project card needs one. Three states, deliberately:

    key absent    this server cannot say — the project document could not be
                  read, or an attached run could not be resolved. A PARTIAL
                  set is not the set, and a list of the four compounds this
                  server can see looks exactly like a complete list of four.
    []            computed, and nothing has been run yet.
    [names...]    a SET: deduplicated, alphabetical, ordered by NOTHING else.
                  Cross-compound ordering by predicted effect was withdrawn
                  2026-08-21 and the prohibition is absolute.

Caching
-------

Nothing in this module is cached, so nothing in it can go stale. The only
memo on this path belongs to ``ProjectStore``, which keys a run's composition
on that document's ``(st_mtime_ns, st_size)`` — and this module bypasses even
that by handing the store the live JobManager instead.

The cold cost that remained was not per-request work at all but two one-time
module imports, and those are warmed on a daemon thread at app-build time
rather than cached. See :func:`_warm` for the measurement and for what that
thread deliberately does NOT warm.
"""
from __future__ import annotations

import inspect
import threading
from typing import Callable, Dict, List, Optional, Tuple

from fastapi import APIRouter, Request

from .cache import (
    all_sections,
    log_internal,
    # The availability half of ``describe_section``'s dispatch. Imported
    # rather than reimplemented so the RULE has one home; see the module
    # docstring for the drift guard that keeps this honest.
    _ingested_available,
    _visium_available,
    _xenium_available,
)
from .jobs import jsonable
from .projects import ProjectError

__all__ = ["build_workspace_router", "workspace_payload", "sections_available",
           "MAX_RECENT_RUNS", "WORKSPACE_CONTRACT_VERSION"]

#: The contract's cap on the recent-runs strip (``docs/WORKSPACE_CONTRACT.md``
#: §1: "across ALL projects, newest first, <= 12").
MAX_RECENT_RUNS = 12

#: Matches ``GET /api/account``'s ``contract`` field. One version for one API.
WORKSPACE_CONTRACT_VERSION = "v1"


# --------------------------------------------------------------------------- #
# Sections — the cheap half
# --------------------------------------------------------------------------- #
def _available(rec: Dict) -> Tuple[bool, Optional[str]]:
    """``describe_section``'s availability dispatch, and only that.

    The branch order is copied EXACTLY from ``cache.describe_section``:
    ingested first (an ingested record can carry any ``source``), then xenium,
    then visium as the fallthrough. Reordering these is not a style choice —
    an ingested Xenium directory would be routed to the built-in Xenium path,
    which looks for a ``cache_dir`` layout it does not have.
    """
    if rec.get("ingested"):
        return _ingested_available(rec)
    if rec.get("source") == "xenium":
        return _xenium_available(rec)
    return _visium_available(rec)


def sections_available() -> int:
    """How many registered sections this server can actually run.

    ``available`` is the same predicate ``GET /api/sections`` publishes: the
    cached bytes are present, and (for an ingested directory) the platform is
    one ``reconstruct_real`` accepts. A section that is registered but whose
    data is not on this disk is NOT counted — the home page's number is what
    you can run, not what the registry has heard of.

    Never raises: a section whose availability check itself blows up is counted
    as unavailable, because it is.
    """
    n = 0
    for rec in all_sections().values():
        try:
            ok, _reason = _available(rec)
        except Exception:  # noqa: BLE001 - one broken record must not 500 the home page
            continue
        if ok:
            n += 1
    return n


def drugs_available() -> int:
    """How many compounds a ``RunSpec``'s ``drug`` field may name.

    ``_known_drugs()`` is the same list ``POST /api/validate`` refuses against,
    so this count and the picker cannot disagree. Imported lazily because
    ``api.py`` imports this module, not the other way round.
    """
    from .api import _known_drugs

    try:
        return len(_known_drugs())
    except Exception:  # noqa: BLE001
        return 0


# --------------------------------------------------------------------------- #
# Run cards
# --------------------------------------------------------------------------- #
def _specimen_label(specimen: Dict) -> str:
    """The specimen's label, or a refusal that says what to do about it.

    A specimen is given a label at ``add_specimen`` time and the store will not
    let you blank one. A label can still fail to resolve: a project document
    edited by hand, restored from an older layout, or written by a version that
    did not have the field. That must render as a sentence, never as ``None``
    formatted into a card and never as a ``TypeError`` that takes the whole
    home page down with it.
    """
    label = specimen.get("label")
    if isinstance(label, str) and label.strip():
        return label.strip()
    sid = specimen.get("id")
    section = specimen.get("section")
    where = f"specimen {sid!r}" if sid else "an unidentified specimen"
    if isinstance(section, str) and section.strip():
        return (f"unnamed ({where}, section {section}) — open the project and "
                f"rename the specimen to label it")
    return (f"unnamed ({where}, no section recorded) — open the project; this "
            f"specimen record is incomplete and cannot be run")


def _run_card(job, project_id: str, project_name, specimen: Dict) -> Dict:
    """One flattened row of the cross-project recent-runs strip.

    Deliberately not a ``RunResult`` (contract §1). Every field is read
    straight off the :class:`~app.server.jobs.Job` — no ``result_view()``, no
    report body, no disk.
    """
    spec = job.spec if isinstance(job.spec, dict) else {}
    return {
        "run_id": job.run_id,
        "project_id": project_id,
        "project_name": project_name,
        "specimen_id": specimen.get("id"),
        "specimen_label": _specimen_label(specimen),
        "section": spec.get("section"),
        "drug": spec.get("drug"),
        "status": job.status,
        "created_at": job.created_at,
        "finished_at": job.finished_at,
        # The same expression ``Job.result_view`` uses, so the strip and the
        # run detail page cannot report two different runtimes. A job that
        # never started has no runtime, and 0.0 would be a claim that it ran
        # instantly.
        "runtime_s": round(job.elapsed_s, 3) if job.started_at else None,
    }


def _card_sort_key(card: Dict) -> float:
    """Newest first. A finished run sorts by when it finished; one still
    queued or running has no finish time and sorts by when it was queued."""
    for key in ("finished_at", "created_at"):
        v = card.get(key)
        if isinstance(v, (int, float)):
            return float(v)
    return float("-inf")


# --------------------------------------------------------------------------- #
# Thumbnail — NOT built here
# --------------------------------------------------------------------------- #
#
# ``ProjectStore`` owns contract §4 and publishes ``thumbnail`` on every row of
# ``store.list()``. This module does not build a second one: two implementations
# of "what may be drawn on a card" is exactly how a rule that says "never
# fabricate a picture of tissue" gets weakened by accident, in the copy nobody
# remembered to update.
#
# What IS supplied here is the resolver. ``store.list(resolve_run=...)`` maps a
# run id to its run view, and the store's docstring is explicit that "passed in
# (the API has a live JobManager) the manager is the truth and nothing is cached
# here". Omitting it makes the store read ``app/runs/<id>/result.json`` from the
# process-wide runs directory — which is both a disk read this endpoint does not
# need and, in a test or a second store, the WRONG runs.
def _run_resolver(jobs) -> Callable[[str], Optional[Dict]]:
    """A run id -> the four fields ``projects.run_composition`` reads.

    Deliberately NOT ``Job.result_view()``: that copies the whole run body,
    backfills the delivery law and recomputes the claim-set verdict on every
    call, and the thumbnail scan can touch every run of a project. These four
    come straight off the ``Job`` object with no disk and no recomputation, and
    they are the complete set ``run_composition`` and ``_thumbnail``'s refusal
    path look at: ``status``, ``report``, ``spec.section`` and ``run_id``.
    """
    def resolve(run_id: str) -> Optional[Dict]:
        job = jobs.get(run_id) if isinstance(run_id, str) else None
        if job is None:
            return None
        body = job.body if isinstance(job.body, dict) else {}
        return {
            "run_id": job.run_id,
            "status": job.status,
            "spec": job.spec if isinstance(job.spec, dict) else {},
            "report": body.get("report"),
        }
    return resolve


def _list_projects(store, resolve: Callable[[str], Optional[Dict]]) -> List[Dict]:
    """``store.list()``, handing it the live job manager when it can take one.

    ``resolve_run`` is being added to ``ProjectStore.list`` by the agent that
    owns that file. Checked rather than assumed, and rather than wrapped in a
    bare ``except TypeError`` — which would also swallow a TypeError raised
    from INSIDE the store and silently retry it down the slower path.
    """
    try:
        accepts = "resolve_run" in inspect.signature(store.list).parameters
    except (TypeError, ValueError):  # a store whose list() cannot be inspected
        accepts = False
    return store.list(resolve_run=resolve) if accepts else store.list()


# --------------------------------------------------------------------------- #
# The payload
# --------------------------------------------------------------------------- #
def _project_sort_key(row: Dict) -> float:
    for key in ("last_opened_at", "updated_at", "created_at"):
        v = row.get(key)
        if isinstance(v, (int, float)):
            return float(v)
    return float("-inf")


def workspace_payload(jobs, store) -> Dict:
    """Build the contract's §1 object. Pure: no request, no globals, no cache.

    ``jobs`` is a :class:`~app.server.jobs.JobManager`, ``store`` a
    :class:`~app.server.projects.ProjectStore`. Passed in rather than imported
    so the tests can build a workspace out of an isolated pair — the same
    reason ``build_account_router`` takes its limits as an argument.
    """
    summaries = _list_projects(store, _run_resolver(jobs))

    cards: List[Dict] = []
    compounds_by_project: Dict[str, List[str]] = {}
    n_specimens_total = 0
    # A SET, not a counter, and that is the whole point. See the `totals`
    # block at the bottom of this function for why runs are counted this way
    # and specimens deliberately are not.
    run_ids_total: set = set()

    for row in summaries:
        pid = row.get("id")
        name = row.get("name")

        # `store.list()` returns this shape for a project whose document is on
        # disk but unreadable. It has no specimens to join and no runs to show;
        # it is still a project and still counted, and the `unreadable` key
        # travels untouched so the card can say so. It gets NO `compounds`
        # key: absent means "this server could not tell you", and that is
        # exactly the truth about a document it could not parse.
        if row.get("unreadable"):
            continue

        # Specimens ARE summed. A specimen is not a shared object: it lives
        # inside the project document, `ProjectStore` has no specimen store,
        # and `add_specimen` allocates its id as max+1 of THAT project's ids --
        # so `s_00001` exists in nearly every project and two projects'
        # `s_00001` are two different patients. Contract §1 draws the line in
        # the duplicate row itself: specimens are COPIED, runs are not.
        # Deduplicating these by id would report one specimen for a workspace
        # holding two.
        n_specimens_total += int(row.get("n_specimens") or 0)

        try:
            # refresh_capability=False: see the module docstring. A capability
            # report re-reads the section's panel gene list per specimen and
            # nothing on this page renders one.
            doc = store.get(pid, refresh_capability=False)
        except ProjectError:
            # ORDINARY. The workspace is a snapshot of a store that two other
            # routes are mutating: `list()` and `get()` are separate reads
            # under separate acquisitions of the store lock, so a DELETE
            # landing between them is expected traffic, not a fault. Nothing
            # is logged, because a log line per deleted project is a log
            # nobody reads.
            continue
        except Exception as exc:  # noqa: BLE001 - one bad document is not a 500
            # NOT ordinary, and the difference from the branch above is the
            # whole reason there are two. A fault here silently removes this
            # project's runs from the strip; swallowing it without a record
            # means the page is wrong and nothing anywhere says so. Logged
            # server-side (never to the browser: `log_internal` is the same
            # path `redact_exception` writes, so a fault has exactly one
            # stderr record in exactly one format), and the page still
            # renders.
            log_internal(exc, f"workspace: specimens of project {pid!r} "
                              f"could not be read")
            continue

        compounds: set = set()
        compounds_complete = True

        for specimen in (doc.get("specimens") or []):
            if not isinstance(specimen, dict):
                compounds_complete = False
                continue
            for run_id in (specimen.get("run_ids") or []):
                # Counted HERE, off the document, and before the job lookup:
                # a dangling id is still a run this workspace references, and
                # `totals.runs` has always reported it (see
                # tests/test_workspace.py::test_a_dangling_run_id_is_skipped_
                # not_fatal). Only `str` is counted -- a non-string entry in a
                # hand-edited document is not a run id, cannot resolve to a
                # job, and may not even be hashable, and a `set.add` TypeError
                # here would take the whole home page down.
                if isinstance(run_id, str):
                    run_ids_total.add(run_id)
                job = jobs.get(run_id) if isinstance(run_id, str) else None
                if job is None:
                    # A dangling id: the association survives in the project
                    # document, the run does not survive in app/runs/. It has
                    # no timestamp, so it cannot be placed in a strip ordered
                    # by time at all. Skipped here; still counted by
                    # `totals.runs`, which reports what the cards report.
                    compounds_complete = False
                    continue
                cards.append(_run_card(job, pid, name, specimen))
                drug = (job.spec or {}).get("drug") if isinstance(
                    job.spec, dict) else None
                if isinstance(drug, str) and drug.strip():
                    compounds.add(drug.strip())
                else:
                    compounds_complete = False

        # A PARTIAL set is not the set. If even one attached run could not be
        # read, this server does not know which compounds this project has run
        # and says so by omitting the key -- the card renders "compounds not
        # reported by this server", which is true, rather than a list that
        # looks complete and is not.
        if compounds_complete:
            compounds_by_project[pid] = sorted(
                compounds, key=lambda c: (c.casefold(), c))

    cards.sort(key=_card_sort_key, reverse=True)

    projects = []
    for row in summaries:
        out = dict(row)
        # Contract §3 and §4, both owned by app/server/projects.py. Defaulted
        # rather than computed: if the store publishes them its values travel
        # untouched, and a row that lacks them (the `unreadable` row does)
        # gets null -- which renders as "not opened yet", NEVER as an old date
        # and never as now.
        out.setdefault("last_opened_at", None)
        out.setdefault("thumbnail", None)
        pid = row.get("id")
        # THREE states, and the middle one is the reason this is a setdefault
        # and not an assignment:
        #   key absent  -> this server cannot say (unreadable project, or an
        #                  attached run it could not resolve)
        #   []          -> computed, and no compound has been run yet
        #   [names...]  -> a SET. Deduplicated, alphabetical, and ordered by
        #                  NOTHING else. Cross-compound ordering by predicted
        #                  effect was withdrawn 2026-08-21 and the prohibition
        #                  is absolute: no metric, no "top", no badge, and no
        #                  sort that could be mistaken for one.
        if pid in compounds_by_project:
            out["compounds"] = compounds_by_project[pid]
        projects.append(out)

    # Contract §3: last_opened_at, falling back to updated_at when it is null.
    # With no project ever opened this reduces exactly to §1's "newest-updated
    # first".
    projects.sort(key=_project_sort_key, reverse=True)

    return {
        "projects": projects,
        "recent_runs": cards[:MAX_RECENT_RUNS],
        "sections_available": sections_available(),
        "drugs_available": drugs_available(),
        # Workspace quantities, not server-wide ones. The server's total run
        # count (including runs attached to no project) is GET /api/health's
        # `runs` block.
        #
        # `runs` IS NOT THE SUM OF THE CARDS, and it must not be. Contract §2:
        # a project "owns nothing but ids", so one run in app/runs/ may be
        # referenced by any number of projects -- contract §1 makes that
        # ordinary on purpose, since `POST /api/projects/{id}/duplicate` copies
        # run ids and re-executes nothing. Summing each row's `n_runs` counted
        # a shared run once per project holding it: MEASURED at 4 on a server
        # holding 2 after a single duplicate, and it grows with every copy.
        # So this counts DISTINCT ids.
        #
        # It costs nothing extra. The ids come out of the `store.get()` this
        # function already performs per project for the run cards and the
        # compound set; no read, no stat and no lock acquisition was added.
        # What it does NOT include is a project whose document could not be
        # read between `list()` and `get()` -- a mid-request DELETE, or the
        # logged fault above. That project has no readable ids, so it
        # contributes none, exactly as it already contributes no run cards and
        # no compound set. Its row and its own `n_runs` still render.
        #
        # `specimens` stays a SUM, for the reason stated where it is
        # accumulated: specimen ids are per-project and collide by
        # construction, so deduplicating them would UNDERCOUNT.
        "totals": {
            "projects": len(summaries),
            "runs": len(run_ids_total),
            "specimens": n_specimens_total,
        },
        # `stub` is a CLIENT-side state: app/web/js/api.js sets `api.stub` when
        # the server cannot be reached and serves js/mock.js instead. A payload
        # that came FROM the server is by construction not stub data, so this
        # is false here and true only in the mock. The field exists so the home
        # page reads one key regardless of which of the two answered it.
        "server": {"contract": WORKSPACE_CONTRACT_VERSION, "stub": False},
    }


# --------------------------------------------------------------------------- #
# Warming
# --------------------------------------------------------------------------- #
def _warm() -> None:
    """Pay this endpoint's one-time IMPORT costs off the request path.

    MEASURED. Neither of these is per-request work; both are the first-touch
    cost of a module, and the first request to touch it pays all of it:

        all_sections()   first call 337.6 ms, every call after it 0.00 ms
                         -- `app.server.ingest_registry` -> `tmesim.data.ingest`
        _known_drugs()   first call 169.7 ms, every call after it 0.03 ms
                         -- `tmesim.molecular.library`

    Together they are 592.6 ms of the 592.6 ms cold / 48.4 ms warm spread this
    endpoint measures on a dev box, and roughly 15x that on the 0.1-CPU
    instance this deploys to. Moving them to a daemon thread at app-build time
    is the same trick `create_app`'s lifespan already plays with `drug_rows()`,
    for the same reason.

    It is a HEAD START, NOT A GUARANTEE: a request that arrives before the
    thread finishes simply pays the import itself, exactly as it does today. On
    a real deploy the browser must fetch `/`, parse it and load its modules
    before it can call this route, so the head start is seconds long.

    What it deliberately does NOT warm is `list_sections()`. That would write
    `app/runs/_section_meta.json` and make a later `GET /api/sections` fast --
    at the cost of spending a 0.1-CPU instance's entire CPU on a file read at
    the precise moment the home page is trying to paint. The Studio's cold
    start is not this page's problem to solve by stealing its CPU.
    """
    try:
        all_sections()
    except Exception:  # noqa: BLE001 - warming is best-effort, like _warm_drugs
        pass
    try:
        drugs_available()
    except Exception:  # noqa: BLE001
        pass


# --------------------------------------------------------------------------- #
# Router
# --------------------------------------------------------------------------- #
async def _optional_body(request: Request) -> Dict:
    """The JSON object in the body, or ``{}`` when there is no body at all.

    Both routes below take a body that is entirely optional — ``/open`` takes
    no parameters and ``/duplicate`` takes at most a name — so an ABSENT body
    is a well-formed request, not a malformed one. A body that is present and
    is not a JSON object still fails exactly as it does on every other route:
    through ``_json_body``, so the 400 a client sees is the same sentence in
    the same shape whichever route produced it.
    """
    from .api import _json_body

    raw = await request.body()
    if not raw.strip():
        return {}
    return await _json_body(request)


def build_workspace_router(*, jobs, store) -> APIRouter:
    """Mount ``GET /api/workspace`` and the home page's two project mutations.

    ``POST /api/projects/{id}/open`` and ``POST /api/projects/{id}/duplicate``
    live here rather than beside the other eleven project routes in ``api.py``
    for one reason: they exist because the workspace exists. Nothing in the
    Studio opens or duplicates a project. The mounting block in ``api.py``
    names them so the next reader looking for a project route finds them.

    Every route here is gated by the same middleware as every other ``/api/``
    route: the gate in ``create_app`` tests
    ``get_route_path(request.scope).startswith("/api/")``, so a router mounted
    here is authenticated by construction and cannot be forgotten. Asserted
    rather than assumed — see
    ``tests/test_workspace.py::test_every_workspace_route_is_gated``.
    """
    router = APIRouter()
    threading.Thread(target=_warm, name="tmesim-warm-workspace",
                     daemon=True).start()

    # The SAME translation every other project route uses, so a refusal the
    # store authored reaches the client with its status and its wording
    # unchanged. Bound once at app-build time: ``api.py`` imports this module
    # from inside ``create_app``, so by the time this runs ``api`` is fully
    # loaded and there is no cycle and no per-request import.
    from .api import _project_errors

    @router.get("/api/workspace")
    def get_workspace() -> Dict:
        """Everything the workspace home page needs, in one round trip.

        See ``docs/WORKSPACE_CONTRACT.md`` §1 for the frozen shape and this
        module's docstring for what is deliberately NOT read to build it.
        """
        return jsonable(workspace_payload(jobs, store))

    @router.post("/api/projects/{project_id}/open")
    async def post_project_open(project_id: str, request: Request) -> Dict:
        """Stamp ``last_opened_at``: a human looked at this project.

        Idempotent — calling it twice is not an error, it is a later
        timestamp. It deliberately does NOT move ``updated_at``; that
        distinction is contract §3 and the entire reason the field exists, and
        the store returns both so a caller can see that it did not move.

        A missing project is a 404, like every other project route, and it
        carries the store's own sentence. The Studio fires this on boot and
        swallows the rejection on purpose: recency is an ANNOTATION, so
        failing to record that a project was opened must never stop it being
        opened. What that swallow must not become is this route lying — a 200
        for an id that does not exist would make "recorded" and "there is
        nothing to record" the same answer to every client, including the ones
        that do not swallow.
        """
        await _optional_body(request)   # reject a malformed body, ignore none
        with _project_errors():
            return jsonable(store.mark_opened(project_id))

    @router.post("/api/projects/{project_id}/duplicate", status_code=201)
    async def post_project_duplicate(project_id: str,
                                     request: Request) -> Dict:
        """Copy a project's specimens and their run REFERENCES. Body ``{name?}``.

        **No run is re-executed and no run is copied.** ``app/runs/`` is the
        single store, a specimen holds ids into it, and after this call two
        specimens in two projects legitimately point at the same run id. That
        is what makes the operation instant, and it is also why deleting
        either project still deletes no run (contract §2).
        """
        body = await _optional_body(request)
        with _project_errors():
            return jsonable(store.duplicate(project_id, body.get("name")))

    return router
