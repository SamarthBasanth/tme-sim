"""A stored run must never show a withdrawn claim as if it were current.

The defect
----------
``app/runs/<id>/result.json`` stores each run's caveats AS THEY WERE WHEN IT
RAN, and that is correct: it is the execution record, and editing it to agree
with a later claim would falsify the one artefact this repo is most careful
about. But the Studio rendered them verbatim, so opening ``r_00052`` today
showed:

    "The CROSS-DRUG RANKING is invariant across all three routes and is the
     defensible output"

-- a claim withdrawn on 2026-08-21 -- with nothing saying it had been. That
sentence is in 65 stored run records.

Why this is not a date cutoff
-----------------------------
"Banner anything created before 2026-08-21" works exactly once. This was the
SECOND claim retirement in a week (the grid-convergence headline was the
first), so the mechanism has to survive the third without a code change.

So the record carries a VERSION. ``app/server/jobs.py`` stamps
:data:`CLAIM_SET_VERSION` into ``provenance.claim_set`` when a run executes;
:func:`claim_set_status` compares that stamp against the current version on
every read and banners anything older or absent. The next withdrawal is a
bumped integer and an appended :data:`CLAIM_SET_LOG` entry, and every record
written before it -- including ones written between now and then -- banners
itself.

What these tests hold down
--------------------------
1. the version record is well-formed and cannot be bumped without saying what
   was withdrawn (:class:`TestTheClaimSetRecord`);
2. the verdict is keyed to the VERSION and never to a date
   (:class:`TestTheVerdictIsKeyedToTheVersion`);
3. a stored record is never rewritten and its caveats are never hidden
   (:class:`TestTheRecordIsNeverRewritten`);
4. a run executed now is stamped and is NOT bannered
   (:class:`TestARunExecutedNowIsNotBannered`);
5. the Studio renders the banner above the caveats and still renders every one
   of them (:class:`TestTheStudioMarksAndDoesNotHide`);
6. this warning and ``report.py``'s "RESTORED FROM A WITHDRAWN PANEL DOCUMENT"
   have not drifted into two vocabularies for one situation
   (:class:`TestOneVocabulary`).
"""
from __future__ import annotations

import json
import os
import time

import pytest

from app.server.jobs import (
    CLAIM_SET_LOG,
    CLAIM_SET_VERSION,
    JobManager,
    claim_set_stamp,
    claim_set_status,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(REPO, "app", "web")

#: The exact sentence that made this necessary, as ``r_00052`` stores it.
RETIRED = ("The CROSS-DRUG RANKING is invariant across all three routes and is "
           "the defensible output")


def _read(*parts: str) -> str:
    with open(os.path.join(REPO, *parts), encoding="utf-8") as fh:
        return fh.read()


def _python_sources():
    """Every Python source that can put text on the product surface.

    ``app/`` and ``tmesim/`` only: ``tests/``, ``experiments/``, ``archive/``
    and ``scripts/`` do not render to a user, and ``archive/`` is deliberately
    full of superseded claims.
    """
    for top in ("app", "tmesim"):
        for root, dirs, files in os.walk(os.path.join(REPO, top)):
            dirs[:] = [d for d in dirs
                       if d not in ("__pycache__", "runs", "reports",
                                    "uploads", "projects")]
            for name in sorted(files):
                if name.endswith(".py"):
                    yield os.path.join(root, name)


def _live_string_literals(path: str):
    """(lineno, value) for every string literal that is NOT prose.

    Excluded: module/class/function docstrings, and the ``CLAIM_SET_LOG``
    assignment itself -- which necessarily contains the retired sentences it
    is the record of. Everything else is a string the product can emit.

    Adjacent literals are already concatenated by the parser, which is the
    whole reason this reads an AST: a caveat is written as three quoted
    fragments on three lines and no line contains the sentence.
    """
    import ast

    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=path)

    exempt = set()

    def _mark_docstring(node):
        body = getattr(node, "body", None)
        if not body:
            return
        first = body[0]
        if (isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)):
            exempt.add(id(first.value))

    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            _mark_docstring(node)
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "CLAIM_SET_LOG"
                for t in node.targets):
            for sub in ast.walk(node.value):
                exempt.add(id(sub))

    for node in ast.walk(tree):
        if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                and id(node) not in exempt):
            yield node.lineno, node.value


def _legacy_record(run_id: str = "r_09999") -> dict:
    """A run document exactly as one written before the stamp existed.

    Not a paraphrase: ``provenance`` has no ``claim_set`` key at all, and the
    caveat is the retired sentence in the words ``r_00052`` uses.
    """
    return {
        "run_id": run_id,
        "spec": {"section": "xenium_319", "dx_um": 12.0,
                 "drug": "temozolomide", "total_time_s": 600.0},
        "status": "done",
        "stage": "done",
        "progress": 1.0,
        "created_at": 1785293834.777318,
        "started_at": 1785293834.9,
        "finished_at": 1785293900.0,
        "report": {"n_cells": 40887},
        "metrics": {"kill_fraction": 0.5},
        "provenance": {
            "alpha": 148.60815765445545,
            "delivery_law": {"caveat": "recorded", "recorded": True},
            "caveats": [
                "ABSOLUTE KILL MUST NOT BE QUOTED FROM THIS SECTION UNDER ANY "
                "STATE ROUTE. " + RETIRED + "; kill_fraction and any kill "
                "ceiling on this section carry an unresolved uncertainty of "
                "that size.",
                "This tissue is an untreated resection with no outcome label. "
                "NOT FOR CLINICAL USE.",
            ],
        },
        "log_tail": ["00:00:00 queued"],
    }


def _store(root: str, doc: dict) -> str:
    d = os.path.join(root, doc["run_id"])
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "result.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
    return path


# --------------------------------------------------------------------------- #
# 1. The version record
# --------------------------------------------------------------------------- #
class TestTheClaimSetRecord:
    def test_the_log_is_monotone_and_starts_at_one(self):
        versions = [e["version"] for e in CLAIM_SET_LOG]
        assert versions == list(range(1, len(versions) + 1)), (
            "CLAIM_SET_LOG must be one entry per version, oldest first, "
            f"numbered from 1; got {versions}")

    def test_a_bump_without_a_note_fails_here(self):
        """The whole promise of the scheme: a version bump is not enough.

        If someone raises ``CLAIM_SET_VERSION`` to retire a claim and does not
        say WHAT was retired, every stale record gets a banner that can only
        say "something changed". This is the assertion that stops that.
        """
        assert CLAIM_SET_LOG[-1]["version"] == CLAIM_SET_VERSION, (
            f"CLAIM_SET_VERSION is {CLAIM_SET_VERSION} but the newest "
            f"CLAIM_SET_LOG entry is {CLAIM_SET_LOG[-1]['version']}. Bumping "
            "the version means appending the entry that says what it withdrew.")

    @pytest.mark.parametrize("entry", CLAIM_SET_LOG,
                             ids=[str(e["version"]) for e in CLAIM_SET_LOG])
    def test_every_entry_says_what_changed_and_where_it_was_measured(self, entry):
        assert len(entry["withdrew"]) > 120, (
            "`withdrew` is rendered to a reader who is about to be shown the "
            "withdrawn sentence. One clause is not enough.")
        assert entry["sources"], "a withdrawal with no source is an assertion"
        for src in entry["sources"]:
            assert os.path.exists(os.path.join(REPO, src)), (
                f"claim set {entry['version']} cites {src}, which does not "
                "exist. A citation that cannot be opened is not a citation.")
        # ISO date, so the note can say WHEN without anything downstream
        # parsing prose for one.
        time.strptime(entry["adopted"], "%Y-%m-%d")

    def test_the_stamp_carries_nothing_that_can_go_stale(self):
        """A stamp says what the run was written against. Nothing else.

        Anything in the stamp that describes the CURRENT claim set would be a
        judgement frozen at write time -- true then, false after the next
        withdrawal -- which is the exact defect this feature exists to fix.
        """
        assert set(claim_set_stamp()) == {"version", "adopted"}
        assert claim_set_stamp()["version"] == CLAIM_SET_VERSION


# --------------------------------------------------------------------------- #
# 2. The verdict
# --------------------------------------------------------------------------- #
class TestTheVerdictIsKeyedToTheVersion:
    def test_an_unstamped_record_is_superseded_and_says_it_is_unrecorded(self):
        st = claim_set_status(None)
        assert st["superseded"] is True
        assert st["recorded"] is False
        assert st["version"] is None
        # An unstamped record predates every recorded claim set, so it is owed
        # every note -- not "version 0 of a series", which it never was.
        assert [c["version"] for c in st["changes"]] == \
            sorted((e["version"] for e in CLAIM_SET_LOG), reverse=True)
        assert st["warning"] and st["headline"] and st["detail"]
        assert "not recoverable" in st["detail"]

    def test_a_record_stamped_current_is_not_bannered(self):
        st = claim_set_status(claim_set_stamp(), ["an ordinary caveat"])
        assert st["superseded"] is False
        assert st["contradicted"] is False
        assert st["recorded"] is True
        assert st["changes"] == []
        assert st["warning"] is None
        assert st["headline"] is None

    def test_an_older_stamp_is_told_only_what_changed_after_it(self):
        st = claim_set_status({"version": 0, "adopted": "2026-01-01"})
        assert st["superseded"] is True
        assert st["recorded"] is True
        assert st["version"] == 0
        assert all(c["version"] > 0 for c in st["changes"])
        assert len(st["changes"]) == len(CLAIM_SET_LOG)
        # It knows which set it was written against, so it says so rather than
        # claiming the era is unrecoverable.
        assert "claim set 0" in st["detail"]

    def test_the_verdict_does_not_put_the_retired_sentence_on_the_wire(self):
        """The one thing a "this claim is withdrawn" payload must not carry.

        ``retired_text`` exists for the source scan below. Shipping it inside
        the verdict would put the withdrawn sentence in the API response that
        exists to say it is withdrawn, one careless render away from being back
        on a screen under a banner disowning it.
        """
        st = claim_set_status(None)
        assert st["changes"]
        for change in st["changes"]:
            assert set(change) == {"version", "adopted", "withdrew", "sources"}
        blob = json.dumps(st)
        for entry in CLAIM_SET_LOG:
            for sentence in entry.get("retired_text") or ():
                assert sentence not in blob

    def test_the_newest_entry_leads_the_change_list(self):
        st = claim_set_status({"version": 0})
        got = [c["version"] for c in st["changes"]]
        assert got == sorted(got, reverse=True), (
            "the most recent withdrawal is the one a stale reader most needs "
            "first")

    def test_a_stamp_from_the_future_is_not_bannered(self):
        """A record newer than this build is not stale; it is unfamiliar.

        Bannering it would tell a reader that a claim they are looking at has
        been withdrawn when in fact this build has not caught up yet.
        """
        st = claim_set_status({"version": CLAIM_SET_VERSION + 5})
        assert st["superseded"] is False
        assert st["changes"] == []

    @pytest.mark.parametrize("junk", [{}, {"version": None}, {"version": "1"},
                                      {"version": 1.0}, "1", 1, []])
    def test_an_unreadable_stamp_fails_towards_the_banner(self, junk):
        """Fail safe: if the stamp cannot be read, say so rather than assume.

        A silently-accepted stamp is a record that stops bannering itself for
        a reason no one can see.
        """
        st = claim_set_status(junk)
        assert st["superseded"] is True
        assert st["recorded"] is False

    def test_the_verdict_never_reads_a_timestamp(self):
        """The point of the scheme, asserted directly.

        Two records with the SAME creation time get opposite verdicts, decided
        only by the stamp; and a record created long before the current claim
        set was adopted is NOT bannered if it carries the current stamp. Any
        reintroduction of a `created_at <= cutoff` rule fails one of these.
        """
        stamped_old_file = claim_set_status(claim_set_stamp(), ["fine"])
        unstamped_new_file = claim_set_status(None, ["fine"])
        assert stamped_old_file["superseded"] is False
        assert unstamped_new_file["superseded"] is True
        # And the function is not given a time it COULD consult. It takes the
        # stamp and the record's own caveats -- both properties of the record,
        # neither a clock. A `created_at`/`before`/`cutoff` parameter would be
        # the date rule coming back in through the signature.
        import inspect as _inspect
        params = set(_inspect.signature(claim_set_status).parameters)
        assert params == {"stamp", "caveats"}, params
        assert not any(w in p for p in params
                       for w in ("date", "time", "cutoff", "created", "before"))


# --------------------------------------------------------------------------- #
# 2b. The finding an age comparison cannot make
# --------------------------------------------------------------------------- #
class TestARecordCanContradictItsOwnStamp:
    """A stamp says what the build DECLARED, not what the build DID.

    ``execute_run`` kept emitting "The CROSS-DRUG RANKING is invariant across
    all three routes and is the defensible output" for a day after that claim
    was withdrawn. Any run executed in that window, once stamping existed,
    would be stamped claim set 1 AND carry a sentence claim set 1 retired --
    and the era check would call it CLEAN, because the version it compares is
    exactly the thing that is lying.

    So the verdict makes two findings. SUPERSEDED is about age. CONTRADICTED
    is about compliance: the record's own caveats disagree with the claim set
    it names. Neither edits the record; both banner.
    """

    def test_a_current_stamp_does_not_excuse_a_retired_sentence(self):
        cav = ["ABSOLUTE KILL MUST NOT BE QUOTED. " + RETIRED + "; and so on."]
        st = claim_set_status(claim_set_stamp(), cav)
        assert st["superseded"] is False, "it is not old"
        assert st["contradicted"] is True, "it is wrong"
        assert st["contradicts"] and st["contradicts"][0]["version"] == 1
        assert "CONTRADICTS THE CLAIM SET IT IS STAMPED WITH" in st["headline"]

    def test_the_contradiction_leads_when_both_findings_hold(self):
        """"This record disagrees with the claim set it names" is a stranger
        fact than "this record is old", and a reader told only the second will
        read the stamp as a reason to trust the text."""
        st = claim_set_status({"version": 1}, [RETIRED])
        # Force both: bump the notional current set by checking a stamp that is
        # simultaneously behind and non-compliant is impossible at n=1, so this
        # asserts the ordering rule directly on the branch that decides it.
        assert st["contradicted"] is True
        assert st["headline"].startswith("THIS RECORD CONTRADICTS")

    def test_an_unstamped_record_is_old_rather_than_non_compliant(self):
        """It never claimed to satisfy anything, so it cannot have failed to.

        Reporting the same record as both would say it is broken when it is
        merely from before the rule existed.
        """
        st = claim_set_status(None, [RETIRED])
        assert st["superseded"] is True
        assert st["contradicted"] is False
        assert st["contradicts"] == []
        assert "SUPERSEDED CLAIM SET" in st["headline"]

    def test_a_sentence_retired_after_the_stamp_is_age_not_contradiction(self):
        """Only claim sets the record CLAIMS to satisfy are checked.

        A record stamped 0 carrying a sentence claim set 1 later retired did
        nothing wrong; it is simply old. Getting this backwards would flag
        every historical record as non-compliant.
        """
        st = claim_set_status({"version": 0}, [RETIRED])
        assert st["superseded"] is True
        assert st["contradicted"] is False

    def test_a_clean_current_record_says_nothing(self):
        st = claim_set_status(claim_set_stamp(),
                              ["WHAT SURVIVES IS THE MAP: within this one "
                               "compound on this one section."])
        assert st["superseded"] is False and st["contradicted"] is False
        assert st["warning"] is None

    def test_r_00090_is_the_live_specimen_and_is_marked(self):
        """``app/runs/r_00090`` IS RETAINED DELIBERATELY. DO NOT DELETE IT.

        If you are tidying ``app/runs/``: this record looks wrong because it
        IS wrong, and that is why it is kept. It must survive any future
        cleanup of that directory.

        It was submitted through the Studio during this change, minutes before
        the emitter in ``execute_run`` was corrected -- so it is a genuine
        artefact of a build that stamped a claim set it was not honouring, and
        it is the only one in existence. It is what proved that the era check
        alone returns CLEAN for exactly the case this feature was built to
        catch. ``r_00091`` is its clean twin: same spec, corrected build.

        Deleting it would be the same move this whole feature refuses --
        editing the record so the guard looks like it was always there. It is
        skipped rather than failed where the store is absent, because
        ``app/runs/`` is gitignored and a fresh checkout has no store at all.
        """
        path = os.path.join(REPO, "app", "runs", "r_00090", "result.json")
        if not os.path.exists(path):
            pytest.skip("app/runs/r_00090 is not present (gitignored store)")
        mgr = JobManager(root=os.path.join(REPO, "app", "runs"),
                         start_worker=False)
        view = mgr.get("r_00090").result_view()
        prov = view["provenance"]
        cs = prov["claim_set"]

        assert any(RETIRED in c for c in prov["caveats"]), (
            "r_00090 is only a specimen while it still carries the retired "
            "sentence; if it no longer does, the record was edited")
        assert cs["version"] == CLAIM_SET_VERSION
        assert cs["superseded"] is False, (
            "and this is the point: the era check calls it clean")
        assert cs["contradicted"] is True, (
            "the compliance check is what catches it")
        assert cs["warning"] and "CONTRADICTS" in cs["warning"]

        # Its clean twin, same spec, corrected build: no finding at all.
        twin = os.path.join(REPO, "app", "runs", "r_00091", "result.json")
        if os.path.exists(twin):
            tcs = mgr.get("r_00091").result_view()["provenance"]["claim_set"]
            assert tcs["superseded"] is False and tcs["contradicted"] is False


# --------------------------------------------------------------------------- #
# 3. The stored record
# --------------------------------------------------------------------------- #
class TestTheRecordIsNeverRewritten:
    def test_reading_a_legacy_run_does_not_touch_its_file(self, tmp_path):
        root = str(tmp_path)
        path = _store(root, _legacy_record())
        before = open(path, "rb").read()

        mgr = JobManager(root=root, start_worker=False)
        job = mgr.get("r_09999")
        assert job is not None
        for _ in range(3):
            view = job.result_view()
        assert view["provenance"]["claim_set"]["superseded"] is True

        assert open(path, "rb").read() == before, (
            "reading a run must not rewrite it. app/runs/ is the execution "
            "record; the verdict is computed for the reader, not stored.")

    def test_the_banner_does_not_hide_the_stored_caveats(self, tmp_path):
        root = str(tmp_path)
        stored = _legacy_record()
        _store(root, stored)
        mgr = JobManager(root=root, start_worker=False)

        served = mgr.get("r_09999").result_view()["provenance"]["caveats"]
        assert served == stored["provenance"]["caveats"], (
            "every caveat, in order, word for word. A reader is entitled to "
            "see what the run said; marking it historical is honest and "
            "deleting it is not.")
        assert any(RETIRED in c for c in served), (
            "the withdrawn sentence is still on screen -- marked, not removed")

    def test_the_persisted_document_carries_the_stamp_and_no_verdict(self, tmp_path):
        """``_persist`` writes ``result_view(annotate=False)``.

        A stored ``"superseded": false`` would be true the moment it was
        written and false forever after the next withdrawal -- this feature's
        own defect, one layer down. So the disk gets the stamp only.
        """
        from app.server.jobs import Job

        mgr = JobManager(root=str(tmp_path), start_worker=False)
        job = Job("r_08888", {"section": "xenium_319", "drug": "temozolomide"})
        job.status = "done"
        job.started_at = time.time()
        job.finished_at = job.started_at + 1
        job.body = {"provenance": {"caveats": ["a"],
                                   "delivery_law": {"caveat": "x"},
                                   "claim_set": claim_set_stamp()}}
        mgr._persist(job, None, status="done")

        with open(os.path.join(str(tmp_path), "r_08888", "result.json")) as fh:
            doc = json.load(fh)
        assert doc["provenance"]["claim_set"] == {
            "version": CLAIM_SET_VERSION,
            "adopted": CLAIM_SET_LOG[-1]["adopted"],
        }
        assert "superseded" not in doc["provenance"]["claim_set"]

        # ...and the same job, read through the API view, does carry one.
        assert job.result_view()["provenance"]["claim_set"]["superseded"] is False

    def test_a_run_with_no_provenance_is_left_alone(self, tmp_path):
        """A failed or queued run has no caveats, so there is nothing to mark."""
        from app.server.jobs import Job

        JobManager(root=str(tmp_path), start_worker=False)
        job = Job("r_07777", {"section": "xenium_319"})
        job.status = "failed"
        job.error = "boom"
        view = job.result_view()
        assert view["provenance"] is None

    def test_the_real_store_banners_r_00052(self):
        """The run named in the bug report, against the real local store.

        ``app/runs/`` is gitignored, so this skips where the store is absent.
        It is here because the synthetic fixture above proves the mechanism and
        this proves it fires on the actual record.
        """
        path = os.path.join(REPO, "app", "runs", "r_00052", "result.json")
        if not os.path.exists(path):
            pytest.skip("app/runs/r_00052 is not present (gitignored store)")
        with open(path) as fh:
            doc = json.load(fh)
        prov = doc.get("provenance") or {}
        assert any(RETIRED in c for c in prov.get("caveats") or []), (
            "r_00052 no longer carries the retired sentence -- if the record "
            "was edited, that is a much worse bug than the one this fixes")
        st = claim_set_status(prov.get("claim_set"))
        assert st["superseded"] is True
        assert st["changes"], "and it is told what was withdrawn"


# --------------------------------------------------------------------------- #
# 4. A run executed now
# --------------------------------------------------------------------------- #
class TestARunExecutedNowIsNotBannered:
    def test_no_live_string_anywhere_repeats_a_retired_claim(self):
        """The stamp would otherwise certify a withdrawn claim as current.

        ``execute_run`` was, until this change, STILL appending "The CROSS-DRUG
        RANKING is invariant across all three routes and is the defensible
        output" to every run on a hypoxia-blind panel. ``api.py``'s submit
        warning, ``app/web/js/mock.js`` and ``report.py``'s limits were all
        corrected when the ranking was withdrawn; this emitter was missed, so
        for a day after the withdrawal every run wrote the withdrawn sentence
        into a FRESH record. That is worse than the five that were caught: a
        stale string in a warning misleads one reader, a stale string in an
        emitter manufactures new artefacts carrying it -- and this feature
        would then have stamped those artefacts current.

        Three deliberate choices:

        * **By AST, not by line.** These caveats are built from adjacent string
          literals, so the sentence spans three source lines and a line scan
          would never see it. The parser concatenates them; this reads what the
          parser produced.
        * **Comments and docstrings are exempt.** This repo corrects a stale
          claim by QUOTING it under a banner rather than deleting it, so the
          sentence legitimately survives in prose. A scan that cannot tell a
          live string from a quoted corpse either fires on every correction or
          forces people to stop quoting what they retract. (Same reasoning as
          ``tests/test_claim_guard.py``, reached a different way: that module
          declares claims opt-in, this one distinguishes code from prose.)
        * **Keyed to CLAIM_SET_LOG, not to this one sentence.** The next
          withdrawal gets this guard by adding ``retired_text`` to its log
          entry -- no new test.

        WHAT IT COVERS
            Every ``.py`` under ``app/`` and ``tmesim/`` -- every Python path
            that can put text on the product surface, including the engine
            modules that reach a user without going through this package
            (``tmesim/data/ingest.py`` rides into ``/api/sections[].caveats``,
            ``tmesim/protocol/resolve.py`` renders as an editor diagnostic).

        WHAT IT DOES NOT, stated plainly because a guard mistaken for complete
        is worse than no guard:

        * **f-strings.** A sentence broken by an interpolation is stored by the
          parser in pieces and is invisible here BY CONSTRUCTION. Closing it
          needs partial evaluation, not parsing. This is the hardest limit of
          the three and it is not attempted.
        * **JavaScript.** ``app/web/js/mock.js`` writes caveats, so a retired
          claim could survive there. Lower severity rather than equal: mock.js
          runs only in STUB MODE, when the backend did not answer, and every
          number it produces is already labelled ``synthetic`` and refused as a
          result. Dev-facing, not customer-facing -- but not zero.
        * **``results/*.md``.** The prose claim surface. Out of scope here;
          ``tests/test_claim_guard.py`` is the module that guards documents.

        The end-to-end assertion in
        ``test_app_server.py::test_result_carries_the_current_claim_set``
        catches any origin regardless of language or string construction, but
        only for the caveats that one run on one section happens to exercise.
        Neither is complete. They fail in DIFFERENT directions -- the scan
        misses constructions, the run misses branches -- which is the whole
        reason both are here, and why neither should be read as coverage.
        """
        retired = [(e["version"], t) for e in CLAIM_SET_LOG
                   for t in e.get("retired_text") or ()]
        assert retired, "claim set 1 retired a sentence; say which"

        checked = 0
        for path in _python_sources():
            checked += 1
            for lineno, text in _live_string_literals(path):
                for version, sentence in retired:
                    assert sentence not in text, (
                        f"{os.path.relpath(path, REPO)}:{lineno} emits a claim "
                        f"withdrawn by claim set {version} as a LIVE string:\n"
                        f"  {sentence}\n"
                        "Correct it in place (this repo replaces a stale claim "
                        "and quotes what it replaced); a comment or docstring "
                        "quoting it is fine, a string the product emits is not.")
        assert checked > 50, f"only {checked} sources scanned; the glob is wrong"

    def test_execute_run_stamps_the_provenance_it_writes(self):
        """Source-level, so it holds without a cached section.

        The end-to-end proof is
        ``test_app_server.py::test_result_carries_the_current_claim_set``.
        """
        src = _read("app", "server", "jobs.py")
        assert '"claim_set": claim_set_stamp(),' in src, (
            "execute_run must stamp the provenance dict it builds, beside the "
            "caveats the stamp is about")


# --------------------------------------------------------------------------- #
# 5. The Studio
# --------------------------------------------------------------------------- #
class TestTheStudioMarksAndDoesNotHide:
    def test_the_notice_is_built_in_one_place(self):
        src = _read("app", "web", "js", "ui", "provenance.js")
        assert "export function claimSetNotice(" in src
        assert "export function claimSetEraLabel(" in src
        # It renders nothing at all for a current, compliant record: an era
        # label on the only era there is says nothing, and a banner would be
        # false. Both findings gate both functions -- a record that contradicts
        # its own stamp must not slip through because it is not OLD.
        assert src.count(
            "if (!cs || !(cs.superseded || cs.contradicted)) return null;") == 2

    def test_the_banner_is_rendered_before_the_caveats_it_qualifies(self):
        src = _read("app", "web", "js", "ui", "inspect.js")
        banner = src.index("caveatsEl.appendChild(claimSetNotice(pv))")
        loop = src.index("cav.forEach((c, i) => caveatsEl.appendChild(")
        assert banner < loop, (
            "a reader who is about to be shown a withdrawn sentence has to be "
            "told before they read it")

    def test_every_stored_caveat_is_still_rendered(self):
        """The banner ADDS. It must never filter, slice or truncate the list."""
        src = _read("app", "web", "js", "ui", "inspect.js")
        assert ("cav.forEach((c, i) => caveatsEl.appendChild("
                "caveatNode(String(c), i + 1)));") in src, (
            "the render loop over the run's own caveats is unconditional and "
            "unfiltered; if it moved, re-check that it still is")
        assert "const cav = (pv && pv.caveats) || [];" in src
        # No filtering of the caveat list anywhere in the panel.
        assert "cav.filter" not in src

    def test_the_panel_says_it_before_the_panel_is_read(self):
        """Twice: the index strip above the list, and the section header --
        which is the only one still on screen once the panel is scrolled.

        And it names WHICH finding. A contradicted record labelled
        "superseded" would be a one-word summary that is false: it is not old,
        and the two findings call for different actions (read it as history
        vs. re-run it).
        """
        src = _read("app", "web", "js", "ui", "inspect.js")
        assert src.count("h('span.stale', staleWord)") == 2
        assert "tallyEl.appendChild(h('span.stale', staleWord))" in src
        assert ("const staleWord = cs && cs.contradicted "
                "? 'contradicted' : 'superseded';") in src

    def test_the_comparison_panel_marks_both_sides(self):
        """Two runs can sit on different claim sets, and then a caveat that
        shows up under "A only" is a change in the product's claims rather than
        a difference between the runs."""
        src = _read("app", "web", "js", "ui", "history.js")
        assert "claimSetNotice" in src
        assert "{ compact: true, label: side }" in src

    def test_the_banner_has_styling_in_both_themes(self):
        css = _read("app", "web", "css", "components.css")
        for sel in (".caveat.crit.cs-superseded", ".cs-change", ".cs-change-hd",
                    ".hr-lbl.cs-era"):
            assert sel in css, f"{sel} is rendered with no rule behind it"
        # Every colour is a token, so light and dark are covered by tokens.css
        # rather than by a second hand-written palette here.
        block = css[css.index("--- superseded claim set"):css.index(".cv-tally {")]
        assert "#" not in block and "rgb(" not in block and "hsl(" not in block, (
            "a literal colour here would be right in one theme only")


# --------------------------------------------------------------------------- #
# 5b. The report front — the copy that leaves the building
# --------------------------------------------------------------------------- #
class TestAReportDocumentIsMarkedToo:
    """``report.html`` is the file a customer KEEPS. It was the worst case.

    ``ReportManager._restore`` warned on a document with ``drugs[]`` -- a
    PANEL, the shape the ranking withdrawal changed. A SINGLE-COMPOUND document
    from the same era has no ``drugs[]``, took no warning at all, and
    re-rendered its stored caveats verbatim into ``report.html``. Everything
    else fixed here lives behind a server; this copy has already been
    downloaded.

    The condition is WIDENED, not moved: a panel document still takes its own
    warning, which says something about the document's SHAPE that a claim-set
    stamp cannot say.
    """

    @staticmethod
    def _doc(report_id: str, **over) -> dict:
        doc = {
            "report_id": report_id,
            "section": "xenium_319",
            "drug": "temozolomide",          # SINGLE compound: no `drugs[]`
            "status": "done",
            "spec_base": {"section": "xenium_319"},
            "spec_overrides": {},
            "arms": {"real": "r_00001", "counterfactual": "r_00002"},
            "warnings": [],
            "caveats": ["ABSOLUTE KILL MUST NOT BE QUOTED. " + RETIRED + "."],
        }
        doc.update(over)
        return doc

    @staticmethod
    def _restore(tmp_path, doc):
        from app.server.report import ReportManager

        d = os.path.join(str(tmp_path), "reports", doc["report_id"])
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "report.json"), "w") as fh:
            json.dump(doc, fh)
        jobs = JobManager(root=os.path.join(str(tmp_path), "runs"),
                          start_worker=False)
        mgr = ReportManager(jobs, root=os.path.join(str(tmp_path), "reports"))
        return mgr.get(doc["report_id"])

    def test_a_single_compound_document_was_exempt_and_is_not_any_more(self, tmp_path):
        """The case that was silently getting nothing. No ``drugs[]``."""
        doc = self._doc("rep_09001")
        assert "drugs" not in doc, "the whole point: this is not a panel"
        rep = self._restore(tmp_path, doc)

        assert not any("RESTORED FROM A WITHDRAWN PANEL" in w
                       for w in rep.warnings), (
            "it is not a panel document, so the panel warning must not claim "
            "it is -- this is a second warning, not a widened lie")
        assert any("SUPERSEDED CLAIM SET" in w for w in rep.warnings), (
            "a single-compound document from the retired era now says so")
        # And the caveats are restored WHOLE, retired sentence included.
        assert rep.body["caveats"] == doc["caveats"]

    def test_a_panel_document_still_takes_its_own_warning(self, tmp_path):
        """The failure mode of a widened condition is that it MOVES the thing
        it was meant to add to. Both fire, and they say different things."""
        doc = self._doc("rep_09002", drug=None,
                        drugs=["temozolomide", "carmustine"])
        rep = self._restore(tmp_path, doc)

        panel = [w for w in rep.warnings if "RESTORED FROM A WITHDRAWN PANEL" in w]
        claim = [w for w in rep.warnings if "SUPERSEDED CLAIM SET" in w]
        assert len(panel) == 1, "the pre-existing warning is untouched"
        assert len(claim) == 1, "and the new one is additional"
        assert "2-compound RANKING" in panel[0]

    def test_a_document_stamped_current_takes_no_claim_set_warning(self, tmp_path):
        from app.server.jobs import claim_set_stamp as _stamp

        doc = self._doc("rep_09003", claim_set=_stamp(),
                        caveats=["WHAT SURVIVES IS THE MAP."])
        rep = self._restore(tmp_path, doc)
        assert not any("SUPERSEDED CLAIM SET" in w or "CONTRADICTS" in w
                       for w in rep.warnings), rep.warnings

    def test_a_document_that_contradicts_its_stamp_is_caught(self, tmp_path):
        """Same two findings as a run: age, and compliance."""
        from app.server.jobs import claim_set_stamp as _stamp

        rep = self._restore(tmp_path, self._doc("rep_09004", claim_set=_stamp()))
        assert any("CONTRADICTS THE CLAIM SET" in w for w in rep.warnings)

    def test_an_assembled_document_carries_the_stamp(self):
        """``attribute_compound`` stamps what it assembles, beside the caveats
        it assembled -- so a document written today is distinguishable, to a
        reader months from now, from one written before the withdrawal."""
        src = _read("app", "server", "report.py")
        assert '"claim_set": claim_set_stamp(),' in src
        assert '"caveats",\n' in src or '"caveats"' in src
        # and it survives a restore rather than being dropped on the way back.
        assert '"claim_set", "window_s", "map_claim"' in src


# --------------------------------------------------------------------------- #
# 6. One vocabulary
# --------------------------------------------------------------------------- #
class TestOneVocabulary:
    def test_reports_already_solved_this_and_still_do(self):
        src = _read("app", "server", "report.py")
        assert "RESTORED FROM A WITHDRAWN PANEL DOCUMENT" in src, (
            "the report front's warning is the thing this was modelled on; if "
            "it moved, this module's wording has to move with it")

    def test_a_set_is_superseded_and_a_claim_is_withdrawn(self):
        """Two words, two meanings, used consistently on both fronts.

        A claim SET is superseded (a newer one exists). An individual claim is
        WITHDRAWN (it was measured and it does not hold). Letting these drift
        into synonyms is how one situation ends up with two vocabularies.
        """
        st = claim_set_status(None)
        assert "SUPERSEDED CLAIM SET" in st["headline"]
        assert all("withdraw" in e["withdrew"].lower() or
                   "WITHDRAWN" in e["withdrew"] for e in CLAIM_SET_LOG)
        js = _read("app", "web", "js", "ui", "provenance.js")
        assert "WITHDRAWN" in js and "SUPERSEDED" in js

    def test_both_fronts_cite_the_same_measurement_for_the_same_withdrawal(self):
        """``report.py`` names two files when it declines to re-render a panel
        ranking. Claim set 1 withdrew that same ranking. Different citations
        for one withdrawal would mean one of them is stale."""
        report_src = _read("app", "server", "report.py")
        ranking = [e for e in CLAIM_SET_LOG if "RANKING" in e["withdrew"]]
        assert ranking, "claim set 1 is the ranking withdrawal"
        for src in ranking[0]["sources"]:
            assert src in report_src, (
                f"claim set {ranking[0]['version']} cites {src}; report.py's "
                "withdrawal of the same claim does not")
