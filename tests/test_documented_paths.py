"""Every URL this repository tells a person to open must answer.

WHY THIS FILE EXISTS
--------------------
On 2026-10-10 ``/ingest`` -- the bring-your-own-tissue page, printed on the
startup banner and named in the README, QUICKSTART, DEPLOYMENT sec. 0, the
workspace contract, the home page and the project grid -- was served by a
~190-line HTML document built as a Python string inside ``app/server/api.py``.
``app/server/browse.py`` had replaced it with a directory browser at
``/tissue``, and the note left behind said that deleting the old route would let
"the static mount serve the new page".

**It would not.** ``StaticFiles(html=True)`` resolves a directory to its
``index.html`` and nothing else: it does not try ``<path>.html``, so with the
route gone ``GET /ingest`` is a 404 while ``GET /ingest.html`` is a 200. That
was measured before the cleanup, not assumed. A one-line deletion would have
404ed the one URL the product's own banner prints, and nothing in this
repository would have noticed -- there was no test over the set of paths the
documents promise.

So this file is that test. It is deliberately NOT a test of the ingest page: it
is a test of the PROMISE, over every page URL at once, because the failure mode
is structural (a route deleted, a router mounted after the catch-all, a page
renamed) and it recurs. ``app/server/api.py``'s own mounting comments say the
same thing three times in three places -- "a line added after it leaves the
route registered and never matched -- the kind of break that shows up as a 404
nobody can explain".

THE TWO HALVES, AND WHY NEITHER IS ENOUGH ALONE
-----------------------------------------------
:func:`test_every_documented_path_answers` is the guard. A table of paths that
answered is worth nothing on its own, though: the way a check like this dies is
that somebody renames a URL, updates the table, and the table stops describing
anything the documents say. So :func:`test_every_documented_path_is_actually_
documented` asserts that each row's ``needle`` really appears in the file named
as documenting it. Rename a URL and you have to change both the code and the
document -- which is the decision this is here to force -- instead of only the
table.
"""
from __future__ import annotations

import os
from typing import NamedTuple, Tuple

import pytest
from starlette.testclient import TestClient

from app.server.api import create_app
from app.server.jobs import JobManager

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


class Promised(NamedTuple):
    """A URL a document tells a person to open."""

    #: The path to request.
    path: str
    #: The file that promises it, repo-relative.
    doc: str
    #: The literal text to find in ``doc``. Usually the path; spelled
    #: separately for ``/`` and for a path that only appears inside a table
    #: cell, where the bare string would match any line.
    needle: str
    #: What a person expects to be looking at.
    what: str


#: Every page URL this repository publishes, and where it publishes it.
#:
#: Pages only. The ``/api/`` surface is gated, is documented in
#: ``API_CONTRACT.md`` and is covered route by route elsewhere; what has no
#: owner is the set of paths a HUMAN is handed -- a banner line, a link in a
#: runbook, a bookmark made from either. ``/healthz`` and ``/api/health`` are in
#: because a platform health check is configured from the deployment page and a
#: 404 there kills a container that is serving correctly, which has happened
#: here twice.
PROMISED: Tuple[Promised, ...] = (
    Promised("/", "docs/WORKSPACE_CONTRACT.md", "| `/` |",
             "the workspace home"),
    Promised("/studio", "docs/WORKSPACE_CONTRACT.md", "`/studio`",
             "the Studio"),
    Promised("/index.html", "docs/WORKSPACE_CONTRACT.md", "index.html",
             "the Studio at its original URL; every old bookmark"),
    Promised("/ingest", "app/server/__main__.py", "/ingest",
             "bring your own tissue -- THE STARTUP BANNER PRINTS THIS ONE"),
    Promised("/ingest", "README.md", "/ingest",
             "bring your own tissue, from the README"),
    Promised("/ingest", "docs/QUICKSTART.md", "/ingest",
             "bring your own tissue, from the quickstart"),
    Promised("/ingest", "docs/DEPLOYMENT.md", "/ingest",
             "bring your own tissue, from the deployment page"),
    Promised("/tissue", "API_CONTRACT.md", "`/tissue`",
             "the same page under a name that says what it is"),
    Promised("/login", "docs/AUTH_CONTRACT.md", "`/login`", "sign in"),
    Promised("/healthz", "docs/DEPLOYMENT.md", "`/healthz`",
             "liveness, for a platform health checker"),
    Promised("/api/health", "docs/DEPLOYMENT.md", "`/api/health`",
             "the richer health route"),
    Promised("/docs", "docs/DEPLOYMENT.md", "`/docs`", "interactive OpenAPI"),
    Promised("/redoc", "docs/DEPLOYMENT.md", "`/redoc`", "the other schema UI"),
    Promised("/openapi.json", "docs/DEPLOYMENT.md", "`/openapi.json`",
             "the machine-readable schema"),
)


@pytest.fixture(scope="module")
def served(tmp_path_factory):
    """The real app, built the way ``python -m app.server`` builds it.

    A module-scoped client over a throwaway run root: these are GETs of pages
    and nothing here writes, so one app is enough and building it once keeps
    the file fast.
    """
    root = tmp_path_factory.mktemp("promised_runs")
    app = create_app(JobManager(root=str(root)))
    with TestClient(app) as client:
        yield client


@pytest.mark.parametrize("row", PROMISED, ids=[
    f"{r.path}@{os.path.basename(r.doc)}" for r in PROMISED])
def test_every_documented_path_answers(served, row: Promised):
    """A documented URL that 404s is the defect, not the cleanup."""
    resp = served.get(row.path)
    assert resp.status_code != 404, (
        f"{row.path} is a 404. {row.doc} tells a person to open it "
        f"({row.what}). Either the route was removed or it is registered "
        f"AFTER the StaticFiles mount at '/', which swallows every path under "
        f"it -- see app/server/api.py::_mount_web. Note that "
        f"StaticFiles(html=True) does NOT serve an extensionless path: it "
        f"resolves a DIRECTORY to its index.html and does not try "
        f"'<path>.html', so deleting a page route does not leave the static "
        f"mount covering it.")
    assert resp.status_code == 200, (
        f"{row.path} answered {resp.status_code}; {row.doc} promises a page "
        f"there ({row.what})")


@pytest.mark.parametrize("row", PROMISED, ids=[
    f"{r.path}@{os.path.basename(r.doc)}" for r in PROMISED])
def test_every_documented_path_is_actually_documented(row: Promised):
    """The non-vacuity half. See this module's docstring.

    Without this, renaming a URL and updating the table below would leave a
    green check over a set of paths no document mentions.
    """
    full = os.path.join(_REPO, row.doc)
    assert os.path.exists(full), f"{row.doc} is gone; this row is stale"
    with open(full, encoding="utf-8") as fh:
        text = fh.read()
    assert row.needle in text, (
        f"{row.doc} no longer contains {row.needle!r}, so this row is no "
        f"longer a promise anybody made. Delete the row, or restore the "
        f"sentence -- do not change the needle to something that happens to "
        f"match.")


def test_the_pages_answer_anonymously_with_the_gate_on(tmp_path, monkeypatch):
    """A page must load in order to prompt for a credential.

    The gate exempts non-``/api/`` GETs precisely so this holds, and it is the
    reason ``/ingest`` can be on a banner that also prints the key. Checked on
    an app built WITH a key, because that is the deployment where getting it
    wrong is invisible locally: every one of these answers 200 on loopback
    whether or not the exemption is right.

    ``/api/health`` is in the list and ``/api/sections`` is the control -- a
    sweep that found everything open would have been testing an ungated app.
    """
    monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "a-key-no-request-here-sends")
    app = create_app(JobManager(root=str(tmp_path / "gated_runs")))
    with TestClient(app) as client:
        for path in ("/", "/studio", "/index.html", "/ingest", "/tissue",
                     "/login", "/healthz", "/api/health", "/docs", "/redoc",
                     "/openapi.json"):
            resp = client.get(path)
            assert resp.status_code == 200, (
                f"{path} answered {resp.status_code} to an anonymous caller. "
                f"A page that cannot load cannot prompt for a credential.")
            assert "a-key-no-request-here-sends" not in resp.text, path
        assert client.get("/api/sections").status_code == 401, (
            "the control failed: this app is not gated at all, so the sweep "
            "above proved nothing")
