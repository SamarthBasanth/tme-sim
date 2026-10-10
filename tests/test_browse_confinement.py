"""``GET /api/browse`` — the directory browser, and the allowlist that holds it.

WHY THIS FILE IS LONG
---------------------
Registering tissue means naming a directory ON THE SERVER'S OWN MACHINE: a
Xenium bundle is gigabytes of matrices read in place, and a section is
registered and persisted by path. The usability fix for "type that path from
memory" is therefore a server-side directory browser, and a route that lists
directories from a client-supplied path is a traversal and a disclosure
primitive the moment it is built carelessly. ``app/server/projects.py::_dir``
already did exactly that here: it joined a request-supplied id onto a path and
``%2e%2e`` turned it into a read-and-delete primitive.

So the confinement is not a side concern of this feature, it IS the feature,
and it is checked the way ``tests/test_scan_confinement.py`` checks the scan
route: **every escape below is attempted over HTTP and shown WORKING first.**
Nothing here asserts that a helper returned ``False`` — a guard that was
accidentally a no-op would sail through a test like that and fail these.

THE TRAP IN THE TEST HARNESS ITSELF
-----------------------------------
``httpx`` — which ``TestClient`` is built on — normalises ``..`` out of a URL
PATH on the client, before a byte is sent. A traversal test written as
``client.get("/api/browse/a/../b")`` therefore measures httpx's URL normaliser
and says NOTHING about the server:

    >>> import httpx
    >>> str(httpx.URL("http://t/api/browse/a/../b"))
    'http://t/api/browse/b'

That is pinned below in :func:`test_httpx_eats_a_dotdot_in_the_url_path` so the
next person does not spend an afternoon on a green test that tests a client
library. Two consequences, both acted on:

* This route takes the path as a QUERY parameter, where httpx leaves ``..``
  alone (also pinned below) — so the traversal tests in this file are real.
* The traversal is ALSO driven through a raw ASGI scope with a literal ``..``
  in it (:func:`_raw_get`), which has no client library in the path at all.

WHAT A LISTING IS ALLOWED TO CONTAIN
------------------------------------
Names and a type. Not sizes, not timestamps, not owners — those are functions
of a file's content or of a stat the caller is not entitled to — and not paths
on the child entries, so the one absolute path in a response is the directory
the caller itself asked about. Pinned in :class:`TestWhatAListingDiscloses`.
"""
from __future__ import annotations

import asyncio
import os
from typing import Dict, List, Optional
from urllib.parse import quote

import pytest

from app.server import cache as cache_mod
from app.server import browse as browse_mod
from tests.test_ingest import _SYNTH_MARKER_PANEL, _write_synthetic_xenium


# --------------------------------------------------------------------------- #
# The mount, and the one thing that can go wrong with it
# --------------------------------------------------------------------------- #
def _mount_browse(app) -> None:
    """``app.include_router(build_browse_router())``, in api.py's position.

    ``include_router`` APPENDS, and ``create_app``'s LAST act is
    ``_mount_web(app)``, which mounts ``StaticFiles`` at ``/`` — a mount that
    matches every path, so anything registered after it is unreachable.
    ``api.py`` puts all nine of its ``include_router`` calls before that line
    and says so in ``_mount_web``'s docstring; this helper reproduces that
    position for tests that must not edit ``api.py``.
    """
    from starlette.routing import Mount

    router = browse_mod.build_browse_router()
    for i, route in enumerate(app.router.routes):
        if isinstance(route, Mount) and route.path in ("", "/"):
            app.router.routes[i:i] = list(router.routes)
            return
    app.include_router(router)


@pytest.fixture
def client(tmp_path, monkeypatch):
    """The real app plus this router, with run history out of the repo."""
    from fastapi.testclient import TestClient

    from app.server import ingest_registry
    from app.server.api import create_app
    from app.server.jobs import JobManager

    runs = tmp_path / "runs"
    runs.mkdir()
    monkeypatch.setattr(cache_mod, "runs_dir", lambda: str(runs))
    monkeypatch.setattr(ingest_registry, "_store", None)
    monkeypatch.setattr(cache_mod, "_meta_cache", None)
    # Ambient configuration must not decide what these tests measure.
    monkeypatch.delenv(cache_mod.DATA_ROOTS_VAR, raising=False)
    monkeypatch.setattr(browse_mod, "_UNCONFINED_WARNED", True)
    app = create_app(JobManager(root=str(runs / "jobs")))
    _mount_browse(app)
    with TestClient(app) as c:
        c._app = app  # type: ignore[attr-defined]
        yield c
    monkeypatch.setattr(ingest_registry, "_store", None)


@pytest.fixture
def estate(tmp_path):
    """A root to allow, a sample inside it, and a sample OUTSIDE it.

    ``outside/`` holds the thing a tester with the key must not be able to
    enumerate, and the disclosure is the NAME: a specimen identifier is exactly
    what travels in a directory name.
    """
    allowed = tmp_path / "allowed"
    allowed.mkdir()
    _write_synthetic_xenium(str(allowed / "Study_A_S001"), _SYNTH_MARKER_PANEL)
    (allowed / "notes.txt").write_text("protocol notes\n")
    (allowed / ".hidden_dir").mkdir()

    outside = tmp_path / "outside"
    outside.mkdir()
    _write_synthetic_xenium(str(outside / "AT10_BRA_5_FO_1"), _SYNTH_MARKER_PANEL)

    # A sibling whose absolute path STARTS WITH the root's string. A prefix
    # test on the raw text calls this "inside"; it is not.
    sibling = tmp_path / "allowed_evil"
    sibling.mkdir()
    _write_synthetic_xenium(str(sibling / "Not_Yours"), _SYNTH_MARKER_PANEL)

    return {
        "root": str(allowed),
        "sample": str(allowed / "Study_A_S001"),
        "outside": str(outside),
        "sibling": str(sibling),
        "tmp": str(tmp_path),
    }


def _browse(client, path=None, **params):
    q = dict(params)
    if path is not None:
        q["path"] = str(path)
    return client.get("/api/browse", params=q)


def _names(resp) -> List[str]:
    return [e["name"] for e in resp.json()["entries"]]


# --------------------------------------------------------------------------- #
# A raw ASGI driver: no client library, no URL normaliser, literal bytes
# --------------------------------------------------------------------------- #
def _raw_get(app, raw_path: str, raw_query: bytes = b"",
             headers: Optional[List] = None) -> Dict:
    """Drive the ASGI app directly with the bytes we choose.

    ``TestClient`` cannot express a literal ``..`` in a URL path (see this
    module's docstring), so for the traversal probes there is no client at all:
    the scope is built here and handed straight to the application, which is
    exactly what a hostile peer on the wire can do. ``path`` and ``raw_path``
    are set independently and deliberately — starlette's ``get_route_path``,
    which both the router and the API-key gate dispatch on, reads ``path``.
    """
    received: List[Dict] = []

    async def receive():
        return {"type": "http.request", "body": b"", "more_body": False}

    async def send(message):
        received.append(message)

    scope = {
        "type": "http", "asgi": {"version": "3.0", "spec_version": "2.3"},
        "http_version": "1.1", "method": "GET", "scheme": "http",
        "path": raw_path, "raw_path": raw_path.encode("utf-8", "surrogateescape"),
        "root_path": "", "query_string": raw_query,
        "headers": list(headers or [(b"host", b"testserver")]),
        "client": ("127.0.0.1", 54321), "server": ("testserver", 80),
    }
    asyncio.get_event_loop_policy().new_event_loop().run_until_complete(
        app(scope, receive, send))
    status = next(m["status"] for m in received
                  if m["type"] == "http.response.start")
    body = b"".join(m.get("body", b"") for m in received
                    if m["type"] == "http.response.body")
    return {"status": status, "text": body.decode("utf-8", "replace")}


def test_httpx_eats_a_dotdot_in_the_url_path():
    """THE HARNESS FINDING. Read this before writing a traversal test here.

    ``httpx`` normalises ``..`` out of a URL PATH before sending, so a test
    written as ``client.get("/x/a/../b")`` asserts a property of httpx and
    nothing whatsoever about the server. Every traversal below therefore goes
    either through a QUERY parameter (where httpx leaves it alone — also pinned
    here) or through :func:`_raw_get`, which has no client in the path.
    """
    import httpx

    assert str(httpx.URL("http://t/api/browse/a/../b")) == "http://t/api/browse/b"
    # The query string is NOT normalised, which is what makes this route's
    # shape testable through TestClient at all.
    url = httpx.URL("http://t/api/browse", params={"path": "/a/../b"})
    assert b".." in url.query
    assert httpx.URL("http://t/api/browse?path=%2e%2e/x").params["path"] == "../x"


def test_a_literal_dotdot_reaches_the_server_through_a_raw_scope(client, estate):
    """Anti-vacuity for :func:`_raw_get`: the ``..`` really does arrive.

    If the raw driver quietly normalised too, every traversal test built on it
    would be measuring the same nothing the httpx ones would have.
    """
    escape = os.path.join(estate["root"], "..", "outside")
    out = _raw_get(client._app, "/api/browse",
                   b"path=" + quote(escape, safe="").encode())
    assert out["status"] == 200, out["text"]
    # Unconfined, so the traversal RESOLVES -- which is the proof it was not
    # stripped on the way in.
    assert "AT10_BRA_5_FO_1" in out["text"]


# --------------------------------------------------------------------------- #
# The default: unset means unconfined, and it says so
# --------------------------------------------------------------------------- #
def test_unset_leaves_browsing_unconfined(client, estate):
    """Backwards compatibility with the documented workflow is why unset is
    unconfined: a bring-your-own-tissue export is almost never inside this
    repo, and ``cache.py``'s block comment makes that an explicit choice.
    ``__main__.py`` confines to ``<repo>/data_cache`` the moment the bind is
    not loopback, so a REACHABLE instance is confined (pinned below)."""
    resp = _browse(client, estate["outside"])
    assert resp.status_code == 200
    assert resp.json()["confined"] is False
    assert _names(resp) == ["AT10_BRA_5_FO_1"]


def test_the_first_unconfined_browse_warns_once_naming_the_variable(
        client, estate, monkeypatch, capsys):
    """An unconfined default is a decision; a SILENT one is not defensible for
    a capability that lists any directory this process can read.

    Separate from ``cache._warn_unconfined``, which is worded for ``POST
    /api/sections/scan`` and fires at most once per process: if that one had
    already gone off, an operator would never be told the BROWSER is unconfined
    too. The capability named has to be the one exercised.
    """
    monkeypatch.setattr(browse_mod, "_UNCONFINED_WARNED", False)
    capsys.readouterr()

    assert _browse(client, estate["outside"]).status_code == 200
    first = capsys.readouterr().err
    assert "/api/browse" in first and "UNCONFINED" in first
    assert cache_mod.DATA_ROOTS_VAR in first

    # Once per process, not once per click: a browser makes one of these calls
    # per navigation, and a per-call warning is log spam an operator learns to
    # filter out.
    assert _browse(client, estate["outside"]).status_code == 200
    assert cache_mod.DATA_ROOTS_VAR not in capsys.readouterr().err


def test_a_non_loopback_bind_confines_this_route_too(monkeypatch):
    """The reachable case, at the one place that decides it.

    ``__main__.py`` sets ``TMESIM_STUDIO_DATA_ROOTS`` to ``<repo>/data_cache``
    for a non-loopback bind when the operator expressed no preference. This
    route reads the SAME variable through the SAME function, so there is
    nothing for it to opt into and nothing it can forget — which is the whole
    reason it does not have a confinement mechanism of its own.
    """
    repo = cache_mod.repo_root()
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR,
                       os.path.join(repo, "data_cache"))
    roots = cache_mod.data_roots()
    assert roots == [os.path.realpath(os.path.join(repo, "data_cache"))]
    # And the browser's starting points are exactly that, not the four
    # convenience anchors it offers when unconfined.
    assert [p["path"] for p in browse_mod._places(roots)] == roots


# --------------------------------------------------------------------------- #
# Real escapes, each shown working first
# --------------------------------------------------------------------------- #
def test_dotdot_traversal_out_of_a_root_is_refused(client, estate, monkeypatch):
    """``<root>/../outside`` — the classic, and the one a string prefix misses."""
    escape = os.path.join(estate["root"], "..", "outside")

    before = _browse(client, escape)
    assert before.status_code == 200
    assert _names(before) == ["AT10_BRA_5_FO_1"]

    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    after = _browse(client, escape)
    assert after.status_code == 403
    assert "AT10_BRA_5_FO_1" not in after.text
    assert '"entries"' not in after.text


def test_dotdot_traversal_is_refused_through_a_raw_scope_too(
        client, estate, monkeypatch):
    """The same escape with no client library anywhere in the path.

    This is the test the httpx finding exists for: it proves the refusal is the
    SERVER's and not a normaliser's.
    """
    escape = os.path.join(estate["root"], "..", "outside")
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    out = _raw_get(client._app, "/api/browse",
                   b"path=" + quote(escape, safe="").encode())
    assert out["status"] == 403, out["text"]
    assert "AT10_BRA_5_FO_1" not in out["text"]


def test_a_percent_encoded_dotdot_is_refused(client, estate, monkeypatch):
    """``%2e%2e`` — the encoding that turned ``projects.py::_dir`` into a
    read-and-delete primitive. It decodes to ``..`` before the handler sees it,
    which is exactly why the check is on the RESOLVED object and not on text.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    encoded = estate["root"] + "/%2e%2e/outside"
    out = _raw_get(client._app, "/api/browse",
                   b"path=" + quote(encoded, safe="%").encode())
    assert out["status"] == 403, out["text"]
    assert "AT10_BRA_5_FO_1" not in out["text"]


def test_a_double_encoded_dotdot_is_not_a_traversal_and_is_not_a_500(
        client, estate, monkeypatch):
    """``%252e%252e`` decodes ONCE, to the literal characters ``%2e%2e``.

    The failure to look for is a handler that decodes twice (or a filesystem
    layer that does), which would turn this back into ``..``. It must stay a
    directory name that does not exist — answered, never a 500.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    out = _raw_get(client._app, "/api/browse",
                   b"path=" + quote(estate["root"] + "/%252e%252e/outside",
                                    safe="%").encode())
    assert out["status"] in (403, 404), out["text"]
    assert "AT10_BRA_5_FO_1" not in out["text"]


def test_an_overlong_utf8_dotdot_is_not_a_dotdot(client, estate, monkeypatch):
    """``%c0%ae`` is the overlong UTF-8 encoding of ``.`` that defeated a
    generation of path filters. It must not become ``.`` here, and must not
    crash the decoder either."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    for seq in ("%c0%ae%c0%ae", "%e0%80%ae%e0%80%ae"):
        out = _raw_get(client._app, "/api/browse",
                       b"path=" + quote(estate["root"] + "/" + seq + "/outside",
                                        safe="%").encode())
        assert out["status"] in (403, 404), (seq, out["text"])
        assert "AT10_BRA_5_FO_1" not in out["text"]


def test_a_nul_byte_is_answered_and_echoes_nothing(client, estate, monkeypatch):
    """``os.path.realpath`` raises ``ValueError: embedded null byte``.

    Uncaught that is a 500 on a route anyone with the key can reach, and the
    traceback is the disclosure. It is answered with a 400 that carries no part
    of the input: a byte string that broke the path API is not something to
    echo back.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    for target in (estate["root"] + "/%00", "%00", estate["root"] + "/x%00.csv"):
        out = _raw_get(client._app, "/api/browse",
                       b"path=" + quote(target, safe="%").encode())
        assert out["status"] == 400, (target, out["text"])
        assert "\x00" not in out["text"]
        assert estate["root"] not in out["text"]


def test_an_overlong_path_is_answered_not_crashed(client, estate, monkeypatch):
    """Longer than PATH_MAX. ``realpath`` raises ``OSError(ENAMETOOLONG)`` on
    some platforms and simply returns the string on others; both are answers,
    neither is a 500."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    deep = estate["root"] + ("/" + "a" * 250) * 40
    resp = _browse(client, deep)
    assert resp.status_code in (400, 403, 404), resp.text


def test_a_tilde_is_expanded_before_the_check(client, estate, monkeypatch):
    """``~`` is what an operator types by hand, and it must not sidestep the
    allowlist. The refusal names the caller's own ``~/...`` text and never the
    home directory it resolved to."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    monkeypatch.setenv("HOME", estate["outside"])
    resp = _browse(client, "~/AT10_BRA_5_FO_1")
    assert resp.status_code == 403
    assert resp.json()["detail"]["path"] == "~/AT10_BRA_5_FO_1"
    assert estate["outside"] not in resp.text


def test_a_relative_path_resolves_against_the_servers_cwd(
        client, estate, monkeypatch):
    """A relative path is the form ``docs/DEPLOYMENT.md`` §1d's own example
    uses. It resolves, so it must be checked the same way."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    monkeypatch.chdir(estate["root"])
    assert _browse(client, "Study_A_S001").status_code == 200
    assert _browse(client, os.path.join("..", "outside")).status_code == 403


def test_a_symlink_out_of_a_root_is_refused(client, estate, monkeypatch):
    """A link INSIDE the allowlist pointing out of it.

    The path TEXT is entirely inside the root, so only resolving the object
    catches this. This is the escape that a ``startswith`` check and an
    ``abspath`` check both wave through.
    """
    link = os.path.join(estate["root"], "shortcut")
    os.symlink(estate["outside"], link)

    before = _browse(client, link)
    assert before.status_code == 200
    assert _names(before) == ["AT10_BRA_5_FO_1"]

    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    after = _browse(client, link)
    assert after.status_code == 403
    assert "AT10_BRA_5_FO_1" not in after.text


def test_dotdot_through_a_symlink_is_judged_on_the_real_object(
        client, estate, monkeypatch):
    """``<root>/link/../…`` — where a TEXTUAL reading and reality disagree.

    ``os.path.abspath`` collapses ``..`` textually, so for this path it reports
    a directory inside the root while the object the path really designates is
    outside it. A check built on ``abspath`` rather than on ``realpath`` gets
    this wrong, and ``cache.py``'s own comment records that its first version
    did.
    """
    os.symlink(estate["outside"], os.path.join(estate["root"], "hop"))
    escape = os.path.join(estate["root"], "hop", "..", "outside",
                          "AT10_BRA_5_FO_1")

    # The textual reading of that path is inside the root...
    assert os.path.abspath(escape).startswith(estate["root"] + os.sep)
    # ...and the real one is not.
    assert not os.path.realpath(escape).startswith(estate["root"] + os.sep)

    assert _browse(client, escape).status_code == 200  # unconfined: answers
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    refused = _browse(client, escape)
    assert refused.status_code == 403
    assert "cells.csv.gz" not in refused.text

    # The same divergence the other way: text that reads OUTSIDE but really
    # designates a sample INSIDE is allowed, and what gets listed is the real
    # directory rather than the text.
    os.symlink(estate["sample"], os.path.join(estate["outside"], "back"))
    inbound = os.path.join(estate["root"], "..", "outside", "back")
    assert not os.path.abspath(inbound).startswith(estate["root"] + os.sep)
    ok = _browse(client, inbound)
    assert ok.status_code == 200, ok.text
    assert ok.json()["path"] == os.path.realpath(estate["sample"])


def test_the_root_itself_is_allowed_and_a_prefix_sibling_is_not(
        client, estate, monkeypatch):
    """The exact-root boundary, in both directions.

    ``<tmp>/allowed`` contains itself and must NOT contain
    ``<tmp>/allowed_evil``, whose absolute path starts with the same string.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])

    at_root = _browse(client, estate["root"])
    assert at_root.status_code == 200
    assert "Study_A_S001" in _names(at_root)
    assert _browse(client, estate["sample"]).status_code == 200

    sibling = _browse(client, estate["sibling"])
    assert sibling.status_code == 403
    assert "Not_Yours" not in sibling.text
    # And the parent of the root is not the root.
    assert _browse(client, estate["tmp"]).status_code == 403


def test_a_set_but_unusable_value_fails_closed(client, estate, monkeypatch):
    """``-e TMESIM_STUDIO_DATA_ROOTS=`` must not mean "unconfined".

    An operator who set the variable asked for confinement; a value that parses
    to no root confines browsing to nothing — including the sample that would
    otherwise be legal — and the listing-less response says which state it is
    in, so a client can distinguish a misconfiguration from an empty disk.
    """
    for value in ("", "   ", os.pathsep * 3):
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, value)
        assert _browse(client, estate["sample"]).status_code == 403, value
        start = _browse(client)
        assert start.status_code == 200
        body = start.json()
        assert body["confined"] is True and body["places"] == []


# --------------------------------------------------------------------------- #
# The symlink that is merely LISTED: no stat of anything outside the roots
# --------------------------------------------------------------------------- #
def test_listing_a_symlink_does_not_stat_what_it_points_at(
        client, estate, monkeypatch):
    """``is_dir()`` FOLLOWS by default, and that is the bug.

    For a link pointing out of the allowlist, following it is a ``stat`` of
    something outside the roots: it answers "does that path exist, and is it a
    directory?" — the question the allowlist exists to refuse. So a symlink is
    reported as ``kind: "link"`` and nothing is resolved. The server decides
    when it is OPENED, through the same ``confine_scan_path`` as everything
    else, which is the second half of this test.
    """
    os.symlink(estate["outside"], os.path.join(estate["root"], "escape"))
    os.symlink(os.path.join(estate["tmp"], "nothing_here"),
               os.path.join(estate["root"], "dangling"))
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])

    resp = _browse(client, estate["root"])
    assert resp.status_code == 200, resp.text
    kinds = {e["name"]: e["kind"] for e in resp.json()["entries"]}
    assert kinds["escape"] == "link"
    # A DANGLING link is the clean proof nothing was followed: a stat would
    # have raised, and a handler that guessed would have called it a file.
    assert kinds["dangling"] == "link"
    # Nothing about either target is in the body.
    assert "outside" not in resp.text and "nothing_here" not in resp.text

    assert _browse(client, os.path.join(estate["root"], "escape")
                   ).status_code == 403
    assert _browse(client, os.path.join(estate["root"], "dangling")
                   ).status_code == 403


# --------------------------------------------------------------------------- #
# The refusal is a diagnostic, and it is not an information leak
# --------------------------------------------------------------------------- #
def test_the_refusal_names_the_path_the_variable_and_the_way_out(
        client, estate, monkeypatch):
    """Same ``{code, message, path}`` shape ``/api/sections/scan`` refuses
    with, because it IS that refusal — ``cache.PathNotAllowed.as_detail``."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    detail = _browse(client, estate["outside"]).json()["detail"]

    assert detail["code"] == "path_outside_data_roots"
    assert detail["path"] == estate["outside"]
    msg = detail["message"]
    assert estate["outside"] in msg
    assert "resolved" in msg
    assert cache_mod.DATA_ROOTS_VAR in msg
    assert "restart" in msg


def test_the_refusal_never_discloses_the_allowed_roots(client, estate,
                                                       monkeypatch):
    """A refused caller learns that it was refused, not what would be accepted.

    This is the property the ROOT DISCLOSURE on the places endpoint must not
    have widened: a browser needs to be told where it may start, and a refusal
    still must not become an oracle for the allowlist.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR,
                       os.pathsep.join([estate["root"], estate["sibling"]]))
    body = _browse(client, estate["outside"]).text
    assert estate["root"] not in body
    assert estate["sibling"] not in body


def test_the_scan_routes_refusal_is_still_not_an_oracle_for_the_roots(
        client, estate, monkeypatch):
    """Nothing about this feature may widen the route that existed before it.

    ``GET /api/browse`` discloses the roots deliberately (you cannot navigate
    to a place you are not told exists). ``POST /api/sections/scan`` does not,
    and that is checked HERE as well as in ``test_scan_confinement.py``,
    because the tempting way to implement the browser is to make the scan
    refusal helpful by naming the allowlist.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR,
                       os.pathsep.join([estate["root"], estate["sibling"]]))
    body = client.post("/api/sections/scan",
                       json={"path": estate["outside"]})
    assert body.status_code == 403
    assert estate["root"] not in body.text
    assert estate["sibling"] not in body.text


def test_a_path_that_is_not_a_directory_is_a_404_naming_only_your_own_text(
        client, estate, monkeypatch):
    """The ``~/x`` rule again, on the other refusal: the message must not carry
    the path the server RESOLVED, only the one the caller sent."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    monkeypatch.setenv("HOME", estate["root"])
    resp = _browse(client, "~/notes.txt")
    assert resp.status_code in (403, 404), resp.text
    assert estate["root"] not in resp.text

    inside = _browse(client, os.path.join(estate["root"], "notes.txt"))
    assert inside.status_code == 404
    assert inside.json()["detail"]["code"] in ("not_a_directory",
                                               "no_such_directory")


# --------------------------------------------------------------------------- #
# What a listing discloses, and what it must never
# --------------------------------------------------------------------------- #
class TestWhatAListingDiscloses:
    """Names and a type. The rest is a stat the caller is not entitled to."""

    def test_an_entry_carries_a_name_and_a_kind_and_nothing_else(
            self, client, estate, monkeypatch):
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        resp = _browse(client, estate["sample"])
        assert resp.status_code == 200, resp.text
        entries = resp.json()["entries"]
        assert entries
        for e in entries:
            assert set(e) == {"name", "kind"}, e

    def test_no_stat_field_appears_anywhere_in_a_listing(
            self, client, estate, monkeypatch):
        """A size is a function of CONTENT. A browser that reports one has
        quietly become a read primitive over every directory it can walk.

        Checked over the KEYS of the whole payload, at every depth, rather
        than over the response text: a substring sweep over the text is
        defeated by any path that happens to contain the word, which is how
        this test first failed against the temporary directory pytest named
        after it.
        """
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        body = _browse(client, estate["sample"]).json()

        def keys(node):
            if isinstance(node, dict):
                for k, v in node.items():
                    yield k
                    yield from keys(v)
            elif isinstance(node, list):
                for item in node:
                    yield from keys(item)

        found = {k.lower() for k in keys(body)}
        for forbidden in ("size", "st_size", "bytes", "mtime", "modified",
                          "ctime", "atime", "uid", "gid", "mode", "owner",
                          "permissions"):
            assert forbidden not in found, forbidden

    def test_only_the_directory_the_caller_asked_about_carries_a_path(
            self, client, estate, monkeypatch):
        """Every absolute path in the body is either the caller's own resolved
        directory, one of its breadcrumbs, or an allowed root. Child entries
        have none at all, so the client composes the next path from a name."""
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        body = _browse(client, estate["root"]).json()
        allowed_paths = ({body["path"]}
                         | {c["path"] for c in body["crumbs"]}
                         | {p["path"] for p in body["places"]})
        assert allowed_paths == {os.path.realpath(estate["root"])}
        assert all("path" not in e for e in body["entries"])

    def test_dotfiles_are_hidden_until_asked_for(self, client, estate,
                                                 monkeypatch):
        """``.DS_Store`` and ``.git`` between a person and their sample is
        noise. It is a default, not a boundary — the flag turns it off."""
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        assert ".hidden_dir" not in _names(_browse(client, estate["root"]))
        assert ".hidden_dir" in _names(
            _browse(client, estate["root"], show_hidden="1"))

    def test_directories_sort_before_files(self, client, estate, monkeypatch):
        """A file list where the way onward is interleaved with the contents
        is a list you have to read rather than scan."""
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        kinds = [e["kind"] for e in _browse(client, estate["root"]).json()["entries"]]
        assert kinds == sorted(kinds, key=lambda k: k != "dir")

    def test_the_filter_runs_on_the_server_before_the_cap(
            self, client, estate, monkeypatch):
        """A client-side filter over a truncated listing finds nothing in the
        part that was cut off, which is the part a big directory is mostly made
        of."""
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        assert _names(_browse(client, estate["root"], q="study")) == [
            "Study_A_S001"]
        assert _names(_browse(client, estate["root"], q="zzz")) == []

    def test_a_huge_directory_is_capped_and_says_so(self, client, estate,
                                                    monkeypatch, tmp_path):
        """A navigation click must not be able to produce an unbounded
        response, and the response must not pretend the directory is small."""
        big = tmp_path / "allowed" / "many"
        big.mkdir()
        for i in range(browse_mod.MAX_ENTRIES + 25):
            (big / f"f{i:05d}.csv").write_text("x")
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        body = _browse(client, str(big)).json()
        assert len(body["entries"]) == browse_mod.MAX_ENTRIES
        assert body["truncated"] is True
        assert body["n_entries"] == browse_mod.MAX_ENTRIES + 25
        # And the per-child detection walk is skipped rather than run over a
        # directory too big to list.
        assert body["diagnosis"] is None


# --------------------------------------------------------------------------- #
# "Versatile": the diagnosis the registry already produced, surfaced
# --------------------------------------------------------------------------- #
class TestTheDiagnosisIsSurfaced:
    """A bare "not ingestable" where the engine had a sentence is the defect."""

    def test_a_recognised_bundle_is_marked_by_name_in_its_parent(
            self, client, estate, monkeypatch):
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        samples = _browse(client, estate["root"]).json()["diagnosis"]["samples"]
        assert "Study_A_S001" in samples
        assert samples["Study_A_S001"]["platform"] == "xenium"
        assert samples["Study_A_S001"]["confidence"] >= 0.5
        assert samples["Study_A_S001"]["registered_as"] is None

    def test_a_recognised_bundle_reports_itself_as_readable(
            self, client, estate, monkeypatch):
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        d = _browse(client, estate["sample"]).json()["diagnosis"]
        assert d["platform"] == "xenium" and d["ingestable"] is True
        assert d["layout"] == "xenium_bundle_extracted"

    def test_an_unrecognised_directory_says_WHAT_IT_EXPECTED(
            self, client, estate, monkeypatch, tmp_path):
        """The sentence that was reaching a console and no screen."""
        empty = tmp_path / "allowed" / "not_a_sample"
        empty.mkdir()
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        d = _browse(client, str(empty)).json()["diagnosis"]
        assert d["ingestable"] is False
        joined = " ".join(d["issues"])
        assert "no recognised layout" in joined
        # Not a bare failure: it names every layout it would have accepted.
        for expected in ("Xenium", "Space Ranger", "MEX"):
            assert expected in joined

    def test_a_half_built_bundle_names_the_missing_file(
            self, client, estate, monkeypatch, tmp_path):
        """The commonest real failure is an interrupted copy, and "no
        cell-feature matrix" is actively misleading about it — the directory is
        right there. ``detect`` already names the missing FILE; this proves the
        naming survives to the wire."""
        broken = tmp_path / "allowed" / "half_copied"
        _write_synthetic_xenium(str(broken), _SYNTH_MARKER_PANEL)
        mdir = broken / "cell_feature_matrix"
        for name in os.listdir(mdir):
            if name.startswith("barcodes"):
                os.remove(mdir / name)
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        d = _browse(client, str(broken)).json()["diagnosis"]
        assert d["ingestable"] is False
        joined = " ".join(d["issues"])
        assert "INCOMPLETE" in joined and "barcodes.tsv" in joined

    def test_the_diagnosis_carries_no_absolute_path(
            self, client, estate, monkeypatch):
        """``detect``'s container NOTE is written as f"{path} holds no sample
        of its own...", with ``path`` absolute.

        Browsing is a read of NAMES, and a diagnosis that smuggled a resolved
        absolute path back out in prose would defeat that — a caller that sent
        a symlink would be told where it points. The stripped form is also the
        better sentence, which is why this is not a redaction that costs
        anything.
        """
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        d = _browse(client, estate["root"]).json()["diagnosis"]
        joined = " ".join(d["issues"] + d["notes"])
        assert "holds no sample of its own" in joined   # the answer survived
        assert "this directory holds no sample" in joined
        assert estate["root"] not in joined             # the path did not
        assert estate["tmp"] not in joined

    def test_the_stripper_handles_the_subdirectory_form_detect_also_writes(
            self, monkeypatch, tmp_path):
        """``_detect_visium`` writes f"{sp}: no scalefactors_json.json" where
        ``sp`` is ``<path>/spatial`` — a path UNDER the directory, not the
        directory itself.

        That branch of ``_detect_visium`` is currently unreachable (``
        _is_spatial_dir`` already demands the file the issue reports missing),
        so this drives :func:`ingest_registry.diagnose` over the shape
        directly rather than claiming an end-to-end path that does not exist.
        The shape is in the engine and one edit from firing; the stripper has
        to be right for it now.
        """
        from app.server import ingest_registry

        here = str(tmp_path)
        monkeypatch.setattr(
            "tmesim.data.ingest.detect",
            lambda path: {
                "path": here, "platform": "visium", "layout": "spaceranger_outs",
                "confidence": 0.5, "ingestable": False, "files": {},
                "issues": [f"{os.path.join(here, 'spatial')}: no "
                           "scalefactors_json.json - spot_diameter_fullres is "
                           "the only way to convert pixels to micrometres."],
                "notes": [f"{here} is the sample root"],
                "contained_samples": [],
            })
        d = ingest_registry.diagnose(here)
        assert d["issues"] == [
            "spatial: no scalefactors_json.json - spot_diameter_fullres is "
            "the only way to convert pixels to micrometres."]
        assert d["notes"] == ["this directory is the sample root"]
        assert here not in " ".join(d["issues"] + d["notes"])

    def test_a_directory_already_registered_says_so(
            self, client, estate, monkeypatch):
        """Keyed on ``realpath`` at both ends, so the same directory reached
        through a symlink is recognised as the section it already is — the
        difference between a browser that offers to register a sample twice and
        one that tells you it is done."""
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        reg = client.post("/api/sections/scan", json={
            "path": estate["sample"], "register": True,
            "section_id": "study_a_s001"})
        assert reg.status_code == 200 and reg.json()["ok"], reg.text
        try:
            d = _browse(client, estate["sample"]).json()["diagnosis"]
            assert d["registered_as"] == "study_a_s001"
            parent = _browse(client, estate["root"]).json()["diagnosis"]
            assert parent["samples"]["Study_A_S001"]["registered_as"] == \
                "study_a_s001"

            # Reached through a link inside the root: still the same section.
            os.symlink(estate["sample"], os.path.join(estate["root"], "same"))
            via = _browse(client, os.path.join(estate["root"], "same"))
            assert via.json()["diagnosis"]["registered_as"] == "study_a_s001"
        finally:
            client.delete("/api/sections/study_a_s001")

    def test_a_container_of_samples_is_not_itself_offered_as_one(
            self, client, estate, monkeypatch):
        """``<root>`` holds a sample but is not one. The old form answered that
        with a failure; the browser answers it by pointing at the children."""
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
        d = _browse(client, estate["root"]).json()["diagnosis"]
        assert d["ingestable"] is False
        assert d["platform"] == "multi_sample"
        assert list(d["samples"]) == ["Study_A_S001"]


# --------------------------------------------------------------------------- #
# Bring your own tissue still works — the product's main input
# --------------------------------------------------------------------------- #
def test_browse_then_register_in_one_action_under_an_allowed_root(
        client, estate, monkeypatch):
    """The whole flow, with the guard on, the way the page drives it.

    Note what registration is NOT: a new route. The browser finds the
    directory and ``POST /api/sections/scan {register: true}`` — the endpoint
    the CLI and the original form already used, with the confinement it already
    had — does the work. This feature adds a way to FIND a directory and not
    one new thing you can do to it.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])

    start = _browse(client).json()
    assert [p["path"] for p in start["places"]] == [
        os.path.realpath(estate["root"])]
    assert start["entries"] == [] and start["path"] is None

    listing = _browse(client, start["places"][0]["path"]).json()
    name = next(e["name"] for e in listing["entries"]
                if listing["diagnosis"]["samples"].get(e["name"]))
    child = os.path.join(listing["path"], name)   # composed from a NAME

    reg = client.post("/api/sections/scan", json={
        "path": child, "register": True, "section_id": "study_a_s001",
        "label": "Study A S001"})
    assert reg.status_code == 200, reg.text
    assert reg.json()["ok"], reg.json().get("error")
    try:
        row = next(r for r in client.get("/api/sections").json()
                   if r["id"] == "study_a_s001")
        assert row["available"] and row["ingested"]
        assert row["cache_dir"] == os.path.realpath(estate["sample"])
    finally:
        client.delete("/api/sections/study_a_s001")


def test_the_crumbs_never_offer_a_path_the_server_would_refuse(
        client, estate, monkeypatch):
    """Breadcrumbs are clipped at the containing root.

    Offering ``/`` on a confined instance is a row of links that all 403 — an
    interface advertising a capability the server will refuse. Every crumb is
    re-requested here, because "they look clipped" is not the property; "they
    all work" is.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    body = _browse(client, estate["sample"]).json()
    assert [c["name"] for c in body["crumbs"]] == ["allowed", "Study_A_S001"]
    for crumb in body["crumbs"]:
        assert _browse(client, crumb["path"]).status_code == 200, crumb

    # Unconfined there is nothing to clip, so the whole ancestry is offered.
    monkeypatch.delenv(cache_mod.DATA_ROOTS_VAR)
    assert _browse(client, estate["sample"]).json()["crumbs"][0]["name"] == "/"


# --------------------------------------------------------------------------- #
# The mount
# --------------------------------------------------------------------------- #
def test_the_router_mounts_with_one_line_before_the_static_catch_all(tmp_path,
                                                                     monkeypatch):
    """``app.include_router(build_browse_router())``, and WHERE it goes.

    ``include_router`` appends; ``create_app``'s last act mounts
    ``StaticFiles`` at ``/``, which matches every path. A line added AFTER that
    one is not a smaller mistake — the route exists and is simply never
    matched, which fails silently. Both halves are pinned so the next reader
    knows the ordering is the mechanism and not a style preference.
    """
    from fastapi.testclient import TestClient

    from app.server import ingest_registry
    from app.server.api import create_app
    from app.server.jobs import JobManager

    runs = tmp_path / "runs"
    runs.mkdir()
    monkeypatch.setattr(cache_mod, "runs_dir", lambda: str(runs))
    monkeypatch.setattr(ingest_registry, "_store", None)
    monkeypatch.delenv(cache_mod.DATA_ROOTS_VAR, raising=False)
    monkeypatch.setattr(browse_mod, "_UNCONFINED_WARNED", True)

    # THE "TOO LATE" HALF IS BUILT ON A BARE APP, NOT ON create_app.
    #
    # It used to append the router to a real `create_app()` and assert 404.
    # That stopped being observable the moment api.py mounted the router
    # properly: create_app now serves /api/browse itself, so a second
    # appended copy is simply never reached and the assertion saw 200.
    #
    # The premise had disappeared, not the mechanism. So this half now
    # reproduces the mechanism directly — a catch-all mount, then an appended
    # router — which is what the ordering claim is actually about, and which
    # stays true however api.py is organised later.
    import os as _os

    from fastapi import FastAPI
    from starlette.staticfiles import StaticFiles

    web = _os.path.join(
        _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
        "app", "web")
    too_late = FastAPI()
    too_late.mount("/", StaticFiles(directory=web, html=True), name="web")
    too_late.include_router(browse_mod.build_browse_router())
    with TestClient(too_late) as c:
        assert c.get("/api/browse").status_code == 404, (
            "a router appended AFTER a StaticFiles catch-all at / was "
            "matched; the ordering claim in api.py is then wrong and the "
            "comment there should be corrected rather than trusted")

    # And the real app, which must mount it before that catch-all.
    in_place = create_app(JobManager(root=str(runs / "b")))
    with TestClient(in_place) as c:
        assert c.get("/api/browse").status_code == 200, (
            "create_app does not serve /api/browse — either the "
            "include_router line is missing or it sits after _mount_web")


def test_the_page_is_served_and_names_no_transfer_of_tissue(client):
    """``/tissue`` serves the page, and ``/ingest`` is untouched.

    The original self-contained form is registered on the app directly, before
    any router, so it keeps winning its path and every bookmark to it survives.

    The page also must not describe registering tissue as sending it anywhere.
    The bytes are read in place and a section is kept by path; a transfer
    metaphor on this screen would be a lie about what the product does, and a
    5 GB one if anybody believed it.
    """
    page = client.get("/tissue")
    assert page.status_code == 200
    assert "text/html" in page.headers["content-type"]
    body = page.text
    assert "/api/browse" in body or "js/ingest/main.js" in body
    assert client.get("/ingest").status_code == 200

    web = os.path.join(cache_mod.repo_root(), "app", "web")
    sources = [os.path.join(web, "ingest.html")]
    js = os.path.join(web, "js", "ingest")
    sources += [os.path.join(js, n) for n in sorted(os.listdir(js))
                if n.endswith(".js")]
    for path in sources:
        with open(path) as fh:
            text = fh.read()
        assert "upload" not in text.lower(), path
        assert "drag" not in text.lower(), path


# --------------------------------------------------------------------------- #
# The page, against a REAL server
#
# A static-file test over app/web would exercise module loading and the first
# render, and would be blind to the thing that actually matters: whether the
# listing a person clicks through is the server's. So this runs the real app on
# a real port with the allowlist SET -- the configuration a reachable instance
# has -- and drives the page the way a visitor does, from the starting point
# the server offers to a registered section.
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def live(tmp_path_factory):
    """The real app on a real port, confined, with its state out of the repo."""
    import threading
    import time

    uvicorn = pytest.importorskip("uvicorn")
    pytest.importorskip("playwright.sync_api",
                        reason="the tissue page test needs playwright")

    from app.server.api import create_app
    from app.server.jobs import JobManager

    state = tmp_path_factory.mktemp("state")
    tissue = tmp_path_factory.mktemp("tissue")
    _write_synthetic_xenium(str(tissue / "Study_B_S002"), _SYNTH_MARKER_PANEL)

    # Module-scoped, so monkeypatch is unavailable; the environment is saved
    # and restored by hand rather than left changed for the rest of the run.
    keys = ("TMESIM_STATE_DIR", cache_mod.DATA_ROOTS_VAR,
            "TMESIM_STUDIO_API_KEY", "TMESIM_STUDIO_ALLOWED_HOSTS")
    saved = {k: os.environ.get(k) for k in keys}
    os.environ["TMESIM_STATE_DIR"] = str(state)
    os.environ[cache_mod.DATA_ROOTS_VAR] = str(tissue)
    os.environ.pop("TMESIM_STUDIO_API_KEY", None)
    os.environ.pop("TMESIM_STUDIO_ALLOWED_HOSTS", None)

    manager = JobManager(root=str(state / "jobs"))
    app = create_app(manager)
    _mount_browse(app)
    server = uvicorn.Server(uvicorn.Config(
        app, host="127.0.0.1", port=0, log_level="error"))
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()
    deadline = time.time() + 30
    while not getattr(server, "started", False) and time.time() < deadline:
        time.sleep(0.05)
    assert server.started, "the test server did not come up"
    port = server.servers[0].sockets[0].getsockname()[1]
    try:
        yield {"base": f"http://127.0.0.1:{port}", "root": str(tissue)}
    finally:
        server.should_exit = True
        thread.join(timeout=15)
        manager.shutdown()
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def test_the_page_walks_the_server_to_a_registered_section(live):
    """The whole interaction, in a browser, with no path typed anywhere.

    This is the brief: a person who could previously only register tissue by
    recalling an absolute path now reaches the same outcome by clicking. The
    assertions follow that sequence, and the last one is that the browser
    reported no page error — which is also the guard against a reversed `h()`
    selector, because `h('div.a#b')` THROWS and takes its whole subtree with
    it, so a pane that renders nothing looks exactly like a page still
    loading.
    """
    from playwright.sync_api import sync_playwright

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        errors = []
        console = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        # Captured so a failure below reports what the page SAID rather than
        # only that a selector never appeared.
        page.on("console", lambda m: console.append(f"{m.type}: {m.text}"))
        try:
            page.goto(live["base"] + "/tissue", wait_until="domcontentloaded")

            # 1. The server offers its allowed root as a starting point. The
            #    allowlist is SET here, so this is the confined configuration.
            start = page.locator(".tb-place").first
            start.wait_for(timeout=15000)
            # Lower-cased before the comparison: the label is set in small
            # caps by `text-transform`, so `inner_text` returns it uppercased.
            assert "allowed on this server" in page.inner_text(
                "#tb-places").lower()

            # 2. Clicking it lists the directory, and the child that is a
            #    recognised bundle is MARKED as one.
            start.click()
            row = page.locator(".tb-row", has_text="Study_B_S002")
            row.wait_for(timeout=15000)
            assert "xenium" in row.inner_text().lower()

            # 3. Opening it says what it is before anything is pressed.
            row.click()
            page.wait_for_selector("#rp-register:not([disabled])", timeout=15000)
            verdict = page.inner_text("#rp-verdict").lower()
            assert "can be read" in verdict
            assert "xenium_bundle_extracted" in verdict

            # 4. One action registers it.
            page.click("#rp-register")
            try:
                page.wait_for_selector(".rp-ok", timeout=120000)
            except Exception:  # noqa: BLE001 - re-raised with the screen's text
                raise AssertionError(
                    "registration did not report success. The panel said: "
                    + page.inner_text("#rp-out") + " | console: "
                    + " || ".join(console))
            assert "added as" in page.inner_text(".rp-ok").lower()

            # 5. ...and the server agrees.
            sections = page.request.get(live["base"] + "/api/sections").json()
            assert "study_b_s002" in [s["id"] for s in sections]

            # 6. Re-reading the directory now says it is already a section,
            #    and the action is withheld rather than offered twice.
            page.click("#tb-reload")
            page.wait_for_selector("#rp-register[disabled]", timeout=15000)
            assert "already a section" in page.inner_text("#rp-verdict").lower()
        finally:
            browser.close()
    assert errors == [], errors


def test_the_page_names_no_transfer_and_renders_a_refusal_verbatim(live):
    """Two properties that are easy to lose and invisible in a screenshot.

    The page must not describe registering tissue as sending it anywhere — the
    bytes are read in place and a section is kept by path, so a transfer
    metaphor would be a lie about the product and a multi-gigabyte one if
    anybody acted on it. And a confinement refusal must arrive on screen in the
    SERVER's words: that 403 names the path sent, says the check was on the
    resolved object, names the variable and says how to add a root. Rewording
    it in the browser is how a UI comes to disagree with its own backend.
    """
    from playwright.sync_api import sync_playwright

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        try:
            page.goto(live["base"] + "/tissue", wait_until="domcontentloaded")
            page.locator(".tb-place").first.wait_for(timeout=15000)
            body = page.inner_text("body").lower()
            for word in ("upload", "drag and drop", "choose file"):
                assert word not in body, word

            # Type a path outside the allowlist: the server refuses and its
            # sentence is what appears.
            page.click("#tb-typed-wrap > summary")
            page.fill("#tb-typed", "/etc")
            page.click("#tb-typed-go")
            page.wait_for_selector(".tb-refusal", timeout=15000)
            shown = page.inner_text(".tb-refusal")
            assert "/etc" in shown
            assert cache_mod.DATA_ROOTS_VAR in shown
            assert "resolved" in shown
            # ...and it still does not disclose the allowlist.
            assert live["root"] not in shown
        finally:
            browser.close()
    assert errors == [], errors
