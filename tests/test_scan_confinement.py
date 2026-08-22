"""``POST /api/sections/scan`` confined to allowlisted roots.

The route takes a filesystem path in the request body and inspects whatever
directory it names. The API key AUTHENTICATES that and does not CONFINE it, so
until ``TMESIM_STUDIO_DATA_ROOTS`` existed the key was equivalent to read
access on the host as the server's uid — which ``docs/DEPLOYMENT.md`` §2 has
always said out loud. That trade was fine for one operator on a laptop and is
not fine for an external field tester.

**Every escape here is a real one, attempted over HTTP.** Nothing below
asserts that a helper returned ``False``: each test first shows what the
unconfined route discloses about the target directory, then shows the same
request refused with the roots configured. A guard that was accidentally a
no-op would pass a "helper returns False" test and fail these.

The four escape shapes covered: ``..`` traversal out of a root, a symlink
inside a root pointing out of it, a ``..`` that traverses THROUGH a symlink
(where a textual ``abspath`` check and the real object disagree), and a
sibling directory whose path merely starts with the root's string. Plus the
two behaviours that must survive the guard: unset stays unconfined, and a
legitimate bring-your-own-tissue ingest under a root still works end to end.
"""
from __future__ import annotations

import os

import pytest

from app.server import cache as cache_mod
from tests.test_ingest import _write_synthetic_xenium, _SYNTH_MARKER_PANEL


@pytest.fixture
def client(tmp_path, monkeypatch):
    """The app, with run history isolated from the repo and no roots set."""
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
    app = create_app(JobManager(root=str(runs / "jobs")))
    with TestClient(app) as c:
        yield c
    monkeypatch.setattr(ingest_registry, "_store", None)


@pytest.fixture
def estate(tmp_path):
    """A root to allow, a sample inside it, and a sample OUTSIDE it.

    ``outside/`` is the thing a tester holding the key must not be able to
    enumerate. It holds a sample directory whose NAME is the disclosure — a
    specimen identifier is exactly the kind of thing §0 warns travels in a
    directory name.
    """
    allowed = tmp_path / "allowed"
    (allowed).mkdir()
    _write_synthetic_xenium(str(allowed / "Study_A_S001"), _SYNTH_MARKER_PANEL)

    outside = tmp_path / "outside"
    outside.mkdir()
    _write_synthetic_xenium(str(outside / "AT10_BRA_5_FO_1"), _SYNTH_MARKER_PANEL)

    # A sibling whose absolute path starts with the root's string. A prefix
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


def _scan(client, path, **body):
    return client.post("/api/sections/scan", json=dict(path=str(path), **body))


def _names(resp):
    """The child sample names an unconfined scan discloses for a container."""
    return [c["name"] for c in resp.json()["detection"]["contained_samples"]]


# --------------------------------------------------------------------------- #
# The default: unset means unconfined, and it says so
# --------------------------------------------------------------------------- #
def test_unset_leaves_scanning_unconfined_exactly_as_before(client, estate):
    """Backwards compatibility is the whole reason unset is not confined.

    Bring-your-own-tissue points the scanner at wherever the operator's export
    lives, which is almost never inside this repo. With the variable unset the
    route reads an arbitrary directory, as it always has.
    """
    resp = _scan(client, estate["outside"])
    assert resp.status_code == 200
    assert _names(resp) == ["AT10_BRA_5_FO_1"]


def test_the_first_unconfined_scan_warns_once_naming_the_variable(
        client, estate, monkeypatch, capsys):
    """The unconfined default is deliberate; a SILENT unconfined default is not.

    This repo's precedent is the unauthenticated loopback bind: keep the
    single-user workflow working, and make the risk loud (``__main__.py``'s
    ``_no_auth_warning``). Same shape here, at the moment the capability is
    actually used.
    """
    monkeypatch.setattr(cache_mod, "_UNCONFINED_WARNED", False)
    capsys.readouterr()

    assert _scan(client, estate["outside"]).status_code == 200
    first = capsys.readouterr().err
    assert "UNCONFINED" in first
    assert cache_mod.DATA_ROOTS_VAR in first
    assert estate["outside"] in first

    # Once per process, not once per scan: this must not become log spam that
    # an operator learns to filter out.
    assert _scan(client, estate["outside"]).status_code == 200
    assert cache_mod.DATA_ROOTS_VAR not in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# Real escapes, each shown working first
# --------------------------------------------------------------------------- #
def test_dotdot_traversal_out_of_a_root_is_refused(client, estate, monkeypatch):
    """``<root>/../outside`` — the classic, and the one a string prefix misses."""
    escape = os.path.join(estate["root"], "..", "outside")

    # BEFORE: unconfined, the traversal works and enumerates the directory.
    before = _scan(client, escape)
    assert before.status_code == 200
    assert _names(before) == ["AT10_BRA_5_FO_1"]

    # AFTER: with the root allowlisted, the same request is refused...
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    after = _scan(client, escape)
    assert after.status_code == 403
    # ...and discloses nothing about what is out there.
    assert "AT10_BRA_5_FO_1" not in after.text
    assert "contained_samples" not in after.text


def test_a_symlink_out_of_a_root_is_refused(client, estate, monkeypatch, tmp_path):
    """A link INSIDE the allowlist pointing out of it.

    The path text is entirely inside the root, so only resolving the object
    catches this.
    """
    link = os.path.join(estate["root"], "shortcut")
    os.symlink(estate["outside"], link)

    before = _scan(client, link)
    assert before.status_code == 200
    assert _names(before) == ["AT10_BRA_5_FO_1"]

    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    after = _scan(client, link)
    assert after.status_code == 403
    assert "AT10_BRA_5_FO_1" not in after.text


def test_dotdot_through_a_symlink_is_judged_on_the_real_object(
        client, estate, monkeypatch):
    """``<root>/link/../…`` — where a TEXTUAL reading and reality disagree.

    ``os.path.abspath`` collapses ``..`` textually, so for this path it reports
    a directory inside the root while the object the path really designates is
    outside it. A check built on ``abspath`` — not just on the raw string —
    still gets this wrong. The guard resolves instead, so it judges the object,
    and the route is then handed that resolved path so the directory checked
    and the directory opened cannot be two different things.
    """
    os.symlink(estate["outside"], os.path.join(estate["root"], "hop"))
    escape = os.path.join(estate["root"], "hop", "..", "outside",
                          "AT10_BRA_5_FO_1")

    # The textual reading of that path is inside the root...
    assert os.path.abspath(escape).startswith(estate["root"] + os.sep)
    # ...and the real one is not.
    assert not os.path.realpath(escape).startswith(estate["root"] + os.sep)

    # Unconfined the route answers about SOME directory (which one depends on
    # whose reading you take — that is the whole problem).
    assert _scan(client, escape).status_code == 200

    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    assert _scan(client, escape).status_code == 403

    # And the same divergence in the other direction: a path that reads as
    # OUTSIDE textually but really designates a sample inside the root is
    # allowed, and what gets scanned is the real sample rather than the text.
    os.symlink(estate["sample"], os.path.join(estate["outside"], "back"))
    inbound = os.path.join(estate["root"], "..", "outside", "back")
    assert not os.path.abspath(inbound).startswith(estate["root"] + os.sep)
    ok = _scan(client, inbound)
    assert ok.status_code == 200
    assert ok.json()["detection"]["path"] == os.path.realpath(estate["sample"])


def test_the_root_itself_is_allowed_and_a_prefix_sibling_is_not(
        client, estate, monkeypatch):
    """The exact-root boundary, in both directions.

    ``<tmp>/allowed`` must contain itself and must NOT contain
    ``<tmp>/allowed_evil``, whose absolute path starts with the same string.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])

    at_the_root = _scan(client, estate["root"])
    assert at_the_root.status_code == 200
    assert _names(at_the_root) == ["Study_A_S001"]

    assert _scan(client, estate["sample"]).status_code == 200

    sibling = _scan(client, estate["sibling"])
    assert sibling.status_code == 403
    assert "Not_Yours" not in sibling.text
    # And the parent of the root is not the root.
    assert _scan(client, estate["tmp"]).status_code == 403


def test_the_registering_branch_is_confined_too(client, estate, monkeypatch):
    """``register: true`` takes the same ``path`` and must not be a way round."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    resp = _scan(client, os.path.join(estate["outside"], "AT10_BRA_5_FO_1"),
                 register=True, section_id="stolen")
    assert resp.status_code == 403
    assert "stolen" not in [r["id"] for r in client.get("/api/sections").json()]


# --------------------------------------------------------------------------- #
# The refusal is a diagnostic, and it is not an information leak
# --------------------------------------------------------------------------- #
def test_the_refusal_names_the_path_the_variable_and_the_way_out(
        client, estate, monkeypatch):
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    detail = _scan(client, estate["outside"]).json()["detail"]

    assert detail["code"] == "path_outside_data_roots"
    assert detail["path"] == estate["outside"]
    msg = detail["message"]
    assert estate["outside"] in msg                  # names the path
    # names what went wrong -- and not merely because the fixture directory
    # happens to be called "outside".
    assert "is outside the directories this server allows" in msg
    assert "resolved" in msg                         # says the check is on the object
    assert cache_mod.DATA_ROOTS_VAR in msg           # names the control
    assert "add its root" in msg                     # names the way out
    assert "restart" in msg


def test_the_resolved_path_reaches_the_console_and_not_the_client(estate,
                                                                  monkeypatch):
    """``~/x`` must not come back as the operator's home directory.

    The refusal has to name a path or it is not a diagnostic; it names the
    caller's OWN text. The resolved path is on the exception for the server
    side and is not in the body.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    monkeypatch.setenv("HOME", estate["outside"])
    with pytest.raises(cache_mod.PathNotAllowed) as exc:
        cache_mod.confine_scan_path("~/AT10_BRA_5_FO_1")
    assert exc.value.path == "~/AT10_BRA_5_FO_1"
    assert exc.value.resolved == os.path.realpath(
        os.path.join(estate["outside"], "AT10_BRA_5_FO_1"))
    assert exc.value.resolved not in str(exc.value.as_detail())


def test_the_refusal_never_discloses_the_allowed_roots(client, estate,
                                                       monkeypatch):
    """A caller learns that it was refused, not what would be accepted."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR,
                       os.pathsep.join([estate["root"], estate["sibling"]]))
    body = _scan(client, estate["outside"]).text
    assert estate["root"] not in body
    assert estate["sibling"] not in body


def test_a_set_but_unusable_value_fails_closed(client, estate, monkeypatch):
    """``-e TMESIM_STUDIO_DATA_ROOTS=`` must not mean "unconfined".

    An operator who set this variable asked for confinement. A value that
    parses to no root confines scanning to nothing — including the sample that
    would otherwise be legal — and says why.
    """
    for value in ("", "   ", os.pathsep * 3):
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, value)
        resp = _scan(client, estate["sample"])
        assert resp.status_code == 403, value
        assert cache_mod.DATA_ROOTS_VAR in resp.json()["detail"]["message"]


# --------------------------------------------------------------------------- #
# Bring your own tissue still works — the product's main input
# --------------------------------------------------------------------------- #
def test_a_legitimate_ingest_under_an_allowed_root_still_works(
        client, estate, monkeypatch):
    """Scan, register, and appear in the section list, with the guard on."""
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])

    scan = _scan(client, estate["sample"])
    assert scan.status_code == 200, scan.text
    assert scan.json()["detection"]["platform"] == "xenium"
    assert scan.json()["capability"] is not None

    reg = _scan(client, estate["sample"], register=True,
                section_id="study_a_s001", label="Study A S001")
    assert reg.status_code == 200, reg.text
    body = reg.json()
    assert body["ok"], body.get("error")

    row = next(r for r in client.get("/api/sections").json()
               if r["id"] == "study_a_s001")
    assert row["available"] and row["ingested"]
    assert row["cache_dir"] == os.path.realpath(estate["sample"])
    client.delete("/api/sections/study_a_s001")


def test_a_relative_path_and_a_tilde_are_resolved_before_the_check(
        client, estate, monkeypatch):
    """Both forms the documented workflow uses go through the same resolution.

    ``docs/DEPLOYMENT.md`` §1d's own example passes a path relative to the
    server's working directory; ``~`` is what an operator types by hand.
    Neither may sidestep the allowlist, and a relative path that lands inside
    a root must still work.
    """
    monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, estate["root"])
    monkeypatch.chdir(estate["root"])
    assert _scan(client, "Study_A_S001").status_code == 200
    assert _scan(client, os.path.join("..", "outside")).status_code == 403

    monkeypatch.setenv("HOME", estate["outside"])
    assert _scan(client, "~/AT10_BRA_5_FO_1").status_code == 403


def test_the_roots_variable_itself_resolves_symlinks_and_the_tilde(estate,
                                                                   monkeypatch):
    """A root given as ``~/x`` or through a symlink still matches real paths.

    Otherwise an operator who configured the allowlist the way they type paths
    would find every legitimate scan refused.
    """
    link = os.path.join(estate["tmp"], "root_link")
    os.symlink(estate["root"], link)
    monkeypatch.setenv("HOME", estate["tmp"])

    assert cache_mod.data_roots(link) == [os.path.realpath(estate["root"])]
    assert cache_mod.data_roots("~/allowed") == [os.path.realpath(estate["root"])]
    # Unset is the one value that means "no confinement at all".
    monkeypatch.delenv(cache_mod.DATA_ROOTS_VAR, raising=False)
    assert cache_mod.data_roots() is None
