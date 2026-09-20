"""The project store's free-text fields are BOUNDED, and ``/api/workspace``
is bounded because of it.

WHAT THIS EXISTS BECAUSE OF
---------------------------
``ProjectStore.create``/``update``/``duplicate``/``add_specimen``/
``update_specimen`` checked their free text for *emptiness* and nothing else.
A name could be empty-rejected at one character and accepted at twenty
million. Measured against a live server on 2026-09-20:

    POST /api/projects   {"name": "pwn", "notes": "N" * 20_971_520}   -> 201
    GET  /api/workspace                                    163 B -> 20,972,002 B

``_summary`` puts ``notes`` in every row of ``store.list()``, and
``GET /api/workspace`` is the SINGLE boot call for the home page. So one
accepted request makes the front door of the product permanently unusable for
everyone, and it survives a restart because it is on disk in ``project.json``.

The 4 MB request-body cap did not stop it: ``api.py::_cap_request_body`` reads
``Content-Length``, and a ``Transfer-Encoding: chunked`` request does not have
one. That is a second defect, in a file this test does not own. This test pins
the store's own half, which holds whether or not the body cap is ever fixed:
**the store decides what it is willing to persist.**

WHY A CHARACTER CAP AND NOT A BYTE CAP
--------------------------------------
The field is text and the limit has to be explicable to the person typing into
it ("a name is 200 characters"), not to whoever is counting UTF-8 bytes
underneath. ``len()`` on a ``str`` is the number the UI can also enforce.
"""
from __future__ import annotations

import pytest

from app.server.projects import (
    MAX_NAME_CHARS,
    MAX_NOTES_CHARS,
    ProjectError,
    ProjectStore,
)


@pytest.fixture
def store(tmp_path):
    return ProjectStore(root=str(tmp_path / "projects"),
                        runs_root=str(tmp_path / "runs"))


# --------------------------------------------------------------------------- #
# The limits are real numbers, not None
# --------------------------------------------------------------------------- #
def test_the_limits_are_finite_and_sane():
    """A cap of None, 0 or 10**9 would make every test below vacuous."""
    assert isinstance(MAX_NAME_CHARS, int) and 0 < MAX_NAME_CHARS <= 10_000
    assert isinstance(MAX_NOTES_CHARS, int) and 0 < MAX_NOTES_CHARS <= 1_000_000


# --------------------------------------------------------------------------- #
# create
# --------------------------------------------------------------------------- #
def test_create_refuses_an_oversize_name(store):
    with pytest.raises(ProjectError) as ei:
        store.create("x" * (MAX_NAME_CHARS + 1))
    assert ei.value.status == 400
    # The refusal must name the field and the limit, or it cannot be acted on.
    assert "name" in ei.value.message
    assert str(MAX_NAME_CHARS) in ei.value.message


def test_create_refuses_oversize_notes(store):
    with pytest.raises(ProjectError) as ei:
        store.create("fine", notes="N" * (MAX_NOTES_CHARS + 1))
    assert ei.value.status == 400
    assert "notes" in ei.value.message
    assert str(MAX_NOTES_CHARS) in ei.value.message


def test_create_accepts_exactly_the_limit(store):
    """The cap must not be off by one against the user who fills the field."""
    doc = store.create("x" * MAX_NAME_CHARS, notes="N" * MAX_NOTES_CHARS)
    assert len(doc["name"]) == MAX_NAME_CHARS
    assert len(doc["notes"]) == MAX_NOTES_CHARS


# --------------------------------------------------------------------------- #
# update / duplicate
# --------------------------------------------------------------------------- #
def test_update_refuses_an_oversize_name(store):
    p = store.create("ok")
    with pytest.raises(ProjectError):
        store.update(p["id"], name="x" * (MAX_NAME_CHARS + 1))
    # and the stored document is untouched
    assert store.get(p["id"], refresh_capability=False)["name"] == "ok"


def test_update_refuses_oversize_notes(store):
    p = store.create("ok", notes="short")
    with pytest.raises(ProjectError):
        store.update(p["id"], notes="N" * (MAX_NOTES_CHARS + 1))
    assert store.get(p["id"], refresh_capability=False)["notes"] == "short"


def test_duplicate_refuses_an_oversize_name(store):
    p = store.create("ok")
    with pytest.raises(ProjectError):
        store.duplicate(p["id"], name="x" * (MAX_NAME_CHARS + 1))


def test_duplicate_of_a_long_named_project_stays_under_the_cap(store):
    """``duplicate`` appends " (copy)", so a name AT the cap would exceed it.

    Duplicating is not a text-entry surface -- the user supplied no name --
    so this must not be a refusal the user cannot act on. It must fit.
    """
    p = store.create("x" * MAX_NAME_CHARS)
    dup = store.duplicate(p["id"])
    assert len(dup["name"]) <= MAX_NAME_CHARS


# --------------------------------------------------------------------------- #
# specimens -- the same free text, one level down
# --------------------------------------------------------------------------- #
def test_update_specimen_refuses_oversize_label_and_notes(store, monkeypatch):
    import app.server.projects as P

    monkeypatch.setattr(P, "all_sections", lambda: {"sec_a"})
    monkeypatch.setattr(P, "capability_report", lambda sid: {"label": "A"})
    p = store.create("ok")
    spec = store.add_specimen(p["id"], "sec_a")

    with pytest.raises(ProjectError):
        store.update_specimen(p["id"], spec["id"],
                              label="x" * (MAX_NAME_CHARS + 1))
    with pytest.raises(ProjectError):
        store.update_specimen(p["id"], spec["id"],
                              notes="N" * (MAX_NOTES_CHARS + 1))


def test_add_specimen_refuses_oversize_label_and_notes(store, monkeypatch):
    import app.server.projects as P

    monkeypatch.setattr(P, "all_sections", lambda: {"sec_a"})
    monkeypatch.setattr(P, "capability_report", lambda sid: {"label": "A"})
    p = store.create("ok")

    with pytest.raises(ProjectError):
        store.add_specimen(p["id"], "sec_a", label="x" * (MAX_NAME_CHARS + 1))
    with pytest.raises(ProjectError):
        store.add_specimen(p["id"], "sec_a", notes="N" * (MAX_NOTES_CHARS + 1))


# --------------------------------------------------------------------------- #
# THE PROPERTY THE CAPS EXIST FOR
# --------------------------------------------------------------------------- #
def test_one_project_cannot_make_the_workspace_boot_call_unbounded(store):
    """``store.list()`` is what ``GET /api/workspace`` returns.

    This is the end-to-end statement, written against the store rather than
    the HTTP route so it holds no matter which router calls it. Without the
    caps this row is ~20 MB; the assertion is deliberately loose (a generous
    multiple of the text caps) so it pins the ORDER OF MAGNITUDE and does not
    fail the next time a field is added to ``_summary``.
    """
    import json

    for i in range(3):
        try:
            store.create(f"p{i}", notes="N" * (MAX_NOTES_CHARS + 10_000_000))
        except ProjectError:
            store.create(f"p{i}", notes="N" * MAX_NOTES_CHARS)

    rows = store.list(thumbnails=False)
    assert len(rows) == 3
    budget = 3 * (MAX_NAME_CHARS + MAX_NOTES_CHARS) + 8_192
    assert len(json.dumps(rows)) <= budget, (
        "GET /api/workspace grows without bound with what a caller stored: "
        f"{len(json.dumps(rows))} B for 3 projects, budget {budget} B"
    )
