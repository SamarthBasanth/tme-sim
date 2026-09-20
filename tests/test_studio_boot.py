"""The Studio has to start when a catalogue call does not.

MEASURED ON THE LIVE DEPLOYMENT 2026-09-20. `GET /api/drugs` did not answer
inside 420 s while `/healthz` answered the same container in 0.2 s and the
same image answered in 0.012 s on a laptop. Building that catalogue runs 80
Brownian diffusivity simulations (5.2 s of 6.5 s locally, minutes on a
throttled host), and `boot()` awaited it in the same `Promise.all` as two
calls that return immediately.

So the Studio rendered `SECTION — RUN — CELLS —` and booted no further: no
project context, no runs, no deep link, no route back to the workspace. One
slow endpoint and the product did not start. It also made contract sec. 8 look
unimplemented, because `POST /api/projects/{id}/open` is called after that
await and never ran.

This file serves the REAL `app/web` in a REAL browser and holds `/api/drugs`
open forever. Nothing is stubbed about the page itself.
"""
from __future__ import annotations

import functools
import http.server
import json
import os
import socketserver
import threading

import pytest

pytest.importorskip("playwright.sync_api")
from playwright.sync_api import sync_playwright  # noqa: E402

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(_ROOT, "app", "web")


@pytest.fixture(scope="module")
def web_server():
    handler = functools.partial(http.server.SimpleHTTPRequestHandler,
                                directory=WEB)

    class Quiet(socketserver.TCPServer):
        allow_reuse_address = True

        def handle_error(self, request, client_address):  # noqa: D102
            pass

    with Quiet(("127.0.0.1", 0), handler) as httpd:
        port = httpd.server_address[1]
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        try:
            yield f"http://127.0.0.1:{port}"
        finally:
            httpd.shutdown()


_SECTIONS = [{
    "id": "xenium_319", "label": "Xenium GBM · 319-panel", "n_cells": 40887,
    "extent_um": [8087.8, 5065.4], "panel": "hBrain_254g", "available": True,
    "source": "xenium", "sample": None, "cache_dir": "data_cache/xenium_gbm",
    "assay": "Xenium In Situ", "accession": "10x", "has_outcome_label": False,
    "reason": None, "remedy": [],
}]
_DEFAULTS = {
    "section": "xenium_319", "dx_um": 21.0, "drug": "temozolomide",
    "sensitivity_mode": "uniform", "vascular_mode": "discrete",
    "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": None,
    "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0, "n_doses": 1,
           "dose_interval_s": 86400.0},
    "total_time_s": 21600.0,
}


def _open_studio(server, *, stall_drugs: bool):
    """Open `/index.html` with the catalogue answered, optionally stalling
    `/api/drugs` the way the live host did."""
    p = sync_playwright().start()
    browser = p.chromium.launch()
    pg = browser.new_page(viewport={"width": 1440, "height": 900})
    errors = []
    pg.on("pageerror", lambda e: errors.append(str(e)))

    def route_api(route):
        url = route.request.url
        if "/api/drugs" in url:
            if stall_drugs:
                return          # held open forever — never fulfilled
            return route.fulfill(status=200, content_type="application/json",
                                 body=json.dumps([]))
        if "/api/sections" in url:
            return route.fulfill(status=200, content_type="application/json",
                                 body=json.dumps(_SECTIONS))
        if "/api/defaults" in url:
            return route.fulfill(status=200, content_type="application/json",
                                 body=json.dumps(_DEFAULTS))
        return route.fulfill(status=200, content_type="application/json",
                             body="[]")

    pg.route("**/api/**", route_api)
    pg.goto(f"{server}/index.html", wait_until="domcontentloaded")
    pg.wait_for_timeout(2500)
    return p, browser, pg, errors


class TestBootSurvivesASlowCatalogueCall:
    def test_the_section_header_fills_when_drugs_never_answers(
            self, web_server):
        """THE LIVE FAILURE, reproduced. The header read em-dashes because
        boot() was still awaiting the compound list."""
        p, browser, pg, errors = _open_studio(web_server, stall_drugs=True)
        try:
            body = pg.inner_text("body")
            assert "xenium_319" in body, (
                "the Studio did not get as far as painting the section header "
                "while /api/drugs was unanswered — boot is still waiting on "
                "the compound list.\n" + body[:400])
            assert not errors, errors[:2]
        finally:
            browser.close()
            p.stop()

    def test_the_route_back_exists_when_drugs_never_answers(self, web_server):
        """The exit link is mounted after the catalogue await. Blocking there
        left a user on a page with no way out AND no content."""
        p, browser, pg, errors = _open_studio(web_server, stall_drugs=True)
        try:
            hrefs = pg.eval_on_selector_all(
                "a[href]",
                "els => els.filter(e => e.offsetParent !== null)"
                "        .map(e => e.getAttribute('href'))")
            assert "/" in hrefs, (
                "no visible link back to the workspace while /api/drugs was "
                f"unanswered; visible links were {hrefs}")
        finally:
            browser.close()
            p.stop()


class TestTheRouteBackIsNotConditionalOnTheUrl:
    def test_studio_without_a_project_still_offers_the_workspace(
            self, web_server):
        """Until 2026-09-20 `← workspace` was mounted ONLY by the project
        context, so `/studio` with no `?project=` had no exit at all — and the
        user most likely to want the workspace is the one who arrived without
        going through it."""
        p, browser, pg, errors = _open_studio(web_server, stall_drugs=False)
        try:
            links = pg.eval_on_selector_all(
                "a[href='/']",
                "els => els.filter(e => e.offsetParent !== null)"
                "        .map(e => (e.innerText || '').trim())")
            assert links, "no visible link to the workspace home at /studio"
            assert any("workspace" in t.lower() for t in links), links
            assert not errors, errors[:2]
        finally:
            browser.close()
            p.stop()


class TestRestoreDoesNotDuplicateARun:
    """`_restore` reads each run's id from INSIDE its result.json, not from the
    directory name, so two directories can declare the same run_id — a copied
    run directory, a restored backup, a hand-edited document.

    `_jobs` is a dict and de-duplicated silently. `_order` was a list and did
    not, and `list()` maps `_order` back through `_jobs`, so it returned THE
    SAME Job OBJECT TWICE.

    Found 2026-09-20 while checking whether cross_section.py's `n_runs` could
    over-count the way workspace.py's did. It could not — but a list-based
    accumulator downstream reported 2 runs for a compound holding 1, so every
    consumer was quietly compensating for this.
    """

    def test_two_directories_with_one_id_yield_one_job(self, tmp_path):
        from app.server.jobs import JobManager

        root = tmp_path / "runs"
        root.mkdir()
        for d in ("r_00001", "r_00002"):
            rd = root / d
            rd.mkdir()
            (rd / "result.json").write_text(json.dumps({
                # BOTH declare r_00001 — the collision this guards.
                "run_id": "r_00001",
                "status": "done",
                "spec": {"section": "xenium_319", "drug": "temozolomide"},
            }))

        mgr = JobManager(str(root))
        listed = mgr.list()
        ids = [j.run_id for j in listed]
        assert ids.count("r_00001") == 1, (
            "JobManager.list() returned the same run more than once: "
            f"{ids}. Every downstream count that trusts this list is then "
            "wrong, and only consumers that happen to de-duplicate are safe.")
        assert len({id(j) for j in listed}) == len(listed), (
            "the same Job object appears twice in list()")

    def test_an_ordinary_pair_of_runs_is_untouched(self, tmp_path):
        """The de-duplication must not lose a genuinely distinct run."""
        from app.server.jobs import JobManager

        root = tmp_path / "runs"
        root.mkdir()
        for d in ("r_00001", "r_00002"):
            rd = root / d
            rd.mkdir()
            (rd / "result.json").write_text(json.dumps({
                "run_id": d, "status": "done",
                "spec": {"section": "xenium_319", "drug": "temozolomide"},
            }))

        ids = sorted(j.run_id for j in JobManager(str(root)).list())
        assert ids == ["r_00001", "r_00002"], ids
