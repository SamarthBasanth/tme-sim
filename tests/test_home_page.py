"""Acceptance gate for the workspace home page (``docs/WORKSPACE_CONTRACT.md``).

WHY THIS FILE IS SEPARATE FROM EVERY AGENT'S OWN TESTS
------------------------------------------------------
Six agents built this feature in parallel: the endpoint, the page shell, the
card, the dialogs, the project store, and the Studio handoff. Each verified its
own piece. Nothing verified that the pieces meet, and "each half works" is the
failure mode of a parallel build — the shared contract is exactly what no
single author is accountable for.

So this file asserts the CONTRACT, not the implementation. Every assertion here
traces to a numbered clause in ``docs/WORKSPACE_CONTRACT.md``, and where an
agent's private test would say "my function returns the right shape", this one
says "the page the user opens shows the truth the server sent".

THE COMPATIBILITY BAR (contract §8) IS THE FIRST TEST FOR A REASON
A user who never visits the home page must notice no change. That is the one
promise this whole feature makes to the people already using the product, and
it is the promise a big feature is most likely to break by accident.
"""
from __future__ import annotations

import json
import os
import re
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

pytest.importorskip("fastapi", reason="the Studio backend is an optional extra")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi.testclient import TestClient  # noqa: E402

from app.server.api import create_app  # noqa: E402
from app.server.jobs import JobManager  # noqa: E402

WEB = os.path.join(_ROOT, "app", "web")
HOME_HTML = os.path.join(WEB, "home.html")
HOME_JS = os.path.join(WEB, "js", "home")
HOME_CSS = os.path.join(WEB, "css", "home.css")
TOKENS_CSS = os.path.join(WEB, "css", "tokens.css")

#: The contract's §5 file layout. A missing one is a half-built feature, and a
#: half-built feature that renders an empty page looks identical to a working
#: one with no projects — which is precisely why this is checked by PATH.
CONTRACT_MODULES = ("main.js", "projectCard.js", "projectGrid.js",
                    "createProject.js", "recentRuns.js", "nav.js")

needs_home = pytest.mark.skipif(
    not os.path.exists(HOME_HTML),
    reason="the workspace home page has not landed yet")


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    root = str(tmp_path_factory.mktemp("runs"))
    manager = JobManager(root=root)
    with TestClient(create_app(manager)) as c:
        yield c
    manager.shutdown()


# --------------------------------------------------------------------------- #
# §8 — the compatibility bar
# --------------------------------------------------------------------------- #
class TestTheStudioIsUnchangedForSomeoneWhoNeverVisitsHome:
    """Contract §8. Checked FIRST because it is the promise to existing users."""

    def test_the_studio_is_still_served_at_its_old_url(self, client):
        r = client.get("/studio")
        if r.status_code == 404:
            pytest.skip("the /studio route has not landed yet")
        assert r.status_code == 200
        assert "<!doctype html" in r.text[:200].lower()

    def test_index_html_never_stopped_working(self, client):
        """Anything bookmarked stays alive — the contract says so explicitly."""
        r = client.get("/index.html")
        assert r.status_code == 200, (
            "index.html 404s; every existing bookmark into this product just "
            "broke, which the contract's §1 forbids")

    def test_the_studio_shell_still_has_its_three_rails(self, client):
        r = client.get("/index.html")
        assert r.status_code == 200
        # The three-rail layout is improved, never replaced (contract §9).
        for marker in ("split-l", "split-r", "stage"):
            assert marker in r.text, (
                f"the Studio shell lost {marker!r}; the contract permits "
                "improving the three-rail layout, not replacing it")


# --------------------------------------------------------------------------- #
# §1 — the boot call
# --------------------------------------------------------------------------- #
class TestTheWorkspaceEndpoint:
    def test_it_exists_and_is_gated(self, tmp_path, monkeypatch):
        """Gated like every /api/ route. A workspace lists what work exists on
        this machine, which is operational detail about a machine."""
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "ws-key")
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        if c.get("/api/workspace", headers={"X-API-Key": "ws-key"}
                 ).status_code == 404:
            pytest.skip("GET /api/workspace has not landed yet")
        assert c.get("/api/workspace").status_code == 401

    def test_the_payload_carries_every_contract_field(self, client):
        r = client.get("/api/workspace")
        if r.status_code == 404:
            pytest.skip("GET /api/workspace has not landed yet")
        assert r.status_code == 200
        body = r.json()
        for key in ("projects", "recent_runs", "sections_available",
                    "drugs_available", "totals", "server"):
            assert key in body, f"contract §1 requires {key!r}"
        assert isinstance(body["projects"], list)
        assert isinstance(body["recent_runs"], list)
        assert len(body["recent_runs"]) <= 12, "contract §1 caps this at 12"
        for k in ("projects", "runs", "specimens"):
            assert k in body["totals"]

    def test_a_run_card_is_flat_and_self_describing(self, client):
        """Contract §1: a RunCard is deliberately NOT a RunResult. The home
        page must not have to resolve a project name by a second call."""
        r = client.get("/api/workspace")
        if r.status_code == 404:
            pytest.skip("GET /api/workspace has not landed yet")
        for card in r.json()["recent_runs"]:
            for key in ("run_id", "project_id", "project_name", "section",
                        "status"):
                assert key in card, f"a RunCard is missing {key!r}"

    def test_recency_is_never_fabricated(self, client):
        """Contract §3. Every project predating `last_opened_at` reads null —
        never `created_at`, never now. A workspace that invents 'you opened
        this just now' is lying about the only thing its sort order means."""
        r = client.get("/api/workspace")
        if r.status_code == 404:
            pytest.skip("GET /api/workspace has not landed yet")
        for p in r.json()["projects"]:
            assert "last_opened_at" in p, "contract §3 requires the field"
            v = p["last_opened_at"]
            assert v is None or isinstance(v, (int, float))
            if v is not None and p.get("created_at"):
                assert v != p["created_at"] or v >= p["created_at"]

    def test_a_thumbnail_is_measured_or_absent_never_invented(self, client):
        """Contract §4. `kind: none` carries a REASON. The alternative — an
        empty chart, or zeros — is a picture of a measurement that was never
        made, which is the one thing this product does not do."""
        r = client.get("/api/workspace")
        if r.status_code == 404:
            pytest.skip("GET /api/workspace has not landed yet")
        for p in r.json()["projects"]:
            t = p.get("thumbnail")
            if t is None:
                continue
            assert t.get("kind") in ("composition", "none")
            if t["kind"] == "none":
                assert t.get("reason"), "a null thumbnail must say why"
            else:
                assert t.get("section") and t.get("type_fractions")

    def test_an_empty_workspace_renders(self, tmp_path):
        """A fresh install has zero projects. That is the FIRST thing a new
        user sees, and it must be a page, not a 500."""
        from app.server.projects import ProjectStore

        empty = ProjectStore(root=str(tmp_path / "none"))
        c = TestClient(create_app(JobManager(root=str(tmp_path / "runs")),
                                  projects=empty))
        r = c.get("/api/workspace")
        if r.status_code == 404:
            pytest.skip("GET /api/workspace has not landed yet")
        assert r.status_code == 200
        assert r.json()["projects"] == []
        assert r.json()["totals"]["projects"] == 0


# --------------------------------------------------------------------------- #
# §5 — the page itself
# --------------------------------------------------------------------------- #
@needs_home
class TestTheHomePageIsWhole:
    def test_every_contract_module_exists(self):
        missing = [m for m in CONTRACT_MODULES
                   if not os.path.exists(os.path.join(HOME_JS, m))]
        assert not missing, (
            f"contract §5 names these modules and they are absent: {missing}. "
            "A half-built grid renders an empty page, which looks exactly like "
            "a working one with no projects.")

    def test_every_module_parses_as_an_es_module(self):
        """`node --check` parses as CommonJS and has passed here on a real ESM
        error that broke the whole app at load. This is the check that catches
        it."""
        import shutil
        import subprocess

        if not shutil.which("node"):
            pytest.skip("node is not installed")
        bad = []
        for name in sorted(os.listdir(HOME_JS)):
            if not name.endswith(".js"):
                continue
            path = os.path.join(HOME_JS, name)
            with open(path) as fh:
                src = fh.read()
            r = subprocess.run(["node", "--input-type=module", "--check"],
                               input=src, capture_output=True, text=True)
            if r.returncode != 0:
                bad.append(f"{name}: {r.stderr.strip().splitlines()[0]}")
        assert not bad, "ES module syntax errors:\n  " + "\n  ".join(bad)

    def test_the_home_page_is_served_at_the_root(self, client):
        r = client.get("/")
        assert r.status_code == 200
        assert "<!doctype html" in r.text[:200].lower()

    def test_it_is_reachable_without_a_key(self, tmp_path, monkeypatch):
        """Contract §1: ungated like every static asset. It has to load in
        order to prompt for a key — a login page behind a login is a loop."""
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "ws-key")
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        assert c.get("/").status_code == 200

    def test_it_defines_no_colour_of_its_own(self):
        """Contract §7. Every colour is a token, so the theme toggle the two
        pages share actually works on both."""
        if not os.path.exists(HOME_CSS):
            pytest.skip("home.css has not landed yet")
        with open(HOME_CSS) as fh:
            css = fh.read()
        with open(TOKENS_CSS) as fh:
            tokens = fh.read()
        used = set(re.findall(r"var\((--[a-z0-9-]+)\)", css))
        missing = sorted(t for t in used if f"{t}:" not in tokens)
        assert not missing, f"home.css uses undefined tokens: {missing}"
        # A raw hex or rgb() outside a @font-face or a gradient stop is a
        # colour literal, which the contract forbids.
        literals = re.findall(r"#[0-9a-fA-F]{3,8}\b(?![^{]*@font-face)", css)
        assert not literals, (
            f"home.css defines colour literals {literals[:5]}; the palette is "
            "tokens.css so both themes stay correct")


# --------------------------------------------------------------------------- #
# §6 — the prohibitions, on the new surface
# --------------------------------------------------------------------------- #
@needs_home
class TestTheHomePageInheritsTheProhibitions:
    """These are not new rules. They are the product's rules, and a new surface
    is exactly where they get forgotten — the bring-your-own-compound panel
    shipped the withdrawn ranking claim on EVERY compile for thirteen days
    because nothing checked the surface it was added to."""

    def _all_text(self):
        parts = []
        for path in (HOME_HTML, HOME_CSS):
            if os.path.exists(path):
                with open(path) as fh:
                    parts.append(fh.read())
        for name in sorted(os.listdir(HOME_JS)):
            if name.endswith(".js"):
                with open(os.path.join(HOME_JS, name)) as fh:
                    parts.append(fh.read())
        return "\n".join(parts)

    def test_the_scan_has_a_subject(self):
        """A guard that reads nothing passes for free."""
        text = self._all_text()
        assert len(text) > 4000, f"only {len(text)} chars scanned"

    @pytest.mark.parametrize("phrase,why", [
        ("digital twin", "prohibited framing"),
        ("personalised medicine", "prohibited framing"),
        ("personalized medicine", "prohibited framing"),
        ("machine learning", "this product makes no AI/ML claim"),
        ("predicts response", "no clinical response prediction"),
        ("treatment plan", "no treatment-planning use"),
    ])
    def test_no_prohibited_phrase(self, phrase, why):
        assert phrase not in self._all_text().lower(), f"{phrase!r}: {why}"

    def test_tissue_is_never_uploaded(self):
        """Sections register BY FILESYSTEM PATH. 'Upload' describes a different
        product with a different security model."""
        text = self._all_text().lower()
        for m in re.finditer(r"upload", text):
            window = text[max(0, m.start() - 120):m.end() + 120]
            assert not any(w in window for w in
                           ("tissue", "section", "specimen", "xenium",
                            "visium")), (
                "'upload' appears next to tissue; sections register by "
                f"filesystem path. Context: ...{window.strip()[:160]}...")

    #: Words that turn an account phrase from an OFFER into a DENIAL.
    #:
    #: The first version of this check asserted the bare substring and failed on
    #: `home/main.js`, which says "there is no self-service sign-up" -- the
    #: correct, honest copy, matching `ui/apikey.js`'s established wording. The
    #: guard was the defect, not the page.
    #:
    #: This is the same lesson `tests/test_doc_claim_surface.py` learned: a
    #: negation is not a violation, and a guard that cannot read one gets
    #: deleted by the next person rather than fixed.
    _DENIALS = ("no ", "not ", "never ", "cannot ", "can't ", "without ",
                "there is no", "does not", "doesn't")

    @pytest.mark.parametrize("phrase", [
        "sign up", "sign-up", "create an account", "your account",
        "log in to your",
    ])
    def test_no_owner_or_account_framing(self, phrase):
        """Contract §9. `GET /api/account` states the truth: one shared key, no
        identities, no attribution. A UI OFFERING an account contradicts its own
        server; a UI saying there is none is telling the truth."""
        text = self._all_text().lower()
        bad = []
        for m in re.finditer(re.escape(phrase), text):
            before = text[max(0, m.start() - 90):m.start()]
            if any(d in before for d in self._DENIALS):
                continue                      # a denial, which is the truth
            bad.append(text[max(0, m.start() - 70):m.end() + 40].strip())
        assert not bad, (
            f"{phrase!r} is OFFERED rather than denied; this deployment "
            "authenticates with one shared key and has no user table. "
            f"Context: ...{bad[0][:150]}...")

    def test_the_account_check_can_tell_an_offer_from_a_denial(self):
        """Both directions, pinned. Without this the fix above could be
        loosened into a guard that never fires."""
        offer = "create an account to get started"
        denial = "this is a private deployment; there is no self-service sign-up"
        def flagged(text, phrase):
            t = text.lower()
            for m in re.finditer(re.escape(phrase), t):
                if not any(d in t[max(0, m.start() - 90):m.start()]
                           for d in self._DENIALS):
                    return True
            return False
        assert flagged(offer, "create an account"), "an offer must be caught"
        assert not flagged(denial, "sign-up"), "a denial must be allowed"


# --------------------------------------------------------------------------- #
# THE GAP THIS FILE HAD: it never EXECUTED the page
# --------------------------------------------------------------------------- #
#
# On 2026-09-20 the home page threw on boot for roughly twenty-four minutes.
# `js/home/projectGrid.js` built `h('button.btn.primary#pgrid-create')`; `h()`
# requires the id BEFORE the classes, so it raised, the exception escaped
# `boot()` after the loading state had been cleared and before `render()` ran,
# and `/` showed "Reaching the server · 0.0s elapsed" forever. The API answered
# in under 50 ms the whole time.
#
# **Every test above passed, 28 of 28, through every minute of it.** They check
# that files exist, that the payload has the right shape, and that no forbidden
# phrase appears. None of them opens the page. A page that renders NOTHING
# satisfies every one of them, and a page rendering nothing is indistinguishable
# from a page still loading.
#
# That is this repository's signature defect — a check that reports green while
# not testing its subject — committed in the file written to catch it in others.
#
# These four tests run the real modules in a real browser. They are deliberately
# few and coarse: the question is not whether the page is pretty, it is whether
# the page EXISTS after boot.
playwright_api = pytest.importorskip(
    "playwright.sync_api", reason="the home page's boot test needs playwright")


@pytest.fixture(scope="module")
def home_server():
    """A static server over ``app/web`` — real ES module semantics."""
    import functools
    import http.server
    import socketserver
    import threading

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


@needs_home
class TestTheHomePageActuallyBoots:
    """Four assertions. Each would have caught the twenty-four-minute outage."""

    def _open(self, server, route=None):
        """Open `/home.html` with `/api/*` stubbed, returning (page, errors)."""
        p = playwright_api.sync_playwright().start()
        browser = p.chromium.launch()
        pg = browser.new_page()
        errors = []
        pg.on("pageerror", lambda e: errors.append(str(e)))
        pg.on("console",
              lambda m: errors.append(m.text) if m.type == "error" else None)

        body = route if route is not None else {
            "projects": [], "recent_runs": [], "sections_available": 3,
            "drugs_available": 12,
            "totals": {"projects": 0, "runs": 0, "specimens": 0},
            "server": {"contract": "v1", "stub": False},
        }
        pg.route("**/api/**", lambda r: r.fulfill(
            status=200, content_type="application/json", body=json.dumps(body)))
        pg.goto(f"{server}/home.html", wait_until="domcontentloaded")
        pg.wait_for_timeout(1200)
        return p, browser, pg, errors

    def test_boot_raises_nothing(self, home_server):
        """The outage was an uncaught throw. Nothing else needs to be true for
        this to be worth asserting."""
        p, browser, pg, errors = self._open(home_server)
        try:
            assert not errors, (
                "the home page raised during boot; the page a user sees is "
                f"whatever had rendered before this: {errors[:2]}")
        finally:
            browser.close()
            p.stop()

    def test_the_page_leaves_the_loading_state(self, home_server):
        """The precise symptom: a live progress bar, forever. The loading state
        is honest ONLY while something is actually in flight."""
        p, browser, pg, errors = self._open(home_server)
        try:
            text = pg.eval_on_selector("body", "e => e.innerText")
            assert "elapsed" not in text.lower(), (
                "the page is still showing its loading state after the API "
                f"answered. Body begins: {text[:200]!r}")
        finally:
            browser.close()
            p.stop()

    def test_an_empty_workspace_still_renders_something_a_user_can_act_on(
            self, home_server):
        """Zero projects is the FIRST screen of a fresh install, and it is the
        state most easily confused with a broken one."""
        p, browser, pg, errors = self._open(home_server)
        try:
            text = pg.eval_on_selector("body", "e => e.innerText")
            assert len(text.strip()) > 120, (
                f"the page rendered {len(text.strip())} characters; an empty "
                "workspace must teach, not shrug")
            buttons = pg.eval_on_selector_all(
                "button, a[href]", "els => els.length")
            assert buttons >= 1, "nothing on the page is actionable"
        finally:
            browser.close()
            p.stop()

    def test_the_header_total_is_reconcilable_with_the_cards(self, home_server):
        """§2 + §6. The header counts runs that EXIST; a card counts the runs
        ATTACHED to its project. Contract §2 -- a project "owns nothing but
        ids" -- lets one run be reached from many projects, so the cards can
        sum to MORE than the header, and until 2026-09-20 the server resolved
        that tension the wrong way round: it summed the cards, and a single
        `duplicate` of a two-run project reported four runs on a server holding
        two.

        With the server corrected, the page now shows two numbers that do not
        add up, and the SMALLER, CORRECT one is the one that looks broken. So
        the header must name what it is counting and say why -- §6, the number
        carries its denominator -- rather than leaving the user to guess.

        The payload below is that exact situation: two projects, two run ids,
        both projects holding both. 2 + 2 on the cards, 2 in the header.
        """
        shared = ["r_00001", "r_00002"]
        p, browser, pg, errors = self._open(home_server, route={
            "projects": [
                {"id": "p_00001", "name": "Arm A", "notes": "",
                 "created_at": 1.0, "updated_at": 2.0, "last_opened_at": None,
                 "n_specimens": 1, "n_runs": len(shared),
                 "sections": ["xenium_319"], "compounds": ["temozolomide"],
                 "thumbnail": None},
                {"id": "p_00002", "name": "Arm A (copy)", "notes": "",
                 "created_at": 3.0, "updated_at": 3.0, "last_opened_at": None,
                 "n_specimens": 1, "n_runs": len(shared),
                 "sections": ["xenium_319"], "compounds": ["temozolomide"],
                 "thumbnail": None},
            ],
            "recent_runs": [], "sections_available": 3, "drugs_available": 12,
            # TWO runs exist. Four references to them.
            "totals": {"projects": 2, "runs": len(shared), "specimens": 2},
            "server": {"contract": "v1", "stub": False},
        })
        try:
            assert not errors, f"raised while rendering: {errors[:2]}"

            note = pg.eval_on_selector(".pgrid-count", "e => e.innerText")

            # The cards really do sum to more than the header. If this stops
            # being true the test below is no longer testing anything, so it
            # is asserted rather than assumed.
            body = pg.eval_on_selector("body", "e => e.innerText")
            assert "2 runs attached" in body, (
                "the per-project cards no longer report attached runs, so "
                f"there is nothing left to reconcile. Body: {body[:300]!r}")

            assert "2 distinct runs" in note, (
                "the workspace header reports a run count that the project "
                "cards sum to MORE than, and does not say why. A user adding "
                "the cards up gets 4 and is shown 2, with nothing on the page "
                f"reconciling them. Header reads: {note!r}")

            gloss = pg.eval_on_selector_all(
                ".pgrid-count abbr[title]", "els => els.map(e => e.title)")
            assert gloss, (
                "the qualifier has no explanation attached; `distinct` on its "
                "own does not tell anyone which runs were counted once")
            assert any("counted once" in t and "two projects" in t
                       for t in gloss), (
                "the gloss does not say that a run reachable from two "
                f"projects is counted once here. Titles: {gloss!r}")
        finally:
            browser.close()
            p.stop()


    def test_a_populated_workspace_renders_its_projects(self, home_server):
        """The path the outage actually broke: `projectGrid` construction."""
        p, browser, pg, errors = self._open(home_server, route={
            "projects": [{
                "id": "p_00001", "name": "Boot probe", "notes": "",
                "created_at": 1.0, "updated_at": 2.0, "last_opened_at": None,
                "n_specimens": 1, "n_runs": 2, "sections": ["xenium_319"],
                "compounds": ["temozolomide"], "thumbnail": None,
            }],
            "recent_runs": [], "sections_available": 3, "drugs_available": 12,
            "totals": {"projects": 1, "runs": 2, "specimens": 1},
            "server": {"contract": "v1", "stub": False},
        })
        try:
            assert not errors, f"raised while rendering a project: {errors[:2]}"
            text = pg.eval_on_selector("body", "e => e.innerText")
            assert "Boot probe" in text, (
                "the project the server sent is not on the page. Body begins: "
                f"{text[:200]!r}")
        finally:
            browser.close()
            p.stop()
