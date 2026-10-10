"""The login portal's front end (``docs/AUTH_CONTRACT.md``).

WHY EVERY TEST HERE OPENS A BROWSER
-----------------------------------
On 2026-09-20 the workspace home page threw on boot and stayed dead for
twenty-four minutes while its acceptance gate reported 28 of 28 green. Every
one of those checks read files, shapes and forbidden phrases; not one of them
opened the page, and a page that renders NOTHING satisfies all of them. So the
assertions below run the real modules in a real browser against a real HTTP
server, and the few that are static say so in their own names.

WHY THE SERVER IS REAL AND NOT ``page.route``
---------------------------------------------
Two of the facts under test are about what goes ON THE WIRE -- whether the
session cookie is attached to an ordinary API call, and whether a stored
``X-API-Key`` is still sent next to it (contract §0). An intercepted request is
the browser's intention; a recorded request is the evidence. The fixture below
is therefore a genuine HTTP server that serves ``app/web`` and answers
``/api/*`` itself, recording every method, path, cookie and header it was
given.

WHAT IS BEING GUARDED, in one line each:
  * §1  the routes and the bodies this front end posts to them
  * §0  the key path still works, because deleting it breaks every workflow,
        render.yaml, and the documented curl flows
  * §5  a session partitions NOTHING -- the clause this feature breaks by
        omission rather than by error, and the one worth the most tests
  * the replacement itself: a person meets a sign-in, never a key field
"""
from __future__ import annotations

import functools
import http.server
import json
import os
import re
import socketserver
import subprocess
import threading
import urllib.parse

import pytest

pytest.importorskip("playwright.sync_api", reason="the portal's tests open a browser")
from playwright.sync_api import sync_playwright  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
WEB = os.path.join(_ROOT, "app", "web")
AUTH_JS = os.path.join(WEB, "js", "auth")
LOGIN_HTML = os.path.join(WEB, "login.html")

#: Everything this agent wrote or edited. Used by the static checks, which name
#: themselves static so nobody mistakes one for a test of the page.
MY_FILES = (
    os.path.join(WEB, "login.html"),
    os.path.join(WEB, "css", "auth.css"),
    os.path.join(AUTH_JS, "shared.js"),
    os.path.join(AUTH_JS, "gate.js"),
    os.path.join(AUTH_JS, "loginPage.js"),
    os.path.join(WEB, "js", "api.js"),
)

needs_portal = pytest.mark.skipif(
    not os.path.exists(LOGIN_HTML), reason="the login portal has not landed yet")


# --------------------------------------------------------------------------- #
# the fixture: a real server that serves app/web AND records every /api/ call
# --------------------------------------------------------------------------- #
_DEFAULT_STATE = {
    # GET /api/auth/me. §1: this route NEVER 401s.
    "me": {"authenticated": False},
    # When true, /api/auth/me denies any request that carries an X-API-Key and
    # accepts the bare one. That is the shape of a server which validates the
    # key header before it looks for a cookie, and it is what hides a good
    # session behind a stale key -- the migration case for every browser that
    # has ever used this product.
    "key_poisons_me": False,
    "login_status": 200,
    "signup_status": 200,
    "workspace_status": 401,
    "sections_status": 401,
    # Every other /api/ route. It has to be settable, because api.js's answer
    # to a 401 on an ORDINARY route is now to navigate to /login -- so a test
    # that wants a signed-in Studio has to let the Studio's whole boot through,
    # not just the one route it is looking at. Leaving this at 401 made four
    # tests assert against a page that had correctly redirected away.
    "default_status": 401,
}

#: The Studio's boot reads these three before it paints anything. Shapes copied
#: from tests/test_studio_boot.py, trimmed to what boot touches.
_DEFAULTS = {
    "section": "xenium_319", "dx_um": 21.0, "drug": "temozolomide",
    "sensitivity_mode": "uniform", "vascular_mode": "discrete",
    "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": None,
    "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0, "n_doses": 1,
           "dose_interval_s": 86400.0},
    "total_time_s": 21600.0,
}
_SECTIONS = [{
    "id": "xenium_319", "label": "Xenium GBM · 319-panel", "n_cells": 40887,
    "extent_um": [8087.8, 5065.4], "panel": "hBrain_254g", "available": True,
    "source": "xenium", "sample": None, "cache_dir": "data_cache/xenium_gbm",
    "assay": "Xenium In Situ", "accession": "10x", "has_outcome_label": False,
    "reason": None, "remedy": [],
}]


class _Portal:
    def __init__(self, base, server):
        self.base = base
        self._server = server

    @property
    def seen(self):
        return self._server.seen

    def reset(self, **state):
        self._server.seen.clear()
        self._server.state = dict(_DEFAULT_STATE)
        self._server.state.update(state)

    def calls(self, path, method=None):
        return [r for r in self.seen
                if r["path"] == path and (method is None or r["method"] == method)]


def _respond(state, rec):
    path, key = rec["path"], rec["key"]
    if path == "/api/auth/me":
        if state["key_poisons_me"] and key:
            return 200, {"authenticated": False}
        return 200, dict(state["me"])
    if path == "/api/auth/login":
        if state["login_status"] == 200:
            return 200, {"authenticated": True, "email": "op@example.org",
                         "via": "session"}
        return state["login_status"], {"detail": "those credentials were refused"}
    if path == "/api/auth/signup":
        if state["signup_status"] == 200:
            return 200, {"authenticated": True, "email": "first@example.org",
                         "via": "session"}
        return state["signup_status"], {"detail": "that token is not this instance"}
    if path == "/api/auth/logout":
        return 200, {"ok": True}
    if path == "/api/workspace":
        if state["workspace_status"] == 200:
            return 200, {"projects": [], "recent_runs": [],
                         "sections_available": 0, "drugs_available": 0,
                         "totals": {"projects": 0, "runs": 0, "specimens": 0},
                         "server": {"contract": "v1", "stub": False}}
        return state["workspace_status"], {"detail": "no accepted credential"}
    if path == "/api/sections":
        if state["sections_status"] == 200:
            return 200, list(_SECTIONS)
        return state["sections_status"], {"detail": "no accepted credential"}
    if path == "/api/defaults":
        if state["default_status"] == 200:
            return 200, dict(_DEFAULTS)
        return state["default_status"], {"detail": "no accepted credential"}
    if path == "/api/account":
        return 200, {"contract": "v1",
                     "auth": {"required": True, "header": "X-API-Key",
                              "scheme": "shared-secret", "matches": bool(key),
                              "caller_fingerprint": "ab12cd34",
                              "key_fingerprint": "ab12cd34", "source": "env"},
                     "usage": {}, "billing": {"status": "not-enabled",
                                              "blocked_on": []}}
    if path == "/api/limits":
        return 200, {"limits": {}}
    if state["default_status"] == 200:
        return 200, []
    return state["default_status"], {"detail": "no accepted credential"}


class _Handler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *a):   # noqa: D102 - the harness is not the subject
        pass

    def _api(self):
        path = urllib.parse.urlparse(self.path).path
        if not path.startswith("/api/"):
            return False
        n = int(self.headers.get("content-length") or 0)
        rec = {
            "method": self.command,
            "path": path,
            "key": self.headers.get("x-api-key"),
            "cookie": self.headers.get("cookie"),
            "body": self.rfile.read(n).decode("utf-8") if n else "",
        }
        self.server.seen.append(rec)
        status, payload = _respond(self.server.state, rec)
        raw = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("content-type", "application/json")
        self.send_header("content-length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)
        return True

    def do_GET(self):           # noqa: N802, D102
        # `/login` and `/studio` as the real app/server mounts them, so the
        # Studio's redirect to /login can be followed all the way to the form.
        path = urllib.parse.urlparse(self.path).path
        if path in ("/login", "/login/"):
            self.path = "/login.html"
        elif path in ("/studio", "/studio/"):
            self.path = "/index.html"
        elif path == "/":
            self.path = "/home.html"
        if not self._api():
            super().do_GET()

    def do_POST(self):          # noqa: N802, D102
        if not self._api():
            self.send_error(405)

    def do_DELETE(self):        # noqa: N802, D102
        if not self._api():
            self.send_error(405)


@pytest.fixture(scope="module")
def portal():
    handler = functools.partial(_Handler, directory=WEB)

    class Server(socketserver.ThreadingTCPServer):
        allow_reuse_address = True
        daemon_threads = True

        def handle_error(self, request, client_address):  # noqa: D102
            pass

    with Server(("127.0.0.1", 0), handler) as httpd:
        httpd.seen = []
        httpd.state = dict(_DEFAULT_STATE)
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        port = httpd.server_address[1]
        try:
            yield _Portal(f"http://127.0.0.1:{port}", httpd)
        finally:
            httpd.shutdown()


@pytest.fixture(scope="module")
def browser():
    with sync_playwright() as p:
        b = p.chromium.launch()
        try:
            yield b
        finally:
            b.close()


class _Page:
    """A page plus the errors it raised. `errors` is the first assertion of
    almost every test here: a page that threw is a page whose remaining
    content is whatever had rendered before the throw."""

    def __init__(self, pg, errors):
        self.pg = pg
        self.errors = errors

    def text(self):
        return self.pg.inner_text("body")


def _open(browser, url, cookies=(), local_storage=None):
    ctx = browser.new_context()
    if cookies:
        ctx.add_cookies(list(cookies))
    pg = ctx.new_page()
    errors: list[str] = []
    pg.on("pageerror", lambda e: errors.append("pageerror: " + str(e)))
    if local_storage:
        # Written before the module runs, because api.js reads the stored key
        # at module top. An IIFE and not a bare arrow function: add_init_script
        # EVALUATES what it is given, so an un-called function expression is a
        # script that does nothing at all -- and a test whose precondition
        # silently did not happen is a test that passes for the wrong reason.
        pg.add_init_script(
            "(() => { try { const o = %s;"
            " for (const k in o) localStorage.setItem(k, o[k]); }"
            " catch (e) {} })()" % json.dumps(local_storage))
    pg.goto(url, wait_until="domcontentloaded")
    pg.wait_for_timeout(1200)
    return _Page(pg, errors), ctx


# --------------------------------------------------------------------------- #
# the replacement itself: a person meets a sign-in, never a key field
# --------------------------------------------------------------------------- #
@needs_portal
class TestTheLoginPageIsASignInAndNotAKeyPrompt:

    def test_it_boots_without_raising(self, browser, portal):
        """The page is executed. Nothing else below means anything if the
        module threw on load, because every one of them would then be reading
        the server-delivered placeholder."""
        portal.reset()
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            assert not page.errors, page.errors[:2]
            assert page.pg.eval_on_selector_all("#auth-form", "e => e.length") == 1, (
                "the sign-in form is not in the document after boot; the page "
                f"shows: {page.text()[:300]!r}")
        finally:
            ctx.close()

    def test_it_asks_for_an_email_and_a_password(self, browser, portal):
        portal.reset()
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            kinds = page.pg.eval_on_selector_all(
                "#auth-form input", "els => els.map(e => [e.id, e.type])")
            assert ["auth-email", "email"] in kinds, kinds
            assert ["auth-password", "password"] in kinds, kinds
        finally:
            ctx.close()

    def test_no_field_on_it_asks_for_the_api_key(self, browser, portal):
        """THE WHOLE POINT. What a visitor used to meet was an input whose
        placeholder read `X-API-Key value` -- a prompt for a shared secret,
        shown to a person who had no way to obtain one. In sign-in mode there
        must be nothing on the page asking for a key.

        The bootstrap token in CLAIM mode is a different thing and is allowed:
        §4 requires it, it is asked of the operator rather than a visitor, and
        it is behind a second click."""
        portal.reset()
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            prompts = page.pg.eval_on_selector_all(
                "input",
                "els => els.map(e => [e.placeholder || '', "
                "(e.labels && e.labels[0] ? e.labels[0].textContent : '')])")
            flat = " ".join(" ".join(p) for p in prompts).lower()
            assert "api-key" not in flat and "api key" not in flat, (
                f"a field on the sign-in form asks for the API key: {prompts}")
            assert "#auth-token" not in page.pg.content(), \
                "the bootstrap-token field is on the default sign-in form"
        finally:
            ctx.close()

    def test_signing_in_posts_the_contract_route_and_shape(self, browser, portal):
        """§1: `POST /api/auth/login`, e-mail and password."""
        portal.reset()
        page, ctx = _open(browser, f"{portal.base}/login?next=%2Fstudio")
        try:
            page.pg.fill("#auth-email", "op@example.org")
            page.pg.fill("#auth-password", "a long enough passphrase")
            page.pg.click("#auth-submit")
            page.pg.wait_for_timeout(900)
            posts = portal.calls("/api/auth/login", "POST")
            assert len(posts) == 1, (
                "exactly one POST /api/auth/login was expected; saw "
                f"{[(r['method'], r['path']) for r in portal.seen]}")
            body = json.loads(posts[0]["body"])
            assert body.get("email") == "op@example.org", body
            assert body.get("password") == "a long enough passphrase", body
            assert not page.errors, page.errors[:2]
        finally:
            ctx.close()

    def test_a_successful_sign_in_goes_where_next_said(self, browser, portal):
        """`?next=` is how the in-app gate gets the reader back to the page
        that refused them. A sign-in that lands somewhere else makes the gate
        a dead end."""
        portal.reset(workspace_status=200, sections_status=200,
                     default_status=200)
        page, ctx = _open(browser, f"{portal.base}/login?next=%2Fstudio")
        try:
            page.pg.fill("#auth-email", "op@example.org")
            page.pg.fill("#auth-password", "a long enough passphrase")
            page.pg.click("#auth-submit")
            page.pg.wait_for_timeout(1500)
            assert "/studio" in page.pg.url, page.pg.url
        finally:
            ctx.close()

    def test_a_hostile_next_is_not_followed(self, browser, portal):
        """`//evil.example` passes a naive startsWith('/') and is a
        protocol-relative URL the browser follows clean off this origin. A
        sign-in page that forwards to somebody else's copy of itself is the
        classic open redirect."""
        portal.reset(workspace_status=200, sections_status=200,
                     default_status=200)
        page, ctx = _open(browser, f"{portal.base}/login?next=%2F%2Fevil.example%2Fx")
        try:
            page.pg.fill("#auth-email", "op@example.org")
            page.pg.fill("#auth-password", "a long enough passphrase")
            page.pg.click("#auth-submit")
            page.pg.wait_for_timeout(1200)
            assert "evil.example" not in page.pg.url, page.pg.url
            assert page.pg.url.startswith(portal.base), page.pg.url
        finally:
            ctx.close()

    def test_a_refusal_stays_on_the_page_and_names_the_reason(self, browser, portal):
        """A wrong password is a 401 on /api/auth/login. api.js's answer to a
        401 is "go to /login" -- and acting on that HERE would navigate the
        reader away from the form they just submitted, to the page they are
        already on, losing the one thing they need to see."""
        portal.reset(login_status=401)
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            # A marker that only survives if this document is never replaced.
            # Checking the URL alone cannot see the failure: a re-route to
            # /login lands on the same URL it started from, so the one
            # observable difference is that the document was torn down.
            page.pg.evaluate("() => { window.__same_document = true; }")
            page.pg.fill("#auth-email", "op@example.org")
            page.pg.fill("#auth-password", "the wrong one entirely")
            page.pg.click("#auth-submit")
            page.pg.wait_for_timeout(900)
            assert "/login" in page.pg.url, (
                f"a failed sign-in navigated away, to {page.pg.url}")
            assert page.pg.evaluate("() => !!window.__same_document"), (
                "the document was reloaded, so the reader lost the reason "
                "their sign-in was refused")
            said = page.pg.inner_text("#auth-status")
            assert "refused" in said.lower(), said
            assert page.pg.eval_on_selector_all("#auth-form", "e => e.length") == 1, \
                "the form is gone, so there is no way to try again"
            assert not page.errors, page.errors[:2]
        finally:
            ctx.close()

    def test_throttling_is_reported_as_the_throttle_and_not_a_fault(
            self, browser, portal):
        """§7 requires failed attempts to be throttled. A 429 reported as a
        server error sends the reader to check the deployment instead of
        waiting."""
        portal.reset(login_status=429)
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            page.pg.fill("#auth-email", "op@example.org")
            page.pg.fill("#auth-password", "the wrong one entirely")
            page.pg.click("#auth-submit")
            page.pg.wait_for_timeout(900)
            said = page.pg.inner_text("#auth-status").lower()
            assert "throttle" in said, said
        finally:
            ctx.close()


@needs_portal
class TestTheFirstRunClaim:
    """§4: sign-up is not open. It takes the operator's own bootstrap token."""

    def test_claim_mode_posts_the_token_to_signup(self, browser, portal):
        portal.reset()
        page, ctx = _open(browser, f"{portal.base}/login?mode=claim")
        try:
            page.pg.fill("#auth-email", "first@example.org")
            page.pg.fill("#auth-password", "a sufficiently long passphrase")
            page.pg.fill("#auth-token", "BOOTSTRAP-TOKEN-ON-THE-WIRE")
            page.pg.click("#auth-submit")
            page.pg.wait_for_timeout(900)
            posts = portal.calls("/api/auth/signup", "POST")
            assert len(posts) == 1, [(r["method"], r["path"]) for r in portal.seen]
            rec = posts[0]
            body = json.loads(rec["body"])
            # §1/§4 pin WHAT the credential is and not where it rides, so it
            # rides both ways and the server reads whichever it implemented.
            assert (rec["key"] == "BOOTSTRAP-TOKEN-ON-THE-WIRE"
                    or body.get("bootstrap_token") == "BOOTSTRAP-TOKEN-ON-THE-WIRE"), (
                "the bootstrap token reached neither the X-API-Key header nor "
                f"the body: header={rec['key']!r} body={body!r}")
            assert body.get("email") == "first@example.org", body
            assert not page.errors, page.errors[:2]
        finally:
            ctx.close()

    def test_claiming_does_not_leave_the_token_stored_as_this_browsers_key(
            self, browser, portal):
        """The token is the operator's deployment secret. Claiming an instance
        from a borrowed laptop must not quietly persist it there."""
        portal.reset()
        page, ctx = _open(browser, f"{portal.base}/login?mode=claim")
        try:
            page.pg.fill("#auth-email", "first@example.org")
            page.pg.fill("#auth-password", "a sufficiently long passphrase")
            page.pg.fill("#auth-token", "BOOTSTRAP-TOKEN-ON-THE-WIRE")
            page.pg.click("#auth-submit")
            page.pg.wait_for_timeout(900)
            stored = page.pg.evaluate(
                "() => [localStorage.getItem('tmesim_api_key'),"
                " sessionStorage.getItem('tmesim_api_key')]")
            assert stored == [None, None], (
                f"the bootstrap token was stored in this browser: {stored}")
        finally:
            ctx.close()

    def test_an_unclaimed_instance_opens_on_the_claim_form(self, browser, portal):
        """When the server says it holds no accounts, the sign-in form cannot
        succeed, and offering it as the default is offering a dead end."""
        portal.reset(me={"authenticated": False, "first_run": True})
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            assert page.pg.eval_on_selector_all("#auth-token", "e => e.length") == 1, (
                "an unclaimed instance did not open on the claim form: "
                f"{page.text()[:300]!r}")
            assert not page.errors, page.errors[:2]
        finally:
            ctx.close()

    def test_the_claim_form_says_sign_up_is_not_open(self, browser, portal):
        """§4's reasoning, not just its mechanism: an internet-reachable
        instance with open registration belongs to whoever finds it first."""
        portal.reset(me={"authenticated": False, "first_run": True})
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            t = page.text().lower()
            assert "not open" in t, t[:400]
            assert "tmesim_studio_api_key" in t, (
                "the claim form does not say WHICH value the token is, so an "
                "operator cannot answer it")
        finally:
            ctx.close()


# --------------------------------------------------------------------------- #
# §5 — a session partitions NOTHING. The clause this feature breaks by
# omission, and the one worth the most tests.
# --------------------------------------------------------------------------- #
@needs_portal
class TestNoSurfaceImpliesAPrivateSpace:

    #: Rendered text, not source. A phrase in a comment is not a phrase a
    #: reader is taught by, and the page is what the reader meets.
    _FORBIDDEN = (
        "your projects", "your project", "your runs", "your workspace",
        "your specimens", "private workspace", "your own workspace",
        "your data",
    )

    def _login_text(self, browser, portal):
        portal.reset()
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            assert not page.errors, page.errors[:2]
            return page.text()
        finally:
            ctx.close()

    def _gate_text(self, browser, portal):
        """The same clause, on the in-app sign-in route. A truth that is only
        on the login page is a truth the signed-out workspace does not tell."""
        portal.reset(workspace_status=401)
        page, ctx = _open(browser, f"{portal.base}/")
        try:
            return page.text()
        finally:
            ctx.close()

    def test_the_scans_have_a_subject(self, browser, portal):
        """Anti-vacuity: a guard that reads an empty page passes for free, and
        this repository has shipped exactly that check before."""
        assert len(self._login_text(browser, portal).strip()) > 900
        assert len(self._gate_text(browser, portal).strip()) > 400

    def test_the_login_page_says_the_store_is_shared(self, browser, portal):
        t = self._login_text(browser, portal).lower()
        assert "way in, not a private space" in t, t[:500]
        assert "every account sees the same" in t, t[:500]

    def test_the_login_page_says_nothing_is_attributed(self, browser, portal):
        """§5 forbids an owner column. The reason it is forbidden is that
        there is no owner: the store records no author, so a reader cannot
        later find out who made a run. That absence has to be stated BEFORE
        they rely on it."""
        t = self._login_text(browser, portal).lower()
        assert "nothing is attributed" in t, t[:600]

    def test_the_in_app_sign_in_route_says_it_too(self, browser, portal):
        t = self._gate_text(browser, portal).lower()
        assert "way in, not a private space" in t, t[:600]

    @pytest.mark.parametrize("phrase", _FORBIDDEN)
    def test_the_login_page_never_implies_privacy(self, browser, portal, phrase):
        t = self._login_text(browser, portal).lower()
        assert phrase not in t, (
            f"{phrase!r} is on the sign-in page over a SHARED store. §5: a UI "
            "that implies privacy teaches a false belief about who can see "
            "what, and that is worse than no login at all.")

    @pytest.mark.parametrize("phrase", _FORBIDDEN)
    def test_the_sign_in_route_never_implies_privacy(self, browser, portal, phrase):
        t = self._gate_text(browser, portal).lower()
        assert phrase not in t, f"{phrase!r} is on the workspace's sign-in route"

    def test_the_privacy_scan_can_actually_fire(self):
        """Pinned in both directions. Without this, narrowing `_FORBIDDEN` to
        nothing would look like an honest UI forever."""
        probe = "everything in your projects is private".lower()
        assert any(p in probe for p in self._FORBIDDEN)
        honest = "an account here is a way in, not a private space".lower()
        assert not any(p in honest for p in self._FORBIDDEN)

    def test_there_is_no_avatar_on_the_sign_in_surfaces(self, browser, portal):
        """§5 names avatars specifically. An avatar is the single strongest
        visual claim that a UI holds per-person things."""
        portal.reset(me={"authenticated": True, "email": "op@example.org",
                         "via": "session"})
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            n = page.pg.eval_on_selector_all(
                "img, svg, [class*=avatar], [class*=gravatar]", "els => els.length")
            assert n == 0, f"{n} image/avatar elements on the signed-in panel"
        finally:
            ctx.close()


# --------------------------------------------------------------------------- #
# api.js: which credential goes on the wire, and what a 401 now means
# --------------------------------------------------------------------------- #
@needs_portal
class TestApiJsCarriesBothCredentials:
    """§0's compatibility bar. The browser cannot CHOOSE between the two --
    `vivum_session` is HttpOnly, so no script can read it or tell that it
    exists -- so both go on the wire and the server decides. These tests read
    what actually arrived, not what the code meant to send."""

    _COOKIE = {"name": "vivum_session", "value": "u_00001.9999999999.deadbeef"}

    def test_the_session_cookie_is_attached_to_an_ordinary_api_call(
            self, browser, portal):
        portal.reset(workspace_status=200, sections_status=200)
        cookie = dict(self._COOKIE, url=portal.base)
        page, ctx = _open(browser, f"{portal.base}/", cookies=[cookie])
        try:
            ws = portal.calls("/api/workspace", "GET")
            assert ws, [(r["method"], r["path"]) for r in portal.seen]
            assert ws[0]["cookie"] and "vivum_session" in ws[0]["cookie"], (
                "GET /api/workspace arrived with no session cookie, so a "
                f"signed-in reader is anonymous to the gate: {ws[0]}")
        finally:
            ctx.close()

    def test_a_stored_key_is_still_sent(self, browser, portal):
        """§0 and §7: deleting or weakening the X-API-Key path is not
        authorised. A reader who only ever uses the key must notice no
        change."""
        portal.reset(workspace_status=200, sections_status=200)
        page, ctx = _open(browser, f"{portal.base}/",
                          local_storage={"tmesim_api_key": "A-REAL-OPERATOR-KEY"})
        try:
            ws = portal.calls("/api/workspace", "GET")
            assert ws, [(r["method"], r["path"]) for r in portal.seen]
            assert ws[0]["key"] == "A-REAL-OPERATOR-KEY", (
                "the stored key did not reach the wire; every workflow, "
                f"render.yaml and documented curl flow is broken: {ws[0]}")
        finally:
            ctx.close()

    def test_both_credentials_ride_together(self, browser, portal):
        """Two credentials, ONE gate. A browser that has both must send both:
        suppressing either one makes this module the thing deciding which
        credential counts, and it has no basis on which to decide."""
        portal.reset(workspace_status=200, sections_status=200)
        cookie = dict(self._COOKIE, url=portal.base)
        page, ctx = _open(browser, f"{portal.base}/", cookies=[cookie],
                          local_storage={"tmesim_api_key": "A-REAL-OPERATOR-KEY"})
        try:
            ws = portal.calls("/api/workspace", "GET")
            assert ws, [(r["method"], r["path"]) for r in portal.seen]
            assert ws[0]["key"] == "A-REAL-OPERATOR-KEY", ws[0]
            assert "vivum_session" in (ws[0]["cookie"] or ""), ws[0]
        finally:
            ctx.close()

    def test_a_stale_key_does_not_mask_a_working_session(self, browser, portal):
        """THE MIGRATION CASE, and it is the common one rather than a corner.

        Every browser that has used this product already holds a key in
        localStorage, and some of those keys are dead. A server that validates
        the header it was given before looking for a cookie refuses the stale
        key and never reaches the perfectly good session beside it. The reader
        then signs in successfully and still sees nothing, which looks like a
        broken login and is really a secret from months ago.
        """
        portal.reset(key_poisons_me=True,
                     me={"authenticated": True, "email": "op@example.org",
                         "via": "session"})
        cookie = dict(self._COOKIE, url=portal.base)
        page, ctx = _open(browser, f"{portal.base}/login", cookies=[cookie],
                          local_storage={"tmesim_api_key": "A-STALE-KEY"})
        try:
            assert not page.errors, page.errors[:2]
            state = page.pg.evaluate(
                "async () => { const m = await import('/js/api.js');"
                " return {auth: m.api.session.authenticated,"
                " sendKey: m.api.sendKey, rejected: m.api.keyRejected,"
                " stored: localStorage.getItem('tmesim_api_key')}; }")
            assert state["auth"] is True, (
                "a good session was hidden behind a stale key: " f"{state}")
            assert state["sendKey"] is False, (
                "the stale key is still being put on the wire, so every gated "
                f"call will keep being refused: {state}")
            assert state["rejected"] is True, state
            # NOT deleted. Refusing a key once is not grounds to destroy an
            # operator's secret, and it may be valid against another host.
            assert state["stored"] == "A-STALE-KEY", state
        finally:
            ctx.close()


@needs_portal
class TestA401IsASignInAndNotAKeyPrompt:

    def test_the_workspace_shows_no_credential_field_at_all(self, browser, portal):
        """THE DEFECT THIS WHOLE TASK IS ABOUT, at the surface a visitor meets
        first. The workspace used to answer a 401 with a full-viewport card
        containing a password input labelled `X-API-Key value`. There must now
        be no field of any kind on the page, and a route to /login instead."""
        portal.reset(workspace_status=401)
        page, ctx = _open(browser, f"{portal.base}/")
        try:
            fields = page.pg.eval_on_selector_all(
                "input[type=password], input[name*=key], input[placeholder*=API]",
                "els => els.map(e => e.placeholder || e.name || e.type)")
            assert fields == [], (
                f"the signed-out workspace still asks for a credential: {fields}")
            hrefs = page.pg.eval_on_selector_all(
                "a[href*='/login']", "els => els.map(e => e.getAttribute('href'))")
            assert hrefs, (
                "the signed-out workspace offers no route to the sign-in page; "
                f"it shows: {page.text()[:400]!r}")
            assert any("next=" in h for h in hrefs), (
                f"the sign-in route does not say where to come back to: {hrefs}")
        finally:
            ctx.close()

    def test_the_workspace_keeps_its_own_chrome(self, browser, portal):
        """It is the landing page. Bouncing the root URL to a login form
        throws away the header, the footer and every honest sentence about
        what this deployment is, before anyone has seen the product."""
        portal.reset(workspace_status=401)
        page, ctx = _open(browser, f"{portal.base}/")
        try:
            assert page.pg.url.rstrip("/") == portal.base.rstrip("/"), (
                f"the workspace navigated away from itself: {page.pg.url}")
            assert "Workspace" in page.text(), page.text()[:300]
        finally:
            ctx.close()

    def test_the_studio_ends_up_at_the_login_page(self, browser, portal):
        """Nothing in the Studio works without a credential, so it does NOT
        keep a half-painted shell; the reader ends up at /login with a way
        back to where they were.

        END TO END, and deliberately NOT a test of the gate. Planting "the
        gate does not navigate" left this green, because the Studio's ordinary
        routes then 401 through `call()` and api.js's own default sends the
        reader to /login anyway. Two mechanisms produce this outcome and that
        redundancy is a feature of the app; it also means this assertion
        cannot attribute the outcome to either one. The test below is the one
        that holds the gate to account.
        """
        portal.reset(sections_status=401, workspace_status=401)
        page, ctx = _open(browser, f"{portal.base}/studio")
        try:
            page.pg.wait_for_timeout(1500)
            assert "/login" in page.pg.url, (
                f"a refused Studio did not route to the sign-in page: {page.pg.url}")
            assert "next=" in page.pg.url, (
                f"the Studio did not say where to come back to: {page.pg.url}")
        finally:
            ctx.close()

    def test_the_gate_itself_routes_and_asks_for_nothing(self, browser, portal):
        """The gate, on its own, with a probe that is nothing but a refusal.

        Run ON the login page, where `requireSignIn()` correctly declines to
        navigate (it is already there) -- which is what makes the gate's own
        output observable instead of being replaced by a page load. What must
        be true: it paints something, that something routes to /login, and it
        contains no field of any kind. `ui/apikey.js`'s card would fail the
        last of those three on its first line.
        """
        portal.reset()
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            out = page.pg.evaluate(
                "async () => {"
                " const g = await import('/js/auth/gate.js');"
                " const done = await Promise.race(["
                "   g.ensureAuthenticated(async () => ({authRequired: true}))"
                "     .then(() => 'resolved'),"
                "   new Promise(r => setTimeout(() => r('pending'), 400))]);"
                " const card = document.querySelector('.gate-overlay');"
                " return {done,"
                "   painted: !!card,"
                "   inputs: card ? card.querySelectorAll('input, textarea, select').length : -1,"
                "   logins: card ? Array.from(card.querySelectorAll(\"a[href*='/login']\"))"
                "       .map(a => a.getAttribute('href')) : [],"
                "   text: card ? card.innerText : ''};"
                " }")
            assert out["painted"], (
                "the gate rendered nothing at all for a refused reader, so "
                "the page is blank with no route off it")
            assert out["inputs"] == 0, (
                "the gate asks a person to type a credential; it is supposed "
                f"to route them to the sign-in page ({out['inputs']} fields)")
            assert out["logins"], (
                f"the gate offers no route to /login: {out['text'][:300]!r}")
            assert out["done"] == "pending", (
                "the gate resolved instead of stopping boot, so the Studio "
                "paints itself over a page that is navigating away")
            assert "api-key" not in out["text"].lower().replace(" ", "-"), out["text"][:300]
        finally:
            ctx.close()


# --------------------------------------------------------------------------- #
# logout
# --------------------------------------------------------------------------- #
@needs_portal
class TestLogout:

    def test_it_posts_the_contract_route_and_forgets_both_credentials(
            self, browser, portal):
        """A browser holding BOTH credentials that clears only the cookie is a
        browser that can still read everything while the UI says it is signed
        out. Of the two ways to resolve that, the UI saying less is worse than
        the browser holding less -- a key is re-pasteable, a false belief
        about who can see what is not."""
        portal.reset(me={"authenticated": True, "email": "op@example.org"})
        page, ctx = _open(browser, f"{portal.base}/login",
                          local_storage={"tmesim_api_key": "A-REAL-OPERATOR-KEY"})
        try:
            out = page.pg.evaluate(
                "async () => { const m = await import('/js/api.js');"
                " const r = await m.logout();"
                " return {r, stored: localStorage.getItem('tmesim_api_key'),"
                " sess: sessionStorage.getItem('tmesim_api_key'),"
                " auth: m.api.session.authenticated}; }")
            assert portal.calls("/api/auth/logout", "POST"), \
                [(r["method"], r["path"]) for r in portal.seen]
            assert out["r"]["sessionCleared"] is True, out
            assert out["r"]["keyForgotten"] is True, out
            assert out["stored"] is None and out["sess"] is None, out
            assert out["auth"] is False, out
        finally:
            ctx.close()

    def test_the_login_page_offers_a_sign_out_when_already_signed_in(
            self, browser, portal):
        portal.reset(me={"authenticated": True, "email": "op@example.org",
                         "via": "session"})
        page, ctx = _open(browser, f"{portal.base}/login")
        try:
            assert not page.errors, page.errors[:2]
            assert page.pg.eval_on_selector_all(
                "#auth-signout", "e => e.length") == 1, page.text()[:400]
            page.pg.click("#auth-signout")
            page.pg.wait_for_timeout(700)
            assert portal.calls("/api/auth/logout", "POST"), \
                [(r["method"], r["path"]) for r in portal.seen]
            said = page.pg.inner_text("#auth-status").lower()
            assert "signed out" in said or "session is closed" in said, said
        finally:
            ctx.close()

    def test_the_account_control_reaches_a_sign_out(self, browser, portal):
        """Required to be reachable from the EXISTING `account` control, which
        is `js/ui/account.js`'s `#btn-account` -- not a new menu beside it."""
        portal.reset(sections_status=200, workspace_status=200,
                     default_status=200,
                     me={"authenticated": True, "email": "op@example.org",
                         "via": "session"})
        page, ctx = _open(browser, f"{portal.base}/studio")
        try:
            page.pg.wait_for_timeout(1200)
            assert page.pg.eval_on_selector_all(
                "#btn-account", "e => e.length") == 1, "no account control"
            page.pg.click("#btn-account")
            page.pg.wait_for_timeout(1000)
            text = page.pg.inner_text("#acct-dlg")
            assert "op@example.org" in text, text[:500]
            assert "Sign out" in text, text[:500]
            # §5, said where a reader is at the exact moment of assuming the
            # opposite: they have just clicked a thing called "account".
            assert "PARTITIONS NOTHING" in text.upper(), text[:700]
        finally:
            ctx.close()


# --------------------------------------------------------------------------- #
# static. Named static so nobody mistakes one for a test of the page.
# --------------------------------------------------------------------------- #
class TestStaticGuards:

    def test_static_every_module_parses_as_an_es_module(self):
        """`node --check` alone silently accepts a duplicate `export const`;
        only `--input-type=module` reads the file as what the browser reads."""
        for path in MY_FILES:
            if not path.endswith(".js"):
                continue
            with open(path, "rb") as fh:
                src = fh.read()
            r = subprocess.run(["node", "--input-type=module", "--check"],
                               input=src, capture_output=True)
            assert r.returncode == 0, (
                f"{os.path.relpath(path, _ROOT)} is not a valid ES module:\n"
                + r.stderr.decode()[:600])

    def test_static_the_api_key_path_is_not_deleted(self):
        """§0 and §7. `.github/workflows/*`, `render.yaml`, the startup banner
        and the documented curl flows all depend on this header, and §7 lists
        "deleting or weakening the X-API-Key path" as not authorised."""
        src = open(os.path.join(WEB, "js", "api.js"), encoding="utf-8").read()
        assert "tmesim_api_key" in src, "api.js no longer reads the stored key"
        assert "x-api-key" in src, "api.js no longer sends the key header"
        assert "credentials: 'same-origin'" in src, \
            "api.js does not attach the session cookie"

    def test_static_no_page_imports_the_key_gate_any_more(self):
        """The replacement has to be complete. `ui/apikey.js` is the module
        whose card asked a visitor for `X-API-Key value`; if either page still
        imports it, the prompt is one 401 away from being back."""
        offenders = []
        for rel in ("js/main.js", "js/home/main.js"):
            src = open(os.path.join(WEB, rel), encoding="utf-8").read()
            for m in re.finditer(r"^\s*import\s[^\n]*apikey\.js", src, re.M):
                offenders.append(f"{rel}: {m.group(0).strip()}")
        assert not offenders, offenders

    def test_static_no_credential_is_hard_coded(self):
        """Never a password, key or secret in a tracked file."""
        pat = re.compile(
            r"(password|api[_-]?key|token|secret)\s*[:=]\s*['\"][^'\"]{8,}['\"]",
            re.I)
        hits = []
        for path in MY_FILES:
            src = open(path, encoding="utf-8").read()
            for m in pat.finditer(src):
                frag = m.group(0)
                # A field NAME, a storage KEY and a placeholder are not values.
                if any(w in frag.lower() for w in
                       ("tmesim_api_key", "bootstrap_token':", "x-api-key",
                        "current-password", "new-password", "$vivum_api_key")):
                    continue
                hits.append(f"{os.path.relpath(path, _ROOT)}: {frag[:70]}")
        assert not hits, hits

    def test_static_the_hard_coded_credential_scan_can_fire(self):
        pat = re.compile(
            r"(password|api[_-]?key|token|secret)\s*[:=]\s*['\"][^'\"]{8,}['\"]",
            re.I)
        assert pat.search("const token = 'hunter2hunter2'")
        assert not pat.search("const token = t")

    @pytest.mark.parametrize("phrase", [
        "digital twin", "personalised medicine", "personalized medicine",
        "machine learning", "predicts response", "treatment plan",
        "diagnostic",
    ])
    def test_static_no_prohibited_framing_on_the_new_surfaces(self, phrase):
        """§8, carried over unchanged and absolute. A new surface is exactly
        where the product's rules get forgotten."""
        for path in MY_FILES:
            src = open(path, encoding="utf-8").read().lower()
            assert phrase not in src, f"{phrase!r} in {os.path.relpath(path, _ROOT)}"

    def test_static_tissue_is_never_uploaded(self):
        """§8: tissue registers BY FILESYSTEM PATH. 'Upload' describes a
        different product with a different security model."""
        for path in MY_FILES:
            src = open(path, encoding="utf-8").read().lower()
            for m in re.finditer(r"upload", src):
                window = src[max(0, m.start() - 120):m.end() + 120]
                assert not any(w in window for w in
                               ("tissue", "section", "specimen", "xenium")), (
                    f"'upload' next to tissue in {os.path.relpath(path, _ROOT)}")

    def test_static_the_framing_scan_has_a_subject(self):
        total = sum(len(open(p, encoding="utf-8").read()) for p in MY_FILES)
        assert total > 20000, f"only {total} chars scanned"
