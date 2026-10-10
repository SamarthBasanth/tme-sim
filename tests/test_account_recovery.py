"""Account recovery: the operator CLI, and the change-password form.

THE TWO GAPS THE LOGIN PORTAL LEFT OPEN
---------------------------------------
``docs/AUTH_CONTRACT.md`` §7 forbids email of any kind and says a forgotten
password is reset by the operator -- and then there was no mechanism for the
operator to do it, so the only route was hand-editing ``users.json``, which
means hand-computing a password hash. §1 lists ``POST /api/auth/password``,
which shipped working and with no front end, because the contract pins the
route and not its body shape. This file guards both of the things that closed
those: ``python -m app.server account ...`` and ``js/auth/password.js``.

EVERY CHECK HERE WAS PLANTED FIRST
----------------------------------
This repository has found seven checks in one week that passed without ever
testing their subject, so a green assertion is not evidence until the defect it
claims to catch has been seen to make it red. Each of these was written by
planting the defect, confirming the plant was REACHED (not just parsed),
watching this file go red, and reverting. The four the brief named, and what
caught each:

====================================  ====================================
plant                                 caught by
====================================  ====================================
the CLI writing the plaintext         ``test_the_cli_writes_a_hash_and_never
password into the store               _the_password``
the reset accepting a password        ``test_a_password_below_the_floor_is
below the floor                       _refused_and_nothing_is_written``
the route accepting a wrong           ``test_a_wrong_current_password_is
current password                      _refused``
the form sending the new password     ``test_the_new_password_never_appears
in a URL                              _in_any_url``
====================================  ====================================

WHY THE FORM'S TESTS OPEN A BROWSER
-----------------------------------
Same reason ``tests/test_login_portal.py`` does, and that file says it best: on
2026-09-20 the home page threw on boot and stayed dead for twenty-four minutes
while 28 of 28 checks reported green, because not one of them opened the page.
A page that renders NOTHING satisfies every static check. So the form is
assembled by the real modules in a real browser against a real HTTP server
that records what went on the wire -- which is also the only way to assert
that a password did not end up in a URL, because that is a fact about the
request and not about the code's intent.

COST
----
No reconstruction, no PDE, no real tissue. The expense is key derivation, which
is the point of a password KDF: one hash is ~155 ms on this interpreter, so
accounts are created with ONE pre-computed hash wherever the KDF is not the
subject.
"""
from __future__ import annotations

import functools
import getpass
import http.server
import io
import json
import os
import re
import socketserver
import stat
import subprocess
import sys
import threading
import urllib.parse

import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app.server import auth as auth_mod                       # noqa: E402
from app.server import __main__ as server_main                # noqa: E402
from app.server import cache as cache_mod                     # noqa: E402
from app.server import compounds as compounds_mod             # noqa: E402
from app.server import ingest_registry                        # noqa: E402
from app.server.api import create_app                         # noqa: E402
from app.server.jobs import JobManager                        # noqa: E402
from app.server.projects import ProjectStore                  # noqa: E402

from starlette.testclient import TestClient                   # noqa: E402

#: The bootstrap token (§4) for the gated app.
API_KEY = "k" * 48
#: Twelve characters exactly, which is §3's floor -- so the boundary is
#: exercised by the ordinary happy path rather than by one test nobody reads.
PASSWORD = "correct-horse"
NEW_PASSWORD = "a-brand-new-passphrase"
EMAIL = "Operator@Example.ORG"
EMAIL_CASEFOLDED = "operator@example.org"
LOOPBACK = ("127.0.0.1", 51000)

WEB = os.path.join(_REPO, "app", "web")
PASSWORD_JS = os.path.join(WEB, "js", "auth", "password.js")
ACCOUNT_JS = os.path.join(WEB, "js", "ui", "account.js")

#: Mine. Used by the static checks, which name themselves static so nobody
#: mistakes one for a test of the page.
MY_FILES = (PASSWORD_JS, ACCOUNT_JS,
            os.path.join(_REPO, "app", "server", "__main__.py"),
            os.path.join(_REPO, "app", "server", "auth.py"))


# --------------------------------------------------------------------------- #
# Fixtures: a real store, and a real gated app over it
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="session")
def known_hash() -> str:
    """``hash_password(PASSWORD)``, computed once for the whole session."""
    return auth_mod.hash_password(PASSWORD)


@pytest.fixture
def store_root(tmp_path, monkeypatch):
    """A private ``TMESIM_STATE_DIR``, which is what carries the store (§6)."""
    for var in ("TMESIM_STUDIO_DATA_ROOTS", "TMESIM_STUDIO_EXPOSE_TRACEBACKS",
                "TMESIM_STUDIO_ALLOWED_ORIGINS", "TMESIM_STUDIO_ALLOWED_HOSTS",
                "TMESIM_SESSION_SECRET"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("TMESIM_STATE_DIR", str(tmp_path / "state"))
    return tmp_path


@pytest.fixture
def store(store_root):
    return auth_mod.UserStore()


@pytest.fixture
def gated(store_root, monkeypatch):
    """The real app, gate on, peer on loopback -- so the CLI and a live server
    can be driven against ONE store and the interaction asserted."""
    monkeypatch.setenv("TMESIM_STUDIO_API_KEY", API_KEY)
    runs = store_root / "runs"
    runs.mkdir(exist_ok=True)
    monkeypatch.setattr(cache_mod, "runs_dir", lambda: str(runs))
    monkeypatch.setattr(ingest_registry, "_store", None)
    monkeypatch.setattr(
        compounds_mod, "COMPOUND_STORE",
        compounds_mod.CompoundStore(root=str(store_root / "compounds")))
    manager = JobManager(root=str(runs / "jobs"), start_worker=False)
    app = create_app(manager, projects=ProjectStore(root=str(store_root / "p")))
    with TestClient(app, client=LOOPBACK) as client:
        yield client
    manager.shutdown()
    monkeypatch.setattr(ingest_registry, "_store", None)


def _cli(*argv) -> int:
    """Run one operator command in-process and return its exit code."""
    return server_main.main(list(argv))


def _stdin(monkeypatch, text: str):
    """A non-tty stdin holding ``text`` -- the ``--stdin`` path."""
    monkeypatch.setattr(sys, "stdin", io.StringIO(text))


def _tty_prompts(monkeypatch, *answers):
    """A tty stdin, and ``getpass`` answering ``answers`` in order."""
    class _Tty(io.StringIO):
        def isatty(self):
            return True
    monkeypatch.setattr(sys, "stdin", _Tty())
    queue = list(answers)
    monkeypatch.setattr(getpass, "getpass",
                        lambda prompt="": queue.pop(0) if queue else "")


# --------------------------------------------------------------------------- #
# 1. The CLI exists, and is the operator's -- not a route
# --------------------------------------------------------------------------- #
class TestTheOperatorPathExists:
    """§7 named this mechanism and did not provide one."""

    def test_the_four_commands_the_operator_needs_are_there(self):
        parser = server_main.build_parser()
        for argv in (["account", "list"],
                     ["account", "reset-password", "--email", EMAIL],
                     ["account", "disable", "--email", EMAIL],
                     ["account", "enable", "--email", EMAIL]):
            args = parser.parse_args(argv)
            assert args.command == "account", argv
            assert callable(getattr(args, "account_func", None)), argv

    def test_no_subcommand_still_starts_the_server(self):
        """BACKWARDS COMPATIBILITY, and not a nicety: tmesim/cli.py::cmd_studio
        calls main() with an argv of bare flags, and so do the Dockerfile,
        render.yaml and both documented deploy flows. A required `serve`
        subcommand would have broken every one of them."""
        args = server_main.build_parser().parse_args(
            ["--host", "0.0.0.0", "--port", "9000", "--api-key", "x"])
        assert getattr(args, "command", None) is None
        assert (args.host, args.port, args.api_key) == ("0.0.0.0", 9000, "x")

    def test_there_is_no_password_flag_on_any_subcommand(self):
        """A password on an argv is readable by every user on the host via
        `ps`, is written to the shell's history file, and is printed into any
        orchestration tool's job log. None of that is undone by changing the
        password later, because the value is already written down somewhere
        nobody is watching. So the flag must not exist -- not merely be
        discouraged in a docstring."""
        parser = server_main.build_parser()
        for argv in (["account", "reset-password", "--email", EMAIL,
                      "--password", "whatever"],
                     ["account", "reset-password", "--email", EMAIL,
                      "--new-password", "whatever"]):
            with pytest.raises(SystemExit):
                parser.parse_args(argv)

    def test_no_http_route_resets_a_password(self, gated):
        """The authority for these commands is filesystem ownership, which is
        not something an HTTP request can present. So there must be no route
        that does this: the only candidate authority the product has is
        TMESIM_STUDIO_API_KEY, and §4 already spends that as the sign-up
        bootstrap token -- making it a master credential over identities is
        more than §4 gave it."""
        paths = {getattr(r, "path", "") for r in gated.app.routes}
        for banned in ("/api/auth/reset", "/api/auth/reset-password",
                       "/api/accounts", "/api/auth/accounts",
                       "/api/auth/disable", "/api/admin"):
            assert banned not in paths, f"{banned} is served"
        # And the operator commands are not reachable from the app at all.
        for name in ("cmd_account_reset_password", "cmd_account_disable",
                     "cmd_account_list"):
            assert hasattr(server_main, name)
        # Read as IMPORTS and not as text: api.py mentions `__main__.py` in a
        # comment about where TMESIM_STUDIO_ALLOWED_HOSTS is set, which a
        # substring check flags and which is not an import of anything.
        import ast
        import app.server.api as api_mod
        tree = ast.parse(open(api_mod.__file__, encoding="utf-8").read())
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported |= {a.name for a in node.names}
            elif isinstance(node, ast.ImportFrom):
                imported.add(node.module or "")
                imported |= {a.name for a in node.names}
        assert "__main__" not in imported, (
            "api.py imports the operator CLI module, which would put the "
            f"reset path one import away from the HTTP surface: {imported}")

    def test_the_cli_module_opens_no_socket_and_serves_nothing(self):
        """An operator command must not be able to become a service by
        accident. Read structurally: no account code path reaches uvicorn."""
        import ast
        import inspect
        tree = ast.parse(inspect.getsource(server_main))
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            if not node.name.startswith(("cmd_account", "_account",
                                         "_read_new_password",
                                         "_private_store_problem")):
                continue
            body = ast.dump(node)
            for banned in ("uvicorn", "create_app", "APIRouter"):
                assert banned not in body, (
                    f"{node.name} references {banned}")


# --------------------------------------------------------------------------- #
# 2. The reset: what it writes
# --------------------------------------------------------------------------- #
class TestTheResetWritesAHash:

    def test_the_cli_writes_a_hash_and_never_the_password(
            self, store, known_hash, monkeypatch, capsys):
        """THE PLANT: ``store.set_password(row["id"], password)`` -- the
        plaintext straight into the store, which is the single worst thing a
        credential tool can do and is one deleted call away at all times."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        _stdin(monkeypatch, NEW_PASSWORD + "\n")
        assert _cli("account", "reset-password", "--email", EMAIL,
                    "--stdin") == server_main.EXIT_OK

        raw = open(store.path, encoding="utf-8").read()
        assert NEW_PASSWORD not in raw, (
            "the new password is in the store in plaintext")
        assert PASSWORD not in raw
        row = store.by_email(EMAIL_CASEFOLDED)
        stored = row["password"]
        assert stored != NEW_PASSWORD
        assert "$" in stored and stored.split("$")[0] in (
            "scrypt", "pbkdf2-sha256"), stored.split("$")[0]
        # The hash is useless if the running interpreter cannot verify it.
        assert auth_mod.verify_password(NEW_PASSWORD, stored) is True
        assert auth_mod.verify_password(PASSWORD, stored) is False
        # Nor may it reach the terminal, which scrolls into a log, a screen
        # share and a bug report.
        out = capsys.readouterr()
        assert NEW_PASSWORD not in out.out and NEW_PASSWORD not in out.err
        assert stored not in out.out and stored not in out.err

    def test_it_does_not_assume_scrypt(self, store, known_hash, monkeypatch):
        """``hashlib.scrypt`` does not exist on this interpreter (macOS links
        LibreSSL, which has no ``EVP_PBE_scrypt``). A CLI that assumed it would
        either crash or write a hash this interpreter cannot verify -- so the
        account would be locked out BY THE RECOVERY TOOL. Forced, rather than
        relying on this machine happening to lack it."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        monkeypatch.setattr(auth_mod, "scrypt_available", lambda: False)
        _stdin(monkeypatch, NEW_PASSWORD + "\n")
        assert _cli("account", "reset-password", "--email", EMAIL,
                    "--stdin") == server_main.EXIT_OK
        stored = store.by_email(EMAIL_CASEFOLDED)["password"]
        assert stored.startswith("pbkdf2-sha256$"), stored.split("$")[0]
        assert auth_mod.verify_password(NEW_PASSWORD, stored) is True

    @pytest.mark.parametrize("password,why", [
        ("short", "below twelve characters"),
        ("", "empty"),
        ("elevenchar", "eleven characters"),
    ])
    def test_a_password_below_the_floor_is_refused_and_nothing_is_written(
            self, store, known_hash, monkeypatch, capsys, password, why):
        """THE PLANT: drop the ``password_problem`` call from the reset, so the
        CLI mints a password the product's own sign-in rules would refuse.

        §3's floor is twelve characters with no composition rules, and the
        rule must come from ``auth.password_problem`` rather than be re-typed
        here -- a CLI with its own floor is a CLI that can write a credential
        the routes reject."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        _stdin(monkeypatch, password + "\n")
        rc = _cli("account", "reset-password", "--email", EMAIL, "--stdin")
        assert rc == server_main.EXIT_REFUSED, why
        assert store.by_email(EMAIL_CASEFOLDED)["password"] == known_hash, (
            "the store was written despite the refusal")
        said = capsys.readouterr()
        assert "nothing was written" in said.err.lower(), said.err

    def test_twelve_characters_is_accepted(self, store, known_hash,
                                           monkeypatch):
        """The floor is a floor and not a fence: the boundary value passes, or
        the refusal above is not testing the boundary."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        _stdin(monkeypatch, "twelvechars!" + "\n")
        assert len("twelvechars!") == auth_mod.PASSWORD_MIN_CHARS
        assert _cli("account", "reset-password", "--email", EMAIL,
                    "--stdin") == server_main.EXIT_OK

    def test_an_unknown_address_is_refused_without_creating_one(
            self, store, known_hash, monkeypatch):
        """A reset must never be a sign-up. §4 governs account creation and
        requires the bootstrap token; a CLI that created an account on a typo
        would be a second way in that §4 never authorised."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        _stdin(monkeypatch, NEW_PASSWORD + "\n")
        assert _cli("account", "reset-password", "--email",
                    "typo@example.org", "--stdin") == server_main.EXIT_REFUSED
        assert store.count() == 1
        assert store.by_email("typo@example.org") is None

    def test_the_address_is_matched_casefolded(self, store, known_hash,
                                               monkeypatch):
        """§6 stores and matches casefolded, so an operator typing the address
        in the case they remember must reach the account."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        _stdin(monkeypatch, NEW_PASSWORD + "\n")
        assert _cli("account", "reset-password", "--email",
                    "OPERATOR@EXAMPLE.ORG", "--stdin") == server_main.EXIT_OK
        assert auth_mod.verify_password(
            NEW_PASSWORD, store.by_email(EMAIL_CASEFOLDED)["password"])

    def test_the_store_stays_owner_only_after_a_cli_write(
            self, store, known_hash, monkeypatch):
        """§6: ``0600``. Written with the default umask it would land ``0644``,
        world-readable, which for a file of password hashes and the session
        signing key is the whole game."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        _stdin(monkeypatch, NEW_PASSWORD + "\n")
        assert _cli("account", "reset-password", "--email", EMAIL,
                    "--stdin") == server_main.EXIT_OK
        mode = stat.S_IMODE(os.stat(store.path).st_mode)
        assert mode == 0o600, oct(mode)

    def test_a_world_readable_store_is_refused_rather_than_written_into(
            self, store, known_hash, monkeypatch, capsys):
        """A fresh hash must not be written into a file other users on this
        host can already read. Refusing rather than silently tightening it,
        because a store that is readable NOW has already been exposed and the
        one person who needs to know is the one running this command."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        os.chmod(store.path, 0o644)
        _stdin(monkeypatch, NEW_PASSWORD + "\n")
        rc = _cli("account", "reset-password", "--email", EMAIL, "--stdin")
        assert rc == server_main.EXIT_CONFIG
        assert store.by_email(EMAIL_CASEFOLDED)["password"] == known_hash
        assert "chmod 600" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# 3. The reset: how the password gets in
# --------------------------------------------------------------------------- #
class TestThePasswordIsNeverEchoed:

    def test_the_prompt_asks_twice_and_a_mismatch_writes_nothing(
            self, store, known_hash, monkeypatch, capsys):
        """A typo in a value nobody can see locks out the account this command
        exists to recover."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        _tty_prompts(monkeypatch, NEW_PASSWORD, NEW_PASSWORD + "x")
        rc = _cli("account", "reset-password", "--email", EMAIL)
        assert rc == server_main.EXIT_REFUSED
        assert store.by_email(EMAIL_CASEFOLDED)["password"] == known_hash
        assert "did not match" in capsys.readouterr().err

    def test_a_matching_pair_from_the_prompt_is_accepted(
            self, store, known_hash, monkeypatch):
        store.create(EMAIL_CASEFOLDED, known_hash)
        _tty_prompts(monkeypatch, NEW_PASSWORD, NEW_PASSWORD)
        assert _cli("account", "reset-password",
                    "--email", EMAIL) == server_main.EXIT_OK
        assert auth_mod.verify_password(
            NEW_PASSWORD, store.by_email(EMAIL_CASEFOLDED)["password"])

    def test_it_refuses_rather_than_echo_when_echo_cannot_be_turned_off(
            self, store, known_hash, monkeypatch, capsys):
        """THE NON-OBVIOUS ONE. ``getpass`` does not fail when it cannot
        control the terminal -- it falls back to ``input()`` and merely WARNS,
        which means it echoes the password onto the screen and into the
        scrollback. That is the exact outcome the prompt exists to prevent, so
        the warning is promoted to an error and the command refuses."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        class _Tty(io.StringIO):
            def isatty(self):
                return True
        monkeypatch.setattr(sys, "stdin", _Tty())

        def _warns(prompt=""):
            import warnings
            warnings.warn("Can not control echo on the terminal",
                          getpass.GetPassWarning)
            return NEW_PASSWORD
        monkeypatch.setattr(getpass, "getpass", _warns)
        rc = _cli("account", "reset-password", "--email", EMAIL)
        assert rc == server_main.EXIT_REFUSED
        assert store.by_email(EMAIL_CASEFOLDED)["password"] == known_hash
        assert "echo" in capsys.readouterr().err.lower()

    def test_stdin_on_a_terminal_is_refused(self, store, known_hash,
                                            monkeypatch, capsys):
        """``--stdin`` with a terminal attached would echo every character
        typed, which is what the flag exists to avoid."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        class _Tty(io.StringIO):
            def isatty(self):
                return True
        monkeypatch.setattr(sys, "stdin", _Tty())
        rc = _cli("account", "reset-password", "--email", EMAIL, "--stdin")
        assert rc == server_main.EXIT_REFUSED
        assert "echo" in capsys.readouterr().err.lower()

    def test_no_terminal_and_no_stdin_flag_names_the_flag(
            self, store, known_hash, monkeypatch, capsys):
        store.create(EMAIL_CASEFOLDED, known_hash)
        _stdin(monkeypatch, "")
        rc = _cli("account", "reset-password", "--email", EMAIL)
        assert rc == server_main.EXIT_REFUSED
        assert "--stdin" in capsys.readouterr().err

    def test_only_the_trailing_newline_is_stripped_from_a_piped_password(
            self, store, known_hash, monkeypatch):
        """``.strip()`` here would set a DIFFERENT password from the one that
        was piped in, and the operator would learn that at the sign-in that
        refuses it, with no way to tell why. A password may legitimately begin
        or end with a space."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        spaced = "  padded passphrase  "
        _stdin(monkeypatch, spaced + "\n")
        assert _cli("account", "reset-password", "--email", EMAIL,
                    "--stdin") == server_main.EXIT_OK
        stored = store.by_email(EMAIL_CASEFOLDED)["password"]
        assert auth_mod.verify_password(spaced, stored) is True
        assert auth_mod.verify_password(spaced.strip(), stored) is False


# --------------------------------------------------------------------------- #
# 4. What a reset does and does NOT end
# --------------------------------------------------------------------------- #
class TestWhatAResetEnds:
    """§2 says there is no per-session revocation. This is the command where
    that is most likely to be assumed away, because every other product an
    operator has used signs the other browsers out on a reset."""

    def _signed_in(self, client, store, known_hash):
        store.create(EMAIL_CASEFOLDED, known_hash)
        res = client.post("/api/auth/login",
                          json={"email": EMAIL, "password": PASSWORD})
        assert res.status_code == 200, res.text
        return client.cookies.get(auth_mod.COOKIE_NAME)

    def test_a_reset_does_not_end_an_existing_session(
            self, gated, store, known_hash, monkeypatch):
        """THE HONEST CLAIM, ASSERTED. If this ever starts failing, the CLI's
        output is a lie and must change with it."""
        token = self._signed_in(gated, store, known_hash)
        assert gated.get("/api/workspace").status_code == 200

        _stdin(monkeypatch, NEW_PASSWORD + "\n")
        assert _cli("account", "reset-password", "--email", EMAIL,
                    "--stdin") == server_main.EXIT_OK

        # Same cookie, same running server: still in. This is not a defect to
        # be fixed here -- it is what a stateless signed session means (§2),
        # and the command says so in its own output.
        assert gated.get("/api/workspace",
                         cookies={auth_mod.COOKIE_NAME: token}
                         ).status_code == 200
        # The old password, however, is gone.
        gated.cookies.clear()
        assert gated.post("/api/auth/login",
                          json={"email": EMAIL,
                                "password": PASSWORD}).status_code == 401

    def test_the_reset_output_says_nobody_was_signed_out(self):
        """The sentence is the deliverable. An operator who resets a password
        during an incident, believes the intruder is out and goes to bed is the
        failure this text exists to prevent -- so it is asserted, not
        reviewed."""
        said = server_main._reset_consequences(EMAIL_CASEFOLDED).lower()
        assert "nobody was signed out" in said
        assert "still signed in" in said
        assert "14 days" in said
        assert "account disable" in said
        assert "rotate-session-key" in said
        # And it must not claim the key path changed (§0).
        assert "x-api-key" in said and "untouched" in said

    def test_disable_ends_a_live_session_at_once(
            self, gated, store, known_hash):
        """The ONE immediate lever over a single account, and the reason it
        works is that ``user_for_request`` re-reads the row on EVERY request
        rather than trusting the cookie."""
        token = self._signed_in(gated, store, known_hash)
        assert gated.get("/api/workspace").status_code == 200

        assert _cli("account", "disable", "--email",
                    EMAIL) == server_main.EXIT_OK
        assert store.by_email(EMAIL_CASEFOLDED)["disabled"] is True
        res = gated.get("/api/workspace",
                        cookies={auth_mod.COOKIE_NAME: token})
        assert res.status_code == 401, res.text

    def test_enable_lets_the_same_cookie_back_in(
            self, gated, store, known_hash):
        """Asserted because the CLI claims it, and it is the surprising half:
        disable/enable is NOT a way to force a re-login."""
        token = self._signed_in(gated, store, known_hash)
        assert _cli("account", "disable", "--email", EMAIL) == 0
        assert gated.get("/api/workspace",
                         cookies={auth_mod.COOKIE_NAME: token}
                         ).status_code == 401
        assert _cli("account", "enable", "--email", EMAIL) == 0
        assert gated.get("/api/workspace",
                         cookies={auth_mod.COOKIE_NAME: token}
                         ).status_code == 200

    def test_disable_does_not_change_the_password(self, store, known_hash):
        store.create(EMAIL_CASEFOLDED, known_hash)
        assert _cli("account", "disable", "--email", EMAIL) == 0
        assert store.by_email(EMAIL_CASEFOLDED)["password"] == known_hash

    def test_a_disabled_account_is_refused_at_sign_in_without_saying_why(
            self, gated, store, known_hash):
        """"that address is disabled here" is the same existence disclosure the
        login route goes to some trouble not to make, in a politer voice."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        assert _cli("account", "disable", "--email", EMAIL) == 0
        res = gated.post("/api/auth/login",
                         json={"email": EMAIL, "password": PASSWORD})
        assert res.status_code == 401
        body = json.dumps(res.json()).lower()
        for leak in ("disabled", "disable", "suspended", "locked"):
            assert leak not in body, body

    def test_the_disable_output_states_both_halves(self):
        off = server_main._disable_consequences(EMAIL_CASEFOLDED, True).lower()
        assert "in effect now" in off
        assert "password is unchanged" in off
        on = server_main._disable_consequences(EMAIL_CASEFOLDED, False).lower()
        assert "also back in" in on, (
            "enable's output does not warn that old cookies come back")

    def test_rotating_the_session_key_invalidates_every_session(
            self, gated, store, known_hash):
        """The only revocation §2 permits, and it is all-or-nothing. Asserted
        in both halves: the key really does change, AND the running process
        keeps working until it restarts (which is why the command says so)."""
        token = self._signed_in(gated, store, known_hash)
        live = gated.app.state.auth
        assert auth_mod.read_session(token, live.secret()) is not None

        assert _cli("account", "rotate-session-key",
                    "--yes") == server_main.EXIT_OK

        # A fresh process reading the same store: the cookie is dead.
        after = auth_mod.AuthService(store=auth_mod.UserStore())
        assert auth_mod.read_session(token, after.secret()) is None
        # The RUNNING one cached its key at first use, on purpose (§2), so it
        # still honours the cookie -- hence "restart the server".
        assert auth_mod.read_session(token, live.secret()) is not None
        # Nothing about the account changed.
        assert store.by_email(EMAIL_CASEFOLDED)["password"] == known_hash

    def test_rotating_needs_an_explicit_yes(self, store):
        assert _cli("account",
                    "rotate-session-key") == server_main.EXIT_REFUSED

    def test_rotating_refuses_when_the_env_var_owns_the_key(
            self, store, monkeypatch, capsys):
        """§2: ``TMESIM_SESSION_SECRET`` wins, and the store's instance_secret
        is then never consulted. Rotating it would succeed, print a success and
        end nothing -- so the command refuses rather than report a thing it did
        not do."""
        monkeypatch.setenv("TMESIM_SESSION_SECRET", "x" * 32)
        before = store.instance_secret()
        rc = _cli("account", "rotate-session-key", "--yes")
        assert rc == server_main.EXIT_REFUSED
        assert store.instance_secret() == before, "it rotated anyway"
        assert "TMESIM_SESSION_SECRET" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# 5. `list` -- and what it may not print
# --------------------------------------------------------------------------- #
class TestList:

    def test_it_prints_the_accounts_and_no_hash(self, store, known_hash,
                                                capsys):
        """The rows come from ``operator_user``, which cannot carry a hash --
        rather than from the raw rows with the caller trusted to omit it. A
        CLI's output goes into a terminal that scrolls into a log, a screen
        share and a bug report."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        store.create("second@example.org", known_hash)
        assert _cli("account", "list") == server_main.EXIT_OK
        out = capsys.readouterr().out
        assert EMAIL_CASEFOLDED in out and "second@example.org" in out
        assert "u_00001" in out and "u_00002" in out
        assert known_hash not in out
        for fragment in known_hash.split("$"):
            if len(fragment) > 8:
                assert fragment not in out, "part of the hash was printed"

    def test_the_operator_projection_carries_no_password_field(
            self, store, known_hash):
        store.create(EMAIL_CASEFOLDED, known_hash)
        for row in store.list_accounts():
            assert "password" not in row, row
            assert set(row) == {"id", "email", "created_at", "last_login_at",
                                "disabled"}, row

    def test_disabled_is_visible_to_the_operator_and_not_to_a_client(
            self, gated, store, known_hash, capsys):
        """The asymmetry is the point: the operator has to see which accounts
        are off, and `/api/auth/me` must not say -- which is why this is a
        second allowlist rather than a flag on ``public_user``."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        assert _cli("account", "disable", "--email", EMAIL) == 0
        capsys.readouterr()
        assert _cli("account", "list") == 0
        assert "DISABLED" in capsys.readouterr().out

        gated.post("/api/auth/login", json={"email": EMAIL,
                                            "password": PASSWORD})
        body = gated.get("/api/auth/me").json()
        assert "disabled" not in json.dumps(body), body

    def test_it_names_the_store_it_is_reading(self, store, known_hash,
                                              capsys):
        """THE LIKELIEST OPERATIONAL MISTAKE with this CLI is not a dangerous
        flag -- it is editing the wrong store and believing it worked. §6 falls
        back to app/auth/ when TMESIM_STATE_DIR is unset, so an operator whose
        server runs from systemd (variable set) and whose shell does not
        (variable unset) resets a password in a file no server will read."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        assert _cli("account", "list") == 0
        out = capsys.readouterr().out
        assert store.path in out

    def test_an_unset_state_dir_warns_that_the_server_may_differ(
            self, monkeypatch, tmp_path, capsys):
        monkeypatch.delenv("TMESIM_STATE_DIR", raising=False)
        monkeypatch.setattr(auth_mod, "auth_dir",
                            lambda: str(tmp_path / "fallback"))
        assert _cli("account", "list") == 0
        out = capsys.readouterr().out
        assert "UNSET" in out
        assert "DIFFERENT STORE" in out.upper()


# --------------------------------------------------------------------------- #
# 6. The change-password ROUTE. The form is built against this shape.
# --------------------------------------------------------------------------- #
class TestTheRouteTheFormIsBuiltAgainst:
    """The contract pins the route and not the body, which is exactly why the
    front end was left out of the first pass. These assertions are the shape
    the form was then written against -- so if the route's fields move, the
    form's test fails HERE, naming the cause, instead of in a browser."""

    def _signed_in(self, client, store, known_hash):
        store.create(EMAIL_CASEFOLDED, known_hash)
        assert client.post("/api/auth/login",
                           json={"email": EMAIL,
                                 "password": PASSWORD}).status_code == 200

    def test_the_body_is_old_password_and_new_password(
            self, gated, store, known_hash):
        self._signed_in(gated, store, known_hash)
        res = gated.post("/api/auth/password",
                         json={"old_password": PASSWORD,
                               "new_password": NEW_PASSWORD})
        assert res.status_code == 200, res.text
        assert res.json()["ok"] is True

    @pytest.mark.parametrize("body", [
        {"old": PASSWORD, "new": NEW_PASSWORD},
        {"current_password": PASSWORD, "new_password": NEW_PASSWORD},
        {"password": PASSWORD, "new_password": NEW_PASSWORD},
        {"new_password": NEW_PASSWORD},
        {"old_password": PASSWORD},
    ])
    def test_any_other_field_naming_is_a_400(self, gated, store, known_hash,
                                             body):
        """The guess the front end refused to make. Each of these is a
        plausible name for the same field and every one of them is refused, so
        a form built on a guess would have been a form that never worked."""
        self._signed_in(gated, store, known_hash)
        res = gated.post("/api/auth/password", json=body)
        assert res.status_code == 400, (body, res.status_code)

    def test_a_wrong_current_password_is_refused(self, gated, store,
                                                 known_hash):
        """THE PLANT: drop the ``verify_password(old_password, ...)`` check, so
        anyone holding a session can set a new password without knowing the
        current one -- which turns a stolen cookie into permanent ownership of
        the account."""
        self._signed_in(gated, store, known_hash)
        res = gated.post("/api/auth/password",
                         json={"old_password": "not-the-password",
                               "new_password": NEW_PASSWORD})
        assert res.status_code == 401, res.text
        # And nothing changed: the old password still signs in.
        gated.cookies.clear()
        assert gated.post("/api/auth/login",
                          json={"email": EMAIL,
                                "password": PASSWORD}).status_code == 200

    def test_the_refusal_envelope_is_what_the_form_reads(
            self, gated, store, known_hash):
        """``{detail: {message, reason?}}`` -- an OBJECT, which is why
        password.js reads the envelope itself instead of going through
        api.js::detailOf (that helper handles a string or an array detail and
        returns '' for an object, which would drop the server's sentence and
        leave the reader with "401 Unauthorized")."""
        self._signed_in(gated, store, known_hash)
        res = gated.post("/api/auth/password",
                         json={"old_password": "not-the-password",
                               "new_password": NEW_PASSWORD})
        detail = res.json()["detail"]
        assert isinstance(detail, dict), detail
        assert isinstance(detail["message"], str) and detail["message"]
        import app.web as _unused  # noqa: F401 - path sanity only
        # The exact thing the form displays.
        assert "current password" in detail["message"]

    def test_the_success_body_carries_the_sentence_the_form_renders(
            self, gated, store, known_hash):
        """The form renders ``other_sessions_note`` VERBATIM rather than
        re-typing it, so the two cannot drift. If the field goes away the form
        falls back -- but this is what keeps it from having to."""
        self._signed_in(gated, store, known_hash)
        body = gated.post("/api/auth/password",
                          json={"old_password": PASSWORD,
                                "new_password": NEW_PASSWORD}).json()
        assert body["other_sessions_invalidated"] is False
        assert isinstance(body["other_sessions_note"], str)
        assert body["other_sessions_note"]

    def test_the_key_alone_cannot_change_a_password(self, gated, store,
                                                   known_hash):
        """§4 gave the key a bootstrap-token role, not authority over
        identities. A key-only caller has no own password to change."""
        store.create(EMAIL_CASEFOLDED, known_hash)
        res = gated.post("/api/auth/password",
                         headers={"X-API-Key": API_KEY},
                         json={"old_password": PASSWORD,
                               "new_password": NEW_PASSWORD})
        assert res.status_code == 401, res.text
        assert store.by_email(EMAIL_CASEFOLDED)["password"] == known_hash


# --------------------------------------------------------------------------- #
# the browser fixture: a real server that serves app/web AND records every
# request target, query string included
# --------------------------------------------------------------------------- #
playwright_api = pytest.importorskip(
    "playwright.sync_api", reason="the form's tests open a browser")
from playwright.sync_api import sync_playwright  # noqa: E402

_DEFAULT_STATE = {
    "me": {"authenticated": True, "email": EMAIL_CASEFOLDED, "via": "session"},
    "password_status": 200,
    "password_body": {
        "ok": True,
        "other_sessions_invalidated": False,
        "other_sessions_note": "THE SERVER WROTE THIS SENTENCE.",
    },
    "default_status": 200,
}

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

    def targets(self):
        """EVERY raw request target this server was given, query string
        included -- not just the /api/ ones, and not the parsed path.

        This is the evidence for "the password is not in a URL". A form whose
        submit is not prevented navigates to the CURRENT page with the field
        values in the query string, and the current page is `/studio`, not an
        API route. Recording only /api/ paths would have made that defect
        invisible to a test that looked for it."""
        return [r["target"] for r in self.seen]


def _respond(state, rec):
    path = rec["path"]
    if path == "/api/auth/me":
        return 200, dict(state["me"])
    if path == "/api/auth/password":
        return state["password_status"], dict(state["password_body"])
    if path == "/api/auth/logout":
        return 200, {"ok": True}
    if path == "/api/workspace":
        return 200, {"projects": [], "recent_runs": [],
                     "sections_available": 0, "drugs_available": 0,
                     "totals": {"projects": 0, "runs": 0, "specimens": 0},
                     "server": {"contract": "v1", "stub": False}}
    if path == "/api/sections":
        return 200, list(_SECTIONS)
    if path == "/api/defaults":
        return 200, dict(_DEFAULTS)
    if path == "/api/account":
        return 200, {"contract": "v1",
                     "auth": {"required": True, "header": "X-API-Key",
                              "scheme": "shared-secret", "matches": True,
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

    def _record(self):
        parsed = urllib.parse.urlparse(self.path)
        n = int(self.headers.get("content-length") or 0)
        rec = {
            "method": self.command,
            # THE RAW TARGET, query string included. See _Portal.targets.
            "target": self.path,
            "path": parsed.path,
            "query": parsed.query,
            "key": self.headers.get("x-api-key"),
            "cookie": self.headers.get("cookie"),
            "body": self.rfile.read(n).decode("utf-8") if n else "",
        }
        self.server.seen.append(rec)
        return rec

    def _api(self, rec):
        if not rec["path"].startswith("/api/"):
            return False
        status, payload = _respond(self.server.state, rec)
        raw = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("content-type", "application/json")
        self.send_header("content-length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)
        return True

    def do_GET(self):           # noqa: N802, D102
        rec = self._record()
        path = rec["path"]
        if path in ("/login", "/login/"):
            self.path = "/login.html"
        elif path in ("/studio", "/studio/"):
            self.path = "/index.html"
        elif path == "/":
            self.path = "/home.html"
        else:
            self.path = path
        if not self._api(rec):
            super().do_GET()

    def do_POST(self):          # noqa: N802, D102
        if not self._api(self._record()):
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


def _open_account_panel(browser, portal):
    """The Studio, with the account panel open. The form lives behind the
    EXISTING `account` control (`#btn-account`), beside logout."""
    ctx = browser.new_context()
    pg = ctx.new_page()
    errors: list = []
    pg.on("pageerror", lambda e: errors.append("pageerror: " + str(e)))
    pg.goto(f"{portal.base}/studio", wait_until="domcontentloaded")
    pg.wait_for_timeout(1400)
    page = _Page(pg, errors)
    assert not page.errors, page.errors[:2]
    assert pg.eval_on_selector_all("#btn-account", "e => e.length") == 1, (
        "no account control on the Studio")
    pg.click("#btn-account")
    pg.wait_for_timeout(900)
    return page, ctx


# --------------------------------------------------------------------------- #
# 7. The form, executed in a browser
# --------------------------------------------------------------------------- #
class TestTheChangePasswordForm:

    def test_it_is_behind_the_account_control_beside_logout(
            self, browser, portal):
        """Where the product owner asked for it: the existing `account`
        control, not a new menu beside it."""
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            assert not page.errors, page.errors[:2]
            dlg = page.pg.inner_text("#acct-dlg")
            assert "Sign out" in dlg, dlg[:400]
            n = page.pg.eval_on_selector_all(
                "#acct-dlg form.acct-form input[type=password]",
                "e => e.length")
            assert n >= 3, (
                f"the three password fields are not in the panel (found {n}); "
                f"the panel shows: {dlg[:400]!r}")
        finally:
            ctx.close()

    def test_it_asks_for_the_current_password_and_the_new_one_twice(
            self, browser, portal):
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            kinds = page.pg.eval_on_selector_all(
                "#acct-dlg input[type=password]",
                "els => els.map(e => [e.id, e.autocomplete])")
            ids = [k[0] for k in kinds]
            assert "pwd-current" in ids, kinds
            assert "pwd-new" in ids, kinds
            assert "pwd-again" in ids, kinds
            got = dict(kinds)
            # Correct autocomplete tokens, so a password manager offers the
            # right thing instead of pasting the old password into the new
            # field.
            assert got["pwd-current"] == "current-password", got
            assert got["pwd-new"] == "new-password", got
            assert got["pwd-again"] == "new-password", got
        finally:
            ctx.close()

    def test_no_password_input_carries_a_name_attribute(
            self, browser, portal):
        """Guard two of three against a password in a URL. A `name` is what a
        form submission serialises into a query string; with none, even an
        unprevented submit has nothing to put there."""
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            names = page.pg.eval_on_selector_all(
                "#acct-dlg input[type=password]",
                "els => els.map(e => e.getAttribute('name'))")
            assert names and all(n is None for n in names), names
        finally:
            ctx.close()

    def test_it_posts_the_route_with_the_body_the_code_defines(
            self, browser, portal):
        """The shape read off ``auth.py::auth_password``: a POST to
        /api/auth/password carrying old_password and new_password."""
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            page.pg.fill("#pwd-current", PASSWORD)
            page.pg.fill("#pwd-new", NEW_PASSWORD)
            page.pg.fill("#pwd-again", NEW_PASSWORD)
            page.pg.click("#acct-dlg form.acct-form button[type=submit]")
            page.pg.wait_for_timeout(900)

            posts = portal.calls("/api/auth/password", "POST")
            assert posts, [(r["method"], r["target"]) for r in portal.seen]
            body = json.loads(posts[-1]["body"])
            assert set(body) == {"old_password", "new_password"}, body
            assert body["old_password"] == PASSWORD
            assert body["new_password"] == NEW_PASSWORD
            assert not page.errors, page.errors[:2]
        finally:
            ctx.close()

    def test_the_new_password_never_appears_in_any_url(
            self, browser, portal):
        """THE PLANT: send the two fields as a query string (or let the <form>
        submit itself, which navigates to /studio with every value in the
        query). Either puts a live password in the address bar, in the
        browser's history, in this server's access log and in any proxy
        between them -- from code that looks entirely ordinary.

        Checked against EVERY request target the server saw, not just the API
        ones, because the form's own unprevented submit goes to /studio."""
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            page.pg.fill("#pwd-current", PASSWORD)
            page.pg.fill("#pwd-new", NEW_PASSWORD)
            page.pg.fill("#pwd-again", NEW_PASSWORD)
            page.pg.click("#acct-dlg form.acct-form button[type=submit]")
            page.pg.wait_for_timeout(1000)

            offenders = [t for t in portal.targets()
                         if NEW_PASSWORD in t or PASSWORD in t
                         or "password=" in t]
            assert not offenders, offenders
            # The address bar, which is where a reader would SEE it.
            assert NEW_PASSWORD not in page.pg.url, page.pg.url
            assert PASSWORD not in page.pg.url, page.pg.url
            # And the request that did carry it was a POST with a body.
            posts = portal.calls("/api/auth/password", "POST")
            assert posts and NEW_PASSWORD in posts[-1]["body"]
            assert posts[-1]["query"] == "", posts[-1]["target"]
        finally:
            ctx.close()

    def test_the_form_does_not_navigate(self, browser, portal):
        """Guard one of three: preventDefault() is the first statement of the
        submit handler. If it is gone the page leaves /studio, and everything
        the reader was looking at goes with it."""
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            before = page.pg.url
            page.pg.fill("#pwd-current", PASSWORD)
            page.pg.fill("#pwd-new", NEW_PASSWORD)
            page.pg.fill("#pwd-again", NEW_PASSWORD)
            page.pg.click("#acct-dlg form.acct-form button[type=submit]")
            page.pg.wait_for_timeout(900)
            assert page.pg.url == before, (before, page.pg.url)
        finally:
            ctx.close()

    def test_it_renders_the_servers_own_sentence_about_other_sessions(
            self, browser, portal):
        """The server ships ``other_sessions_note`` precisely so this surface
        does not re-type it. Reporting a password change as "signed out
        everywhere" would be the comfortable lie the route refuses to tell."""
        portal.reset(password_body={
            "ok": True, "other_sessions_invalidated": False,
            "other_sessions_note": "ONLY-THE-SERVER-SAYS-THIS-SENTENCE."})
        page, ctx = _open_account_panel(browser, portal)
        try:
            page.pg.fill("#pwd-current", PASSWORD)
            page.pg.fill("#pwd-new", NEW_PASSWORD)
            page.pg.fill("#pwd-again", NEW_PASSWORD)
            page.pg.click("#acct-dlg form.acct-form button[type=submit]")
            page.pg.wait_for_timeout(900)
            said = page.pg.inner_text("#acct-dlg")
            assert "ONLY-THE-SERVER-SAYS-THIS-SENTENCE." in said, said[-800:]
            assert "NOTHING WAS SIGNED OUT" in said.upper(), said[-800:]
        finally:
            ctx.close()

    def test_the_fields_are_cleared_on_success(self, browser, portal):
        """A filled password field left on screen is a password on screen."""
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            page.pg.fill("#pwd-current", PASSWORD)
            page.pg.fill("#pwd-new", NEW_PASSWORD)
            page.pg.fill("#pwd-again", NEW_PASSWORD)
            page.pg.click("#acct-dlg form.acct-form button[type=submit]")
            page.pg.wait_for_timeout(900)
            # Scoped to the form: the Access tab also holds a type=password
            # input for replacing the X-API-Key, which is not this form's and
            # must not be cleared by it.
            values = page.pg.eval_on_selector_all(
                "#acct-dlg form.acct-form input[type=password]",
                "els => els.map(e => e.value)")
            assert values == ["", "", ""], values
        finally:
            ctx.close()

    def test_a_refusal_shows_the_servers_own_words(self, browser, portal):
        """The route's refusal envelope is an object, which api.js::detailOf
        drops. If this reads "401" instead of the sentence, the form is going
        through the wrong reader."""
        portal.reset(password_status=401, password_body={
            "detail": {"message": "your current password is incorrect"}})
        page, ctx = _open_account_panel(browser, portal)
        try:
            page.pg.fill("#pwd-current", "wrong-but-long-enough")
            page.pg.fill("#pwd-new", NEW_PASSWORD)
            page.pg.fill("#pwd-again", NEW_PASSWORD)
            page.pg.click("#acct-dlg form.acct-form button[type=submit]")
            page.pg.wait_for_timeout(900)
            said = page.pg.inner_text("#acct-dlg")
            assert "your current password is incorrect" in said, said[-600:]
            assert not page.errors, page.errors[:2]
        finally:
            ctx.close()

    def test_a_mismatch_and_a_short_password_are_not_sent_at_all(
            self, browser, portal):
        """The three mistakes that need no server to diagnose. The server
        remains the authority -- this never decides a password is GOOD."""
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            page.pg.fill("#pwd-current", PASSWORD)
            page.pg.fill("#pwd-new", NEW_PASSWORD)
            page.pg.fill("#pwd-again", NEW_PASSWORD + "x")
            page.pg.click("#acct-dlg form.acct-form button[type=submit]")
            page.pg.wait_for_timeout(600)
            assert not portal.calls("/api/auth/password", "POST"), "sent anyway"
            assert "do not match" in page.pg.inner_text("#acct-dlg")

            page.pg.fill("#pwd-new", "short")
            page.pg.fill("#pwd-again", "short")
            page.pg.click("#acct-dlg form.acct-form button[type=submit]")
            page.pg.wait_for_timeout(600)
            assert not portal.calls("/api/auth/password", "POST"), "sent anyway"
            assert "12 characters" in page.pg.inner_text("#acct-dlg")
        finally:
            ctx.close()

    def test_the_form_sends_no_api_key(self, browser, portal):
        """§1 and the route's docstring: session only. The key identifies a
        program, and a program has no own password."""
        portal.reset()
        ctx = browser.new_context()
        pg = ctx.new_page()
        pg.add_init_script(
            "(() => { try { localStorage.setItem('tmesim_api_key',"
            " 'A-REAL-OPERATOR-KEY'); } catch (e) {} })()")
        try:
            pg.goto(f"{portal.base}/studio", wait_until="domcontentloaded")
            pg.wait_for_timeout(1400)
            pg.click("#btn-account")
            pg.wait_for_timeout(900)
            pg.fill("#pwd-current", PASSWORD)
            pg.fill("#pwd-new", NEW_PASSWORD)
            pg.fill("#pwd-again", NEW_PASSWORD)
            pg.click("#acct-dlg form.acct-form button[type=submit]")
            pg.wait_for_timeout(900)
            posts = portal.calls("/api/auth/password", "POST")
            assert posts, [(r["method"], r["target"]) for r in portal.seen]
            assert posts[-1]["key"] is None, posts[-1]["key"]
        finally:
            ctx.close()

    def test_no_form_is_offered_to_a_browser_with_no_session(
            self, browser, portal):
        """The route is session-only, so drawing the form for a key-only
        browser would offer a control the server will always refuse."""
        portal.reset(me={"authenticated": False})
        ctx = browser.new_context()
        pg = ctx.new_page()
        try:
            pg.goto(f"{portal.base}/studio", wait_until="domcontentloaded")
            pg.wait_for_timeout(1400)
            if pg.eval_on_selector_all("#btn-account", "e => e.length") != 1:
                pytest.skip("the Studio redirected a signed-out browser away")
            pg.click("#btn-account")
            pg.wait_for_timeout(900)
            assert pg.eval_on_selector_all(
                "#acct-dlg input#pwd-new", "e => e.length") == 0
        finally:
            ctx.close()


# --------------------------------------------------------------------------- #
# 8. §5 -- nothing here implies a private space
# --------------------------------------------------------------------------- #
class TestNothingImpliesAPrivateSpace:

    def test_static_the_forbidden_wording_is_absent(self):
        """§5: no "your projects", no owner column, no avatars. It breaks by
        OMISSION rather than by error, so the phrases are checked literally.

        COMMENTS ARE STRIPPED FIRST, and that is not a loophole -- §5 is about
        what the product SAYS to a reader, and both of these files discuss the
        banned wording at length in order to forbid it. A check that could not
        tell "no avatars" in a comment from an avatar on screen would make the
        only correct way to pass it to stop explaining the rule."""
        banned = ("your projects", "your workspace", "my projects",
                  "private workspace", "owner", "avatar", "gravatar",
                  "your data", "your runs", "your password")
        for path in (PASSWORD_JS, ACCOUNT_JS):
            raw = open(path, encoding="utf-8").read()
            code = re.sub(r"/\*.*?\*/", " ", raw, flags=re.S)
            code = "\n".join(re.sub(r"(^|\s)//.*$", "", line)
                             for line in code.splitlines()).lower()
            for phrase in banned:
                assert phrase not in code, f"{path} says {phrase!r}"

    def test_that_check_can_actually_fire(self):
        """Anti-vacuity. Stripping comments is one regex away from stripping
        everything, and a check over an empty string passes forever."""
        import textwrap
        sample = textwrap.dedent("""
            // no avatar and no owner column here
            /* and none here either */
            const label = 'Your projects';
        """)
        code = re.sub(r"/\*.*?\*/", " ", sample, flags=re.S)
        code = "\n".join(re.sub(r"(^|\s)//.*$", "", line)
                         for line in code.splitlines()).lower()
        assert "your projects" in code, code
        assert "avatar" not in code and "owner" not in code, code

    def test_the_form_says_a_password_changes_the_way_in_and_nothing_else(
            self, browser, portal):
        portal.reset()
        page, ctx = _open_account_panel(browser, portal)
        try:
            said = page.pg.inner_text("#acct-dlg").lower()
            assert "shared" in said
            assert "every account reaches the same one" in said or \
                   "every account that can sign in here sees the same one" in said
        finally:
            ctx.close()


# --------------------------------------------------------------------------- #
# 9. static guards
# --------------------------------------------------------------------------- #
class TestStaticGuards:

    def test_static_every_module_parses_as_an_es_module(self):
        """`node --check` alone silently accepts a duplicate `export const`;
        only `--input-type=module` reads the file as what the browser reads."""
        for path in (PASSWORD_JS, ACCOUNT_JS):
            with open(path, "rb") as fh:
                src = fh.read()
            out = subprocess.run(
                ["node", "--input-type=module", "--check"],
                input=src, capture_output=True)
            assert out.returncode == 0, (
                path + "\n" + out.stderr.decode()[:600])

    def test_static_no_h_selector_puts_the_id_after_the_classes(self):
        """`h('div.cls#id')` THROWS and takes the whole surrounding render with
        it; a panel that builds nothing looks exactly like one still loading.
        This file's own copy of the repo-wide guard, because the form is built
        inside a dialog nobody looks at until they click."""
        bad = re.compile(r"""h\(\s*['"][a-zA-Z0-9-]*\.[^'"#]*#[^'"]*['"]""")
        for path in (PASSWORD_JS, ACCOUNT_JS):
            for i, line in enumerate(open(path, encoding="utf-8"), 1):
                assert not bad.search(line), f"{path}:{i}: {line.strip()}"

    def test_static_the_secret_never_goes_into_a_url_in_the_source(self):
        """The two shapes that put a password in a query string, read
        syntactically: a template or a concatenation that appends a password
        variable to the path, and a `name` attribute on a password input."""
        text = open(PASSWORD_JS, encoding="utf-8").read()
        for shape in ("?old_password=", "?new_password=",
                      "old_password=' +", 'old_password=" +',
                      "new_password=' +", 'new_password=" +',
                      "${oldPassword}", "${newPassword}"):
            assert shape not in text, f"password.js builds a URL with {shape!r}"
        # preventDefault must be the first statement of the submit handler.
        handler = text.split("onsubmit:", 1)
        assert len(handler) == 2, "the submit handler moved"
        after = handler[1]
        assert after.index("event.preventDefault()") < after.index("run()"), (
            "run() is called before preventDefault(); a throw inside it would "
            "let the browser submit these fields as a GET query string")

    def test_static_the_request_field_names_match_the_route(self):
        """One source for the pairing. If someone renames a field in auth.py
        this fails HERE, naming the cause, rather than in a browser."""
        server = open(os.path.join(_REPO, "app", "server", "auth.py"),
                      encoding="utf-8").read()
        found = set(re.findall(r'_text_field\(body, "(\w+_password)"\)',
                               server))
        assert found == {"old_password", "new_password"}, found
        client = open(PASSWORD_JS, encoding="utf-8").read()
        for name in found:
            assert name in client, f"password.js does not send {name}"

    def test_static_the_cli_never_prints_a_hash_or_a_password(self):
        """Read syntactically: no print/format in app/server/__main__.py may
        reference a password or hash variable."""
        import ast
        import inspect
        tree = ast.parse(inspect.getsource(server_main))
        secretish = {"password", "new_password", "first", "second",
                     "password_hash", "hashed", "stored"}
        offenders = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if not (isinstance(func, ast.Name) and func.id == "print"):
                continue
            for arg in node.args:
                for inner in ast.walk(arg):
                    if isinstance(inner, ast.Name) and inner.id in secretish:
                        offenders.append(
                            f"line {node.lineno}: print() references "
                            f"{inner.id!r}")
        assert not offenders, offenders

    def test_static_the_account_store_is_gitignored(self):
        """§6's fallback path when TMESIM_STATE_DIR is unset -- the documented
        local workflow -- is app/auth/, inside the tree. That file holds every
        password hash AND the instance_secret the session signing key is
        derived from (§2), so a committed copy hands over both the hashes to
        attack offline and the key to forge a cookie for any account. Every
        sibling Studio store was ignored and this one was not."""
        out = subprocess.run(
            ["git", "check-ignore", "-q", "app/auth/users.json"],
            cwd=_REPO, capture_output=True)
        assert out.returncode == 0, (
            "app/auth/users.json is not gitignored; the account store would "
            "be staged by `git add -A`")
        tracked = subprocess.run(
            ["git", "ls-files", "app/auth"], cwd=_REPO,
            capture_output=True, text=True).stdout.strip()
        assert not tracked, f"the account store is already tracked: {tracked}"

    def test_static_no_email_path_was_introduced(self):
        """§7: no verification mail, no reset link, no SMTP. The way it arrives
        is one "helpful" import, and a password-reset feature is exactly where
        somebody reaches for it."""
        import ast
        for path in MY_FILES:
            if not path.endswith(".py"):
                continue
            tree = ast.parse(open(path, encoding="utf-8").read())
            imported = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported |= {a.name.split(".")[0] for a in node.names}
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module.split(".")[0])
            for banned in ("smtplib", "email", "mailbox", "ssl_email"):
                assert banned not in imported, f"{path} imports {banned}"
        for path in (PASSWORD_JS,):
            text = open(path, encoding="utf-8").read().lower()
            for banned in ("mailto:", "smtp", "reset link", "verification "
                           "email"):
                assert banned not in text, f"{path} mentions {banned!r}"

    def test_static_no_roles_were_introduced(self):
        """§7: a role model with one role is just a column that lies later. A
        disable/enable CLI is where an `is_admin` arrives."""
        store_src = open(os.path.join(_REPO, "app", "server", "auth.py"),
                         encoding="utf-8").read()
        for banned in ('"role"', '"roles"', '"permissions"', '"is_admin"',
                       '"scopes"', '"tier"'):
            assert banned not in store_src, f"auth.py stores {banned}"

    def test_static_the_contract_amendment_is_appended_not_rewritten(self):
        """The contract is FROZEN. An amendment is appended and dated; the
        frozen text above it is byte-identical, which the quoted facts in
        tests/test_auth_portal.py::test_the_contract_is_unedited also pin."""
        text = open(os.path.join(_REPO, "docs", "AUTH_CONTRACT.md"),
                    encoding="utf-8").read()
        assert "**Frozen 2026-10-09.**" in text
        assert "## 9. Amendment" in text, "the amendment is not there"
        head, _, amendment = text.partition("## 9. Amendment")
        # Everything §8 said is still above the amendment, unchanged.
        assert "## 7. What this contract does NOT authorise" in head
        assert "Deleting or weakening the `X-API-Key` path." in head
        assert "python -m app.server account" in amendment
        assert "old_password" in amendment and "new_password" in amendment
