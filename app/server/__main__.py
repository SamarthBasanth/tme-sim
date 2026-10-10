"""``./.venv/bin/python -m app.server`` — start the tmesim Studio server, and
the operator commands that are not reachable over HTTP.

TWO MODES, AND NO SUBCOMMAND STILL MEANS "SERVE"
------------------------------------------------
``python -m app.server [--host ...]``        start the server, exactly as before
``python -m app.server account <command>``   the operator's account commands

The subcommand is OPTIONAL and the server flags stayed on the top-level parser
for a reason that is not style: ``tmesim/cli.py::cmd_studio`` builds an argv of
``--host/--port/--log-level`` and calls :func:`main` directly, and
``docs/DEPLOYMENT.md``, the ``Dockerfile`` and ``render.yaml`` all invoke the
module with bare flags. Requiring a ``serve`` subcommand would have broken
every one of them to make a help text tidier.

WHY THE ACCOUNT COMMANDS ARE HERE AND NOT ON A ROUTE
----------------------------------------------------
``docs/AUTH_CONTRACT.md`` §7 forbids email of any kind -- no verification mail,
no reset link, no SMTP -- and says a forgotten password is reset by the
operator. It did not say HOW, and until now the only way was to hand-edit
``users.json``, which means hand-computing a password hash. These commands are
that missing mechanism.

They are a CLI and not an endpoint because the credential they hand out is a
password for an account that can do everything the API key can do (§7: no
roles). An HTTP route that resets somebody's password needs an authority above
every account to authorise it, and the only candidate this product has is
``TMESIM_STUDIO_API_KEY`` -- which §4 already spends as the sign-up bootstrap
token, is printed in a hosting dashboard, and is pasted into CI. Making that
value a master credential over identities is more than §4 gave it. The
filesystem already expresses the authority correctly: ``users.json`` is
``0600`` under a ``0700`` directory (§6), so "can edit the account store" and
"is the operator on this host" are the same sentence. These commands inherit
that and add nothing.

So there is no route, no port, no token to steal, and nothing here is
importable by ``api.py``. Asserted by ``tests/test_account_recovery.py``.

Binds to 127.0.0.1 by default: this is a local research workstation application,
not a network service, and the data it exposes is a reconstruction of a
de-identified tissue section. It performs no diagnosis and is not built for
identifiable patient data. Pass ``--host`` explicitly if you really mean to
expose it -- doing so with no ``--api-key`` and no ``TMESIM_STUDIO_API_KEY`` set
auto-generates one and prints it once, rather than serving unauthenticated on a
non-loopback address. See docs/DEPLOYMENT.md (§0 for the de-identified-tissue
scope).

That auto-generation is keyed on the BIND ADDRESS, which is not the same thing
as who can reach the process: a loopback bind fronted by an ngrok/Cloudflare
tunnel or a reverse proxy is world-reachable while this process still sees
``127.0.0.1`` and generates nothing. The default is deliberately left alone --
unauthenticated loopback is the documented single-user workflow -- but every
start with no key now prints ``_no_auth_warning()`` so the operator cannot
reach that situation without having been told.
"""
from __future__ import annotations

import argparse
import getpass
import os
import secrets
import stat
import sys
import time

_LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1"}

#: Exit codes. 0/1 are the ordinary pair; 78 is ``EX_CONFIG`` from
#: ``sysexits.h`` and is already what this module returns for an unusable
#: ``TMESIM_STATE_DIR``, so a misconfigured store answers with the same number
#: as a misconfigured state directory rather than a second convention.
EXIT_OK = 0
EXIT_REFUSED = 1
EXIT_CONFIG = 78


def _ensure_repo_on_path() -> str:
    """Put the repository root on ``sys.path`` so ``import tmesim`` works from
    any working directory."""
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    if root not in sys.path:
        sys.path.insert(0, root)
    return root


def _no_auth_warning(host: str) -> str:
    """The startup banner block printed whenever no API key is in force.

    SECURITY: the auto-generate rule above reasons about the bind address, and
    the bind address is not reachability. Tunnelling a loopback service
    (``ngrok http 8000``, ``cloudflared``, an SSH ``-R`` forward, a reverse
    proxy with no auth of its own) is the ordinary way people demo software,
    and it is exactly the case where the two differ: ``is_loopback`` is true,
    no key is generated, and the server is on the internet unauthenticated.
    Nothing in the process can detect that, so the only honest fix is to say
    it out loud on every unauthenticated start rather than to guess.

    Kept as a function so the wording is asserted once, in
    tests/test_app_server.py, instead of being duplicated per call site.
    """
    return "\n".join((
        "",
        "  !!  WARNING: NO API KEY -- every /api/ route is UNAUTHENTICATED.",
        "  !!",
        f"  !!  Binding {host} is NOT protection. This process knows only what",
        "  !!  it bound to, not who can reach it. Put an ngrok/Cloudflare",
        "  !!  tunnel, an SSH -R forward or a reverse proxy in front of a",
        "  !!  loopback bind and it is reachable from the internet -- while",
        "  !!  this process still sees 127.0.0.1 and auto-generates no key,",
        "  !!  because that rule reads the bind address, not reachability.",
        "  !!  If anything at all forwards traffic here, you are exposed.",
        "  !!",
        "  !!  Unauthenticated, POST /api/sections/scan takes an arbitrary",
        "  !!  filesystem path and reports on any directory this process can",
        "  !!  read -- that is filesystem reconnaissance, not just data access.",
        "  !!",
        "  !!  Fix: restart with TMESIM_STUDIO_API_KEY=<a-real-secret> set, or",
        "  !!  pass --api-key <a-real-secret>. See docs/DEPLOYMENT.md §1e.",
        "  !!  A key authenticates /api/sections/scan; it does not confine the",
        "  !!  path it will read. Only give it to people you'd give a shell to.",
        "",
    ))


# --------------------------------------------------------------------------- #
# `account` — the operator commands (docs/AUTH_CONTRACT.md §7)
# --------------------------------------------------------------------------- #
def _auth_module():
    """``app.server.auth``, or ``None`` if the optional extra is missing.

    Imported HERE rather than at module top because ``auth.py`` imports
    fastapi: a checkout without ``pip install -e '.[app]'`` must still be able
    to run ``python -m app.server --help`` and get a usage message instead of
    an ImportError traceback. Same trade, and the same wording,
    ``tmesim/cli.py::cmd_studio`` makes for the server itself.
    """
    _ensure_repo_on_path()
    try:
        from . import auth as auth_mod
    except ImportError:
        # Not bound, and not reported in detail. The one actionable fact is
        # which extra is missing, and the caller prints that.
        return None
    return auth_mod


def _store_location_lines(store) -> list:
    """Where this command is about to write, and the one way that goes wrong.

    PRINTED BY EVERY ``account`` COMMAND, INCLUDING THE READ-ONLY ONE, because
    the likeliest operational mistake with this CLI is not a dangerous flag --
    it is editing the wrong store and believing it worked. §6 puts the store
    under ``TMESIM_STATE_DIR`` and falls back to ``app/auth/`` when that is
    unset, so an operator who starts the server from systemd or Docker (with
    the variable set) and then runs this command from an interactive shell
    (without it) resets a password in a file no server will ever read. The
    command succeeds, the operator tells someone their new password, and the
    sign-in refuses it. Naming the path on every invocation is cheap and is
    the only thing that makes that visible.
    """
    root = os.environ.get("TMESIM_STATE_DIR", "").strip()
    lines = [f"  store       {store.path}"]
    if root:
        lines.append(f"  state dir   TMESIM_STATE_DIR={root}")
    else:
        lines.append("  state dir   TMESIM_STATE_DIR is UNSET, so this is the "
                     "in-repo fallback (app/auth/).")
        lines.append("              IF THE SERVER RUNS WITH THAT VARIABLE SET "
                     "IT IS USING A DIFFERENT STORE, and this command will "
                     "not affect it.")
        lines.append("              Export the same value you start the "
                     "server with, then run this again.")
    return lines


def _private_store_problem(path: str):
    """Why this process must not write a password hash here, or ``None``.

    THE AUTHORITY FOR THESE COMMANDS IS FILESYSTEM OWNERSHIP (see the module
    docstring), so it is worth checking that the filesystem actually says what
    it is assumed to say rather than trusting it:

    * **owned by somebody else** -- then this is not the operator's store and
      the right answer is to stop, not to write into another user's file and
      possibly leave it owned by them with our contents in it.
    * **readable or writable by group or other** -- then a fresh hash is about
      to be written into a file other users on this host can read offline at
      their leisure, and the session signing key sits in the same file (§2).
      ``UserStore._write`` does re-chmod to ``0600``, so this is not a
      permanent state -- but a store that is world-readable NOW has already
      been exposed, and silently tightening it would hide that from the one
      person who needs to know.

    Both refuse rather than warn, for writes. A missing file is not a problem:
    no store yet is the ordinary first-run state (§4), and ``UserStore``
    creates it ``0600`` under a ``0700`` directory.
    """
    try:
        info = os.stat(path)
    except FileNotFoundError:
        return None
    except OSError:
        return (f"the account store at {path} exists but cannot be examined. "
                f"Check the permissions on it and on its directory.")
    if hasattr(os, "geteuid") and info.st_uid != os.geteuid():
        return (f"the account store at {path} is owned by uid "
                f"{info.st_uid}, and this process is uid {os.geteuid()}. "
                f"These commands are the operator's, and the operator is "
                f"whoever owns the store -- run this as that user (the same "
                f"one the server runs as) rather than writing into their "
                f"file.")
    leaked = stat.S_IMODE(info.st_mode) & 0o077
    if leaked:
        return (f"the account store at {path} is mode "
                f"{stat.S_IMODE(info.st_mode):04o}, which lets other users on "
                f"this host read it. It holds every password hash and the "
                f"session signing key (docs/AUTH_CONTRACT.md §2, §6), so it "
                f"must be 0600 before anything new is written into it. Fix "
                f"it with:  chmod 600 {path}  (and chmod 700 on its "
                f"directory), then run this again. Refusing rather than "
                f"quietly tightening it, because a store that is readable now "
                f"has already been readable for as long as it has existed.")
    return None


def _stamp(value) -> str:
    """A timestamp an operator can read, or ``never``."""
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return "never"
    if value <= 0:
        return "never"
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(value))


def _read_new_password(from_stdin: bool, email: str):
    """``(password, why_not)``. Exactly one of the two is meaningful.

    THERE IS NO ``--password`` FLAG AND THERE WILL NOT BE ONE. A password on
    an argv is readable by every user on the host for as long as the process
    lives (``ps -ef`` shows the full command line on Linux and macOS), it is
    appended to the shell's history file where it outlives the reset by
    months, and in any orchestration tool it is printed into the job log. None
    of those is recoverable by changing the password later, because the value
    is already written down somewhere nobody is watching.

    So: a no-echo prompt, or an explicit pipe.

    * **prompt** -- :func:`getpass.getpass` writes its prompt to the terminal
      device (``/dev/tty``) and not to stdout, and turns echo off. Asked
      twice, because a typo in a value nobody can see locks out the account
      this command exists to recover.
    * ``--stdin`` -- one line, for a password manager or a deploy tool. The
      trailing newline is removed AND NOTHING ELSE: a password may legitimately
      start or end with a space, and ``.strip()`` would set a different
      password from the one that was piped in, which would then fail at
      sign-in with no explanation.

    THE FALLBACK IS REFUSED, WHICH IS THE NON-OBVIOUS PART. When
    :mod:`getpass` cannot control the terminal it does not fail -- it falls
    back to :func:`input` and warns, which means it ECHOES the password onto
    the screen and into the scrollback. That is the exact outcome this
    function exists to prevent, so the warning is promoted to an error here
    and the command refuses instead.
    """
    if from_stdin:
        if sys.stdin.isatty():
            return None, ("--stdin was passed but stdin is a terminal, so "
                          "anything typed would be echoed to the screen and "
                          "kept in the scrollback. Drop --stdin to be "
                          "prompted without echo, or pipe the password in.")
        line = sys.stdin.readline()
        if not line:
            return None, ("--stdin was passed and stdin was empty. Nothing "
                          "was written.")
        return (line[:-1] if line.endswith("\n") else line), ""
    if not sys.stdin.isatty():
        return None, ("there is no terminal here to prompt on without echo. "
                      "Pipe the password in and pass --stdin, e.g.  "
                      "printf '%s' \"$NEW\" | python -m app.server account "
                      "reset-password --email ... --stdin")
    import warnings
    try:
        with warnings.catch_warnings():
            # Promote getpass's "could not turn off echo" warning to an
            # exception. Left as a warning it is a line of noise above a
            # password that is now on screen.
            warnings.simplefilter("error", getpass.GetPassWarning)
            first = getpass.getpass(f"New password for {email} (not echoed): ")
            second = getpass.getpass("Repeat it: ")
    except getpass.GetPassWarning:
        return None, ("this terminal will not turn echo off, so the password "
                      "would have been printed as it was typed. Refusing. "
                      "Pipe it in with --stdin instead.")
    except (EOFError, KeyboardInterrupt):
        return None, "cancelled. Nothing was written."
    if first != second:
        # A plain `!=`, and deliberately so: both halves were typed by the
        # person running this command one second apart, neither is being
        # compared against a stored secret, and there is no remote attacker to
        # time. compare_digest here would be cargo cult -- §2's rule is about
        # comparing a PRESENTED secret with a STORED one, which this is not.
        return None, "the two entries did not match. Nothing was written."
    return first, ""


def _reset_consequences(email: str) -> str:
    """What a reset ends and what it does not. One copy, so a test can pin it.

    §2 says there is no per-session revocation, and this is the command where
    that fact is most likely to be assumed away: every other product an
    operator has used signs the other browsers out on a password reset. This
    one cannot, and saying so is the difference between an operator who
    disables the account during an incident and one who resets the password,
    believes the intruder is out, and goes to bed.
    """
    return "\n".join((
        f"  Reset the password for {email}.",
        "",
        "  WHAT THIS DID",
        "    The new password works on the next sign-in. The old one does not.",
        "",
        "  WHAT THIS DID NOT DO -- NOBODY WAS SIGNED OUT",
        "    A session here is a signed token with no server-side record "
        "(docs/AUTH_CONTRACT.md §2),",
        "    so there is nothing to revoke and this product has no "
        "per-session revocation at all.",
        "    Every vivum_session cookie already issued for this account keeps "
        "opening the gate",
        "    until its own expiry, up to 14 days from when it was issued. If "
        "the reason for this",
        "    reset is that somebody else had the password, THEY ARE STILL "
        "SIGNED IN.",
        "",
        "  THE TWO THINGS THAT DO END A LIVE SESSION",
        "    python -m app.server account disable --email <address>",
        "        Ends this one account, at once: `disabled` is re-read on "
        "every request, so the",
        "        next one a held cookie makes is refused. `enable` lets those "
        "same cookies back",
        "        in if they have not expired, so it is not a clean undo.",
        "    python -m app.server account rotate-session-key --yes   (then "
        "restart the server)",
        "        Ends EVERY session on this deployment, every account. There "
        "is no narrower",
        "        version of this, which is the honest shape of a stateless "
        "session.",
        "",
        "  UNTOUCHED",
        "    X-API-Key. A program holding the key never used this password "
        "and still does not",
        "    need one (§0), so nothing in CI or in a documented curl flow "
        "changed here.",
    ))


def _disable_consequences(email: str, disabled: bool) -> str:
    """What ``disable``/``enable`` actually changes. One copy; see above."""
    if disabled:
        return "\n".join((
            f"  Disabled {email}.",
            "",
            "    In effect NOW, including for browsers already signed in: the "
            "account row is",
            "    re-read on every request, so the next request any held "
            "cookie makes is refused.",
            "    This is the only immediate lever this product has over one "
            "account.",
            "    Sign-in is refused with the same words an unknown address "
            "gets -- the server will",
            "    not tell a caller that an address is disabled here, because "
            "that is the same",
            "    disclosure as confirming the account exists.",
            "    The password is unchanged. `enable` restores access, "
            "INCLUDING to any cookie that",
            "    has not yet expired, so if a session is what you are shutting "
            "out, reset the",
            "    password before enabling it again -- or rotate the session "
            "key.",
            "    X-API-Key is untouched (§0): disabling an account does not "
            "close the key path.",
        ))
    return "\n".join((
        f"  Enabled {email}.",
        "",
        "    Sign-in works again with the account's existing password, which "
        "was never changed.",
        "    ALSO BACK IN: any vivum_session cookie issued before the account "
        "was disabled, if it",
        "    has not reached its own 14-day expiry. Disabling and re-enabling "
        "is therefore not a",
        "    way to force a re-login. If that is what you want, reset the "
        "password (which still",
        "    ends no session) or rotate the session key and restart.",
    ))


def _account_table(rows: list) -> list:
    """``account list``'s table, as lines. No hash can appear in it: the rows
    are ``operator_user`` projections, which do not carry one."""
    head = ("id", "address", "created", "last sign-in", "state")
    body = [(str(r.get("id") or "?"), str(r.get("email") or "?"),
             _stamp(r.get("created_at")), _stamp(r.get("last_login_at")),
             "DISABLED" if r.get("disabled") else "active") for r in rows]
    widths = [max(len(head[i]), *(len(b[i]) for b in body))
              for i in range(len(head))]
    out = ["  " + "  ".join(h.ljust(widths[i]) for i, h in enumerate(head)),
           "  " + "  ".join("-" * w for w in widths)]
    for b in body:
        out.append("  " + "  ".join(c.ljust(widths[i])
                                    for i, c in enumerate(b)))
    return out


def cmd_account_list(args, auth_mod) -> int:
    """Every account. Reads; writes nothing."""
    store = auth_mod.UserStore()
    for line in _store_location_lines(store):
        print(line)
    # A WARNING here and a refusal for the writing commands. This command
    # prints no secret, so refusing it would make the store's own permissions
    # harder to diagnose -- which is the opposite of useful when the message
    # is about permissions.
    problem = _private_store_problem(store.path)
    if problem is not None:
        print(f"  WARNING: {problem}", file=sys.stderr)
    rows = auth_mod.UserStore().list_accounts()
    print(f"  accounts    {len(rows)}")
    print("")
    if not rows:
        print("  No accounts in this store.")
        print("  The first one is claimed through POST /api/auth/signup with "
              "this deployment's")
        print("  TMESIM_STUDIO_API_KEY as the bootstrap token (§4) -- these "
              "commands deliberately")
        print("  do not create accounts, so there is exactly one way in and "
              "§4 governs it.")
        return EXIT_OK
    for line in _account_table(rows):
        print(line)
    return EXIT_OK


def cmd_account_reset_password(args, auth_mod) -> int:
    """§7's operator reset, which had no mechanism until now."""
    store = auth_mod.UserStore()
    for line in _store_location_lines(store):
        print(line)
    problem = _private_store_problem(store.path)
    if problem is not None:
        print(f"REFUSING: {problem}", file=sys.stderr)
        return EXIT_CONFIG

    email = auth_mod.normalise_email(args.email)
    row = store.by_email(email)
    if row is None:
        # NOT the login route's deliberately ambiguous refusal. That one
        # exists so an anonymous caller cannot ask who has an account here;
        # this caller already owns the store and can read every address in it
        # with `account list`, so withholding the answer would protect nothing
        # and would make a typo in an address indistinguishable from a bug.
        print(f"REFUSED: there is no account for {email} in this store. "
              f"`python -m app.server account list` prints every address. "
              f"Nothing was written.", file=sys.stderr)
        return EXIT_REFUSED

    password, why = _read_new_password(bool(args.stdin), email)
    if password is None:
        print(f"REFUSED: {why}", file=sys.stderr)
        return EXIT_REFUSED
    problem = auth_mod.password_problem(password)
    if problem is not None:
        # The SAME rule and the same sentence the sign-up and change-password
        # routes use (§3), from the same function -- a CLI with its own floor
        # is a CLI that can mint a password the product refuses to accept.
        print(f"REFUSED: {problem}. Nothing was written.", file=sys.stderr)
        return EXIT_REFUSED

    # hash_password, not a hand-rolled one: on this interpreter scrypt does
    # not exist (macOS links LibreSSL) and it returns a pbkdf2-sha256 string
    # which the running server verifies, because the algorithm names itself in
    # the stored string. A CLI that assumed scrypt would write a hash this
    # interpreter cannot verify, and the account would be locked out BY THE
    # RECOVERY TOOL.
    store.set_password(row["id"], auth_mod.hash_password(password))
    # Nothing below prints, logs or echoes the password or the hash. The only
    # output is what changed.
    print("")
    print(_reset_consequences(email))
    return EXIT_OK


def _set_disabled(args, auth_mod, disabled: bool) -> int:
    store = auth_mod.UserStore()
    for line in _store_location_lines(store):
        print(line)
    problem = _private_store_problem(store.path)
    if problem is not None:
        print(f"REFUSING: {problem}", file=sys.stderr)
        return EXIT_CONFIG
    email = auth_mod.normalise_email(args.email)
    row = store.by_email(email)
    if row is None:
        print(f"REFUSED: there is no account for {email} in this store. "
              f"Nothing was written.", file=sys.stderr)
        return EXIT_REFUSED
    if bool(row.get("disabled")) == disabled:
        print("")
        print(f"  {email} is already "
              f"{'disabled' if disabled else 'enabled'}. Nothing was written.")
        return EXIT_OK
    store.set_disabled(row["id"], disabled)
    print("")
    print(_disable_consequences(email, disabled))
    return EXIT_OK


def cmd_account_disable(args, auth_mod) -> int:
    return _set_disabled(args, auth_mod, True)


def cmd_account_enable(args, auth_mod) -> int:
    return _set_disabled(args, auth_mod, False)


def cmd_account_rotate_session_key(args, auth_mod) -> int:
    """Invalidate every session on this deployment (§2's only revocation)."""
    if not args.yes:
        print("REFUSED: this signs out every account on this deployment, not "
              "one. There is no per-session and no per-account version of it "
              "-- §2 sessions are signed tokens with no server-side record, "
              "so the only thing that can be invalidated is the key they were "
              "all signed with. Pass --yes if that is what you mean.",
              file=sys.stderr)
        return EXIT_REFUSED

    configured = os.environ.get("TMESIM_SESSION_SECRET", "").strip()
    if configured:
        # Refusing rather than rotating-and-warning. Rotating would succeed,
        # print a success, and end nothing at all, because §2 says
        # TMESIM_SESSION_SECRET wins and the store's instance_secret is then
        # never consulted. A command that cannot do what its name says must
        # not report that it did.
        print("REFUSED: TMESIM_SESSION_SECRET is set in this environment, so "
              "the session signing key comes from that variable and the "
              "store's instance_secret is never consulted (§2). Rotating it "
              "here would end no session and would print a success. To "
              "invalidate every session on a deployment configured that way, "
              "set TMESIM_SESSION_SECRET to a new value (at least 16 "
              "characters) and restart the server.", file=sys.stderr)
        return EXIT_REFUSED

    store = auth_mod.UserStore()
    for line in _store_location_lines(store):
        print(line)
    problem = _private_store_problem(store.path)
    if problem is not None:
        print(f"REFUSING: {problem}", file=sys.stderr)
        return EXIT_CONFIG
    store.rotate_instance_secret()
    print("")
    print("  Rotated this store's instance_secret.")
    print("")
    print("    NOT IN EFFECT YET. The running server cached the old signing "
          "key at first use")
    print("    (§2, deliberately -- the same property TMESIM_STUDIO_API_KEY "
          "has), so RESTART THE")
    print("    SERVER. Until you do, every existing session keeps working and "
          "nothing has changed.")
    print("    After the restart every vivum_session cookie on this "
          "deployment is invalid and")
    print("    everyone signs in again. Passwords and accounts are untouched: "
          "this rotated a")
    print("    signing key, not a credential.")
    print("    X-API-Key is untouched too (§0) -- scripts and CI never had a "
          "session.")
    return EXIT_OK


def _account_main(args) -> int:
    """Dispatch one ``account`` subcommand, and turn a refusal into an exit code."""
    auth_mod = _auth_module()
    if auth_mod is None:
        print("[account] the account store needs the optional 'app' extra "
              "(fastapi), which is not installed.", file=sys.stderr)
        print("[account] install it with:  pip install -e '.[app]'",
              file=sys.stderr)
        return EXIT_CONFIG

    func = getattr(args, "account_func", None)
    if func is None:
        # `python -m app.server account` with no subcommand. argparse leaves
        # the sub-namespace empty rather than erroring, because an optional
        # subparser cannot be `required` without also breaking the bare
        # `python -m app.server` that starts the server.
        print("REFUSED: `account` needs a subcommand. They are: list, "
              "reset-password, disable, enable, rotate-session-key. Try "
              "`python -m app.server account --help`.", file=sys.stderr)
        return EXIT_REFUSED

    from .cache import state_dir_problem
    # The same check, and the same exit code, that starting the server makes.
    # An unusable TMESIM_STATE_DIR is exactly as fatal here: auth_dir() would
    # otherwise raise a PermissionError traceback out of os.makedirs with the
    # real cause (the running uid is not the mount's owner) named nowhere.
    problem = state_dir_problem()
    if problem is not None:
        print(f"REFUSING: {problem}", file=sys.stderr)
        return EXIT_CONFIG

    try:
        return int(func(args, auth_mod))
    except auth_mod.AuthError as refusal:
        # The sentence is pulled onto a LOCAL before it is formatted, which is
        # the shape `auth.py::_refusal_response` documents and for the same
        # reason: tests/test_response_redaction.py flags any stringification
        # of a caught exception name, and that sweep covers every module under
        # app/server/ including this one. `refusal.message` is a sentence
        # written for a human by the code that raised it -- the store's own
        # refusals -- and is not exception text.
        said = refusal.message
        print(f"REFUSED: {said}", file=sys.stderr)
        return EXIT_REFUSED
    except OSError:
        # Not stringified: a filesystem error here is about a path the
        # operator supplied via TMESIM_STATE_DIR, and the actionable fact is
        # which store it was, which is already on screen.
        print("REFUSED: the account store could not be written. Check the "
              "permissions on the path printed above, and that it is on a "
              "writable filesystem.", file=sys.stderr)
        return EXIT_CONFIG


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m app.server",
        description="Vivum Studio backend (API_CONTRACT.md v1).",
    )
    p.add_argument("--host", default="127.0.0.1",
                   help="bind address (default 127.0.0.1, local only)")
    # PORT is the convention every container host uses to tell a process where
    # to listen -- Render, Heroku, Cloud Run and Fly all set it and route to it.
    # Ignoring it does not produce an error: the host routes to the port IT
    # assigned, nothing is listening there, and every request HANGS until it
    # times out. DNS resolves, the build is green, the dashboard looks alive,
    # and the service answers nothing. That is what a hosted deploy of this
    # image did before 2026-08-22.
    #
    # Precedence is explicit > host-assigned > local default, so `--port 9000`
    # still wins and an unset PORT still gives 8000 exactly as before.
    p.add_argument("--port", type=int,
                   default=int(os.environ.get("PORT") or 8000),
                   help="bind port (default: $PORT if the host sets one, else "
                        "8000)")
    p.add_argument("--reload", action="store_true",
                   help="uvicorn autoreload (development)")
    p.add_argument("--log-level", default="info",
                   choices=("critical", "error", "warning", "info", "debug"))
    p.add_argument("--cache-size", type=int, default=None,
                   help="max cached reconstructions (default 4; each holds a "
                        "full cell population and several float64 grids)")
    p.add_argument("--api-key", default=None,
                   help="require this value in every request's X-API-Key "
                        "header on /api/* routes. If omitted and --host is "
                        "non-loopback, one is generated and printed once. "
                        "Pass --api-key '' to explicitly disable (not "
                        "recommended off loopback).")

    # OPTIONAL, and the default (no subcommand at all) is still "serve" --
    # see the module docstring. argparse leaves `args.command` as None, which
    # is what main() dispatches on.
    sub = p.add_subparsers(dest="command", metavar="{account}")

    acct = sub.add_parser(
        "account",
        help="operator account commands (not reachable over HTTP)",
        description="Operator account commands for the login portal "
                    "(docs/AUTH_CONTRACT.md §7). Runs as the user that owns "
                    "the account store; there is no route that does any of "
                    "this.")
    # metavar lists them rather than letting argparse print the tuple, so
    # `reset-password` is visible in the one-line usage an operator reads
    # when they have forgotten the name of the command they need.
    asub = acct.add_subparsers(
        dest="account_command",
        metavar="{list,reset-password,disable,enable,rotate-session-key}")

    a_list = asub.add_parser(
        "list", help="every account: id, address, created, last sign-in, "
                     "whether it is disabled")
    a_list.set_defaults(account_func=cmd_account_list)

    a_reset = asub.add_parser(
        "reset-password",
        help="set a new password for one account (§7's operator reset)",
        description="Set a new password for one account. The password is "
                    "read from a no-echo prompt, or from stdin with "
                    "--stdin. THERE IS DELIBERATELY NO --password FLAG: "
                    "anything on a command line is visible in `ps` to every "
                    "user on the host and is written to the shell's history "
                    "file.")
    a_reset.add_argument("--email", required=True,
                         help="the account to reset, as it appears in "
                              "`account list`")
    a_reset.add_argument(
        "--stdin", action="store_true",
        help="read the new password from stdin (one line, the trailing "
             "newline stripped) instead of prompting. For a password manager "
             "or an orchestration tool: `pass show x | python -m app.server "
             "account reset-password --email a@b.co --stdin`. Refused when "
             "stdin is a terminal, because typing into it would echo.")
    a_reset.set_defaults(account_func=cmd_account_reset_password)

    a_dis = asub.add_parser(
        "disable",
        help="turn an account off. The ONLY thing here that ends a live "
             "session, because `disabled` is re-read on every request")
    a_dis.add_argument("--email", required=True)
    a_dis.set_defaults(account_func=cmd_account_disable)

    a_en = asub.add_parser(
        "enable",
        help="turn an account back on. Note that this also lets any "
             "unexpired cookie it already had back in")
    a_en.add_argument("--email", required=True)
    a_en.set_defaults(account_func=cmd_account_enable)

    a_rot = asub.add_parser(
        "rotate-session-key",
        help="invalidate EVERY session on this deployment (all accounts). "
             "Takes effect on the next server start",
        description="Replace the account store's `instance_secret`, which is "
                    "what the session signing key is derived from "
                    "(docs/AUTH_CONTRACT.md §2). This is the only way to end "
                    "a session that has already been issued, because §2 "
                    "sessions carry no server-side record. It ends them for "
                    "every account at once, not one -- there is no "
                    "per-session revocation in this product.")
    a_rot.add_argument(
        "--yes", action="store_true",
        help="required. Confirms that every signed-in browser on this "
             "deployment is to be signed out.")
    a_rot.set_defaults(account_func=cmd_account_rotate_session_key)

    return p


def main(argv=None) -> int:
    """Parse, then either serve or run an operator command.

    NO SUBCOMMAND MEANS SERVE, which keeps every existing caller working:
    ``tmesim/cli.py::cmd_studio`` calls this with an argv of bare flags, and
    so do the Dockerfile, render.yaml and both documented deploy flows.
    """
    args = build_parser().parse_args(argv)
    if getattr(args, "command", None) == "account":
        return _account_main(args)
    return _serve(args)


def _serve(args) -> int:
    """Start the server. Byte-for-byte what ``main()`` has always done."""
    root = _ensure_repo_on_path()

    if args.cache_size is not None:
        os.environ["TMESIM_RECON_CACHE_SIZE"] = str(args.cache_size)

    is_loopback = args.host in _LOOPBACK_HOSTS
    api_key = args.api_key
    generated = False
    if api_key is None:
        # --api-key wasn't passed. Respect an already-configured environment
        # (Docker ENV, systemd, a shell export) rather than silently clobber
        # it -- only decide a value ourselves when the environment has none.
        api_key = os.environ.get("TMESIM_STUDIO_API_KEY", "")
        if not api_key and not is_loopback:
            api_key = secrets.token_urlsafe(24)
            generated = True
    os.environ["TMESIM_STUDIO_API_KEY"] = api_key
    # Record WHICH of the three sources won, for /api/account to report.
    #
    # The banner below already works this out to print a label, but the banner
    # is scrollback: an operator reading Studio's Account panel three weeks
    # later cannot see it. Without this, account.py would have to infer the
    # source from the key's presence alone and could not distinguish "the
    # operator chose this secret" from "the process invented one at boot and
    # printed it exactly once" -- which is precisely the distinction someone
    # who has LOST the key needs. Unset means the variable was already in the
    # environment when the app was built (Docker ENV, systemd, a shell export)
    # and this module never ran.
    if api_key:
        os.environ["TMESIM_STUDIO_API_KEY_SOURCE"] = (
            "generated" if generated
            else "flag" if args.api_key is not None
            else "environment")

    # SECURITY: on the default loopback bind the only legitimate Host values
    # are the loopback names, so pin them -- that blocks DNS rebinding, where
    # a hostile domain resolves to 127.0.0.1 and thereby becomes same-origin
    # with this server (CORS cannot help there). Left unset for a non-loopback
    # bind, where we do not know the operator's hostname and refusing it would
    # break the deployment; set TMESIM_STUDIO_ALLOWED_HOSTS yourself there.
    #
    # IPv6 IS EXCLUDED, and not as an oversight. Starlette's TrustedHostMiddleware
    # matches on ``headers["host"].split(":")[0]`` to strip the port. That is
    # correct for a name or an IPv4 literal and wrong for an IPv6 one: a browser
    # sends ``Host: [::1]:8000``, which splits to ``[``. No allowlist entry can
    # ever equal ``[``, so pinning here would make ``--host ::1`` reject 100 % of
    # its own traffic with "400 Invalid host header" -- a supported flag broken
    # by its own protection. Enabling the middleware on a bind it cannot express
    # buys no safety and costs the whole server, so IPv6 loopback keeps the
    # non-loopback treatment: unset, and the operator may set it explicitly.
    _ipv6_bind = ":" in args.host
    if is_loopback and not _ipv6_bind and "TMESIM_STUDIO_ALLOWED_HOSTS" not in os.environ:
        os.environ["TMESIM_STUDIO_ALLOWED_HOSTS"] = "localhost,127.0.0.1"

    import uvicorn

    from .cache import (DATA_ROOTS_VAR, RECON_CACHE, runs_dir,
                        state_dir_problem)

    # SECURITY: same trade as the API key immediately above, for the same
    # reason and by the same rule. POST /api/sections/scan resolves and reads
    # an operator-supplied path; unset, it is UNCONFINED, which is the right
    # default on loopback (bring-your-own-tissue points at directories outside
    # this repo, and a data_cache/ default would break every existing local
    # workflow) and the wrong one the moment the process is reachable.
    #
    # So: off loopback, and ONLY when the operator has expressed no preference,
    # confine to the repo's own data_cache. Presence is tested, not truthiness
    # -- an explicitly empty value means "confine to nothing", which
    # confine_scan_path treats as fail-closed, and that is a deliberate choice
    # this must not overwrite.
    roots_generated = False
    if not is_loopback and DATA_ROOTS_VAR not in os.environ:
        os.environ[DATA_ROOTS_VAR] = os.path.join(root, "data_cache")
        roots_generated = True

    # BEFORE the banner, because the banner is what used to crash. The first
    # use of TMESIM_STATE_DIR is `runs_dir()` three lines down, and an
    # unwritable state directory turned that into a PermissionError traceback
    # out of os.makedirs -- a dead service, at boot, with the actual cause
    # (running uid is not the mount's owner) named nowhere in the output.
    #
    # Exit 78 is EX_CONFIG from sysexits.h: the configuration is wrong, as
    # opposed to a crash. A platform that distinguishes them can; one that
    # does not still gets the sentence on stderr.
    problem = state_dir_problem()
    if problem is not None:
        print(f"REFUSING TO START: {problem}", file=sys.stderr)
        sys.stderr.flush()
        return 78

    print(f"Vivum Studio backend")
    print(f"  repo        {root}")
    print(f"  run history {runs_dir()}")
    print(f"  recon cache maxsize={RECON_CACHE.maxsize}")
    print(f"  serving     http://{args.host}:{args.port}")
    # Bring-your-own-tissue is the whole value proposition and it was
    # reachable only by hand-POSTing JSON (results/first_user_readiness.md).
    # Printing the page here is how a first user finds out it exists.
    print(f"  your tissue http://{args.host}:{args.port}/ingest")
    if api_key:
        # Name the source that actually won. There are three, and this line is
        # exactly where an operator looks to work out which: it used to print
        # "set via --api-key" for a key that came from the environment, which
        # is the one situation where the answer is not obvious (an exported
        # TMESIM_STUDIO_API_KEY and an --api-key flag can disagree, and the
        # flag wins) and the one place the banner could have said so.
        if generated:
            label = "generated (not persisted -- save it now)"
        elif args.api_key is not None:
            label = "set via --api-key"
        else:
            label = "set via TMESIM_STUDIO_API_KEY"
        print(f"  api key     {api_key}  [{label}]")
        print(f"              every request needs header: X-API-Key: {api_key}")
    else:
        # No key, for ANY reason -- loopback default, `--api-key ''`, or an
        # empty TMESIM_STUDIO_API_KEY in the environment. The loopback case
        # used to print nothing at all, which is the defect: the risk was
        # invisible precisely where it is easiest to hit (see
        # _no_auth_warning). The default is unchanged -- still no key, still
        # serves -- only now it says so.
        print("  api key     NONE -- no authentication on any /api/ route")
        if not is_loopback:
            print(f"              bound to {args.host}, which is not a "
                  "loopback address")
        print(_no_auth_warning(args.host), end="")

    # Reported for the same reason the key is: the operator's one chance to
    # notice that a capability is or is not confined is the line they are
    # already reading. Deliberately OUTSIDE the key's if/else -- the scan path
    # is governed independently of authentication, and nesting it under either
    # branch would make it invisible in exactly one of the two configurations.
    _roots = os.environ.get(DATA_ROOTS_VAR)
    if _roots is not None:
        _label = ("defaulted for this non-loopback bind" if roots_generated
                  else f"set via {DATA_ROOTS_VAR}")
        print(f"  tissue roots {_roots or 'NONE -- every scan is refused'}  "
              f"[{_label}]")
        if roots_generated:
            print("              tissue outside that directory will be "
                  f"REFUSED; set {DATA_ROOTS_VAR} to include yours")

    # OUTSIDE the key's if/else, for the reason the tissue-roots line above
    # says: this is about the MACHINE, not about authentication, and nesting it
    # under either branch would hide it in exactly one configuration. Placing
    # it INSIDE cost an hour once already -- dropped between the `if` body and
    # its `else`, it parses cleanly and silently rebinds the `else`, so a
    # correctly secured server announces "no authentication on any /api/
    # route". Two statements at the wrong indent, no syntax error, opposite
    # meaning.
    #
    # Silent on a correctly configured host: it fires only when this server's
    # grid ceiling still permits the spec measured to OOM-kill a container of
    # this host's size. An operator who sees this line has a real problem, not
    # boilerplate.
    from .account import oom_warning
    from .api import _MAX_GRID_CELLS, _MAX_GRID_CELLS_SOURCE

    _oom = oom_warning(_MAX_GRID_CELLS, _MAX_GRID_CELLS_SOURCE == "default")
    if _oom:
        print(_oom)

    if args.reload:
        uvicorn.run("app.server.api:app", host=args.host, port=args.port,
                    reload=True, log_level=args.log_level)
        return 0

    from .api import create_app

    uvicorn.run(create_app(), host=args.host, port=args.port,
                log_level=args.log_level)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
