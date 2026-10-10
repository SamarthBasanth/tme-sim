"""auth.py — the login portal's server side.

``docs/AUTH_CONTRACT.md`` IS THE SPECIFICATION. It was frozen before any code
so that this module and the front end could be built against it in parallel.
Every section reference below (§1-§8) points there, and anywhere this file
departs from it the departure is written down in the docstring that departs,
not left for a reader to discover by comparing behaviour with prose.

THE ONE THING TO UNDERSTAND BEFORE EDITING
------------------------------------------
This is the only authentication in the product. There are exactly two
credentials and one gate (§0):

* a **browser session** — the ``vivum_session`` cookie this module issues;
* **``X-API-Key``** — the pre-existing shared secret, which every CI workflow,
  every documented ``curl`` flow and ``render.yaml`` already depend on.

The gate in ``api.py`` accepts EITHER. Nothing here removes, weakens or
replaces the key path; a reader who only ever uses the key notices no change,
and that is the contract's compatibility bar.

WHAT A SESSION IS NOT (§5)
--------------------------
It identifies the person at the keyboard and it partitions NOTHING. Projects,
runs and specimens stay shared, ``GET /api/workspace`` returns the same
workspace to every authenticated user, and no payload here says "your"
anything. The half-built version is the dangerous one: a UI that says "your
projects" over a shared store teaches a false belief about who can see what.
``/api/auth/me`` therefore ships ``shared_workspace: true`` and a sentence
saying so, so the client cannot invent a private space the server does not
have.

CRYPTOGRAPHIC PRIMITIVES — ALL FROM THE STANDARD LIBRARY, NONE INVENTED
-----------------------------------------------------------------------
* passwords: :func:`hashlib.scrypt`, parameters embedded in the stored string
  (§3), with the one documented substitution in :func:`hash_password`'s
  docstring for an interpreter whose SSL library does not provide scrypt;
* sessions: ``hmac.new(secret, body, sha256)`` (§2), and every comparison of a
  secret or a signature in this file goes through
  :func:`secrets.compare_digest`. A ``==`` on either leaks a prefix by timing.
  The gate in ``api.py`` already held that standard for the key; this matches
  it rather than lowering it.

OUT OF SCOPE, DELIBERATELY (§7)
-------------------------------
No email of any kind — no verification, no reset link, no SMTP. No OAuth or
SSO. No roles or permissions: every account can do exactly what the key can
do, because a role model with one role is a column that lies later. No
per-user data partitioning. Nothing weakens the ``X-API-Key`` path.
"""
from __future__ import annotations

import base64
import binascii
import hashlib
import hmac
import ipaddress
import json
import os
import re
import secrets
import threading
import time
from typing import Dict, List, Optional, Tuple

from fastapi import APIRouter, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response

from .cache import (
    PRIVATE_FILE_MODE,
    own_only,
    private_makedirs,
    redact_exception,
    repo_root,
    state_dir,
)

# --------------------------------------------------------------------------- #
# §2 — the session cookie
# --------------------------------------------------------------------------- #

#: Cookie name, fixed by §2 so the front end can rely on it.
COOKIE_NAME = "vivum_session"

#: 14 days (§2). There is deliberately NO refresh-on-use in v1: re-login is
#: cheap, and a sliding window is a thing to get wrong quietly -- the bug it
#: invites is a session that never actually expires because every poll of
#: ``/api/auth/me`` renews it.
SESSION_TTL_S = 14 * 24 * 60 * 60

#: Domain separation for the signing key DERIVED from the store's
#: ``instance_secret``. Not a secret. It exists so the random hex sitting in
#: ``users.json`` is never itself the HMAC key: if that file is ever reused for
#: any other keyed purpose, the two keys are unrelated by construction.
_SESSION_KEY_LABEL = b"vivum-session-signing-key-v1"

#: Floor on an operator-supplied ``TMESIM_SESSION_SECRET``. A short secret is a
#: forgeable session cookie, which is total compromise of every account, and it
#: fails SILENTLY -- there is no symptom until someone tries. So a value too
#: short to be a secret stops the process instead, with the fix named. The
#: precedent is ``TMESIM_MAX_GRID_CELLS``, where a typo also refuses to boot
#: rather than falling back to something dangerous.
SESSION_SECRET_MIN_CHARS = 16

# --------------------------------------------------------------------------- #
# §3 — passwords
# --------------------------------------------------------------------------- #

#: scrypt cost parameters (§3). They live IN the stored string, so raising them
#: later does not invalidate a single existing hash -- an old hash verifies
#: against the parameters it was written with and can be re-hashed on next
#: login whenever someone decides to.
SCRYPT_N = 2 ** 15
SCRYPT_R = 8
SCRYPT_P = 1

#: OpenSSL refuses a scrypt call whose working set exceeds ``maxmem``, and its
#: own default sits a hair above what ``n=2**15, r=8`` needs (128*n*r = 32 MiB)
#: -- close enough that a build with a tighter default would raise instead of
#: hashing. Stated explicitly at double the requirement so the call cannot
#: start failing because of a library default nobody here chose. It is a
#: resource guard, NOT a KDF parameter: it does not enter the derived key, so
#: it has no business in the stored string.
SCRYPT_MAXMEM = 128 * SCRYPT_N * SCRYPT_R * 2

#: 16 bytes from ``secrets.token_bytes`` (§3).
SALT_BYTES = 16

#: Derived-key length. 32 bytes is a 256-bit tag; the length is recoverable
#: from the stored string, so :func:`verify_password` re-derives at exactly the
#: length it finds rather than assuming this constant, and changing it here
#: cannot break an existing hash.
DK_LEN = 32

#: Minimum password length (§3). No composition rules: they push people toward
#: ``Password1!`` and a length floor does more.
PASSWORD_MIN_CHARS = 12

#: Refused BEFORE hashing (§3). scrypt on a megabyte of input is a denial of
#: service, and the check therefore has to come first -- which is why
#: :func:`verify_password` repeats it rather than trusting the sign-up path to
#: have bounded what is now in the store.
PASSWORD_MAX_BYTES = 1024

#: PBKDF2 iterations for the fallback described in :func:`hash_password`.
#: Measured 155 ms on the development interpreter, which is the right order for
#: a login and is itself a rate limiter.
PBKDF2_ITERATIONS = 600_000

#: RFC 5321's limit on a forward path. A bound, not a validator.
EMAIL_MAX_CHARS = 254

#: ``u_`` plus a zero-padded counter (§6). It is a key in a JSON document, not
#: a path segment, but it is also the first field of a SIGNED session body, and
#: pinning its shape means a token can be rejected on syntax before any
#: cryptography runs.
USER_ID_RE = re.compile(r"^u_\d{5,12}$")

#: Deliberately minimal (see :func:`normalise_email`). It is not RFC 5322 and
#: does not pretend to be.
_EMAIL_RE = re.compile(r"^[^\s@,;:<>\"'\\]+@[^\s@,;:<>\"'\\]+\.[^\s@,;:<>\"'\\]+$")

#: §6, verbatim. A store written by a future version says so in its own first
#: field, and this refuses it rather than guessing at a shape it does not know.
STORE_VERSION = 1


class AuthError(Exception):
    """A DELIBERATE refusal, with the HTTP status it means.

    Carries a sentence written FOR THE CALLER and rendered verbatim, exactly as
    ``ProjectError`` and ``IngestError`` are. None of these messages is ever
    built from an exception's text or from a filesystem path -- see
    ``tests/test_response_redaction.py``, whose sweep covers this module like
    every other file under ``app/server/``.
    """

    def __init__(self, message: str, *, status: int = 400,
                 reason: Optional[str] = None,
                 retry_after_s: Optional[int] = None):
        super().__init__(message)
        self.message = message
        self.status = int(status)
        self.reason = reason
        self.retry_after_s = retry_after_s

    def payload(self) -> Dict:
        detail: Dict = {"message": self.message}
        if self.reason:
            detail["reason"] = self.reason
        if self.retry_after_s is not None:
            detail["retry_after_s"] = int(self.retry_after_s)
        return {"detail": detail}


# --------------------------------------------------------------------------- #
# Passwords
# --------------------------------------------------------------------------- #
def _b64(raw: bytes) -> str:
    return base64.b64encode(raw).decode("ascii")


def _unb64(text: str) -> bytes:
    return base64.b64decode(text.encode("ascii"), validate=True)


def scrypt_available() -> bool:
    """Does this interpreter's SSL library provide scrypt?

    ``hashlib.scrypt`` is present only when Python was built against an SSL
    library that implements ``EVP_PBE_scrypt``. OpenSSL does; **LibreSSL does
    not**, and macOS's system Python links LibreSSL 2.8.3, so the attribute is
    simply absent there. See :func:`hash_password` for what that costs.
    """
    return hasattr(hashlib, "scrypt")


def password_problem(password: str) -> Optional[str]:
    """Why this password is unacceptable, as one sentence, or ``None`` (§3).

    Returns a sentence rather than raising so the sign-up and change-password
    routes share one set of rules and one set of words, and so the length floor
    the client enforces can be read off ``/api/auth/me`` instead of being
    re-typed in JavaScript and drifting.
    """
    if not isinstance(password, str) or not password:
        return "a password is required"
    # BEFORE the length-in-characters test, because this is the denial of
    # service and the other is a policy: a 10 MB "password" must be refused
    # without ever reaching a KDF.
    if len(password.encode("utf-8")) > PASSWORD_MAX_BYTES:
        return (f"a password above {PASSWORD_MAX_BYTES} bytes is refused "
                f"before hashing")
    if len(password) < PASSWORD_MIN_CHARS:
        return (f"a password needs at least {PASSWORD_MIN_CHARS} characters; "
                f"there are no composition rules, so length is the whole "
                f"requirement")
    return None


def hash_password(password: str) -> str:
    """``scrypt$<n>$<r>$<p>$<salt_b64>$<hash_b64>`` (§3).

    THE ONE DEVIATION FROM §3 IN THIS FILE, AND IT IS A SUBSTITUTION RATHER
    THAN A WEAKENING OF THE SPECIFIED PATH.

    §3 names ``hashlib.scrypt``, and that is what runs wherever it exists --
    which includes every deployment target: ``Dockerfile`` builds on
    ``python:3.11-slim`` and both GitHub workflows pin Python 3.11, all of them
    OpenSSL, all of them with scrypt. The development interpreter in this
    repository is macOS's Python 3.9.6 against **LibreSSL 2.8.3**, where
    ``hashlib.scrypt`` does not exist at all (``AttributeError`` on the
    attribute, not a failure at call time).

    The alternatives were: require scrypt and have this module fail to import
    on the development machine, which makes the login portal untestable where
    it is written; add ``bcrypt`` or ``argon2``, which §3 forbids for a good
    reason (``constraints.txt`` is a numerically load-bearing pin set and a
    login form is a poor reason to add to it); or write scrypt by hand, which
    is inventing a primitive. So:

    * scrypt present -> ``scrypt$32768$8$1$...``, exactly §3;
    * scrypt absent  -> ``pbkdf2-sha256$600000$...`` via
      :func:`hashlib.pbkdf2_hmac`, which is always present.

    The algorithm NAMES ITSELF in the stored string, so the two coexist in one
    store and :func:`verify_password` dispatches on what it finds. §3 put the
    parameters in the string so they could be raised later; putting the
    algorithm there costs nothing and is what makes this recoverable.

    **PBKDF2-HMAC-SHA256 IS NOT AS GOOD AS SCRYPT AND THIS DOES NOT PRETEND IT
    IS.** It has no memory hardness, so an attacker with GPUs gets a much
    better rate per unit cost against a stolen store. That is why it is the
    fallback and not the default, why the active algorithm is reported on
    ``/api/account`` rather than left implicit, and why this paragraph is here.

    One asymmetry worth knowing: a ``pbkdf2-sha256`` hash verifies on any
    interpreter, a ``scrypt`` hash only on one that has scrypt. Moving a store
    from a Linux deployment to this laptop therefore refuses those logins --
    it fails CLOSED, and ``/api/account`` names the algorithm so the cause is
    visible rather than mysterious.
    """
    problem = password_problem(password)
    if problem is not None:
        raise AuthError(problem, status=400)
    salt = secrets.token_bytes(SALT_BYTES)
    raw = password.encode("utf-8")
    if scrypt_available():
        dk = hashlib.scrypt(raw, salt=salt, n=SCRYPT_N, r=SCRYPT_R,
                            p=SCRYPT_P, dklen=DK_LEN, maxmem=SCRYPT_MAXMEM)
        return (f"scrypt${SCRYPT_N}${SCRYPT_R}${SCRYPT_P}"
                f"${_b64(salt)}${_b64(dk)}")
    dk = hashlib.pbkdf2_hmac("sha256", raw, salt, PBKDF2_ITERATIONS,
                             dklen=DK_LEN)
    return f"pbkdf2-sha256${PBKDF2_ITERATIONS}${_b64(salt)}${_b64(dk)}"


def verify_password(password: str, stored: str) -> bool:
    """Does ``password`` match ``stored``? Constant-time, and never raises.

    Total by construction: a truncated, hand-edited, future-algorithm or
    otherwise unparsable stored string returns ``False``. §7 says a forgotten
    password is reset by the operator editing the store, so a hand-edited
    string is an EXPECTED input here, and the correct answer to one that cannot
    be parsed is "this password does not match" -- not a 500 out of the login
    route, and not an exception whose text could reach the response body.

    The length bound from §3 is re-applied here and not merely at sign-up. The
    store is not the only way a password string arrives: this function is what
    an unauthenticated POST reaches, so the input it must refuse cheaply is the
    one that comes off the wire.
    """
    if not isinstance(password, str) or not isinstance(stored, str):
        return False
    raw = password.encode("utf-8")
    if not raw or len(raw) > PASSWORD_MAX_BYTES:
        return False
    parts = stored.split("$")
    try:
        if parts[0] == "scrypt":
            if len(parts) != 6 or not scrypt_available():
                return False
            n, r, p = int(parts[1]), int(parts[2]), int(parts[3])
            salt, expected = _unb64(parts[4]), _unb64(parts[5])
            # Cost parameters come out of the store, so they are bounded
            # before they are handed to a KDF: a hand-edited `n` of 2**40 is
            # an out-of-memory, not a login failure.
            if not (1 < n <= 2 ** 20 and 0 < r <= 64 and 0 < p <= 16):
                return False
            if not expected:
                return False
            dk = hashlib.scrypt(raw, salt=salt, n=n, r=r, p=p,
                                dklen=len(expected),
                                maxmem=max(SCRYPT_MAXMEM, 128 * n * r * 2))
        elif parts[0] == "pbkdf2-sha256":
            # pbkdf2-sha256$<iterations>$<salt_b64>$<hash_b64> -- four fields,
            # one fewer than scrypt's because there is one cost parameter
            # instead of three.
            if len(parts) != 4:
                return False
            iterations = int(parts[1])
            salt, expected = _unb64(parts[2]), _unb64(parts[3])
            if not (0 < iterations <= 10_000_000) or not expected:
                return False
            dk = hashlib.pbkdf2_hmac("sha256", raw, salt, iterations,
                                     dklen=len(expected))
        else:
            return False
    except (ValueError, TypeError, IndexError, binascii.Error, MemoryError):
        # No exception text escapes, and none is logged as a fault: an
        # unparsable hash is an expected state (§7 invites the operator to edit
        # this file by hand), not an internal error to hand a reference id for.
        return False
    # NEVER `==`: a bytes comparison short-circuits on the first differing
    # byte, which turns a derived key into an oracle you can walk.
    return secrets.compare_digest(dk, expected)


#: A real hash of a value nobody knows, built once per process on first use.
#:
#: WHY: ``POST /api/auth/login`` for an address that does not exist must cost
#: what one for an address that does costs. Returning early on an unknown
#: email makes the route an existence oracle measurable from outside -- the
#: text of the two refusals being identical (which it is) does not help if one
#: of them comes back in 1 ms and the other in 155 ms.
_DUMMY_HASH: Optional[str] = None
_DUMMY_LOCK = threading.Lock()


def _burn_a_verification(password: str) -> None:
    """Do the work a real verification would, and discard the answer."""
    global _DUMMY_HASH
    if _DUMMY_HASH is None:
        with _DUMMY_LOCK:
            if _DUMMY_HASH is None:
                _DUMMY_HASH = hash_password(secrets.token_urlsafe(24))
    verify_password(password, _DUMMY_HASH)


# --------------------------------------------------------------------------- #
# Email
# --------------------------------------------------------------------------- #
def normalise_email(email: str) -> str:
    """Casefolded and stripped (§6), or :class:`AuthError` if it is not usable.

    DELIBERATELY MINIMAL. Validating an address properly means RFC 5322, and
    a regex that half-implements it rejects real addresses -- which, with §7
    forbidding any email path, nobody can appeal. Nothing here sends mail, so
    the only properties this needs are the ones that keep the value safe to
    store in JSON and to match on: one ``@``, a dot in the domain, no
    whitespace or control characters, a length bound.

    Casefolded rather than lowercased: ``str.casefold`` is the Unicode-aware
    operation, and an address differing only in case must be the same account
    -- otherwise ``Sam@x.com`` and ``sam@x.com`` are two accounts and the
    second sign-up silently succeeds.
    """
    if not isinstance(email, str):
        raise AuthError("an email address is required", status=400)
    value = email.strip().casefold()
    if not value:
        raise AuthError("an email address is required", status=400)
    if len(value) > EMAIL_MAX_CHARS:
        raise AuthError(
            f"an email address above {EMAIL_MAX_CHARS} characters is refused",
            status=400)
    if any(ord(ch) < 0x20 or ord(ch) == 0x7F for ch in value):
        raise AuthError("that email address contains a control character",
                        status=400)
    if value.count("@") != 1 or not _EMAIL_RE.match(value):
        raise AuthError(
            "that does not look like an email address",
            status=400,
            reason=("the check is deliberately loose -- one @, a dot in the "
                    "domain, no spaces -- because nothing in this deployment "
                    "sends mail and a strict pattern would reject real "
                    "addresses nobody could appeal"))
    return value


# --------------------------------------------------------------------------- #
# §6 — storage
# --------------------------------------------------------------------------- #
def auth_dir() -> str:
    """``<state_dir>/auth``, else ``app/auth`` (§6).

    The same shape as ``projects.projects_dir``, and for the same reason: one
    persistent disk mounts at ``TMESIM_STATE_DIR`` and every store becomes a
    child of it, so accounts survive a redeploy on a host that has a disk --
    and a host WITHOUT one loses accounts exactly as it loses projects, which
    the home page already discloses.

    Read per call, never cached, so exporting ``TMESIM_STATE_DIR`` after this
    module loaded still takes effect and a test can point it at ``tmp_path``.
    """
    root = state_dir()
    if root:
        return private_makedirs(os.path.join(root, "auth"))
    return private_makedirs(os.path.join(repo_root(), "app", "auth"))


class UserStore:
    """``users.json`` (§6): one file, one lock, same discipline as ``ProjectStore``.

    Not a database, for the reason the project store is not one: this holds a
    handful of rows, the write rate is a login a day, and a plain JSON file is
    inspectable, diffable and -- per §7, which has no password-reset path --
    repairable by hand with an editor.

    Owner-only throughout: ``private_makedirs`` for the directory and
    ``PRIVATE_FILE_MODE`` on the file (§6). Written with the default umask it
    would land ``0644``, world-readable, which for a file of password hashes
    and the session signing key is the whole game. Writes go through a temp
    file and ``os.replace`` so an interrupted write cannot leave a half-parsed
    account store -- a corrupt ``users.json`` locks every user out at once.
    """

    def __init__(self, root: Optional[str] = None):
        self.root = root or auth_dir()
        private_makedirs(self.root)
        self.path = os.path.join(self.root, "users.json")
        self._lock = threading.RLock()

    # -- reading ------------------------------------------------------------ #
    def _blank(self) -> Dict:
        return {"version": STORE_VERSION, "instance_secret": None,
                "users": []}

    def _read(self) -> Dict:
        """The document, or a blank one. Call under ``self._lock``."""
        try:
            with open(self.path, encoding="utf-8") as fh:
                doc = json.load(fh)
        except FileNotFoundError:
            # NOT an error: no file means no accounts yet, which is the
            # first-run state §4 exists to serve.
            return self._blank()
        except (OSError, ValueError) as exc:
            # A store that cannot be read must not look like a store with no
            # accounts in it -- that would turn a permissions mistake into an
            # open sign-up. Refuse loudly, with the detail logged server-side
            # under a reference and NOT in the response body.
            raise AuthError(
                redact_exception(exc, "the account store could not be read"),
                status=503,
                reason=("logins are refused until this is fixed, rather than "
                        "treated as 'no accounts exist' -- which would make "
                        "an unreadable store look like a fresh install")
            ) from exc
        if not isinstance(doc, dict):
            raise AuthError(
                "the account store is not a JSON object",
                status=503,
                reason="see docs/AUTH_CONTRACT.md §6 for the expected shape")
        version = doc.get("version")
        if version != STORE_VERSION:
            raise AuthError(
                f"the account store declares version {version!r}; this server "
                f"understands version {STORE_VERSION}",
                status=503,
                reason=("refusing rather than guessing at a shape written by "
                        "a different version"))
        users = doc.get("users")
        doc["users"] = users if isinstance(users, list) else []
        return doc

    def _write(self, doc: Dict) -> None:
        """Persist atomically, owner-only. Call under ``self._lock``."""
        private_makedirs(self.root)
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, indent=1, sort_keys=False)
        own_only(tmp, PRIVATE_FILE_MODE)
        os.replace(tmp, self.path)
        # The mode travels with the temp file through os.replace, but repeat it
        # on the destination: a pre-existing users.json from before this code
        # (or from a hand edit) keeps whatever mode it had otherwise.
        own_only(self.path, PRIVATE_FILE_MODE)

    # -- the session signing material (§2) ---------------------------------- #
    def instance_secret(self) -> str:
        """The store's own random secret, generated once on first write (§2).

        32 bytes of ``secrets.token_hex``. This is the fallback source for the
        session signing key, and §2 is explicit about what it must NOT be: a
        value derived from ``TMESIM_STUDIO_API_KEY``. Rotating the API key
        would then silently invalidate every session, and the two secrets have
        different lifetimes on purpose -- the key is printed in a dashboard and
        pasted into CI, a session key is never seen by anyone.
        """
        with self._lock:
            doc = self._read()
            existing = doc.get("instance_secret")
            if isinstance(existing, str) and len(existing) >= 32:
                return existing
            fresh = secrets.token_hex(32)
            doc["instance_secret"] = fresh
            self._write(doc)
            return fresh

    # -- queries ------------------------------------------------------------ #
    def count(self) -> int:
        with self._lock:
            return len(self._read()["users"])

    def by_email(self, email: str) -> Optional[Dict]:
        """Matched casefolded (§6), on a value already through
        :func:`normalise_email`."""
        wanted = (email or "").strip().casefold()
        if not wanted:
            return None
        with self._lock:
            for row in self._read()["users"]:
                if isinstance(row, dict) and row.get("email") == wanted:
                    return dict(row)
        return None

    def by_id(self, user_id: str) -> Optional[Dict]:
        if not USER_ID_RE.match(user_id or ""):
            return None
        with self._lock:
            for row in self._read()["users"]:
                if isinstance(row, dict) and row.get("id") == user_id:
                    return dict(row)
        return None

    # -- mutations ---------------------------------------------------------- #
    def create(self, email: str, password_hash: str) -> Dict:
        """Add an account. The uniqueness check and the write are ONE critical
        section, so two simultaneous sign-ups for one address cannot both pass
        the check and both write."""
        with self._lock:
            doc = self._read()
            for row in doc["users"]:
                if isinstance(row, dict) and row.get("email") == email:
                    raise AuthError(
                        "an account already exists for that email address",
                        status=409)
            if doc.get("instance_secret") is None:
                doc["instance_secret"] = secrets.token_hex(32)
            row = {
                "id": self._next_id(doc),
                "email": email,
                "password": password_hash,
                "created_at": time.time(),
                "last_login_at": None,
                "disabled": False,
            }
            doc["users"].append(row)
            self._write(doc)
            return dict(row)

    def _next_id(self, doc: Dict) -> str:
        """Next free ``u_NNNNN``. Call under ``self._lock``.

        Derived from the highest id present rather than from ``len(users)``: a
        row deleted by hand (§7's reset path) would otherwise hand out an id
        that is already in a live session cookie, and that cookie would then
        authenticate as the NEW account.
        """
        highest = 0
        for row in doc["users"]:
            raw = row.get("id", "") if isinstance(row, dict) else ""
            if USER_ID_RE.match(raw):
                highest = max(highest, int(raw[2:]))
        return f"u_{highest + 1:05d}"

    def set_password(self, user_id: str, password_hash: str) -> None:
        with self._lock:
            doc = self._read()
            for row in doc["users"]:
                if isinstance(row, dict) and row.get("id") == user_id:
                    row["password"] = password_hash
                    self._write(doc)
                    return
        raise AuthError("no such account", status=404)

    def mark_login(self, user_id: str) -> None:
        """Stamp ``last_login_at``. Best-effort: a successful login must not
        become a failed one because the stamp could not be written."""
        try:
            with self._lock:
                doc = self._read()
                for row in doc["users"]:
                    if isinstance(row, dict) and row.get("id") == user_id:
                        row["last_login_at"] = time.time()
                        self._write(doc)
                        return
        except (AuthError, OSError):
            return


def public_user(row: Optional[Dict]) -> Optional[Dict]:
    """The fields a client may see. Never the hash, by construction.

    An allowlist and not a ``del row["password"]``: a field added to the store
    later is absent from a payload until someone adds it here, which is the
    direction a mistake should fail in.
    """
    if not row:
        return None
    return {
        "id": row.get("id"),
        "email": row.get("email"),
        "created_at": row.get("created_at"),
        "last_login_at": row.get("last_login_at"),
    }


# --------------------------------------------------------------------------- #
# §2 — issuing and reading a session
# --------------------------------------------------------------------------- #
def issue_session(user_id: str, secret: bytes, *,
                  now: Optional[float] = None,
                  ttl_s: int = SESSION_TTL_S) -> Tuple[str, int]:
    """``<user_id>.<expiry_unix>.<hmac>`` (§2), and the expiry.

    The HMAC is SHA-256 over ``"<user_id>.<expiry>"`` -- the signed body is
    exactly the two fields that precede it, so there is no gap between what
    the cookie claims and what the signature covers. In particular the EXPIRY
    IS INSIDE THE MAC: a client that edits it invalidates the signature, which
    is what makes a stateless session safe to expire without a server-side
    session table.
    """
    expiry = int((time.time() if now is None else now) + ttl_s)
    body = f"{user_id}.{expiry}"
    sig = hmac.new(secret, body.encode("ascii"), hashlib.sha256).hexdigest()
    return f"{body}.{sig}", expiry


def read_session(token: Optional[str], secret: bytes, *,
                 now: Optional[float] = None) -> Optional[str]:
    """The ``user_id`` this token proves, or ``None``. Never raises.

    Order is deliberate: parse, then VERIFY THE SIGNATURE, then check the
    expiry. Nothing in the token is acted on before the HMAC has said the
    server wrote it, and the expiry test is therefore a test on trusted data.

    ``secrets.compare_digest`` on the signature, and this is the line §2 names
    explicitly. A ``==`` here short-circuits on the first differing hex
    character, so the time a rejection takes tells an attacker how much of the
    signature they have right -- 64 hex characters forged in 64 * 16 guesses
    instead of 16**64. The gate in ``api.py`` has always held this standard for
    the key; this is the same standard, not a lower one.
    """
    if not token or not isinstance(token, str):
        return None
    # A cookie value is attacker-controlled: bound it before any work.
    if len(token) > 512:
        return None
    parts = token.split(".")
    if len(parts) != 3:
        return None
    user_id, expiry_raw, sig = parts
    # Syntax first, so a malformed token costs a regex and not an HMAC. The
    # id shape also guarantees it contains no "." -- which is what makes the
    # three-way split unambiguous and the signed body unforgeable by shifting
    # a dot from one field into another.
    if not USER_ID_RE.match(user_id) or not expiry_raw.isdigit():
        return None
    if len(expiry_raw) > 12:
        return None
    body = f"{user_id}.{expiry_raw}"
    expected = hmac.new(secret, body.encode("ascii"), hashlib.sha256).hexdigest()
    if not secrets.compare_digest(sig, expected):
        return None
    if int(expiry_raw) <= (time.time() if now is None else now):
        return None
    return user_id


def request_is_loopback(request: Request) -> bool:
    """Did this request arrive on loopback? Decides ``Secure`` (§2).

    §2 asks for ``Secure`` whenever the request did NOT arrive on loopback,
    and both halves of that are real defects: a hard ``Secure`` breaks local
    development over http, because the browser accepts the ``Set-Cookie`` and
    then never sends it back, so login "succeeds" and nothing is logged in.
    Omitting it on a public host puts a session cookie on the wire in clear.

    SO THIS FAILS CLOSED: ``Secure`` unless the request is *provably* a direct
    loopback plaintext connection. Three signals, in order:

    * **any forwarding header** (``X-Forwarded-For``, ``X-Forwarded-Proto``,
      ``Forwarded``) means a proxy is in front. This matters because behind a
      reverse proxy the PEER IS THE PROXY, frequently ``127.0.0.1`` -- naive
      peer-address detection would drop ``Secure`` on exactly the public
      deployment that needs it most. Treated as not-loopback.
    * **``https`` scheme** -> not loopback, so a local TLS setup still gets
      ``Secure``.
    * otherwise the **peer address**, parsed with :mod:`ipaddress` rather than
      string-matched: ``127.0.0.1`` is not the only loopback address, and
      ``127.0.0.1.evil.com`` must not look like one.

    A client CAN forge ``X-Forwarded-For`` on a loopback request and get a
    ``Secure`` cookie it then cannot send back over http. That costs the
    forger their own login and nobody else anything, which is the right side
    of this trade to be wrong on.
    """
    for header in ("x-forwarded-for", "x-forwarded-proto", "forwarded"):
        if request.headers.get(header):
            return False
    if request.url.scheme == "https":
        return False
    client = request.client
    host = getattr(client, "host", None) if client else None
    if not host:
        # No peer address at all (an ASGI transport that does not supply one).
        # Unknown is not loopback.
        return False
    try:
        return ipaddress.ip_address(host.strip("[]")).is_loopback
    except ValueError:
        return False


def set_session_cookie(response: Response, token: str, *,
                       loopback: bool, max_age: int) -> None:
    """The cookie flags, in one place (§2).

    ``HttpOnly`` so a script cannot read the session -- which is the whole
    reason this is a cookie and not a token in ``localStorage``.
    ``SameSite=Lax`` so a cross-site POST does not carry it, which is this
    design's CSRF defence: the gate accepts the cookie on every ``/api/``
    route, so without ``Lax`` any page the operator visited could drive the
    whole API as them. ``Path=/`` because the gate it satisfies covers the
    whole app. ``Secure`` per :func:`request_is_loopback`.
    """
    response.set_cookie(
        COOKIE_NAME, token,
        max_age=max_age,
        path="/",
        httponly=True,
        samesite="lax",
        secure=not loopback,
    )


def clear_session_cookie(response: Response, *, loopback: bool) -> None:
    """Expire the cookie with the SAME flags it was set with.

    Not optional tidiness: a browser matches a deletion against name, path and
    the security attributes, so clearing with different flags can leave the
    original cookie in place and a logout that silently does nothing.
    """
    response.delete_cookie(
        COOKIE_NAME,
        path="/",
        httponly=True,
        samesite="lax",
        secure=not loopback,
    )


# --------------------------------------------------------------------------- #
# §7 — throttling failed logins
# --------------------------------------------------------------------------- #
class LoginThrottle:
    """Failed-login backoff, keyed so it cannot be aimed at a victim.

    THE PROBLEM WITH THE OBVIOUS DESIGN. A counter keyed on the EMAIL ALONE is
    a lockout primitive handed to the internet: anyone who knows an address
    posts five wrong passwords and the real owner is locked out. The defence
    becomes the attack, and the victim cannot do anything about it -- which is
    worse than the brute force it was meant to stop, because brute force
    against a 12-character minimum behind scrypt is already expensive.

    SO NOTHING HERE IS KEYED ON THE EMAIL ALONE. Two counters, both over a
    sliding window:

    * **(email, peer address)** -- the primary one. An attacker hammering
      ``victim@example.com`` from their own address throttles themselves
      against that address, and the victim logging in from anywhere else is
      completely unaffected. There is no input an attacker can supply that
      slows down a request they did not send.
    * **(peer address)**, across every address tried -- the secondary one,
      against credential stuffing that walks a list of emails from one place.
      A success deliberately does NOT clear this one; otherwise one known-good
      login resets the budget for the guessing.

    THE PEER ADDRESS IS ``request.client.host`` AND NEVER A HEADER. Keying on
    ``X-Forwarded-For`` would be friendlier behind a proxy and would also make
    the throttle decorative, since the attacker chooses that value and can put
    a fresh one on every request. The cost of that choice is stated rather than
    hidden: **behind a reverse proxy every request shares the proxy's address**,
    so the per-peer counter degrades into one global budget. That is why the
    per-peer limit is loose and the (email, peer) limit is the tight one -- a
    tight global counter would be a denial of service against the deployment's
    own users. :data:`CLIENT_MAX` and :data:`PAIR_MAX` are class attributes so
    a deployment that needs different numbers can set them without a new
    environment variable (and without a new line on the deployment page that
    could then go stale).

    WHAT THIS DOES NOT STOP, stated plainly: a distributed attack from many
    addresses against one account. Nothing keyed on the caller can, and the
    thing that actually defends that case is §3 -- a 12-character minimum and a
    memory-hard KDF, i.e. cost per guess. This reduces the cheap cases and is
    not presented as more.

    In-memory and per-process, so it resets on restart and is not shared
    between workers. For a backoff that is acceptable and the alternative is a
    second store; what is NOT acceptable is unbounded growth, so the table is
    capped and evicts least-recently-failed.
    """

    #: Failures from one peer against one address before refusal.
    PAIR_MAX = 5
    #: Failures from one peer across all addresses before refusal. Loose on
    #: purpose; see the class docstring on proxies.
    CLIENT_MAX = 50
    #: Sliding window. Failures older than this are forgotten, so a refusal
    #: always clears itself without an operator doing anything.
    WINDOW_S = 15 * 60
    #: Distinct keys tracked. A bound, so an attacker cycling addresses fills
    #: a fixed-size table rather than the heap.
    MAX_TRACKED = 4096

    def __init__(self):
        self._lock = threading.Lock()
        self._fails: Dict[Tuple[str, str], List[float]] = {}

    @staticmethod
    def _peer(request: Request) -> str:
        client = request.client
        host = getattr(client, "host", None) if client else None
        return host or "unknown"

    def _prune(self, key: Tuple[str, str], now: float) -> List[float]:
        """Call under ``self._lock``."""
        kept = [t for t in self._fails.get(key, ()) if t > now - self.WINDOW_S]
        if kept:
            self._fails[key] = kept
        else:
            self._fails.pop(key, None)
        return kept

    def _keys(self, email: str, request: Request) -> Tuple[Tuple[str, str],
                                                           Tuple[str, str]]:
        peer = self._peer(request)
        return (("pair", f"{peer}|{(email or '')[:EMAIL_MAX_CHARS]}"),
                ("peer", peer))

    def retry_after_s(self, email: str, request: Request, *,
                      now: Optional[float] = None) -> Optional[int]:
        """Seconds until this caller may try again, or ``None``.

        Checked BEFORE the password is verified, which is the point: a
        throttled attempt must cost a dict lookup, not a KDF. Otherwise the
        throttle still spends the server's CPU on every guess and the defence
        pays the attacker's bill.
        """
        now = time.time() if now is None else now
        pair_key, peer_key = self._keys(email, request)
        with self._lock:
            waits = []
            for key, limit in ((pair_key, self.PAIR_MAX),
                               (peer_key, self.CLIENT_MAX)):
                kept = self._prune(key, now)
                if len(kept) >= limit:
                    # The oldest failure still in the window is what has to
                    # age out before there is room for another attempt.
                    waits.append(int(kept[0] + self.WINDOW_S - now) + 1)
            return max(waits) if waits else None

    def record_failure(self, email: str, request: Request, *,
                       now: Optional[float] = None) -> None:
        now = time.time() if now is None else now
        pair_key, peer_key = self._keys(email, request)
        with self._lock:
            for key in (pair_key, peer_key):
                kept = self._prune(key, now)
                kept.append(now)
                self._fails[key] = kept
            while len(self._fails) > self.MAX_TRACKED:
                # Evict whichever key failed least recently. Bounded work and
                # it keeps the keys that are currently under attack.
                oldest = min(self._fails, key=lambda k: self._fails[k][-1])
                self._fails.pop(oldest, None)

    def clear_pair(self, email: str, request: Request) -> None:
        """Forget this (email, peer)'s failures after a SUCCESS.

        The per-peer counter is deliberately left alone -- see the class
        docstring.
        """
        pair_key, _peer_key = self._keys(email, request)
        with self._lock:
            self._fails.pop(pair_key, None)


# --------------------------------------------------------------------------- #
# The service
# --------------------------------------------------------------------------- #
class AuthService:
    """The store, the signing key and the throttle, as one object.

    Built once by ``api.py::create_app`` and hung on ``app.state.auth``, so the
    gate and the routes share one instance -- a second instance would mean a
    second throttle table, and a throttle you can reset by hitting a different
    route is not a throttle.
    """

    def __init__(self, store: Optional[UserStore] = None):
        self.store = store or UserStore()
        self.throttle = LoginThrottle()
        self._secret: Optional[bytes] = None
        self._secret_lock = threading.Lock()

    def secret(self) -> bytes:
        """The session signing key (§2).

        ``TMESIM_SESSION_SECRET`` if set, else derived from the store's own
        random ``instance_secret``. **NOT derived from
        ``TMESIM_STUDIO_API_KEY``**, which §2 forbids and gives the reason for:
        rotating the key would silently invalidate every session, and the two
        secrets have different lifetimes on purpose. There is no code path
        here that reads that variable.

        Cached per process after first use, like the other build-time
        variables: the key cannot be rotated under a running process, which is
        exactly the property ``TMESIM_STUDIO_API_KEY`` already has. The cache
        also means the common request path does not read a file to check a
        cookie.
        """
        if self._secret is not None:
            return self._secret
        with self._secret_lock:
            if self._secret is None:
                self._secret = self._derive_secret()
            return self._secret

    def _derive_secret(self) -> bytes:
        configured = os.environ.get("TMESIM_SESSION_SECRET", "").strip()
        if configured:
            if len(configured) < SESSION_SECRET_MIN_CHARS:
                raise AuthError(
                    f"TMESIM_SESSION_SECRET is {len(configured)} characters; "
                    f"at least {SESSION_SECRET_MIN_CHARS} are required",
                    status=500,
                    reason=("a guessable session secret is a forgeable session "
                            "cookie, which is every account at once, and it "
                            "fails with no symptom -- so this refuses instead. "
                            "Unset the variable to use the store's own random "
                            "secret, or set a long one "
                            "(python -c 'import secrets;"
                            "print(secrets.token_urlsafe(32))')"))
            # Used as the HMAC key directly, which is what §2 says it is.
            return configured.encode("utf-8")
        # Derived, with domain separation, so the hex in users.json is never
        # itself a key (§2).
        material = self.store.instance_secret().encode("ascii")
        return hmac.new(_SESSION_KEY_LABEL, material, hashlib.sha256).digest()

    def secret_source(self) -> str:
        return ("TMESIM_SESSION_SECRET"
                if os.environ.get("TMESIM_SESSION_SECRET", "").strip()
                else "the account store's instance_secret")

    # -- identification ----------------------------------------------------- #
    def user_for_request(self, request: Request) -> Optional[Dict]:
        """The account this request's cookie proves, or ``None``.

        Total: every failure mode -- no cookie, a tampered cookie, an expired
        cookie, a deleted account, a disabled account, an unreadable store --
        returns ``None``. This is what the GATE calls, so the only safe
        direction to fail in is "not authenticated", and the caller then falls
        back to ``X-API-Key`` exactly as before.

        A disabled or deleted account is checked on EVERY request rather than
        only at login. The session is stateless, so there is no server-side
        record to revoke; looking the account up each time is what makes
        ``disabled: true`` take effect at all, and what stops a cookie
        outliving the account it names.
        """
        try:
            token = request.cookies.get(COOKIE_NAME)
            if not token:
                return None
            user_id = read_session(token, self.secret())
            if user_id is None:
                return None
            row = self.store.by_id(user_id)
            if not row or row.get("disabled"):
                return None
            return row
        except AuthError:
            return None
        except OSError:
            return None


# --------------------------------------------------------------------------- #
# §1 — the routes
# --------------------------------------------------------------------------- #
#: The shared sentence that stops the client inventing a private space (§5).
SHARED_WORKSPACE_NOTICE = (
    "An account is a way in, not a private space. Projects, runs and "
    "specimens are shared by everyone who can sign in to this deployment, and "
    "signing in does not create a space of your own."
)

#: The ONE answer a failed login gives, whatever actually went wrong (§7 and
#: the redaction rules). Built once, at module level, so the call sites
#: physically cannot drift into two different sentences -- which is how an
#: existence oracle gets reintroduced by somebody being helpful.
_LOGIN_REFUSAL = (
    "email or password is incorrect",
    ("the same answer comes back whether no account has that address or the "
     "password was wrong. Distinguishing them would turn this route into a "
     "way to ask who has an account here"),
)


def _login_refused() -> AuthError:
    return AuthError(_LOGIN_REFUSAL[0], status=401, reason=_LOGIN_REFUSAL[1])


async def _json_object(request: Request) -> Dict:
    """The body as a JSON object, or a 400 that quotes nothing back.

    Deliberately NOT ``api.py::_json_body``, which answers with
    ``f"...{exc}"`` -- a json decoder's message is harmless enough there, but
    this is the one unauthenticated POST surface in the product and the rule
    it has to keep is that no exception text reaches it. The decoder's own
    words add nothing a client can act on beyond "that was not JSON".
    """
    try:
        body = await request.json()
    except Exception:  # noqa: BLE001 - any decode failure is one 400
        raise AuthError(
            "the request body is not valid JSON",
            status=400,
            reason="send a JSON object with content-type: application/json")
    if not isinstance(body, dict):
        raise AuthError("the request body must be a JSON object", status=400)
    return body


def _text_field(body: Dict, name: str) -> str:
    value = body.get(name)
    if not isinstance(value, str) or not value:
        raise AuthError(f"{name} is required", status=400)
    return value


def _same_secret(presented: str, configured: str) -> bool:
    """``compare_digest`` on two secrets of any encoding.

    Encoded to bytes first because ``secrets.compare_digest`` raises
    ``TypeError`` on a ``str`` holding a non-ASCII character -- so a non-ASCII
    token would otherwise be a 500 out of the sign-up route instead of a
    refusal. Never ``==``, for the reason the gate in ``api.py`` gives.
    """
    return secrets.compare_digest(
        (presented or "").encode("utf-8"), (configured or "").encode("utf-8"))


def _refusal_response(refusal: AuthError) -> JSONResponse:
    """One place where an :class:`AuthError` becomes a response.

    ``Retry-After`` is a real header and not decoration: a client that honours
    it backs off correctly without the front end reimplementing the window.

    THE SECONDS ARE BOUND TO A LOCAL BEFORE BEING FORMATTED, and that shape is
    deliberate -- do not "simplify" it back to
    ``str(int(refusal.retry_after_s))`` inline.
    ``tests/test_response_redaction.py`` binds every parameter annotated with
    an ``Error`` type and flags any ``str()`` over it, because a helper that
    takes a caught exception and formats it is how that leak pattern
    regenerates one call away from the handler. The check is syntactic and
    cannot see that this field is an integer, so the inline version fails the
    sweep. Pulling the int out first makes what is being formatted an integer
    in the source as well as in fact -- which is the honest shape anyway, and
    leaves the sweep watching this module for the thing it is actually for.
    """
    retry_after_s = refusal.retry_after_s
    headers = ({"Retry-After": str(int(retry_after_s))}
               if retry_after_s is not None else None)
    return JSONResponse(refusal.payload(), status_code=refusal.status,
                        headers=headers)


def build_auth_router(auth: AuthService, api_key: str = "") -> APIRouter:
    """The six routes of §1.

    ``api_key`` is the configured ``TMESIM_STUDIO_API_KEY``, PASSED IN rather
    than read from the environment here -- the same reasoning
    ``account.build_account_router`` gives for its limits. ``api.py`` reads it
    once at build time to decide whether to install the gate at all, so
    re-reading it per request here could have sign-up accept a bootstrap token
    the gate is not enforcing. One read, one value, no drift.
    """
    router = APIRouter()
    configured_key = (api_key or "").strip()

    # ------------------------------------------------------------------ #
    # GET /login -- the page. No credential (§1).
    # ------------------------------------------------------------------ #
    @router.get("/login", include_in_schema=False)
    def login_page():
        """Serve ``app/web/login.html``.

        The page itself is the front end's file, not this module's. This route
        exists so ``/login`` works without an extension and so it is
        registered BEFORE the static mount at ``/`` -- starlette matches in
        registration order, and a ``StaticFiles`` mount at ``/`` swallows
        everything under it.

        Ungated by necessity: a key gate in front of the sign-in page is a
        loop. The fallback is shown only when the file is not on disk, so a
        backend deployed ahead of the front end says what is missing instead
        of 404ing like a broken server.
        """
        path = os.path.join(repo_root(), "app", "web", "login.html")
        if os.path.exists(path):
            return FileResponse(path, media_type="text/html",
                                headers={"Cache-Control": "no-cache"})
        return HTMLResponse(_LOGIN_PLACEHOLDER)

    # ------------------------------------------------------------------ #
    # GET /api/auth/me -- session OR key OR nobody. NEVER 401 (§1).
    # ------------------------------------------------------------------ #
    @router.get("/api/auth/me")
    def auth_me(request: Request) -> JSONResponse:
        """Who am I, and what should the page render?

        **This route never returns 401**, which §1 states as a property rather
        than an accident: the page calls it to decide what to draw, so a 401
        here would mean the sign-in page could not load in order to let you
        sign in.

        It also carries the facts the client would otherwise hard-code -- the
        password floor, whether any account exists yet, whether sign-up is
        possible at all -- so there is one source for them instead of a copy
        in JavaScript that goes stale.
        """
        user = auth.user_for_request(request)
        via = None
        if user is not None:
            via = "session"
        elif configured_key and _same_secret(
                request.headers.get("x-api-key", ""), configured_key):
            # §1: the key satisfies this route too. A program holding the key
            # is authenticated and has no account, and saying so is more
            # honest than reporting it as nobody.
            via = "key"
        try:
            accounts = auth.store.count()
        except AuthError:
            # An unreadable store must not 500 the route the page boots on.
            # None is "cannot tell", which is not the same claim as zero.
            accounts = None
        return JSONResponse({
            "authenticated": via is not None,
            "via": via,
            "user": public_user(user),
            "accounts_exist": None if accounts is None else accounts > 0,
            "signup_requires_bootstrap_token": True,
            "signup_available": bool(configured_key),
            "password_min_length": PASSWORD_MIN_CHARS,
            "session": {
                "cookie": COOKIE_NAME,
                "lifetime_days": SESSION_TTL_S // 86400,
                "refresh_on_use": False,
            },
            # §5, in the payload and not only in the UI's memory.
            "shared_workspace": True,
            "notice": SHARED_WORKSPACE_NOTICE,
        })

    # ------------------------------------------------------------------ #
    # POST /api/auth/login -- no credential; sets the cookie (§1).
    # ------------------------------------------------------------------ #
    @router.post("/api/auth/login")
    async def auth_login(request: Request) -> JSONResponse:
        """email + password -> the session cookie.

        THE ORDER OF OPERATIONS IS THE SECURITY OF THIS ROUTE:

        1. parse and normalise -- cheap, and bounds the password's size before
           any KDF touches it (§3);
        2. ask the throttle (§7) -- a refusal here must cost a dict lookup and
           not 150 ms of key derivation, or the throttle spends the server's
           CPU on the attacker's behalf;
        3. verify -- and verify SOMETHING even when the address is unknown, so
           the two outcomes take the same time as well as returning the same
           words;
        4. on success only, clear this (email, peer)'s failures and issue.
        """
        try:
            body = await _json_object(request)
            email = normalise_email(body.get("email"))
            password = _text_field(body, "password")

            wait = auth.throttle.retry_after_s(email, request)
            if wait is not None:
                raise AuthError(
                    "too many failed sign-in attempts from this address",
                    status=429,
                    reason=("the limit is per (address, caller) and per "
                            "caller, never per address alone -- so it cannot "
                            "be used to lock somebody else out"),
                    retry_after_s=wait)

            row = auth.store.by_email(email)
            if row is None or row.get("disabled"):
                # Constant work, then the one refusal. A disabled account is
                # answered identically to an unknown one: "this address is
                # disabled here" is the same disclosure in a politer voice.
                _burn_a_verification(password)
                auth.throttle.record_failure(email, request)
                raise _login_refused()
            if not verify_password(password, row.get("password", "")):
                auth.throttle.record_failure(email, request)
                raise _login_refused()

            auth.throttle.clear_pair(email, request)
            auth.store.mark_login(row["id"])
            token, expiry = issue_session(row["id"], auth.secret())
            fresh = auth.store.by_id(row["id"]) or row
            response = JSONResponse({
                "authenticated": True,
                "via": "session",
                "user": public_user(fresh),
                "expires_at_unix": expiry,
                "shared_workspace": True,
                "notice": SHARED_WORKSPACE_NOTICE,
            })
            set_session_cookie(
                response, token,
                loopback=request_is_loopback(request),
                max_age=max(0, expiry - int(time.time())))
            return response
        except AuthError as refusal:
            return _refusal_response(refusal)

    # ------------------------------------------------------------------ #
    # POST /api/auth/logout -- clears the cookie (§1).
    # ------------------------------------------------------------------ #
    @router.post("/api/auth/logout")
    def auth_logout(request: Request) -> JSONResponse:
        """Clear the cookie, and do it idempotently.

        A DELIBERATE READING OF §1, which lists this route's credential as
        "session". A logout called WITHOUT a valid session is answered 200 and
        the cookie is cleared anyway, rather than 401: the caller's goal is "I
        am not signed in", refusing cannot make that less true, and a 401 here
        would leave a stale or tampered cookie sitting in the browser with no
        way for the page to remove it (``HttpOnly`` means script cannot).
        ``was_authenticated`` reports which case it was, so the client is not
        guessing.
        """
        user = auth.user_for_request(request)
        response = JSONResponse({
            "authenticated": False,
            "was_authenticated": user is not None,
        })
        clear_session_cookie(response, loopback=request_is_loopback(request))
        return response

    # ------------------------------------------------------------------ #
    # POST /api/auth/signup -- bootstrap token (§4).
    # ------------------------------------------------------------------ #
    @router.post("/api/auth/signup")
    async def auth_signup(request: Request) -> JSONResponse:
        """Create an account. Requires the operator's bootstrap token (§4).

        WHY SIGN-UP IS NOT OPEN: an internet-reachable instance with open
        sign-up is claimed by whoever finds it first. The token is the
        operator's existing ``TMESIM_STUDIO_API_KEY`` -- deliberate reuse and
        NOT a new mechanism, because the operator already has that value, it
        already proves "I run this server", and §9 of the workspace contract
        forbids a new API-key mechanism. **The requirement does not lapse once
        the first account exists** (§4).

        No bootstrap token configured means NO SIGN-UP. Falling back to open
        sign-up when the variable is empty would make the one configuration
        that cannot authenticate anybody the one that lets anybody in. On such
        an instance the gate is not installed either, so there is nothing an
        account would get you that is not already open.

        **A successful sign-up does NOT sign you in.** The caller here is the
        operator, who may be creating an account for somebody else; handing
        them a session as that person is a surprise, and the browser case
        costs one extra POST. ``next`` names it.
        """
        try:
            body = await _json_object(request)
            if not configured_key:
                raise AuthError(
                    "sign-up is unavailable on this deployment",
                    status=503,
                    reason=("sign-up requires the operator's bootstrap token, "
                            "and this server has no TMESIM_STUDIO_API_KEY set "
                            "-- so there is no token that could authorise it. "
                            "Set that variable and restart"))
            # Throttled on the CALLER, with a synthetic key belonging to no
            # account -- so guessing at the token can never fill a real
            # address's counter and lock its owner out.
            wait = auth.throttle.retry_after_s(_SIGNUP_THROTTLE_KEY, request)
            if wait is not None:
                raise AuthError(
                    "too many failed sign-up attempts from this address",
                    status=429, retry_after_s=wait)
            presented = (request.headers.get("x-api-key")
                         or body.get("bootstrap_token") or "")
            if not isinstance(presented, str) or not _same_secret(
                    presented, configured_key):
                auth.throttle.record_failure(_SIGNUP_THROTTLE_KEY, request)
                raise AuthError(
                    "sign-up requires the operator's bootstrap token",
                    status=403,
                    reason=("send it as the X-API-Key header or as "
                            "bootstrap_token in the body. It is this "
                            "deployment's TMESIM_STUDIO_API_KEY, which the "
                            "startup banner prints and the hosting dashboard "
                            "shows"))
            auth.throttle.clear_pair(_SIGNUP_THROTTLE_KEY, request)
            email = normalise_email(body.get("email"))
            password = _text_field(body, "password")
            problem = password_problem(password)
            if problem is not None:
                raise AuthError(problem, status=400)
            row = auth.store.create(email, hash_password(password))
            return JSONResponse({
                "user": public_user(row),
                "authenticated": False,
                "next": "POST /api/auth/login",
                "shared_workspace": True,
                "notice": SHARED_WORKSPACE_NOTICE,
            }, status_code=201)
        except AuthError as refusal:
            return _refusal_response(refusal)

    # ------------------------------------------------------------------ #
    # POST /api/auth/password -- session (§1).
    # ------------------------------------------------------------------ #
    @router.post("/api/auth/password")
    async def auth_password(request: Request) -> JSONResponse:
        """Change your own password: old + new.

        SESSION ONLY, and not the key. ``X-API-Key`` identifies a program,
        which has no "own password" to change -- letting the key set any
        account's password would make it a master credential over identities,
        which is more than §4 gave it.

        This route is in the gate's exempt set (like every ``/api/auth/*``
        path), so the 401 below is the handler's own and not the middleware's.
        That is the trap in exempting a prefix: the routes under it must still
        enforce something themselves.
        """
        try:
            user = auth.user_for_request(request)
            if user is None:
                raise AuthError(
                    "this request carries no valid session",
                    status=401,
                    reason=("sign in first. Changing a password needs the "
                            "session for the account being changed; the "
                            "X-API-Key credential identifies a program and "
                            "has no account of its own"))
            body = await _json_object(request)
            old_password = _text_field(body, "old_password")
            new_password = _text_field(body, "new_password")

            email = user.get("email", "")
            wait = auth.throttle.retry_after_s(email, request)
            if wait is not None:
                raise AuthError(
                    "too many failed attempts from this address",
                    status=429, retry_after_s=wait)
            if not verify_password(old_password, user.get("password", "")):
                auth.throttle.record_failure(email, request)
                # Naming the wrong field is not a disclosure here: the caller
                # already proved they hold this account's session, so they
                # learn nothing they did not already have.
                raise AuthError("your current password is incorrect",
                                status=401)
            problem = password_problem(new_password)
            if problem is not None:
                raise AuthError(problem, status=400)
            if secrets.compare_digest(old_password.encode("utf-8"),
                                      new_password.encode("utf-8")):
                raise AuthError(
                    "the new password is the same as the current one",
                    status=400)
            auth.throttle.clear_pair(email, request)
            auth.store.set_password(user["id"], hash_password(new_password))
            return JSONResponse({
                "ok": True,
                # THE HONEST FIELD. A session is a signed string, not a row in
                # a table, so there is nothing to revoke: another browser
                # holding a cookie for this account keeps working until it
                # expires. Reporting this as `true` would be the comfortable
                # lie, and a client that believed it would tell the user they
                # had been signed out everywhere.
                "other_sessions_invalidated": False,
                "other_sessions_note": (
                    "Sessions are signed tokens with no server-side record, "
                    "so changing a password cannot end the ones already "
                    "issued. They expire on their own within "
                    f"{SESSION_TTL_S // 86400} days. To end them now, set a "
                    "new TMESIM_SESSION_SECRET and restart the server, which "
                    "invalidates every session on this deployment at once."),
            })
        except AuthError as refusal:
            return _refusal_response(refusal)

    return router


#: The throttle key sign-up uses. A value :func:`normalise_email` can never
#: produce, so a bootstrap-token guess cannot land in any real account's
#: counter.
_SIGNUP_THROTTLE_KEY = "\x00signup"


#: Shown at ``/login`` when ``app/web/login.html`` is not on disk. The same
#: approach as ``api.py::_PLACEHOLDER``: say what is missing rather than 404.
#: It states §5's point in its own words, because that sentence has to be true
#: of whatever is at this URL -- including this.
_LOGIN_PLACEHOLDER = """<!doctype html><meta charset="utf-8">
<title>Vivum &mdash; sign in</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
 body{background:#13161a;color:#c8ccd2;font:14px/1.6 -apple-system,
   BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;margin:0;
   padding:48px 16px}
 .b{border:1px solid #2a2e33;padding:18px 20px;max-width:620px;margin:0 auto}
 h1{font-size:15px;font-weight:600;letter-spacing:.04em;text-transform:uppercase;
    margin:0 0 14px;color:#f0f2f4}
 code{color:#7fb6e8} p{margin:0 0 12px}
</style>
<div class="b">
<h1>Sign in &mdash; page not built</h1>
<p>The sign-in API is serving. <code>app/web/login.html</code> does not exist
yet, so this placeholder is shown instead.</p>
<p>The routes are <code>POST /api/auth/login</code>,
<code>POST /api/auth/logout</code>, <code>POST /api/auth/signup</code>,
<code>POST /api/auth/password</code> and <code>GET /api/auth/me</code>.</p>
<p>An account is a way in, not a private space: projects, runs and specimens
are shared by everyone who can sign in to this deployment.</p>
</div>
"""
