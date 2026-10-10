"""The login portal: ``docs/AUTH_CONTRACT.md`` turned into assertions.

WHAT THIS FILE IS FOR
---------------------
``app/server/auth.py`` is the only authentication in the product. Every check
here was written by planting the defect it is supposed to catch, confirming the
plant was REACHED, watching this file go red, and reverting -- because this
session has already found nine checks that passed without ever testing their
subject, and an auth test that passes vacuously is worse than no auth test at
all: it is a reason not to look again.

The six plants the brief named, and the test that caught each:

====================================  ====================================
plant                                 caught by
====================================  ====================================
``compare_digest`` -> ``==``          ``TestNoTimingOracle``
an expired session accepted           ``test_an_expired_session_is_refused``
                                      and ``test_an_expired_cookie_does_not
                                      _open_the_gate``
a tampered cookie accepted            ``test_a_tampered_signature_is_
                                      refused`` and its gate twin
sign-up without the token             ``test_signup_without_the_token_is_
                                      refused``
the hash stored in plaintext          ``TestPasswordStorage``
a login error revealing existence     ``test_unknown_email_and_wrong_
                                      password_are_byte_identical``
====================================  ====================================

``TestNoTimingOracle`` deserves a word, because it is the one check here that
reads source rather than behaviour, and that is not laziness. Replacing
``secrets.compare_digest(sig, expected)`` with ``sig == expected`` changes NO
observable answer -- both accept exactly the same tokens. The whole difference
is how long a rejection takes, and the leak is that a short-circuiting compare
tells an attacker how many leading characters of a forged signature are right,
which turns 16**64 guesses into about 64*16. A behavioural test cannot see
that without timing measurements far too noisy to put in a suite, so the
property is asserted where it is decidable: in the syntax.

COST
----
No reconstruction, no PDE, no real tissue. The expense is key derivation,
which is the point of a password KDF: one hash is ~155 ms on the development
interpreter, so accounts are created by inserting ONE pre-computed hash
(``known_hash``, session-scoped) wherever the KDF is not the subject.
"""
from __future__ import annotations

import ast
import json
import os
import re
import stat
import sys
import time

import pytest
from starlette.testclient import TestClient

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app.server import auth as auth_mod                       # noqa: E402
from app.server import cache as cache_mod                     # noqa: E402
from app.server import compounds as compounds_mod             # noqa: E402
from app.server import ingest_registry                        # noqa: E402
from app.server.api import create_app                         # noqa: E402
from app.server.jobs import JobManager                        # noqa: E402
from app.server.projects import ProjectStore                  # noqa: E402

#: The bootstrap token (§4) for every gated app in this file.
API_KEY = "k" * 48

#: Twelve characters exactly, which is the floor (§3) -- so the boundary is
#: exercised by the ordinary happy path rather than by one test nobody reads.
PASSWORD = "correct-horse"
OTHER_PASSWORD = "another-valid-passphrase"
EMAIL = "Operator@Example.COM"
EMAIL_CASEFOLDED = "operator@example.com"

#: Loopback, so ``request_is_loopback`` is true and the cookie is issued
#: WITHOUT ``Secure`` -- which is what lets httpx's cookie jar send it back
#: over http. A test client at the default ("testclient", 50000) gets a
#: ``Secure`` cookie that http can never return, and every session test would
#: fail for a reason that has nothing to do with sessions. ``Secure`` itself is
#: asserted in :class:`TestCookieFlags`, off the raw header.
LOOPBACK = ("127.0.0.1", 51000)
OTHER_PEER = ("203.0.113.7", 51001)


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="session")
def known_hash() -> str:
    """``hash_password(PASSWORD)``, computed once for the whole session."""
    return auth_mod.hash_password(PASSWORD)


def _isolate(tmp_path, monkeypatch, *, api_key: str):
    """Point every store at ``tmp_path`` and build the app.

    ``TMESIM_STATE_DIR`` is what carries the account store (§6), so setting it
    here is also the test that ``auth_dir()`` honours it -- see
    :func:`test_the_store_lives_under_the_state_dir`.
    """
    for var in ("TMESIM_STUDIO_DATA_ROOTS", "TMESIM_STUDIO_EXPOSE_TRACEBACKS",
                "TMESIM_STUDIO_ALLOWED_ORIGINS", "TMESIM_STUDIO_ALLOWED_HOSTS",
                "TMESIM_SESSION_SECRET"):
        monkeypatch.delenv(var, raising=False)
    if api_key:
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", api_key)
    else:
        monkeypatch.delenv("TMESIM_STUDIO_API_KEY", raising=False)
    monkeypatch.setenv("TMESIM_STATE_DIR", str(tmp_path / "state"))

    runs = tmp_path / "runs"
    runs.mkdir(exist_ok=True)
    monkeypatch.setattr(cache_mod, "runs_dir", lambda: str(runs))
    monkeypatch.setattr(ingest_registry, "_store", None)
    monkeypatch.setattr(
        compounds_mod, "COMPOUND_STORE",
        compounds_mod.CompoundStore(root=str(tmp_path / "compounds")))
    manager = JobManager(root=str(runs / "jobs"), start_worker=False)
    app = create_app(manager, projects=ProjectStore(root=str(tmp_path / "p")))
    return app, manager


@pytest.fixture
def gated(tmp_path, monkeypatch):
    """The app with the gate ON, and a client whose peer is loopback."""
    app, manager = _isolate(tmp_path, monkeypatch, api_key=API_KEY)
    with TestClient(app, client=LOOPBACK) as client:
        yield client
    manager.shutdown()
    monkeypatch.setattr(ingest_registry, "_store", None)


@pytest.fixture
def ungated(tmp_path, monkeypatch):
    """The app with no key configured -- the local workflow (§0's
    compatibility bar)."""
    app, manager = _isolate(tmp_path, monkeypatch, api_key="")
    with TestClient(app, client=LOOPBACK) as client:
        yield client
    manager.shutdown()
    monkeypatch.setattr(ingest_registry, "_store", None)


def _service(client: TestClient):
    return client.app.state.auth


def _add_account(client: TestClient, email: str, password_hash: str) -> dict:
    """Insert an account with a pre-computed hash, skipping the KDF."""
    return _service(client).store.create(
        auth_mod.normalise_email(email), password_hash)


def _sign_in(client: TestClient, email: str = EMAIL,
             password: str = PASSWORD):
    return client.post("/api/auth/login",
                       json={"email": email, "password": password})


# --------------------------------------------------------------------------- #
# §3 — passwords
# --------------------------------------------------------------------------- #
class TestPasswordStorage:
    """The stored string, its parameters, and that it is not the password.

    PLANTED to confirm this class is not vacuous, in the version of the
    defect that does NOT announce itself: ``hash_password`` returning
    ``password`` unchanged AND ``verify_password`` growing a plaintext branch,
    so ``users.json`` held ``correct-horse`` in clear while **every login,
    gate and route test in this file stayed green**. Four checks here went
    red. That is the argument for asserting on the STORED STRING rather than
    trusting that a working login implies a hashed one.
    """

    def test_the_stored_string_names_its_algorithm_and_parameters(self):
        """§3: the parameters live IN the string so they can be raised later
        without invalidating existing hashes."""
        stored = auth_mod.hash_password(PASSWORD)
        fields = stored.split("$")
        if auth_mod.scrypt_available():
            assert fields[0] == "scrypt"
            assert fields[1:4] == [str(auth_mod.SCRYPT_N),
                                   str(auth_mod.SCRYPT_R),
                                   str(auth_mod.SCRYPT_P)], stored
            assert int(fields[1]) == 2 ** 15, "§3 pins n = 2**15"
            assert (int(fields[2]), int(fields[3])) == (8, 1), "§3 pins r, p"
        else:
            # The documented substitution for an interpreter whose SSL library
            # has no scrypt. It must still NAME itself, which is the property
            # that lets one store hold both.
            assert fields[0] == "pbkdf2-sha256", stored
            assert int(fields[1]) == auth_mod.PBKDF2_ITERATIONS
        assert len(fields) in (4, 6), stored

    def test_the_password_is_not_recoverable_from_the_store(self, tmp_path,
                                                            monkeypatch):
        """The whole point. Not the password, not reversibly encoded, not a
        bare digest of it."""
        import base64
        import hashlib

        stored = auth_mod.hash_password(PASSWORD)
        assert PASSWORD not in stored
        raw = PASSWORD.encode("utf-8")
        for encoded in (base64.b64encode(raw).decode(),
                        base64.b16encode(raw).decode(),
                        base64.b32encode(raw).decode(),
                        raw.hex(),
                        hashlib.sha256(raw).hexdigest(),
                        hashlib.md5(raw).hexdigest()):
            assert encoded not in stored, (
                "the stored string contains a reversible or unsalted encoding "
                "of the password")

    def test_the_salt_is_fresh_every_time(self):
        """16 random bytes per hash (§3), so two accounts with one password do
        not share a stored string and a rainbow table buys nothing."""
        a, b = auth_mod.hash_password(PASSWORD), auth_mod.hash_password(PASSWORD)
        assert a != b
        salt_a, salt_b = a.split("$")[-2], b.split("$")[-2]
        import base64
        assert len(base64.b64decode(salt_a)) == auth_mod.SALT_BYTES
        assert salt_a != salt_b

    def test_it_round_trips_and_a_wrong_password_does_not(self):
        stored = auth_mod.hash_password(PASSWORD)
        assert auth_mod.verify_password(PASSWORD, stored) is True
        assert auth_mod.verify_password(PASSWORD + "x", stored) is False
        assert auth_mod.verify_password(PASSWORD.upper(), stored) is False
        assert auth_mod.verify_password("", stored) is False

    def test_a_plaintext_store_entry_does_not_verify(self):
        """The plant, pinned as a property.

        If a future edit stored the password instead of a hash,
        ``verify_password`` must not accept it -- the format has no "plain"
        branch, so a plaintext entry is unparsable and refused. That makes the
        plaintext defect fail CLOSED (nobody can log in) rather than open.
        """
        assert auth_mod.verify_password(PASSWORD, PASSWORD) is False
        assert auth_mod.verify_password(PASSWORD, "plain$" + PASSWORD) is False

    @pytest.mark.parametrize("stored", [
        "", "scrypt", "scrypt$$$$$", "scrypt$32768$8$1$notb64$notb64",
        "scrypt$32768$8$1$AAAA", "pbkdf2-sha256$0$AAAA$AAAA",
        "pbkdf2-sha256$600000$AAAA$", "bcrypt$12$whatever",
        "scrypt$1099511627776$8$1$AAAA$AAAA",
    ])
    def test_a_malformed_stored_string_is_refused_not_raised(self, stored):
        """§7 invites the operator to edit this file by hand, so a broken
        entry is an EXPECTED input. The right answer is False -- not a 500 out
        of the login route, and not an exception whose text could reach the
        response body. The last case is a hand-edited cost parameter that
        would otherwise be an out-of-memory, not a login failure.
        """
        assert auth_mod.verify_password(PASSWORD, stored) is False

    @pytest.mark.parametrize("password,why", [
        ("", "empty"),
        ("short", "five characters"),
        ("elevenchars", "eleven characters, one below the floor"),
    ])
    def test_a_password_below_the_floor_is_refused(self, password, why):
        assert auth_mod.password_problem(password) is not None, why
        with pytest.raises(auth_mod.AuthError):
            auth_mod.hash_password(password)

    def test_twelve_characters_is_accepted(self):
        """The floor is 12 and not 13. A boundary asserted in both
        directions, because an off-by-one here silently rejects a password the
        documentation promises."""
        assert len("a" * 12) == auth_mod.PASSWORD_MIN_CHARS
        assert auth_mod.password_problem("a" * 12) is None

    def test_an_oversized_password_is_refused_before_hashing(self):
        """§3: scrypt on a megabyte of input is a denial of service, so the
        bound is checked BEFORE the KDF -- and on the verify path too, which
        is the one an unauthenticated POST reaches."""
        huge = "a" * (auth_mod.PASSWORD_MAX_BYTES + 1)
        assert auth_mod.password_problem(huge) is not None
        with pytest.raises(auth_mod.AuthError):
            auth_mod.hash_password(huge)
        stored = auth_mod.hash_password(PASSWORD)
        started = time.time()
        assert auth_mod.verify_password(huge, stored) is False
        assert time.time() - started < 0.05, (
            "an oversized password reached the KDF on the verify path")

    def test_the_byte_bound_is_bytes_and_not_characters(self):
        """A multi-byte character must not slip past a character count: 400
        astral characters are 1,600 bytes."""
        assert auth_mod.password_problem("\U0001F300" * 400) is not None


# --------------------------------------------------------------------------- #
# §2 — sessions
# --------------------------------------------------------------------------- #
class TestSessionTokens:
    SECRET = b"s" * 32

    def test_the_shape_is_the_contract(self):
        """§2: ``<user_id>.<expiry_unix>.<hmac>``, HMAC SHA-256 over
        ``"<user_id>.<expiry>"``."""
        import hashlib
        import hmac

        token, expiry = auth_mod.issue_session("u_00001", self.SECRET)
        user_id, expiry_raw, sig = token.split(".")
        assert user_id == "u_00001"
        assert int(expiry_raw) == expiry
        assert re.fullmatch(r"[0-9a-f]{64}", sig), "SHA-256, hex"
        assert sig == hmac.new(self.SECRET, b"u_00001." + expiry_raw.encode(),
                               hashlib.sha256).hexdigest(), (
            "the signature is not an HMAC-SHA256 over exactly the two fields "
            "that precede it")

    def test_the_lifetime_is_fourteen_days(self):
        assert auth_mod.SESSION_TTL_S == 14 * 24 * 3600
        _token, expiry = auth_mod.issue_session("u_00001", self.SECRET)
        assert 13.9 * 86400 < expiry - time.time() < 14.1 * 86400

    def test_a_valid_token_reads_back(self):
        token, _ = auth_mod.issue_session("u_00042", self.SECRET)
        assert auth_mod.read_session(token, self.SECRET) == "u_00042"

    def test_an_expired_session_is_refused(self):
        """PLANTED (``int(expiry_raw) <= now`` -> ``False``): red here."""
        token, _ = auth_mod.issue_session("u_00001", self.SECRET, ttl_s=-1)
        assert auth_mod.read_session(token, self.SECRET) is None

    def test_expiry_is_refused_exactly_at_the_boundary(self):
        """``<=``, not ``<``: a token whose expiry equals now is spent."""
        now = 1_000_000.0
        token, expiry = auth_mod.issue_session(
            "u_00001", self.SECRET, now=now, ttl_s=100)
        assert auth_mod.read_session(
            token, self.SECRET, now=expiry - 1) == "u_00001"
        assert auth_mod.read_session(token, self.SECRET, now=expiry) is None

    def test_a_tampered_signature_is_refused(self):
        """One character of the signature flipped. PLANTED
        (``compare_digest`` -> ``True``): red here."""
        token, _ = auth_mod.issue_session("u_00001", self.SECRET)
        head, sig = token.rsplit(".", 1)
        for index in (0, 31, 63):
            flipped = list(sig)
            flipped[index] = "0" if sig[index] != "0" else "1"
            assert auth_mod.read_session(
                f"{head}.{''.join(flipped)}", self.SECRET) is None, (
                f"a signature with character {index} flipped was accepted")

    def test_an_extended_expiry_is_refused(self):
        """The expiry is INSIDE the MAC, which is what makes a stateless
        session safe to expire: editing it invalidates the signature."""
        token, _ = auth_mod.issue_session("u_00001", self.SECRET, ttl_s=-1)
        user_id, _expiry, sig = token.split(".")
        forged = f"{user_id}.{int(time.time()) + 99999}.{sig}"
        assert auth_mod.read_session(forged, self.SECRET) is None

    def test_a_token_signed_with_another_key_is_refused(self):
        token, _ = auth_mod.issue_session("u_00001", b"d" * 32)
        assert auth_mod.read_session(token, self.SECRET) is None

    @pytest.mark.parametrize("token", [
        "", "u_00001", "u_00001.1", "u_00001.1.2.3", "..",
        "u_00001.notanumber.abc", "../../etc/passwd.1.abc",
        "u_1.1.abc", "u_00001.-1.abc", "u_00001." + "9" * 40 + ".abc",
        "x" * 600,
    ])
    def test_a_malformed_token_is_refused_not_raised(self, token):
        assert auth_mod.read_session(token, self.SECRET) is None

    def test_a_user_id_shaped_to_smuggle_a_dot_is_refused(self):
        """The id shape guarantees it holds no ``.``, which is what makes the
        three-way split unambiguous.

        If an id could contain a dot, one signed body could be re-cut into a
        different (id, expiry) pair under the same signature -- so the shape
        is enforced on the way in as well as on the way out. A token for a
        dotted id simply cannot authenticate: it splits into four fields.
        """
        assert auth_mod.USER_ID_RE.match("u_00001.123") is None
        smuggled, _ = auth_mod.issue_session("u_00001.123", self.SECRET)
        assert len(smuggled.split(".")) == 4
        assert auth_mod.read_session(smuggled, self.SECRET) is None
        # And the id field is validated before the HMAC runs, so a correctly
        # signed token for a malformed id is still refused.
        for bad in ("u_1", "u_" + "9" * 20, "U_00001", "admin", ""):
            token, _ = auth_mod.issue_session(bad, self.SECRET)
            assert auth_mod.read_session(token, self.SECRET) is None, bad


class TestSessionSecret:
    """§2's three rules about where the signing key comes from."""

    def test_it_is_not_derived_from_the_api_key(self, tmp_path, monkeypatch):
        """§2 forbids this and gives the reason: rotating the API key would
        otherwise sign out every browser, and the two secrets have different
        lifetimes on purpose.

        Asserted behaviourally -- issue a session, rotate the key, and the
        session must still verify -- and then structurally, because a
        behavioural test alone would pass if the derivation merely ignored the
        key by accident.
        """
        monkeypatch.delenv("TMESIM_SESSION_SECRET", raising=False)
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "first-key-value-aaaa")
        service = auth_mod.AuthService(
            auth_mod.UserStore(root=str(tmp_path / "auth")))
        before = service.secret()
        token, _ = auth_mod.issue_session("u_00001", before)

        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "second-key-value-bbbb")
        fresh = auth_mod.AuthService(
            auth_mod.UserStore(root=str(tmp_path / "auth")))
        assert fresh.secret() == before, (
            "the signing key moved when TMESIM_STUDIO_API_KEY was rotated")
        assert auth_mod.read_session(token, fresh.secret()) == "u_00001"

        source = open(os.path.join(_REPO, "app", "server", "auth.py"),
                      encoding="utf-8").read()
        derive = source.split("def _derive_secret")[1].split("\n    def ")[0]
        assert "TMESIM_STUDIO_API_KEY" not in derive, (
            "_derive_secret reads the API key; §2 forbids deriving the "
            "session secret from it")

    def test_the_env_var_wins_when_set(self, tmp_path, monkeypatch):
        monkeypatch.setenv("TMESIM_SESSION_SECRET", "x" * 40)
        service = auth_mod.AuthService(
            auth_mod.UserStore(root=str(tmp_path / "auth")))
        assert service.secret() == ("x" * 40).encode("utf-8")
        assert service.secret_source() == "TMESIM_SESSION_SECRET"

    def test_otherwise_it_derives_from_the_stores_own_random_secret(
            self, tmp_path, monkeypatch):
        monkeypatch.delenv("TMESIM_SESSION_SECRET", raising=False)
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        service = auth_mod.AuthService(store)
        secret = service.secret()
        assert len(secret) == 32
        instance = store.instance_secret()
        assert len(instance) >= 64, "32 bytes of token_hex"
        # Derived, not used raw: the hex in users.json is never itself a key.
        assert secret != instance.encode("ascii")
        assert secret != bytes.fromhex(instance)
        assert "instance_secret" in service.secret_source()

    def test_two_instances_invent_different_secrets(self, tmp_path):
        a = auth_mod.UserStore(root=str(tmp_path / "a")).instance_secret()
        b = auth_mod.UserStore(root=str(tmp_path / "b")).instance_secret()
        assert a != b

    def test_a_short_configured_secret_refuses_rather_than_warns(
            self, tmp_path, monkeypatch):
        """A guessable session secret is a forgeable cookie, i.e. every
        account at once, and it has no symptom. So it is a refusal."""
        monkeypatch.setenv("TMESIM_SESSION_SECRET", "tooshort")
        service = auth_mod.AuthService(
            auth_mod.UserStore(root=str(tmp_path / "auth")))
        with pytest.raises(auth_mod.AuthError) as caught:
            service.secret()
        assert "TMESIM_SESSION_SECRET" in str(caught.value)
        assert auth_mod.SESSION_SECRET_MIN_CHARS == 16


class TestNoTimingOracle:
    """Every secret and signature comparison in ``auth.py`` is
    ``compare_digest``.

    THIS IS THE ONE SOURCE-READING CHECK IN THIS FILE, and the module
    docstring explains why it has to be: ``==`` accepts exactly the tokens
    ``compare_digest`` accepts, so no behavioural assertion can tell them
    apart. What differs is the time a rejection takes, and a short-circuiting
    compare over a 64-character hex signature turns 16**64 guesses into about
    64*16.

    PLANTED, both ways: ``read_session``'s ``compare_digest`` replaced with
    ``sig == expected``, and ``_same_secret``'s with ``presented ==
    configured``. Every behavioural test in this file stayed GREEN; this class
    went red on each. That is the whole argument for its existence.
    """

    #: A comparison involving any of these names is comparing a secret.
    SECRET_NAMES = frozenset({
        "sig", "signature", "expected", "secret", "token", "presented",
        "configured", "password", "dk", "digest", "api_key", "stored",
        "old_password", "new_password",
    })

    def _tree(self):
        path = os.path.join(_REPO, "app", "server", "auth.py")
        with open(path, encoding="utf-8") as fh:
            return ast.parse(fh.read()), fh.name

    def test_no_equality_comparison_of_a_secret(self):
        tree, path = self._tree()
        offenders = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Compare):
                continue
            if not any(isinstance(op, (ast.Eq, ast.NotEq))
                       for op in node.ops):
                continue
            names = {n.id for n in ast.walk(node)
                     if isinstance(n, ast.Name)}
            names |= {n.attr for n in ast.walk(node)
                      if isinstance(n, ast.Attribute)}
            hit = names & self.SECRET_NAMES
            if hit:
                offenders.append(
                    f"{path}:{node.lineno}  {ast.unparse(node)}  "
                    f"(names {sorted(hit)})")
        assert not offenders, (
            "a secret or signature is compared with ==/!=, which leaks its "
            "prefix by timing. Use secrets.compare_digest:\n"
            + "\n".join(offenders))

    def test_the_signature_check_uses_compare_digest(self):
        """Named function by named function, so deleting the comparison
        outright is caught too -- not only replacing it."""
        tree, _path = self._tree()
        wanted = {"read_session", "_same_secret", "verify_password"}
        seen = set()
        for node in ast.walk(tree):
            if (isinstance(node, ast.FunctionDef)
                    and node.name in wanted):
                calls = {ast.unparse(n.func) for n in ast.walk(node)
                         if isinstance(n, ast.Call)}
                assert any(c.endswith("compare_digest") for c in calls), (
                    f"{node.name} does not call secrets.compare_digest")
                seen.add(node.name)
        assert seen == wanted, f"functions not found: {sorted(wanted - seen)}"

    def test_the_gate_still_uses_compare_digest_for_the_key(self):
        """The standard this file matches rather than lowers. If the gate's
        own comparison regressed, the session work would be defending a door
        with a hole next to it."""
        path = os.path.join(_REPO, "app", "server", "api.py")
        with open(path, encoding="utf-8") as fh:
            source = fh.read()
        # Sliced to the end of the function rather than a fixed byte count: a
        # character window silently stopped covering the comparison once the
        # gate grew a comment block, which is a check quietly ceasing to test
        # its subject -- the exact failure this file exists to avoid.
        gate = source.split("async def _require_api_key")[1]
        gate = gate.split("\n    jobs = ")[0]
        assert "x-api-key" in gate, "the slice no longer covers the gate body"
        assert "secrets.compare_digest(" in gate
        for leak in ("== api_key", "api_key ==", "!= api_key", "api_key !="):
            assert leak not in gate, f"the gate compares the key with {leak!r}"


# --------------------------------------------------------------------------- #
# §6 — storage
# --------------------------------------------------------------------------- #
class TestStore:
    def test_the_store_lives_under_the_state_dir(self, tmp_path, monkeypatch):
        """§6: ``<state_dir>/auth/users.json``, so accounts survive a redeploy
        on a host that has a disk."""
        monkeypatch.setenv("TMESIM_STATE_DIR", str(tmp_path / "state"))
        assert auth_mod.auth_dir() == str(tmp_path / "state" / "auth")
        store = auth_mod.UserStore()
        assert store.path == str(tmp_path / "state" / "auth" / "users.json")

    def test_it_falls_back_to_app_auth_when_unset(self, monkeypatch):
        """``app/auth`` with no state dir, mirroring ``projects_dir``.

        ``private_makedirs`` is stubbed to the identity because the real one
        would CREATE ``app/auth/`` inside the developer's checkout as a side
        effect of asserting a path. A test that litters the repository it is
        testing gets noticed and then disabled, and an empty directory is
        invisible to ``git status``, so nobody would notice this one until it
        shipped in an image.
        """
        monkeypatch.delenv("TMESIM_STATE_DIR", raising=False)
        monkeypatch.setattr(auth_mod, "private_makedirs", lambda path: path)
        assert auth_mod.auth_dir().endswith(os.path.join("app", "auth"))

    def test_it_is_read_per_call_not_cached_at_import(self, tmp_path,
                                                      monkeypatch):
        monkeypatch.setenv("TMESIM_STATE_DIR", str(tmp_path / "one"))
        first = auth_mod.auth_dir()
        monkeypatch.setenv("TMESIM_STATE_DIR", str(tmp_path / "two"))
        assert auth_mod.auth_dir() != first

    def test_the_file_and_directory_are_owner_only(self, tmp_path,
                                                   known_hash):
        """§6: ``0700``/``0600`` via ``private_makedirs`` and
        ``PRIVATE_FILE_MODE``. Written with the default umask this file of
        password hashes and the session signing key would land 0644 --
        readable by every local account on the machine."""
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        store.create("a@b.co", known_hash)
        assert stat.S_IMODE(os.stat(store.path).st_mode) == 0o600
        assert stat.S_IMODE(os.stat(store.root).st_mode) == 0o700
        assert auth_mod.PRIVATE_FILE_MODE == 0o600

    def test_the_document_shape_is_the_contract(self, tmp_path, known_hash):
        """§6's JSON, field for field."""
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        store.create("a@b.co", known_hash)
        with open(store.path, encoding="utf-8") as fh:
            doc = json.load(fh)
        assert doc["version"] == 1
        assert re.fullmatch(r"[0-9a-f]{64}", doc["instance_secret"])
        assert set(doc["users"][0]) == {
            "id", "email", "password", "created_at", "last_login_at",
            "disabled"}
        row = doc["users"][0]
        assert row["id"] == "u_00001"
        assert row["last_login_at"] is None and row["disabled"] is False

    def test_email_is_stored_and_matched_casefolded(self, tmp_path,
                                                    known_hash):
        """§6. Without this, ``Sam@x.com`` and ``sam@x.com`` are two accounts
        and the second sign-up silently succeeds."""
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        store.create(auth_mod.normalise_email("Op@Example.COM"), known_hash)
        assert store.by_email("op@example.com")["email"] == "op@example.com"
        assert store.by_email("OP@EXAMPLE.COM") is not None
        assert store.by_email("  Op@Example.com  ") is not None

    def test_a_duplicate_email_is_refused(self, tmp_path, known_hash):
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        store.create("a@b.co", known_hash)
        with pytest.raises(auth_mod.AuthError) as caught:
            store.create("a@b.co", known_hash)
        assert caught.value.status == 409

    def test_ids_come_from_the_highest_present_not_the_count(self, tmp_path,
                                                             known_hash):
        """A row deleted by hand (§7's reset path) must not cause an id to be
        REUSED -- a live cookie naming it would then authenticate as the new
        account."""
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        store.create("a@b.co", known_hash)
        store.create("c@d.co", known_hash)
        doc = store._read()
        doc["users"] = [r for r in doc["users"] if r["id"] != "u_00001"]
        store._write(doc)
        assert store.create("e@f.co", known_hash)["id"] == "u_00003"

    def test_an_unreadable_store_refuses_rather_than_looking_empty(
            self, tmp_path):
        """A permissions mistake must not look like a fresh install -- that
        would turn it into open sign-up."""
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        with open(store.path, "w", encoding="utf-8") as fh:
            fh.write("{not json")
        with pytest.raises(auth_mod.AuthError) as caught:
            store.count()
        assert caught.value.status == 503

    def test_a_future_version_is_refused(self, tmp_path):
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        with open(store.path, "w", encoding="utf-8") as fh:
            json.dump({"version": 99, "users": []}, fh)
        with pytest.raises(auth_mod.AuthError):
            store.count()

    def test_a_store_fault_text_carries_no_path(self, tmp_path):
        """The redaction rule, on the one error path here that wraps an
        unexpected exception."""
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        with open(store.path, "w", encoding="utf-8") as fh:
            fh.write("{not json")
        with pytest.raises(auth_mod.AuthError) as caught:
            store.count()
        message = caught.value.message
        assert str(tmp_path) not in message
        assert "users.json" not in message
        assert "reference" in message

    def test_the_public_projection_never_carries_the_hash(self, tmp_path,
                                                          known_hash):
        store = auth_mod.UserStore(root=str(tmp_path / "auth"))
        row = store.create("a@b.co", known_hash)
        public = auth_mod.public_user(row)
        assert "password" not in public
        assert known_hash not in json.dumps(public)
        assert set(public) == {"id", "email", "created_at", "last_login_at"}


class TestEmailNormalisation:
    @pytest.mark.parametrize("value", [
        "", "   ", "noatsign", "no@domain", "two@at@signs.com",
        "spaces in@example.com", "a@b.co\nx", "a@b.co\x00",
        "a" * 250 + "@example.com",
    ])
    def test_a_bad_address_is_refused(self, value):
        with pytest.raises(auth_mod.AuthError):
            auth_mod.normalise_email(value)

    @pytest.mark.parametrize("value", [
        "a@b.co", "Operator@Example.COM", "  x@y.org  ", "a+tag@b.co.uk",
    ])
    def test_a_usable_address_survives(self, value):
        assert "@" in auth_mod.normalise_email(value)
        assert auth_mod.normalise_email(value) == value.strip().casefold()


# --------------------------------------------------------------------------- #
# §2 — cookie flags
# --------------------------------------------------------------------------- #
class TestCookieFlags:
    """§2: ``HttpOnly``, ``SameSite=Lax``, ``Path=/``, and ``Secure`` only
    when the request did not arrive on loopback.

    Read off the raw ``set-cookie`` header rather than the client's jar,
    because the jar is exactly what hides a wrong ``Secure``.
    """

    def _set_cookie(self, client, known_hash, **kwargs):
        _add_account(client, EMAIL, known_hash)
        response = client.post("/api/auth/login",
                               json={"email": EMAIL, "password": PASSWORD},
                               **kwargs)
        assert response.status_code == 200, response.text
        header = response.headers["set-cookie"]
        assert header.startswith(f"{auth_mod.COOKIE_NAME}=")
        return header

    def test_the_always_on_flags(self, gated, known_hash):
        header = self._set_cookie(gated, known_hash)
        lowered = header.lower()
        assert "httponly" in lowered, "a script could read the session"
        assert "samesite=lax" in lowered, (
            "without SameSite=Lax a cross-site POST carries the cookie and "
            "the whole API is CSRF-able as the signed-in user")
        assert "path=/" in lowered
        assert f"max-age={auth_mod.SESSION_TTL_S}" in lowered.replace(" ", "")

    def test_no_secure_flag_on_loopback(self, gated, known_hash):
        """A hard ``Secure`` breaks local development over http: the browser
        accepts the Set-Cookie and then never sends it back, so login
        "succeeds" and nothing is signed in."""
        header = self._set_cookie(gated, known_hash)
        assert "secure" not in header.lower(), header

    def test_secure_when_a_proxy_is_in_front(self, gated, known_hash):
        """Behind a reverse proxy THE PEER IS THE PROXY, frequently
        127.0.0.1. Naive peer-address detection would drop ``Secure`` on
        exactly the public deployment that needs it most -- so a forwarding
        header means not-loopback even though this request's peer IS
        loopback."""
        header = self._set_cookie(
            gated, known_hash, headers={"x-forwarded-for": "198.51.100.4"})
        assert "secure" in header.lower(), header

    def test_secure_on_a_non_loopback_peer(self, tmp_path, monkeypatch,
                                           known_hash):
        app, manager = _isolate(tmp_path, monkeypatch, api_key=API_KEY)
        try:
            with TestClient(app, client=OTHER_PEER) as client:
                header = self._set_cookie(client, known_hash)
                assert "secure" in header.lower(), header
        finally:
            manager.shutdown()

    def test_secure_under_https(self, tmp_path, monkeypatch, known_hash):
        app, manager = _isolate(tmp_path, monkeypatch, api_key=API_KEY)
        try:
            with TestClient(app, client=LOOPBACK,
                            base_url="https://testserver") as client:
                header = self._set_cookie(client, known_hash)
                assert "secure" in header.lower(), header
        finally:
            manager.shutdown()

    def test_logout_clears_with_the_same_flags(self, gated, known_hash):
        """A browser matches a deletion on name, path AND the security
        attributes, so clearing with different flags can leave the original
        cookie in place -- a logout that silently does nothing."""
        self._set_cookie(gated, known_hash)
        header = gated.post("/api/auth/logout").headers["set-cookie"]
        lowered = header.lower()
        assert "path=/" in lowered
        assert "httponly" in lowered and "samesite=lax" in lowered
        assert ("max-age=0" in lowered.replace(" ", "")
                or "expires=thu, 01 jan 1970" in lowered)


# --------------------------------------------------------------------------- #
# §1 — the routes
# --------------------------------------------------------------------------- #
class TestMe:
    def test_it_never_401s_without_a_credential(self, gated):
        """§1 states this as a property: the page asks this to decide what to
        render, so a 401 would mean the sign-in page cannot load in order to
        let you sign in."""
        response = gated.get("/api/auth/me")
        assert response.status_code == 200, response.text
        body = response.json()
        assert body["authenticated"] is False
        assert body["via"] is None and body["user"] is None
        assert body["accounts_exist"] is False

    def test_it_reports_the_key_as_a_credential(self, gated):
        body = gated.get("/api/auth/me",
                         headers={"x-api-key": API_KEY}).json()
        assert body["authenticated"] is True and body["via"] == "key"
        assert body["user"] is None, "the key is a program, not an account"

    def test_it_reports_the_session(self, gated, known_hash):
        _add_account(gated, EMAIL, known_hash)
        assert _sign_in(gated).status_code == 200
        body = gated.get("/api/auth/me").json()
        assert body["authenticated"] is True and body["via"] == "session"
        assert body["user"]["email"] == EMAIL_CASEFOLDED
        # The HASH must not be here. Scanned for the stored string and for a
        # literal "password" KEY, rather than for the substring "password" --
        # which `password_min_length` contains legitimately, and which would
        # make this assertion fail for a reason that is not a leak.
        serialised = json.dumps(body)
        assert known_hash not in serialised
        assert '"password"' not in serialised
        assert "scrypt$" not in serialised and "pbkdf2" not in serialised

    def test_it_carries_the_facts_the_client_would_otherwise_hard_code(
            self, gated):
        body = gated.get("/api/auth/me").json()
        assert body["password_min_length"] == auth_mod.PASSWORD_MIN_CHARS
        assert body["signup_requires_bootstrap_token"] is True
        assert body["session"]["cookie"] == "vivum_session"
        assert body["session"]["lifetime_days"] == 14
        assert body["session"]["refresh_on_use"] is False

    def test_it_states_that_nothing_is_partitioned(self, gated):
        """§5, in the payload and not only in the UI's memory."""
        body = gated.get("/api/auth/me").json()
        assert body["shared_workspace"] is True
        assert "not a private space" in body["notice"]


class TestSignup:
    def test_signup_without_the_token_is_refused(self, gated):
        """§4. PLANTED: the ``if not _same_secret(...)`` guard replaced with
        ``if False:``, so no token was required at all. A stranger then
        claimed the instance with ``201`` and ``u_00001``; red here and in
        the six parametrised cases below."""
        response = gated.post("/api/auth/signup",
                              json={"email": "x@y.co", "password": PASSWORD})
        assert response.status_code == 403, response.text
        assert gated.get("/api/auth/me").json()["accounts_exist"] is False

    @pytest.mark.parametrize("token", [
        "", "wrong", API_KEY[:-1], API_KEY + "x", API_KEY.upper(),
    ])
    def test_a_near_miss_token_is_refused(self, gated, token):
        response = gated.post("/api/auth/signup",
                              json={"email": "x@y.co", "password": PASSWORD,
                                    "bootstrap_token": token})
        assert response.status_code in (403, 429), response.text

    def test_the_token_is_accepted_in_the_header_or_the_body(self, gated):
        first = gated.post("/api/auth/signup",
                           json={"email": "a@y.co", "password": PASSWORD},
                           headers={"x-api-key": API_KEY})
        assert first.status_code == 201, first.text
        second = gated.post("/api/auth/signup",
                            json={"email": "b@y.co", "password": PASSWORD,
                                  "bootstrap_token": API_KEY})
        assert second.status_code == 201, second.text

    def test_the_token_is_still_required_after_the_first_account(
            self, gated, known_hash):
        """§4 says so explicitly: there is no first-run window that closes."""
        _add_account(gated, EMAIL, known_hash)
        response = gated.post("/api/auth/signup",
                              json={"email": "second@y.co",
                                    "password": PASSWORD})
        assert response.status_code == 403, response.text

    def test_the_created_account_is_returned_without_its_hash(self, gated):
        body = gated.post("/api/auth/signup",
                          json={"email": "New@Y.CO", "password": PASSWORD},
                          headers={"x-api-key": API_KEY}).json()
        assert body["user"]["email"] == "new@y.co", "casefolded (§6)"
        assert body["user"]["id"] == "u_00001"
        assert "password" not in json.dumps(body)

    def test_signup_does_not_sign_you_in(self, gated):
        """The operator may be creating an account for somebody else; handing
        them a session as that person is a surprise."""
        response = gated.post("/api/auth/signup",
                              json={"email": "x@y.co", "password": PASSWORD},
                              headers={"x-api-key": API_KEY})
        assert response.json()["authenticated"] is False
        assert "set-cookie" not in response.headers
        assert gated.get("/api/auth/me").json()["authenticated"] is False

    def test_a_duplicate_is_refused(self, gated):
        for _ in range(1):
            gated.post("/api/auth/signup",
                       json={"email": "x@y.co", "password": PASSWORD},
                       headers={"x-api-key": API_KEY})
        again = gated.post("/api/auth/signup",
                           json={"email": "X@Y.CO", "password": PASSWORD},
                           headers={"x-api-key": API_KEY})
        assert again.status_code == 409, again.text

    def test_a_short_password_is_refused(self, gated):
        response = gated.post("/api/auth/signup",
                              json={"email": "x@y.co", "password": "short"},
                              headers={"x-api-key": API_KEY})
        assert response.status_code == 400, response.text
        assert "12" in response.json()["detail"]["message"]

    def test_signup_is_unavailable_with_no_key_configured(self, ungated):
        """Falling open here would make the one configuration that cannot
        authenticate anybody the one that lets anybody in."""
        response = ungated.post("/api/auth/signup",
                                json={"email": "x@y.co",
                                      "password": PASSWORD})
        assert response.status_code == 503, response.text
        assert "bootstrap token" in response.json()["detail"]["reason"]


class TestLogin:
    def test_a_good_password_signs_you_in(self, gated, known_hash):
        _add_account(gated, EMAIL, known_hash)
        response = _sign_in(gated)
        assert response.status_code == 200, response.text
        body = response.json()
        assert body["authenticated"] is True
        assert body["user"]["email"] == EMAIL_CASEFOLDED
        assert body["expires_at_unix"] > time.time()
        assert auth_mod.COOKIE_NAME in gated.cookies

    def test_the_email_is_matched_casefolded(self, gated, known_hash):
        _add_account(gated, EMAIL, known_hash)
        assert _sign_in(gated, email="OPERATOR@EXAMPLE.COM").status_code == 200

    def test_last_login_is_stamped(self, gated, known_hash):
        _add_account(gated, EMAIL, known_hash)
        assert _sign_in(gated).json()["user"]["last_login_at"] is not None

    def test_unknown_email_and_wrong_password_are_byte_identical(
            self, gated, known_hash):
        """THE EXISTENCE ORACLE. PLANTED (the unknown-email branch given its
        own "no account with that address" message): red here.

        Body, status and headers compared whole -- a difference anywhere in
        the response is a way to ask who has an account on this deployment.
        """
        _add_account(gated, EMAIL, known_hash)
        wrong_password = _sign_in(gated, password="wrong-but-long-enough")
        unknown_email = _sign_in(gated, email="nobody@example.com",
                                 password="wrong-but-long-enough")
        assert wrong_password.status_code == 401
        assert unknown_email.status_code == 401
        assert wrong_password.json() == unknown_email.json(), (
            "the two refusals differ, so this route can be asked whether an "
            "address has an account")
        assert "set-cookie" not in wrong_password.headers
        assert "set-cookie" not in unknown_email.headers
        # Byte-identity above is the property that matters; this scans the
        # MESSAGE for a wording that would name which half failed. Only the
        # message, because the `reason` field explains the policy and
        # necessarily uses the words "no account has that address" -- it is a
        # constant that is identical in both responses, so it discloses
        # nothing, and scanning it would be scanning the explanation for the
        # thing it explains.
        message = wrong_password.json()["detail"]["message"].lower()
        for leak in ("no such", "not found", "unknown", "no account",
                     "does not exist", "incorrect password", "wrong password",
                     "disabled"):
            assert leak not in message, f"the refusal message says {leak!r}"
        assert message == "email or password is incorrect"

    def test_a_disabled_account_is_refused_identically(self, gated,
                                                       known_hash):
        """"This address is disabled here" is the same disclosure in a
        politer voice."""
        _add_account(gated, EMAIL, known_hash)
        service = _service(gated)
        doc = service.store._read()
        doc["users"][0]["disabled"] = True
        service.store._write(doc)
        disabled = _sign_in(gated)
        unknown = _sign_in(gated, email="nobody@example.com")
        assert disabled.status_code == 401
        assert disabled.json() == unknown.json()

    def test_an_unknown_email_costs_the_same_work(self, gated, known_hash):
        """Identical WORDS do not help if one refusal returns in 1 ms and the
        other in 155 ms. The unknown-email path burns a real verification.

        The bound is deliberately loose -- this asserts the same order of
        magnitude, not a tight timing, because a tight one would be flaky on
        a loaded machine and would get deleted.
        """
        _add_account(gated, EMAIL, known_hash)
        started = time.time()
        _sign_in(gated, password="wrong-but-long-enough")
        real = time.time() - started
        started = time.time()
        _sign_in(gated, email="nobody@example.com",
                 password="wrong-but-long-enough")
        unknown = time.time() - started
        assert unknown > real / 4, (
            f"an unknown address was refused in {unknown:.3f}s against "
            f"{real:.3f}s for a wrong password -- the gap is an existence "
            f"oracle measurable from outside")

    @pytest.mark.parametrize("body,status", [
        ({}, 400),
        ({"email": "x@y.co"}, 400),
        ({"password": PASSWORD}, 400),
        ({"email": "notanemail", "password": PASSWORD}, 400),
        ({"email": "x@y.co", "password": 12345}, 400),
        ({"email": None, "password": None}, 400),
    ])
    def test_a_malformed_body_is_a_400_not_a_500(self, gated, body, status):
        assert _sign_in.__name__  # keep the helper referenced
        assert gated.post("/api/auth/login", json=body).status_code == status

    def test_a_non_json_body_is_a_400_quoting_nothing(self, gated):
        response = gated.post("/api/auth/login", content=b"{not json",
                              headers={"content-type": "application/json"})
        assert response.status_code == 400
        text = json.dumps(response.json())
        for leak in ("JSONDecode", "Expecting", "line 1", "char 1"):
            assert leak not in text, f"the decoder's own text reached the wire"


class TestLogout:
    def test_it_clears_the_session(self, gated, known_hash):
        _add_account(gated, EMAIL, known_hash)
        _sign_in(gated)
        assert gated.get("/api/auth/me").json()["authenticated"] is True
        response = gated.post("/api/auth/logout")
        assert response.status_code == 200
        assert response.json() == {"authenticated": False,
                                   "was_authenticated": True}
        assert gated.get("/api/auth/me").json()["authenticated"] is False

    def test_it_is_idempotent(self, gated):
        """A 401 here would leave a stale cookie in the browser with no way
        for the page to remove it -- ``HttpOnly`` means script cannot."""
        response = gated.post("/api/auth/logout")
        assert response.status_code == 200
        assert response.json()["was_authenticated"] is False


class TestPasswordChange:
    def _signed_in(self, client, known_hash):
        _add_account(client, EMAIL, known_hash)
        assert _sign_in(client).status_code == 200

    def test_it_requires_a_session(self, gated):
        response = gated.post("/api/auth/password",
                              json={"old_password": PASSWORD,
                                    "new_password": OTHER_PASSWORD})
        assert response.status_code == 401, response.text

    def test_the_api_key_is_not_enough(self, gated, known_hash):
        """The key identifies a PROGRAM, which has no own password. Letting it
        set any account's password would make it a master credential over
        identities -- more than §4 gave it."""
        _add_account(gated, EMAIL, known_hash)
        response = gated.post("/api/auth/password",
                              json={"old_password": PASSWORD,
                                    "new_password": OTHER_PASSWORD},
                              headers={"x-api-key": API_KEY})
        assert response.status_code == 401, response.text

    def test_a_wrong_current_password_is_refused(self, gated, known_hash):
        self._signed_in(gated, known_hash)
        response = gated.post("/api/auth/password",
                              json={"old_password": "wrong-but-long",
                                    "new_password": OTHER_PASSWORD})
        assert response.status_code == 401, response.text

    def test_it_changes_the_password(self, gated, known_hash):
        self._signed_in(gated, known_hash)
        response = gated.post("/api/auth/password",
                              json={"old_password": PASSWORD,
                                    "new_password": OTHER_PASSWORD})
        assert response.status_code == 200, response.text
        gated.post("/api/auth/logout")
        assert _sign_in(gated, password=PASSWORD).status_code == 401
        assert _sign_in(gated, password=OTHER_PASSWORD).status_code == 200

    def test_the_new_password_must_clear_the_floor(self, gated, known_hash):
        self._signed_in(gated, known_hash)
        response = gated.post("/api/auth/password",
                              json={"old_password": PASSWORD,
                                    "new_password": "short"})
        assert response.status_code == 400

    def test_it_says_plainly_that_other_sessions_survive(self, gated,
                                                         known_hash):
        """THE HONEST FIELD. A session is a signed string, not a row in a
        table, so there is nothing to revoke. Reporting this as ``true`` would
        be the comfortable lie, and a client believing it would tell the user
        they had been signed out everywhere."""
        self._signed_in(gated, known_hash)
        body = gated.post("/api/auth/password",
                          json={"old_password": PASSWORD,
                                "new_password": OTHER_PASSWORD}).json()
        assert body["other_sessions_invalidated"] is False
        assert "TMESIM_SESSION_SECRET" in body["other_sessions_note"]

    def test_another_browsers_cookie_really_does_survive(self, gated,
                                                         known_hash):
        """The claim above, measured rather than asserted. If this ever starts
        failing, the field is now a lie in the other direction and must be
        flipped to ``true``."""
        self._signed_in(gated, known_hash)
        stolen = gated.cookies[auth_mod.COOKIE_NAME]
        gated.post("/api/auth/password",
                   json={"old_password": PASSWORD,
                         "new_password": OTHER_PASSWORD})
        gated.cookies.clear()
        gated.cookies.set(auth_mod.COOKIE_NAME, stolen)
        assert gated.get("/api/auth/me").json()["authenticated"] is True


# --------------------------------------------------------------------------- #
# §7 — the throttle
# --------------------------------------------------------------------------- #
class TestThrottle:
    def test_failures_from_one_caller_are_throttled(self, gated, known_hash):
        _add_account(gated, EMAIL, known_hash)
        limit = auth_mod.LoginThrottle.PAIR_MAX
        for attempt in range(limit):
            response = _sign_in(gated, password="wrong-but-long-enough")
            assert response.status_code == 401, (attempt, response.text)
        blocked = _sign_in(gated, password="wrong-but-long-enough")
        assert blocked.status_code == 429, blocked.text
        assert int(blocked.headers["retry-after"]) > 0
        assert blocked.json()["detail"]["retry_after_s"] > 0

    def test_the_right_password_is_also_refused_while_throttled(
            self, gated, known_hash):
        """Otherwise the throttle is a speed bump an attacker drives over the
        moment they guess right."""
        _add_account(gated, EMAIL, known_hash)
        for _ in range(auth_mod.LoginThrottle.PAIR_MAX):
            _sign_in(gated, password="wrong-but-long-enough")
        assert _sign_in(gated).status_code == 429

    def test_a_victim_elsewhere_is_not_locked_out(self, tmp_path, monkeypatch,
                                                  known_hash):
        """THE POINT OF THE WHOLE DESIGN, and the reason nothing here is keyed
        on the email alone.

        An attacker fills the budget against ``EMAIL`` from their own address.
        The account's real owner, at a different address, must sign in
        normally -- otherwise the throttle is a lockout primitive handed to
        anyone who knows an email address, and the defence has become the
        attack.
        """
        app, manager = _isolate(tmp_path, monkeypatch, api_key=API_KEY)
        try:
            attacker = TestClient(app, client=("198.51.100.9", 4444))
            victim = TestClient(app, client=LOOPBACK)
            with attacker, victim:
                _add_account(victim, EMAIL, known_hash)
                for _ in range(auth_mod.LoginThrottle.PAIR_MAX + 2):
                    attacker.post("/api/auth/login",
                                  json={"email": EMAIL,
                                        "password": "wrong-but-long"})
                assert attacker.post(
                    "/api/auth/login",
                    json={"email": EMAIL, "password": PASSWORD}
                ).status_code == 429, "the attacker is not throttled at all"
                assert _sign_in(victim).status_code == 200, (
                    "the account owner was locked out by somebody else's "
                    "failed attempts -- the throttle is keyed on the email")
        finally:
            manager.shutdown()

    def test_a_success_clears_the_pair(self, gated, known_hash):
        _add_account(gated, EMAIL, known_hash)
        for _ in range(auth_mod.LoginThrottle.PAIR_MAX - 1):
            _sign_in(gated, password="wrong-but-long-enough")
        assert _sign_in(gated).status_code == 200
        for _ in range(auth_mod.LoginThrottle.PAIR_MAX - 1):
            assert _sign_in(
                gated, password="wrong-but-long-enough").status_code == 401

    def test_the_window_rolls_without_an_operator(self, gated, known_hash):
        """A refusal that needs a human to clear it is an outage."""
        throttle = _service(gated).throttle
        request = _FakeRequest("127.0.0.1")
        now = 1_000_000.0
        for _ in range(throttle.PAIR_MAX):
            throttle.record_failure("a@b.co", request, now=now)
        assert throttle.retry_after_s("a@b.co", request, now=now) is not None
        later = now + throttle.WINDOW_S + 1
        assert throttle.retry_after_s("a@b.co", request, now=later) is None

    def test_a_refusal_costs_no_key_derivation(self, gated, known_hash):
        """Checked BEFORE the password is verified -- otherwise the throttle
        still spends the server's CPU on every guess and the defence pays the
        attacker's bill."""
        _add_account(gated, EMAIL, known_hash)
        for _ in range(auth_mod.LoginThrottle.PAIR_MAX):
            _sign_in(gated, password="wrong-but-long-enough")
        started = time.time()
        assert _sign_in(gated, password="wrong-but-long").status_code == 429
        assert time.time() - started < 0.05, (
            "a throttled attempt still ran the KDF")

    def test_the_table_is_bounded(self, gated):
        """An attacker cycling addresses must fill a fixed-size table rather
        than the heap."""
        throttle = _service(gated).throttle
        for index in range(throttle.MAX_TRACKED + 200):
            throttle.record_failure(f"user{index}@example.com",
                                    _FakeRequest(f"10.0.{index % 250}.1"))
        assert len(throttle._fails) <= throttle.MAX_TRACKED

    def test_signup_throttling_cannot_touch_a_real_account(self, gated,
                                                           known_hash):
        """Guessing at the bootstrap token must not fill any account's
        counter -- otherwise sign-up becomes the lockout tool that §7's
        keying was designed to prevent."""
        _add_account(gated, EMAIL, known_hash)
        for _ in range(auth_mod.LoginThrottle.PAIR_MAX + 1):
            gated.post("/api/auth/signup",
                       json={"email": EMAIL, "password": PASSWORD,
                             "bootstrap_token": "wrong"})
        assert _sign_in(gated).status_code == 200, (
            "failed sign-ups for this address blocked its owner's login")


class _FakeRequest:
    """The two attributes the throttle reads. Enough to drive it with an
    injected clock, which is what makes the window test deterministic instead
    of a sleep."""

    def __init__(self, host: str):
        self.client = type("C", (), {"host": host})()
        self.headers = {}


# --------------------------------------------------------------------------- #
# §0 — the gate takes either credential
# --------------------------------------------------------------------------- #
class TestGate:
    def test_the_key_still_works(self, gated):
        """§0's compatibility bar: a reader who only ever uses the key
        notices no change."""
        assert gated.get("/api/runs",
                         headers={"x-api-key": API_KEY}).status_code == 200

    def test_neither_credential_is_refused(self, gated):
        response = gated.get("/api/runs")
        assert response.status_code == 401, response.text

    def test_a_session_opens_the_gate(self, gated, known_hash):
        _add_account(gated, EMAIL, known_hash)
        assert _sign_in(gated).status_code == 200
        response = gated.get("/api/runs")
        assert response.status_code == 200, (
            "a signed-in browser could not reach a data route")
        assert gated.get("/api/workspace").status_code == 200
        assert gated.get("/api/account").status_code == 200

    def test_a_tampered_cookie_does_not_open_the_gate(self, gated,
                                                      known_hash):
        """PLANTED (``compare_digest`` -> ``True`` in ``read_session``): red
        here, which is the gate-level twin of the unit check."""
        _add_account(gated, EMAIL, known_hash)
        _sign_in(gated)
        good = gated.cookies[auth_mod.COOKIE_NAME]
        head, sig = good.rsplit(".", 1)
        flipped = ("0" if sig[0] != "0" else "1") + sig[1:]
        gated.cookies.clear()
        gated.cookies.set(auth_mod.COOKIE_NAME, f"{head}.{flipped}")
        assert gated.get("/api/runs").status_code == 401
        assert gated.get("/api/auth/me").json()["authenticated"] is False

    def test_an_expired_cookie_does_not_open_the_gate(self, gated,
                                                      known_hash):
        """PLANTED (the expiry test removed from ``read_session``): red
        here."""
        row = _add_account(gated, EMAIL, known_hash)
        stale, _ = auth_mod.issue_session(
            row["id"], _service(gated).secret(), ttl_s=-60)
        gated.cookies.set(auth_mod.COOKIE_NAME, stale)
        assert gated.get("/api/runs").status_code == 401
        assert gated.get("/api/auth/me").json()["authenticated"] is False

    def test_a_cookie_for_a_deleted_account_does_not_open_the_gate(
            self, gated, known_hash):
        """The session is stateless, so looking the account up on every
        request is what stops a cookie outliving the account it names."""
        _add_account(gated, EMAIL, known_hash)
        _sign_in(gated)
        assert gated.get("/api/runs").status_code == 200
        service = _service(gated)
        doc = service.store._read()
        doc["users"] = []
        service.store._write(doc)
        assert gated.get("/api/runs").status_code == 401

    def test_a_cookie_for_a_disabled_account_does_not_open_the_gate(
            self, gated, known_hash):
        """``disabled: true`` (§6) has to take effect somewhere, and with no
        server-side session record the only place is per request."""
        _add_account(gated, EMAIL, known_hash)
        _sign_in(gated)
        service = _service(gated)
        doc = service.store._read()
        doc["users"][0]["disabled"] = True
        service.store._write(doc)
        assert gated.get("/api/runs").status_code == 401

    def test_a_cookie_signed_with_another_deployments_key_is_refused(
            self, gated, known_hash):
        row = _add_account(gated, EMAIL, known_hash)
        forged, _ = auth_mod.issue_session(row["id"], b"someone-elses-key" * 2)
        gated.cookies.set(auth_mod.COOKIE_NAME, forged)
        assert gated.get("/api/runs").status_code == 401

    @pytest.mark.parametrize("path", [
        "/login", "/api/auth/me",
    ])
    def test_the_portal_is_reachable_without_a_credential(self, gated, path):
        """§1: these are the only new ungated paths. A credential gate in
        front of the page that obtains the credential is a loop."""
        assert gated.get(path).status_code == 200, path

    def test_the_login_post_is_reachable_without_a_credential(self, gated):
        response = gated.post("/api/auth/login",
                              json={"email": "x@y.co", "password": "x" * 12})
        assert response.status_code != 401 or (
            response.json()["detail"]["message"]
            != "missing or incorrect X-API-Key header"), (
            "the gate refused the login route, which nobody can then pass")

    def test_the_liveness_exemption_is_intact(self, gated):
        """Deliberate, and the reason is recorded in ``api.py``: a platform
        probe cannot carry a key, and pointing one at a gated route killed a
        healthy container twice."""
        for path in ("/healthz", "/api/health"):
            assert gated.get(path).status_code == 200, path
            assert gated.head(path).status_code == 200, path
        anonymous = gated.get("/api/health").json()
        assert anonymous == {"ok": True}, (
            "the liveness route disclosed operational detail to an "
            "unauthenticated caller")
        authenticated = gated.get(
            "/api/health", headers={"x-api-key": API_KEY}).json()
        assert "runs_dir" in authenticated

    def test_a_post_to_healthz_still_needs_a_credential(self, gated):
        """The exemption is GET/HEAD only, so the fail-closed rule is
        untouched."""
        assert gated.post("/healthz").status_code in (401, 405)

    def test_a_host_header_that_smuggles_a_path_does_not_skip_the_gate(
            self, gated):
        """GHSA-86qp-5c8j-p5mr, re-asserted because the gate's decision
        changed.

        ``request.url`` is rebuilt from the CLIENT-SUPPLIED Host header, so
        ``Host: 127.0.0.1:8000/../`` yields ``url.path='/..//api/runs'`` -- the
        ``startswith('/api/')`` test fails, the gate is skipped, and the router
        still dispatches. The fix is to decide on
        ``get_route_path(request.scope)``, and THE NEW EXEMPT SET MUST BE
        COMPUTED ON THAT SAME STRING: a ``/api/auth/`` prefix test against
        ``request.url.path`` would have reopened the hole in a new place.
        """
        for host in ("127.0.0.1:8000/../", "evil.example/..", "x/../",
                     "127.0.0.1:8000/..%2f", "x/../api/auth/",
                     "x/..//api/auth/me",
                     # THESE FOUR WERE ADDED BY A PLANT, and they are the
                     # point of re-asserting this at all.
                     #
                     # Planting `is_auth_portal = request.url.path...` while
                     # leaving route_path correct was NOT caught by the six
                     # hosts above: starlette splits the netloc at the first
                     # "/", so `Host: x/../` yields a path beginning "/.." and
                     # the `/api/auth/` prefix test fails honestly. The shape
                     # that attacks the NEW exemption is a host whose tail IS
                     # the exempt prefix -- `Host: evil.example/api/auth/x`
                     # with `GET /api/runs` gives
                     # url.path='/api/auth/x/api/runs', which starts with
                     # "/api/auth/" while the router still dispatches
                     # /api/runs. That is the original CVE reopened in a new
                     # place, and nothing in the suite looked for it.
                     "evil.example/api/auth/x",
                     "127.0.0.1:8000/api/auth/",
                     "x/api/auth/me",
                     "x/login"):
            response = gated.get("/api/runs", headers={"Host": host})
            assert response.status_code == 401, (
                f"Host: {host!r} skipped the gate")

    def test_the_auth_prefix_cannot_be_smuggled_onto_a_data_route(self,
                                                                 gated):
        """The exemption is a PREFIX, so the interesting attack is making a
        gated path look like it is under that prefix."""
        for path in ("/api/auth/../runs", "/api/auth/%2e%2e/runs",
                     "/api/auth/..%2Fruns"):
            response = gated.get(path)
            assert response.status_code in (401, 404, 307, 405), (
                f"{path} reached a data route unauthenticated "
                f"({response.status_code})")
            if response.status_code == 200:
                pytest.fail(f"{path} returned 200")

    def test_the_gate_is_off_when_no_key_is_configured(self, ungated):
        """Unchanged local behaviour (§0). The login portal must not have
        quietly introduced a gate where there was none."""
        assert ungated.get("/api/runs").status_code == 200
        assert ungated.get("/api/auth/me").json()["authenticated"] is False


# --------------------------------------------------------------------------- #
# §5 — a session partitions nothing
# --------------------------------------------------------------------------- #
class TestNothingIsPartitioned:
    def test_two_accounts_see_the_same_workspace(self, tmp_path, monkeypatch,
                                                 known_hash):
        """§5, measured. The half-built version is the dangerous one: a UI
        that says "your projects" over a shared store teaches a false belief
        about who can see what."""
        app, manager = _isolate(tmp_path, monkeypatch, api_key=API_KEY)
        try:
            one = TestClient(app, client=LOOPBACK)
            two = TestClient(app, client=("127.0.0.2", 51002))
            with one, two:
                service = app.state.auth
                service.store.create("a@x.co", known_hash)
                service.store.create("b@x.co", known_hash)
                assert one.post("/api/auth/login",
                                json={"email": "a@x.co",
                                      "password": PASSWORD}).status_code == 200
                assert two.post("/api/auth/login",
                                json={"email": "b@x.co",
                                      "password": PASSWORD}).status_code == 200
                made = one.post("/api/projects",
                                json={"name": "shared study"})
                assert made.status_code in (200, 201), made.text
                project_id = made.json()["id"]
                theirs = two.get("/api/workspace").json()
                ids = [p["id"] for p in theirs.get("projects", [])]
                assert project_id in ids, (
                    "the second account cannot see the first's project -- the "
                    "stores have been partitioned, which §5 forbids")
        finally:
            manager.shutdown()

    def test_no_payload_claims_the_data_is_yours(self, gated, known_hash):
        """§5: no "your projects", no owner column. Checked on the payloads a
        signed-in client actually reads."""
        _add_account(gated, EMAIL, known_hash)
        _sign_in(gated)
        for path in ("/api/auth/me", "/api/account", "/api/workspace"):
            text = json.dumps(gated.get(path).json()).lower()
            for claim in ("your project", "your run", "your specimen",
                          "my project", "owned by", "owner_id"):
                assert claim not in text, f"{path} says {claim!r}"

    def test_the_account_panel_still_tells_the_truth_about_attribution(
            self, gated, known_hash):
        """The field a UI would read before captioning a run with somebody's
        name. Its REASON changed when sessions landed and its value did
        not."""
        _add_account(gated, EMAIL, known_hash)
        _sign_in(gated)
        body = gated.get("/api/account").json()
        attribution = body["usage"]["attribution"]
        assert attribution.startswith("none")
        assert "shared" in attribution
        session = body["auth"]["session"]
        assert session["supported"] is True
        assert session["signed_in"] is True
        assert session["partitions_data"] is False
        assert body["auth"]["self_service_signup"] is False
        assert "browser-session" in body["auth"]["schemes"]

    def test_the_account_panel_reports_a_session_caller_without_a_key(
            self, gated, known_hash):
        """``matches`` was documented as "reported rather than assumed,
        because the day this route stops being gated is the day the field
        starts mattering". The day arrived: a signed-in browser reaches this
        route with no key at all."""
        _add_account(gated, EMAIL, known_hash)
        _sign_in(gated)
        body = gated.get("/api/account").json()
        assert body["auth"]["matches"] is False
        assert body["auth"]["caller_fingerprint"] is None
        assert body["auth"]["session"]["signed_in_as"] == EMAIL_CASEFOLDED


# --------------------------------------------------------------------------- #
# §7 — what is NOT here
# --------------------------------------------------------------------------- #
def test_no_email_path_exists_anywhere_in_the_module():
    """§7: no verification mail, no reset link, no SMTP config. Half-built
    email is worse than none, and the way it arrives is one "helpful" import.
    """
    path = os.path.join(_REPO, "app", "server", "auth.py")
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    for banned in ("smtplib", "email", "ssl_email", "mailbox"):
        assert banned not in imported, f"auth.py imports {banned}"


def test_no_roles_or_permissions_field_exists(tmp_path, known_hash):
    """§7: a role model with one role is just a column that lies later."""
    store = auth_mod.UserStore(root=str(tmp_path / "auth"))
    row = store.create("a@b.co", known_hash)
    for banned in ("role", "roles", "permissions", "scopes", "is_admin",
                   "admin", "tier", "plan"):
        assert banned not in row, f"the stored user carries a {banned!r} field"


def test_the_contract_is_unedited():
    """This file's specification is FROZEN. If it moved, these assertions are
    measuring the wrong thing and someone needs to know."""
    path = os.path.join(_REPO, "docs", "AUTH_CONTRACT.md")
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    assert "**Frozen 2026-10-09.**" in text
    # The five facts this file is built on, quoted from the contract so a
    # silent edit to either side shows up as a failure here.
    assert "Name `vivum_session`" in text
    assert "`<user_id>.<expiry_unix>.<hmac>`" in text
    assert "scrypt$<n>$<r>$<p>$<salt_b64>$<hash_b64>" in text
    assert "`n=2**15, r=8, p=1`" in text
    assert "Minimum length 12 characters" in text
    assert "must\n  NOT be derived from `TMESIM_STUDIO_API_KEY`" in text
    assert "Compared with `secrets.compare_digest`, always. Never `==`." in text


def test_the_deployment_page_documents_the_new_variable():
    """``tests/test_app_server.py::TestEnvVarsAreDocumented`` already enforces
    presence; this pins that the entry says the two things an operator has to
    know, rather than merely containing the string."""
    path = os.path.join(_REPO, "docs", "DEPLOYMENT.md")
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    assert "### `TMESIM_SESSION_SECRET`" in text
    assert "TMESIM_STUDIO_API_KEY" in text.split(
        "### `TMESIM_SESSION_SECRET`")[1][:3000], (
        "the entry does not say it is not the API key")
    assert "16 characters" in text
