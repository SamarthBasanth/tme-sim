"""What the backend does when a client sends something it did not expect.

Every route in ``app/server/api.py`` (plus the ``/api/protocol`` and
``/api/account`` routers it mounts) driven with what a careless or hostile
client actually sends: no body, the wrong JSON *shape*, the wrong *type* in a
RunSpec field, a filesystem path where a section id belongs, a percent-encoded
``..`` in an identifier, and every route with the API-key gate armed.

**What this file is for.** The three failures it exists to catch are, in
descending order of what they would cost:

1. **A 500.** A refusal is a product surface: it tells a researcher what they
   sent that was wrong. ``Internal Server Error`` tells them nothing, tells the
   operator nothing without console access, and on ``POST /api/validate`` —
   whose entire job is to turn bad input into a readable message — it is a
   contradiction in terms. ``validate_spec`` already carries a comment saying
   exactly that (``app/server/api.py``, the ``lookup_section`` note).
2. **A leak in the refusal.** ``traceback.format_exc()`` carries the repo root
   and the operator's OS username; an unredacted ``str(exc)`` carries the
   module layout. ``tests/test_response_redaction.py`` watches the *source* for
   the pattern; this file watches the *wire*, with a real failed run in the
   manager and a sweep over the routes that report it.
3. **The gate leaking round the edges.** With ``TMESIM_STUDIO_API_KEY`` set,
   every ``/api/*`` route must refuse without the header, and the two liveness
   routes that are deliberately exempt must be incapable of saying anything
   beyond "the process answers".

**What this file deliberately does not do.** It asserts only behaviour that is
correct. Input classes that currently produce a 500, or that are accepted when
they should be refused, are reported rather than pinned here — a test that
asserts a defect makes the defect permanent, and a test that passes without
exercising its subject is worse than no test at all.

**Cost.** No reconstruction, no PDE, no real tissue. The one failed run is
injected (``start_worker=False`` plus a hand-set status, the pattern
``test_app_server.py::TestTracebacksAreNotServedByDefault`` uses), and every
store is rooted in a ``tmp_path`` so run history, projects and compounds never
touch ``app/runs/``, ``app/projects/`` or ``app/compounds/``.
"""
from __future__ import annotations

import json
import os
import sys

import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

fastapi = pytest.importorskip("fastapi", reason="the app/server backend needs fastapi")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi.testclient import TestClient  # noqa: E402

from app.server import cache as cache_mod  # noqa: E402
from app.server import compounds as compounds_mod  # noqa: E402
from app.server import ingest_registry  # noqa: E402
from app.server.api import create_app, default_spec, validate_spec  # noqa: E402
from app.server.cache import repo_root  # noqa: E402
from app.server.jobs import JobManager  # noqa: E402
from app.server.projects import ProjectStore  # noqa: E402

#: Strings that must never appear in a response body. ``repo_root()`` and the
#: interpreter's own prefix are the two absolute paths a stack trace or an
#: unredacted ``str(exc)`` drags in; the rest are the shape of a Python
#: traceback, which is never a thing a browser should be shown.
LEAK_MARKERS = (
    "Traceback (most recent call last)",
    'File "',
    "site-packages",
    ".venv/",
)


def assert_no_leak(resp, *, extra=()):
    """A response body carries no traceback, no host path and no secret."""
    body = resp.text
    for marker in LEAK_MARKERS:
        assert marker not in body, f"{marker!r} in {resp.request.url}: {body[:400]}"
    assert repo_root() not in body, f"repo root in {resp.request.url}: {body[:400]}"
    for secret in extra:
        assert secret not in body, f"secret in {resp.request.url}: {body[:400]}"


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture
def client(tmp_path, monkeypatch):
    """The app with EVERY store rooted in tmp_path and no ambient config.

    Run history, projects, compounds and the ingested-section registry are all
    redirected, so a hostile POST in this file can never register a section in,
    or delete a compound from, the developer's real ``app/`` tree. The three
    environment variables the app reads are cleared, because ambient
    configuration must not decide what these tests measure.
    """
    for var in ("TMESIM_STUDIO_API_KEY", "TMESIM_STUDIO_DATA_ROOTS",
                "TMESIM_STUDIO_EXPOSE_TRACEBACKS",
                "TMESIM_STUDIO_ALLOWED_ORIGINS", "TMESIM_STUDIO_ALLOWED_HOSTS"):
        monkeypatch.delenv(var, raising=False)

    runs = tmp_path / "runs"
    runs.mkdir()
    monkeypatch.setattr(cache_mod, "runs_dir", lambda: str(runs))
    monkeypatch.setattr(ingest_registry, "_store", None)
    monkeypatch.setattr(
        compounds_mod, "COMPOUND_STORE",
        compounds_mod.CompoundStore(root=str(tmp_path / "compounds")))

    manager = JobManager(root=str(runs / "jobs"), start_worker=False)
    store = ProjectStore(root=str(tmp_path / "projects"))
    app = create_app(manager, projects=store)
    with TestClient(app) as c:
        c._manager = manager      # type: ignore[attr-defined]
        c._projects = store       # type: ignore[attr-defined]
        yield c
    manager.shutdown()
    monkeypatch.setattr(ingest_registry, "_store", None)


@pytest.fixture(scope="session", autouse=True)
def _warm_the_drug_library():
    """Compile the drug library once, before the first app is built.

    ``create_app`` starts a background ``_warm_drugs`` thread per app, and this
    file builds one app per test. Cold, they all miss the memo at once and
    dozens of concurrent RDKit compiles contend for the GIL — the file took
    nine minutes, nearly all of it that. Warming the module-level memo first
    makes every one of those threads a no-op.
    """
    from app.server.api import drug_rows

    drug_rows()


@pytest.fixture
def spec():
    """A RunSpec that ``/api/validate`` accepts, to mutate one field at a time."""
    return dict(default_spec())


#: Every route that takes a JSON **object** body. ``/api/validate`` is excluded
#: on purpose: it answers rather than refuses (see its own test below).
OBJECT_BODY_ROUTES = [
    "/api/sections/scan",
    "/api/compounds/compile",
    "/api/compounds",
    "/api/runs",
    "/api/projects",
    "/api/reports",
]


# --------------------------------------------------------------------------- #
# 1. Shape abuse — the body is not the shape the route was written for
# --------------------------------------------------------------------------- #
class TestBodyShape:
    @pytest.mark.parametrize("route", OBJECT_BODY_ROUTES)
    @pytest.mark.parametrize("body,kind", [
        ([1, 2, 3], "list"), ("hello", "str"), (5, "int"), (True, "bool"),
    ])
    def test_a_json_scalar_or_array_where_an_object_belongs_is_a_400(
            self, client, route, body, kind):
        """A body of the wrong JSON *shape* must be a 400, never a 500.

        Every one of these routes reaches into the body with ``.get()``. On a
        list that is ``AttributeError``, on a string it silently succeeds for
        the wrong reason, and on an int it is ``AttributeError`` again — all of
        which surface as ``Internal Server Error``. Operationally: a UI bug
        that posts an array instead of an object would produce a blank error
        toast and a log line nobody can act on, instead of a sentence naming
        the shape that arrived.
        """
        resp = client.post(route, json=body)
        assert resp.status_code == 400, resp.text
        assert_no_leak(resp)

    @pytest.mark.parametrize("route", OBJECT_BODY_ROUTES)
    @pytest.mark.parametrize("raw,kind", [
        (b"", "empty"),
        (b"not json at all", "text"),
        (b'{"unclosed": ', "truncated"),
        (b"[" * 2000 + b"]" * 2000, "2000-deep nesting"),
    ])
    def test_an_unparseable_body_is_a_400_that_says_so(
            self, client, route, raw, kind):
        """Malformed JSON, including a nesting depth that exhausts the decoder's
        own recursion limit, is a 400 and not a crash.

        The 2000-deep case is the one that matters: ``json.loads`` raises
        ``RecursionError``, which is not an exception most handlers think to
        catch. If it escaped, one 4 KB request would return a 500 — and on a
        server where the recursion limit is shared with the worker thread, a
        deliberate near-limit body is a way to make *other* work fail.
        """
        resp = client.post(route, content=raw,
                           headers={"content-type": "application/json"})
        assert resp.status_code == 400, resp.text
        assert_no_leak(resp)

    @pytest.mark.parametrize("route", OBJECT_BODY_ROUTES + ["/api/validate"])
    def test_no_body_at_all_is_a_refusal_not_a_crash(self, client, route):
        """A POST with no body — what a `curl -X POST` typo sends.

        This is the single most common hand-driven mistake against this API and
        it must produce a sentence, not a stack.
        """
        resp = client.post(route)
        assert resp.status_code in (200, 400), resp.text
        assert_no_leak(resp)

    def test_duplicate_json_keys_resolve_to_the_last_and_never_error(
            self, client):
        """RFC 8259 leaves duplicate keys undefined; Python takes the last.

        Worth pinning because a client that builds JSON by string concatenation
        (the ``/ingest`` page's form, a shell script) can emit a duplicate key
        by accident. The requirement is not *which* value wins — it is that the
        request is answered deterministically rather than rejected or crashed.
        """
        resp = client.post("/api/projects",
                           content=b'{"name":"first","name":"second"}',
                           headers={"content-type": "application/json"})
        assert resp.status_code == 201, resp.text
        assert resp.json()["name"] == "second"

    def test_an_unrecognised_run_spec_field_is_a_warning_not_a_refusal(
            self, client, spec):
        """A field this build does not know must WARN and run, not refuse.

        The Studio's spec editor and this server version out of step by one
        release is the normal case, not the exceptional one. Refusing would
        make every forward-compatible client unusable; silently dropping would
        let a user believe a knob took effect. The contract is the third thing:
        ignored, and said out loud.
        """
        spec["not_a_real_knob"] = 7
        body = client.post("/api/validate", json=spec).json()
        assert body["ok"] is True, body["errors"]
        assert any("not_a_real_knob" in w for w in body["warnings"]), \
            body["warnings"]

    def test_post_validate_answers_a_broken_body_rather_than_refusing_it(
            self, client):
        """``/api/validate`` is the code editor's endpoint and must always
        answer 200 with a verdict.

        Every other route refuses malformed JSON with a 400. This one cannot:
        the editor posts on every keystroke, and half-typed JSON is the *normal*
        state of its buffer. A 4xx here makes the editor's error panel go blank
        exactly when the user most needs it to say what is wrong.
        """
        for raw in (b"", b"{", b"not json", b"[" * 2000 + b"]" * 2000):
            resp = client.post("/api/validate", content=raw,
                               headers={"content-type": "application/json"})
            assert resp.status_code == 200, (raw[:40], resp.text)
            body = resp.json()
            assert body["ok"] is False and body["errors"]
            assert_no_leak(resp)

    @pytest.mark.parametrize("body,kind", [
        ([1], "list"), ("x", "str"), (5, "int"), (True, "bool"),
    ])
    def test_post_validate_names_the_shape_that_arrived(self, client, body, kind):
        """The verdict on a non-object RunSpec names the type it got.

        ``validate_spec`` is called directly by the CLI as well as over HTTP,
        so this is the one message a user sees from two surfaces; it has to be
        specific enough to act on without a second request.
        """
        verdict = client.post("/api/validate", json=body).json()
        assert verdict["ok"] is False
        assert any(kind in e for e in verdict["errors"]), verdict["errors"]


# --------------------------------------------------------------------------- #
# 2. Type and range abuse inside a RunSpec
# --------------------------------------------------------------------------- #
class TestRunSpecTypes:
    """``validate_spec`` is the ONLY thing between a JSON body and a solver
    that allocates float64 grids. Every refusal here is a refusal to allocate.
    """

    @pytest.mark.parametrize("field", ["dx_um", "total_time_s"])
    @pytest.mark.parametrize("literal", ["NaN", "Infinity", "-Infinity"])
    def test_a_non_finite_number_is_refused_before_it_reaches_the_solver(
            self, client, field, literal):
        """``NaN`` and ``Infinity`` are not JSON, but Python's decoder accepts
        them and every JS client emits them for an empty numeric input.

        A NaN ``dx_um`` propagates silently: the grid extent becomes NaN, every
        field is NaN, and the run COMPLETES with a kill_fraction of NaN that
        renders as a blank number in the UI. That is the worst outcome
        available — a result that looks like an answer. It has to die at the
        boundary.
        """
        raw = ('{"section":"xenium_319","dx_um":12,"total_time_s":3600,'
               '"drug":"temozolomide","pk":{"plasma_peak":50,'
               '"half_life_s":14400,"n_doses":1,"dose_interval_s":86400}}')
        doc = json.loads(raw)
        body = json.dumps(doc).replace(
            f'"{field}": {doc[field]}', f'"{field}": {literal}')
        resp = client.post("/api/validate", content=body.encode(),
                           headers={"content-type": "application/json"})
        verdict = resp.json()
        assert verdict["ok"] is False, verdict
        assert any(field in e and "finite" in e for e in verdict["errors"]), \
            verdict["errors"]

        run = client.post("/api/runs", content=body.encode(),
                          headers={"content-type": "application/json"})
        assert run.status_code == 400, run.text
        assert_no_leak(run)

    @pytest.mark.parametrize("value,why", [
        ("12", "a string, not a number"),
        (-12.0, "negative"),
        (0, "zero"),
        (True, "a bool, which is an int in Python"),
        ([12], "an array"),
        ({"v": 12}, "an object"),
        (None, "null"),
    ])
    def test_dx_um_must_be_a_finite_positive_number(self, client, spec,
                                                    value, why):
        """Grid memory scales as 1/dx^2, so ``dx_um`` is the field a hostile
        client reaches for first.

        ``"12"`` is the interesting one: a spec editor that does not coerce its
        text input sends the string. Accepting it via a lenient ``float()``
        would mean the type of the field depends on which client sent it, and
        the ``dx <= 0`` guard below would then be comparing a string.
        """
        spec["dx_um"] = value
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"] is False, (why, verdict)
        assert any("dx_um" in e for e in verdict["errors"]), verdict["errors"]
        assert client.post("/api/runs", json=spec).status_code == 400

    def test_a_dx_below_the_half_micron_floor_is_an_error_not_a_warning(
            self, client, spec):
        """``ok = not errors``, and ``POST /api/runs`` refuses only when ``ok``
        is false — so a WARNING here would not block the allocation.

        ``dx_um = 1e-308`` on a 5000 um section asks for ~1e622 voxels. The
        refusal has to be an error, in the ``errors`` list, or the ceiling is
        decorative.
        """
        spec["dx_um"] = 1e-308
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"] is False
        assert any("dx_um" in e for e in verdict["errors"]), verdict["errors"]
        assert client.post("/api/runs", json=spec).status_code == 400

    @pytest.mark.parametrize("value", ["3600", -1, 0, 1e-308, [1], {"v": 1}, True])
    def test_total_time_s_must_be_a_positive_window(self, client, spec, value):
        """A zero or sub-timestep window makes every exposure metric exactly 0.

        Not an error the engine can detect after the fact: 0 is a legal
        ``kill_fraction``. A user reading "this drug kills nothing" from a
        window they set to 0 by accident is a wrong conclusion the product
        handed them, so the window is checked against the solver's own dt at
        the boundary.
        """
        spec["total_time_s"] = value
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"] is False, verdict
        assert any("total_time_s" in e for e in verdict["errors"]), \
            verdict["errors"]
        assert client.post("/api/runs", json=spec).status_code == 400

    @pytest.mark.parametrize("value,why", [
        ("", "empty string"),
        ("a" * 100_000, "100 kB of id"),
        ("../../../../etc/passwd", "a filesystem path"),
        ("xenium_319\x00", "a trailing NUL byte"),
        (["xenium_319"], "an array"),
        ({"id": "xenium_319"}, "an object"),
        (1, "an integer"),
        ("xenium_З19", "a Cyrillic З where a 3 belongs"),
    ])
    def test_section_is_only_ever_resolved_as_a_known_string_id(
            self, client, spec, value, why):
        """``section`` names a registry key. It must never be treated as a path.

        ``lookup_section`` is a dict ``get``, so a list or dict argument raises
        ``TypeError: unhashable type`` — a 500 out of the endpoint whose job is
        to explain bad input. The traversal and NUL cases matter for a
        different reason: the id reaches the reconstruction cache and the
        ingested-section registry, and neither should ever be asked about a
        path. The Cyrillic case is the one a human cannot see: it must be
        refused as unknown rather than silently matched.
        """
        spec["section"] = value
        resp = client.post("/api/validate", json=spec)
        assert resp.status_code == 200, resp.text
        verdict = resp.json()
        assert verdict["ok"] is False, (why, verdict)
        assert any(e.startswith("section:") for e in verdict["errors"]), \
            verdict["errors"]
        assert_no_leak(resp)

        run = client.post("/api/runs", json=spec)
        assert run.status_code == 400, run.text
        assert_no_leak(run)

    @pytest.mark.parametrize("value,why", [
        ("notadrug", "unknown name"),
        ("temozolomide\x00", "a known name with a NUL appended"),
        ("a" * 100_000, "100 kB of name"),
        ({"name": "temozolomide"}, "an object"),
        (["temozolomide"], "an array"),
        (12, "an integer"),
    ])
    def test_an_unusable_drug_is_refused_and_never_guessed(
            self, client, spec, value, why):
        """A drug name selects a calibrated parameter set. A near-miss must not
        resolve.

        The failure this prevents is a fuzzy or prefix match: a run labelled
        ``temozolomide\\x00`` that silently simulated temozolomide would produce
        a result attributed to a compound the user did not name.
        """
        spec["drug"] = value
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"] is False, (why, verdict)
        assert any(e.startswith("drug:") for e in verdict["errors"]), \
            verdict["errors"]
        assert client.post("/api/runs", json=spec).status_code == 400

    def test_a_null_drug_is_a_legal_reconstruction_only_run_that_says_so(
            self, client, spec):
        """``drug: null`` is not an error — it is the reconstruction-only mode.

        Pinned because it is the one null in the spec that must NOT be
        refused, and because the warning that explains it was inverted for a
        month (it fired on every run that HAD a drug and was silent on the one
        case it exists for).
        """
        spec["drug"] = None
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"] is True, verdict["errors"]
        assert any("drug is null" in w for w in verdict["warnings"]), \
            verdict["warnings"]

    @pytest.mark.parametrize("pk,why", [
        ("x", "a string"),
        ([1], "an array"),
        (None, "null"),
        ({"n_doses": "3"}, "a stringly-typed dose count"),
        ({"n_doses": 1.5}, "a fractional dose count"),
        ({"n_doses": 0}, "zero doses"),
        ({"plasma_peak": 50, "half_life_s": 0, "n_doses": 1,
          "dose_interval_s": 1}, "a zero half-life the profile divides by"),
        ({"plasma_peak": 0, "half_life_s": 100, "n_doses": 1,
          "dose_interval_s": 1}, "a zero peak that delivers no drug"),
    ])
    def test_the_pk_block_is_type_checked_field_by_field(
            self, client, spec, pk, why):
        """``plasma_profile`` divides by ``half_life`` and multiplies by
        ``plasma_peak``; neither may arrive as 0 or as a string.

        A zero half-life is a division by zero inside the solver — an inf
        plasma curve, and a run that finishes with a kill_fraction of 1.0 that
        means nothing. A zero peak is the mirror image: 0.0 everywhere, which
        reads as "this drug does not work".
        """
        spec["pk"] = pk
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"] is False, (why, verdict)
        assert any(e.startswith("pk") for e in verdict["errors"]), \
            verdict["errors"]

    @pytest.mark.parametrize("field,value", [
        ("sensitivity_mode", {"a": 1}), ("sensitivity_mode", "nonsense"),
        ("vascular_mode", ["kde"]), ("vascular_mode", 7),
        ("oxygen_mode", 7), ("oxygen_mode", None),
        ("state_mode", {"a": 1}),
        ("resistance_off", "yes"), ("attribute_mechanism", "yes"),
        ("time_resolved", {"a": 1}),
        ("mgmt", [1]), ("mgmt", {"a": 1}), ("mgmt", -5),
        ("timelapse_frames", "x"), ("timelapse_frames", 1.5),
        ("timelapse_frames", -1), ("timelapse_frames", 1),
    ])
    def test_every_enumerated_and_boolean_knob_is_checked_by_value_not_by_truth(
            self, client, spec, field, value):
        """A mode is a closed set and a flag is a bool. Neither may be coerced.

        ``resistance_off`` is the one with teeth: it selects the
        COUNTERFACTUAL arm whose kill_fraction is the denominator of rho, not a
        therapy result. Python would read the string ``"no"`` as true, so a
        client sending a form field verbatim would silently get the
        counterfactual and quote it as a prediction.
        """
        spec[field] = value
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"] is False, verdict
        assert any(field in e for e in verdict["errors"]), verdict["errors"]

    def test_a_timelapse_frame_count_above_the_ceiling_is_capped_not_refused(
            self, client, spec):
        """``timelapse_frames`` far above ``MAX_TIMELAPSE_FRAMES`` is a warning
        and a cap, because the payload writer caps it too.

        This one is deliberately NOT a refusal: the ceiling is enforced where
        the memory is spent (``payloads.encode_timelapse``), so the validator's
        job is to tell the user the number they will actually get. A refusal
        here would be a second, drifting copy of the ceiling.
        """
        spec["timelapse_frames"] = 10 ** 9
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"] is True, verdict["errors"]
        assert any("MAX_TIMELAPSE_FRAMES" in w for w in verdict["warnings"]), \
            verdict["warnings"]

    def test_a_refused_spec_never_names_a_host_path_or_a_module(
            self, client, spec):
        """Sweep: the refusal text is written for a researcher, not a debugger.

        ``validate_spec`` interpolates the caller's own values into its
        messages. What it must never interpolate is anything from the server —
        the repo root, a module name, a stack. This is the wire-side companion
        to ``tests/test_response_redaction.py``, which watches the source.
        """
        for field, value in [("section", "../../etc"), ("dx_um", "x"),
                             ("drug", {"a": 1}), ("pk", "x"),
                             ("total_time_s", -1), ("mgmt", [1])]:
            bad = dict(spec)
            bad[field] = value
            for route in ("/api/validate", "/api/runs"):
                assert_no_leak(client.post(route, json=bad))


# --------------------------------------------------------------------------- #
# 3. POST /api/sections/scan — the one route that takes a filesystem path
# --------------------------------------------------------------------------- #
#: Paths that are not sample directories. Each one is a thing a real client
#: sends by accident, or a probe a hostile one sends on purpose.
HOSTILE_PATHS = [
    ("/dev/zero", "a character device"),
    ("/dev/stdin", "a pipe the reader would block on"),
    ("/proc/self/environ", "the process's own environment"),
    ("/tmp/a\x00b", "an embedded NUL byte"),
    ("/" + "a" * 100_000, "a 100 kB path"),
    ("../" * 4000, "4000 levels of traversal"),
    ("~", "a tilde the shell would have expanded"),
    ("//etc", "a double-slash root"),
    ("file:///etc", "a URL where a path belongs"),
]

#: The subset the UNCONFINED test drives. ``~`` and a 4000-deep ``..`` resolve
#: to the home directory and to ``/``, and an unconfined scan really does walk
#: them — correct behaviour, and not something a test should do to the machine
#: it runs on. Both stay in the CONFINED list above, where they are refused
#: before anything is opened, which is the assertion that matters for them.
UNCONFINED_PROBE_PATHS = [(p, w) for p, w in HOSTILE_PATHS
                          if p not in ("~", "../" * 4000)]


class TestScanPathHandling:
    """``tests/test_scan_confinement.py`` proves the allowlist holds against
    four escape shapes. This class goes where it does not: paths that are not
    directories at all, and the disclosure properties of the refusal.
    """

    @pytest.mark.parametrize("path,why", UNCONFINED_PROBE_PATHS)
    def test_a_path_that_is_not_a_sample_directory_is_an_ANSWER_not_a_crash(
            self, client, path, why):
        """Unconfined, an unreadable or nonsensical path returns a report
        saying so.

        The route's contract is "never 500s on a bad directory: an
        unrecognised or broken sample is an ANSWER". That is load-bearing for
        bring-your-own-tissue: the ``/ingest`` page shows ``detection.issues``
        verbatim, so a user who typed the wrong path must read what was wrong
        with it rather than a blank failure. ``/dev/stdin`` is the case that
        would hang rather than crash if the scanner opened before it stat'd.
        """
        resp = client.post("/api/sections/scan", json={"path": path})
        assert resp.status_code == 200, resp.text
        body = resp.json()
        assert body["detection"]["issues"], body
        # The echoed path is the CALLER'S path resolved, and nothing else: the
        # report must not accrete anything from the server on the way through.
        assert body["detection"]["path"] == os.path.abspath(
            os.path.expanduser(path)), body["detection"]["path"]
        if path.startswith("/"):
            assert_no_leak(resp)
        else:
            # A RELATIVE path resolves against the server's working directory,
            # so the echo necessarily contains it. That is what an UNCONFINED
            # scan is — read access to the whole host — and the confined case
            # below refuses these outright rather than resolving them.
            for marker in LEAK_MARKERS:
                assert marker not in resp.text

    @pytest.mark.parametrize("path,why", HOSTILE_PATHS)
    def test_the_same_paths_are_refused_when_roots_are_configured(
            self, client, tmp_path, monkeypatch, path, why):
        """With the allowlist on, none of them is even attempted.

        The point is that confinement is decided BEFORE anything is opened, so
        a path that would have hung or errored never reaches the scanner. A
        403 here and a 200 above is the whole difference between a laptop
        deployment and one with a field tester holding the key.
        """
        root = tmp_path / "allowed"
        root.mkdir()
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR, str(root))
        resp = client.post("/api/sections/scan", json={"path": path})
        assert resp.status_code == 403, resp.text
        assert resp.json()["detail"]["code"] == "path_outside_data_roots"
        assert_no_leak(resp)

    def test_the_refusal_discloses_neither_the_roots_nor_the_resolved_path(
            self, client, tmp_path, monkeypatch):
        """A caller learns that it was refused. Not what would be accepted, and
        not where its own text landed on this host.

        Both halves are disclosures. Naming the roots hands an unauthorised
        caller the map of the operator's data layout in one request; echoing
        the RESOLVED path turns the route into a ``realpath`` oracle — send
        ``~``, read back the operator's home directory, send ``.``, read back
        the server's working directory. The refusal names the caller's OWN
        text, which is the only string in the exchange the caller already had.
        """
        allowed = tmp_path / "allowed"
        allowed.mkdir()
        other = tmp_path / "second_root"
        other.mkdir()
        monkeypatch.setenv(cache_mod.DATA_ROOTS_VAR,
                           os.pathsep.join([str(allowed), str(other)]))
        monkeypatch.setenv("HOME", str(tmp_path / "home"))

        # Every one of these RESOLVES to something different from its own text,
        # which is what makes the echo testable: if the refusal named the
        # resolved object, the resolved string would appear and the sent one
        # would not.
        for sent in ("~/private", ".", "../elsewhere", "~"):
            resolved = os.path.realpath(os.path.expanduser(sent))
            assert resolved != sent, sent
            resp = client.post("/api/sections/scan", json={"path": sent})
            assert resp.status_code == 403, (sent, resp.text)
            body = resp.text
            assert str(allowed) not in body, sent
            assert str(other) not in body, sent
            assert resolved not in body, sent
            # It IS a diagnostic: it names the caller's text and the control.
            detail = resp.json()["detail"]
            assert detail["path"] == sent
            assert cache_mod.DATA_ROOTS_VAR in detail["message"]

    def test_a_body_with_no_usable_path_is_a_400_naming_what_was_wanted(
            self, client):
        """No ``path`` key, an empty one, or a null is a 400, not a scan of
        the server's working directory.

        The failure mode being closed is a falsy path becoming ``""`` and then
        resolving to CWD — which would make the *absence* of a field a request
        to enumerate the deployment directory.
        """
        for body in ({}, {"path": ""}, {"path": None}, {"path": 0},
                     {"register": True}):
            resp = client.post("/api/sections/scan", json=body)
            assert resp.status_code == 400, (body, resp.text)
            assert "path" in resp.json()["detail"]
            assert_no_leak(resp)

    @pytest.mark.parametrize("section_id", [
        "../../evil", "a\x00b", "../..", "Study/A", "..%2f..", "/etc/passwd",
    ])
    def test_a_hostile_section_id_at_registration_is_sanitised_not_joined(
            self, client, tmp_path, section_id):
        """``register: true`` takes a caller-chosen id. It becomes a registry
        KEY, never a path component.

        If it were joined to a directory, ``../../evil`` would write the
        section record outside the registry. The id is normalised to
        ``[a-z0-9_]`` first, so the traversal is not refused — it is spelled
        away, which is the right answer for a field a user types.
        """
        sample = tmp_path / "sample"
        sample.mkdir()
        before = sorted(p.name for p in tmp_path.iterdir())
        resp = client.post("/api/sections/scan",
                           json={"path": str(sample), "register": True,
                                 "section_id": section_id})
        assert resp.status_code == 200, resp.text
        body = resp.json()
        sid = body["section_id"]
        # The property, not one particular spelling: whatever id survives
        # normalisation cannot address anything outside the registry.
        assert sid, body
        assert "/" not in sid and "\\" not in sid, sid
        assert ".." not in sid and "\x00" not in sid, sid
        assert sid == sid.lower(), sid
        # An empty directory is not a sample: the registration is refused on
        # its merits, and nothing was created anywhere.
        assert body["ok"] is False
        assert sorted(p.name for p in tmp_path.iterdir()) == before, \
            "registration created something outside the registry"
        assert_no_leak(resp)


# --------------------------------------------------------------------------- #
# 4. Identifiers in the URL path
# --------------------------------------------------------------------------- #
#: Every id template in the app that takes a caller-supplied identifier and
#: resolves it against an in-memory table.
IN_MEMORY_ID_ROUTES = [
    ("GET", "/api/runs/{i}"),
    ("GET", "/api/runs/{i}/status"),
    ("GET", "/api/runs/{i}/traceback"),
    ("DELETE", "/api/runs/{i}"),
    ("GET", "/api/runs/{i}/geometry"),
    ("GET", "/api/runs/{i}/mechanism"),
    ("GET", "/api/runs/{i}/timelapse"),
    ("GET", "/api/reports/{i}"),
    ("GET", "/api/reports/{i}/html"),
    ("DELETE", "/api/reports/{i}"),
    ("GET", "/api/protocol/runs/{i}"),
    ("GET", "/api/protocol/sweeps/{i}"),
    ("DELETE", "/api/protocol/sweeps/{i}"),
]

#: Identifiers that a router or a store must not mistake for something else.
#: ``%2e%2e`` survives to the handler as ``..``; the rest are refused by the
#: router itself because a decoded ``/`` breaks single-segment matching.
HOSTILE_IDS = [
    ("%2e%2e", "an encoded .."),
    ("%2e%2e%2f%2e%2e%2fsecret", "an encoded ../../secret"),
    ("%252e%252e%252fsecret", "a double-encoded traversal"),
    ("%2Fetc%2Fpasswd", "an encoded absolute path"),
    ("..%5cwindows", "a backslash traversal"),
    ("a%00b", "an embedded NUL byte"),
    ("a%0d%0aX-Injected:%201", "a CRLF header injection"),
    ("p_0000١", "an Arabic-Indic 1 where an ASCII 1 belongs"),
]


class TestIdentifierHandling:
    @pytest.mark.parametrize("ident,why", HOSTILE_IDS)
    @pytest.mark.parametrize("method,template", IN_MEMORY_ID_ROUTES)
    def test_run_and_report_ids_are_resolved_in_memory_and_never_on_disk(
            self, client, method, template, ident, why):
        """A run or report id is a key in a dict, so an unknown one is a 404
        before any path is built from it.

        This is the design that makes the whole ``/api/runs`` family immune to
        the identifier classes above: ``JobManager.get`` and
        ``ReportManager.get`` are ``dict.get``, and the payload routes call
        ``_require`` *before* ``payload_path``. If a future refactor made an id
        lookup lazy — stat the directory, then load — every entry in this
        parametrisation becomes a filesystem probe with a caller-chosen name.
        """
        resp = client.request(method, template.format(i=ident))
        # 404 from the handler (the id was not in the table) or from the
        # router (the decoded id spans a "/" and matches no single-segment
        # route). A DELETE whose decoded path matches nothing falls through to
        # the catch-all web mount, which serves GET only, so starlette answers
        # 405 -- also a refusal without a handler, which is the property here.
        assert resp.status_code in (404, 405), (resp.status_code, resp.text)
        assert_no_leak(resp)

    @pytest.mark.parametrize("ident,why", HOSTILE_IDS)
    def test_a_hostile_id_never_reaches_a_file_outside_the_report_store(
            self, client, tmp_path, ident, why):
        """The same, proved against a real file planted where ``..`` points.

        A status assertion alone would pass against a store that read the file
        and then happened to 404 for another reason. This plants a
        ``report.json`` one level above the report root — the exact object a
        one-level traversal reaches — and asserts its contents never appear in
        a response and that it survives a DELETE.
        """
        reports_root = client.app.state.reports.root
        canary = os.path.join(os.path.dirname(reports_root), "report.json")
        os.makedirs(os.path.dirname(canary), exist_ok=True)
        with open(canary, "w") as fh:
            fh.write('{"report_id":"UP","status":"done","CANARY":"do-not-serve"}')

        for method, template in [("GET", "/api/reports/{i}"),
                                 ("GET", "/api/reports/{i}/html"),
                                 ("DELETE", "/api/reports/{i}")]:
            resp = client.request(method, template.format(i=ident))
            assert "CANARY" not in resp.text, (method, ident, resp.text[:200])
        assert os.path.exists(canary), f"{ident} deleted a file outside the store"

    @pytest.mark.parametrize("name", [
        "%2e%2e", "%2e%2e%2f%2e%2e%2fetc%2fpasswd", "a%00b", "not_a_field",
        "DRUG_AUC",
    ])
    def test_an_unknown_field_name_is_a_404_that_lists_the_legal_set(
            self, client, name):
        """``/api/runs/{id}/field/{name}`` builds a filename from ``name``.

        The check is an allowlist membership test that runs BEFORE the run is
        even looked up, which is why a traversal in this position can never
        become a path. The uppercase case pins that the allowlist is exact
        rather than case-folded — a near-miss must be a 404 naming the real
        names, not a silent miss that looks like an absent payload.
        """
        resp = client.get(f"/api/runs/r_00001/field/{name}")
        assert resp.status_code == 404, resp.text
        if "%2f" not in name.lower():
            # A single-segment name reaches the handler, and the handler names
            # the legal set. A name that decodes to several segments never
            # matches the route at all, which is a refusal one layer earlier.
            assert "unknown field" in resp.text, resp.text
        assert_no_leak(resp)

    def test_a_binary_payload_route_refuses_an_unknown_run_before_touching_disk(
            self, client):
        """404 for an unknown run on every binary route.

        These four routes are the only ones that return a ``FileResponse``, so
        they are the only place where a caller-supplied string could become a
        served file. The run must be resolved first, always.
        """
        for suffix in ("geometry", "mechanism", "timelapse", "field/drug_auc"):
            resp = client.get(f"/api/runs/no_such_run/{suffix}")
            assert resp.status_code == 404, (suffix, resp.text)
            assert_no_leak(resp)


# --------------------------------------------------------------------------- #
# 5. The gate
# --------------------------------------------------------------------------- #
API_KEY = "hostile-input-test-key-0123456789abcdef"


@pytest.fixture
def gated(tmp_path, monkeypatch):
    """The same app with ``TMESIM_STUDIO_API_KEY`` set.

    A separate fixture because the gate is installed by ``create_app`` from the
    environment at build time; setting the variable on the shared client would
    do nothing.
    """
    for var in ("TMESIM_STUDIO_DATA_ROOTS", "TMESIM_STUDIO_EXPOSE_TRACEBACKS",
                "TMESIM_STUDIO_ALLOWED_ORIGINS", "TMESIM_STUDIO_ALLOWED_HOSTS"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("TMESIM_STUDIO_API_KEY", API_KEY)

    runs = tmp_path / "runs"
    runs.mkdir()
    monkeypatch.setattr(cache_mod, "runs_dir", lambda: str(runs))
    monkeypatch.setattr(ingest_registry, "_store", None)
    monkeypatch.setattr(
        compounds_mod, "COMPOUND_STORE",
        compounds_mod.CompoundStore(root=str(tmp_path / "compounds")))
    manager = JobManager(root=str(runs / "jobs"), start_worker=False)
    app = create_app(manager, projects=ProjectStore(root=str(tmp_path / "p")))
    with TestClient(app) as c:
        yield c
    manager.shutdown()
    monkeypatch.setattr(ingest_registry, "_store", None)


def _concrete_api_urls(app):
    """Every ``/api/*`` route with its path parameters filled in."""
    out = []
    for route in app.routes:
        path = getattr(route, "path", None)
        methods = sorted(getattr(route, "methods", None) or [])
        if not path or not path.startswith("/api/"):
            continue
        url = path
        for name in ("run_id", "section_id", "project_id", "specimen_id",
                     "report_id", "compound_id", "sweep_id", "name"):
            url = url.replace("{" + name + "}", "zz")
        assert "{" not in url, f"unfilled path parameter in {path}"
        for method in methods:
            if method in ("HEAD", "OPTIONS"):
                continue
            out.append((method, url))
    return sorted(out)


class TestTheGate:
    def test_every_api_route_refuses_without_the_header(self, gated):
        """Enumerated from ``app.routes``, so a NEW route joins this test the
        day it is added.

        A hand-written list is the failure mode: it is correct on the day it is
        written and silently incomplete forever after. The gate itself is
        fail-closed by design, and this is the assertion that keeps it honest —
        one route added outside ``/api/`` that mutates state, or one
        ``include_router`` with the wrong prefix, and this goes red.
        """
        urls = _concrete_api_urls(gated.app)
        assert len(urls) > 40, f"only found {len(urls)} routes to check"
        ungated = []
        for method, url in urls:
            resp = gated.request(
                method, url, json={} if method in ("POST", "PATCH", "PUT") else None)
            if resp.status_code != 401:
                ungated.append((method, url, resp.status_code))
        # /api/health is the ONE documented exemption, and it withholds in the
        # handler instead -- see the next test.
        assert [u[:2] for u in ungated] == [("GET", "/api/health")], ungated

    def test_the_liveness_routes_tell_an_anonymous_caller_only_that_it_answers(
            self, gated):
        """``/healthz`` and an anonymous ``GET /api/health`` return exactly
        ``{"ok": true}`` — no runs directory, no counts, no cache statistics.

        Both are exempt from the gate because a platform health checker cannot
        send a header, and pointing one at a 401 gets a healthy container
        killed. The exemption is only safe while the BODY is empty of
        information: the runs directory is an absolute host path, and the run
        counts and cache statistics are a usage profile of the deployment. If
        this test goes red, an unauthenticated stranger is reading operational
        detail off a route that exists to say "yes".
        """
        for url in ("/healthz", "/api/health"):
            resp = gated.get(url)
            assert resp.status_code == 200, resp.text
            assert resp.json() == {"ok": True}, resp.text
            assert_no_leak(resp, extra=(API_KEY,))

    def test_the_diagnostics_on_api_health_still_require_the_key(self, gated):
        """With the key, ``/api/health`` reports what an operator needs.

        The other half of the contract above: the route is exempt from the
        GATE, the DATA is not. A version of this that withheld unconditionally
        would be safe and useless.
        """
        resp = gated.get("/api/health", headers={"x-api-key": API_KEY})
        assert resp.status_code == 200
        body = resp.json()
        for key in ("runs_dir", "runs", "projects_dir", "reconstruction_cache"):
            assert key in body, body

    @pytest.mark.parametrize("header,why", [
        ({}, "no header at all"),
        ({"x-api-key": ""}, "an empty header"),
        ({"x-api-key": "wrong"}, "a wrong key"),
        ({"x-api-key": API_KEY[:-1]}, "the key with its last byte removed"),
        ({"x-api-key": API_KEY + "x"}, "the key with a byte appended"),
        ({"x-api-key": API_KEY.upper()}, "the key case-folded"),
    ])
    def test_a_near_miss_key_is_refused_exactly_like_no_key(
            self, gated, header, why):
        """A prefix of the key must not get further than a wrong one.

        The comparison is ``secrets.compare_digest`` precisely so that the
        truncated-key case cannot be distinguished by timing; this test pins
        the observable half of that — the response is byte-identical, so a
        caller cannot binary-search the key one character at a time by reading
        status codes or messages.
        """
        resp = gated.get("/api/runs", headers=header)
        assert resp.status_code == 401, resp.text
        assert resp.json() == {"detail": "missing or incorrect X-API-Key header"}
        assert_no_leak(resp, extra=(API_KEY,))

    def test_a_host_header_that_smuggles_a_path_does_not_skip_the_gate(
            self, gated):
        """GHSA-86qp-5c8j-p5mr, which was live in this file's history.

        ``request.url`` is rebuilt by concatenating the CLIENT-SUPPLIED Host
        header with the path and re-parsing, so ``Host: 127.0.0.1:8000/../``
        yields ``url.path='/..//api/runs'`` — the ``startswith('/api/')`` test
        fails, the gate is skipped, and the router still dispatches the
        endpoint. The fix was to decide on ``get_route_path(request.scope)``,
        the same string the router matches on, so the gate and the dispatch
        cannot physically disagree.
        """
        for host in ("127.0.0.1:8000/../", "evil.example/..",
                     "x/../", "127.0.0.1:8000/..%2f"):
            resp = gated.get("/api/runs", headers={"Host": host})
            assert resp.status_code == 401, (host, resp.status_code, resp.text)
            assert_no_leak(resp, extra=(API_KEY,))

    def test_no_refusal_on_any_route_ever_echoes_the_key(self, gated):
        """The server must never put its own configured key in a body.

        The sweep is over every ``/api/*`` route with NO header, so every body
        is a refusal the middleware wrote. A refusal that quoted the expected
        value — "expected ``abc123``, got ``''``" is a real shape of mistake —
        would hand the key to the caller who did not have it, which is the
        whole gate given away in one 401.
        """
        for method, url in _concrete_api_urls(gated.app):
            resp = gated.request(
                method, url,
                json={} if method in ("POST", "PATCH", "PUT") else None)
            assert API_KEY not in resp.text, (method, url, resp.text[:300])

    def test_the_key_is_reported_only_as_a_non_reversing_fingerprint(
            self, gated):
        """``/api/account`` identifies the key without transmitting it.

        An operator with two saved keys needs to know WHICH one this server
        holds; that is a legitimate need and the reason the field exists. It is
        met with an HMAC prefix, not the secret. If this went red, the account
        page — the one screen an operator is most likely to screenshot into a
        support thread — would be printing the shared secret.
        """
        resp = gated.get("/api/account", headers={"x-api-key": API_KEY})
        assert resp.status_code == 200, resp.text
        assert API_KEY not in resp.text
        body = resp.json()
        fp = body["auth"]["key_fingerprint"]
        assert fp and fp not in API_KEY and API_KEY not in fp
        assert len(fp) < len(API_KEY)

    def test_the_ungated_documentation_routes_disclose_no_operational_detail(
            self, gated):
        """``/openapi.json``, ``/docs``, ``/redoc`` and ``/ingest`` answer
        without the key. None of them may carry state or a secret.

        They are outside ``/api/`` and serve GET, so the gate's fail-closed
        rule lets them through by design — the front end has to be able to
        load. What that exemption must never be worth is a host path, a run
        count, the runs directory or the key: a schema is a map of the API, and
        an ungated map is a decision; an ungated *deployment state* is a leak.
        """
        for url in ("/openapi.json", "/docs", "/redoc", "/ingest"):
            resp = gated.get(url)
            assert resp.status_code == 200, (url, resp.status_code)
            assert API_KEY not in resp.text, url
            assert repo_root() not in resp.text, url
            for marker in ("runs_dir", "projects_dir", "reconstruction_cache"):
                assert marker not in resp.text, (url, marker)


# --------------------------------------------------------------------------- #
# 6. The error surface
# --------------------------------------------------------------------------- #
class TestErrorSurface:
    def _failed_run(self, client):
        """A run in the ``failed`` state with a traceback that carries the repo
        root — injected, so this costs nothing and the stack is deterministic.
        """
        job = client._manager.submit(dict(default_spec()))
        job.status = "failed"
        job.error = "ValueError: boom"
        job.traceback = (
            'Traceback (most recent call last):\n'
            f'  File "{repo_root()}/tmesim/x.py", line 1, in solve\n'
            'ValueError: boom\n')
        return job

    def test_a_failed_run_reports_its_diagnosis_without_its_stack(self, client):
        """Every surface that mentions a failed run stays clean by default.

        ``traceback.format_exc()`` always carries absolute host paths — the
        repo root, the operator's OS username and the internal module layout.
        The typed diagnosis is the product; the stack is a console artefact.
        The sweep is the point: it is not enough for the dedicated traceback
        route to withhold, because ``result_view`` and ``status_view`` are
        built from the same job object.
        """
        job = self._failed_run(client)
        for url in (f"/api/runs/{job.run_id}",
                    f"/api/runs/{job.run_id}/status",
                    f"/api/runs/{job.run_id}/traceback",
                    f"/api/runs/{job.run_id}/geometry",
                    f"/api/runs/{job.run_id}/mechanism",
                    f"/api/runs/{job.run_id}/timelapse",
                    f"/api/runs/{job.run_id}/field/drug_auc",
                    "/api/runs"):
            resp = client.get(url)
            assert_no_leak(resp)

        body = client.get(f"/api/runs/{job.run_id}/traceback").json()
        assert body["traceback"] is None
        assert body["error"] == "ValueError: boom"
        assert "TMESIM_STUDIO_EXPOSE_TRACEBACKS" in body["traceback_withheld"]

    def test_expose_tracebacks_changes_exactly_one_route(
            self, client, monkeypatch):
        """The opt-in is scoped to ``/api/runs/{id}/traceback`` and nothing else.

        An environment variable that widens the disclosure surface has to have
        a knowable blast radius. An operator who sets this to debug one failing
        run must not discover afterwards that ``/api/runs``, the report
        document or ``/api/health`` also started carrying stacks. The
        assertion is the *difference*: every other body is byte-identical
        before and after.
        """
        job = self._failed_run(client)
        others = [f"/api/runs/{job.run_id}", f"/api/runs/{job.run_id}/status",
                  "/api/runs", "/api/health", "/api/sections", "/api/limits",
                  "/api/reports"]
        before = {u: client.get(u).text for u in others}

        monkeypatch.setenv("TMESIM_STUDIO_EXPOSE_TRACEBACKS", "1")
        for url in others:
            assert client.get(url).text == before[url], url
            assert_no_leak(client.get(url))

        served = client.get(f"/api/runs/{job.run_id}/traceback").json()["traceback"]
        assert "ValueError: boom" in served

    def test_the_served_stack_has_the_repo_root_replaced(
            self, client, monkeypatch):
        """Even opted in, the absolute repo path is scrubbed to ``<repo>``.

        The opt-in buys the operator the shape of the failure — which module,
        which line, which exception. It does not buy a stranger holding the key
        the deployment's directory layout or the operator's username, which is
        what a raw ``format_exc`` would hand over.
        """
        job = self._failed_run(client)
        monkeypatch.setenv("TMESIM_STUDIO_EXPOSE_TRACEBACKS", "1")
        resp = client.get(f"/api/runs/{job.run_id}/traceback")
        assert "<repo>/tmesim/x.py" in resp.json()["traceback"]
        assert repo_root() not in resp.text

    def test_a_404_names_the_id_it_was_given_and_nothing_else(self, client):
        """A "not found" must not describe where it looked.

        ``no such run 'x'`` is a diagnostic. ``no such run:
        /Users/…/app/runs/x/result.json does not exist`` is a directory
        listing, and it is the single easiest way for a store to start leaking
        its layout one 404 at a time.
        """
        for url in ("/api/runs/nope", "/api/reports/nope",
                    "/api/protocol/runs/nope", "/api/compounds/nope"):
            resp = client.get(url)
            assert resp.status_code == 404, (url, resp.text)
            assert_no_leak(resp)

    def test_a_refusal_from_every_post_route_stays_clean(self, client):
        """The whole refusal surface in one sweep, over every shape class.

        This is the regression net for the class of change that reintroduces
        ``f"{type(exc).__name__}: {exc}"`` on a client-facing path — the house
        pattern ``tests/test_response_redaction.py`` was written to watch. That
        file reads the source; this one reads the wire, so a leak that arrives
        through a helper it does not scan still gets caught.
        """
        bodies = [
            {}, [1], "x", 5, True,
            {"path": "/nonexistent-\x00-path"},
            {"name": None}, {"section": None, "drug": None},
            {"section": "../../etc/passwd", "drug": "temozolomide"},
            {"smiles": "not-a-smiles", "name": "x"},
            {"source": "\x00\x01\x02"},
        ]
        routes = OBJECT_BODY_ROUTES + [
            "/api/validate", "/api/protocol/validate", "/api/protocol/plan"]
        for route in routes:
            for body in bodies:
                resp = client.post(route, json=body)
                assert resp.status_code < 500, (route, body, resp.text[:200])
                assert_no_leak(resp)


# --------------------------------------------------------------------------- #
# Regression guards for the defects this file's first pass FOUND
# --------------------------------------------------------------------------- #
#
# The first version of this file asserted survivors only, and deliberately left
# the defects unpinned rather than freeze the broken behaviour. They are fixed
# now, so these pin the fix.
class TestTraversalInStoreIdentifiers:
    """`%2e%2e` in `project_id` / `compound_id` was a read AND delete primitive.

    Starlette percent-decodes before route matching, so `%2e%2e` arrives as a
    literal `..` inside ONE path segment and the router's `[^/]+` never sees a
    slash. Joined onto the store root it resolved one level up. Measured with a
    planted canary: GET served the file above the root, and DELETE removed it.

    Bounded to one level and to the fixed filename `project.json` /
    `compound.json`, and on a real deployment those do not exist — so the blast
    radius on the day it was found was nil. That is not the same as safe: it was
    a delete primitive one root-reconfiguration from mattering, and both stores
    already DEFINED the regex that stops it and applied it only when listing.

    The canary is the point of this test. A status-code assertion passes if the
    read is refused for any reason at all; only a file that is still on disk
    afterwards proves the delete did not happen.
    """

    def _client_with_canary(self, tmp_path):
        from app.server.projects import ProjectStore

        root = tmp_path / "projects"
        root.mkdir()
        canary = tmp_path / "project.json"          # ONE LEVEL ABOVE the root
        canary.write_text(json.dumps({"id": "CANARY", "name": "unreachable"}))
        app = create_app(JobManager(root=str(tmp_path / "runs")),
                         projects=ProjectStore(root=str(root)))
        return TestClient(app), canary

    @pytest.mark.parametrize("method,path", [
        ("GET", "/api/projects/%2e%2e"),
        ("GET", "/api/projects/%2e%2e/compare"),
        ("DELETE", "/api/projects/%2e%2e"),
        ("GET", "/api/projects/abc%00def"),
    ])
    def test_the_canary_above_the_store_root_is_untouched(
            self, tmp_path, method, path):
        c, canary = self._client_with_canary(tmp_path)
        r = c.request(method, path)
        assert r.status_code == 404, (
            f"{method} {path} -> {r.status_code}; a malformed id is a lookup "
            "miss, not a traversal and not a crash")
        assert "CANARY" not in r.text, "the file above the store root was served"
        assert canary.exists(), "the file above the store root was DELETED"

    @pytest.mark.parametrize("method,path", [
        ("GET", "/api/compounds/%2e%2e"),
        ("DELETE", "/api/compounds/%2e%2e"),
        ("GET", "/api/compounds/abc%00def"),
    ])
    def test_the_compound_store_refuses_the_same_shapes(
            self, tmp_path, method, path):
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        assert c.request(method, path).status_code == 404

    def test_a_wellformed_id_that_is_simply_absent_is_also_404(self, tmp_path):
        """The refusal must not be reachable ONLY through the regex — an
        ordinary missing project has to keep behaving."""
        c, _ = self._client_with_canary(tmp_path)
        assert c.get("/api/projects/p_99999").status_code == 404


class TestOversizedIntegersAreRefusedNotCrashed:
    """A JSON integer above `sys.float_info.max` was a 500 in three sites.

    Python parses it exactly, as an `int`, so every `isinstance(v, (int, float))`
    check passed and `float(v)` was where it raised. The boundary is 309 digits;
    308 was a clean 200. It reached `_as_number`, the inline `float(nd) !=
    int(nd)` wholeness checks, and `_num` — and patching only the first moved
    the crash one line rather than removing it.

    A 500 out of `/api/validate` is the failure that endpoint's own docstring
    argues against: its whole job is turning bad input into a readable message.
    """

    BIG = int("9" * 309)

    def _spec(self, **over):
        spec = {"section": "xenium_319", "dx_um": 12.0, "total_time_s": 3600.0,
                "drug": "temozolomide",
                "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
                       "n_doses": 1, "dose_interval_s": 86400.0}}
        pk = over.pop("pk_field", None)
        if pk:
            spec["pk"][pk] = self.BIG
        spec.update(over)
        return spec

    @pytest.mark.parametrize("desc,over", [
        ("dx_um", {"dx_um": BIG}),
        ("total_time_s", {"total_time_s": BIG}),
        ("mgmt", {"mgmt": BIG}),
        ("pk.plasma_peak", {"pk_field": "plasma_peak"}),
        ("pk.half_life_s", {"pk_field": "half_life_s"}),
        ("pk.n_doses", {"pk_field": "n_doses"}),
        ("pk.dose_interval_s", {"pk_field": "dose_interval_s"}),
    ])
    def test_validate_answers_rather_than_crashing(self, tmp_path, desc, over):
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        r = c.post("/api/validate", json=self._spec(**over))
        assert r.status_code == 200, f"{desc}: {r.status_code}, body {r.text[:200]}"
        assert r.json()["ok"] is False, f"{desc} was accepted"

    def test_an_oversized_timelapse_frames_is_clamped_not_refused(
            self, tmp_path):
        """`timelapse_frames` is the one field here that must NOT become an
        error, and it is worth saying why rather than leaving it off the list.

        Above MAX_TIMELAPSE_FRAMES the value is CLAMPED, with a warning, and
        the run reports the count it actually used — bounded by construction,
        so there is no allocation to refuse. That is a different and better
        answer than `pk.n_doses`, which was HONOURED at any size. The only
        thing this field owed was not to crash on an integer float() cannot
        take.
        """
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        r = c.post("/api/validate", json=self._spec(timelapse_frames=self.BIG))
        assert r.status_code == 200
        body = r.json()
        assert body["ok"] is True, "clamping is not a refusal"
        assert any("capped" in w for w in body["warnings"]), \
            "clamping silently would be worse than refusing"

    def test_the_refusal_never_prints_the_integer(self, tmp_path):
        """`str()` of an int over ~4,300 digits raises on Python 3.11+, so a
        message that repr'd the value would re-create the crash inside the
        handler for it."""
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        huge = int("9" * 5000)
        r = c.post("/api/validate", json=self._spec(dx_um=huge))
        assert r.status_code == 200
        assert "9999999999" not in r.text

    def test_a_merely_large_float_is_still_handled(self, tmp_path):
        """308 digits fits in a float and must keep its existing behaviour —
        the fix is for the values that cannot be converted, not a new ceiling."""
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        r = c.post("/api/validate", json=self._spec(dx_um=float("1e308")))
        assert r.status_code == 200


class TestUnboundedWorkIsRefused:
    """`pk.n_doses > 100` was a WARNING, and `ok = not errors`.

    `plasma_profile` allocates `np.arange(n_doses)` plus broadcast arrays of
    shape `(n_doses, len(times))`. `POST /api/runs` accepted `n_doses = 10**15`
    with a 202. Same class as the OOM the grid ceiling exists for: one ordinary
    POST that takes the process down for everyone. The grid was capped; the dose
    count was not.
    """

    def _spec(self, **over):
        spec = {"section": "xenium_319", "dx_um": 30.0, "total_time_s": 3600.0,
                "drug": "temozolomide",
                "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
                       "n_doses": 1, "dose_interval_s": 86400.0}}
        for k, v in over.items():
            if k.startswith("pk."):
                spec["pk"][k[3:]] = v
            else:
                spec[k] = v
        return spec

    @pytest.mark.parametrize("desc,over", [
        ("n_doses 1e9", {"pk.n_doses": 10 ** 9}),
        ("n_doses 1e15", {"pk.n_doses": 10 ** 15}),
        ("total_time_s 1e308", {"total_time_s": 1e308}),
    ])
    def test_runs_refuses_rather_than_accepting(self, tmp_path, desc, over):
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        assert c.post("/api/runs", json=self._spec(**over)).status_code == 400, desc

    def test_a_real_schedule_is_untouched(self, tmp_path):
        """A full Stupp course is 30-42 administrations. The ceiling must be
        nowhere near anything a user would legitimately ask for."""
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        v = c.post("/api/validate", json=self._spec(**{"pk.n_doses": 42}))
        assert v.status_code == 200 and v.json()["ok"] is True

    def test_the_ceilings_are_published(self, tmp_path):
        """A client that cannot see a limit finds it by tripping it — the same
        argument that put max_grid_cells on this route."""
        import os as _os

        key = "limits-key"
        _os.environ["TMESIM_STUDIO_API_KEY"] = key
        try:
            c = TestClient(create_app(JobManager(root=str(tmp_path))))
            lim = c.get("/api/limits", headers={"X-API-Key": key}).json()
            assert lim["max_n_doses"] > 0
            assert lim["max_total_time_s"] > 0
        finally:
            _os.environ.pop("TMESIM_STUDIO_API_KEY", None)


class TestProtocolLimitArgumentsAreValidated:
    """`int(body["max_runs"])` and `float(body["max_window_s"])`, bare, ABOVE
    the try/except that wraps make_plan — so the protocol module's own "never a
    bare traceback" contract did not cover its own argument parsing."""

    @pytest.mark.parametrize("body", [
        {"source": "", "max_runs": "abc"},
        {"source": "", "max_runs": [1]},
        {"source": "", "max_runs": {"a": 1}},
        {"source": "", "max_window_s": "abc"},
        {"source": "", "max_window_s": {"a": 1}},
        {"source": "", "max_runs": int("9" * 400)},
        {"source": "", "max_window_s": int("9" * 400)},
    ])
    def test_neither_route_500s(self, tmp_path, body):
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        assert c.post("/api/protocol/plan", json=body).status_code < 500
        assert c.post("/api/protocol/run", json=body).status_code < 500


class TestScanPathMustBeAString:
    """`if not body.get("path")` is a truthiness test; the next line was
    `str(body["path"])`. A dict became the literal path `{'a': 1}`, resolved
    against the server's working directory and echoed back — a type error
    turned into a disclosure of the server's CWD."""

    @pytest.mark.parametrize("value", [{"a": 1}, [1, 2], 5, True, "   "])
    def test_a_non_string_path_is_a_400(self, tmp_path, value):
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        r = c.post("/api/sections/scan", json={"path": value})
        assert r.status_code == 400
        assert os.getcwd() not in r.text


class TestLivenessAnswersHeadProbes:
    """`HEAD /healthz` returned 404 while the gate already exempted HEAD.

    FastAPI registers only the method you name, so `@app.get` left HEAD
    unrouted — and `api.py`'s key middleware had been written as
    `request.method in ("GET", "HEAD")` on this path, i.e. the gate expected
    HEAD probes the router then refused. A platform health check configured for
    HEAD would have marked the service down while it was serving perfectly.

    That is the exact failure this route exists because of: a health check
    pointed at something it could not satisfy killed a healthy container here
    twice. Nothing in this repo probes with HEAD, which is why it took an input
    sweep to notice.
    """

    @pytest.mark.parametrize("method", ["GET", "HEAD"])
    @pytest.mark.parametrize("path", ["/healthz", "/api/health"])
    def test_a_probe_gets_200_with_or_without_a_key(
            self, tmp_path, monkeypatch, method, path):
        for key in (None, "probe-key"):
            if key is None:
                monkeypatch.delenv("TMESIM_STUDIO_API_KEY", raising=False)
            else:
                monkeypatch.setenv("TMESIM_STUDIO_API_KEY", key)
            c = TestClient(create_app(JobManager(root=str(tmp_path / str(key)))))
            r = c.request(method, path)
            assert r.status_code == 200, (
                f"{method} {path} with key={key!r} -> {r.status_code}; a "
                "platform probe cannot send a header and cannot retry its way "
                "out of a 404")

    def test_head_does_not_widen_what_an_anonymous_caller_learns(
            self, tmp_path, monkeypatch):
        """Adding a method must not add a disclosure. GET already withholds the
        runs directory, counts and cache stats from an unauthenticated caller;
        HEAD returns no body at all, so the property to check is that the GET
        body is still the bare liveness answer."""
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "probe-key")
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        body = c.get("/api/health").json()
        assert body == {"ok": True}
        assert "runs_dir" not in body and "reconstruction_cache" not in body


class TestTheBodyCapCoversEveryRouter:
    """`protocol_api.py` capped bodies at 4 MB. Nothing else did.

    That module's own docstring carries the argument — "a limit that fires
    after the allocation is not a limit" — written after a ~500 MB POST cost
    ~2.8 GB resident and OOM'd the process. The limit was then applied to one
    router out of two.

    MEASURED before the fix: a 64 MB body was buffered AND parsed by
    /api/validate, /api/runs, /api/reports, /api/projects, /api/compounds and
    /api/sections/scan, in 0.07-0.45 s each, while /api/protocol/plan refused
    the same body in 0.009 s.
    """

    ROUTES = ("/api/validate", "/api/runs", "/api/reports", "/api/projects",
              "/api/compounds", "/api/sections/scan", "/api/protocol/plan")

    @pytest.mark.parametrize("route", ROUTES)
    def test_an_oversized_body_is_refused_by_every_post_route(
            self, tmp_path, route):
        from app.server.api import MAX_REQUEST_BODY_BYTES

        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        big = b'{"x":"' + b"a" * (MAX_REQUEST_BODY_BYTES + 1024) + b'"}'
        r = c.post(route, content=big,
                   headers={"content-type": "application/json"})
        assert r.status_code == 413, (
            f"{route} -> {r.status_code}; it buffered a "
            f"{len(big)} B body instead of refusing it")

    def test_a_body_under_the_cap_is_untouched(self, tmp_path):
        """The cap must not become a second, quieter validator."""
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        spec = {"section": "xenium_319", "dx_um": 30.0, "total_time_s": 3600.0,
                "drug": "temozolomide",
                "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
                       "n_doses": 1, "dose_interval_s": 86400.0}}
        r = c.post("/api/validate", json=spec)
        assert r.status_code == 200 and r.json()["ok"] is True

    def test_the_gate_still_answers_first(self, tmp_path, monkeypatch):
        """Middleware ORDER, checked rather than assumed.

        If the cap ran ahead of the key gate, an anonymous caller would learn
        the body limit before authenticating — a small disclosure, and an
        avoidable one. It runs behind: an unauthenticated request is a 401
        whatever its size, and only a caller holding the key sees the 413.

        Nothing is buffered either way. The gate reads headers only, so an
        oversized anonymous body is refused without being read.
        """
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "gate-key")
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        big = b'{"x":"' + b"a" * (5 * 1024 * 1024) + b'"}'
        hdr = {"content-type": "application/json"}
        assert c.post("/api/validate", content=big, headers=hdr
                      ).status_code == 401
        assert c.post("/api/validate", json={"a": 1}).status_code == 401
        assert c.post("/api/validate", content=big,
                      headers={**hdr, "X-API-Key": "gate-key"}
                      ).status_code == 413

    def test_a_liveness_probe_is_never_refused_by_the_cap(self, tmp_path):
        """GET/HEAD are skipped rather than trusted. Refusing a probe would be
        a new way to fail a health check, which this file has done twice."""
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        for method in ("GET", "HEAD"):
            assert c.request(method, "/healthz").status_code == 200


class TestPathNormalisationCannotReachAGatedRoute:
    """The gate and the router must resolve the SAME path, always.

    `api.py` already defends the Host-header variant of this (GHSA-86qp-5c8j-
    p5mr) by deciding on `get_route_path(request.scope)` — the string
    `starlette.routing.Route.matches()` uses — rather than `request.url.path`,
    which starlette rebuilds from the client-supplied Host header. That fix is
    what makes gate and dispatch physically unable to disagree.

    This pins the property from the other direction: no amount of `..`, `.`,
    doubled slashes or case folding may let an anonymous caller REACH a data
    route. Checked by content, not only status — a 200 that returns the
    Studio's `index.html` is the static mount doing its job, and a 200 that
    returns a JSON array is a breach.

    `/api/projects/../..` DOES return 200 with `index.html`, and that is
    correct: after normalisation it is a request for a static asset, and static
    assets are deliberately ungated so the page can load and prompt for a key.
    Recorded here because it looks alarming in a log and is not.
    """

    PROBES = ["/x/../api/projects", "/css/../api/projects", "/./api/projects",
              "/api/./projects", "/api//projects", "/a/b/../../api/projects",
              "/%2e/api/projects", "/api/projects/./", "/API/PROJECTS",
              "/api/projects/../..", "/api/%2e%2e/projects"]

    @pytest.mark.parametrize("probe", PROBES)
    def test_no_probe_reaches_the_project_store(
            self, tmp_path, monkeypatch, probe):
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "gate-key")
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        r = c.get(probe)
        if r.status_code != 200:
            return                      # 401 or 404: never reached a handler
        body = r.text.lstrip()
        assert not body.startswith(("[", "{")), (
            f"{probe} -> 200 with a JSON body: path normalisation reached a "
            f"gated route without the key. Body: {body[:120]}")
        assert body.startswith("<!doctype"), (
            f"{probe} -> 200 with neither JSON nor the SPA shell: {body[:120]}")

    def test_the_control_still_works(self, tmp_path, monkeypatch):
        """If the key stopped working the test above would pass vacuously —
        every probe would 401 for the wrong reason."""
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "gate-key")
        c = TestClient(create_app(JobManager(root=str(tmp_path))))
        r = c.get("/api/projects", headers={"X-API-Key": "gate-key"})
        assert r.status_code == 200
        assert r.text.lstrip().startswith("["), \
            "the control did not reach the store, so the probes prove nothing"
