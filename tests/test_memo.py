"""Guards on ``app/server/memo.py`` -- the exact memo and the exact algebra.

The thing being guarded is not "the cache works". It is that the cache CANNOT
serve a number the current engine would not produce. A cache that is merely
usually right is worse than no cache here: this repository asserts bit-exact
stored references and quotes figures to ten significant figures, so a silently
stale hit would corrupt that and leave nothing behind to find it by.

Every class below names the specific way a wrong number could reach a reader if
the property failed, and every one of them was CONFIRMED to fail against a
planted defect before being kept. The plants, and what each turned red:

* the solver constants dropped from the key composition
  -> ``TestKeyComposition::test_a_solver_constant_changes_the_key``
* the engine fingerprint dropped from the key composition
  -> ``TestKeyComposition::test_the_engine_fingerprint_changes_the_key``,
     ``TestStaleness::test_a_changed_engine_turns_a_hit_into_a_miss``
* ``canonical`` dropping null values, so two distinct specs collide
  -> ``TestCanonical::test_distinct_specs_do_not_collide``
* ``n_doses`` added to ``SCALED_AXES``, so an unmeasured knob gets scaled
  -> ``TestScalingRefusals::test_n_doses_is_never_scaled``
* ``kill_fraction`` added to ``READOUT_RULES``
  -> ``TestScalingRefusals::test_kill_fraction_is_never_served``
* ``payload_status`` trusting the record instead of the filesystem
  -> ``TestPayloadHonesty::test_a_pruned_payload_is_not_present``,
     ``TestPayloadHonesty::test_a_bare_rm_is_reported_as_missing``

Cost: the classes that need a stored run share ONE module-scoped fixture that
solves three specs on the cheapest cached section at ``dx = 30 um`` -- a
reconstruction and three transport solves, a few seconds. Everything else
parses or hashes and touches no tissue.
"""
from __future__ import annotations

import json
import math
import os
import time

import pytest

from app.server import cache as cache_mod
from app.server import jobs as jobs_mod
from app.server import memo as memo_mod
from app.server.jobs import JobManager
from app.server.memo import (
    READOUT_RULES,
    SCALED_AXES,
    Unmemoisable,
    canonical,
    memo_lookup,
    payload_status,
    plasma_scalars,
    run_keys,
    scale_lookup,
    spec_digest,
)

WINDOW_S = 3600.0
DX_UM = 30.0


def _cheapest_section():
    """The first cached section, preferring a small Visium one.

    Deliberately not a fixed id: this suite must run wherever ANY tissue is
    cached, and nothing it asserts is a property of a particular specimen.
    """
    available = [s["id"] for s in cache_mod.list_sections() if s["available"]]
    if not available:
        return None
    visium = [s for s in available if s.startswith("gbmspace_")]
    return (visium or available)[0]


SECTION = _cheapest_section()
requires_section = pytest.mark.skipif(
    SECTION is None, reason="no section is cached in data_cache/")


def _spec(**overrides):
    spec = {
        "section": SECTION,
        "dx_um": DX_UM,
        "drug": "temozolomide",
        "sensitivity_mode": "uniform",
        "vascular_mode": "discrete",
        "oxygen_mode": "pde",
        "state_mode": "oxygen",
        "mgmt": None,
        "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
               "n_doses": 1, "dose_interval_s": 86400.0},
        "total_time_s": WINDOW_S,
        "attribute_mechanism": False,
        "timelapse_frames": 0,
    }
    pk = overrides.pop("pk", None)
    spec.update(overrides)
    if pk:
        spec["pk"] = dict(spec["pk"], **pk)
    return spec


#: A dose the ``store`` fixture does NOT solve, so a scaling test cannot be
#: answered by the exact memo instead and quietly test nothing.
UNSOLVED_PEAK = 150.0


def _matching(mgr, spec, field="run_key"):
    """Every stored run the lookup would consider for ``spec``.

    Mutation tests have to contaminate ALL of them: the fixture deliberately
    solves one spec twice, so patching one copy leaves the other answering and
    the test passes while proving nothing.
    """
    want = run_keys(spec)[field]
    return [j for j in mgr.list()
            if j.status == "done" and isinstance(j.body, dict)
            and ((j.body.get("provenance") or {}).get("memo") or {}
                 ).get(field) == want]


def _numbers(body, prefix=""):
    """Every finite numeric leaf of a run body, by dotted path.

    The unit of comparison for "indistinguishable from a fresh solve".
    """
    out = {}
    if isinstance(body, dict):
        for k, v in body.items():
            if str(k).startswith("_"):
                continue
            out.update(_numbers(v, f"{prefix}.{k}" if prefix else str(k)))
    elif isinstance(body, (list, tuple)):
        for i, v in enumerate(body):
            out.update(_numbers(v, f"{prefix}[{i}]"))
    elif isinstance(body, bool) or body is None or isinstance(body, str):
        pass
    elif isinstance(body, (int, float)) and math.isfinite(float(body)):
        out[prefix] = float(body)
    return out


# --------------------------------------------------------------------------- #
# A solved store, shared by every class that needs one.
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def store(tmp_path_factory):
    """A JobManager holding three real runs plus the timings they cost.

    ``base`` and ``repeat`` are the SAME spec solved twice -- the second one is
    what a memo hit has to be indistinguishable from, and solving it is also
    the only way to know the engine is deterministic at all. ``dose`` is the
    same spec at 4x the plasma peak, which is what the algebra is checked
    against.
    """
    if SECTION is None:
        pytest.skip("no section is cached in data_cache/")
    root = str(tmp_path_factory.mktemp("memo_runs"))
    mgr = JobManager(root=root)
    specs = {
        "base": _spec(),
        "repeat": _spec(),
        "dose": _spec(pk={"plasma_peak": 200.0}),
    }
    runs, timings = {}, {}
    try:
        for name, spec in specs.items():
            t0 = time.time()
            job = mgr.submit(spec)
            deadline = time.time() + 900.0
            while job.status not in ("done", "failed", "cancelled"):
                if time.time() > deadline:
                    pytest.fail(f"run {name} did not finish")
                time.sleep(0.05)
            assert job.status == "done", (name, job.error)
            timings[name] = time.time() - t0
            runs[name] = job.run_id
        yield {"mgr": mgr, "specs": specs, "runs": runs, "timings": timings}
    finally:
        mgr.shutdown()


# --------------------------------------------------------------------------- #
# 0. Canonicalisation
#
# How a wrong number gets out if this fails: two specs that are DIFFERENT runs
# reduce to one key, and the second one is answered with the first one's
# numbers. That is the only failure mode in this file that produces a
# confidently wrong answer rather than a missed optimisation, so it is first.
# --------------------------------------------------------------------------- #
class TestCanonical:
    def test_key_order_does_not_change_the_key(self):
        a = {"section": "s", "dx_um": 12.0, "pk": {"plasma_peak": 1.0,
                                                   "half_life_s": 2.0}}
        b = {"pk": {"half_life_s": 2.0, "plasma_peak": 1.0},
             "dx_um": 12.0, "section": "s"}
        assert spec_digest(a) == spec_digest(b)

    def test_int_and_float_spelling_of_one_number_is_one_key(self):
        # ``execute_run`` coerces every numeric spec field with ``float()``, so
        # 12 and 12.0 ARE the same run and must not be two cache entries.
        assert spec_digest({"dx_um": 12}) == spec_digest({"dx_um": 12.0})

    def test_float_formatting_does_not_change_the_key(self):
        # Two spellings that parse to the same double.
        a = json.loads('{"total_time_s": 0.1}')
        b = json.loads('{"total_time_s": 0.10000000000000001}')
        assert a["total_time_s"] == b["total_time_s"]
        assert spec_digest(a) == spec_digest(b)

    @pytest.mark.parametrize("left,right,why", [
        ({"timelapse_frames": None}, {},
         "an explicit null and an absent field are DIFFERENT runs: "
         "execute_run reads spec.get('timelapse_frames', 24) and then "
         "int(x) if x is not None else 0, so null means 0 frames and "
         "absent means 24"),
        ({"resistance_off": True}, {"resistance_off": 1},
         "a bool and an int must not merge: nothing guarantees every "
         "consumer of a spec field coerces the same way"),
        ({"mgmt": None}, {"mgmt": 0.0},
         "mgmt=None is 'not supplied' and mgmt=0 is a supplied zero; "
         "execute_run branches on `mgmt is not None`"),
        ({"drug": "temozolomide"}, {"drug": "Temozolomide"},
         "compound names are resolved case-sensitively upstream"),
        ({"pk": {"plasma_peak": 50.0}}, {"pk": {}},
         "an omitted PK field resolves to a default; the memo does NOT "
         "resolve defaults, and treating the two as one key would make "
         "that resolution part of the cache"),
        ({"a": 1.0}, {"a": 1.0, "b": None},
         "an extra field, even a null one, is a different spec"),
        ({"x": [1.0, 2.0]}, {"x": [2.0, 1.0]},
         "list order is meaningful"),
        ({"x": "1"}, {"x": 1.0},
         "a string and a number must not merge"),
        ({"x": "1.0"}, {"x": 1.0},
         "and not even when the string is the number's own repr -- which is "
         "what an untagged canonicaliser produces"),
        ({"x": True}, {"x": "True"},
         "nor a bool and its text"),
    ])
    def test_distinct_specs_do_not_collide(self, left, right, why):
        assert spec_digest(left) != spec_digest(right), why

    def test_an_unrepresentable_number_is_refused_not_approximated(self):
        with pytest.raises(Unmemoisable):
            canonical({"timelapse_frames": 2 ** 53 + 1})
        with pytest.raises(Unmemoisable):
            canonical({"dx_um": float("nan")})
        with pytest.raises(Unmemoisable):
            canonical({"dx_um": float("inf")})

    def test_two_keys_with_the_same_text_are_refused(self):
        # Not reachable over HTTP (JSON object keys are strings) and refused
        # anyway, because canonicalising it would silently merge them.
        with pytest.raises(Unmemoisable):
            canonical({1: "a", "1": "b"})

    def test_an_unkeyable_spec_gets_no_memo_stamp(self):
        assert jobs_mod._memo_stamp({"dx_um": float("nan")}) is None


# --------------------------------------------------------------------------- #
# 1. What the key is made of
#
# How a wrong number gets out if this fails: the engine, or a solver constant
# the RunSpec does not carry, changes -- and the key does not, so yesterday's
# answer is served for today's engine with nothing marking it.
# --------------------------------------------------------------------------- #
@requires_section
class TestKeyComposition:
    @pytest.mark.parametrize("const,value", [
        ("PDE_DT_S", 5.0),
        ("PDE_SAVE_EVERY", 250),
        ("PDE_CONVECTION", False),
        ("SIM_SEED", 1),
        ("DEFAULT_TIMELAPSE_FRAMES", 12),
        ("TIMELAPSE_BUDGET_BYTES", 1024),
    ])
    def test_a_solver_constant_changes_the_key(self, monkeypatch, const, value):
        """A frozen RunSpec carries no solver settings, so they must be in the
        key or two runs computed with different numerics share one entry."""
        spec = _spec()
        before = run_keys(spec)
        monkeypatch.setattr(jobs_mod, const, value)
        after = run_keys(spec)
        assert after["run_key"] != before["run_key"], (
            f"{const} changed and the run key did not")
        assert after["shape_key"] != before["shape_key"], (
            f"{const} changed and the shape key did not")

    def test_the_engine_fingerprint_changes_the_key(self, monkeypatch):
        spec = _spec()
        before = run_keys(spec)
        stale = dict(memo_mod.engine_fingerprint())
        stale["digest"] = "0" * 64
        monkeypatch.setattr(memo_mod, "_FINGERPRINT", stale)
        after = run_keys(spec)
        assert after["run_key"] != before["run_key"]
        assert after["shape_key"] != before["shape_key"]

    def test_the_claim_set_version_changes_the_key(self, monkeypatch):
        """A withdrawal must make every memo a miss, so a repeat re-solves
        into caveats written under the CURRENT claim set rather than being
        handed older ones under a banner."""
        spec = _spec()
        before = run_keys(spec)
        monkeypatch.setattr(jobs_mod, "CLAIM_SET_VERSION",
                            jobs_mod.CLAIM_SET_VERSION + 1)
        assert run_keys(spec)["run_key"] != before["run_key"]

    def test_the_tissue_fingerprint_changes_the_key(self, monkeypatch):
        spec = _spec()
        before = run_keys(spec)
        monkeypatch.setattr(memo_mod, "_fingerprint",
                            lambda paths: "something-else")
        assert run_keys(spec)["run_key"] != before["run_key"]

    def test_an_unregistered_section_does_not_inherit_a_key(self):
        assert memo_mod.tissue_fingerprint("no-such-section") != \
            memo_mod.tissue_fingerprint(SECTION)

    def test_the_memo_format_version_is_in_the_key(self, monkeypatch):
        spec = _spec()
        before = run_keys(spec)
        monkeypatch.setattr(memo_mod, "MEMO_FORMAT_VERSION",
                            memo_mod.MEMO_FORMAT_VERSION + 1)
        assert run_keys(spec)["run_key"] != before["run_key"]

    def test_the_two_scaled_axes_are_the_only_shape_key_freedom(self):
        """``shape_key`` must collapse exactly the two measured axes and
        nothing else, or the algebra gets applied across a knob nobody
        measured."""
        base = run_keys(_spec())["shape_key"]
        for axis, value in (("plasma_peak", 123.0), ("half_life_s", 999.0)):
            assert run_keys(_spec(pk={axis: value}))["shape_key"] == base, axis
        for axis, value in (("n_doses", 2), ("dose_interval_s", 1.0)):
            assert run_keys(_spec(pk={axis: value}))["shape_key"] != base, axis
        for field, value in (("total_time_s", 21600.0), ("dx_um", 12.0),
                             ("drug", "doxorubicin"), ("state_mode", "markers"),
                             ("time_resolved", True), ("resistance_off", True)):
            assert run_keys(_spec(**{field: value}))["shape_key"] != base, field

    def test_every_engine_entrypoint_is_inside_the_hashed_set(self):
        """The fingerprint hashes ``tmesim/`` whole plus a NAMED list of server
        modules. A new server module in the run chain that nobody adds to that
        list is a hole in the key, so the list is pinned against the
        entrypoints a run records about itself."""
        declared = set(memo_mod._ENGINE_FILES)
        entrypoints = [
            "reconstruct.real_pipeline.reconstruct_real",
            "reconstruct.real_pipeline.apply_drug_sensitivity",
            "validation.calibration.calibrate_alpha",
            "validation.calibration.calibrated_drugparams",
            "molecular.binding_site_barrier.apply_barrier_to_drugparams",
            "sim.engine.simulate_response",
            "validation.vessel_profile.direct_vessel_distance_um",
            "validation.delivery_metrics.compute_all",
            "app.server.reach.compute_reach",
        ]
        for name in entrypoints:
            if name.startswith("app.server."):
                module = "app/server/" + name.split(".")[2] + ".py"
                assert module in declared, (
                    f"{name} reaches a run's numbers and {module} is not in "
                    "memo._ENGINE_FILES, so editing it would not invalidate a "
                    "single memo")
            else:
                # Everything else is under tmesim/, which is hashed whole.
                assert "tmesim" in memo_mod._ENGINE_TREES

    def test_the_memo_module_does_not_import_fastapi(self):
        """``jobs.py`` imports ``memo`` to stamp every run, and fastapi is in
        the OPTIONAL ``[app]`` extra. A module-level fastapi import here would
        make the engine depend on the web stack."""
        src = open(os.path.join(memo_mod._repo_root(),
                                "app", "server", "memo.py")).read()
        assert "fastapi" not in src, (
            "app/server/memo.py imports fastapi; the router belongs in "
            "app/server/memo_api.py")

    def test_the_declared_engine_files_exist(self):
        root = memo_mod._repo_root()
        missing = [p for p in memo_mod._ENGINE_FILES
                   if not os.path.isfile(os.path.join(root, p))]
        assert not missing, (
            f"memo._ENGINE_FILES names files that are gone: {missing}. The "
            "fingerprint records them as 'absent', which is honest, but a "
            "stale list means something else moved unwatched.")

    def test_no_host_path_reaches_the_returned_key_components(self):
        """The key is built from absolute host paths; the response is not
        allowed to carry them."""
        text = json.dumps(run_keys(_spec())["components"])
        for probe in (memo_mod._repo_root(), os.path.expanduser("~")):
            assert probe not in text, probe


# --------------------------------------------------------------------------- #
# 2. The central claim: a hit is a re-solve
# --------------------------------------------------------------------------- #
@requires_section
class TestExactMemo:
    def test_the_engine_is_deterministic_on_a_repeated_spec(self, store):
        """If this fails nothing else in this file means anything: the memo's
        promise is only interesting because a re-solve would agree."""
        mgr = store["mgr"]
        a = _numbers(mgr.get(store["runs"]["base"]).body)
        b = _numbers(mgr.get(store["runs"]["repeat"]).body)
        shared = set(a) & set(b)
        assert len(shared) > 100, "almost nothing was compared"
        drift = {k: (a[k], b[k]) for k in shared
                 if a[k] != b[k] and not k.endswith(("_time_s", "_s"))
                 and "time" not in k}
        assert not drift, f"the same spec solved twice disagreed: {drift}"

    def test_a_hit_is_bit_identical_to_a_fresh_solve(self, store):
        mgr = store["mgr"]
        hit = memo_lookup(mgr, store["specs"]["base"])
        assert hit["mode"] == "exact_memo", hit.get("reason")
        fresh = _numbers(mgr.get(store["runs"]["repeat"]).body)
        served = _numbers(hit["result"])
        compared = 0
        for key, value in fresh.items():
            if "time" in key or key.endswith("_s"):
                continue      # wall clocks are not results
            if key not in served:
                continue
            compared += 1
            assert served[key] == value, (
                f"{key}: memo served {served[key]!r}, a fresh solve of the "
                f"same spec produced {value!r}")
        assert compared > 100, f"only {compared} numbers were compared"

    def test_a_hit_declares_that_it_is_not_a_solve(self, store):
        answer = memo_mod.answer_spec(store["mgr"], store["specs"]["base"])
        assert answer["mode"] == "exact_memo"
        assert answer["is_solve"] is False
        assert "NOT A SOLVE" in answer["disclosure"]
        # ``base`` and ``repeat`` are the same spec solved twice, so either is
        # a correct answer; the lookup takes the newest.
        assert answer["source_run_id"] in (store["runs"]["base"],
                                           store["runs"]["repeat"])
        assert answer["source_run_id"] in answer["disclosure"], (
            "a hit must name the execution it is recalling, or a client "
            "cannot tell that its run_id belongs to an older run")

    def test_a_hit_is_faster_than_the_solve_it_replaces(self, store):
        mgr = store["mgr"]
        t0 = time.time()
        for _ in range(5):
            hit = memo_lookup(mgr, store["specs"]["base"])
        elapsed = (time.time() - t0) / 5
        assert hit["mode"] == "exact_memo"
        assert elapsed < store["timings"]["repeat"] / 10.0, (
            f"a memo hit took {elapsed:.3f}s against {store['timings']['repeat']:.3f}s "
            "for the solve it replaces")

    def test_a_different_spec_does_not_hit(self, store):
        for spec in (_spec(dx_um=12.0), _spec(drug="doxorubicin"),
                     _spec(total_time_s=21600.0), _spec(pk={"n_doses": 2})):
            out = memo_lookup(store["mgr"], spec)
            assert out["mode"] == "miss", spec

    def test_a_key_match_alone_is_not_enough(self, store):
        """The key is a sha256 and the stored run carries its own spec, so the
        spec half of the match is checked for EQUALITY too. Forced here by
        stamping a DIFFERENT spec's run with this spec's key: the key says yes
        and the specs say no, and no must win."""
        mgr, spec = store["mgr"], store["specs"]["base"]
        other = _spec(pk={"plasma_peak": 200.0})
        victim = _matching(mgr, other)
        assert victim, "the 4x-dose run is not in the store"
        job = victim[0]
        body = job.body
        prov = dict(body["provenance"])
        prov["memo"] = dict(prov["memo"], run_key=run_keys(spec)["run_key"])
        job.body = dict(body, provenance=prov)
        # Hide the genuine copies so only the forged one could answer.
        genuine = [(j, j.status) for j in _matching(mgr, spec)]
        for real, _ in genuine:
            real.status = "hidden-for-this-test"
        try:
            out = memo_lookup(mgr, spec)
            assert out["mode"] == "miss", (
                "a run whose stamped key matched but whose spec did not was "
                "served")
        finally:
            for real, status in genuine:
                real.status = status
            job.body = body

    def test_a_run_that_did_not_finish_is_never_served(self, store):
        """A failed or cancelled run has no answer to recall. The gate is
        ``status == 'done'``, asserted rather than assumed because a failed
        run still carries a spec and could still carry a stamp."""
        mgr, spec = store["mgr"], store["specs"]["base"]
        jobs_matching = _matching(mgr, spec)
        assert jobs_matching
        saved = [(j, j.status) for j in jobs_matching]
        for job in jobs_matching:
            job.status = "failed"
        try:
            assert memo_lookup(mgr, spec)["mode"] == "miss"
        finally:
            for job, status in saved:
                job.status = status

    def test_a_protocol_sweep_is_never_a_memo_candidate(self, store):
        """``execute_run`` dispatches a ``kind='protocol'`` spec to
        ``execute_protocol_job`` BEFORE the provenance block, so a TPL sweep
        never gets a stamp and can never be recalled. Fail-closed by
        construction rather than by a check -- pinned here so the dispatch
        cannot move above the stamp without this failing."""
        import inspect

        src = inspect.getsource(jobs_mod.execute_run)
        dispatch = src.index("execute_protocol_job(spec, ctx, cache=cache)")
        stamp = src.index("_memo_stamp(spec)")
        assert dispatch < stamp, (
            "the protocol dispatch no longer returns before the memo stamp, "
            "so a TPL sweep could now be stamped and recalled as if it were "
            "a RunSpec")

    def test_an_empty_store_misses(self, tmp_path):
        mgr = JobManager(root=str(tmp_path / "empty"), start_worker=False)
        try:
            out = memo_lookup(mgr, _spec())
            assert out["mode"] == "miss"
            assert "POST /api/runs" in out["reason"]
        finally:
            mgr.shutdown()


# --------------------------------------------------------------------------- #
# 3. Staleness
#
# How a wrong number gets out if this fails: the engine changes and the memo
# keeps answering. This is the whole reason the module exists.
# --------------------------------------------------------------------------- #
@requires_section
class TestStaleness:
    def test_a_changed_engine_turns_a_hit_into_a_miss(self, store, monkeypatch):
        mgr, spec = store["mgr"], store["specs"]["base"]
        assert memo_lookup(mgr, spec)["mode"] == "exact_memo"
        stale = dict(memo_mod.engine_fingerprint())
        stale["digest"] = "f" * 64
        monkeypatch.setattr(memo_mod, "_FINGERPRINT", stale)
        out = memo_lookup(mgr, spec)
        assert out["mode"] == "miss", (
            "the engine fingerprint moved and the memo still answered")

    def test_a_changed_alpha_turns_a_hit_into_a_miss(self, store, monkeypatch):
        """``calibrate_alpha()`` is in the fingerprint because a source hash
        cannot see a changed calibration INPUT."""
        mgr, spec = store["mgr"], store["specs"]["base"]
        assert memo_lookup(mgr, spec)["mode"] == "exact_memo"
        monkeypatch.setattr(memo_mod, "_FINGERPRINT", None)
        monkeypatch.setattr("tmesim.validation.calibration.calibrate_alpha",
                            lambda *a, **k: 19.0)
        try:
            assert memo_lookup(mgr, spec)["mode"] == "miss"
        finally:
            memo_mod._FINGERPRINT = None

    def test_a_changed_solver_constant_turns_a_hit_into_a_miss(
            self, store, monkeypatch):
        mgr, spec = store["mgr"], store["specs"]["base"]
        assert memo_lookup(mgr, spec)["mode"] == "exact_memo"
        monkeypatch.setattr(jobs_mod, "PDE_DT_S", 5.0)
        assert memo_lookup(mgr, spec)["mode"] == "miss"

    def test_a_run_with_no_memo_stamp_is_never_a_candidate(self, store):
        """Every run written before this field existed. Its engine state was
        not recorded, so it cannot be shown to match -- and a memo that
        assumed it did would be guessing."""
        mgr, spec = store["mgr"], store["specs"]["base"]
        jobs_matching = _matching(mgr, spec)
        assert jobs_matching, "nothing matched, so nothing was tested"
        saved = [(j, j.body) for j in jobs_matching]
        for job in jobs_matching:
            prov = dict(job.body["provenance"])
            prov.pop("memo")
            job.body = dict(job.body, provenance=prov)
        try:
            assert memo_lookup(mgr, spec)["mode"] == "miss"
        finally:
            for job, body in saved:
                job.body = body

    def test_a_run_stamped_by_an_older_key_format_is_never_a_candidate(
            self, store):
        mgr, spec = store["mgr"], store["specs"]["base"]
        jobs_matching = _matching(mgr, spec)
        assert jobs_matching, "nothing matched, so nothing was tested"
        saved = [(j, j.body) for j in jobs_matching]
        for job in jobs_matching:
            prov = dict(job.body["provenance"])
            prov["memo"] = dict(prov["memo"],
                                format=prov["memo"]["format"] - 1)
            job.body = dict(job.body, provenance=prov)
        try:
            assert memo_lookup(mgr, spec)["mode"] == "miss"
        finally:
            for job, body in saved:
                job.body = body


# --------------------------------------------------------------------------- #
# 4. The claim set
#
# How a wrong number gets out if this fails: a withdrawn claim is laundered.
# ``/api/runs`` banners a stale record; a memo that copied it without the
# banner would strip the one thing marking it.
# --------------------------------------------------------------------------- #
@requires_section
class TestClaimSet:
    def test_a_withdrawal_makes_the_memo_a_miss(self, store, monkeypatch):
        mgr, spec = store["mgr"], store["specs"]["base"]
        assert memo_lookup(mgr, spec)["mode"] == "exact_memo"
        monkeypatch.setattr(jobs_mod, "CLAIM_SET_VERSION",
                            jobs_mod.CLAIM_SET_VERSION + 1)
        monkeypatch.setattr(jobs_mod, "CLAIM_SET_LOG", jobs_mod.CLAIM_SET_LOG + (
            {"version": jobs_mod.CLAIM_SET_VERSION, "adopted": "2099-01-01",
             "withdrew": "a test withdrawal", "sources": (), "retired_text": ()},
        ))
        out = memo_lookup(mgr, spec)
        assert out["mode"] == "miss", (
            "a claim was withdrawn and the memo kept serving caveats written "
            "before it")

    def test_a_contradicted_record_is_refused_even_at_the_same_version(
            self, store, monkeypatch):
        """The case a version comparison cannot see: a record stamped with the
        CURRENT claim set whose caveats still carry a sentence that claim set
        retired. ``r_00090`` in ``app/runs/`` is the real specimen."""
        mgr, spec = store["mgr"], store["specs"]["base"]
        retired = jobs_mod.CLAIM_SET_LOG[-1].get("retired_text") or ()
        if not retired:
            pytest.skip("the current claim set retired no literal text")
        jobs_matching = _matching(mgr, spec)
        assert jobs_matching, "nothing matched, so nothing was tested"
        saved = [(j, j.body) for j in jobs_matching]
        for job in jobs_matching:
            prov = dict(job.body["provenance"])
            prov["caveats"] = list(prov["caveats"]) + [f"...{retired[0]}..."]
            job.body = dict(job.body, provenance=prov)
        try:
            out = memo_lookup(mgr, spec)
            assert out["mode"] == "miss"
            assert "CONTRADICT" in out["reason"].upper()
        finally:
            for job, body in saved:
                job.body = body

    def test_one_clean_copy_still_answers_when_another_is_contradicted(
            self, store):
        """The flip side, and the reason the gate is per-candidate: a
        contradicted copy must not take a clean one down with it."""
        mgr, spec = store["mgr"], store["specs"]["base"]
        retired = jobs_mod.CLAIM_SET_LOG[-1].get("retired_text") or ()
        jobs_matching = _matching(mgr, spec)
        if not retired or len(jobs_matching) < 2:
            pytest.skip("needs two stored copies and a retired sentence")
        job = jobs_matching[0]
        body = job.body
        prov = dict(body["provenance"])
        prov["caveats"] = list(prov["caveats"]) + [f"...{retired[0]}..."]
        job.body = dict(body, provenance=prov)
        try:
            out = memo_lookup(mgr, spec)
            assert out["mode"] == "exact_memo"
            assert out["source_run_id"] != job.run_id
        finally:
            job.body = body

    def test_a_hit_carries_the_live_claim_set_verdict(self, store):
        """The verdict is a fact about the READER'S moment, so it has to be
        computed on this read and not copied from the record."""
        hit = memo_lookup(store["mgr"], store["specs"]["base"])
        status = hit["result"]["provenance"]["claim_set"]
        assert status["current"] == jobs_mod.CLAIM_SET_VERSION
        assert status["superseded"] is False
        assert status["contradicted"] is False


# --------------------------------------------------------------------------- #
# 5. Payload honesty
#
# How a wrong number gets out if this fails: the memo hands back a run whose
# binary payload is gone while reporting it as present. The caller skips the
# solve, asks for the file, and gets a 410 for a run the memo said was
# complete. ``scripts/prune_run_payloads.py`` exists because ``rm`` makes the
# server explain confidently that the payload was never captured.
# --------------------------------------------------------------------------- #
@requires_section
class TestPayloadHonesty:
    def test_a_complete_run_reports_its_payloads_present(self, store):
        mgr, run_id = store["mgr"], store["runs"]["base"]
        status = payload_status(mgr, run_id, mgr.get(run_id).body)
        assert status["complete"] is True
        assert status["payloads"]["geometry"]["state"] == "present"
        # This fixture asks for no time-lapse, so there is nothing to lose.
        assert status["payloads"]["timelapse"]["state"] == "not_produced"

    def test_a_pruned_payload_is_not_present(self, store, tmp_path):
        """Exactly what ``scripts/prune_run_payloads.py`` leaves behind.

        Applied to EVERY copy of the spec, because one surviving copy would
        legitimately answer and the assertion would then prove nothing.
        """
        mgr, spec = store["mgr"], store["specs"]["base"]
        jobs_matching = _matching(mgr, spec)
        assert jobs_matching, "nothing matched, so nothing was tested"
        saved = []
        for job in jobs_matching:
            path = mgr.payload_path(job.run_id, "geometry")
            backup = tmp_path / f"{job.run_id}-geometry.bin"
            with open(path, "rb") as fh:
                backup.write_bytes(fh.read())
            saved.append((job, job.body, path, backup))
            os.remove(path)
            job.body = dict(job.body, pruned_payloads={
                "geometry": {"at": "2026-08-30",
                             "bytes": backup.stat().st_size}})
        try:
            job = jobs_matching[0]
            status = payload_status(mgr, job.run_id, job.body)
            rec = status["payloads"]["geometry"]
            assert rec["state"] == "pruned"
            assert rec["available"] is False
            assert status["complete"] is False
            assert "deleted on 2026-08-30" in rec["reason"]
            hit = memo_lookup(mgr, spec)
            assert hit["mode"] == "exact_memo"
            assert hit["payloads"]["complete"] is False, (
                "the memo answered with a run whose payload is gone and "
                "called the answer complete")
            assert hit["payloads"]["payloads"]["geometry"]["state"] == "pruned"
            # And the run document beside it must agree -- ``result_view``
            # reconciles a pruned block at read time, so a memo hit cannot
            # claim the file is there in one field and gone in another.
            assert "pruned_payloads" in hit["result"]
        finally:
            for job, body, path, backup in saved:
                with open(path, "wb") as fh:
                    fh.write(backup.read_bytes())
                job.body = body

    def test_a_complete_copy_is_preferred_over_a_pruned_one(self, store, tmp_path):
        """Two copies of one spec, one pruned: the memo must hand back the one
        that can still serve its payloads."""
        mgr, spec = store["mgr"], store["specs"]["base"]
        jobs_matching = _matching(mgr, spec)
        if len(jobs_matching) < 2:
            pytest.skip("needs two stored copies of one spec")
        victim = jobs_matching[0]
        path = mgr.payload_path(victim.run_id, "geometry")
        backup = tmp_path / "victim.bin"
        with open(path, "rb") as fh:
            backup.write_bytes(fh.read())
        body = victim.body
        os.remove(path)
        victim.body = dict(body, pruned_payloads={
            "geometry": {"at": "2026-08-30", "bytes": 1}})
        try:
            hit = memo_lookup(mgr, spec)
            assert hit["mode"] == "exact_memo"
            assert hit["source_run_id"] != victim.run_id
            assert hit["payloads"]["complete"] is True
        finally:
            with open(path, "wb") as fh:
                fh.write(backup.read_bytes())
            victim.body = body

    def test_a_bare_rm_is_reported_as_missing(self, store, tmp_path):
        """No prune record at all -- a bare ``rm``, a partial copy, a failed
        restore. The record says nothing, so only the filesystem can answer."""
        mgr, run_id = store["mgr"], store["runs"]["base"]
        path = mgr.payload_path(run_id, "geometry")
        backup = tmp_path / "geometry2.bin"
        with open(path, "rb") as fh:
            backup.write_bytes(fh.read())
        os.remove(path)
        try:
            rec = payload_status(mgr, run_id,
                                 mgr.get(run_id).body)["payloads"]["geometry"]
            assert rec["state"] == "missing", (
                "a payload that is simply gone was not reported gone")
            assert rec["available"] is False
            assert "no prune recorded" in rec["reason"]
        finally:
            with open(path, "wb") as fh:
                fh.write(backup.read_bytes())

    def test_a_scaled_answer_offers_no_payload_at_all(self, store):
        out = scale_lookup(store["mgr"], _spec(pk={"plasma_peak": UNSOLVED_PEAK}))
        assert out["mode"] == "algebraic_scale", out.get("reason")
        assert "payloads" not in out
        assert "source_run_id" in out, (
            "a scaled answer must name the run it was derived from")


# --------------------------------------------------------------------------- #
# 6. The algebra, checked against a real re-solve
# --------------------------------------------------------------------------- #
@requires_section
class TestAlgebra:
    def test_the_plasma_scalar_matches_what_the_run_recorded(self, store):
        """The scale factor is a ratio of two numbers this module computes.
        The denominator must be the number the stored run actually used, bit
        for bit, or the factor is a guess."""
        mgr = store["mgr"]
        for name, run_id in store["runs"].items():
            spec = mgr.get(run_id).spec
            computed = plasma_scalars(spec["pk"], spec["total_time_s"])
            recorded = mgr.get(run_id).body["extra_metrics"]["plasma_mean_uM"]
            assert computed["plasma_mean"] == recorded, name

    def test_a_dose_change_is_answered_without_solving(self, store):
        out = scale_lookup(store["mgr"], _spec(pk={"plasma_peak": 200.0}))
        assert out["mode"] == "algebraic_scale", out.get("reason")
        assert out["derivation"]["axes_changed"] == ["plasma_peak"]
        assert out["derivation"]["plasma_peak_ratio"] == 4.0
        assert out["derivation"]["bit_exact"] is True, (
            "a 4x dose change IS a power of two, so the scaling is exact in "
            "floating point as well as in algebra")

    def test_a_scaled_dose_matches_the_solve_it_replaces(self, store):
        """The one measurement that decides whether any of this is allowed."""
        mgr = store["mgr"]
        out = scale_lookup(mgr, store["specs"]["dose"])
        assert out["mode"] == "algebraic_scale", out.get("reason")
        solved = mgr.get(store["runs"]["dose"]).body
        checked = 0
        for block, leaves in out["readouts"].items():
            for leaf, rec in leaves.items():
                truth = (solved.get(block) or {}).get(leaf)
                if rec["value"] is None or truth is None:
                    continue
                checked += 1
                if rec.get("basis") == "invariant" and rec["max_measured_ulp"] == 0:
                    assert rec["value"] == float(truth), (
                        f"{block}.{leaf} is declared a bit-identical "
                        f"invariant and is not")
                    continue
                allowed = max(rec["max_measured_ulp"], 1.0) * math.ulp(
                    max(abs(rec["value"]), abs(float(truth))))
                assert abs(rec["value"] - float(truth)) <= allowed, (
                    f"{block}.{leaf}: scaled {rec['value']!r} against solved "
                    f"{truth!r}, outside the {rec['max_measured_ulp']} ulp "
                    "this rule declares")
        assert checked >= 15, f"only {checked} readouts were checked"

    def test_a_4x_dose_is_bit_identical_not_merely_close(self, store):
        """4 is a power of two, so there is no round-off to allow for. The
        surrogate study's "exactly 4.0" is true BECAUSE of this, and a 1.5x
        change does not get the same guarantee -- which is why the response
        reports ``bit_exact`` per request rather than claiming it once."""
        mgr = store["mgr"]
        out = scale_lookup(mgr, store["specs"]["dose"])
        solved = mgr.get(store["runs"]["dose"]).body
        for block, leaves in out["readouts"].items():
            for leaf, rec in leaves.items():
                truth = (solved.get(block) or {}).get(leaf)
                if rec["value"] is None or truth is None:
                    continue
                assert rec["value"] == float(truth), (
                    f"{block}.{leaf}: a power-of-two dose change should be "
                    f"bit-exact, got {rec['value']!r} vs {truth!r}")

    def test_the_reach_block_is_copied_because_it_has_no_pk_input(self, store):
        mgr = store["mgr"]
        out = scale_lookup(mgr, store["specs"]["dose"])
        solved = mgr.get(store["runs"]["dose"]).body["reach"]
        assert _numbers(out["reach"]) == _numbers(solved), (
            "reach.py solves the steady-state field from a UNIT source, so "
            "the whole block is PK-independent by construction")
        assert "UNIT source" in out["reach_basis"]

    def test_the_report_and_compound_blocks_are_copied(self, store):
        mgr = store["mgr"]
        out = scale_lookup(mgr, store["specs"]["dose"])
        solved = mgr.get(store["runs"]["dose"]).body
        assert _numbers(out["report"]) == _numbers(solved["report"])
        assert _numbers(out["drug_params"]) == _numbers(solved["drug_params"])

    def test_the_verdict_stability_census_travels_with_the_numbers(self, store):
        """A delivery number read without it could be one of the seven
        compounds in twelve whose delivery/dose verdict INVERTS under the
        transvascular term the transport law lacks, with nothing saying so."""
        mgr = store["mgr"]
        out = scale_lookup(mgr, store["specs"]["dose"])
        solved = mgr.get(store["runs"]["dose"]).body["extra_metrics"]
        key = "extra_metrics.verdict_stability"
        if solved.get("verdict_stability") is None:
            pytest.skip("the census did not score this compound")
        assert out["copied_verbatim"][key] == solved["verdict_stability"]
        assert key not in out["refused"]

    def test_a_scaled_answer_declares_that_it_is_not_a_solve(self, store):
        answer = memo_mod.answer_spec(store["mgr"],
                                      _spec(pk={"plasma_peak": UNSOLVED_PEAK}),
                                      allow_scaling=True)
        assert answer["mode"] == "algebraic_scale", answer.get("reason")
        assert answer["is_solve"] is False
        assert "NOT A SOLVE AND NOT A STORED RUN" in answer["disclosure"]
        assert "rearranged by an exact algebraic identity" in answer["disclosure"]

    def test_scaling_is_opt_in(self, store):
        answer = memo_mod.answer_spec(store["mgr"],
                                      _spec(pk={"plasma_peak": UNSOLVED_PEAK}))
        assert answer["mode"] == "miss", (
            "a caller who asked for an exact memo was handed an algebraic "
            "rearrangement without asking")
        assert "allow_scaling=true" in answer["reason"]

    def test_a_half_life_change_is_not_claimed_bit_exact(self, store):
        out = scale_lookup(store["mgr"], _spec(pk={"half_life_s": 57600.0}))
        assert out["mode"] == "algebraic_scale", out.get("reason")
        assert out["derivation"]["axes_changed"] == ["half_life_s"], (
            "the candidate differing along ONE axis must be preferred over "
            "the one differing along two")
        assert out["derivation"]["bit_exact"] is False, (
            "the plasma-average ratio for a half-life change is not a power "
            "of two, so the scaled magnitudes are not bit-identical and must "
            "not be advertised as such")
        assert "EXACT IN EXACT ARITHMETIC" in out["derivation"]["exactness"]
        assert out["derivation"]["plasma_peak_ratio"] == 1.0

    def test_a_disagreeing_stored_plasma_mean_refuses(self, store):
        """If the PK model, the time grid or the default resolution moves under
        a stored record, the factor would be wrong with nothing downstream able
        to see it. Refuse instead."""
        mgr = store["mgr"]
        query = _spec(pk={"plasma_peak": UNSOLVED_PEAK})
        candidates = _matching(mgr, query, field="shape_key")
        assert candidates, "nothing matched, so nothing was tested"
        saved = [(j, j.body) for j in candidates]
        for job in candidates:
            xm = job.body["extra_metrics"]
            job.body = dict(job.body, extra_metrics=dict(
                xm, plasma_mean_uM=float(xm["plasma_mean_uM"]) * 2))
        try:
            out = scale_lookup(mgr, query)
            assert out["mode"] == "miss"
            assert "does not match" in out["reason"]
        finally:
            for job, body in saved:
                job.body = body


# --------------------------------------------------------------------------- #
# 7. What the algebra refuses
#
# How a wrong number gets out if this fails: a readout nobody measured, or one
# measured NOT to scale, is multiplied by a factor anyway. This is the class
# that stops the identity being widened past its evidence.
# --------------------------------------------------------------------------- #
@requires_section
class TestScalingRefusals:
    @pytest.mark.parametrize("spec_kwargs,knob", [
        ({"pk": {"n_doses": 2}}, "n_doses"),
        ({"pk": {"dose_interval_s": 1800.0}}, "dose_interval_s"),
        ({"total_time_s": 21600.0}, "total_time_s"),
        ({"dx_um": 12.0}, "dx_um"),
        ({"drug": "doxorubicin"}, "drug"),
        ({"state_mode": "markers"}, "state_mode"),
        ({"resistance_off": True}, "resistance_off"),
    ])
    def test_n_doses_is_never_scaled(self, store, spec_kwargs, knob):
        """``n_doses`` and ``dose_interval_s`` reach the solver only through
        ``plasma_mean``, so the derivation says they scale too -- and every row
        of the measurement this module rests on carries ``n_doses = 1``. A
        derived-but-unmeasured identity is not an exact answer, so they are
        refused. ``total_time_s`` is not even derivable: it moves the
        integration window as well as the plasma average."""
        out = scale_lookup(store["mgr"], _spec(**spec_kwargs))
        assert out["mode"] == "miss", (
            f"{knob} changed and an algebraic answer was returned for it")

    def test_time_resolved_transport_is_refused(self, store):
        """On that path the plasma curve drives the solver step by step, so a
        half-life change reshapes the source in time instead of rescaling
        it."""
        out = scale_lookup(store["mgr"],
                           _spec(time_resolved=True, pk={"plasma_peak": 200.0}))
        assert out["mode"] == "miss"
        assert "time_resolved" in out["reason"]

    def test_kill_fraction_is_never_served(self, store):
        mgr = store["mgr"]
        out = scale_lookup(mgr, store["specs"]["dose"])
        assert out["mode"] == "algebraic_scale"
        for block, leaves in out["readouts"].items():
            assert "kill_fraction" not in leaves, (
                "an absolute kill fraction is never a defensible output of "
                "this product and it is not a Hill-linear quantity either")
        assert "metrics.kill_fraction" in out["refused"]
        assert "kill" not in json.dumps(out["readouts"])

    def test_no_rule_names_a_kill_or_coverage_readout(self):
        for block, rules in READOUT_RULES.items():
            for leaf in rules:
                lowered = leaf.lower()
                assert "kill" not in lowered, (block, leaf)
                assert "coverage" not in lowered, (block, leaf)
                assert "viable" not in lowered, (block, leaf)

    @pytest.mark.parametrize("leaf", [
        "kill_fraction", "therapeutic_coverage_1x", "therapeutic_coverage_0_5x",
        "therapeutic_coverage_2x", "viable_tumor_fraction",
        "penetration_length_um",
    ])
    def test_every_measured_non_scalable_readout_is_refused_with_a_reason(
            self, store, leaf):
        out = scale_lookup(store["mgr"], store["specs"]["dose"])
        hits = [k for k in out["refused"] if k.endswith("." + leaf)]
        assert hits, f"{leaf} was neither served nor refused"
        for key in hits:
            assert "REFUSED" in out["refused"][key]
            assert "Re-solve" in out["refused"][key]

    def test_an_unmeasured_readout_is_refused_rather_than_guessed(
            self, store, monkeypatch):
        """The fail-closed half. A metric added to ``execute_run`` is
        unmeasured until somebody measures it, and the scaler must refuse it
        instead of assuming it behaves like its neighbours."""
        mgr = store["mgr"]
        query = _spec(pk={"plasma_peak": UNSOLVED_PEAK})
        candidates = _matching(mgr, query, field="shape_key")
        assert candidates, "nothing matched, so nothing was tested"
        saved = [(j, j.body) for j in candidates]
        for job in candidates:
            job.body = dict(job.body, extra_metrics=dict(
                job.body["extra_metrics"], a_brand_new_metric=1.0))
        try:
            out = scale_lookup(mgr, query)
            assert out["mode"] == "algebraic_scale", out.get("reason")
            assert "a_brand_new_metric" not in out["readouts"]["extra_metrics"]
            assert "has not been measured" in \
                out["refused"]["extra_metrics.a_brand_new_metric"]
        finally:
            for job, body in saved:
                job.body = body

    def test_a_stored_null_is_not_turned_into_a_number(self, store):
        mgr = store["mgr"]
        query = _spec(pk={"plasma_peak": UNSOLVED_PEAK})
        candidates = _matching(mgr, query, field="shape_key")
        assert candidates, "nothing matched, so nothing was tested"
        saved = [(j, j.body) for j in candidates]
        for job in candidates:
            job.body = dict(job.body, extra_metrics=dict(
                job.body["extra_metrics"], E_p50_uM=None))
        try:
            out = scale_lookup(mgr, query)
            rec = out["readouts"]["extra_metrics"]["E_p50_uM"]
            assert rec["value"] is None
            assert rec["basis"] == "stored_null"
        finally:
            for job, body in saved:
                job.body = body

    def test_the_scaled_axes_list_is_exactly_the_measured_pair(self):
        assert SCALED_AXES == ("plasma_peak", "half_life_s"), (
            "this list is the entire scope of the algebraic answer. Adding to "
            "it widens a claim past the measurement in "
            "experiments/exact_memo/probe_scaling.py.")

    def test_every_rule_carries_a_measured_bound_and_a_reason(self):
        for block, rules in READOUT_RULES.items():
            for leaf, rule in rules.items():
                a, b, ulp, note = rule
                assert isinstance(a, int) and isinstance(b, int), (block, leaf)
                assert ulp >= 0.0, (block, leaf)
                assert len(note) > 15, (
                    f"{block}.{leaf} has no written reason for its exponents")


# --------------------------------------------------------------------------- #
# 8. The route
#
# Mounted here exactly as it is mounted in ``create_app`` -- one line -- so
# this suite exercises the real HTTP surface without ``api.py`` having to carry
# a second copy of the wiring.
# --------------------------------------------------------------------------- #
@requires_section
class TestRoute:
    @pytest.fixture(scope="class")
    def client(self, store):
        from fastapi.testclient import TestClient

        from app.server.api import create_app
        from app.server.memo_api import build_memo_router

        app = create_app(store["mgr"])
        # THE ONE LINE, and the PLACEMENT is load-bearing exactly as the other
        # routers' mounting block in ``api.py`` says: ``include_router``
        # appends, and ``_mount_web`` mounts StaticFiles at "/" which matches
        # every path, so a router added after it is registered and never
        # reached. ``create_app`` has already mounted the web root by the time
        # it returns, so this fixture inserts the route at the FRONT -- which
        # is also the evidence that appending would not work.
        router = build_memo_router(jobs=store["mgr"])
        app.router.routes[:0] = router.routes
        with TestClient(app) as c:
            yield c

    def test_appending_after_create_app_would_never_be_reached(self, store):
        """The failure the fixture above avoids, asserted rather than assumed:
        this is why the mount line goes inside ``create_app``, ahead of
        ``_mount_web``."""
        import os as _os

        from fastapi import FastAPI
        from fastapi.testclient import TestClient
        from starlette.staticfiles import StaticFiles

        from app.server.api import create_app
        from app.server.memo_api import build_memo_router

        # BUILT ON A BARE APP, NOT ON create_app, and the reason is worth
        # keeping. This half used to append the router to a real create_app()
        # and assert a non-200. That stopped being observable the moment
        # api.py mounted the router properly: create_app now serves /api/memo
        # itself, so a second appended copy is never reached and the
        # assertion saw the FIRST copy's 200.
        #
        # The premise had disappeared; the mechanism had not. This reproduces
        # the mechanism directly -- a StaticFiles catch-all at `/`, then an
        # appended router -- which is what the placement claim is about and
        # stays true however api.py is organised later. Third test in this
        # repository to need this correction, after the browse and
        # cross-section ones.
        web = _os.path.join(
            _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
            "app", "web")
        too_late = FastAPI()
        too_late.mount("/", StaticFiles(directory=web, html=True), name="web")
        too_late.include_router(build_memo_router(jobs=store["mgr"]))
        with TestClient(too_late) as c:
            r = c.post("/api/memo", json=store["specs"]["base"])
        assert r.status_code != 200, (
            "a router appended after a StaticFiles catch-all at / was "
            "reached; the placement note in api.py's mounting block would "
            "then be wrong")

        # And the real app must serve it, which is the other half of the claim.
        real = create_app(store["mgr"])
        with TestClient(real) as c:
            assert c.post("/api/memo",
                          json=store["specs"]["base"]).status_code == 200, (
                "create_app does not serve /api/memo -- the include_router "
                "line is missing or sits after _mount_web")

    def test_a_bare_runspec_body_hits(self, client, store):
        r = client.post("/api/memo", json=store["specs"]["base"])
        assert r.status_code == 200, r.text
        body = r.json()
        assert body["mode"] == "exact_memo"
        assert body["is_solve"] is False
        assert body["scaling_requested"] is False

    def test_a_wrapped_body_can_ask_for_scaling(self, client):
        r = client.post("/api/memo", json={
            "spec": _spec(pk={"plasma_peak": UNSOLVED_PEAK}),
            "allow_scaling": True})
        assert r.status_code == 200, r.text
        body = r.json()
        assert body["mode"] == "algebraic_scale", body.get("reason")
        assert body["is_solve"] is False
        assert "NOT A SOLVE" in body["disclosure"]

    def test_an_invalid_spec_is_answered_not_crashed(self, client):
        r = client.post("/api/memo", json={"section": "no-such-section",
                                           "dx_um": 12.0})
        assert r.status_code == 200, r.text
        body = r.json()
        assert body["mode"] == "miss"
        assert body["validate"]["ok"] is False

    def test_a_non_object_body_is_refused_in_words(self, client):
        r = client.post("/api/memo", json={"spec": 7})
        assert r.status_code == 200, r.text
        assert r.json()["mode"] == "miss"

    def test_the_response_carries_no_host_path(self, client, store):
        r = client.post("/api/memo", json=store["specs"]["base"])
        for probe in (memo_mod._repo_root(), os.path.expanduser("~")):
            assert probe not in r.text, probe

    def test_the_route_is_under_api_so_the_gate_covers_it(self, client):
        paths = {getattr(r, "path", None) for r in client.app.router.routes}
        assert "/api/memo" in paths
