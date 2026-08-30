"""Tests for the Vivum Studio backend (``app/server/``).

Four things are checked, in rising order of what they would cost if they broke:

1. **Binary payload round-trip** — the exact byte layout API_CONTRACT.md freezes
   (magic, little-endian header fields, typed-array-ready bodies). The front end
   slices these into typed arrays; a silent layout change would render garbage.
2. **Job lifecycle** — submit -> queued -> running -> done, live status with a
   log tail, cancel, and history that survives a restart.
3. **Reconstruction cache** — the same
   ``(section, dx, vascular_mode, oxygen_mode, state_mode)`` is reconstructed once and
   reused, and different keys do not collide.
4. **CLI parity** — a run through the API produces numbers IDENTICAL (exact
   equality, not tolerance) to the same run driven straight through the
   production entry points the way ``experiments/patient_report.py`` drives
   them. This is the test that stops the API drifting into a second, subtly
   different pipeline.

The end-to-end tests use the smallest real section (``xenium_319``, 40,887
measured cells) at a coarse ``dx`` and a short window so the whole file runs in
well under a minute.
"""
from __future__ import annotations

import json
import os
import re
import struct
import sys
import time

import numpy as np
import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

fastapi = pytest.importorskip("fastapi", reason="the app/server backend needs fastapi")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi.testclient import TestClient  # noqa: E402

from app.server import cache as cache_mod  # noqa: E402
from app.server.cache import repo_root  # noqa: E402
from app.server.api import create_app, default_spec, validate_spec  # noqa: E402
from app.server.jobs import PDE_DT_S, PDE_SAVE_EVERY, SIM_SEED, JobManager  # noqa: E402
from app.server.payloads import (  # noqa: E402
    FIELD_HEADER_STRUCT,
    GEOMETRY_HEADER_STRUCT,
    MAGIC_FIELD,
    MAGIC_GEOMETRY,
    PAYLOAD_VERSION,
    decode_field,
    decode_geometry,
    encode_field,
    encode_geometry,
)

# Small, fast, real: the 40,887-cell brain section on a coarse grid.
#: Every environment variable ``app.server.__main__.main()`` WRITES into
#: ``os.environ``. The three helpers that call ``main()`` in-process must
#: isolate all of them, and they share this tuple SPECIFICALLY so a variable
#: cannot be added to main() and then remembered at only one of the three
#: sites -- which is how this has now failed twice:
#:
#:   TMESIM_STUDIO_ALLOWED_HOSTS leaked into every app built by every later
#:   test module and 400'd them all.
#:   TMESIM_STUDIO_DATA_ROOTS (added to main() 2026-08-22, for the non-loopback
#:   default) leaked into tests/test_ingest.py and produced eight 403
#:   path_outside_data_roots failures, because a scan confined to
#:   <repo>/data_cache refuses pytest's tmp_path. It passed in isolation and
#:   failed only in a full run, which reads exactly like a flaky test and is
#:   not one.
#:
#: If you add a variable to main(), add it here.
MAIN_WRITES_ENV = (
    "TMESIM_STUDIO_API_KEY",
    # Added 2026-08-30 with /api/account, which reports WHICH of main()'s three
    # key sources won. Omitting it from this tuple leaked a stale "flag" or
    # "generated" into every later app in the session -- the same failure the
    # note above describes, one variable later.
    "TMESIM_STUDIO_API_KEY_SOURCE",
    "TMESIM_STUDIO_ALLOWED_HOSTS",
    "TMESIM_STUDIO_DATA_ROOTS",
)

SECTION = "xenium_319"
DX_UM = 30.0
WINDOW_S = 600.0
DRUG = "temozolomide"


def _spec(**overrides):
    spec = {
        "section": SECTION,
        "dx_um": DX_UM,
        "drug": DRUG,
        "sensitivity_mode": "uniform",
        "vascular_mode": "discrete",
        "oxygen_mode": "pde",
        "mgmt": None,
        "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
               "n_doses": 1, "dose_interval_s": 86400.0},
        "total_time_s": WINDOW_S,
    }
    spec.update(overrides)
    return spec


def _section_available(section_id: str = SECTION) -> bool:
    return any(s["id"] == section_id and s["available"]
               for s in cache_mod.list_sections())


requires_section = pytest.mark.skipif(
    not _section_available(),
    reason=f"section {SECTION!r} is not cached in data_cache/",
)


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    """A TestClient whose run history is isolated from the real ``app/runs/``."""
    root = str(tmp_path_factory.mktemp("runs"))
    manager = JobManager(root=root)
    app = create_app(manager)
    with TestClient(app) as c:
        c._manager = manager  # type: ignore[attr-defined]
        c._root = root  # type: ignore[attr-defined]
        yield c
    manager.shutdown()


def _await(client, run_id, timeout=300.0):
    """Poll ``/status`` the way the UI does until the run reaches a terminal
    state, then return the final status body."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = client.get(f"/api/runs/{run_id}/status").json()
        if body["status"] in ("done", "failed", "cancelled"):
            return body
        time.sleep(0.2)
    raise AssertionError(f"run {run_id} did not finish within {timeout}s")


def _run(client, spec, timeout=300.0):
    posted = client.post("/api/runs", json=spec)
    assert posted.status_code == 202, posted.text
    run_id = posted.json()["run_id"]
    status = _await(client, run_id, timeout)
    if status["status"] != "done":
        tb = client.get(f"/api/runs/{run_id}/traceback").json()
        raise AssertionError(
            f"run {run_id} ended {status['status']}: {tb.get('error')}\n"
            f"{tb.get('traceback')}"
        )
    return client.get(f"/api/runs/{run_id}").json()


# --------------------------------------------------------------------------- #
# 1. Binary payloads
# --------------------------------------------------------------------------- #
class TestPayloadRoundTrip:
    def test_geometry_round_trip(self):
        rng = np.random.default_rng(0)
        n = 5000
        x = rng.uniform(0, 8120.0, n)
        y = rng.uniform(0, 5100.0, n)
        types = rng.integers(0, 5, n)

        blob = encode_geometry(x, y, types, (8120.0, 5100.0))
        # magic | version | n | ex | ey | x f32[n] | y f32[n] | type u8[n]
        assert len(blob) == GEOMETRY_HEADER_STRUCT.size + n * 4 + n * 4 + n

        out = decode_geometry(blob)
        assert out["magic"] == MAGIC_GEOMETRY
        assert out["version"] == PAYLOAD_VERSION
        assert out["n_cells"] == n
        assert out["extent_um"] == (8120.0, 5100.0)
        # float32 is the wire type, so compare at float32 precision.
        assert np.array_equal(out["x"], x.astype("<f4"))
        assert np.array_equal(out["y"], y.astype("<f4"))
        assert np.array_equal(out["cell_type"], types.astype(np.uint8))

    def test_geometry_header_is_little_endian_and_exact(self):
        blob = encode_geometry([1.5], [2.5], [3], (10.0, 20.0))
        assert blob[:4] == b"TMSG"
        assert struct.unpack_from("<I", blob, 4)[0] == PAYLOAD_VERSION
        assert struct.unpack_from("<I", blob, 8)[0] == 1
        assert struct.unpack_from("<f", blob, 12)[0] == 10.0
        assert struct.unpack_from("<f", blob, 16)[0] == 20.0
        assert struct.unpack_from("<f", blob, 20)[0] == 1.5   # x[0]
        assert struct.unpack_from("<f", blob, 24)[0] == 2.5   # y[0]
        assert blob[28] == 3                                   # type[0]
        assert len(blob) == 29
        # A big-endian reader would see a different (wrong) count.
        assert struct.unpack_from(">I", blob, 8)[0] != 1

    def test_field_round_trip_and_row_major(self):
        rows, cols = 37, 53
        values = np.arange(rows * cols, dtype=np.float64).reshape(rows, cols)
        blob = encode_field(values, dx_um=12.0)
        assert len(blob) == FIELD_HEADER_STRUCT.size + rows * cols * 4

        out = decode_field(blob)
        assert out["magic"] == MAGIC_FIELD
        assert (out["rows"], out["cols"]) == (rows, cols)
        assert out["dx_um"] == 12.0
        assert out["min"] == 0.0
        assert out["max"] == float(rows * cols - 1)
        assert np.array_equal(out["values"], values.astype("<f4"))
        # Row-major: element (r, c) sits at r*cols + c in the flat body.
        flat = np.frombuffer(blob, dtype="<f4",
                             offset=FIELD_HEADER_STRUCT.size)
        assert flat[5 * cols + 7] == values[5, 7]

    def test_field_range_ignores_non_finite_and_never_invents_one(self):
        values = np.array([[1.0, np.nan], [np.inf, 3.0]])
        out = decode_field(encode_field(values, dx_um=1.0))
        assert out["min"] == 1.0 and out["max"] == 3.0

        allnan = np.full((2, 2), np.nan)
        out2 = decode_field(encode_field(allnan, dx_um=1.0))
        assert np.isnan(out2["min"]) and np.isnan(out2["max"])

    def test_malformed_payloads_raise(self):
        blob = encode_field(np.zeros((2, 2)), dx_um=1.0)
        with pytest.raises(ValueError, match="bad magic"):
            decode_geometry(blob)
        with pytest.raises(ValueError, match="expected"):
            decode_field(blob[:-4])
        with pytest.raises(ValueError, match="uint8"):
            encode_geometry([0.0], [0.0], [999], (1.0, 1.0))


# --------------------------------------------------------------------------- #
# 2. Static endpoints + validation
# --------------------------------------------------------------------------- #
class TestStaticEndpoints:
    def test_sections_enumerate_without_crashing(self, client):
        rows = client.get("/api/sections").json()
        ids = {r["id"] for r in rows}
        assert {"xenium_319", "xenium_io", "gse324268"} <= ids
        assert len([r for r in rows if r["id"].startswith("visium")]) == 4
        for r in rows:
            assert set(r) >= {"id", "label", "n_cells", "extent_um", "panel",
                              "available"}
            if not r["available"]:
                assert r["reason"], f"{r['id']} unavailable with no reason"
                assert r["n_cells"] is None and r["extent_um"] is None
            if r["n_cells"] is None:
                assert r.get("n_cells_reason")

    def test_drugs_carry_measured_ic50_and_its_source(self, client):
        rows = client.get("/api/drugs").json()
        by = {r["name"]: r for r in rows}
        assert {"temozolomide", "doxorubicin", "bevacizumab"} <= set(by)
        for r in rows:
            assert set(r) >= {"name", "mw", "logp", "ic50_uM", "ic50_source",
                              "bbb_penetrance"}
            assert r["ic50_source"], f"{r['name']} has no IC50 provenance"

        tmz = by["temozolomide"]
        # The real measured value, cross-checked against potency.resolve_ic50.
        from tmesim.molecular.potency import resolve_ic50

        truth = resolve_ic50("temozolomide")
        assert tmz["ic50_uM"] == pytest.approx(truth["ic50_uM"], rel=1e-12)
        assert truth["confidence"] == "measured"
        assert tmz["ic50_source"] == f"GDSC measured (n={truth['n_cell_lines']})"

    def test_defaults_are_a_valid_runspec(self, client):
        spec = client.get("/api/defaults").json()
        assert set(spec) == {"section", "dx_um", "drug", "sensitivity_mode",
                             "vascular_mode", "oxygen_mode", "state_mode",
                             "mgmt", "pk", "total_time_s"}
        verdict = client.post("/api/validate", json=spec).json()
        assert verdict["ok"], verdict["errors"]

    def test_engine_mode_tuples_match_the_registry(self):
        from tmesim.physics.oxygen_proxy import OXYGEN_MODES as ENGINE_OXY
        from tmesim.reconstruct.real_fields import VASCULAR_MODES as ENGINE_VASC
        from tmesim.reconstruct.cell_state import STATE_MODES as ENGINE_STATE
        from tmesim.reconstruct.real_pipeline import (
            STATE_MODE_DEFAULT,
            SENSITIVITY_MODES as ENGINE_SENS,
        )

        assert set(cache_mod.OXYGEN_MODES) == set(ENGINE_OXY)
        assert set(cache_mod.VASCULAR_MODES) == set(ENGINE_VASC)
        assert set(cache_mod.SENSITIVITY_MODES) == set(ENGINE_SENS)
        assert set(cache_mod.STATE_MODES) == set(ENGINE_STATE)
        assert cache_mod._STATE_MODE_DEFAULT == STATE_MODE_DEFAULT


class TestIngestPage:
    """``GET /ingest`` -- the server-rendered bring-your-own-tissue form.

    ``results/first_user_readiness.md`` found the ingest path fully working
    through the API and the CLI but exposed NOWHERE a user could click, so a
    researcher who would not hand-POST JSON could not load their own sample.
    This page is that surface, served from the backend so it needs no front-end
    build. These tests pin the two things that make it useful: it is reachable,
    and it drives the real endpoint rather than a parallel implementation.
    """

    def test_page_is_served(self, client):
        r = client.get("/ingest")
        assert r.status_code == 200
        assert "text/html" in r.headers["content-type"]

    def test_page_drives_the_real_scan_endpoint(self, client):
        body = client.get("/ingest").text
        # The whole point: one fetch, to the endpoint curl would have used.
        assert "/api/sections/scan" in body
        assert body.count("fetch(") == 1
        # Carries its own assets -- no external host, and no dependency on
        # app/web/ having been built.
        assert "http://" not in body and "https://" not in body
        assert "<script" in body and "src=" not in body

    def test_page_offers_the_api_key_a_remote_deployment_needs(self, client):
        # Non-loopback deployments always have a key (see __main__), and /api/*
        # is gated while this page is not -- so it must be able to send one.
        assert "X-API-Key" in client.get("/ingest").text

    def test_scan_answers_a_bad_directory_instead_of_500ing(self, client, tmp_path):
        empty = tmp_path / "not_a_sample"
        empty.mkdir()
        r = client.post("/api/sections/scan", json={"path": str(empty)})
        assert r.status_code == 200
        assert r.json()["error"], "an unusable directory must come back with a reason"

    def test_scan_requires_a_path(self, client):
        assert client.post("/api/sections/scan", json={}).status_code == 400


class TestValidate:
    def test_unknown_section_and_drug_are_actionable(self):
        v = validate_spec(_spec(section="not_a_section", drug="aspirin"))
        assert not v["ok"]
        joined = " ".join(v["errors"])
        assert "not_a_section" in joined and "Known ids" in joined
        assert "aspirin" in joined and "temozolomide" in joined

    def test_bad_modes_are_rejected(self):
        v = validate_spec(_spec(sensitivity_mode="magic",
                                vascular_mode="swirl", oxygen_mode="vibes"))
        assert not v["ok"]
        assert len(v["errors"]) >= 3

    @requires_section
    def test_markers_oxygen_rejected_on_a_panel_without_hypoxia_genes(self):
        # The 319-gene brain panel carries none of HYPOXIA_MARKER_GENES, so
        # hypoxia is UNMEASURED there and the engine would raise mid-run.
        v = validate_spec(_spec(oxygen_mode="markers"))
        assert not v["ok"]
        assert any("markers" in e and "hypoxia" in e for e in v["errors"])

    def test_unknown_state_mode_is_rejected(self):
        v = validate_spec(_spec(state_mode="vibes"))
        assert not v["ok"]
        assert any("state_mode" in e and "vibes" in e for e in v["errors"])

    def test_default_state_mode_is_accepted_and_silent_about_the_route(self):
        v = validate_spec(_spec(state_mode="oxygen"))
        assert v["ok"], v["errors"]
        # The default must not warn that it "changes" anything: it is
        # bit-identical to the pre-option pipeline.
        assert not any("changes cells.state" in w for w in v["warnings"])

    def test_marker_state_rejected_on_visium(self):
        v = validate_spec(_spec(section="visium_gbm1", state_mode="markers",
                                dx_um=100.0))
        assert not v["ok"]
        assert any("state_mode='markers'" in e and "Xenium only" in e
                   for e in v["errors"])

    @requires_section
    def test_marker_state_warns_loudly_on_the_319_panel(self):
        """results/cell_state_impact.md: the 319-gene flagship panel carries NO
        hypoxia gene, so HYPOXIC is UNMEASURED there, and proliferation detection
        is 13x below the published Ki-67 range so the per-cell call is refused.
        The route must not silently produce a number."""
        v = validate_spec(_spec(state_mode="markers"))
        joined = " ".join(v["warnings"])
        assert "UNMEASURED" in joined
        assert "hypoxia-response transcript" in joined
        # It must say what the fallback is and that its LEVEL is not an index.
        assert "NEIGHBOURHOOD" in joined or "neighbourhood" in joined
        assert "1.42%" in joined and "19-40%" in joined
        # And it must never be silent: some warning has to name state_mode.
        assert any("state_mode" in w for w in v["warnings"])

    @requires_section
    def test_absolute_kill_is_warned_against_on_that_panel_under_every_route(self):
        """§5: on a panel with no hypoxia gene NO state route clears the bar, so
        an absolute kill must not be quoted from it under ANY route."""
        for mode in ("oxygen", "vessel_distance", "markers"):
            v = validate_spec(_spec(state_mode=mode))
            assert any("ABSOLUTE KILL" in w for w in v["warnings"]), mode
            hit = next(w for w in v["warnings"] if "ABSOLUTE KILL" in w)
            # CONVERTED 2026-08-21, not weakened. This line was
            # `assert "RANKING" in hit  # what IS defensible`. Compound ordering
            # was measured and withdrawn (results/pk_ranking_robustness.md,
            # results/capacity_ranking_robustness.md), so a warning that still
            # named it would be pointing a reader at a retired output. The
            # CONTRACT is unchanged and is the reason the assertion existed: a
            # warning that says what is NOT defensible must also say what IS,
            # or it leaves the reader with no usable output at all. Only the
            # answer moved -- from the cross-drug ranking to the within-map
            # relative reading -- so the assertion moved with it and got
            # stricter rather than looser.
            assert "MAP" in hit, mode                    # what IS defensible
            assert "RELATIVE" in hit, mode               # ... and in what sense
            assert "another compound" in hit, mode       # ... and what is not
            # The withdrawal is CITED at the point of use, not just recorded in
            # a results file, so a reader meeting this warning can check it.
            assert "pk_ranking_robustness.md" in hit, mode
            assert "capacity_ranking_robustness.md" in hit, mode
            # Unchanged from the original assertion: the measured size of the
            # swing this warning exists to report.
            assert "0.3252" in hit           # the measured size of the swing

    def test_absolute_kill_warning_is_not_raised_without_a_drug(self):
        v = validate_spec(_spec(drug=None))
        assert not any("ABSOLUTE KILL" in w for w in v["warnings"])

    def test_markers_oxygen_rejected_on_visium(self):
        v = validate_spec(_spec(section="visium_gbm1", oxygen_mode="markers",
                                dx_um=100.0))
        assert not v["ok"]
        assert any("Xenium only" in e for e in v["errors"])

    def test_out_of_range_pk_is_rejected(self):
        v = validate_spec(_spec(pk={"plasma_peak": 0.0, "half_life_s": -1.0,
                                    "n_doses": 0, "dose_interval_s": 0.0}))
        assert not v["ok"]
        joined = " ".join(v["errors"])
        assert "plasma_peak" in joined
        assert "half_life_s" in joined
        assert "n_doses" in joined

    def test_total_time_below_the_solver_step_is_rejected(self):
        v = validate_spec(_spec(total_time_s=PDE_DT_S / 2.0))
        assert not v["ok"]
        assert any("dt" in e for e in v["errors"])

    def test_non_numeric_dx_is_rejected_not_coerced(self):
        v = validate_spec(_spec(dx_um="twelve"))
        assert not v["ok"]
        assert any("dx_um" in e for e in v["errors"])

    @requires_section
    def test_cfl_is_flagged_once_a_speed_has_been_measured(self, client):
        """After a real run the peak Darcy speed is known for that key, so a dx
        below ``dt * max|v|`` must be a hard error naming the numbers.

        RE-PINNED 2026-08-14 — ``results/area_conservation_promotion.md``.

        This test used to derive both probe meshes from THIS key's own speed. It
        therefore passed only when no *other* cached mesh for the section carried
        more than twice this key's speed — because when the probe dx is not
        itself a cached key, ``validate_spec`` falls back to tier 2 and uses the
        **maximum** speed recorded for the section at ANY mesh
        (``_cfl_check``, and ``solve_ifp`` reads the vascular field, not the
        mesh). ``app/runs/_cfl_speeds.json`` is git-ignored machine-local state
        that accumulates across sessions and already held a ``dx = 12 um`` key at
        0.6336 um/s against this ``dx = 30 um`` key's 0.2762, so the margin was
        1.09x and the test was one small change away from red on any developer
        machine with a finer mesh cached.

        The vessel-source promotion (``geometry='binary'`` -> ``'area'``) lowered
        the peak Darcy speed — the same source total spread over 3x the support
        gives a smoother pressure field — and tipped it over. **The promotion
        exposed the fragility; it did not create it.**

        The fix derives the probe meshes from the speed the CFL rule will
        ACTUALLY use, so the test pins the rule rather than the contents of a
        developer's cache directory. **No assertion was removed or weakened** —
        the error case, its message, and the safe case are all still asserted.
        """
        _run(client, _spec())  # measures and records the speed
        speed = cache_mod.RECON_CACHE.known_speed(
            cache_mod.make_key(SECTION, DX_UM, "discrete", "pde"))
        assert speed and speed > 0

        # The speed tier 2 will quote for a dx that is not itself a cached key.
        recorded = cache_mod.RECON_CACHE.speeds_for_section(SECTION, "discrete")
        ref_speed = max(recorded.values())
        assert ref_speed >= speed, (
            "this key's own speed must be among those recorded for the section")

        # dx small enough that dt*max|v|/dx > 1 for the speed actually quoted.
        bad_dx = 0.5 * ref_speed * PDE_DT_S
        v = validate_spec(_spec(dx_um=bad_dx))
        assert not v["ok"]
        assert any("CFL" in e and "dt*max|v|/dx" in e for e in v["errors"])
        # And the safe dx is honoured.
        assert validate_spec(_spec(dx_um=2.0 * ref_speed * PDE_DT_S))["ok"]

    def test_post_runs_rejects_an_invalid_spec_without_enqueueing(self, client):
        before = len(client.get("/api/runs").json())
        r = client.post("/api/runs", json=_spec(section="nope"))
        assert r.status_code == 400
        assert "errors" in r.json()["detail"]
        assert len(client.get("/api/runs").json()) == before


# --------------------------------------------------------------------------- #
# 3. Job lifecycle, payload endpoints, cache
# --------------------------------------------------------------------------- #
@requires_section
class TestJobLifecycle:
    def test_run_reports_stage_progress_and_a_log_tail(self, client):
        posted = client.post("/api/runs", json=_spec())
        assert posted.status_code == 202
        run_id = posted.json()["run_id"]
        assert run_id.startswith("r_")

        seen_stages = set()
        deadline = time.time() + 300
        while time.time() < deadline:
            st = client.get(f"/api/runs/{run_id}/status").json()
            assert set(st) >= {"status", "progress", "stage", "elapsed_s",
                               "log_tail"}
            assert 0.0 <= st["progress"] <= 1.0
            assert isinstance(st["log_tail"], list) and st["log_tail"]
            seen_stages.add(st["stage"])
            if st["status"] in ("done", "failed", "cancelled"):
                break
            time.sleep(0.05)

        assert st["status"] == "done", st
        assert st["progress"] == 1.0
        assert st["elapsed_s"] > 0
        assert "done" in seen_stages

    def test_result_carries_mandatory_provenance(self, client):
        res = _run(client, _spec())
        assert res["status"] == "done"
        assert res["runtime_s"] > 0
        assert res["report"]["n_cells"] > 0

        m = res["metrics"]
        assert set(m) >= {"kill_fraction", "mean_tumor_auc",
                          "therapeutic_coverage_1x", "penetration_length_um"}
        for name, value in m.items():
            if value is None:
                assert res["metric_reasons"].get(name), (
                    f"{name} is null with no reason")

        p = res["provenance"]
        # RE-PINNED at the two-route promotion (results/two_route_promotion.md).
        # SUPERSEDED: 148.6 +/- 0.5. alpha is fitted in closed form from
        # L = sqrt(D_eff/k_eff) on doxorubicin, whose compiled D fell 6.616x when
        # the two-route cell term became the shipping default, so alpha followed.
        # It is a GAUGE COORDINATE (results/alpha_identifiability.md), so this
        # pin exists to catch silent drift, not to assert a physical constant.
        assert p["alpha"] == pytest.approx(22.2235, abs=0.5)
        assert p["ic50_uM"] > 0
        assert p["ic50_source"].startswith("GDSC measured")
        assert p["sensitivity_mode"] == "uniform"
        assert p["sensitivity_upper_bound"] is True
        assert p["sensitivity_note"]
        assert len(p["caveats"]) >= 4
        assert any("UPPER BOUND" in c for c in p["caveats"])
        assert any("NOT FOR CLINICAL USE" in c for c in p["caveats"])

    def test_result_carries_the_current_claim_set(self, client):
        """A run executed NOW is stamped, and is therefore not bannered.

        The other half of tests/test_claim_set.py, which proves that a record
        WITHOUT the stamp (all 65 runs already in app/runs/, r_00052 among
        them) is bannered. This is the half that proves the banner is keyed to
        the claim set and not to "is this run old" -- a run that is seconds old
        and a run from last month get the same verdict if they carry the same
        stamp, and opposite verdicts if they do not.
        """
        from app.server.jobs import (CLAIM_SET_LOG, CLAIM_SET_VERSION,
                                     claim_set_stamp)

        res = _run(client, _spec())
        cs = res["provenance"]["claim_set"]
        assert cs["version"] == CLAIM_SET_VERSION
        assert cs["recorded"] is True
        assert cs["superseded"] is False, (
            "a run executed by this build is written against this build's "
            "claims; bannering it would be false")
        assert cs["warning"] is None and cs["changes"] == []

        # The document ON DISK carries the stamp only. A stored verdict would
        # be true when written and false after the next withdrawal.
        with open(os.path.join(client._root, res["run_id"],
                               "result.json")) as fh:
            stored = json.load(fh)
        assert stored["provenance"]["claim_set"] == claim_set_stamp()
        assert "superseded" not in stored["provenance"]["claim_set"]
        assert stored["provenance"]["claim_set"]["adopted"] == \
            CLAIM_SET_LOG[-1]["adopted"]

        # A stamp is only worth having if it is true. This section is
        # hypoxia-blind, so it takes the "ABSOLUTE KILL MUST NOT BE QUOTED"
        # caveat -- which, until the claim set landed, still closed by naming
        # the CROSS-DRUG RANKING as "the defensible output". That claim was
        # withdrawn on 2026-08-21, so a run stamped claim set 1 must not carry
        # it: the stamp would be certifying a withdrawn claim as current.
        cav = res["provenance"]["caveats"]
        assert any("MUST NOT BE QUOTED" in c for c in cav), (
            "the caveat this guards is missing; re-point the assertion below")
        assert not any(
            "RANKING is invariant across all three routes and is the "
            "defensible output" in c for c in cav)
        assert any("WHAT SURVIVES IS THE MAP" in c for c in cav)

    def test_reconstruction_only_run_nulls_drug_metrics_with_reasons(self, client):
        res = _run(client, _spec(drug=None))
        assert res["status"] == "done"
        assert res["report"]["n_cells"] > 0
        for name, value in res["metrics"].items():
            assert value is None
            assert "drug" in res["metric_reasons"][name]
        assert res["provenance"]["ic50_uM"] is None
        assert res["provenance"]["ic50_reason"]

        # ...and the drug_auc field says why it does not exist.
        run_id = res["run_id"]
        r = client.get(f"/api/runs/{run_id}/field/drug_auc")
        assert r.status_code == 404
        assert "drug=null" in r.json()["detail"]["reason"]
        # The geometry and the three reconstruction fields still exist.
        assert client.get(f"/api/runs/{run_id}/geometry").status_code == 200
        for name in ("vascular_density", "ecm_density", "oxygen"):
            assert client.get(
                f"/api/runs/{run_id}/field/{name}").status_code == 200

    def test_binary_endpoints_match_the_reconstruction(self, client):
        res = _run(client, _spec())
        run_id = res["run_id"]

        geo = client.get(f"/api/runs/{run_id}/geometry")
        assert geo.status_code == 200
        assert geo.headers["content-type"] == "application/octet-stream"
        g = decode_geometry(geo.content)
        assert g["n_cells"] == res["report"]["n_cells"]
        assert g["extent_um"] == pytest.approx(
            tuple(res["report"]["extent_um"]), rel=1e-5)

        rows, cols = res["report"]["grid"]
        for name in ("vascular_density", "ecm_density", "oxygen", "drug_auc"):
            r = client.get(f"/api/runs/{run_id}/field/{name}")
            assert r.status_code == 200, (name, r.text)
            f = decode_field(r.content)
            assert (f["rows"], f["cols"]) == (rows, cols)
            assert f["dx_um"] == pytest.approx(res["report"]["dx_um"])
            assert f["values"].shape == (rows, cols)

        assert client.get(f"/api/runs/{run_id}/field/nope").status_code == 404
        assert client.get("/api/runs/r_99999/geometry").status_code == 404

    def test_cancel_a_queued_run(self, client):
        # Two at once: the worker is single-threaded, so the second is queued
        # and can be cancelled before it ever starts.
        first = client.post("/api/runs", json=_spec()).json()["run_id"]
        second = client.post(
            "/api/runs", json=_spec(drug="lomustine")).json()["run_id"]
        r = client.delete(f"/api/runs/{second}")
        assert r.status_code == 200 and r.json()["cancelled"] is True

        body = _await(client, second)
        assert body["status"] == "cancelled"
        assert client.get(f"/api/runs/{second}").json()["error"]
        _await(client, first)

        # Cancelling a finished run is a no-op that says so.
        again = client.delete(f"/api/runs/{second}").json()
        assert again["cancelled"] is False and "already" in again["reason"]

    def _failed_job(self, client, monkeypatch, exc, suffix):
        import app.server.jobs as jobs_mod

        def boom(spec, ctx, cache=None):
            raise exc

        manager = JobManager(root=str(client._root) + suffix)
        monkeypatch.setattr(jobs_mod, "execute_run", boom)
        job = manager.submit(_spec())
        deadline = time.time() + 30
        while job.status not in ("failed", "done", "cancelled"):
            if time.time() > deadline:
                raise AssertionError("job never finished")
            time.sleep(0.05)
        manager.shutdown()
        return job

    def test_failed_runs_surface_their_error(self, client, monkeypatch):
        """An UNEXPECTED failure: redacted in the served field, kept in full
        where only an operator can reach it.

        This test used to assert the engine's raw text in ``job.error``. That
        assertion encoded a leak rather than a contract: ``job.error`` is
        served in full by ``GET /api/runs/{id}``, and a real run can fail with
        an exception whose message was built from an absolute host path --
        ``visium._locate`` resolves to ``<repo>/data_cache/visium_gbm/...`` and
        ``visium.py`` raises ``FileNotFoundError(f"no matrix.mtx(.gz) in
        {matrix_dir}")``. The assertion is not dropped, it MOVES: the detail
        must still exist, on the channel that is not the response body.
        """
        job = self._failed_job(
            client, monkeypatch, RuntimeError("synthetic engine failure"),
            "_fail")

        assert job.status == "failed"
        # The served field says what happened and how to have it looked up...
        assert "synthetic engine failure" not in job.error
        assert "RuntimeError" not in job.error
        ref = re.search(r"reference ([0-9a-f]{8})", job.error)
        assert ref, f"a failed run must be quotable; got {job.error!r}"
        # ...and the detail survives behind the opt-in traceback endpoint.
        assert job.traceback
        assert "synthetic engine failure" in job.traceback
        assert "RuntimeError" in job.traceback

    def test_a_deliberate_refusal_reaches_the_client_word_for_word(
            self, client, monkeypatch):
        """The other half of the same branch, and the one that breaks quietly.

        ``IngestError`` is the product's informative "no" -- it names every
        layout the ingester accepts. Both the CLI and the Studio render it
        verbatim, so redaction must not touch it. If this ever starts
        returning a reference number instead, the Studio has begun disagreeing
        with ``vivum ingest``.
        """
        from tmesim.data.ingest import IngestError

        refusal = ("no recognised layout. Expected one of: a Xenium bundle "
                   "(cells.csv.gz + cell_feature_matrix[.tar.gz] + "
                   "gene_panel.json), a Space Ranger output, or a MEX triple.")
        job = self._failed_job(client, monkeypatch, IngestError(refusal),
                               "_refused")

        assert job.status == "failed"
        assert job.error == refusal
        assert "unexpected internal failure" not in job.error
        assert "reference" not in job.error
        # The traceback is still captured for whoever wants it.
        assert job.traceback

    def test_the_refusal_names_are_the_real_classes(self):
        """``jobs.RUN_REFUSALS`` matches by NAME up the MRO, so that an error
        path never runs a heavy import (``CompoundRefused`` lives behind
        RDKit). A name that no longer corresponds to a real class would fail
        open -- redacting a deliberate refusal -- so pin the names."""
        from app.server.compounds import CompoundError
        from app.server.jobs import RUN_REFUSALS, is_run_refusal
        from app.server.projects import ProjectError
        from tmesim.data.ingest import IngestError
        from tmesim.protocol.errors import TplError

        for klass in (IngestError, TplError, CompoundError, ProjectError):
            assert klass.__name__ in RUN_REFUSALS
            assert is_run_refusal(klass("x"))
        assert "CompoundRefused" in RUN_REFUSALS
        assert not is_run_refusal(RuntimeError("x"))
        assert not is_run_refusal(FileNotFoundError("/private/path"))

    def test_history_survives_a_restart(self, client):
        res = _run(client, _spec())
        run_id = res["run_id"]

        # A brand-new manager over the same directory = a server restart.
        restarted = JobManager(root=client._root, start_worker=False)
        try:
            job = restarted.get(run_id)
            assert job is not None
            assert job.status == "done"
            assert job.result_view()["metrics"] == res["metrics"]
            assert os.path.exists(restarted.payload_path(run_id, "geometry"))
            assert os.path.exists(
                restarted.payload_path(run_id, "field_drug_auc"))
            # Run ids keep counting up rather than colliding.
            new_id = restarted.submit(_spec()).run_id
            assert new_id > run_id
        finally:
            restarted.shutdown()


@requires_section
class TestReconstructionCache:
    def test_same_key_is_reconstructed_once(self, client):
        cache_mod.RECON_CACHE.clear()
        before = cache_mod.RECON_CACHE.stats()

        first = _run(client, _spec())
        mid = cache_mod.RECON_CACHE.stats()
        assert mid["misses"] == before["misses"] + 1
        assert first["report"]["reconstruction_cache_hit"] is False

        # Same geometry key, different drug -> must reuse the reconstruction.
        second = _run(client, _spec(drug="lomustine"))
        after = cache_mod.RECON_CACHE.stats()
        assert after["misses"] == mid["misses"], "reconstruction rebuilt"
        assert after["hits"] == mid["hits"] + 1
        assert second["report"]["reconstruction_cache_hit"] is True
        # Same tissue, so the geometry report is identical...
        assert second["report"]["n_cells"] == first["report"]["n_cells"]
        assert second["report"]["grid"] == first["report"]["grid"]
        # ...but the drug physics is not.
        assert second["metrics"]["kill_fraction"] != \
            first["metrics"]["kill_fraction"]

    def test_a_different_dx_is_a_different_key(self, client):
        before = cache_mod.RECON_CACHE.stats()
        _run(client, _spec(dx_um=40.0))
        after = cache_mod.RECON_CACHE.stats()
        assert after["misses"] == before["misses"] + 1

    def test_key_is_the_five_contract_fields(self):
        k = cache_mod.make_key("s", 12.0, "Discrete", "PDE")
        assert k == ("s", 12.0, "discrete", "pde", "oxygen")
        assert cache_mod.make_key("s", 12.0, "discrete", "pde") == \
            cache_mod.make_key("s", 12.000000000001, "discrete", "pde")
        assert cache_mod.make_key("s", 12.0, "kde", "pde") != k

    def test_state_mode_is_part_of_the_key(self):
        """``state_mode`` rewrites ``cells.state``, so a reconstruction built
        under one route must never be served for another. Reusing it would
        report the wrong band occupancies — the largest single term in an
        absolute kill (results/cell_state_impact.md section 3.2)."""
        base = cache_mod.make_key("s", 12.0, "discrete", "pde", "oxygen")
        assert cache_mod.make_key("s", 12.0, "discrete", "pde", "markers") != base
        assert cache_mod.make_key(
            "s", 12.0, "discrete", "pde", "vessel_distance") != base
        assert cache_mod.make_key("s", 12.0, "discrete", "pde", "OXYGEN") == base

    def test_a_run_never_inherits_the_previous_run_sensitivity(self, client):
        """The cached microenvironment is shared, and
        ``apply_drug_sensitivity`` writes to it. Every run must start from the
        pristine reconstruction, so a determinant-mode run cannot contaminate a
        later uniform-mode one."""
        _run(client, _spec(sensitivity_mode="determinant", drug="doxorubicin"))
        after = _run(client, _spec(sensitivity_mode="uniform", drug=DRUG))
        sens = after["report"]["sensitivity"]
        assert sens["mode"] == "uniform"
        assert sens["mean_all_cells"] == pytest.approx(1.0)
        assert sens["is_uniform"] is True


# --------------------------------------------------------------------------- #
# 4. CLI parity — the anti-drift test
# --------------------------------------------------------------------------- #
@requires_section
def test_api_run_is_numerically_identical_to_the_cli(client):
    """Drive the production entry points directly, exactly as
    ``experiments/patient_report.py`` does, and require EXACT equality with the
    same run submitted through the API.

    Exact, not approximate: the engine is deterministic given ``cfg.seed``, so
    any difference at all means the API has forked the pipeline.
    """
    from tmesim.core import PDEConfig, PKConfig, SimConfig
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.library import get_molecule
    from tmesim.reconstruct.real_pipeline import (
        apply_drug_sensitivity,
        reconstruct_real,
    )
    from tmesim.sim.engine import simulate_response
    from tmesim.validation import delivery_metrics
    from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

    spec = _spec()

    # ---- the CLI path, written out in full ------------------------------- #
    micro = reconstruct_real(source="xenium", sample=None, dx_um=DX_UM,
                             oxygen_mode="pde")
    apply_drug_sensitivity(micro, DRUG, mode="uniform")
    alpha = calibrate_alpha()
    dp = apply_barrier_to_drugparams(
        calibrated_drugparams(DRUG, alpha), get_molecule(DRUG))
    cfg = SimConfig(
        drug=dp,
        pk=PKConfig(plasma_peak=spec["pk"]["plasma_peak"],
                    half_life=spec["pk"]["half_life_s"],
                    n_doses=spec["pk"]["n_doses"],
                    dose_interval=spec["pk"]["dose_interval_s"]),
        pde=PDEConfig(dt=PDE_DT_S, total_time=WINDOW_S, convection=True,
                      save_every=PDE_SAVE_EVERY),
        seed=SIM_SEED,
    )
    cli = simulate_response(micro, cfg)
    cli_metrics = delivery_metrics.compute_all(micro, cli, dp)
    cli_head = cli_metrics["headline"]

    # ---- the API path ------------------------------------------------------ #
    res = _run(client, spec)
    api = res["metrics"]

    assert api["kill_fraction"] == cli.kill_fraction
    assert api["mean_tumor_auc"] == cli_head["mean_tumor_auc_uM_s"]
    assert api["therapeutic_coverage_1x"] == \
        cli_head["therapeutic_coverage_1x_ic50"]
    assert api["penetration_length_um"] == cli_head["penetration_length_L_um"]

    # The provenance must describe that same compiled drug, not a lookalike.
    assert res["provenance"]["alpha"] == alpha
    assert res["provenance"]["ic50_uM"] == dp.ic50
    assert res["drug_params"]["bbb_penetrance"] == dp.bbb_penetrance
    assert res["drug_params"]["binding_1_per_s"] == dp.binding
    assert res["drug_params"]["diffusivity_um2_s"] == dp.diffusivity

    # ...and so must the reconstruction report and the extra metrics.
    assert res["report"]["n_cells"] == micro.cells.n
    assert res["extra_metrics"]["penetration_index"] == \
        cli.extras["penetration_index"]
    assert res["extra_metrics"]["plasma_mean_uM"] == cli.extras["plasma_mean"]


class TestApiKeyAuth:
    """The optional TMESIM_STUDIO_API_KEY gate on /api/* (app/server/api.py).

    Off by default (unset env var), so every test above this class runs
    against a completely unauthenticated app, exactly as before this gate was
    added -- these tests build their OWN app with the env var explicitly set,
    isolated via monkeypatch, rather than touching the module-scoped `client`
    fixture other tests share.
    """

    def _make_client(self, monkeypatch, key, root):
        if key is None:
            monkeypatch.delenv("TMESIM_STUDIO_API_KEY", raising=False)
        else:
            monkeypatch.setenv("TMESIM_STUDIO_API_KEY", key)
        manager = JobManager(root=str(root))
        app = create_app(manager)
        return TestClient(app), manager

    def test_a_platform_health_check_gets_200_and_nothing_else(
            self, monkeypatch, tmp_path):
        """Liveness answers an anonymous probe; diagnostics do not.

        WHY THIS EXISTS. A platform health checker cannot send an API key.
        Render was configured to probe ``/api/health``, got 401 every 10 s, and
        killed a container that was serving correctly on a port it had already
        detected -- the log read "401 Unauthorized" a hundred times and then
        "Timed Out". Adding ``/healthz`` did not fix it, because the path to
        probe lives in a dashboard, not in this repo.

        So both liveness routes answer an anonymous GET, and ``/api/health``
        WITHHOLDS ITS PAYLOAD instead: the runs directory, run counts, projects
        root, cache statistics and solver timestep are operational detail about
        a machine and still need the key. Exempting a route from the gate and
        exempting its DATA are different things, and only the first is done
        here.
        """
        client, manager = self._make_client(monkeypatch, "secret-123", tmp_path)
        try:
            with client:
                # Liveness, unauthenticated: 200, and literally nothing else.
                for path in ("/healthz", "/api/health"):
                    r = client.get(path)
                    assert r.status_code == 200, path
                    assert r.json() == {"ok": True}, (
                        f"{path} leaked detail to an anonymous caller: "
                        f"{r.json()}")

                # With the key, /api/health is the full diagnostic again.
                r = client.get("/api/health", headers={"X-API-Key": "secret-123"})
                assert r.status_code == 200
                body = r.json()
                for k in ("runs_dir", "runs", "reconstruction_cache"):
                    assert k in body, f"authenticated /api/health lost {k}"

                # And nothing else opened up. A real data route and a
                # state-changing one both still refuse.
                assert client.get("/api/sections").status_code == 401
                assert client.post(
                    "/api/sections/scan", json={"path": "/etc"}
                ).status_code == 401
                # POST to a liveness path is not a read, so it is not exempt.
                assert client.post("/healthz").status_code == 401
        finally:
            manager.shutdown()

    def test_no_key_configured_leaves_api_open(self, monkeypatch, tmp_path):
        client, manager = self._make_client(monkeypatch, None, tmp_path)
        try:
            with client:
                assert client.get("/api/sections").status_code == 200
        finally:
            manager.shutdown()

    def test_key_configured_rejects_missing_header(self, monkeypatch, tmp_path):
        client, manager = self._make_client(monkeypatch, "secret-123", tmp_path)
        try:
            with client:
                res = client.get("/api/sections")
                assert res.status_code == 401
                assert "X-API-Key" in res.json()["detail"]
        finally:
            manager.shutdown()

    def test_key_configured_rejects_wrong_header(self, monkeypatch, tmp_path):
        client, manager = self._make_client(monkeypatch, "secret-123", tmp_path)
        try:
            with client:
                res = client.get("/api/sections", headers={"x-api-key": "wrong"})
                assert res.status_code == 401
        finally:
            manager.shutdown()

    def test_key_configured_accepts_correct_header(self, monkeypatch, tmp_path):
        client, manager = self._make_client(monkeypatch, "secret-123", tmp_path)
        try:
            with client:
                res = client.get("/api/sections", headers={"x-api-key": "secret-123"})
                assert res.status_code == 200
        finally:
            manager.shutdown()

    def test_key_configured_does_not_gate_static_files(self, monkeypatch, tmp_path):
        """The frontend itself must still load with no header, so it can
        render a prompt for the key rather than failing to appear at all."""
        client, manager = self._make_client(monkeypatch, "secret-123", tmp_path)
        try:
            with client:
                res = client.get("/")
                assert res.status_code == 200
        finally:
            manager.shutdown()

    # -- GHSA-86qp-5c8j-p5mr: Host-header path poisoning ------------------- #
    #
    # starlette builds request.url by concatenating the client-supplied Host
    # header with the path and re-parsing it, so a "/" inside the Host value
    # terminates the authority and the remainder lands in url.path:
    #
    #     Host: x/../   ->   request.url.path == "/..//api/compounds"
    #
    # The gate used to test `request.url.path.startswith("/api/")`, which is
    # then False, while the ROUTER still dispatches on scope["path"], which is
    # unchanged. Gate skipped, endpoint runs, compound library returned with no
    # key. Verified live against starlette 0.49.3 before the fix.
    @pytest.mark.parametrize("host", [
        "x/../",                 # bare poisoned authority
        "testserver/../",        # ditto, using the host TestClient really uses
        "testserver:80/../",     # survives TrustedHostMiddleware's ":" split
        "x/..%2f",
        "x/y/",
    ])
    def test_a_poisoned_host_header_does_not_bypass_the_key(
            self, monkeypatch, tmp_path, host):
        client, manager = self._make_client(monkeypatch, "secret-123", tmp_path)
        try:
            with client:
                res = client.get("/api/sections", headers={"Host": host})
                assert res.status_code == 401, (
                    f"Host: {host!r} reached the endpoint without a key "
                    f"(status {res.status_code})")
        finally:
            manager.shutdown()

    def test_a_poisoned_host_header_does_not_bypass_the_key_on_a_mutation(
            self, monkeypatch, tmp_path):
        """DELETE matters more than GET: the bypass reached the destructive
        routes too. A 401 here, never a 404/200 from the endpoint itself."""
        client, manager = self._make_client(monkeypatch, "secret-123", tmp_path)
        try:
            with client:
                res = client.delete("/api/compounds/nope",
                                    headers={"Host": "x/../"})
                assert res.status_code == 401
        finally:
            manager.shutdown()

    def test_the_gate_uses_the_same_path_the_router_dispatches_on(
            self, monkeypatch, tmp_path):
        """The invariant behind the fix. A gate and a router that derive the
        path differently can always be made to disagree; deriving it the same
        way is what makes that impossible, so assert it structurally."""
        import inspect

        from app.server import api as api_mod

        src = inspect.getsource(api_mod.create_app)
        assert "get_route_path(request.scope)" in src
        gate = src.split("_require_api_key", 1)[1].split("jobs =", 1)[0]
        # Comments in the gate explain WHY url.path is banned, so strip them
        # before looking for actual uses of it.
        code = "\n".join(line for line in gate.splitlines()
                         if not line.lstrip().startswith("#"))
        assert "request.url.path" not in code, (
            "the auth gate must never make a decision on request.url.path -- "
            "the Host header can poison it (GHSA-86qp-5c8j-p5mr)")

    def test_a_non_get_route_outside_api_still_needs_the_key(
            self, monkeypatch, tmp_path):
        """Fail closed. Static assets are read-only, so anything that is not a
        plain GET/HEAD is authenticated regardless of its path -- a future
        state-changing route added outside /api/ is covered by default."""
        client, manager = self._make_client(monkeypatch, "secret-123", tmp_path)
        try:
            with client:
                assert client.post("/some/new/route").status_code == 401
        finally:
            manager.shutdown()


class TestCorsIsNotWildcard:
    """C2: `allow_origins=["*"]` meant any page the operator browsed could
    read AND delete their compound library through 127.0.0.1 -- a browser is a
    local process, so binding to loopback is no defence at all."""

    def _make_client(self, monkeypatch, tmp_path, origins=None):
        monkeypatch.delenv("TMESIM_STUDIO_API_KEY", raising=False)
        if origins is None:
            monkeypatch.delenv("TMESIM_STUDIO_ALLOWED_ORIGINS", raising=False)
        else:
            monkeypatch.setenv("TMESIM_STUDIO_ALLOWED_ORIGINS", origins)
        manager = JobManager(root=str(tmp_path))
        return TestClient(create_app(manager)), manager

    def test_no_cross_origin_grant_by_default(self, monkeypatch, tmp_path):
        client, manager = self._make_client(monkeypatch, tmp_path)
        try:
            with client:
                res = client.get("/api/sections",
                                 headers={"Origin": "https://evil.example"})
                assert "access-control-allow-origin" not in res.headers
        finally:
            manager.shutdown()

    def test_preflight_does_not_approve_delete_from_any_origin(
            self, monkeypatch, tmp_path):
        client, manager = self._make_client(monkeypatch, tmp_path)
        try:
            with client:
                res = client.options(
                    "/api/compounds/anything",
                    headers={"Origin": "https://evil.example",
                             "Access-Control-Request-Method": "DELETE"})
                assert res.headers.get("access-control-allow-origin") is None
        finally:
            manager.shutdown()

    def test_an_allowlisted_origin_is_granted_exactly_itself(
            self, monkeypatch, tmp_path):
        client, manager = self._make_client(
            monkeypatch, tmp_path, origins="http://localhost:5173")
        try:
            with client:
                res = client.get("/api/sections",
                                 headers={"Origin": "http://localhost:5173"})
                assert res.headers["access-control-allow-origin"] == \
                    "http://localhost:5173"
                other = client.get("/api/sections",
                                   headers={"Origin": "https://evil.example"})
                assert other.headers.get("access-control-allow-origin") is None
        finally:
            manager.shutdown()

    def test_a_wildcard_in_the_env_var_is_ignored_not_honoured(
            self, monkeypatch, tmp_path):
        """An operator reaching for the old behaviour must not get it back by
        writing `*` in the env var."""
        client, manager = self._make_client(monkeypatch, tmp_path, origins="*")
        try:
            with client:
                res = client.get("/api/sections",
                                 headers={"Origin": "https://evil.example"})
                assert res.headers.get("access-control-allow-origin") is None
        finally:
            manager.shutdown()


class TestTracebacksAreNotServedByDefault:
    """H4: traceback.format_exc() always carries absolute host paths -- the
    repo root and the operator's OS username."""

    def _failed_job(self, manager):
        job = manager.submit({"section": "x", "dx_um": 12.0})
        job.status = "failed"
        job.error = "ValueError: boom"
        job.traceback = (f'Traceback (most recent call last):\n'
                         f'  File "{repo_root()}/tmesim/x.py", line 1\n'
                         f'ValueError: boom\n')
        return job

    def test_traceback_is_withheld_and_the_typed_error_is_not(
            self, monkeypatch, tmp_path):
        monkeypatch.delenv("TMESIM_STUDIO_API_KEY", raising=False)
        monkeypatch.delenv("TMESIM_STUDIO_EXPOSE_TRACEBACKS", raising=False)
        manager = JobManager(root=str(tmp_path), start_worker=False)
        try:
            job = self._failed_job(manager)
            with TestClient(create_app(manager)) as client:
                res = client.get(f"/api/runs/{job.run_id}/traceback")
                assert res.status_code == 200
                body = res.json()
                assert body["traceback"] is None
                assert body["error"] == "ValueError: boom"
                assert "TMESIM_STUDIO_EXPOSE_TRACEBACKS" in \
                    body["traceback_withheld"]
                assert repo_root() not in res.text
        finally:
            manager.shutdown()

    def test_opt_in_serves_the_stack_with_the_repo_root_scrubbed(
            self, monkeypatch, tmp_path):
        monkeypatch.delenv("TMESIM_STUDIO_API_KEY", raising=False)
        monkeypatch.setenv("TMESIM_STUDIO_EXPOSE_TRACEBACKS", "1")
        manager = JobManager(root=str(tmp_path), start_worker=False)
        try:
            job = self._failed_job(manager)
            with TestClient(create_app(manager)) as client:
                res = client.get(f"/api/runs/{job.run_id}/traceback")
                assert "ValueError: boom" in res.json()["traceback"]
                assert "<repo>/tmesim/x.py" in res.json()["traceback"]
                assert repo_root() not in res.text
        finally:
            manager.shutdown()


class TestGridSizeIsBounded:
    """H5: grid memory scales as 1/dx^2 and dx_um was bounded only by
    `must be > 0`. Everything else was a WARNING, and warnings do not block --
    `ok = not errors`, and POST /api/runs only refuses when ok is False. One
    JSON POST could ask for ~20 GB."""

    @pytest.mark.parametrize("dx", [0.001, 0.01, 0.1, 0.49])
    def test_a_dx_below_the_floor_is_an_error_not_a_warning(self, dx):
        spec = dict(default_spec())
        spec["dx_um"] = dx
        verdict = validate_spec(spec)
        assert verdict["ok"] is False
        assert any("dx_um" in e for e in verdict["errors"])

    def test_a_grid_above_the_ceiling_is_refused_by_post_runs(self, client):
        spec = dict(default_spec())
        spec["dx_um"] = 0.5          # legal dx, but ~1.6e8 voxels on this section
        res = client.post("/api/runs", json=spec)
        assert res.status_code == 400
        assert "ceiling" in json.dumps(res.json())

    def test_the_default_spec_is_still_accepted(self, client):
        res = client.post("/api/validate", json=default_spec())
        assert res.json()["ok"] is True


class TestUnexpectedExceptionsDoNotLeakInternals:
    """Three routes answered an unexpected exception with
    ``f"{type(exc).__name__}: {exc}"`` straight into the response body. For
    the ingest routes that wrapped the caller's own path; for
    ``list_sections`` -- which every client hits on load -- it wrapped the
    whole section registry, so a missing ``data_cache/`` download surfaced a
    SERVER-side absolute path the user never supplied."""

    def test_enumeration_failure_gives_a_reference_not_a_path(self,
                                                              monkeypatch):
        boom = FileNotFoundError(
            f"{repo_root()}/data_cache/xenium_gbm/cells.csv.gz is missing")

        def explode(sid):
            raise boom

        monkeypatch.setattr(cache_mod, "describe_section", explode)
        rows = cache_mod.list_sections()
        assert rows, "list_sections must still enumerate, never raise"
        for row in rows:
            assert row["available"] is False
            reason = row["reason"]
            assert repo_root() not in reason
            assert "FileNotFoundError" not in reason
            # ...but the operator can still find the real failure.
            assert "reference" in reason

    def test_the_real_detail_goes_to_stderr(self, monkeypatch, capsys):
        def explode(sid):
            raise FileNotFoundError(f"{repo_root()}/data_cache/gone")

        monkeypatch.setattr(cache_mod, "describe_section", explode)
        reason = cache_mod.list_sections()[0]["reason"]
        err = capsys.readouterr().err
        ref = reason.split("reference ")[1].split()[0].rstrip(".,")
        assert f"[vivum:{ref}]" in err
        assert "FileNotFoundError" in err
        assert repo_root() in err

    def test_a_deliberate_refusal_is_never_redacted(self, tmp_path):
        """The distinction that matters. The UI renders these verbatim so it
        can never disagree with the CLI about a sample's limits."""
        from app.server import ingest_registry

        out = ingest_registry.scan(str(tmp_path))
        assert "no recognised layout" in out["error"]
        assert "reference" not in out["error"]


class TestAutoGeneratedApiKey:
    """app/server/__main__.py's decision of whether to require an API key,
    exercised directly against ``main()`` with uvicorn.run mocked out (no
    real socket bind) -- distinct from tests/test_cli.py's TestStudioDelegates
    -To-AppServer, which only checks argv construction and mocks main() itself
    away entirely."""

    @pytest.fixture(autouse=True)
    def _isolate_env_main_writes(self):
        """``main()`` WRITES into ``os.environ``. Several tests here call it
        without monkeypatch having recorded the variables first, and
        ``monkeypatch.delenv(raising=False)`` records nothing when the
        variable is absent -- so main()'s value survived the test and leaked
        into every app built by every later module. Caught live:
        TMESIM_STUDIO_ALLOWED_HOSTS leaking made TrustedHostMiddleware 400
        every TestClient request in test_ingest/test_projects/test_protocol.
        """
        names = ("TMESIM_STUDIO_API_KEY", "TMESIM_STUDIO_ALLOWED_HOSTS")
        saved = {n: os.environ.get(n) for n in names}
        try:
            yield
        finally:
            for n, v in saved.items():
                if v is None:
                    os.environ.pop(n, None)
                else:
                    os.environ[n] = v

    def _run_main(self, monkeypatch, argv):
        import app.server.__main__ as server_main

        # main() does `import uvicorn` locally (so importing __main__ doesn't
        # require the 'app' extra), so patch the real module in sys.modules
        # rather than an attribute on server_main.
        import sys
        import types
        fake_uvicorn = types.ModuleType("uvicorn")
        fake_uvicorn.run = lambda *a, **k: None
        monkeypatch.setitem(sys.modules, "uvicorn", fake_uvicorn)
        # main() WRITES these into os.environ. A bare delenv(raising=False)
        # records nothing when the variable is absent, so monkeypatch's undo
        # would leave main()'s value behind for the rest of the session --
        # which leaked TMESIM_STUDIO_ALLOWED_HOSTS into every app built by
        # every later test module and 400'd them all. setenv first so the
        # restore is always recorded.
        for var in MAIN_WRITES_ENV:
            monkeypatch.setenv(var, os.environ.get(var, ""))
            monkeypatch.delenv(var, raising=False)
        rc = server_main.main(argv)
        return rc, os.environ.get("TMESIM_STUDIO_API_KEY", "")

    def test_loopback_bind_pins_the_host_header_to_the_loopback_names(
            self, monkeypatch):
        """DNS rebinding: a hostile domain that resolves to 127.0.0.1 becomes
        same-origin with this server, which no CORS policy can stop. On the
        default loopback bind the only legitimate Host values are known."""
        self._run_main(monkeypatch, ["--host", "127.0.0.1"])
        assert os.environ["TMESIM_STUDIO_ALLOWED_HOSTS"] == "localhost,127.0.0.1"

    def test_a_non_loopback_bind_does_not_guess_the_operators_hostname(
            self, monkeypatch):
        """We do not know what name a remote deployment is reached by, and
        refusing it would break the deployment outright. Left to the operator."""
        self._run_main(monkeypatch, ["--host", "0.0.0.0"])
        assert os.environ.get("TMESIM_STUDIO_ALLOWED_HOSTS", "") == ""

    def test_an_explicit_allowed_hosts_is_never_clobbered(self, monkeypatch):
        """A Docker ENV / systemd / shell export wins, exactly as it already
        does for TMESIM_STUDIO_API_KEY."""
        import sys
        import types

        import app.server.__main__ as server_main

        fake_uvicorn = types.ModuleType("uvicorn")
        fake_uvicorn.run = lambda *a, **k: None
        monkeypatch.setitem(sys.modules, "uvicorn", fake_uvicorn)
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "")
        monkeypatch.setenv("TMESIM_STUDIO_ALLOWED_HOSTS", "studio.internal")
        assert server_main.main(["--host", "127.0.0.1"]) == 0
        assert os.environ["TMESIM_STUDIO_ALLOWED_HOSTS"] == "studio.internal"

    def test_loopback_host_sets_no_key_by_default(self, monkeypatch, capsys):
        rc, key = self._run_main(monkeypatch, ["--host", "127.0.0.1"])
        assert rc == 0
        assert key == ""
        out = capsys.readouterr().out
        # The DEFAULT is unchanged: loopback still comes up with no key and
        # still serves. What changed is that it now says so -- this used to
        # assert the string "api key" was absent entirely, i.e. that the
        # banner was silent about it, which was the defect (see
        # _no_auth_warning). The real invariant is that no SECRET is printed,
        # because there is none to print.
        assert "X-API-Key:" not in out
        assert "generated" not in out
        assert "  api key     NONE" in out

    def test_loopback_with_no_key_warns_that_a_tunnel_defeats_the_bind(
            self, monkeypatch, capsys):
        """The auto-generate rule reads the bind address, and a loopback
        service behind an ngrok/Cloudflare tunnel or a reverse proxy is
        world-reachable with `is_loopback` still true -- so no key is
        generated and nothing used to be said. Making that loud is the whole
        fix; the default deliberately did NOT change."""
        import app.server.__main__ as server_main

        rc, key = self._run_main(monkeypatch, ["--host", "127.0.0.1"])
        assert rc == 0
        assert key == ""
        out = capsys.readouterr().out
        assert server_main._no_auth_warning("127.0.0.1") in out
        assert "UNAUTHENTICATED" in out
        # names the failure mode, the exposed primitive, and the fix
        assert "tunnel" in out
        assert "reverse proxy" in out
        assert "/api/sections/scan" in out
        assert "TMESIM_STUDIO_API_KEY" in out

    def test_the_no_key_warning_is_absent_once_a_key_is_in_force(
            self, monkeypatch, capsys):
        """It must be a signal, not banner noise: an operator who has done the
        right thing never sees it."""
        rc, key = self._run_main(
            monkeypatch, ["--host", "127.0.0.1", "--api-key", "a-real-secret"])
        assert rc == 0
        assert key == "a-real-secret"
        out = capsys.readouterr().out
        assert "UNAUTHENTICATED" not in out
        assert "NO API KEY" not in out

    def test_nonloopback_host_autogenerates_a_key(self, monkeypatch, capsys):
        rc, key = self._run_main(monkeypatch, ["--host", "0.0.0.0"])
        assert rc == 0
        assert len(key) > 20  # secrets.token_urlsafe(24) -- a real random secret
        out = capsys.readouterr().out
        assert key in out
        assert "generated" in out

    def test_nonloopback_host_with_explicit_key_uses_it_verbatim(self, monkeypatch, capsys):
        rc, key = self._run_main(monkeypatch, ["--host", "0.0.0.0", "--api-key", "my-fixed-key"])
        assert rc == 0
        assert key == "my-fixed-key"
        assert "my-fixed-key" in capsys.readouterr().out

    def test_nonloopback_host_with_explicit_empty_key_disables_auth(self, monkeypatch, capsys):
        rc, key = self._run_main(monkeypatch, ["--host", "0.0.0.0", "--api-key", ""])
        assert rc == 0
        assert key == ""
        assert "NONE" in capsys.readouterr().out

    def test_preexisting_env_key_is_respected_not_clobbered(self, monkeypatch, capsys):
        """Regression test: an earlier version of this logic unconditionally
        set TMESIM_STUDIO_API_KEY based on --host/--api-key alone, silently
        erasing a key someone had already exported (Docker ENV, systemd, a
        shell export) whenever --host defaulted to loopback and --api-key
        wasn't passed. Caught live: `TMESIM_STUDIO_API_KEY=x python -m
        app.server` served completely unauthenticated despite the env var."""
        import app.server.__main__ as server_main
        import sys
        import types

        fake_uvicorn = types.ModuleType("uvicorn")
        fake_uvicorn.run = lambda *a, **k: None
        monkeypatch.setitem(sys.modules, "uvicorn", fake_uvicorn)
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "pre-existing-from-environment")

        rc = server_main.main([])  # default --host 127.0.0.1, no --api-key
        assert rc == 0
        assert os.environ["TMESIM_STUDIO_API_KEY"] == "pre-existing-from-environment"


class TestReconCacheSizeIsValidated:
    """``TMESIM_RECON_CACHE_SIZE`` used to be
    ``int(os.environ.get("TMESIM_RECON_CACHE_SIZE", "4"))`` inline at import,
    which failed in both directions at once: an empty or non-numeric value
    killed the process at import with a bare ``ValueError: invalid literal for
    int() with base 10: ''`` and a traceback naming ``app/server/cache.py``
    (in Docker, a restart loop with no diagnostic), while ``0`` and ``-1`` were
    accepted silently and disable the cache -- the eviction loop drops the
    entry it has just inserted, so every run rebuilds a 1 s to ~90 s
    reconstruction while ``maxsize=0`` is reported as healthy.

    The contract now: unset is the documented default and is silent, a usable
    value is used verbatim, and anything else is a diagnostic on stderr naming
    the variable, the value and the valid range, followed by the default. See
    docs/DEPLOYMENT.md 4.
    """

    def test_unset_is_the_documented_default_and_says_nothing(self, capsys):
        assert cache_mod.recon_cache_size(None) == cache_mod.RECON_CACHE_SIZE_DEFAULT
        assert cache_mod.RECON_CACHE_SIZE_DEFAULT == 4
        assert capsys.readouterr().err == ""

    @pytest.mark.parametrize("raw,expected", [("1", 1), ("4", 4), ("32", 32),
                                              ("  7 ", 7)])
    def test_a_usable_value_is_used_verbatim_and_says_nothing(
            self, raw, expected, capsys):
        assert cache_mod.recon_cache_size(raw) == expected
        assert capsys.readouterr().err == ""

    @pytest.mark.parametrize("raw", ["", "   ", "abc", "4.5", "0", "-1", "-100"])
    def test_an_unusable_value_diagnoses_and_falls_back(self, raw, capsys):
        """The whole point: it must not raise, and it must not be silent."""
        assert cache_mod.recon_cache_size(raw) == 4
        err = capsys.readouterr().err
        # Names the variable, the offending value and the valid range, so the
        # operator can fix it from the message alone.
        assert "TMESIM_RECON_CACHE_SIZE" in err
        assert repr(raw) in err
        assert ">= 1" in err
        assert "Continuing with 4" in err

    def test_the_empty_string_no_longer_kills_the_process_at_import(
            self, monkeypatch, capsys):
        """`-e TMESIM_RECON_CACHE_SIZE=` on a docker run is enough to produce
        the old ValueError, and a container that cannot import its own cache
        module restarts forever."""
        monkeypatch.setenv("TMESIM_RECON_CACHE_SIZE", "")
        size = cache_mod.recon_cache_size()   # must not raise
        assert size == 4
        assert "not usable" in capsys.readouterr().err

    def test_zero_is_refused_rather_than_silently_disabling_the_cache(self):
        """A capacity of 0 is not a small cache, it is a cache that evicts
        every entry the moment it is built."""
        assert cache_mod.recon_cache_size("0") == 4

    def test_the_live_cache_was_built_through_the_validator(self):
        assert cache_mod.RECON_CACHE.maxsize >= 1


class TestApiKeyBannerNamesItsSource:
    """The banner's ``[...]`` label is precisely the line an operator reads to
    work out which of the three sources won. It printed ``[set via --api-key]``
    for a key that came from ``TMESIM_STUDIO_API_KEY`` with no flag passed --
    naming the one source that had NOT supplied it, in the one situation where
    the answer is not obvious."""

    def _banner(self, monkeypatch, capsys, argv, env_key=None):
        import sys as _sys
        import types

        import app.server.__main__ as server_main

        fake_uvicorn = types.ModuleType("uvicorn")
        fake_uvicorn.run = lambda *a, **k: None
        monkeypatch.setitem(_sys.modules, "uvicorn", fake_uvicorn)
        for var in MAIN_WRITES_ENV:
            monkeypatch.setenv(var, os.environ.get(var, ""))
            monkeypatch.delenv(var, raising=False)
        if env_key is not None:
            monkeypatch.setenv("TMESIM_STUDIO_API_KEY", env_key)
        assert server_main.main(argv) == 0
        return capsys.readouterr().out

    def test_a_key_from_the_environment_names_the_environment(
            self, monkeypatch, capsys):
        out = self._banner(monkeypatch, capsys, ["--host", "127.0.0.1"],
                           env_key="from-the-environment")
        assert "from-the-environment  [set via TMESIM_STUDIO_API_KEY]" in out
        assert "set via --api-key" not in out

    def test_an_explicit_flag_still_names_the_flag_and_still_wins(
            self, monkeypatch, capsys):
        """--api-key beats an exported key; the banner must say which value is
        actually in force AND where it came from."""
        out = self._banner(
            monkeypatch, capsys,
            ["--host", "127.0.0.1", "--api-key", "from-the-flag"],
            env_key="from-the-environment")
        assert "from-the-flag  [set via --api-key]" in out
        assert "from-the-environment" not in out

    def test_a_generated_key_still_says_generated(self, monkeypatch, capsys):
        out = self._banner(monkeypatch, capsys, ["--host", "0.0.0.0"])
        assert "[generated (not persisted -- save it now)]" in out


class TestIpv6LoopbackIsExcludedFromTheHostAllowlist:
    """docs/DEPLOYMENT.md 4 used to say an ``::1`` bind got the automatic
    loopback allowlist and therefore 400'd every request. That was true, and
    was a real bug: Starlette matches on ``host.split(":")[0]``, so a browser's
    ``Host: [::1]:8000`` splits to ``[`` and no allowlist entry can ever match
    it -- the middleware rejected 100 % of its own traffic.

    IPv6 loopback is now excluded from the fill-in, so it serves. The cost,
    which the doc must keep stating, is that it then has NO DNS-rebinding
    protection, unlike the IPv4 loopback default."""

    def _run(self, monkeypatch, host):
        import sys as _sys
        import types

        import app.server.__main__ as server_main

        fake_uvicorn = types.ModuleType("uvicorn")
        fake_uvicorn.run = lambda *a, **k: None
        monkeypatch.setitem(_sys.modules, "uvicorn", fake_uvicorn)
        for var in MAIN_WRITES_ENV:
            monkeypatch.setenv(var, os.environ.get(var, ""))
            monkeypatch.delenv(var, raising=False)
        assert server_main.main(["--host", host]) == 0
        return os.environ.get("TMESIM_STUDIO_ALLOWED_HOSTS", "")

    def test_ipv6_loopback_gets_no_allowlist(self, monkeypatch):
        assert self._run(monkeypatch, "::1") == ""

    def test_ipv4_loopback_still_gets_one(self, monkeypatch):
        assert self._run(monkeypatch, "127.0.0.1") == "localhost,127.0.0.1"

    def test_an_ipv6_bind_still_honours_an_explicit_allowlist(self, monkeypatch):
        """The documented workaround: set it yourself, to a NAME. A bracketed
        literal still cannot match, for the same split(':') reason."""
        import sys as _sys
        import types

        import app.server.__main__ as server_main

        fake_uvicorn = types.ModuleType("uvicorn")
        fake_uvicorn.run = lambda *a, **k: None
        monkeypatch.setitem(_sys.modules, "uvicorn", fake_uvicorn)
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "")
        monkeypatch.setenv("TMESIM_STUDIO_ALLOWED_HOSTS", "localhost")
        assert server_main.main(["--host", "::1"]) == 0
        assert os.environ["TMESIM_STUDIO_ALLOWED_HOSTS"] == "localhost"


class TestAccountAndLimits:
    """``/api/limits`` and ``/api/account`` — the deployment's own surface.

    These two routes exist because everything else in this API answers a
    question about TISSUE, and a user also has questions about the INSTANCE:
    which key am I sending, what will this server refuse, what has it done, and
    is anyone charging me. See ``app/server/account.py``.

    The tests below are mostly about the two ways this surface could lie:

    * the ceiling the spec editor shows a user before a run could drift away
      from the ceiling ``validate_spec`` actually refuses with, leaving the
      editor confidently green on a spec the server rejects. The client cannot
      call ``validate_spec`` on every keystroke, so it reimplements the
      formula, so there ARE two copies and the only question is whether
      anything checks them against each other. ``test_limits_formula_matches_
      validate`` is that check.
    * the billing block could read as though billing existed. It is
      scaffolding for an integration that has not been built, and a scaffold
      that does not say so is a claim.
    """

    _K = "account-test-key"
    _H = {"X-API-Key": _K}

    def _make_client(self, monkeypatch, root, key=_K):
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", key)
        # __main__.py records the provenance; a test that builds the app
        # directly has to clear it or it inherits whatever ran last.
        monkeypatch.delenv("TMESIM_STUDIO_API_KEY_SOURCE", raising=False)
        monkeypatch.delenv("TMESIM_BILLING_PROVIDER", raising=False)
        manager = JobManager(root=str(root))
        return TestClient(create_app(manager)), manager

    # -- the gate ---------------------------------------------------------- #
    def test_both_routes_are_gated(self, monkeypatch, tmp_path):
        """Neither is a liveness route. A fingerprint, a usage count and a
        refusal threshold are operational detail about a machine."""
        c, mgr = self._make_client(monkeypatch, tmp_path)
        try:
            assert c.get("/api/limits").status_code == 401
            assert c.get("/api/account").status_code == 401
            assert c.get("/api/limits", headers=self._H).status_code == 200
            assert c.get("/api/account", headers=self._H).status_code == 200
        finally:
            mgr.shutdown()

    # -- the formula, checked against the refusal it predicts --------------- #
    def test_limits_formula_matches_validate(self, monkeypatch, tmp_path):
        """The published formula must predict the actual refusal boundary.

        Method: lower the ceiling to a value whose critical ``dx`` lands in a
        range the other dx rules do not touch, compute the boundary from the
        formula ``/api/limits`` publishes, then ask ``/api/validate`` on either
        side of it. A spec 2% coarser than the boundary must NOT be refused for
        grid size; a spec 2% finer must be.

        Asserted on the GRID error specifically, not on ``ok``: both of these
        specs also trip the coarse-mesh warnings, and a test that just watched
        ``ok`` flip would pass for the wrong reason the moment an unrelated
        rule started erroring.
        """
        import math

        import app.server.api as api_mod

        monkeypatch.setattr(api_mod, "_MAX_GRID_CELLS", 100_000)
        # The source is RESOLVED at import alongside the value now, so patching
        # the value alone would leave the payload reporting "default" for a
        # ceiling that is nothing of the kind.
        monkeypatch.setattr(api_mod, "_MAX_GRID_CELLS_SOURCE",
                            "TMESIM_MAX_GRID_CELLS")
        c, mgr = self._make_client(monkeypatch, tmp_path)
        try:
            limits = c.get("/api/limits", headers=self._H).json()
            assert limits["max_grid_cells"] == 100_000
            assert limits["max_grid_cells_source"] == "TMESIM_MAX_GRID_CELLS"

            # The client reimplements this string. Asserting the literal means
            # editing the published formula forces editing this test, which
            # forces a human to look at whether the client still agrees.
            assert limits["grid_cells_formula"] == \
                "(extent_um[0] / dx_um) * (extent_um[1] / dx_um)"

            row = next(r for r in c.get("/api/sections", headers=self._H).json()
                       if r["id"] == SECTION)
            ex, ey = row["extent_um"]

            def cells(dx):  # the published formula, implemented here
                return (ex / dx) * (ey / dx)

            ceiling = limits["max_grid_cells"]
            dx_crit = math.sqrt(ex * ey / ceiling)
            dx_ok, dx_bad = dx_crit * 1.02, dx_crit * 0.98
            # The premise of the test, checked rather than assumed.
            assert cells(dx_ok) < ceiling < cells(dx_bad)
            assert dx_ok > limits["min_dx_um"] and dx_bad > limits["min_dx_um"]

            def grid_errors(dx):
                body = {**default_spec(), "section": SECTION, "dx_um": dx}
                verdict = c.post("/api/validate", json=body,
                                 headers=self._H).json()
                return [e for e in verdict["errors"] if "ceiling" in e]

            assert grid_errors(dx_ok) == [], (
                "the server refused a spec the published formula says is "
                "under the ceiling")
            bad = grid_errors(dx_bad)
            assert bad, ("the server accepted a spec the published formula "
                         "says is over the ceiling")
            assert f"{ceiling:,}" in bad[0]
        finally:
            mgr.shutdown()

    # -- the key ----------------------------------------------------------- #
    def test_account_never_echoes_the_key(self, monkeypatch, tmp_path):
        """A fingerprint is an identifier. The secret itself must not appear
        anywhere in the payload, at any nesting depth."""
        c, mgr = self._make_client(monkeypatch, tmp_path)
        try:
            raw = c.get("/api/account", headers=self._H).text
            assert self._K not in raw
            body = c.get("/api/account", headers=self._H).json()
            assert body["auth"]["key_fingerprint"]
            assert body["auth"]["key_fingerprint"] == \
                body["auth"]["caller_fingerprint"]
            assert body["auth"]["matches"] is True
        finally:
            mgr.shutdown()

    def test_fingerprint_distinguishes_two_keys(self):
        """The panel's whole job is telling a user WHICH saved key this
        browser is sending, so two keys must not share a fingerprint."""
        from app.server.account import key_fingerprint

        assert key_fingerprint("alpha") != key_fingerprint("beta")
        assert key_fingerprint("alpha") == key_fingerprint("alpha")
        assert key_fingerprint("") is None
        assert key_fingerprint("   ") is None
        # Not a truncated bare digest of the secret: a bare sha256 of a
        # low-entropy operator-chosen key is a dictionary lookup.
        import hashlib
        assert key_fingerprint("alpha") != \
            hashlib.sha256(b"alpha").hexdigest()[:12]

    def test_key_source_is_reported_and_unknown_stays_unknown(
            self, monkeypatch, tmp_path):
        """An operator who has LOST the key needs to know whether they chose
        it or the process invented one at boot and printed it once."""
        c, mgr = self._make_client(monkeypatch, tmp_path)
        try:
            # Nothing recorded it -> it was already in the environment.
            assert c.get("/api/account", headers=self._H
                         ).json()["auth"]["source"] == "environment"
            for recorded, reported in (("generated", "generated-at-startup"),
                                       ("flag", "--api-key"),
                                       ("environment", "environment"),
                                       ("nonsense", "unknown")):
                monkeypatch.setenv("TMESIM_STUDIO_API_KEY_SOURCE", recorded)
                assert c.get("/api/account", headers=self._H
                             ).json()["auth"]["source"] == reported
        finally:
            mgr.shutdown()

    def test_main_records_which_key_source_won(self, monkeypatch, tmp_path):
        """The banner works this out to print a label; the panel needs it
        weeks later, when the banner is scrollback nobody has."""
        import sys as _sys
        import types

        import app.server.__main__ as server_main

        fake_uvicorn = types.ModuleType("uvicorn")
        fake_uvicorn.run = lambda *a, **k: None
        monkeypatch.setitem(_sys.modules, "uvicorn", fake_uvicorn)
        monkeypatch.setenv("TMESIM_STUDIO_DATA_ROOTS", str(tmp_path))
        # main() writes these into os.environ DIRECTLY, and a bare
        # delenv(raising=False) records NOTHING when the variable is absent --
        # so main()'s value survives teardown and leaks into every app built
        # afterwards. setenv first, so the restore is always recorded. This is
        # the same loop _run_main uses over the same shared tuple; a
        # hand-rolled second copy here is what let the leak in the first time.
        for var in MAIN_WRITES_ENV:
            monkeypatch.setenv(var, os.environ.get(var, ""))
            monkeypatch.delenv(var, raising=False)

        monkeypatch.delenv("TMESIM_STUDIO_API_KEY", raising=False)
        monkeypatch.delenv("TMESIM_STUDIO_API_KEY_SOURCE", raising=False)
        assert server_main.main(["--host", "0.0.0.0"]) == 0
        assert os.environ["TMESIM_STUDIO_API_KEY_SOURCE"] == "generated"

        monkeypatch.delenv("TMESIM_STUDIO_API_KEY_SOURCE", raising=False)
        assert server_main.main(["--api-key", "chosen-on-the-flag"]) == 0
        assert os.environ["TMESIM_STUDIO_API_KEY_SOURCE"] == "flag"

        monkeypatch.delenv("TMESIM_STUDIO_API_KEY_SOURCE", raising=False)
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "already-exported")
        assert server_main.main([]) == 0
        assert os.environ["TMESIM_STUDIO_API_KEY_SOURCE"] == "environment"

    # -- usage and billing -------------------------------------------------- #
    def test_usage_says_it_is_an_inventory_not_a_meter(
            self, monkeypatch, tmp_path):
        """Deleting a run decrements these counts. Anything that reads like a
        meter here would be the first half of a billing claim."""
        c, mgr = self._make_client(monkeypatch, tmp_path)
        try:
            usage = c.get("/api/account", headers=self._H).json()["usage"]
            assert "not a meter" in usage["kind"]
            assert usage["attribution"].startswith("none")
            assert usage["runs_total"] == sum(usage["runs_by_status"].values())
        finally:
            mgr.shutdown()

    def test_billing_is_never_enabled_and_says_why(
            self, monkeypatch, tmp_path):
        """Naming a provider must not be mistakable for having built one.

        With TMESIM_BILLING_PROVIDER set and no implementation mounted, the
        block reports configured-but-not-implemented and ``enabled`` stays
        False -- a state an operator will notice, rather than a silent
        half-integration.
        """
        c, mgr = self._make_client(monkeypatch, tmp_path)
        try:
            b = c.get("/api/account", headers=self._H).json()["billing"]
            assert b["enabled"] is False
            assert b["status"] == "not-configured"
            assert b["provider"] is None
            assert b["blocked_on"], "the preconditions must be stated"
            assert any("shared key" in r for r in b["blocked_on"])

            monkeypatch.setenv("TMESIM_BILLING_PROVIDER", "some-psp")
            b = c.get("/api/account", headers=self._H).json()["billing"]
            assert b["enabled"] is False, \
                "naming a provider must never flip billing on"
            assert b["status"] == "configured-but-not-implemented"
            assert b["provider"] == "some-psp"
        finally:
            mgr.shutdown()


class TestPrunedPayloads:
    """A payload deleted to reclaim disk must not be reported as one that was
    never produced.

    ``scripts/prune_run_payloads.py`` removes ``timelapse.bin`` from old runs
    (measured 2026-08-30: 1,164 MB across 73 runs, the largest thing this
    repository generates). Without a record of the prune, ``_binary`` falls
    through to "no drug time-lapse was captured for this run (timelapse_frames
    was 0, or this run predates the feature)" -- a confident sentence in which
    every clause is false for a run that captured one and had it pruned.

    That is this project's signature defect, correct reasoning pointed at the
    wrong object, and it is the ONLY reason the prune writes a record at all.
    So the record has to survive a restart and the endpoint has to prefer it.
    """

    def _restored(self, tmp_path, doc):
        (tmp_path / "r_00001").mkdir()
        (tmp_path / "r_00001" / "result.json").write_text(json.dumps(doc))
        manager = JobManager(root=str(tmp_path))
        return TestClient(create_app(manager)), manager

    def _doc(self, **extra):
        doc = {"run_id": "r_00001", "status": "done",
               "spec": {"section": SECTION, "dx_um": 30.0, "drug": "cisplatin"},
               "timelapse": {"n_frames": 12}}
        doc.update(extra)
        return doc

    def test_a_pruned_payload_is_410_and_says_it_existed(self, tmp_path):
        c, mgr = self._restored(tmp_path, self._doc(pruned_payloads={
            "timelapse": {"at": "2026-08-30", "bytes": 27476084,
                          "reason": "Pruned to reclaim disk."}}))
        try:
            r = c.get("/api/runs/r_00001/timelapse")
            assert r.status_code == 410, "a pruned payload is Gone, not Not Found"
            reason = r.json()["detail"]["reason"]
            assert "DID produce" in reason
            assert "2026-08-30" in reason
            assert "27,476,084" in reason
            # The false explanation must NOT be the one that comes out.
            assert "timelapse_frames was 0" not in reason
            assert "predates the feature" not in reason
        finally:
            mgr.shutdown()

    def test_without_a_prune_record_the_old_reason_still_applies(self, tmp_path):
        """The prune branch must not swallow the genuine never-produced case:
        a run that really captured nothing still gets its own reason."""
        c, mgr = self._restored(tmp_path, self._doc(
            timelapse={"reason": "timelapse_frames was 0"}))
        try:
            r = c.get("/api/runs/r_00001/timelapse")
            assert r.status_code == 404
            assert "timelapse_frames was 0" in r.json()["detail"]["reason"]
        finally:
            mgr.shutdown()

    def test_a_pruned_run_does_not_also_report_the_payload_as_available(
            self, tmp_path):
        """The run document and the payload endpoint must not disagree.

        MEASURED on the live store 2026-08-30: 51 of 54 pruned runs reported
        `timelapse.available: true, n_frames: 24, reason: null` NEXT TO their
        own `pruned_payloads` record, because `result_view` copies `self.body`
        wholesale. So `GET /api/runs/{id}` asserted, in one payload, both that
        the file was there and that it had been deleted -- and a client that
        trusts the document never asks for the payload and never sees the
        correct 410. Found by a claim-surface audit, not by the prune's own
        tests.
        """
        c, mgr = self._restored(tmp_path, self._doc(
            timelapse={"available": True, "n_frames": 24, "reason": None},
            pruned_payloads={"timelapse": {"at": "2026-08-30",
                                           "bytes": 27476084}}))
        try:
            tl = c.get("/api/runs/r_00001").json()["timelapse"]
            assert tl["available"] is False, \
                "the document still claims a deleted payload is available"
            assert tl["pruned"] is True
            assert tl["reason"] and "deleted on 2026-08-30" in tl["reason"]
            assert "27,476,084" in tl["reason"]
            # The run DID capture 24 frames. That fact is history and stays.
            assert tl["n_frames"] == 24
        finally:
            mgr.shutdown()

    def test_an_unpruned_run_is_left_alone(self, tmp_path):
        c, mgr = self._restored(tmp_path, self._doc(
            timelapse={"available": True, "n_frames": 24, "reason": None}))
        try:
            tl = c.get("/api/runs/r_00001").json()["timelapse"]
            assert tl["available"] is True
            assert "pruned" not in tl
        finally:
            mgr.shutdown()

    def test_the_prune_record_survives_a_restart(self, tmp_path):
        """``JobManager._restore`` copies a whitelist of keys out of
        result.json. A key missing from that whitelist is silently dropped, and
        the endpoint then falls back to the false reason on the NEXT boot --
        green today, lying tomorrow."""
        (tmp_path / "r_00001").mkdir()
        (tmp_path / "r_00001" / "result.json").write_text(json.dumps(
            self._doc(pruned_payloads={"timelapse": {"at": "2026-08-30",
                                                     "bytes": 1}})))
        mgr = JobManager(root=str(tmp_path))
        try:
            job = mgr.get("r_00001")
            assert job is not None
            assert "pruned_payloads" in (job.body or {}), \
                "the prune record did not survive _restore's key whitelist"
        finally:
            mgr.shutdown()


class TestHostMemoryAndOomExposure:
    """``app/server/account.py``: what the host grants, and whether this
    server's grid ceiling still permits the spec measured to kill it.

    The measurement behind this is TWO POINTS AT ONE INSTANCE SIZE (46,070
    voxels peaks at 220 MiB, 286,200 voxels OOM-kills, both on 512 MB). That
    supports a COMPARISON and not a law -- 220 MiB over 46,070 voxels is ~5 KiB
    a voxel against the 8 bytes a float64 field costs, so the peak is dominated
    by the cell population and not the grid. The first test below is what stops
    a later change from quietly fitting a line through those two points.
    """

    def test_no_ceiling_is_ever_derived_from_detected_memory(
            self, monkeypatch):
        """Detecting 128 MB must not lower the ceiling by itself.

        A derived ceiling would be a bytes-per-voxel law fitted to two points
        of a quantity that is mostly not the independent variable: a confident
        number with nothing under it. The code reports the exposure and leaves
        the ceiling to the operator, and this pins that.
        """
        import importlib

        import app.server.api as api_mod

        monkeypatch.delenv("TMESIM_MAX_GRID_CELLS", raising=False)
        monkeypatch.setenv("TMESIM_HOST_MEMORY_MB", "128")
        importlib.reload(api_mod)
        try:
            assert api_mod._MAX_GRID_CELLS == api_mod._DEFAULT_MAX_GRID_CELLS, (
                "the ceiling moved with detected memory — that is a derived "
                "law fitted to two measured points, which this repository "
                "does not do")
        finally:
            monkeypatch.delenv("TMESIM_HOST_MEMORY_MB", raising=False)
            importlib.reload(api_mod)

    def test_the_warning_fires_only_when_the_ceiling_permits_the_kill(
            self, monkeypatch):
        from app.server.account import MEASURED_OOM_VOXELS, oom_warning

        monkeypatch.setenv("TMESIM_HOST_MEMORY_MB", "512")
        # default ceiling on a small box: the spec that killed a container is
        # still allowed through, so this must be loud.
        loud = oom_warning(40_000_000, True)
        assert loud and "OOM-KILLED" in loud
        assert f"{MEASURED_OOM_VOXELS:,}" in loud
        assert "TMESIM_MAX_GRID_CELLS" in loud
        # a ceiling BELOW the measured kill already refuses it: silent.
        assert oom_warning(MEASURED_OOM_VOXELS - 1, False) is None
        # exactly at it still permits it (the refusal is `>`), so still loud.
        assert oom_warning(MEASURED_OOM_VOXELS, False) is not None

    def test_undetected_memory_is_not_reported_as_unlimited(self, monkeypatch):
        """A host can be capped by something this process cannot see. Unknown
        must read as unknown, never as headroom."""
        from app.server.account import memory_limit, oom_exposure

        monkeypatch.delenv("TMESIM_HOST_MEMORY_MB", raising=False)
        monkeypatch.setattr("app.server.account._CGROUP_V2", "/nonexistent/v2")
        monkeypatch.setattr("app.server.account._CGROUP_V1", "/nonexistent/v1")
        m = memory_limit()
        assert m["bytes"] is None
        assert m["source"] == "not declared"
        ex = oom_exposure(40_000_000, True)
        assert ex["host_at_or_below_measured_instance"] is None, \
            "unknown memory must be None, not False"

    def test_cgroup_v2_max_and_a_v1_sentinel_both_read_as_no_limit(
            self, monkeypatch, tmp_path):
        from app.server.account import memory_limit

        monkeypatch.delenv("TMESIM_HOST_MEMORY_MB", raising=False)
        v2 = tmp_path / "memory.max"
        v1 = tmp_path / "limit_in_bytes"
        monkeypatch.setattr("app.server.account._CGROUP_V2", str(v2))
        monkeypatch.setattr("app.server.account._CGROUP_V1", str(v1))

        v2.write_text("max\n")                      # cgroup v2 "no limit"
        v1.write_text("9223372036854771712\n")      # cgroup v1 sentinel
        assert memory_limit()["bytes"] is None

        v2.write_text("536870912\n")
        got = memory_limit()
        assert got["bytes"] == 536870912
        assert "cgroup v2" in got["source"]

    def test_limits_publishes_the_exposure(self, monkeypatch, tmp_path):
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "k")
        manager = JobManager(root=str(tmp_path))
        c = TestClient(create_app(manager))
        try:
            ex = c.get("/api/limits", headers={"X-API-Key": "k"}
                       ).json()["oom_exposure"]
            assert ex["measured"]["oom_voxels"] == 286_200
            assert ex["measured"]["instance_mb"] == 512
            assert "NOT a bytes-per-voxel law" in ex["measured"]["note"]
            assert isinstance(ex["ceiling_permits_measured_oom"], bool)
        finally:
            manager.shutdown()


class TestBannerElseBinding:
    """A regression guard for a two-line edit that parses and inverts meaning.

    Adding a block between ``if api_key:``'s body and its ``else:`` at the
    OUTER indent rebinds the ``else`` to the new ``if``. Python accepts it. The
    result is that a correctly secured server prints "api key NONE -- no
    authentication on any /api/ route", which is the single most alarming line
    the banner can emit and would be false.

    This happened twice in one session while adding the OOM warning below the
    key block, so it gets a test rather than a comment.
    """

    def test_a_secured_server_never_announces_no_authentication(
            self, monkeypatch, tmp_path, capsys):
        import sys as _sys
        import types

        import app.server.__main__ as server_main

        fake_uvicorn = types.ModuleType("uvicorn")
        fake_uvicorn.run = lambda *a, **k: None
        monkeypatch.setitem(_sys.modules, "uvicorn", fake_uvicorn)
        monkeypatch.setenv("TMESIM_STUDIO_DATA_ROOTS", str(tmp_path))
        for var in MAIN_WRITES_ENV:
            monkeypatch.setenv(var, os.environ.get(var, ""))
            monkeypatch.delenv(var, raising=False)
        monkeypatch.setenv("TMESIM_STUDIO_API_KEY", "a-real-secret")
        # PIN THE OOM WARNING OFF. Without this the test passes for the wrong
        # reason: with the bug planted the banner becomes
        # `if _oom: print(_oom) else: print("api key NONE...")`, and on a host
        # with an unconfigured ceiling `_oom` is truthy, so the else never runs
        # and the assertion below never sees the damage. Caught by planting the
        # real bug and watching this test stay green. A ceiling below the
        # measured kill silences the warning, which is exactly the state that
        # exercises the else.
        monkeypatch.setenv("TMESIM_MAX_GRID_CELLS", "100000")
        import importlib

        import app.server.api as api_mod

        # try/finally, not a trailing reload: _MAX_GRID_CELLS is read at
        # IMPORT, so a failed assertion below would otherwise leave every later
        # test in this session running against a 100,000 ceiling.
        importlib.reload(api_mod)
        try:
            assert server_main.main(["--host", "0.0.0.0"]) == 0
            out = capsys.readouterr().out
            assert "a-real-secret" in out, "the banner should report the key"
            assert "no authentication on any /api/ route" not in out, (
                "a server WITH a key announced that it has none — the `else` "
                "of `if api_key:` has been rebound to a later `if`")
        finally:
            monkeypatch.delenv("TMESIM_MAX_GRID_CELLS", raising=False)
            importlib.reload(api_mod)

    def test_the_else_belongs_to_the_key_check_not_a_later_if(self):
        """Read structurally, from the AST, so it cannot be satisfied by a
        coincidence of runtime state."""
        import ast
        import inspect

        import app.server.__main__ as server_main

        src = inspect.getsource(server_main)
        tree = ast.parse(src)
        owners = []
        for node in ast.walk(tree):
            if isinstance(node, ast.If) and node.orelse:
                seg = ast.get_source_segment(src, node.orelse[0]) or ""
                if "no authentication on any /api/ route" in seg:
                    owners.append(ast.get_source_segment(src, node.test))
        assert owners == ["api_key"], (
            f"the no-authentication `else` hangs off {owners!r}, not `api_key`")


class TestCeilingResolution:
    """``_resolve_max_grid_cells`` — where the grid ceiling comes from.

    Four branches, and the interesting one is what `auto` does NOT do. It
    adopts 46,070 voxels: the largest grid MEASURED to survive a 512 MB box.
    That is applying a measurement, not fitting one — the only two points that
    exist are "46,070 survived" and "286,200 died", so the only honest
    automatic ceiling is the one that survived, and it is knowingly
    conservative about everything between them.

    It is OPT-IN. The ceiling does not move on its own when a small cgroup
    limit is detected, because a server that silently refuses a spec the
    operator expects to run, on the strength of a number it read out of /sys,
    is a debugging session nobody asked for. Detection always drives the
    WARNING; it drives the CEILING only by name.
    """

    def _resolve(self, monkeypatch, value, mem=None):
        import app.server.api as api_mod

        if value is None:
            monkeypatch.delenv("TMESIM_MAX_GRID_CELLS", raising=False)
        else:
            monkeypatch.setenv("TMESIM_MAX_GRID_CELLS", value)
        if mem is None:
            monkeypatch.delenv("TMESIM_HOST_MEMORY_MB", raising=False)
        else:
            monkeypatch.setenv("TMESIM_HOST_MEMORY_MB", mem)
        return api_mod._resolve_max_grid_cells()

    def test_unset_is_the_default(self, monkeypatch):
        import app.server.api as api_mod
        assert self._resolve(monkeypatch, None) == \
            (api_mod._DEFAULT_MAX_GRID_CELLS, "default")

    def test_an_explicit_number_is_reported_as_explicit(self, monkeypatch):
        """Even when it happens to equal the default. The source used to be
        INFERRED from the value, so `TMESIM_MAX_GRID_CELLS=40000000` reported
        itself as "default" — an operator who set it deliberately was told the
        server had ignored them."""
        import app.server.api as api_mod
        assert self._resolve(monkeypatch, "123456") == \
            (123456, "TMESIM_MAX_GRID_CELLS")
        assert self._resolve(
            monkeypatch, str(api_mod._DEFAULT_MAX_GRID_CELLS)) == \
            (api_mod._DEFAULT_MAX_GRID_CELLS, "TMESIM_MAX_GRID_CELLS")

    def test_auto_on_a_small_host_adopts_the_measured_survivor(
            self, monkeypatch):
        from app.server.account import MEASURED_OK_VOXELS, MEASURED_OOM_VOXELS

        cells, source = self._resolve(monkeypatch, "auto", mem="512")
        assert cells == MEASURED_OK_VOXELS
        assert source.startswith("auto (")
        # The whole point: below the count that killed a box of this size.
        assert cells < MEASURED_OOM_VOXELS

    def test_auto_on_a_large_or_undeclared_host_says_which_branch_ran(
            self, monkeypatch):
        """`auto` must never hide the fact that it did nothing."""
        import app.server.api as api_mod

        cells, source = self._resolve(monkeypatch, "auto", mem="16000")
        assert cells == api_mod._DEFAULT_MAX_GRID_CELLS
        assert "no small-host limit detected" in source

        cells, source = self._resolve(monkeypatch, "auto", mem=None)
        assert cells == api_mod._DEFAULT_MAX_GRID_CELLS
        assert "no small-host limit detected" in source

    def test_a_typo_refuses_to_start(self, monkeypatch):
        """The dangerous failure is a typo falling back to 40 M on a small
        host, which permits the spec measured to kill it. Loud, not silent."""
        with pytest.raises(SystemExit) as exc:
            self._resolve(monkeypatch, "1O0000")   # letter O, not zero
        assert "neither an integer nor 'auto'" in str(exc.value)
        assert "40,000,000" in str(exc.value)

    def test_detection_alone_never_moves_the_ceiling(self, monkeypatch):
        """The opt-in guarantee, stated as a test: a tiny detected host with
        nothing set must leave the ceiling exactly where it was."""
        import app.server.api as api_mod
        assert self._resolve(monkeypatch, None, mem="128") == \
            (api_mod._DEFAULT_MAX_GRID_CELLS, "default")


class TestEnvVarsAreDocumented:
    """CI for one page of the claim surface.

    ``docs/DEPLOYMENT.md`` §4 opened with "Six. There is no config file; these
    are the whole surface" while the code read TEN. ``TMESIM_MAX_GRID_CELLS``
    had been the difference between a refusal and a dead container for weeks
    and was not on the page at all. Nobody wrote a wrong sentence: the sentence
    was true when it was written and the code moved underneath it, which is
    this project's most persistent defect and the one thing `results/` and
    `docs/` have no CI for.

    This is a small, exact piece of that CI. It does not check that the prose
    is RIGHT — only that no variable the server reads is absent from the page
    that claims to list them all.
    """

    _DOC = os.path.join(_REPO, "docs", "DEPLOYMENT.md")

    #: Scanned trees. `app/server` alone was NOT enough: a claim-surface audit
    #: on 2026-08-30 counted 12 TMESIM_* in the repository against the 10 this
    #: found, and the difference included TMESIM_PROTOCOL_MAX_WINDOW_S, which
    #: lives in tmesim/protocol/dosing.py and is reachable from a running
    #: server through the /api/protocol/* routes. A knob an operator can turn
    #: on a deployment belongs on the deployment page wherever its code sits.
    _SRC_TREES = (("app", "server"), ("tmesim", "protocol"))

    #: Variables that are NOT deployment knobs, each with the reason. An
    #: exclusion needs a justification here rather than a silent narrowing of
    #: the scan, so adding a variable anywhere forces a decision instead of
    #: quietly widening the gap between the page and the code.
    _NOT_DEPLOYMENT = {
        "TMESIM_BAKED": "build-time only (scripts/build_studio_export.py); "
                        "never read by a running server",
    }

    def _vars_in_code(self):
        found = set()
        for parts in self._SRC_TREES:
            d = os.path.join(_REPO, *parts)
            for root, _dirs, files in os.walk(d):
                for name in sorted(files):
                    if not name.endswith(".py"):
                        continue
                    with open(os.path.join(root, name)) as fh:
                        found.update(re.findall(r"TMESIM_[A-Z_]+", fh.read()))
        return found - set(self._NOT_DEPLOYMENT)

    def test_every_env_var_the_server_reads_is_on_the_page(self):
        with open(self._DOC) as fh:
            doc = fh.read()
        missing = sorted(v for v in self._vars_in_code() if v not in doc)
        assert not missing, (
            "docs/DEPLOYMENT.md claims to list the whole environment surface "
            f"and does not mention: {missing}")

    def test_the_stated_count_matches_the_code(self):
        """The number in the opening sentence is itself a claim."""
        with open(self._DOC) as fh:
            doc = fh.read()
        m = re.search(r"## 4\. Environment variables\s+(\w+)\.", doc)
        assert m, "the section no longer opens with a count"
        words = {"four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8,
                 "nine": 9, "ten": 10, "eleven": 11, "twelve": 12}
        stated = words.get(m.group(1).lower())
        assert stated is not None, f"unparsable count {m.group(1)!r}"
        actual = len(self._vars_in_code())
        assert stated == actual, (
            f"docs/DEPLOYMENT.md §4 says {m.group(1)} ({stated}); the server "
            f"reads {actual}: {sorted(self._vars_in_code())}")


class TestPruneCitationProtection:
    """``scripts/prune_run_payloads.py`` must not delete what a record cites.

    Its whole justification is "a reader following a citation must find what
    the citation points at". On 2026-08-30 it deleted `r_00028`, which
    `docs/RUNBOOK_DEPLOY.md:368` cites as part of the range
    ``r_00027`-`r_00029``. The plain id scan matched the two ENDPOINTS, found
    "something" on that line, and looked like it had worked. The middle of a
    range is named nowhere.
    """

    def _fn(self):
        import importlib.util

        path = os.path.join(_REPO, "scripts", "prune_run_payloads.py")
        spec = importlib.util.spec_from_file_location("_prune", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_a_range_in_an_isolated_tree_protects_its_middle(self, tmp_path):
        """Against a TEMP repo, not the real one.

        The first version of this test asserted `r_00028 in cited_run_ids()`
        over the live tree, and it passed with the range logic DELETED --
        because `results/claim_surface_audit_2026-08-30.md`, the audit that
        found the bug, names `r_00028` in prose, and `results/` is a citation
        root. The test was green for a reason with nothing to do with ranges.
        Caught by planting the original bug; it would otherwise have shipped as
        a guard that guards nothing.
        """
        mod = self._fn()
        docs = tmp_path / "docs"
        docs.mkdir()
        (docs / "note.md").write_text(
            "three runs, cited as a range: `r_00027`-`r_00029`, all "
            "temozolomide.\n")
        cited = mod.cited_run_ids(repo=str(tmp_path))
        assert cited == {"r_00027", "r_00028", "r_00029"}, (
            "a range cites its middle, and the middle is named nowhere on the "
            f"line; got {sorted(cited)}")

    def test_the_real_citation_that_was_missed_is_protected_now(self):
        """The live tree, kept as a second check -- weaker than the one above
        (see its docstring) but it is the actual record that broke."""
        mod = self._fn()
        cited = mod.cited_run_ids()
        assert {"r_00027", "r_00028", "r_00029"} <= cited

    def test_the_middle_of_a_range_is_not_findable_by_the_plain_id_scan(self):
        """The premise, checked rather than assumed: if RUN_ID_RE could already
        see `r_00028` in that line, the range logic would be redundant and this
        whole class would be theatre."""
        mod = self._fn()
        line = b"three runs: `r_00027`-`r_00029`, all temozolomide."
        plain = {m.decode() for m in mod.RUN_ID_RE.findall(line)}
        assert plain == {"r_00027", "r_00029"}
        assert "r_00028" not in plain

    @pytest.mark.parametrize("text,expect", [
        (b"`r_00027`-`r_00029`", {"r_00027", "r_00028", "r_00029"}),
        ("r_00040–r_00042".encode(), {"r_00040", "r_00041", "r_00042"}),
        ("r_00050—r_00051".encode(), {"r_00050", "r_00051"}),
        (b"r_00060 to r_00061", {"r_00060", "r_00061"}),
        (b"r_00070--r_00071", {"r_00070", "r_00071"}),
    ])
    def test_every_dash_a_human_types_is_read_as_a_range(self, text, expect):
        mod = self._fn()
        got = set()
        for lo, hi in mod.RANGE_RE.findall(text):
            got.update(f"r_{n:05d}" for n in range(int(lo), int(hi) + 1))
        assert got == expect

    def test_a_descending_range_still_expands(self):
        mod = self._fn()
        pairs = mod.RANGE_RE.findall(b"r_00009-r_00007")
        assert pairs, "a descending range must still be recognised"

    def test_an_absurd_range_is_not_treated_as_a_citation(self):
        """Two unrelated ids on one line must not protect the whole store.
        MAX_RANGE is what stops a false positive costing every run's disk."""
        mod = self._fn()
        assert mod.MAX_RANGE < 100
