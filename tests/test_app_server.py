"""Tests for the tmesim Studio backend (``app/server/``).

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

import os
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
            assert "RANKING" in hit          # what IS defensible
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
        below ``dt * max|v|`` must be a hard error naming the numbers."""
        _run(client, _spec())  # measures and records the speed
        speed = cache_mod.RECON_CACHE.known_speed(
            cache_mod.make_key(SECTION, DX_UM, "discrete", "pde"))
        assert speed and speed > 0

        # dx small enough that dt*max|v|/dx > 1 for this exact key.
        bad_dx = 0.5 * speed * PDE_DT_S
        v = validate_spec(_spec(dx_um=bad_dx))
        assert not v["ok"]
        assert any("CFL" in e and "dt*max|v|/dx" in e for e in v["errors"])
        # And the safe dx is honoured.
        assert validate_spec(_spec(dx_um=2.0 * speed * PDE_DT_S))["ok"]

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

    def test_failed_runs_surface_their_error(self, client, monkeypatch):
        import app.server.jobs as jobs_mod

        def boom(spec, ctx, cache=None):
            raise RuntimeError("synthetic engine failure")

        manager = JobManager(root=str(client._root) + "_fail")
        monkeypatch.setattr(jobs_mod, "execute_run", boom)
        job = manager.submit(_spec())
        deadline = time.time() + 30
        while job.status not in ("failed", "done", "cancelled"):
            if time.time() > deadline:
                raise AssertionError("job never finished")
            time.sleep(0.05)
        assert job.status == "failed"
        assert "synthetic engine failure" in job.error
        assert job.traceback
        manager.shutdown()

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


class TestAutoGeneratedApiKey:
    """app/server/__main__.py's decision of whether to require an API key,
    exercised directly against ``main()`` with uvicorn.run mocked out (no
    real socket bind) -- distinct from tests/test_cli.py's TestStudioDelegates
    -To-AppServer, which only checks argv construction and mocks main() itself
    away entirely."""

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
        monkeypatch.delenv("TMESIM_STUDIO_API_KEY", raising=False)
        rc = server_main.main(argv)
        return rc, os.environ.get("TMESIM_STUDIO_API_KEY", "")

    def test_loopback_host_sets_no_key_by_default(self, monkeypatch, capsys):
        rc, key = self._run_main(monkeypatch, ["--host", "127.0.0.1"])
        assert rc == 0
        assert key == ""
        assert "api key" not in capsys.readouterr().out

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
