"""Tests for the per-cell mechanism attribution and the drug time-lapse.

The load-bearing claim of this feature is that it is a FAITHFUL read-out of the
kill model rather than a plausible-looking story told alongside it. Three things
have to hold for that claim to survive, and each has a test here:

1. **Bit-exactness.** The labels rebuild ``kill_fraction`` with ``==``, not
   ``allclose``, on a real reconstruction — and the engine refuses to return an
   attribution that does not.
2. **Inertness.** Asking for the attribution, or for time-lapse frames, must not
   move a single bit of the simulation. Both are compared against the same run
   with the options off, on ``kill_fraction`` and on the whole AUC field.
3. **No rounding into the nearest bucket.** A survival the bands do not explain
   is labelled ``UNATTRIBUTED`` and counted. The test drives the engine into
   exactly that state (``max_kill_rate < 1``, a suppression the bands do not
   describe) and asserts the count lands there instead of in CHANCE.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from tmesim.cells.pkpd import apply_drug_kill  # noqa: E402
from tmesim.core import (  # noqa: E402
    CellState,
    CellType,
    DrugParams,
    PDConfig,
    PDEConfig,
    SimConfig,
)
from tmesim.sim.engine import simulate_response  # noqa: E402
from tmesim.sim.mechanism import (  # noqa: E402
    ALREADY_DEAD,
    KILLED,
    LABELS,
    NOT_TUMOR,
    SURVIVED_RESISTANT,
    SURVIVED_STOCHASTIC,
    SURVIVED_UNDERDOSED,
    SURVIVED_UNREACHED,
    UNATTRIBUTED,
    MechanismReconciliationError,
    MechanismThresholds,
    attribute,
    classify,
    reconcile,
)

PDE = PDEConfig(dt=10.0, total_time=1800.0, convection=True, save_every=100)


# --------------------------------------------------------------------------- #
# 1. Faithfulness on a real reconstruction
# --------------------------------------------------------------------------- #
def test_attribution_reconciles_bit_exactly_with_kill_fraction(
        real_micro, real_drug):
    """The KILLED count must rebuild kill_fraction bit-for-bit, not nearly."""
    cfg = SimConfig(drug=real_drug, pde=PDE, seed=0)
    res = simulate_response(real_micro, cfg, attribute=True)
    rec = res.mechanism["reconciliation"]

    assert rec["bit_exact"] is True, rec
    # The load-bearing comparison, restated here so the test fails on the claim
    # itself and not only on the flag the module set.
    assert rec["kill_fraction_from_labels"] == res.kill_fraction
    assert rec["abs_difference"] == 0.0

    labels = res.mechanism["labels"]
    assert labels.dtype == np.uint8
    assert labels.size == real_micro.cells.n
    # Every cell is labelled, and every tumour cell gets a tumour label.
    tumor = np.asarray(real_micro.cells.cell_type) == int(CellType.TUMOR)
    assert np.all(labels[~tumor] == NOT_TUMOR)
    assert np.all(labels[tumor] != NOT_TUMOR)


def test_attribution_is_bit_identical_to_a_run_without_it(real_micro, real_drug):
    """Asking WHY must not change WHAT. Both the kill and the exposure field."""
    base = simulate_response(real_micro, SimConfig(drug=real_drug, pde=PDE, seed=0))
    with_attr = simulate_response(
        real_micro, SimConfig(drug=real_drug, pde=PDE, seed=0), attribute=True)

    assert with_attr.kill_fraction == base.kill_fraction
    assert with_attr.viable_tumor_fraction == base.viable_tumor_fraction
    assert np.array_equal(with_attr.transport.auc, base.transport.auc)
    assert np.array_equal(with_attr.viable_history, base.viable_history)
    assert base.mechanism is None       # off by default
    assert with_attr.mechanism is not None


def test_survival_breakdown_sums_to_one_and_counts_add_up(real_micro, real_drug):
    m = simulate_response(
        real_micro, SimConfig(drug=real_drug, pde=PDE, seed=0),
        attribute=True).mechanism
    counts = m["counts"]

    assert sum(counts.values()) == m["n_cells"]
    assert (counts["KILLED"] + counts["ALREADY_DEAD"]
            + m["n_survived"]) == m["n_tumor"]
    if m["n_survived"]:
        assert abs(sum(m["survival_breakdown"].values()) - 1.0) < 1e-12
    assert set(m["label_names"]) == set(LABELS)


# --------------------------------------------------------------------------- #
# 2. The bands themselves
# --------------------------------------------------------------------------- #
def _diag(ratio, modulation, p_kill, alive=None, killed=None):
    """A hand-built diagnostics dict with the five arrays classify() reads."""
    n = len(ratio)
    return {
        "ratio": np.asarray(ratio, dtype=float),
        "modulation": np.asarray(modulation, dtype=float),
        "p_kill": np.asarray(p_kill, dtype=float),
        "alive": np.ones(n, bool) if alive is None else np.asarray(alive, bool),
        "killed": np.zeros(n, bool) if killed is None else np.asarray(killed, bool),
        "hill": 1.0,
    }


def test_bands_place_each_cell_in_the_documented_category():
    th = MechanismThresholds().resolve(PDConfig(), hill=1.0)
    assert th.resistant_modulation == 0.5          # max(0.35, 0.5), derived
    assert th.stochastic_min_p_kill == 0.25        # f(1.0) * 0.5, derived

    #            unreached  underdosed  resistant  stochastic  killed  pre-dead
    ratio =      [0.01,     0.5,        4.0,       4.0,        4.0,    4.0]
    modulation = [1.0,      1.0,        0.35,      1.0,        1.0,    1.0]
    p_kill =     [0.01,     0.33,       0.28,      0.80,       0.80,   0.0]
    diag = _diag(ratio, modulation, p_kill,
                 alive=[1, 1, 1, 1, 1, 0], killed=[0, 0, 0, 0, 1, 0])
    ctype = np.full(6, int(CellType.TUMOR))

    got = classify(ctype, diag, th)
    assert list(got) == [SURVIVED_UNREACHED, SURVIVED_UNDERDOSED,
                         SURVIVED_RESISTANT, SURVIVED_STOCHASTIC,
                         KILLED, ALREADY_DEAD]


def test_exposure_binds_before_resistance_and_the_overlap_is_reported():
    """A cell short of drug AND resistant is labelled by its exposure band —
    and the count of such cells is surfaced rather than buried."""
    th = MechanismThresholds()
    diag = _diag(ratio=[0.01, 0.5, 4.0], modulation=[0.35, 0.35, 0.35],
                 p_kill=[0.003, 0.11, 0.28])
    ctype = np.full(3, int(CellType.TUMOR))
    out = attribute(ctype, diag, kill_fraction=0.0, pd=PDConfig(), thresholds=th)

    assert list(out["labels"]) == [SURVIVED_UNREACHED, SURVIVED_UNDERDOSED,
                                   SURVIVED_RESISTANT]
    assert out["joint_limited"]["count"] == 2
    assert out["joint_limited"]["fraction_of_exposure_limited"] == 1.0


def test_unexplained_survival_is_unattributed_not_reassigned():
    """A ceiling the bands do not model must NOT be called bad luck.

    A cell at exactly the adequate boundary with no resistance would have
    ``p_kill = 0.5`` under the production ceiling — honest chance. Drop
    ``max_kill_rate`` to 0.4 and the same cell's kill probability falls to 0.2,
    below the derived CHANCE floor of 0.25, while exposure and biology both
    still read as adequate. Calling that "bad luck" would be a lie: it had a one
    in five chance. The ceiling is a suppression the four bands do not describe,
    so the cell must land in UNATTRIBUTED and be counted there.
    """
    pd_cfg = PDConfig(max_kill_rate=0.4)
    th = MechanismThresholds().resolve(pd_cfg, hill=1.0)
    assert th.stochastic_min_p_kill == 0.25
    diag = _diag(ratio=[1.0], modulation=[1.0], p_kill=[0.4 * 0.5])
    out = attribute(np.array([int(CellType.TUMOR)]), diag,
                    kill_fraction=0.0, pd=pd_cfg)

    assert list(out["labels"]) == [UNATTRIBUTED]
    assert out["counts"]["UNATTRIBUTED"] == 1
    assert out["counts"]["SURVIVED_STOCHASTIC"] == 0
    assert out["survival_breakdown"]["unattributed"] == 1.0
    # With the production ceiling the identical cell is honest chance, so the
    # category is not a catch-all that swallows the adequate-exposure branch.
    assert list(attribute(np.array([int(CellType.TUMOR)]),
                          _diag([1.0], [1.0], [0.5]), 0.0,
                          pd=PDConfig())["labels"]) == [SURVIVED_STOCHASTIC]


def test_non_finite_exposure_is_unattributed():
    out = attribute(np.array([int(CellType.TUMOR)]),
                    _diag([np.nan], [1.0], [0.8]), 0.0, pd=PDConfig())
    assert list(out["labels"]) == [UNATTRIBUTED]


def test_thresholds_follow_the_pdconfig_the_run_used():
    """resistant_modulation is derived, so a counterfactual resistance-off run
    moves the band with it instead of judging against a stale constant."""
    off = MechanismThresholds().resolve(
        PDConfig(hypoxia_resistance=1.0, quiescent_resistance=1.0), hill=1.2)
    assert off.resistant_modulation == 1.0
    lenient = MechanismThresholds().resolve(
        PDConfig(hypoxia_resistance=0.1, quiescent_resistance=0.2), hill=1.0)
    assert lenient.resistant_modulation == 0.2
    assert lenient.stochastic_min_p_kill == pytest.approx(0.5 * 0.2)


def test_a_lying_attribution_raises_instead_of_returning():
    diag = _diag(ratio=[4.0], modulation=[1.0], p_kill=[0.8])
    with pytest.raises(MechanismReconciliationError):
        attribute(np.array([int(CellType.TUMOR)]), diag,
                  kill_fraction=0.5, pd=PDConfig())   # nothing was killed


def test_reconcile_mirrors_viable_tumor_fraction_arithmetic():
    labels = np.array([KILLED] * 3 + [SURVIVED_STOCHASTIC] * 7
                      + [NOT_TUMOR] * 5, dtype=np.uint8)
    rec = reconcile(labels, 1.0 - 7 / 10)
    assert rec["bit_exact"] is True
    assert rec["n_tumor"] == 10 and rec["n_killed"] == 3


# --------------------------------------------------------------------------- #
# 3. pkpd's diagnostics hook is inert
# --------------------------------------------------------------------------- #
def test_diagnostics_hook_changes_nothing_about_the_kill(real_micro, real_drug):
    """apply_drug_kill must behave identically with and without diagnostics."""
    from tmesim.core import Microenvironment, TransportResult

    def run(diag):
        micro = Microenvironment(
            domain=real_micro.domain, fields=real_micro.fields,
            cells=real_micro.cells.copy(), patient_id=real_micro.patient_id)
        auc = np.linspace(0.0, 4.0 * real_drug.ic50 * 600.0,
                          micro.domain.nx * micro.domain.ny)
        transport = TransportResult(
            auc=auc.reshape(micro.domain.shape),
            final_concentration=micro.domain.zeros(),
            times=np.array([0.0, 600.0]),
            mean_concentration_history=np.array([0.0, 0.0]))
        p = apply_drug_kill(micro, transport, PDConfig(), real_drug,
                            rng=np.random.default_rng(7), diagnostics=diag)
        return p, micro.cells.state.copy()

    p_plain, state_plain = run(None)
    diag = {}
    p_diag, state_diag = run(diag)

    assert np.array_equal(p_plain, p_diag)
    assert np.array_equal(state_plain, state_diag)
    # and the dict really carries the arrays the kill decided on
    for key in ("E", "ratio", "p_hill", "modulation", "p_kill", "draws",
                "alive", "killed"):
        assert key in diag, key
    assert np.array_equal(diag["p_kill"], p_diag)
    assert np.array_equal(
        diag["killed"], diag["alive"] & (diag["draws"] < diag["p_kill"]))


# --------------------------------------------------------------------------- #
# 4. Time-lapse
# --------------------------------------------------------------------------- #
def test_frame_capture_does_not_perturb_the_solve(real_micro, real_drug):
    base = simulate_response(real_micro, SimConfig(drug=real_drug, pde=PDE, seed=0))
    filmed = simulate_response(real_micro, SimConfig(drug=real_drug, pde=PDE, seed=0),
                               capture_frames=8)

    assert filmed.kill_fraction == base.kill_fraction
    assert np.array_equal(filmed.transport.auc, base.transport.auc)
    assert np.array_equal(filmed.transport.final_concentration,
                          base.transport.final_concentration)
    assert base.transport.frames is None


def test_frames_span_the_window_and_start_from_the_zero_initial_condition(
        real_micro, real_drug):
    res = simulate_response(real_micro, SimConfig(drug=real_drug, pde=PDE, seed=0),
                            capture_frames=6)
    frames, times = res.transport.frames, res.transport.frame_times

    assert frames.shape[0] == 6 == times.size
    assert frames.shape[1:] == real_micro.domain.shape
    assert times[0] == 0.0
    assert times[-1] == pytest.approx(PDE.total_time)
    assert np.all(np.diff(times) > 0)
    # t=0 is the initial condition: no drug anywhere.
    assert np.all(frames[0] == 0.0)
    # the last captured frame IS the solver's final concentration
    assert np.array_equal(frames[-1], res.transport.final_concentration)
    # drug enters: total tissue drug is non-decreasing early on
    assert frames[-1].sum() > frames[1].sum()


def test_capture_frames_is_clipped_to_the_number_of_steps(real_micro, real_drug):
    """Asking for more frames than there are steps returns every step, not
    duplicates of the last one."""
    short = PDEConfig(dt=10.0, total_time=40.0, convection=True, save_every=1)
    res = simulate_response(real_micro, SimConfig(drug=real_drug, pde=short, seed=0),
                            capture_frames=30)
    assert res.transport.frames.shape[0] == 5     # steps 0..4
    assert np.unique(res.transport.frame_times).size == 5


def test_capture_frames_below_two_is_refused():
    from tmesim.physics.diffusion import simulate_drug_transport

    with pytest.raises(ValueError, match="at least a first and a last frame"):
        simulate_drug_transport(
            _tiny_micro(), DrugParams(name="x"), _tiny_micro().domain.zeros(),
            PDEConfig(dt=1.0, total_time=10.0, convection=False, save_every=1),
            capture_frames=1)


def _tiny_micro():
    from tmesim.core import CellPopulation, Domain, Microenvironment, TissueFields

    domain = Domain(nx=4, ny=4, dx=10.0)
    return Microenvironment(
        domain=domain,
        fields=TissueFields(ecm_density=domain.zeros(),
                            vascular_density=domain.zeros()),
        cells=CellPopulation(
            x=np.array([5.0]), y=np.array([5.0]),
            cell_type=np.array([int(CellType.TUMOR)]),
            state=np.array([int(CellState.VIABLE)]),
            sensitivity=np.array([1.0])),
    )


# --------------------------------------------------------------------------- #
# 5. The wire: TMSM / TMSL payloads and the two routes
# --------------------------------------------------------------------------- #
fastapi = pytest.importorskip("fastapi", reason="the app/server backend needs fastapi")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi.testclient import TestClient  # noqa: E402

from app.server.api import create_app, validate_spec  # noqa: E402
from app.server.jobs import JobManager  # noqa: E402
from app.server.payloads import (  # noqa: E402
    MAGIC_MECHANISM,
    MAGIC_TIMELAPSE,
    MAX_TIMELAPSE_FRAMES,
    decode_mechanism,
    decode_timelapse,
    encode_field,
    encode_mechanism,
    encode_timelapse,
)

_SECTION = "xenium_319"
_SPEC = {
    "section": _SECTION, "dx_um": 30.0, "drug": "temozolomide",
    "sensitivity_mode": "uniform", "vascular_mode": "discrete",
    "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": None,
    "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
           "n_doses": 1, "dose_interval_s": 86400.0},
    "total_time_s": 600.0, "timelapse_frames": 4,
}


def test_tmsm_round_trips_and_refuses_an_out_of_range_code():
    labels = np.array([0, 1, 2, 3, 4, 5, 6, 7], dtype=np.uint8)
    agg = {"label_names": list(LABELS), "counts": {"KILLED": 1}}
    blob = encode_mechanism(labels, agg)
    assert blob[:4] == MAGIC_MECHANISM

    back = decode_mechanism(blob)
    assert back["n_cells"] == 8
    assert back["n_labels"] == len(LABELS)
    assert np.array_equal(back["labels"], labels)
    assert back["aggregate"] == agg
    # A truncated payload must raise, not render garbage.
    with pytest.raises(ValueError, match="expected"):
        decode_mechanism(blob[:-1])
    with pytest.raises(ValueError, match="uint8"):
        encode_mechanism(np.array([300]), agg)


def test_tmsl_frames_are_byte_identical_to_the_tmsf_field_payload():
    """The container extends the existing field payload; it does not replace it.

    That is the whole reason the front end can decode a frame with the TMSF
    reader it already ships, so it is asserted at the byte level.
    """
    rng = np.random.default_rng(0)
    frames = rng.random((3, 5, 7))
    blob = encode_timelapse(frames, [0.0, 10.0, 20.0], 12.0, interval_s=10.0)
    assert blob[:4] == MAGIC_TIMELAPSE

    # header 20 B, then per frame: 8 B prefix + a whole TMSF block
    solo = encode_field(frames[1], 12.0)
    assert solo in blob

    back = decode_timelapse(blob)
    assert back["n_frames"] == 3
    assert back["interval_s"] == 10.0 and back["window_s"] == 20.0
    assert [f["frame"] for f in back["frames"]] == [0, 1, 2]
    assert [f["time_s"] for f in back["frames"]] == [0.0, 10.0, 20.0]
    assert np.allclose(back["frames"][2]["values"], frames[2], atol=1e-6)

    with pytest.raises(ValueError, match="MAX_TIMELAPSE_FRAMES"):
        encode_timelapse(np.zeros((MAX_TIMELAPSE_FRAMES + 1, 2, 2)),
                         np.arange(MAX_TIMELAPSE_FRAMES + 1), 1.0)


def test_validate_accepts_the_new_fields_and_warns_about_the_static_plateau():
    ok = validate_spec(dict(_SPEC))
    assert ok["ok"] is True
    assert any("never CLEAR" in w for w in ok["warnings"])

    bad = validate_spec(dict(_SPEC, timelapse_frames=1))
    assert bad["ok"] is False
    assert any("at least a first and a last frame" in e for e in bad["errors"])

    tr = validate_spec(dict(_SPEC, time_resolved=True))
    assert tr["ok"] is True
    assert any("SCHEDULE effect" in w for w in tr["warnings"])
    assert not any("never CLEAR" in w for w in tr["warnings"])


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    jobs = JobManager(root=str(tmp_path_factory.mktemp("mech_runs")))
    with TestClient(create_app(jobs)) as c:
        yield c
    jobs.shutdown()


def _await(client, run_id, timeout=600.0):
    import time as _t
    t0 = _t.time()
    while _t.time() - t0 < timeout:
        st = client.get(f"/api/runs/{run_id}/status").json()
        if st["status"] in ("done", "failed", "cancelled"):
            return st
        _t.sleep(0.4)
    raise AssertionError(f"run {run_id} did not finish in {timeout}s")


def test_routes_serve_a_reconciling_attribution_and_a_playable_timelapse(client):
    run_id = client.post("/api/runs", json=_SPEC).json()["run_id"]
    st = _await(client, run_id)
    assert st["status"] == "done", client.get(
        f"/api/runs/{run_id}/traceback").json()

    run = client.get(f"/api/runs/{run_id}").json()

    # --- mechanism ------------------------------------------------------- #
    r = client.get(f"/api/runs/{run_id}/mechanism")
    assert r.status_code == 200
    assert r.headers["content-type"] == "application/octet-stream"
    mech = decode_mechanism(r.content)
    agg = mech["aggregate"]

    geom = client.get(f"/api/runs/{run_id}/geometry")
    n_cells = int(np.frombuffer(geom.content, dtype="<u4", count=1, offset=8)[0])
    assert mech["n_cells"] == n_cells, "labels are not aligned to the geometry"

    # THE claim: the labels reproduce the engine's own kill_fraction exactly.
    assert agg["reconciliation"]["bit_exact"] is True
    assert (agg["reconciliation"]["kill_fraction_from_engine"]
            == run["metrics"]["kill_fraction"])
    assert agg["counts"]["KILLED"] + agg["counts"]["ALREADY_DEAD"] \
        + agg["n_survived"] == agg["n_tumor"]
    assert run["mechanism"]["counts"] == agg["counts"]
    assert "labels" not in run["mechanism"], "the JSON must not carry 850k labels"

    # --- time-lapse ------------------------------------------------------- #
    r = client.get(f"/api/runs/{run_id}/timelapse")
    assert r.status_code == 200
    tl = decode_timelapse(r.content)
    assert tl["n_frames"] == 4
    assert [f["frame"] for f in tl["frames"]] == [0, 1, 2, 3]
    assert tl["frames"][0]["time_s"] == 0.0
    assert tl["frames"][-1]["time_s"] == pytest.approx(_SPEC["total_time_s"])
    assert np.all(tl["frames"][0]["values"] == 0.0)      # zero initial condition
    assert tl["frames"][-1]["values"].max() > 0.0        # drug arrived
    assert run["timelapse"]["available"] is True
    assert run["timelapse"]["n_frames"] == 4
    # the static-source caveat must travel with the frames
    assert "quasi-steady state" in run["timelapse"]["plateau_note"]


def test_a_run_without_them_says_why_rather_than_404ing_blankly(client):
    spec = dict(_SPEC, drug=None, timelapse_frames=0)
    run_id = client.post("/api/runs", json=spec).json()["run_id"]
    assert _await(client, run_id)["status"] == "done"

    r = client.get(f"/api/runs/{run_id}/mechanism")
    assert r.status_code == 404
    assert "reconstruction-only" in r.json()["detail"]["reason"]

    r = client.get(f"/api/runs/{run_id}/timelapse")
    assert r.status_code == 404
    assert "reconstruction-only" in r.json()["detail"]["reason"]

    run = client.get(f"/api/runs/{run_id}").json()
    assert run["mechanism"] is None
    assert run["timelapse"]["available"] is False
