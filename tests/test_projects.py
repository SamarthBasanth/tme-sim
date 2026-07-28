"""Tests for the multi-patient layer (``app/server/projects.py``).

Three things, in rising order of what they would cost if they broke:

1. **Persistence** — a project is a JSON document under ``app/projects/`` and
   must survive a process restart with its specimens, its capability blocks and
   its run associations intact.
2. **Association** — a specimen holds run IDS into ``app/runs/``. Runs are never
   copied, never deleted by a project operation, and a run computed on ANOTHER
   section can never be filed under a specimen (that would silently corrupt every
   comparison built on the project).
3. **The comparison's capability-aware flagging** — the part that carries the
   scientific claim. Rank agreement is the primary result; an absolute-kill
   delta is refused unless BOTH panels can support a measured cell state AND both
   runs used one. This is ``results/cell_state_impact.md`` §5 (a panel with no
   hypoxia-response transcript cannot support an absolute kill under any route,
   because kill moves 0.3252 absolute / 1.77x with the route alone) applied to
   the two-patient case ``results/heldout_patient.md`` §3.3/§3.5 measured: the
   ranking replicated, the absolute kills did not.

The comparison tests build their run payloads by hand rather than simulating, so
the FLAGGING LOGIC is tested exhaustively and in milliseconds. The end-to-end
path (real sections, real runs, real ranking) is exercised by
``test_compare_on_two_real_sections``, which is skipped unless both sections are
cached.
"""
from __future__ import annotations

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
from app.server.api import create_app  # noqa: E402
from app.server.jobs import JobManager  # noqa: E402
from app.server.projects import (  # noqa: E402
    ProjectError,
    ProjectStore,
    capability_report,
    compare_project,
    spec_signature,
)

SECTION_A = "xenium_319"      # 319-gene brain panel: NO hypoxia-response gene
SECTION_B = "xenium_io"       # 480-gene IO panel: carries VEGFA


def _available(section_id: str) -> bool:
    return any(s["id"] == section_id and s["available"]
               for s in cache_mod.list_sections())


requires_both = pytest.mark.skipif(
    not (_available(SECTION_A) and _available(SECTION_B)),
    reason=f"sections {SECTION_A!r} and {SECTION_B!r} are not both cached",
)


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture
def store(tmp_path):
    return ProjectStore(root=str(tmp_path / "projects"))


@pytest.fixture
def client(tmp_path):
    """A TestClient whose runs AND projects are isolated from the real dirs."""
    manager = JobManager(root=str(tmp_path / "runs"), start_worker=False)
    app = create_app(manager, ProjectStore(root=str(tmp_path / "projects")))
    with TestClient(app) as c:
        c._manager = manager  # type: ignore[attr-defined]
        yield c
    manager.shutdown()


def _fake_run(run_id, section, drug, kill, *, dx=30.0, state_mode="oxygen",
              state_measured=False, created_at=1000.0, status="done",
              report=None, **metrics):
    """A RunResult view shaped exactly like ``Job.result_view()`` returns one."""
    return {
        "run_id": run_id,
        "status": status,
        "created_at": created_at,
        "spec": {
            "section": section, "dx_um": dx, "drug": drug,
            "sensitivity_mode": "uniform", "vascular_mode": "discrete",
            "oxygen_mode": "pde", "state_mode": state_mode, "mgmt": None,
            "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
                   "n_doses": 1, "dose_interval_s": 86400.0},
            "total_time_s": 600.0,
        },
        "metrics": {
            "kill_fraction": kill,
            "mean_tumor_auc": metrics.get("mean_tumor_auc", 1.0e6),
            "therapeutic_coverage_1x": metrics.get("therapeutic_coverage_1x", 0.5),
            "penetration_length_um": metrics.get("penetration_length_um", 200.0),
        },
        "report": report if report is not None else {
            "n_cells": 40887,
            "celltype_fractions": {"TUMOR": 0.738, "ENDOTHELIAL": 0.0974},
            "celltype_counts": {"TUMOR": 30174, "ENDOTHELIAL": 3984},
            "extent_um": [8100.0, 5088.0],
            "tissue_coverage_fraction": 0.677,
            "mean_vascular_density": 0.025,
            "mean_ecm_density": 0.0759,
            "mean_oxygen": 0.2955,
            "ecm_markers_used": ["COL1A2", "POSTN"],
            "cell_state": {
                "mode": state_mode,
                "source": "measured markers" if state_measured else "simulated O2",
                "is_measured": state_measured,
                "unmeasured_states": [] if state_measured else None,
                "fractions_tumor": {"VIABLE": 0.1, "QUIESCENT": 0.9,
                                    "HYPOXIC": 0.0},
            },
        },
        "provenance": {"state_is_measured": state_measured,
                       "state_mode": state_mode},
    }


def _capability(section, *, hypoxia, panel="panel-X"):
    """A capability block with only the fields the comparison reads."""
    return {
        "section": section, "panel": panel, "genes_known": True,
        "hypoxia_present": ["VEGFA"] if hypoxia else [],
        "hypoxia_missing": [] if hypoxia else ["VEGFA", "HIF1A", "CA9"],
        "supports_cell_state": bool(hypoxia),
        "absolute_kill_supported": bool(hypoxia),
        "absolute_kill_reason": (
            "carries VEGFA" if hypoxia else
            "carries none of the hypoxia-response transcripts, so NO state route "
            "clears the bar (results/cell_state_impact.md §5)"),
        "state_routes": {
            "oxygen": {"supported": True, "is_measured": False, "reason": ""},
            "vessel_distance": {"supported": True, "is_measured": False,
                                "reason": ""},
            "markers": {"supported": bool(hypoxia), "is_measured": bool(hypoxia),
                        "reason": ""},
        },
        "available": True,
    }


def _project(runs_a, runs_b, *, hypoxia_a=False, hypoxia_b=True):
    """A two-specimen project plus a resolver over the supplied run views."""
    by_id = {r["run_id"]: r for r in (runs_a + runs_b)}
    project = {
        "id": "p_00001", "name": "T",
        "specimens": [
            {"id": "s_00001", "label": "A", "section": SECTION_A,
             "run_ids": [r["run_id"] for r in runs_a],
             "capability": _capability(SECTION_A, hypoxia=hypoxia_a,
                                       panel="brain_319")},
            {"id": "s_00002", "label": "B", "section": SECTION_B,
             "run_ids": [r["run_id"] for r in runs_b],
             "capability": _capability(SECTION_B, hypoxia=hypoxia_b,
                                       panel="io_480")},
        ],
    }
    return project, (lambda rid: by_id.get(rid))


# --------------------------------------------------------------------------- #
# 1. Persistence
# --------------------------------------------------------------------------- #
class TestPersistence:
    def test_project_survives_a_new_store_instance(self, store, tmp_path):
        p = store.create("Study 1", notes="two patients")
        assert p["id"] == "p_00001"

        reopened = ProjectStore(root=store.root)
        got = reopened.get(p["id"], refresh_capability=False)
        assert got["name"] == "Study 1"
        assert got["notes"] == "two patients"
        assert got["specimens"] == []

        # The document is on disk where the contract says it is.
        assert os.path.exists(
            os.path.join(store.root, p["id"], "project.json"))

    def test_ids_do_not_collide_after_a_restart(self, store):
        store.create("one")
        store.create("two")
        reopened = ProjectStore(root=store.root)
        third = reopened.create("three")
        assert third["id"] == "p_00003"
        assert {p["id"] for p in reopened.list()} == {"p_00001", "p_00002",
                                                      "p_00003"}

    @requires_both
    def test_specimen_and_capability_survive_a_restart(self, store):
        p = store.create("Study")
        s = store.add_specimen(p["id"], SECTION_A, label="Patient A")
        assert s["id"] == "s_00001"
        assert s["capability"]["section"] == SECTION_A

        reopened = ProjectStore(root=store.root)
        got = reopened.get(p["id"], refresh_capability=False)
        (spec,) = got["specimens"]
        assert spec["label"] == "Patient A"
        assert spec["section"] == SECTION_A
        assert spec["capability"]["absolute_kill_supported"] is False

    def test_update_and_delete_round_trip(self, store):
        p = store.create("old name")
        store.update(p["id"], name="new name", notes="n")
        assert ProjectStore(root=store.root).get(p["id"])["name"] == "new name"
        store.delete(p["id"])
        with pytest.raises(ProjectError) as exc:
            ProjectStore(root=store.root).get(p["id"])
        assert exc.value.status == 404

    def test_unknown_project_is_404_and_empty_name_is_400(self, store):
        with pytest.raises(ProjectError) as exc:
            store.get("p_99999")
        assert exc.value.status == 404
        with pytest.raises(ProjectError) as exc:
            store.create("   ")
        assert exc.value.status == 400


# --------------------------------------------------------------------------- #
# 2. Association
# --------------------------------------------------------------------------- #
@requires_both
class TestAssociation:
    def test_runs_are_associated_not_duplicated(self, store):
        p = store.create("S")
        s = store.add_specimen(p["id"], SECTION_A)
        run = _fake_run("r_00001", SECTION_A, "temozolomide", 0.1)

        out = store.attach_run(p["id"], s["id"], "r_00001", run)
        assert out["attached"] is True

        doc = ProjectStore(root=store.root).get(p["id"], refresh_capability=False)
        (spec,) = doc["specimens"]
        # Only the ID is persisted — no copy of the result body anywhere.
        assert spec["run_ids"] == ["r_00001"]
        blob = open(os.path.join(store.root, p["id"], "project.json")).read()
        assert "kill_fraction" not in blob

    def test_attaching_twice_is_idempotent(self, store):
        p = store.create("S")
        s = store.add_specimen(p["id"], SECTION_A)
        run = _fake_run("r_00001", SECTION_A, "temozolomide", 0.1)
        store.attach_run(p["id"], s["id"], "r_00001", run)
        again = store.attach_run(p["id"], s["id"], "r_00001", run)
        assert again["attached"] is False
        assert again["n_runs"] == 1

    def test_a_run_from_another_section_is_refused(self, store):
        """The one association error that would silently corrupt a comparison."""
        p = store.create("S")
        s = store.add_specimen(p["id"], SECTION_A)
        foreign = _fake_run("r_00009", SECTION_B, "cisplatin", 0.3)
        with pytest.raises(ProjectError) as exc:
            store.attach_run(p["id"], s["id"], "r_00009", foreign)
        assert exc.value.status == 409
        assert SECTION_B in exc.value.message and SECTION_A in exc.value.message
        assert store.get(p["id"])["specimens"][0]["run_ids"] == []

    def test_attaching_a_nonexistent_run_is_404(self, store):
        p = store.create("S")
        s = store.add_specimen(p["id"], SECTION_A)
        with pytest.raises(ProjectError) as exc:
            store.attach_run(p["id"], s["id"], "r_99999", None)
        assert exc.value.status == 404

    def test_detach_and_delete_retain_the_run(self, store):
        p = store.create("S")
        s = store.add_specimen(p["id"], SECTION_A)
        run = _fake_run("r_00001", SECTION_A, "temozolomide", 0.1)
        store.attach_run(p["id"], s["id"], "r_00001", run)

        out = store.detach_run(p["id"], s["id"], "r_00001")
        assert out["detached"] is True
        assert store.get(p["id"])["specimens"][0]["run_ids"] == []
        with pytest.raises(ProjectError) as exc:
            store.detach_run(p["id"], s["id"], "r_00001")
        assert exc.value.status == 404

        store.attach_run(p["id"], s["id"], "r_00001", run)
        removed = store.remove_specimen(p["id"], s["id"])
        assert removed["runs_retained"] == ["r_00001"]

        p2 = store.create("S2")
        s2 = store.add_specimen(p2["id"], SECTION_A)
        store.attach_run(p2["id"], s2["id"], "r_00001", run)
        deleted = store.delete(p2["id"])
        assert deleted["runs_retained"] == ["r_00001"]

    def test_one_specimen_per_section(self, store):
        p = store.create("S")
        store.add_specimen(p["id"], SECTION_A)
        with pytest.raises(ProjectError) as exc:
            store.add_specimen(p["id"], SECTION_A)
        assert exc.value.status == 409

    def test_unknown_section_is_rejected(self, store):
        p = store.create("S")
        with pytest.raises(ProjectError) as exc:
            store.add_specimen(p["id"], "not_a_section")
        assert exc.value.status == 400


# --------------------------------------------------------------------------- #
# 3. Capability
# --------------------------------------------------------------------------- #
@requires_both
class TestCapability:
    def test_the_319_panel_cannot_support_an_absolute_kill(self):
        """results/cell_state_impact.md §5, read off the panel's own gene list."""
        cap = capability_report(SECTION_A)
        assert cap["genes_known"] is True
        assert cap["hypoxia_present"] == []
        assert cap["absolute_kill_supported"] is False
        assert cap["supports_cell_state"] is False
        assert "hypoxia" in cap["absolute_kill_reason"].lower()
        assert "0.3252" in cap["absolute_kill_reason"]

    def test_the_io_panel_can(self):
        cap = capability_report(SECTION_B)
        assert "VEGFA" in cap["hypoxia_present"]
        assert cap["absolute_kill_supported"] is True
        assert cap["supports_cell_state"] is True
        assert cap["state_routes"]["markers"]["supported"] is True

    def test_the_per_cell_detection_rate_is_never_claimed_from_the_panel(self):
        """A gene being ON a panel says nothing about how often it fires — the
        distinction that disqualifies the 319 panel's proliferation call."""
        cap = capability_report(SECTION_A)
        assert "MKI67" in cap["proliferation_present"]
        assert cap["state_routes"]["markers"]["per_cell_call_knowable"] is False

    def test_unknown_section_raises(self):
        with pytest.raises(ProjectError):
            capability_report("nope")


# --------------------------------------------------------------------------- #
# 4. Matching
# --------------------------------------------------------------------------- #
class TestMatching:
    def test_signature_ignores_section_and_drug_only(self):
        a = _fake_run("r1", SECTION_A, "temozolomide", 0.1)["spec"]
        b = _fake_run("r2", SECTION_B, "cisplatin", 0.4)["spec"]
        assert spec_signature(a) == spec_signature(b)

    def test_signature_separates_a_different_mesh_or_state_route(self):
        a = _fake_run("r1", SECTION_A, "temozolomide", 0.1)["spec"]
        assert spec_signature(a) != spec_signature(
            _fake_run("r2", SECTION_A, "temozolomide", 0.1, dx=12.0)["spec"])
        assert spec_signature(a) != spec_signature(
            _fake_run("r3", SECTION_A, "temozolomide", 0.1,
                      state_mode="markers")["spec"])

    def test_the_counterfactual_arm_is_never_matched_against_a_real_one(self):
        """resistance_off=True switches this tumour's resistance OFF; that arm
        is a what-if, not a prediction, and must not be differenced against a
        normal arm of the same drug on another patient."""
        a = _fake_run("r1", SECTION_A, "temozolomide", 0.1)["spec"]
        counterfactual = dict(a, resistance_off=True)
        assert spec_signature(a) != spec_signature(counterfactual)
        # ...but "flag off" and "flag not sent" are the SAME experiment.
        assert spec_signature(a) == spec_signature(dict(a, resistance_off=False))

    def test_an_unknown_spec_field_fails_closed(self):
        a = _fake_run("r1", SECTION_A, "temozolomide", 0.1)["spec"]
        assert spec_signature(a) != spec_signature(dict(a, future_knob="x"))
        assert spec_signature(a) == spec_signature(dict(a, future_knob=None))
        assert spec_signature(a) == spec_signature(
            dict(a, pk=dict(a["pk"], future_pk_knob=False)))
        assert spec_signature(a) != spec_signature(
            dict(a, pk=dict(a["pk"], future_pk_knob=3.0)))

    def test_runs_at_a_different_mesh_are_not_matched(self):
        runs_a = [_fake_run("a1", SECTION_A, "cisplatin", 0.37),
                  _fake_run("a2", SECTION_A, "temozolomide", 0.01)]
        runs_b = [_fake_run("b1", SECTION_B, "cisplatin", 0.31, dx=12.0),
                  _fake_run("b2", SECTION_B, "temozolomide", 0.02, dx=12.0)]
        project, resolve = _project(runs_a, runs_b)
        doc = compare_project(project, resolve)
        assert doc["ranking"]["agrees"] is None
        assert "NOT EVALUABLE" in doc["ranking"]["verdict"]
        assert doc["kill"] == []


# --------------------------------------------------------------------------- #
# 5. The comparison — ranking is the primary result
# --------------------------------------------------------------------------- #
class TestRanking:
    #: The ordering results/heldout_patient.md §3.3 found identical on both
    #: sections, restricted to the three drugs this study runs.
    _A = [("cisplatin", 0.3699), ("doxorubicin", 0.2165), ("temozolomide", 0.0147)]
    _B = [("cisplatin", 0.3162), ("doxorubicin", 0.2284), ("temozolomide", 0.0274)]

    def _doc(self, a=None, b=None, **kw):
        a = a if a is not None else self._A
        b = b if b is not None else self._B
        runs_a = [_fake_run(f"a{i}", SECTION_A, d, k)
                  for i, (d, k) in enumerate(a)]
        runs_b = [_fake_run(f"b{i}", SECTION_B, d, k)
                  for i, (d, k) in enumerate(b)]
        project, resolve = _project(runs_a, runs_b, **kw)
        return compare_project(project, resolve)

    def test_ranking_agreement_is_detected_and_is_the_headline(self):
        doc = self._doc()
        rk = doc["ranking"]
        assert rk["agrees"] is True
        assert rk["drugs"] == ["cisplatin", "doxorubicin", "temozolomide"]
        assert (rk["by_specimen"]["s_00001"]["order"]
                == rk["by_specimen"]["s_00002"]["order"]
                == ["cisplatin", "doxorubicin", "temozolomide"])
        assert rk["pairwise"][0]["kendall_tau_b"] == 1.0
        assert rk["pairwise"][0]["swapped"] == []
        assert doc["headline"].startswith("RANKING AGREES")
        assert doc["reading_order"][1].startswith("ranking")

    def test_a_swap_is_detected_and_named(self):
        swapped = [("doxorubicin", 0.40), ("cisplatin", 0.30),
                   ("temozolomide", 0.03)]
        rk = self._doc(b=swapped)["ranking"]
        assert rk["agrees"] is False
        pw = rk["pairwise"][0]
        assert pw["agrees"] is False
        assert sorted(pw["swapped"][0]) == ["cisplatin", "doxorubicin"]
        assert pw["discordant_pairs"] == 1
        assert "DISAGREES" in rk["verdict"]

    def test_ranking_uses_only_the_drugs_both_specimens_ran(self):
        rk = self._doc(b=self._B + [("lomustine", 0.22)])["ranking"]
        assert "lomustine" not in rk["drugs"]
        assert rk["agrees"] is True

    def test_one_drug_is_not_a_ranking(self):
        rk = self._doc(a=[("cisplatin", 0.37)], b=[("cisplatin", 0.31)])["ranking"]
        assert rk["agrees"] is None
        assert "at least 2" in rk["verdict"]

    def test_a_queued_run_is_not_ranked(self):
        runs_a = [_fake_run("a0", SECTION_A, "cisplatin", 0.37),
                  _fake_run("a1", SECTION_A, "temozolomide", None,
                            status="queued")]
        runs_b = [_fake_run("b0", SECTION_B, "cisplatin", 0.31),
                  _fake_run("b1", SECTION_B, "temozolomide", 0.03)]
        project, resolve = _project(runs_a, runs_b)
        doc = compare_project(project, resolve)
        assert doc["ranking"]["drugs"] == ["cisplatin"]
        assert doc["specimens"][0]["runs_not_done"][0]["status"] == "queued"

    def test_a_missing_run_is_reported_not_swallowed(self):
        project, resolve = _project(
            [_fake_run("a0", SECTION_A, "cisplatin", 0.37)],
            [_fake_run("b0", SECTION_B, "cisplatin", 0.31)])
        project["specimens"][0]["run_ids"].append("r_gone")
        doc = compare_project(project, resolve)
        assert doc["specimens"][0]["runs_missing"] == ["r_gone"]
        assert doc["problems"] and "r_gone" in doc["problems"][0]


# --------------------------------------------------------------------------- #
# 6. The comparison — capability-aware flagging of the ABSOLUTE kill
# --------------------------------------------------------------------------- #
class TestKillComparability:
    def _doc(self, *, state="oxygen", measured_a=False, measured_b=False, **kw):
        """Two specimens x two drugs at an IDENTICAL RunSpec (``state`` applies
        to both sides, or the runs would not match at all)."""
        runs_a = [_fake_run(f"a{i}", SECTION_A, d, k, state_mode=state,
                            state_measured=measured_a)
                  for i, (d, k) in enumerate([("cisplatin", 0.37),
                                              ("temozolomide", 0.01)])]
        runs_b = [_fake_run(f"b{i}", SECTION_B, d, k, state_mode=state,
                            state_measured=measured_b)
                  for i, (d, k) in enumerate([("cisplatin", 0.31),
                                              ("temozolomide", 0.02)])]
        project, resolve = _project(runs_a, runs_b, **kw)
        return compare_project(project, resolve)

    def test_a_panel_without_hypoxia_blocks_the_absolute_delta(self):
        """THE case: the 319 panel cannot support an absolute kill, so no
        cross-specimen kill difference involving it is readable."""
        doc = self._doc()  # A has no hypoxia gene, B does
        assert doc["kill"], "the drugs are matched and should produce rows"
        for row in doc["kill"]:
            assert row["comparable"] is False
            joined = " ".join(row["blockers"])
            assert "hypoxia" in joined.lower()
            assert "cell_state_impact" in joined
        assert "NOT-COMPARABLE" in doc["headline"]
        # ...while the numbers themselves are still present, with their delta.
        row = next(r for r in doc["kill"] if r["drug"] == "cisplatin")
        kf = row["metrics"]["kill_fraction"]
        assert kf["values"] == {"s_00001": 0.37, "s_00002": 0.31}
        assert kf["delta"]["s_00002"] == pytest.approx(0.31 - 0.37)

    def test_the_unsupported_specimen_is_labelled(self):
        doc = self._doc()
        a = next(s for s in doc["specimens"] if s["id"] == "s_00001")
        b = next(s for s in doc["specimens"] if s["id"] == "s_00002")
        assert "ABSOLUTE KILL NOT SUPPORTED BY THIS PANEL" in a["labels"]
        assert a["absolute_kill_supported"] is False
        assert "UNSUPPORTED STATE ROUTE(S): markers" in a["labels"]
        assert b["labels"] == []
        assert b["absolute_kill_supported"] is True

    def test_two_supporting_panels_on_an_unmeasured_route_are_still_blocked(self):
        """Both panels COULD measure a state; neither run did. That is not a
        patient difference either — it is the default route on both sides."""
        doc = self._doc(hypoxia_a=True)
        for row in doc["kill"]:
            assert row["comparable"] is False
            assert any("no specimen used a MEASURED cell-state route" in b
                       for b in row["blockers"])

    def test_measured_against_assumed_is_refused(self):
        """Both sides asked for state_mode='markers' and both panels clear the
        hypoxia bar, but only one run came back with a MEASURED label. That
        delta is the state rule, not the patient."""
        doc = self._doc(hypoxia_a=True, state="markers",
                        measured_a=True, measured_b=False)
        row = next(r for r in doc["kill"] if r["drug"] == "cisplatin")
        assert row["comparable"] is False
        assert any("MEASURED on" in b and "ASSUMED on" in b
                   for b in row["blockers"])
        assert row["state_measured"] == {"s_00001": True, "s_00002": False}

    def test_both_measured_and_both_supported_is_comparable_with_caveats(self):
        doc = self._doc(hypoxia_a=True, state="markers",
                        measured_a=True, measured_b=True)
        for row in doc["kill"]:
            assert row["comparable"] is True
            assert row["blockers"] == []
            joined = " ".join(row["caveats"])
            assert "DIFFERENT panels" in joined      # panel confound survives
            assert "UPPER BOUND" in joined           # uniform sensitivity
            assert "uncalibrated" in joined          # oxygen level
            assert "NOT FOR CLINICAL USE" in joined
        assert "NOT-COMPARABLE" not in doc["headline"]


# --------------------------------------------------------------------------- #
# 7. The comparison — composition
# --------------------------------------------------------------------------- #
class TestComposition:
    def _rows(self, **kw):
        runs_a = [_fake_run("a0", SECTION_A, "cisplatin", 0.37)]
        runs_b = [_fake_run("b0", SECTION_B, "cisplatin", 0.31, report={
            "n_cells": 816769,
            "celltype_fractions": {"TUMOR": 0.7909, "ENDOTHELIAL": 0.0418},
            "celltype_counts": {"TUMOR": 646, "ENDOTHELIAL": 34179},
            "extent_um": [11472.0, 10860.0],
            "tissue_coverage_fraction": 0.7183,
            "mean_vascular_density": 0.0579,
            "mean_ecm_density": 0.4571,
            "mean_oxygen": 0.1232,
            "ecm_markers_used": ["VCAN", "DCN", "FBLN1"],
            "cell_state": {"mode": "oxygen", "source": "simulated O2",
                           "is_measured": False, "unmeasured_states": None,
                           "fractions_tumor": {"VIABLE": 0.0, "QUIESCENT": 0.03,
                                               "HYPOXIC": 0.97}},
        })]
        project, resolve = _project(runs_a, runs_b, **kw)
        doc = compare_project(project, resolve)
        return {r["key"]: r for r in doc["composition"]["rows"]}, doc

    def test_the_vascular_difference_is_shown(self):
        rows, doc = self._rows()
        v = rows["mean_vascular_density"]
        assert v["values"]["s_00001"] == pytest.approx(0.025)
        assert v["values"]["s_00002"] == pytest.approx(0.0579)
        assert v["ratio"]["s_00002"] == pytest.approx(0.0579 / 0.025)
        assert doc["composition"]["from_run"] == {"s_00001": "a0",
                                                  "s_00002": "b0"}

    def test_endothelial_fraction_and_areal_density_both_appear(self):
        """heldout_patient.md §1.1: the fraction FELL while the areal density
        ROSE on exactly these two sections. Showing only one inverts it."""
        rows, _ = self._rows()
        frac = rows["frac_endothelial"]
        dens = rows["endothelial_per_mm2"]
        assert frac["ratio"]["s_00002"] < 1.0
        assert dens["ratio"]["s_00002"] > 1.0
        assert frac["comparable"] is False   # different panels, different gates

    def test_ecm_and_oxygen_are_never_comparable(self):
        rows, _ = self._rows()
        assert rows["mean_ecm_density"]["comparable"] is False
        assert "no common absolute scale" in rows["mean_ecm_density"]["reason"]
        assert rows["ecm_markers_used"]["values"]["s_00002"] == [
            "VCAN", "DCN", "FBLN1"]
        assert rows["mean_oxygen"]["comparable"] is False
        assert rows["mean_oxygen"]["provenance"]["s_00001"]["kind"] == "assumed"

    def test_cell_counts_and_density_are_comparable(self):
        rows, _ = self._rows()
        assert rows["n_cells"]["comparable"] is True
        assert rows["n_cells"]["provenance"]["s_00001"]["kind"] == "measured"
        assert rows["cell_density_mm2"]["comparable"] is True
        assert rows["cell_density_mm2"]["values"]["s_00002"] > \
            rows["cell_density_mm2"]["values"]["s_00001"]

    def test_an_unmeasured_state_band_is_null_not_zero(self):
        runs_a = [_fake_run("a0", SECTION_A, "cisplatin", 0.37)]
        runs_a[0]["report"]["cell_state"] = {
            "mode": "markers", "source": "measured markers", "is_measured": True,
            "unmeasured_states": ["HYPOXIC"],
            "fractions_tumor": {"VIABLE": 0.298, "QUIESCENT": 0.702,
                                "HYPOXIC": 0.0},
        }
        runs_b = [_fake_run("b0", SECTION_B, "cisplatin", 0.31,
                            state_mode="markers", state_measured=True)]
        project, resolve = _project(runs_a, runs_b)
        rows = {r["key"]: r
                for r in compare_project(project, resolve)["composition"]["rows"]}
        hyp = rows["state_hypoxic"]
        assert hyp["values"]["s_00001"] is None, \
            "an UNMEASURED band must be null, never a measured-looking 0.0"
        assert hyp["provenance"]["s_00001"]["kind"] == "unknown"
        assert "computed and found to be zero" in hyp["provenance"]["s_00001"]["why"]
        assert hyp["comparable"] is False
        assert rows["state_viable"]["provenance"]["s_00001"]["kind"] == "measured"


# --------------------------------------------------------------------------- #
# 8. Routes
# --------------------------------------------------------------------------- #
@requires_both
class TestRoutes:
    def test_crud_and_association_over_http(self, client):
        assert client.get("/api/projects").json() == []

        r = client.post("/api/projects", json={"name": "S", "notes": "n"})
        assert r.status_code == 201, r.text
        pid = r.json()["id"]
        assert client.get("/api/projects").json()[0]["id"] == pid

        r = client.post(f"/api/projects/{pid}/specimens",
                        json={"section": SECTION_A, "label": "A"})
        assert r.status_code == 201, r.text
        sid = r.json()["id"]
        assert r.json()["capability"]["absolute_kill_supported"] is False

        r = client.post(f"/api/projects/{pid}/specimens",
                        json={"section": SECTION_B})
        assert r.status_code == 201
        assert r.json()["capability"]["absolute_kill_supported"] is True

        # duplicate section -> 409
        assert client.post(f"/api/projects/{pid}/specimens",
                           json={"section": SECTION_A}).status_code == 409
        # unknown section -> 400
        assert client.post(f"/api/projects/{pid}/specimens",
                           json={"section": "zzz"}).status_code == 400
        # unknown project -> 404
        assert client.get("/api/projects/p_99999").status_code == 404

        # attach a real run submitted through /api/runs (the worker is stopped,
        # so it stays queued — association does not require a finished run).
        posted = client.post("/api/runs", json={
            "section": SECTION_A, "dx_um": 30.0, "drug": "temozolomide",
            "sensitivity_mode": "uniform", "vascular_mode": "discrete",
            "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": None,
            "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
                   "n_doses": 1, "dose_interval_s": 86400.0},
            "total_time_s": 600.0})
        assert posted.status_code == 202, posted.text
        run_id = posted.json()["run_id"]

        att = client.post(f"/api/projects/{pid}/specimens/{sid}/runs",
                          json={"run_id": run_id})
        assert att.status_code == 201, att.text
        assert att.json()["attached"] is True
        assert client.post(f"/api/projects/{pid}/specimens/{sid}/runs",
                           json={"run_id": "r_99999"}).status_code == 404
        assert client.post(f"/api/projects/{pid}/specimens/{sid}/runs",
                           json={}).status_code == 400

        got = client.get(f"/api/projects/{pid}/specimens/{sid}").json()
        assert got["run_ids"] == [run_id]
        assert got["runs"][0]["run_id"] == run_id

        det = client.delete(
            f"/api/projects/{pid}/specimens/{sid}/runs/{run_id}")
        assert det.status_code == 200 and det.json()["detached"] is True
        # the run itself is untouched
        assert client.get(f"/api/runs/{run_id}").status_code == 200

        patched = client.patch(f"/api/projects/{pid}", json={"name": "S2"})
        assert patched.json()["name"] == "S2"
        assert client.delete(f"/api/projects/{pid}").json()["deleted"] is True
        assert client.get(f"/api/projects/{pid}").status_code == 404

    def test_compare_route_and_bad_metric(self, client):
        pid = client.post("/api/projects", json={"name": "S"}).json()["id"]
        client.post(f"/api/projects/{pid}/specimens", json={"section": SECTION_A})
        doc = client.get(f"/api/projects/{pid}/compare")
        assert doc.status_code == 200, doc.text
        body = doc.json()
        assert body["ranking"]["agrees"] is None
        assert body["specimens"][0]["labels"]
        assert client.get(f"/api/projects/{pid}/compare",
                          params={"metric": "nope"}).status_code == 400

    def test_section_capability_route(self, client):
        a = client.get(f"/api/sections/{SECTION_A}/capability").json()
        b = client.get(f"/api/sections/{SECTION_B}/capability").json()
        assert a["absolute_kill_supported"] is False
        assert b["absolute_kill_supported"] is True
        assert client.get("/api/sections/zzz/capability").status_code == 400

    def test_health_reports_the_project_store(self, client):
        h = client.get("/api/health").json()
        assert "projects_dir" in h and h["projects"] == 0
