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
   scientific claim. Every comparison is SAME-COMPOUND: one compound's own value
   on specimen A against that same compound's value on specimen B, plus the
   tissue contrast between them. An absolute-kill delta is refused unless BOTH
   panels can support a measured cell state AND both runs used one. This is
   ``results/cell_state_impact.md`` §5 (a panel with no hypoxia-response
   transcript cannot support an absolute kill under any route, because kill
   moves 0.3252 absolute / 1.77x with the route alone) applied to the
   two-specimen case ``results/heldout_patient.md`` §3.5 measured: the absolute
   kills did not travel.

**Rank agreement used to be item 3's primary result and is WITHDRAWN.** The
tests below no longer assert that a cross-compound ordering agrees, because
``compare_project`` no longer computes one: ``results/pk_ranking_robustness.md``
measures per-compound plasma exposure spanning 113.1x across the compounds the
GBM panel ranks with the published ordering breaking at a spread of 1.53x, and
``results/capacity_ranking_robustness.md`` moves the breast top-1 pick under this
repo's own structural capacity rule. Every assertion that guarded the ranking was
CONVERTED to the same-compound contract rather than dropped — the matching rules
those tests pinned down (which runs pair, which are excluded, which compounds are
common to both specimens) are all still load-bearing.

The comparison tests build their run payloads by hand rather than simulating, so
the FLAGGING LOGIC is tested exhaustively and in milliseconds. The end-to-end
path (real sections, real runs) is exercised over HTTP by ``TestRoutes``, which
is skipped unless both sections are cached.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import sys
import time

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
from app.server import projects as projects_mod  # noqa: E402
from app.server.projects import (  # noqa: E402
    ProjectError,
    ProjectStore,
    capability_report,
    compare_project,
    run_composition,
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
    # runs_root is pinned to the tmp dir, NOT left to default to the real
    # app/runs/: the thumbnail reads run documents off disk, and a test that
    # silently read this developer's 96 real runs would pass or fail on what
    # happens to be cached in this checkout.
    return ProjectStore(root=str(tmp_path / "projects"),
                        runs_root=str(tmp_path / "runs"))


@pytest.fixture
def client(tmp_path):
    """A TestClient whose runs AND projects are isolated from the real dirs."""
    manager = JobManager(root=str(tmp_path / "runs"), start_worker=False)
    app = create_app(manager, ProjectStore(root=str(tmp_path / "projects"),
                                           runs_root=str(tmp_path / "runs")))
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
        """Converted from an assertion on ``ranking.agrees is None``: the same
        pair of runs is still unmatched, and the document still says so — now on
        the headline, which is the only verdict left."""
        runs_a = [_fake_run("a1", SECTION_A, "cisplatin", 0.37),
                  _fake_run("a2", SECTION_A, "temozolomide", 0.01)]
        runs_b = [_fake_run("b1", SECTION_B, "cisplatin", 0.31, dx=12.0),
                  _fake_run("b2", SECTION_B, "temozolomide", 0.02, dx=12.0)]
        project, resolve = _project(runs_a, runs_b)
        doc = compare_project(project, resolve)
        assert doc["headline_tag"] == "NOT EVALUABLE"
        assert "NOT EVALUABLE" in doc["headline"]
        assert doc["kill"] == []
        assert "ranking" not in doc


# --------------------------------------------------------------------------- #
# 5. The comparison — the SAME compound across specimens
#
# Was "5. The comparison — ranking is the primary result". Every test in this
# class guarded a cross-compound ordering; each one has been converted to the
# contract that replaced it (the same compound on two specimens), because the
# MATCHING rules they pinned down did not change when the ordering went away.
# --------------------------------------------------------------------------- #
class TestSameCompoundAcrossSpecimens:
    #: Three compounds run on both sections. These were the values whose
    #: ORDERING results/heldout_patient.md §3.3 found identical; the ordering is
    #: no longer computed or asserted, and the values are kept only because they
    #: are real and exercise the matching.
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

    def test_the_matched_compounds_are_the_headline(self):
        """Converted from ``test_ranking_agreement_is_detected_and_is_the_headline``.
        One assertion per line of the old contract, restated on the new one."""
        doc = self._doc()
        # was: rk["agrees"] is True -> there is no agreement verdict to give
        assert "ranking" not in doc
        # ...and the tag is about COMPARABILITY, not agreement. Specimen A's
        # panel carries no hypoxia transcript in this fixture, so the absolute
        # levels are blocked while the matching itself succeeded.
        assert doc["headline_tag"] == "SAME COMPOUND · NOT COMPARABLE"
        # was: rk["drugs"] == [...] -> one row per compound matched on both
        assert [r["drug"] for r in doc["kill"]] == [
            "cisplatin", "doxorubicin", "temozolomide"]
        # was: the two specimens' ORDERS are equal -> every row holds the
        # compound fixed and carries both specimens, which is the whole contract
        for row in doc["kill"]:
            assert row["compound"] == row["drug"]
            assert row["varies"] == "specimen"
            assert set(row["metrics"]["kill_fraction"]["values"]) == {
                "s_00001", "s_00002"}
        # was: kendall_tau_b == 1.0 -> the matched-set count, which is what that
        # statistic was computed over
        assert doc["matched_sets"][0]["n_matched_drugs"] == 3
        assert doc["matched_sets"][0]["is_primary"] is True
        # was: swapped == [] -> no cross-compound quantity is emitted at all
        assert all("swapped" not in r and "kendall_tau_b" not in r
                   for r in doc["kill"])
        # was: headline startswith "RANKING AGREES"
        assert doc["headline"].startswith("3 compounds matched across 2 specimens")
        assert "No ordering of compounds is computed" in doc["headline"]
        # was: reading_order[1].startswith("ranking")
        assert doc["reading_order"][1].startswith("kill")
        assert not any("ranking" in entry for entry in doc["reading_order"])

    def test_reordering_the_values_changes_no_verdict_because_there_is_none(self):
        """Converted from ``test_a_swap_is_detected_and_named``. The same input
        that used to be reported as a rank SWAP now produces exactly what it
        should: each compound compared against itself across the two specimens,
        and no statement about how the compounds sit relative to each other."""
        swapped = [("doxorubicin", 0.40), ("cisplatin", 0.30),
                   ("temozolomide", 0.03)]
        doc = self._doc(b=swapped)
        # was: rk["agrees"] is False
        assert "ranking" not in doc
        assert doc["headline_tag"] == "SAME COMPOUND · NOT COMPARABLE"
        rows = {r["drug"]: r for r in doc["kill"]}
        # was: pw["agrees"] is False -> the two compounds that "swapped" moved in
        # OPPOSITE cross-specimen directions, which is a same-compound fact and
        # is still reported, per compound, as its own ratio
        assert rows["doxorubicin"]["metrics"]["kill_fraction"]["ratio"]["s_00002"] \
            == pytest.approx(0.40 / 0.2165)
        assert rows["cisplatin"]["metrics"]["kill_fraction"]["ratio"]["s_00002"] \
            == pytest.approx(0.30 / 0.3699)
        # was: sorted(pw["swapped"][0]) == ["cisplatin", "doxorubicin"]
        assert sorted(rows) == ["cisplatin", "doxorubicin", "temozolomide"]
        # was: pw["discordant_pairs"] == 1 -> no pair of DIFFERENT compounds is
        # placed on one axis anywhere in the document
        assert all(len({r["compound"]}) == 1 for r in doc["kill"])
        # was: "DISAGREES" in rk["verdict"]
        assert "DISAGREE" not in doc["headline"]
        assert "RANKING" not in doc["headline"]

    def test_only_the_compounds_both_specimens_ran_are_compared(self):
        """Converted from ``test_ranking_uses_only_the_drugs_both_specimens_ran``.
        Same matching rule, now guarding the same-compound rows."""
        doc = self._doc(b=self._B + [("lomustine", 0.22)])
        drugs = [r["drug"] for r in doc["kill"]]
        assert "lomustine" not in drugs
        assert drugs == ["cisplatin", "doxorubicin", "temozolomide"]

    def test_one_compound_on_both_specimens_is_still_a_comparison(self):
        """Converted from ``test_one_drug_is_not_a_ranking``. One compound was
        the case a ranking could NOT be taken over; it is exactly the case the
        product is now for, so the assertion flips from refusal to acceptance."""
        doc = self._doc(a=[("cisplatin", 0.37)], b=[("cisplatin", 0.31)])
        # was: rk["agrees"] is None
        assert doc["headline_tag"] != "NOT EVALUABLE"
        assert [r["drug"] for r in doc["kill"]] == ["cisplatin"]
        # was: "at least 2" in rk["verdict"]
        assert "1 compound matched across 2 specimens" in doc["headline"]

    def test_a_queued_run_is_excluded_from_the_comparison(self):
        """Converted from ``test_a_queued_run_is_not_ranked``. A run that has
        not finished has no metrics to compare, so it is dropped from the
        matched set and named in ``runs_not_done`` — unchanged behaviour, new
        name for the block it is excluded from."""
        runs_a = [_fake_run("a0", SECTION_A, "cisplatin", 0.37),
                  _fake_run("a1", SECTION_A, "temozolomide", None,
                            status="queued")]
        runs_b = [_fake_run("b0", SECTION_B, "cisplatin", 0.31),
                  _fake_run("b1", SECTION_B, "temozolomide", 0.03)]
        project, resolve = _project(runs_a, runs_b)
        doc = compare_project(project, resolve)
        assert [r["drug"] for r in doc["kill"]] == ["cisplatin"]
        assert doc["specimens"][0]["runs_not_done"][0]["status"] == "queued"

    def test_the_withdrawn_ordering_is_cited_where_it_was_withdrawn(self):
        """New, and the reason the class above can be trusted: the document
        states what it compares and names both measurements that took the
        ordering away, so the removal is on screen rather than only in git."""
        doc = self._doc()
        assert "SAME compound across specimens" in doc["compares"]
        assert "results/pk_ranking_robustness.md" in doc["compares"]
        assert "results/capacity_ranking_robustness.md" in doc["compares"]
        claim = " ".join(doc["kill"][0]["caveats"])
        assert "RELATIVE exposure between regions is the defensible output" in claim
        assert "comparison BETWEEN compounds is not" in claim
        # the two invariances, at their two DIFFERENT strengths
        assert "1.41e-14" in claim and "5,200x span" in claim
        assert "1e-9" in claim and "only where the field has converged" in claim
        # ...and the limit that stops it being read as a general guarantee
        assert "NEITHER INVARIANCE PROTECTS AGAINST ERROR IN THE TRANSPORT LAW" in claim
        assert "results/descriptor_falsification.md" in claim
        assert "A LIST" in claim and "UNSWEPT" in claim

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
        assert doc["headline_tag"] == "SAME COMPOUND · MATCHED"
        assert doc["headline_level"] == "ok"


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
        # was: body["ranking"]["agrees"] is None — one specimen, no run, so
        # nothing is matched and the route says so on the headline instead.
        assert "ranking" not in body
        assert body["headline_tag"] == "NOT EVALUABLE"
        assert body["specimens"][0]["labels"]
        assert client.get(f"/api/projects/{pid}/compare",
                          params={"metric": "nope"}).status_code == 400
        # ``metric`` still selects a column, so a valid one is still accepted and
        # is echoed back on the document it produced.
        ok = client.get(f"/api/projects/{pid}/compare",
                        params={"metric": "penetration_length_um"})
        assert ok.status_code == 200, ok.text
        assert ok.json()["metric"] == "penetration_length_um"

    def test_section_capability_route(self, client):
        a = client.get(f"/api/sections/{SECTION_A}/capability").json()
        b = client.get(f"/api/sections/{SECTION_B}/capability").json()
        assert a["absolute_kill_supported"] is False
        assert b["absolute_kill_supported"] is True
        assert client.get("/api/sections/zzz/capability").status_code == 400

    def test_health_reports_the_project_store(self, client):
        h = client.get("/api/health").json()
        assert "projects_dir" in h and h["projects"] == 0


# --------------------------------------------------------------------------- #
# 6. The workspace additions — last_opened_at, duplicate, thumbnail
#
# docs/WORKSPACE_CONTRACT.md §1, §3 and §4. Three things the home page needs and
# three different ways they could lie:
#
#   §3 last_opened_at — by claiming an open that never happened. Every project
#      in this store predates the field, and the only honest value for a
#      document that never recorded one is null. Backfilling it with now() would
#      put every old project at the top of a recency sort; backfilling it with
#      created_at would silently redefine "opened" as "created".
#   §1 duplicate — by copying, re-running, or disturbing a run. A project owns
#      nothing but ids (§2); a duplicate holds the SAME ids, and app/runs/ must
#      come out of the operation byte-identical.
#   §4 thumbnail — by drawing a chart nobody measured. An empty stacked bar
#      reads as "all one cell type" and a zero reads as "measured, and it was
#      zero". Neither has been computed when no run has finished, so the only
#      other output is a refusal carrying the true reason.
#
# These tests build run documents on disk by hand rather than simulating,
# because the thumbnail reads app/runs/<id>/result.json and that file format —
# not the JobManager's memory — is what it has to survive.
# --------------------------------------------------------------------------- #
#: Written into a run document and asserted as LITERALS below. They are NOT read
#: back out of the module or out of the payload they are checking: a test whose
#: expected value is derived from the thing under test moves with the bug.
_THUMB_FRACTIONS = {"TUMOR": 0.5823855993347519, "STROMAL": 0.11776359233986353,
                    "IMMUNE": 0.19480519480519481,
                    "ENDOTHELIAL": 0.10504561352018979, "NECROTIC": 0.0}
_THUMB_N_CELLS = 40887


def _write_run_doc(runs_root, run_id, *, status="done", section=SECTION_A,
                   report="default", drug="cisplatin"):
    """Persist one run document exactly where ``JobManager._persist`` puts it."""
    d = os.path.join(runs_root, run_id)
    os.makedirs(d, exist_ok=True)
    if report == "default":
        report = {
            "n_cells": _THUMB_N_CELLS,
            "celltype_fractions": dict(_THUMB_FRACTIONS),
            "celltype_counts": {"TUMOR": 23812, "STROMAL": 4815, "IMMUNE": 7965,
                                "ENDOTHELIAL": 4295, "NECROTIC": 0},
            "source": "xenium",
        }
    doc = {"run_id": run_id, "status": status, "created_at": 1000.0,
           "spec": {"section": section, "drug": drug, "dx_um": 30.0},
           "metrics": {"kill_fraction": 0.1}, "report": report}
    with open(os.path.join(d, "result.json"), "w") as fh:
        json.dump(doc, fh, indent=1)
    return doc


def _seed_project(store, name, specimens, *, legacy=False):
    """Create a project and give it ``specimens`` by writing project.json.

    ``legacy=True`` writes the document a pre-2026-09-20 build would have
    written: no ``last_opened_at`` key at all. That is the migration case, and
    it cannot be produced through the store's own writers any more.
    """
    p = store.create(name)
    path = os.path.join(store.root, p["id"], "project.json")
    with open(path) as fh:
        doc = json.load(fh)
    doc["specimens"] = specimens
    if legacy:
        doc.pop("last_opened_at", None)
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=1)
    return p["id"]


def _fingerprint(path):
    """(bytes, mtime_ns) — what "this file was not touched" actually means."""
    st = os.stat(path)
    with open(path, "rb") as fh:
        return fh.read(), st.st_mtime_ns


class TestLastOpenedAt:
    """§3. ``updated_at`` is the document; ``last_opened_at`` is a human."""

    def test_a_project_that_predates_the_field_reads_as_null(self, store):
        pid = _seed_project(store, "legacy study", [], legacy=True)
        on_disk = json.load(open(os.path.join(store.root, pid, "project.json")))
        assert "last_opened_at" not in on_disk, "the legacy fixture must be legacy"

        doc = store.get(pid, refresh_capability=False)
        assert doc["last_opened_at"] is None
        (row,) = store.list()
        assert row["last_opened_at"] is None

    def test_the_migration_is_never_now_and_never_created_at(self, store):
        """The two wrong defaults, named, because both look harmless."""
        pid = _seed_project(store, "legacy study", [], legacy=True)
        created = store.get(pid, refresh_capability=False)["created_at"]
        before = time.time()
        value = store.get(pid, refresh_capability=False)["last_opened_at"]
        after = time.time()
        assert value is None
        # not "now"
        assert not (isinstance(value, float) and before <= value <= after)
        # not created_at
        assert value != created

    def test_reading_a_legacy_project_rewrites_no_document(self, store):
        """A read must not write. Backfilling on read would rewrite the whole
        store on the first home-page load and move every mtime."""
        pid = _seed_project(store, "legacy study", [], legacy=True)
        path = os.path.join(store.root, pid, "project.json")
        before = _fingerprint(path)

        store.get(pid, refresh_capability=False)
        store.list()
        store.list()
        store.thumbnail(pid)

        assert _fingerprint(path) == before
        assert b"last_opened_at" not in before[0]
        assert b"last_opened_at" not in _fingerprint(path)[0]

    def test_mark_opened_survives_the_round_trip(self, store):
        pid = _seed_project(store, "S", [], legacy=True)
        out = store.mark_opened(pid)
        assert out["ok"] is True
        stamped = out["last_opened_at"]
        assert isinstance(stamped, float)

        reopened = ProjectStore(root=store.root, runs_root=store.runs_root)
        assert reopened.get(pid, refresh_capability=False)["last_opened_at"] \
            == stamped
        assert reopened.list()[0]["last_opened_at"] == stamped

    def test_mark_opened_does_not_move_updated_at(self, store):
        """THE distinction. ``updated_at`` says a background job attached a run;
        it cannot say that if looking at a project counts as touching it."""
        pid = _seed_project(store, "S", [])
        was = store.get(pid, refresh_capability=False)["updated_at"]
        time.sleep(0.01)
        out = store.mark_opened(pid)
        now = store.get(pid, refresh_capability=False)

        assert now["updated_at"] == was
        assert out["updated_at"] == was
        assert now["last_opened_at"] > was
        # ...and a real document change still does move it.
        store.update(pid, notes="edited")
        assert store.get(pid, refresh_capability=False)["updated_at"] > was

    def test_mark_opened_is_idempotent_and_404s_on_a_stranger(self, store):
        pid = _seed_project(store, "S", [])
        first = store.mark_opened(pid)["last_opened_at"]
        time.sleep(0.01)
        second = store.mark_opened(pid)["last_opened_at"]
        assert second > first
        assert set(store.get(pid, refresh_capability=False)) == set(
            store.get(pid, refresh_capability=False))
        with pytest.raises(ProjectError) as exc:
            store.mark_opened("p_99999")
        assert exc.value.status == 404

    def test_the_summary_carries_the_contract_field_names(self, store):
        """``ProjectSummary`` is a hard interface: another agent is building
        GET /api/workspace against these keys."""
        _seed_project(store, "S", [])
        (row,) = store.list()
        for key in ("id", "name", "notes", "created_at", "updated_at",
                    "last_opened_at", "n_specimens", "n_runs", "sections",
                    "thumbnail"):
            assert key in row, f"ProjectSummary lost {key!r}"


class TestDuplicate:
    """§1. Specimens and run REFERENCES. Nothing in app/runs/ moves."""

    @staticmethod
    def _two_specimens():
        return [
            {"id": "s_00001", "label": "Patient A", "section": SECTION_A,
             "notes": "n", "created_at": 500.0,
             "run_ids": ["r_00001", "r_00002"],
             "capability": {"section": SECTION_A, "panel": "brain_319"}},
            {"id": "s_00002", "label": "Patient B", "section": SECTION_B,
             "notes": "", "created_at": 501.0, "run_ids": ["r_00003"],
             "capability": {"section": SECTION_B, "panel": "io_480"}},
        ]

    def test_specimens_and_run_references_are_copied(self, store):
        pid = _seed_project(store, "Original", self._two_specimens())
        dup = store.duplicate(pid)

        assert dup["id"] != pid
        assert dup["name"] == "Original (copy)"
        assert dup["duplicated_from"] == pid
        assert dup["last_opened_at"] is None
        assert [s["id"] for s in dup["specimens"]] == ["s_00001", "s_00002"]
        assert [s["label"] for s in dup["specimens"]] == ["Patient A", "Patient B"]
        assert [s["section"] for s in dup["specimens"]] == [SECTION_A, SECTION_B]
        assert [s["run_ids"] for s in dup["specimens"]] == \
            [["r_00001", "r_00002"], ["r_00003"]]

        # ...and it survives the round trip.
        reopened = ProjectStore(root=store.root, runs_root=store.runs_root)
        got = reopened.get(dup["id"], refresh_capability=False)
        assert [s["run_ids"] for s in got["specimens"]] == \
            [["r_00001", "r_00002"], ["r_00003"]]

    def test_the_original_document_is_not_touched(self, store):
        pid = _seed_project(store, "Original", self._two_specimens())
        path = os.path.join(store.root, pid, "project.json")
        before = _fingerprint(path)

        store.duplicate(pid, name="A copy")

        assert _fingerprint(path) == before, \
            "duplicate() wrote to the project it was copying FROM"
        src = store.get(pid, refresh_capability=False)
        assert src["name"] == "Original"
        assert src["last_opened_at"] is None
        assert "duplicated_from" not in src

    def test_no_run_is_re_executed_or_copied(self, store):
        """app/runs/ is the single store. A duplicate is ids, not results."""
        runs = store.runs_root
        for rid in ("r_00001", "r_00002", "r_00003"):
            _write_run_doc(runs, rid, section=SECTION_A)
        before = {name: _fingerprint(os.path.join(runs, name, "result.json"))
                  for name in sorted(os.listdir(runs))}

        pid = _seed_project(store, "Original", self._two_specimens())
        store.duplicate(pid)

        after = {name: _fingerprint(os.path.join(runs, name, "result.json"))
                 for name in sorted(os.listdir(runs))}
        assert after == before, "duplicate() created, copied or rewrote a run"
        assert sorted(after) == ["r_00001", "r_00002", "r_00003"]

    def test_both_projects_resolve_the_same_run(self, store):
        """Two specimens in two projects legitimately point at ONE run id."""
        _write_run_doc(store.runs_root, "r_00001", section=SECTION_A)
        pid = _seed_project(store, "Original", self._two_specimens())
        dup = store.duplicate(pid)["id"]

        seen = []

        def resolver(run_id):
            seen.append(run_id)
            return {"run_id": run_id, "status": "done",
                    "spec": {"section": SECTION_A},
                    "report": {"n_cells": _THUMB_N_CELLS,
                               "celltype_fractions": dict(_THUMB_FRACTIONS)}}

        a = store.thumbnail(pid, resolver)
        b = store.thumbnail(dup, resolver)
        assert a["run_id"] == b["run_id"] == "r_00001"
        assert a["type_fractions"] == b["type_fractions"]
        assert seen == ["r_00001", "r_00001"]

        # The same id in both documents, read back off disk.
        ids = lambda p: [s["run_ids"] for s in  # noqa: E731
                         store.get(p, refresh_capability=False)["specimens"]]
        assert ids(pid) == ids(dup) == [["r_00001", "r_00002"], ["r_00003"]]

    def test_an_explicit_name_and_the_error_cases(self, store):
        pid = _seed_project(store, "Original", [])
        assert store.duplicate(pid, name="  Renamed  ")["name"] == "Renamed"
        with pytest.raises(ProjectError) as exc:
            store.duplicate(pid, name="   ")
        assert exc.value.status == 400
        with pytest.raises(ProjectError) as exc:
            store.duplicate("p_99999")
        assert exc.value.status == 404


class TestThumbnail:
    """§4. A measured composition, or a refusal that says why. Never a shape."""

    _SPECIMEN = [{"id": "s_00001", "label": "A", "section": SECTION_A,
                  "run_ids": ["r_00001"], "capability": {}}]

    def test_a_finished_run_supplies_the_composition_it_recorded(self, store):
        _write_run_doc(store.runs_root, "r_00001", section=SECTION_A)
        pid = _seed_project(store, "S", self._SPECIMEN)

        t = store.thumbnail(pid)
        assert t["kind"] == "composition"
        assert t["section"] == SECTION_A
        # Literals, not values read back out of the payload under test.
        assert t["n_cells"] == 40887
        assert t["type_fractions"]["TUMOR"] == 0.5823855993347519
        assert t["type_fractions"]["ENDOTHELIAL"] == 0.10504561352018979
        assert t["type_fractions"]["NECROTIC"] == 0.0
        assert set(t["type_fractions"]) == {"TUMOR", "STROMAL", "IMMUNE",
                                            "ENDOTHELIAL", "NECROTIC"}
        assert store.list()[0]["thumbnail"] == t

    @pytest.mark.parametrize("status", ["queued", "running", "error",
                                        "cancelled"])
    def test_an_unfinished_run_is_never_drawn(self, store, status):
        """A failed run still carries a report — the reconstruction succeeds
        before the solve fails. Drawing it presents a run that did not finish
        as one that did."""
        _write_run_doc(store.runs_root, "r_00001", status=status,
                       section=SECTION_A)
        pid = _seed_project(store, "S", self._SPECIMEN)

        t = store.thumbnail(pid)
        assert t["kind"] == "none"
        assert status in t["reason"]
        assert "type_fractions" not in t and "n_cells" not in t

    def test_no_run_at_all_refuses_with_a_true_reason(self, store):
        empty = _seed_project(store, "no specimens", [])
        t = store.thumbnail(empty)
        assert t["kind"] == "none" and "no specimens" in t["reason"]

        no_runs = _seed_project(store, "no runs", [
            {"id": "s_00001", "label": "A", "section": SECTION_A,
             "run_ids": [], "capability": {}}])
        t = store.thumbnail(no_runs)
        assert t["kind"] == "none" and "no run is attached" in t["reason"]
        assert "1 specimen" in t["reason"]

    def test_a_missing_run_is_reported_not_drawn_as_empty(self, store):
        pid = _seed_project(store, "S", self._SPECIMEN)   # r_00001 never written
        t = store.thumbnail(pid)
        assert t["kind"] == "none"
        assert "not in app/runs/" in t["reason"]

    def test_a_finished_run_with_no_composition_says_so(self, store):
        """Zeros are not the answer to "nothing was measured"."""
        _write_run_doc(store.runs_root, "r_00001", report={"n_cells": 40887})
        pid = _seed_project(store, "S", self._SPECIMEN)
        t = store.thumbnail(pid)
        assert t["kind"] == "none"
        assert "celltype_fractions" in t["reason"]
        assert "r_00001" in t["reason"]

    def test_the_first_usable_run_wins_and_the_scan_stops_there(self, store):
        """Cost: list() runs on every home-page load. The scan must not read
        every run of every project to draw one card."""
        specimens = [{"id": "s_00001", "label": "A", "section": SECTION_A,
                      "run_ids": ["r_00001", "r_00002", "r_00003"],
                      "capability": {}}]
        pid = _seed_project(store, "S", specimens)
        seen = []

        def resolver(run_id):
            seen.append(run_id)
            return {"run_id": run_id, "status": "done",
                    "spec": {"section": SECTION_A},
                    "report": {"n_cells": _THUMB_N_CELLS,
                               "celltype_fractions": dict(_THUMB_FRACTIONS)}}

        t = store.thumbnail(pid, resolver)
        assert t["run_id"] == "r_00001"
        assert seen == ["r_00001"], \
            f"the scan read {len(seen)} run(s) to draw one card: {seen}"

    def test_the_cache_does_not_outlive_the_run_document(self, store):
        """A run that finishes AFTER the first home-page load must appear. A
        thumbnail cached against the PROJECT would not: a run finishing does
        not touch project.json."""
        _write_run_doc(store.runs_root, "r_00001", status="running")
        pid = _seed_project(store, "S", self._SPECIMEN)
        assert store.thumbnail(pid)["kind"] == "none"

        time.sleep(0.01)
        _write_run_doc(store.runs_root, "r_00001", status="done")
        assert store.thumbnail(pid)["kind"] == "composition"

        # ...and a run removed from app/runs/ stops being drawn.
        shutil.rmtree(os.path.join(store.runs_root, "r_00001"))
        gone = store.thumbnail(pid)
        assert gone["kind"] == "none" and "not in app/runs/" in gone["reason"]

    def test_run_composition_refuses_the_shapes_it_cannot_measure(self):
        base = {"run_id": "r_1", "status": "done", "spec": {"section": "sec"},
                "report": {"n_cells": 10, "celltype_fractions": {"TUMOR": 1.0}}}
        assert run_composition(base)["type_fractions"] == {"TUMOR": 1.0}
        assert run_composition(dict(base, status="error")) is None
        assert run_composition(dict(base, report={})) is None
        assert run_composition(dict(base, report={"n_cells": 10,
                                                  "celltype_fractions": {}})) is None
        assert run_composition(dict(base, report={
            "celltype_fractions": {"TUMOR": 1.0}})) is None      # no n_cells
        assert run_composition(dict(base, report={
            "n_cells": 10, "celltype_fractions": {"TUMOR": "0.5"}})) is None
        assert run_composition(dict(base, report={
            "n_cells": 10, "celltype_fractions": {"TUMOR": float("nan")}})) is None
        assert run_composition(None) is None

    def test_a_run_id_is_never_a_path(self, store):
        """SECURITY. ``%2e%2e`` was a read-and-delete traversal through
        ProjectStore._dir (fixed f9d443d). The thumbnail opens
        app/runs/<run_id>/result.json, so a run_id is a path segment there too,
        and a project document is a file an operator can edit."""
        # The runs root MUST exist, or `<root>/../result.json` fails to stat for
        # the boring reason that `<root>` is not there and this test would pass
        # while never reaching the traversal it is about.
        os.makedirs(store.runs_root, exist_ok=True)
        outside = os.path.dirname(store.runs_root.rstrip(os.sep))
        with open(os.path.join(outside, "result.json"), "w") as fh:
            json.dump({"run_id": "canary", "status": "done",
                       "spec": {"section": "SECRET_SECTION"},
                       "report": {"n_cells": 1,
                                  "celltype_fractions": {"CANARY": 1.0}}}, fh)

        for bad in ("..", "../..", "r_0000", "", "r_00001/../..",
                    "/etc/passwd", "r_00001\x00"):
            pid = _seed_project(store, f"bad {bad!r}", [
                {"id": "s_00001", "label": "A", "section": SECTION_A,
                 "run_ids": [bad], "capability": {}}])
            t = store.thumbnail(pid)
            blob = json.dumps(t)
            assert "CANARY" not in blob, f"{bad!r} read a file outside app/runs/"
            assert "SECRET_SECTION" not in blob
            assert t["kind"] == "none", f"{bad!r} produced a thumbnail"
            assert store._run_path(bad) is None, f"{bad!r} was turned into a path"


# --------------------------------------------------------------------------- #
# A specimen id is not a path segment — and that is why there is no
# SPECIMEN_ID_RE
# --------------------------------------------------------------------------- #
class TestASpecimenIdReachesNothingThatNeedsAGuard:
    """SECURITY, and the standing evidence behind a DELETED regex.

    ``app/server/projects.py`` defined ``SPECIMEN_ID_RE`` for weeks and applied
    it to nothing. That is the exact shape of the traversal this file already
    suffered — ``PROJECT_ID_RE`` guarded only the directory listing while
    ``_dir`` joined a request-supplied id onto the store root, giving
    ``%2e%2e`` a read-AND-DELETE primitive (fixed f9d443d). The regex was
    deleted rather than wired, because a specimen id reaches nothing a regex
    could protect.

    A deletion is only safe while that stays true, so the premise is PINNED
    here rather than asserted in a comment: the day someone joins a specimen id
    to a path, opens it, or makes it a dict key, this class goes red and the id
    needs the treatment ``run_id`` gets in ``ProjectStore._run_path``.
    """

    #: Every shape the traversal that bit ``_dir`` came in, plus the
    #: key-shadowing shapes a dict lookup would have been vulnerable to.
    HOSTILE = ("..", "../..", "%2e%2e", "/etc/passwd", "s_00001/../..",
               "s_00001\x00", "....//", "__proto__", "constructor",
               "specimens", "id", "s_" + "9" * 200, "", "S_00001", "s_1")

    @staticmethod
    def _seed(store):
        return _seed_project(store, "traversal", [
            {"id": "s_00001", "label": "A", "section": SECTION_A,
             "run_ids": ["r_00001"], "capability": {}}])

    def test_no_specimen_id_ever_reaches_os_path_join(self, store):
        """The load-bearing one. ``os.path.join`` is instrumented AS THIS MODULE
        SEES IT, so any future line that builds a path out of a specimen id is
        caught whether or not the resulting path escapes anything."""
        pid = self._seed(store)
        marker = "MARKER_SPECIMEN_ID_../../.."
        seen = []
        real_join = projects_mod.os.path.join

        def spy(*parts):
            seen.append(parts)
            return real_join(*parts)

        projects_mod.os.path.join = spy
        try:
            for call in (
                lambda: ProjectStore._find_specimen(
                    store.get(pid, refresh_capability=False), marker),
                lambda: store.update_specimen(pid, marker, label="x"),
                lambda: store.remove_specimen(pid, marker),
                lambda: store.attach_run(pid, marker, "r_00001",
                                         {"spec": {"section": SECTION_A}}),
                lambda: store.detach_run(pid, marker, "r_00001"),
            ):
                with pytest.raises(ProjectError) as exc:
                    call()
                assert exc.value.status == 404
        finally:
            projects_mod.os.path.join = real_join

        # The premise, checked rather than assumed: the instrumentation must
        # have been REACHED, or this test passes by observing nothing.
        assert seen, "os.path.join was never called; the spy was not installed"
        hits = [p for p in seen if any(marker in str(x) for x in p)]
        assert hits == [], (
            f"a specimen id was joined into a filesystem path: {hits}")

    def test_a_hostile_specimen_id_neither_reads_nor_deletes_a_neighbour(
            self, store):
        """The canary sits one level above the store root — exactly the object
        ``%2e%2e`` reached through ``_dir``."""
        os.makedirs(store.root, exist_ok=True)
        outside = os.path.dirname(store.root.rstrip(os.sep))
        canary = os.path.join(outside, "project.json")
        with open(canary, "w") as fh:
            json.dump({"id": "CANARY", "name": "CANARY_SECRET",
                       "specimens": [{"id": "s_00001", "label": "SECRET_LABEL",
                                      "section": "SECRET_SECTION"}]}, fh)

        pid = self._seed(store)
        for bad in self.HOSTILE:
            for call in (
                lambda b=bad: store.update_specimen(pid, b, label="x"),
                lambda b=bad: store.remove_specimen(pid, b),
                lambda b=bad: store.detach_run(pid, b, "r_00001"),
            ):
                with pytest.raises(ProjectError) as exc:
                    call()
                blob = str(exc.value)
                assert "SECRET_LABEL" not in blob, f"{bad!r} read the canary"
                assert "CANARY_SECRET" not in blob
            assert os.path.exists(canary), f"{bad!r} DELETED the canary"

        # ...and the real document is untouched: still one specimen, still its
        # run association.
        doc = store.get(pid, refresh_capability=False)
        assert [s["id"] for s in doc["specimens"]] == ["s_00001"]
        assert doc["specimens"][0]["run_ids"] == ["r_00001"]
        # ...and nothing named after a hostile id was created under the root.
        assert sorted(os.listdir(store.root)) == [pid]

    def test_the_regex_is_gone_and_did_not_come_back_unwired(self):
        """An unused guard that looks like protection is its own hazard: the
        next reader assumes it is doing something. If a ``SPECIMEN_ID_RE`` is
        ever reintroduced it must be APPLIED, not merely defined."""
        assert not hasattr(projects_mod, "SPECIMEN_ID_RE"), (
            "SPECIMEN_ID_RE is back. If it is needed it must be wired to the "
            "boundary that needs it -- and note ^s_\\d{5}$ is exactly-five "
            "digits, which would brick the 100,000th specimen for the reason "
            "RUN_ID_RE documents with \\d{5,}.")
        with open(os.path.join(_REPO, "app", "server", "projects.py")) as fh:
            src = fh.read()
        for name in ("PROJECT_ID_RE", "RUN_ID_RE"):
            # Every id regex this module defines is APPLIED somewhere in it:
            # one definition, at least one call site, ignoring prose mentions.
            applied = re.findall(rf"{name}\s*\.\s*(?:match|search|findall)\(",
                                 src)
            assert applied, f"{name} is defined and never applied"


class TestASpecimenIdIsNotAPathOverHTTP:
    """The same question at the boundary that actually took the traversal.

    ``httpx`` normalises ``..`` out of a request path CLIENT-side, so a
    ``TestClient`` call is not evidence about this at all — it turns
    ``/api/projects/p_00001/specimens/..`` into ``DELETE /api/projects/p_00001``
    and deletes the project, which looks exactly like a traversal and is not
    one. The literal segment therefore has to be pushed through a raw ASGI
    scope.
    """

    @staticmethod
    def _raw(app, path, method="GET", body=b""):
        """Drive the ASGI app with a literal, un-normalised ``raw_path``.

        ``body`` is sent with a JSON content type because PATCH parses its body
        BEFORE it looks the specimen up: with an empty body it answers 400 and
        never reaches the id at all, which would be a PATCH row that tests
        nothing.
        """
        import asyncio

        got = []

        async def go():
            headers = [(b"host", b"testserver")]
            if body:
                headers += [(b"content-type", b"application/json"),
                            (b"content-length", str(len(body)).encode())]
            scope = {"type": "http", "asgi": {"version": "3.0"},
                     "http_version": "1.1", "method": method, "scheme": "http",
                     "path": path, "raw_path": path.encode(), "root_path": "",
                     "query_string": b"", "headers": headers,
                     "client": ("testclient", 123),
                     "server": ("testserver", 80)}

            async def receive():
                return {"type": "http.request", "body": body,
                        "more_body": False}

            async def send(msg):
                got.append(msg)

            await app(scope, receive, send)

        asyncio.run(go())
        status = next((m["status"] for m in got
                       if m["type"] == "http.response.start"), None)
        body = b"".join(m.get("body", b"") for m in got
                        if m["type"] == "http.response.body")
        return status, body.decode("utf-8", "replace")

    def test_a_literal_dotdot_in_the_specimen_slot_is_a_404(self, tmp_path):
        from app.server.api import create_app

        store = ProjectStore(root=str(tmp_path / "projects"),
                             runs_root=str(tmp_path / "runs"))
        pid = _seed_project(store, "S", [
            {"id": "s_00001", "label": "A", "section": SECTION_A,
             "run_ids": [], "capability": {}}])
        doc_path = os.path.join(store.root, pid, "project.json")

        manager = JobManager(root=str(tmp_path / "runs"), start_worker=False)
        app = create_app(manager, store)
        try:
            # The premise: the same route with a REAL id answers 200, so a 404
            # below is the id being rejected and not the route being absent.
            ok, _ = self._raw(app, f"/api/projects/{pid}/specimens/s_00001")
            assert ok == 200, "the positive control did not route; test is void"

            patch_body = json.dumps({"label": "x"}).encode()
            for bad in ("..", "%2e%2e", "%2fetc%2fpasswd", "s_00001%00",
                        "....//", "__proto__"):
                for method in ("GET", "PATCH", "DELETE"):
                    status, body = self._raw(
                        app, f"/api/projects/{pid}/specimens/{bad}", method,
                        patch_body if method == "PATCH" else b"")
                    assert status in (404, 405), \
                        f"{method} {bad!r} -> {status}: {body[:200]}"
                assert os.path.exists(doc_path), \
                    f"{bad!r} deleted the project document"
        finally:
            manager.shutdown()

        assert [s["id"] for s in store.get(
            pid, refresh_capability=False)["specimens"]] == ["s_00001"]
