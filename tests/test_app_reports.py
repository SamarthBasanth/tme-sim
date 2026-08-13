"""Therapy reports — ``app/server/report.py`` and the ``/api/reports`` routes.

Three things are under test and they are not the same thing:

1. **The science is shared, not forked.** The attribution the Studio runs must be
   the *identical function object* ``experiments/patient_report.py`` runs.
2. **The limits are carried, prominently.** Every non-negotiable caveat this
   repository has measured has to be in the document, with its number, above the
   ranking rather than under it — and the NOT FOR CLINICAL USE banner has to be
   at the top AND the bottom.
3. **The document is self-contained.** A rendered report must make no network
   request of any kind: it has to read the same off a USB stick as it does in
   the product.

The end-to-end arm runs the real engine on the real cached section, small and
coarse, and asserts the ranking is assembled from the two real arms per drug.
"""
from __future__ import annotations

import importlib.util
import os
import re
import sys
import time

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

pytest.importorskip("fastapi", reason="the Studio backend is an optional extra")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi.testclient import TestClient  # noqa: E402

from app.server import cache as cache_mod  # noqa: E402
from app.server.api import create_app, validate_report_request  # noqa: E402
from app.server.jobs import JobManager  # noqa: E402
from app.server.report import (  # noqa: E402
    DEFAULT_PANEL,
    NOT_FOR_CLINICAL_USE,
    build_limits,
    render_html,
)

SECTION = "xenium_319"
DX_UM = 30.0
WINDOW_S = 600.0
#: Two drugs is the smallest panel that HAS a ranking and a panel mean.
DRUGS = ["temozolomide", "cisplatin"]


def _section_available(section_id: str = SECTION) -> bool:
    return any(s["id"] == section_id and s["available"]
               for s in cache_mod.list_sections())


requires_section = pytest.mark.skipif(
    not _section_available(),
    reason=f"section {SECTION!r} is not cached in data_cache/",
)


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    """A TestClient whose runs AND reports are isolated from the real dirs."""
    root = str(tmp_path_factory.mktemp("study"))
    manager = JobManager(root=os.path.join(root, "runs"))
    app = create_app(manager)
    with TestClient(app) as c:
        c._manager = manager  # type: ignore[attr-defined]
        yield c
    manager.shutdown()


def _await_report(client, report_id, timeout=1800.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = client.get(f"/api/reports/{report_id}").json()
        if body["status"] in ("done", "failed", "cancelled"):
            return body
        time.sleep(0.5)
    raise AssertionError(f"report {report_id} did not finish within {timeout}s")


# --------------------------------------------------------------------------- #
# 1. The science is imported, not re-implemented
# --------------------------------------------------------------------------- #
def _load_patient_report():
    path = os.path.join(_ROOT, "experiments", "patient_report.py")
    spec = importlib.util.spec_from_file_location("_patient_report", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestSharedAttribution:
    def test_cli_and_studio_call_the_same_function_object(self):
        """Not "equivalent code" — the SAME object. Two implementations of an
        attribution would be two sciences."""
        pytest.importorskip("matplotlib",
                            reason="experiments/patient_report.py imports viz")
        from tmesim.report import attribution

        cli = _load_patient_report()
        assert cli.attribute is attribution.attribute
        assert cli._gap_reason is attribution.gap_reason
        assert cli._TERMS is attribution.TERMS
        assert cli._TERM_LABEL is attribution.TERM_LABEL

    def test_delivery_plus_potency_is_the_drive_exactly(self):
        """The decomposition's own identity, on numbers it did not choose."""
        import math

        from tmesim.report.attribution import attribute

        rows = [
            {"name": "a", "tumor_auc": 1.69e6, "ic50": 24.6, "kill": 0.37,
             "kill_no_resistance": 0.74, "L_um": 160.0},
            {"name": "b", "tumor_auc": 1.08e4, "ic50": 0.5, "kill": 0.22,
             "kill_no_resistance": 0.45, "L_um": 73.0},
            {"name": "c", "tumor_auc": 7.06e5, "ic50": 606.0, "kill": 0.015,
             "kill_no_resistance": 0.030, "L_um": 296.0},
        ]
        attribute(rows, window_s=3600.0 * 6)
        log_drive = [math.log10(r["drive"]) for r in rows]
        mean = sum(log_drive) / len(log_drive)
        for r, ld in zip(rows, log_drive):
            assert r["delivery_decades"] + r["potency_decades"] == \
                pytest.approx(ld - mean, abs=1e-12)


# --------------------------------------------------------------------------- #
# 2. The request contract
# --------------------------------------------------------------------------- #
class TestReportRequest:
    def test_defaults_to_the_five_drug_panel(self):
        v = validate_report_request({"section": SECTION})
        assert v["ok"], v["errors"]
        assert v["drugs"] == list(DEFAULT_PANEL)
        assert v["base_spec"]["drug"] is None

    def test_sensitivity_mode_cannot_be_overridden(self):
        """'determinant' is not a ranking on a panel that measures some drugs'
        determinants and not others; a report must refuse to rank in it."""
        v = validate_report_request({
            "section": SECTION,
            "spec_overrides": {"sensitivity_mode": "determinant", "dx_um": 20.0},
        })
        assert v["ok"], v["errors"]
        assert v["base_spec"]["sensitivity_mode"] == "uniform"
        assert v["base_spec"]["dx_um"] == 20.0     # a NON-owned field passes through
        assert any("sensitivity_mode" in w and "IGNORED" in w
                   for w in v["warnings"])

    def test_unknown_drug_and_section_are_actionable(self):
        v = validate_report_request({"section": SECTION, "drugs": ["aspirin"]})
        assert not v["ok"]
        assert any("aspirin" in e and "Available" in e for e in v["errors"])
        v = validate_report_request({"section": "nope", "drugs": DRUGS})
        assert not v["ok"]
        assert any("nope" in e for e in v["errors"])

    def test_a_one_drug_panel_is_warned_about(self):
        v = validate_report_request({"section": SECTION, "drugs": ["cisplatin"]})
        assert v["ok"]
        assert any("CROSS-drug" in w or "no ranking" in w for w in v["warnings"])

    def test_the_spec_the_arms_will_run_is_validated_up_front(self):
        """A report must never enqueue 2N runs /api/validate would refuse."""
        v = validate_report_request({
            "section": SECTION, "drugs": DRUGS,
            "spec_overrides": {"dx_um": -1.0},
        })
        assert not v["ok"]
        assert any("dx_um" in e for e in v["errors"])

    def test_post_rejects_an_invalid_request_without_enqueueing(self, client):
        before = len(client.get("/api/runs").json())
        r = client.post("/api/reports",
                        json={"section": SECTION, "drugs": ["aspirin"]})
        assert r.status_code == 400
        assert len(client.get("/api/runs").json()) == before
        assert client.get("/api/reports/rep_99999").status_code == 404


# --------------------------------------------------------------------------- #
# 3. The limits block — the non-negotiable content, with its numbers
# --------------------------------------------------------------------------- #
class TestLimits:
    @pytest.fixture(scope="class")
    def limits(self):
        spec = {"dx_um": 12.0, "oxygen_mode": "pde", "state_mode": "oxygen"}
        return build_limits(SECTION, spec, DEFAULT_PANEL)

    def test_every_mandatory_limit_is_present(self, limits):
        ids = [x["id"] for x in limits]
        for required in ("not_for_clinical_use", "ranking_is_the_output",
                         "cell_state_support", "kill_is_an_interval",
                         "uniform_upper_bound", "oxygen_pde_falsified",
                         "no_outcome_label"):
            assert required in ids, f"missing limit {required!r}"

    def test_the_banner_is_the_first_limit(self, limits):
        assert limits[0]["id"] == "not_for_clinical_use"
        assert "NOT FOR CLINICAL USE" in limits[0]["body"]

    def test_the_ranking_is_named_as_the_defensible_output(self, limits):
        lim = next(x for x in limits if x["id"] == "ranking_is_the_output")
        assert "0.3252" in lim["body"] and "1.77x" in lim["body"]
        assert "cell_state_impact.md" in lim["source"]

    @requires_section
    def test_the_319_panel_is_disqualified_from_an_absolute_kill(self, limits):
        """results/cell_state_impact.md: no hypoxia gene, 1.42% proliferation."""
        lim = next(x for x in limits if x["id"] == "cell_state_support")
        assert lim["applies"] is True
        assert lim["severity"] == "critical"
        assert "1.42%" in lim["body"]
        assert "19-40%" in lim["body"]
        assert "UNMEASURED" in lim["body"]
        assert "VEGFA" in lim["body"]   # the hypoxia transcripts it lacks

    def test_absolute_kill_is_an_interval_with_the_measured_span(self, limits):
        lim = next(x for x in limits if x["id"] == "kill_is_an_interval")
        assert "0.0048-0.1442" in lim["body"]
        assert "30x" in lim["body"]
        assert "kill_calibration.md" in lim["source"]

    def test_uniform_is_declared_an_upper_bound(self, limits):
        lim = next(x for x in limits if x["id"] == "uniform_upper_bound")
        assert "UPPER BOUND" in lim["body"]
        assert "MGMT" in lim["body"]

    def test_the_oxygen_field_is_declared_weak_in_the_corrected_direction(self, limits):
        """The oxygen limit is PINNED, in the corrected direction.

        This assertion used to require ``"FALSIFIED"`` in the title. That claim
        was an inverted inequality — the solve out-predicts the baselines it was
        said to lose to (results/oxygen_demotion.md). The limit is kept and still
        pinned, because a reversal that turns into an overclaim in the opposite
        direction is the same error with its sign flipped: what must reach the
        reader is that the winner is WEAK, and that the unfitted-``consumption``
        criticism survives the reversal untouched.
        """
        lim = next(x for x in limits if x["id"] == "oxygen_pde_falsified")
        assert "oxygen_demotion.md" in lim["source"]
        assert "WEAK" in lim["title"]
        body = lim["body"]
        # The corrected ordering, with the numbers it rests on.
        assert "out-predicting" in body
        assert "-0.1523" in body and "-0.1204" in body
        # The limit on it — this must never read as a validation.
        assert "best of several weak predictors" in body
        assert "never 'the oxygen solve is validated'" in body
        # The criticism that survives the reversal.
        assert "never been fitted" in body
        assert "invariant to monotone rescaling" in body
        # The withdrawn claim survives as history only, never as an assertion.
        assert "has been withdrawn" in body
        assert "FALSIFIED" not in lim["title"]

    def test_no_outcome_label_anywhere(self, limits):
        lim = next(x for x in limits if x["id"] == "no_outcome_label")
        assert "UNKNOWN" in lim["body"]
        assert "no outcome" in lim["body"].lower()

    def test_the_grid_convergence_claim_is_retracted_not_asserted(self, limits):
        """The resolution limit is PINNED, in the corrected direction.

        It used to assert, at the renderer's LOWEST severity ("info"), that
        "kill_fraction is grid-converged for dx <= 20 um", citing
        results/numerical_rigor.md §2 — a claim that file's own correction block
        withdraws and that results/grid_convergence_breast.md and
        results/grid_convergence_gbm.md both contradict. It also had the sign of
        the penetration trend backwards. This was the only place in the product
        where a caveat surface was less honest than the documentation, so the
        limit is kept and pinned rather than deleted: what must reach the reader
        is that NO dx has been shown to be converged, and that the coarse-grid
        error runs the other way.
        """
        lim = next(x for x in limits if x["id"] == "resolution")
        assert lim["applies"] is True
        # never reassurance again
        assert lim["severity"] != "info"
        assert "WITHDRAWN" in lim["title"]
        body = lim["body"]
        assert "NO dx AT WHICH THIS PIPELINE HAS BEEN SHOWN TO BE " \
               "GRID-CONVERGED" in body.replace("THERE IS ", "")
        # the retracted sentence survives as history, never as an assertion
        assert "BOTH HALVES ARE WITHDRAWN" in body
        # the two ladders that contradict it, with their measured numbers
        assert "0.1484" in body and "0.1537" in body
        assert "FAIL -> PASS" in body
        assert "NON-MONOTONE" in body
        assert "0.96" in body and "-3.33" in body
        # the corrected direction of the penetration trend
        assert "0.1110" in body and "0.3062" in body
        assert "UNDER-states" in body
        # ...and the overcorrection guard: the solver itself is not impeached
        assert "4e-14" in body
        for src in ("numerical_rigor.md", "grid_convergence_breast.md",
                    "grid_convergence_gbm.md"):
            assert src in lim["source"]


# --------------------------------------------------------------------------- #
# 4. The rendered document
# --------------------------------------------------------------------------- #
_EXTERNAL = re.compile(
    r"""(?ix)
    https?://           |   # any absolute URL
    \burl\s*\(          |   # a CSS resource reference
    @import             |   # a stylesheet import
    <img\b              |   # an image
    <script\b           |   # a script (inline or otherwise)
    <link\b             |   # a stylesheet / preload link
    <iframe\b           |   # a nested browsing context
    \bsrcset\b          |
    \bfetch\s*\(
    """
)


def _doc(status="done"):
    """A minimal but structurally complete report document."""
    return {
        "report_id": "rep_00001", "status": status, "section": SECTION,
        "drugs": DRUGS, "run_ids": ["r_1", "r_2"], "window_s": 3600.0 * 6,
        "sensitivity_mode": "uniform", "generator": "app/server/report.py",
        "attribution_implementation": "tmesim.report.attribution.attribute",
        "spec_base": {"dx_um": 12.0, "vascular_mode": "discrete",
                      "oxygen_mode": "pde", "state_mode": "oxygen"},
        "section_detail": {"label": "Xenium GBM", "panel": "hBrain_254g",
                           "has_outcome_label": False},
        "microenvironment": {"n_cells": 40887, "n_tumor": 30174,
                             "mean_ecm_density": 0.076, "grid": [424, 675]},
        "limits": build_limits(SECTION, {"dx_um": 12.0}, DRUGS),
        "caveats": ["a caveat carried by an arm"],
        "warnings": [],
        "attribution": {"method": "delivery_i = ...",
                        "check": {"spearman_total_vs_kill": 1.0,
                                  "faithful": "faithful",
                                  "spearman_delivery_vs_kill": 1.0,
                                  "spearman_potency_vs_kill": -1.0}},
        "ranking": [
            {"rank": 1, "name": "cisplatin", "kill": 0.3699,
             "kill_no_resistance": 0.7377, "rho": 0.5014, "hypoxia_cost": 0.4986,
             "tumor_auc": 1.69e6, "c_tumor_uM": 78.2, "ic50": 24.6,
             "L_um": 160.0, "d_eff": 225.0, "bbb": 0.858, "pen": 0.42,
             "coverage_1x": 0.9955,
             "mw": 298.0, "diffusivity": 236.0, "binding": 0.0087,
             "clearance": 7.6e-05, "hill": 1.2,
             "delivery_decades": 0.98, "potency_decades": -0.06,
             "hypoxia_decades": -0.01, "total_decades": 0.91,
             "exposure_pct": 100.0, "potency_pct": 50.0,
             "dominant": "DELIVERY (advantage)",
             "gap_reason": "Versus **temozolomide** immediately below it ...",
             "provenance": {
                 "run_id": "r_1", "counterfactual_run_id": "r_2",
                 "alpha": 148.6, "alpha_source": "calibrate_alpha()",
                 "ic50": {"ic50_uM": 24.6, "confidence": "measured",
                          "source": "GDSC measured (n=34 glioma cell lines)",
                          "n_cell_lines": 34, "note": "median across 34 lines",
                          "provenance_kind": "measured"},
                 "sensitivity_mode": "uniform", "sensitivity_upper_bound": True,
                 "sensitivity_mean": 1.0, "sensitivity_determinant_gene": None,
                 "mgmt": {"state": "unmodelled - MGMT not on panel",
                          "detail": "MGMT is not on this section's panel."},
                 "seed": 0, "state_mode": "oxygen", "caveats": []}},
        ],
    }


class TestRenderedDocument:
    @pytest.fixture(scope="class")
    def html(self):
        return render_html(_doc())

    def test_it_makes_no_external_request_of_any_kind(self, html):
        """Self-contained is a hard requirement: the file has to read the same
        off a USB stick, from an attachment, or inside an air-gapped review."""
        hits = _EXTERNAL.findall(html)
        assert not hits, f"rendered report references something external: {hits}"

    def test_it_inlines_the_studio_design_tokens(self, html):
        """Product-styled, and from the product's OWN tokens.css — not a
        lookalike palette that will drift the first time a token changes."""
        tokens = open(os.path.join(_ROOT, "app", "web", "css",
                                   "tokens.css")).read()
        marker = "--accent:      hsl(191 72% 54%);"
        assert marker in tokens          # the file still says what we think
        assert marker in html            # ...and the document carries it
        assert "--font-mono" in html and "tabular-nums" in html

    def test_the_banner_is_at_the_top_and_the_bottom(self, html):
        wheres = re.findall(r'class="nfcu"[^>]*data-where="(\w+)"', html)
        assert wheres == ["top", "bottom"], wheres
        assert html.count(NOT_FOR_CLINICAL_USE.split(".")[0]) >= 3
        assert "NOT FOR CLINICAL USE" in html[:html.index("</title>")]

    def test_the_limits_come_before_the_ranking(self, html):
        """Prominent, not a footnote."""
        limits_at = html.index("Limits — read these before any number")
        ranking_at = html.index("Ranked panel")
        first_kill = html.index("0.3699")
        assert limits_at < ranking_at < first_kill

    def test_every_mandatory_limit_number_survives_rendering(self, html):
        for needle in ("0.0048-0.1442", "30x", "0.3252", "1.77x", "1.42%",
                       "19-40%", "UPPER BOUND", "FALSIFIED"):
            assert needle in html, f"{needle!r} was lost in rendering"

    def test_per_drug_provenance_is_rendered(self, html):
        assert "GDSC measured (n=34 glioma cell lines)" in html
        assert "unmodelled - MGMT not on panel" in html
        assert "148.6" in html                      # alpha
        assert "mode=uniform" in html
        assert "upper_bound=True" in html

    def test_the_counterfactual_is_named_as_the_source_of_rho(self, html):
        assert "COUNTERFACTUAL" in html
        assert "r_2" in html                        # the resistance-off arm

    def test_the_hero_carries_the_condition_the_same_page_states(self, html):
        """The hero names the top of a ranking, which is the defensible output.
        It may NOT name it as settled: limit 2b on the same page says the
        ordering is law-dependent, and the delivery law is precisely what
        overturns this pick (cisplatin -> vincristine)."""
        hero = html[html.index('<div class="ch">'):html.index("2 &middot;")] \
            if "2 &middot;" in html else html
        # the ranking headline is still there — it is the claim
        assert "Top of this ranking" in html
        assert "cisplatin" in html
        # ...and it now states the law it is conditional on, in the hero itself
        assert "S = PS*Cp" in html or "S = PS*(Cp - C)" in html
        assert "LAW-DEPENDENT" in html
        assert "vincristine" in html
        assert "SET, not a\n     winner" in html or "SET, not a winner" in html
        # the pre-existing honest lines are NOT weakened
        assert "This is a model output, not a recommendation" in html
        assert "ORDER is the claim" in html
        assert "UPPER BOUND" in hero

    def test_the_hero_names_the_law_stored_with_the_document(self):
        """Not a constant: the law that produced THIS document."""
        doc = _doc()
        doc["delivery_law"] = {
            "law": "S = PS*(Cp - C)", "alternative_law": "S = PS*Cp",
            "title": "DELIVERY LAW S = PS*(Cp - C) (exchange)",
            "headline": "x", "caveat": "y", "reference": "z",
        }
        html = render_html(doc)
        assert "Top of this ranking, under S = PS*(Cp - C)" in " ".join(
            html.split())

    def test_the_retired_metric_lost_its_column_to_the_computed_one(self, html):
        """`penetration_index` is retired (results/delivery_metrics.md §4) — a
        scale-invariant shape statistic that reports the section's vascular
        geometry, shared by every drug. `therapeutic_coverage_1x` was already
        computed and stored for every row and was never rendered. Swapped."""
        head = html[html.index("<th class='num'>predicted kill</th>"):]
        head = head[:head.index("</thead>")]
        assert "coverage @1x IC50" in head
        # the retired metric is gone as a COLUMN; it survives only in the
        # footnote that explains that it was retired, which is history
        assert "pen. index" not in head
        assert "0.9955" in html                 # the row's coverage, rendered
        # the reader is told why, on the page, with the citation
        assert "delivery_metrics.md" in html
        assert "scale-invariant" in html

    def test_a_report_that_is_not_done_still_renders_its_limits(self):
        html = render_html(_doc(status="running"))
        assert "0.0048-0.1442" in html
        assert not _EXTERNAL.findall(html)


# --------------------------------------------------------------------------- #
# 5. End to end on the real engine
# --------------------------------------------------------------------------- #
@requires_section
class TestReportEndToEnd:
    @pytest.fixture(scope="class")
    def report(self, client):
        posted = client.post("/api/reports", json={
            "section": SECTION, "drugs": DRUGS,
            "spec_overrides": {"dx_um": DX_UM, "total_time_s": WINDOW_S},
        })
        assert posted.status_code == 202, posted.text
        body = posted.json()
        assert len(body["run_ids"]) == 2 * len(DRUGS)
        done = _await_report(client, body["report_id"])
        assert done["status"] == "done", done.get("error")
        return done

    def test_it_enqueues_two_arms_per_drug_on_the_shared_queue(self, client, report):
        """Not a private simulation path: every arm is a normal run."""
        runs = {r["run_id"]: r for r in client.get("/api/runs").json()}
        for drug in DRUGS:
            arms = report["arms"][drug]
            real, cf = runs[arms["real"]], runs[arms["counterfactual"]]
            assert real["spec"]["drug"] == drug
            assert real["spec"]["resistance_off"] is False
            assert cf["spec"]["resistance_off"] is True
            # Identical in every other respect — that is what makes rho a ratio.
            assert {k: v for k, v in real["spec"].items() if k != "resistance_off"} \
                == {k: v for k, v in cf["spec"].items() if k != "resistance_off"}
            assert cf["provenance"]["pd_hypoxia_resistance"] == 1.0
            assert cf["provenance"]["pd_quiescent_resistance"] == 1.0

    def test_the_reconstruction_is_built_once_and_shared(self, client, report):
        """The whole reason a report goes THROUGH the queue and not around it."""
        runs = {r["run_id"]: r for r in client.get("/api/runs").json()}
        hits = [runs[rid]["report"]["reconstruction_cache_hit"]
                for rid in report["run_ids"]]
        assert hits.count(False) <= 1, "a report rebuilt the reconstruction"
        assert hits.count(True) >= len(report["run_ids"]) - 1

    def test_the_counterfactual_never_kills_less_than_the_real_arm(self, report):
        """Switching resistance OFF cannot lower the kill. If it does, rho is
        not a resistance ratio and the attribution term is meaningless."""
        for row in report["ranking"]:
            assert row["kill_no_resistance"] >= row["kill"] - 1e-12
            assert 0.0 < row["rho"] <= 1.0 + 1e-9
            assert row["rho"] == pytest.approx(
                row["kill"] / row["kill_no_resistance"], rel=1e-9)

    def test_the_ranking_is_sorted_by_predicted_kill(self, report):
        kills = [r["kill"] for r in report["ranking"]]
        assert kills == sorted(kills, reverse=True)
        assert [r["rank"] for r in report["ranking"]] == \
            list(range(1, len(DRUGS) + 1))
        assert report["top_ranked"]["drug"] == report["ranking"][0]["name"]

    def test_the_attribution_reproduces_the_simulator_ranking(self, report):
        check = report["attribution"]["check"]
        assert check["spearman_total_vs_kill"] == pytest.approx(1.0)
        assert check["faithful"] == "faithful"
        assert report["attribution"]["implementation"] == \
            "tmesim.report.attribution.attribute"

    def test_every_row_carries_its_ic50_provenance_and_mgmt_status(self, report):
        for row in report["ranking"]:
            prov = row["provenance"]
            assert prov["ic50"]["source"]
            assert prov["ic50"]["confidence"] in (
                "measured", "literature", "default")
            assert prov["sensitivity_mode"] == "uniform"
            assert prov["sensitivity_upper_bound"] is True
            assert prov["mgmt"]["state"]
            assert prov["alpha"] and prov["alpha_source"]
        tmz = next(r for r in report["ranking"] if r["name"] == "temozolomide")
        # MGMT is the determinant for TMZ and is on neither Xenium panel.
        assert tmz["provenance"]["mgmt"]["state"] == \
            "unmodelled - MGMT not on panel"

    def test_the_tumour_auc_comes_from_the_arm_and_is_not_recomputed(
            self, client, report):
        runs = {r["run_id"]: r for r in client.get("/api/runs").json()}
        for row in report["ranking"]:
            arm = runs[row["provenance"]["run_id"]]
            assert row["tumor_auc"] == arm["metrics"]["mean_tumor_auc"]
            assert row["kill"] == arm["metrics"]["kill_fraction"]
            assert row["ic50"] == arm["provenance"]["ic50_uM"]

    def test_the_html_route_serves_a_self_contained_document(self, client, report):
        r = client.get(f"/api/reports/{report['report_id']}/html")
        assert r.status_code == 200
        assert r.headers["content-type"].startswith("text/html")
        assert not _EXTERNAL.findall(r.text)
        assert "NOT FOR CLINICAL USE" in r.text
        assert report["ranking"][0]["name"] in r.text

    def test_it_is_listed_and_persisted(self, client, report):
        rows = client.get("/api/reports").json()
        row = next(x for x in rows if x["report_id"] == report["report_id"])
        assert row["status"] == "done"
        assert row["n_runs"] == 2 * len(DRUGS)
        assert row["not_for_clinical_use"] is True
        path = os.path.join(client._manager.root, "..", "reports",
                            report["report_id"], "report.html")
        assert os.path.exists(os.path.abspath(path))

    def test_the_limits_travel_with_the_json_too(self, report):
        ids = [x["id"] for x in report["limits"]]
        assert "kill_is_an_interval" in ids
        assert "ranking_is_the_output" in ids
        assert report["not_for_clinical_use"] is True
