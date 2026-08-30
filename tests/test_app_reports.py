"""The single-compound map document — ``app/server/report.py`` and ``/api/reports``.

Three things are under test and they are not the same thing:

1. **The attribution exists at N = 1, and it is absolute.** The report's own
   decomposition must be referenced to quantities from the SAME run, not to a
   panel mean. The panel-centred function in ``tmesim.report.attribution`` is
   identically zero on one compound; it is untouched, still used by
   ``experiments/patient_report.py``, and deliberately NOT what the Studio calls.
2. **The limits are carried, prominently.** Every non-negotiable caveat this
   repository has measured has to be in the document, with its number, above the
   map rather than under it — and the NOT FOR CLINICAL USE banner has to be at
   the top AND the bottom.
3. **The document is self-contained.** A rendered report must make no network
   request of any kind: it has to read the same off a USB stick as it does in
   the product.

The end-to-end arm runs the real engine on the real cached section, small and
coarse, and asserts the document is assembled from the two real arms of ONE
compound.

WHAT CHANGED ON 2026-08-21, BECAUSE THESE TESTS USED TO ASSERT THE OPPOSITE.
This file used to pin a ranked five-compound panel, a ``rank`` per row, a sort by
predicted kill, and a Spearman check that the attribution "reproduces the
simulator ranking". The cross-compound ranking was withdrawn:
``results/pk_ranking_robustness.md`` sec. 0 measures per-compound plasma exposure
spanning 113.1x across the compounds the GBM panel ranked while the order breaks
at 1.53x, and ``results/capacity_ranking_robustness.md`` sec. 0 changes the
breast top-1 under this repo's own structural capacity rule. No assertion was
deleted in this conversion; each was re-pointed at the single-compound contract
that replaced the one it pinned.
"""
from __future__ import annotations

import importlib.util
import math
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
    MAP_CLAIM,
    NOT_FOR_CLINICAL_USE,
    TERMS,
    attribute_compound,
    build_limits,
    render_html,
)

SECTION = "xenium_319"
DX_UM = 30.0
WINDOW_S = 600.0
#: ONE compound is the only legal shape. It used to be a panel of two, described
#: here as "the smallest panel that HAS a ranking and a panel mean" — the very
#: property the withdrawn decomposition depended on and the new one does not.
DRUG = "cisplatin"


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
# 1. The decomposition: absolute, single-compound, and NOT the panel-centred one
# --------------------------------------------------------------------------- #
def _load_patient_report():
    path = os.path.join(_ROOT, "experiments", "patient_report.py")
    spec = importlib.util.spec_from_file_location("_patient_report", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestAttribution:
    def test_the_studio_no_longer_shares_the_panel_centred_attribution(self):
        """This test used to assert that the Studio and the CLI call the SAME
        function object, on the reasoning that two implementations of an
        attribution would be two sciences. That reasoning is intact and the CLI
        half of it is still pinned below — what changed is that the shared
        function is PANEL-CENTRED and therefore identically zero at N = 1, so a
        single-compound document cannot use it and must not pretend to.
        """
        pytest.importorskip("matplotlib",
                            reason="experiments/patient_report.py imports viz")
        from tmesim.report import attribution

        from app.server import report as studio

        # The CLI's four bindings, unchanged: the panel-centred science stayed
        # exactly where it was and this pass did not fork or edit it.
        cli = _load_patient_report()
        assert cli.attribute is attribution.attribute
        assert cli._gap_reason is attribution.gap_reason
        assert cli._TERMS is attribution.TERMS
        assert cli._TERM_LABEL is attribution.TERM_LABEL

        # ...and the Studio's is a DIFFERENT, absolute function, by design.
        assert studio.attribute_compound is not attribution.attribute
        assert studio.TERMS != attribution.TERMS
        src = open(os.path.join(_ROOT, "app", "server", "report.py")).read()
        assert "tmesim.report.attribution.attribute" not in src.split(
            "ATTRIBUTION_RATIONALE")[0].split("Why the attribution")[0], \
            "report.py must not call the panel-centred attribution"

        # The reason is stated in the module, not only here.
        assert "identically zero" in studio.ATTRIBUTION_RATIONALE
        assert "PANEL-CENTRED" in studio.ATTRIBUTION_RATIONALE

    def test_delivery_plus_potency_is_the_drive_exactly(self):
        """The decomposition's own identity, on numbers it did not choose.

        Re-pointed from the panel-centred form
        (``delivery_i + potency_i == log10 drive_i - mean_j log10 drive_j``) to
        the absolute one. The identity survives the change of decomposition; only
        the right-hand side loses its panel mean.
        """
        rows = [
            {"tumor_auc": 1.69e6, "ic50": 24.6, "c_plasma_uM": 100.0,
             "kill": 0.37, "kill_no_resistance": 0.74, "L_um": 160.0},
            {"tumor_auc": 1.08e4, "ic50": 0.5, "c_plasma_uM": 3.2,
             "kill": 0.22, "kill_no_resistance": 0.45, "L_um": 73.0},
            {"tumor_auc": 7.06e5, "ic50": 606.0, "c_plasma_uM": 880.0,
             "kill": 0.015, "kill_no_resistance": 0.030, "L_um": 296.0},
        ]
        for r in rows:
            block = attribute_compound(r, window_s=3600.0 * 6)
            assert block["null_reason"] is None
            assert r["delivery_decades"] + r["dose_decades"] == \
                pytest.approx(math.log10(r["drive"]), abs=1e-12)

    def test_the_exact_identity_holds_over_a_wide_numeric_sweep(self):
        """DELIVERY + DOSE == log10(drive), swept rather than spot-checked.

        Measured worst deviation over this sweep: 1.78e-15, i.e. floating-point
        round-off and nothing else.
        """
        import random

        random.seed(20260821)
        worst = 0.0
        for _ in range(2000):
            row = {
                "tumor_auc": 10 ** random.uniform(2, 8),
                "ic50": 10 ** random.uniform(-3, 3),
                "c_plasma_uM": 10 ** random.uniform(-3, 3),
                "kill": 0.3, "kill_no_resistance": 0.6,
            }
            T = 10 ** random.uniform(1, 6)
            assert attribute_compound(row, window_s=T)["null_reason"] is None
            worst = max(worst, abs(row["delivery_decades"] + row["dose_decades"]
                                   - math.log10(row["drive"])))
        assert worst < 1e-12, worst

    def test_delivery_does_not_move_when_the_plasma_level_is_rescaled(self):
        """DELIVERY is INVARIANT to the dose, and that is the point of it.

        ``C_tumor`` is exactly linear in the plasma source — worst relative
        deviation 1.41e-14 over plasma factors 0.005 to 26, a 5,200x span, at the
        shipping configuration (results/pk_ranking_robustness.md sec. 1.2) — so
        rescaling the dose multiplies ``C_tumor`` and ``C_plasma`` together and
        the factor cancels from ``log10(C_tumor / C_plasma)``.

        RESISTANCE does not move either (it is a ratio of two kills of the same
        compound at the same dose). The whole dose uncertainty therefore lands in
        DOSE, which is exactly the quarantine the decomposition exists to
        provide — so this test also asserts that DOSE DOES move, by the full
        log10 of the factor. A DOSE term that did not move would mean the
        uncertainty had gone somewhere unnamed.
        """
        base = {"tumor_auc": 1.69e6, "ic50": 24.6, "c_plasma_uM": 100.0,
                "kill": 0.3699, "kill_no_resistance": 0.7377}
        attribute_compound(base, window_s=3600.0 * 6)

        for factor in (0.005, 0.033, 0.66, 1.0, 26.0, 5200.0):
            scaled = dict(base, tumor_auc=1.69e6 * factor,
                          c_plasma_uM=100.0 * factor)
            attribute_compound(scaled, window_s=3600.0 * 6)
            assert scaled["delivery_decades"] == pytest.approx(
                base["delivery_decades"], abs=1e-12), factor
            assert scaled["resistance_decades"] == pytest.approx(
                base["resistance_decades"], abs=1e-12), factor
            # ...and the uncertainty is in DOSE, by exactly log10(factor).
            assert scaled["dose_decades"] - base["dose_decades"] == \
                pytest.approx(math.log10(factor), abs=1e-12), factor

    def test_one_compound_produces_NON_ZERO_terms(self):
        """The regression the panel-centred form would fail BY CONSTRUCTION.

        ``delivery_i = log10 C_i - mean_j log10 C_j`` over a single row is
        exactly zero, and so are the other two terms. Any future change that
        reintroduces a panel mean fails here.
        """
        row = {"tumor_auc": 1.69e6, "ic50": 24.6, "c_plasma_uM": 100.0,
               "kill": 0.3699, "kill_no_resistance": 0.7377}
        block = attribute_compound(row, window_s=3600.0 * 6)
        assert block["null_reason"] is None
        for key in TERMS:
            assert row[key] is not None, key
            assert math.isfinite(row[key]), key
            assert abs(row[key]) > 1e-6, \
                f"{key} is zero on a single compound — a panel mean is back"
        # And each carries its documented sign on this row.
        assert row["delivery_decades"] < 0      # exposure below plasma
        assert row["dose_decades"] > 0          # plasma above IC50
        assert row["resistance_decades"] < 0    # resistance costs kill

    def test_a_missing_plasma_level_REFUSES_rather_than_substituting(self):
        """C_plasma is read from the run and never invented."""
        for bad in (None, float("nan"), 0.0, -1.0):
            row = {"tumor_auc": 1.69e6, "ic50": 24.6, "c_plasma_uM": bad,
                   "kill": 0.37, "kill_no_resistance": 0.74}
            block = attribute_compound(row, window_s=3600.0 * 6)
            assert block["null_reason"], bad
            assert "plasma_mean_uM" in block["null_reason"]
            assert "NOT substituted" in block["null_reason"]
            for key in TERMS:
                assert row[key] is None, (bad, key)

    def test_a_missing_counterfactual_refuses_only_the_resistance_term(self):
        """DELIVERY and DOSE read no kill, so they survive a missing rho."""
        row = {"tumor_auc": 1.69e6, "ic50": 24.6, "c_plasma_uM": 100.0,
               "kill": 0.37, "kill_no_resistance": 0.0}
        block = attribute_compound(row, window_s=3600.0 * 6)
        assert block["null_reason"] is None
        assert block["resistance_null_reason"]
        assert "never assumed" in block["resistance_null_reason"]
        assert row["resistance_decades"] is None
        assert row["delivery_decades"] is not None
        assert row["dose_decades"] is not None


# --------------------------------------------------------------------------- #
# 2. The request contract
# --------------------------------------------------------------------------- #
class TestReportRequest:
    def test_a_report_is_one_compound_and_a_panel_is_refused(self):
        """Converted from ``test_defaults_to_the_five_drug_panel``.

        That test asserted ``v["drugs"] == list(DEFAULT_PANEL)`` — a report with
        no compound named defaulted to ranking five. There is no default panel
        and no panel at all: a request naming none is an error, and a request
        naming ``drugs[]`` is refused rather than truncated to its first element.
        The one assertion that survives verbatim is that ``base_spec.drug`` is
        None, because the base spec is still the shared, drug-free RunSpec the
        manager clones per arm.
        """
        v = validate_report_request({"section": SECTION})
        assert not v["ok"]
        assert any("drug" in e and "required" in e for e in v["errors"])

        v = validate_report_request({"section": SECTION, "drug": DRUG})
        assert v["ok"], v["errors"]
        assert v["drug"] == DRUG
        assert v["base_spec"]["drug"] is None
        assert v["map_claim"] == MAP_CLAIM

        # A panel is a 400 with the reason, never a silent truncation.
        v = validate_report_request({"section": SECTION,
                                     "drugs": ["cisplatin", "temozolomide"]})
        assert not v["ok"]
        err = " ".join(v["errors"])
        assert "113.1x" in err and "1.53x" in err
        assert "pk_ranking_robustness.md" in err
        assert "capacity_ranking_robustness.md" in err
        assert "not truncated" in err

    def test_sensitivity_mode_cannot_be_overridden(self):
        """'determinant' would move with measurement coverage rather than with
        pharmacology, so a report must refuse to run in it.

        The assertion is unchanged; only the REASON in the warning moved, from
        "the only mode in which a cross-drug ranking is fair" to the one that
        outlives the ranking.
        """
        v = validate_report_request({
            "section": SECTION, "drug": DRUG,
            "spec_overrides": {"sensitivity_mode": "determinant", "dx_um": 20.0},
        })
        assert v["ok"], v["errors"]
        assert v["base_spec"]["sensitivity_mode"] == "uniform"
        assert v["base_spec"]["dx_um"] == 20.0     # a NON-owned field passes through
        assert any("sensitivity_mode" in w and "IGNORED" in w
                   for w in v["warnings"])
        assert any("UPPER BOUND" in w and "measurement coverage" in w
                   for w in v["warnings"])

    def test_unknown_drug_and_section_are_actionable(self):
        v = validate_report_request({"section": SECTION, "drug": "aspirin"})
        assert not v["ok"]
        assert any("aspirin" in e and "Available" in e for e in v["errors"])
        v = validate_report_request({"section": "nope", "drug": DRUG})
        assert not v["ok"]
        assert any("nope" in e for e in v["errors"])

    def test_one_compound_is_the_ONLY_legal_shape(self):
        """Inverted from ``test_a_one_drug_panel_is_warned_about``.

        One compound used to earn a warning that "a one-drug panel has no
        ranking and no panel mean, so every attribution term is 0 by
        construction". Both halves are now true of the PRODUCT rather than of a
        degenerate request: there is no ranking, and the decomposition no longer
        has a panel mean to be zero against. One compound is the contract, so it
        warns about nothing.
        """
        v = validate_report_request({"section": SECTION, "drug": "cisplatin"})
        assert v["ok"], v["errors"]
        assert v["drug"] == "cisplatin"
        assert not any("CROSS-drug" in w or "no ranking" in w
                       for w in v["warnings"]), v["warnings"]

    def test_the_spec_the_arms_will_run_is_validated_up_front(self):
        """A report must never enqueue runs /api/validate would refuse."""
        v = validate_report_request({
            "section": SECTION, "drug": DRUG,
            "spec_overrides": {"dx_um": -1.0},
        })
        assert not v["ok"]
        assert any("dx_um" in e for e in v["errors"])

    def test_post_rejects_an_invalid_request_without_enqueueing(self, client):
        before = len(client.get("/api/runs").json())
        r = client.post("/api/reports",
                        json={"section": SECTION, "drug": "aspirin"})
        assert r.status_code == 400
        assert len(client.get("/api/runs").json()) == before
        # ...and so does a panel, which is the shape this endpoint retired.
        r = client.post("/api/reports",
                        json={"section": SECTION, "drugs": [DRUG]})
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
        return build_limits(SECTION, spec, DRUG)

    def test_every_mandatory_limit_is_present(self, limits):
        ids = [x["id"] for x in limits]
        for required in ("not_for_clinical_use", "map_is_the_output",
                         "transvascular_delivery_law",
                         "diffusivity_two_route_partial",
                         "cell_state_support", "kill_is_an_interval",
                         "uniform_upper_bound", "oxygen_pde_falsified",
                         "no_outcome_label", "one_anchor",
                         "pk_simplification"):
            assert required in ids, f"missing limit {required!r}"
        # The withdrawn id must not linger beside its replacement.
        assert "ranking_is_the_output" not in ids

    def test_the_banner_is_the_first_limit(self, limits):
        assert limits[0]["id"] == "not_for_clinical_use"
        assert "NOT FOR CLINICAL USE" in limits[0]["body"]

    def test_the_map_is_named_as_the_defensible_output(self, limits):
        """Converted from ``test_the_ranking_is_named_as_the_defensible_output``.

        The two pinned numbers (0.3252 absolute, 1.77x relative) and the
        cell_state_impact.md citation are UNCHANGED: they were always about the
        absolute kill's dependence on the cell-state route, which withdrawing the
        ranking does not touch. What is added is the claim that replaced the
        ordering, and the two invariances at their true and DIFFERENT strengths.
        """
        lim = next(x for x in limits if x["id"] == "map_is_the_output")
        # unchanged, still true, still about absolute kill
        assert "0.3252" in lim["body"] and "1.77x" in lim["body"]
        assert "cell_state_impact.md" in lim["source"]

        body = lim["body"]
        # the claim, verbatim
        assert MAP_CLAIM in body
        # the withdrawal, with the numbers that forced it
        assert "WITHDRAWN" in body
        assert "113.1x" in body and "1.53x" in body
        assert "1307.5x" in body
        # The ORDERING ITSELF IS NOT PRINTED. This assertion used to be
        # `assert "vinorelbine -> doxorubicin" in body`, a leftover the
        # conversion named in this docstring missed: the whole point of the
        # withdrawal is that this document does not put two compound names in
        # an order, and a limits paragraph that does it "as evidence" still
        # prints an ordering to a reader who is skimming. The FORCE of the
        # caveat is kept -- the top-1 changes, and the file that names the pair
        # is cited -- which is what the next two assertions check.
        assert "vinorelbine" not in body.lower()
        assert "top-1 CHANGES" in body
        assert "capacity_ranking_robustness.md" in body
        # invariance 1 — EXACT, with its measurement
        assert "1.41e-14" in body
        assert "0.005 to 26" in body and "5,200x" in body
        # invariance 2 — exact ONLY for a converged shape, never unconditional
        assert "CONVERGED GRADIENT SHAPE ONLY" in body
        assert "1e-9" in body
        assert "AUC_R(x,T) = R * AUC_1(x,T/R)" in body
        assert "residual" in body
        # what neither protects against
        assert "SHAPE of the gradient" in body
        assert "descriptor_falsification.md" in body
        # ...and that the two are a LIST, explicitly
        assert "LIST OF THE AXES THAT WERE SWEPT" in body
        assert "UNTESTED, not covered" in body
        for src in ("pk_ranking_robustness.md", "capacity_ranking_robustness.md",
                    "descriptor_falsification.md"):
            assert src in lim["source"]

    def test_the_delivery_law_limit_survived_and_now_qualifies_the_MAP(self, limits):
        """The re-parented disclosure. Its title used to be suffixed "— the
        ranking above is LAW-DEPENDENT". The limit is about the SOURCE TERM, not
        about the ranking, so it must survive the ranking with a title that
        describes what it now qualifies: the SHAPE of the map."""
        lim = next(x for x in limits
                   if x["id"] == "transvascular_delivery_law")
        assert lim["severity"] == "critical"
        assert "LAW-DEPENDENT" in lim["title"]
        assert "MAP" in lim["title"]
        assert "ranking above" not in lim["title"]
        assert "boundary condition" in lim["body"]
        assert "relative exposure between two regions" in lim["body"]

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
        """The assertions are unchanged. The RATIONALE inside the limit is what
        moved: it used to justify ``uniform`` as "the only setting in which a
        cross-drug ranking is fair", which is an argument about a ranking. The
        reason that outlives it is pinned here alongside the original two."""
        lim = next(x for x in limits if x["id"] == "uniform_upper_bound")
        assert "UPPER BOUND" in lim["body"]
        assert "MGMT" in lim["body"]
        # the reason that survives the withdrawal
        assert "no-resistance UPPER BOUND is the honest reading" in lim["body"]
        assert "an estimate would not be" in lim["body"]
        assert "MEASUREMENT-COVERAGE" in lim["body"]
        # ...and the superseded reason survives as history, not as an argument
        assert "THE REASON IT IS FORCED HAS CHANGED" in lim["body"]

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
        # Re-pointed: the sentence that used to say "the ranking does not rest on
        # the oxygen field" now names the TWO things on the page that DO read it.
        assert "EXPOSURE MAP" in body
        assert "RESISTANCE decade" in body
        assert "'viable' tumour-cell band" in body

    def test_no_outcome_label_anywhere(self, limits):
        lim = next(x for x in limits if x["id"] == "no_outcome_label")
        assert "UNKNOWN" in lim["body"]
        assert "no outcome" in lim["body"].lower()
        # Re-pointed from "the ranking below has never been scored" — and the
        # limit says explicitly that it did not weaken in the move.
        assert "THE MAP BELOW HAS THEREFORE NEVER BEEN SCORED" in lim["body"]
        assert "did not weaken this limit by one word" in lim["body"]

    def test_the_single_anchor_is_named_as_a_MAP_SHAPE_risk(self, limits):
        """Re-pointed: alpha sets how steeply exposure falls away from a vessel,
        so it re-shapes the map rather than scaling it — which is why it is not
        divided out by either invariance in the map limit."""
        lim = next(x for x in limits if x["id"] == "one_anchor")
        assert "doxorubicin" in lim["body"] and "mitoxantrone" in lim["body"]
        assert "MAP-SHAPE RISK" in lim["body"]
        assert "do not divide it out" in lim["body"]

    def test_the_pk_simplification_names_the_invariance_it_buys(self, limits):
        """Re-pointed from "it preserves total delivered drug and the relative
        ranking". There is no ranking; what the single plasma level actually
        preserves is the RATIO between two points of the map, because it enters
        as one multiplicative factor on the source term."""
        lim = next(x for x in limits if x["id"] == "pk_simplification")
        assert "ONE MULTIPLICATIVE FACTOR" in lim["body"]
        assert "1.41e-14" in lim["body"]
        assert "no schedule optimisation" in lim["body"]

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


def _compound_row():
    """One compound's row, with its decades computed by the real function.

    Built through :func:`attribute_compound` rather than hand-written, so the
    fixture cannot drift out of its own identity.
    """
    row = {
        "name": DRUG, "kill": 0.3699, "kill_no_resistance": 0.7377,
        "tumor_auc": 1.69e6, "ic50": 24.6, "c_plasma_uM": 100.0,
        "plasma_peak_uM": 140.0,
        "L_um": 160.0, "d_eff": 225.0, "bbb": 0.858, "pen": 0.42,
        "coverage_0_5x": 0.9991, "coverage_1x": 0.9955, "coverage_2x": 0.9702,
        "median_tumor_auc": 1.51e6, "mean_auc": 1.2e6,
        "E_p10_uM": 21.0, "E_p50_uM": 70.0, "E_p90_uM": 132.0,
        "viable_tumor_fraction": 0.88, "hypoxic_tumor_fraction": 0.12,
        "penetration_length_fit_um": 151.0,
        "mw": 298.0, "diffusivity": 236.0, "binding": 0.0087,
        "clearance": 7.6e-05, "hill": 1.2,
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
            "seed": 0, "state_mode": "oxygen", "caveats": []},
    }
    block = attribute_compound(row, window_s=3600.0 * 6)
    return row, block


def _doc(status="done"):
    """A minimal but structurally complete report document."""
    row, attribution = _compound_row()
    return {
        "report_id": "rep_00001", "status": status, "section": SECTION,
        "drug": DRUG, "run_ids": ["r_1", "r_2"], "window_s": 3600.0 * 6,
        "sensitivity_mode": "uniform", "generator": "app/server/report.py",
        "attribution_implementation": "app.server.report.attribute_compound",
        "map_claim": MAP_CLAIM,
        "spec_base": {"dx_um": 12.0, "vascular_mode": "discrete",
                      "oxygen_mode": "pde", "state_mode": "oxygen"},
        "section_detail": {"label": "Xenium GBM", "panel": "hBrain_254g",
                           "has_outcome_label": False},
        "microenvironment": {"n_cells": 40887, "n_tumor": 30174,
                             "mean_ecm_density": 0.076, "grid": [424, 675]},
        "limits": build_limits(SECTION, {"dx_um": 12.0}, DRUG),
        "caveats": ["a caveat carried by an arm"],
        "warnings": [],
        "attribution": attribution,
        "compound": row,
        "mechanism": {
            "null_reason": None,
            "unreached_fraction_of_tumor": 0.4123,
            "unreached_fraction_of_survivors": 0.6544,
            "survival_breakdown": {"delivery": 0.6544, "dose": 0.2011,
                                   "resistance": 0.1201, "chance": 0.0244},
            "counts": {"SURVIVED_UNREACHED": 12440}, "n_tumor": 30174,
            "n_survived": 19010, "killed_fraction_of_tumor": 0.3699,
            "limiter_vocabulary": {"delivery_decades": "DELIVERY"},
            "definition": "SURVIVED_UNREACHED is a TUMOR cell that survived ...",
        },
        "reach": {
            "null_reason": None,
            "headline_under_reached_fraction": 0.2871, "headline_x": 0.1,
            "curve": {"x": [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0],
                      "under_reached_fraction": [0.02, 0.05, 0.16, 0.2871,
                                                 0.44, 0.71, 0.90]},
            "falloff": {"distance_um": [5.0, 15.0, 25.0, 35.0],
                        "normalised_mean_auc": [1.0, 0.71, 0.48, 0.31],
                        "count": [911, 1420, 1188, 743],
                        "bin_um": 10.0, "min_bin_count": 30},
            "E_ref_uM_s": 2.44e6,
            "E_ref_source": "delivery_metrics.penetration_length_fit(...) "
                            "['drug']['perivascular_mean'] — REUSED verbatim",
            "n_viable": 26553, "n_tumor": 30174,
            "viable_definition": "a TUMOR cell whose local oxygen ...",
            "convention": "x = 0.1 is a REPORTING CONVENTION, not a physical "
                          "constant.",
            "scale_invariance": "PLASMA/DOSE SCALE: exact ... 1.41e-14 ... "
                                "TISSUE CAPACITY: exact for a CONVERGED "
                                "gradient shape ... 1e-9.",
        },
        "reach_reason": None,
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

    def test_the_limits_come_before_the_map(self, html):
        """Prominent, not a footnote. Converted from
        ``test_the_limits_come_before_the_ranking``: the middle landmark was the
        "Ranked panel" heading and is now the single-compound reading, and the
        first number is still the run's own kill."""
        limits_at = html.index("Limits — read these before any number")
        map_at = html.index("The single-compound reading")
        first_kill = html.index("0.3699")
        assert limits_at < map_at < first_kill

    def test_every_mandatory_limit_number_survives_rendering(self, html):
        for needle in ("0.0048-0.1442", "30x", "0.3252", "1.77x", "1.42%",
                       "19-40%", "UPPER BOUND", "FALSIFIED",
                       # the numbers the map claim now rests on
                       "113.1x", "1.53x", "1.41e-14", "1307.5x", "1e-9"):
            assert needle in html, f"{needle!r} was lost in rendering"

    def test_the_map_claim_is_on_the_page_verbatim(self, html):
        """The one sentence the whole document exists to support."""
        assert MAP_CLAIM in html
        assert html.count(MAP_CLAIM) >= 2      # the limit AND the hero

    def test_per_drug_provenance_is_rendered(self, html):
        assert "GDSC measured (n=34 glioma cell lines)" in html
        assert "unmodelled - MGMT not on panel" in html
        assert "148.6" in html                      # alpha
        assert "mode=uniform" in html
        assert "upper_bound=True" in html

    def test_the_provenance_chip_legend_still_renders(self, html):
        """It was rendered beside the ranking table and must not have gone with
        it: a colour a reader has to infer carries no information."""
        assert "a number somebody observed in an experiment" in html
        assert "a number this simulation produced" in html
        assert "not measured on this tissue" in html

    def test_the_counterfactual_is_named_as_the_source_of_rho(self, html):
        assert "COUNTERFACTUAL" in html
        assert "r_2" in html                        # the resistance-off arm

    def test_the_hero_leads_with_the_per_cell_unreached_fraction(self, html):
        """Converted from ``test_the_hero_carries_the_condition_the_same_page_states``.

        The hero used to name the top of a ranking and to carry the condition
        that ordering was true under. It now leads with SURVIVED_UNREACHED — the
        strongest thing in the product — and carries the conditions the MAP is
        true under. Every honest line the old hero carried is still asserted
        here; none was dropped in the conversion.
        """
        hero = html[html.index('<div class="hero">'):]
        hero = hero[:hero.index("</section>")]

        # the headline IS the per-cell answer, not a supporting figure
        assert "41.2% of this tumour's cells were never reached" in hero
        assert "SURVIVED_UNREACHED" in hero
        assert DRUG in hero
        # the claim, in the hero, verbatim
        assert MAP_CLAIM in hero
        # the pre-existing honest lines are NOT weakened
        assert "This is a model output, not a recommendation" in hero
        assert "UPPER BOUND" in hero
        assert "LAW-DEPENDENT" in hero
        # ...and the law it is conditional on is still stated IN the hero
        assert "S = PS*Cp" in hero or "S = PS*(Cp - C)" in hero
        # the two invariances, at their true and DIFFERENT strengths
        assert "1.41e-14" in hero and "5,200x" in hero
        assert "CONVERGED GRADIENT SHAPE ONLY" in hero
        assert "1e-9" in hero
        assert "UNTESTED, not covered" in hero
        # The two hero disclosures whose SUBJECT was the ranking are re-pointed,
        # not dropped: the law's leverage is still quoted with the swap that
        # measured it, and the missing potency readout is re-aimed at the terms
        # that actually read an IC50.
        assert "vincristine" in hero
        assert "47x" in hero
        assert "EVIDENCE OF THE LAW'S LEVERAGE" in hero
        assert "SET, not a winner" in hero
        assert "0.587" in hero
        assert "ex-vivo potency readout" in hero
        assert "not the exposure map, which reads no potency at all" in hero
        # the condition is not a footnote: both the law limit and the map limit
        # it cites are reachable from the hero on this page
        assert 'href="#lim-transvascular_delivery_law"' in hero
        assert 'id="lim-transvascular_delivery_law"' in html
        assert 'href="#lim-map_is_the_output"' in hero
        assert 'id="lim-map_is_the_output"' in html

    def test_the_hero_names_the_law_stored_with_the_document(self):
        """Not a constant: the law that produced THIS document."""
        doc = _doc()
        doc["delivery_law"] = {
            "law": "S = PS*(Cp - C)", "alternative_law": "S = PS*Cp",
            "title": "DELIVERY LAW S = PS*(Cp - C) (exchange)",
            "headline": "x", "caveat": "y", "reference": "z",
        }
        html = render_html(doc)
        hero = html[html.index('<div class="hero">'):]
        hero = hero[:hero.index("</section>")]
        flat = " ".join(hero.split())
        assert "produced under <span class=\"mono\">S = PS*(Cp - C)</span>" in flat
        assert "S = PS*Cp" in flat        # the alternative it is contrasted with

    def test_the_coverage_ladder_replaced_the_single_quoted_rung(self, html):
        """Converted from
        ``test_the_retired_metric_lost_its_column_to_the_computed_one``.

        That test pinned the swap of ``pen. index`` for ``therapeutic_coverage``
        at 1x. Both halves survive and both are extended. The retired metric is
        still on the page only as a labelled TISSUE-GEOMETRY statistic with its
        citation, never as a delivery number; and 1x IC50 is no longer quoted
        alone — all three rungs the run already computed are rendered, with the
        note that coverage is the ABSOLUTE quantity on the page.
        """
        assert "0.9955" in html                 # the 1x rung, still rendered
        assert "0.9991" in html and "0.9702" in html   # ...and the ladder
        for rung in ("0.5x IC50", "1x IC50", "2x IC50"):
            assert rung in html
        # the reader is told why the index was retired, on the page, with its
        # citation — a methodology disclosure, not a ranking artefact
        assert "delivery_metrics.md" in html
        assert "scale-invariant" in html
        assert "TISSUE-GEOMETRY statistic" in html
        assert "never as a delivery number" in html
        # ...and that coverage is the one thing scaling the field would scale
        assert "Coverage is an ABSOLUTE quantity" in html

    def test_the_reach_curve_is_rendered_in_full_not_just_its_headline(self, html):
        """The threshold must read as the convention it is, so the whole ladder
        is printed beside the point, and the falloff with distance from vessel
        with it."""
        for x in ("0.01", "0.02", "0.05", "0.1", "0.2", "0.5"):
            assert x in html
        assert "0.2871" in html                    # F(0.1), the headline
        assert "REPORTING CONVENTION" in html
        # the exposure falloff with distance from the nearest vessel
        assert "Exposure falloff with distance from the nearest vessel" in html
        assert "0.7100" in html and "0.4800" in html
        # and the scale-invariance sentence travels with it
        assert "TISSUE CAPACITY: exact for a CONVERGED" in html

    def test_a_missing_reach_block_degrades_with_a_stated_reason(self):
        """Absent is not zero. A section that says nothing is
        indistinguishable from a section whose answer was zero."""
        doc = _doc()
        doc["reach"] = None
        doc["reach_reason"] = "this run predates app/server/reach.py"
        html = render_html(doc)
        assert "The scale-free under-reached curve is NOT available" in html
        assert "this run predates app/server/reach.py" in html
        assert not _EXTERNAL.findall(html)

    def test_the_three_decades_are_rendered_with_what_they_mean(self, html):
        assert "DELIVERY" in html and "DOSE" in html and "RESISTANCE" in html
        assert "decades lost blood -&gt; interstitium" in html or \
               "decades LOST getting from blood" in html
        assert "DELIVERY + DOSE" in html
        assert "log10(drive)" in html
        # the rationale for the replacement is ON the page, not only in code
        assert "identically zero" in html
        assert "quarantine" in html.lower()

    def test_no_rank_survives_anywhere_in_the_document(self, html):
        """The ranking is withdrawn, so no rank ordinal may leak back in."""
        assert "Top of this ranking" not in html
        assert "Ranked panel" not in html
        assert "Rank gap" not in html
        assert "percentile of this panel" not in html
        assert "class='num rank'" not in html

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
            "section": SECTION, "drug": DRUG,
            "spec_overrides": {"dx_um": DX_UM, "total_time_s": WINDOW_S},
        })
        assert posted.status_code == 202, posted.text
        body = posted.json()
        assert len(body["run_ids"]) == 2
        done = _await_report(client, body["report_id"])
        assert done["status"] == "done", done.get("error")
        return done

    def test_it_enqueues_two_arms_of_one_compound_on_the_shared_queue(
            self, client, report):
        """Not a private simulation path: both arms are normal runs.

        Converted from ``test_it_enqueues_two_arms_per_drug_...``: still two
        arms, now of the one compound. Every assertion is unchanged; only the
        loop over a panel is gone.
        """
        runs = {r["run_id"]: r for r in client.get("/api/runs").json()}
        arms = report["arms"]
        assert set(arms) == {"real", "counterfactual"}
        real, cf = runs[arms["real"]], runs[arms["counterfactual"]]
        assert real["spec"]["drug"] == DRUG
        assert cf["spec"]["drug"] == DRUG
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
        not a resistance ratio and the RESISTANCE decade is meaningless."""
        row = report["compound"]
        assert row["kill_no_resistance"] >= row["kill"] - 1e-12
        assert 0.0 < row["rho"] <= 1.0 + 1e-9
        assert row["rho"] == pytest.approx(
            row["kill"] / row["kill_no_resistance"], rel=1e-9)

    def test_the_document_is_about_exactly_one_compound(self, report):
        """Converted from ``test_the_ranking_is_sorted_by_predicted_kill``.

        That test asserted the rows were sorted by kill, that ``rank`` ran
        1..N, and that ``top_ranked.drug`` matched row 0. There is no order, no
        rank and no ``top_ranked``. What replaces all three is the single-subject
        invariant: the document names one compound, the compound block IS that
        compound, and the summary headline agrees with it — including that the
        summary carries no kill and no rank.
        """
        assert report["drug"] == DRUG
        assert report["compound"]["name"] == DRUG
        assert "ranking" not in report
        assert "top_ranked" not in report
        assert report["compound"].get("rank") is None
        assert report["headline"]["drug"] == DRUG
        assert "kill" not in report["headline"]
        assert report["map_claim"] == MAP_CLAIM

    def test_the_attribution_is_absolute_and_non_zero_on_real_numbers(self, report):
        """Converted from ``test_the_attribution_reproduces_the_simulator_ranking``.

        That test read a Spearman rank correlation between the net attribution
        and the simulated kill across the panel and required it to be 1.0. There
        is no panel and therefore no rank correlation. What replaces it is the
        property the new decomposition actually has and the old one did not: the
        terms are non-zero at N = 1 and the exact identity holds on numbers the
        engine produced rather than on numbers a test chose.
        """
        row = report["compound"]
        att = report["attribution"]
        assert att["implementation"] == "app.server.report.attribute_compound"
        assert att["null_reason"] is None
        assert report["attribution_implementation"] == \
            "app.server.report.attribute_compound"

        for key in TERMS:
            assert row[key] is not None, key
            assert math.isfinite(row[key]), key
            assert abs(row[key]) > 1e-9, f"{key} is zero — a panel mean is back"

        assert row["delivery_decades"] + row["dose_decades"] == \
            pytest.approx(math.log10(row["drive"]), abs=1e-12)
        # DELIVERY is a ratio to the run's OWN plasma mean, not to a panel.
        assert row["delivery_decades"] == pytest.approx(
            math.log10(row["c_tumor_uM"] / row["c_plasma_uM"]), abs=1e-12)
        assert row["resistance_decades"] == pytest.approx(
            math.log10(row["rho"]), abs=1e-12)

    def test_the_compound_carries_its_ic50_provenance_and_mgmt_status(self, report):
        prov = report["compound"]["provenance"]
        assert prov["ic50"]["source"]
        assert prov["ic50"]["confidence"] in ("measured", "literature", "default")
        assert prov["sensitivity_mode"] == "uniform"
        assert prov["sensitivity_upper_bound"] is True
        assert prov["mgmt"]["state"]
        assert prov["alpha"] and prov["alpha_source"]

    def test_the_tumour_auc_and_the_plasma_mean_come_from_the_arm(
            self, client, report):
        """Neither is recomputed here. ``plasma_mean_uM`` in particular is READ
        from the run's extra_metrics, because it is the reference of the DELIVERY
        term and inventing it would invent the term."""
        runs = {r["run_id"]: r for r in client.get("/api/runs").json()}
        row = report["compound"]
        arm = runs[row["provenance"]["run_id"]]
        assert row["tumor_auc"] == arm["metrics"]["mean_tumor_auc"]
        assert row["kill"] == arm["metrics"]["kill_fraction"]
        assert row["ic50"] == arm["provenance"]["ic50_uM"]
        assert row["c_plasma_uM"] == arm["extra_metrics"]["plasma_mean_uM"]
        assert row["coverage_1x"] == arm["metrics"]["therapeutic_coverage_1x"]

    def test_the_reach_block_is_consumed_from_the_run_not_recomputed(
            self, client, report):
        """app/server/reach.py owns this curve; the document consumes it."""
        runs = {r["run_id"]: r for r in client.get("/api/runs").json()}
        arm = runs[report["compound"]["provenance"]["run_id"]]
        if arm.get("reach") is None:
            assert report["reach_reason"], \
                "a run with no reach block must produce a stated reason"
            return
        assert report["reach"] == arm["reach"]

    def test_the_per_cell_unreached_fraction_is_the_headline(self, report):
        mech = report["mechanism"]
        if mech.get("null_reason"):
            assert report["headline"]["unreached_null_reason"]
            return
        frac = mech["unreached_fraction_of_tumor"]
        assert 0.0 <= frac <= 1.0
        assert report["headline"]["unreached_fraction_of_tumor"] == frac
        assert mech["counts"]["SURVIVED_UNREACHED"] / mech["n_tumor"] == \
            pytest.approx(frac, rel=1e-9)

    def test_the_html_route_serves_a_self_contained_document(self, client, report):
        r = client.get(f"/api/reports/{report['report_id']}/html")
        assert r.status_code == 200
        assert r.headers["content-type"].startswith("text/html")
        assert not _EXTERNAL.findall(r.text)
        assert "NOT FOR CLINICAL USE" in r.text
        assert report["compound"]["name"] in r.text
        assert MAP_CLAIM in r.text

    def test_it_is_listed_and_persisted(self, client, report):
        rows = client.get("/api/reports").json()
        row = next(x for x in rows if x["report_id"] == report["report_id"])
        assert row["status"] == "done"
        assert row["n_runs"] == 2
        assert row["not_for_clinical_use"] is True
        assert row["drug"] == DRUG
        # The delivery law travels with the SUMMARY, re-parented off the ranking
        # it used to qualify: a list row is read without opening the document.
        assert row["delivery_law"]["law"]
        assert row["map_claim"] == MAP_CLAIM
        path = os.path.join(client._manager.root, "..", "reports",
                            report["report_id"], "report.html")
        assert os.path.exists(os.path.abspath(path))

    def test_the_limits_travel_with_the_json_too(self, report):
        ids = [x["id"] for x in report["limits"]]
        assert "kill_is_an_interval" in ids
        assert "map_is_the_output" in ids
        assert "transvascular_delivery_law" in ids
        assert "ranking_is_the_output" not in ids
        assert report["not_for_clinical_use"] is True
