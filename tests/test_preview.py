"""The fast preview is a COARSE SOLVE, and these are the ways it can lie.

Four things make ``app/server/preview.py`` honest rather than fast, and each
one is a planted defect away from being gone. Every class below was confirmed
to FAIL against that defect before it was kept:

1. the preview is refused on a spot section, where ``dx_um`` is clamped to
   100 um and the "preview" would be the identical solve compared with itself;
2. the preview is refused for the FITTED PENETRATION LENGTH, which disagrees by
   up to 16.0405 %, and for the PENETRATION INDEX, whose dx agreement is
   unreconciled across three of this repository's own records;
3. the MEASURED agreement travels with the number, and is recomputed from the
   measured pairs rather than copied into source;
4. the number is never rendered without the word PREVIEW beside it, at no
   smaller than its own size.

The browser class at the bottom executes the REAL ``app/web`` page in Chromium.
A page that returns 200 is not a page that booted, and a panel that parses is
not a panel that painted.
"""
from __future__ import annotations

import functools
import http.server
import json
import os
import re
import socketserver
import threading

import pytest

from app.server.preview import (
    AGREEMENT_SOURCE,
    MEASURED_COMPOUNDS,
    MEASURED_SECTIONS,
    PREVIEW_DX_UM,
    PREVIEW_LABEL,
    PREVIEWABLE_READOUT,
    REFUSED_READOUTS,
    SPOT_PREVIEW_REFUSAL,
    TARGET_DX_UM,
    PreviewError,
    _audit_payload,
    measured_agreement,
    preview_payload,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(_ROOT, "app", "web")

BASE_SPEC = {
    "section": "xenium_319",
    "dx_um": 12.0,
    "drug": "temozolomide",
    "sensitivity_mode": "uniform",
    "vascular_mode": "discrete",
    "oxygen_mode": "pde",
    "state_mode": "oxygen",
    "mgmt": None,
    "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0, "n_doses": 1,
           "dose_interval_s": 86400.0},
    "total_time_s": 21600.0,
}


class _Jobs:
    """Stands in for JobManager. The preview only ever reads ``.cache``."""


def _refused(**over):
    spec = dict(BASE_SPEC)
    spec.update(over)
    with pytest.raises(PreviewError) as ei:
        preview_payload(_Jobs(), spec)
    return ei.value


# --------------------------------------------------------------------------- #
# 1. The agreement is a measurement, not a literal
# --------------------------------------------------------------------------- #
class TestTheAgreementIsDerivedFromTheMeasuredPairs:
    """The error bar is recomputed from ``coarse_baseline.json`` on every
    request. A second copy of a measured number in source is this
    repository's recurring defect (``results/codon_claim_surface_drift``)."""

    def test_the_source_file_is_read_and_not_mirrored(self):
        ag = measured_agreement()
        assert ag is not None, (
            f"{AGREEMENT_SOURCE} did not resolve; the preview cannot be served "
            "without it and must refuse rather than fall back to a literal")
        assert ag["source"] == AGREEMENT_SOURCE

    def test_only_pairs_whose_grid_actually_differs_are_counted(self):
        """18 of the 30 stored pairs compared a run WITH ITSELF.

        On a Visium section ``_reconstruct_visium`` clamps dx_um below 50 um to
        100 um, so dx = 30 and dx = 12 produce the identical grid. Counting
        those is how the study behind this feature first reported a perfect
        0.000 % result on all four readouts.
        """
        ag = measured_agreement()
        assert ag["n_pairs"] == 12, (
            f"expected the 12 genuine pairs, got {ag['n_pairs']}. If this is "
            "30 the tautologies are being counted and the error bar is a "
            "measurement of nothing.")
        assert ag["pairs_discarded_same_grid"] == 18
        assert ag["n_sections"] == 2 and ag["n_compounds"] == 6

    def test_the_figures_match_the_published_record(self):
        """``results/surrogate_feasibility_2026-10-09.md`` sec. 6, to 4 dp.

        Including the PERCENTILE CONVENTION: nearest-rank, the same as
        ``experiments/surrogate_feasibility/derived.py``. Linear interpolation
        over 12 values puts the p90 at 1.1757 % and would put this feature's
        error bar at odds with the record it cites.
        """
        ag = measured_agreement()
        assert round(ag["median_rel_err_pct"], 4) == 0.7133
        assert round(ag["p90_rel_err_pct"], 4) == 1.2001
        assert round(ag["max_rel_err_pct"], 4) == 1.7620
        assert round(ag["max_abs_err"], 4) == 0.0909
        assert round(ag["cost_ratio"], 3) == 6.780
        assert round(ag["coarse_median_s"], 2) == 8.67
        assert round(ag["target_median_s"], 2) == 58.77

    def test_the_error_bar_is_small_against_a_named_comparator(self):
        """0.0909 against a 0.5357 across-donor range is 17 %, and against the
        readout's own 0.202-wide nuisance band is 45 %. Both comparators are
        carried with the figure so the reader does not have to find them."""
        ag = measured_agreement()
        assert ag["across_donor_range"] == 0.5357238964376265
        assert ag["nuisance_band_width"] == 0.202
        assert "results/" in ag["across_donor_range_source"]
        assert "results/" in ag["nuisance_band_source"]
        assert ag["max_abs_err"] / ag["across_donor_range"] < 0.20

    def test_no_comparator_is_an_ordering_between_compounds(self):
        """``results/delivery_spread_noise.md`` also reports a compound
        SEPARATION. Cross-compound ranking was withdrawn 2026-08-21, so the
        within-compound span is quoted from it and the separation is not, and
        the compound list is a SET."""
        ag = measured_agreement()
        assert ag["compounds"] == sorted(ag["compounds"]), (
            "the compound list is not alphabetical, so its order carries "
            "information it must not carry")
        # The two `*_source` strings legitimately contain the word "ordering",
        # because what they say is that there ISN'T one. So the scan is over
        # the DATA, not over the prose that disclaims it.
        data = {k: v for k, v in ag.items() if not k.endswith("_source")}
        blob = json.dumps(data).lower()
        for word in ("rank", "best", "better", "worse", "top-1", "winner"):
            assert word not in blob, (
                f"{word!r} in the agreement data; the compound list is a set "
                "and no ordering between compounds is available here")
        # And no per-compound figure at all: a number beside each compound is a
        # ranking however it got there.
        assert all(isinstance(c, str) for c in ag["compounds"])


# --------------------------------------------------------------------------- #
# 2. The refusals
# --------------------------------------------------------------------------- #
class TestSpotSectionsAreRefused:
    """Not because the statistics are weak there -- because ``dx_um`` is not a
    lever at all, so the preview is a NO-OP that would report 0.000 %."""

    def test_a_visium_section_is_refused(self):
        exc = _refused(section="visium_gbm2")
        assert exc.status == 422
        assert exc.message == SPOT_PREVIEW_REFUSAL

    def test_the_refusal_names_the_clamp_that_makes_it_true(self):
        """A refusal defended by a false reason is one relaxation away from
        being dropped, which is the lesson SPOT_SPATIAL_REFUSAL records about
        itself. So this one names the line of code."""
        assert "_reconstruct_visium" in SPOT_PREVIEW_REFUSAL
        assert ">= 50.0 else 100.0" in SPOT_PREVIEW_REFUSAL
        assert "100 um" in SPOT_PREVIEW_REFUSAL

    def test_the_clamp_the_refusal_cites_is_really_in_the_pipeline(self):
        """The refusal quotes a line of ``real_pipeline.py``. If that line
        changes, the preview's reason for refusing Visium has changed and this
        test is the thing that says so."""
        src = os.path.join(_ROOT, "tmesim", "reconstruct", "real_pipeline.py")
        with open(src) as fh:
            text = fh.read()
        assert "visium_dx = float(dx_um) if dx_um and dx_um >= 50.0 else 100.0" in text, (
            "the Visium dx clamp this module's refusal is built on is no "
            "longer in real_pipeline.py. Re-measure before relaxing the "
            "refusal: if dx_um has become a lever on a spot section, a preview "
            "there may be possible, and if it has become a different clamp the "
            "refusal's wording is now wrong.")

    def test_a_spot_section_is_refused_before_anything_is_solved(self):
        """The refusal is cheap. A preview that reconstructed first and then
        refused would cost the thing it is refusing to cost."""
        exc = _refused(section="visium_gbm2")
        assert exc.detail.get("platform") == "visium"
        assert exc.detail.get("clamped_dx_um") == 100.0


class TestTheFittedPenetrationLengthIsRefused:
    def test_it_is_refused_by_name(self):
        exc = _refused(readout="penetration_length_um")
        assert exc.status == 422
        assert "16.0405" in exc.message and "17.3158" in exc.message

    def test_the_refusal_compares_it_to_the_noise_this_repo_records(self):
        """16 % is the same order as the 19.4 % a pure placement reseed moves
        it by. A preview whose error is the size of the readout's own noise is
        a different answer, not a preview."""
        assert "19.4" in REFUSED_READOUTS["penetration_length_um"]

    def test_it_is_not_the_previewable_readout(self):
        assert PREVIEWABLE_READOUT == "E_p90_over_p10"
        assert "penetration_length_um" not in (PREVIEWABLE_READOUT,)


class TestThePenetrationIndexIsRefused:
    """A DEVIATION from sec. 8 of the feasibility record, on evidence that
    record did not assemble: its dx agreement is unreconciled across three of
    this repository's own measurements, in two directions."""

    def test_it_is_refused_by_name(self):
        exc = _refused(readout="penetration_index")
        assert exc.status == 422
        assert "0.2247" in exc.message and "0.3035" in exc.message
        assert "0.396" in exc.message and "0.196" in exc.message

    def test_all_three_disagreeing_ladders_are_still_in_the_repo(self):
        """The refusal rests on three records disagreeing. If they stop
        disagreeing, the refusal should be revisited -- and this test is what
        makes that a decision rather than a drift."""
        nr = os.path.join(_ROOT, "results", "numerical_rigor.md")
        ct = os.path.join(_ROOT, "results", "cross_tissue.md")
        with open(nr) as fh:
            nr_text = fh.read()
        with open(ct) as fh:
            ct_text = fh.read()
        assert "0.2247" in nr_text and "0.3035" in nr_text, (
            "numerical_rigor.md sec. 2's rising penetration-index ladder is "
            "gone; re-derive the preview's refusal of this readout")
        assert "0.396" in ct_text and "0.196" in ct_text, (
            "cross_tissue.md's FALLING penetration-index ladder is gone; "
            "re-derive the preview's refusal of this readout")
        assert "Neither ladder should be quoted" in nr_text

    def test_the_disagreement_sits_in_the_rung_a_preview_would_use(self):
        """The third measurement, live from the data file, and the sharpest
        reason to refuse this readout.

        Doxorubicin on ``xenium_319`` at the SAME two grids the published
        ladder used (170x271 and 424x675):

          dx = 12   0.313970  vs the ladder's 0.3035  ->   3.4 % apart
          dx = 30   0.310702  vs the ladder's 0.2247  ->  38.3 % apart

        The two records AGREE at the shipped resolution and disagree by a
        factor of nearly 1.4 at the coarse one. The coarse rung is the rung a
        preview is made of, so the disagreement is not somewhere off to the
        side of this feature -- it is exactly underneath it.
        """
        path = os.path.join(_ROOT, AGREEMENT_SOURCE)
        with open(path) as fh:
            cb = json.load(fh)
        p = next(x for x in cb["pairs"]
                 if x["section"] == "xenium_319" and x["drug"] == "doxorubicin")
        assert p["coarse"]["grid"] == [170, 271]
        assert p["fine"]["grid"] == [424, 675]

        rel = abs(p["coarse"]["penetration_index"]
                  - p["fine"]["penetration_index"]) / p["fine"]["penetration_index"]
        assert rel < 0.02, rel

        fine_gap = abs(p["fine"]["penetration_index"] - 0.3035) / 0.3035
        coarse_gap = abs(p["coarse"]["penetration_index"] - 0.2247) / 0.2247
        assert fine_gap < 0.05, (
            f"the fine rung no longer agrees with numerical_rigor.md sec. 2 "
            f"({fine_gap:.1%} apart); the shape of the disagreement this "
            "refusal rests on has changed")
        assert coarse_gap > 0.25, (
            f"the COARSE rung now agrees with numerical_rigor.md sec. 2 "
            f"({coarse_gap:.1%} apart). If the two records have been "
            "reconciled at dx = 30 um, the penetration index may be "
            "previewable after all -- re-measure and revisit the refusal "
            "deliberately rather than leaving this test red.")
        assert coarse_gap > fine_gap * 5

    def test_dx_is_the_largest_driver_of_this_index_per_the_shipped_caveat(self):
        """Every run already carries that sentence. Previewing a readout across
        a change in its own largest driver needs a settled figure."""
        with open(os.path.join(_ROOT, "app", "server", "jobs.py")) as fh:
            text = fh.read()
        assert "LARGEST single driver of the penetration index" in text


class TestAbsoluteLevelsAreRefused:
    def test_kill_is_refused_whatever_it_is_called(self):
        for name in ("kill_fraction", "KILL_FRACTION", "predicted_kill",
                     "kill"):
            exc = _refused(readout=name)
            assert exc.status == 422
            assert "ABSOLUTE KILL FRACTION IS NEVER A DEFENSIBLE OUTPUT" in exc.message

    def test_the_best_agreeing_readout_is_still_refused(self):
        """``median_auc_over_peak`` agrees to 0.3215 % -- better than the
        readout that IS previewed -- and is refused anyway, because agreement
        is not defensibility."""
        exc = _refused(readout="median_auc_over_peak")
        assert "0.3215" in exc.message
        assert "agreement is not defensibility" in exc.message

    def test_the_reach_ladder_is_refused(self):
        exc = _refused(readout="reach_headline_F")
        assert "27.5" in exc.message

    def test_an_unknown_readout_is_refused_rather_than_attempted(self):
        exc = _refused(readout="whatever_comes_next")
        assert "NO MEASURED" in exc.message
        assert PREVIEWABLE_READOUT in exc.message


class TestOnlyTheMeasuredEnvelopeIsServed:
    """The measured agreement travels with the number, so a pair nobody
    measured has no number to show."""

    def test_an_unmeasured_section_is_refused(self):
        exc = _refused(section="xenium_io")
        assert exc.status == 422
        assert "has not been measured" in exc.message
        assert exc.detail["measured_sections"] == list(MEASURED_SECTIONS)

    def test_an_unmeasured_compound_is_refused(self):
        exc = _refused(drug="not_a_compound")
        assert "has not been measured" in exc.message
        assert exc.detail["measured_compounds"] == list(MEASURED_COMPOUNDS)

    def test_a_target_dx_nobody_paired_is_refused(self):
        """``api.py::_default_dx_um`` returns 21 um on a memory-constrained
        host. There is no 30/21 pair, so there is no preview there -- and that
        is a measurement to take, not a tolerance to widen."""
        exc = _refused(dx_um=21.0)
        assert "NO MEASURED PAIR" in exc.message
        assert exc.detail["measured_target_dx_um"] == TARGET_DX_UM

    def test_a_run_already_at_the_preview_resolution_is_refused(self):
        exc = _refused(dx_um=PREVIEW_DX_UM)
        assert exc.status == 409
        assert "ALREADY AT dx = 30 um" in exc.message

    def test_every_refusal_says_what_to_do_next(self):
        """WORKSPACE_CONTRACT sec. 6: a refusal says what to do next. Never a
        bare dash."""
        cases = [
            {"section": "visium_gbm2"},
            {"section": "xenium_io"},
            {"drug": "not_a_compound"},
            {"dx_um": 21.0},
            {"dx_um": PREVIEW_DX_UM},
            {"readout": "penetration_length_um"},
            {"readout": "penetration_index"},
            {"readout": "median_auc_over_peak"},
            {"readout": "kill_fraction"},
        ]
        for over in cases:
            msg = _refused(**over).message.lower()
            assert any(w in msg for w in (
                "run the real solve", "wait for the run", "submit the run",
                "add the pair", "add the pairs", "settle the ladder",
                "settle the coarse rung", "measure the pairs",
                "read the curve", "post /api/runs", "solve it at both",
                "what to do")), f"{over}: no next action in {msg[:200]!r}"


# --------------------------------------------------------------------------- #
# 3. The number never travels alone
# --------------------------------------------------------------------------- #
class TestThePayloadCannotShedItsLabelOrItsErrorBar:
    def test_a_real_preview_carries_both(self):
        pytest.importorskip("numpy")
        row = _section_row("xenium_319")
        if row is None or not row.get("available"):
            pytest.skip("xenium_319 is not reconstructable on this machine")
        body = preview_payload(_Jobs(), dict(BASE_SPEC))
        assert body["is_preview"] is True
        assert body["quotable"] is False
        assert body["preview_label"] == PREVIEW_LABEL
        assert body["agreement"]["n_pairs"] == 12
        assert body["readout"] == PREVIEWABLE_READOUT
        assert body["preview_dx_um"] == PREVIEW_DX_UM
        assert body["target_dx_um"] == TARGET_DX_UM
        assert body["value"] is not None
        # The engine's own limits, carried rather than summarised.
        assert body["caveats"], "the engine's caveats were dropped"
        assert any("dx = 30 um" in c for c in body["caveats"]), (
            "the RESOLUTION caveat does not state the PREVIEW's own grid")

    def test_the_label_says_do_not_quote_and_names_both_grids(self):
        assert "DO NOT QUOTE THIS NUMBER" in PREVIEW_LABEL
        assert "30 um" in PREVIEW_LABEL and "12 um" in PREVIEW_LABEL
        assert "no run_id" in PREVIEW_LABEL

    def test_a_body_without_the_label_is_refused_by_the_audit(self):
        good = {"is_preview": True, "quotable": False,
                "preview_label": PREVIEW_LABEL, "value": 2.0,
                "agreement": measured_agreement()}
        _audit_payload(dict(good))
        for drop in ("preview_label", "agreement"):
            bad = dict(good)
            bad.pop(drop)
            with pytest.raises(PreviewError) as ei:
                _audit_payload(bad)
            assert ei.value.status == 500, drop

    def test_a_body_claiming_to_be_quotable_is_refused(self):
        bad = {"is_preview": True, "quotable": True,
               "preview_label": PREVIEW_LABEL, "value": 2.0,
               "agreement": measured_agreement()}
        with pytest.raises(PreviewError):
            _audit_payload(bad)

    def test_a_refused_readout_cannot_ride_along_in_the_payload(self):
        bad = {"is_preview": True, "quotable": False,
               "preview_label": PREVIEW_LABEL, "value": 2.0,
               "agreement": measured_agreement(),
               "extras": {"penetration_length_um": 277.7}}
        with pytest.raises(PreviewError) as ei:
            _audit_payload(bad)
        assert "penetration_length" in ei.value.message

    def test_an_absolute_level_cannot_ride_along_in_the_payload(self):
        for key in ("kill_fraction", "median_tumor_auc", "E_p90_uM",
                    "therapeutic_coverage_1x"):
            bad = {"is_preview": True, "quotable": False,
                   "preview_label": PREVIEW_LABEL, "value": 2.0,
                   "agreement": measured_agreement(), key: 1.0}
            with pytest.raises(PreviewError):
                _audit_payload(bad)

    def test_no_preview_is_served_without_the_measurement(self, monkeypatch):
        """If the data file goes, THE FEATURE GOES -- not the error bar.

        A server that kept serving the number after losing the statement of how
        far it may be from the real solve would be serving a wrong number
        quickly, which is the one outcome this whole module is written against.
        """
        import app.server.preview as mod

        monkeypatch.setattr(mod, "_agreement_path",
                            lambda: "/nonexistent/coarse_baseline.json")
        assert mod.measured_agreement() is None
        with pytest.raises(PreviewError) as ei:
            mod.preview_payload(_Jobs(), dict(BASE_SPEC))
        assert ei.value.status == 503
        assert "THE PREVIEW IS UNAVAILABLE" in ei.value.message


def _section_row(section):
    from app.server.cache import list_sections

    return next((r for r in list_sections() if r.get("id") == section), None)


# --------------------------------------------------------------------------- #
# 4. A preview is not a run
# --------------------------------------------------------------------------- #
class TestAPreviewIsNotARun:
    def test_nothing_is_persisted_and_no_run_id_is_minted(self):
        row = _section_row("xenium_319")
        if row is None or not row.get("available"):
            pytest.skip("xenium_319 is not reconstructable on this machine")
        body = preview_payload(_Jobs(), dict(BASE_SPEC))
        assert "run_id" not in body, (
            "a preview with a run_id is a run: it can be opened, attached to a "
            "specimen, exported and reported, which are the four routes by "
            "which a fast number becomes a quoted one")
        assert "_payloads" not in body

    def test_the_module_never_enqueues(self):
        """Read from the source, because the behaviour is an absence."""
        with open(os.path.join(_ROOT, "app", "server", "preview.py")) as fh:
            src = fh.read()
        for call in ("jobs.submit", ".submit(", "JobManager("):
            assert call not in src, (
                f"{call!r} in preview.py: a preview must not enter the run "
                "store")
        assert 'result.pop("_payloads", None)' in src, (
            "the run payloads are no longer dropped; a preview that writes "
            "blobs is a run")


class TestTheRouterMountsAndAnswers:
    """The one line that mounts this feature, exercised.

    ``app/server/api.py`` is owned by other work this round, so the mount is
    reported rather than made; this builds the same router onto a bare app so
    that "it mounts in one line" is a tested claim and not an assurance.
    """

    def _client(self):
        fastapi = pytest.importorskip("fastapi")
        pytest.importorskip("fastapi.testclient")
        from fastapi.testclient import TestClient

        from app.server.preview import build_preview_router

        app = fastapi.FastAPI()
        # THE MOUNT LINE, verbatim apart from the local `app`/`jobs` names.
        app.include_router(build_preview_router(jobs=_Jobs()))
        return TestClient(app)

    def test_the_agreement_route_answers_without_solving_anything(self):
        r = self._client().get("/api/preview/agreement")
        assert r.status_code == 200, r.text
        body = r.json()
        assert body["available"] is True
        assert body["agreement"]["n_pairs"] == 12
        assert body["previewable"] == [PREVIEWABLE_READOUT]
        assert body["preview_label"] == PREVIEW_LABEL
        # The refusals are readable BEFORE spending the 8.7 s, so the absence
        # of a measurement is discoverable rather than only encountered.
        names = {r_["readout"] for r_ in body["refused"]}
        assert {"penetration_length_um", "penetration_index",
                "median_auc_over_peak"} <= names
        assert "_reconstruct_visium" in body["spot_refusal"]

    def test_a_refusal_comes_back_as_a_structured_detail(self):
        r = self._client().post("/api/preview", json=dict(
            BASE_SPEC, section="visium_gbm2"))
        assert r.status_code == 422, r.text
        detail = r.json()["detail"]
        assert detail["error"] == SPOT_PREVIEW_REFUSAL
        assert detail["platform"] == "visium"

    def test_a_non_object_body_is_a_400_not_a_stack(self):
        r = self._client().post("/api/preview", json=[1, 2, 3])
        assert r.status_code == 400, r.text

    def test_the_route_is_under_api_so_the_gate_covers_it(self):
        """``create_app`` authenticates by path prefix: it tests
        ``get_route_path(request.scope).startswith("/api/")``. A router mounted
        there is gated by construction and nothing about it opts in."""
        from app.server.preview import build_preview_router

        router = build_preview_router(jobs=_Jobs())
        paths = {r_.path for r_ in router.routes}
        assert paths == {"/api/preview", "/api/preview/agreement"}, paths
        assert all(p.startswith("/api/") for p in paths)


# --------------------------------------------------------------------------- #
# 5. The page, in a browser
# --------------------------------------------------------------------------- #
pytest.importorskip("playwright.sync_api")
from playwright.sync_api import sync_playwright  # noqa: E402

_SECTIONS = [{
    "id": "xenium_319", "label": "Xenium GBM - 319-panel", "n_cells": 40887,
    "extent_um": [8087.8, 5065.4], "panel": "hBrain_254g", "available": True,
    "source": "xenium", "sample": None, "cache_dir": "data_cache/xenium_gbm",
    "assay": "Xenium In Situ", "accession": "10x", "has_outcome_label": False,
    "reason": None, "remedy": [], "spatial_spread_valid": True, "caveats": [],
}]
_DEFAULTS = dict(BASE_SPEC)

#: What `POST /api/preview` answers in the happy path. Shaped exactly like the
#: real payload, with the real label and a real agreement block.
_PREVIEW_OK = {
    "is_preview": True,
    "quotable": False,
    "preview_label": PREVIEW_LABEL,
    "readout": PREVIEWABLE_READOUT,
    "readout_name": "delivery spread  P90 / P10",
    "value": 1.9935720609199155,
    "unit": "x",
    "null_reason": None,
    "n_tumor": 23812,
    "derivation": "extra_metrics.E_p90_uM / extra_metrics.E_p10_uM from a "
                  "dx = 30 um solve, over the tumour cell set.",
    "agreement": {
        "readout": PREVIEWABLE_READOUT,
        "preview_dx_um": 30.0, "target_dx_um": 12.0,
        "n_pairs": 12, "n_sections": 2, "n_compounds": 6,
        "sections": ["xenium_319", "xenium_breast"],
        "compounds": ["cisplatin", "doxorubicin", "lomustine",
                      "methotrexate", "temozolomide", "vincristine"],
        "median_rel_err_pct": 0.7132559626060712,
        "p90_rel_err_pct": 1.2001473762142536,
        "max_rel_err_pct": 1.7619996350758897,
        "median_abs_err": 0.024797022150895964,
        "max_abs_err": 0.09089332263845673,
        "across_donor_range": 0.5357238964376265,
        "across_donor_range_source": "results/surrogate_feasibility_2026-10-09.md sec. 5.2",
        "nuisance_band_width": 0.202,
        "nuisance_band_source": "results/delivery_spread_noise.md sec. 0",
        "coarse_median_s": 8.667764902114868,
        "target_median_s": 58.77110505104065,
        "cost_ratio": 6.780422140510635,
        "loadavg_1m_min": 2.37451171875, "loadavg_1m_max": 9.7802734375,
        "n_cpu": [10], "pairs_discarded_same_grid": 18,
        "source": AGREEMENT_SOURCE,
        "record": "results/surrogate_feasibility_2026-10-09.md sec. 6",
    },
    "preview_dx_um": 30.0, "target_dx_um": 12.0,
    "preview_grid": [170, 271], "preview_wall_clock_s": 8.4,
    "section": "xenium_319", "drug": "temozolomide",
    "refused": [{"readout": k, "why": v} for k, v in sorted(REFUSED_READOUTS.items())],
    "caveats": ["RESOLUTION: this run used dx = 30 um. THERE IS NO dx AT WHICH "
                "THIS PIPELINE HAS BEEN SHOWN TO BE GRID-CONVERGED."],
    "claim_set": {"version": 1, "adopted": "2026-08-21"},
    "superseded_by": "the real run.",
}


@pytest.fixture(scope="module")
def web_server():
    handler = functools.partial(http.server.SimpleHTTPRequestHandler,
                                directory=WEB)

    class Quiet(socketserver.TCPServer):
        allow_reuse_address = True

        def handle_error(self, request, client_address):  # noqa: D102
            pass

    with Quiet(("127.0.0.1", 0), handler) as httpd:
        port = httpd.server_address[1]
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        try:
            yield f"http://127.0.0.1:{port}"
        finally:
            httpd.shutdown()


def _open_studio(server, *, preview):
    """Open the REAL Studio page with /api/preview answered by ``preview``.

    ``preview`` is either a dict (200) or an ``(status, detail)`` tuple, which
    is how the server answers every refusal.
    """
    p = sync_playwright().start()
    browser = p.chromium.launch()
    pg = browser.new_page(viewport={"width": 1440, "height": 900})
    errors = []
    pg.on("pageerror", lambda e: errors.append(str(e)))
    console = []
    pg.on("console", lambda m: console.append(f"{m.type}: {m.text}"))

    def route_api(route):
        url = route.request.url
        if "/api/preview" in url:
            if isinstance(preview, tuple):
                status, detail = preview
                return route.fulfill(
                    status=status, content_type="application/json",
                    body=json.dumps({"detail": detail}))
            return route.fulfill(status=200, content_type="application/json",
                                 body=json.dumps(preview))
        if "/api/sections" in url:
            return route.fulfill(status=200, content_type="application/json",
                                 body=json.dumps(_SECTIONS))
        if "/api/defaults" in url:
            return route.fulfill(status=200, content_type="application/json",
                                 body=json.dumps(_DEFAULTS))
        if "/api/validate" in url:
            return route.fulfill(status=200, content_type="application/json",
                                 body=json.dumps({"ok": True, "errors": [],
                                                  "warnings": []}))
        if "/api/drugs" in url:
            return route.fulfill(status=200, content_type="application/json",
                                 body=json.dumps([]))
        if "/api/runs" in url and route.request.method == "POST":
            return route.fulfill(status=202, content_type="application/json",
                                 body=json.dumps({"run_id": "r_00042"}))
        return route.fulfill(status=200, content_type="application/json",
                             body="[]")

    pg.route("**/api/**", route_api)
    pg.goto(f"{server}/index.html", wait_until="domcontentloaded")
    pg.wait_for_timeout(2200)
    return p, browser, pg, errors, console


def _click_preview(pg):
    pg.wait_for_selector("#sect-preview", timeout=8000)
    btn = pg.query_selector("#sect-preview button.btn")
    assert btn is not None, "the Preview control did not render in the rail"
    btn.click()
    pg.wait_for_timeout(900)


class TestThePreviewPanelPaintsInARealBrowser:
    """A panel that parses is not a panel that painted. `h()` throws on a
    reversed selector and the throw takes the whole subtree with it, so a
    missing preview looks exactly like a rail still loading."""

    def test_the_panel_is_mounted_in_the_commit_group(self, web_server):
        p, browser, pg, errors, _ = _open_studio(web_server, preview=_PREVIEW_OK)
        try:
            assert not errors, errors[:3]
            el = pg.query_selector("#sect-preview")
            assert el is not None, "no #sect-preview section in the Setup rail"
            group = pg.eval_on_selector(
                "#sect-preview", "el => el.closest('.rail-group').id")
            assert group == "rg-job", (
                f"the preview is mounted in {group!r}; it belongs in the "
                "Commit group, beside the button it is advising on")
            # And NOT in the Read rail, which is where results live.
            rail = pg.eval_on_selector(
                "#sect-preview", "el => el.closest('.rail').className")
            assert "rail-left" in rail, rail
        finally:
            browser.close()
            p.stop()

    def test_the_number_renders_with_the_word_preview_beside_it(self, web_server):
        p, browser, pg, errors, _ = _open_studio(web_server, preview=_PREVIEW_OK)
        try:
            _click_preview(pg)
            val = pg.query_selector("#preview-value")
            assert val is not None, (
                "no preview value rendered.\n"
                + pg.inner_text("#sect-preview")[:500])
            assert "1.99" in val.inner_text()
            body = pg.inner_text("#sect-preview")
            assert "PREVIEW" in body
            assert "NOT A RESULT" in body
            assert "DO NOT QUOTE THIS NUMBER" in body
            assert not errors, errors[:3]
        finally:
            browser.close()
            p.stop()

    def test_the_measured_agreement_is_on_screen_with_the_number(self, web_server):
        """Constraint: the measured agreement travels WITH the number. A
        preview that does not state how far it may be from the real solve is
        just a wrong number shown quickly."""
        p, browser, pg, errors, _ = _open_studio(web_server, preview=_PREVIEW_OK)
        try:
            _click_preview(pg)
            body = pg.inner_text("#sect-preview")
            for figure in ("0.7133", "1.2001", "1.7620", "0.0909", "12"):
                assert figure in body, (
                    f"the measured agreement is missing {figure!r} from the "
                    f"panel.\n{body[:600]}")
            assert "6.78" in body, "the cost ratio is not shown"
            assert AGREEMENT_SOURCE in body, "the error bar names no source"
        finally:
            browser.close()
            p.stop()

    def test_the_caveat_is_never_smaller_than_the_number(self, web_server):
        """docs/WORKSPACE_CONTRACT.md sec. 6, measured in the browser rather
        than read off a stylesheet."""
        p, browser, pg, errors, console = _open_studio(
            web_server, preview=_PREVIEW_OK)
        try:
            _click_preview(pg)
            sizes = pg.evaluate(
                """() => {
                  const v = document.querySelector('#preview-value');
                  const keys = Array.from(document.querySelectorAll(
                      '#sect-preview .caveat.sev-2 .cv-key'));
                  const flagged = document.querySelector(
                      '#sect-preview [data-caveat-floor]');
                  return {
                    value: v ? parseFloat(getComputedStyle(v).fontSize) : null,
                    keys: keys.map(k => parseFloat(getComputedStyle(k).fontSize)),
                    flag: flagged ? flagged.dataset.caveatFloor : null,
                  };
                }""")
            assert sizes["value"], sizes
            assert sizes["keys"], "no hard-stop subject line beside the number"
            smallest = min(sizes["keys"])
            assert smallest >= sizes["value"], (
                f"the caveat subject is {smallest}px against a "
                f"{sizes['value']}px number. docs/WORKSPACE_CONTRACT.md sec. 6: "
                "a caveat is never set smaller than the number it qualifies.")
            assert sizes["flag"] == "ok", sizes
            assert not [c for c in console if "WORKSPACE_CONTRACT" in c], console

    # The panel's own in-browser assertion must agree with this test's
    # measurement; if it does not, one of the two is reading the wrong node.
        finally:
            browser.close()
            p.stop()

    def test_the_error_bar_itself_is_as_large_as_the_number(self, web_server):
        """The disagreement figure is the other half of the number, so it is a
        hard stop too -- not a footnote in small print."""
        p, browser, pg, errors, _ = _open_studio(web_server, preview=_PREVIEW_OK)
        try:
            _click_preview(pg)
            txt = pg.evaluate(
                """() => Array.from(document.querySelectorAll(
                     '#sect-preview .caveat.sev-2 .cv-subj'))
                     .map(e => e.innerText)""")
            assert any("MAY DIFFER FROM THE REAL SOLVE" in t for t in txt), txt
        finally:
            browser.close()
            p.stop()


class TestARefusalRendersAsARefusal:
    def test_the_visium_refusal_renders_in_full_with_no_number(self, web_server):
        p, browser, pg, errors, _ = _open_studio(
            web_server,
            preview=(422, {"error": SPOT_PREVIEW_REFUSAL,
                           "platform": "visium",
                           "remedy": "run the real solve"}))
        try:
            _click_preview(pg)
            body = pg.inner_text("#sect-preview")
            assert pg.query_selector("#preview-value") is None, (
                "a refusal left a number on screen")
            assert "NO PREVIEW ON A SPOT SECTION" in body
            assert "_reconstruct_visium" in body
            assert "run the real solve" in body.lower()
            assert not errors, errors[:3]
        finally:
            browser.close()
            p.stop()

    def test_the_penetration_length_refusal_renders(self, web_server):
        p, browser, pg, errors, _ = _open_studio(
            web_server,
            preview=(422, {"error": REFUSED_READOUTS["penetration_length_um"]}))
        try:
            _click_preview(pg)
            body = pg.inner_text("#sect-preview")
            assert "16.0405" in body and "19.4" in body
            assert pg.query_selector("#preview-value") is None
        finally:
            browser.close()
            p.stop()

    def test_the_unavailable_measurement_renders_as_no_feature(self, web_server):
        p, browser, pg, errors, _ = _open_studio(
            web_server,
            preview=(503, {"error": "THE PREVIEW IS UNAVAILABLE ON THIS SERVER "
                                    "because the measured agreement between "
                                    "the coarse and the shipped grid is not "
                                    "present."}))
        try:
            _click_preview(pg)
            body = pg.inner_text("#sect-preview")
            assert "THE PREVIEW IS UNAVAILABLE" in body
            assert pg.query_selector("#preview-value") is None, (
                "a number was shown with no error bar beside it")
        finally:
            browser.close()
            p.stop()


class TestThePreviewIsReplacedByTheRealSolve:
    def test_submitting_a_run_drops_the_preview_number(self, web_server):
        """Not kept beside, not annotated. One readout with two live values is
        an invitation to quote the convenient one."""
        p, browser, pg, errors, _ = _open_studio(web_server, preview=_PREVIEW_OK)
        try:
            _click_preview(pg)
            assert pg.query_selector("#preview-value") is not None
            pg.wait_for_timeout(600)
            btn = pg.query_selector("#btn-submit-run")
            assert btn is not None
            if btn.is_disabled():
                pytest.skip("Submit is gated in this fixture; the supersede "
                            "path is exercised by the unit below")
            btn.click()
            pg.wait_for_timeout(500)
            assert pg.query_selector("#preview-value") is None, (
                "the preview number survived the real run's submission")
            body = pg.inner_text("#sect-preview")
            assert "Superseded" in body
        finally:
            browser.close()
            p.stop()

    def test_the_submit_handler_supersedes_before_posting(self):
        """Read from the source: the drop happens BEFORE the POST, so there is
        no window in which a real run is in flight and a preview number is
        still on screen."""
        with open(os.path.join(_ROOT, "app", "web", "js", "ui", "spec.js")) as fh:
            src = fh.read()
        m = re.search(r"previewPanel\.supersede\(\);\s*onSubmit\(", src)
        assert m, (
            "the Submit handler no longer drops the preview before posting; a "
            "preview number beside a running real solve is two live values for "
            "one readout")


class TestThePreviewSurfaceMakesNoProhibitedClaim:
    def test_no_ai_ml_or_clinical_framing(self):
        for path in (os.path.join(_ROOT, "app", "server", "preview.py"),
                     os.path.join(_ROOT, "app", "web", "js", "ui", "preview.js")):
            with open(path) as fh:
                text = fh.read().lower()
            for phrase in ("digital twin", "personalised medicine",
                           "personalized medicine", "machine learning",
                           "predicts response", "treatment plan",
                           "diagnos"):
                assert phrase not in text, f"{phrase!r} in {path}"
            for m in re.finditer(r"upload", text):
                window = text[max(0, m.start() - 120):m.end() + 120]
                assert not any(w in window for w in
                               ("tissue", "section", "specimen", "xenium",
                                "visium")), window[:200]

    def test_the_js_module_parses_as_a_module(self):
        """Every JS module in this app must pass
        `node --input-type=module --check`."""
        import shutil
        import subprocess

        node = shutil.which("node")
        if node is None:
            pytest.skip("node is not installed")
        path = os.path.join(_ROOT, "app", "web", "js", "ui", "preview.js")
        with open(path, "rb") as fh:
            src = fh.read()
        r = subprocess.run([node, "--input-type=module", "--check"],
                           input=src, capture_output=True)
        assert r.returncode == 0, r.stderr.decode()[:600]

    def test_no_h_selector_puts_classes_before_the_id(self):
        """`h('div.a#id')` THROWS and takes the whole subtree with it."""
        path = os.path.join(_ROOT, "app", "web", "js", "ui", "preview.js")
        with open(path) as fh:
            src = fh.read()
        bad = re.findall(r"""h\(\s*['"]([a-zA-Z0-9-]*(?:\.[^.#'"\s]+)+#[^'"\s]+)['"]""",
                         src)
        assert not bad, bad
