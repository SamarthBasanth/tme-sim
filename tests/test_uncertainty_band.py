"""The interval as a first-class output — and the guards that keep it honest.

What is under test here is not arithmetic. The arithmetic is four lines. What is
under test is that the feature cannot quietly become the thing it exists to
prevent:

1. **The two kinds must stay apart.** The potency bar is BIOLOGICAL (patients
   differ; no better model narrows it) and the cell-state bar is EPISTEMIC (our
   ignorance; one Ki-67 stain deletes it). A payload that merged them would tell
   a reader the opposite of the truth about both. ``test_plant_*`` merges them
   and watches the guard go red.
2. **A band must never be fabricated.** A compound with no recorded span, a
   quoted round-number literature range, a section the span was not measured on
   — each must REFUSE, not interpolate. Each refusal is planted and checked.
3. **A band must be a band around THIS run.** The re-evaluation is gated on
   reproducing the run's own ``kill_fraction`` to within the engine's Bernoulli
   sampling scale. The gate is planted by corrupting the exposure.
4. **Absolute kill fraction is still not a defensible output.** Every interval
   carries ``quotable = False``. The interval is the size of the refusal.

Every guard here is first exercised on the HEALTHY object and asserted green —
so that when the planted version goes red we know the check reached its subject
and did not pass by accident. That is the failure mode this repository has had
seven times in one week (``results/`` 2026-08 vacuous-green audit) and it is
cheaper to prevent than to find.

The run is REAL: the 40,887-cell open 10x Xenium GBM section, reconstructed and
simulated through the production entry points at a coarse grid and a short
window so the file runs in seconds.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import time

import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from tmesim.core.config import PDEConfig  # noqa: E402
from tmesim.sim.engine import SimConfig, simulate_response  # noqa: E402
from tmesim.validation.calibration import calibrate_alpha  # noqa: E402
from tmesim.validation.uncertainty import (  # noqa: E402
    BIOLOGICAL,
    EPISTEMIC,
    F_STATE_BY_ROUTE,
    POTENCY_SPAN,
    REPRODUCTION_SIGMA_GATE,
    STATE_SPAN_FACTOR,
    UNCHARACTERISED_METRICS,
    FinishedRun,
    assert_kinds_separated,
    finished_run_from_result,
    kill_ceiling_interval,
    kill_fraction_interval,
    run_interval_report,
)

#: The frozen single-anchor calibration constant. Asserted at the TOP and the
#: BOTTOM of this file: a band is a re-evaluation of the run's own Hill, so an
#: alpha that moved under the feature would move every endpoint silently.
ALPHA_FROZEN = 22.223477711172563

#: Same short, real configuration tests/test_mechanism.py uses.
PDE = PDEConfig(dt=10.0, total_time=1800.0, convection=True, save_every=100)

FLAGSHIP_SECTION = "xenium_319"


def test_calibrate_alpha_frozen_at_start():
    assert calibrate_alpha() == ALPHA_FROZEN, calibrate_alpha()


@pytest.fixture(scope="module")
def real_run(real_micro, real_drug) -> FinishedRun:
    """A REAL finished run, reduced to the per-cell arrays the band needs.

    ``finished_run_from_result`` is the default-on in-process path: it reads the
    exposure and oxygen fields the engine just returned and re-derives nothing
    the engine did not already decide.
    """
    res = simulate_response(real_micro, SimConfig(drug=real_drug, pde=PDE, seed=0))
    return finished_run_from_result(
        real_micro, res, real_drug,
        compound="temozolomide", section_id=FLAGSHIP_SECTION, run_id="test")


# --------------------------------------------------------------------------- #
# 1. The transcription is checked against the RECORD, not against itself
# --------------------------------------------------------------------------- #
def _read_result(name: str) -> str:
    with open(os.path.join(_REPO, "results", name)) as fh:
        return fh.read()


def test_potency_span_matches_the_dated_record():
    """The stored IC50 endpoints must be the ones results/ actually records.

    A span transcribed once and then drifting from its source is this project's
    recurring defect (results/codon_claim-surface-drift). So the endpoints are
    grepped out of the record rather than trusted.
    """
    doc = _read_result("kill_calibration.md")
    tmz = POTENCY_SPAN["temozolomide"]
    assert "151.1-1882" in doc, "the temozolomide IC50 band moved in the record"
    assert (tmz.lo, tmz.hi, tmz.point) == (151.1, 1882.0, 362.1)
    cis = POTENCY_SPAN["cisplatin"]
    assert "8.929-84.7" in doc
    assert (cis.lo, cis.hi, cis.point) == (8.929, 84.7, 24.55)
    # The combined temozolomide band the whole feature is named after.
    assert "0.0048-0.1442" in doc


def test_cell_state_span_matches_the_dated_record():
    doc = _read_result("cell_state_impact.md")
    for v in ("0.5502", "0.6490", "0.9550"):
        assert v in doc, f"F_state {v} is no longer in the record"
    # The recorded F_state values must equal the recorded band occupancies run
    # through the recorded multipliers. If they ever stop agreeing, one of the
    # two tables in the record is stale and the span is meaningless.
    #
    # TOLERANCE, from the record's OWN precision and nothing else: the
    # occupancies are quoted to two decimals of a percent (+-5e-5 each) and
    # F_state to four decimals (+-5e-5), so the arithmetic can agree to at best
    # 1*5e-5 + 0.5*5e-5 + 0.35*5e-5 + 5e-5 ~ 1.4e-4. The vessel-distance row
    # lands 6e-5 out for exactly that reason (0.95506 recomputed vs 0.9550
    # recorded) and a tighter tolerance would be testing the rounding.
    _ROUNDING = 2e-4
    assert F_STATE_BY_ROUTE["oxygen"] == pytest.approx(
        0.1004 + 0.5 * 0.8996, abs=_ROUNDING)
    assert F_STATE_BY_ROUTE["markers"] == pytest.approx(
        0.2980 + 0.5 * 0.7020, abs=_ROUNDING)
    assert F_STATE_BY_ROUTE["vessel_distance"] == pytest.approx(
        0.9108 + 0.5 * 0.0867 + 0.35 * 0.0026, abs=_ROUNDING)
    # results/kill_calibration.md sec 5.1's own arithmetic: 0.0326 x factor
    # is its recorded state-route endpoint, 0.0566.
    assert round(0.0326 * STATE_SPAN_FACTOR[1], 4) == 0.0566


# --------------------------------------------------------------------------- #
# 2. The reproduction gate — and the plant that fires it
# --------------------------------------------------------------------------- #
def test_reproduction_gate_is_reached_and_green(real_run):
    """GUARD REACHED, healthy case: the arrays reproduce the run's own kill."""
    iv = kill_fraction_interval(real_run)
    rep = iv["reproduction"]
    assert rep["n_tumor"] > 0
    assert rep["passed"] is True, rep
    assert abs(rep["z"]) <= REPRODUCTION_SIGMA_GATE, rep
    # The residual must be sampling noise, not agreement by construction: a
    # gate that can only pass is not a gate.
    assert rep["bernoulli_sigma"] > 0.0
    assert iv["status"] == "BANDED"


def test_plant_corrupted_exposure_refuses_the_band(real_run):
    """PLANT: double every cell's exposure. The run's kill no longer describes
    these arrays, and the band must be REFUSED rather than drawn around the
    wrong number."""
    planted = FinishedRun(
        compound=real_run.compound, section_id=real_run.section_id,
        kill_fraction=real_run.kill_fraction,
        E_uM=real_run.E_uM * 2.0,            # <-- the defect
        modulation=real_run.modulation,
        ic50_uM=real_run.ic50_uM, hill=real_run.hill,
        max_kill_rate=real_run.max_kill_rate,
        state_mode=real_run.state_mode, dx_um=real_run.dx_um)
    iv = kill_fraction_interval(planted)
    assert iv["reproduction"]["passed"] is False, iv["reproduction"]
    assert abs(iv["reproduction"]["z"]) > REPRODUCTION_SIGMA_GATE
    assert iv["status"] == "REFUSED"
    assert iv["band"] is None
    assert "do not reproduce" in iv["reason"]


# --------------------------------------------------------------------------- #
# 3. The two kinds — and the plant that merges them
# --------------------------------------------------------------------------- #
def test_both_kinds_are_present_and_mean_opposite_things(real_run):
    iv = kill_fraction_interval(real_run)
    comps = iv["components"]
    assert set(comps) == {"potency", "cell_state_route"}
    assert comps["potency"]["kind"] == BIOLOGICAL
    assert comps["cell_state_route"]["kind"] == EPISTEMIC
    assert sorted(iv["covers_kinds"]) == [BIOLOGICAL, EPISTEMIC]

    # THE distinction, asserted rather than merely written down: a better model
    # removes one bar and cannot touch the other.
    assert comps["potency"]["does_a_better_model_shrink_it"] == "no"
    assert comps["cell_state_route"]["does_a_better_model_shrink_it"] == "yes"
    assert comps["potency"]["what_shrinks_it"] != \
        comps["cell_state_route"]["what_shrinks_it"]
    assert "this patient" in comps["potency"]["what_shrinks_it"].lower()
    assert "ki-67" in comps["cell_state_route"]["what_shrinks_it"].lower()

    # Each component names the dated record it came out of.
    for c in comps.values():
        assert c["source_document"].startswith("results/")
        assert c["source_section"]


def test_guard_is_reached_and_green_on_a_healthy_interval(real_run):
    """GUARD REACHED: assert_kinds_separated passes on the real payload, so the
    red result in the next test is a reaction to the plant and not to the guard
    never having been exercised."""
    assert_kinds_separated(kill_fraction_interval(real_run))
    assert_kinds_separated(kill_ceiling_interval(real_run))


def test_plant_merging_the_two_kinds_is_rejected(real_run):
    """PLANT: collapse the two components into one undifferentiated bar — the
    single most damaging thing this feature could do."""
    iv = kill_fraction_interval(real_run)
    lo = min(c["lo"] for c in iv["components"].values())
    hi = max(c["hi"] for c in iv["components"].values())
    iv["components"] = {"uncertainty": {"lo": lo, "hi": hi}}   # <-- the defect
    with pytest.raises(AssertionError, match="unlabelled or unknown kinds"):
        assert_kinds_separated(iv)


def test_plant_dropping_the_kind_label_is_rejected(real_run):
    """PLANT: keep both bars but strip what makes them different."""
    iv = kill_fraction_interval(real_run)
    iv["components"]["potency"].pop("what_shrinks_it")          # <-- the defect
    with pytest.raises(AssertionError, match="what_shrinks_it"):
        assert_kinds_separated(iv)


def test_plant_marking_an_absolute_kill_quotable_is_rejected(real_run):
    """PLANT: flip the refusal. An absolute-kill interval that presents itself
    as quotable is the prohibition this whole feature is built around."""
    iv = kill_fraction_interval(real_run)
    iv["quotable"] = True                                       # <-- the defect
    with pytest.raises(AssertionError, match="quotable=False"):
        assert_kinds_separated(iv)


# --------------------------------------------------------------------------- #
# 4. Refusals — never a fabricated band
# --------------------------------------------------------------------------- #
def _with_compound(run: FinishedRun, compound: str,
                   section_id: str = None) -> FinishedRun:
    return FinishedRun(
        compound=compound,
        section_id=section_id if section_id is not None else run.section_id,
        kill_fraction=run.kill_fraction, E_uM=run.E_uM,
        modulation=run.modulation, ic50_uM=run.ic50_uM, hill=run.hill,
        max_kill_rate=run.max_kill_rate, state_mode=run.state_mode,
        dx_um=run.dx_um)


def test_compound_with_no_recorded_span_refuses_the_biological_bar(real_run):
    """PLANT: ask for a band on a compound nothing has ever measured a
    between-patient spread for. It must refuse the biological bar and SAY the
    remaining band is narrower than the real uncertainty."""
    iv = kill_fraction_interval(_with_compound(real_run, "lomustine"))
    assert iv["covers_kinds"] == [EPISTEMIC]
    assert iv["status"] == "BANDED_PARTIAL"
    assert "NARROWER" in iv["combination"]
    refusal = next(r for r in iv["refusals"]
                   if r["metric"] == "kill_fraction.biological")
    assert refusal["band"] is None
    assert "no recorded potency span" in refusal["reason"]
    assert "lomustine" in refusal["reason"]


def test_doxorubicin_quoted_range_is_refused_not_banded(real_run):
    """A quoted 0.1-1 uM literature range has no population, no n and no
    quantile. results/kill_calibration.md sec 4 says it is NOT a measured anchor
    in this disease, so it must not become an error bar."""
    assert POTENCY_SPAN["doxorubicin"].bandable is False
    iv = kill_fraction_interval(_with_compound(real_run, "doxorubicin"))
    assert "potency" not in iv["components"]
    refusal = next(r for r in iv["refusals"]
                   if r["metric"] == "kill_fraction.biological")
    assert "quoted_literature_range" in refusal["reason"]


def test_cisplatin_band_names_its_population_as_cell_lines(real_run):
    """Cisplatin's spread IS measured — across 34 immortalised cell lines, not
    across patients. It may be banded; it may not be called a between-patient
    bar."""
    iv = kill_fraction_interval(_with_compound(real_run, "cisplatin"))
    pot = iv["components"]["potency"]
    assert "CELL LINES" in pot["population"]
    assert "not patients" in pot["population"].lower()
    assert "NOT PATIENTS" in pot["caveat"]


def test_cell_state_span_does_not_transfer_to_another_section(real_run):
    """Three sections on three panels is a sample size of three and nothing
    transfers between them (results/cell_state_impact.md headline)."""
    iv = kill_fraction_interval(_with_compound(real_run, "temozolomide",
                                               section_id="xenium_io"))
    assert "cell_state_route" not in iv["components"]
    refusal = next(r for r in iv["refusals"]
                   if r["metric"] == "kill_fraction.epistemic")
    assert "does not transfer" in refusal["reason"]
    assert iv["covers_kinds"] == [BIOLOGICAL]


def test_no_span_at_all_gives_no_band(real_run):
    iv = kill_fraction_interval(_with_compound(real_run, "lomustine",
                                               section_id="xenium_io"))
    assert iv["status"] == "UNCHARACTERISED"
    assert iv["band"] is None
    assert iv["components"] == {}
    assert len(iv["refusals"]) == 2


@pytest.mark.parametrize("metric", sorted(UNCHARACTERISED_METRICS))
def test_uncharacterised_metrics_report_that_and_not_a_band(real_run, metric):
    rep = run_interval_report(real_run, point_metrics={metric: 1.23})
    rec = rep["uncharacterised"][metric]
    assert rec["band"] is None
    assert rec["status"] == "UNCHARACTERISED"
    assert rec["quotable"] is False
    assert rec["point"] == 1.23          # the point is still served, unadorned
    assert "no measured span" in rec["reason"].lower()


# --------------------------------------------------------------------------- #
# 5. What the band says about itself
# --------------------------------------------------------------------------- #
def test_absolute_kill_is_never_quotable(real_run):
    for iv in (kill_fraction_interval(real_run), kill_ceiling_interval(real_run)):
        assert iv["quotable"] is False
        assert "not a defensible output" in iv["not_quotable_reason"] or \
            "stop quoting" in iv["not_quotable_reason"]


def test_band_brackets_the_point_and_is_wide(real_run):
    iv = kill_fraction_interval(real_run)
    band, point = iv["band"], iv["point"]
    assert band["lo"] < point < band["hi"], (band, point)
    # The finding the feature exists to surface: this is not a tight interval.
    assert band["span_x"] > 5.0, band


def test_state_bar_is_one_sided_upward(real_run):
    """The shipping route sits at the BOTTOM of the measured route span, so the
    epistemic bar only ever widens upward. Getting this backwards would put the
    prediction's lower edge somewhere no route measured."""
    c = kill_fraction_interval(real_run)["components"]["cell_state_route"]
    assert c["one_sided"] is True
    assert c["factor_lo"] == 1.0
    assert c["lo"] == pytest.approx(real_run.kill_fraction)
    assert c["hi"] > c["lo"]


def test_potency_inversion_is_the_right_way_round(real_run):
    """A LOWER IC50 is MORE potent and must give a HIGHER kill. This inversion
    is the easiest thing in the feature to get backwards and the hardest to
    notice once it is."""
    from tmesim.validation.uncertainty import expected_kill_at_ic50
    span = POTENCY_SPAN["temozolomide"]
    assert expected_kill_at_ic50(real_run, span.lo) > \
        expected_kill_at_ic50(real_run, span.hi)


def test_kill_ceiling_carries_epistemic_only(real_run):
    """The ceiling saturates the Hill term, so intrinsic potency drops out of it
    entirely: a metric with an epistemic bar and a STRUCTURALLY ZERO biological
    one. It is the clearest demonstration that the two are different objects."""
    iv = kill_ceiling_interval(real_run)
    assert iv["covers_kinds"] == [EPISTEMIC]
    assert set(iv["components"]) == {"cell_state_route"}
    refusal = next(r for r in iv["refusals"]
                   if r["metric"] == "kill_ceiling.biological")
    assert "structural absence, not a data gap" in refusal["reason"]
    assert iv["band"]["lo"] == pytest.approx(F_STATE_BY_ROUTE["oxygen"])
    assert iv["band"]["hi"] == pytest.approx(F_STATE_BY_ROUTE["vessel_distance"])


def test_ceiling_discloses_when_its_own_point_falls_outside_the_band(real_run):
    """The ceiling band's endpoints are the RECORDED per-route measurement; the
    run's own F_state is recomputed from the run's own cells. When they differ
    the point can fall outside its own band, and that must be DISCLOSED rather
    than fixed by widening a measured endpoint to swallow it."""
    iv = kill_ceiling_interval(real_run)
    assert "point_inside_band" in iv
    if not iv["point_inside_band"]:
        assert "not moved to contain it" in iv["point_outside_band_disclosure"]
    # Either way the endpoints are the record's, untouched.
    assert iv["band"]["lo"] == pytest.approx(F_STATE_BY_ROUTE["oxygen"])
    assert iv["band"]["hi"] == pytest.approx(F_STATE_BY_ROUTE["vessel_distance"])


def test_report_states_what_it_does_not_cover(real_run):
    rep = run_interval_report(real_run)
    blob = " ".join(rep["not_covered"]).upper()
    for term in ("MODEL-FORM ERROR", "OUTCOME VALIDATION", "MGMT",
                 "CROSS-COMPOUND"):
        assert term in blob
    assert "REFUSED" in rep["cross_compound"]


def test_no_entry_point_takes_two_compounds():
    """Cross-compound ranking was withdrawn 2026-08-21, absolutely and in any
    form. Structural, not a promise: every entry point takes ONE FinishedRun,
    which carries ONE compound."""
    import inspect

    from tmesim.validation import uncertainty as U

    for fn in (U.kill_fraction_interval, U.kill_ceiling_interval,
               U.run_interval_report, U.expected_kill_at_ic50):
        params = inspect.signature(fn).parameters
        assert sum(1 for p in params.values()
                   if p.annotation in (FinishedRun, "FinishedRun")) <= 1
        assert not any("compound" in n for n in params), fn.__name__


# --------------------------------------------------------------------------- #
# 6. Cost — measured, not guessed (PRE-REGISTERED: default-on iff < 1% of a run)
# --------------------------------------------------------------------------- #
#: MEASURED, in this checkout, over the 86 stored runs in ``app/runs`` that
#: carry a runtime: median 39.91 s (min 1.20, max 240.06). PR-2 asks for < 1% of
#: the run the band annotates, so against the shipping default that is a budget
#: of ~0.40 s. The band's cost is FLAT (it solves nothing), so the assertion
#: below is on the absolute cost, with a 15x regression margin over the ~3 ms it
#: actually takes -- anything that made it re-solve or re-reconstruct would blow
#: straight through it.
BAND_COST_CEILING_S = 0.05


def test_band_is_cheap_enough_to_be_default_on(real_micro, real_drug):
    """PRE-REGISTERED PR-2: default-on iff the band costs < 1% of its run.

    REPORTED, NOT TUNED: this test's own run is the smallest one the suite can
    build (dx = 80 um, a 30-minute window) and takes ~0.23 s, against which the
    band's flat ~3 ms is ~1.2% -- PR-2 FAILS on a run that small, and that is
    stated in results/uncertainty_band_2026-09-20.md rather than hidden by
    picking a kinder denominator. On the shipping configuration the same 3 ms is
    0.011% of a 40 s run and PR-2 passes with three orders of magnitude to
    spare. The band is default-on-cheap for every run anyone actually makes and
    is not for a run shorter than about half a second.
    """
    t0 = time.perf_counter()
    res = simulate_response(real_micro, SimConfig(drug=real_drug, pde=PDE, seed=0))
    run_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    run = finished_run_from_result(real_micro, res, real_drug,
                                   compound="temozolomide",
                                   section_id=FLAGSHIP_SECTION)
    rep = run_interval_report(run)
    band_s = time.perf_counter() - t0

    assert rep["intervals"]["kill_fraction"]["status"] == "BANDED"
    assert band_s < BAND_COST_CEILING_S, (
        f"the band must stay milliseconds-flat: measured {band_s:.4f}s "
        f"(ceiling {BAND_COST_CEILING_S}s). This run took {run_s:.2f}s, i.e. "
        f"{100 * band_s / run_s:.3f}% -- see the docstring for why the ratio on "
        f"a 0.2 s run is not the PR-2 test.")


# --------------------------------------------------------------------------- #
# 7. The route
# --------------------------------------------------------------------------- #
fastapi = pytest.importorskip("fastapi", reason="the app/server backend needs fastapi")
pytest.importorskip("httpx", reason="starlette's TestClient needs httpx")

from fastapi import FastAPI  # noqa: E402
from fastapi.testclient import TestClient  # noqa: E402

from app.server.jobs import JobManager  # noqa: E402
from app.server.uncertainty_api import build_uncertainty_router  # noqa: E402


def _client(root):
    jobs = JobManager(root=str(root))
    app = FastAPI()
    app.include_router(build_uncertainty_router(jobs))
    return TestClient(app), jobs


def test_spans_route_serves_the_registry_and_the_distinction(tmp_path):
    client, jobs = _client(tmp_path)
    try:
        body = client.get("/api/uncertainty/spans").json()
        assert set(body["kinds"]) == {BIOLOGICAL, EPISTEMIC}
        assert body["kinds"][BIOLOGICAL]["does_a_better_model_shrink_it"] == "no"
        assert body["kinds"][EPISTEMIC]["does_a_better_model_shrink_it"] == "yes"
        assert "worse than no feature" in body["the_distinction"]
        assert body["potency_spans"]["doxorubicin"]["bandable"] is False
        assert body["potency_spans"]["temozolomide"]["bandable"] is True
        assert body["cell_state_factor"]["lo"] == 1.0
        assert "REFUSED" in body["cross_compound"]
        # No second copy of a span on the wire.
        assert body["f_state_by_route"] == F_STATE_BY_ROUTE
    finally:
        jobs.shutdown()


def test_unknown_run_is_404_and_a_traversal_never_reaches_disk(tmp_path):
    client, jobs = _client(tmp_path)
    try:
        assert client.get("/api/runs/r_99999/uncertainty").status_code == 404
        for bad in ("..%2F..%2Fetc%2Fpasswd", "a%2Fb", "."):
            r = client.get(f"/api/runs/{bad}/uncertainty")
            assert r.status_code == 404, (bad, r.status_code)
    finally:
        jobs.shutdown()


def _copy_a_finished_run(tmp_path):
    """Copy one REAL finished run from this checkout into an isolated store."""
    src_root = os.path.join(_REPO, "app", "runs")
    if not os.path.isdir(src_root):
        return None
    for name in sorted(os.listdir(src_root), reverse=True):
        d = os.path.join(src_root, name)
        doc_path = os.path.join(d, "result.json")
        if not os.path.isfile(doc_path):
            continue
        with open(doc_path) as fh:
            doc = json.load(fh)
        if doc.get("status") != "done" or not (doc.get("spec") or {}).get("drug"):
            continue
        if not os.path.isfile(os.path.join(d, "field_drug_auc.bin")):
            continue
        shutil.copytree(d, os.path.join(str(tmp_path), name))
        return name
    return None


def test_run_route_bands_a_real_finished_run(tmp_path):
    run_id = _copy_a_finished_run(tmp_path)
    if run_id is None:
        pytest.skip("no finished run with a drug payload in app/runs/")
    client, jobs = _client(tmp_path)
    try:
        body = client.get(f"/api/runs/{run_id}/uncertainty").json()
        assert body["resimulated"] is False
        kf = body["intervals"]["kill_fraction"]
        assert kf["reproduction"]["passed"] is True, kf["reproduction"]
        assert kf["quotable"] is False
        assert kf["status"] in ("BANDED", "BANDED_PARTIAL")
        # INCLUSIVE, and a degenerate side must be EXPLAINED rather than
        # assumed away.
        #
        # This asserted strict interiority and passed only because
        # `_copy_a_finished_run` takes the newest run and the newest run
        # happened to be temozolomide, which has both components. On 2026-09-20
        # a doxorubicin run became the newest and it went red at
        # `0.4769 < 0.4769` -- correctly. Doxorubicin's potency span is REFUSED
        # (its 0.1-1 uM is a quoted range with no population, n or quantile),
        # so the only surviving component is the epistemic cell-state route,
        # which is one-sided upward. `lo == point` is then the honest band, and
        # `BANDED_PARTIAL` is the status that says so.
        #
        # Widening the endpoints to make a strict inequality true would be
        # inventing uncertainty that was never characterised -- the exact thing
        # this module refuses to do elsewhere.
        lo, hi, point = kf["band"]["lo"], kf["band"]["hi"], kf["point"]
        assert lo <= point <= hi, (
            f"the point estimate {point} is outside its own band "
            f"[{lo}, {hi}]")
        assert lo < hi, "a band with no width is not a band"
        if lo == point or hi == point:
            assert kf["status"] == "BANDED_PARTIAL", (
                "a one-sided band must be labelled BANDED_PARTIAL, so a "
                f"reader knows a component was refused; got {kf['status']!r}")
        assert_kinds_separated(kf)
        # The whole payload must survive the JSON round trip the client makes.
        json.dumps(body)
    finally:
        jobs.shutdown()


def test_run_route_refuses_a_run_whose_payloads_are_gone(tmp_path):
    """PLANT: delete the exposure field the band re-evaluates. The route must
    refuse with the missing payload NAMED — not substitute another run's field,
    and not 500."""
    run_id = _copy_a_finished_run(tmp_path)
    if run_id is None:
        pytest.skip("no finished run with a drug payload in app/runs/")
    os.remove(os.path.join(str(tmp_path), run_id, "field_drug_auc.bin"))
    client, jobs = _client(tmp_path)
    try:
        r = client.get(f"/api/runs/{run_id}/uncertainty")
        assert r.status_code == 200           # a refusal is a result
        body = r.json()
        assert body["status"] == "REFUSED"
        assert body["intervals"] == {}
        assert body["missing_payloads"] == ["auc"]
    finally:
        jobs.shutdown()


def test_calibrate_alpha_frozen_at_end():
    assert calibrate_alpha() == ALPHA_FROZEN, calibrate_alpha()
