"""THE FIRST TEN MINUTES — regressions found by using the Studio as a stranger.

These are UI tests, and they exist because the defects they pin were not
findable by reading the code. Every one of them was found by opening the Studio
with no prior knowledge and reading it top to bottom, which is the one thing
nobody who works on it can do twice.

What is pinned here:

1.  **Absolute kill is not the first thing on the rail.** Until 2026-08-22 the
    READ rail opened on ``Kill fraction`` rendered at the largest type size the
    panel has, four sections above the run's own caveat reading "ABSOLUTE KILL
    MUST NOT BE QUOTED FROM THIS SECTION UNDER ANY STATE ROUTE". A reader who
    has never seen this product reads top-down and quotes the one number the
    project forbids. This is the project's own hard rule failing inside its own
    UI, so it gets a test rather than a comment.

2.  **The long caveat is restructured losslessly.** ``scale_invariance`` is one
    ~5,900-character string and it is split on its own axis markers so each
    swept axis gets its own paragraph. A restructuring that DROPS text would be
    a caveat silently shortened, which is the one thing this repo may never do
    to a disclosure. The test asserts character-for-character reassembly
    against the REAL backend strings, both modes, not against a fixture.

3.  **The finding is stated in counts with a visible denominator.** The claim
    "99.8% under-dosed" has been mis-stated inside this project by attaching it
    to the wrong denominator (survivors, n=21,140, versus viable tumour cells,
    n=2,008 off the reach ladder). A percentage on screen without its
    denominator beside it is how that happens, so the denominator is pinned.

The browser is driven against a plain static server over ``app/web``, NOT
against a running Codon backend: these assertions are about rendering, and a
test that needs the physics stack to check a heading is a test that will be
skipped. Run payloads are FIXTURES built in this file — see ``_run_fixture``.
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

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(REPO, "app", "web")

playwright_api = pytest.importorskip(
    "playwright.sync_api", reason="first-session UI tests need playwright"
)


# --------------------------------------------------------------------------- #
# a static server over app/web, so the ES modules load with real module
# semantics (file:// blocks them) and nothing needs the Codon backend
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def web_server():
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=WEB)

    class Quiet(socketserver.TCPServer):
        allow_reuse_address = True

        def handle_error(self, request, client_address):  # noqa: D102
            pass

    with Quiet(("127.0.0.1", 0), handler) as httpd:
        port = httpd.server_address[1]
        t = threading.Thread(target=httpd.serve_forever, daemon=True)
        t.start()
        try:
            yield f"http://127.0.0.1:{port}"
        finally:
            httpd.shutdown()


@pytest.fixture(scope="module")
def page(web_server):
    with playwright_api.sync_playwright() as p:
        browser = p.chromium.launch()
        pg = browser.new_page()
        errors: list[str] = []
        pg.on("pageerror", lambda e: errors.append(str(e)))
        # Any page under the served root works; the modules are imported by URL.
        pg.goto(f"{web_server}/index.html", wait_until="domcontentloaded")
        pg.wait_for_timeout(300)
        pg._codon_errors = errors  # type: ignore[attr-defined]
        yield pg
        browser.close()


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
def _run_fixture(**over):
    """A RunResult shaped like the real thing.

    The numbers are the flagship section's, so the assertions below read like
    the screen a first-time user actually meets. They are a FIXTURE and not a
    claim: nothing in this file is a source of truth about the physics, and
    every figure here is also carried by the run record it was copied from.
    """
    run = {
        "run_id": "r_fixture",
        "status": "done",
        "spec": {"section": "xenium_319", "drug": "temozolomide",
                 "sensitivity_mode": "uniform"},
        "metrics": {
            "kill_fraction": 0.11221232991768859,
            "mean_tumor_auc": 4225543.987491742,
            "therapeutic_coverage_1x": 0.0,
            "penetration_length_um": 278.4919817884945,
        },
        # Every figure below is r_00090's, the shipped default RunSpec on the
        # flagship section: plasma_peak = 50 uM against a measured temozolomide
        # IC50 of 606 uM, so the ratio the verdict turns on is 0.083 and
        # therapeutic_coverage_1x is 0.0 above.
        "extra_metrics": {"E_p10_uM": 135.42888240722527,
                          "E_p50_uM": 185.38644552397483,
                          "E_p90_uM": 270.1005830234793,
                          "plasma_peak_uM": 50.0,
                          "plasma_mean_uM": 31.09499320754238,
                          "median_E_over_ic50": 0.30589466787932706,
                          "p10_E_over_ic50": 0.22346279356151275,
                          "therapeutic_coverage_0_5x": 0.049092894338988745,
                          "therapeutic_coverage_2x": 0.0,
                          # Supplied by app/server/jobs.py so the panel reads
                          # rather than divides. The PEAK ratio is the bound
                          # (C <= Cp); the mean one is context and inflates
                          # with the window length, not with the violation.
                          "plasma_peak_over_ic50": 0.08250189678505047,
                          "tissue_median_over_plasma_peak": 3.7077289104794966,
                          "tissue_median_over_plasma_mean": 5.961629551941236},
        "drug_params": {"name": "temozolomide", "ic50_uM": 606.0466722391783,
                        "bbb_penetrance": 0.12537158358322562},
        "report": {"source": "xenium", "n_cells": 40887, "n_tumor": 23812},
        "provenance": {
            "caveats": [], "claim_set": None,
            "ic50_source": "GDSC measured (n=51)",
            "ic50_n_cell_lines": 51,
            "ic50_confidence": "measured",
            # The SHIPPED default: the one-way leak, which is unbounded above.
            "delivery_law": {"transvascular_exchange": False,
                             "law": "S = PS*Cp",
                             "alternative_law": "S = PS*(Cp - C)"},
        },
        "mechanism": {
            "n_tumor": 23812,
            "n_survived": 21140,
            "counts": {
                "KILLED": 2672, "SURVIVED_UNREACHED": 46,
                "SURVIVED_UNDERDOSED": 21094, "SURVIVED_RESISTANT": 0,
                "SURVIVED_STOCHASTIC": 0, "UNATTRIBUTED": 0,
                "ALREADY_DEAD": 0, "NOT_TUMOR": 17075,
            },
            "survival_breakdown": {
                "delivery": 0.0021759697256386, "dose": 0.9978240302743614,
                "resistance": 0.0, "chance": 0.0, "unattributed": 0.0,
            },
            "survival_breakdown_note": "fractions of SURVIVING tumour cells, "
                                       "so they sum to 1.",
            "exposure_shortfall_fold": 3.2690991540738192,
            "exposure_axis_note": "UNREACHED and UNDERDOSED are two bands on "
                                  "ONE axis, E/IC50.",
            "thresholds": {"unreached_ratio": 0.1, "adequate_ratio": 1.0},
            "joint_limited": {
                "count": 19603,
                "fraction_of_exposure_limited": 0.9272942289498581,
                # The backend's own wording. Neither panel used to print it,
                # so both paraphrased a string nobody re-read.
                "note": "exposure-limited survivals that are ALSO at or below "
                        "the resistance threshold.",
            },
            "reconciliation": {"bit_exact": True},
        },
    }
    run.update(over)
    return run


def _backend_notes():
    """The REAL ``scale_invariance`` strings, both modes, off the backend.

    Imported rather than pasted: a copy would keep passing after the backend
    reworded the axis markers the split relies on, which is precisely the
    regression this test is for.
    """
    from app.server import reach as R  # noqa: PLC0415

    return {"finite_window": R.SCALE_INVARIANCE_NOTE,
            "steady_state": R.STEADY_SCALE_INVARIANCE_NOTE}


# --------------------------------------------------------------------------- #
# 1. the long caveat is restructured, never shortened
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mode", ["finite_window", "steady_state"])
def test_scale_invariance_split_is_lossless(page, web_server, mode):
    """``splitAxes`` reassembles to the original, character for character.

    This is the whole licence for touching a caveat's layout. Restructure,
    group, progressively disclose — never remove. If the reassembly differs by
    one character the restructuring has eaten a disclosure and this fails.
    """
    note = _backend_notes()[mode]
    out = page.evaluate(
        """async ([base, s]) => {
             const m = await import(base + '/js/ui/reach.js');
             const parts = m.splitAxes(s);
             return {
               n: parts.length,
               joined: parts.map(p => (p.head || '') + p.body).join(''),
               heads: parts.map(p => p.head).filter(Boolean),
             };
           }""",
        [web_server, note],
    )
    assert out["joined"] == note, (
        f"splitAxes lost or altered text in {mode} mode: "
        f"{len(note)} chars in, {len(out['joined'])} chars out"
    )
    # And it actually split — a lossless no-op would pass the check above while
    # leaving the wall of prose exactly as unreadable as it was.
    assert out["n"] >= 4, (
        f"{mode}: expected a lead-in plus three swept axes, got {out['n']} "
        f"parts (heads={out['heads']})"
    )
    joined_heads = " ".join(out["heads"])
    for axis in ("PLASMA/DOSE SCALE", "TISSUE CAPACITY", "TRANSVASCULAR PENETRANCE"):
        assert axis in joined_heads, f"{mode}: axis {axis!r} was not given its own head"


def test_split_axes_passes_through_unmarked_text(page, web_server):
    """No markers -> one part carrying everything. The degradation has to be
    'render exactly as before', never 'render some of it'."""
    s = "A caveat with no numbered axes in it at all. It must survive whole."
    out = page.evaluate(
        """async ([base, s]) => {
             const m = await import(base + '/js/ui/reach.js');
             const parts = m.splitAxes(s);
             return { n: parts.length, joined: parts.map(p => (p.head || '') + p.body).join('') };
           }""",
        [web_server, s],
    )
    assert out["n"] == 1
    assert out["joined"] == s


# --------------------------------------------------------------------------- #
# 2. the finding, in counts, with its denominator on screen
# --------------------------------------------------------------------------- #
_RENDER_JS = """async ([base, mod, run]) => {
  const m = await import(base + '/js/ui/' + mod);
  const p = mod === 'finding.js' ? m.createFindingPanel() : m.createReachPanel();
  p.render(run);
  document.body.appendChild(p.el);
  const t = p.el.innerText;
  p.el.remove();
  return t;
}"""


@pytest.fixture()
def finding_text(page, web_server):
    return page.evaluate(_RENDER_JS, [web_server, "finding.js", _run_fixture()])


def test_finding_states_counts_and_denominator(finding_text):
    """Both counts AND the denominator they are shares of.

    "99.8% under-dosed" has already been mis-attached to the wrong denominator
    inside this project — survivors (21,140) versus viable tumour cells (2,008)
    off the reach ladder. A percentage whose denominator is not beside it on
    screen is how a reader repeats that, so the denominator is not optional.
    """
    assert "21,140" in finding_text, "the survivor denominator is not on screen"
    assert "23,812" in finding_text, "the tumour-cell total is not on screen"
    assert "21,094" in finding_text, "the under-dosed COUNT is not on screen"
    assert "46" in finding_text, "the unreached COUNT is not on screen"
    assert "99.8%" in finding_text
    assert "0.2%" in finding_text


def test_finding_names_what_would_refute_it(finding_text):
    """The claim replaced a compound ranking BECAUSE it is falsifiable. A
    falsifiable claim that never names its own test is being sold as an
    unfalsifiable one."""
    assert "WHAT WOULD REFUTE THIS" in finding_text
    lowered = finding_text.lower()
    assert "refute" in lowered
    # It must name a MEASUREMENT, not an outcome.
    assert any(w in lowered for w in ("mass spectrometry", "microdialysis"))


def test_finding_carries_the_dose_counterfactual_limit(finding_text):
    """"Arrived but too weak" must not be readable as "more drug would work".
    92.7% of the exposure-limited survivors are also below the resistance
    threshold, and that limit belongs beside the number, not in a block."""
    assert "19,603" in finding_text
    assert "92.7%" in finding_text
    assert "more drug would" in finding_text.lower()


def test_finding_does_not_headline_absolute_kill(finding_text):
    """The kill count may appear as the denominator arithmetic. It may not be
    the panel's headline, and it must arrive already carrying its bound."""
    assert "DO NOT" in finding_text or "not a result" in finding_text.lower()
    # The bound travels with it, on the same panel, in the same breath.
    assert "0.3252" in finding_text and "1.77x" in finding_text
    # And the panel leads on the survivor arithmetic, not on a kill rate: the
    # tumour total is stated before the killed count anywhere on the panel.
    assert finding_text.index("23,812") < finding_text.index("2,672")


def test_finding_refuses_a_drugless_run(page, web_server):
    """A reconstruction-only run has no delivery-versus-potency finding. It
    must say so and print nothing, rather than rendering an empty split."""
    run = _run_fixture(spec={"section": "xenium_319", "drug": None},
                       mechanism=None)
    text = page.evaluate(_RENDER_JS, [web_server, "finding.js", run])
    assert "no compound" in text.lower()
    assert "99.8" not in text and "21,094" not in text


# --------------------------------------------------------------------------- #
# 3. the hard rule: absolute kill is not where the eye lands
# --------------------------------------------------------------------------- #
def test_kill_fraction_carries_an_inline_hard_stop():
    """The metric table declares a hard stop for absolute kill.

    Read out of the source rather than the DOM because it is a declaration:
    inspect.js renders whatever carries one, so pinning the declaration pins
    every render site at once.
    """
    src = open(os.path.join(WEB, "js", "ui", "provenance.js"), encoding="utf-8").read()
    block = src[src.index("kill_fraction: {"):]
    block = block[: block.index("},")]
    assert "hard_stop:" in block, (
        "kill_fraction lost its inline hard stop. Absolute kill is not a "
        "defensible output of this project and must not be quotable off the "
        "rail without its bound attached."
    )
    assert "upper bound" in block


def test_no_metric_is_rendered_as_the_lead():
    """No metric gets the largest type size on the rail.

    The first key of METRICS is ``kill_fraction``, so ANY "emphasise the first
    metric" rule puts an absolute kill percentage at the top of the READ rail
    at --fs-xl. That is the defect this file was opened for. The rule is gone;
    this stops it coming back under a different name.
    """
    src = open(os.path.join(WEB, "js", "ui", "inspect.js"), encoding="utf-8").read()
    # Strip comments: this repo documents a withdrawn decision by quoting it.
    code = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    code = re.sub(r"^\s*//.*$", "", code, flags=re.M)
    assert "'.lead'" not in code and '".lead"' not in code
    assert "'div.metric' + (first" not in code
    assert re.search(r"\.metric\.lead", code) is None


def test_finding_section_precedes_metrics_on_the_rail():
    """The rail opens on the finding, not on four absolute numbers."""
    src = open(os.path.join(WEB, "js", "ui", "inspect.js"), encoding="utf-8").read()
    body = src[src.index("const railBody ="):]
    i_find = body.index("'The finding'")
    i_metrics = body.index("'Metrics'")
    i_caveats = body.index("'Caveats'")
    assert i_find < i_metrics, "Metrics is rendered above the finding"
    assert i_find < i_caveats


# --------------------------------------------------------------------------- #
# 4. the two objects are named as two questions, where the map starts
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "mode,expect",
    [("steady_state", "GEOMETRY PERMIT"), ("finite_window", "WINDOW ACHIEVE")],
)
def test_reach_panel_names_which_question_before_the_evidence(
    page, web_server, mode, expect
):
    """A reader must be told which of the two questions they are looking at the
    answer to BEFORE the falloff curve, the ladder and four capacity caveats,
    all of which mean different things in the two modes."""
    from app.server import reach as R  # noqa: PLC0415

    run = _run_fixture()
    run["reach"] = {
        "null_reason": None,
        "headline_x": 0.1,
        "headline_under_reached_fraction": 0.00049800796812749,
        "curve": {"x": [0.01, 0.1, 1.0],
                  "under_reached_fraction": [0.0005, 0.0005, 0.3934]},
        "falloff": {"distance_um": [5.0, 55.0], "normalised_mean_auc": [1.0, 0.6],
                    "bin_um": 10.0},
        "E_ref_uM_s": 4708982.107545392,
        "E_ref_units": "uM" if mode == "steady_state" else "uM*s",
        "E_ref_source": "penetration_length_fit(...)['drug']['perivascular_mean']",
        "window_mode": mode,
        "window_note": "note",
        "units": {"E_ref": "uM" if mode == "steady_state" else "uM*s"},
        "convention": R.CONVENTION_NOTE,
        "scale_invariance": _backend_notes()[mode],
    }
    text = page.evaluate(_RENDER_JS, [web_server, "reach.js", run])
    assert expect in text, f"{mode}: the panel does not name its own question"
    assert "different questions" in text.lower()
    # BEFORE the evidence it governs, not after it. Anchored on the delivery
    # spread, which is the first number of that half.
    assert text.index(expect) < text.index("DELIVERY SPREAD"), (
        "the object banner is rendered below the delivery evidence it governs"
    )
    # ...and before the capacity caveats, which say opposite things in the two
    # modes and are unreadable without knowing which mode you are in.
    assert text.index(expect) < text.index("CAPACITY")


def test_reach_panel_does_not_call_a_concentration_an_exposure(page, web_server):
    """In steady-state mode the reference is a CONCENTRATION (uM); the payload
    field is still named ``E_ref_uM_s``, which is simply wrong there. The panel
    reads the units, never the field name."""
    from app.server import reach as R  # noqa: PLC0415

    run = _run_fixture()
    run["reach"] = {
        "null_reason": None, "headline_x": 0.1,
        "curve": {"x": [0.1], "under_reached_fraction": [0.0005]},
        "falloff": {"distance_um": [5.0], "normalised_mean_auc": [1.0], "bin_um": 10.0},
        "E_ref_uM_s": 12.5, "E_ref_units": "uM",
        "E_ref_source": "penetration_length_fit(C_ss, ...)",
        "window_mode": "steady_state", "units": {"E_ref": "uM"},
        "convention": R.CONVENTION_NOTE,
        "scale_invariance": R.STEADY_SCALE_INVARIANCE_NOTE,
    }
    text = page.evaluate(_RENDER_JS, [web_server, "reach.js", run])
    assert "settled drug concentration" in text
    assert "reference concentration" in text, (
        "the reference is glossed as an 'exposure' in a mode where it is a "
        "concentration"
    )
    # The value is printed with its own units, read off E_ref_units.
    assert re.search(r"1\.250e1\s*uM\b", text), text[:400]
    # And the exposure wording must NOT survive into this mode.
    assert "reference exposure" not in text


def test_no_console_errors_on_module_load(page):
    assert page._codon_errors == []  # type: ignore[attr-defined]


# --------------------------------------------------------------------------- #
# 5. the sign of the verdict is not silently an input
#
# THE DEFECT. The shipped default RunSpec carries plasma_peak = 50 uM and
# temozolomide's measured IC50 is 606 uM, so plasma_peak / IC50 = 0.083 and
# therapeutic_coverage_1x comes back 0.0. The rail then reported the compound as
# DOSE-LIMITED without once printing the dose it was given. The first question a
# pharmacologist asks is the right one — "you told it the peak was 12x under
# IC50 and it told you the dose was too low; what did the tissue contribute?" —
# and until 2026-08-22 no panel answered it.
#
# js/ui/reach.js had answered the structurally identical question for
# bevacizumab since it was written (bbb_penetrance = 0 gives an identically zero
# source term, and the panel says THIS IS AN INPUT, NOT A FINDING). These pin
# the same treatment on the dose axis.
# --------------------------------------------------------------------------- #
@pytest.fixture()
def reach_run():
    """A fixture with a reach block, for the panel that restates the verdict."""
    from app.server import reach as R  # noqa: PLC0415

    run = _run_fixture()
    run["reach"] = {
        "null_reason": None, "headline_x": 0.1,
        "curve": {"x": [0.1, 1.0], "under_reached_fraction": [0.0005, 0.3934]},
        "falloff": {"distance_um": [5.0, 55.0], "normalised_mean_auc": [1.0, 0.6],
                    "bin_um": 10.0},
        "E_ref_uM_s": 4708982.107545392, "E_ref_units": "uM*s",
        "E_ref_source": "penetration_length_fit(...)['drug']['perivascular_mean']",
        "window_mode": "finite_window", "units": {"E_ref": "uM*s"},
        "convention": R.CONVENTION_NOTE,
        "scale_invariance": R.SCALE_INVARIANCE_NOTE,
    }
    return run


def _render(page, web_server, mod, run):
    return page.evaluate(_RENDER_JS, [web_server, mod, run])


def test_the_ratio_itself_is_on_screen_beside_the_verdict(page, web_server):
    """The RATIO, with BOTH of its terms, so the reader can do the division.

    Not a grade, not a band, not a warning triangle: the two numbers and their
    quotient. A reader who disagrees with this panel has to be able to disagree
    with an arithmetic they can see.
    """
    text = _render(page, web_server, "finding.js", _run_fixture())
    assert "plasma_peak / IC50" in text, "the ratio is not named on screen"
    assert "50 uM" in text, "the entered plasma peak is not printed"
    assert "606.0 uM" in text, "the IC50 the ratio is taken against is not printed"
    assert "0.083" in text, "the ratio itself is not printed"
    # And the half that is NOT an input is attributed, or the whole thing reads
    # as a ratio between two numbers the operator chose.
    assert "GDSC" in text


def test_a_sub_ic50_dose_is_named_as_an_input(page, web_server):
    """Below 1, the panel says so in the same words it already uses for
    bbb_penetrance = 0: this part of the verdict was handed to the model."""
    text = _render(page, web_server, "finding.js", _run_fixture())
    assert "PART OF THIS VERDICT IS AN INPUT" in text
    assert "RunSpec" in text
    # The fold, so "below 1" is not the only quantity a reader carries away.
    assert "12.1x" in text


def test_dose_input_disclosure_invents_no_threshold():
    """No tuned cut anywhere in the module.

    A threshold chosen now would be chosen KNOWING what the shipped default
    does, which is post-hoc selection wearing a convention's clothes — this
    product removed a single-threshold headline for exactly that reason. The
    only comparison permitted is against 1.0, which is not a choice: it is the
    definitional point where the entered peak equals the half-effect
    concentration.
    """
    src = open(os.path.join(WEB, "js", "ui", "doseinput.js"), encoding="utf-8").read()
    code = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    code = re.sub(r"^\s*//.*$", "", code, flags=re.M)
    cmps = re.findall(r"(?:ratio|foldUnder|coverage1x|medianOverIc50|p10OverIc50)"
                      r"\s*(?:<=?|>=?)\s*([0-9.]+)", code)
    assert cmps, "the ratio is never compared to anything — has the branch gone?"
    assert set(cmps) == {"1"}, (
        f"doseinput.js compares a reported quantity against {sorted(set(cmps))}. "
        "Only 1 is permitted: any other number is a threshold this module chose "
        "after seeing what the data does."
    )


@pytest.mark.parametrize(
    "law,expect,forbid",
    [
        # gradient exchange: C_ss <= Cp by construction, so a sub-IC50 peak makes
        # the verdict arithmetic and the panel must say the RunSpec settled it.
        ({"transvascular_exchange": True, "law": "S = PS*(Cp - C)"},
         "settled by the RunSpec", "UNBOUNDED ABOVE"),
        # one-way leak (the SHIPPED default): unbounded above, so the arithmetic
        # did NOT force the answer and the panel must not claim it did.
        ({"transvascular_exchange": False, "law": "S = PS*Cp"},
         "did not force the answer", "settled by the RunSpec"),
        # absent: a third state. Neither is assumed.
        (None, "CANNOT be read off it", "settled by the RunSpec"),
    ],
)
def test_what_follows_from_a_sub_ic50_dose_is_read_off_the_delivery_law(
    page, web_server, law, expect, forbid
):
    """The consequence of a ratio below 1 is NOT law-independent.

    Under ``S = PS*(Cp - C)`` tissue is bounded by plasma, so coverage at 1x
    IC50 is 0.0 by arithmetic and the verdict's sign is the spec's. Under the
    shipped one-way leak ``S = PS*Cp`` the steady state is unbounded above, so
    it is not — and on the flagship the map really does run above plasma. A
    panel that stated either consequence unconditionally would be wrong half
    the time, so the branch is pinned in all three states.
    """
    run = _run_fixture()
    prov = dict(run["provenance"])
    if law is None:
        prov.pop("delivery_law", None)
    else:
        prov["delivery_law"] = law
    run["provenance"] = prov
    text = _render(page, web_server, "finding.js", run)
    assert expect in text, f"delivery_law={law}: panel does not say {expect!r}"
    assert forbid not in text, (
        f"delivery_law={law}: panel says {forbid!r}, which is only true under "
        "the other law"
    )


def test_the_unbounded_law_is_quantified_not_merely_named(page, web_server):
    """Under the leak the map exceeds the plasma that fed it, and the panel puts
    a number on it rather than leaving 'unbounded above' as an abstraction.

    185.4 uM in the median tumour cell against a 31.1 uM plasma mean over the
    same window. That is the source term's known defect, NOT the tissue
    concentrating drug, and the panel has to say which — otherwise the single
    most striking number on it reads as a biological finding.
    """
    text = _render(page, web_server, "finding.js", _run_fixture())
    # THE BOUND IS AGAINST THE PEAK. C <= Cp is a statement about the highest
    # concentration the blood ever reached. The plasma MEAN is the intuitive
    # comparison and the wrong one: it falls as the window lengthens, so the
    # ratio inflates with the schedule rather than with the violation. Measured
    # on a short window the mean reading is 0.848 and the peak reading 0.836 —
    # both under 1 while the law is exactly as unbounded — so a panel reasoning
    # off the mean would report the bound satisfied there and violated here.
    assert "C <= Cp" in text, "the bound is not named"
    assert "185.4 uM" in text and "3.71x the peak" in text
    assert "active uptake" in text, (
        "the excess is printed without saying it is the law rather than the "
        "tissue"
    )
    # The mean is allowed, as context, and must be marked as not the bound.
    assert "5.96x the plasma MEAN" in text
    assert "Do not read that as the bound" in text


def test_a_dose_at_or_above_ic50_is_reported_as_not_forcing_the_verdict(
    page, web_server
):
    """The complement, stated just as plainly.

    Without this branch the block is a standing apology for the shipped default
    rather than a measurement of one run, and a reader learns nothing from its
    absence.
    """
    run = _run_fixture()
    # Both, as the backend writes them: it derives the ratio from the same peak
    # it reports. Setting only one is a payload that cannot exist.
    run["extra_metrics"] = dict(run["extra_metrics"], plasma_peak_uM=700.0,
                                plasma_peak_over_ic50=700.0 / 606.0466722391783)
    text = _render(page, web_server, "finding.js", run)
    assert "DOES NOT FORCE THIS VERDICT" in text
    assert "PART OF THIS VERDICT IS AN INPUT" not in text


def test_no_ratio_is_shown_when_the_payload_does_not_carry_one(page, web_server):
    """A half-known ratio is not a weaker disclosure, it is a different and
    unfounded one. Nothing is substituted for the missing half."""
    run = _run_fixture()
    run["drug_params"] = {"name": "temozolomide"}  # no ic50_uM
    text = _render(page, web_server, "finding.js", run)
    assert "plasma_peak / IC50" not in text
    assert "0.083" not in text
    # ...and the panel still renders everything else it always did.
    assert "21,094" in text


def test_the_ratio_reaches_the_panel_that_prints_coverage(page, web_server, reach_run):
    """`therapeutic_coverage_1x` is printed at the BOTTOM of the reach panel,
    three screens below the verdict. "0.0% of the tumour reached IC50" is a
    complete sentence and reads as a result, so the ratio that constrains it has
    to be attached where it is printed and not only where the verdict is."""
    text = _render(page, web_server, "reach.js", reach_run)
    assert "PART OF THIS VERDICT IS AN INPUT" in text, "not beside the verdict"
    assert "THIS LADDER STARTS BELOW ITS OWN TOP RUNG BY INPUT" in text
    # beside the coverage table, not above the whole panel
    assert text.index("THIS LADDER STARTS") > text.index("0.5x IC50")


# --------------------------------------------------------------------------- #
# 6. deduplicating the AUTHORING never weakens the DISCLOSURE
# --------------------------------------------------------------------------- #
def test_joint_limited_caveat_is_identical_on_both_panels(page, web_server, reach_run):
    """The verdict is stated twice on the rail, so its limits are too. They used
    to be WRITTEN twice, and had already drifted: finding.js carried the cell
    count and reach.js did not, so the same caveat was materially weaker on one
    of the two screens that states the claim it qualifies.

    Both now render js/ui/verdict.js. If one of them loses the count again, or
    the bodies diverge, this fails.
    """
    finding = _render(page, web_server, "finding.js", _run_fixture())
    reach = _render(page, web_server, "reach.js", reach_run)
    body = ("also at or below the modelled resistance threshold")
    for name, text in (("finding", finding), ("reach", reach)):
        assert "19,603" in text, f"{name}: the joint-limited COUNT is missing"
        assert "92.7%" in text, f"{name}: the joint-limited share is missing"
        assert body in text, f"{name}: the joint-limited body has drifted"
        assert "more drug would fix it" in text, f"{name}: the consequence is gone"
    # And the backend's own wording, which neither panel used to print.
    assert "backend ·" in finding and "backend ·" in reach


# --------------------------------------------------------------------------- #
# 7. the rail is navigable by heading
# --------------------------------------------------------------------------- #
def test_every_rail_section_title_is_a_heading():
    """Measured 2026-08-22: ZERO elements matching h1..h6 in the READ rail, over
    ~30,000 characters of disclosure in one scroller. A screen-reader user had
    no heading list, so the only route to "Delivery or dose" was to read from
    the top — and this rail's whole argument is that limits sit beside the
    numbers they qualify, which only pays off if you can reach the place.

    Pinned on the SOURCE because it is a construction rule: any new section
    written as `h('summary', 'Title'` puts a bare text node back in.
    """
    src = open(os.path.join(WEB, "js", "ui", "inspect.js"), encoding="utf-8").read()
    bare = re.findall(r"h\('summary',\s*'[^']+'", src)
    assert not bare, (
        f"{len(bare)} rail section(s) title themselves with a bare string "
        f"instead of sectTitle(): {bare[:3]}"
    )
    assert src.count("sectTitle(") >= 14


def test_the_reach_panel_has_headings_inside_it(page, web_server, reach_run):
    """The longest panel on the rail: four in-panel labels over ~14,000
    characters. They were styled <div>s, so there was no way to jump from the
    verdict to the potency table."""
    tags = page.evaluate(
        """async ([base, run]) => {
             const m = await import(base + '/js/ui/reach.js');
             const p = m.createReachPanel();
             p.render(run);
             return [...p.el.querySelectorAll('h1,h2,h3,h4,h5,h6')]
               .map(e => e.tagName + ':' + e.textContent.trim());
           }""",
        [web_server, reach_run],
    )
    assert len(tags) >= 4, f"reach panel exposes no heading structure: {tags}"
    assert all(t.startswith("H4:") for t in tags), (
        f"in-panel labels must be h4 under the h3 section title: {tags}"
    )
    assert any("potency evidence" in t for t in tags)


def test_the_ratio_is_read_from_the_payload_not_divided_in_the_browser(
    page, web_server
):
    """All three ratios are supplied by app/server/jobs.py so they carry
    provenance, can be asserted against and land in result.json. A number
    computed in a browser has none of those properties.

    Pinned by handing the panel a SUPPLIED ratio that disagrees with the two
    concentrations beside it. If the client is still dividing, it prints 0.083.
    """
    run = _run_fixture()
    run["extra_metrics"] = dict(run["extra_metrics"],
                                plasma_peak_over_ic50=0.222)
    text = _render(page, web_server, "finding.js", run)
    assert "= 0.222" in text, "the panel divided instead of reading the field"
    assert "Derived in this browser" not in text


def test_a_record_predating_the_field_says_its_ratio_was_derived(page, web_server):
    """~90 run records here were written before those fields existed. Refusing
    to disclose on them would hide the ratio exactly where a reader is
    re-reading an old result; dividing silently would present a browser
    arithmetic as a recorded one. It divides, and it says so."""
    run = _run_fixture()
    em = dict(run["extra_metrics"])
    for k in ("plasma_peak_over_ic50", "tissue_median_over_plasma_peak",
              "tissue_median_over_plasma_mean"):
        em.pop(k)
    run["extra_metrics"] = em
    text = _render(page, web_server, "finding.js", run)
    assert "0.083" in text, "an old record loses the disclosure entirely"
    assert "Derived in this browser" in text
    assert "result.json" in text


def test_the_bound_is_read_off_the_peak_where_the_two_readings_disagree(
    page, web_server
):
    """The case that separates the two ratios, which the flagship cannot.

    On the shipped 6 h default the peak reading (3.71) and the mean reading
    (5.96) are BOTH above 1, so a panel wired to the wrong one says the same
    thing and no assertion on that run can tell them apart. This payload
    straddles: the median tumour cell sits at 40 uM, under the 50 uM peak the
    blood actually reached (0.800) but over the 31.1 uM mean the plasma
    averaged across the window (1.286). It is not a contrived shape — it is
    what any window long enough for plasma to decay produces.

    C <= Cp is a statement about the HIGHEST concentration the blood ever
    reached, so this run does NOT break the bound. A panel reasoning off the
    mean would report that it does, and would then report the bound satisfied
    on a short window and violated on a long one with nothing but the dosing
    schedule between them.
    """
    run = _run_fixture()
    run["extra_metrics"] = dict(
        run["extra_metrics"],
        E_p50_uM=40.0,
        tissue_median_over_plasma_peak=40.0 / 50.0,           # 0.800
        tissue_median_over_plasma_mean=40.0 / 31.09499320754238,  # 1.286
    )
    text = _render(page, web_server, "finding.js", run)
    assert "0.80x the peak" in text, "the bound is not read against the peak"
    assert "That is the bound broken" not in text, (
        "the panel calls the bound broken on a run whose tissue never exceeds "
        "the plasma peak — it is reading the mean, which inflates with the "
        "window length rather than with the violation"
    )
    # It must NOT read as reassurance either: the law has no ceiling whichever
    # side of 1 one run happens to land on.
    assert "NOT evidence the law respects it" in text
    # The mean is still shown, still marked as not the bound.
    assert "1.29x the plasma MEAN" in text
    assert "Do not read that as the bound" in text
