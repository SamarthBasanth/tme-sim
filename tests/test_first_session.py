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
        "extra_metrics": {"E_p10_uM": 135.42888240722527,
                          "E_p90_uM": 270.1005830234793},
        "report": {"source": "xenium", "n_cells": 40887, "n_tumor": 23812},
        "provenance": {"caveats": [], "claim_set": None},
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
            "joint_limited": {"count": 19603,
                              "fraction_of_exposure_limited": 0.9272942289498581},
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
