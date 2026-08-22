// reach.js — DELIVERY OR DOSE: the discrimination this product now sells.
//
// One compound, one section, one question: when this compound failed on this
// tissue, was it because the drug never arrived, or because it arrived and was
// not potent enough at the level it arrived at? Every number here comes from
// the backend (`app/server/reach.py` and the per-cell mechanism attribution);
// nothing on this panel is computed in the client and nothing is a constant
// written here.
//
// WHY THIS AND NOT A RANKING
// --------------------------
// This panel replaced a ranked drug panel, and the argument for the trade is
// not that ranking was unpopular — it is that the new claim can be WRONG in a
// way the old one could not:
//
//   "This compound reached the tissue and was still 3.3-fold short of IC50" is
//   checkable against a bench measurement on the same section. "Cisplatin
//   outranks temozolomide" was only ever checkable against an outcome nobody
//   has.
//
// The product got smaller and acquired a falsification criterion. Every number
// on this panel is therefore written to be checkable rather than persuasive:
// the shortfall is quoted as a fold-change against a measurable concentration,
// the delivery spread is a ratio of two percentiles of one field, and the
// verdict is read off per-cell labels that reconcile bit-exactly with
// kill_fraction rather than being asserted alongside it.
//
// WHY THIS PANEL EXISTS, AND WHY ITS HEADLINE IS A RATIO
// -----------------------------------------------------
// The obvious headline — "what fraction of the tumour reached 1x IC50" — is
// `therapeutic_coverage_1x`, and it is NOT usable as the headline. It compares
// an exposure to an IC50, so it reads BOTH the plasma level and the potency,
// and those are the two per-compound quantities measured to be uncertain enough
// to break the cross-compound ordering this product used to publish
// (results/pk_ranking_robustness.md). Its own implementation says so:
// tmesim/validation/delivery_metrics.py calls it "absolute -- scaling the field
// scales the answer". Shipping it as the headline would rebuild the withdrawn
// failure one compound at a time.
//
// What survives is a RATIO TAKEN INSIDE ONE MAP: this region got N times less
// than that region. A per-compound plasma factor is an exact linear rescale of
// the whole field, so it moves every point together and cancels from a ratio
// between points. The headline is therefore the fraction of viable tumour cells
// whose exposure is below a fixed fraction of that same map's own perivascular
// reference — no IC50, no dose, no other compound.
//
// Coverage is NOT deleted. It is the number a pharmacologist will ask for, and
// hiding it would be its own dishonesty. It is DEMOTED, and it is shown across
// the dose ladder rather than as a single value, so the reader can watch it
// move instead of being handed one figure that looks settled.
//
// THE CURVE IS THE POINT, NOT THE NUMBER ON IT. F(x) is drawn across the whole
// ladder precisely so that the 0.1 threshold reads as the reporting convention
// it is. A single number at one threshold acquires false authority — that is
// exactly how `1x IC50` came to be defended as though it meant something.

import { h, clear } from '../util/dom.js';
import { pct, num, sci, NULL_GLYPH } from '../util/fmt.js';
import { doseInputBlock, doseFraming } from './doseinput.js';
import { jointLimitedCaveat } from './verdict.js';

/* Chart geometry. The rail is ~347 px wide, so these are drawn into a viewBox
 * and scaled by CSS rather than sized in pixels. */
const VB_W = 300;
const VB_H = 84;
const PAD_L = 30;
const PAD_R = 8;
const PAD_T = 8;
const PAD_B = 18;

const PLOT_W = VB_W - PAD_L - PAD_R;
const PLOT_H = VB_H - PAD_T - PAD_B;

function svgEl(name, attrs, ...kids) {
  const el = document.createElementNS('http://www.w3.org/2000/svg', name);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v !== null && v !== undefined) el.setAttribute(k, String(v));
  }
  for (const k of kids) if (k) el.appendChild(k);
  return el;
}

function svgTitle(text) {
  const t = document.createElementNS('http://www.w3.org/2000/svg', 'title');
  t.textContent = text;
  return t;
}

/**
 * One line chart, single series.
 *
 * Single series, so there is no legend — the section title names it (the
 * dataviz rule: a legend box for one series is noise). The mark is a 2 px line
 * with 8 px-diameter markers; grid and axis are recessive `--line-1`; every
 * label wears a TEXT token rather than the series colour, so identity is
 * carried by the mark and never by coloured type. Each point carries a native
 * `<title>`, which is the hover layer, and the same numbers are printed as a
 * table underneath — that table is the accessible view, and it is also why no
 * value is ever invented to fill a gap.
 */
function lineChart({ xs, ys, xLabel, yLabel, fmtX, fmtY, markIndex, yMax }) {
  const n = xs.length;
  const top = Number.isFinite(yMax) ? yMax : Math.max(1, ...ys.filter(Number.isFinite));
  const px = (i) => PAD_L + (n === 1 ? PLOT_W / 2 : (PLOT_W * i) / (n - 1));
  const py = (v) => PAD_T + PLOT_H - (PLOT_H * (Number.isFinite(v) ? v / top : 0));

  const kids = [];
  // Recessive frame: a baseline and a top rule, nothing more. No gridlines at
  // this size — at 84 px tall they would out-weigh the data.
  kids.push(svgEl('line', {
    x1: PAD_L, y1: PAD_T + PLOT_H, x2: PAD_L + PLOT_W, y2: PAD_T + PLOT_H,
    stroke: 'var(--line-1)', 'stroke-width': 1,
  }));
  kids.push(svgEl('line', {
    x1: PAD_L, y1: PAD_T, x2: PAD_L, y2: PAD_T + PLOT_H,
    stroke: 'var(--line-1)', 'stroke-width': 1,
  }));

  const pts = [];
  for (let i = 0; i < n; i += 1) if (Number.isFinite(ys[i])) pts.push(`${px(i)},${py(ys[i])}`);
  if (pts.length > 1) {
    kids.push(svgEl('polyline', {
      points: pts.join(' '), fill: 'none', stroke: 'var(--accent)',
      'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round',
    }));
  }
  for (let i = 0; i < n; i += 1) {
    if (!Number.isFinite(ys[i])) continue;
    const isMark = i === markIndex;
    const dot = svgEl('circle', {
      cx: px(i), cy: py(ys[i]), r: isMark ? 4 : 2.6,
      fill: isMark ? 'var(--accent)' : 'var(--bg-1)',
      stroke: 'var(--accent)', 'stroke-width': 2,
    }, svgTitle(`${xLabel} ${fmtX(xs[i])} — ${yLabel} ${fmtY(ys[i])}`));
    kids.push(dot);
  }

  // Axis extremes only. Selective labels: the two ends of x, the top of y.
  const tick = (x, y, text, anchor) => svgEl('text', {
    x, y, 'text-anchor': anchor || 'middle',
    fill: 'var(--fg-2)', 'font-size': 8, 'font-family': 'var(--font-mono)',
  }, null);
  const t1 = tick(PAD_L, VB_H - 5, '', 'start'); t1.textContent = fmtX(xs[0]);
  const t2 = tick(PAD_L + PLOT_W, VB_H - 5, '', 'end'); t2.textContent = fmtX(xs[n - 1]);
  const t3 = tick(PAD_L - 4, PAD_T + 4, '', 'end'); t3.textContent = fmtY(top);
  const t4 = tick(PAD_L - 4, PAD_T + PLOT_H, '', 'end'); t4.textContent = fmtY(0);
  kids.push(t1, t2, t3, t4);

  // Size in CSS, not in SVG attributes: `height="auto"` is not a valid SVG
  // length and the browser logs it as an error on every render. The viewBox
  // carries the aspect ratio, so width:100% + height:auto in CSS scales it
  // correctly at any rail width without an attribute the parser rejects.
  return svgEl('svg', {
    viewBox: `0 0 ${VB_W} ${VB_H}`,
    role: 'img', preserveAspectRatio: 'xMidYMid meet',
    style: 'display:block;width:100%;height:auto',
    'aria-label': `${yLabel} against ${xLabel}`,
  }, ...kids);
}

/** A prose paragraph at the disclosure floor. Every sentence on this panel goes
 *  through here so none of them can drift below it.
 *
 *  Variadic, so a paragraph can carry glossed terms (see `term()`) inline
 *  rather than being forced to a single flat string. A lone string argument
 *  still behaves exactly as before. */
function prose(...kids) {
  return h('p', {
    style: {
      margin: '0 0 var(--sp-4)',
      fontSize: 'var(--fs-disclosure)',
      lineHeight: 'var(--lh-body)',
      color: 'var(--fg-1)',
    },
  }, ...kids);
}

/* An in-panel label, and it is an `<h4>` rather than a styled `<div>`.
 *
 * This panel is the longest thing on the rail and it has four of these —
 * delivery evidence, the two charts, potency evidence. They were divs, so the
 * only heading structure inside ~14,000 characters was the section title above
 * them, and a screen-reader user had no way to jump from the verdict to the
 * potency table. h4 because the section title is h3 (js/ui/rail.js::sectTitle).
 * Every visual property is still set here, so nothing moves. */
function microLabel(text) {
  return h('h4', {
    style: {
      fontSize: 'var(--fs-micro)', color: 'var(--fg-2)',
      fontFamily: 'var(--font-mono)', letterSpacing: 'var(--tr-wide)',
      textTransform: 'uppercase', margin: 'var(--sp-4) 0 var(--sp-2)',
      fontWeight: '400',
    },
  }, text);
}

/* ------------------------------------------------- INVENTED VOCABULARY ---
 * Words this project made up, and what they mean to someone who has never
 * read this repo.
 *
 * A first-time reader was measured against this panel and stalled on all six.
 * They are not obscure by accident — each one names a quantity that had no
 * name before, which is exactly the kind of term a reader cannot look up and
 * cannot guess. "Under-reached fraction" in particular reads as though it must
 * mean "the fraction the drug failed to reach", and it does not: it is a count
 * against a threshold that is itself a chosen convention.
 *
 * Rendered as a `<abbr>` with a dotted underline: the term keeps its exact
 * spelling (so it still matches the results files and the API payload) and
 * acquires a hover gloss and a screen-reader expansion. NOTHING is replaced —
 * a reader who already knows the vocabulary sees the same words in the same
 * places. This adds a reading, it does not substitute one. */
const GLOSSARY = {
  'under-reached fraction':
    'The share of viable tumour cells whose exposure sits below a chosen '
    + 'fraction of the best-perfused tissue in this same map. It is a count '
    + 'against a threshold, NOT a measurement of where the drug failed to go — '
    + 'the threshold is a reporting convention this project picked, and the '
    + 'whole ladder of thresholds is printed so the choice stays visible.',
  'ladder':
    'The full table of thresholds, rather than one. Printed because reading a '
    + 'single threshold as "the answer" is the mistake that made "1x IC50" look '
    + 'settled when it was not.',
  'rung':
    'One row of that table — one threshold. Different rungs have very '
    + 'different robustness, so "the ladder moves X" and "the rung you are '
    + 'using moves X" are different statements.',
  'delivery spread':
    'How uneven delivery is inside this one section: the exposure of the '
    + 'best-perfused tenth of tumour cells divided by that of the worst-'
    + 'perfused tenth. Both numbers come from the same map, so the dose, the '
    + 'potency and the vessel-wall permeability all cancel.',
  // Two entries, because the reference is NOT the same physical quantity in
  // the two modes: the finite-window map normalises by an exposure (uM*s), the
  // steady map by a settled concentration (uM). The payload field is called
  // `E_ref_uM_s` in both, which is why the word on screen is chosen from
  // `E_ref_units` and never from the field name.
  'reference exposure':
    'The exposure in the tissue immediately next to a vessel, in this same '
    + 'run. Every ratio on this panel is some other exposure divided by this '
    + 'one, which is what makes those ratios carry no dose and no potency.',
  'reference concentration':
    'The settled drug concentration in the tissue immediately next to a '
    + 'vessel, in this same run. Every ratio on this panel is some other '
    + 'concentration divided by this one, which is what makes those ratios '
    + 'carry no dose and no potency.',
  'exposure':
    'The time-integral of drug concentration at a point — concentration '
    + 'multiplied by time, in uM.s. It is what the transport solve produces; '
    + 'it is not a concentration and not a dose.',
};

/** The term, verbatim, carrying its plain-language gloss. */
function term(word) {
  const gloss = GLOSSARY[String(word).toLowerCase()];
  if (!gloss) return word;
  return h('abbr.gloss', { title: gloss }, word);
}

/* The axis markers the backend already writes into `scale_invariance`:
 * "(1) PLASMA/DOSE SCALE:", "(2) TISSUE CAPACITY:", "(3) TRANSVASCULAR
 * PENETRANCE:". Deliberately narrow — a digit in parens, then SHOUTED words,
 * then a colon — so it cannot fire on ordinary prose.
 *
 * The em dash is in the class because the two modes do not word the marker
 * identically: the finite-window note writes "(2) TISSUE CAPACITY:" and the
 * steady-state note writes "(2) TISSUE CAPACITY — A STATEMENT ABOUT THE MAP:".
 * Without it the steady note splits into two parts instead of three and the
 * capacity argument — the whole reason the product moved to a steady map —
 * stays buried inside the plasma paragraph. Nothing is lost either way; the
 * lossless contract on splitAxes() holds for any input. */
const AXIS_MARK = /\(\d+\)\s+[A-Z][A-Z0-9 /,'–—-]*[A-Z]\s*:/g;

/**
 * Split one long caveat into `{head, body}` paragraphs on its own axis
 * markers, LOSSLESSLY.
 *
 * The contract, and the only thing that makes this safe to do to a caveat:
 *
 *     splitAxes(s).map(p => (p.head || '') + p.body).join('') === s
 *
 * for every input, including one with no markers at all (which returns a
 * single part carrying the whole string). Nothing is dropped, nothing is
 * reordered, nothing is reworded — the marker is promoted from mid-sentence
 * to its own line and that is the entire transformation.
 *
 * @param {string} s
 * @returns {{head: string|null, body: string}[]}
 */
export function splitAxes(s) {
  const text = String(s);
  const marks = [];
  AXIS_MARK.lastIndex = 0;
  let m;
  while ((m = AXIS_MARK.exec(text)) !== null) marks.push({ at: m.index, mark: m[0] });
  if (!marks.length) return [{ head: null, body: text }];

  const out = [];
  if (marks[0].at > 0) out.push({ head: null, body: text.slice(0, marks[0].at) });
  for (let i = 0; i < marks.length; i += 1) {
    const start = marks[i].at + marks[i].mark.length;
    const end = i + 1 < marks.length ? marks[i + 1].at : text.length;
    out.push({ head: marks[i].mark, body: text.slice(start, end) });
  }
  return out;
}

// Leading zero KEPT deliberately. Stripped, `0.1x` renders as `.1x`, which in a
// mono column sits one glyph away from `1x` — and the difference between them
// is a factor of ten in the quantity this whole section is about.
const fmtRatio = (v) => String(v);
const fmtFrac = (v) => (Number.isFinite(v) ? pct(v, 1) + '%' : NULL_GLYPH);

/**
 * THE VERDICT — is this a delivery problem or a dose problem — and the three
 * things that qualify it: the exposure shortfall, how much of the answer was
 * an input, and the joint-limited counterfactual.
 *
 * EXTRACTED FROM render() 2026-08-22, verbatim. It depends on `run` alone,
 * which is now readable off the signature rather than by tracing ~90 lines
 * upward through a closure. See renderFalloff below for the character-exact
 * losslessness check that licences moving a block of caveats at all.
 */
function renderVerdict(el, run) {
  /* ---- THE VERDICT: is this a delivery problem or a dose problem? ----- *
   * This is the answer the panel exists to give, and it is the first thing
   * on it. It is read straight off the per-cell mechanism attribution, which
   * is per-compound by construction and reconciles bit-exactly with
   * kill_fraction, so nothing here is a second opinion about the same
   * numbers — it is a restatement of the engine's own labels.
   *
   * Why the verdict and not the reach number leads: measured on the flagship
   * Xenium GBM section at dx = 12 um, the under-reached fraction at 0.1x is
   * 0.05% — the drug gets essentially everywhere — while 99.8% of surviving
   * tumour cells are DOSE-limited. A panel that led with "0.0% under-reached"
   * would have buried its own finding. */
  const mech = run.mechanism || null;
  const bd = (mech && mech.survival_breakdown) || null;
  if (bd) {
    const entries = Object.entries(bd).filter(([, v]) => Number.isFinite(v) && v > 0);
    entries.sort((a, b) => b[1] - a[1]);
    const [topKey, topVal] = entries[0] || [null, null];
    const LABEL = {
      delivery: 'DELIVERY-LIMITED', dose: 'DOSE-LIMITED',
      resistance: 'RESISTANCE-LIMITED', chance: 'CHANCE',
      unattributed: 'UNATTRIBUTED',
    };
    const MEANING = {
      delivery: 'the drug did not arrive in a pharmacologically meaningful '
        + 'amount. This is the case the exposure map below is evidence FOR.',
      dose: 'the drug arrived and stayed below IC50. The exposure map below '
        + 'is what rules delivery OUT — it is the reason this reads as a '
        + 'potency finding rather than an unresolved one.',
      resistance: 'exposure was adequate and the cell state suppressed the '
        + 'kill anyway.',
      chance: 'exposure and biology were both adequate; the cell drew above '
        + 'its kill probability.',
      unattributed: 'not explained by the model’s own bands. A non-zero share '
        + 'here is itself the finding.',
    };
    if (topKey) {
      el.appendChild(h('div', {
        style: { display: 'flex', alignItems: 'baseline', gap: 'var(--sp-4)',
                 margin: '0 0 var(--sp-3)' },
      },
      h('span', { style: { fontSize: 'var(--fs-xl)', fontFamily: 'var(--font-mono)',
                           color: 'var(--fg-0)', fontWeight: '600' } },
      fmtFrac(topVal)),
      h('span', { style: { fontSize: 'var(--fs-disclosure)', color: 'var(--fg-1)',
                           lineHeight: 'var(--lh-body)' } },
      h('b', LABEL[topKey] || topKey.toUpperCase()), ' — ',
      `${fmtFrac(topVal)} of SURVIVING tumour cells. `,
      MEANING[topKey] || '')));
    }
    const fold = mech.exposure_shortfall_fold;
    if (Number.isFinite(fold)) {
      el.appendChild(prose(
        `The median tumour cell would need ${num(fold, 2)}x its exposure to `
        + 'reach IC50. Above 1 the median cell is short of drug; below 1 it is '
        + 'already past the half-effect point and something else is holding '
        + 'the kill down. This is an ABSOLUTE quantity — it reads the potency '
        + '— and is stated here because the verdict it supports would be '
        + 'unreadable without it, not because it is protected by anything in '
        + 'the invariances below.'));
    }
    /* ---- HOW MUCH OF THIS VERDICT WAS AN INPUT? ---------------------
     * Compact form, directly beside the verdict, because that is the number
     * it qualifies. The full separation — what the map added versus what the
     * spec fixed — is under `The finding`, which is the section immediately
     * above this one on the same rail and is never collapsed.
     *
     * The panel above this one already refuses the structurally identical
     * case for bbb_penetrance = 0 (see the null_reason branch). This is that
     * refusal extended to the dose axis, and it too picks no threshold: it
     * prints plasma_peak / IC50 and says what follows by arithmetic below 1.
     * js/ui/doseinput.js carries the argument. */
    el.append(...doseInputBlock(run, {
      compact: true,
      pointTo: 'What the map added on top of that ratio — the spread across '
        + 'this tissue — is separated out under The finding, above.',
    }));

    // AUTHORED IN js/ui/verdict.js. This panel used to paraphrase the same
    // caveat independently and WITHOUT the cell count that finding.js
    // printed, so the identical limit was weaker here than three sections
    // up. The lead is still this panel's own wording; the body is not.
    const jlCaveat = jointLimitedCaveat(run, {
      lead: 'DOSE IS NOT A CLEAN COUNTERFACTUAL HERE. ',
    });
    if (jlCaveat) el.appendChild(jlCaveat);
  }
}

/**
 * Delivery heterogeneity — P90/P10 of tumour-cell exposure — and the caveat
 * that it cannot see a permeability barrier.
 *
 * EXTRACTED FROM render() 2026-08-22, verbatim; depends on `run` alone.
 */
function renderSpread(el, run) {
  /* ---- delivery heterogeneity, the headline of this half ------------- *
   * P90/P10 over tumour-cell exposure. It is a ratio of two percentiles of
   * ONE field, so it is scale-free for exactly the same reason the reach
   * ratio is, and it carries signal even where the reach headline does not.
   * That is not hypothetical: measured on the flagship Xenium GBM section at
   * dx = 12 um, the under-reached fraction at 0.1x is 0.05% — the drug gets
   * essentially everywhere — while this ratio is a legible 2.0x.
   *
   * It also needs NO CONVENTION. There was a single-threshold headline here
   * (F at 0.1x) and it was removed rather than retuned: any threshold chosen
   * now would be chosen knowing what the field does, which is post-hoc
   * selection wearing a convention's clothes. The full F(x) ladder is still
   * published below as the evidence; what is gone is the practice of reading
   * one rung of it as the answer.
   *
   * And its denominator is GUARANTEED grid-independent rather than merely
   * observed to be. `exposure_percentiles` is computed over the TUMOUR cell
   * set, which tmesim/validation/delivery_metrics.py::tumor_cell_exposure
   * documents as "a FIXED set of measured Xenium cells, identical at every dx
   * — which is what makes these statistics grid-comparable". That is a
   * structural property: the cells are measured, so no grid can move them.
   *
   * CORRECTION 2026-08-22. This passage previously argued the stronger case
   * that the viable band carries "~10x grid sensitivity", citing "89% of
   * tumour cells at dx = 30 um" against 2,008 at dx = 12. The dx = 12 figure
   * is right and the dx = 30 figure was wrong by a factor of ~10.5. Measured
   * through reconstruct_real on this section, temozolomide, 6 h, convection
   * ON: 2,008 of 23,812 at dx = 12 um (8.43%), 2,007 at dx = 15, 2,019 at
   * dx = 30 (8.48%). Across a 2.5x change in grid the viable band moves by
   * 12 cells out of ~2,010 — 0.6%, not 10x.
   *
   * So the contrast this paragraph used to draw does not exist, and the
   * reason to prefer this ratio is the weaker but true one above: the viable
   * band is a COMPUTED quantity that happens to be stable on this section,
   * carrying no guarantee that it is stable on another, whereas a measured
   * cell set cannot move by construction. Stability observed on one section
   * is not stability. */
  const ex0 = run.extra_metrics || {};
  const p10 = ex0.E_p10_uM;
  const p90 = ex0.E_p90_uM;
  if (Number.isFinite(p10) && Number.isFinite(p90) && p10 > 0) {
    el.appendChild(h('div', {
      style: { display: 'flex', alignItems: 'baseline', gap: 'var(--sp-4)',
               margin: '0 0 var(--sp-3)' },
    },
    h('span', { style: { fontSize: 'var(--fs-lg)', fontFamily: 'var(--font-mono)',
                         color: 'var(--fg-0)' } }, num(p90 / p10, 2) + 'x'),
    h('span', { style: { fontSize: 'var(--fs-disclosure)', color: 'var(--fg-1)',
                         lineHeight: 'var(--lh-body)' } },
    h('b', term('DELIVERY SPREAD')), ' — the best-perfused tenth of tumour cells '
      + 'receives this many times the exposure of the worst-perfused tenth. A '
      + 'ratio of two percentiles of one field: no dose, no potency, no '
      + 'threshold to choose, and computed over the tumour cell set rather '
      + 'than the viable band, so it does not move with the grid the way a '
      + '"fraction of viable tumour" does.')));
    // THIS MAP CANNOT SEE A PERMEABILITY BARRIER, and it is the first thing a
    // reader will assume it does. bbb_penetrance is a pure global multiplier
    // on the source term, so it divides out of every within-map ratio —
    // verified bit-identical at bbb = 0.1 and bbb = 1.0. Scale-freedom is why
    // this figure works, and the barrier is a scale; the two cannot both be
    // had. It also explains the bevacizumab case structurally rather than as
    // an edge case: bbb = 0 is the one point where the quantity this metric
    // is blind to becomes the entire answer.
    el.appendChild(h('div.caveat.crit', h('span',
      h('b', 'THIS IS NOT A PERMEABILITY FINDING. '),
      'The reach ratio is EXACTLY invariant to transvascular penetrance — '
      + 'verified bit-identical at bbb_penetrance 0.1 and 1.0 — because '
      + 'penetrance multiplies every point of the map equally and therefore '
      + 'divides out of a ratio between points. That is the same property '
      + 'that makes this number immune to the dose, and it cannot be had '
      + 'selectively. So the delivery evidence here depends on the compound '
      + 'ONLY through its penetration length, L = sqrt(D_eff/k). It says how '
      + 'STEEPLY exposure falls away from a vessel; it says nothing about how '
      + 'much crosses the vessel wall in the first place. Do not read it as a '
      + 'brain-penetration or barrier result.')));
  }
}

/**
 * The exposure-versus-vessel-distance figure, and the three limits on it.
 *
 * EXTRACTED FROM render() 2026-08-22, verbatim. render() was ~660 lines in one
 * closure and every block in it read as though it might depend on the ones
 * above; this one depends on `reach` and nothing else, which is now visible in
 * its signature instead of having to be established by reading upward.
 *
 * The move is text-identical and was checked that way: the panel's rendered
 * innerText is compared character-for-character across eight payload shapes
 * before and after. That check is the only licence for touching a function
 * whose output is mostly caveats — a refactor that drops a sentence is a
 * disclosure silently shortened, and this repo may not do that.
 */
function renderFalloff(el, reach) {
  /* ---- falloff ------------------------------------------------------- */
  const fall = reach.falloff || {};
  const fd = fall.distance_um || [];
  const fv = fall.normalised_mean_auc || [];
  if (fd.length && fv.length === fd.length) {
    el.appendChild(microLabel('exposure vs vessel distance'));
    el.appendChild(lineChart({
      xs: fd, ys: fv,
      xLabel: 'distance', yLabel: 'exposure / reference',
      fmtX: (v) => num(v, 0) + 'um', fmtY: (v) => num(v, 2),
      markIndex: 0, yMax: Math.max(1, ...fv.filter(Number.isFinite)),
    }));
    el.appendChild(prose(
      'Each bin is the mean exposure over tissue at that distance from the '
      + 'nearest endothelial cell, divided by the perivascular reference, so '
      + 'the first bin is 1.00 by construction. This curve is the SHAPE of '
      + 'the map, and the shape is the one thing a within-map ratio does not '
      + 'protect: it is set by the effective diffusivity, whose descriptor '
      + 'inputs have never been scored '
      + '(results/descriptor_falsification.md). Read it as the best available '
      + 'estimate of a gradient, not as a measured profile.'));
    // THE SCOPE LIMIT ON THIS HALF OF THE PANEL, and it is not symmetric.
    // The obvious objection to any delivery map is "isn't this just showing
    // me where my vessels aren't?" — the field inherits the shape of the
    // vascular source wherever the penetration length is SHORT compared with
    // vessel spacing. Measured against the KDE source null on this section,
    // the excess over the null is ~7.9x for a long-penetration compound and
    // only 1.34x for doxorubicin, a short-L one. So this curve is strong
    // evidence for compounds that penetrate far and WEAK evidence for
    // compounds that do not — and the weak case is the one where a customer
    // is most likely to suspect they are looking at their own angiogram.
    // Stated here rather than only in a results file, because a scope limit
    // that lives away from the figure it limits does not get read.
    // THE CAPACITY RESIDUAL, and it is the largest limit on this figure —
    // larger than the descriptor gap and larger than the source-null one.
    // Measured on the flagship section: the normalised falloff moves 0.3772
    // ABSOLUTE at capacity R = 39 (the best-available constant), on a curve
    // whose entire range is 1.0 to 0.51, with the fitted penetration length
    // collapsing 279.7 -> 76.7 um. Those three figures reproduced to four
    // significant figures on re-measurement (0.37723, 1.0 -> 0.51284,
    // 279.715 -> 76.740) and are the solid part of this block.
    // What did NOT survive is the asymmetry that used to be drawn from them
    // -- the ladder is not the safe counterpart, it moves 0.4021 over the
    // same range. The real asymmetry is BETWEEN RUNGS, not between the
    // ladder and the curve. See the dated correction below.
    /* CORRECTED 2026-08-22. This caveat previously read "THE SHAPE IS NOT
     * CAPACITY-INVARIANT — the ladder above nearly is, and this is not",
     * and printed the two figures side by side on that basis. The asymmetry
     * does not exist: measured on the flagship at R = 39 the falloff moves
     * 0.3772 and the whole 7-rung ladder moves 0.4021. They move by the SAME
     * amount. What is nearly immune is the bottom four rungs (x <= 0.1,
     * which include the headline); the top rung F(1.0) goes 0.4056 -> 0.8077
     * and is plotted on the same chart the caveat was reassuring about.
     * It also asserted the 6 h window "has not converged", citing 0.01302
     * decaying as 281/T. That was measuring the ESTIMATOR: any AUC-based
     * ratio carries an O(1/T) tail however converged the field is, and the
     * field is converged to 1.112e-05 at 6 h.
     * Sources: results/reach_invariances_measure.md, steady_state_reach.md. */
    el.appendChild(h('div.caveat.crit', h('span',
      reach.window_mode === 'steady_state'
        ? [h('b', 'CAPACITY CANNOT REACH THIS FIGURE. '),
           'A tissue-capacity term R multiplies the TIME DERIVATIVE and '
           + 'nothing else, and a converged steady state has no time '
           + 'derivative — so R is absent from this map by inspection, not '
           + 'small. Measured across R = 1 to 4184: this curve moves at most '
           + '4.5e-14 and the ladder above moves exactly 0.000e+00, which is '
           + 'float64 round-off and an unmoved threshold count, not a '
           + 'residual. This is the one limit that a finite-window reading '
           + 'of the same field DOES carry — there, at R = 39, this curve '
           + 'moves 0.3772 and the ladder 0.4021.']
        : [h('b', 'THIS SHAPE AND THE LADDER ABOVE ARE BOTH '
             + 'CAPACITY-SENSITIVE, BY ABOUT THE SAME AMOUNT. '),
           'The engine assigns every compound a tissue capacity of 1.0; '
           + 'measured values span 3.2 to 4184 mL/g (n = 33, best available '
           + 'constant 39). Capacity rescales TIME, so over a finite window '
           + 'it changes the shape the field has reached. At R = 39 this '
           + 'curve moves 0.3772 absolute on a curve whose entire range is '
           + '1.0 to 0.51, the fitted penetration length collapses from '
           + '279.7 to 76.7 um, and the ladder above moves 0.4021 — the same '
           + 'size. Do not read one as solid because the other is flagged. '
           + 'The exception is the BOTTOM of the ladder: F(0.1) and below, '
           + 'which include the headline, are bit-identical out to R = 80. '
           + 'The movement is almost entirely in the top rung F(1.0), '
           + '0.4056 to 0.8077, which is on this chart.'])));
    el.appendChild(h('div.caveat', h('span',
      h('b', 'HOW MUCH OF THIS IS THE VASCULATURE? '),
      'Where a compound’s penetration length is short compared with the '
      + 'spacing between vessels, the exposure field largely inherits the '
      + 'SHAPE of the vascular source, and this curve then says more about '
      + 'the tissue’s vessels than about the compound. Scored against the '
      + 'source null on this section, the excess over vasculature-alone runs '
      + 'about 7.9x for a long-penetration compound but only 1.34x for '
      + 'doxorubicin, which penetrates poorly. Read the falloff as compound-'
      + 'specific evidence in proportion to the penetration length reported '
      + 'in Metrics, and treat a short one as close to a vessel map.')));
  }
}

/**
 * The demoted absolute reading: therapeutic coverage, as a ladder.
 *
 * EXTRACTED FROM render() 2026-08-22, verbatim — see renderFalloff above for
 * the losslessness check that licences the move. Depends on `run` alone.
 */
function renderPotency(el, run) {
  /* ---- the demoted absolute reading ---------------------------------- */
  const ex = run.extra_metrics || {};
  const cov = [
    ['0.5x IC50', ex.therapeutic_coverage_0_5x],
    ['1x IC50', run.metrics ? run.metrics.therapeutic_coverage_1x : null],
    ['2x IC50', ex.therapeutic_coverage_2x],
  ].filter(([, v]) => v !== undefined);
  if (cov.length) {
    el.appendChild(microLabel('potency evidence'));
    el.appendChild(h('div.tblwrap', h('table.tbl.zebra',
      h('thead', h('tr', h('th', 'threshold'), h('th.right', 'of tumour above it'))),
      h('tbody', cov.map(([k, v]) => h('tr',
        h('td.mono', k),
        h('td.right.num.mono', fmtFrac(v))))))));
    el.appendChild(h('div.caveat', h('span',
      h('b', 'ABSOLUTE, and shown as a ladder for that reason. '),
      'Coverage compares an exposure to an IC50, so unlike everything above '
      + 'it, it reads both the plasma level and the potency — the two '
      + 'per-compound quantities whose spread withdrew the cross-compound '
      + 'ordering (results/pk_ranking_robustness.md). Scaling the exposure '
      + 'field scales this answer, which is exactly what the ratios above are '
      + 'immune to. It is shown across three thresholds rather than at one so '
      + 'that its movement is visible; a single coverage number looks settled '
      + 'and is not.')));
    /* AND WHERE THE LADDER STARTS FROM. This table is the one place on the
     * rail where `therapeutic_coverage_1x` is printed, and it is the number a
     * reader is most likely to quote — "0.0% of the tumour reached IC50" is a
     * complete sentence and reads as a result. On the shipped default it is
     * also the number most constrained by the spec: the entered plasma peak
     * is under IC50, so under a plasma-bounded source term every row of this
     * table at 1x and above is 0.0% for every tissue.
     *
     * One sentence, not the block: the ratio is stated in full beside the
     * verdict at the top of this panel and again under The finding. What this
     * adds is the pointer at the exact figure it governs, which is three
     * screens below both of them. */
    const framing = doseFraming(run);
    if (framing && framing.ratio < 1) {
      el.appendChild(h('div.caveat.crit', h('span',
        h('b', 'THIS LADDER STARTS BELOW ITS OWN TOP RUNG BY INPUT. '),
        `plasma_peak / IC50 for this run is ${num(framing.ratio, 3)} `
        + `(${num(framing.peak_uM, 3)} uM entered against an IC50 of `
        + `${num(framing.ic50_uM, 3)} uM), so the entered peak is `
        + `${num(framing.foldUnder, 3)}x under the half-effect concentration `
        + 'before any transport happens. Do not read the 1x and 2x rows as '
        + 'measurements of this tissue without that ratio attached — stated '
        + 'in full beside the verdict at the top of this panel.')));
    }
  }
}

export function createReachPanel() {
  const el = h('div');

  function renderEmpty(msg) {
    clear(el);
    const p = prose(msg);
    p.classList.add('dim');
    el.appendChild(p);
  }

  /**
   * @param {object|null} run a RunResult, or null.
   */
  function render(run) {
    if (!run) {
      renderEmpty('No run selected. The reach map appears once a run reaches '
        + 'status=done.');
      return;
    }
    const reach = run.reach || null;
    if (!reach) {
      renderEmpty(`Run ${run.run_id} carried no reach block. That is a backend `
        + 'omission rather than a result — no number is substituted for it.');
      return;
    }
    clear(el);

    /* ---- the refusal path. A refusal is a result, and it is stated in full,
     *      never collapsed into a dash.
     *
     * WHAT THE MODEL WAS TOLD vs WHAT THE MODEL FOUND. This distinction is the
     * whole content of the branch below, and it is the panel's worst failure
     * mode if it is missed. A compound parameterised with bbb_penetrance = 0
     * has an identically zero source term, so its exposure field is zero
     * everywhere, the perivascular reference is zero, and every reach ratio
     * would be a division by it — the map correctly refuses. But the mechanism
     * attribution simultaneously reports 100% DELIVERY-limited, and a reader
     * meeting a bare divide-by-zero message next to a confident 100% would
     * reasonably conclude the model DISCOVERED that the drug cannot reach the
     * tumour. It did not. It was told, in the compound's own parameters, and
     * then correctly propagated that input to its spatial consequence. That is
     * a useful thing for the panel to do and a dishonest thing for it to
     * imply, so it says which it is.
     *
     * Measured example (bevacizumab, Xenium GBM, dx = 12 um, 6 h): every one of
     * 23,812 tumour cells is SURVIVED_UNREACHED, mean tumour AUC is exactly
     * 0.0, and bbb_penetrance is 0.0 by parameterisation. --------------- */
    if (reach.null_reason) {
      const dp = run.drug_params || {};
      const bbb = dp.bbb_penetrance;
      const toldNotFound = Number.isFinite(bbb) && bbb === 0;
      el.appendChild(h('div.caveat.crit', h('span',
        h('b', 'NO REACH MAP. '), String(reach.null_reason))));
      if (toldNotFound) {
        el.appendChild(h('div.caveat.crit', h('span',
          h('b', 'THIS IS AN INPUT, NOT A FINDING. '),
          `This compound carries bbb_penetrance = 0, so the transport solve was `
          + 'given a source term of exactly zero and no drug enters the tissue '
          + 'anywhere. A verdict of DELIVERY-limited here is the model '
          + 'correctly propagating something it was TOLD to its spatial '
          + 'consequence — it is not something the model FOUND. There is no '
          + 'reach map because there is nothing to take a ratio of, and the '
          + 'delivery evidence on this panel can neither support nor challenge '
          + 'the verdict. To get a discrimination the model actually made, run '
          + 'a compound with a non-zero penetrance.')));
      }
      if (reach.convention) el.appendChild(prose(String(reach.convention)));
      return;
    }

    const x = reach.headline_x;

    renderVerdict(el, run);

    /* ---- DELIVERY EVIDENCE — the half that rules delivery in or out ----- */
    el.appendChild(microLabel('delivery evidence'));

    /* ---- WHICH QUESTION THIS MAP ANSWERS -------------------------------
     * MOVED TO THE TOP OF THIS HALF, 2026-08-22, and reworded from a label
     * into a question.
     *
     * The payload has always known whether the map came from a converged
     * steady state or from a finite exposure window; when it was first put on
     * screen it went in at the BOTTOM of the panel, under every caveat it
     * governs. Measured on a first session: a reader reaches the falloff
     * curve, the reach ladder and four capacity caveats — all of which mean
     * different things in the two modes — roughly 600 words before being told
     * which mode they are in. By then the reading is already formed.
     *
     * It leads this half instead. It does NOT lead the whole panel, and that
     * placement is deliberate: the verdict above (delivery vs dose) is read
     * off the per-cell mechanism attribution, which is not this object and
     * does not change with it. A banner over the verdict saying "FINITE
     * WINDOW" would attach a qualifier to a number it does not qualify.
     *
     * Phrased as the question rather than the object name. "steady_state" is
     * an implementation word; "what does this tissue's geometry permit" is
     * the thing a bench scientist is actually being handed, and the two modes
     * are only distinguishable to a newcomer if the difference is stated as
     * two different questions rather than two different nouns.
     *
     * Three states, deliberately, because there are three. A record written
     * before `window_mode` existed has the key ABSENT, which is NOT the same
     * as finite_window: treating a missing key as the historical default is
     * how a record ends up labelled with a claim nobody verified it makes.
     * Absent is reported as unknown and says so. */
    const wmode = reach.window_mode;
    if (wmode === 'steady_state' || wmode === 'finite_window') {
      const steady = wmode === 'steady_state';
      el.appendChild(h('div.caveat.info', h('span',
        h('b', steady
          ? 'THIS ANSWERS: WHAT DOES THIS TISSUE’S GEOMETRY PERMIT? '
          : 'THIS ANSWERS: WHAT DID THIS EXPOSURE WINDOW ACHIEVE? '),
        steady
          ? 'The map below is the reach this geometry permits once transport '
            + 'has settled — vessel spacing and this compound\'s transport, with '
            + 'no dosing schedule in it. No exposure window enters the '
            + 'computation, so nothing below is a statement about how long the '
            + 'drug was present. '
          : 'The map below is the exposure this run actually accumulated over '
            + 'its own window, so it is a property of the schedule and the '
            + 'window length as well as of the tissue. It is NOT the same '
            + 'object as "the reach this geometry permits" — a longer window '
            + 'or a different schedule moves it. ',
        h('b', 'These are different questions'),
        ', and an answer to one is not an answer to the other.',
        reach.window_note ? ' ' + String(reach.window_note) : '')));
    } else {
      el.appendChild(h('div.caveat', h('span',
        h('b', 'WHICH QUESTION THIS ANSWERS IS NOT RECORDED. '),
        'This run predates the field that says whether the map below is the '
        + 'reach the geometry permits (settled transport, no window) or the '
        + 'reach one finite exposure window achieved (entangled with the '
        + 'schedule). Those are different questions, and which one this is '
        + 'cannot be read off the record. It is not being assumed to be '
        + 'either.')));
    }

    /* SCOPE: this half means what it claims ONLY where the cell positions were
     * MEASURED. It is not a data-provenance footnote, it is a statement about
     * this specific figure, so it sits on top of it.
     *
     * The chain is: the exposure field's SOURCE is a density over ENDOTHELIAL
     * CELL POSITIONS; the falloff bins by distance to the nearest endothelial
     * cell; the perivascular reference is that curve's first bin; the delivery
     * spread is a percentile ratio over the same field. On a spot platform
     * (Visium and the GBmap sections) the pipeline does not measure where cells
     * sit — reconstruct.real_visium synthesises 12 pseudo-cells per spot at the
     * spot centroid with uniform jitter — so every one of those four
     * quantities is a statement about the SYNTHESIS, not about the tissue's
     * vascular architecture.
     *
     * THE ERROR CLASS, NAMED, BECAUSE IT WILL RECUR: do not confuse WHERE THE
     * PDE IS DISCRETISED with WHERE ITS SOURCE TERM COMES FROM. They feel like
     * one question — "is this computed on the grid or on the cells?" — and they
     * are not, because a grid solve inherits everything its source knows. An
     * earlier version of this very caveat argued that the exposure map survived
     * a positional refusal "because it is solved on the grid", which is true of
     * the discretisation and irrelevant to the claim. Anyone reasoning about
     * which outputs survive a positional refusal will reach for the same
     * distinction; reach for it in the right direction.
     *
     * A POSITIVE ALLOWLIST, deliberately. An unrecognised source is treated as
     * unverified rather than assumed measured, so a new adapter cannot acquire
     * the strong reading by default. As of 2026-08-21 four of the thirty-eight
     * shipped sections carry measured per-cell positions. */
    const MEASURED_POSITION_SOURCES = ['xenium', 'codex'];
    const src = ((run.report || {}).source || '').toLowerCase();
    if (!MEASURED_POSITION_SOURCES.includes(src)) {
      el.appendChild(h('div.caveat.crit', h('span',
        h('b', 'POSITIONS ON THIS SECTION ARE NOT MEASURED. '),
        `Its source is "${src || 'unknown'}", which is not one of the platforms `
        + 'this build recognises as carrying per-cell positions '
        + `(${MEASURED_POSITION_SOURCES.join(', ')}). On a spot platform the `
        + 'pipeline synthesises pseudo-cells at spot centroids, and the '
        + 'exposure field’s source term is a density over ENDOTHELIAL CELL '
        + 'positions — so the delivery spread, the falloff curve and the '
        + 'perivascular reference below all measure that synthesis rather than '
        + 'this tissue’s vasculature. Read none of them as tissue architecture. '
        + 'The potency evidence and the per-cell verdict do not depend on where '
        + 'the cells sit and are unaffected by this.')));
    }

    renderSpread(el, run);

    /* THE REFERENCE EXPOSURE, in that order: what it is, then where it came
     * from. This paragraph used to OPEN on the backend's `E_ref_source`
     * string, which is a Python expression —
     * `delivery_metrics.penetration_length_fit(...)['drug']['perivascular_mean']`
     * — printed as the first eleven words a bench scientist reads about the
     * quantity every ratio on this panel is divided by. The expression is not
     * removed: it is the traceable provenance and it is the thing an engineer
     * greps for. It is DEMOTED to a mono provenance line under the sentence
     * that says, in words, what the quantity is. */
    /* WHICH QUANTITY THE REFERENCE IS — and it is NOT the same quantity in the
     * two modes. The finite-window map normalises by an EXPOSURE (AUC, uM*s);
     * the steady map normalises by a CONCENTRATION (uM). The payload says so in
     * `E_ref_units`, and the field it arrives in is called `E_ref_uM_s` in both
     * cases, which is a name that is simply wrong in steady mode. Calling a
     * concentration an exposure on screen would be this panel repeating that
     * mistake in words, where a reader can act on it — so the word is read off
     * the units, never assumed. */
    const refUnits = reach.E_ref_units || (reach.units && reach.units.E_ref) || 'uM*s';
    const refIsConc = /^\s*uM\s*$/i.test(String(refUnits));
    el.appendChild(prose(
      'Every ratio on this panel is ',
      refIsConc ? 'some concentration' : 'some exposure',
      ' divided by the ',
      term(refIsConc ? 'reference concentration' : 'reference exposure'), ': the ',
      refIsConc ? 'settled drug concentration' : term('exposure'),
      ` in the tissue immediately beside a vessel — the first `
      + `${(reach.falloff && reach.falloff.bin_um) || 10} um distance bin from the `
      + 'nearest endothelial cell — in this same run. Because both halves of '
      + 'every ratio come from one map, the dose and the potency divide out and '
      + 'no number here carries either.'));
    if (reach.E_ref_source) {
      el.appendChild(h('div', {
        style: {
          fontSize: 'var(--fs-micro)', fontFamily: 'var(--font-mono)',
          color: 'var(--fg-2)', lineHeight: 'var(--lh-body)',
          margin: '0 0 var(--sp-4)', overflowWrap: 'anywhere',
        },
        title: 'The backend\'s own description of where this value is read '
          + 'from, verbatim. Quoted so the number can be traced to the call '
          + 'that produced it.',
      }, `source · ${String(reach.E_ref_source)}`));
    }
    if (Number.isFinite(reach.E_ref_uM_s)) {
      el.appendChild(prose(
        `On this run it is ${sci(reach.E_ref_uM_s, 3)} ${refUnits}. That absolute `
        + 'value is shown so the ratios can be unwound, not because it is '
        + 'defensible on its own — it carries the plasma level and the delivery '
        + 'law, and the ratios above are exactly what it is divided out of.'));
    }

    /* ---- the curve ----------------------------------------------------- */
    const curve = reach.curve || {};
    const cx = curve.x || [];
    const cy = curve.under_reached_fraction || [];
    if (cx.length && cy.length === cx.length) {
      el.appendChild(microLabel('under-reached fraction vs threshold'));
      /* WHAT THE READER IS LOOKING AT, BEFORE THE CHART.
       * "Under-reached fraction" is a term this project invented, and it reads
       * to a newcomer as "the fraction the drug did not reach" — a measurement.
       * It is not one: it is a count of cells below a threshold, and the
       * threshold is a convention. The chart used to arrive with only its own
       * name above it, which is precisely the case where an invented word gets
       * read as the meaning it merely sounds like. */
      el.appendChild(prose(
        'Read this as: "how much of the viable tumour sits below X times the ',
        term(refIsConc ? 'reference concentration' : 'reference exposure'),
        '?" Each row of the ', term('ladder'),
        ' below is one threshold X. The curve is the point, not any one ',
        term('rung'), ' of it — the thresholds are a reporting convention, '
        + 'nothing in the tissue or the compound picks one out.'));
      el.appendChild(lineChart({
        xs: cx, ys: cy,
        xLabel: 'threshold', yLabel: 'under-reached',
        fmtX: (v) => fmtRatio(v) + 'x', fmtY: (v) => pct(v, 0) + '%',
        markIndex: cx.indexOf(x), yMax: 1,
      }));
      // The table view. It is the accessible form of the chart above, and it is
      // also the reason nothing on this panel needs to be estimated off a
      // pixel: every plotted value is printed.
      const tbl = h('table.tbl.zebra',
        h('thead', h('tr',
          h('th', 'below'), h('th.right', 'of viable tumour'))),
        // The headline threshold is marked in the table so the reader can see
        // WHICH row the big number above was picked off. It is emphasis, not a
        // verdict — the point of printing the whole ladder is that the choice
        // of row is a convention.
        h('tbody', cx.map((v, i) => h('tr',
          { title: v === x ? 'the threshold the headline above is read at' : null },
          h('td.mono' + (v === x ? '.strong' : ''), fmtRatio(v) + 'x'),
          h('td.right.num.mono' + (v === x ? '.strong' : ''), fmtFrac(cy[i]))))));
      el.appendChild(h('div.tblwrap', { style: { margin: 'var(--sp-3) 0 0' } }, tbl));
    }

    if (reach.convention) el.appendChild(prose(String(reach.convention)));

    /* THE LADDER'S BOUND, NAMED. A bound that says where it stops holding is
     * worth more than one that implies it holds generally, so this states the
     * capacity range over which the ladder is usable AND what happens past it.
     * The principle stands; the figures it was originally written with do not
     * -- see the dated correction below. */
    /* CORRECTED 2026-08-22. All three shipped figures here were refuted by
     * direct measurement (results/reach_invariances_measure.md): the claim was
     * "<= 0.0011 over R = 1 to 39, exactly flat to R = 60, up to 0.80 over
     * 3.2-4184", and NO arm and NO reading reproduces any of the three. The
     * defect was not a bad measurement but a claim that outgrew one: read as
     * the HEADLINE RUNG the numbers were conservative and safe to rely on;
     * read as written -- "the whole ladder" -- all three understate, 0.0011 by
     * a factor of 365. Stated per-rung below, because a single max over rungs
     * reports the worst rung and the headline is the best one. */
    el.appendChild(h('div.caveat', h('span',
      reach.window_mode === 'steady_state'
        ? [h('b', 'THIS LADDER DOES NOT DEPEND ON TISSUE CAPACITY. '),
           'The engine assigns every compound a capacity of 1.0, and measured '
           + 'values span 3.2 to 4184 mL/g (n = 33). On a converged steady '
           + 'state that assignment does not matter: capacity multiplies the '
           + 'time derivative, which this map does not have. Measured across '
           + 'the full 3.2 to 4184 range, every rung moves exactly 0.000e+00 — '
           + 'it counts cells below a threshold and no cell crossed one.']
        : [h('b', 'THE BOTTOM OF THIS LADDER SURVIVES TISSUE CAPACITY; THE TOP '
             + 'OF IT DOES NOT. '),
           'The engine assigns every compound a capacity of 1.0; measured '
           + 'values span 3.2 to 4184 mL/g (n = 33, best available constant '
           + '39). Over R = 1 to 39 the headline rung F(0.1) and every rung '
           + 'below it are bit-identical — 0.000000, and flat out to R = 80. '
           + 'The WHOLE ladder moves 0.4021 over that same range, because the '
           + 'top rung F(1.0) carries essentially all of it. Across the full '
           + 'measured range the headline rung moves 0.7409 and the whole '
           + 'ladder 0.9427, both worse than the cross-compound ordering this '
           + 'product withdrew. Read the rung you are using, not the ladder.'])));

    renderFalloff(el, reach);

    /* ---- the invariance statement, from the backend --------------------
     * RESTRUCTURED, NOT REDUCED, 2026-08-22.
     *
     * This is one string and it is ~1,100 words — longer than everything else
     * on this panel put together. It was rendered as a single undifferentiated
     * `<span>`, and on a first session it was the block a reader bounced off
     * hardest: it opens on an argument, contains three numbered axes, two
     * dated self-corrections about what a previous version of itself got
     * wrong, a digression about an estimator versus a field, and closes on the
     * one sentence that most needs to be read — "an axis not named here was
     * not tested."
     *
     * A caveat nobody reads protects nobody. That is an argument about
     * PRESENTATION and never about content, so: not one word is removed, not
     * one word is reordered, and nothing is collapsed or clamped. The string
     * is split ONLY on the axis markers it already writes into itself —
     * "(1) PLASMA/DOSE SCALE:", "(2) TISSUE CAPACITY:", "(3) TRANSVASCULAR
     * PENETRANCE:" — so each swept axis becomes its own paragraph under its
     * own name, in the original order, with the lead-in and the closing
     * sentence kept as their own paragraphs around them.
     *
     * If the backend rewords those markers the regex stops matching and the
     * whole string renders as one block again, exactly as it does today. That
     * is the intended failure: degrade to the current behaviour, never to a
     * partial caveat. `splitAxes` returns a single-element array on no match
     * and the join is asserted lossless in tests/test_first_session.py. */
    if (reach.scale_invariance) {
      const full = String(reach.scale_invariance);
      const parts = splitAxes(full);
      el.appendChild(h('div.caveat.info', { title: full },
        parts.map(({ head, body }) => h('div', { style: { margin: '0 0 var(--sp-3)' } },
          head ? h('div', {
            style: {
              fontSize: 'var(--fs-micro)', color: 'var(--fg-2)',
              fontFamily: 'var(--font-mono)', letterSpacing: 'var(--tr-wide)',
              margin: '0 0 var(--sp-1)',
            },
          }, head) : null,
          h('span', body)))));
    }
    if (reach.penetration_fit_error) {
      el.appendChild(h('div.caveat', h('span',
        h('b', 'Falloff fit did not converge: '), String(reach.penetration_fit_error),
        ' The reference exposure above is a BIN MEAN and does not read a fitted '
        + 'parameter, so it is unaffected; the penetration LENGTH is the number '
        + 'this refusal applies to.')));
    }

    renderPotency(el, run);
  }

  return { el, render };
}
