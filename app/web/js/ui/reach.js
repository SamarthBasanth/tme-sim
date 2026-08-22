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
import { pct, num, NULL_GLYPH } from '../util/fmt.js';

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
 *  through here so none of them can drift below it. */
function prose(text, cls) {
  return h('p' + (cls ? '.' + cls : ''), {
    style: {
      margin: '0 0 var(--sp-4)',
      fontSize: 'var(--fs-disclosure)',
      lineHeight: 'var(--lh-body)',
      color: 'var(--fg-1)',
    },
  }, text);
}

function microLabel(text) {
  return h('div', {
    style: {
      fontSize: 'var(--fs-micro)', color: 'var(--fg-2)',
      fontFamily: 'var(--font-mono)', letterSpacing: 'var(--tr-wide)',
      textTransform: 'uppercase', margin: 'var(--sp-4) 0 var(--sp-2)',
    },
  }, text);
}

// Leading zero KEPT deliberately. Stripped, `0.1x` renders as `.1x`, which in a
// mono column sits one glyph away from `1x` — and the difference between them
// is a factor of ten in the quantity this whole section is about.
const fmtRatio = (v) => String(v);
const fmtFrac = (v) => (Number.isFinite(v) ? pct(v, 1) + '%' : NULL_GLYPH);

export function createReachPanel() {
  const el = h('div');

  function renderEmpty(msg) {
    clear(el);
    el.appendChild(prose(msg, 'dim'));
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
      const joint = mech.joint_limited;
      if (joint && Number.isFinite(joint.fraction_of_exposure_limited)) {
        el.appendChild(h('div.caveat', h('span',
          h('b', 'DOSE IS NOT A CLEAN COUNTERFACTUAL HERE. '),
          `${fmtFrac(joint.fraction_of_exposure_limited)} of the exposure-limited `
          + 'survivals are ALSO at or below the resistance threshold. They are '
          + 'labelled by their exposure band because below IC50 the exposure is '
          + 'the binding constraint, but raising the dose alone would not have '
          + 'killed them. Read the verdict as "exposure is the first thing in '
          + 'the way", not as "more drug would fix it".')));
      }
    }

    /* ---- DELIVERY EVIDENCE — the half that rules delivery in or out ----- */
    el.appendChild(microLabel('delivery evidence'));

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
      h('b', 'DELIVERY SPREAD'), ' — the best-perfused tenth of tumour cells '
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

    el.appendChild(prose(
      `The reference is ${reach.E_ref_source || 'the perivascular mean of this run'}`
      + `. It is an exposure, in ${(reach.units && reach.units.E_ref) || 'uM.s'}, `
      + `measured over tissue in the first ${(reach.falloff && reach.falloff.bin_um) || 10} um `
      + 'distance bin from the nearest endothelial cell. Every number in this '
      + 'section is that exposure divided into another one, so the ratio carries '
      + 'no dose and no potency.'));

    /* ---- the curve ----------------------------------------------------- */
    const curve = reach.curve || {};
    const cx = curve.x || [];
    const cy = curve.under_reached_fraction || [];
    if (cx.length && cy.length === cx.length) {
      el.appendChild(microLabel('under-reached fraction vs threshold'));
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
     * Measured on the flagship section: <= 0.0011 movement over R = 1 to 39,
     * flat to R = 60, and up to 0.80 over the full measured 3.2-4184 range. */
    el.appendChild(h('div.caveat', h('span',
      h('b', 'THIS LADDER IS USABLE FOR TISSUE CAPACITY R ≲ 60, AND NOT BEYOND '
        + 'IT. '),
      'The engine assigns every compound a capacity of 1.0. Over R = 1 to 39 — '
      + '1 to the best available constant from 33 measured values — the whole '
      + 'ladder moves by at most 0.0011, and it is exactly flat out to R = 60. '
      + 'Past that it degrades fast: across the full measured range of 3.2 to '
      + '4184 mL/g it moves by up to 0.80, which is worse than the '
      + 'cross-compound ordering this product withdrew. The bound is stated '
      + 'rather than implied because a reader who assumes it holds generally '
      + 'would be making the same mistake the withdrawn ordering made.')));

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
      // collapsing 279.7 -> 76.7 um. The under-reached ladder is nearly immune
      // over the same range (<= 0.0011); the SHAPE is not. This is exactly the
      // asymmetry the panel must not blur, because the two figures sit adjacent
      // and look equally solid.
      el.appendChild(h('div.caveat.crit', h('span',
        h('b', 'THE SHAPE IS NOT CAPACITY-INVARIANT — the ladder above nearly '
          + 'is, and this is not. '),
        'The engine assigns every compound a tissue capacity of 1.0; measured '
        + 'values span 3.2 to 4184 mL/g (n = 33, best available constant 39). A '
        + 'capacity term rescales time rather than level, so it cancels from a '
        + 'ratio only where the field has CONVERGED, and at the shipping 6 h '
        + 'window it has not — the gap to the infinite-window limit is 0.01302 '
        + 'and decays as 281/T with no plateau, so about 78 h would be needed '
        + 'to reach 0.001. At R = 39 this curve moves 0.3772 absolute and the '
        + 'fitted penetration length collapses from 279.7 to 76.7 um. Read the '
        + 'SHAPE as carrying an unquantified capacity term of that size. The '
        + 'under-reached ladder above moves at most 0.0011 over the same range '
        + 'and is flat to R = 60, which is why the two are not presented as '
        + 'equally solid.')));
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

    /* ---- WHICH OBJECT THIS IS ------------------------------------------
     * The payload has always known whether the map came from a converged
     * steady state or from a finite exposure window; until now nothing put it
     * on the screen. Those are different claims — "the reach this geometry
     * permits" versus "the reach this exposure achieved" — and a panel that
     * renders them identically invites the reader to carry a conclusion from
     * one to the other.
     *
     * Three states, deliberately, because there are three. A record written
     * before `window_mode` existed has the key ABSENT, which is NOT the same
     * as finite_window: treating a missing key as the historical default is
     * how a record ends up labelled with a claim nobody verified it makes.
     * Absent is reported as unknown and says so. */
    const wmode = reach.window_mode;
    if (wmode === 'steady_state' || wmode === 'finite_window') {
      const label = wmode === 'steady_state' ? 'STEADY STATE' : 'FINITE WINDOW';
      const gist = wmode === 'steady_state'
        ? 'the reach this geometry permits once transport has converged — no '
          + 'exposure window enters the computation'
        : 'the reach achieved over one finite exposure window, which is a '
          + 'property of that window as well as of the geometry';
      el.appendChild(h('div.caveat.info', h('span',
        h('b', label + ' — '), gist,
        reach.window_note ? '. ' + String(reach.window_note) : '.')));
    } else {
      el.appendChild(h('div.caveat', h('span',
        h('b', 'Object not recorded. '),
        'This run predates the field that says whether the map is a converged '
        + 'steady state or a finite-window exposure, so which one it is cannot '
        + 'be read off the record. It is not being assumed to be either.')));
    }

    /* ---- the invariance statement, verbatim from the backend ----------- */
    if (reach.scale_invariance) {
      el.appendChild(h('div.caveat.info', h('span', String(reach.scale_invariance))));
    }
    if (reach.penetration_fit_error) {
      el.appendChild(h('div.caveat', h('span',
        h('b', 'Falloff fit did not converge: '), String(reach.penetration_fit_error),
        ' The reference exposure above is a BIN MEAN and does not read a fitted '
        + 'parameter, so it is unaffected; the penetration LENGTH is the number '
        + 'this refusal applies to.')));
    }

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
    }
  }

  return { el, render };
}
