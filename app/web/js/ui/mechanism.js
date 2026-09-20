// ui/mechanism.js — the two pieces this feature adds:
//
//   1. the MECHANISM BREAKDOWN: every attribution category with its count and
//      its share, the delivery / dose / resistance / chance roll-up as one
//      stacked bar, and the bit-exact reconciliation against the engine's own
//      kill_fraction. Rows toggle their category in the viewer.
//   2. the TIME-LAPSE SCRUBBER: play / pause / step / seek over the captured
//      drug-concentration frames.
//
// (1) used to be a floating card parked in the stage's top-right corner. It is
// now DOCKED into the READ rail — see the note above createMechanismLegend.
// (2) is still a stage overlay: it is a transport control for the picture, and
// it has to be under the picture it scrubs.
//
// Design rule being honoured (API_CONTRACT.md): every category carries a
// visible label AND a number. The categorical palette clears its CVD target on
// the dark surround but three slots sit under 3:1 on the light one, so
// colour-alone identity is not available and is not used.

import { h, fill } from '../util/dom.js';
import { num, int, pct, NULL_GLYPH } from '../util/fmt.js';
import { MECH, MECH_META, LIMITER_LABEL, dominantReason } from '../mechanism.js';

/* ==================================================== breakdown ======== */

// WHY THIS IS NO LONGER A STAGE OVERLAY
// -------------------------------------
// docs/UI_ARCHITECTURE.md 2.3 booked the legend's 300px as a permanent cost
// against the specimen and named the escape hatch as "a collapse control on
// the card, not a mode coupling". Both alternatives were wrong in the same
// way: they treated a 300x380 px card sitting ON TOP of the tissue as a width
// problem. It is not. The stage is the product's primary output and the card
// OCCLUDED the top-right quarter of it in every colour mode, including the two
// modes it is not the key for.
//
// A mode coupling was rejected because "tying it to mech mode would hide a
// result behind a display setting" — and that reasoning is still right, which
// is exactly why the answer is not to hide it but to MOVE it. The panel is
// read once per run (what was killed, and what limited the survivors) and then
// used episodically (toggling a category to isolate it spatially). Once-per-
// run reference numbers belong in the READ rail — that rail's whole job is
// "what does the run say" — and specifically in RESULT, immediately under
// METRICS, because this breakdown IS the decomposition of the `Kill fraction`
// metric printed two lines above it and reconciles with it bit-exactly.
//
// There is precedent for the interactive half too: `type` colouring has no
// stage key either. Its per-channel visibility toggles live in the rail under
// VIEW > Channels. Mechanism now matches it.
//
// Nothing is hidden by this: the rail section is `open` by default, RESULT is
// the group the rail opens on, and selecting the `mech` colour source scrolls
// the rail to it (main.js). What the stage gets back is a quarter of itself.
export function createMechanismLegend(viewer, { onToggle } = {}) {
  const title = h('div.m-title');
  const headline = h('span.m-headline');
  const filter = h('div.m-filter');
  const rows = h('div.m-rows');
  const bar = h('div.m-bar');
  const note = h('div.m-note');
  const body = h('div.m-body', filter, rows, bar, note);

  const card = h('div.mechleg.mechleg-docked',
    { role: 'figure', 'aria-label': 'Mechanism attribution breakdown' },
    title, body);

  let agg = null;
  /** Row element per mechanism code, so a full repaint can put the keyboard
   *  back where it was. A control that loses focus when you activate it sends
   *  a keyboard user to the top of a 2,200 px rail. */
  const rowEls = [];

  /** Every category the viewer is currently drawing / not drawing. */
  function hiddenCodes() {
    const out = [];
    for (let i = 0; i < MECH_META.length; i += 1) {
      if (viewer.labelVis && !viewer.labelVis[i]) out.push(i);
    }
    return out;
  }

  function showAllCategories() {
    for (let i = 0; i < viewer.labelVis.length; i += 1) viewer.setLabelVisible(i, true);
    render();
    if (onToggle) onToggle(null, true);
  }

  /** Draw ONLY `code`. Same gesture, same words and same keys as the channel
   *  solo in js/ui/inspect.js — two controls that do the same thing get the
   *  same label (docs/UI_ARCHITECTURE.md section 5). */
  function isolate(code) {
    for (let i = 0; i < viewer.labelVis.length; i += 1) {
      viewer.setLabelVisible(i, i === code);
    }
    render();
    if (rowEls[code] && typeof rowEls[code].focus === 'function') {
      rowEls[code].focus();
    }
    if (onToggle) onToggle(code, true);
  }

  /* --- WHAT THE VIEWER IS ACTUALLY DRAWING, AND HOW TO UNDO IT -----------
   *
   * The round trip this panel is half of: clicking a survival row on THE
   * FINDING colours the tissue by mechanism, hides every OTHER category, and
   * — because main.js routes it through setColourMode — scrolls the READ rail
   * down to this section. So a reader who clicks one row at the top of the
   * rail arrives here, several screens from what they clicked, looking at a
   * picture that is no longer showing most of its cells.
   *
   * What told them so was `opacity: 0.36` on the rows that had gone, and a
   * `show all` ghost button in the title strip that appeared only when
   * something was hidden. Two problems, and the first is a house rule:
   *
   *   - COLOUR (here, dimming) WAS THE ONLY CHANNEL saying a category is
   *     hidden. docs/UI_ARCHITECTURE.md 4.5: never. Every row now says the
   *     word `hidden` and carries `aria-pressed`, so the state is in the text
   *     and in the accessibility tree, not only in the alpha.
   *   - THE UNDO DID NOT SAY WHAT IT WOULD UNDO. "show all" beside a card with
   *     six dimmed rows is an instruction, not an explanation of the state the
   *     reader has arrived in.
   *
   * This band states the filter in a sentence, in counts, names the one
   * gesture that clears it, and is a real <button> so the whole thing is on
   * the keyboard path with the global focus ring. It exists only while a
   * filter is applied, so nothing new appears in the ordinary case.
   *
   * It is separate from render() because a row toggle has to repaint it
   * WITHOUT rebuilding the rows: a full repaint replaces the element the
   * keyboard is on, which drops focus to the top of a 2,200 px rail.
   */
  function renderFilter() {
    const counts = (agg && agg.counts) || {};
    const nTumor = (agg && agg.n_tumor) || 0;
    const hidden = hiddenCodes();
    const anyHidden = hidden.length > 0;

    // No 'mechanism' word in the title: the section's own <summary> says it
    // one line above, and a card that repeats its container's title is how a
    // rail ends up with two headings per section.
    fill(title,
      headline,
      h('span.grow'),
      h('span.mono', {
        style: { color: anyHidden ? 'var(--accent)' : 'var(--fg-2)' },
        title: anyHidden
          ? 'How many mechanism categories the specimen is currently drawing. '
            + 'It is a property of the PICTURE, not of the run: every count '
            + 'below is the whole run either way.'
          : 'The specimen is drawing every mechanism category.',
      }, anyHidden
        ? `${MECH_META.length - hidden.length}/${MECH_META.length} shown`
        : 'all shown'));

    // Emptied AND unstyled: `fill` removes the children, not the inline
    // background and rule a previous render set, so a cleared filter would
    // otherwise leave an empty accent band behind it.
    fill(filter);
    filter.removeAttribute('style');
    if (!anyHidden) return;

    const shownCodes = [];
    for (let i = 0; i < MECH_META.length; i += 1) {
      if (hidden.indexOf(i) === -1) shownCodes.push(i);
    }
    const shownNames = shownCodes.map((c) => MECH_META[c].short).join(', ');
    const shownCells = shownCodes.reduce(
      (s, c) => s + (counts[MECH_META[c].name] || 0), 0);
    /* Styled from JS: this module owns no stylesheet, every value below is a
     * token, and no colour literal appears. `--accent` is correct rather than
     * semantic here — a filter is a SELECTION, which is exactly what the one
     * accent is reserved for; an amber band would say the run is in a warning
     * state, and it is not. */
    filter.style.background = 'var(--accent-bg)';
    filter.style.borderLeft = '3px solid var(--accent)';
    filter.style.padding = 'var(--sp-3) var(--sp-5)';
    filter.style.margin = '0 0 var(--sp-4)';
    fill(filter,
      h('div.m-filter-hd', {
        style: {
          fontFamily: 'var(--font-mono)', fontSize: 'var(--fs-md)',
          fontWeight: '600', color: 'var(--accent)',
          lineHeight: 'var(--lh-tight)', overflowWrap: 'anywhere',
        },
      }, shownCodes.length === 1
        ? `SHOWING ONLY ${shownNames}`
        : `SHOWING ${shownCodes.length} OF ${MECH_META.length} CATEGORIES`),
      h('div.m-filter-bd', {
        style: {
          fontSize: 'var(--fs-disclosure)', color: 'var(--fg-1)',
          lineHeight: 'var(--lh-body)', margin: 'var(--sp-2) 0 var(--sp-3)',
        },
      },
      shownCodes.length
        ? `The specimen is drawing ${int(shownCells)} of ${int(nTumor)} `
          + `tumour cells — ${shownNames}. `
        : 'The specimen is drawing no tumour cells at all: every category is '
          + 'hidden. ',
      'The counts and shares below are the whole run and do not change with '
      + 'this filter; only the picture does.'),
      h('button.btn.ghost.tiny', {
        type: 'button',
        'aria-label': `Show all ${MECH_META.length} mechanism categories in `
          + 'the specimen again',
        title: 'Draw every mechanism category in the specimen again. This '
          + 'undoes an isolate — from a row on The finding, from a row below, '
          + 'or from a Selected cell — in one action.',
        onclick: showAllCategories,
      }, 'show all categories'));
  }

  function render() {
    if (!agg) {
      // An empty state, not a disappearance. A reconstruction-only run has
      // nothing to attribute and the section says so where the numbers would
      // be — a section that silently vanishes reads as a bug.
      fill(title);
      fill(filter);
      fill(rows, h('div.m-empty',
        'This run carries no per-cell attribution. Submit a run with '
        + 'attribute_mechanism (the default) to get one.'));
      fill(bar);
      bar.style.display = 'none';
      fill(note);
      return;
    }

    const counts = agg.counts || {};
    const nTumor = agg.n_tumor || 0;
    const killed = counts.KILLED || 0;
    const nSurv = Number.isFinite(agg.n_survived) ? agg.n_survived : null;

    /* THE HEADLINE OF THIS SECTION IS NO LONGER A KILL PERCENTAGE. 2026-09-20.
     *
     * It was `"11.2% killed"`, in `.m-headline`, first line of the section —
     * i.e. an ABSOLUTE KILL FRACTION as the headline of a panel, with no bound
     * anywhere on it. That is the one thing this project forbids, and both of
     * the other surfaces that print the number already refuse it: js/ui/
     * inspect.js renders `DO NOT QUOTE THIS NUMBER` under the metric, and
     * js/ui/finding.js prints the count only as the denominator arithmetic
     * with its 0.3252-absolute state-route bound attached.
     *
     * The number is NOT removed — it is a real output and it is the
     * denominator the survivor shares below are complements of, so hiding it
     * would be its own dishonesty. What changes is that the headline slot now
     * carries what this section is FOR: how many cells survived, which is what
     * every row and every segment below is a share of. The kill count, its
     * denominator and its bound go into the note block at the foot of the
     * panel, where the other limits on these numbers already live. */
    headline.textContent = nSurv !== null
      ? `${int(nSurv)} survivors of ${int(nTumor)} tumour cells`
      : `${int(nTumor)} tumour cells`;

    renderFilter();

    /* --- one row per category ---------------------------------------- */
    const order = [
      MECH.KILLED,
      MECH.SURVIVED_UNREACHED, MECH.SURVIVED_UNDERDOSED,
      MECH.SURVIVED_RESISTANT, MECH.SURVIVED_STOCHASTIC,
      MECH.UNATTRIBUTED, MECH.ALREADY_DEAD, MECH.NOT_TUMOR,
    ];
    const kids = [];
    rowEls.length = 0;
    for (const code of order) {
      const meta = MECH_META[code];
      const n = counts[meta.name] || 0;
      // The denominator differs by category and must not be fudged: the
      // outcome categories are shares of the tumour population, the survival
      // categories are shares of the survivors, and NOT_TUMOR is neither.
      const denom = meta.limiter ? (agg.n_survived || 0) : (code === MECH.NOT_TUMOR ? 0 : nTumor);
      if (code === MECH.SURVIVED_UNREACHED) kids.push(h('div.m-sep'));
      const on = !!(viewer.labelVis && viewer.labelVis[code]);
      /* WHAT THE PERCENTAGE IN THIS ROW IS A PERCENTAGE OF — and, for the one
       * row that has no denominator, why its cell holds a dash.
       *
       * Three different denominators share one column here (tumour cells,
       * survivors, none), the difference was stated only in a source comment,
       * and NOT_TUMOR rendered a bare em dash in a column of percentages —
       * the shape of a missing number, a zero and a refusal at once. The
       * column is 42 px wide and cannot hold a sentence, so the sentence goes
       * where there is room for it: on the row, for the pointer and for the
       * accessibility tree. */
      const share = denom > 0
        ? `${int(n)} of ${int(denom)} `
          + `${meta.limiter ? 'surviving tumour cells' : 'tumour cells'} `
          + `(${pct(n / denom, 1)}%).`
        : (code === MECH.NOT_TUMOR
          ? `${int(n)} cells that are not tumour. No share is shown because `
            + 'this category is not a share of the tumour or of the '
            + 'survivors — there is no denominator to divide by, so the '
            + 'column holds a dash rather than a number.'
          : `${int(n)} cells. No share is shown: this run reports no `
            + 'denominator for this category, and none is assumed.');
      const desc = ((agg.label_descriptions || {})[meta.name] || meta.name)
        + `\n\n${share}`;
      /* THE ROW SAYS ITS OWN STATE, IN A WORD AND IN ARIA.
       *
       * It was `role="button"` with no `aria-pressed`, so a screen reader
       * announced a button and not whether it was on — and `.off` is
       * `opacity: 0.36`, which makes the visible state colour-alone. Both are
       * fixed by the same two additions, and the `hidden` marker below is the
       * non-colour channel docs/UI_ARCHITECTURE.md 4.5 requires. */
      const row = h('div.m-row', {
        role: 'button',
        tabindex: '0',
        'aria-pressed': String(on),
        title: `${desc}\n\nClick, or Enter, to ${on ? 'hide' : 'show'} this `
          + 'category in the specimen. Alt-click, or Alt+Enter, shows only '
          + 'this one. Neither changes a count below — the numbers are the '
          + 'whole run either way.',
      },
        h('span.m-sw', { style: { background: `var(${meta.var})` } }),
        h('span.m-name', {
          style: meta.alarm && n > 0 ? { color: 'var(--crit)' } : null,
        }, meta.short),
        h('span.m-n', int(n)),
        h('span.m-p', denom > 0 ? pct(n / denom, 1) + '%' : NULL_GLYPH),
      );
      if (!n) row.classList.add('zero');
      if (!on) row.classList.add('off');
      // The word, beside the dimming, so "hidden from the picture" is legible
      // without reading an alpha value. Removed and re-added by flip() rather
      // than left in the DOM at zero opacity.
      const marker = h('span.m-hidden', {
        style: {
          fontFamily: 'var(--font-mono)', fontSize: 'var(--fs-micro)',
          color: 'var(--fg-2)', letterSpacing: 'var(--tr-wide)',
        },
      }, 'hidden');
      if (!on) row.appendChild(marker);
      /* Updated IN PLACE, and the filter band is repainted on its own. The
       * obvious implementation — call render() — replaces the element the
       * keyboard is standing on, which sends focus to the body and the reader
       * to the top of a 2,200 px rail every time they toggle a row. */
      const flip = () => {
        const next = !viewer.labelVis[code];
        viewer.setLabelVisible(code, next);
        row.classList.toggle('off', !next);
        row.setAttribute('aria-pressed', String(next));
        row.title = `${desc}\n\nClick, or Enter, to ${next ? 'hide' : 'show'} `
          + 'this category in the specimen. Alt-click, or Alt+Enter, shows '
          + 'only this one. Neither changes a count below — the numbers are '
          + 'the whole run either way.';
        if (!next) row.appendChild(marker);
        else if (marker.parentNode) marker.parentNode.removeChild(marker);
        renderFilter();
        if (onToggle) onToggle(code, next);
      };
      row.addEventListener('click', (e) => {
        if (e.altKey) { e.preventDefault(); isolate(code); return; }
        flip();
      });
      row.addEventListener('keydown', (e) => {
        if (e.key !== 'Enter' && e.key !== ' ') return;
        e.preventDefault();
        if (e.altKey) isolate(code);
        else flip();
      });
      rowEls[code] = row;
      kids.push(row);
    }
    /* Both routes to the isolate gesture, named on the panel rather than left
     * to be discovered. Same shape and almost the same words as the channel
     * solo hint in js/ui/inspect.js, because it is the same gesture doing the
     * same thing to a different category set — and until now the mechanism
     * rows carried no hint at all, so "click a row to hide it" was a gesture
     * a reader could only find by trying it. */
    kids.push(h('div.m-hint', {
      style: {
        fontSize: 'var(--fs-disclosure)', color: 'var(--fg-2)',
        lineHeight: 'var(--lh-body)', margin: 'var(--sp-3) 0 0',
      },
      title: 'Isolate draws exactly one mechanism category in the specimen '
        + 'and hides the others. It is the same action a survival row on The '
        + 'finding performs, and "show all categories" undoes it in one '
        + 'click. It changes the picture only; every count on this panel is '
        + 'the whole run.',
    }, 'click a row to hide it · isolate one: alt-click · alt+enter'));
    fill(rows, ...kids);

    /* --- the roll-up: what fraction of survival is what --------------- */
    const b = agg.survival_breakdown || {};
    const segs = [];
    for (const code of [MECH.SURVIVED_UNREACHED, MECH.SURVIVED_UNDERDOSED,
      MECH.SURVIVED_RESISTANT, MECH.SURVIVED_STOCHASTIC, MECH.UNATTRIBUTED]) {
      const meta = MECH_META[code];
      const v = b[meta.limiter] || 0;
      if (v <= 0) continue;
      segs.push(h('span', {
        style: { width: (v * 100).toFixed(3) + '%', background: `var(${meta.var})` },
        title: `${LIMITER_LABEL[meta.limiter]} ${(v * 100).toFixed(1)}% of survivors`,
      }));
    }
    fill(bar, ...segs);
    bar.style.display = segs.length ? 'flex' : 'none';

    /* --- the sentence, and the audit --------------------------------- */
    const dom = dominantReason(agg);
    const rec = agg.reconciliation || {};
    const joint = agg.joint_limited || {};
    const kids2 = [];
    /* THE SENTENCE LEADS ON THE DECOMPOSITION, NOT ON THE KILL RATE.
     *
     * It opened `<b>11.2%</b> of tumour cells killed.` — a bold absolute kill
     * percentage, with no denominator and no bound, as the first words of the
     * note under a bar chart of survivors. The share is what this panel is
     * about; the kill count is its complement, and it is stated below as a
     * count over its denominator with the limit that governs it, exactly as
     * js/ui/finding.js already states it. Nothing is removed. */
    kids2.push(h('div',
      dom
        ? h('span', 'Of the survivors, the largest share is ',
          h('b', { style: { color: `var(${MECH_META.find((m) => m.limiter === dom.limiter).var})` } },
            dom.label), '-limited (', pct(dom.fraction, 1), '% of ',
          nSurv !== null ? int(nSurv) : NULL_GLYPH, ' survivors).')
        : h('span', 'No survivors to attribute.')));

    if (joint.count > 0) {
      kids2.push(h('div', { style: { marginTop: '2px' } },
        int(joint.count), ' exposure-limited cells are ALSO below the ',
        'resistance threshold — raising the dose alone would not have killed ',
        'them. They are counted under their exposure band.'));
    }
    if ((counts.UNATTRIBUTED || 0) > 0) {
      kids2.push(h('div.m-bad', { style: { marginTop: '2px' } },
        int(counts.UNATTRIBUTED), ' survivals are UNATTRIBUTED — not explained ',
        'by the bands, and deliberately NOT reassigned to the nearest bucket.'));
    }
    const th = agg.thresholds || {};
    kids2.push(h('div', { style: { marginTop: '2px', color: 'var(--fg-3)' } },
      `bands: E/IC50 < ${num(th.unreached_ratio, 2)} unreached · `,
      `< ${num(th.adequate_ratio, 2)} underdosed · `,
      `modulation <= ${num(th.resistant_modulation, 2)} resistant`));
    kids2.push(h('div', { class: rec.bit_exact ? 'm-ok' : 'm-bad' },
      rec.bit_exact
        ? `reconciles bit-exact with kill_fraction = ${num(rec.kill_fraction_from_engine, 6)}`
        : `DOES NOT RECONCILE (delta ${rec.abs_difference})`));
    fill(note, ...kids2);

    /* THE KILL COUNT, WITH ITS DENOMINATOR AND ITS BOUND, where the other
     * limits on this panel live — and never as the headline.
     *
     * Same treatment js/ui/finding.js gives it, in the same words, for the
     * same reason: this panel prints a kill number, so a reader must not be
     * able to read that number without meeting the bound on it. `.m-bad` is
     * not used — this is not a failure of the run, it is a limit on how one of
     * its outputs may be quoted — so it takes the caveat severity the rest of
     * the rail uses for exactly that. */
    note.appendChild(h('div.caveat.crit', { style: { marginTop: 'var(--sp-4)' } },
      h('span',
        h('b', 'THE KILL COUNT IS THE DENOMINATOR HERE, NOT A RESULT. '),
        `${int(killed)} of ${int(nTumor)} tumour cells are labelled KILLED, `
        + 'which is what leaves the survivor count this panel decomposes. Do '
        + 'not quote it as a kill rate. It is an upper bound with no intrinsic '
        + 'resistance modelled, and on a section that cannot measure hypoxia '
        + 'it moves by up to 0.3252 absolute (1.77x) purely with which '
        + 'cell-state route is chosen, everything else held identical. The '
        + 'shares above are ratios among survivors and do not inherit that '
        + 'movement. The full statement is in Caveats, below, unedited.')));
  }

  return {
    node: card,
    /** @param {object|null} aggregate the TMSM JSON block */
    set(aggregate) { agg = aggregate; render(); },
    refresh: render,
    get aggregate() { return agg; },
  };
}

/* ==================================================== scrubber ========= */

/** Frames per second when the time-lapse is playing. */
const PLAY_FPS = 8;

/**
 * Pick the frame the series should OPEN on.
 *
 * The series starts at t=0, which by construction is empty — no drug has
 * entered the tissue yet. Opening there is defensible arithmetic and terrible
 * evidence: the first thing a visitor sees on a view sold as "where does the
 * drug go" is a frame with no drug in it, and the panel underneath has to
 * apologise for it in orange.
 *
 * So open on the frame that carries the MOST DRUG — argmax of the mean over
 * the sampled grid. Two reasons for that statistic rather than "the last
 * frame":
 *
 *   - with a static plasma source the concentration climbs to a quasi-steady
 *     state and never clears, so the mean is monotone and the argmax IS the
 *     last frame: maximum penetration, which is the honest headline;
 *   - with time_resolved=true the drug enters AND clears, so the last frame is
 *     empty for the same reason the first one is. The argmax lands on peak
 *     exposure instead, and the same rule keeps working.
 *
 * `frame.max` alone would not do: it peaks the moment the vessel-adjacent
 * voxel saturates, which can be long before the drug has gone anywhere.
 *
 * Ties resolve to the LATEST frame, so a perfectly flat series still opens on
 * the end of the window rather than at t=0.
 *
 * Returning 0 is a legitimate answer (a series that never carries any drug),
 * and then nothing has been misrepresented: t=0 is genuinely what there is.
 */
export function openingFrameIndex(series) {
  if (!series || !series.n_frames) return 0;
  let best = 0, bestMean = -Infinity;
  for (let k = 0; k < series.frames.length; k++) {
    const v = series.frames[k].values;
    if (!v || !v.length) continue;
    let s = 0;
    for (let j = 0; j < v.length; j++) s += v[j];
    const mean = s / v.length;
    if (mean >= bestMean) { bestMean = mean; best = k; }
  }
  return bestMean > 0 ? best : 0;
}

export function createTimelapseScrubber({ onFrame }) {
  // Glyph-only buttons: the character itself is not a usable accessible name
  // (a screen reader announces the Unicode character, e.g. "black right-
  // pointing triangle", not "play"), and `title` alone is not dependably
  // exposed as one either. aria-label carries the real name and tracks the
  // action the button performs, per docs/UI_ARCHITECTURE.md section 5 ("a
  // toggle states the state it produces").
  const btnPlay = h('button.btn.ghost.sm.s-btn', { title: 'Play / pause (space)', 'aria-label': 'Play' }, '▶');
  const btnStop = h('button.btn.ghost.sm.s-btn', {
    title: 'Jump to the first frame — t = 0, before any drug has entered the tissue',
    'aria-label': 'Back to first frame, t = 0',
  }, '⏮');
  const range = h('input', {
    type: 'range', min: '0', max: '0', step: '1', value: '0',
    'aria-label': 'Time-lapse frame',
  });
  const read = h('div.s-read', NULL_GLYPH);
  const meta = h('div.s-meta');
  const card = h('div.ov.ov-timelapse.hidden',
    h('div.ov-card.scrub', btnStop, btnPlay, range, read, meta));

  let series = null;      // the decoded TMSL
  let info = null;        // RunResult.timelapse metadata
  let i = 0;
  let timer = 0;

  function paint() {
    if (!series) return;
    const f = series.frames[i];
    read.textContent = '';
    // The series no longer opens at t=0, so the clock is the single most
    // load-bearing character on this control: a viewer who thinks the frame in
    // front of them is t=0 misreads the penetration by the whole window. It is
    // therefore set apart — its own type size, its own colour, `t =` spelled
    // out — rather than being one more grey token in a run-on line.
    read.appendChild(h('span.s-t', `t = ${fmtTime(f.time_s)}`));
    read.appendChild(h('span.unit', ` · frame ${i + 1}/${series.n_frames}`));
    read.appendChild(h('span.unit', ` · max ${num(f.max, 3)} uM`));
    // t=0 is a legitimate place to stand, and it is empty by construction.
    // Say which of the two it is rather than letting a blank field read as a
    // broken overlay.
    if (f.max <= 0) read.appendChild(h('span.s-zero', ' · no drug yet'));
    range.value = String(i);
    if (onFrame) onFrame(f, i, series);
  }

  function go(next) {
    if (!series) return;
    i = Math.max(0, Math.min(series.n_frames - 1, next));
    paint();
  }

  function play() {
    if (!series || timer) return;
    if (i >= series.n_frames - 1) i = 0;
    btnPlay.textContent = '❚❚';
    btnPlay.setAttribute('aria-label', 'Pause');
    timer = setInterval(() => {
      if (i >= series.n_frames - 1) { pause(); return; }
      go(i + 1);
    }, 1000 / PLAY_FPS);
  }

  function pause() {
    if (timer) clearInterval(timer);
    timer = 0;
    btnPlay.textContent = '▶';
    btnPlay.setAttribute('aria-label', 'Play');
  }

  btnPlay.addEventListener('click', () => (timer ? pause() : play()));
  btnStop.addEventListener('click', () => { pause(); go(0); });
  range.addEventListener('input', () => { pause(); go(Number(range.value)); });

  return {
    node: card,
    get playing() { return !!timer; },
    toggle() { return timer ? pause() : play(); },
    pause,
    step(d) { pause(); go(i + d); },
    clear() {
      pause();
      series = null; info = null; i = 0;
      card.classList.add('hidden');
    },
    /**
     * @param {object} decoded  the decoded TMSL payload
     * @param {object} metadata RunResult.timelapse (carries the caveats)
     */
    set(decoded, metadata) {
      pause();
      series = decoded; info = metadata || {}; i = 0;
      if (!series || !series.n_frames) { card.classList.add('hidden'); return; }
      card.classList.remove('hidden');
      range.max = String(series.n_frames - 1);
      i = openingFrameIndex(series);
      const bits = [
        `${series.n_frames} frames`,
        `every ${fmtTime(info.realised_interval_s || series.interval_s)}`,
        `over ${fmtTime(series.window_s)}`,
        `${(series.bytes / 1e6).toFixed(1)} MB`,
        `drug concentration (uM), ${info.grid ? info.grid.join('x') : '?'} grid`,
      ];
      fill(meta, h('span', bits.join(' · ')));
      // State the default out loud, next to the control that set it. The
      // series is showing a chosen frame, not its beginning, and the operator
      // is told so in words before they read a number off the picture.
      if (i > 0) {
        meta.appendChild(h('span.s-open',
          ` — opened at t = ${fmtTime(series.frames[i].time_s)} `
          + `(frame ${i + 1}/${series.n_frames}), the frame carrying the most `
          + 'drug. It is NOT the start of the window: press ⏮ for t = 0, '
          + 'before dosing.'));
      }
      if (info.plateau_note) {
        meta.appendChild(h('span.warn', ' — static plasma source: the drug is ' +
          'seen to enter and spread toward a quasi-steady state, but never to ' +
          'clear. Set time_resolved=true to watch it fall.'));
      }
      if (info.reason) {
        meta.appendChild(h('span.warn', ' — ' + info.reason));
      }
      paint();
    },
  };
}

/**
 * The one clock format for time-lapse instants. Exported because the stage's
 * field legend now prints the same value (main.js) and the two must not be
 * able to disagree about what time the picture is.
 */
export function fmtTime(s) {
  if (!Number.isFinite(s)) return NULL_GLYPH;
  if (s < 90) return `${num(s, 0)} s`;
  if (s < 5400) return `${num(s / 60, 1)} min`;
  return `${num(s / 3600, 2)} h`;
}
