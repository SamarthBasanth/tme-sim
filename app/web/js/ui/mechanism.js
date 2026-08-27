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
  const rows = h('div.m-rows');
  const bar = h('div.m-bar');
  const note = h('div.m-note');
  const body = h('div.m-body', rows, bar, note);

  const card = h('div.mechleg.mechleg-docked',
    { role: 'figure', 'aria-label': 'Mechanism attribution breakdown' },
    title, body);

  let agg = null;

  function render() {
    if (!agg) {
      // An empty state, not a disappearance. A reconstruction-only run has
      // nothing to attribute and the section says so where the numbers would
      // be — a section that silently vanishes reads as a bug.
      fill(title);
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

    headline.textContent = `${pct(killed / (nTumor || 1), 1)}% killed`;

    // No 'mechanism' word here any more: the section's own <summary> says it
    // one line above, and a card that repeats its container's title is how a
    // rail ends up with two headings per section.
    // SHOW ALL, and it only exists when something is hidden.
    //
    // Rows here toggle one category at a time, which was fine while toggling
    // was the only way to hide anything -- you turned one off, you turned it
    // back on. Then THE FINDING gained a click that ISOLATES a category
    // (finding.js -> main.js onIsolateMechanism), which turns six rows off at
    // once. Getting back meant clicking all six. A one-click action that takes
    // six clicks to undo is a trap, and the undo belongs on the legend that
    // owns the visibility state rather than on the panel that triggered it.
    //
    // Rendered only when at least one category is hidden, so it is absent in
    // the ordinary case and cannot read as a control the reader has to think
    // about before they have hidden anything.
    const anyHidden = Array.prototype.some.call(
      viewer.labelVis, (v, i) => i < MECH_META.length && !v);
    const showAll = h('button.btn.ghost.tiny', {
      type: 'button',
      title: 'Show every category again',
      onclick: () => {
        for (let i = 0; i < viewer.labelVis.length; i += 1) {
          viewer.setLabelVisible(i, true);
        }
        render();
        if (onToggle) onToggle(null, true);
      },
    }, 'show all');

    fill(title,
      headline,
      h('span.grow'),
      anyHidden ? showAll : null,
      h('span.mono', { style: { color: 'var(--fg-2)' } },
        `${int(nTumor)} tumour cells`));

    /* --- one row per category ---------------------------------------- */
    const order = [
      MECH.KILLED,
      MECH.SURVIVED_UNREACHED, MECH.SURVIVED_UNDERDOSED,
      MECH.SURVIVED_RESISTANT, MECH.SURVIVED_STOCHASTIC,
      MECH.UNATTRIBUTED, MECH.ALREADY_DEAD, MECH.NOT_TUMOR,
    ];
    const kids = [];
    for (const code of order) {
      const meta = MECH_META[code];
      const n = counts[meta.name] || 0;
      // The denominator differs by category and must not be fudged: the
      // outcome categories are shares of the tumour population, the survival
      // categories are shares of the survivors, and NOT_TUMOR is neither.
      const denom = meta.limiter ? (agg.n_survived || 0) : (code === MECH.NOT_TUMOR ? 0 : nTumor);
      if (code === MECH.SURVIVED_UNREACHED) kids.push(h('div.m-sep'));
      const row = h('div.m-row', {
        role: 'button',
        tabindex: '0',
        title: (agg.label_descriptions || {})[meta.name] || meta.name,
      },
        h('span.m-sw', { style: { background: `var(${meta.var})` } }),
        h('span.m-name', {
          style: meta.alarm && n > 0 ? { color: 'var(--crit)' } : null,
        }, meta.short),
        h('span.m-n', int(n)),
        h('span.m-p', denom > 0 ? pct(n / denom, 1) + '%' : NULL_GLYPH),
      );
      if (!n) row.classList.add('zero');
      if (!viewer.labelVis[code]) row.classList.add('off');
      const flip = () => {
        const on = !viewer.labelVis[code];
        viewer.setLabelVisible(code, on);
        row.classList.toggle('off', !on);
        if (onToggle) onToggle(code, on);
      };
      row.addEventListener('click', flip);
      row.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); flip(); }
      });
      kids.push(row);
    }
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
    kids2.push(h('div',
      h('b', pct(killed / (nTumor || 1), 1) + '%'), ' of tumour cells killed. ',
      dom
        ? h('span', 'Of the survivors, the largest share is ',
          h('b', { style: { color: `var(${MECH_META.find((m) => m.limiter === dom.limiter).var})` } },
            dom.label), '-limited (', pct(dom.fraction, 1), '%).')
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
