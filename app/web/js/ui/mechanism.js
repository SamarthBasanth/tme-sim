// ui/mechanism.js — the two stage overlays this feature adds:
//
//   1. the MECHANISM LEGEND: every attribution category with its count and its
//      share, the delivery / dose / resistance / chance roll-up as one stacked
//      bar, and the bit-exact reconciliation against the engine's own
//      kill_fraction. Rows toggle their category in the viewer.
//   2. the TIME-LAPSE SCRUBBER: play / pause / step / seek over the captured
//      drug-concentration frames.
//
// Both are new files and both are additive: they are appended to #canvas-wrap
// alongside the existing overlays and never touch them.
//
// Design rule being honoured (API_CONTRACT.md): every category carries a
// visible label AND a number. The categorical palette clears its CVD target on
// the dark surround but three slots sit under 3:1 on the light one, so
// colour-alone identity is not available and is not used.

import { h, fill } from '../util/dom.js';
import { num, int, pct, NULL_GLYPH } from '../util/fmt.js';
import { MECH, MECH_META, LIMITER_LABEL, dominantReason } from '../mechanism.js';

/* ====================================================== legend ========= */

export function createMechanismLegend(viewer, { onToggle } = {}) {
  const title = h('div.m-title');
  const rows = h('div.m-rows');
  const bar = h('div.m-bar');
  const note = h('div.m-note');
  const card = h('div.ov.ov-mech.ov-card.mechleg.hidden.interactive',
    { role: 'figure', 'aria-label': 'Mechanism attribution legend' },
    title, rows, bar, note);

  let agg = null;

  function render() {
    if (!agg) { card.classList.add('hidden'); return; }
    card.classList.remove('hidden');

    const counts = agg.counts || {};
    const nTumor = agg.n_tumor || 0;
    const killed = counts.KILLED || 0;

    fill(title,
      h('span', 'mechanism'),
      h('span.grow'),
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
    setHidden(on) { card.classList.toggle('hidden', !!on || !agg); },
  };
}

/* ==================================================== scrubber ========= */

/** Frames per second when the time-lapse is playing. */
const PLAY_FPS = 8;

export function createTimelapseScrubber({ onFrame }) {
  const btnPlay = h('button.btn.ghost.sm.s-btn', { title: 'Play / pause (space)' }, '▶');
  const btnStop = h('button.btn.ghost.sm.s-btn', { title: 'Back to first frame' }, '⏮');
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
    read.appendChild(h('span', `t ${fmtTime(f.time_s)}`));
    read.appendChild(h('span.unit', ` · frame ${i + 1}/${series.n_frames}`));
    read.appendChild(h('span.unit', ` · max ${num(f.max, 3)} uM`));
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
    timer = setInterval(() => {
      if (i >= series.n_frames - 1) { pause(); return; }
      go(i + 1);
    }, 1000 / PLAY_FPS);
  }

  function pause() {
    if (timer) clearInterval(timer);
    timer = 0;
    btnPlay.textContent = '▶';
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
      const bits = [
        `${series.n_frames} frames`,
        `every ${fmtTime(info.realised_interval_s || series.interval_s)}`,
        `over ${fmtTime(series.window_s)}`,
        `${(series.bytes / 1e6).toFixed(1)} MB`,
        `drug concentration (uM), ${info.grid ? info.grid.join('x') : '?'} grid`,
      ];
      fill(meta, h('span', bits.join(' · ')));
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

function fmtTime(s) {
  if (!Number.isFinite(s)) return NULL_GLYPH;
  if (s < 90) return `${num(s, 0)} s`;
  if (s < 5400) return `${num(s / 60, 1)} min`;
  return `${num(s / 3600, 2)} h`;
}
