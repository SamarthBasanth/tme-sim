// overlays.js — the things drawn over the canvas: scale bar, field legend,
// hover read-out, and the render statistics.
//
// The scale bar is recomputed from the live camera on every frame, so it is
// correct at any zoom — it is not a static graphic with a caption.

import { h, fill } from '../util/dom.js';
import { nice125, num, int, auto, pct, NULL_GLYPH } from '../util/fmt.js';
import { gradientCss } from '../gl/ramps.js';
import { fieldMeta } from '../fields.js';
import { CELL_TYPES } from '../gl/viewer.js';
import { chip, fieldProvenance } from './provenance.js';

const CHAN_VARS = ['--ch-tumor', '--ch-stromal', '--ch-immune', '--ch-endothelial', '--ch-necrotic'];

export function createOverlays(viewer) {
  /* ------------------------------------------------------- scale bar --- */
  const barEl = h('div.bar');
  const barTxt = h('div.txt', '—');
  const scalebar = h('div.ov.ov-scalebar.ov-card.scalebar', { 'aria-hidden': 'true' }, barEl, barTxt);

  function updateScalebar() {
    const zoom = viewer.cam.zoom;            // css px per um
    if (!(zoom > 0)) return;
    const targetPx = 110;
    let lenUm = nice125(targetPx / zoom);
    let px = lenUm * zoom;
    // keep the bar inside a sane pixel band
    let guard = 0;
    while (px < 46 && guard++ < 8) { lenUm = nice125(lenUm * 2.2); px = lenUm * zoom; }
    while (px > 190 && guard++ < 16) { lenUm = nice125(lenUm / 2.2); px = lenUm * zoom; }
    barEl.style.width = px.toFixed(1) + 'px';
    const len = lenUm >= 1000
      ? `${num(lenUm / 1000, lenUm % 1000 === 0 ? 0 : 2)} mm`
      : `${num(lenUm, lenUm < 1 ? 2 : 0)} um`;
    // The scale bar is exact; the DOTS are not, once a cell is smaller than a
    // pixel. Say so on the same element rather than letting the bar imply that
    // everything in the frame is to scale. See MIN_POINT_CSS_PX in gl/viewer.js.
    const over = viewer.pointOversize;
    barTxt.textContent = over > 1.1 ? `${len} · cells ${num(over, 1)}× scale` : len;
    scalebar.title = over > 1.1
      ? `Scale bar is exact. Cells are ${num(over, 1)}× life size: a ${num(viewer.cellDiamUm, 0)} um `
        + `cell is ${num(viewer.trueDiamCssPx, 2)} px at this zoom, below the 1.75 px legibility `
        + 'floor. Apparent density is inflated by roughly the square of that. Zoom in to see them to scale.'
      : 'Scale bar is exact. Cells are drawn to scale at this zoom.';
  }

  /* ---------------------------------------------------------- legend --- */
  const legTitle = h('div.l-title');
  const legRamp = h('div.l-ramp');
  const legTicks = h('div.l-ticks');
  const legNote = h('div.l-note');
  const legend = h('div.ov.ov-legend.ov-card.legend.hidden', { role: 'figure', 'aria-label': 'Field legend' },
    legTitle, legRamp, legTicks, legNote);

  let legendField = null, legendName = null;

  // The legend labels the WINDOW, not the payload range, because the window is
  // what the colours on screen actually mean. When the two differ the legend
  // says so and prints the full served range underneath — a narrowed window
  // that is not declared is a lie about the image.
  function updateLegend(name, field) {
    legendName = name || null;
    legendField = field || null;
    if (!legendName || !legendField) { legend.classList.add('hidden'); return; }
    const meta = fieldMeta(legendName);
    legend.classList.remove('hidden');
    fill(legTitle, h('span', meta.label), h('span.unit', meta.unit));
    legRamp.style.backgroundImage = gradientCss(viewer.rampName, 24, viewer.fieldAlphaRamp);

    const [lo, hi] = viewer.windowRange();
    fill(legTicks,
      h('span', auto(lo)),
      h('span', auto(lo + (hi - lo) * 0.5)),
      h('span', auto(hi)));

    const windowed = Math.abs(lo - legendField.min) > 1e-9 * Math.abs(legendField.max - legendField.min || 1)
      || Math.abs(hi - legendField.max) > 1e-9 * Math.abs(legendField.max - legendField.min || 1);
    // Run-aware, not the static table: these pixels ARE the field, so if the
    // run's own metrics escalate to ASSUMED the ramp beside them must too.
    const fp = fieldProvenance(legendName);
    fill(legNote,
      h('span', `${legendField.rows}×${legendField.cols} @ ${num(legendField.dx, 1)} um`),
      ' ',
      chip(fp.kind, fp.why),
      // The overlay's opacity carries the value as well as its hue, so the
      // legend has to say that the bottom of the ramp is nothing at all —
      // otherwise a transparent region reads as "no data" when it is a
      // measured zero. And when the window floor is above the data floor,
      // "invisible" additionally covers everything clamped below it.
      viewer.fieldAlphaRamp
        ? h(windowed ? 'div.l-win' : 'div.l-enc', {
          title: 'The overlay fades to fully transparent at the bottom of the window: '
            + 'opacity is linear in the value, with no contrast stretch. Bare background '
            + 'therefore means a value at or below the window floor, not missing data.',
        }, windowed ? '0 opacity at window floor' : 'transparent = 0')
        : null,
      windowed
        ? h('div.l-win', { title: 'The colour ramp is mapped across the window, not the full data range. Values outside it are clamped — and, with the opacity ramp on, drawn as nothing.' },
          'windowed · data ', auto(legendField.min), ' → ', auto(legendField.max))
        : null);
  }

  /** Re-label the legend after the window moved, without refetching anything. */
  function refreshWindow() { if (legendName) updateLegend(legendName, legendField); }

  /* --------------------------------------------------------- read-out -- */
  const roBody = h('div');
  const readout = h('div.ov.ov-readout.ov-card.readout', { role: 'status', 'aria-live': 'off' }, roBody);

  function updateReadout(info, fieldName, pinned) {
    if (!info) {
      fill(roBody, h('div.r-empty', 'hover a cell'));
      readout.style.borderColor = '';
      return;
    }
    const meta = fieldName ? fieldMeta(fieldName) : null;
    const comp = info.composition;
    const total = comp.total;
    const rows = [];
    for (let i = 0; i < CELL_TYPES.length; i++) {
      const c = comp.counts[i];
      if (!c) continue;
      rows.push(h('div', {
        style: {
          display: 'grid', gridTemplateColumns: '8px 1fr auto auto', gap: '4px',
          alignItems: 'center', fontSize: 'var(--fs-micro)', lineHeight: '1.45',
        },
      },
        h('span', { style: { width: '7px', height: '7px', background: `var(${CHAN_VARS[i]})` } }),
        h('span', { style: { color: 'var(--fg-2)' } }, CELL_TYPES[i]),
        h('span.mono', int(c)),
        h('span.mono', { style: { color: 'var(--fg-2)', minWidth: '34px', textAlign: 'right' } },
          total > 0 ? pct(c / total, 0) + '%' : NULL_GLYPH),
      ));
    }

    fill(roBody,
      h('div.r-head',
        h('span.r-swatch', { style: { background: `var(${CHAN_VARS[info.type] || '--fg-2'})` } }),
        h('b', info.typeName),
        h('span.grow'),
        h('span.mono.dim', '#' + info.index),
        pinned ? h('span.chip.st-running', 'pinned') : null),
      h('dl.kv',
        h('dt', 'x'), h('dd', num(info.x, 1), h('span.unit', 'um')),
        h('dt', 'y'), h('dd', num(info.y, 1), h('span.unit', 'um')),
        meta ? h('dt', meta.short) : null,
        meta ? h('dd', info.fieldValue === null || info.fieldValue === undefined
          ? h('span.dim', NULL_GLYPH)
          : h('span', auto(info.fieldValue)), h('span.unit', meta.unit)) : null,
      ),
      h('div.hr-lbl', { style: { margin: '4px 0 1px' } },
        `local composition r=${num(info.compositionRadiusUm, 0)} um`),
      total < 0
        ? h('div.r-empty', 'radius spans too much of the section')
        : h('div', rows.length ? rows : h('div.r-empty', 'no neighbours in radius')),
      h('div', { style: { fontSize: 'var(--fs-micro)', color: 'var(--fg-3)', marginTop: '2px' } },
        total >= 0 ? `${int(total)} cells in radius` : ''),
    );
    readout.style.borderColor = pinned ? 'var(--accent)' : '';
  }

  /* ------------------------------------------------------ empty state -- */
  // NOT `.interactive`: the notice spans the full width of the top strip at
  // z-index 3 and contains nothing clickable, so pointer-events:auto only
  // meant it swallowed drags that started in that band.
  const empty = h('div.ov.ov-empty',
    h('div.notice',
      h('h3', 'No geometry loaded'),
      h('p', 'Submit a RunSpec, or open a completed run from the history dock, to load its cell geometry over ',
        h('code.badge-code', 'GET /api/runs/{id}/geometry'), '.'),
      h('p.notice-keys',
        h('kbd', 'wheel'), ' zoom · ', h('kbd', 'drag'), ' pan · ',
        h('kbd', 'click'), ' pin a cell · ', h('kbd', 'F'), ' fit · ',
        h('kbd', '?'), ' all keys'),
    ));

  return {
    nodes: [empty, scalebar, legend, readout],
    updateScalebar, updateLegend, updateReadout, refreshWindow,
    setEmpty(on) { empty.classList.toggle('hidden', !on); },
    refreshRamp() {
      legRamp.style.backgroundImage = gradientCss(viewer.rampName, 24, viewer.fieldAlphaRamp);
    },
  };
}
