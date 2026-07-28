// wl.js — WINDOW / LEVEL over the field overlay, with the value histogram
// that the window is set against.
//
// This is the control every diagnostic viewer has and the reason it has it:
// the served field is a float grid whose interesting structure is almost
// never spread evenly across [min, max]. A drug-AUC grid is a spike at zero
// with a thin perivascular tail; mapped edge-to-edge it renders as a flat
// wash. Narrowing the window onto the tail is not "enhancement" — the data is
// unchanged, only the transfer function moves, and both ends of the window
// are stated numerically at all times so the reader always knows what the
// colour they are looking at is worth.
//
// Values clamped by the window are counted and reported (`clipped lo / hi`),
// because a window that saturates 40% of the region is a window that is
// hiding something, and the operator has to be told.

import { h } from '../util/dom.js';
import { auto, pct, clamp } from '../util/fmt.js';
import { buildHistogram, quantile, rebin } from '../analysis.js';
import { sample } from '../gl/ramps.js';
import { fieldMeta } from '../fields.js';

const H = 80;            // histogram height in css px
const PAD_B = 17;        // ramp strip (4 px) + the numeric axis row below it
const GRAB = 7;          // px within which a pointer grabs a window edge

function cssv(el, name) {
  return getComputedStyle(el).getPropertyValue(name).trim() || '#888';
}

/** Canvas text in the SAME mono stack as the DOM. Canvas has no
 *  font-variant-numeric, but every face in --font-mono is fixed-advance, so
 *  digits still sit in columns. */
function monoFont(el, px) {
  return `${px}px ${cssv(el, '--font-mono') || 'monospace'}`;
}

export function createWindowLevel({ viewer, onChange }) {
  let field = null;
  let name = null;
  let hist = null;
  let histSrc = null;      // the values buffer `hist` was built from
  let brush = null;        // a single value marked on the axis (hovered cell)
  let roiHist = null;      // {hist, label} overlaid outline
  let display = null;      // rebinned counts at the current pixel width
  let logScale = true;

  const cv = h('canvas.wl-canvas', { 'aria-hidden': 'true', title:
    'Field value histogram. Drag an edge to change the window width, drag inside to move the level, wheel to widen or narrow.' });
  const ctx = cv.getContext('2d');

  const centreInp = h('input.inp.num', {
    type: 'number', step: 'any', id: 'wl-centre',
    'aria-label': 'Window centre in field units',
    oninput: () => commit(),
  });
  const widthInp = h('input.inp.num', {
    type: 'number', step: 'any', min: '0', id: 'wl-width',
    'aria-label': 'Window width in field units',
    oninput: () => commit(),
  });
  const unitC = h('span.unit', '');
  const unitW = h('span.unit', '');

  const loOut = h('span.mono', '—');
  const hiOut = h('span.mono', '—');
  const clipOut = h('span.mono', '—');

  const btnFull = h('button.btn.sm', {
    title: 'Window = the full served range (TMSF min..max). Shift+W.',
    onclick: () => { viewer.resetWindow(); sync(); },
  }, 'full');
  const btnP98 = h('button.btn.sm', {
    title: 'Window = the 2nd to 98th percentile of this field, from a 1024-bin histogram CDF. W.',
    onclick: () => autoWindow(0.02, 0.98),
  }, 'p2–p98');
  const btnP995 = h('button.btn.sm', {
    title: 'Window = the 0.5th to 99.5th percentile — keeps more of the tail than p2–p98.',
    onclick: () => autoWindow(0.005, 0.995),
  }, 'p0.5–p99.5');
  // The label states the STATE, not an ambiguous noun. It read `log count` with
  // aria-pressed=true, from which a reader could not tell whether pressing it
  // turns log on or off — the only toggle in this panel that did not restate
  // its own state.
  const btnLog = h('button.btn.sm', {
    'aria-pressed': 'true',
    title: 'Histogram count axis. Log is the default because a field with a large zero spike hides its own tail on a linear count axis.',
    onclick: () => {
      logScale = !logScale;
      btnLog.setAttribute('aria-pressed', String(logScale));
      btnLog.textContent = logScale ? 'count: log' : 'count: linear';
      draw();
    },
  }, 'count: log');

  const root = h('div.wl',
    h('div.wl-plot', cv),
    h('div.fgrid.wl-grid',
      h('label.lbl', { for: 'wl-centre' }, 'Centre'),
      h('div.with-unit', centreInp, unitC),
      h('label.lbl', { for: 'wl-width' }, 'Width'),
      h('div.with-unit', widthInp, unitW),
    ),
    h('div.btn-row.wl-btns', btnFull, btnP98, btnP995, h('span.grow'), btnLog),
    h('dl.kv.wl-kv',
      h('dt', 'window'),
      h('dd', loOut, h('span.dim', ' → '), hiOut),
      h('dt', { title: 'Fraction of the grid clamped to the low and high ends of the window. A window that saturates a large fraction is hiding structure.' }, 'clipped'),
      h('dd', clipOut, h('span.unit', '%')),
    ),
  );

  /* ------------------------------------------------------------ state --- */

  function setField(fieldName, f) {
    name = fieldName;
    field = f || null;
    hist = null;
    roiHist = null;
    brush = null;
    histSrc = null;
    if (field) {
      hist = buildHistogram(field.values, field.min, field.max);
      histSrc = field.values;
      const meta = fieldMeta(name || '');
      unitC.textContent = meta.unit;
      unitW.textContent = meta.unit;
    } else {
      unitC.textContent = '';
      unitW.textContent = '';
    }
    root.classList.toggle('empty', !field);
    sync();
  }

  function autoWindow(qLo, qHi) {
    if (!hist) return;
    const lo = quantile(hist, qLo), hi = quantile(hist, qHi);
    if (lo === null || hi === null || !(hi > lo)) { viewer.resetWindow(); sync(); return; }
    viewer.setWindow((lo + hi) / 2, hi - lo);
    sync();
  }

  function commit() {
    const c = parseFloat(centreInp.value);
    const w = parseFloat(widthInp.value);
    if (!Number.isFinite(c) || !Number.isFinite(w) || w <= 0) return;
    viewer.setWindow(c, w);
    sync(true);
  }

  /** Repaint everything that depends on the window, and tell the app. */
  function sync(fromInput = false) {
    if (!field) {
      loOut.textContent = '—'; hiOut.textContent = '—'; clipOut.textContent = '—';
      centreInp.value = ''; widthInp.value = '';
      centreInp.disabled = widthInp.disabled = true;
      btnFull.disabled = btnP98.disabled = btnP995.disabled = true;
      draw();
      return;
    }
    centreInp.disabled = widthInp.disabled = false;
    btnFull.disabled = btnP98.disabled = btnP995.disabled = false;
    const [lo, hi] = viewer.windowRange();
    if (!fromInput) {
      centreInp.value = fmtIn((lo + hi) / 2);
      widthInp.value = fmtIn(hi - lo);
    }
    loOut.textContent = auto(lo);
    hiOut.textContent = auto(hi);
    clipOut.textContent = clipText(lo, hi);
    draw();
    if (onChange) onChange();
  }

  // Six significant figures in the editable boxes. More is false precision on
  // a grid whose step is tens of micrometres, and it overflows the field.
  function fmtIn(v) {
    if (!Number.isFinite(v)) return '';
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e6)) return v.toExponential(4);
    return String(Number(v.toPrecision(6)));
  }

  function clipText(lo, hi) {
    if (!hist || !hist.total) return '—';
    const w = (hist.max - hist.min) / hist.n;
    let below = 0, above = 0;
    for (let i = 0; i < hist.n; i++) {
      const v = hist.min + (i + 0.5) * w;
      if (v < lo) below += hist.bins[i];
      else if (v > hi) above += hist.bins[i];
    }
    return `${pct(below / hist.total, 1)} lo · ${pct(above / hist.total, 1)} hi`;
  }

  /* ------------------------------------------------------------- paint -- */

  function draw() {
    const rect = cv.getBoundingClientRect();
    const w = Math.max(40, Math.round(rect.width || root.clientWidth || 300));
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(H * dpr)) {
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(H * dpr);
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, H);

    const bg = cssv(root, '--bg-void');
    const line = cssv(root, '--line-1');
    const fg2 = cssv(root, '--fg-2');
    const fg3 = cssv(root, '--fg-3');
    const fg0 = cssv(root, '--fg-0');
    const accent = cssv(root, '--accent');
    const focus = cssv(root, '--focus');

    ctx.fillStyle = bg;
    ctx.fillRect(0, 0, w, H);

    // The histogram must never outlive the pixels it describes: if the values
    // buffer has been swapped underneath us (the time-lapse path calls
    // setFieldValues, which reuses the texture), rebuild before drawing rather
    // than label a new frame with the previous frame's distribution.
    if (field && field.values !== histSrc) {
      hist = buildHistogram(field.values, field.min, field.max);
      histSrc = field.values;
    }

    const plotH = H - PAD_B;
    if (!hist || !hist.total) {
      ctx.fillStyle = fg3;
      ctx.font = monoFont(root, 9.5);
      ctx.fillText('no field overlay loaded', 6, plotH / 2 + 3);
      ctx.strokeStyle = line;
      ctx.strokeRect(0.5, 0.5, w - 1, H - 1);
      return;
    }

    const cols = Math.max(24, Math.min(w, 512));
    display = rebin(hist, cols);
    let peak = 0;
    for (let i = 0; i < cols; i++) if (display[i] > peak) peak = display[i];
    const scale = (v) => (logScale ? Math.log1p(v) / Math.log1p(peak || 1) : v / (peak || 1));

    const [lo, hi] = viewer.windowRange();
    const x = (v) => ((v - hist.min) / (hist.max - hist.min)) * w;
    const cw = w / cols;

    // bars: dim outside the window, bright inside
    for (let i = 0; i < cols; i++) {
      const v0 = hist.min + (i / cols) * (hist.max - hist.min);
      const v1 = hist.min + ((i + 1) / cols) * (hist.max - hist.min);
      const inside = v1 >= lo && v0 <= hi;
      const bh = Math.max(display[i] > 0 ? 1 : 0, scale(display[i]) * (plotH - 3));
      ctx.fillStyle = inside ? fg2 : fg3;
      ctx.globalAlpha = inside ? 0.95 : 0.38;
      ctx.fillRect(i * cw, plotH - bh, Math.max(0.7, cw - 0.25), bh);
    }
    ctx.globalAlpha = 1;

    // ROI histogram, as an outline over the same axis
    if (roiHist && roiHist.hist && roiHist.hist.total) {
      const r = rebin(roiHist.hist, cols);
      let rp = 0;
      for (let i = 0; i < cols; i++) if (r[i] > rp) rp = r[i];
      // The ROI curve is normalised to its OWN peak: it is a shape comparison
      // against the whole field, not a count comparison. Said in the caption.
      ctx.beginPath();
      for (let i = 0; i < cols; i++) {
        const vv = logScale ? Math.log1p(r[i]) / Math.log1p(rp || 1) : r[i] / (rp || 1);
        const bh = vv * (plotH - 3);
        const px = i * cw + cw / 2;
        if (i === 0) ctx.moveTo(px, plotH - bh); else ctx.lineTo(px, plotH - bh);
      }
      ctx.strokeStyle = accent;
      ctx.lineWidth = 1;
      ctx.stroke();
    }

    // window shading + edges
    const xl = clamp(x(lo), 0, w), xh = clamp(x(hi), 0, w);
    ctx.strokeStyle = accent;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(Math.round(xl) + 0.5, 0); ctx.lineTo(Math.round(xl) + 0.5, plotH);
    ctx.moveTo(Math.round(xh) - 0.5, 0); ctx.lineTo(Math.round(xh) - 0.5, plotH);
    ctx.stroke();
    // grab handles, so the affordance is visible and not just hoverable
    ctx.fillStyle = accent;
    ctx.fillRect(Math.round(xl) - 1.5, 0, 4, 5);
    ctx.fillRect(Math.round(xh) - 2.5, 0, 4, 5);

    // the ramp, painted only across the window: this strip IS the legend
    const stripY = plotH + 1;
    const rampName = viewer.rampName || 'viridis';
    const stripW = Math.max(1, xh - xl);
    for (let px = 0; px < stripW; px++) {
      ctx.fillStyle = sample(rampName, px / Math.max(1, stripW - 1));
      ctx.fillRect(xl + px, stripY, 1, 4);
    }
    ctx.fillStyle = fg3;
    if (xl > 0) ctx.fillRect(0, stripY, xl, 4);
    if (xh < w) ctx.fillRect(xh, stripY, w - xh, 4);

    // brush: the value of the cell under the cursor / pin
    if (brush !== null && Number.isFinite(brush)) {
      const bx = clamp(x(brush), 0, w);
      ctx.strokeStyle = focus;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(Math.round(bx) + 0.5, 0);
      ctx.lineTo(Math.round(bx) + 0.5, plotH + 5);
      ctx.stroke();
      ctx.fillStyle = focus;
      ctx.beginPath();
      ctx.moveTo(bx, 0); ctx.lineTo(bx - 3, 5); ctx.lineTo(bx + 3, 5);
      ctx.closePath(); ctx.fill();
    }

    // Axis row, BELOW the ramp strip: the data ends in muted type, the window
    // centre in full contrast, because the centre is the number the operator
    // is actually setting.
    ctx.font = monoFont(root, 9);
    ctx.fillStyle = fg2;
    ctx.textBaseline = 'bottom';
    ctx.fillText(auto(hist.min), 1, H - 1);
    const hiTxt = auto(hist.max);
    ctx.fillText(hiTxt, w - ctx.measureText(hiTxt).width - 1, H - 1);
    ctx.fillStyle = fg0;
    const mid = auto((lo + hi) / 2);
    const midW = ctx.measureText(mid).width;
    const mx = clamp((xl + xh) / 2 - midW / 2,
      ctx.measureText(auto(hist.min)).width + 6,
      w - ctx.measureText(hiTxt).width - midW - 6);
    if (xh - xl > midW + 12 && mx > 0) ctx.fillText(mid, mx, H - 1);

    ctx.strokeStyle = line;
    ctx.strokeRect(0.5, 0.5, w - 1, H - 1);
  }

  /* -------------------------------------------------------- interaction - */

  let drag = null;   // {mode:'lo'|'hi'|'move', x0, lo0, hi0}

  cv.addEventListener('pointerdown', (e) => {
    if (!hist) return;
    const r = cv.getBoundingClientRect();
    const px = e.clientX - r.left;
    const [lo, hi] = viewer.windowRange();
    const toX = (v) => ((v - hist.min) / (hist.max - hist.min)) * r.width;
    const xl = toX(lo), xh = toX(hi);
    let mode = 'move';
    if (Math.abs(px - xl) <= GRAB) mode = 'lo';
    else if (Math.abs(px - xh) <= GRAB) mode = 'hi';
    else if (px < xl || px > xh) mode = 'move';
    drag = { mode, x0: px, lo0: lo, hi0: hi, w: r.width };
    cv.setPointerCapture(e.pointerId);
    e.preventDefault();
  });

  cv.addEventListener('pointermove', (e) => {
    const r = cv.getBoundingClientRect();
    if (!drag) {
      if (!hist) return;
      const [lo, hi] = viewer.windowRange();
      const toX = (v) => ((v - hist.min) / (hist.max - hist.min)) * r.width;
      const px = e.clientX - r.left;
      const near = Math.abs(px - toX(lo)) <= GRAB || Math.abs(px - toX(hi)) <= GRAB;
      cv.style.cursor = near ? 'ew-resize' : 'grab';
      return;
    }
    const px = e.clientX - r.left;
    const span = hist.max - hist.min;
    const dv = ((px - drag.x0) / drag.w) * span;
    if (drag.mode === 'lo') {
      const lo = Math.min(drag.lo0 + dv, drag.hi0 - span * 1e-4);
      viewer.setWindow((lo + drag.hi0) / 2, drag.hi0 - lo);
    } else if (drag.mode === 'hi') {
      const hi = Math.max(drag.hi0 + dv, drag.lo0 + span * 1e-4);
      viewer.setWindow((drag.lo0 + hi) / 2, hi - drag.lo0);
    } else {
      viewer.setWindow((drag.lo0 + drag.hi0) / 2 + dv, drag.hi0 - drag.lo0);
    }
    sync();
  });

  const endDrag = (e) => {
    if (!drag) return;
    drag = null;
    try { cv.releasePointerCapture(e.pointerId); } catch { /* already gone */ }
  };
  cv.addEventListener('pointerup', endDrag);
  cv.addEventListener('pointercancel', endDrag);

  cv.addEventListener('wheel', (e) => {
    if (!hist) return;
    e.preventDefault();
    const [lo, hi] = viewer.windowRange();
    const f = Math.exp(e.deltaY * 0.0012);
    const c = (lo + hi) / 2;
    viewer.setWindow(c, Math.max((hist.max - hist.min) * 1e-4, (hi - lo) * f));
    sync();
  }, { passive: false });

  const ro = new ResizeObserver(() => draw());
  ro.observe(cv);

  draw();

  return {
    el: root,
    setField,
    sync,
    draw,
    /** Mark one value on the axis — the hovered or pinned cell. */
    setBrush(v) {
      const next = (v === null || v === undefined || !Number.isFinite(v)) ? null : v;
      if (next === brush) return;
      brush = next;
      draw();
    },
    /** Overlay an ROI's own value distribution as an outline. */
    setRoiHist(hOrNull) { roiHist = hOrNull; draw(); },
    get hist() { return hist; },
    autoWindow,
  };
}
