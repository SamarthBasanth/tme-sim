// tools.js — the measurement tools: LINE PROFILE, REGION (ROI) STATISTICS and
// POINT-TO-POINT DISTANCE, plus the 2D annotation layer they are drawn on.
//
// The annotation layer is a plain 2D canvas over the WebGL canvas, not more
// GL geometry. Two reasons: the cell pass is the frame budget and nothing
// should be added to it, and an annotation redraw only happens when the camera
// or the annotation changes — a few dozen strokes, off the hot path entirely.
//
// What each tool is FOR, on this data:
//   profile  — the perivascular gradient. Drag a line out from a vessel and
//              the drug-AUC decay along it is the transport model's central
//              claim, read directly instead of inferred from a colour.
//   region   — composition and exposure for one piece of tissue. The
//              section-level kill fraction is an average over a section that
//              is not homogeneous; a region says what it is an average OF.
//   distance — because a length on screen means nothing without um.

import { h, clear, fill } from '../util/dom.js';
import { auto, num, int, pct, um, NULL_GLYPH, clamp } from '../util/fmt.js';
import { CELL_TYPES } from '../gl/viewer.js';
import { sampleLine, rectFieldStats, rectCellStats, buildHistogram } from '../analysis.js';
import { fieldMeta } from '../fields.js';
import { chip } from './provenance.js';

export const TOOLS = ['pick', 'profile', 'region', 'distance'];
const CHAN_VARS = ['--ch-tumor', '--ch-stromal', '--ch-immune', '--ch-endothelial', '--ch-necrotic'];
const MAX_ROI = 4;
const PLOT_H = 104;

const cssv = (el, n) => getComputedStyle(el).getPropertyValue(n).trim() || '#888';
/** Canvas text in the same mono stack the DOM uses. */
const monoFont = (el, px) => `${px}px ${cssv(el, '--font-mono') || 'monospace'}`;

export function createTools({ viewer, getFieldName, onChange, onRoiChange }) {
  let tool = 'pick';
  let profile = null;         // {x0,y0,x1,y1}
  let profileData = null;
  let profileMark = null;     // fraction along the line, from the plot hover
  let distance = null;        // {x0,y0,x1,y1}
  let rois = [];
  let roiSeq = 0;
  let pending = null;         // geometry being dragged out right now

  /* ================================================== annotation layer == */

  const annot = h('canvas.ov.ov-annot', { 'aria-hidden': 'true' });
  const actx = annot.getContext('2d');

  let hadContent = false;

  function drawAnnot() {
    const wrap = annot.parentElement;
    if (!wrap) return;
    // Free when there is nothing annotated: no resize check, no clear, no
    // upload. This matters because the layer is repainted once per presented
    // frame while the camera is moving, and the common case is an empty layer.
    const has = !!(profile || distance || pending || rois.length);
    if (!has && !hadContent) return;
    hadContent = has;
    const w = wrap.clientWidth, hgt = wrap.clientHeight;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    if (annot.width !== Math.round(w * dpr) || annot.height !== Math.round(hgt * dpr)) {
      annot.width = Math.round(w * dpr);
      annot.height = Math.round(hgt * dpr);
      annot.style.width = w + 'px';
      annot.style.height = hgt + 'px';
    }
    actx.setTransform(dpr, 0, 0, dpr, 0, 0);
    actx.clearRect(0, 0, w, hgt);
    if (!viewer.geom) return;

    const accent = cssv(annot, '--accent');
    const focus = cssv(annot, '--focus');
    const bg = cssv(annot, '--bg-1');
    const S = (x, y) => viewer.worldToScreen(x, y);

    actx.font = monoFont(annot, 10);
    actx.textBaseline = 'middle';

    const tag = (x, y, text, colour) => {
      const tw = actx.measureText(text).width;
      actx.fillStyle = bg;
      actx.globalAlpha = 0.88;
      actx.fillRect(x, y - 7, tw + 6, 14);
      actx.globalAlpha = 1;
      actx.strokeStyle = colour;
      actx.lineWidth = 1;
      actx.strokeRect(x + 0.5, y - 6.5, tw + 5, 13);
      actx.fillStyle = colour;
      actx.fillText(text, x + 3, y + 0.5);
    };

    /* --- regions ---------------------------------------------------- */
    for (const r of rois) {
      const a = S(r.x0, r.y0), b = S(r.x1, r.y1);
      actx.strokeStyle = accent;
      actx.lineWidth = 1;
      actx.strokeRect(Math.round(Math.min(a.x, b.x)) + 0.5, Math.round(Math.min(a.y, b.y)) + 0.5,
        Math.abs(b.x - a.x), Math.abs(b.y - a.y));
      tag(Math.min(a.x, b.x), Math.min(a.y, b.y) - 9, r.label, accent);
    }

    /* --- line profile ------------------------------------------------ */
    if (profile) {
      const a = S(profile.x0, profile.y0), b = S(profile.x1, profile.y1);
      actx.strokeStyle = accent;
      actx.lineWidth = 1.25;
      actx.beginPath(); actx.moveTo(a.x, a.y); actx.lineTo(b.x, b.y); actx.stroke();
      capMark(a, b, accent);
      capMark(b, a, accent);
      if (profileData) {
        tag((a.x + b.x) / 2 + 6, (a.y + b.y) / 2 - 10,
          `${um(profileData.lengthUm).v} ${um(profileData.lengthUm).u}`, accent);
      }
      if (profileMark !== null) {
        const px = a.x + (b.x - a.x) * profileMark, py = a.y + (b.y - a.y) * profileMark;
        actx.fillStyle = focus;
        actx.beginPath(); actx.arc(px, py, 3.5, 0, Math.PI * 2); actx.fill();
      }
    }

    /* --- distance ----------------------------------------------------- */
    if (distance) {
      const a = S(distance.x0, distance.y0), b = S(distance.x1, distance.y1);
      actx.strokeStyle = focus;
      actx.lineWidth = 1.25;
      actx.setLineDash([4, 3]);
      actx.beginPath(); actx.moveTo(a.x, a.y); actx.lineTo(b.x, b.y); actx.stroke();
      actx.setLineDash([]);
      capMark(a, b, focus);
      capMark(b, a, focus);
      const L = Math.hypot(distance.x1 - distance.x0, distance.y1 - distance.y0);
      const u = um(L);
      tag((a.x + b.x) / 2 + 6, (a.y + b.y) / 2 - 10, `${u.v} ${u.u}`, focus);
    }

    /* --- in-progress -------------------------------------------------- */
    if (pending) {
      const a = S(pending.x0, pending.y0), b = S(pending.x1, pending.y1);
      actx.strokeStyle = accent;
      actx.lineWidth = 1;
      actx.setLineDash([3, 3]);
      if (pending.kind === 'region') {
        actx.strokeRect(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.abs(b.x - a.x), Math.abs(b.y - a.y));
        const wu = Math.abs(pending.x1 - pending.x0), hu = Math.abs(pending.y1 - pending.y0);
        actx.setLineDash([]);
        tag(Math.max(a.x, b.x) + 5, Math.max(a.y, b.y),
          `${num(wu, 0)} × ${num(hu, 0)} um`, accent);
      } else {
        actx.beginPath(); actx.moveTo(a.x, a.y); actx.lineTo(b.x, b.y); actx.stroke();
        actx.setLineDash([]);
        const L = Math.hypot(pending.x1 - pending.x0, pending.y1 - pending.y0);
        const u = um(L);
        tag((a.x + b.x) / 2 + 6, (a.y + b.y) / 2 - 10, `${u.v} ${u.u}`, accent);
      }
      actx.setLineDash([]);
    }

    function capMark(p, other, colour) {
      const dx = other.x - p.x, dy = other.y - p.y;
      const L = Math.hypot(dx, dy) || 1;
      const nx = -dy / L * 5, ny = dx / L * 5;
      actx.strokeStyle = colour;
      actx.lineWidth = 1;
      actx.beginPath();
      actx.moveTo(p.x - nx, p.y - ny);
      actx.lineTo(p.x + nx, p.y + ny);
      actx.stroke();
    }
  }

  /* ====================================================== profile panel = */

  const plot = h('canvas.plot-canvas', { 'aria-hidden': 'true' });
  const pctx = plot.getContext('2d');
  const profStats = h('dl.kv');
  const profHint = h('div.tool-hint');
  const profClear = h('button.btn.ghost.sm', {
    onclick: () => { profile = null; profileData = null; profileMark = null; paintProfile(); drawAnnot(); },
    title: 'Discard the profile line',
  }, 'clear');
  const profileEl = h('div.tool-panel',
    h('div.frow.tool-bar',
      h('span.grow'),
      profClear),
    h('div.plot-wrap', plot),
    profStats,
    profHint,
  );

  function computeProfile() {
    if (!profile || !viewer.field) { profileData = null; return; }
    const L = Math.hypot(profile.x1 - profile.x0, profile.y1 - profile.y0);
    // One sample per screen pixel, bounded — more than that is inventing
    // resolution the grid does not have, less is throwing it away.
    const n = clamp(Math.round(L * viewer.cam.zoom), 32, 512);
    profileData = sampleLine(viewer, profile.x0, profile.y0, profile.x1, profile.y1, n);
  }

  function paintProfile() {
    const has = !!(profile && profileData && profileData.valid);
    profileEl.classList.toggle('empty', !profile);
    profClear.disabled = !profile;
    const fname = getFieldName();

    if (!profile) {
      profHint.textContent = 'No profile. Press L, then drag a line across the tissue.';
      clear(profStats);
      paintPlot();
      return;
    }
    if (!fname) {
      profHint.textContent = 'A profile needs a field overlay — pick one under Field overlay.';
      clear(profStats);
      paintPlot();
      return;
    }
    const meta = fieldMeta(fname);
    const d = profileData;
    const L = um(d ? d.lengthUm : 0);
    fill(profStats,
      h('dt', 'length'), h('dd', L.v, h('span.unit', L.u)),
      h('dt', 'samples'), h('dd', int(d ? d.n : 0), h('span.unit', 'n'),
        h('span.dim.step', `@ ${num(d ? d.stepUm : 0, 2)} um`)),
      h('dt', 'min'), h('dd', has ? auto(d.min) : NULL_GLYPH, h('span.unit', meta.unit)),
      h('dt', 'max'), h('dd', has ? auto(d.max) : NULL_GLYPH, h('span.unit', meta.unit)),
      h('dt', 'mean'), h('dd', has ? auto(d.mean) : NULL_GLYPH, h('span.unit', meta.unit)),
      h('dt', 'source'), h('dd', chip(meta.provenance,
        `${meta.source}. Sampled bilinearly from the served ${fname} grid along the drawn line; no smoothing, no fitting.`)),
    );
    profHint.textContent = d && d.outside
      ? `${d.outside} of ${d.n} samples fall outside the grid and are not plotted.`
      : 'Hover the plot to mark that position on the tissue.';
    paintPlot();
  }

  function paintPlot() {
    const rect = plot.getBoundingClientRect();
    const w = Math.max(60, Math.round(rect.width || 300));
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    if (plot.width !== Math.round(w * dpr) || plot.height !== Math.round(PLOT_H * dpr)) {
      plot.width = Math.round(w * dpr);
      plot.height = Math.round(PLOT_H * dpr);
    }
    pctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const bg = cssv(plot, '--bg-void');
    const line = cssv(plot, '--line-1');
    const fg0 = cssv(plot, '--fg-0');
    const fg3 = cssv(plot, '--fg-3');
    const accent = cssv(plot, '--accent');
    const focus = cssv(plot, '--focus');
    pctx.fillStyle = bg;
    pctx.fillRect(0, 0, w, PLOT_H);
    pctx.strokeStyle = line;
    pctx.strokeRect(0.5, 0.5, w - 1, PLOT_H - 1);
    pctx.font = monoFont(plot, 9);

    const d = profileData;
    if (!d || !d.valid) {
      pctx.fillStyle = fg3;
      pctx.textBaseline = 'middle';
      pctx.fillText(profile ? 'no field value along this line' : 'no profile drawn', 6, PLOT_H / 2);
      return;
    }

    const padL = 4, padR = 4, padT = 6, padB = 12;
    const iw = w - padL - padR, ih = PLOT_H - padT - padB;
    let lo = d.min, hi = d.max;
    if (hi <= lo) { hi = lo + 1e-9; }
    const px = (t) => padL + t * iw;
    const py = (v) => padT + (1 - (v - lo) / (hi - lo)) * ih;

    // the window band: where the current window sits relative to this profile
    const [wl, wh] = viewer.windowRange();
    const yTop = py(Math.min(hi, Math.max(lo, wh)));
    const yBot = py(Math.min(hi, Math.max(lo, wl)));
    if (yBot > yTop) {
      pctx.fillStyle = accent;
      pctx.globalAlpha = 0.10;
      pctx.fillRect(padL, yTop, iw, yBot - yTop);
      pctx.globalAlpha = 1;
    }

    pctx.strokeStyle = fg0;
    pctx.lineWidth = 1.1;
    pctx.beginPath();
    let started = false;
    for (let i = 0; i < d.n; i++) {
      const v = d.v[i];
      if (!Number.isFinite(v)) { started = false; continue; }
      const X = px(i / (d.n - 1)), Y = py(v);
      if (!started) { pctx.moveTo(X, Y); started = true; } else pctx.lineTo(X, Y);
    }
    pctx.stroke();

    if (profileMark !== null) {
      const X = px(profileMark);
      pctx.strokeStyle = focus;
      pctx.beginPath(); pctx.moveTo(X + 0.5, padT); pctx.lineTo(X + 0.5, padT + ih); pctx.stroke();
      const i = clamp(Math.round(profileMark * (d.n - 1)), 0, d.n - 1);
      const v = d.v[i];
      if (Number.isFinite(v)) {
        pctx.fillStyle = focus;
        pctx.beginPath(); pctx.arc(X, py(v), 2.5, 0, Math.PI * 2); pctx.fill();
        const unit = getFieldName() ? fieldMeta(getFieldName()).unit : '';
        const txt = `${num(d.d[i], 0)} um · ${auto(v)} ${unit}`;
        const tw = pctx.measureText(txt).width;
        pctx.fillStyle = cssv(plot, '--bg-1');
        pctx.fillRect(clamp(X + 4, 0, w - tw - 6), padT, tw + 5, 12);
        pctx.fillStyle = fg0;
        pctx.textBaseline = 'top';
        pctx.fillText(txt, clamp(X + 6, 2, w - tw - 4), padT + 1.5);
      }
    }

    pctx.fillStyle = fg3;
    pctx.textBaseline = 'bottom';
    pctx.fillText('0', padL, PLOT_H - 2);
    const endTxt = `${num(d.lengthUm, 0)} um`;
    pctx.fillText(endTxt, w - padR - pctx.measureText(endTxt).width, PLOT_H - 2);
    pctx.textBaseline = 'top';
    const fn = getFieldName();
    pctx.fillText(`${auto(hi)}${fn ? ' ' + fieldMeta(fn).unit : ''}`, padL + 1, 1);
  }

  plot.addEventListener('pointermove', (e) => {
    if (!profileData) return;
    const r = plot.getBoundingClientRect();
    profileMark = clamp((e.clientX - r.left - 4) / Math.max(1, r.width - 8), 0, 1);
    paintPlot();
    drawAnnot();
  });
  plot.addEventListener('pointerleave', () => {
    if (profileMark === null) return;
    profileMark = null; paintPlot(); drawAnnot();
  });

  /* ======================================================= region panel = */

  const roiTblWrap = h('div.tblwrap');
  const roiHint = h('div.tool-hint');
  const roiClear = h('button.btn.ghost.sm', {
    onclick: () => { rois = []; roiSeq = 0; paintRois(); drawAnnot(); },
    title: 'Discard every region',
  }, 'clear all');
  const roiCount = h('span.mono.dim', '0');
  const regionEl = h('div.tool-panel',
    h('div.frow.tool-bar', roiCount, h('span.unit', `of ${MAX_ROI}`), h('span.grow'), roiClear),
    roiTblWrap, roiHint,
  );

  function computeRoi(r) {
    const rect = {
      x0: Math.min(r.x0, r.x1), x1: Math.max(r.x0, r.x1),
      y0: Math.min(r.y0, r.y1), y1: Math.max(r.y0, r.y1),
    };
    r.cells = rectCellStats(viewer, rect);
    r.fieldStats = viewer.field ? rectFieldStats(viewer.field, viewer.geom.ex, viewer.geom.ey, rect) : null;
    r.areaUm2 = (rect.x1 - rect.x0) * (rect.y1 - rect.y0);
    // The region's own exposure distribution, binned against the FIELD-WIDE
    // axis so it can be overlaid on the window/level histogram and compared
    // like for like. Same samples the statistics came from — the region is
    // walked once.
    r.hist = (viewer.field && r.cells && r.cells.values && r.cells.values.length)
      ? buildHistogram(r.cells.values, viewer.field.min, viewer.field.max)
      : null;
  }

  function paintRois() {
    regionEl.classList.toggle('empty', !rois.length);
    roiClear.disabled = !rois.length;
    roiCount.textContent = String(rois.length);
    clear(roiTblWrap);
    if (!rois.length) {
      roiHint.textContent = 'No region. Press R, then drag a rectangle over the tissue.';
      if (onRoiChange) onRoiChange(null);
      return;
    }
    roiHint.textContent = rois.length >= MAX_ROI
      ? `${MAX_ROI} regions is the limit — clear one to draw another.`
      : 'Composition is measured; every field statistic is computed on the served grid.';

    const fname = getFieldName();
    const meta = fname ? fieldMeta(fname) : null;
    const head = h('tr', h('th', 'region'), rois.map((r) => h('th.right', r.label)), h('th', ''));
    const rows = [];

    const row = (label, get, unitTxt, prov, title) => {
      rows.push(h('tr',
        h('td', { title: title || '' }, label),
        rois.map((r) => h('td.right.num', get(r), unitTxt ? h('span.unit', unitTxt) : null)),
        h('td', prov ? chip(prov.kind, prov.why) : null),
      ));
    };

    const measured = { kind: 'measured', why: 'Cell centroids and phenotype calls come from the segmented spatial assay; counting them inside a rectangle adds no model.' };
    const computedGrid = meta
      ? { kind: meta.provenance, why: `${meta.source}. Area statistics are over the served grid nodes inside the region, so they inherit the run dx — a region narrower than one dx reports nothing rather than an interpolated value.` }
      : null;
    const computedCell = meta
      ? { kind: meta.provenance, why: `${meta.source}, bilinearly interpolated at each measured cell centroid. An interpolation of a computed field is computed, never measured.` }
      : null;

    row('area', (r) => num(r.areaUm2 / 1e6, 4), 'mm2',
      { kind: 'computed', why: 'Rectangle drawn by the operator, in the section coordinate frame from the TMSG extent.' });
    row('cells', (r) => int(r.cells ? r.cells.total : null), 'n', measured);
    row('density', (r) => num(r.cells && r.areaUm2 ? r.cells.total / (r.areaUm2 / 1e6) : null, 0), 'n/mm2', measured);

    CELL_TYPES.forEach((t, i) => {
      rows.push(h('tr',
        h('td',
          h('span.roi-sw', { style: { background: `var(${CHAN_VARS[i]})` } }),
          h('span', t.toLowerCase())),
        rois.map((r) => h('td.right.num',
          r.cells && r.cells.total ? pct(r.cells.counts[i] / r.cells.total, 1) : NULL_GLYPH,
          h('span.unit', '%'))),
        h('td', i === 0 ? chip(measured.kind, measured.why) : null),
      ));
    });

    if (meta) {
      rows.push(h('tr.roi-sep', h('td', { colspan: String(rois.length + 2) }, `${meta.label} · ${meta.unit}`)));
      const fs = (r) => (r.fieldStats && r.fieldStats.n ? r.fieldStats : null);
      row('area mean', (r) => (fs(r) ? auto(fs(r).mean) : NULL_GLYPH), meta.unit, computedGrid,
        'Mean over the grid nodes inside the region.');
      row('area sd', (r) => (fs(r) ? auto(fs(r).sd) : NULL_GLYPH), meta.unit, null);
      row('p10 / p90', (r) => (fs(r) ? `${auto(fs(r).p10)} / ${auto(fs(r).p90)}` : NULL_GLYPH), meta.unit, null,
        'Quantiles from a 512-bin histogram CDF of the region, not a sort: worst case one bin of error.');
      row('grid nodes', (r) => (fs(r) ? int(fs(r).n) : NULL_GLYPH), 'n', null,
        'How many served grid nodes the region contains. A small count is a weak statistic.');
      const ex = (r) => (r.cells && r.cells.exposure ? r.cells.exposure : null);
      row('cell exposure μ', (r) => (ex(r) ? auto(ex(r).mean) : NULL_GLYPH), meta.unit, computedCell,
        'Field interpolated at every cell centroid in the region, then averaged.');
      row('tumour exposure μ', (r) => (ex(r) && ex(r).tumorN ? auto(ex(r).tumorMean) : NULL_GLYPH), meta.unit, computedCell,
        'Same, restricted to cells whose measured phenotype call is TUMOR.');
      row('tumour cells', (r) => (ex(r) ? int(ex(r).tumorN) : NULL_GLYPH), 'n', measured);
    }

    rows.push(h('tr',
      h('td', ''),
      rois.map((r) => h('td.right', h('button.btn.ghost.sm', {
        onclick: () => { rois = rois.filter((x) => x !== r); paintRois(); drawAnnot(); },
        title: `Remove ${r.label}`,
      }, 'remove'))),
      h('td', ''),
    ));

    roiTblWrap.appendChild(h('table.tbl.roi-tbl', h('thead', head), h('tbody', rows)));
    if (onRoiChange) onRoiChange(rois.length ? rois[rois.length - 1] : null);
  }

  /* ========================================================= interaction = */

  /** @returns {boolean} true when the tool consumed the gesture. */
  function pointerDown(world) {
    if (tool === 'pick' || !viewer.geom) return false;
    pending = { kind: tool === 'region' ? 'region' : 'line', x0: world.x, y0: world.y, x1: world.x, y1: world.y };
    return true;
  }

  function pointerMove(world) {
    if (!pending) return false;
    pending.x1 = world.x; pending.y1 = world.y;
    drawAnnot();
    return true;
  }

  function pointerUp(world) {
    if (!pending) return false;
    pending.x1 = world.x; pending.y1 = world.y;
    const dx = pending.x1 - pending.x0, dy = pending.y1 - pending.y0;
    const tiny = Math.hypot(dx, dy) * viewer.cam.zoom < 6;
    const geom = { x0: pending.x0, y0: pending.y0, x1: pending.x1, y1: pending.y1 };
    pending = null;
    if (tiny) { drawAnnot(); return true; }

    if (tool === 'profile') {
      profile = geom;
      computeProfile();
      paintProfile();
    } else if (tool === 'distance') {
      distance = geom;
    } else if (tool === 'region') {
      if (rois.length >= MAX_ROI) rois.shift();
      const r = { ...geom, label: 'R' + (++roiSeq) };
      computeRoi(r);
      rois.push(r);
      paintRois();
    }
    drawAnnot();
    if (onChange) onChange();
    return true;
  }

  function recompute() {
    if (profile) computeProfile();
    for (const r of rois) computeRoi(r);
    paintProfile();
    paintRois();
    drawAnnot();
  }

  function clearActive() {
    let did = false;
    if (pending) { pending = null; did = true; }
    if (tool === 'profile' && profile) { profile = null; profileData = null; profileMark = null; paintProfile(); did = true; }
    if (tool === 'distance' && distance) { distance = null; did = true; }
    if (tool === 'region' && rois.length) { rois.pop(); paintRois(); did = true; }
    if (did) drawAnnot();
    return did;
  }

  function clearAll() {
    profile = null; profileData = null; profileMark = null;
    distance = null; rois = []; roiSeq = 0; pending = null;
    paintProfile(); paintRois(); drawAnnot();
  }

  const ro = new ResizeObserver(() => { paintPlot(); });
  ro.observe(plot);

  paintProfile();
  paintRois();

  return {
    annot,
    profileEl,
    regionEl,
    draw: drawAnnot,
    recompute,
    clearActive,
    clearAll,
    pointerDown, pointerMove, pointerUp,
    get tool() { return tool; },
    setTool(t) {
      if (!TOOLS.includes(t)) return;
      tool = t;
      pending = null;
      drawAnnot();
    },
    get rois() { return rois; },
    get hasProfile() { return !!profile; },
    repaintPanels() { paintProfile(); paintRois(); },
    /** Only the profile plot — the one thing here that depends on the window.
     *  This runs on every frame of a window drag, so it must not rebuild the
     *  region table: that is a DOM rebuild under the operator's pointer, for a
     *  table of numbers the window does not change. */
    repaintPlot() { paintPlot(); },
  };
}
