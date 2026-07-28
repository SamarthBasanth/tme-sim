// main.js — boot and wiring.

import { h, $, raf } from './util/dom.js';
import { createStore, INITIAL } from './state.js';
import * as api from './api.js';
import { Viewer, CELL_TYPES } from './gl/viewer.js';
import { createSpecPanel } from './ui/spec.js';
import { createInspectPanel } from './ui/inspect.js';
import { createHistoryPanel } from './ui/history.js';
import { createOverlays } from './ui/overlays.js';
import { makeSplitter } from './ui/splitter.js';
import { toast, notifyError } from './ui/toast.js';
import { num, int, bytes, NULL_GLYPH } from './util/fmt.js';

const store = createStore({ ...INITIAL });

/* ============================================================ theme ==== */

const THEME_KEY = 'tmesim.theme';
function initTheme() {
  let t = null;
  try { t = localStorage.getItem(THEME_KEY); } catch { /* private mode */ }
  if (t !== 'light' && t !== 'dark') {
    t = window.matchMedia('(prefers-color-scheme: light)').matches ? 'light' : 'dark';
  }
  applyTheme(t);
}
function applyTheme(t) {
  document.documentElement.setAttribute('data-theme', t);
  const b = $('#btn-theme');
  if (b) { b.textContent = 'theme: ' + t; b.setAttribute('aria-pressed', String(t === 'light')); }
  try { localStorage.setItem(THEME_KEY, t); } catch { /* ignore */ }
  store.set({ theme: t });
}
function toggleTheme() {
  applyTheme(document.documentElement.getAttribute('data-theme') === 'light' ? 'dark' : 'light');
}

/* ============================================================= boot ==== */

let viewer = null;
let overlays = null;
let specPanel = null;
let inspect = null;
let history = null;

async function boot() {
  initTheme();

  /* --- viewer --------------------------------------------------------- */
  const canvas = $('#gl');
  try {
    viewer = new Viewer(canvas);
  } catch (err) {
    $('#canvas-wrap').appendChild(h('div.ov.ov-empty',
      h('div.empty',
        h('h3', 'WebGL2 unavailable'),
        h('p', 'This workstation requires WebGL2 for the cell viewer. ' +
          'The run control, history, and provenance panels still work; the viewport does not.'),
        h('p.mono.dim', { style: { marginTop: '6px' } }, String(err && err.message)))));
    notifyError('WebGL2 context could not be created', err);
  }

  if (viewer) {
    overlays = createOverlays(viewer);
    for (const n of overlays.nodes) $('#canvas-wrap').appendChild(n);
    overlays.setEmpty(true);
  }

  /* --- panels --------------------------------------------------------- */
  specPanel = createSpecPanel({ store, onSubmit: submitRun });
  inspect = createInspectPanel({
    viewer: viewer || stubViewer(),
    onChange: () => { if (viewer) viewer.dirty = true; updateColourLabel(); },
    onBenchmark: runBenchmark,
    onResetView: () => { if (viewer) { viewer.fit(); tickOverlays(); } },
  });
  history = createHistoryPanel({ onOpen: openRun, onRefresh: refreshRuns });

  const body = $('#body');
  body.insertBefore(specPanel.el, body.firstChild);
  body.appendChild(inspect.el);
  body.appendChild(history.el);

  inspect.onFieldSelect = selectField;

  /* --- splitters ------------------------------------------------------ */
  makeSplitter($('#split-l'), { axis: 'x', varName: '--rail-l', min: 240, max: 560, onChange: onLayout });
  makeSplitter($('#split-r'), { axis: 'x', varName: '--rail-r', min: 250, max: 560, invert: true, onChange: onLayout });
  makeSplitter($('#split-h'), { axis: 'y', varName: '--dock-h', min: 90, max: 560, invert: true, onChange: onLayout });

  /* --- chrome --------------------------------------------------------- */
  $('#btn-theme').addEventListener('click', () => {
    toggleTheme();
    if (viewer) { viewer.refreshTheme(); }
    if (overlays) overlays.refreshRamp();
  });
  const dlg = $('#help-dlg');
  $('#btn-help').addEventListener('click', () => dlg.showModal());
  $('#help-close').addEventListener('click', () => dlg.close());

  $('#tb-fit').addEventListener('click', () => { if (viewer) { viewer.fit(); tickOverlays(); } });
  $('#tb-zoom-in').addEventListener('click', () => zoomCentre(1.35));
  $('#tb-zoom-out').addEventListener('click', () => zoomCentre(1 / 1.35));
  $('#tb-pin-clear').addEventListener('click', clearPin);

  if (viewer) wireCanvas(canvas);
  wireKeys();

  window.addEventListener('resize', onLayout);
  new ResizeObserver(onLayout).observe($('#canvas-wrap'));

  api.api.onModeChange = paintBackendState;

  /* --- data ----------------------------------------------------------- */
  const probe = await api.probe();
  paintBackendState(api.api);
  if (!probe.live) {
    toast('warn',
      'Backend unreachable — running on STUB DATA',
      'Every number on screen is synthesised by js/mock.js against API_CONTRACT.md. Start the server with ./.venv/bin/python -m app.server and reload.',
      0);
  }

  try {
    const [sections, drugs, defaults] = await Promise.all([
      api.getSections(), api.getDrugs(), api.getDefaults(),
    ]);
    store.set({ sections, drugs });
    specPanel.setSections(sections);
    specPanel.setDrugs(drugs);
    specPanel.loadSpec(defaults);
    paintSectionHeader(defaults.section);
  } catch (err) {
    notifyError('Failed to load catalogue (/api/sections, /api/drugs, /api/defaults)', err);
  }

  await refreshRuns();

  // Open the newest completed run so the workstation is not empty on arrival.
  const done = store.get().runs.find((r) => r.status === 'done');
  if (done) openRun(done.run_id);

  if (viewer) viewer.start();
  setInterval(paintStatusBar, 250);
  setInterval(refreshRuns, 6000);
}

function stubViewer() {
  // Keeps the inspect panel constructible when WebGL2 is missing.
  return {
    visible: [1, 1, 1, 1, 1], rampName: 'viridis', dirty: false,
    setVisible() {}, setColorMode() {}, setRamp() {}, fit() {},
  };
}

/* ========================================================== chrome ===== */

function paintBackendState(a) {
  const led = $('#hdr-led');
  const txt = $('#hdr-backend');
  if (!led || !txt) return;
  if (a.stub) {
    led.className = 'led warn';
    txt.textContent = 'STUB DATA · ' + (a.stubReason || 'backend unreachable');
    txt.style.color = 'var(--warn)';
    txt.title = 'The FastAPI backend did not answer. Responses are synthesised locally by js/mock.js against API_CONTRACT.md and must not be read as results.';
  } else if (a.reachable) {
    led.className = 'led ok';
    txt.textContent = 'live · ' + location.host;
    txt.style.color = 'var(--fg-1)';
    txt.title = 'Connected to the FastAPI backend.';
  } else {
    led.className = 'led';
    txt.textContent = 'probing…';
  }
}

function paintSectionHeader(sectionId) {
  const s = (store.get().sections || []).find((x) => x.id === sectionId);
  $('#hdr-section').textContent = sectionId || NULL_GLYPH;
  $('#hdr-section-meta').textContent = s
    ? `${int(s.n_cells)} cells · ${num(s.extent_um[0] / 1000, 2)}×${num(s.extent_um[1] / 1000, 2)} mm · ${s.panel}`
    : '';
}

function paintStatusBar() {
  const s = viewer ? viewer.stats : null;
  const set = (id, v) => { const el = $(id); if (el) el.textContent = v; };
  if (s) {
    set('#sb-mode', s.mode);
    set('#sb-fps', s.fps ? s.fps.toFixed(1) : NULL_GLYPH);
    set('#sb-frame', s.frameMs ? s.frameMs.toFixed(2) : NULL_GLYPH);
    set('#sb-drawn', s.total ? `${int(s.drawn)} / ${int(s.total)}` : NULL_GLYPH);
    set('#sb-vbo', s.gpuBytes ? bytes(s.gpuBytes) : NULL_GLYPH);
    set('#sb-index', s.indexBytes ? `${bytes(s.indexBytes)} · ${num(s.indexMs, 0)} ms` : NULL_GLYPH);
    const fpsEl = $('#sb-fps');
    if (fpsEl) fpsEl.style.color = s.fps && s.fps < 30 ? 'var(--warn)' : 'var(--fg-1)';
  }
  set('#sb-latency', api.api.lastLatencyMs === null ? NULL_GLYPH : num(api.api.lastLatencyMs, 1));
}

function updateColourLabel() {
  const el = $('#tb-colour');
  if (!el || !viewer) return;
  el.textContent = viewer.colorMode === 1 ? `field: ${store.get().fieldName || 'none'}` : 'cell type';
}

const tickOverlays = raf(() => {
  if (!viewer || !overlays) return;
  overlays.updateScalebar();
  const el = $('#tb-zoom'); if (el) el.textContent = num(viewer.cam.zoom, 4);
  const el2 = $('#tb-scale'); if (el2) el2.textContent = num(viewer.umPerPx, 3);
});

function onLayout() {
  if (viewer) { viewer.dirty = true; }
  tickOverlays();
}

/* ======================================================= interaction === */

function zoomCentre(f) {
  if (!viewer) return;
  const { w, h: hh } = viewer.cssSize;
  viewer.zoomAt(w / 2, hh / 2, f);
  tickOverlays();
}

function clearPin() {
  if (!viewer) return;
  viewer.selected = -1;
  store.set({ pinned: null });
  viewer.dirty = true;
  overlays.updateReadout(viewer.hovered >= 0 ? viewer.cellInfo(viewer.hovered) : null, store.get().fieldName, false);
}

function wireCanvas(canvas) {
  let panning = false, moved = 0, lastX = 0, lastY = 0, pid = null;

  canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    const r = canvas.getBoundingClientRect();
    const f = Math.exp(-e.deltaY * (e.deltaMode === 1 ? 0.02 : 0.0018));
    viewer.zoomAt(e.clientX - r.left, e.clientY - r.top, Math.max(0.3, Math.min(3, f)));
    tickOverlays();
  }, { passive: false });

  canvas.addEventListener('pointerdown', (e) => {
    if (e.button !== 0) return;
    panning = true; moved = 0; lastX = e.clientX; lastY = e.clientY; pid = e.pointerId;
    canvas.setPointerCapture(pid);
    canvas.classList.add('panning');
    canvas.focus();
  });

  const hoverUpdate = raf(() => {
    if (!viewer.geom) return;
    const info = viewer.hovered >= 0 ? viewer.cellInfo(viewer.hovered) : null;
    const pinnedInfo = store.get().pinned;
    overlays.updateReadout(
      info || (pinnedInfo ? viewer.cellInfo(pinnedInfo) : null),
      store.get().fieldName,
      !info && !!pinnedInfo);
  });

  canvas.addEventListener('pointermove', (e) => {
    const r = canvas.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    if (panning) {
      moved += Math.abs(e.clientX - lastX) + Math.abs(e.clientY - lastY);
      viewer.panByScreen(e.clientX - lastX, e.clientY - lastY);
      lastX = e.clientX; lastY = e.clientY;
      tickOverlays();
      return;
    }
    const w = viewer.screenToWorld(px, py);
    const cur = $('#tb-cursor');
    if (cur) cur.textContent = `${num(w.x, 1)}, ${num(w.y, 1)}`;
    const idx = viewer.pick(px, py);
    if (idx !== viewer.hovered) {
      viewer.hovered = idx;
      viewer.dirty = true;
      hoverUpdate();
    }
  });

  const endPan = (e) => {
    if (!panning) return;
    panning = false;
    canvas.classList.remove('panning');
    try { canvas.releasePointerCapture(pid); } catch { /* already released */ }
    if (moved < 4) {
      const r = canvas.getBoundingClientRect();
      const idx = viewer.pick(e.clientX - r.left, e.clientY - r.top);
      viewer.selected = idx;
      store.set({ pinned: idx >= 0 ? idx : null });
      viewer.dirty = true;
      overlays.updateReadout(idx >= 0 ? viewer.cellInfo(idx) : null, store.get().fieldName, idx >= 0);
    }
  };
  canvas.addEventListener('pointerup', endPan);
  canvas.addEventListener('pointercancel', () => { panning = false; canvas.classList.remove('panning'); });
  canvas.addEventListener('pointerleave', () => {
    if (panning) return;
    viewer.hovered = -1;
    viewer.dirty = true;
    const p = store.get().pinned;
    overlays.updateReadout(p !== null ? viewer.cellInfo(p) : null, store.get().fieldName, p !== null);
  });
}

function wireKeys() {
  window.addEventListener('keydown', (e) => {
    const t = e.target;
    const typing = t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.tagName === 'SELECT' || t.isContentEditable);

    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
      e.preventDefault(); submitRun(store.get().spec); return;
    }
    if ((e.ctrlKey || e.metaKey) && (e.key === 'l' || e.key === 'L')) {
      e.preventDefault(); specPanel.focusEditor(); return;
    }
    if (typing) return;
    if (e.ctrlKey || e.metaKey || e.altKey) return;

    const step = e.shiftKey ? 120 : 30;
    switch (e.key) {
      case 'f': case 'F': if (viewer) { viewer.fit(); tickOverlays(); } break;
      case '+': case '=': zoomCentre(1.3); break;
      case '-': case '_': zoomCentre(1 / 1.3); break;
      case 'ArrowLeft': if (viewer) { viewer.panByScreen(step, 0); tickOverlays(); e.preventDefault(); } break;
      case 'ArrowRight': if (viewer) { viewer.panByScreen(-step, 0); tickOverlays(); e.preventDefault(); } break;
      case 'ArrowUp': if (viewer) { viewer.panByScreen(0, step); tickOverlays(); e.preventDefault(); } break;
      case 'ArrowDown': if (viewer) { viewer.panByScreen(0, -step); tickOverlays(); e.preventDefault(); } break;
      case 'Escape': clearPin(); break;
      case 't': case 'T':
        toggleTheme(); if (viewer) viewer.refreshTheme(); if (overlays) overlays.refreshRamp(); break;
      case 'c': case 'C':
        if (viewer) { viewer.setColorMode(viewer.colorMode ? 0 : 1); updateColourLabel(); } break;
      case '0':
        if (viewer) { CELL_TYPES.forEach((_, i) => viewer.setVisible(i, true)); inspect.buildChannels(viewer.counts, viewer.geom ? viewer.geom.n : 0); } break;
      default: {
        const n = Number(e.key);
        if (Number.isInteger(n) && n >= 1 && n <= CELL_TYPES.length && viewer) {
          viewer.setVisible(n - 1, !viewer.visible[n - 1]);
          inspect.buildChannels(viewer.counts, viewer.geom ? viewer.geom.n : 0);
          // re-apply visibility to the freshly built checkboxes
          CELL_TYPES.forEach((_, i) => {
            const cb = inspect.el.querySelectorAll('.chan input')[i];
            if (cb) { cb.checked = !!viewer.visible[i]; cb.closest('.chan').classList.toggle('off', !viewer.visible[i]); }
          });
        }
      }
    }
  });
}

/* ============================================================ data ===== */

async function refreshRuns() {
  try {
    const runs = await api.listRuns();
    store.set({ runs });
    history.setRuns(runs);
    const active = store.get().activeRunId;
    if (active) {
      const r = runs.find((x) => x.run_id === active);
      if (r && r.status === 'done' && (!store.get().activeRun || store.get().activeRun.status !== 'done')) {
        openRun(active);
      }
    }
  } catch (err) {
    const el = $('#sb-msg');
    if (el) el.textContent = 'GET /api/runs failed: ' + (err.message || err);
  }
}

async function openRun(id) {
  history.setActive(id);
  store.set({ activeRunId: id });
  $('#hdr-run').textContent = id;
  let run = null;
  try {
    run = await api.getRun(id);
  } catch (err) {
    notifyError(`GET /api/runs/${id} failed`, err);
    return;
  }
  store.set({ activeRun: run });
  inspect.renderRun(run);
  const st = $('#hdr-run-status');
  while (st.firstChild) st.removeChild(st.firstChild);
  st.appendChild(h('span.chip.st-' + run.status, run.status));
  if (run.status === 'failed' && run.error) {
    toast('err', `Run ${id} failed`, String(run.error), 0);
  }
  paintSectionHeader(run.spec && run.spec.section);

  if (!viewer) return;
  if (run.status !== 'done') {
    viewer.clearGeometry();
    overlays.setEmpty(true);
    $('#hdr-cells').textContent = NULL_GLYPH;
    return;
  }
  await loadGeometry(id);
  if (store.get().fieldName) await selectField(store.get().fieldName);
}

async function loadGeometry(id) {
  const msg = $('#sb-msg');
  msg.textContent = `GET /api/runs/${id}/geometry …`;
  try {
    const t0 = performance.now();
    const g = await api.getGeometry(id);
    const fetchMs = performance.now() - t0;
    const t1 = performance.now();
    const { indexMs } = viewer.setGeometry(g);
    const uploadMs = performance.now() - t1;
    viewer.fit();
    overlays.setEmpty(false);
    inspect.buildChannels(viewer.counts, g.n);
    $('#hdr-cells').textContent = int(g.n);
    store.set({ geometryMeta: { n: g.n, bytes: g.bytes, fetchMs, uploadMs, indexMs } });
    msg.textContent = `geometry ${int(g.n)} cells · ${bytes(g.bytes)} · fetch+decode ${num(fetchMs, 1)} ms · upload+index ${num(uploadMs, 1)} ms`;
    tickOverlays();
  } catch (err) {
    msg.textContent = 'geometry load failed: ' + (err.message || err);
    notifyError(`GET /api/runs/${id}/geometry failed`, err);
    viewer.clearGeometry();
    overlays.setEmpty(true);
  }
}

async function selectField(name) {
  store.set({ fieldName: name || null });
  updateColourLabel();
  if (!viewer) return;
  if (!name) {
    viewer.setField(null);
    store.set({ field: null });
    overlays.updateLegend(null, null);
    inspect.setFieldInfo(null, null);
    return;
  }
  const id = store.get().activeRunId;
  if (!id) {
    toast('warn', 'No run selected', 'A field overlay is a property of a run. Open a completed run first.');
    inspect.selectField(null);
    store.set({ fieldName: null });
    return;
  }
  const msg = $('#sb-msg');
  msg.textContent = `GET /api/runs/${id}/field/${name} …`;
  try {
    const f = await api.getField(id, name);
    inspect.setRampFor(name);
    viewer.setField(f);
    store.set({ field: f });
    overlays.updateLegend(name, f);
    inspect.setFieldInfo(name, f);
    overlays.refreshRamp();
    msg.textContent = `field ${name} · ${f.rows}×${f.cols} @ ${num(f.dx, 2)} um · ${bytes(f.bytes)}`;
  } catch (err) {
    msg.textContent = `field ${name} failed: ` + (err.message || err);
    notifyError(`GET /api/runs/${id}/field/${name} failed`, err);
    inspect.selectField(null);
    store.set({ fieldName: null });
    viewer.setField(null);
    overlays.updateLegend(null, null);
  }
}

async function submitRun(spec) {
  if (!spec) return;
  specPanel.setSubmitting(true);
  try {
    const res = await api.postRun(spec);
    const id = res && (res.run_id || res.id);
    if (!id) throw new api.ApiError('POST /api/runs returned no run_id', { status: 502, payload: res });
    toast('ok', 'Run queued', id, 3500);
    await refreshRuns();
    history.setActive(id);
    store.set({ activeRunId: id });
    $('#hdr-run').textContent = id;
    specPanel.track(id, async (doneId, status) => {
      await refreshRuns();
      if (status === 'done') openRun(doneId);
      else if (status === 'failed') {
        const r = await api.getRun(doneId).catch(() => null);
        toast('err', `Run ${doneId} failed`, (r && r.error) || 'no error message returned', 0);
        if (r) { store.set({ activeRun: r }); inspect.renderRun(r); }
      } else {
        toast('warn', `Run ${doneId} ${status}`, '', 4000);
      }
    });
  } catch (err) {
    notifyError('POST /api/runs failed', err);
  } finally {
    specPanel.setSubmitting(false);
  }
}

/* ======================================================= benchmark ===== */

const BENCH_N = 850000;

async function runBenchmark(outEl) {
  if (!viewer) return;
  outEl.textContent = `loading ${int(BENCH_N)}-cell synthetic TMSG payload…`;
  await new Promise((r) => setTimeout(r, 30));
  const res = await benchmark850k((s) => { outEl.textContent = s; });
  outEl.innerHTML = '';
  outEl.appendChild(h('div', `n = ${int(res.n)} cells, ${res.canvas[0]}×${res.canvas[1]} px @ dpr ${res.dpr}`));
  outEl.appendChild(h('div', `presented ${res.medianMs.toFixed(2)} ms/f  →  `,
    h('b', { style: { color: res.fpsMedian >= 30 ? 'var(--ok)' : 'var(--crit)' } }, res.fpsMedian.toFixed(1) + ' fps'),
    ' (vsync-capped)'));
  outEl.appendChild(h('div', `p95 ${res.p95Ms.toFixed(2)} ms → ${res.fpsP5.toFixed(1)} fps floor`));
  outEl.appendChild(h('div', `gpu-bound ${res.gpuMedianMs.toFixed(2)} ms/f → ${res.gpuFpsMedian.toFixed(1)} fps ceiling`));
  outEl.appendChild(h('div', `decode ${res.decodeMs.toFixed(1)} ms · upload+index ${res.uploadMs.toFixed(1)} ms · vbo ${bytes(res.gpuBytes)}`));
  toast(res.fpsMedian >= 30 ? 'ok' : 'err',
    `850k benchmark: ${res.fpsMedian.toFixed(1)} fps median`,
    `${res.medianMs.toFixed(2)} ms/frame, p95 ${res.p95Ms.toFixed(2)} ms. Viewer now holds the synthetic phantom, not a run.`,
    9000);
  $('#hdr-section').textContent = 'phantom_850k';
  $('#hdr-section-meta').textContent = 'synthetic render-stress payload — not a specimen';
  return res;
}

/**
 * The measured number. Builds a real 850k TMSG buffer, runs it through the
 * production decode -> quantise -> upload path, then times full-detail
 * rendering with gl.finish() per frame while the camera sweeps.
 */
async function benchmark850k(progress = () => {}) {
  const t0 = performance.now();
  const buf = api.synthTMSG(BENCH_N, 24000, 18000);
  const g = api.decodeTMSG(buf);
  const decodeMs = performance.now() - t0;
  progress(`decoded ${int(g.n)} cells (${bytes(buf.byteLength)}) in ${decodeMs.toFixed(1)} ms — uploading…`);
  await new Promise((r) => setTimeout(r, 20));

  const t1 = performance.now();
  viewer.setGeometry(g);
  const uploadMs = performance.now() - t1;
  viewer.fit();
  if (overlays) { overlays.setEmpty(false); }
  if (inspect) inspect.buildChannels(viewer.counts, g.n);
  $('#hdr-cells').textContent = int(g.n);
  progress(`uploaded in ${uploadMs.toFixed(1)} ms — rendering 180 frames…`);
  await new Promise((r) => setTimeout(r, 20));

  const res = await viewer.benchmark({ frames: 180, pan: true });
  res.decodeMs = decodeMs;
  res.uploadMs = uploadMs;
  res.gpuBytes = viewer.stats.gpuBytes;
  res.indexMs = viewer.stats.indexMs;
  res.indexBytes = viewer.stats.indexBytes;
  return res;
}

/* --- exposed for headless perf runs ------------------------------------ */
window.__tmesim = {
  get viewer() { return viewer; },
  get store() { return store; },
  api,
  benchmark850k,
  runBenchmark,
};

boot().catch((err) => {
  // eslint-disable-next-line no-console
  console.error('boot failed', err);
  notifyError('Boot failed', err);
});
