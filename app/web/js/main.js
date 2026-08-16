// main.js — boot and wiring.

import { h, $, raf, wireEdgeShadow } from './util/dom.js';
import { createStore, INITIAL } from './state.js';
import * as api from './api.js';
import { Viewer, CELL_TYPES, COLOR_MODE } from './gl/viewer.js';
import { createSpecPanel } from './ui/spec.js';
import { createInspectPanel } from './ui/inspect.js';
import { createHistoryPanel } from './ui/history.js';
import { createReportsPanel } from './ui/reports.js';
import { createProjectsPanel } from './ui/projects.js';
import { createOverlays } from './ui/overlays.js';
import { createMechanismLegend, createTimelapseScrubber, fmtTime as fmtClock } from './ui/mechanism.js';
import { createWindowLevel } from './ui/wl.js';
import { createTools, TOOLS } from './ui/tools.js';
import { makeSplitter } from './ui/splitter.js';
import { toast, notifyError, notifyRunFailed } from './ui/toast.js';
import { fieldMeta } from './fields.js';
import { setSynthetic, setActiveRun } from './ui/provenance.js';
import { num, int, auto, bytes, dur, NULL_GLYPH } from './util/fmt.js';
import { ensureAuthenticated } from './ui/apikey.js';
import { shouldShowWelcome, createWelcomeOverlay } from './ui/welcome.js';

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
/**
 * Flip the theme and repaint everything that caches a colour.
 *
 * ONE routine, because there used to be two: the #btn-theme click handler and
 * the `T` key case were parallel copies of the same five calls, and they had
 * drifted — the click path called mechLegend.refresh() and the key path did
 * not, so pressing T left the mechanism legend painted in the outgoing theme's
 * colours until something else forced a re-render.
 */
function applyThemeChange() {
  toggleTheme();
  if (viewer) viewer.refreshTheme();
  if (overlays) overlays.refreshRamp();
  if (mechLegend) mechLegend.refresh();
  if (wl) wl.draw();
  if (tools) { tools.repaintPanels(); tools.draw(); }
}

/* ============================================================= boot ==== */

let viewer = null;
let overlays = null;
let mechLegend = null;
let scrubber = null;
// True while the viewer's field texture is holding a TIME-LAPSE frame rather
// than a served /field/{name} overlay. It is tracked separately from
// store.fieldName because that name is the SERVER field selector: writing
// 'drug_concentration' into it would make the next openRun ask the backend for
// a field it does not serve.
let timelapseField = false;
let specPanel = null;
let inspect = null;
let history = null;
let reportsPanel = null;
let projectsPanel = null;
let wl = null;
let tools = null;

async function boot() {
  initTheme();

  /* --- viewer --------------------------------------------------------- */
  const canvas = $('#gl');
  try {
    viewer = new Viewer(canvas);
  } catch (err) {
    $('#canvas-wrap').appendChild(h('div.ov.ov-empty',
      h('div.notice',
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

    // The mechanism breakdown is built here (it needs the viewer, whose label
    // visibility its rows drive) but MOUNTED in the READ rail, below —
    // `mechanismEl` on createInspectPanel, the same seam windowLevelEl and the
    // three tool panels already use. It used to be appended to #canvas-wrap as
    // a floating card and it occluded the stage's top-right quarter in every
    // colour mode; see js/ui/mechanism.js for the full argument.
    mechLegend = createMechanismLegend(viewer, {
      onToggle: () => { viewer.dirty = true; },
    });
    // The scrubber stays on the stage: it is a transport control for the
    // picture, so it has to be under the picture it scrubs.
    scrubber = createTimelapseScrubber({ onFrame: onTimelapseFrame });
    $('#canvas-wrap').appendChild(scrubber.node);

    wl = createWindowLevel({ viewer, onChange: onWindowChange });
    tools = createTools({
      viewer,
      getFieldName: () => store.get().fieldName,
      onChange: () => { paintToolState(); },
      onRoiChange: (roi) => { if (wl) wl.setRoiHist(roi && roi.hist ? { hist: roi.hist, label: roi.label } : null); },
    });
    $('#canvas-wrap').appendChild(tools.annot);
    // One repaint per presented frame keeps the annotation layer glued to the
    // tissue through every camera path — inertia, LOD relax, keyboard pan —
    // without adding anything to the GL pass. tickOverlays is rAF-coalesced,
    // and the layer early-outs when there is nothing drawn on it.
    viewer.onFrame = () => { if (viewer.interacting) tickOverlays(); };
  }

  /* --- panels --------------------------------------------------------- */
  specPanel = createSpecPanel({
    store,
    onSubmit: submitRun,
    // A user-compiled compound joins (or leaves) the drug catalogue, so every
    // panel that lists drugs is refetched — patching only the picker would
    // leave the reports panel offering a compound that no longer exists.
    onDrugsChanged: () => refreshDrugs(),
    // Ingesting or removing a section changes the section catalogue the same
    // way: the picker, the projects panel and the protocol validator all read
    // it, so it is refetched once here rather than patched in three places.
    onSectionsChanged: (sectionId) => refreshSections(sectionId),
  });
  inspect = createInspectPanel({
    viewer: viewer || stubViewer(),
    onChange: () => { if (viewer) { viewer.dirty = true; } syncColourControl(); if (wl) wl.draw(); },
    onBenchmark: runBenchmark,
    onResetView: () => { if (viewer) { viewer.fit(); tickOverlays(); } },
    windowLevelEl: wl ? wl.el : null,
    profileEl: tools ? tools.profileEl : null,
    regionEl: tools ? tools.regionEl : null,
    distanceEl: tools ? tools.distanceEl : null,
    mechanismEl: mechLegend ? mechLegend.node : null,
  });
  // The tracked job's poll payload is already in the store (spec.js `poll`);
  // this is the only subscriber, and it exists so the top bar tells the truth
  // about the run it names for as long as that run is in flight.
  store.on('status', (s) => {
    const st = store.get();
    if (!s || !st.trackedRunId || st.trackedRunId !== st.activeRunId) return;
    paintRunChip(s.status, s.elapsed_s);
  });

  history = createHistoryPanel({ onOpen: openRun, onRefresh: refreshRuns });
  // Additive: the reports panel mounts as a third dock column and brings its
  // own stylesheet, so nothing above it has to change.
  reportsPanel = createReportsPanel({
    store,
    getSection: () => (store.get().spec && store.get().spec.section)
      || (store.get().activeRun && store.get().activeRun.spec
        && store.get().activeRun.spec.section) || null,
  });
  history.el.classList.add('with-reports');
  history.el.appendChild(h('div.dock-div'));
  history.el.appendChild(reportsPanel.el);
  // Additive: the project / specimen browser and the cross-specimen comparison
  // live in their own modal with their own stylesheet, and inject their own
  // toolbar button — the three-rail layout is untouched.
  projectsPanel = createProjectsPanel({ store });

  const body = $('#body');
  body.insertBefore(specPanel.el, body.firstChild);
  body.appendChild(inspect.el);
  body.appendChild(history.el);

  inspect.onFieldSelect = selectField;

  /* --- splitters ------------------------------------------------------ */
  makeSplitter($('#split-l'), {
    axis: 'x', varName: '--rail-l', min: 240, max: 560, onChange: onLayout,
    ariaLabel: 'Resize the setup rail. Arrow keys resize, Shift+Arrow resizes further, Home/End snap to the limits.',
  });
  makeSplitter($('#split-r'), {
    axis: 'x', varName: '--rail-r', min: 250, max: 560, invert: true, onChange: onLayout,
    ariaLabel: 'Resize the read rail. Arrow keys resize, Shift+Arrow resizes further, Home/End snap to the limits.',
  });
  makeSplitter($('#split-h'), {
    axis: 'y', varName: '--dock-h', min: 90, max: 560, invert: true, onChange: onLayout,
    ariaLabel: 'Resize the run history dock. Arrow keys resize, Shift+Arrow resizes further, Home/End snap to the limits.',
  });

  /* --- chrome --------------------------------------------------------- */
  $('#btn-theme').addEventListener('click', applyThemeChange);
  const dlg = $('#help-dlg');
  $('#btn-help').addEventListener('click', () => dlg.showModal());
  $('#help-close').addEventListener('click', () => dlg.close());

  $('#tb-fit').addEventListener('click', () => { if (viewer) { viewer.fit(); tickOverlays(); } });
  $('#tb-zoom-in').addEventListener('click', () => zoomCentre(1.35));
  $('#tb-zoom-out').addEventListener('click', () => zoomCentre(1 / 1.35));
  $('#tb-pin-clear').addEventListener('click', clearPin);

  for (const b of $('#tb-tools').children) {
    b.addEventListener('click', () => setTool(b.dataset.tool));
  }
  for (const b of $('#tb-colour').children) {
    b.addEventListener('click', () => setColourMode(Number(b.dataset.mode)));
  }
  // Below ~1150px the toolbar no longer fits its own tools + colour-by +
  // pin-clear and falls back to its own horizontal scroll (layout.css); mark
  // the edges so that is visible instead of reading as controls that vanished.
  wireEdgeShadow($('.stage-toolbar'));

  if (viewer) { wireCanvas(canvas); setTool('pick'); }
  syncColourControl();
  wireKeys();

  window.addEventListener('resize', onLayout);
  new ResizeObserver(onLayout).observe($('#canvas-wrap'));

  api.api.onModeChange = paintBackendState;

  /* --- data ----------------------------------------------------------- */
  // Blocks here (showing the access-key gate) if the backend is reachable
  // but requires a key this browser doesn't have or has wrong. A no-op when
  // no key is configured (the ordinary local/no-auth case) or when the
  // backend is genuinely unreachable (probe() below falls back to stub, as
  // before this gate existed).
  await ensureAuthenticated(api.probe);

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
    reportsPanel.setDrugs(drugs);
    specPanel.loadSpec(defaults);
    paintSectionHeader(defaults.section);
  } catch (err) {
    notifyError('Failed to load catalogue (/api/sections, /api/drugs, /api/defaults)', err);
  }

  specPanel.refreshCompounds();
  specPanel.refreshTissue();

  await refreshRuns();

  // Open the newest completed SINGLE-SPEC run so the workstation is not empty
  // on arrival. Protocol (TPL) runs are deliberately excluded here: they have
  // no single `/api/runs/{id}/geometry` payload the way a regular RunSpec run
  // does (a sweep can be many points, not one reconstruction), so opening one
  // here 404s on geometry load and leaves the viewer looking broken with only
  // a small corner toast explaining why -- confirmed live: the newest run in
  // this workstation's own history was a protocol run and did exactly that.
  const done = store.get().runs.find((r) => r.status === 'done' && r.spec && r.spec.kind !== 'protocol');
  if (done) openRun(done.run_id);

  // AFTER the auto-open, deliberately. A run still in flight has no geometry
  // and no metrics yet, so opening the last completed one first gives the
  // operator something real to look at while they wait; re-attaching second
  // then leaves the top bar naming the run that is actually happening rather
  // than the one on screen — which is the contradiction paintRunChip exists to
  // prevent, and whose tooltip already says the metrics belong to the previous
  // run.
  reattachLiveRun();

  if (viewer) viewer.start();

  const openTour = () => {
    const welcome = createWelcomeOverlay({
      onGoToProtocol: () => { const b = $('#btn-protocol'); if (b) b.click(); },
    });
    document.body.appendChild(welcome.node);
  };
  if (shouldShowWelcome()) openTour();
  const tourBtn = $('#btn-tour');
  if (tourBtn) tourBtn.addEventListener('click', openTour);

  setInterval(paintStatusBar, 250);
  setInterval(refreshRuns, 6000);
  // Reports poll fast while arms are in flight, slowly when idle. A report is
  // 2 runs per drug, so it is minutes long — a fixed fast poll would be noise.
  await reportsPanel.refresh();
  (async function pollReports() {
    for (;;) {
      await new Promise((r) => setTimeout(r, reportsPanel.busy() ? 2500 : 15000));
      await reportsPanel.refresh();
    }
  }());
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
  // The STUB state is a property of the whole session, not of the status bar.
  // Two things are set here and NOTHING else has to know about the stub:
  //   * data-mode on <html> drives the shell-wide chrome in css/components.css
  //     (striped frame, viewport tag, coloured section header) — every panel is
  //     inside it, so there is no glance that misses it;
  //   * the provenance latch re-labels EVERY chip in EVERY panel as SYNTHETIC,
  //     which is what stops a fabricated cell count wearing a green "measured".
  // Both are re-evaluated on every mode change, so a mid-session degradation is
  // caught too, not just one that was present at boot.
  const stubWhy = 'The FastAPI backend did not answer (' + (a.stubReason || 'backend unreachable')
    + '). Every response is being synthesised locally by js/mock.js from a seeded PRNG '
    + 'against the shape of API_CONTRACT.md. There is no tissue, no solve and no measurement behind any of it.';
  document.documentElement.setAttribute('data-mode', a.stub ? 'stub' : 'live');
  setSynthetic(!!a.stub, stubWhy);
  // Re-render whatever is already on screen under the new labelling.
  if (a.stub && inspect && store.get().activeRun) inspect.renderRun(store.get().activeRun);
  if (a.stub && history) history.setRuns(store.get().runs || []);

  const led = $('#hdr-led');
  const txt = $('#hdr-backend');
  if (!led || !txt) return;
  if (a.stub) {
    led.className = 'led warn';
    txt.textContent = 'STUB DATA · ' + (a.stubReason || 'backend unreachable');
    txt.style.color = 'var(--synthetic)';
    txt.title = stubWhy;
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
  if (reportsPanel) reportsPanel.setSection();
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
    set('#sb-lod', s.total ? (s.lod * 100).toFixed(0) : NULL_GLYPH);
    set('#sb-vbo', s.gpuBytes ? bytes(s.gpuBytes) : NULL_GLYPH);
    set('#sb-index', s.indexBytes ? `${bytes(s.indexBytes)} · ${num(s.indexMs, 0)} ms` : NULL_GLYPH);
    const fpsEl = $('#sb-fps');
    if (fpsEl) fpsEl.style.color = s.fps && s.fps < 30 ? 'var(--warn)' : 'var(--fg-1)';
  }
  set('#sb-latency', api.api.lastLatencyMs === null ? NULL_GLYPH : num(api.api.lastLatencyMs, 1));
}

/* ----------------------------------------------------------- controls -- */

function setColourMode(m) {
  if (!viewer) return;
  if (m === COLOR_MODE.FIELD && !store.get().fieldName && !timelapseField) {
    toast('warn', 'No field overlay', 'Colouring cells by field value needs a field. Pick one under Field overlay.');
    return;
  }
  if (m === COLOR_MODE.LABEL && !viewer.labels) {
    const r = store.get().activeRun;
    toast('warn', 'No mechanism attribution for this run',
      (r && r.mechanism_reason)
      || 'This run carries no per-cell attribution. Submit a run with '
         + 'attribute_mechanism (the default) to get one.');
    return;
  }
  viewer.setColorMode(m);
  syncColourControl();
  // Selecting `mech` promotes the breakdown from "the run's finding" to "the
  // key for what is on screen right now". It lives in the rail rather than on
  // the stage, so bring it to the operator instead of making them hunt for the
  // legend of the picture they are looking at.
  if (m === COLOR_MODE.LABEL && inspect && inspect.revealMechanism) {
    inspect.revealMechanism();
  }
}

/** Cycle 'colour by' through the modes this run can actually offer. */
function cycleColourMode() {
  if (!viewer) return;
  const available = [COLOR_MODE.TYPE];
  if (store.get().fieldName || timelapseField) available.push(COLOR_MODE.FIELD);
  if (viewer.labels) available.push(COLOR_MODE.LABEL);
  const at = available.indexOf(viewer.colorMode);
  setColourMode(available[(at + 1) % available.length]);
}

function syncColourControl() {
  const m = viewer ? viewer.colorMode : 0;
  // ONE availability rule, applied to both surfaces. The rail used to carry a
  // two-state checkbox with its own idea of the same setting; it now carries
  // the same three-state control, so the two cannot disagree about which mode
  // is even selectable.
  const disabled = [
    false,
    !store.get().fieldName && !timelapseField,
    !(viewer && viewer.labels),
  ];
  for (const b of $('#tb-colour').children) {
    const mode = Number(b.dataset.mode);
    b.setAttribute('aria-pressed', String(mode === m));
    b.disabled = !!disabled[mode];
  }
  if (inspect) inspect.setColorMode(m, disabled);
}

function setTool(t) {
  if (!tools || !TOOLS.includes(t)) return;
  tools.setTool(t);
  for (const b of $('#tb-tools').children) {
    b.setAttribute('aria-pressed', String(b.dataset.tool === t));
  }
  const canvas = $('#gl');
  canvas.classList.toggle('tooling', t !== 'pick');
  paintToolState();
}

function paintToolState() {
  const msg = $('#sb-msg');
  if (!msg || !tools) return;
  if (tools.tool === 'pick') {
    // Do not leave the previous tool's instruction sitting in the status
    // strip after the operator has left that tool.
    if (msg.dataset.tool) { msg.textContent = ''; delete msg.dataset.tool; }
    return;
  }
  msg.dataset.tool = tools.tool;
  msg.textContent = {
    profile: 'line profile — drag across the tissue; middle-drag pans',
    region: 'region — drag a rectangle; middle-drag pans',
    distance: 'distance — drag between two points; middle-drag pans',
  }[tools.tool] || '';
}

/** Anything that changes the window: relabel the legend, the plot band, the
 *  read-out. The pixels themselves are re-mapped by the shader. */
function onWindowChange() {
  if (overlays) overlays.refreshWindow();
  if (tools) tools.repaintPlot();
}

const tickOverlays = raf(() => {
  if (!viewer || !overlays) return;
  overlays.updateScalebar();
  const el2 = $('#tb-scale'); if (el2) el2.textContent = num(viewer.umPerPx, 3);
  if (tools) tools.draw();
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
  $('#tb-pin-clear').disabled = true;
  overlays.updateReadout(viewer.hovered >= 0 ? viewer.cellInfo(viewer.hovered) : null, store.get().fieldName, false);
}

/* Window / level by dragging on the image, the way every diagnostic viewer
   does it: horizontal moves the LEVEL (window centre), vertical changes the
   WIDTH. Right-drag is the convention; Alt+left-drag is the same gesture for a
   pointer with one button. Both are scaled to the full data range so the feel
   is identical whether the field spans 0–1 or 0–1e7. */
function applyWindowDrag(base, dxPx, dyPx, widthPx) {
  const [flo, fhi] = viewer.fullRange();
  const span = fhi - flo || 1;
  const centre = base.c + (dxPx / Math.max(120, widthPx)) * span * 1.6;
  const width = Math.max(span * 1e-4, base.w * Math.exp(dyPx * 0.006));
  viewer.setWindow(centre, width);
  if (wl) wl.sync();
}

function wireCanvas(canvas) {
  let mode = null;                // null | 'pan' | 'tool' | 'window'
  let moved = 0, lastX = 0, lastY = 0, pid = null;
  let base = null;                // window at gesture start

  // The canvas fills the fixed centre stage — it is never inside a scrolling
  // container (html/body are overflow:hidden by contract) — so its rect only
  // changes on resize, not on every pointer event. pointermove alone can fire
  // past 100/s while panning; reading layout that often for a rect that is
  // almost always unchanged was measurable, avoidable work. Cache it and
  // refresh only when the element itself actually resizes.
  let rect = canvas.getBoundingClientRect();
  new ResizeObserver(() => { rect = canvas.getBoundingClientRect(); }).observe(canvas);

  // Same reasoning for the status-bar read-outs this handler paints into on
  // every hover: three querySelector calls per mousemove, for elements that
  // never change identity, is three DOM walks the hot path does not need.
  const tbCursor = $('#tb-cursor');
  const tbValue = $('#tb-value');
  const tbValueUnit = $('#tb-value-unit');

  function paintCursorValue(world) {
    if (!tbValue) return;
    const name = store.get().fieldName;
    if (!name || !viewer.field) { tbValue.textContent = NULL_GLYPH; if (tbValueUnit) tbValueUnit.textContent = ''; return; }
    const v = viewer.sampleField(world.x, world.y);
    tbValue.textContent = v === null || !Number.isFinite(v) ? NULL_GLYPH : auto(v);
    if (tbValueUnit) tbValueUnit.textContent = fieldMeta(name).unit;
  }

  canvas.addEventListener('contextmenu', (e) => e.preventDefault());

  canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    const r = rect;
    const f = Math.exp(-e.deltaY * (e.deltaMode === 1 ? 0.02 : 0.0018));
    viewer.zoomAt(e.clientX - r.left, e.clientY - r.top, Math.max(0.3, Math.min(3, f)));
    tickOverlays();
  }, { passive: false });

  canvas.addEventListener('pointerdown', (e) => {
    const r = rect;
    const px = e.clientX - r.left, py = e.clientY - r.top;
    canvas.focus();
    moved = 0; lastX = e.clientX; lastY = e.clientY;

    const windowing = (e.button === 2 || (e.button === 0 && e.altKey)) && !!viewer.field;
    if (windowing) {
      const [lo, hi] = viewer.windowRange();
      base = { c: (lo + hi) / 2, w: hi - lo, x: e.clientX, y: e.clientY };
      mode = 'window';
      canvas.classList.add('windowing');
    } else if (e.button === 1 || (e.button === 0 && tools.tool === 'pick')) {
      mode = 'pan';
      canvas.classList.add('panning');
    } else if (e.button === 0) {
      mode = 'tool';
      tools.pointerDown(viewer.screenToWorld(px, py));
    } else {
      return;
    }
    pid = e.pointerId;
    canvas.setPointerCapture(pid);
    e.preventDefault();
  });

  const hoverUpdate = raf(() => {
    if (!viewer.geom) return;
    const info = viewer.hovered >= 0 ? viewer.cellInfo(viewer.hovered) : null;
    const pinnedInfo = store.get().pinned;
    const shown = info || (pinnedInfo !== null ? viewer.cellInfo(pinnedInfo) : null);
    overlays.updateReadout(shown, store.get().fieldName, !info && pinnedInfo !== null);
    // Brushing: the value of the cell under the pointer is marked on the
    // histogram, so "where does this cell sit in the distribution" is one
    // glance rather than an inference from a colour.
    if (wl) wl.setBrush(shown ? shown.fieldValue : null);
  });

  canvas.addEventListener('pointermove', (e) => {
    const r = rect;
    const px = e.clientX - r.left, py = e.clientY - r.top;
    const dx = e.clientX - lastX, dy = e.clientY - lastY;

    if (mode === 'pan') {
      moved += Math.abs(dx) + Math.abs(dy);
      viewer.panByScreen(dx, dy);
      lastX = e.clientX; lastY = e.clientY;
      tickOverlays();
      return;
    }
    if (mode === 'window') {
      applyWindowDrag(base, e.clientX - base.x, e.clientY - base.y, r.width);
      return;
    }
    if (mode === 'tool') {
      tools.pointerMove(viewer.screenToWorld(px, py));
      return;
    }

    const w = viewer.screenToWorld(px, py);
    if (tbCursor) tbCursor.textContent = `${num(w.x, 1)}, ${num(w.y, 1)}`;
    paintCursorValue(w);
    const idx = viewer.pick(px, py);
    if (idx !== viewer.hovered) {
      viewer.hovered = idx;
      viewer.dirty = true;
      hoverUpdate();
    }
  });

  const endGesture = (e) => {
    if (!mode) return;
    const r = rect;
    const was = mode;
    mode = null;
    canvas.classList.remove('panning', 'windowing');
    try { canvas.releasePointerCapture(pid); } catch { /* already released */ }
    if (was === 'tool') {
      tools.pointerUp(viewer.screenToWorld(e.clientX - r.left, e.clientY - r.top));
    } else if (was === 'pan' && moved < 4) {
      const idx = viewer.pick(e.clientX - r.left, e.clientY - r.top);
      viewer.selected = idx;
      store.set({ pinned: idx >= 0 ? idx : null });
      viewer.dirty = true;
      $('#tb-pin-clear').disabled = idx < 0;
      const info = idx >= 0 ? viewer.cellInfo(idx) : null;
      overlays.updateReadout(info, store.get().fieldName, idx >= 0);
      if (wl) wl.setBrush(info ? info.fieldValue : null);
    }
  };
  canvas.addEventListener('pointerup', endGesture);
  canvas.addEventListener('pointercancel', () => {
    mode = null;
    canvas.classList.remove('panning', 'windowing');
  });
  canvas.addEventListener('pointerleave', () => {
    if (mode) return;
    viewer.hovered = -1;
    viewer.dirty = true;
    const p = store.get().pinned;
    const info = p !== null ? viewer.cellInfo(p) : null;
    overlays.updateReadout(info, store.get().fieldName, p !== null);
    if (wl) wl.setBrush(info ? info.fieldValue : null);
    if (tbValue) tbValue.textContent = NULL_GLYPH;
  });
}

function wireKeys() {
  window.addEventListener('keydown', (e) => {
    const t = e.target;
    const typing = t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.tagName === 'SELECT' || t.isContentEditable);

    // Something nearer the event already claimed this key. The splitters
    // (js/ui/splitter.js) bind Arrow/Home/End on the separator elements and
    // preventDefault; this listener is on window and never saw that, so with a
    // splitter focused an arrow key BOTH resized the pane and panned the
    // viewer. defaultPrevented is the general form of that fix — any future
    // element-level binding is covered without either side knowing about the
    // other.
    if (e.defaultPrevented) return;

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
      case 'Escape':
        // The protocol workspace's own Esc-closes-it handler
        // (js/ui/protocol.js) never sees this keypress when the TPL editor
        // has focus: the editor claims Escape first, to blur itself back out
        // of the Tab-trap it needs for indentation, and stops the event from
        // bubbling to the workspace root. That leaves the workspace open and
        // this fallback — reached once focus has landed back on <body> — as
        // the only remaining path to the "Esc closes it" contract every
        // workspace is required to honour (docs/UI_ARCHITECTURE.md §2.7).
        // Checked first so a stray Esc after leaving the editor cannot also
        // fall through to the viewer commands below and mutate the specimen
        // hidden behind the open workspace.
        if (window.__tmesimProtocol && window.__tmesimProtocol.visible) {
          window.__tmesimProtocol.hide();
          break;
        }
        if (e.shiftKey && tools) { tools.clearAll(); break; }
        if (tools && tools.clearActive()) break;
        clearPin();
        break;
      case '?': $('#help-dlg').showModal(); break;
      case 't': case 'T':
        // Routed through the same function as #btn-theme. The two used to be
        // parallel copies of five calls, and the copies had drifted: the click
        // path refreshed the mechanism legend and the key path did not, so
        // pressing T left the legend painted in the previous theme's colours.
        applyThemeChange();
        break;
      case 'v': case 'V': setTool('pick'); break;
      case 'l': setTool('profile'); break;
      case 'L': setTool('profile'); break;
      case 'r': case 'R': setTool('region'); break;
      case 'm': case 'M': setTool('distance'); break;
      case 'w': if (wl) wl.autoWindow(0.02, 0.98); break;
      // `wl` is null when WebGL2 is unavailable. Every neighbouring branch
      // null-checks it; this one dereferenced it and would throw.
      case 'W': if (viewer && viewer.field) { viewer.resetWindow(); if (wl) wl.sync(); } break;
      // The scrubber's play button has always advertised `(space)` in its
      // tooltip and nothing registered it — scrubber.toggle() was dead code.
      // Bound only while a time-lapse is actually loaded, so space keeps its
      // native meaning everywhere else.
      case ' ':
        if (scrubber && scrubber.node && !scrubber.node.classList.contains('hidden')
            && !(t && (t.tagName === 'BUTTON' || t.tagName === 'A'))) {
          e.preventDefault();
          scrubber.toggle();
        }
        break;
      case 'c': case 'C': cycleColourMode(); break;
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

/**
 * Refetch the drug catalogue and repaint every panel that lists drugs.
 *
 * Called when the operator compiles or removes one of their own compounds
 * (js/ui/compound.js). `GET /api/drugs` is the single source for both the run
 * picker and the report panel's panel-selector, so both are repainted from the
 * one response rather than each keeping its own idea of what exists.
 */
async function refreshDrugs() {
  try {
    const drugs = await api.getDrugs();
    store.set({ drugs });
    specPanel.setDrugs(drugs);
    reportsPanel.setDrugs(drugs);
    const spec = store.get().spec;
    if (spec) specPanel.loadSpec(spec);   // keep the current selection selected
  } catch (err) {
    notifyError('GET /api/drugs failed after a compound change', err);
  }
}

/**
 * Refetch GET /api/sections after the operator ingested or removed one of their
 * own.
 *
 * `selectId` is passed ONLY on a successful ingest, and moving the selection is
 * then the point: the operator just deliberately loaded that specimen, and
 * landing on it is what makes the picker's not-validated line appear against
 * the section they are about to run. On a REMOVAL it is null and the selection
 * is left alone — repainting a list must not silently move a run onto a
 * different specimen.
 */
async function refreshSections(selectId) {
  try {
    const sections = await api.getSections();
    store.set({ sections });
    specPanel.setSections(sections);
    const spec = store.get().spec;
    const target = selectId && sections.some((s) => s.id === selectId && s.available)
      ? selectId : null;
    if (spec) specPanel.loadSpec(target ? { ...spec, section: target } : spec);
    if (target) paintSectionHeader(target);
  } catch (err) {
    notifyError('GET /api/sections failed after a tissue change', err);
  }
}

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

/**
 * Paint the top bar's run status chip.
 *
 * The top bar is the one readout that is never scrolled away, and for the whole
 * of a run it stated the opposite of the truth. `submitRun` set `#hdr-run` to
 * the new id and nothing touched the chip, which was still the *previous*
 * run's — so from the moment Submit was pressed until the job finished, the
 * header read `RUN r_00077  DONE` while the run history two panels below read
 * `r_00077  RUNNING`. Measured on a 36 s reconstruction: the contradiction was
 * on screen for its entire duration. Worse than merely wrong, it was wrong in
 * the reassuring direction, and the metrics beside it belonged to the previous
 * run, so a user who trusted it would have read the last run's numbers as this
 * one's.
 *
 * `elapsed` is appended while a job is in flight so the header also answers
 * "how long has this been going", which previously required scrolling the
 * Setup rail down to the Active job section.
 */
function paintRunChip(status, elapsedS) {
  const st = $('#hdr-run-status');
  if (!st) return;
  while (st.firstChild) st.removeChild(st.firstChild);
  if (!status) return;
  const live = status === 'queued' || status === 'running';
  st.appendChild(h('span.chip.st-' + status, {
    title: live
      ? 'This run is still on the queue. The metrics on the right belong to the '
        + 'last completed run, not to this one, until it finishes.'
      : '',
  }, status));
  if (live && Number.isFinite(elapsedS)) {
    st.appendChild(h('span.mono.dim', { style: { marginLeft: 'var(--sp-3)' } },
      dur(elapsedS)));
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
  // A field overlay's provenance depends on the RunSpec that produced it (a KDE
  // vascular field is not the same claim as a discrete one), but the legend and
  // the ROI tables do not hold a RunResult. Park it where they can reach it.
  setActiveRun(run);
  inspect.renderRun(run);
  paintRunChip(run.status);
  // No `&& run.error` guard. A failed run with no message is still a failed
  // run, and notifyRunFailed() has a card for exactly that case — the guard
  // meant the one failure that says least about itself announced itself least.
  if (run.status === 'failed') {
    notifyRunFailed(id, run.error);
  }
  paintSectionHeader(run.spec && run.spec.section);

  if (!viewer) return;
  if (run.status !== 'done') {
    clearResultSurfaces(run);
    return;
  }
  await loadGeometry(id);
  await loadMechanism(id, run);
  await loadTimelapse(id, run);
  if (store.get().fieldName) await selectField(store.get().fieldName);
}

/**
 * Take down everything on screen that was DERIVED FROM A RUN, because the run
 * now in the header did not produce it.
 *
 * WHY THIS EXISTS, AND WHY IT IS ONE FUNCTION
 * -------------------------------------------
 * paintRunChip's comment above already names this defect in its earlier form:
 * "the metrics beside it belonged to the previous run, so a user who trusted
 * it would have read the last run's numbers as this one's." That was fixed for
 * the CHIP. It was not fixed for the stage. A run submitted from the Setup
 * rail and watched to failure left the header reading `r_00085 FAILED` over a
 * full stage — the previous run's cells, its mechanism attribution at "11.2%
 * killed", its drug-concentration ramp labelled 0 → 571.1 uM, its time-lapse
 * scrubber at frame 19/19 — none of it produced by r_00085, all of it sitting
 * in the position where r_00085's answer belongs. A screenshot of that travels
 * without the FAILED chip attached to the number, which is the whole failure
 * mode this product exists to not have: the engine is fine and the surface
 * lies.
 *
 * Two paths reached that state and disagreed. `openRun` cleared some of it;
 * `trackRun`'s failure branch cleared none of it. Neither cleared all of it —
 * openRun left the field texture, the field legend with its ticks, the frame
 * stamp, the window/level histogram, the ROI and profile read-outs and the
 * status line, every one of them a number belonging to a run that is not this
 * one. So the two paths do not merely agree now; they agree on a list that was
 * rebuilt from what is actually on screen. A run's terminal state looks the
 * same whether you watched it fail or opened it an hour later.
 *
 * NOT cleared: `store.fieldName`. That is the operator's CHOICE of overlay,
 * not a value produced by a run — openRun's done path re-selects it, which is
 * how opening a completed run keeps showing the field you were looking at. The
 * DATA behind it is dropped here; only the request survives.
 */
function clearResultSurfaces(run) {
  // The time-lapse first: it owns the field texture while it is playing, and a
  // paused scrubber left holding a frame would put it straight back.
  timelapseField = false;
  if (scrubber) scrubber.clear();
  if (overlays) overlays.setFieldFrame(null);

  if (viewer) {
    viewer.clearGeometry();        // also clears the mechanism labels
    viewer.setField(null);
    if (viewer.colorMode === 1) viewer.setColorMode(0);
  }
  store.set({ field: null });

  // The hover/pin card names a CELL INDEX in a geometry that no longer exists
  // — "ENDOTHELIAL #13915 · 23 cells in radius", a per-cell composition read
  // off the previous run's tissue. It survived the first version of this
  // function and was caught only by looking at the screenshot.
  clearPin();
  if (wl) wl.setBrush(null);

  if (mechLegend) mechLegend.set(null);
  if (overlays) overlays.updateLegend(null, null);
  if (inspect) inspect.setFieldInfo(null, null);
  if (wl) { wl.setField(null, null); wl.setRoiHist(null); }
  // Annotations, not just their numbers: a ROI outline drawn on the previous
  // run's tissue is a claim about geometry that is no longer loaded.
  if (tools) tools.clearAll();
  syncColourControl();

  $('#hdr-cells').textContent = NULL_GLYPH;
  const msg = $('#sb-msg');
  if (msg && run) {
    msg.textContent = `run ${run.run_id} ${run.status} · stage cleared: `
      + 'no geometry, no field, no attribution';
  }
  if (overlays) overlays.setEmpty(true, noResultNotice(run));
}

/**
 * What the emptied stage says about itself.
 *
 * Blank is better than stale; blank with a reason is better than blank. Every
 * word here is either ours (procedural) or the run's own `status`, which the
 * server sent. The failure's REASON is not restated — it is the server's
 * sentence, it is long, and it is already rendered in full in two places that
 * are not a canvas overlay. This says where to read it.
 */
function noResultNotice(run) {
  const id = (run && run.run_id) || 'This run';
  const st = (run && run.status) || '';
  if (st === 'failed') {
    return {
      head: 'This run produced nothing to draw',
      body: `${id} failed. Anything that was on this stage was produced by a `
        + 'different run, so it has been cleared rather than left standing '
        + 'under this one’s name. What the server said about the failure is '
        + 'in Result → Metrics, on the right.',
    };
  }
  if (st === 'cancelled') {
    return {
      head: 'This run was stopped before it produced anything',
      body: `${id} was cancelled. Anything that was on this stage was produced `
        + 'by a different run, so it has been cleared rather than left '
        + 'standing under this one’s name.',
    };
  }
  return {
    head: 'This run has not produced anything yet',
    body: `${id} is ${st || 'not finished'}. The stage is empty on purpose: `
      + 'the cells and fields that were here belong to a different run, and '
      + 'this one has produced nothing to replace them with.',
  };
}

/* ------------------------------------------------- mechanism attribution -- */

/**
 * Fetch the per-cell mechanism attribution and hand it to the viewer.
 *
 * A missing attribution is a normal state, not an error: a reconstruction-only
 * run has no kill to attribute. The backend returns 404 with the run's OWN
 * recorded reason, which is surfaced verbatim rather than replaced with a
 * generic message.
 */
async function loadMechanism(id, run) {
  if (!viewer) return;
  viewer.clearLabels();
  if (mechLegend) mechLegend.set(null);
  syncColourControl();

  const msg = $('#sb-msg');
  try {
    const m = await api.getMechanism(id);
    viewer.setLabels(m.labels);
    if (mechLegend) mechLegend.set(m.aggregate);
    syncColourControl();
    const rec = m.aggregate && m.aggregate.reconciliation;
    if (rec && rec.bit_exact === false) {
      // The backend refuses to emit a non-reconciling attribution, so this can
      // only mean the payload was tampered with in transit. Say so loudly.
      toast('err', 'Mechanism attribution does not reconcile',
        `labels give kill_fraction ${rec.kill_fraction_from_labels}, the engine `
        + `reported ${rec.kill_fraction_from_engine}. Do not read this map.`, 0);
    }
    msg.textContent = `mechanism ${int(m.n_cells)} labels · ${bytes(m.bytes)} · `
      + `${m.aggregate.counts.KILLED} killed / ${m.aggregate.n_tumor} tumour`;
  } catch (err) {
    const reason = (err && err.payload && err.payload.detail
      && err.payload.detail.reason)
      || (run && run.mechanism_reason) || (err && err.message);
    msg.textContent = 'no mechanism attribution: ' + reason;
    // A run that simply has no attribution is an ordinary state and the line
    // above is enough. A payload that does not LINE UP with the geometry is
    // not: it would be a map of the wrong cells, so it is escalated rather
    // than left in the status bar. setLabels throws before installing
    // anything, so nothing misaligned is ever on screen when this fires.
    if (err instanceof Error && /refusing to render a misaligned/.test(err.message)) {
      notifyError('Mechanism payload does not match the geometry', err);
    }
  }
}

/* ------------------------------------------------------------- time-lapse -- */

async function loadTimelapse(id, run) {
  if (!scrubber) return;
  scrubber.clear();
  timelapseField = false;
  if (overlays) overlays.setFieldFrame(null);
  store.set({ timelapse: null });
  const info = (run && run.timelapse) || null;
  if (info && info.available === false) return;   // the run says there is none
  try {
    const tl = await api.getTimelapse(id);
    scrubber.set(tl, info || {});
    store.set({ timelapse: tl });
  } catch {
    // A run without a time-lapse is ordinary; the reason already rides in
    // RunResult.timelapse.reason and is shown in the run's provenance panel.
    scrubber.clear();
  }
}

/**
 * Paint one time-lapse frame. The frames are drug CONCENTRATION, which is a
 * different quantity from the /field/drug_auc overlay, so the field is
 * (re)installed as its own named field the first time and thereafter only its
 * pixels are replaced — keeping the operator's window fixed across the series
 * so the animation shows the drug moving, not the ramp re-normalising.
 */
function onTimelapseFrame(frame, i, series) {
  if (!viewer) return;
  const f = viewer.field;
  const same = f && f.name === 'drug_concentration'
    && f.rows === frame.rows && f.cols === frame.cols;
  if (same) {
    viewer.setFieldValues(frame.values, frame.min, frame.max, true);
  } else {
    // Window across the WHOLE series, not this frame, so scrubbing is
    // comparable frame to frame.
    let lo = Infinity, hi = -Infinity;
    for (const fr of series.frames) {
      if (Number.isFinite(fr.min)) lo = Math.min(lo, fr.min);
      if (Number.isFinite(fr.max)) hi = Math.max(hi, fr.max);
    }
    viewer.setField({
      name: 'drug_concentration',
      rows: frame.rows, cols: frame.cols, dx: frame.dx,
      min: Number.isFinite(lo) ? lo : 0,
      max: Number.isFinite(hi) ? hi : 1,
      values: frame.values,
    });
    timelapseField = true;
    store.set({ field: viewer.field });
    if (overlays) overlays.updateLegend('drug_concentration', viewer.field);
    syncColourControl();
  }
  // Stamp the stage legend with the instant it is labelling. The series opens
  // on the frame carrying the most drug, not on t=0, so the key to these
  // pixels has to name the time — otherwise the picture is a penetration
  // claim with no clock on it.
  if (overlays) {
    overlays.setFieldFrame(
      `t = ${fmtClock(frame.time_s)} · frame ${i + 1}/${series.n_frames}`);
  }
  const el = $('#sb-msg');
  if (el) {
    el.textContent = `time-lapse frame ${i + 1}/${series.n_frames} · `
      + `t = ${num(frame.time_s, 0)} s · max ${num(frame.max, 4)} uM`;
  }
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
    // Annotations are section coordinates. A different geometry means they
    // point at nothing, so they go rather than silently re-projecting onto
    // tissue they were not drawn on.
    if (tools) tools.clearAll();
    $('#tb-pin-clear').disabled = true;
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
  // The operator asked for a served overlay: hand the field texture back and
  // stop the time-lapse rather than letting the two fight over it.
  if (timelapseField) { timelapseField = false; if (scrubber) scrubber.pause(); }
  // A served field is not an instant of anything — drop the time stamp with
  // the frame it belonged to, rather than leaving it to label the next field.
  if (overlays) overlays.setFieldFrame(null);
  store.set({ fieldName: name || null });
  if (!viewer) return;
  if (!name) {
    viewer.setField(null);
    store.set({ field: null });
    overlays.updateLegend(null, null);
    inspect.setFieldInfo(null, null);
    if (wl) wl.setField(null, null);
    if (tools) tools.recompute();
    if (viewer.colorMode === 1) viewer.setColorMode(0);
    syncColourControl();
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
    viewer.setField(f);           // installs the full range as the window
    store.set({ field: f });
    overlays.updateLegend(name, f);
    inspect.setFieldInfo(name, f);
    overlays.refreshRamp();
    if (wl) wl.setField(name, f);
    if (tools) tools.recompute();
    syncColourControl();
    msg.textContent = `field ${name} · ${f.rows}×${f.cols} @ ${num(f.dx, 2)} um · ${bytes(f.bytes)}`;
  } catch (err) {
    msg.textContent = `field ${name} failed: ` + (err.message || err);
    notifyError(`GET /api/runs/${id}/field/${name} failed`, err);
    inspect.selectField(null);
    store.set({ fieldName: null });
    viewer.setField(null);
    overlays.updateLegend(null, null);
    if (wl) wl.setField(null, null);
    if (tools) tools.recompute();
    syncColourControl();
  }
}

/**
 * Point the Setup rail's job monitor at a run and own what happens when it
 * settles.
 *
 * Extracted from submitRun because it is no longer only submit that needs it:
 * `boot()` re-attaches to a run that was already in flight when the page
 * loaded (see reattachLiveRun), and a re-attached job must finish EXACTLY the
 * way a submitted one does — same refresh, same openRun, same failure toast.
 * Two copies of this callback would drift, and the branch that would drift
 * unnoticed is the failure branch, which is the one nobody exercises twice.
 */
function trackRun(id) {
  specPanel.track(id, async (doneId, status) => {
    await refreshRuns();
    if (status === 'done') openRun(doneId);
    else if (status === 'failed') {
      // Paint the chip BEFORE handing over, because openRun's own GET can
      // fail: a run that has failed must not be left reading `running`
      // because the request that would have said so failed too.
      paintRunChip('failed');
      // Then the SAME path a user takes when they open the run later. This
      // branch used to be a reimplementation of a subset of openRun — fetch,
      // toast, re-render the inspector — and the part it did not reimplement
      // was the part that takes the previous run's cells, fields and
      // attribution off the screen. So a failure you watched happen and a
      // failure you opened afterwards did not look alike, and the one you
      // watched was the one that lied. See clearResultSurfaces.
      await openRun(doneId);
    } else {
      paintRunChip(status);
      toast('warn', `Run ${doneId} ${status}`, '', 4000);
    }
  });
}

/**
 * Pick up a run that was still in flight when this page loaded, and re-attach
 * the job monitor to it.
 *
 * WHY: the job lives in the SERVER process, not in this tab — `POST /api/runs`
 * returns an id and the work continues whether or not a browser is watching.
 * The monitor did not know that. It was armed only by submitRun, so a reload,
 * a closed laptop lid, or a second tab left the operator with `Active job:
 * idle`, a disabled Cancel button and no progress, while the reconstruction
 * carried on invisibly. On a run measured in minutes — which is the normal
 * case, and the exact case where someone walks away — the product looked like
 * it had lost their work.
 *
 * No client-side bookkeeping is used to find it. The server's own run list is
 * the authority, so this also recovers a run submitted from a DIFFERENT tab or
 * before a browser crash, which a stored id would not. Protocol (TPL) runs are
 * excluded for the same reason the auto-open below excludes them: they are not
 * a single RunSpec reconstruction and this monitor is not their read-out.
 *
 * @returns {object|null} the run it attached to, if any.
 */
function reattachLiveRun() {
  const live = store.get().runs.find((r) => (r.status === 'queued' || r.status === 'running')
    && r.spec && r.spec.kind !== 'protocol');
  if (!live) return null;
  history.setActive(live.run_id);
  store.set({ activeRunId: live.run_id });
  $('#hdr-run').textContent = live.run_id;
  paintRunChip(live.status, NaN);
  trackRun(live.run_id);
  toast('info', `Still running: ${live.run_id}`,
    'This run was already in flight when the page loaded, so the job monitor '
    + 'in the Setup rail has picked it back up. It is running in the server '
    + 'process, not in this tab — that is why closing the tab did not stop it. '
    + 'The metrics on the right belong to the last completed run until this one '
    + 'finishes.', 14000);
  return live;
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
    // Say `queued` NOW rather than leaving the previous run's chip standing
    // until the first poll returns. See paintRunChip.
    paintRunChip('queued', 0);
    trackRun(id);
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
  get projects() { return projectsPanel; },
  get tools() { return tools; },
  get wl() { return wl; },
  api,
  benchmark850k,
  runBenchmark,
  // Load a run into the viewer + inspector from outside the shell. The protocol
  // workspace (js/ui/protocol.js) calls this when a sweep cell is clicked.
  openRun,
};

boot().catch((err) => {
  // eslint-disable-next-line no-console
  console.error('boot failed', err);
  notifyError('Boot failed', err);
});
