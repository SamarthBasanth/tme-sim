// smoke.mjs — headless sanity check of results/studio_export.html.
//
// There is no browser here, so this does not prove the pixels. What it does
// prove, and what would otherwise be untestable before shipping a 9 MB file:
//
//   1. the whole document parses as JS/JSON (the bundler did not corrupt it),
//   2. the baked payload round-trips and its binaries decode as valid
//      TMSG / TMSF with the cell counts and grid the RunResults claim,
//   3. every module in the bundle instantiates against a stub DOM, i.e. the
//      import rewriting produced a working dependency graph,
//   4. the caveats survived verbatim — including the falsification note.
//
// Run: node scripts/studio_export/smoke.mjs [file]

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const FILE = path.resolve(process.argv[2]
  || path.join(HERE, '..', '..', 'results', 'studio_export.html'));
const html = fs.readFileSync(FILE, 'utf8');

let failures = 0;
const ok = (cond, msg) => {
  if (cond) { console.log('  ok   ' + msg); } else { failures++; console.log('  FAIL ' + msg); }
};

/* ---------------------------------------------------------- 1. payload -- */
console.log('payload  <- ' + FILE);
const m = /<script id="tmesim-baked" type="application\/json">([\s\S]*?)<\/script>/.exec(html);
ok(!!m, 'baked JSON block present');
const B = JSON.parse(m[1].replace(/<\\\//g, '</'));
ok(B.runs.length === 3, `3 runs baked (${B.runs.map((r) => r.run_id).join(', ')})`);
ok(B.runs.every((r) => r.status === 'done'), 'every baked run is status=done');
ok(B.sections.length > 0 && B.drugs.length > 0, 'sections + drugs catalogue baked');

const b64 = (s) => Buffer.from(s, 'base64');
const dv = (buf) => new DataView(buf.buffer, buf.byteOffset, buf.byteLength);

for (const r of B.runs) {
  const p = B.payloads[r.run_id];
  const g = b64(B.blobs[p.geometry]);
  const d = dv(g);
  const magic = g.subarray(0, 4).toString('latin1');
  const n = d.getUint32(8, true);
  ok(magic === 'TMSG', `${r.run_id} geometry magic TMSG`);
  ok(n === r.report.n_cells, `${r.run_id} geometry n=${n} matches report.n_cells`);
  ok(g.length === 20 + n * 9, `${r.run_id} geometry length exact for n=${n}`);

  for (const name of ['vascular_density', 'ecm_density', 'oxygen', 'drug_auc']) {
    const f = b64(B.blobs[p[name]]);
    const fd = dv(f);
    const rows = fd.getUint32(8, true), cols = fd.getUint32(12, true);
    const dx = fd.getFloat32(16, true);
    const okMagic = f.subarray(0, 4).toString('latin1') === 'TMSF';
    ok(okMagic && f.length === 28 + rows * cols * 4,
      `${r.run_id} field ${name}: TMSF ${rows}x${cols} @ ${dx} um, length exact`);
    ok(rows === r.report.grid[0] && cols === r.report.grid[1],
      `${r.run_id} field ${name} grid matches report.grid`);
    ok(Number.isFinite(fd.getFloat32(20, true)) && Number.isFinite(fd.getFloat32(24, true)),
      `${r.run_id} field ${name} has a finite min/max header`);
  }
}


/* -------------------------------------------------------- 2b. mechanism -- */
console.log('mechanism attribution');
for (const r of B.runs) {
  const p = B.payloads[r.run_id];
  if (!p.mechanism) { ok(false, `${r.run_id} has NO mechanism payload baked`); continue; }
  const m = b64(B.blobs[p.mechanism]);
  const d = dv(m);
  ok(m.subarray(0, 4).toString('latin1') === 'TMSM', `${r.run_id} mechanism magic TMSM`);
  const n = d.getUint32(8, true);
  const nLabels = d.getUint32(12, true);
  const jsonLen = d.getUint32(16, true);
  ok(m.length === 20 + jsonLen + n, `${r.run_id} TMSM length exact (${n} labels, ${jsonLen} B json)`);

  // 1:1 with the geometry payload's cell order — the whole premise of painting
  // the labels straight into the existing vertex buffer.
  const g = b64(B.blobs[p.geometry]);
  const gn = dv(g).getUint32(8, true);
  ok(n === gn && n === r.report.n_cells,
    `${r.run_id} TMSM n_cells=${n} aligns 1:1 with TMSG n_cells=${gn} and report.n_cells`);

  const agg = JSON.parse(m.subarray(20, 20 + jsonLen).toString('utf8'));
  const labels = m.subarray(20 + jsonLen);
  ok(labels.length === n, `${r.run_id} label array is ${labels.length} bytes, one per cell`);
  ok(nLabels === (agg.label_names || []).length || nLabels > 0,
    `${r.run_id} declares ${nLabels} label codes`);

  // Independently recount the labels here and check the aggregate is not
  // just JSON that says what it likes.
  const NOT_TUMOR = 7, KILLED = 0, ALREADY_DEAD = 6;
  const hist = new Array(8).fill(0);
  for (let i = 0; i < n; i++) hist[labels[i]]++;
  const names = agg.label_names || ['KILLED', 'SURVIVED_UNREACHED', 'SURVIVED_UNDERDOSED',
    'SURVIVED_RESISTANT', 'SURVIVED_STOCHASTIC', 'UNATTRIBUTED', 'ALREADY_DEAD', 'NOT_TUMOR'];
  const mismatch = names.map((nm, i) => [nm, hist[i], (agg.counts || {})[nm] || 0])
    .filter(([, a, b2]) => a !== b2);
  ok(mismatch.length === 0,
    `${r.run_id} recounted label histogram matches aggregate.counts`
    + (mismatch.length ? ' -> ' + JSON.stringify(mismatch) : ''));

  // THE check: rebuild kill_fraction from the labels with the engine's own
  // integer arithmetic and compare with === , not a tolerance.
  const nTumor = n - hist[NOT_TUMOR];
  const nDead = hist[KILLED] + hist[ALREADY_DEAD];
  const rebuilt = nTumor === 0 ? 0 : 1.0 - (nTumor - nDead) / nTumor;
  const engine = r.metrics.kill_fraction;
  const rec = agg.reconciliation || {};
  ok(rec.bit_exact === true, `${r.run_id} aggregate declares bit_exact=${rec.bit_exact}`);
  ok(rebuilt === engine,
    `${r.run_id} kill rebuilt from labels === engine kill_fraction  `
    + `(${rebuilt} vs ${engine}, delta ${Math.abs(rebuilt - engine)})`);
  ok(rec.kill_fraction_from_labels === engine && rec.kill_fraction_from_engine === engine,
    `${r.run_id} recorded reconciliation agrees with RunResult.metrics.kill_fraction`);
}

/* ---------------------------------------------------------- 2. caveats -- */
console.log('provenance');
const CANARY = 'oxygen (pde): FALSIFIED AS VALUE-ADD';
for (const r of B.runs) {
  const cav = r.provenance.caveats || [];
  ok(cav.some((c) => c.startsWith(CANARY)), `${r.run_id} carries "${CANARY}"`);
  ok(cav.some((c) => /NOT FOR CLINICAL USE/.test(c)), `${r.run_id} carries the clinical-use caveat`);
  ok(cav.some((c) => /UPPER BOUND/.test(c)), `${r.run_id} carries the upper-bound caveat`);
  ok(!!r.provenance.ic50_source && r.provenance.alpha > 0, `${r.run_id} has alpha + IC50 provenance`);
}
ok(html.includes(CANARY), 'the falsification caveat is literally in the file');

/* ------------------------------------------------------- 3. the bundle -- */
console.log('bundle');
const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map((x) => x[1]);
ok(scripts.length === 2, `${scripts.length} inline <script> blocks`);
const bundleSrc = scripts[scripts.length - 1];
ok(!/^\s*(import|export)\s/m.test(bundleSrc), 'no leftover ESM import/export in the bundle');

/* A stub DOM just rich enough for the modules to construct. Anything the code
   actually needs and this lacks shows up as a thrown error, which is the
   point. */
function makeEl(tag) {
  const el = {
    tagName: (tag || 'div').toUpperCase(),
    children: [], childNodes: [], style: {}, dataset: {}, nodeType: 1,
    classList: {
      _s: new Set(),
      add(...c) { c.forEach((x) => this._s.add(x)); },
      remove(...c) { c.forEach((x) => this._s.delete(x)); },
      toggle(c, on) { if (on === undefined) { this._s.has(c) ? this._s.delete(c) : this._s.add(c); } else if (on) this._s.add(c); else this._s.delete(c); },
      contains(c) { return this._s.has(c); },
    },
    get firstChild() { return this.children[0] || null; },
    get textContent() { return this._t || ''; },
    set textContent(v) { this._t = String(v); this.children.length = 0; },
    appendChild(c) { this.children.push(c); return c; },
    insertBefore(c) { this.children.unshift(c); return c; },
    removeChild(c) { const i = this.children.indexOf(c); if (i >= 0) this.children.splice(i, 1); return c; },
    remove() {},
    setAttribute(k, v) { this[k] = v; },
    getAttribute(k) { return this[k] === undefined ? null : String(this[k]); },
    removeAttribute(k) { delete this[k]; },
    addEventListener() {}, removeEventListener() {},
    querySelector() { return makeEl('div'); },
    querySelectorAll() { return []; },
    closest() { return makeEl('div'); },
    getBoundingClientRect() { return { left: 0, top: 0, width: 900, height: 600 }; },
    getContext() { return null; },   // forces the "WebGL2 unavailable" branch
    focus() {}, showModal() {}, close() {},
    scrollTop: 0, scrollHeight: 0, value: '', checked: false, disabled: false,
    offsetWidth: 900, offsetHeight: 600, clientWidth: 900, clientHeight: 600,
  };
  return el;
}
const byId = new Map();
const doc = {
  documentElement: makeEl('html'),
  body: makeEl('body'),
  head: makeEl('head'),
  title: '',
  createElement: (t) => makeEl(t),
  createTextNode: (t) => ({ nodeType: 3, textContent: String(t) }),
  createDocumentFragment: () => makeEl('fragment'),
  getElementById: (id) => { if (!byId.has(id)) byId.set(id, makeEl('div')); return byId.get(id); },
  querySelector: () => makeEl('div'),
  querySelectorAll: () => [],
  addEventListener() {},
};
const win = {
  document: doc,
  location: { host: '', href: 'file:///studio_export.html' },
  matchMedia: () => ({ matches: false, addEventListener() {} }),
  localStorage: { getItem: () => null, setItem() {} },
  addEventListener() {},
  requestAnimationFrame: () => 0,
  performance: { now: () => 0 },
  ResizeObserver: class { observe() {} },
  MutationObserver: class { observe() {} disconnect() {} },
  getComputedStyle: () => ({ getPropertyValue: () => '300px' }),
  Blob: class { constructor(a) { this.size = String(a && a[0] || '').length; } },
  navigator: { clipboard: { writeText: async () => {} } },
  setInterval: () => 0, clearInterval: () => {}, setTimeout: () => 0, clearTimeout: () => {},
  atob: (s) => Buffer.from(s, 'base64').toString('latin1'),
  console,
};
win.window = win;
win.self = win;
win.__TMESIM_BAKED = B;

const errors = [];
const origErr = console.error;
console.error = (...a) => errors.push(a.map(String).join(' '));
try {
  const fn = new Function(
    'window', 'document', 'location', 'performance', 'localStorage', 'matchMedia',
    'requestAnimationFrame', 'ResizeObserver', 'MutationObserver', 'getComputedStyle',
    'Blob', 'navigator', 'setInterval', 'clearInterval', 'setTimeout', 'clearTimeout',
    'atob', 'self',
    bundleSrc);
  fn(win, doc, win.location, win.performance, win.localStorage, win.matchMedia,
    win.requestAnimationFrame, win.ResizeObserver, win.MutationObserver,
    win.getComputedStyle, win.Blob, win.navigator, win.setInterval, win.clearInterval,
    win.setTimeout, win.clearTimeout, win.atob, win);
  ok(true, 'bundle evaluates and every module instantiates against a stub DOM');
  ok(!!win.__tmesim && !!win.__tmesim.api, 'window.__tmesim wired (api + store exported)');
  const a = win.__tmesim.api;
  ok(a.api.staticExport === true && a.api.stub === false,
    'api layer reports staticExport=true, stub=false (mock.js is not bundled)');
  const g = a.decodeTMSG(Buffer.from(B.blobs[B.payloads[B.runs[0].run_id].geometry], 'base64').buffer.slice(0));
  ok(g.n === B.runs[0].report.n_cells, `production decodeTMSG returns n=${g.n}`);
} catch (e) {
  ok(false, 'bundle evaluates: ' + e.stack.split('\n').slice(0, 3).join(' | '));
} finally {
  console.error = origErr;
}
if (errors.length) console.log('  note console.error during boot: ' + errors.join(' ;; ').slice(0, 400));

/* ---------------------------------------------------- 4. externalities -- */
console.log('self-containment');
const hrefs = [...html.matchAll(/\b(?:src|href)\s*=\s*"([^"]*)"/g)].map((x) => x[1]);
const external = hrefs.filter((u) => /^(https?:|\/\/)/i.test(u));
ok(external.length === 0, `0 external src/href (${hrefs.length} attrs checked)`);
ok(!/\bfetch\s*\(/.test(html), 'no fetch( anywhere in the document');
ok(!/XMLHttpRequest|EventSource|WebSocket|importScripts/.test(html), 'no XHR / SSE / WebSocket / importScripts');
ok(!/@import|url\((?!['"]?data:)/.test(html), 'no CSS @import and no non-data url()');
const size = fs.statSync(FILE).size;
ok(size < 12 * 1024 * 1024, `size ${(size / 1048576).toFixed(2)} MiB < 12 MiB`);

console.log(failures ? `\n${failures} FAILURE(S)` : '\nall checks passed');
process.exit(failures ? 1 : 0);
