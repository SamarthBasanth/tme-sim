// mock.js — a local stand-in for app/server, used ONLY when the real backend
// is unreachable. It implements the frozen API_CONTRACT surface exactly,
// including the little-endian TMSG / TMSF binary payloads, so the production
// decode path in api.js is the same code in both modes. Nothing here is ever
// presented to the user as a real result: api.js flips a global STUB flag and
// the status bar says so in the clearest terms available.
//
// The numbers below are shaped to be *plausible*, not authoritative. Every
// metric it emits is tagged `assumed` in its provenance block.

/* ------------------------------------------------------------- rng ------ */
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const gauss = (r) => {
  let u = 0, v = 0;
  while (u === 0) u = r();
  while (v === 0) v = r();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
};

/* --------------------------------------------------------- catalogue ---- */

export const SECTIONS = [
  { id: 'xenium_319', label: 'Xenium GBM · 319-panel', n_cells: 40887, extent_um: [8120.0, 5100.0], panel: 'hBrain_254g', available: true },
  { id: 'xenium_gbm07', label: 'Xenium GBM tumor-07 (GSE324268)', n_cells: 162340, extent_um: [10400.0, 7200.0], panel: 'hBrain_266g', available: true },
  { id: 'visium_gbm_a1', label: 'Visium GBM A1 · spot-deconvolved', n_cells: 4992, extent_um: [6500.0, 6500.0], panel: 'whole-tx (deconv)', available: true },
  { id: 'ivygap_ct_1', label: 'Ivy GAP · cellular tumour ROI', n_cells: 18240, extent_um: [4200.0, 3100.0], panel: 'anatomic/ISH', available: true },
  { id: 'phantom_850k', label: 'QA phantom · 850k cells (render stress)', n_cells: 850000, extent_um: [24000.0, 18000.0], panel: 'synthetic', available: true },
  { id: 'htan_gbm_ms', label: 'HTAN GBM multi-section (controlled access)', n_cells: 0, extent_um: [0.0, 0.0], panel: 'CODEX 52-plex', available: false },
];

export const DRUGS = [
  { name: 'temozolomide', mw: 194.15, logp: -2.8, ic50_uM: 606.0, ic50_source: 'GDSC measured (n=51 glioma lines)', bbb_penetrance: 0.35 },
  { name: 'lomustine', mw: 233.7, logp: 2.83, ic50_uM: 28.4, ic50_source: 'GDSC measured (n=38)', bbb_penetrance: 0.72 },
  { name: 'carmustine', mw: 214.05, logp: 1.53, ic50_uM: 41.2, ic50_source: 'GDSC measured (n=29)', bbb_penetrance: 0.68 },
  { name: 'doxorubicin', mw: 543.52, logp: 1.27, ic50_uM: 0.51, ic50_source: 'GDSC measured (n=64)', bbb_penetrance: 0.04 },
  { name: 'cisplatin', mw: 300.05, logp: -2.19, ic50_uM: 8.7, ic50_source: 'GDSC measured (n=57)', bbb_penetrance: 0.12 },
  { name: 'etoposide', mw: 588.56, logp: 0.6, ic50_uM: 2.3, ic50_source: 'GDSC measured (n=61)', bbb_penetrance: 0.05 },
];

export const DEFAULT_SPEC = {
  section: 'xenium_319',
  dx_um: 12.0,
  drug: 'temozolomide',
  sensitivity_mode: 'uniform',
  vascular_mode: 'discrete',
  oxygen_mode: 'pde',
  state_mode: 'oxygen',
  mgmt: null,
  pk: { plasma_peak: 1.0, half_life_s: 14400, n_doses: 1, dose_interval_s: 86400 },
  total_time_s: 600.0,
};

/* -------------------------------------------------- geometry synthesis -- */

const geomCache = new Map();

/**
 * A GBM-shaped point cloud: necrotic centre, dense cellular tumour shell,
 * infiltrative margin, tortuous vessels with endothelium along them, and
 * perivascular immune clustering. Enough structure that the channel toggles
 * and field overlays exercise real visual cases.
 */
export function buildGeometry(sectionId) {
  if (geomCache.has(sectionId)) return geomCache.get(sectionId);
  const sec = SECTIONS.find((s) => s.id === sectionId) || SECTIONS[0];
  const n = sec.n_cells;
  const [EX, EY] = sec.extent_um;
  const r = mulberry32(hash(sectionId));

  // --- vessel skeleton -------------------------------------------------
  const nVessels = Math.max(8, Math.round(Math.sqrt(n) / 6));
  const vx = [], vy = [];
  for (let v = 0; v < nVessels; v++) {
    let x = r() * EX, y = r() * EY;
    let th = r() * Math.PI * 2;
    const steps = 60 + Math.floor(r() * 140);
    const step = Math.min(EX, EY) / 90;
    for (let s = 0; s < steps; s++) {
      th += gauss(r) * 0.22;
      x += Math.cos(th) * step;
      y += Math.sin(th) * step;
      if (x < 0 || x > EX) { th = Math.PI - th; x = Math.max(0, Math.min(EX, x)); }
      if (y < 0 || y > EY) { th = -th; y = Math.max(0, Math.min(EY, y)); }
      vx.push(x); vy.push(y);
    }
  }
  const nv = vx.length;

  const cx = EX * 0.5, cy = EY * 0.5;
  const coreR = Math.min(EX, EY) * 0.30;
  const necR = coreR * 0.34;

  const x = new Float32Array(n), y = new Float32Array(n), type = new Uint8Array(n);
  // TUMOR 0, STROMAL 1, IMMUNE 2, ENDOTHELIAL 3, NECROTIC 4
  for (let i = 0; i < n; i++) {
    const u = r();
    let t, px, py;
    if (u < 0.075) {                 // endothelium hugs the vessel skeleton
      t = 3;
      const k = (r() * nv) | 0;
      px = vx[k] + gauss(r) * 9;
      py = vy[k] + gauss(r) * 9;
    } else if (u < 0.135) {          // necrotic core
      t = 4;
      const a = r() * Math.PI * 2, rad = Math.sqrt(r()) * necR;
      px = cx + Math.cos(a) * rad * 1.25;
      py = cy + Math.sin(a) * rad;
    } else if (u < 0.275) {          // immune: perivascular + diffuse
      t = 2;
      if (r() < 0.62) {
        const k = (r() * nv) | 0;
        px = vx[k] + gauss(r) * 46;
        py = vy[k] + gauss(r) * 46;
      } else { px = r() * EX; py = r() * EY; }
    } else if (u < 0.44) {           // stroma / reactive astrocytes, peripheral
      t = 1;
      const a = r() * Math.PI * 2;
      const rad = coreR * (1.05 + Math.abs(gauss(r)) * 0.9);
      px = cx + Math.cos(a) * rad * 1.4 + gauss(r) * 120;
      py = cy + Math.sin(a) * rad + gauss(r) * 120;
    } else {                         // tumour: dense shell + infiltrative tail
      t = 0;
      const a = r() * Math.PI * 2;
      let rad;
      if (r() < 0.72) rad = necR + Math.abs(gauss(r)) * coreR * 0.55;
      else rad = coreR * (1 + Math.abs(gauss(r)) * 1.5);
      px = cx + Math.cos(a) * rad * 1.32 + gauss(r) * 30;
      py = cy + Math.sin(a) * rad + gauss(r) * 30;
    }
    x[i] = Math.max(0, Math.min(EX, px));
    y[i] = Math.max(0, Math.min(EY, py));
    type[i] = t;
  }

  const g = { n, extent: [EX, EY], x, y, type, vx, vy, cx, cy, coreR, necR };
  geomCache.set(sectionId, g);
  return g;
}

function hash(s) {
  let h = 2166136261;
  for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
  return h >>> 0;
}

/* ------------------------------------------------------ binary encode --- */

// Header per contract: magic(4) | version u32 | n_cells u32 | extent_x f32 |
// extent_y f32  ==  20 bytes, then x f32[n], y f32[n], type u8[n].
export function encodeTMSG(g) {
  const n = g.n;
  const header = 20;
  const buf = new ArrayBuffer(header + n * 4 * 2 + n);
  const dv = new DataView(buf);
  const u8 = new Uint8Array(buf);
  u8[0] = 0x54; u8[1] = 0x4D; u8[2] = 0x53; u8[3] = 0x47;
  dv.setUint32(4, 1, true);
  dv.setUint32(8, n, true);
  dv.setFloat32(12, g.extent[0], true);
  dv.setFloat32(16, g.extent[1], true);
  new Float32Array(buf, header, n).set(g.x);
  new Float32Array(buf, header + n * 4, n).set(g.y);
  u8.set(g.type, header + n * 8);
  return buf;
}

export function encodeTMSF(f) {
  const nvals = f.rows * f.cols;
  const header = 4 + 4 + 4 + 4 + 4 + 4 + 4;             // 28 B
  const buf = new ArrayBuffer(header + nvals * 4);
  const dv = new DataView(buf);
  const u8 = new Uint8Array(buf);
  u8[0] = 0x54; u8[1] = 0x4D; u8[2] = 0x53; u8[3] = 0x46;
  dv.setUint32(4, 1, true);
  dv.setUint32(8, f.rows, true);
  dv.setUint32(12, f.cols, true);
  dv.setFloat32(16, f.dx, true);
  dv.setFloat32(20, f.min, true);
  dv.setFloat32(24, f.max, true);
  new Float32Array(buf, header, nvals).set(f.values);
  return buf;
}

/* ----------------------------------------------------------- fields ---- */

const fieldCache = new Map();

export function buildField(sectionId, name, dx) {
  const key = `${sectionId}|${name}|${dx}`;
  if (fieldCache.has(key)) return fieldCache.get(key);
  const g = buildGeometry(sectionId);
  const [EX, EY] = g.extent;
  const cols = Math.max(8, Math.min(768, Math.round(EX / dx)));
  const rows = Math.max(8, Math.min(768, Math.round(EY / dx)));
  const sx = cols / EX, sy = rows / EY;
  const vals = new Float32Array(rows * cols);

  // vascular density: stamp a small kernel around every vessel node
  const vasc = new Float32Array(rows * cols);
  const rad = Math.max(1, Math.round(70 * sx));
  for (let k = 0; k < g.vx.length; k++) {
    const ci = Math.round(g.vx[k] * sx), ri = Math.round(g.vy[k] * sy);
    for (let dr = -rad; dr <= rad; dr++) {
      const rr = ri + dr; if (rr < 0 || rr >= rows) continue;
      for (let dc = -rad; dc <= rad; dc++) {
        const cc = ci + dc; if (cc < 0 || cc >= cols) continue;
        const d2 = (dr * dr + dc * dc) / (rad * rad);
        if (d2 > 1) continue;
        vasc[rr * cols + cc] += (1 - d2) * 0.16;
      }
    }
  }
  let vmax = 0; for (let i = 0; i < vasc.length; i++) if (vasc[i] > vmax) vmax = vasc[i];
  if (vmax > 0) for (let i = 0; i < vasc.length; i++) vasc[i] = Math.min(1, vasc[i] / vmax);

  const cxg = g.cx * sx, cyg = g.cy * sy, coreG = g.coreR * sx;

  for (let rI = 0; rI < rows; rI++) {
    for (let cI = 0; cI < cols; cI++) {
      const i = rI * cols + cI;
      const dxx = (cI - cxg) / 1.32, dyy = (rI - cyg);
      const dist = Math.sqrt(dxx * dxx + dyy * dyy) / Math.max(1e-6, coreG);
      const v = vasc[i];
      switch (name) {
        case 'vascular_density':
          vals[i] = v * (dist < 1 ? 0.45 + 0.55 * dist : 1);   // sparse in core
          break;
        case 'ecm_density':
          vals[i] = Math.max(0, Math.min(1, 0.22 + 0.5 * Math.exp(-Math.pow(dist - 1.15, 2) / 0.32) + 0.12 * v));
          break;
        case 'oxygen': {
          const supply = v * (dist < 1 ? 0.45 + 0.55 * dist : 1);
          const demand = Math.exp(-Math.pow(dist, 2) / 1.6);
          vals[i] = Math.max(0.005, Math.min(1, 0.06 + 1.35 * supply - 0.55 * demand));
          break;
        }
        case 'drug_auc': {
          const supply = v * (dist < 1 ? 0.45 + 0.55 * dist : 1);
          const ecm = 0.22 + 0.5 * Math.exp(-Math.pow(dist - 1.15, 2) / 0.32);
          vals[i] = 3.1e6 * (0.08 + supply) * (1.25 - 0.65 * ecm);
          break;
        }
        default:
          vals[i] = v;
      }
    }
  }
  let mn = Infinity, mx = -Infinity;
  for (let i = 0; i < vals.length; i++) { if (vals[i] < mn) mn = vals[i]; if (vals[i] > mx) mx = vals[i]; }
  const f = { name, rows, cols, dx: EX / cols, min: mn, max: mx, values: vals };
  fieldCache.set(key, f);
  return f;
}

/* -------------------------------------------------------- validation --- */

const SPEC_KEYS = ['section', 'dx_um', 'drug', 'sensitivity_mode', 'vascular_mode',
  'oxygen_mode', 'state_mode', 'mgmt', 'pk', 'total_time_s'];
const PK_KEYS = ['plasma_peak', 'half_life_s', 'n_doses', 'dose_interval_s'];

export function validateSpec(spec) {
  const errors = [], warnings = [];
  const E = (path, msg) => errors.push({ path, message: msg });
  const W = (path, msg) => warnings.push({ path, message: msg });

  if (spec === null || typeof spec !== 'object' || Array.isArray(spec)) {
    return { ok: false, errors: [{ path: '', message: 'RunSpec must be a JSON object' }], warnings: [] };
  }
  for (const k of Object.keys(spec)) {
    if (!SPEC_KEYS.includes(k)) E(k, `unknown field "${k}" — not in RunSpec`);
  }
  for (const k of SPEC_KEYS) {
    if (!(k in spec)) E(k, `missing required field "${k}"`);
  }

  const sec = SECTIONS.find((s) => s.id === spec.section);
  if (spec.section !== undefined) {
    if (typeof spec.section !== 'string') E('section', 'section must be a string id');
    else if (!sec) E('section', `unknown section "${spec.section}" — see GET /api/sections`);
    else if (!sec.available) E('section', `section "${spec.section}" is not available in this deployment (controlled access)`);
  }

  if (spec.dx_um !== undefined) {
    if (typeof spec.dx_um !== 'number' || !Number.isFinite(spec.dx_um)) E('dx_um', 'dx_um must be a finite number (um)');
    else if (spec.dx_um < 2) E('dx_um', 'dx_um < 2 um — grid would exceed memory budget');
    else if (spec.dx_um > 100) E('dx_um', 'dx_um > 100 um — coarser than the cell pitch, transport is meaningless');
    else if (spec.dx_um < 6) W('dx_um', 'dx_um < 6 um: reconstruction cost grows as dx^-2; expect minutes');
    else if (spec.dx_um > 25) W('dx_um', 'dx_um > 25 um under-resolves the perivascular gradient');
  }

  if (spec.drug !== undefined && spec.drug !== null) {
    if (typeof spec.drug !== 'string') E('drug', 'drug must be a string name or null');
    else if (!DRUGS.some((d) => d.name === spec.drug)) E('drug', `unknown drug "${spec.drug}" — see GET /api/drugs`);
  }

  const enums = {
    sensitivity_mode: ['uniform', 'determinant'],
    vascular_mode: ['discrete', 'kde'],
    oxygen_mode: ['pde', 'vessel_distance', 'markers'],
    state_mode: ['oxygen', 'vessel_distance', 'markers'],
  };
  for (const [k, allowed] of Object.entries(enums)) {
    if (spec[k] !== undefined && !allowed.includes(spec[k])) {
      E(k, `${k} must be one of ${allowed.map((v) => `"${v}"`).join(' | ')}`);
    }
  }
  if (spec.vascular_mode === 'kde') W('vascular_mode', 'kde vessel field did not converge under grid refinement; discrete is the production default');
  if (spec.oxygen_mode === 'markers') W('oxygen_mode', 'marker-derived oxygen is a SURROGATE, not a measured pO2 — every downstream number inherits that');
  if (spec.state_mode && spec.state_mode !== 'oxygen') {
    W('state_mode', `state_mode=${spec.state_mode} relabels cells.state on the RECONSTRUCTION only; sim.engine re-bands from its own oxygen solve, so it moves the reported band occupancies, not kill_fraction`);
  }
  if (spec.state_mode === 'markers' && spec.section === 'xenium_319') {
    W('state_mode', 'state_mode=markers on the 319-gene panel: it carries NO hypoxia-response transcript, so HYPOXIC is UNMEASURED and assigned to no cell — that is NOT "computed and found to be zero". Proliferation detection is 13x below the published Ki-67 range, so the per-cell call falls back to a neighbourhood ENRICHMENT rule whose level is not an index.');
  }
  if (spec.section === 'xenium_319' && spec.drug) {
    W('state_mode', 'ABSOLUTE KILL SHOULD NOT BE QUOTED FROM THIS SECTION UNDER ANY state_mode: no state route clears the bar on this panel and kill moves up to 0.3252 absolute (1.77x) across routes. What is defensible is the MAP: within this one compound on this one section, the RELATIVE exposure between regions, because a per-compound plasma factor is an exact linear rescale of the whole field and cancels from a ratio between points. The absolute level does not survive, and neither does a comparison against another compound.');
  }
  if (spec.sensitivity_mode === 'determinant' && spec.drug && !['temozolomide', 'lomustine', 'carmustine'].includes(spec.drug)) {
    W('sensitivity_mode', `no determinant gene on panel for "${spec.drug}"; will fall back to uniform 1.0`);
  }

  if (spec.mgmt !== undefined && spec.mgmt !== null) {
    if (typeof spec.mgmt !== 'number') E('mgmt', 'mgmt must be a number (log2 CPM) or null');
    else if (spec.mgmt < -2 || spec.mgmt > 12) E('mgmt', 'mgmt outside the calibrated range [-2, 12] log2 CPM');
    if (spec.drug !== 'temozolomide') W('mgmt', 'MGMT modulates alkylator response only; it is ignored for this drug');
  }

  if (spec.pk !== undefined) {
    if (spec.pk === null || typeof spec.pk !== 'object') E('pk', 'pk must be an object');
    else {
      for (const k of Object.keys(spec.pk)) if (!PK_KEYS.includes(k)) E(`pk.${k}`, `unknown pk field "${k}"`);
      for (const k of PK_KEYS) if (!(k in spec.pk)) E(`pk.${k}`, `missing pk field "${k}"`);
      if (typeof spec.pk.plasma_peak === 'number' && spec.pk.plasma_peak <= 0) E('pk.plasma_peak', 'plasma_peak must be > 0 uM');
      if (typeof spec.pk.half_life_s === 'number' && spec.pk.half_life_s <= 0) E('pk.half_life_s', 'half_life_s must be > 0 s');
      if (typeof spec.pk.n_doses === 'number' && (spec.pk.n_doses < 1 || spec.pk.n_doses % 1 !== 0)) E('pk.n_doses', 'n_doses must be an integer >= 1');
      if (typeof spec.pk.n_doses === 'number' && spec.pk.n_doses > 30) E('pk.n_doses', 'n_doses > 30 exceeds the scheduling budget');
      if (typeof spec.pk.dose_interval_s === 'number' && spec.pk.dose_interval_s <= 0) E('pk.dose_interval_s', 'dose_interval_s must be > 0 s');
    }
  }

  if (spec.total_time_s !== undefined) {
    if (typeof spec.total_time_s !== 'number' || !Number.isFinite(spec.total_time_s)) E('total_time_s', 'total_time_s must be a finite number (s)');
    else if (spec.total_time_s <= 0) E('total_time_s', 'total_time_s must be > 0 s');
    else if (spec.total_time_s > 2.6e6) E('total_time_s', 'total_time_s > 30 d exceeds the integration budget');
    else if (spec.pk && typeof spec.pk.half_life_s === 'number' && spec.total_time_s < 2 * spec.pk.half_life_s) {
      W('total_time_s', 'integration window shorter than 2 plasma half-lives: AUC will be truncated and kill_fraction understated');
    }
  }
  if (spec.drug === null) W('drug', 'reconstruction-only run: no PK/PD stage, drug metrics will be null');

  return { ok: errors.length === 0, errors, warnings };
}

/* ------------------------------------------- verdict stability --------- */
//
// A MIRROR OF app/server/jobs.py::verdict_stability, AND IT HAS TO BE ONE.
//
// THE FINDING prints, directly under its verdict sentence, whether that
// verdict survives the transvascular exchange term the transport law does not
// have. When extra_metrics.verdict_stability is absent the panel falls back to
// "this run's compound was not one of the twelve scored" -- which on the
// hosted demo told every visitor that TEMOZOLOMIDE, the compound the Studio
// opens on and one the census scored as REVERSING, had not been looked at.
// The panel's most important caveat read as an absence of evidence when the
// evidence exists and says the opposite.
//
// Three tiers, same as the server: the delivery/dose verdict inverts (5 of
// 12), it holds while the four-way dominant limiter inverts (3 of 12), or both
// hold (4 of 12). tests/test_verdict_stability.py asserts THESE lists and the
// Python ones against results/verdict_inversion_census_data.json, so the copy
// here cannot drift from the copy there or from the run that measured either.
const VERDICT_REVERSING = ['carmustine', 'cisplatin', 'lomustine', 'procarbazine', 'temozolomide'];
const VERDICT_LIMITER_ONLY = ['irinotecan', 'palbociclib', 'ribociclib'];
const VERDICT_STABLE = ['doxorubicin', 'etoposide', 'methotrexate', 'vincristine'];
const VERDICT_SECTION = 'xenium_319';

function verdictStability(drugName) {
  if (!drugName) return null;
  const key = String(drugName).trim().toLowerCase();
  const reverses = VERDICT_REVERSING.includes(key);
  const limiterOnly = VERDICT_LIMITER_ONLY.includes(key);
  const holds = VERDICT_STABLE.includes(key);
  // Absent from all three means the census did not score it, and absence is
  // not evidence of stability — so it stays null and the panel says so.
  if (!reverses && !limiterOnly && !holds) return null;
  const verdict = reverses ? 'reverses' : limiterOnly ? 'holds_with_caveat' : 'holds';
  const note = reverses
    ? `${key} REVERSES its delivery/dose verdict when the missing `
      + 'transvascular exchange term is switched on, so the verdict above is '
      + 'one of the unstable ones.'
    : limiterOnly
      ? `${key} KEEPS its delivery/dose verdict when the missing `
        + 'transvascular exchange term is switched on -- but its four-way '
        + 'dominant limiter does NOT, so the label under the sentence changes '
        + 'even though the sentence does not.'
      : `${key} KEEPS its delivery/dose verdict when the missing `
        + 'transvascular exchange term is switched on, and keeps most of its '
        + 'exposure under the split.';
  return {
    compound: key,
    stable: verdict === 'holds',
    verdict,
    delivery_dose_verdict_inverts: reverses,
    dominant_limiter_inverts: reverses || limiterOnly,
    measured_on: VERDICT_SECTION,
    source: 'results/verdict_inversion_census.md',
    note,
  };
}

/* ------------------------------------------- mechanism attribution ----- */
//
// WHY THE STUB HAS TO CARRY THIS AT ALL
// -------------------------------------
// THE FINDING (js/ui/finding.js) leads the READ rail and states the product's
// central claim. It reads `run.mechanism`, and this file did not synthesise
// one -- so on every hosted instance whose backend is asleep, and on every
// offline read, the product's headline panel rendered as an apology and told
// the reader it was "a backend omission". There is no backend in stub mode.
// The panel, its click-to-isolate, the legend's SHOW ALL and the verdict's own
// reversal caveat were all unreachable on the demo that most people see.
//
// The reach block one section down already carries the same argument in its
// own comment -- "a stub that omits it would let the panel's empty state ship
// untested" -- and THE FINDING sits above reach in the rail.
//
// WHAT MAKES IT HONEST RATHER THAN JUST PRESENT
// ---------------------------------------------
// 1. THE LABEL ARRAY IS THE SOURCE OF TRUTH. Labels are assigned per cell
//    first; every count, every share, and kill_fraction itself are then TALLIED
//    off that array. Nothing is planned and then hoped for, so the legend's
//    reconciliation cannot fail the way it can when a float is picked and
//    integers are fitted to it afterwards. Same rule the reach ladder already
//    follows ("the headline is PICKED OFF the ladder").
// 2. LEVEL AND SHAPE COME FROM THE RUN, NOT FROM A DIE. The exposure falls off
//    exponentially with distance from the vessel skeleton at this run's OWN
//    penetration_length_um; its level is solved so the share at or above IC50
//    equals therapeutic_coverage_1x, and the kill rate so the mean kill
//    probability equals kill_fraction. So the delivery/dose split moves with
//    the COMPOUND -- a short-penetration compound leaves UNREACHED cells out in
//    the infiltrative margin, a freely-penetrating one leaves UNDERDOSED cells
//    everywhere -- rather than being one fixed picture wearing different
//    labels, and the two metrics printed beside the breakdown are readings of
//    the same label array rather than independent draws.
// 3. IT IS STILL SYNTHETIC AND SAYS SO. js/ui/finding.js and js/ui/mechanism.js
//    both check isSynthetic() and mark the panel; the reconciliation's own
//    `method` string says that in stub mode both sides of the comparison are
//    the same arithmetic and it audits nothing.
//
// Codes, band thresholds and wording are mirrored from tmesim/sim/mechanism.py
// and tmesim/core/config.py; that module is the authority and this is a
// stand-in that has to stay contract-shaped.

const MECH_LABELS = [
  'KILLED', 'SURVIVED_UNREACHED', 'SURVIVED_UNDERDOSED', 'SURVIVED_RESISTANT',
  'SURVIVED_STOCHASTIC', 'UNATTRIBUTED', 'ALREADY_DEAD', 'NOT_TUMOR',
];
const M_KILLED = 0, M_UNREACHED = 1, M_UNDERDOSED = 2, M_RESISTANT = 3,
  M_STOCHASTIC = 4, M_UNATTRIBUTED = 5, M_NOT_TUMOR = 7;

const MECH_LIMITER = {
  SURVIVED_UNREACHED: 'delivery',
  SURVIVED_UNDERDOSED: 'dose',
  SURVIVED_RESISTANT: 'resistance',
  SURVIVED_STOCHASTIC: 'chance',
  UNATTRIBUTED: 'unattributed',
};

// Verbatim from tmesim.sim.mechanism.LABEL_DESCRIPTION: the legend shows these
// as row tooltips and a stub that paraphrased them would teach the wrong
// vocabulary.
const MECH_DESCRIPTION = {
  KILLED: "killed by the drug (the engine's own Bernoulli outcome)",
  SURVIVED_UNREACHED: 'DELIVERY-limited: local exposure E/IC50 below the '
    + 'unreached band — the drug never arrived in a pharmacologically '
    + 'meaningful amount',
  SURVIVED_UNDERDOSED: 'DOSE-limited: drug arrived but E stayed below IC50, so '
    + 'the Hill term never reached half of its ceiling',
  SURVIVED_RESISTANT: 'RESISTANCE-limited: exposure was adequate, but '
    + "sensitivity x state suppressed the kill probability by at least the "
    + "model's own mildest resistance penalty",
  SURVIVED_STOCHASTIC: 'CHANCE: exposure and biology both adequate; the cell '
    + 'drew above its kill probability',
  UNATTRIBUTED: 'not explained by the bands above — reported, never reassigned',
  ALREADY_DEAD: 'tumour cell already DEAD before the drug was applied; counted '
    + 'in kill_fraction but not killed by the drug',
  NOT_TUMOR: 'non-tumour cell; kill_fraction is a tumour-only quantity',
};

// hill 1.2 is the shipped per-compound value in tmesim/core/config.py; the two
// resistance floors are PDConfig's, and resistant_modulation resolves to
// max(hypoxia 0.35, quiescent 0.5). stochastic_min_p_kill is then
// f(adequate_ratio) * resistant_modulation = 0.5 * 0.5 at hill 1.2.
const MECH_HILL = 1.2;
const MECH_THRESHOLDS = {
  unreached_ratio: 0.1,
  adequate_ratio: 1.0,
  resistant_modulation: 0.5,
};

const hillP = (ratio) => {
  if (!(ratio > 0)) return 0;
  const r = Math.pow(ratio, MECH_HILL);
  return r / (1 + r);
};

/* Distance from every grid cell to the nearest vessel node, in um.
 *
 * Two-pass chamfer over the same grid buildField uses, seeded at the vessel
 * skeleton, so the attribution and the AUC overlay are shapes of the SAME
 * vasculature. A per-cell scan against every vessel node would be O(n * nv) and
 * the phantom section is 850k cells. */
const distCache = new Map();

function vesselDistance(sectionId, dx) {
  const key = `${sectionId}|${dx}`;
  if (distCache.has(key)) return distCache.get(key);
  const g = buildGeometry(sectionId);
  const [EX, EY] = g.extent;
  const cols = Math.max(8, Math.min(768, Math.round(EX / dx)));
  const rows = Math.max(8, Math.min(768, Math.round(EY / dx)));
  const cw = EX / cols, ch = EY / rows;            // cell size in um
  const BIG = 1e9;
  const d = new Float32Array(rows * cols).fill(BIG);
  for (let k = 0; k < g.vx.length; k++) {
    const ci = Math.min(cols - 1, Math.max(0, Math.round(g.vx[k] / cw)));
    const ri = Math.min(rows - 1, Math.max(0, Math.round(g.vy[k] / ch)));
    d[ri * cols + ci] = 0;
  }
  const diag = Math.hypot(cw, ch);
  const relax = (i, j, w) => { const v = d[j] + w; if (v < d[i]) d[i] = v; };
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      const i = r * cols + c;
      if (r > 0) relax(i, i - cols, ch);
      if (c > 0) relax(i, i - 1, cw);
      if (r > 0 && c > 0) relax(i, i - cols - 1, diag);
      if (r > 0 && c < cols - 1) relax(i, i - cols + 1, diag);
    }
  }
  for (let r = rows - 1; r >= 0; r--) {
    for (let c = cols - 1; c >= 0; c--) {
      const i = r * cols + c;
      if (r < rows - 1) relax(i, i + cols, ch);
      if (c < cols - 1) relax(i, i + 1, cw);
      if (r < rows - 1 && c < cols - 1) relax(i, i + cols + 1, diag);
      if (r < rows - 1 && c > 0) relax(i, i + cols - 1, diag);
    }
  }
  const out = { rows, cols, cw, ch, dist: d };
  distCache.set(key, out);
  return out;
}

/**
 * Per-cell labels plus the aggregate block, in the shape
 * tmesim.sim.mechanism.attribute() returns (minus `labels`, which travels as
 * the binary half of TMSM).
 *
 * @returns {null|{labels: Uint8Array, aggregate: object, kill_fraction: number}}
 */
function mechanismFor(spec, rng, targets) {
  if (!spec.drug) return null;
  const g = buildGeometry(spec.section);
  const vd = vesselDistance(spec.section, spec.dx_um);
  const d = DRUGS.find((x) => x.name === spec.drug) || DRUGS[0];

  // TWO SCALARS, TWO TARGETS, AND NEITHER IS INVENTED HERE.
  //
  //   SHAPE  the exposure falls off exponentially with distance from the
  //          vessel skeleton at this run's own penetration_length_um. That is
  //          what makes the delivery/dose split a property of the COMPOUND: a
  //          poorly-penetrating one leaves UNREACHED cells out in the
  //          infiltrative margin, a freely-penetrating one leaves UNDERDOSED
  //          cells everywhere.
  //   LEVEL  solved so the share of tumour cells at or above IC50 equals
  //          therapeutic_coverage_1x -- the number the METRICS block prints
  //          two lines under this breakdown, and whose own definition is
  //          "fraction of tumour cells whose exposure reaches 1x IC50". A
  //          panel that showed a coverage figure its own attribution
  //          contradicted would be self-refuting on one screen.
  //   RATE   solved so the mean kill probability equals kill_fraction, the
  //          engine's max_kill_rate standing in for everything between an
  //          exposure and a dead cell that this stub does not model.
  //
  // Both solves are bisections on a monotone quantity, so neither needs a
  // derivative and neither can wander: the level only ever raises coverage,
  // the rate only ever raises kill.
  const L = Math.max(20, targets.penetration_length_um || 120);
  const coverTarget = Math.max(0, Math.min(0.999, targets.therapeutic_coverage_1x || 0));
  const killTarget = Math.max(0, Math.min(0.999, targets.kill_fraction || 0));

  const [EX, EY] = g.extent;
  const sx = vd.cols / EX, sy = vd.rows / EY;
  const falloff = (i2) => {
    const c = Math.min(vd.cols - 1, Math.max(0, Math.round(g.x[i2] * sx)));
    const r = Math.min(vd.rows - 1, Math.max(0, Math.round(g.y[i2] * sy)));
    return Math.exp(-vd.dist[r * vd.cols + c] / L);
  };

  const tumour = [];
  for (let i2 = 0; i2 < g.n; i2++) if (g.type[i2] === 0) tumour.push(i2);
  // Both solves are for ONE scalar each, so they run on a subsample: the
  // phantom section is 850k cells and a 60-step bisection over all of them
  // buys no precision the tally would keep.
  const stride = Math.max(1, Math.ceil(tumour.length / 20000));
  const sample = [];
  for (let k = 0; k < tumour.length; k += stride) sample.push(falloff(tumour[k]));

  const bisect = (fn, want) => {
    let lo = -8, hi = 8;
    for (let it = 0; it < 60; it++) {
      const mid = (lo + hi) / 2;
      if (fn(Math.pow(10, mid)) < want) lo = mid; else hi = mid;
    }
    return Math.pow(10, (lo + hi) / 2);
  };
  const coverAt = (lv) => {
    if (!sample.length) return 0;
    let n2 = 0;
    for (let k = 0; k < sample.length; k++) if (lv * sample[k] >= MECH_THRESHOLDS.adequate_ratio) n2 += 1;
    return n2 / sample.length;
  };
  const level = coverTarget > 0 ? bisect(coverAt, coverTarget) : Math.pow(10, -8);

  const killAt = (rate) => {
    if (!sample.length) return 0;
    let acc = 0;
    for (let k = 0; k < sample.length; k++) acc += hillP(level * sample[k]) * Math.min(1, rate);
    return acc / sample.length;
  };
  const maxKillRate = Math.min(1, bisect(killAt, killTarget));

  // DERIVED from this run, exactly as tmesim.sim.mechanism.MechanismThresholds
  // .resolve() derives it: f(adequate_ratio) * resistant_modulation, times the
  // run's own kill rate. Hardcoding the max_kill_rate = 1 value (0.25) would
  // put every adequately-exposed survivor of a low-rate run into UNATTRIBUTED,
  // which the legend paints as an alarm — and it would be the stub's own
  // arithmetic setting the alarm off, not anything about the run.
  const stochasticFloor = hillP(MECH_THRESHOLDS.adequate_ratio)
    * MECH_THRESHOLDS.resistant_modulation * maxKillRate;

  // sensitivity x state. `uniform` sets every cell to 1.0 — which is why no
  // cell can be RESISTANCE-limited under it, exactly as in the real rule
  // (modulation must be BELOW 1.0, not merely at or under the threshold).
  const uniform = spec.sensitivity_mode === 'uniform';

  /* --- the labels, then the tally ------------------------------------- */
  const labels = new Uint8Array(g.n).fill(M_NOT_TUMOR);
  const ratios = new Float64Array(tumour.length);
  const mods = new Float64Array(tumour.length);
  let nAdequate = 0;
  for (let k = 0; k < tumour.length; k++) {
    const i2 = tumour[k];
    const ratio = level * falloff(i2);
    ratios[k] = ratio;
    if (ratio >= MECH_THRESHOLDS.adequate_ratio) nAdequate += 1;
    const modulation = uniform ? 1.0 : Math.max(0.08, Math.min(1, 0.791 + gauss(rng) * 0.22));
    mods[k] = modulation;
    const p = hillP(ratio) * modulation * maxKillRate;
    if (rng() < p) { labels[i2] = M_KILLED; continue; }
    if (ratio < MECH_THRESHOLDS.unreached_ratio) labels[i2] = M_UNREACHED;
    else if (ratio < MECH_THRESHOLDS.adequate_ratio) labels[i2] = M_UNDERDOSED;
    else if (modulation < 1 && modulation <= MECH_THRESHOLDS.resistant_modulation) labels[i2] = M_RESISTANT;
    else if (p >= stochasticFloor) labels[i2] = M_STOCHASTIC;
    else labels[i2] = M_UNATTRIBUTED;
  }

  const counts = {};
  MECH_LABELS.forEach((nm) => { counts[nm] = 0; });
  for (let i = 0; i < labels.length; i++) counts[MECH_LABELS[labels[i]]] += 1;

  // The engine's own arithmetic, on the tally above, in the same order:
  // integer counts, one true division, one subtraction from 1.0. The value
  // that ships as metrics.kill_fraction IS this one, so the reconciliation is
  // an identity here rather than a check — see `method` below.
  const nTumor = counts.KILLED + counts.SURVIVED_UNREACHED
    + counts.SURVIVED_UNDERDOSED + counts.SURVIVED_RESISTANT
    + counts.SURVIVED_STOCHASTIC + counts.UNATTRIBUTED + counts.ALREADY_DEAD;
  const nDead = counts.KILLED + counts.ALREADY_DEAD;
  const viable = nTumor === 0 ? 0.0 : (nTumor - nDead) / nTumor;
  const killFraction = 1.0 - viable;

  const nSurvived = nTumor - nDead;
  const survival = {};
  Object.entries(MECH_LIMITER).forEach(([nm, lim]) => {
    survival[lim] = nSurvived ? counts[nm] / nSurvived : 0.0;
  });

  const sorted = Float64Array.from(ratios).sort();
  const q = (f) => (sorted.length
    ? sorted[Math.min(sorted.length - 1, Math.max(0, Math.round(f * (sorted.length - 1))))]
    : null);
  const median = q(0.5);

  // Exposure-limited survivals that would ALSO have survived a dose increase.
  // Zero under `uniform` by construction, because nothing is suppressed there.
  let joint = 0, expLimited = 0;
  for (let k = 0; k < tumour.length; k++) {
    const lab = labels[tumour[k]];
    if (lab !== M_UNREACHED && lab !== M_UNDERDOSED) continue;
    expLimited += 1;
    // The same predicate the RESISTANT band uses, applied to a cell the
    // exposure band already claimed: below 1.0 AND at or under the threshold.
    if (mods[k] < 1 && mods[k] <= MECH_THRESHOLDS.resistant_modulation) joint += 1;
  }

  const aggregate = {
    label_names: MECH_LABELS.slice(),
    label_descriptions: { ...MECH_DESCRIPTION },
    counts,
    n_cells: g.n,
    n_tumor: nTumor,
    n_survived: nSurvived,
    killed_fraction_of_tumor: nTumor ? counts.KILLED / nTumor : 0.0,
    survival_breakdown: survival,
    survival_breakdown_note: 'fractions of SURVIVING tumour cells, so they sum '
      + "to 1. 'delivery' = SURVIVED_UNREACHED, 'dose' = SURVIVED_UNDERDOSED, "
      + "'resistance' = SURVIVED_RESISTANT, 'chance' = SURVIVED_STOCHASTIC.",
    joint_limited: {
      count: joint,
      fraction_of_exposure_limited: expLimited ? joint / expLimited : 0.0,
      note: 'exposure-limited survivals that are ALSO at or below the '
        + 'resistance threshold. They are labelled by their exposure band '
        + '(documented precedence: below IC50 the exposure is the binding '
        + 'constraint), but raising the dose alone would not have killed them. '
        + 'Reported so the delivery/dose share is not read as a clean '
        + 'counterfactual.',
    },
    exposure_ratio_percentiles_tumor: {
      p05: q(0.05), p50: median, p95: q(0.95),
      max: sorted.length ? sorted[sorted.length - 1] : null,
    },
    exposure_shortfall_fold: median > 0 ? MECH_THRESHOLDS.adequate_ratio / median : null,
    exposure_shortfall_note: 'adequate_ratio / median(E/IC50) over tumour '
      + 'cells: the fold change in exposure the median tumour cell would need '
      + 'to reach IC50. Above 1 the median cell is short of drug; below 1 it '
      + 'is already past the half-effect point and something else is holding '
      + 'the kill down.',
    exposure_axis_note: 'UNREACHED and UNDERDOSED are two bands on ONE axis, '
      + 'E/IC50, not two independent mechanisms. The kill model reaches p_hill '
      + 'through that ratio alone, so it cannot separate "too little drug '
      + 'arrived" from "the drug is too weak" — and neither can this '
      + 'attribution. Read the split as a statement about how far short the '
      + 'exposure fell, not about which upstream term caused it. SYNTHETIC '
      + 'RUN — these labels were fabricated by js/mock.js and are not a '
      + 'result.',
    thresholds: { ...MECH_THRESHOLDS, stochastic_min_p_kill: stochasticFloor },
    reconciliation: {
      n_tumor: nTumor,
      n_killed: counts.KILLED,
      n_already_dead: counts.ALREADY_DEAD,
      kill_fraction_from_labels: killFraction,
      kill_fraction_from_engine: killFraction,
      bit_exact: true,
      abs_difference: 0,
      // Said in the payload as well as on screen. In a real run the two sides
      // are computed independently and == is a genuine audit. Here there is no
      // engine: the same tally produced both, so this line is an identity and
      // establishes nothing about anything.
      method: 'STUB (js/mock.js): both sides are the same tally of the same '
        + 'fabricated label array, so this comparison is an identity and NOT '
        + 'an audit. The real check is 1.0 - (n_tumor - (KILLED + '
        + 'ALREADY_DEAD)) / n_tumor against the engine, compared with == (not '
        + 'allclose).',
    },
    synthetic: true,
  };
  return {
    labels,
    aggregate,
    kill_fraction: killFraction,
    // Tallied, not the target: the share of tumour cells the label pass found
    // at or above IC50. The bisection lands close, the tally is what ships.
    coverage_1x: nTumor ? nAdequate / nTumor : 0,
  };
}

// TMSM: magic(4) | version u32 | n_cells u32 | n_labels u32 | json_len u32
//       | aggregate json | label u8[n_cells]
export function encodeTMSM(labels, aggregate) {
  const body = new TextEncoder().encode(JSON.stringify(aggregate));
  const header = 20;
  const buf = new ArrayBuffer(header + body.length + labels.length);
  const dv = new DataView(buf);
  const u8 = new Uint8Array(buf);
  u8[0] = 0x54; u8[1] = 0x4D; u8[2] = 0x53; u8[3] = 0x4D;
  dv.setUint32(4, 1, true);
  dv.setUint32(8, labels.length, true);
  dv.setUint32(12, (aggregate.label_names || []).length, true);
  dv.setUint32(16, body.length, true);
  u8.set(body, header);
  u8.set(labels, header + body.length);
  return buf;
}

/* ---------------------------------------------------------- run store -- */

let seq = 7;
const runs = [];

/* The values here are TARGETS for the attribution, not the shipped answer.
 *
 * mechanismFor() calibrates its exposure field to `kill_fraction` and
 * `therapeutic_coverage_1x`, then tick() overwrites both with what its label
 * tally actually produced. The legend reconciles the per-cell labels against
 * kill_fraction with `==`, and the panel prints coverage two lines under a
 * breakdown of the same cells: a number drawn here and labels drawn there
 * would agree only by accident, and the demo would ship both the red DOES NOT
 * RECONCILE line the real backend refuses to emit and a coverage figure its
 * own attribution contradicts. */
function metricsFor(spec, rng) {
  if (!spec.drug) {
    return {
      kill_fraction: null, mean_tumor_auc: null,
      therapeutic_coverage_1x: null, penetration_length_um: null,
    };
  }
  const d = DRUGS.find((x) => x.name === spec.drug) || DRUGS[0];
  const bbb = d.bbb_penetrance;
  const pot = 1 / (1 + d.ic50_uM / 40);
  const mg = spec.mgmt === null || spec.mgmt === undefined ? 1 : 1 / (1 + Math.exp((spec.mgmt - 5.4) * 0.55));
  const kill = Math.max(0.001, Math.min(0.97, bbb * pot * mg * (0.85 + rng() * 0.3)));
  const auc = spec.pk.plasma_peak * spec.pk.half_life_s * spec.pk.n_doses * bbb * 1.9e3 * (0.9 + rng() * 0.2);
  return {
    kill_fraction: kill,
    mean_tumor_auc: auc,
    therapeutic_coverage_1x: Math.max(0.002, Math.min(0.99, kill * (1.9 + rng() * 0.7))),
    penetration_length_um: 60 + bbb * 380 * (0.85 + rng() * 0.3),
  };
}

/* The `reach` block, in the shape app/server/reach.py returns.
 *
 * SYNTHETIC, like everything else in this file, and it must stay contract-
 * shaped for the same reason the rest of it does: the Reach section is now the
 * product's primary answer, and a stub that omits it would let the panel's
 * empty state ship untested. Two properties of the real thing are reproduced
 * exactly rather than approximated, because they are the properties the panel
 * asserts in prose: the curve is monotone non-decreasing in the threshold, and
 * the first falloff bin is 1.0 by construction. The headline is PICKED OFF the
 * ladder rather than computed separately, so the number and the chart cannot
 * disagree here any more than they can in the backend. */
const REACH_LADDER = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];
const REACH_HEADLINE_X = 0.1;

function reachFor(spec, rng) {
  if (!spec.drug) {
    return {
      null_reason: 'requires a drug; this RunSpec has drug=null '
        + '(reconstruction-only run), so there is no exposure field and no '
        + 'reach ratio to take',
      headline_under_reached_fraction: null, headline_x: REACH_HEADLINE_X,
      curve: null, falloff: null, E_ref_uM_s: null,
    };
  }
  const d = DRUGS.find((x) => x.name === spec.drug) || DRUGS[0];
  // A decay length drives both curves, so they stay consistent with each other.
  const L = 60 + d.bbb_penetrance * 380 * (0.85 + rng() * 0.3);
  const bins = [];
  for (let i = 0; i < 14; i += 1) bins.push(i * 10 + 5);
  const fall = bins.map((dist) => Math.exp(-(dist - bins[0]) / L));
  // F(x) = fraction of tissue whose normalised exposure is below x, read off
  // the same exponential — monotone in x by construction.
  const F = REACH_LADDER.map((x) => {
    const dCut = -L * Math.log(Math.min(1, x));
    const frac = 1 - Math.exp(-dCut / (L * 2.2));
    return Math.max(0, Math.min(1, 1 - frac));
  }).map((v, i, arr) => (i ? Math.max(v, arr[i - 1]) : v));
  return {
    null_reason: null,
    headline_under_reached_fraction: F[REACH_LADDER.indexOf(REACH_HEADLINE_X)],
    headline_x: REACH_HEADLINE_X,
    curve: { x: REACH_LADDER.slice(), under_reached_fraction: F },
    falloff: {
      distance_um: bins, normalised_mean_auc: fall,
      count: bins.map(() => 200), bin_um: 10.0, min_bin_count: 30,
    },
    E_ref_uM_s: 1e5 * (0.7 + rng() * 0.6),
    E_ref_source: 'SYNTHETIC — js/mock.js fabricated this from a seeded PRNG '
      + 'against the shape of app/server/reach.py. The real value is the '
      + 'perivascular mean of penetration_length_fit(...)["drug"]',
    n_viable: 20000, n_tumor: 24000,
    viable_definition: 'SYNTHETIC placeholder for the engine’s viable '
      + 'tumour band; the real definition comes from delivery_metrics.',
    convention: 'The 0.1 threshold is a REPORTING CONVENTION, not a physical '
      + 'constant. The whole curve is published so that the choice of '
      + 'threshold is visible as a choice.',
    scale_invariance: 'SYNTHETIC RUN — this reach block was fabricated by '
      + 'js/mock.js and is not a result. In a real run this field carries the '
      + 'measured invariances of the within-map ratio and must be read there.',
    penetration_fit_error: null,
    units: { E_ref: 'uM.s', reach_ratio: 'dimensionless',
      under_reached_fraction: 'fraction of viable tumour cells',
      distance_um: 'um', normalised_mean_auc: 'dimensionless' },
  };
}

function reportFor(spec) {
  const sec = SECTIONS.find((s) => s.id === spec.section) || SECTIONS[0];
  const g = buildGeometry(sec.id);
  const counts = [0, 0, 0, 0, 0];
  for (let i = 0; i < g.n; i++) counts[g.type[i]]++;
  const names = ['TUMOR', 'STROMAL', 'IMMUNE', 'ENDOTHELIAL', 'NECROTIC'];
  const celltype_counts = {}, celltype_fractions = {};
  names.forEach((nm, i) => { celltype_counts[nm] = counts[i]; celltype_fractions[nm] = counts[i] / g.n; });
  const cols = Math.round(sec.extent_um[0] / spec.dx_um), rows = Math.round(sec.extent_um[1] / spec.dx_um);
  return {
    patient_id: sec.id,
    source: sec.panel,
    n_cells: g.n,
    n_tumor: counts[0],
    grid: [rows, cols],
    dx_um: spec.dx_um,
    extent_um: sec.extent_um,
    celltype_counts,
    celltype_fractions,
    mean_vascular_density: 0.184,
    mean_ecm_density: 0.412,
    mean_oxygen: 0.371,
    oxygen: {
      mode: spec.oxygen_mode,
      source: spec.oxygen_mode === 'markers' ? 'HIF1A/CA9/VEGFA marker composite' :
        spec.oxygen_mode === 'vessel_distance' ? 'distance transform from discrete vessels' :
          'steady-state reaction-diffusion solve',
      is_surrogate: spec.oxygen_mode === 'markers',
      caveat: spec.oxygen_mode === 'markers'
        ? 'Marker composite is a SURROGATE for pO2. It has never been calibrated against a measured oxygen electrode on this tissue.'
        : null,
      detail: `oxygen_mode=${spec.oxygen_mode}`,
    },
    oxygen_mode: spec.oxygen_mode,
    cell_state: {
      mode: spec.state_mode || 'oxygen',
      source: spec.state_mode === 'markers'
        ? 'measured proliferation transcripts; hypoxia UNMEASURED on this panel'
        : spec.state_mode === 'vessel_distance'
          ? 'update_cell_states banding exp(-d/L) in distance to measured endothelial cells'
          : 'update_cell_states banding the SIMULATED oxygen field',
      is_measured: spec.state_mode === 'markers',
      unmeasured_states: spec.state_mode === 'markers' && spec.section === 'xenium_319' ? ['HYPOXIC'] : [],
      caveat: (spec.state_mode || 'oxygen') === 'oxygen'
        ? 'NOT A MEASUREMENT: the VIABLE fraction this produces moves with the tissue-mask convention and the mesh, neither of which is a measurement.'
        : 'sim.engine re-bands cells.state from its own oxygen solve, so this route did not drive kill_fraction.',
      detail: `state_mode=${spec.state_mode || 'oxygen'}`,
    },
    state_mode: spec.state_mode || 'oxygen',
    ecm_markers_used: ['COL1A1', 'FN1', 'SPARC'],
    tissue_mask_mode: 'alpha-shape of detected cells',
    tissue_coverage_fraction: 0.6412,
    sensitivity: {
      mode: spec.sensitivity_mode,
      upper_bound: spec.sensitivity_mode === 'uniform',
      is_uniform: spec.sensitivity_mode === 'uniform',
      determinant_gene: spec.sensitivity_mode === 'determinant' ? 'MGMT' : null,
      kind: spec.sensitivity_mode === 'uniform' ? 'none (uniform 1.0)' : 'per-cell determinant expression',
      mean_all_cells: spec.sensitivity_mode === 'uniform' ? 1.0 : 0.826,
      mean_tumor: spec.sensitivity_mode === 'uniform' ? 1.0 : 0.791,
      mgmt_applied: spec.mgmt !== null && spec.mgmt !== undefined,
      mgmt_scalar: spec.mgmt !== null && spec.mgmt !== undefined ? 1 / (1 + Math.exp((spec.mgmt - 5.4) * 0.55)) : null,
      mgmt_input: spec.mgmt ?? null,
      caveat: 'ABSOLUTE kill_fraction is CONDITIONAL on this sensitivity assumption. It is deterministic and documented, but it is NOT a measured drug response — this section has no outcome label.',
    },
  };
}

function provenanceFor(spec, rng) {
  const d = DRUGS.find((x) => x.name === spec.drug);
  const caveats = [
    'This section carries no measured drug-response label. Every absolute kill number is a model output conditional on the sensitivity assumption in force, not a prediction that has been scored.',
    'Vascular geometry is reconstructed from endothelial marker positions, not from an imaged lumen. Vessel radius is assumed, not measured.',
  ];
  if (spec.oxygen_mode === 'markers') {
    caveats.push('oxygen_mode=markers: the oxygen field is a transcriptional surrogate. It has not been calibrated against a measured pO2 on this tissue.');
  }
  if (spec.state_mode && spec.state_mode !== 'oxygen') {
    caveats.push(`state_mode=${spec.state_mode}: sim.engine re-bands cells.state from its own oxygen solve before the kill, so this route changed the REPORTED band occupancies and their provenance, not this run's kill_fraction.`);
  }
  if (spec.section === 'xenium_319') {
    caveats.push('ABSOLUTE KILL MUST NOT BE QUOTED FROM THIS SECTION UNDER ANY STATE ROUTE: its panel measures no hypoxia, so no route clears the bar and kill moves up to 0.3252 absolute across routes. What is defensible is the MAP: within this one compound on this one section, the RELATIVE exposure between regions. The absolute level is not, and a comparison against another compound is not — that claim was measured and withdrawn (results/pk_ranking_robustness.md, results/capacity_ranking_robustness.md).');
  }
  if (spec.vascular_mode === 'kde') {
    caveats.push('vascular_mode=kde did not converge under grid refinement in the held-out check; discrete vessels are the production default.');
  }
  if (spec.sensitivity_mode === 'uniform') {
    caveats.push('sensitivity_mode=uniform sets every cell to 1.0, so kill_fraction is a NO-RESISTANCE UPPER BOUND, not an expected response.');
  }
  if (spec.mgmt !== null && spec.mgmt !== undefined) {
    caveats.push('MGMT enters as a single per-patient scalar. The spatial MGMT signal was measured and found to be noise; do not read spatial MGMT structure into this run.');
  }
  return {
    alpha: 148.6 + rng() * 9,
    ic50_uM: d ? d.ic50_uM : null,
    ic50_source: d ? d.ic50_source : null,
    sensitivity_note: spec.sensitivity_mode === 'uniform'
      ? 'uniform 1.0 across all cells (drug-agnostic reconstruction default) — every kill number is a no-resistance upper bound on this one map, not an expected response and not a number to set beside another compound'
      : 'per-cell determinant expression scaled to [0,1]; the resulting kill is conditional on a determinant this panel may not measure for this compound',
    caveats,
  };
}

const STAGES = [
  ['loading section', 0.06],
  ['segment + phenotype', 0.16],
  ['vessel reconstruction', 0.3],
  ['ECM / tissue fields', 0.42],
  ['oxygen solve', 0.58],
  ['IFP + velocity', 0.66],
  ['drug transport PDE', 0.9],
  ['pharmacodynamic kill', 0.97],
  ['metrics + provenance', 1.0],
];

export function createRun(spec) {
  const v = validateSpec(spec);
  if (!v.ok) {
    const err = new Error('spec did not validate');
    err.status = 422;
    err.payload = { detail: v.errors };
    throw err;
  }
  seq += 1;
  const id = 'r_' + String(seq).padStart(5, '0');
  const rng = mulberry32(hash(id + JSON.stringify(spec)));
  const durationMs = 9000 + rng() * 9000;
  const run = {
    run_id: id,
    spec: JSON.parse(JSON.stringify(spec)),
    status: 'queued',
    report: null,
    metrics: null,
    provenance: null,
    mechanism: null,
    mechanism_reason: null,
    extra_metrics: null,
    runtime_s: null,
    error: null,
    created_at: Date.now() / 1000,
    _t0: performance.now(),
    _dur: durationMs,
    _rng: rng,
    _cancelled: false,
    _log: ['[queue] accepted ' + id, `[queue] section=${spec.section} dx=${spec.dx_um}um drug=${spec.drug ?? 'none'}`],
  };
  runs.unshift(run);
  return { run_id: id };
}

function tick(run) {
  if (run.status === 'done' || run.status === 'failed' || run.status === 'cancelled') return run;
  const el = performance.now() - run._t0;
  if (el < 700) { run.status = 'queued'; return run; }
  run.status = 'running';
  const p = Math.min(1, (el - 700) / run._dur);
  let stage = STAGES[0][0];
  for (const [nm, upto] of STAGES) { if (p <= upto) { stage = nm; break; } stage = nm; }
  run._progress = p;
  run._stage = stage;
  const wantLines = Math.floor(p * 14);
  while (run._log.length < wantLines + 2) {
    run._log.push(`[${stage}] step ${run._log.length - 1} ok`);
  }
  if (p >= 1) {
    run.status = 'done';
    run.runtime_s = run._dur / 1000 * 4.1;
    run.report = reportFor(run.spec);
    run.metrics = metricsFor(run.spec, run._rng);
    // The attribution calibrates to those two targets and then REPLACES them
    // with what its own tally gave, so the breakdown, the kill fraction and
    // the coverage figure on screen are three readings of one label array
    // rather than three independent draws that happen to sit near each other.
    run._mech = mechanismFor(run.spec, run._rng, run.metrics);
    run.mechanism = run._mech ? run._mech.aggregate : null;
    run.mechanism_reason = run._mech ? null
      : 'reconstruction-only run: no compound was given, so there is no kill '
        + 'to attribute';
    if (run._mech) {
      run.metrics.kill_fraction = run._mech.kill_fraction;
      run.metrics.therapeutic_coverage_1x = run._mech.coverage_1x;
    }
    // A MEASURED result about the compound, not about this run, so it does not
    // depend on anything synthesised above and is not synthesised itself: the
    // same three tiers app/server/jobs.py emits, from the same census.
    run.extra_metrics = { verdict_stability: verdictStability(run.spec.drug) };
    run.reach = reachFor(run.spec, run._rng);
    run.provenance = provenanceFor(run.spec, run._rng);
    run._log.push('[done] wrote results');
  }
  return run;
}

export function listRuns() { runs.forEach(tick); return runs.map(publicRun); }
export function getRun(id) {
  const r = runs.find((x) => x.run_id === id);
  if (!r) { const e = new Error('no such run'); e.status = 404; throw e; }
  return publicRun(tick(r));
}
export function runStatus(id) {
  const r = runs.find((x) => x.run_id === id);
  if (!r) { const e = new Error('no such run'); e.status = 404; throw e; }
  tick(r);
  return {
    status: r.status,
    progress: r.status === 'done' ? 1 : (r._progress || 0),
    stage: r.status === 'done' ? 'complete' : r.status === 'cancelled' ? 'cancelled' : (r._stage || 'queued'),
    elapsed_s: (performance.now() - r._t0) / 1000,
    log_tail: r._log.slice(-12),
  };
}
export function cancelRun(id) {
  const r = runs.find((x) => x.run_id === id);
  if (!r) { const e = new Error('no such run'); e.status = 404; throw e; }
  if (r.status === 'done') { const e = new Error('run already complete'); e.status = 409; throw e; }
  r.status = 'cancelled';
  r._log.push('[cancel] terminated by operator');
  return { run_id: id, status: 'cancelled' };
}

function publicRun(r) {
  return {
    run_id: r.run_id, spec: r.spec, status: r.status, report: r.report,
    metrics: r.metrics, reach: r.reach, provenance: r.provenance,
    mechanism: r.mechanism || null,
    mechanism_reason: r.mechanism_reason || null,
    extra_metrics: r.extra_metrics || null,
    runtime_s: r.runtime_s, error: r.error, created_at: r.created_at,
  };
}

/** Two finished runs so history + comparison have content on first paint.
 *
 * BOTH ARMS ARE THE SAME COMPOUND, and that is deliberate. The seed used to be
 * temozolomide vs lomustine, which made the shipped demo of the comparison
 * panel a two-DRUG comparison — the one reading this product no longer
 * supports. The legitimate use of A/B is a single compound under two
 * conditions, so the demo now varies the dose and the mesh and holds the drug
 * fixed. A shipped example teaches people what the tool is for.
 */
export function seedHistory() {
  if (runs.length) return;
  const a = { ...DEFAULT_SPEC, drug: 'temozolomide' };
  const b = {
    ...DEFAULT_SPEC,
    drug: 'temozolomide',
    dx_um: 15.0,
    pk: { ...(DEFAULT_SPEC.pk || {}) },
  };
  if (b.pk && typeof b.pk.plasma_peak === 'number') b.pk.plasma_peak = b.pk.plasma_peak * 2;
  for (const spec of [a, b]) {
    const { run_id } = createRun(spec);
    const r = runs.find((x) => x.run_id === run_id);
    r._t0 = performance.now() - r._dur - 2000;
    r.created_at = Date.now() / 1000 - (spec === a ? 5400 : 1800);
    tick(r);
  }
}

/* ------------------------------------------------------------- router -- */

/** Route a contract path to a mock response. Returns {json} or {buffer}. */
export function handle(method, path, body) {
  const p = path.replace(/^\/+/, '/').split('?')[0];
  if (method === 'GET' && p === '/api/sections') return { json: SECTIONS };
  if (method === 'GET' && p === '/api/drugs') return { json: DRUGS };
  if (method === 'GET' && p === '/api/defaults') return { json: JSON.parse(JSON.stringify(DEFAULT_SPEC)) };
  if (method === 'POST' && p === '/api/validate') return { json: validateSpec(body) };
  if (method === 'POST' && p === '/api/runs') return { json: createRun(body) };
  if (method === 'GET' && p === '/api/runs') return { json: listRuns() };

  let m = /^\/api\/runs\/([^/]+)$/.exec(p);
  if (m && method === 'GET') return { json: getRun(m[1]) };
  if (m && method === 'DELETE') return { json: cancelRun(m[1]) };

  m = /^\/api\/runs\/([^/]+)\/status$/.exec(p);
  if (m && method === 'GET') return { json: runStatus(m[1]) };

  m = /^\/api\/runs\/([^/]+)\/geometry$/.exec(p);
  if (m && method === 'GET') {
    const r = getRun(m[1]);
    return { buffer: encodeTMSG(buildGeometry(r.spec.section)) };
  }

  m = /^\/api\/runs\/([^/]+)\/mechanism$/.exec(p);
  if (m && method === 'GET') {
    const r = runs.find((x) => x.run_id === m[1]);
    if (!r) { const e = new Error('no such run'); e.status = 404; throw e; }
    tick(r);
    // 404 with the run's OWN reason, the same shape app/server/api.py returns
    // for a run with nothing to attribute — main.js surfaces `detail.reason`
    // verbatim rather than inventing a message.
    if (!r._mech) {
      const e = new Error('no mechanism attribution');
      e.status = 404;
      e.payload = { detail: { reason: r.mechanism_reason
        || 'this run carries no per-cell mechanism attribution' } };
      throw e;
    }
    return { buffer: encodeTMSM(r._mech.labels, r._mech.aggregate) };
  }

  m = /^\/api\/runs\/([^/]+)\/field\/([^/]+)$/.exec(p);
  if (m && method === 'GET') {
    const r = getRun(m[1]);
    return { buffer: encodeTMSF(buildField(r.spec.section, m[2], r.spec.dx_um)) };
  }

  const e = new Error(`mock: no route for ${method} ${p}`);
  e.status = 404;
  throw e;
}
