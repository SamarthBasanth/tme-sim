// api_baked.js — the STATIC EXPORT's replacement for js/api.js.
//
// Same module surface, same return shapes, same binary decode path. The only
// difference is where the bytes come from: this module issues no network
// request of any kind. Every response is read out of the payload baked into
// the file by scripts/build_studio_export.py.
//
// What is real here:
//   * /api/sections, /api/drugs, /api/defaults  — verbatim backend responses
//   * /api/runs, /api/runs/{id}, /status        — verbatim RunResults, produced
//                                                 by the production engine
//   * /geometry, /field/{name}                  — the identical TMSG / TMSF
//                                                 bytes the server would serve
//
// What CANNOT be real, and is refused rather than faked:
//   * POST /api/runs      — there is no engine here
//   * DELETE /api/runs/id — there is no job to cancel
//   * POST /api/validate  — for any spec that was not precomputed
//
// Nothing in this file synthesises a number. If a route has no baked answer it
// raises an ApiError that says so in as many words.

const BAKED = (typeof window !== 'undefined' && window.__TMESIM_BAKED) || null;

export class ApiError extends Error {
  constructor(message, { status = 0, path = '', payload = null, cause = null } = {}) {
    super(message);
    this.name = 'ApiError';
    this.status = status;
    this.path = path;
    this.payload = payload;
    this.cause = cause;
  }
}

export const api = {
  base: '',
  /** false: this is NOT js/mock.js. Nothing here is synthesised. */
  stub: false,
  /** There is no backend to reach. Said plainly rather than left ambiguous. */
  reachable: false,
  staticExport: true,
  lastLatencyMs: null,
  onModeChange: null,
};

if (!BAKED) {
  // Loud, immediate, and impossible to mistake for a working app.
  throw new Error(
    'studio_export: window.__TMESIM_BAKED is missing. This file is a static '
    + 'export and carries its data inline; without it there is nothing to show.',
  );
}

/* ----------------------------------------------------------- payloads --- */

const bufCache = new Map();

/** base64 -> ArrayBuffer, with the byte length the encoder recorded. */
function decodeB64(b64) {
  const bin = atob(b64);
  const n = bin.length;
  const buf = new ArrayBuffer(n);
  const u8 = new Uint8Array(buf);
  for (let i = 0; i < n; i++) u8[i] = bin.charCodeAt(i);
  return buf;
}

function blob(key) {
  if (bufCache.has(key)) return bufCache.get(key);
  const b64 = BAKED.blobs[key];
  if (b64 === undefined) {
    throw new ApiError(
      `baked payload ${String(key)} is not in this export`, { status: 404 });
  }
  const buf = decodeB64(b64);
  bufCache.set(key, buf);
  return buf;
}

const timed = (fn) => {
  const t0 = performance.now();
  const out = fn();
  api.lastLatencyMs = performance.now() - t0;
  return out;
};

function runOr404(id) {
  const r = BAKED.runs.find((x) => x.run_id === id);
  if (!r) {
    throw new ApiError(
      `no such run ${id} in this static export. It carries `
      + `${BAKED.runs.length} precomputed runs: `
      + `${BAKED.runs.map((x) => x.run_id).join(', ')}.`,
      { status: 404, path: `/api/runs/${id}` },
    );
  }
  return r;
}

/* ------------------------------------------------------------- routes --- */

export async function probe() {
  return { live: true, sections: BAKED.sections, staticExport: true };
}

export const getSections = async () => timed(() => BAKED.sections);
export const getDrugs = async () => timed(() => BAKED.drugs);
export const getDefaults = async () => timed(() => JSON.parse(JSON.stringify(BAKED.defaults)));
export const listRuns = async () => timed(() => BAKED.runs);
export const getRun = async (id) => timed(() => runOr404(id));

export const getStatus = async (id) => timed(() => {
  const r = runOr404(id);
  return {
    run_id: r.run_id,
    status: r.status,
    progress: r.progress,
    stage: r.stage,
    elapsed_s: r.runtime_s,
    log_tail: r.log_tail || [],
  };
});

/**
 * The three precomputed specs carry the verdict the real ``POST /api/validate``
 * returned for them, recorded at bake time. Anything else is UNVALIDATED, and
 * is reported as unvalidated rather than waved through.
 */
export async function postValidate(spec) {
  return timed(() => {
    const hit = BAKED.validate[canonical(spec)];
    if (hit) return hit;

    // A front end that carries a RunSpec field this backend build did not have
    // (or the other way round) must not read as "your spec is broken". Drop the
    // fields no baked spec knows, and if THAT matches, hand back the recorded
    // verdict with the same "unrecognised field(s)" warning the real
    // /api/validate emits for them.
    const known = new Set();
    for (const r of BAKED.runs) Object.keys(r.spec || {}).forEach((k) => known.add(k));
    const trimmed = {};
    const extra = [];
    for (const k of Object.keys(spec || {})) {
      if (known.has(k)) trimmed[k] = spec[k]; else extra.push(k);
    }
    const near = extra.length ? BAKED.validate[canonical(trimmed)] : null;
    if (near) {
      return {
        ok: near.ok,
        errors: near.errors.slice(),
        warnings: near.warnings.concat([
          `unrecognised field(s) ${JSON.stringify(extra.sort())}: not present in the `
          + 'RunSpec that was validated when this export was built, so the verdict '
          + 'shown ignores them.',
        ]),
      };
    }
    return {
      ok: false,
      errors: [{
        path: '',
        message:
          'STATIC EXPORT — this RunSpec has NOT been validated. Validation is a '
          + 'backend route (POST /api/validate) that reads the section registry, '
          + 'the drug library and the measured advection CFL of the reconstruction; '
          + 'none of that is in this file. Only the precomputed specs carry a real '
          + 'verdict. Run "./.venv/bin/python -m app.server" to validate an edited spec.',
      }],
      warnings: [],
    };
  });
}

export async function postRun() {
  throw new ApiError(
    'STATIC EXPORT — POST /api/runs is not available. There is no engine behind '
    + 'this page; the runs shown were computed ahead of time by the real backend. '
    + 'To submit a run, start the local server: ./.venv/bin/python -m app.server',
    { status: 501, path: '/api/runs' },
  );
}

export async function cancelRun(id) {
  throw new ApiError(
    `STATIC EXPORT — DELETE /api/runs/${id} is not available. Nothing is running: `
    + 'every run in this export finished before the file was written.',
    { status: 501, path: `/api/runs/${id}` },
  );
}

/* ------------------------------------------- routes with no baked answer -- */
/* Reports and projects are server-side state: a report is rendered by the
   backend and a project is mutable bookkeeping. Neither can exist in a file.
   They refuse with the reason rather than returning an empty list, which would
   read as "you have no projects" instead of "this cannot work here".
   build_studio_export.py fails the build if app/web/js/api.js grows an export
   that is missing here, so this list cannot silently fall behind. */
function noServer(route) {
  return () => {
    throw new ApiError(
      `STATIC EXPORT — ${route} is server-side state and is not available in a `
      + 'single file. Start the local server to use it: '
      + './.venv/bin/python -m app.server',
      { status: 501, path: route },
    );
  };
}

export const postReport = noServer('POST /api/reports');
export const listReports = noServer('GET /api/reports');
export const getReport = noServer('GET /api/reports/{id}');
export const cancelReport = noServer('DELETE /api/reports/{id}');
export const reportHtmlUrl = () => null;   // no URL exists; callers must handle null

export const listProjects = noServer('GET /api/projects');
export const createProject = noServer('POST /api/projects');
export const getProject = noServer('GET /api/projects/{id}');
export const patchProject = noServer('PATCH /api/projects/{id}');
export const deleteProject = noServer('DELETE /api/projects/{id}');
export const addSpecimen = noServer('POST /api/projects/{id}/specimens');
export const getSpecimen = noServer('GET /api/projects/{id}/specimens/{sid}');
export const patchSpecimen = noServer('PATCH /api/projects/{id}/specimens/{sid}');
export const deleteSpecimen = noServer('DELETE /api/projects/{id}/specimens/{sid}');
export const attachRun = noServer('POST .../specimens/{sid}/runs');
export const detachRun = noServer('DELETE .../specimens/{sid}/runs/{run}');
export const compareProject = noServer('GET /api/projects/{id}/compare');
export const getCapability = noServer('GET /api/capability/{section}');

/* --- compound library: mutable server-side state ------------------------- */
/* A stored compound is a row the backend writes to app/compounds/, and
   compiling SMILES is an rdkit call in the server process. Neither can happen
   in a file, and an empty list would read as "your library is empty" rather
   than "this cannot work here". */
export const postCompoundCompile = noServer('POST /api/compounds/compile');
export const postCompound = noServer('POST /api/compounds');
export const listCompounds = noServer('GET /api/compounds');
export const getCompound = noServer('GET /api/compounds/{id}');
export const deleteCompound = noServer('DELETE /api/compounds/{id}');

/* --- bring-your-own tissue ----------------------------------------------- */
/* Sections register BY FILESYSTEM PATH on the server's own machine — nothing
   is ever uploaded — so a scan is the server reading a directory it can see.
   A file opened from file:// has no such directory and no such server, and
   the confinement that bounds the real route (TMESIM_STUDIO_DATA_ROOTS) lives
   in that server's environment. Refused, with the command that does work. */
export const postSectionScan = noServer('POST /api/sections/scan');
export const deleteSection = noServer('DELETE /api/sections/{id}');

/**
 * No-op, deliberately, and NOT a refusal.
 *
 * This is the one name in this list that is not a route: in js/api.js it sets
 * a header for later requests and remembers it in localStorage. This export
 * issues no request, so there is no header to set and no gate to satisfy —
 * `probe()` reports the baked data as reachable and ui/apikey.js never opens.
 *
 * It stays callable rather than throwing because a caller reaching it is not
 * an error the way `postRun()` is, and it deliberately does NOT write to
 * localStorage: persisting a secret from a page that cannot use it would
 * leave a real API key in browser storage for a file that has no server to
 * send it to.
 */
export function setApiKey() { /* no request is ever issued from this file */ }

export async function getGeometry(id) {
  return timed(() => {
    const r = runOr404(id);
    const key = (BAKED.payloads[r.run_id] || {}).geometry;
    if (!key) {
      throw new ApiError(
        `cell geometry was not baked for ${r.run_id}`, { status: 404 });
    }
    return decodeTMSG(blob(key));
  });
}

/**
 * Per-cell mechanism attribution (TMSM) and the drug-concentration time-lapse
 * (TMSL). Both are ordinary run payloads and ARE baked when the run produced
 * them — but a run only produces them if its RunSpec asked, so a run that did
 * not is reported as such rather than as an empty result.
 */
export async function getMechanism(id) {
  return timed(() => decodeTMSM(runPayload(id, 'mechanism', 'mechanism attribution',
    'attribute_mechanism')));
}

export async function getTimelapse(id) {
  return timed(() => decodeTMSL(runPayload(id, 'timelapse', 'drug time-lapse',
    'timelapse_frames')));
}

function runPayload(id, name, label, specKey) {
  const r = runOr404(id);
  const key = (BAKED.payloads[r.run_id] || {})[name];
  if (key) return blob(key);
  const asked = r.spec && r.spec[specKey];
  throw new ApiError(
    `${label} was not produced by run ${r.run_id}: its RunSpec `
    + (asked === undefined
      ? `carries no ${specKey}, so the engine wrote no ${name} payload `
        + '(the run predates the feature or ran with it off)'
      : `has ${specKey} = ${JSON.stringify(asked)}`)
    + '. Nothing is substituted.',
    { status: 404, path: `/api/runs/${id}/${name}` },
  );
}

export async function getField(id, name) {
  return timed(() => {
    const r = runOr404(id);
    const key = (BAKED.payloads[r.run_id] || {})[name];
    if (!key) {
      throw new ApiError(
        `field ${name} was not produced by run ${r.run_id}`
        + (name === 'drug_auc' && !(r.spec && r.spec.drug)
          ? ' (this RunSpec had drug=null, so there is no drug field)'
          : ''),
        { status: 404 });
    }
    const f = decodeTMSF(blob(key));
    f.name = name;
    return f;
  });
}

/** Stable key order so a spec hashes the same however it was assembled. */
function canonical(o) {
  if (o === null || typeof o !== 'object') return JSON.stringify(o);
  if (Array.isArray(o)) return '[' + o.map(canonical).join(',') + ']';
  const keys = Object.keys(o).sort();
  return '{' + keys.map((k) => JSON.stringify(k) + ':' + canonical(o[k])).join(',') + '}';
}

/* ------------------------------------------------------- binary decode -- */
/* Byte-for-byte the decoders from js/api.js — the export must exercise the
   production decode -> quantise -> upload path, not a shortcut around it. */

const MAGIC = (b, s) => String.fromCharCode(b[0], b[1], b[2], b[3]) === s;

export function decodeTMSG(buf) {
  const u8 = new Uint8Array(buf);
  if (buf.byteLength < 20 || !MAGIC(u8, 'TMSG')) {
    throw new ApiError('geometry payload is not TMSG (bad magic) — refusing to render it', { status: 502 });
  }
  const dv = new DataView(buf);
  const version = dv.getUint32(4, true);
  const n = dv.getUint32(8, true);
  const ex = dv.getFloat32(12, true);
  const ey = dv.getFloat32(16, true);
  const need = 20 + n * 8 + n;
  if (buf.byteLength < need) {
    throw new ApiError(`TMSG truncated: ${buf.byteLength} B for n_cells=${n} (need ${need} B)`, { status: 502 });
  }
  return {
    version, n, extent: [ex, ey],
    x: new Float32Array(buf, 20, n),
    y: new Float32Array(buf, 20 + n * 4, n),
    type: new Uint8Array(buf, 20 + n * 8, n),
    bytes: buf.byteLength,
  };
}

export function decodeTMSF(buf) {
  const u8 = new Uint8Array(buf);
  if (buf.byteLength < 28 || !MAGIC(u8, 'TMSF')) {
    throw new ApiError('field payload is not TMSF (bad magic) — refusing to render it', { status: 502 });
  }
  const dv = new DataView(buf);
  const version = dv.getUint32(4, true);
  const rows = dv.getUint32(8, true);
  const cols = dv.getUint32(12, true);
  const dx = dv.getFloat32(16, true);
  const min = dv.getFloat32(20, true);
  const max = dv.getFloat32(24, true);
  const need = 28 + rows * cols * 4;
  if (buf.byteLength < need) {
    throw new ApiError(`TMSF truncated: ${buf.byteLength} B for ${rows}x${cols} (need ${need} B)`, { status: 502 });
  }
  return {
    version, rows, cols, dx, min, max,
    values: new Float32Array(buf, 28, rows * cols),
    bytes: buf.byteLength,
  };
}

/**
 * TMSM: magic(4) | version u32 | n_cells u32 | n_labels u32 | json_len u32
 *       | aggregate json[json_len] (UTF-8) | label u8[n_cells]
 *
 * One byte per cell in geometry order, so the labels go straight into a vertex
 * buffer. The aggregate breakdown rides along as a short JSON block — a legend
 * without its counts is decoration, and a second round trip for four integers
 * is not worth a request.
 */
export function decodeTMSM(buf) {
  const u8 = new Uint8Array(buf);
  if (buf.byteLength < 20 || !MAGIC(u8, 'TMSM')) {
    throw new ApiError('mechanism payload is not TMSM (bad magic) — refusing to render it', { status: 502 });
  }
  const dv = new DataView(buf);
  const version = dv.getUint32(4, true);
  const n = dv.getUint32(8, true);
  const nLabels = dv.getUint32(12, true);
  const jsonLen = dv.getUint32(16, true);
  const need = 20 + jsonLen + n;
  if (buf.byteLength < need) {
    throw new ApiError(`TMSM truncated: ${buf.byteLength} B for n_cells=${n} (need ${need} B)`, { status: 502 });
  }
  let aggregate;
  try {
    aggregate = JSON.parse(new TextDecoder().decode(new Uint8Array(buf, 20, jsonLen)));
  } catch (e) {
    throw new ApiError('TMSM aggregate block is not valid JSON', { status: 502, cause: e });
  }
  return {
    version, n_cells: n, n_labels: nLabels, aggregate,
    labels: new Uint8Array(buf, 20 + jsonLen, n),
    bytes: buf.byteLength,
  };
}

/**
 * TMSL: magic(4) | version u32 | n_frames u32 | interval_s f32 | window_s f32
 *       | n_frames x [ frame u32 | time_s f32 | <a complete TMSF block> ]
 *
 * Each frame is a byte-identical TMSF block, so it is decoded by the reader
 * already in this file — the container adds a frame index and a time stamp and
 * nothing else. Values are Float32Array views over the SAME response buffer, so
 * a 24-frame series costs one allocation, not 24.
 */
export function decodeTMSL(buf) {
  const u8 = new Uint8Array(buf);
  if (buf.byteLength < 20 || !MAGIC(u8, 'TMSL')) {
    throw new ApiError('timelapse payload is not TMSL (bad magic) — refusing to render it', { status: 502 });
  }
  const dv = new DataView(buf);
  const version = dv.getUint32(4, true);
  const nFrames = dv.getUint32(8, true);
  const interval = dv.getFloat32(12, true);
  const window = dv.getFloat32(16, true);

  const frames = [];
  let off = 20;
  for (let i = 0; i < nFrames; i++) {
    if (off + 8 + 28 > buf.byteLength) {
      throw new ApiError(`TMSL truncated at frame ${i} of ${nFrames}`, { status: 502 });
    }
    const index = dv.getUint32(off, true);
    const t = dv.getFloat32(off + 4, true);
    if (index !== i) {
      throw new ApiError(`TMSL frame index out of order: got ${index}, expected ${i}`, { status: 502 });
    }
    off += 8;
    if (!MAGIC(new Uint8Array(buf, off, 4), 'TMSF')) {
      throw new ApiError(`TMSL frame ${i} is not a TMSF block`, { status: 502 });
    }
    const rows = dv.getUint32(off + 8, true);
    const cols = dv.getUint32(off + 12, true);
    const size = 28 + rows * cols * 4;
    if (off + size > buf.byteLength) {
      throw new ApiError(`TMSL frame ${i} truncated (${rows}x${cols})`, { status: 502 });
    }
    frames.push({
      frame: index,
      time_s: t,
      rows, cols,
      dx: dv.getFloat32(off + 16, true),
      min: dv.getFloat32(off + 20, true),
      max: dv.getFloat32(off + 24, true),
      values: new Float32Array(buf, off + 28, rows * cols),
    });
    off += size;
  }
  return { version, n_frames: nFrames, interval_s: interval, window_s: window, frames, bytes: buf.byteLength };
}

/**
 * Client-side synthetic TMSG for the render benchmark. This one is honest to
 * keep offline: it measures THIS machine's GPU on a phantom point cloud and
 * never claims to be a specimen — the header says `phantom_850k` while it is
 * loaded.
 */
export function synthTMSG(n, extentX = 24000, extentY = 18000, seed = 12345) {
  const header = 20;
  const buf = new ArrayBuffer(header + n * 8 + n);
  const u8 = new Uint8Array(buf);
  u8[0] = 0x54; u8[1] = 0x4D; u8[2] = 0x53; u8[3] = 0x47;
  const dv = new DataView(buf);
  dv.setUint32(4, 1, true);
  dv.setUint32(8, n, true);
  dv.setFloat32(12, extentX, true);
  dv.setFloat32(16, extentY, true);
  const xs = new Float32Array(buf, header, n);
  const ys = new Float32Array(buf, header + n * 4, n);
  const ts = new Uint8Array(buf, header + n * 8, n);
  let s = seed >>> 0;
  const rnd = () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; };
  const cx = extentX / 2, cy = extentY / 2, R = Math.min(extentX, extentY) * 0.42;
  for (let i = 0; i < n; i++) {
    const a = rnd() * Math.PI * 2;
    const rr = Math.sqrt(rnd()) * R * (0.6 + 0.6 * rnd());
    xs[i] = Math.max(0, Math.min(extentX, cx + Math.cos(a) * rr * 1.3));
    ys[i] = Math.max(0, Math.min(extentY, cy + Math.sin(a) * rr));
    ts[i] = (rnd() * 5) | 0;
  }
  return buf;
}

/** js/mock.js is deliberately NOT bundled: nothing in this export is synthetic. */
export const mock = null;
