// api.js — the ONLY place that talks to app/server.
//
// Routes consumed (all of API_CONTRACT.md v1):
//   GET    /api/sections
//   POST   /api/sections/scan                 (bring your own tissue)
//   DELETE /api/sections/{id}                 (remove an ingested section)
//   GET    /api/drugs
//   GET    /api/defaults
//   POST   /api/runs
//   GET    /api/runs
//   GET    /api/runs/{id}
//   GET    /api/runs/{id}/status
//   DELETE /api/runs/{id}
//   GET    /api/runs/{id}/geometry          (TMSG binary)
//   GET    /api/runs/{id}/field/{name}      (TMSF binary)
//   POST   /api/validate
//
// Policy:
//  * one automatic retry on a network error or 5xx, with a 400 ms backoff;
//    4xx is never retried (the request itself is wrong).
//  * a failure is surfaced, never swallowed. `error` on a RunResult is shown.
//  * when the backend is unreachable at boot we fall back to js/mock.js and
//    set `api.stub = true`; the shell then labels every surface STUB DATA.

import * as mock from './mock.js';

const RETRY_MS = 400;
const TIMEOUT_MS = 30000;

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

const API_KEY_STORAGE = 'tmesim_api_key';

export const api = {
  base: '',
  stub: false,
  /** null until probe() has run. */
  reachable: null,
  lastLatencyMs: null,
  onModeChange: null,
  /** Set once probe() gets a 401 -- a real server that requires a key we
   * don't have (or have wrong), as distinct from "unreachable" (status 0/5xx),
   * which falls back to stub data instead. main.js shows the access-key gate
   * when this is true and clears it once a working key is stored. */
  authRequired: false,
  apiKey: (() => {
    try { return localStorage.getItem(API_KEY_STORAGE) || ''; }
    catch { return ''; }
  })(),
};

/** Persist a new key (or clear with '') and use it on every subsequent call. */
export function setApiKey(key) {
  api.apiKey = key || '';
  try {
    if (api.apiKey) localStorage.setItem(API_KEY_STORAGE, api.apiKey);
    else localStorage.removeItem(API_KEY_STORAGE);
  } catch { /* localStorage unavailable (private mode, etc.) -- key still works this session */ }
}

function setStub(on, why) {
  if (api.stub === on) return;
  api.stub = on;
  api.stubReason = why || null;
  if (api.onModeChange) api.onModeChange(api);
}

async function rawFetch(method, path, { body, binary, signal, timeoutMs } = {}) {
  const ctrl = new AbortController();
  const timer = setTimeout(() => ctrl.abort(), timeoutMs || TIMEOUT_MS);
  if (signal) signal.addEventListener('abort', () => ctrl.abort(), { once: true });
  const t0 = performance.now();
  try {
    const headers = body !== undefined
      ? { 'content-type': 'application/json', accept: binary ? 'application/octet-stream' : 'application/json' }
      : { accept: binary ? 'application/octet-stream' : 'application/json' };
    if (api.apiKey) headers['x-api-key'] = api.apiKey;
    const res = await fetch(api.base + path, {
      method,
      headers,
      body: body !== undefined ? JSON.stringify(body) : undefined,
      signal: ctrl.signal,
      cache: 'no-store',
    });
    api.lastLatencyMs = performance.now() - t0;
    if (!res.ok) {
      let payload = null, text = '';
      try {
        const ct = res.headers.get('content-type') || '';
        if (ct.includes('json')) payload = await res.json();
        else text = (await res.text()).slice(0, 400);
      } catch { /* body unreadable — the status is still the story */ }
      throw new ApiError(
        detailOf(payload) || text || `${res.status} ${res.statusText}`,
        { status: res.status, path, payload },
      );
    }
    if (binary) return await res.arrayBuffer();
    if (res.status === 204) return null;
    return await res.json();
  } finally {
    clearTimeout(timer);
  }
}

function detailOf(payload) {
  if (!payload) return '';
  if (typeof payload.detail === 'string') return payload.detail;
  if (Array.isArray(payload.detail)) {
    return payload.detail
      .map((d) => (d && d.msg ? `${(d.loc || []).slice(1).join('.')}: ${d.msg}` : (d.message || JSON.stringify(d))))
      .join('; ');
  }
  if (typeof payload.error === 'string') return payload.error;
  return '';
}

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

/** One call with exactly one retry on network / 5xx, then mock fallback. */
async function call(method, path, opts = {}) {
  // `noFallback` opts a route OUT of both the retry and the mock: one attempt,
  // and every failure propagates. Used by the ingest routes, where a retry
  // would restart a minutes-long synchronous read the server is still doing,
  // and where degrading to js/mock.js would answer a question about the
  // operator's own tissue with a seeded PRNG.
  if (opts.noFallback) return rawFetch(method, path, opts);
  if (api.stub) return fromMock(method, path, opts);

  for (let attempt = 0; attempt < 2; attempt++) {
    try {
      const out = await rawFetch(method, path, opts);
      api.reachable = true;
      return out;
    } catch (err) {
      const status = err instanceof ApiError ? err.status : 0;
      const retryable = status === 0 || status >= 500;
      if (retryable && attempt === 0) { await sleep(RETRY_MS); continue; }
      if (status === 0 || status >= 500) {
        // The server is genuinely not answering. Degrade to the stub rather
        // than leaving the operator with a dead workstation, and say so.
        api.reachable = false;
        setStub(true, status === 0 ? 'backend unreachable' : `backend ${status}`);
        return fromMock(method, path, opts);
      }
      throw err;
    }
  }
  throw new ApiError('unreachable', { path });
}

function fromMock(method, path, opts) {
  try {
    const r = mock.handle(method, path, opts.body);
    if (r.buffer) return Promise.resolve(r.buffer);
    return Promise.resolve(r.json);
  } catch (e) {
    return Promise.reject(new ApiError(e.message, { status: e.status || 500, path, payload: e.payload || null, cause: e }));
  }
}

/** Probe the real backend once at boot. */
export async function probe() {
  try {
    const secs = await rawFetch('GET', '/api/sections');
    api.reachable = true;
    api.authRequired = false;
    setStub(false);
    return { live: true, sections: secs };
  } catch (err) {
    // A 401 means a REAL, reachable server that requires a key we don't have
    // (or have wrong) -- fundamentally different from "not answering at
    // all," and must not be papered over with fake stub data: that would
    // silently show an operator plausible-looking numbers that are not this
    // server's real reconstruction. Surface it and let main.js gate the app
    // on a correct key instead.
    if (err instanceof ApiError && err.status === 401) {
      api.reachable = true;
      api.authRequired = true;
      return { live: false, authRequired: true, sections: [] };
    }
    api.reachable = false;
    api.authRequired = false;
    mock.seedHistory();
    setStub(true, err instanceof ApiError && err.status ? `backend ${err.status}` : 'backend unreachable');
    return { live: false, sections: await fromMock('GET', '/api/sections', {}) };
  }
}

/* ------------------------------------------------------------ routes --- */

export const getSections = () => call('GET', '/api/sections');
export const getDrugs = () => call('GET', '/api/drugs');
export const getDefaults = () => call('GET', '/api/defaults');
export const postValidate = (spec) => call('POST', '/api/validate', { body: spec });
export const postRun = (spec) => call('POST', '/api/runs', { body: spec });
export const listRuns = () => call('GET', '/api/runs');
export const getRun = (id) => call('GET', `/api/runs/${encodeURIComponent(id)}`);
export const getStatus = (id) => call('GET', `/api/runs/${encodeURIComponent(id)}/status`);
export const cancelRun = (id) => call('DELETE', `/api/runs/${encodeURIComponent(id)}`);

/* --- bring your own compound -------------------------------------------- */
// A refusal comes back as 422 with `{detail: {code, message, field}}`. That is a
// 4xx, so `call()` throws rather than retrying or falling back to the mock — the
// UI reads `err.payload.detail` and puts the sentence on the field at fault.
// Compiling is never mocked: js/ui/compound.js refuses outright in stub mode
// rather than let a seeded PRNG produce something shaped like chemistry.
export const postCompoundCompile = (body) => call('POST', '/api/compounds/compile', { body });
export const postCompound = (body) => call('POST', '/api/compounds', { body });
export const listCompounds = () => call('GET', '/api/compounds');
export const getCompound = (id) => call('GET', `/api/compounds/${encodeURIComponent(id)}`);
export const deleteCompound = (id) => call('DELETE', `/api/compounds/${encodeURIComponent(id)}`);

/* --- therapy reports ---------------------------------------------------- */
export const postReport = (body) => call('POST', '/api/reports', { body });
export const listReports = () => call('GET', '/api/reports');
export const getReport = (id) => call('GET', `/api/reports/${encodeURIComponent(id)}`);
export const cancelReport = (id) => call('DELETE', `/api/reports/${encodeURIComponent(id)}`);
/** Absolute path of the rendered document — used as an iframe src and as a download. */
export const reportHtmlUrl = (id) => `${api.base}/api/reports/${encodeURIComponent(id)}/html`;

export async function getGeometry(id) {
  const buf = await call('GET', `/api/runs/${encodeURIComponent(id)}/geometry`, { binary: true });
  return decodeTMSG(buf);
}

/**
 * Per-cell mechanism attribution (TMSM). Returns
 * `{n_cells, labels: Uint8Array, aggregate, bytes}` where `labels` is aligned
 * one-to-one with the geometry payload's cell order.
 */
export async function getMechanism(id) {
  const buf = await call('GET', `/api/runs/${encodeURIComponent(id)}/mechanism`, { binary: true });
  return decodeTMSM(buf);
}

/**
 * Drug-concentration time-lapse (TMSL). Returns
 * `{n_frames, interval_s, window_s, frames: [...TMSF], bytes}`.
 */
export async function getTimelapse(id) {
  const buf = await call('GET', `/api/runs/${encodeURIComponent(id)}/timelapse`, { binary: true });
  return decodeTMSL(buf);
}

export async function getField(id, name) {
  const buf = await call('GET', `/api/runs/${encodeURIComponent(id)}/field/${encodeURIComponent(name)}`, { binary: true });
  const f = decodeTMSF(buf);
  f.name = name;
  return f;
}

/* ------------------------------------------------------- binary decode -- */

const MAGIC = (b, s) => String.fromCharCode(b[0], b[1], b[2], b[3]) === s;

/**
 * TMSG: magic(4) | version u32 | n_cells u32 | extent_x f32 | extent_y f32
 *       | x f32[n] | y f32[n] | type u8[n]      (little-endian)
 * Decoded straight into typed-array views over the response buffer — no copy,
 * no per-cell JS object, which is the entire point of the binary payload.
 */
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

/**
 * TMSF: magic(4) | version u32 | rows u32 | cols u32 | dx_um f32
 *       | min f32 | max f32 | values f32[rows*cols]
 */
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

/* -------------------------------------------------- local synth payload -- */

/**
 * Build a real TMSG buffer for an arbitrary cell count, entirely client-side.
 * Used by the perf harness so the 850k benchmark exercises the identical
 * decode -> quantise -> upload path as a live payload.
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

/* ------------------------------------------------- projects / specimens --
 *
 * The multi-patient layer (API_CONTRACT.md "Projects, specimens, and the
 * cross-specimen comparison"). A specimen holds run IDS into app/runs/, so
 * nothing here ever moves a run: attach/detach are pure association.
 *
 * js/mock.js does NOT synthesise these — a study is state the operator built,
 * and inventing one would put fabricated patients on screen. When the backend
 * is unreachable these reject, and ui/projects.js says so instead of drawing
 * a fake project.
 */
const P = (id) => `/api/projects/${encodeURIComponent(id)}`;
const S = (id, sid) => `${P(id)}/specimens/${encodeURIComponent(sid)}`;

export const listProjects = () => call('GET', '/api/projects');
export const createProject = (name, notes = '') =>
  call('POST', '/api/projects', { body: { name, notes } });
export const getProject = (id) => call('GET', P(id));
export const patchProject = (id, patch) => call('PATCH', P(id), { body: patch });
export const deleteProject = (id) => call('DELETE', P(id));

export const addSpecimen = (id, section, label, notes = '') =>
  call('POST', `${P(id)}/specimens`, { body: { section, label, notes } });
export const getSpecimen = (id, sid) => call('GET', S(id, sid));
export const patchSpecimen = (id, sid, patch) => call('PATCH', S(id, sid), { body: patch });
export const deleteSpecimen = (id, sid) => call('DELETE', S(id, sid));

export const attachRun = (id, sid, runId) =>
  call('POST', `${S(id, sid)}/runs`, { body: { run_id: runId } });
export const detachRun = (id, sid, runId) =>
  call('DELETE', `${S(id, sid)}/runs/${encodeURIComponent(runId)}`);

export const compareProject = (id, metric = 'kill_fraction') =>
  call('GET', `${P(id)}/compare?metric=${encodeURIComponent(metric)}`);

export const getCapability = (sectionId) =>
  call('GET', `/api/sections/${encodeURIComponent(sectionId)}/capability`);

/* --- bring your own tissue ---------------------------------------------- */
// One SYNCHRONOUS endpoint, not a job: `POST /api/sections/scan` does the whole
// detect / ingest / validate / capability-report pass inside the request, so
// there is no run id to poll and the browser holds the connection open for the
// duration. A full ingest of a large section is minutes, which is why the caller
// passes its own `timeoutMs` — the module default of 30 s would abort a request
// the server is still working on and leave the operator thinking it failed.
//
// A bad directory is an ANSWER, not an error: the endpoint returns 200 with
// `error` populated. Only a malformed BODY is a 4xx. Neither is retried, and
// neither is ever answered from js/mock.js (see `noFallback`).
export const postSectionScan = (body, { timeoutMs } = {}) =>
  call('POST', '/api/sections/scan', { body, timeoutMs, noFallback: true });

/** Remove an INGESTED section. Built-in sections are refused with a 400. */
export const deleteSection = (sectionId) =>
  call('DELETE', `/api/sections/${encodeURIComponent(sectionId)}`, { noFallback: true });

export { mock };
