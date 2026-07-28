// api.js — the ONLY place that talks to app/server.
//
// Routes consumed (all of API_CONTRACT.md v1):
//   GET    /api/sections
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

export const api = {
  base: '',
  stub: false,
  /** null until probe() has run. */
  reachable: null,
  lastLatencyMs: null,
  onModeChange: null,
};

function setStub(on, why) {
  if (api.stub === on) return;
  api.stub = on;
  api.stubReason = why || null;
  if (api.onModeChange) api.onModeChange(api);
}

async function rawFetch(method, path, { body, binary, signal } = {}) {
  const ctrl = new AbortController();
  const timer = setTimeout(() => ctrl.abort(), TIMEOUT_MS);
  if (signal) signal.addEventListener('abort', () => ctrl.abort(), { once: true });
  const t0 = performance.now();
  try {
    const res = await fetch(api.base + path, {
      method,
      headers: body !== undefined
        ? { 'content-type': 'application/json', accept: binary ? 'application/octet-stream' : 'application/json' }
        : { accept: binary ? 'application/octet-stream' : 'application/json' },
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
    setStub(false);
    return { live: true, sections: secs };
  } catch (err) {
    api.reachable = false;
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

export async function getGeometry(id) {
  const buf = await call('GET', `/api/runs/${encodeURIComponent(id)}/geometry`, { binary: true });
  return decodeTMSG(buf);
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

export { mock };
