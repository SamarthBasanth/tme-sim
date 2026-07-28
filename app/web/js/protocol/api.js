// protocol/api.js — the ONLY place that talks to the protocol routes.
//
// FROZEN INTERFACE consumed by app/web/js/ui/protocol.js. The parser/executor
// side (`tmesim/protocol/` + `app/server/`) implements PROTOCOL_LANGUAGE.md v1;
// this is the wire shape the editor codes to.
//
//   POST /api/protocol/validate   {source}
//     -> { ok, errors[], warnings[], symbols? }
//        diagnostic := { line, col?, end_col?, message, severity?, code?, remedy?,
//                        remedy_line?, remedy_text? }
//        `line` is 1-BASED. `col` / `end_col` are 0-BASED character offsets in
//        that line, half-open. Both optional; without them the whole line is
//        marked. `remedy_text` is TPL the editor may offer to insert verbatim.
//
//   POST /api/protocol/plan       {source}
//     -> { ok, run_count, max_runs, estimated_cost_s, cost_note,
//          axes:[{name, unit, values[]}], runs:[{index, coords[], label, params, spec, est_s}],
//          refusal: {code, message, remedy, remedy_text, line} | null,
//          errors[], warnings[], caveats[], assumptions[] }
//
//   POST /api/protocol/upload     multipart/form-data, field `file`
//     -> { file:{name, size, format}, drug:{name, source, fields:{...}},
//          warnings[], caveats[] }
//        field := { value, unit, provenance, source, note }
//        provenance ∈ measured | computed | user-supplied | flagged-default
//
//   POST /api/protocol/run        {source, confirm_run_count}
//     -> { sweep_id, run_count, runs:[{index, run_id, status}] }
//        The server MUST refuse if confirm_run_count disagrees with its own
//        expansion: the whole point of plan-before-run is that the operator
//        confirmed a specific number.
//
//   GET  /api/protocol/sweeps/{id}   -> sweep payload (see stub.js for shape)
//   DELETE /api/protocol/sweeps/{id} -> cancel
//
// POLICY: exactly one retry on a network error or 5xx, 400 ms backoff. A 404 or
// 501 means the executor has not landed in this deployment yet — we degrade to
// the local stub, flip `protoApi.stub`, and every surface says STUB in words.

import { api as coreApi } from '../api.js';
import * as stub from './stub.js';

const RETRY_MS = 400;
const TIMEOUT_MS = 60000;

export class ProtocolError extends Error {
  constructor(message, { status = 0, path = '', payload = null } = {}) {
    super(message);
    this.name = 'ProtocolError';
    this.status = status;
    this.path = path;
    this.payload = payload;
  }
}

export const protoApi = {
  stub: false,
  reason: null,
  probed: false,
  onModeChange: null,
};

function setStub(on, why) {
  if (protoApi.stub === on && protoApi.reason === why) return;
  protoApi.stub = on;
  protoApi.reason = on ? (why || 'protocol executor unreachable') : null;
  if (protoApi.onModeChange) protoApi.onModeChange(protoApi);
}

function detailOf(payload) {
  if (!payload) return '';
  if (typeof payload.detail === 'string') return payload.detail;
  if (Array.isArray(payload.detail)) {
    return payload.detail.map((d) => (d && d.msg ? `${(d.loc || []).slice(1).join('.')}: ${d.msg}` : JSON.stringify(d))).join('; ');
  }
  if (typeof payload.error === 'string') return payload.error;
  return '';
}

async function rawFetch(method, path, { body, form } = {}) {
  const ctrl = new AbortController();
  const timer = setTimeout(() => ctrl.abort(), TIMEOUT_MS);
  try {
    const init = { method, signal: ctrl.signal, cache: 'no-store', headers: { accept: 'application/json' } };
    if (form) {
      init.body = form;                       // let the browser set the boundary
    } else if (body !== undefined) {
      init.headers['content-type'] = 'application/json';
      init.body = JSON.stringify(body);
    }
    const res = await fetch(coreApi.base + path, init);
    if (!res.ok) {
      let payload = null;
      let text = '';
      try {
        const ct = res.headers.get('content-type') || '';
        if (ct.includes('json')) payload = await res.json();
        else text = (await res.text()).slice(0, 400);
      } catch { /* the status is still the story */ }
      throw new ProtocolError(detailOf(payload) || text || `${res.status} ${res.statusText}`,
        { status: res.status, path, payload });
    }
    if (res.status === 204) return null;
    return await res.json();
  } finally {
    clearTimeout(timer);
  }
}

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

/**
 * One call, exactly one retry on network / 5xx, then a documented degradation.
 * A 4xx that is NOT 404 is a real answer from the executor (a refusal, a bad
 * request) and is raised so the caller can render it — never swallowed.
 */
async function call(method, path, opts = {}, stubFn) {
  if (protoApi.stub || coreApi.stub) {
    if (!protoApi.stub) setStub(true, coreApi.stub ? 'backend unreachable' : 'protocol routes absent');
    return stubFn();
  }
  for (let attempt = 0; attempt < 2; attempt++) {
    try {
      const out = await rawFetch(method, path, opts);
      setStub(false);
      return out;
    } catch (err) {
      const status = err instanceof ProtocolError ? err.status : 0;
      if ((status === 0 || status >= 500) && attempt === 0) { await sleep(RETRY_MS); continue; }
      // 404 = no such route. 405 = the path resolves to the static-file
      // catch-all, which answers GET only — i.e. the executor is not mounted
      // yet. 501 = mounted but unimplemented. All three mean "not deployed",
      // and none of them is an answer about the operator's program.
      if (status === 404 || status === 405 || status === 501) {
        setStub(true, `protocol routes not deployed (HTTP ${status} on ${path})`);
        return stubFn();
      }
      if (status === 0 || status >= 500) {
        setStub(true, status === 0 ? 'protocol executor unreachable' : `protocol executor ${status}`);
        return stubFn();
      }
      throw err;
    }
  }
  throw new ProtocolError('unreachable', { path });
}

/* ------------------------------------------------------------- routes --- */

export const validate = (source) =>
  call('POST', '/api/protocol/validate', { body: { source } }, () => stub.validate(source));

export const plan = (source) =>
  call('POST', '/api/protocol/plan', { body: { source } }, () => stub.plan(source));

export function upload(file) {
  const form = new FormData();
  form.append('file', file, file.name);
  return call('POST', '/api/protocol/upload', { form }, () => stub.upload(file));
}

export const start = (source, confirmRunCount) =>
  call('POST', '/api/protocol/run', { body: { source, confirm_run_count: confirmRunCount } },
    () => stub.start(source, confirmRunCount));

export const sweep = (id) =>
  call('GET', `/api/protocol/sweeps/${encodeURIComponent(id)}`, {}, () => stub.sweep(id));

export const cancel = (id) =>
  call('DELETE', `/api/protocol/sweeps/${encodeURIComponent(id)}`, {}, () => stub.cancel(id));
