// net.js — the only place this page talks to the server.
//
// Deliberately NOT js/api.js. That module is the Studio's client and it has
// one behaviour this page must never have: when the backend does not answer it
// falls back to js/mock.js and serves synthetic data with `api.stub = true`.
// Synthetic sections are the right call for a Studio demo and a catastrophe
// here — this page's entire job is to tell you what is REALLY on the server's
// disk, and a fabricated directory listing is worse than an error message.
// api.js is also not ours to edit.
//
// Routes used, and only these three:
//   GET  /api/browse                     list the names in one directory
//   POST /api/sections/scan              inspect a directory (read only)
//   POST /api/sections/scan {register}   register it as a section
//
// The second and third are the SAME endpoint the original /ingest form and
// `vivum ingest` use. Nothing on this page is a new capability: it is a way to
// find the directory without typing it.

/** The Studio's own storage key, on purpose: a person who has already entered
 *  the access key for /studio must not be asked for it again here, and one
 *  entered here works there. js/api.js owns this name; we only read and write
 *  the same slot. */
const API_KEY_STORAGE = 'tmesim_api_key';

export const net = {
  /** '' when the server needs no key (the loopback default). */
  apiKey: readKey(),
  /** Set once a call comes back 401: a real server asking for a key. */
  authRequired: false,
};

function readKey() {
  try {
    return localStorage.getItem(API_KEY_STORAGE)
      || sessionStorage.getItem(API_KEY_STORAGE) || '';
  } catch {
    // Private mode, or site data blocked. The key still works for this page
    // load; it simply will not survive a reload.
    return '';
  }
}

/** Persist a key (or clear it with ''). `persist` false keeps it to this tab. */
export function setApiKey(key, persist = true) {
  net.apiKey = key || '';
  try {
    if (net.apiKey && persist) {
      localStorage.setItem(API_KEY_STORAGE, net.apiKey);
      sessionStorage.removeItem(API_KEY_STORAGE);
    } else if (net.apiKey) {
      sessionStorage.setItem(API_KEY_STORAGE, net.apiKey);
      localStorage.removeItem(API_KEY_STORAGE);
    } else {
      localStorage.removeItem(API_KEY_STORAGE);
      sessionStorage.removeItem(API_KEY_STORAGE);
    }
  } catch { /* storage unavailable; the in-memory key still works */ }
}

export class ServerRefusal extends Error {
  constructor(message, { status = 0, code = '', payload = null } = {}) {
    super(message);
    this.name = 'ServerRefusal';
    this.status = status;
    this.code = code;
    this.payload = payload;
  }
}

/** The server's own sentence, whichever shape it refused in.
 *
 *  This backend refuses in two shapes: `detail` as a string (FastAPI's
 *  default, and what the body validators use) and `detail` as
 *  `{code, message, ...}` (the confinement 403, every compound refusal). Both
 *  are authored FOR THE READER — the 403 names the path you sent, the
 *  environment variable and how to add a root — so the one thing this must not
 *  do is replace them with wording of its own. It only finds the sentence. */
function refusalFrom(status, body) {
  const detail = body && body.detail !== undefined ? body.detail : body;
  if (detail && typeof detail === 'object' && typeof detail.message === 'string') {
    return new ServerRefusal(detail.message, {
      status, code: detail.code || '', payload: detail,
    });
  }
  if (typeof detail === 'string' && detail) {
    return new ServerRefusal(detail, { status, payload: body });
  }
  return new ServerRefusal(
    'the server refused this request with status ' + status
    + ' and no explanation.', { status, payload: body });
}

async function call(path, { method = 'GET', body = null } = {}) {
  const headers = {};
  if (body !== null) headers['Content-Type'] = 'application/json';
  if (net.apiKey) headers['X-API-Key'] = net.apiKey;

  let resp;
  try {
    resp = await fetch(path, {
      method, headers,
      body: body === null ? undefined : JSON.stringify(body),
    });
  } catch (err) {
    throw new ServerRefusal(
      'the server did not answer. It may have stopped, or this page may have '
      + 'been left open while it restarted.', { status: 0, code: 'unreachable' });
  }

  let payload = null;
  try { payload = await resp.json(); } catch { payload = null; }

  if (resp.status === 401) {
    net.authRequired = true;
    throw new ServerRefusal(
      'this server was started with an access key, so it needs one on every '
      + 'request. Enter it below.', { status: 401, code: 'auth_required' });
  }
  if (!resp.ok) throw refusalFrom(resp.status, payload);
  net.authRequired = false;
  return payload;
}

/** GET /api/browse — one directory's names, or the starting points. */
export function browse({ path = '', q = '', showHidden = false } = {}) {
  const qs = new URLSearchParams();
  if (path) qs.set('path', path);
  if (q) qs.set('q', q);
  if (showHidden) qs.set('show_hidden', '1');
  const tail = qs.toString();
  return call('/api/browse' + (tail ? '?' + tail : ''));
}

/** POST /api/sections/scan — read only. Milliseconds unless `measureRates`. */
export function inspect(path, { measureRates = false } = {}) {
  return call('/api/sections/scan', {
    method: 'POST',
    body: { path, measure_detection_rates: !!measureRates },
  });
}

/** POST /api/sections/scan with register:true — one action, from the browser.
 *
 *  A full ingest runs first and registration is REFUSED, not degraded, if the
 *  sample does not load and validate. `{ok:false, error}` comes back with a
 *  200 for exactly that case, so a caller must read `ok` and not the status. */
export function register(path, { sectionId = '', label = '' } = {}) {
  const body = { path, register: true };
  if (sectionId) body.section_id = sectionId;
  if (label) body.label = label;
  return call('/api/sections/scan', { method: 'POST', body });
}
