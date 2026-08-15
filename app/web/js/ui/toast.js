// toast.js — transient notices. Errors are NEVER auto-dismissed: a failure
// stays on screen until the operator dismisses it.
//
// WHY notifyError() IS NOT A ONE-LINER
// ------------------------------------
// It used to be. It rendered `title` — which at all 30-odd call sites is the
// ROUTE, e.g. `GET /api/runs/r_00078/field/oxygen failed` — as the headline,
// and `err.message + " [HTTP 502]"` underneath. That is an accurate report and
// a useless one: it tells a field tester which line of our code was unhappy and
// nothing about what happened to them or what to do next. It also violates
// docs/UI_ARCHITECTURE.md §5 head-on — "a control is named for what it does,
// not for the endpoint it calls; endpoint names belong in `title`."
//
// The shape here is lifted from js/ui/tissue.js's showRefusal(), which is the
// house style for a failure, and it keeps that panel's load-bearing rule:
//
//   THE SERVER'S SENTENCE AND THE UI'S ADVICE ARE VISUALLY SEPARATE, BECAUSE
//   ONLY THE FIRST IS EVIDENCE.
//
// So an error card is four parts, in this order:
//
//   1. a headline in the operator's terms, chosen from the FAILURE CLASS;
//   2. one line of what to do next, also from the class — procedural only,
//      and it never restates a fact the server did not send;
//   3. the server's own message, verbatim, under a label that says so. It is
//      never clamped, scrolled or truncated: a reason behind a "show more" is
//      a reason that does not get read;
//   4. the route and HTTP status, dim, last. That is the line you paste into a
//      bug report, and it is the only place an internal path appears.
//
// A raw traceback is never rendered. Where the failure is genuinely server-side
// (5xx) the card says where the detail is — the server's own console — rather
// than trying to put a Python stack in a 380 px card.
//
// Call sites did not change. A title that is route-shaped is demoted to (4)
// and the headline comes from the class; a title that is already prose
// ("WebGL2 context could not be created", "Boot failed", "Cannot format: JSON
// does not parse") is kept as the headline, because those sites are describing
// a CLIENT-side fault that no status code classifies.

import { h, clear } from '../util/dom.js';

let host = null;

function ensure() {
  if (!host) {
    host = h('div.toasts', { role: 'log', 'aria-live': 'polite', 'aria-relevant': 'additions' });
    document.body.appendChild(host);
  }
  return host;
}

/**
 * @param {'info'|'ok'|'warn'|'err'} kind
 * @param {string} title
 * @param {string|Node} [detail]
 * @param {number} [ttl] ms; 0 = sticky. Errors default to sticky.
 */
export function toast(kind, title, detail, ttl) {
  const el = h('div.toast' + (kind === 'info' ? '' : '.' + kind),
    h('div.bar'),
    h('div.toast-body',
      h('b', title),
      detail == null ? null : (detail instanceof Node ? detail : h('code', detail))),
    h('button.btn.ghost.icon.sm.x', {
      title: 'Dismiss', 'aria-label': 'Dismiss notification',
      onclick: () => el.remove(),
    }, '×'),
  );
  ensure().appendChild(el);
  const life = ttl === undefined ? (kind === 'err' ? 0 : kind === 'warn' ? 9000 : 4500) : ttl;
  if (life > 0) setTimeout(() => el.remove(), life);
  return el;
}

/* ---------------------------------------------------------- classes ----- */

/**
 * The failure classes this front end can actually produce, each with what
 * happened and what to do about it. Ordered; first match wins.
 *
 * Every remedy is PROCEDURAL — a next action — because the facts belong to the
 * server's own sentence, which is rendered underneath it in full. None of them
 * quotes a number, a timeout, or a limit: those would go stale against
 * js/api.js the moment either changed, and a stale reassurance is worse than
 * none. (js/ui/tissue.js's REMEDIES table is built to the same rule.)
 */
/*
 * `saidLbl` is not decoration. The verbatim block is captioned with WHOSE
 * SENTENCE IT IS, and the three answers are genuinely different: the server
 * found something, the browser gave up, or one of OUR OWN format checks in
 * js/api.js refused a payload. Captioning a decoder's message "What the server
 * said" would attribute this project's words to the backend — the same
 * category of mislabel this whole product exists to avoid, applied to the box
 * whose entire job is to mark which words are evidence.
 */
const CLASSES = [
  {
    test: (e) => e && e.name === 'AbortError',
    head: 'The browser stopped waiting for the server',
    next: 'The request ran past the time this Studio waits for an answer. THE '
      + 'SERVER MAY STILL BE WORKING ON IT — check the run history, or the '
      + 'section picker after an ingest, before you start over. Do not assume '
      + 'nothing happened.',
    saidLbl: 'What the browser reported',
  },
  {
    test: (e) => statusOf(e) === 0,
    head: 'Could not reach the Codon server',
    next: 'The request did not complete: either it never arrived or the answer '
      + 'never came back. Check the server process is still running, and that '
      + 'this page is still pointed at it, then try again.',
    // Nothing came back, so there is nothing the server said. This is the
    // browser's own account of why.
    saidLbl: 'What the browser reported',
  },
  {
    test: (e) => statusOf(e) === 401,
    head: 'This server needs an access key',
    next: 'The key this browser sent was missing or wrong, so nothing was done. '
      // Same wording as ui/apikey.js's gate on purpose: it is the same fact, and
      // the server prints it under exactly that label (app/server/__main__.py).
      + 'Reload the Studio and enter the key printed once at server start, next '
      + 'to "api key".',
    saidLbl: 'What the server said',
  },
  {
    test: (e) => statusOf(e) === 403,
    head: 'The server declined to do that',
    next: 'The request arrived and the server refused to act on it. Its reason '
      + 'is below.',
    saidLbl: 'What the server said',
  },
  {
    test: (e) => statusOf(e) === 404,
    head: 'The server has no record of that',
    next: 'It was deleted, or it was never on this server — a run id from a '
      + 'different Codon instance will not resolve here. Refresh the list you '
      + 'opened it from and see whether it is still there.',
    saidLbl: 'What the server said',
  },
  {
    test: (e) => statusOf(e) === 409,
    head: 'The server is not ready to answer that yet',
    next: 'Something this depends on is still in flight, or it is already in '
      + 'the state you asked for. Its reason is below; nothing was changed.',
    saidLbl: 'What the server said',
  },
  {
    test: (e) => { const s = statusOf(e); return s >= 400 && s < 500; },
    head: 'The server rejected the request',
    next: 'This is a verdict on what was sent, not a fault — nothing was '
      + 'changed. Change what the reason below names, then send it again.',
    saidLbl: 'What the server said',
  },
  {
    // The binary decoders in js/api.js raise ApiError(..., {status: 502}) and
    // never set `path`; every error that came off the wire does. That makes an
    // empty path the reliable discriminator, rather than matching on the
    // wording of the message, which is ours to change.
    test: (e) => statusOf(e) === 502 && !(e && e.path),
    head: 'The server’s answer was not in a form this Studio can read',
    next: 'Nothing was drawn, deliberately: rendering a payload that failed its '
      + 'own format check would put a picture on screen that no longer matches '
      + 'the data behind it. This normally means the page and the backend are '
      + 'from different builds — reload the page, and if it persists restart '
      + 'the server from the same checkout this page was served by.',
    // OUR check, OUR sentence. The server sent bytes and said nothing.
    saidLbl: 'What failed the format check',
  },
  {
    test: (e) => statusOf(e) >= 500,
    head: 'The server hit a fault handling that',
    next: 'This is a fault on the server side, not something wrong with what '
      + 'you sent. The full detail is in the console of the process running '
      + 'Codon, not here — its short message is below.',
    saidLbl: 'What the server said',
  },
];

function statusOf(err) {
  const s = err && err.status;
  return Number.isFinite(s) ? s : 0;
}

/**
 * Is this error about TALKING TO THE SERVER at all?
 *
 * It has to be asked before the table above is consulted, because the table's
 * `status === 0` row would otherwise swallow every client-side fault that has
 * no status — `notifyError('Cannot format: JSON does not parse', syntaxError)`
 * would be answered with "Could not reach the Codon server". That is not a
 * cosmetic mislabel: it sends the operator to restart a backend that is fine.
 *
 * js/api.js guarantees the two shapes below — every transport failure leaves it
 * as an `ApiError` carrying a status (0 for "never completed"), and an abort
 * keeps the DOM's own `AbortError`. Anything else is ours, and degrades to "no
 * remedy line", never to a WRONG one.
 */
const isTransport = (e) => !!e && (e.name === 'ApiError' || e.name === 'AbortError');

/** A title that is really a route: `GET /api/…`, or anything naming one. */
const ROUTE_TITLE = /^(GET|POST|PUT|PATCH|DELETE)\s|\/api\//;

/**
 * The route + status line. Prefers `err.path` (what api.js actually called)
 * and falls back to the caller's route-shaped title.
 */
function techLine(title, err) {
  const bits = [];
  const path = err && err.path;
  if (path) bits.push(path);
  else if (ROUTE_TITLE.test(title)) bits.push(String(title).replace(/\s+failed$/i, ''));
  const s = statusOf(err);
  if (s) bits.push('HTTP ' + s);
  return bits.join('  ·  ');
}

/**
 * Render a failure as an answer.
 *
 * @param {string} title  what the CALLER was doing — usually a route.
 * @param {Error|import('../api.js').ApiError} err
 */
export function notifyError(title, err) {
  const said = (err && (err.message || String(err))) || '';
  const routeShaped = ROUTE_TITLE.test(String(title || ''));
  const hit = isTransport(err)
    ? CLASSES.find((c) => { try { return c.test(err); } catch { return false; } })
    : null;

  // A prose title is a client-side fault the caller already described better
  // than any status code could ("WebGL2 context could not be created"). Keep
  // it. Only a route-shaped title gets replaced by the class headline — and
  // only when there is a class to replace it with.
  const headline = (routeShaped && hit) ? hit.head : (title || 'Something failed');
  const remedy = hit ? hit.next : '';
  // Unclassified: this is our own code failing, so it is captioned neutrally
  // rather than attributed to anyone.
  const saidLbl = hit ? hit.saidLbl : 'Detail';

  const tech = techLine(String(title || ''), err);

  const body = h('div.toast-detail',
    remedy ? h('div.toast-remedy', remedy) : null,
    said
      ? h('div.toast-said',
        h('div.toast-said-lbl', saidLbl),
        // Verbatim and IN FULL — no clamp, no inner scroller, no "show more".
        h('code', said))
      : null,
    tech ? h('div.toast-tech', { title: 'The route and status. This is the line to quote in a bug report.' }, tech) : null,
  );

  // The whole object still goes to the console, where a developer can expand
  // it — that is where a stack belongs, and it is the reason none is rendered.
  // eslint-disable-next-line no-console
  console.error(title, err);
  return toast('err', headline, body);
}

export function clearToasts() { if (host) clear(host); }
