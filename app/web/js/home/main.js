// main.js — the workspace home page's boot (docs/WORKSPACE_CONTRACT.md).
//
// ONE CALL, FOUR STATES.
// ----------------------
// The whole first paint comes from a single `GET /api/workspace`, because four
// round-trips on a cold Render instance is four chances to look broken (§1).
// What that one call can produce is four genuinely different screens, and this
// module's job is to keep them different:
//
//   LOADING      the call is in flight. This instance sleeps, so a first paint
//                is realistically 20-35 s. The page shows a true elapsed count
//                and says what is being waited on, escalating as it goes.
//                Never a spinner that is really an error.
//   EMPTY        the server answered and holds zero projects. The most
//                important screen on the page: it has to teach the product and
//                end on a first run. Owned by projectGrid.js.
//   POPULATED    the server answered with projects. Grid, search, sort, and
//                the cross-project recent-runs strip.
//   UNREACHABLE  nothing answered, or the answer is js/mock.js's. This is NOT
//                "you have no projects" — it is "nothing was answered" — and
//                it is never allowed to look like the empty state.
//
// A fifth branch that is not a state: a 401 means a real, reachable server that
// wants a key, and it gets the SAME gate the Studio uses (ui/apikey.js).
//
// WHY THE SIBLING MODULES ARE IMPORTED DYNAMICALLY.
// ------------------------------------------------
// nav.js, projectCard.js, recentRuns.js and createProject.js are owned by other
// agents and land separately. A static `import` of a module that is not there
// throws at LOAD time and the result is a blank page — the worst possible
// failure for a landing page. Each is therefore behind a dynamic `import()` in
// a try/catch, and a miss degrades to a named, visible placeholder that says
// which file is absent and what still works. The EXPORT NAMES below are not
// guesses: each was read out of the landed module with grep, and the signatures
// they are called with are the ones in those modules' own JSDoc.
//
// This module owns exactly three session-wide things, because two modules
// writing any of them is the failure docs/UI_ARCHITECTURE.md §2.1 records:
// the `data-theme` stamp at boot, the `data-mode` stub latch, and the single
// workspace fetch.

import { h, $, fill, clear } from '../util/dom.js';
import * as api from '../api.js';
import { dur, int } from '../util/fmt.js';
import { createProjectGrid } from './projectGrid.js';
import { createApiKeyGate } from '../ui/apikey.js';

/* ============================================================== theme ==== */

/* The SAME localStorage key the Studio's js/main.js writes, so a user who
   picked light there lands on a light home page (contract §7). nav.js exports
   its own initTheme/toggleTheme against this identical key; this copy exists
   because it runs SYNCHRONOUSLY at module top, before the dynamic import of
   nav.js has resolved, and the theme must be stamped before the first paint
   rather than one microtask into it. Once nav.js mounts, its button owns the
   toggle and this fallback button has been removed with the fallback bar. */
const THEME_KEY = 'tmesim.theme';

function applyTheme(t) {
  document.documentElement.setAttribute('data-theme', t);
  try { localStorage.setItem(THEME_KEY, t); } catch { /* private mode */ }
  const b = $('#home-theme');
  if (b) {
    b.textContent = 'theme: ' + t;
    b.setAttribute('aria-pressed', String(t === 'light'));
  }
}

function initTheme() {
  let t = null;
  try { t = localStorage.getItem(THEME_KEY); } catch { /* private mode */ }
  if (t !== 'light' && t !== 'dark') {
    t = window.matchMedia('(prefers-color-scheme: light)').matches ? 'light' : 'dark';
  }
  applyTheme(t);
}

function toggleTheme() {
  applyTheme(document.documentElement.getAttribute('data-theme') === 'light' ? 'dark' : 'light');
}

/* ================================================== module resolution ==== */

function errText(err) {
  if (!err) return 'no reason given';
  if (err.message) return String(err.message);
  return String(err);
}

/**
 * Load a sibling module and pick a named export.
 * Returns {fn, mod, why} — `fn` null means "use the placeholder", and `why` is
 * a sentence fit to show a reader, not a stack trace.
 */
async function resolve(path, names) {
  let mod = null;
  try {
    mod = await import(path);
  } catch (err) {
    return { fn: null, mod: null, why: `${path} did not load: ${errText(err)}` };
  }
  for (const n of names) {
    if (mod && typeof mod[n] === 'function') return { fn: mod[n], name: n, mod, why: '' };
  }
  const found = mod ? Object.keys(mod).filter((k) => typeof mod[k] === 'function') : [];
  return {
    fn: null,
    mod,
    why: `${path} loaded, but exports no function named ${names.join(' / ')}`
      + (found.length ? ` (it exports: ${found.join(', ')}).` : ' (it exports no functions).'),
  };
}

/** A node out of whatever shape a factory returned. `h()` and
 *  createProjectCard/createRecentRuns return an element directly; the house
 *  `createX()` convention returns an object, and createNav's key for it is
 *  `el`. All three, or null. */
function nodeOf(r) {
  if (!r) return null;
  if (r.nodeType) return r;
  if (r.el && r.el.nodeType) return r.el;
  if (r.node && r.node.nodeType) return r.node;
  return null;
}

function missNote(what, why, stillWorks) {
  return h('div.mount-miss',
    h('b', what + ' '),
    stillWorks,
    why ? h('span.mm-why', why) : null);
}

/* ================================================================ dom ==== */

const shell = $('#home-shell');
const statusHost = $('#home-status');
const bodyHost = $('#home-body');

/* =========================================================== the wait ==== */

/* A COLD START IS THE NORMAL CASE HERE, NOT THE EXCEPTION.
 *
 * Measured against the real deployment, on an instance that had gone to sleep:
 *   GET /healthz         200 in 22.4 s      the cheapest route on the server,
 *                                           returning a literal {"ok": true}.
 *                                           That is the CONTAINER waking, not
 *                                           the API being slow.
 *   GET /api/sections    12.7 s first call, 0.08 s second, at the 0.1 CPU a
 *                                           free instance gets (recorded in
 *                                           api.js's PROBE_TIMEOUT_MS comment)
 * So a first paint for anyone arriving at a link is realistically 20-35 s, and
 * a bare spinner reads as broken somewhere past five of them. Every phase below
 * therefore says what is being waited on and why, and none of them says or
 * implies that anything has failed until something has.
 *
 * THE BUDGET, and why it is not the 120 s PROBE_TIMEOUT_MS. getWorkspace() goes
 * through api.js's call(), which uses TIMEOUT_MS = 30000 per attempt with one
 * automatic retry RETRY_MS = 400 apart; PROBE_TIMEOUT_MS applies only to
 * probe()'s own /api/sections call, which this page does not make. That is NOT
 * a 30 s deadline on the wake, and the difference matters: if the container is
 * still waking when the first attempt aborts at 30 s, the second attempt starts
 * at 30.4 s against a server that is by then finishing its wake, and it has its
 * own 30 s. The wall-clock tolerance for a cold start is therefore ~60 s, which
 * clears the measured 20-35 s with room. This page adds no deadline of its own
 * and will wait out both attempts. */
const PHASES = [
  () => ['Asking this server for its projects and its most recent runs — one call to ',
    h('span.mono', 'GET /api/workspace'), '.'],
  () => ['This is normal on a first visit and nothing has failed. The instance sleeps when '
    + 'nothing is using it, and the first request after a sleep has to wake the container '
    + 'before any route can answer at all. Measured on this deployment: ',
    h('span.mono', '/healthz'), ', the cheapest route on the server, took ',
    h('b', '22.4 s'), ' on a sleeping instance and is immediate once it is awake. Expect '
    + 'roughly 20 to 35 s here, once, and then a fast page for the rest of the session.'],
  () => ['Past 20 s, which is where the measured cold start lands. Still inside the first '
    + 'attempt, which has a 30 s limit. The Studio does not depend on this page — it is '
    + 'linked in the bar above and it reads its own routes.'],
  () => ['The first attempt reached its 30 s limit. js/api.js waits 400 ms and tries once '
    + 'more, and that second attempt is usually the one that lands: the container that was '
    + 'still waking during the first attempt has finished waking by the time the second one '
    + 'arrives. It gets its own 30 s. If both miss — about 60 s in total — this page falls '
    + 'back to example data from js/mock.js and labels every surface as such, rather than '
    + 'show fabricated work without saying so.'],
  () => ['Both attempts have now passed their 30 s limit, so the labelled fallback should '
    + 'be happening. If this message is still here in a few seconds, reload the page: that '
    + 'is a failure this page did not catch, and it is not a statement about your work.'],
];

function phaseOf(s) {
  // The explanation has to be on screen before a spinner starts reading as a
  // fault, which is at about five seconds.
  if (s < 4) return 0;
  if (s < 20) return 1;
  if (s < 30) return 2;
  if (s < 61) return 3;
  return 4;
}

/** Render the live LOADING state. Returns a stop function. */
function startLoading() {
  const t0 = performance.now();
  // aria-hidden: this ticks four times a second, and a live region that
  // re-announces a counter is a live region nobody can use. The WHY below it
  // changes only at a phase boundary, and that is what gets announced.
  const elapsed = h('span.mono.boot-el', { 'aria-hidden': 'true' }, '0.0s elapsed');
  const why = h('p.boot-why');

  fill(statusHost,
    h('div.boot',
      h('div.boot-head',
        h('span.led.busy', { 'aria-hidden': 'true' }),
        h('h2.boot-title', 'Reaching the server'),
        elapsed),
      h('div.progress.indet', { role: 'progressbar', 'aria-label': 'Contacting the server' },
        h('i')),
      why,
      h('div.boot-out',
        h('a.btn.ghost.sm', { href: '/studio' }, 'Open the Studio instead'))));
  statusHost.setAttribute('aria-busy', 'true');

  let phase = -1;
  const tick = () => {
    const s = (performance.now() - t0) / 1000;
    elapsed.textContent = dur(s) + ' elapsed';
    const p = phaseOf(s);
    if (p !== phase) { phase = p; fill(why, PHASES[p]()); }
  };
  tick();
  const timer = setInterval(tick, 250);
  return () => { clearInterval(timer); statusHost.setAttribute('aria-busy', 'false'); };
}

/* ============================================================ banners ==== */

function banner(kind, title, chip, bodyNodes, footNodes) {
  return h('div.hbanner.' + kind,
    h('div.hb-head', h('h2', title), chip || null),
    h('div.hb-body', bodyNodes),
    footNodes && footNodes.length ? h('div.hb-foot', footNodes) : null);
}

function retryBtn() {
  return h('button.btn.primary', { type: 'button', onclick: () => location.reload() },
    'Try again');
}

function studioBtn(label) {
  return h('a.btn', { href: '/studio' }, label || 'Open the Studio');
}

/** api.js fell back to js/mock.js. Labelled loudly, exactly as the Studio
 *  labels it (contract §6). */
function stubBanner() {
  // api.api, not api: `import * as api` is the MODULE namespace, and the live
  // flags (stub, stubReason, reachable, apiKey) hang off the `api` object that
  // module exports. js/main.js reads them as api.api.stub for the same reason.
  const a = api.api || {};
  const why = a.stubReason ? String(a.stubReason) : 'the backend did not answer';
  return banner('stub', 'Example data — this is not this server’s workspace',
    h('span.chip.synthetic', { title: 'Fabricated in this browser by js/mock.js' }, 'synthetic'),
    [
      h('p', 'The workspace call did not get an answer (', h('b', why),
        '), so js/api.js fell back to js/mock.js.'),
      h('p', 'Every project, run and count on this page below this banner was generated in '
        + 'this browser from a seeded pseudo-random number generator. None of it is stored '
        + 'anywhere, none of it came off a tissue section, and opening one of these '
        + 'projects in the Studio will not find it.'),
      h('p', 'If this was a cold start, both 30 s attempts missed a container that was '
        + 'still waking — measured at 22.4 s on this deployment, and occasionally longer. '
        + 'The wake is finished by now, so trying again is the fix and it will be fast.'),
    ],
    [retryBtn(), studioBtn()]);
}

/** The SERVER said its own workspace is stub data. A different claim from the
 *  one above, and it deserves different words. */
function serverStubBanner() {
  return banner('stub', 'This server reports its workspace as stub data',
    h('span.chip.synthetic', { title: 'server.stub === true in the /api/workspace payload' },
      'synthetic'),
    [
      h('p', 'The call reached the server and the server answered, but its payload carries ',
        h('span.mono', '"server": {"stub": true}'),
        ' — it is telling this page that the projects and runs below are placeholders it '
        + 'generated, not work it has stored.'),
    ],
    [studioBtn()]);
}

function contractBanner(got) {
  return banner('warn', 'This server announces a different workspace contract', null,
    [
      h('p', 'This page is built against contract ', h('span.mono', 'v1'),
        ' and the server announced ', h('span.mono', String(got)),
        '. Fields this page reads may be absent or may mean something else. Anything '
        + 'missing is shown as missing below, never filled in with a guess.'),
    ],
    []);
}

/**
 * The backend did not answer AND js/mock.js has no /api/workspace route
 * either, so the call rejected from the fallback rather than from the server.
 *
 * THIS BANNER EXISTS BECAUSE THE OBVIOUS CODE GETS IT WRONG. api.js's call()
 * wraps a network failure as ApiError(status 0), retries, latches stub mode,
 * and then asks js/mock.js — whose router throws `mock: no route for GET
 * /api/workspace` with `status: 404`. So the error that arrives here carries a
 * 404 that NO SERVER EVER SENT, and reporting it as "the server answered 404"
 * sends the reader to look for a missing route on a server that in fact said
 * nothing at all. The stub latch is what tells the two apart, and it is
 * checked before the status is.
 */
function stubRefusedBanner(err) {
  const a = api.api || {};
  const why = a.stubReason ? String(a.stubReason) : 'the backend did not answer';
  const detail = errText(err);
  return banner('crit', 'The workspace call did not come back', null,
    [
      h('p', 'Two things failed here, and only the first is about this server. ',
        h('b', why), ' — js/api.js made two attempts with a 30 s limit each, 400 ms '
        + 'apart, and then latched this session into stub mode.'),
      h('p', 'It then asked js/mock.js for the workspace instead, and js/mock.js carries no ',
        h('span.mono', '/api/workspace'),
        ' route, so there is no fabricated project list to show either. That is the better '
        + 'of the two outcomes: nothing on this page has been invented to fill the gap.'),
      h('p', 'Any 404 in the detail below is the local mock refusing, not this server '
        + 'answering. The server said nothing at all.'),
      h('p', 'If this instance sleeps when idle, a reload is often all it takes: the wake '
        + 'was measured at 22.4 s on this deployment, the first request pays it, and the '
        + 'next one lands immediately.'),
      detail ? h('span.hb-detail', detail) : null,
    ],
    [retryBtn(), studioBtn(), h('a.btn.ghost', { href: '/api/health' }, 'Check server health')]);
}

/** Nothing came back. The message is specific to WHAT failed, because "not
 *  reachable" and "answered 404" send a reader to completely different places. */
function unreachableBanner(err) {
  const status = err && typeof err.status === 'number' ? err.status : null;
  const detail = errText(err);
  const lines = [];

  if (status === 404 || status === 405) {
    lines.push(h('p', 'The server answered ', h('b', String(status)), ' for ',
      h('span.mono', 'GET /api/workspace'),
      '. That route is part of the workspace contract (docs/WORKSPACE_CONTRACT.md §1) and '
      + 'this deployment does not serve it. That is a missing route, not a broken server.'));
    lines.push(h('p', 'Nothing else is affected: the Studio reads its own routes and works '
      + 'exactly as it always has, project store included.'));
  } else if (status === 403) {
    lines.push(h('p', 'The server answered ', h('b', '403'), ' for ',
      h('span.mono', 'GET /api/workspace'),
      '. The key this browser is sending reached a live server and was refused for this '
      + 'route specifically. Whoever runs this instance can say why; there is no '
      + 'self-service sign-up.'));
  } else if (status && status >= 500) {
    lines.push(h('p', 'The server answered ', h('b', String(status)),
      '. js/api.js already retried once, 400 ms apart, and got the same thing. The server '
      + 'is running and failing, rather than absent.'));
  } else if (status === 0 || status === null) {
    lines.push(h('p', 'No HTTP status came back at all: the request never reached a server '
      + 'that answered. js/api.js made two attempts with a 30 s limit each, 400 ms apart. '
      + 'That pattern is DNS, the network between here and the host, or an instance that '
      + 'is not running.'));
    lines.push(h('p', 'If this instance sleeps when idle, a reload is often all it takes: '
      + 'the wake was measured at 22.4 s on this deployment, the first request pays it, '
      + 'and the next one lands immediately.'));
  } else {
    lines.push(h('p', 'The server answered ', h('b', String(status)), ' for ',
      h('span.mono', 'GET /api/workspace'), '.'));
  }

  lines.push(h('p', 'What this page will not do is guess. There is no project list below, '
    + 'because none was received — not because there are none.'));
  if (detail) lines.push(h('span.hb-detail', detail));

  return banner('crit', 'The workspace call did not come back', null, lines,
    [retryBtn(), studioBtn(), h('a.btn.ghost', { href: '/api/health' }, 'Check server health')]);
}

/* =============================================================== gate ==== */

/** A 401 is a live server asking for a key. Same gate as the Studio, driven
 *  here so that one attempt is still exactly one /api/workspace call. */
function showGate(onAuthed) {
  // A 401 is proof of a live server, so the page is not in stub mode and must
  // not wear the synthetic chrome while the gate is up.
  document.documentElement.setAttribute('data-mode', 'live');
  if (api.api) api.api.reachable = true;
  paintBackend();
  fill(statusHost, banner('info', 'This deployment needs an access key', null,
    [h('p', 'The server answered ', h('b', '401'), ' for ', h('span.mono', '/api/workspace'),
      '. It is there, it is answering, and it wants a key before it will say what work it '
      + 'holds.')],
    []));

  const gate = createApiKeyGate({
    onSubmit: async (key, remember) => {
      api.setApiKey(key, remember !== false);
      try {
        const ws = await api.getWorkspace();
        gate.remove();
        onAuthed({ ok: true, ws });
        return true;
      } catch (err) {
        if (err && err.status === 401) {
          // Confirmed wrong by the server. Clear it from both stores so a
          // reload does not silently retry a key already known to be rejected.
          api.setApiKey('');
          return false;
        }
        gate.remove();
        onAuthed({ ok: false, err });
        return true;
      }
    },
  });
  document.body.appendChild(gate.node);
}

/* ============================================================== open ===== */

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

/**
 * Stamp recency, then go.
 *
 * Recency is an ANNOTATION: a failure to record it never blocks the open, and
 * neither does a slow one. The POST gets 600 ms of the reader's time and not a
 * millisecond more, for two reasons — on a warm server it is an 0.08 s-class
 * call and finishes inside that, and on a cold one the Studio stamps the same
 * project itself on arrival (contract §8.1), so waiting here buys nothing that
 * is not about to happen anyway.
 */
async function stampOpen(id) {
  if (!id) return;
  try {
    await Promise.race([api.openProject(id).catch(() => null), sleep(600)]);
  } catch { /* recency is an annotation; never block an open on it */ }
}

async function openProject(p) {
  const id = p && p.id;
  if (!id) return;
  await stampOpen(id);
  location.assign('/studio?project=' + encodeURIComponent(String(id)));
}

/* ============================================================== boot ===== */

initTheme();
const themeBtn = $('#home-theme');
if (themeBtn) themeBtn.addEventListener('click', toggleTheme);

/* nav.js builds its OWN <header class="topbar nav-bar" role="banner"> and
   inserts it as the first child of `mount`. home.html ships a working static
   bar in that same slot, so the static one is DETACHED first — two role=banner
   landmarks is a worse bug than a plain bar, and js/ui/account.js finds its
   slot with $('.topbar-group.last'), which would otherwise resolve to the
   static bar's group. If createNav throws, the static bar goes straight back;
   it is a working bar, not a placeholder. */
let nav = null;
async function mountNav() {
  const r = await resolve('./nav.js', ['createNav']);
  const fallback = $('#home-nav-fallback');
  if (!r.fn || !shell) return;
  const anchor = fallback && fallback.parentNode === shell ? fallback : null;
  try {
    if (anchor) anchor.remove();
    const out = r.fn({ page: 'home', mount: shell, version: 'v1' });
    if (!nodeOf(out)) throw new Error('createNav returned no element');
    nav = out;
  } catch {
    nav = null;
    if (anchor && shell.firstChild) shell.insertBefore(anchor, shell.firstChild);
    else if (anchor) shell.appendChild(anchor);
  }
}

/** Tell the bar's LED what the one call found. nav.js's setBackend reads
 *  api.js's live `api` object and explicitly does NOT set data-mode — that is
 *  this module's job, and two writers is the §2.1 failure. */
function paintBackend() {
  if (nav && typeof nav.setBackend === 'function') {
    try { nav.setBackend(api.api); } catch { /* the bar is not load-bearing */ }
  }
}

let grid = null;
let recentFn = null;
let recentWhy = '';
let dialogs = null;
let lastProjects = [];

async function main() {
  // Fire the one call FIRST, so resolving four local modules never sits in
  // front of the network on a cold instance.
  const wsPromise = api.getWorkspace().then(
    (ws) => ({ ok: true, ws }),
    (err) => ({ ok: false, err }),
  );

  const stopLoading = startLoading();

  // mountNav() is awaited alongside the rest rather than fired and forgotten:
  // render() calls paintBackend(), and a nav that is still resolving at that
  // moment leaves the bar's LED reading "probing…" over a call that has
  // already come back. The network is not delayed by this — wsPromise is
  // already in flight above, and these four are local module loads.
  const [, cardMod, recentMod, createMod, res] = await Promise.all([
    mountNav(),
    resolve('./projectCard.js', ['createProjectCard']),
    resolve('./recentRuns.js', ['createRecentRuns']),
    resolve('./createProject.js', ['createProjectDialogs']),
    wsPromise,
  ]);

  stopLoading();
  wire(cardMod, recentMod, createMod);
  render(res);
}

/** Re-run the one call and repaint. The single place anything that mutates a
 *  project comes back to, which is what createProjectDialogs' `onChanged` is
 *  for. */
async function refresh() {
  try { render({ ok: true, ws: await api.getWorkspace() }); }
  catch (err) { render({ ok: false, err }); }
}

function wire(cardMod, recentMod, createMod) {
  recentFn = recentMod.fn;
  recentWhy = recentMod.why;

  /* createProjectDialogs(defaults) -> {create, rename, edit, duplicate, remove}
     `projects` may be a function, which is what keeps the name-collision check
     honest after a refresh without re-binding. `sections` is deliberately not
     passed: openCreateProject fetches GET /api/sections itself only when it
     opens, so the home page's boot stays one call. */
  if (createMod.fn) {
    try {
      dialogs = createMod.fn({ projects: () => lastProjects, onChanged: refresh });
    } catch (err) {
      dialogs = null;
      recentWhy = recentWhy || '';
      createMod.why = 'js/home/createProject.js threw while binding: ' + errText(err);
    }
  }

  const act = (name) => (dialogs && typeof dialogs[name] === 'function'
    ? (project) => dialogs[name](project)
    : undefined);

  /* createProjectCard(project, handlers) -> <article>, never null.
     handlers: onOpen(project, ev) / onRename / onDuplicate / onDelete. The
     card's own delete confirmation states that runs are NOT deleted before
     onDelete is ever called (contract §2). */
  const renderCard = cardMod.fn
    ? (p) => nodeOf(cardMod.fn(p, {
      onOpen: openProject,
      onRename: act('rename'),
      onDuplicate: act('duplicate'),
      onDelete: act('remove'),
    }))
    : null;

  grid = createProjectGrid({
    renderCard,
    onOpen: openProject,
    onCreate: dialogs && typeof dialogs.create === 'function'
      ? () => dialogs.create()
      : null,
    onRetry: () => location.reload(),
    cardMiss: cardMod.fn ? '' : (cardMod.why || ''),
  });
}

/* createRecentRuns(runs, handlers) -> <section class="rr">, never null. It
   brings its own <h2> and aria-labelledby, so this page must NOT wrap it in a
   second heading; only the placeholder path needs one. */
function recentSection(runs) {
  if (recentFn) {
    try {
      const el = nodeOf(recentFn(runs, {}));
      if (el) return el;                       // guarded: never appendChild(null)
      recentWhy = recentWhy || 'js/home/recentRuns.js returned nothing.';
    } catch (err) {
      recentWhy = 'js/home/recentRuns.js threw: ' + errText(err);
    }
  }
  return h('section#recent.home-sect', { 'aria-labelledby': 'recent-h' },
    h('h2#recent-h.home-h2', 'Recent runs, across every project'),
    missNote(
      'The recent-runs strip is not on this page.',
      recentWhy || 'js/home/recentRuns.js is not available.',
      [
        `This server reported ${int(runs.length)} recent run${runs.length === 1 ? '' : 's'} `
        + 'across all projects (the workspace call returns at most 12). They are not lost: '
        + 'every one of them is in the Studio’s run history, and each project above '
        + 'opens into the Studio with its own runs attached.',
      ]));
}

/* ====================================================== the four states == */

/* LOADING AND EMPTY MUST NEVER BE MISTAKEN FOR EACH OTHER, and the guarantee is
 * structural rather than cosmetic: the grid is not constructed until wire()
 * runs, and wire() runs only after the workspace promise has SETTLED, so
 * #home-body is genuinely empty for the whole wait and the words "no projects"
 * cannot be on screen while the call is still in flight. A 22 s wake that
 * showed "no projects yet" and then a populated grid would have told a new
 * reader something false, and no amount of styling fixes that ordering.
 * On top of that they share no class: the wait is a bordered .boot card with a
 * live elapsed count, an indeterminate bar and a pulsing LED; the empty state
 * is a .home-empty panel with a numbered three-step start. */
function render(res) {
  clear(bodyHost);
  if (!grid) return;

  /* --- UNREACHABLE: nothing answered ---------------------------------- */
  if (!res.ok) {
    const err = res.err;
    if (err && err.status === 401) { showGate(render); return; }

    // ORDER MATTERS. The stub latch is checked BEFORE the status, because once
    // api.js has fallen back, the status on the error is js/mock.js's router
    // refusing — not a server's answer. See stubRefusedBanner.
    const stubbed = !!(api.api && api.api.stub);
    if (!stubbed && api.api && err && typeof err.status === 'number' && err.status > 0) {
      // A real status means a real server answered, which is exactly what
      // api.js's own rawFetch records on any completed response.
      api.api.reachable = true;
    }
    // The Studio's rule verbatim (js/main.js): stub latched -> 'stub'. The
    // session IS in stub mode even though this particular call produced no
    // fabricated data to show, and one rule for this attribute beats a second
    // interpretation of it.
    document.documentElement.setAttribute('data-mode', stubbed ? 'stub' : 'live');
    paintBackend();
    fill(statusHost, stubbed ? stubRefusedBanner(err) : unreachableBanner(err));
    lastProjects = [];
    grid.render({ projects: [], unreachable: true });
    bodyHost.appendChild(grid.node);
    return;
  }

  const ws = (res.ws && typeof res.ws === 'object') ? res.ws : {};
  const server = (ws.server && typeof ws.server === 'object') ? ws.server : {};

  const clientStub = !!(api.api && api.api.stub);
  const serverStub = server.stub === true;

  // Labelled EXACTLY as the Studio labels it (contract §6): js/main.js sets
  // this same attribute, and css/components.css hangs the hatched magenta
  // frame and the "SYNTHETIC · NO BACKEND · js/mock.js · NOT DATA" band off
  // html[data-mode="stub"]. One attribute, not a second divergent treatment.
  document.documentElement.setAttribute('data-mode',
    (clientStub || serverStub) ? 'stub' : 'live');
  paintBackend();

  const banners = [];
  if (clientStub) banners.push(stubBanner());
  else if (serverStub) banners.push(serverStubBanner());
  if (server.contract && server.contract !== 'v1') banners.push(contractBanner(server.contract));

  if (banners.length) fill(statusHost, banners);
  else clear(statusHost);

  const projects = Array.isArray(ws.projects) ? ws.projects : [];
  const runs = Array.isArray(ws.recent_runs) ? ws.recent_runs : [];
  lastProjects = projects;

  /* --- EMPTY / POPULATED ---------------------------------------------- */
  grid.render({
    projects,
    recent_runs: runs,
    sections_available: ws.sections_available,
    drugs_available: ws.drugs_available,
    totals: ws.totals,
    stub: clientStub || serverStub,
    unreachable: false,
  });
  bodyHost.appendChild(grid.node);

  // recent_runs is "across ALL projects", so with zero projects there is
  // nothing for the strip to hold and an empty strip would only dilute the one
  // screen on this page that has to teach the product.
  if (projects.length) bodyHost.appendChild(recentSection(runs));
}

main();
