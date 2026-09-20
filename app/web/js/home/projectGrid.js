// projectGrid.js — the project grid, its search, its sort, and the empty
// state (docs/WORKSPACE_CONTRACT.md §5).
//
// THE SORT IS THE LOAD-BEARING PART OF THIS FILE.
// -----------------------------------------------
// The contract's §3 makes `last_opened_at` and `updated_at` two different
// facts on purpose: `updated_at` moves when the project DOCUMENT changes (a
// specimen added, a run attached, notes edited), `last_opened_at` moves when a
// HUMAN looks at it. Sorting by the first puts a project at the top because a
// background job touched it; sorting by the second puts it at the top because
// you were working on it.
//
// The default here is "last opened, falling back to last modified", which
// means the list is ordered by a mixture of two different facts. A control
// labelled only "recent" would be the control lying about which one it used,
// so the caption under it states, with counts and denominators, how many
// projects are placed by which stamp. `null` is rendered as "not opened yet" —
// never as an old date, never as "now".
//
// WHAT THIS FILE WILL NOT DO
// --------------------------
// It does not order anything by a predicted effect. Cross-compound ranking was
// withdrawn 2026-08-21 and it is not smuggled back in as a sort option, a
// "best" badge or a default. Every sort offered below is over a timestamp, a
// name, or a count of stored objects. No kill fraction appears on this page.
//
// Owns: the <section id="projects"> element and everything inside it. The
// individual card is js/home/projectCard.js's; when that module is missing
// this file draws a REDUCED card that says so rather than leaving a hole.

import { h, fill, clear } from '../util/dom.js';
import { int, ts } from '../util/fmt.js';

const SORTS = [
  { id: 'recency', label: 'Last opened, then last modified' },
  { id: 'updated', label: 'Last modified' },
  { id: 'created', label: 'Created, newest first' },
  { id: 'name',    label: 'Name, A to Z' },
  { id: 'runs',    label: 'Attached runs, most first' },
];
const DEFAULT_SORT = 'recency';

/** A finite number, or null. The backend sends `null` for "we do not know",
 *  and this module never converts that into 0 or into a date. */
function tnum(v) {
  return (typeof v === 'number' && Number.isFinite(v)) ? v : null;
}

function cnum(v) {
  const n = tnum(v);
  return n === null ? null : Math.max(0, Math.round(n));
}

/** Which stamp this project is actually placed by under the default sort. */
function basisOf(p) {
  const opened = tnum(p && p.last_opened_at);
  if (opened !== null) return { key: 'last_opened_at', label: 'last opened', t: opened };
  return { key: 'updated_at', label: 'last modified', t: tnum(p && p.updated_at) };
}

/** Descending, with "we do not know" always last — never first, and never
 *  silently treated as 0, which would sort a missing stamp as 1970. */
function descNulls(a, b) {
  if (a === null && b === null) return 0;
  if (a === null) return 1;
  if (b === null) return -1;
  return b - a;
}

function nameKey(p) {
  const n = p && typeof p.name === 'string' ? p.name.trim() : '';
  return n ? n.toLowerCase() : null;
}

function haystack(p) {
  const bits = [p.id, p.name, p.notes];
  for (const k of ['sections', 'drugs', 'compounds']) {
    if (Array.isArray(p[k])) bits.push(p[k].join(' '));
  }
  return bits.filter((x) => typeof x === 'string' && x).join(' \n ').toLowerCase();
}

/** Every whitespace-separated token must appear somewhere. */
function matches(p, tokens) {
  if (!tokens.length) return true;
  const hay = haystack(p);
  for (const t of tokens) if (!hay.includes(t)) return false;
  return true;
}

function sectionSet(p) {
  const raw = Array.isArray(p && p.sections) ? p.sections : [];
  return Array.from(new Set(raw.filter((s) => typeof s === 'string' && s))).sort();
}

function studioHref(id) {
  return '/studio?project=' + encodeURIComponent(String(id));
}

/**
 * @param {object} [opts]
 * @param {(project: object, ctx: object) => (Node|null)} [opts.renderCard]
 *        js/home/projectCard.js's card factory, injected by main.js. It may be
 *        absent, it may throw, and it may return null — all three are handled
 *        at the call site, because appendChild(null) is a TypeError that kills
 *        the rest of the render.
 * @param {(project: object) => void} [opts.onOpen]  stamp recency, then navigate.
 * @param {() => void} [opts.onCreate]  open js/home/createProject.js's dialog.
 * @param {() => void} [opts.onRetry]   re-attempt the workspace call.
 * @param {string} [opts.cardMiss]      why renderCard is missing, if it is.
 */
export function createProjectGrid(opts = {}) {
  const renderCard = typeof opts.renderCard === 'function' ? opts.renderCard : null;
  const onOpen = typeof opts.onOpen === 'function' ? opts.onOpen : null;
  const onCreate = typeof opts.onCreate === 'function' ? opts.onCreate : null;
  const onRetry = typeof opts.onRetry === 'function' ? opts.onRetry : null;
  const cardMiss = typeof opts.cardMiss === 'string' ? opts.cardMiss : '';

  let projects = [];
  let facts = {};
  let stub = false;
  let unreachable = false;
  let query = '';
  let sortId = DEFAULT_SORT;

  /* --- controls, built ONCE --------------------------------------------
   * Rebuilding the search input on every keystroke takes the caret with it.
   * Only the list region below is re-rendered. */

  const search = h('input#pgrid-q.inp', {
    type: 'search',
    autocomplete: 'off',
    spellcheck: 'false',
    placeholder: 'name, notes, project id, section id',
    oninput: () => { query = search.value; apply(); },
  });

  const sortSel = h('select#pgrid-sort.inp', {
    onchange: () => { sortId = sortSel.value; apply(); },
  }, SORTS.map((s) => h('option', { value: s.id }, s.label)));
  sortSel.value = DEFAULT_SORT;

  const basisNote = h('p.pgrid-basis');
  const countNote = h('p.pgrid-count', { 'aria-live': 'polite' });

  // CREATE LIVES IN THE CONTROLS BAR, NOT ONLY IN THE EMPTY STATE.
  //
  // It was rendered solely inside `div.home-empty`, so the button existed
  // only while a user had ZERO projects. Making a second project from the
  // home page was impossible -- confirmed in a browser by the agent that
  // built the dialogs. The empty state kept its own prominent call to
  // action, which is right for a first-time user; this is the one a
  // returning user needs, and the two are not alternatives.
  const headCreate = onCreate
    ? h('div.pgrid-field.pgrid-new',
        h('button#pgrid-create.btn.primary', {
          type: 'button',
          title: 'Create a project. Nothing is simulated and no run is '
            + 'started.',
          onclick: () => onCreate(),
        }, 'New project'))
    : null;

  const head = h('div.pgrid-head',
    h('div.pgrid-field.grow',
      h('label.lbl', { for: 'pgrid-q' }, 'Search projects'),
      search),
    h('div.pgrid-field',
      h('label.lbl', { for: 'pgrid-sort' }, 'Order by'),
      sortSel),
    headCreate);

  const list = h('ul.pgrid');
  const listRegion = h('div.pgrid-region');

  const heading = h('h2#pgrid-h.home-h2', 'Projects');
  const body = h('div.pgrid-body');
  const node = h('section#projects.home-sect', {
    tabindex: '-1',
    'aria-labelledby': 'pgrid-h',
  }, heading, body);

  /* --- the sort --------------------------------------------------------- */

  function sorted(items) {
    const arr = items.slice();
    const tie = (a, b) => {
      const an = nameKey(a); const bn = nameKey(b);
      if (an !== bn) {
        if (an === null) return 1;
        if (bn === null) return -1;
        return an < bn ? -1 : 1;
      }
      return String(a.id || '') < String(b.id || '') ? -1 : 1;
    };
    let cmp;
    if (sortId === 'updated') {
      cmp = (a, b) => descNulls(tnum(a.updated_at), tnum(b.updated_at));
    } else if (sortId === 'created') {
      cmp = (a, b) => descNulls(tnum(a.created_at), tnum(b.created_at));
    } else if (sortId === 'name') {
      cmp = () => 0;
    } else if (sortId === 'runs') {
      cmp = (a, b) => descNulls(cnum(a.n_runs), cnum(b.n_runs));
    } else {
      cmp = (a, b) => descNulls(basisOf(a).t, basisOf(b).t);
    }
    arr.sort((a, b) => cmp(a, b) || tie(a, b));
    return arr;
  }

  /* --- the caption that says WHICH stamp -------------------------------- */

  function updateBasis() {
    const n = projects.length;
    const opened = projects.filter((p) => tnum(p.last_opened_at) !== null).length;
    const out = [];

    if (sortId === 'recency') {
      out.push(h('b', 'Ordered by last opened'),
        ' — the last time a person opened the project through this workspace. ');
      if (opened === n) {
        out.push(`All ${int(n)} of ${int(n)} projects carry an opened time.`);
      } else if (opened === 0) {
        out.push(`0 of ${int(n)} projects has been opened through this workspace yet, so `
          + `all ${int(n)} are placed by `, h('b', 'last modified'), ' instead — a stamp '
          + 'that moves when a run is attached or a note is edited, not when a person looks '
          + 'at the project.');
      } else {
        out.push(`${int(opened)} of ${int(n)} projects carry an opened time; the other `
          + `${int(n - opened)} are placed by `, h('b', 'last modified'), ' instead — a '
          + 'stamp that moves when a run is attached or a note is edited, not when a person '
          + 'looks at the project.');
      }
    } else if (sortId === 'updated') {
      const have = projects.filter((p) => tnum(p.updated_at) !== null).length;
      out.push(h('b', 'Ordered by last modified'),
        ' — when the project document last changed: a specimen added, a run attached, '
        + 'a note edited. This is not when anyone looked at it. ');
      out.push(have === n
        ? `All ${int(n)} of ${int(n)} projects carry a modified time.`
        : `${int(have)} of ${int(n)} projects carry a modified time; the other `
          + `${int(n - have)} have none and are placed last.`);
    } else if (sortId === 'created') {
      const have = projects.filter((p) => tnum(p.created_at) !== null).length;
      out.push(h('b', 'Ordered by creation time'), ', newest first. ');
      out.push(have === n
        ? `All ${int(n)} of ${int(n)} projects carry a creation time.`
        : `${int(have)} of ${int(n)} projects carry a creation time; the other `
          + `${int(n - have)} have none and are placed last.`);
    } else if (sortId === 'name') {
      const unnamed = projects.filter((p) => nameKey(p) === null).length;
      out.push(h('b', 'Ordered by name'), ', A to Z, case-insensitive. ');
      out.push(unnamed === 0
        ? `All ${int(n)} of ${int(n)} projects are named.`
        : `${int(unnamed)} of ${int(n)} projects have no name and are placed last.`);
    } else if (sortId === 'runs') {
      out.push(h('b', 'Ordered by attached runs'), ', most first. The count is '
        + 'run REFERENCES held by a project’s specimens: a run that two specimens both '
        + 'point at is counted twice, and the runs themselves live in this server’s run '
        + `store, not in the project. All ${int(n)} of ${int(n)} projects are ordered this way.`);
    }

    fill(basisNote, out);
  }

  /* --- one card, with the appendChild(null) guard at the call site ------- */

  function cardFor(p) {
    const basis = basisOf(p);
    const ctx = { basis: basis.key, basisLabel: basis.label, basisTime: basis.t, stub, onOpen };
    let el = null;
    let why = cardMiss;

    if (renderCard) {
      try {
        el = renderCard(p, ctx);
      } catch (err) {
        el = null;
        why = 'js/home/projectCard.js threw while drawing this card: '
          + ((err && err.message) ? err.message : String(err));
      }
      if (el && !el.nodeType) {
        why = 'js/home/projectCard.js returned a ' + (typeof el) + ', not a DOM node.';
        el = null;
      } else if (!el && !why) {
        why = 'js/home/projectCard.js returned nothing for this project.';
      }
    } else if (!why) {
      why = 'js/home/projectCard.js is not available.';
    }

    // NEVER hand a null to appendChild: it is a TypeError, and it takes the
    // rest of the grid down with it.
    return el || reducedCard(p, basis, why);
  }

  function reducedCard(p, basis, why) {
    const specimens = cnum(p.n_specimens);
    const runs = cnum(p.n_runs);
    const secs = sectionSet(p);
    const name = (typeof p.name === 'string' && p.name.trim()) ? p.name.trim() : null;

    const counts = [];
    counts.push(h('span', specimens === null
      ? 'specimen count not reported'
      : `${int(specimens)} specimen${specimens === 1 ? '' : 's'}`));
    counts.push(h('span', runs === null
      ? 'run count not reported'
      : `${int(runs)} run reference${runs === 1 ? '' : 's'}`
        + (specimens === null
          ? ''
          : ` across ${int(specimens)} specimen${specimens === 1 ? '' : 's'}`)));

    const when = basis.t === null
      ? h('div.pc-when', h('span.pc-when-k', 'Last opened: '),
          h('span.pc-never', 'not opened yet, and no modified time either'))
      : h('div.pc-when',
          h('span.pc-when-k', basis.key === 'last_opened_at'
            ? 'Last opened: ' : 'Not opened yet · last modified: '),
          h('span.mono', ts(basis.t)));

    return h('article.pcard-min',
      h('h3', name || h('span.pc-unnamed', 'Unnamed project')),
      h('div.pc-id', String(p.id || 'no project id')),
      h('div.pc-facts', counts),
      when,
      h('div.pc-sections', secs.length
        ? ['Sections, as a set: ', h('span.mono', secs.join(', '))]
        : 'No section is attached to this project yet.'),
      h('div.pc-act',
        h('a.btn.primary', {
          href: studioHref(p.id),
          onclick: (ev) => {
            if (ev.metaKey || ev.ctrlKey || ev.shiftKey || ev.altKey || ev.button !== 0) return;
            if (!onOpen) return;
            ev.preventDefault();
            onOpen(p);
          },
        }, 'Open in the Studio')),
      why ? h('div.mount-miss',
        h('b', 'Reduced card. '),
        'The full project card, including the measured-composition bar, is drawn by '
        + 'js/home/projectCard.js. It is not drawing this one, so what you see above is '
        + 'projectGrid.js’s fallback: the stored facts, and nothing derived from them.',
        h('span.mm-why', why)) : null);
  }

  /* --- the states this section can be in -------------------------------- */

  function noMatchBlock() {
    const n = projects.length;
    return h('div.pgrid-none',
      h('p', 'No project matches ', h('span.pn-q', '“' + query.trim() + '”'), '.'),
      h('p', `All ${int(n)} of ${int(n)} projects on this server are filtered out by it. `
        + 'The search reads the project name, its notes, its project id and the section ids '
        + 'attached to it — every word you type has to appear in one of those.'),
      h('div.he-act',
        h('button.btn.primary', {
          type: 'button',
          onclick: () => { query = ''; search.value = ''; apply(); search.focus(); },
        }, 'Clear the search')));
  }

  /** Reachable server, zero projects. The screen a new reader actually starts
   *  on, so it has to teach the product and end on a first run. */
  function emptyBlock() {
    const sections = cnum(facts.sections_available);
    const drugs = cnum(facts.drugs_available);

    const secLine = sections === null
      ? 'This server did not say how many tissue sections it has registered.'
      : (sections === 0
        ? '0 tissue sections are registered on this server right now, so there is nothing '
          + 'to simulate until one is.'
        : `${int(sections)} tissue section${sections === 1 ? '' : 's'} `
          + `${sections === 1 ? 'is' : 'are'} already registered on this server.`);

    const drugLine = drugs === null
      ? 'This server did not say how many compounds it carries.'
      : (drugs === 0
        ? '0 compounds are loaded on this server, so a run cannot be submitted yet.'
        : `${int(drugs)} compound${drugs === 1 ? '' : 's'} `
          + `${drugs === 1 ? 'is' : 'are'} loaded on this server.`);

    const createAct = onCreate
      ? h('div.he-act',
          h('button.btn.primary', { type: 'button', onclick: () => onCreate() }, 'New project'))
      : h('div.mount-miss',
          h('b', 'There is no create-project button on this page. '),
          'js/home/createProject.js did not load, so this step has no control. The Studio '
          + 'carries the same project store behind its ',
          h('span.mono', 'projects'),
          ' launcher, and a project made there appears here on reload.',
          cardMiss ? h('span.mm-why', cardMiss) : null);

    return h('div.home-empty',
      h('h2', 'No projects on this server yet.'),

      h('p.he-lede',
        h('b', 'Vivum simulates what a compound does inside a measured tumour tissue section'),
        ' — where it reaches, where it does not, and which cells that accounts for. '
        + 'The tissue is a real, de-identified section; the transport, the dose and the '
        + 'per-cell outcome are computed from it.'),

      h('p.he-lede',
        'A ', h('b', 'project'), ' is where that work is kept. It owns nothing but ids: it '
        + 'holds specimens, and each specimen points at runs that live in this server’s '
        + 'run store. That is why deleting a project never deletes its runs, and why two '
        + 'specimens can point at the same run.'),

      h('div.he-limit',
        'Before the first run: everything this product produces is model output. No run '
        + 'here has been scored against a patient outcome, and nothing in it is for '
        + 'clinical use.'),

      h('ol.he-steps',
        h('li.he-step',
          h('span.he-num', '1'),
          h('div.he-body',
            h('h3.home-h3', 'Make a project'),
            h('p', 'A name, and notes if you want them. Nothing is simulated yet — this is '
              + 'just the folder the runs get kept in, so that in a week you can still tell '
              + 'which specimen a number came off.'),
            createAct)),

        h('li.he-step' + (sections === 0 ? '.blocked' : ''),
          h('span.he-num', '2'),
          h('div.he-body',
            h('h3.home-h3', 'Point the server at a tissue section'),
            h('p', 'Tissue is registered by ', h('b', 'filesystem path'),
              ' on the machine running this server. This page never moves image data and '
              + 'never holds any. ', secLine),
            h('div.he-act',
              h('a.btn', { href: '/ingest' }, 'Register a section by path')))),

        h('li.he-step',
          h('span.he-num', '3'),
          h('div.he-body',
            h('h3.home-h3', 'Run a therapy, and read the mechanism'),
            h('p', 'The Studio is where a section, a compound and a RunSpec become a run: '
              + 'the cell viewer, the metrics, the per-cell attribution of why each cell '
              + 'died or survived, and the caveats that come with all of it. ', drugLine),
            h('div.he-act',
              h('a.btn.primary', { href: '/studio' }, 'Open the Studio'))))),

      h('p.he-foot',
        'You do not need a project to run the simulator. ',
        h('a', { href: '/studio' }, 'The Studio'),
        ' works on its own exactly as it always has — a project is what makes a run '
        + 'findable, comparable and still explicable a week later.'));
  }

  /** The server could not be reached. This is NOT "you have no projects": it
   *  is "nothing was answered", and the two must never look alike. */
  function unreachableBlock() {
    return h('div.home-empty',
      h('h2', 'No project list, because this page could not reach the server.'),
      h('p.he-lede',
        'This is not “there are no projects”. It is “nothing answered”. '
        + 'A project list may well exist on this server; this page has no way to know, and '
        + 'it will not guess one. The banner above this section says what failed.'),
      h('p.he-lede',
        'A sleeping instance is the usual cause, and it clears itself: the first request '
        + 'wakes the container and the one after it is fast.'),
      h('div.he-act',
        onRetry
          ? h('button.btn.primary', { type: 'button', onclick: () => onRetry() },
              'Try the workspace call again')
          : null,
        h('a.btn', { href: '/studio' }, 'Open the Studio')));
  }

  /* --- render ------------------------------------------------------------ */

  function apply() {
    if (unreachable) {
      heading.textContent = 'Projects';
      fill(body, unreachableBlock());
      return;
    }
    if (!projects.length) {
      heading.textContent = 'Getting started';
      fill(body, emptyBlock());
      return;
    }

    heading.textContent = 'Projects';
    if (body.firstChild !== head) fill(body, head, basisNote, countNote, listRegion);

    updateBasis();

    const tokens = query.trim().toLowerCase().split(/\s+/).filter(Boolean);
    const shown = sorted(projects.filter((p) => matches(p, tokens)));
    const n = projects.length;

    const tail = [];
    const specimens = cnum(facts.totals && facts.totals.specimens);
    const runs = cnum(facts.totals && facts.totals.runs);
    if (specimens !== null) tail.push(`${int(specimens)} specimen${specimens === 1 ? '' : 's'}`);
    if (runs !== null) tail.push(`${int(runs)} run${runs === 1 ? '' : 's'}`);

    fill(countNote, tokens.length
      ? [h('b', `${int(shown.length)} of ${int(n)}`),
         ` project${n === 1 ? '' : 's'} match this search.`]
      : [h('b', `${int(n)}`),
         ` project${n === 1 ? '' : 's'} on this server`
         + (tail.length ? `, holding ${tail.join(' and ')}.` : '.')]);

    if (!shown.length) { fill(listRegion, noMatchBlock()); return; }

    clear(list);
    for (const p of shown) list.appendChild(h('li', cardFor(p)));
    fill(listRegion, list);
  }

  /**
   * @param {object} ws  the /api/workspace payload, plus the two flags main.js
   *                     resolves from api.js: {stub, unreachable}.
   */
  function render(ws) {
    const w = ws || {};
    projects = Array.isArray(w.projects)
      ? w.projects.filter((p) => p && typeof p === 'object')
      : [];
    facts = w;
    stub = w.stub === true;
    unreachable = w.unreachable === true;
    apply();
  }

  return {
    node,
    render,
    focusSearch: () => { if (body.contains(search)) search.focus(); },
  };
}
