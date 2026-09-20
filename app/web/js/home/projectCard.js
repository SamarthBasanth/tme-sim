// projectCard.js — one project, as a card.
//
// WHAT THE CARD IS. A project is the workspace's unit of kept work, and this
// card is the only place most people will ever read one. It therefore has to
// answer, without being opened: what is this called, what did I write about
// it, how much is in it, what tissue and which compounds it touches, when I
// last looked at it, and what it measured. Everything on it is a fact the
// server reported; nothing on it is inferred, ordered or scored here.
//
// FOUR THINGS THIS CARD IS NOT ALLOWED TO DO, and where each is enforced:
//
//   1. It never ORDERS COMPOUNDS. Cross-compound ranking was withdrawn
//      2026-08-21 (results/pk_ranking_robustness.md,
//      results/capacity_ranking_robustness.md). The compounds a project has
//      run are drawn as a SET — deduplicated, alphabetical, no metric, no
//      "best", no first-among-equals position — and the card says the word
//      "alphabetical" on screen so the order cannot be misread as a result.
//      An alphabetical list read as a ranking is the exact failure the
//      withdrawal is about, so the label is not optional decoration.
//   2. It never puts an ABSOLUTE KILL FRACTION on a card. There is no metric
//      on this card at all. The only numbers are counts of things the project
//      contains, and a composition chart's fractions of a stated cell count.
//   3. The COMPOSITION THUMBNAIL is a chart, never a picture of tissue
//      (WORKSPACE_CONTRACT.md §4). It is drawn as a labelled stacked bar with
//      a written key, at full card width and 12px tall — a shape no one reads
//      as a micrograph — and it carries a sentence saying what it is.
//   4. DELETE states, at the point of deletion, that the runs are not deleted.
//      A project owns nothing but ids (§2), so deleting one detaches runs; it
//      never removes them from app/runs/. The card refuses to let that be
//      discovered afterwards: the second step of the delete is a paragraph,
//      not a tooltip, and it is rendered by THIS module rather than delegated,
//      so no dialog owner can drop it.
//
// RECENCY, §3. `last_opened_at` is "a human looked at it"; `updated_at` is
// "the document changed". They are not interchangeable and this card never
// prints one under the other's name. When `last_opened_at` is null the card
// says "not opened yet" — never an old date, never "now" — and prints the
// document's own update time on a separate line, labelled as what it is. The
// fallback the contract asks for is a SORT fallback, and it is exported as
// `projectRecency().sortAt` for the grid rather than smuggled into the label.
//
// KEYBOARD. Every control is a real <button> or a real <a href>. The card is
// not a clickable <div>: "open" is a link to the Studio that works with no
// JavaScript wired at all, which is also why a middle-click and a ctrl-click
// behave the way a link should. Focus rings come from base.css's global
// :focus-visible rule; nothing here removes an outline.

import { h, clear } from '../util/dom.js';
import { int, ts, pct } from '../util/fmt.js';

/* ------------------------------------------------------------ constants */

/** Notes longer than this are truncated in the card, with the full text one
 *  keypress away. Truncation is only honest if the reader can see that it
 *  happened AND undo it; both are on screen, not in a title. */
const NOTES_CLAMP = 220;

/** Mirrors gl/viewer.js's CELL_TYPES / CELL_TYPE_VARS. Copied, not imported:
 *  gl/viewer.js is the WebGL renderer and the home page must not pull a
 *  renderer in to draw five swatches. ui/inspect.js, ui/tools.js and
 *  ui/overlays.js each keep the same local copy for the same reason.
 *  Verified against app/web/js/gl/viewer.js:28-29. */
const TYPE_VAR = {
  TUMOR: '--ch-tumor',
  STROMAL: '--ch-stromal',
  IMMUNE: '--ch-immune',
  ENDOTHELIAL: '--ch-endothelial',
  NECROTIC: '--ch-necrotic',
};

/** A label the renderer has no channel for still gets a row in the key with
 *  its own name and its own percentage — colour is never the only channel
 *  (UI_ARCHITECTURE.md §4.5) — but it does NOT borrow a cell type's colour,
 *  because drawing an unrecognised label in tumour red asserts an identity
 *  the server never sent. */
const UNKNOWN_TYPE_VAR = '--fg-3';

/* ----------------------------------------------------------------- time */

/** Epoch seconds (what projects.py writes with time.time()) or epoch
 *  milliseconds -> ms. Same >1e11 discrimination util/fmt.js's ts() uses, so
 *  the card and the run table never disagree about what a timestamp means. */
export function epochMs(t) {
  const n = Number(t);
  if (t === null || t === undefined || t === '' || !Number.isFinite(n) || n <= 0) return null;
  return n > 1e11 ? n : n * 1000;
}

const REL_UNITS = [
  [31557600, 'year'], [2629800, 'month'], [86400, 'day'],
  [3600, 'hour'], [60, 'minute'],
];

/**
 * "3 days ago". Null for a timestamp that is not a timestamp.
 *
 * A clock-skewed future stamp is reported as such rather than rounded down to
 * "just now": a home page that says a project was opened a moment ago when the
 * server's clock is a day ahead is lying about the one field §3 exists for.
 */
export function relTime(t, now = Date.now()) {
  const ms = epochMs(t);
  if (ms === null) return null;
  const s = Math.round((now - ms) / 1000);
  if (s < -90) return 'at a timestamp in the future';
  if (s < 60) return 'less than a minute ago';
  for (const [size, name] of REL_UNITS) {
    if (s >= size) {
      const n = Math.floor(s / size);
      return n + ' ' + name + (n === 1 ? '' : 's') + ' ago';
    }
  }
  return 'less than a minute ago';
}

/**
 * The §3 recency rule, in one place, because the card and the grid's sort
 * must not implement it twice and differently.
 *
 * @returns {{kind: 'opened'|'never'|'unknown', at: number|null,
 *            sortAt: number|null, label: string, abs: string|null,
 *            secondary: string|null, title: string}}
 *   `kind`     — what we actually know.
 *   `sortAt`   — the contract's sort key: last_opened_at, falling back to
 *                updated_at, then created_at, and NULL when there is no
 *                timestamp at all. Null rather than 0, because 0 sorts a
 *                project with no stamps as 1970 and puts it at the bottom
 *                for the wrong reason; projectGrid.js's descNulls() already
 *                treats null as "we do not know" and sorts it last on
 *                purpose. A SORT key only; it is never the thing rendered
 *                as "opened".
 *   `label`    — safe to print verbatim. Null last_opened_at renders as
 *                "not opened yet", never as an old date and never as "now".
 *   `secondary`— the document's own update time, labelled as a document
 *                change, or null when there is nothing to say.
 */
export function projectRecency(project, now = Date.now()) {
  const p = project || {};
  const opened = epochMs(p.last_opened_at);
  const updated = epochMs(p.updated_at);
  const created = epochMs(p.created_at);
  const updLine = updated === null ? null
    : 'document updated ' + relTime(updated, now) + ' (' + ts(updated / 1000) + ')';

  if (opened !== null) {
    return {
      kind: 'opened',
      at: opened,
      sortAt: opened,
      label: 'opened ' + relTime(opened, now),
      abs: ts(opened / 1000),
      secondary: updLine,
      title: 'last_opened_at — when a human last opened this project. '
        + 'Distinct from updated_at, which moves whenever the project '
        + 'document changes, including when a background job attaches a run.',
    };
  }
  return {
    kind: updated === null && created === null ? 'unknown' : 'never',
    at: null,
    sortAt: updated !== null ? updated : created,
    label: 'not opened yet',
    abs: null,
    secondary: updLine,
    title: 'last_opened_at is null: this project has not been opened since the '
      + 'server started recording that, so there is no open time to show. The '
      + 'document’s own update time is shown beside it, and is what this '
      + 'card is sorted on.',
  };
}

/* ------------------------------------------------------------- counting */

const nOrNull = (v) => (Number.isFinite(Number(v)) && Number(v) >= 0 ? Math.round(Number(v)) : null);

function specimenCount(p) {
  const n = nOrNull(p.n_specimens);
  if (n !== null) return n;
  return Array.isArray(p.specimens) ? p.specimens.length : null;
}

function runCount(p) {
  const n = nOrNull(p.n_runs);
  if (n !== null) return n;
  if (!Array.isArray(p.specimens)) return null;
  return p.specimens.reduce((a, s) => a + ((s && Array.isArray(s.run_ids)) ? s.run_ids.length : 0), 0);
}

const plural = (n, one, many) => int(n) + ' ' + (n === 1 ? one : many);

/* ------------------------------------------------------ delete sentence */

/**
 * The sentence the delete step must carry, with this project's own counts in
 * it. Exported so createProject.js can put the SAME words on any dialog it
 * builds — two different accounts of what delete does is worse than one.
 *
 * WORKSPACE_CONTRACT.md §2: a project owns nothing but ids.
 */
export function deleteDisclosure(project) {
  const p = project || {};
  const r = runCount(p);
  const s = specimenCount(p);
  const runs = r === null
    ? 'Its attached runs are NOT deleted.'
    : (r === 0
      ? 'This project has 0 runs attached, and deleting a project never '
        + 'deletes runs in any case.'
      : 'Its ' + plural(r, 'attached run', 'attached runs') + ' are NOT deleted.');
  const spec = s === null
    ? 'Its specimens and their run attachments go with the project document'
    : (s === 0
      ? 'It has 0 specimens; the project document is all there is to delete'
      : (s === 1
        ? 'Its 1 specimen, and that specimen\u2019s run attachments, goes with '
          + 'the project document'
        : 'Its ' + int(s) + ' specimens, and their run attachments, go with '
          + 'the project document'));
  return runs + ' A project holds run IDS, not runs: app/runs/ stays the '
    + 'single store for every run, its results and its payloads, and deleting '
    + 'a project detaches runs rather than removing them. ' + spec
    + ', which is what is deleted here. This cannot be undone.';
}

/* ------------------------------------------------- composition thumbnail */

/**
 * WORKSPACE_CONTRACT.md §4, both branches. Always returns an element, never
 * null — appendChild(null) is a TypeError that would kill the rest of the
 * card's render.
 *
 * `{kind: "composition"}` draws `type_fractions` as a stacked bar. It is a
 * CHART OF MEASURED COMPOSITION and it must never read as a rendering of
 * tissue, so: it is a 12px bar the full width of the card, every segment is
 * named in a written key with its own percentage, the cell count it is a
 * fraction OF is printed above it, and a sentence under it says in words that
 * this is a chart and not a picture. A card that merely avoided claiming to
 * be tissue would still be read as tissue; this one says what it is.
 *
 * `{kind: "none"}` renders the server's reason as visible text inside a
 * plainly-empty box — not a grey rectangle that reads as a failed image, and
 * not a hover-only explanation, which is unreachable from a keyboard.
 */
export function compositionThumb(thumb) {
  if (!thumb || typeof thumb !== 'object') {
    return placeholder('No composition chart',
      'This server did not report a thumbnail for this project. Nothing is '
      + 'wrong with the project; the field is absent from the response.');
  }
  if (thumb.kind === 'none') {
    const why = typeof thumb.reason === 'string' && thumb.reason.trim()
      ? thumb.reason.trim()
      : 'The server gave no reason. Open the project to see what it contains.';
    return placeholder('No composition yet', why);
  }
  if (thumb.kind !== 'composition') {
    return placeholder('No composition chart',
      'The server sent a thumbnail of kind "' + String(thumb.kind)
      + '", which this page does not know how to draw. It is not drawn rather '
      + 'than drawn wrongly.');
  }

  const fr = (thumb.type_fractions && typeof thumb.type_fractions === 'object')
    ? thumb.type_fractions : {};
  const rows = Object.keys(fr)
    .map((k) => ({ label: String(k), v: Number(fr[k]) }))
    .filter((r) => Number.isFinite(r.v) && r.v > 0);

  if (!rows.length) {
    return placeholder('No composition yet',
      'The thumbnail is of kind "composition" but carries no cell-type '
      + 'fractions above zero, so there is nothing to chart.');
  }

  // Keys are drawn in the renderer's channel order where they are known, so
  // the same tissue reads the same way here as it does in the Studio; unknown
  // labels follow, alphabetically. This is an order over CELL TYPES, which is
  // a display convention — it is not, and cannot become, an order over
  // compounds or over any metric.
  const known = Object.keys(TYPE_VAR);
  rows.sort((a, b) => {
    const ia = known.indexOf(a.label.toUpperCase());
    const ib = known.indexOf(b.label.toUpperCase());
    if (ia !== ib) return (ia < 0 ? 99 : ia) - (ib < 0 ? 99 : ib);
    return a.label.localeCompare(b.label);
  });

  const sum = rows.reduce((a, r) => a + r.v, 0);
  // Widths are scaled only when the fractions OVERFLOW the bar, and the
  // scaling is stated. When they fall short, the shortfall is left as a
  // visible, named remainder rather than silently stretched to fill: a bar
  // normalised to 100% would report a complete census the server did not send.
  const denom = Math.max(1, sum);
  const nCells = nOrNull(thumb.n_cells);
  const ofWhat = nCells === null
    ? 'of the reconstruction’s cells (count not reported)'
    : 'of ' + int(nCells) + ' cells';

  const seg = (label, v, varName, extraCls) =>
    h('span.pc-comp-seg' + (extraCls || ''), {
      style: { width: (v / denom * 100).toFixed(3) + '%', background: 'var(' + varName + ')' },
      title: label + ' — ' + pct(v, 1) + '% ' + ofWhat,
    });

  const keyRow = (label, v, varName, note) =>
    h('li.pc-comp-key-row',
      h('span.pc-comp-sw', { style: { background: 'var(' + varName + ')' } }),
      h('span.pc-comp-lbl', label),
      h('span.pc-comp-val.mono', pct(v, 1) + '%'),
      note ? h('span.pc-comp-lbl.dim', note) : null);

  const bar = h('div.pc-comp-bar', { 'aria-hidden': 'true' },
    rows.map((r) => seg(r.label, r.v, TYPE_VAR[r.label.toUpperCase()] || UNKNOWN_TYPE_VAR,
      TYPE_VAR[r.label.toUpperCase()] ? '' : '.pc-comp-seg-unknown')));

  const keys = rows.map((r) => keyRow(
    r.label, r.v, TYPE_VAR[r.label.toUpperCase()] || UNKNOWN_TYPE_VAR,
    TYPE_VAR[r.label.toUpperCase()] ? null : '(no channel for this label)'));

  if (sum < 0.995) {
    const rest = 1 - sum;
    bar.appendChild(seg('unlabelled', rest, '--line-1', '.pc-comp-seg-rest'));
    keys.push(keyRow('unlabelled', rest, '--line-1', null));
  }

  const warn = sum > 1.005
    ? h('p.pc-comp-warn', 'The reported fractions sum to ' + pct(sum, 1)
      + '%, not 100%. The bar is scaled to fit; the key shows each fraction as '
      + 'the server sent it.')
    : null;

  return h('figure.pc-comp',
    h('figcaption.pc-comp-head',
      h('span.pc-comp-title', 'Measured composition'),
      thumb.section ? h('span.pc-comp-sect.mono', String(thumb.section)) : null,
      h('span.pc-comp-denom', nCells === null
        ? 'cell count not reported'
        : int(nCells) + ' cells')),
    bar,
    h('ul.pc-comp-key', keys),
    warn,
    h('p.pc-comp-note',
      'Cell-type fractions ' + ofWhat + ', measured from the reconstruction. '
      + 'This is a chart of that measurement. It is not a picture of the '
      + 'tissue and nothing about its shape corresponds to where cells are.'));
}

function placeholder(title, reason) {
  return h('div.pc-thumb-none',
    h('span.pc-thumb-none-title', title),
    h('p.pc-thumb-none-reason', reason));
}

/* -------------------------------------------------------------- the set */

/**
 * The compounds a project has run, drawn as a SET.
 *
 * Alphabetical, deduplicated, and labelled "alphabetical" on screen. There is
 * no metric here, no ordering by one, and no badge on any member: ordering
 * compounds by predicted effect is withdrawn (2026-08-21) and an alphabetical
 * list that does not SAY it is alphabetical is read as a ranking by roughly
 * everyone. The label is load-bearing, not politeness.
 *
 * A missing field and an empty list are different facts and are printed
 * differently: "no compound has been run" is a statement about the project,
 * "this server did not report" is a statement about the response.
 */
function compoundSet(project) {
  const raw = project.compounds !== undefined ? project.compounds
    : (project.drugs !== undefined ? project.drugs : undefined);
  if (raw === undefined || raw === null) {
    return h('div.pc-set',
      h('span.pc-set-note', 'Compounds not reported by this server.'));
  }
  const names = Array.from(new Set(
    (Array.isArray(raw) ? raw : [raw])
      .filter((d) => d !== null && d !== undefined && String(d).trim() !== '')
      .map((d) => String(d).trim()),
  )).sort((a, b) => a.localeCompare(b));

  if (!names.length) {
    return h('div.pc-set', h('span.pc-set-note', 'No compound has been run yet.'));
  }
  return h('div.pc-set',
    h('ul.pc-set-items', names.map((d) => h('li.pc-set-item.mono', d))),
    h('p.pc-set-note', {
      title: 'Cross-compound ranking is withdrawn (2026-08-21): per-compound '
        + 'plasma exposure spans 113x across this panel and the published '
        + 'ordering breaks at a spread of 1.53x. This is the set of compounds '
        + 'this project has run — membership only.',
    }, plural(names.length, 'compound', 'compounds') + ' run, listed '
      + 'alphabetically. This is a set, not an order: nothing here ranks them.'));
}

function sectionSet(project) {
  const raw = Array.isArray(project.sections) ? project.sections : null;
  if (!raw) return null;
  const named = raw.filter((s) => s !== null && s !== undefined && String(s).trim() !== '')
    .map((s) => String(s).trim());
  const missing = raw.length - named.length;
  const names = Array.from(new Set(named)).sort((a, b) => a.localeCompare(b));
  if (!names.length && !missing) return null;
  return h('div.pc-set',
    names.length
      ? h('ul.pc-set-items', names.map((s) => h('li.pc-set-item.mono', s)))
      : null,
    h('p.pc-set-note',
      (names.length ? plural(names.length, 'section', 'sections') + ' touched'
        : 'No section recorded')
      + (missing ? ', and ' + plural(missing, 'specimen has', 'specimens have')
        + ' no section recorded' : '') + '.'));
}

/* ---------------------------------------------------------------- notes */

function notesBlock(project, domId) {
  const text = typeof project.notes === 'string' ? project.notes : '';
  if (!text.trim()) {
    return h('p.pc-notes.pc-notes-empty', { id: domId }, 'No notes.');
  }
  const body = h('span.pc-notes-text', { id: domId });
  if (text.length <= NOTES_CLAMP) {
    body.textContent = text;
    return h('p.pc-notes', body);
  }
  // Cut on a word boundary when one is near enough that the cut does not lose
  // most of a line; otherwise cut hard. Either way the ellipsis is real text
  // and the full note is one button press away — a title attribute is not a
  // way to reach anything from a keyboard.
  const space = text.lastIndexOf(' ', NOTES_CLAMP);
  const head = text.slice(0, space > NOTES_CLAMP * 0.6 ? space : NOTES_CLAMP).trimEnd();
  let open = false;
  body.textContent = head + '…';
  const more = h('button.btn.ghost.sm.pc-notes-more', {
    type: 'button',
    'aria-expanded': 'false',
    'aria-controls': domId,
  }, 'show all ' + int(text.length) + ' characters');
  more.addEventListener('click', () => {
    open = !open;
    body.textContent = open ? text : head + '…';
    more.textContent = open
      ? 'show less (' + int(text.length) + ' characters in full)'
      : 'show all ' + int(text.length) + ' characters';
    more.setAttribute('aria-expanded', open ? 'true' : 'false');
  });
  return h('p.pc-notes', body, ' ', more);
}

/* -------------------------------------------------------------- actions */

export const studioHrefForProject = (id) =>
  '/studio?project=' + encodeURIComponent(String(id));

const safeId = (s) => String(s == null ? 'unknown' : s).replace(/[^A-Za-z0-9_-]/g, '-');

/** Run a handler, keep the card honest about being busy, and never leave a
 *  control stuck disabled because a promise rejected. */
function invoke(card, btn, fn, project, ev) {
  let out;
  try {
    out = fn(project, ev);
  } catch (err) {
    console.error('projectCard: handler threw', err);
    return;
  }
  if (!out || typeof out.then !== 'function') return;
  card.setAttribute('aria-busy', 'true');
  if (btn) btn.disabled = true;
  const done = () => {
    card.removeAttribute('aria-busy');
    if (btn && btn.isConnected) btn.disabled = false;
  };
  out.then(done, (err) => { done(); console.error('projectCard: handler rejected', err); });
}

/* ----------------------------------------------------------------- card */

/**
 * @param {object} project  a ProjectSummary: {id, name, notes, created_at,
 *   updated_at, n_specimens, n_runs, sections[]} plus the contract's
 *   last_opened_at and thumbnail, plus `compounds` (alias `drugs`).
 * @param {object} [handlers]
 *   onOpen(project, ev)      optional. Absent -> the link navigates natively
 *                            to /studio?project=<id>, which is the whole
 *                            behaviour the contract asks for, so the card is
 *                            useful before anything is wired.
 *   onRename(project)        optional; opens the owner's rename dialog.
 *   onDuplicate(project)     optional.
 *   onDelete(project)        optional; called ONLY after this card has shown
 *                            the runs-are-not-deleted paragraph and the
 *                            reader has confirmed against it.
 *   Each may return a promise; the card marks itself aria-busy until it
 *   settles.
 * @returns {HTMLElement} an <article>. Never null.
 */
export function createProjectCard(project, handlers = {}) {
  const p = (project && typeof project === 'object') ? project : {};
  const hs = handlers || {};
  const pid = p.id == null ? '' : String(p.id);
  const uid = safeId(pid || Math.random().toString(36).slice(2));
  const titleId = 'pc-title-' + uid;
  const notesId = 'pc-notes-' + uid;

  const name = (typeof p.name === 'string' && p.name.trim()) ? p.name.trim() : null;

  // A project whose document could not be read comes back from projects.py's
  // list() as {id, name: null, unreadable: "<why>"}. It gets a card, because a
  // project that has silently vanished from the grid is worse than one that
  // explains itself — and no actions, because every one of them would 500.
  if (p.unreadable) {
    return h('article.pc.pc-unreadable', {
      'aria-labelledby': titleId,
      dataset: { projectId: pid },
    },
      h('header.pc-head',
        h('h3.pc-name', { id: titleId }, name || pid || 'unreadable project'),
        pid ? h('span.pc-id.mono', pid) : null),
      h('p.pc-unreadable-why',
        'This project document could not be read, so nothing about it is '
        + 'shown. The server said: ' + String(p.unreadable)
        + (pid ? ' The file is app/projects/' + pid + '/project.json on the '
          + 'server; its runs are unaffected either way, because runs live in '
          + 'app/runs/.' : '')));
  }

  const rec = projectRecency(p);
  const nSpec = specimenCount(p);
  const nRuns = runCount(p);

  const openHref = pid ? studioHrefForProject(pid) : null;
  const openIntercept = (ev) => {
    if (ev.defaultPrevented) return;
    // Leave the browser's own link gestures alone: a new tab is a link's job.
    if (ev.button !== 0 || ev.metaKey || ev.ctrlKey || ev.shiftKey || ev.altKey) return;
    if (typeof hs.onOpen !== 'function') return;   // native navigation
    ev.preventDefault();
    invoke(card, null, hs.onOpen, p, ev);
  };

  const nameEl = openHref
    ? h('a.pc-name-link', {
      href: openHref,
      onclick: openIntercept,
      title: 'Open ' + (name || pid) + ' in the Studio (' + openHref + ')',
    }, name || pid)
    : h('span.pc-name-link', name || 'unnamed project');

  const facts = h('dl.pc-facts',
    h('div.pc-fact',
      h('dt.pc-fact-k.lbl', 'specimens'),
      h('dd.pc-fact-v.mono', nSpec === null ? 'not reported'
        : plural(nSpec, 'specimen', 'specimens'))),
    h('div.pc-fact',
      h('dt.pc-fact-k.lbl', 'runs'),
      h('dd.pc-fact-v.mono', {
        title: 'Runs attached to this project’s specimens. A specimen '
          + 'holds run IDS into app/runs/; the project does not own them.',
      }, nRuns === null ? 'not reported' : plural(nRuns, 'run attached', 'runs attached'))),
    h('div.pc-fact.pc-recency',
      h('dt.pc-fact-k.lbl', 'last opened'),
      h('dd.pc-fact-v', { title: rec.title },
        h('span', { class: rec.kind === 'opened' ? null : 'pc-recency-never' }, rec.label),
        rec.abs ? h('span.pc-recency-abs.mono', rec.abs) : null,
        rec.secondary ? h('span.pc-recency-doc', rec.secondary) : null)));

  const card = h('article.pc', {
    'aria-labelledby': titleId,
    dataset: { projectId: pid },
  },
    h('header.pc-head',
      h('h3.pc-name', { id: titleId }, nameEl),
      pid ? h('span.pc-id.mono', { title: 'project id' }, pid) : null),
    h('div.pc-thumb', compositionThumb(p.thumbnail)),
    facts,
    h('div.pc-sets',
      sectionSet(p),
      compoundSet(p)),
    notesBlock(p, notesId));

  /* -- the action row --------------------------------------------------- */

  const actions = h('footer.pc-actions');
  const openBtn = openHref
    ? h('a.btn.sm.pc-action.pc-open', {
      href: openHref,
      onclick: openIntercept,
      title: 'Open in the Studio with this project’s runs attached',
    }, 'open')
    : null;
  if (openBtn) actions.appendChild(openBtn);

  const missing = [];
  const secondary = (label, key, title) => {
    const wired = typeof hs[key] === 'function';
    const btn = h('button.btn.ghost.sm.pc-action', {
      type: 'button',
      title: wired ? title : 'Unavailable: this page did not pass a ' + key + ' handler.',
      disabled: !wired,
    }, label);
    if (wired) btn.addEventListener('click', () => invoke(card, btn, hs[key], p));
    else missing.push(label);
    actions.appendChild(btn);
    return btn;
  };

  secondary('rename', 'onRename', 'Rename this project');
  secondary('duplicate', 'onDuplicate',
    'Copy this project’s specimens and their run IDS. Runs are neither '
    + 're-executed nor copied — two specimens may point at the same run.');

  /* -- delete: the disclosure is the second step, not a tooltip ---------- */

  const confirmId = 'pc-confirm-' + uid;
  const delBtn = h('button.btn.ghost.sm.pc-action.pc-del', {
    type: 'button',
    'aria-expanded': 'false',
    'aria-controls': confirmId,
    title: typeof hs.onDelete === 'function'
      ? 'Delete this project document. It does not delete its runs; the next '
        + 'step says exactly what is and is not removed.'
      : 'Unavailable: this page did not pass an onDelete handler.',
    disabled: typeof hs.onDelete !== 'function',
  }, 'delete');
  if (typeof hs.onDelete !== 'function') missing.push('delete');
  actions.appendChild(delBtn);

  if (missing.length) {
    actions.appendChild(h('p.pc-missing',
      missing.join(', ') + ' ' + (missing.length === 1 ? 'is' : 'are')
      + ' unavailable on this page: the handler was not wired. Open the '
      + 'project to do it in the Studio.'));
  }
  card.appendChild(actions);

  const confirmTitleId = confirmId + '-t';
  const confirmEl = h('div.pc-confirm', {
    id: confirmId,
    role: 'group',
    'aria-labelledby': confirmTitleId,
    hidden: true,
  },
    h('h4.pc-confirm-title', { id: confirmTitleId },
      'Delete ' + (name || pid || 'this project') + '?'),
    // WORKSPACE_CONTRACT.md §2. This paragraph is the product telling the
    // truth about a destructive action at the moment it is taken. It is not
    // shortened to fit, it is not moved into a title, and it is not delegated
    // to whatever dialog happens to be wired: it is rendered here so that it
    // cannot go missing.
    h('p.pc-confirm-text', deleteDisclosure(p)),
    h('div.pc-confirm-actions',
      h('button.btn.sm.pc-confirm-cancel', {
        type: 'button',
        onclick: () => closeConfirm(true),
      }, 'cancel'),
      h('button.btn.sm.danger.pc-confirm-go', {
        type: 'button',
        onclick: () => {
          closeConfirm(false);
          if (typeof hs.onDelete === 'function') invoke(card, delBtn, hs.onDelete, p);
        },
      }, 'delete the project, keep the runs')));
  card.appendChild(confirmEl);

  function openConfirm() {
    confirmEl.hidden = false;
    delBtn.setAttribute('aria-expanded', 'true');
    // Focus lands on cancel, not on the destructive button: the keyboard path
    // into a delete must not be Tab-Enter-gone.
    const cancel = confirmEl.querySelector('.pc-confirm-cancel');
    if (cancel) cancel.focus();
  }
  function closeConfirm(refocus) {
    confirmEl.hidden = true;
    delBtn.setAttribute('aria-expanded', 'false');
    if (refocus && delBtn.isConnected) delBtn.focus();
  }
  delBtn.addEventListener('click', () => {
    // DELEGATE WHEN THERE IS SOMEWHERE TO DELEGATE TO.
    //
    // This card's inline confirm and `createProject.js`'s delete dialog were
    // both built to §2, independently, and wiring them in series produced
    // THREE presses and the same guarantee stated twice -- found in a browser,
    // not by reading.
    //
    // The dialog wins for one substantive reason, not by seniority: it fetches
    // the project and reports the number of DISTINCT runs retained, while this
    // card only has `n_runs` from the summary, which is
    // `sum(len(run_ids))` -- run REFERENCES. Two specimens may hold the same
    // run id, so the card's count can overstate what survives, and the whole
    // point of the sentence is that the number is true.
    //
    // The inline confirm stays as the FALLBACK for a page that never wired a
    // handler, so the disclosure still cannot go missing -- which was the
    // reason it was written here in the first place.
    if (typeof hs.onDelete === 'function') {
      invoke(card, delBtn, hs.onDelete, p);
      return;
    }
    if (confirmEl.hidden) openConfirm(); else closeConfirm(true);
  });
  // Escape backs out, from anywhere inside the card. Listener is on the card,
  // not on the document, so a removed card takes its listener with it.
  card.addEventListener('keydown', (ev) => {
    if (ev.key === 'Escape' && !confirmEl.hidden) {
      ev.stopPropagation();
      closeConfirm(true);
    }
  });

  return card;
}

/** Replace a card's contents in place, keeping its node identity (and the
 *  scroll position of the grid around it) after a rename or a duplicate. */
export function refreshProjectCard(el, project, handlers) {
  const next = createProjectCard(project, handlers);
  el.className = next.className;
  for (const a of Array.from(next.attributes)) {
    if (a.name !== 'class') el.setAttribute(a.name, a.value);
  }
  clear(el);
  while (next.firstChild) el.appendChild(next.firstChild);
  return el;
}
