// browser.js — walk the server's own filesystem instead of reciting a path.
//
// WHAT THIS REPLACES
// A single <input type="text"> into which the operator typed an absolute path
// from memory. The bytes of a tissue export are read IN PLACE on the machine
// running the server (a Xenium bundle is gigabytes of matrices, and a section
// is registered and persisted by path), so a browser file picker cannot
// express the input at all — it can only offer to send a copy of data that is
// already on the right disk. This pane is the other way to solve it: the
// server lists its own directory names, and the person clicks.
//
// WHAT IT SHOWS, AND WHY THAT IS ALL
// Names and a type. `GET /api/browse` returns nothing else — no sizes, no
// timestamps, no owners, and no paths on the child entries. The one absolute
// path in a response is the directory the caller itself asked about. So this
// pane composes the next path from the current one plus a name, and a symlink
// is shown as a link and resolved by the SERVER when it is clicked, never
// here.

import { h, fill, clear, cls, debounce } from '../util/dom.js';
import { browse, ServerRefusal } from './net.js';

const PLATFORM_LABEL = {
  xenium: '10x Xenium',
  visium: 'Space Ranger / Visium',
  mex: 'MEX + coordinates',
  multi_sample: 'holds samples',
};

/** A recognised-bundle badge, or nothing. `samples` is keyed by NAME. */
function bundleChip(info) {
  if (!info) return null;
  const label = PLATFORM_LABEL[info.platform] || info.platform || 'sample';
  const conf = typeof info.confidence === 'number'
    ? ' ' + Math.round(info.confidence * 100) + '%' : '';
  return h('span.chip.bundle', { title: 'this directory looks like a ' + label
    + ' bundle' + (conf ? ' (detection confidence' + conf + ')' : '') },
    label + conf);
}

function kindGlyph(kind) {
  if (kind === 'dir') return '/';
  if (kind === 'link') return '->';
  if (kind === 'unknown') return '?';
  return '';
}

/**
 * The browser pane.
 *
 * @param {object} opts
 * @param {(state: object) => void} opts.onPlace called whenever the listed
 *        directory changes, with `{path, diagnosis, entries, confined}` —
 *        the right-hand pane renders from that and nothing else.
 * @param {(refusal: ServerRefusal) => void} opts.onRefusal a 401 or a 403 the
 *        page as a whole has to answer (the access key, the allowlist).
 */
export function createBrowser({ onPlace, onRefusal }) {
  const state = {
    path: '', q: '', showHidden: false,
    confined: null, places: [], crumbs: [], entries: [],
    diagnosis: null, truncated: false, nEntries: 0, busy: false,
  };

  const placesEl = h('div#tb-places.tb-places');
  const crumbsEl = h('nav#tb-crumbs.tb-crumbs', { 'aria-label': 'Directory path' });
  const listEl = h('div#tb-list.tb-list', { role: 'list' });
  const statusEl = h('p#tb-status.tb-status', { 'aria-live': 'polite' });

  const filterEl = h('input#tb-filter.inp.sm', {
    type: 'search', spellcheck: 'false',
    placeholder: 'filter these names',
    'aria-label': 'Filter the names in this directory',
  });
  const hiddenEl = h('input#tb-hidden', { type: 'checkbox' });
  const upEl = h('button#tb-up.btn.ghost.sm', {
    type: 'button', title: 'Go to the containing directory',
  }, 'Up');
  const reloadEl = h('button#tb-reload.btn.ghost.sm', {
    type: 'button', title: 'Read this directory again',
  }, 'Reload');

  const typedEl = h('input#tb-typed.inp', {
    type: 'text', spellcheck: 'false',
    placeholder: '/data/xenium_run_2026_08/output-XETG00001',
    'aria-label': 'Go to a directory by path',
  });
  const typedGo = h('button#tb-typed-go.btn.sm', { type: 'button' }, 'Go');

  const el = h('section#tb.tb', { 'aria-label': 'Tissue directories on this server' },
    h('div#tb-bar.tb-bar',
      upEl, reloadEl, filterEl,
      h('label.tb-check', hiddenEl, 'hidden entries')),
    placesEl,
    crumbsEl,
    listEl,
    statusEl,
    h('details#tb-typed-wrap.tb-typed',
      h('summary', 'Go straight to a path'),
      h('p.tb-hint',
        'A path on the machine running this server. Browsing above does the '
        + 'same thing without the typing; this is here for a path you already '
        + 'have in a terminal.'),
      h('div.tb-typed-row', typedEl, typedGo)));

  // --- rendering --------------------------------------------------------- //
  function renderPlaces() {
    const rows = state.places.map((p) => h('button.tb-place', {
      type: 'button', title: p.path,
      'aria-pressed': state.path === p.path ? 'true' : 'false',
      onclick: () => goto(p.path),
    }, h('span.tb-place-name', p.name), h('span.tb-place-why', p.why)));
    if (!rows.length) {
      // Set-but-empty TMESIM_STUDIO_DATA_ROOTS. The server refuses every path
      // in that state by design, and "no starting points" is a configuration
      // answer, not an empty disk.
      fill(placesEl, h('p.tb-hint',
        state.confined
          ? 'This server allows no tissue directory at all: its allowlist is '
            + 'set but names nothing usable, so every path is refused. Whoever '
            + 'runs it has to point TMESIM_STUDIO_DATA_ROOTS at the directory '
            + 'holding your export, and restart.'
          : 'No starting point could be offered. Use the path box below.'));
      return;
    }
    fill(placesEl, h('span.tb-places-label',
      state.confined ? 'Allowed on this server' : 'Start from'), rows);
  }

  function renderCrumbs() {
    if (!state.crumbs.length) { clear(crumbsEl); return; }
    const rows = [];
    state.crumbs.forEach((c, i) => {
      if (i) rows.push(h('span.tb-sep', '/'));
      const last = i === state.crumbs.length - 1;
      rows.push(last
        ? h('span.tb-crumb.here', { 'aria-current': 'location' }, c.name)
        : h('button.tb-crumb', {
          type: 'button', onclick: () => goto(c.path),
        }, c.name));
    });
    fill(crumbsEl, rows);
  }

  function renderList() {
    const samples = (state.diagnosis && state.diagnosis.samples) || {};
    if (!state.entries.length) {
      fill(listEl, h('p.tb-hint', state.q
        ? 'No name here matches ' + JSON.stringify(state.q) + '.'
        : state.path
          ? 'This directory is empty' + (state.showHidden ? '.'
            : ' apart from hidden entries.')
          : 'Pick a starting point above.'));
      return;
    }
    const rows = state.entries.map((e) => {
      const info = samples[e.name] || null;
      const nav = e.kind === 'dir' || e.kind === 'link';
      const bits = [
        h('span.tb-kind', { 'aria-hidden': 'true' }, kindGlyph(e.kind)),
        h('span.tb-name', e.name),
        bundleChip(info),
        info && info.registered_as
          ? h('span.chip.registered', { title: 'already a section on this '
            + 'server, as ' + info.registered_as }, 'registered')
          : null,
        e.kind === 'link'
          ? h('span.tb-note', 'link — the server resolves it when you open it')
          : null,
      ];
      if (!nav) return h('div.tb-row.file', { role: 'listitem' }, bits);
      return h('button.tb-row', {
        type: 'button', role: 'listitem',
        onclick: () => goto(join(state.path, e.name)),
      }, bits);
    });
    fill(listEl, rows);
  }

  function renderStatus(message) {
    const bits = [];
    if (message) bits.push(h('span.tb-status-msg', message));
    if (state.truncated) {
      bits.push(h('span.tb-status-msg',
        'Showing the first ' + state.entries.length + ' of ' + state.nEntries
        + ' matching names. Type in the filter to narrow it — the filter runs '
        + 'on the server, so it searches all of them.'));
    }
    if (!bits.length && state.path && state.nEntries) {
      bits.push(h('span.tb-status-msg', state.nEntries + ' entries'));
    }
    fill(statusEl, bits);
  }

  function render(message) {
    renderPlaces();
    renderCrumbs();
    renderList();
    renderStatus(message);
    upEl.disabled = !state.path || state.crumbs.length < 2;
    reloadEl.disabled = !state.path;
  }

  // --- navigation -------------------------------------------------------- //
  /** Compose a child path from the listed directory and a NAME.
   *
   *  Done here rather than returned by the server because the server returns
   *  no paths on child entries at all — see this file's header. A trailing
   *  separator on the parent (`/`) must not produce `//name`, which resolves
   *  the same on POSIX but reads like a bug in a breadcrumb. */
  function join(dir, name) {
    if (!dir) return name;
    return dir.endsWith('/') ? dir + name : dir + '/' + name;
  }

  async function load(path, { message = '' } = {}) {
    state.busy = true;
    cls(el, 'busy', true);
    try {
      const body = await browse({
        path, q: state.q, showHidden: state.showHidden,
      });
      state.path = body.path || '';
      state.confined = !!body.confined;
      state.places = body.places || [];
      state.crumbs = body.crumbs || [];
      state.entries = body.entries || [];
      state.diagnosis = body.diagnosis || null;
      state.truncated = !!body.truncated;
      state.nEntries = body.n_entries || 0;
      render(message);
      if (onPlace) onPlace(state);
    } catch (err) {
      if (!(err instanceof ServerRefusal)) throw err;
      // A refusal is an ANSWER here, and the server's own wording is the
      // answer: the 403 names the path that was sent, says the check was on
      // the resolved object and not on the text, and names the control. Any
      // rewording would make this page disagree with the server.
      state.entries = [];
      state.diagnosis = null;
      state.truncated = false;
      render(null);
      fill(statusEl, h('span.tb-refusal', err.message));
      if (onRefusal) onRefusal(err);
    } finally {
      state.busy = false;
      cls(el, 'busy', false);
    }
  }

  function goto(path) { return load(path); }

  const refilter = debounce(() => {
    state.q = filterEl.value;
    load(state.path);
  }, 180);

  filterEl.addEventListener('input', refilter);
  hiddenEl.addEventListener('change', () => {
    state.showHidden = hiddenEl.checked;
    load(state.path);
  });
  upEl.addEventListener('click', () => {
    if (state.crumbs.length >= 2) goto(state.crumbs[state.crumbs.length - 2].path);
  });
  reloadEl.addEventListener('click', () => load(state.path));
  typedGo.addEventListener('click', () => {
    const v = typedEl.value.trim();
    if (v) goto(v);
  });
  typedEl.addEventListener('keydown', (ev) => {
    if (ev.key === 'Enter') { ev.preventDefault(); typedGo.click(); }
  });

  return { el, state, goto, reload: () => load(state.path), start: () => load('') };
}
