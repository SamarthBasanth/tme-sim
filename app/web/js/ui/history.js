// history.js — the bottom dock: run history, per-run detail, A/B comparison.
//
// Two explicit radio columns pick the comparison pair. Selecting a run for
// comparison is deliberately NOT the same gesture as opening it in the
// viewer — conflating "look at this" with "compare against this" is how
// people end up reading the wrong pair of numbers.
//
// WHAT THIS PANEL IS FOR, restated after the 2026-09-20 pass. It was a list:
// 96 rows, newest first, no way to reach one except scrolling, and a row for a
// failed run was the same row as a successful one with four em-dashes where
// its numbers would be. It is now a place you can work:
//
//   * FIND       — search + status / section / compound filters, and a count
//                  that says how many of how many are on screen.
//   * READ       — a failed run states its failure ON THE ROW, and its detail
//                  carries the server's sentence in full, the traceback route
//                  and the run's own log tail.
//   * DISTINGUISH— a payload that was PRUNED is not a payload that was never
//                  produced. The server says which, in a 410 carrying the date
//                  and the byte count; both are rendered, and the 410 itself
//                  can be demanded from the server and is printed verbatim.
//   * ACT        — "the same run but dx = 30" is one action, on the row, with
//                  the changed field named before and after submission.
//
// THE SORT IS TIME OR STATUS, NEVER A METRIC. A run history that can be
// ordered by kill fraction is a compound ranking with a scroll bar, and
// ordering compounds by predicted effect was measured and withdrawn on
// 2026-08-21 (results/pk_ranking_robustness.md,
// results/capacity_ranking_robustness.md). There is no sort-by-metric, no
// "best run" badge, and no default that puts a number first.

import { h, clear } from '../util/dom.js';
import * as api from '../api.js';
import { toast, notifyError } from './toast.js';
import { METRICS, metricProvenance, chip, claimSetNotice } from './provenance.js';
import { auto, num, int, pct, sci, dur, ts, bytes, NULL_GLYPH, relPct } from '../util/fmt.js';

function fmtMetric(key, v) {
  if (v === null || v === undefined) return NULL_GLYPH;
  const m = METRICS[key];
  switch (m && m.fmt) {
    case 'pct': return pct(v, 2);
    case 'sci': return sci(v, 3);
    case 'num1': return num(v, 1);
    default: return auto(v);
  }
}
const metricUnit = (key) => (METRICS[key] ? (METRICS[key].fmt === 'pct' ? '%' : METRICS[key].unit) : '');

const SPEC_ROWS = [
  ['section', (s) => s.section],
  ['dx_um', (s) => num(s.dx_um, 1) + ' um'],
  ['drug', (s) => s.drug ?? 'none'],
  ['sensitivity_mode', (s) => s.sensitivity_mode],
  ['vascular_mode', (s) => s.vascular_mode],
  ['oxygen_mode', (s) => s.oxygen_mode],
  ['state_mode', (s) => s.state_mode ?? 'oxygen'],
  ['mgmt', (s) => (s.mgmt === null || s.mgmt === undefined ? 'null' : num(s.mgmt, 2) + ' log2CPM')],
  ['pk.plasma_peak', (s) => num(s.pk && s.pk.plasma_peak, 2) + ' uM'],
  ['pk.half_life_s', (s) => int(s.pk && s.pk.half_life_s) + ' s'],
  ['pk.n_doses', (s) => int(s.pk && s.pk.n_doses)],
  ['pk.dose_interval_s', (s) => int(s.pk && s.pk.dose_interval_s) + ' s'],
  ['total_time_s', (s) => int(s.total_time_s) + ' s'],
];

/** Column count of the run table — the empty-state row spans all of it. */
const COL_COUNT = 12;

/** The five job states, in LIFECYCLE order. This is the only ordering this
 *  panel knows other than time; it is a property of the queue, not of any
 *  result, which is exactly why it is allowed to exist. */
const STATUS_ORDER = ['running', 'queued', 'done', 'failed', 'cancelled'];

const STATUS_TITLE = {
  queued: 'on the queue — the worker has not started it, and no result exists yet',
  running: 'executing now — the numbers on this row are not its own until it is done',
  done: 'finished — the metrics on this row belong to this run',
  failed: 'produced no result. The server’s own sentence is on the row and in full under details',
  cancelled: 'stopped before it finished — no result was written',
};

/** The left-edge stripe that lets the five states be told apart at a glance,
 *  without reading. Colour is never the only channel: every row also carries
 *  its status chip in words, and a non-done row states its state in prose
 *  where its metrics would be. */
const STATUS_INK = {
  queued: 'var(--line-2)',
  running: 'var(--accent)',
  done: 'var(--ok)',
  failed: 'var(--crit)',
  cancelled: 'var(--fg-3)',
};

/** The two binary payloads a run record can account for. `geometry` and the
 *  four fields are not here: a run document carries no block for them, so this
 *  panel would have to guess, and guessing is the defect this whole section
 *  exists to avoid. */
const PAYLOADS = [
  {
    name: 'timelapse',
    label: 'drug time-lapse',
    route: 'GET /api/runs/{id}/timelapse',
    fetch: (id) => api.getTimelapse(id),
    // The time-lapse block states its own availability.
    present: (b) => b.available === true,
    reason: (r, b) => b.reason,
  },
  {
    name: 'mechanism',
    label: 'mechanism attribution',
    route: 'GET /api/runs/{id}/mechanism',
    fetch: (id) => api.getMechanism(id),
    // The mechanism block does NOT carry an `available` flag — its existence
    // IS the attribution (jobs.py `_mechanism_block` returns the block or
    // None, and puts the sentence in `mechanism_reason` beside it). Asking it
    // the time-lapse's question returned false for every run that HAS an
    // attribution, which is the same class of error as the prune 404: a
    // confident sentence pointed at the wrong object.
    present: (b) => !!b && typeof b === 'object',
    reason: (r) => r.mechanism_reason,
  },
];

/* ----------------------------------------------------- re-run lineage --- *
 * A LOCAL note, and it says so wherever it is rendered. The server records no
 * parent/child relationship between runs, so this is written here when the
 * re-run form is used and read back on the next load. It is a POINTER, never
 * evidence: what actually differs between two runs is the spec difference in
 * the Comparison panel, computed from the two RunSpecs themselves. */
const LINEAGE_KEY = 'vivum_rerun_lineage_v1';
const LINEAGE_MAX = 200;

function loadLineage() {
  try {
    const raw = JSON.parse(localStorage.getItem(LINEAGE_KEY) || '{}');
    return raw && typeof raw === 'object' ? raw : {};
  } catch { return {}; }
}

function saveLineage(map) {
  try {
    const ids = Object.keys(map);
    if (ids.length > LINEAGE_MAX) {
      ids.sort((x, y) => (map[x].at || 0) - (map[y].at || 0));
      for (const id of ids.slice(0, ids.length - LINEAGE_MAX)) delete map[id];
    }
    localStorage.setItem(LINEAGE_KEY, JSON.stringify(map));
  } catch { /* storage unavailable — the note is a convenience, not a record */ }
}

/* -------------------------------------------------------- the re-run --- *
 * NUMERIC PARAMETERS ONLY, and the omission is deliberate rather than lazy.
 * Section, compound and the four modes are held exactly as the parent run had
 * them: this form exists for the one comparison this product makes — the same
 * compound on the same specimen with one parameter moved — and the four mode
 * enumerations are owned by the RunSpec editor (js/ui/spec.js), which is where
 * a second copy of them would go stale. */
const RERUN_FIELDS = [
  { k: 'dx_um', path: ['dx_um'], label: 'dx', unit: 'um', step: '1', dp: 1,
    hint: 'grid pitch. Cost is quadratic in this field, and the server refuses above its grid ceiling.' },
  { k: 'total_time_s', path: ['total_time_s'], label: 'total time', unit: 's', step: '60', dp: 0,
    hint: 'simulated duration of the whole run.' },
  { k: 'mgmt', path: ['mgmt'], label: 'MGMT', unit: 'log2CPM', step: '0.1', dp: 2, nullable: true,
    hint: 'leave empty for null — null means "not supplied", and is not the same as zero.' },
  { k: 'pk.plasma_peak', path: ['pk', 'plasma_peak'], label: 'plasma peak', unit: 'uM', step: '1', dp: 2,
    hint: 'peak plasma concentration of one dose.' },
  { k: 'pk.half_life_s', path: ['pk', 'half_life_s'], label: 'half life', unit: 's', step: '600', dp: 0,
    hint: 'plasma elimination half-life.' },
  { k: 'pk.n_doses', path: ['pk', 'n_doses'], label: 'doses', unit: 'doses', step: '1', dp: 0,
    // The unit IS the label here; printing both gives "doses doses".
    noUnitLabel: true,
    hint: 'how many doses in the schedule.' },
  { k: 'pk.dose_interval_s', path: ['pk', 'dose_interval_s'], label: 'dose interval', unit: 's', step: '3600', dp: 0,
    hint: 'time between doses.' },
];

const getPath = (obj, path) => path.reduce((o, k) => (o === null || o === undefined ? o : o[k]), obj);

function setPath(obj, path, value) {
  let o = obj;
  for (let i = 0; i < path.length - 1; i += 1) {
    if (o[path[i]] === null || typeof o[path[i]] !== 'object') o[path[i]] = {};
    o = o[path[i]];
  }
  o[path[path.length - 1]] = value;
  return obj;
}

/** The field's value as the input shows it: '' for null/absent. */
function fieldText(spec, f) {
  const v = getPath(spec || {}, f.path);
  if (v === null || v === undefined || Number.isNaN(v)) return '';
  return String(v);
}

/** One line of prose for a value + its unit, with null spelled out. */
function fieldValue(f, raw) {
  if (raw === '' || raw === null || raw === undefined) return 'null';
  const n = Number(raw);
  if (!Number.isFinite(n)) return String(raw);
  return `${num(n, f.dp)} ${f.unit}`;
}

/** First sentence (or first `max` characters) — the row-sized form of a
 *  message whose full text is one keystroke away in the detail row and in the
 *  cell's own tooltip. Nothing is ever ONLY truncated. */
function firstSentence(text, max = 84) {
  const s = String(text || '').trim().replace(/\s+/g, ' ');
  if (!s) return '';
  const stop = s.indexOf('. ');
  if (stop > 0 && stop < max) return s.slice(0, stop + 1);
  if (s.length <= max) return s;
  // Cut on a word boundary: a mid-word cut reads as a typo, not as a cut.
  const sp = s.lastIndexOf(' ', max);
  return s.slice(0, sp > max * 0.6 ? sp : max).trim() + ' …';
}

/** An ApiError's server-side detail, as lines. FastAPI puts a whole object in
 *  `detail` for the interesting refusals (409 not-yet, 410 pruned), and
 *  api.js's `detailOf()` returns '' for an object — so `err.message` on a 410
 *  is the bare string "410 Gone" and the reason, the date and the byte count
 *  are reachable ONLY here. That is how the prune reason got swallowed once
 *  already. */
function serverDetail(err) {
  const d = err && err.payload && err.payload.detail;
  if (!d || typeof d !== 'object' || Array.isArray(d)) return null;
  return d;
}

export function createHistoryPanel({ onOpen, onRefresh }) {
  let runs = [];
  let activeId = null;
  let aId = null;
  let bId = null;
  /** At most one run's detail is open; keeping it to one is what stops the
   *  dock from becoming a wall of expanded rows in a 214 px column. */
  let expandedId = null;
  let lineage = loadLineage();
  /** Unsubmitted re-run edits, per run id, so a 6-second poll repaint does not
   *  eat what the operator is typing (main.js: setInterval(refreshRuns, 6000)). */
  const drafts = new Map();
  let lastSig = '';
  let lastCmpSig = '';

  const filters = { q: '', status: 'all', section: 'all', drug: 'all', sort: 'newest' };
  const filtersAreDefault = () => (
    !filters.q && filters.status === 'all' && filters.section === 'all'
    && filters.drug === 'all' && filters.sort === 'newest');

  const tbody = h('tbody');
  const countEl = h('span.hint', '0 runs');
  const cmpEl = h('div.panel-body.pad');

  /** The open run's detail, under the list rather than inside the table. It
   *  has its own scroll and a ceiling, so it never pushes the run list out of
   *  the panel. */
  const detailEl = h('div.run-detail', {
    id: 'run-detail-panel',
    style: {
      // `0 1 auto`, NOT `0 0 auto`: the dock is 214 px tall by default, and a
      // detail region that refuses to shrink takes 58% of that and leaves the
      // run list 4 px — a list you cannot see is not a list you can pick from.
      // It shrinks, scrolls inside itself, and the table below keeps a floor
      // of its own (`tableWrap`, two rows plus the header).
      flex: '1 1 auto',
      minHeight: '0',
      overflowY: 'auto',
      display: 'none',
      borderTop: '1px solid var(--line-2)',
      background: 'var(--bg-1)',
    },
  });

/* THE SORT LIVES ON TWO COLUMN HEADERS, and only two.
 *
 * It was a fourth `<select>` in the filter bar until the bar was measured
 * against the real Studio: 100 px of a 213 px dock, half the panel spent on
 * controls. Two headers cost nothing, and they put the refusal where it is
 * actually asked — hovering `kill / %` to sort by it is the moment to say why
 * this product does not. Every other column carries that sentence as its own
 * tooltip; these two carry `aria-sort` and a marker. */
const NOT_SORTABLE = 'This column is NOT sortable, and that is deliberate: ordering runs by an '
  + 'outcome is a ranking of what produced them, and ordering compounds by '
  + 'predicted effect was measured and withdrawn on 2026-08-21 '
  + '(results/pk_ranking_robustness.md, results/capacity_ranking_robustness.md). '
  + 'Sort by time or by status instead.';

  const sortMark = { submitted: h('span.unit', ''), status: h('span.unit', '') };

  function sortableTh(cls, label, key, describe) {
    return h(cls, {
      tabindex: '0',
      title: describe,
      'aria-label': `${label}, sortable column header`,
      style: { cursor: 'pointer' },
      onclick: () => cycleSort(key),
      onkeydown: (e) => {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); cycleSort(key); }
      },
    }, label, sortMark[key]);
  }

  function cycleSort(key) {
    if (key === 'submitted') filters.sort = filters.sort === 'newest' ? 'oldest' : 'newest';
    else filters.sort = filters.sort === 'status' ? 'newest' : 'status';
    render(true);
  }

  function paintSortMarks() {
    const th = sortMark.submitted.parentElement;
    const ts2 = sortMark.status.parentElement;
    sortMark.submitted.textContent = filters.sort === 'newest' ? '\u25be newest'
      : filters.sort === 'oldest' ? '\u25b4 oldest' : '';
    sortMark.status.textContent = filters.sort === 'status' ? '\u25be grouped' : '';
    if (th) {
      th.setAttribute('aria-sort', filters.sort === 'newest' ? 'descending'
        : filters.sort === 'oldest' ? 'ascending' : 'none');
    }
    if (ts2) ts2.setAttribute('aria-sort', filters.sort === 'status' ? 'other' : 'none');
  }

  const table = h('table.tbl.zebra',
    h('thead', h('tr',
      h('th', { style: { width: '22px' }, title: 'comparison slot A' }, 'A'),
      h('th', { style: { width: '22px' }, title: 'comparison slot B' }, 'B'),
      h('th', { title: 'The run id the server assigned. Not sortable: ids are issued in submission order, so sorting by one is sorting by time — use the submitted column, which says so.' }, 'run'),
      sortableTh('th', 'status', 'status',
        'Sort into queue order: running, queued, done, failed, cancelled. A '
        + 'property of the QUEUE, not of any result. Press again to return to time order.'),
      h('th', { title: 'The tissue section this run reconstructed. Filter by it above.' }, 'section'),
      h('th', { title: 'The compound this run simulated. Filter by it above.' }, 'drug'),
      h('th.right', { title: 'Grid pitch of the reconstruction, in micrometres. Not sortable; it is a setting, and the filter bar plus the re-run form are how you move it.' }, 'dx / um'),
      h('th.right', { title: NOT_SORTABLE }, 'kill / %'),
      h('th.right', { title: NOT_SORTABLE }, 'AUC / uM.s'),
      h('th.right', { title: NOT_SORTABLE }, 'pen / um'),
      h('th.right', { title: 'Wall-clock seconds this run took. Not sortable: it is a property of this machine and this queue, not of the run\u2019s result.' }, 'runtime'),
      sortableTh('th', 'submitted', 'submitted',
        'Sort by submission time. Press to flip between newest first and oldest first.'),
    )),
    tbody,
  );

  /** The list's own scroller. Held in a variable because the detail region
   *  below needs it to have a floor it cannot be squeezed past. */
  const tableWrap = h('div.panel-body.tblwrap', { style: { minHeight: '64px' } }, table);

  /* ==================================================== the filter bar === */

  const searchInp = h('input.inp', {
    type: 'search',
    placeholder: 'search runs',
    'aria-label': 'Search runs by id, section, compound or status',
    title: 'Matches the run id, the section, the compound, the status, and the '
      + 'word "pruned" for a run whose binary payload was deleted to reclaim disk.',
    style: { flex: '1 1 130px', minWidth: '96px' },
    oninput: () => { filters.q = searchInp.value.trim().toLowerCase(); render(true); },
  });

  /* A SELECT, not a segmented control, and the reason is measured. Six `.seg`
   * buttons carrying their counts are ~340 px wide; the run-history dock
   * column is 280 px at a 1440 px viewport (layout.css `.dock.with-reports`
   * floors: 280 + 240 + 300), and `.seg` is `overflow: hidden`. Screenshot QA
   * against the real Studio caught the result: the `cancelled` button was
   * clipped away entirely — a filter with 13 runs behind it that could not be
   * clicked and could not be scrolled to. Four uniform selects wrap instead of
   * clipping, and they cost two rows of dock height rather than four. */
  const statusSel = h('select.inp', {
    'aria-label': 'Filter the run list by status',
    title: 'Only runs in this state. A count of 0 means this server has none.',
    style: { flex: '1 1 84px', width: 'auto', minWidth: '76px' },
    onchange: () => { filters.status = statusSel.value; render(true); },
  });

  const sectionSel = h('select.inp', {
    'aria-label': 'Filter the run list by section',
    title: 'Only runs on this tissue section. Sections come from the runs in this list.',
    style: { flex: '1 1 84px', width: 'auto', minWidth: '76px' },
    onchange: () => { filters.section = sectionSel.value; render(true); },
  });

  const drugSel = h('select.inp', {
    'aria-label': 'Filter the run list by compound',
    title: 'Only runs of this compound. Filtering to ONE compound is also what '
      + 'makes the A/B metric delta below legal: a comparison that holds the '
      + 'compound fixed is the comparison this product makes.',
    style: { flex: '1 1 84px', width: 'auto', minWidth: '76px' },
    onchange: () => { filters.drug = drugSel.value; render(true); },
  });

  const clearBtn = h('button.btn.ghost.sm', {
    title: 'Clear the search and every filter, and return the sort to newest first.',
    onclick: () => {
      filters.q = ''; filters.status = 'all'; filters.section = 'all';
      filters.drug = 'all'; filters.sort = 'newest';
      searchInp.value = '';
      render(true);
    },
  }, 'clear filters');

  const filterBar = h('div.hist-filters', {
    style: {
      flex: '0 0 auto',
      display: 'flex',
      flexWrap: 'wrap',
      alignItems: 'center',
      gap: 'var(--sp-3) var(--sp-4)',
      padding: 'var(--sp-3) var(--sp-5)',
      borderBottom: '1px solid var(--line-1)',
      background: 'var(--bg-1)',
    },
  },
  searchInp, statusSel, sectionSel, drugSel, clearBtn);

  function paintStatusOptions() {
    const counts = { all: runs.length };
    for (const s of STATUS_ORDER) counts[s] = 0;
    for (const r of runs) {
      const s = r.status || 'queued';
      counts[s] = (counts[s] || 0) + 1;
    }
    const keep = filters.status;
    clear(statusSel);
    // Every state is listed even at zero: a state that disappears from the
    // list when it is empty makes "are there any failed runs?" unanswerable
    // without already knowing the answer.
    for (const key of ['all', ...STATUS_ORDER]) {
      statusSel.appendChild(h('option', {
        value: key,
        title: key === 'all' ? '' : STATUS_TITLE[key],
      }, key === 'all'
        ? `any status (${counts.all})`
        : `${key} (${counts[key] || 0})`));
    }
    statusSel.value = keep;
    if (statusSel.value !== keep) { filters.status = 'all'; statusSel.value = 'all'; }
    statusSel.title = filters.status === 'all'
      ? `Every run this server has, of any status — ${runs.length} run${runs.length === 1 ? '' : 's'}.`
      : `${counts[filters.status] || 0} of ${runs.length} run${runs.length === 1 ? '' : 's'}: `
        + STATUS_TITLE[filters.status];
  }

  function paintSelect(sel, values, allLabel, current) {
    const keep = current;
    clear(sel);
    sel.appendChild(h('option', { value: 'all' }, allLabel));
    for (const v of values) sel.appendChild(h('option', { value: v }, v));
    sel.value = values.includes(keep) ? keep : 'all';
    if (sel.value !== keep) return 'all';
    return keep;
  }

  function paintFilterOptions() {
    const sections = [...new Set(runs.map((r) => (r.spec && r.spec.section) || '').filter(Boolean))].sort();
    const drugs = [...new Set(runs.map((r) => (r.spec && (r.spec.drug ?? 'none')) || 'none'))].sort();
    filters.section = paintSelect(sectionSel, sections, `all sections (${sections.length})`, filters.section);
    filters.drug = paintSelect(drugSel, drugs, `all compounds (${drugs.length})`, filters.drug);
    paintStatusOptions();
    paintSortMarks();
    clearBtn.disabled = filtersAreDefault();
  }

  function haystack(r) {
    const s = r.spec || {};
    const pruned = Object.keys(r.pruned_payloads || {});
    const lin = lineage[r.run_id];
    return [
      r.run_id, r.status || 'queued', s.section, s.drug ?? 'none',
      pruned.length ? 'pruned ' + pruned.join(' ') : '',
      lin ? 're-run rerun of ' + lin.parent : '',
    ].filter(Boolean).join(' ').toLowerCase();
  }

  function visible() {
    let out = runs.filter((r) => {
      const s = r.spec || {};
      if (filters.status !== 'all' && (r.status || 'queued') !== filters.status) return false;
      if (filters.section !== 'all' && s.section !== filters.section) return false;
      if (filters.drug !== 'all' && (s.drug ?? 'none') !== filters.drug) return false;
      if (filters.q && !haystack(r).includes(filters.q)) return false;
      return true;
    });
    // `runs` arrives newest first (GET /api/runs). Nothing below sorts by a
    // metric, and nothing here may be made to.
    if (filters.sort === 'oldest') out = out.slice().reverse();
    else if (filters.sort === 'status') {
      out = out.slice().sort((x, y) => {
        const ix = STATUS_ORDER.indexOf(x.status || 'queued');
        const iy = STATUS_ORDER.indexOf(y.status || 'queued');
        return (ix < 0 ? 99 : ix) - (iy < 0 ? 99 : iy);
      });
    }
    return out;
  }

  /* ================================================== focus preservation = *
   * The run list repaints every 6 seconds whether or not anything changed
   * (main.js). Before this, that repaint destroyed every row in the table —
   * so a keyboard user's focus fell to <body> twice a minute, and the re-run
   * form below could not have existed at all. Two defences: the signature
   * below skips the repaint entirely when nothing changed, and this restores
   * focus (and the caret) when a repaint does happen. */
  function captureFocus() {
    const el = document.activeElement;
    if (!el || !(tbody.contains(el) || detailEl.contains(el))) return null;
    const row = el.closest ? el.closest('tr[data-run]') : null;
    return {
      fk: (el.dataset && el.dataset.fk) || null,
      act: (el.dataset && el.dataset.act) || null,
      run: row ? row.dataset.run : null,
      sel: el.selectionStart != null ? [el.selectionStart, el.selectionEnd] : null,
      isRow: el === row && !!row,
    };
  }

  function restoreFocus(f) {
    if (!f) return;
    const esc = (s) => (window.CSS && CSS.escape ? CSS.escape(s) : s);
    let el = null;
    if (f.fk) el = detailEl.querySelector(`[data-fk="${esc(f.fk)}"]`);
    else if (f.run && f.isRow) el = tbody.querySelector(`tr[data-run="${esc(f.run)}"]`);
    else if (f.run && f.act) el = tbody.querySelector(`tr[data-run="${esc(f.run)}"] [data-act="${esc(f.act)}"]`)
      || detailEl.querySelector(`[data-act="${esc(f.act)}"]`);
    else if (f.act) el = detailEl.querySelector(`[data-act="${esc(f.act)}"]`);
    if (!el) return;
    try { el.focus({ preventScroll: true }); } catch { el.focus(); }
    if (f.sel && el.setSelectionRange) {
      try { el.setSelectionRange(f.sel[0], f.sel[1]); } catch { /* not a text input */ }
    }
  }

  /** Everything the dock's DOM is a function of. Derived from the WHOLE run
   *  list rather than the filtered one so the filter controls are rebuilt on
   *  the same schedule as the rows — a poll that changes nothing must not
   *  close an open <select> or blow away the focus ring on a status button. */
  function signature() {
    return [
      JSON.stringify(filters), expandedId, activeId, aId, bId,
      runs.map((r) => [r.run_id, r.status, r.progress, r.stage, r.runtime_s,
        r.error ? 1 : 0,
        (r.spec && r.spec.section) || '', (r.spec && r.spec.drug) || '',
        Object.keys(r.pruned_payloads || {}).join('+')].join('~')).join('|'),
    ].join('#');
  }

  /* ========================================================== the rows === */

  function open(id) {
    activeId = id;
    render(true);
    onOpen && onOpen(id);
  }

  function setSlot(slot, id) {
    if (slot === 'A') aId = id; else bId = id;
    render(true);
    renderCompare(true);
  }

  function toggleDetail(id) {
    expandedId = expandedId === id ? null : id;
    render(true);
    const esc = (s) => (window.CSS && CSS.escape ? CSS.escape(s) : s);
    const btn = tbody.querySelector(`tr[data-run="${esc(id)}"] [data-act="details"]`);
    if (btn) { try { btn.focus({ preventScroll: true }); } catch { btn.focus(); } }
    // The detail region has its own scroll and lives under the table, so
    // opening one never moves the list — neither down nor (the bug this
    // replaced) sideways.
    detailEl.scrollTop = 0;
  }

  function focusRow(fromId, step) {
    const rows = Array.from(tbody.querySelectorAll('tr[data-run]'));
    const i = rows.findIndex((el) => el.dataset.run === fromId);
    const next = step === 'first' ? rows[0]
      : step === 'last' ? rows[rows.length - 1]
        : rows[i + step];
    if (next) next.focus();
  }

  function statusCell(r) {
    const st = r.status || 'queued';
    const live = st === 'running';
    return h('td',
      h('span.chip.st-' + st, { title: STATUS_TITLE[st] || st }, st),
      live
        ? h('div.progress', {
          style: { width: '38px', marginTop: 'var(--sp-1)' },
          title: `${pct(r.progress || 0, 0)}% of this run's stages are complete`,
        }, h('i', { style: { width: pct(r.progress || 0, 0) + '%' } }))
        : null);
  }

  /** The three metric cells — or, for a run that has no metrics, ONE cell in
   *  their place that says what happened instead. A failed run used to render
   *  as four em-dashes, which is the least legible thing on this panel and was
   *  attached to the row that most needed reading. */
  function outcomeCells(r) {
    const st = r.status || 'queued';
    const m = r.metrics || {};
    if (st === 'done') {
      return [
        h('td.right.num', fmtMetric('kill_fraction', m.kill_fraction)),
        h('td.right.num', fmtMetric('mean_tumor_auc', m.mean_tumor_auc)),
        h('td.right.num', fmtMetric('penetration_length_um', m.penetration_length_um)),
      ];
    }
    if (st === 'failed') {
      const said = r.error || '';
      return [h('td', {
        colspan: '3',
        title: said || 'The server reported this run as failed and sent no reason with it.',
      },
      h('span', { style: { color: 'var(--crit)' } },
        said ? firstSentence(said) : 'failed — the server sent no reason with it'))];
    }
    if (st === 'cancelled') {
      return [h('td', { colspan: '3' }, h('span.dim',
        `cancelled${r.stage ? ' at stage ' + r.stage : ''} — no result was written`))];
    }
    if (st === 'running') {
      return [h('td', { colspan: '3' }, h('span.dim',
        `${r.stage || 'running'} · ${pct(r.progress || 0, 0)}% of stages`))];
    }
    return [h('td', { colspan: '3' }, h('span.dim', 'queued — no result yet'))];
  }

  function runRow(r) {
    const st = r.status || 'queued';
    const s = r.spec || {};
    const isActive = r.run_id === activeId;
    const isOpen = expandedId === r.run_id;
    const prunedNames = Object.keys(r.pruned_payloads || {});
    const lin = lineage[r.run_id];

    // The stripe goes on the first cell, and is skipped on the active row:
    // components.css already draws `inset 2px 0 0 var(--accent)` there for
    // `aria-selected`, and an inline style would win over it and take the
    // selection cue away.
    const firstStyle = isActive ? {} : { boxShadow: `inset 2px 0 0 ${STATUS_INK[st] || 'var(--line-2)'}` };

    const tr = h('tr.sel-able', {
      tabindex: '0',
      dataset: { run: r.run_id, status: st },
      'aria-selected': String(isActive),
      'aria-label': `run ${r.run_id}, ${st}`,
      onclick: (e) => {
        const t = e.target;
        if (t && (t.tagName === 'INPUT' || t.tagName === 'BUTTON')) return;
        open(r.run_id);
      },
      onkeydown: (e) => {
        if (e.target !== tr) return;   // let the radios and buttons keep their keys
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(r.run_id); }
        else if (e.key === 'ArrowDown') { e.preventDefault(); focusRow(r.run_id, 1); }
        else if (e.key === 'ArrowUp') { e.preventDefault(); focusRow(r.run_id, -1); }
        else if (e.key === 'Home') { e.preventDefault(); focusRow(r.run_id, 'first'); }
        else if (e.key === 'End') { e.preventDefault(); focusRow(r.run_id, 'last'); }
        else if (e.key === 'ArrowRight' && !isOpen) { e.preventDefault(); toggleDetail(r.run_id); }
        else if (e.key === 'ArrowLeft' && isOpen) { e.preventDefault(); toggleDetail(r.run_id); }
      },
    },
    h('td', { style: firstStyle }, h('input', {
      type: 'radio', name: 'cmp-a', checked: aId === r.run_id,
      'aria-label': `set ${r.run_id} as comparison A`,
      onchange: () => { aId = r.run_id; renderCompare(true); },
    })),
    h('td', h('input', {
      type: 'radio', name: 'cmp-b', checked: bId === r.run_id,
      'aria-label': `set ${r.run_id} as comparison B`,
      onchange: () => { bId = r.run_id; renderCompare(true); },
    })),
    // THE DISCLOSURE IS HERE, not in a trailing column. The run table is ~700 px
    // of natural width inside a 280 px dock column (layout.css keeps it at
    // `min-width: max-content` and scrolls), so a control in the LAST column is
    // off screen until you scroll sideways — measured on the real Studio. The
    // leftmost cells are the ones always visible, so the one new affordance
    // this panel adds lives in one.
    h('td.mono',
      h('button.btn.tiny', {
        dataset: { act: 'details' },
        'aria-expanded': String(isOpen),
        'aria-controls': 'run-detail-panel',
        'aria-label': `${isOpen ? 'Hide' : 'Show'} the detail of run ${r.run_id}`,
        style: { marginRight: 'var(--sp-3)' },
        title: isOpen
          ? 'Collapse this run\u2019s detail'
          : 'Payload availability (including what was pruned), the failure in full, '
            + 'and a re-run of this exact spec with one parameter moved',
        onclick: (e) => { e.stopPropagation(); toggleDetail(r.run_id); },
      }, isOpen ? '\u25be' : '\u25b8'),
      r.run_id,
      prunedNames.length
        ? h('span.chip.st-pruned', {
          style: { color: 'var(--warn)', background: 'var(--warn-bg)', marginLeft: 'var(--sp-3)' },
          title: `${prunedNames.join(' + ')} was produced by this run and later deleted to `
            + 'reclaim disk. Open details for the date, the byte count and how to reproduce it.',
        }, 'pruned')
        : null,
      lin
        ? h('button.btn.tiny', {
          dataset: { act: 'lineage' },
          style: { marginLeft: 'var(--sp-3)' },
          title: `Re-run of ${lin.parent}, noted by THIS browser${lin.at ? ' on ' + ts(lin.at) : ''}`
            + `${lin.changes && lin.changes.length ? ' (' + lin.changes.map((c) => c.text).join('; ') + ')' : ''}. `
            + 'Press to set A = ' + lin.parent + ' and B = ' + r.run_id + '; the spec difference '
            + 'in the Comparison panel is then computed from the two RunSpecs themselves.',
          onclick: (e) => { e.stopPropagation(); aId = lin.parent; bId = r.run_id; render(true); renderCompare(true); },
        }, 'from ' + lin.parent)
        : null),
    statusCell(r),
    h('td', { title: s.section }, s.section || NULL_GLYPH),
    h('td', (s.drug ?? 'none') || NULL_GLYPH),
    h('td.right.num', num(s.dx_um, 1)),
    ...outcomeCells(r),
    h('td.right.num', { title: st === 'running' ? 'elapsed so far' : 'wall-clock seconds of this run' },
      dur(r.runtime_s)),
    h('td.mono', { style: { fontSize: 'var(--fs-micro)' } }, ts(r.created_at)),
    );
    return tr;
  }

  /* ======================================================= the detail === */

  const lbl = (text) => h('div.hr-lbl', text);
  const note = (...kids) => h('div', {
    style: {
      fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)', color: 'var(--fg-1)',
    },
  }, ...kids);

  /* --- failure --------------------------------------------------------- */
  function failureBlock(r) {
    const out = h('div', { style: { display: 'flex', flexDirection: 'column', gap: 'var(--sp-3)' } });
    const btn = h('button.btn.ghost.sm', {
      dataset: { act: 'trace' },
      title: 'GET /api/runs/{id}/traceback and GET /api/runs/{id}/status — the '
        + 'typed error, the stack if this server serves one, and the run’s own log tail.',
      onclick: (e) => { e.stopPropagation(); loadFailureDetail(r.run_id, out, btn); },
    }, 'fetch traceback and log');

    return [
      lbl('Why this run produced no result'),
      // The server's sentence, verbatim and in full. No clamp, no "show more":
      // a deliberate refusal names what it wanted, and a redacted internal
      // failure carries the reference id that is the only way to reach the
      // detail behind it.
      h('div.caveat.crit', h('span', r.error
        ? [h('b', 'The server said: '), r.error]
        : 'The server reported this run as failed and sent no reason with it, so this '
          + 'Studio has nothing of its account to show. The failure is recorded '
          + 'server-side: read the console of the process running Vivum.')),
      note('No binary payload exists for this run, and that is a THIRD state '
        + 'again: payloads are written when a run completes, so a request for '
        + 'one answers ', h('b', 'HTTP 409'), ' naming this run’s status — not '
        + '404 (never produced) and not 410 (produced, then deleted).'),
      h('div.btn-row', btn),
      out,
    ];
  }

  async function loadFailureDetail(runId, out, btn) {
    clear(out);
    btn.disabled = true;
    out.appendChild(h('span.hint', 'asking the server…'));
    try {
      if (api.api.stub) {
        clear(out);
        out.appendChild(h('div.caveat', h('span',
          'The traceback route needs the live backend. js/mock.js does not '
          + 'synthesise one, and a synthesised stack would be a fabricated '
          + 'diagnosis of a failure that did not happen. Start the server with '
          + './.venv/bin/python -m app.server.')));
        return;
      }
      const headers = { accept: 'application/json' };
      if (api.api.apiKey) headers['x-api-key'] = api.api.apiKey;
      const res = await fetch(`${api.api.base}/api/runs/${encodeURIComponent(runId)}/traceback`,
        { cache: 'no-store', headers });
      if (!res.ok) throw new Error(`${res.status} ${res.statusText}`);
      const tb = await res.json();
      let status = null;
      try { status = await api.getStatus(runId); } catch { status = null; }

      clear(out);
      if (tb.error) {
        out.appendChild(note(h('b', 'typed error: '), h('code', tb.error)));
      }
      if (tb.traceback) {
        out.appendChild(note('Stack, with this server’s repo root replaced by ',
          h('code', '<repo>'), ':'));
        out.appendChild(h('pre', {
          style: {
            maxHeight: '180px', overflow: 'auto', fontSize: 'var(--fs-micro)',
            background: 'var(--bg-0)', border: '1px solid var(--line-1)',
            padding: 'var(--sp-4)', margin: '0', whiteSpace: 'pre',
          },
        }, tb.traceback));
      } else if (tb.traceback_withheld) {
        // Verbatim: it already names the environment variable and the console.
        out.appendChild(h('div.caveat', h('span', tb.traceback_withheld)));
      }
      const log = (status && status.log_tail) || [];
      if (log.length) {
        out.appendChild(note(`The run’s own log, last ${log.length} line${log.length === 1 ? '' : 's'} `
          + 'of what the worker recorded (GET /api/runs/{id}/status):'));
        out.appendChild(h('pre', {
          style: {
            maxHeight: '140px', overflow: 'auto', fontSize: 'var(--fs-micro)',
            background: 'var(--bg-0)', border: '1px solid var(--line-1)',
            padding: 'var(--sp-4)', margin: '0', whiteSpace: 'pre-wrap',
          },
        }, log.join('\n')));
      } else if (status) {
        out.appendChild(note('This run kept no log lines.'));
      }
    } catch (err) {
      clear(out);
      out.appendChild(h('div.caveat.crit', h('span',
        h('b', 'GET /api/runs/' + runId + '/traceback failed: '),
        (err && err.message) || String(err))));
    } finally {
      btn.disabled = false;
    }
  }

  /* --- payloads -------------------------------------------------------- */
  /**
   * PRUNED is not ABSENT, and this is the only place in the Studio that says
   * so. `scripts/prune_run_payloads.py` deletes a run's `timelapse.bin` to
   * reclaim disk and writes what it removed into the run's own record;
   * `app/server/api.py::_binary` then answers 410 Gone with the date, the byte
   * count and how to get it back, INSTEAD of the 404 whose reason ("no
   * time-lapse was captured for this run") is false for that run. Rendering
   * both as one grey dash would put the false version back on screen.
   */
  function payloadState(r, p) {
    const rec = (r.pruned_payloads || {})[p.name];
    const block = r[p.name];
    if (rec && typeof rec === 'object') {
      return { kind: 'pruned', at: rec.at, nbytes: rec.bytes, reason: rec.reason };
    }
    if (block && typeof block === 'object' && block.pruned) {
      return { kind: 'pruned', at: null, nbytes: null, reason: block.reason };
    }
    const why = p.reason(r, block || {});
    if (block && typeof block === 'object' && p.present(block)) {
      const bits = [];
      if (block.n_frames) {
        bits.push(`${block.n_frames} of ${block.requested_frames ?? block.n_frames} frames requested`);
      }
      if (Number.isFinite(block.payload_bytes)) {
        bits.push(`${int(block.payload_bytes)} bytes (${bytes(block.payload_bytes)})`);
      }
      return { kind: 'available', detail: bits.join(' · ') };
    }
    if (block !== undefined && (why || block !== null)) {
      return { kind: 'absent', reason: why || 'this run recorded no reason' };
    }
    if (why) return { kind: 'absent', reason: why };
    return {
      kind: 'unknown',
      reason: `this run record carries no ${p.name} entry at all, which is what a run `
        + 'that predates the field looks like. Nothing here says it was or was not produced.',
    };
  }

  function payloadRow(r, p) {
    const st = payloadState(r, p);
    const out = h('div', { style: { display: 'flex', flexDirection: 'column', gap: 'var(--sp-3)' } });
    const chipStyle = {
      pruned: { color: 'var(--warn)', background: 'var(--warn-bg)' },
      available: { color: 'var(--ok)', background: 'var(--ok-bg)' },
      absent: { color: 'var(--fg-3)', background: 'transparent' },
      unknown: { color: 'var(--fg-3)', background: 'transparent' },
    }[st.kind];
    const words = {
      pruned: 'pruned', available: 'available', absent: 'not produced', unknown: 'not recorded',
    }[st.kind];

    const head = h('div', {
      style: { display: 'flex', alignItems: 'center', gap: 'var(--sp-4)', flexWrap: 'wrap' },
    },
    h('span', { style: { fontSize: 'var(--fs-xs)', color: 'var(--fg-1)' } }, p.label),
    h('span.chip', { style: chipStyle, title: p.route.replace('{id}', r.run_id) }, words));

    if (st.kind === 'pruned' || st.kind === 'unknown') {
      head.appendChild(h('button.btn.tiny', {
        dataset: { act: 'probe-' + p.name },
        title: `${p.route.replace('{id}', r.run_id)} — ask for it and print what the server answers, verbatim.`,
        onclick: (e) => { e.stopPropagation(); probePayload(r.run_id, p, out); },
      }, 'ask the server for it'));
    }

    const body = st.kind === 'pruned'
      ? note(h('b', 'This run DID produce it. '),
        st.at
          ? `The file was deleted on ${st.at} to reclaim disk: ${int(st.nbytes)} bytes `
            + `(${bytes(st.nbytes)}). `
          : 'The file was deleted to reclaim disk; this record does not carry the date. ',
        'Re-running the RunSpec in this record reproduces it — the spec is complete, '
        + 'and the re-run form below starts from it. A request for it answers ',
        h('b', 'HTTP 410 Gone'), ', not 404: the difference is that it existed.',
        st.reason ? [' The record’s own note: ', h('em', st.reason)] : null)
      : st.kind === 'available'
        ? note('Served by ', h('code', p.route.replace('{id}', r.run_id)),
          st.detail ? ' · ' + st.detail : '')
        : st.kind === 'absent'
          ? note(h('b', 'Never produced by this run, and not deleted. '),
            'A request for it answers ', h('b', 'HTTP 404'), ' with the run’s own '
            + 'reason, which is: ', h('em', st.reason))
          // Not a third phrasing of "absent": this record cannot tell you
          // either way, and saying "never produced" here would be the same
          // confident-sentence-wrong-object mistake as the 404 the prune
          // record exists to prevent.
          : note(h('b', 'This record does not say. '), st.reason,
            ' Asking the server is the only way to find out — it answers from the '
            + 'filesystem, not from this document.');

    return h('div', {
      style: {
        display: 'flex', flexDirection: 'column', gap: 'var(--sp-2)',
        padding: 'var(--sp-3) 0', borderBottom: '1px solid var(--line-0)',
      },
    }, head, body, out);
  }

  async function probePayload(runId, p, out) {
    clear(out);
    out.appendChild(h('span.hint', `asking ${p.route.replace('{id}', runId)}…`));
    try {
      const got = await p.fetch(runId);
      clear(out);
      out.appendChild(h('div.caveat.info', h('span',
        'The server SERVED it. This run record’s prune note is therefore stale '
        + `(${(got && got.bytes) ? int(got.bytes) + ' bytes received' : 'payload received'}) `
        + '— press refresh runs to re-read the list.')));
    } catch (err) {
      const d = serverDetail(err);
      clear(out);
      out.appendChild(h('div.caveat', h('span',
        h('b', `HTTP ${(err && err.status) || '?'} · `),
        d
          ? [d.message || '', d.reason ? ' — ' + d.reason : '',
            d.recoverable ? [' ', h('b', 'To get it back: '), d.recoverable] : null]
          : ((err && err.message) || String(err)))));
    }
  }

  /* --- the re-run composer -------------------------------------------- */
  function draftFor(runId, spec) {
    if (!drafts.has(runId)) {
      const d = {};
      for (const f of RERUN_FIELDS) d[f.k] = fieldText(spec, f);
      drafts.set(runId, d);
    }
    return drafts.get(runId);
  }

  function changeList(r) {
    const d = drafts.get(r.run_id) || {};
    const out = [];
    for (const f of RERUN_FIELDS) {
      const was = fieldText(r.spec, f);
      const now = (d[f.k] ?? '').trim();
      if (now === was) continue;
      if (Number(now) === Number(was) && now !== '' && was !== '') continue;
      out.push({
        field: f.k,
        from: fieldValue(f, was),
        to: fieldValue(f, now),
        text: `${f.k} ${fieldValue(f, was)} → ${fieldValue(f, now)}`,
      });
    }
    return out;
  }

  function specWithDraft(r) {
    const spec = JSON.parse(JSON.stringify(r.spec || {}));
    const d = drafts.get(r.run_id) || {};
    for (const f of RERUN_FIELDS) {
      const raw = (d[f.k] ?? '').trim();
      if (raw === '') {
        if (f.nullable) setPath(spec, f.path, null);
        continue;
      }
      const n = Number(raw);
      if (!Number.isFinite(n)) continue;
      setPath(spec, f.path, n);
    }
    return spec;
  }

  function rerunBlock(r) {
    const d = draftFor(r.run_id, r.spec);
    const summary = h('div');
    const valEl = h('div', { style: { display: 'flex', flexDirection: 'column', gap: 'var(--sp-2)' } });
    const rows = [];

    const grid = h('div', {
      style: {
        display: 'grid',
        gridTemplateColumns: 'minmax(78px, max-content) 86px minmax(0, 1fr)',
        gap: 'var(--sp-2) var(--sp-4)',
        alignItems: 'center',
      },
    });

    const submitBtn = h('button.btn.primary.sm', {
      dataset: { act: 'rerun' },
      title: 'POST /api/validate, then POST /api/runs with this run’s spec and your edits.',
      onclick: (e) => { e.stopPropagation(); submitRerun(r, submitBtn, valEl); },
    }, 'submit re-run');

    function paintChanges() {
      const changed = changeList(r);
      const changedKeys = new Set(changed.map((c) => c.field));
      for (const row of rows) {
        const on = changedKeys.has(row.f.k);
        row.inp.style.borderColor = on ? 'var(--accent)' : '';
        clear(row.was);
        row.was.appendChild(on
          ? h('span', { style: { color: 'var(--accent)', fontSize: 'var(--fs-xs)' } },
            `was ${fieldValue(row.f, fieldText(r.spec, row.f))}`)
          : h('span', { style: { color: 'var(--fg-2)', fontSize: 'var(--fs-xs)' } },
            fieldValue(row.f, fieldText(r.spec, row.f))));
      }
      clear(summary);
      summary.appendChild(changed.length
        ? note(h('b', { style: { color: 'var(--accent)' } },
          `${changed.length} field${changed.length === 1 ? '' : 's'} changed: `),
        changed.map((c) => c.text).join(' · '),
        `. Same compound (${r.spec && (r.spec.drug ?? 'none')}), same section `
          + `(${(r.spec && r.spec.section) || NULL_GLYPH}).`)
        : note('No field changed yet. Submitting now would re-run this spec '
          + 'unaltered — which is a determinism check, not a comparison.'));
      submitBtn.textContent = changed.length ? 'submit re-run' : 'submit identical re-run';
    }

    for (const f of RERUN_FIELDS) {
      const was = h('div');
      const inp = h('input.inp.num', {
        type: 'number',
        step: f.step,
        value: d[f.k] ?? '',
        dataset: { fk: `${r.run_id}:${f.k}` },
        'aria-label': `${f.label} for a re-run of ${r.run_id}, in ${f.unit}`,
        title: f.hint,
        // No stopPropagation: main.js's window key handler already ignores
        // events whose target is an INPUT, and swallowing them here would take
        // Ctrl+Enter (submit the RunSpec) away from anyone whose cursor
        // happens to be in this field.
        oninput: () => { d[f.k] = inp.value; paintChanges(); },
      });
      rows.push({ f, inp, was });
      grid.appendChild(h('label', {
        style: { fontSize: 'var(--fs-xs)', color: 'var(--fg-1)' },
        title: f.hint,
      }, f.label, f.noUnitLabel ? null : h('span.unit', f.unit)));
      grid.appendChild(inp);
      grid.appendChild(was);
    }

    paintChanges();

    return [
      lbl('Re-run this spec with one parameter moved'),
      note('Starts from ', h('b', r.run_id), '’s own RunSpec. ',
        h('b', 'Section, compound and the four modes are held fixed here'),
        ' — change those in the RunSpec editor on the left rail. Holding the '
        + 'compound fixed is what makes the A/B metric delta below legal: a '
        + 'per-compound scale factor cancels out of a comparison that does not '
        + 'change the compound, and does not cancel out of one that does.'),
      grid,
      summary,
      h('div.btn-row',
        submitBtn,
        h('button.btn.ghost.sm', {
          dataset: { act: 'rerun-reset' },
          title: 'Put every field back to the value this run used.',
          onclick: (e) => {
            e.stopPropagation();
            for (const row of rows) {
              d[row.f.k] = fieldText(r.spec, row.f);
              row.inp.value = d[row.f.k];
            }
            clear(valEl);
            paintChanges();
          },
        }, 'reset fields')),
      valEl,
    ];
  }

  async function submitRerun(r, btn, valEl) {
    const changes = changeList(r);
    const spec = specWithDraft(r);
    btn.disabled = true;
    clear(valEl);
    valEl.appendChild(h('span.hint', 'validating…'));
    try {
      let verdict = null;
      try {
        verdict = await api.postValidate(spec);
      } catch {
        // Non-fatal: POST /api/runs validates the same spec server-side and
        // refuses with the same errors. Losing the preview is not losing the
        // guard.
        verdict = null;
      }
      clear(valEl);
      if (verdict && verdict.ok === false) {
        for (const e of verdict.errors || []) {
          valEl.appendChild(h('div.caveat.crit', h('span', h('b', 'refused: '), e)));
        }
        if (!(verdict.errors || []).length) {
          valEl.appendChild(h('div.caveat.crit', h('span',
            'POST /api/validate refused this spec and named no field.')));
        }
        return;
      }
      for (const w of (verdict && verdict.warnings) || []) {
        valEl.appendChild(h('div.caveat', h('span', h('b', 'warning: '), w)));
      }
      const res = await api.postRun(spec);
      const newId = res && res.run_id;
      if (newId) {
        recordLineage(newId, r.run_id, changes);
        aId = r.run_id;
        bId = newId;
      }
      drafts.delete(r.run_id);
      expandedId = null;
      toast('ok', `Run ${newId || '?'} queued`,
        (changes.length
          ? `Re-run of ${r.run_id} with ${changes.map((c) => c.text).join(', ')}. `
          : `Re-run of ${r.run_id} with no field changed. `)
        + `A is ${r.run_id} and B is ${newId}, so the Comparison panel holds their `
        + 'spec difference now and their metric delta when this one finishes.',
        7000);
      render(true);
      renderCompare(true);
      if (onRefresh) await onRefresh();
    } catch (err) {
      notifyError('POST /api/runs failed', err);
    } finally {
      btn.disabled = false;
    }
  }

  function recordLineage(childId, parentId, changes) {
    lineage = { ...lineage, [childId]: { parent: parentId, at: Date.now() / 1000, changes } };
    saveLineage(lineage);
  }

  function detailBody(r) {
    const st = r.status || 'queued';
    const lin = lineage[r.run_id];
    const kids = [];

    if (st === 'failed') kids.push(...failureBlock(r));

    if (st === 'cancelled') {
      kids.push(lbl('What this run left behind'));
      kids.push(h('div.caveat', h('span',
        `This run was cancelled${r.stage ? ' at stage ' + r.stage : ''}. A cancelled run `
        + 'writes no result and no binary payload, so there is nothing of it to read. '
        + 'Its RunSpec is intact and is the starting point of the re-run form below.')));
    }

    if (st === 'queued' || st === 'running') {
      kids.push(lbl('This run is still in flight'));
      kids.push(note(st === 'running'
        ? `Stage ${r.stage || 'unknown'} · ${pct(r.progress || 0, 0)}% of this run’s `
          + `stages · ${dur(r.runtime_s)} elapsed. Its payloads are written when it `
          + 'completes, so a request for one answers HTTP 409 with this status until then.'
        : 'It is on the queue and the worker has not started it. Nothing has been '
          + 'computed, so every metric column on its row is empty rather than zero, '
          + 'and a request for any payload answers HTTP 409 with this status.'));
    }

    if (st === 'done' || Object.keys(r.pruned_payloads || {}).length) {
      kids.push(lbl('Binary payloads'));
      for (const p of PAYLOADS) kids.push(payloadRow(r, p));
    }

    kids.push(...rerunBlock(r));

    if (lin) {
      kids.push(h('div.caveat.info', h('span',
        h('b', `Re-run of ${lin.parent}. `),
        `Recorded by THIS browser${lin.at ? ' on ' + ts(lin.at) : ''} when the re-run form `
        + 'was used'
        + (lin.changes && lin.changes.length
          ? `, changing ${lin.changes.map((c) => c.text).join('; ')}` : ', changing nothing')
        + '. It is a local note and not a server-side link: the server records no '
        + 'run lineage, and a run submitted from another browser or from the CLI '
        + 'carries none. What actually differs between the two is the spec '
        + 'difference in the Comparison panel, computed from the two RunSpecs.')));
    }

    kids.push(h('div.btn-row',
      h('button.btn.ghost.sm', {
        dataset: { act: 'open' },
        title: 'Load this run’s geometry and result into the viewer and the READ rail',
        onclick: (e) => { e.stopPropagation(); open(r.run_id); },
      }, 'open in viewer'),
      h('button.btn.ghost.sm', {
        dataset: { act: 'set-a' },
        title: 'Use this run as comparison A',
        onclick: (e) => { e.stopPropagation(); setSlot('A', r.run_id); },
      }, 'set as A'),
      h('button.btn.ghost.sm', {
        dataset: { act: 'set-b' },
        title: 'Use this run as comparison B',
        onclick: (e) => { e.stopPropagation(); setSlot('B', r.run_id); },
      }, 'set as B')));

    return kids;
  }

  /**
   * Paint the open run's detail into its own region UNDER the list.
   *
   * It was a `<tr colspan=13>` inside the table for one release-candidate
   * afternoon, and screenshot QA against the real Studio killed that: the
   * dock's table is `min-width: max-content` and scrolls sideways in a 280 px
   * column, so a detail cell 520 px wide (a) could not wrap to the column, and
   * (b) took the whole dock's height the moment it opened, leaving the run
   * list — the thing you are choosing FROM — entirely off screen. Below the
   * table it wraps to the real column width, keeps its own scroll, and the
   * list stays visible above it.
   */
  function paintDetail(list) {
    const r = expandedId ? list.find((x) => x.run_id === expandedId) : null;
    clear(detailEl);
    if (!r) {
      // A run filtered out of the list closes its own detail rather than
      // hanging on as a panel about a row nobody can see.
      if (expandedId && !list.some((x) => x.run_id === expandedId)) expandedId = null;
      detailEl.style.display = 'none';
      tableWrap.style.flex = '';
      tableWrap.style.height = '';
      return;
    }
    detailEl.style.display = '';
    // The list keeps a share of the panel with a hard floor (header + 2 rows);
    // the detail takes everything else. Measured at the DEFAULT dock height
    // (--dock-h: 214px): without this the detail resolved to 27 px, which is
    // a disclosure you cannot read, attached to a button that claims to open
    // one. Drag the dock taller and both grow.
    tableWrap.style.flex = '0 0 auto';
    tableWrap.style.height = '30%';
    detailEl.appendChild(h('div.panel-head',
      h('h2', 'Run detail'),
      h('span.mono.strong', r.run_id),
      h('span.chip.st-' + (r.status || 'queued'),
        { title: STATUS_TITLE[r.status || 'queued'] }, r.status || 'queued'),
      h('span.grow'),
      h('button.btn.ghost.sm', {
        dataset: { act: 'close-detail' },
        title: 'Close this run’s detail',
        onclick: () => toggleDetail(r.run_id),
      }, 'close')));
    detailEl.appendChild(h('div', {
      style: {
        display: 'flex', flexDirection: 'column', gap: 'var(--sp-4)',
        padding: 'var(--sp-4) var(--sp-5) var(--sp-6)',
      },
    }, ...detailBody(r)));
  }

  /* ========================================================== render ==== */

  function render(force) {
    const sig = signature();
    if (!force && sig === lastSig) return;
    lastSig = sig;
    // After the guard, never before: `GET /api/runs` is re-read every 6 s
    // (main.js), and rebuilding the three selects on every one of those polls
    // closed an open dropdown and dropped the focus ring twice a minute for
    // no change at all.
    paintFilterOptions();
    const list = visible();

    const focus = captureFocus();
    clear(tbody);

    countEl.textContent = list.length === runs.length
      ? `${runs.length} run${runs.length === 1 ? '' : 's'}`
      : `${list.length} of ${runs.length} run${runs.length === 1 ? '' : 's'}`;
    countEl.title = list.length === runs.length
      ? 'Every run this server has.'
      : 'Filtered. Press clear filters to see every run this server has.';

    if (!list.length) {
      tbody.appendChild(h('tr', h('td', { colspan: String(COL_COUNT), class: 'dim' },
        runs.length
          ? [`No run matches this filter (${runs.length} run${runs.length === 1 ? '' : 's'} in total). `,
            h('button.btn.tiny', { onclick: () => clearBtn.click() }, 'clear filters')]
          : 'No runs yet. Submit a RunSpec from the left rail.')));
      paintDetail(list);
      restoreFocus(focus);
      return;
    }

    for (const r of list) {
      tbody.appendChild(runRow(r));
    }
    paintDetail(list);
    restoreFocus(focus);
  }

  /* The caveat union, lifted out of renderCompare so the cross-compound
   * refusal path can render it too. A comparison that declines to subtract
   * still owes the reader both runs' caveats — the refusal removes a number,
   * never a disclosure. */
  function renderCaveatDiff(A, B) {
    const ca = (A.provenance && A.provenance.caveats) || [];
    const cb = (B.provenance && B.provenance.caveats) || [];
    // Same rule as the inspect rail: a stored caveat is never shown as current
    // if the claim set it was written against has been superseded. It matters
    // MORE here — the two sides can sit on different claim sets, and a caveat
    // that appears under "A only" because B was run after a withdrawal is a
    // difference in the product's claims, not in the two runs.
    for (const [side, run] of [['A', A], ['B', B]]) {
      const n = claimSetNotice(run.provenance, { compact: true, label: side });
      if (n) cmpEl.appendChild(n);
    }
    const only = (x, other) => x.filter((c) => !other.includes(c));
    const onlyA = only(ca, cb), onlyB = only(cb, ca);
    if (!onlyA.length && !onlyB.length) {
      cmpEl.appendChild(h('div.caveat.info', h('span',
        ca.length ? 'Both runs carry the same caveats; the comparison inherits all of them.'
          : 'Neither run reported caveats.')));
      return;
    }
    for (const c of onlyA) cmpEl.appendChild(h('div.caveat', h('span', h('b', 'A only: '), c)));
    for (const c of onlyB) cmpEl.appendChild(h('div.caveat', h('span', h('b', 'B only: '), c)));
  }

  function renderCompare(force) {
    const A = runs.find((r) => r.run_id === aId) || null;
    const B = runs.find((r) => r.run_id === bId) || null;
    const sig = [aId, bId, A && A.status, B && B.status,
      A && A.runtime_s, B && B.runtime_s].join('~');
    if (!force && sig === lastCmpSig) return;
    lastCmpSig = sig;
    clear(cmpEl);

    if (!A || !B) {
      cmpEl.appendChild(h('div.notice', { style: { maxWidth: 'none', textAlign: 'left', padding: '4px 0' } },
        h('p.dim', 'Pick two runs with the ', h('b', 'A'), ' and ', h('b', 'B'),
          ' radio columns to see their metrics side by side with deltas.')));
      return;
    }
    if (A.run_id === B.run_id) {
      cmpEl.appendChild(h('div.caveat.info', h('span', 'A and B are the same run; every delta is zero by construction.')));
    }

    /* --- spec difference, first: a delta is meaningless without it ----- */
    const diffs = SPEC_ROWS
      .map(([k, f]) => {
        let va, vb;
        try { va = String(f(A.spec || {})); } catch { va = NULL_GLYPH; }
        try { vb = String(f(B.spec || {})); } catch { vb = NULL_GLYPH; }
        return { k, va, vb, same: va === vb };
      });
    const changed = diffs.filter((d) => !d.same);

    cmpEl.appendChild(h('div.hr-lbl', `Spec difference — ${changed.length} of ${diffs.length} fields differ`));
    const specTbl = h('table.tbl',
      h('thead', h('tr',
        h('th', 'field'),
        h('th', { title: A.run_id }, 'A · ' + A.run_id),
        h('th', { title: B.run_id }, 'B · ' + B.run_id))),
      h('tbody', (changed.length ? changed : diffs.slice(0, 3)).map((d) => h('tr',
        h('td', d.k),
        h('td.mono', { style: d.same ? null : { color: 'var(--fg-0)' } }, d.va),
        h('td.mono', { style: d.same ? null : { color: 'var(--accent)' } }, d.vb),
      ))),
    );
    // `.cmp-fields` — three columns of single short words in a 240 px dock
    // column. It WRAPS (layout.css); it must never ellipse, because `unifo`
    // and `unifor` are two different modes and a truncated mode name is a
    // wrong answer rather than a short one.
    cmpEl.appendChild(h('div.tblwrap.cmp-fields', specTbl));
    if (!changed.length) {
      cmpEl.appendChild(h('div.caveat.info', h('span',
        'The two specs are identical. Any metric difference is run-to-run nondeterminism, which for a seeded engine should be exactly zero — if it is not, that is a finding.')));
    }

    /* --- metric comparison --------------------------------------------- */
    /* CROSS-COMPOUND DELTAS ARE REFUSED, and the refusal is the panel's whole
     * position on the repositioning. A/B stays — it is the right tool for one
     * compound at two doses, on two meshes, under two routes, and none of that
     * is affected by anything measured this week. What is affected is Δ between
     * two DIFFERENT compounds: that is an ordering of two, and ordering
     * compounds by predicted effect is exactly the claim that was withdrawn on
     * 2026-08-21 (results/pk_ranking_robustness.md,
     * results/capacity_ranking_robustness.md).
     *
     * The refusal is deliberately WHOLE-TABLE rather than per-metric. There is
     * an argument that comparing TRANSPORT between two compounds — penetration
     * depth rather than kill — is defensible, since a decay length is a shape
     * quantity and survives both the plasma and the capacity failures. But it
     * still rests on the unscored descriptor path
     * (results/descriptor_falsification.md), and adjudicating metric-by-metric
     * would smuggle a design decision into a removal pass. Suppress the table,
     * state why, widen later on evidence. Reversible beats clever.
     *
     * Everything else on this panel survives untouched: the spec diff above,
     * the caveat union below, and both runs' own numbers, which are still on
     * screen in the run history table. Nothing is hidden — only the SUBTRACTION
     * between two compounds is declined. */
    const drugA = (A.spec && A.spec.drug) || null;
    const drugB = (B.spec && B.spec.drug) || null;
    const crossCompound = !!(drugA && drugB && drugA !== drugB);
    cmpEl.appendChild(h('div.hr-lbl', 'Metrics'));
    if (crossCompound) {
      cmpEl.appendChild(h('div.caveat.crit', h('span',
        h('b', 'NO CROSS-COMPOUND DELTA. '),
        `A is ${drugA} and B is ${drugB}, so the metric delta between them is `
        + 'an ordering of two compounds by predicted effect. That claim was '
        + 'measured and withdrawn: per-compound plasma exposure across the '
        + 'compounds this product ranked spans 113x, while the published '
        + 'ordering survives a spread of only 1.53x '
        + '(results/pk_ranking_robustness.md), and the breast top-1 pick '
        + 'changes under this repo’s own structural capacity rule '
        + '(results/capacity_ranking_robustness.md). Both runs’ numbers '
        + 'are unchanged and still on screen in the run history above; what is '
        + 'declined is the subtraction between them. Compare one compound '
        + 'against itself — two doses, two meshes, two routes — and '
        + 'this table returns, because a per-compound scale factor cancels out '
        + 'of a comparison that holds the compound fixed. The ',
        h('b', 'details'), ' button on either row re-runs that run’s own spec '
        + 'with one parameter moved, which is that comparison.')));
      cmpEl.appendChild(h('div.hr-lbl', 'Caveats that differ'));
      renderCaveatDiff(A, B);
      return;
    }
    const rows = Object.keys(METRICS).map((key) => {
      const m = METRICS[key];
      const va = A.metrics ? A.metrics[key] : null;
      const vb = B.metrics ? B.metrics[key] : null;
      const pa = metricProvenance(key, A);
      const pb = metricProvenance(key, B);
      const both = Number.isFinite(va) && Number.isFinite(vb);
      const d = both ? vb - va : null;
      const rp = both ? relPct(va, vb) : null;
      let dcls = 'flat';
      if (both && Math.abs(rp || 0) > 1e-9) {
        const improving = m.better === 'higher' ? d > 0 : d < 0;
        dcls = improving ? 'up' : 'down';
      }
      return h('tr',
        h('td', { title: m.src }, m.label),
        h('td.right.num', fmtMetric(key, va), h('span.unit', metricUnit(key))),
        h('td', chip(pa.kind, pa.why)),
        h('td.right.num', fmtMetric(key, vb), h('span.unit', metricUnit(key))),
        h('td', chip(pb.kind, pb.why)),
        h('td.right', both
          ? h('span.delta.' + dcls, fmtMetric(key, Math.abs(d)))
          : h('span.dim', NULL_GLYPH)),
        h('td.right', rp === null
          ? h('span.dim', NULL_GLYPH)
          : h('span.delta.' + dcls, pct(Math.abs(rp), 1) + '%')),
      );
    });
    const mTbl = h('table.tbl',
      h('thead', h('tr',
        h('th', 'metric'),
        h('th.right', 'A'), h('th', 'prov'),
        h('th.right', 'B'), h('th', 'prov'),
        h('th.right', 'Δ abs'), h('th.right', 'Δ rel'))),
      h('tbody', rows),
    );
    // `.cmp-metrics` — seven columns, four of them numeric and two of them
    // provenance chips. This one cannot wrap (a seven-column table folded into
    // 240 px sets its own headings one letter per line), so it keeps its
    // natural width and scrolls, with a permanently visible scrollbar.
    cmpEl.appendChild(h('div.tblwrap.cmp-metrics', mTbl));

    /* --- caveats union -------------------------------------------------- */
    cmpEl.appendChild(h('div.hr-lbl', 'Caveats that differ'));
    renderCaveatDiff(A, B);

    /* --- deltas are conditional ---------------------------------------- */
    const anyAssumed = Object.keys(METRICS).some((k) =>
      metricProvenance(k, A).kind === 'assumed' || metricProvenance(k, B).kind === 'assumed');
    if (anyAssumed) {
      cmpEl.appendChild(h('div.caveat', h('span',
        'At least one metric on each side is ASSUMED rather than computed from measured inputs. The delta between two upper bounds is a delta between two upper bounds — it is not a measured treatment effect.')));
    }
  }

  const root = h('div.dock',
    h('div.panel',
      h('div.panel-head',
        h('h2', 'Run history'),
        countEl,
        h('span.grow'),
        h('button.btn.ghost.sm', {
          onclick: () => onRefresh && onRefresh(),
          // Named for its object, not its verb: the reports panel two
          // columns to the right carried a button with the identical label
          // `refresh` hitting a different endpoint.
          title: 'Re-read the run list. GET /api/runs',
        }, 'refresh runs')),
      filterBar,
      tableWrap,
      detailEl),
    h('div.dock-div'),
    h('div.panel',
      h('div.panel-head', h('h2', 'Comparison'), h('span.grow'),
        h('span.hint', 'A vs B')),
      cmpEl),
  );

  render(true);
  renderCompare(true);

  return {
    el: root,
    setRuns(list) {
      runs = list || [];
      if (!aId && runs.length) aId = runs[Math.min(1, runs.length - 1)].run_id;
      if (!bId && runs.length) bId = runs[0].run_id;
      render();
      renderCompare();
    },
    setActive(id) { activeId = id; render(); },
    get activeId() { return activeId; },
  };
}
