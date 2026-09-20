// crosssection.js — ONE COMPOUND, MANY SECTIONS.
//
// Self-contained the way ui/projects.js is: it injects its own toolbar button,
// its own <dialog> and its own stylesheet, so it adds nothing to the three-rail
// workstation layout and owns no file anyone else edits. Every class it defines
// is prefixed `xs-`; everything else on screen is a shared component class that
// already exists (.btn, .tbl, .chip, .caveat, .hr-lbl, .hint, .mono).
//
// WHICH COMPARISON THIS IS, AND WHICH ONE IT IS NOT
// -------------------------------------------------
// Ordering two DIFFERENT COMPOUNDS by predicted effect was withdrawn on
// 2026-08-21: per-compound plasma exposure spans 113.1x across the compounds
// the GBM panel ranks while the published ordering first breaks at a spread of
// 1.53x (results/pk_ranking_robustness.md), and the breast top-1 pick changes
// under this repo's own structural capacity rule
// (results/capacity_ranking_robustness.md). No flag, no mode and no caveat
// restores it, and neither does an axis label: one grid on one colour scale IS
// a ranking however the compounds got there.
//
// The SAME compound on two pieces of tissue is a different question and a
// legitimate one, and it is the only comparison this view draws. The compound
// is FIXED; the SECTION varies.
//
// WHY THE ORDER ON SCREEN IS WHAT IT IS
// ------------------------------------
//   1. SECTIONS     — what each panel and PLATFORM can support, and the
//                     denominator every fraction below is a fraction OF. A
//                     number is never shown above its own capability verdict,
//                     and a spot-resolution section's per-cell statistics are
//                     refused here before they could be missed below.
//   2. REACH LADDER — one compound's F(x) on each section, on the ladder the
//                     RUNS reported. Where a section is refused, the refusal
//                     is rendered IN PLACE OF the number, with its reason —
//                     never a blank, never a dash, never a zero.
//   3. DISAGREEMENT — where the answer is NOT stable across tissue. This is
//                     the primary result: a compound whose curves CROSS
//                     between two patients' tissue cannot order those two
//                     tissues at all, and that is the most informative thing
//                     this engine currently says. It is placed above the
//                     absolute block for exactly that reason.
//   4. ABSOLUTE     — this compound's own level per section, and NEVER a
//                     headline. Withheld wholesale unless the backend says
//                     comparable: a number shown under a refusal is still a
//                     number shown.
//
// Every verdict, reason, blocker and caveat string comes from the backend
// verbatim (app/server/cross_section.py). This module decides layout and
// emphasis; it never decides comparability and it never composes a caveat.
//
// THE GUARD
// ---------
// `assertOneCompoundAxis` walks the payload recursively and throws if a second
// compound is named anywhere in it. render() calls it FIRST and paints the
// refusal instead of the axis when it throws, so this view cannot put two
// compounds on one axis even if a server hands it two. It is exported and pure
// (no DOM) so tests/test_cross_section.py can execute it under node.

import { h, clear, fill, $ } from '../util/dom.js';
import { auto, num, int, pct, sci, ts, NULL_GLYPH } from '../util/fmt.js';
import { api } from '../api.js';

/** Keys that name a compound, matched at any depth. Mirrors
 *  app/server/cross_section.py::_COMPOUND_KEYS — the same audit on both sides
 *  of the wire, because a guard that only runs on the server is a guard that a
 *  stale build or a hand-edited payload walks straight past. */
const COMPOUND_KEYS = ['compound', 'drug', 'compounds', 'drugs'];

export const WITHDRAWN_CROSS_COMPOUND_AXIS =
  'REFUSED: this view was handed more than one compound for one axis. Ordering '
  + 'two different compounds by predicted effect was WITHDRAWN on 2026-08-21, '
  + 'absolutely — no flag, no mode, no caveat and no axis label restores it, '
  + 'because one grid on one colour scale IS a ranking however the compounds '
  + 'got there (results/pk_ranking_robustness.md: per-compound plasma exposure '
  + 'spans 113.1x while the published ordering breaks at 1.53x; '
  + 'results/capacity_ranking_robustness.md: the breast top-1 pick changes '
  + 'under this repo’s own structural capacity rule). Nothing is rendered.';

function collectCompounds(node, found) {
  if (Array.isArray(node)) {
    for (const item of node) collectCompounds(item, found);
    return found;
  }
  if (node && typeof node === 'object') {
    for (const key of Object.keys(node)) {
      const value = node[key];
      if (COMPOUND_KEYS.includes(key)) {
        if (typeof value === 'string' && value.trim()) found.add(value.trim().toLowerCase());
        else if (Array.isArray(value)) {
          for (const item of value) {
            if (typeof item === 'string' && item.trim()) found.add(item.trim().toLowerCase());
          }
        }
      }
      collectCompounds(value, found);
    }
  }
  return found;
}

/**
 * The mechanical guard. Returns the ONE compound on this axis, or throws.
 *
 * Pure and DOM-free on purpose: it is the one piece of this module a test can
 * execute directly, and it is the piece that must never be skipped.
 */
export function assertOneCompoundAxis(payload) {
  const found = [...collectCompounds(payload, new Set())].sort();
  if (found.length === 1) return found[0];
  const err = new Error(`${WITHDRAWN_CROSS_COMPOUND_AXIS} Compounds found: `
    + `${found.length ? found.join(', ') : '(none)'}.`);
  err.name = 'CrossCompoundAxisRefused';
  err.compounds = found;
  throw err;
}

/* -------------------------------------------------------------- transport */

/** GET a JSON route with the same key header api.js uses.
 *
 *  Deliberately NOT falling back to js/mock.js the way api.js's `call()` does.
 *  There is no mock for this route, and a fabricated cross-section axis is the
 *  one failure this view must never produce: the whole output is "do these two
 *  real tissues agree", and seeded noise would answer it. When the backend is
 *  unreachable the panel says so and draws nothing. */
async function getJSON(path) {
  const headers = { accept: 'application/json' };
  if (api.apiKey) headers['x-api-key'] = api.apiKey;
  const res = await fetch(api.base + path, { headers, cache: 'no-store' });
  if (!res.ok) {
    let detail = `${res.status} ${res.statusText}`;
    try {
      const body = await res.json();
      const d = body && body.detail;
      if (typeof d === 'string') detail = d;
      else if (d && typeof d.error === 'string') detail = d.error;
    } catch { /* the status is still the story */ }
    const err = new Error(detail);
    err.status = res.status;
    throw err;
  }
  return res.json();
}

/* ------------------------------------------------------------ formatting */

const fracText = (v) => (typeof v === 'number' && Number.isFinite(v) ? pct(v, 3) : NULL_GLYPH);
const countText = (v) => (typeof v === 'number' && Number.isFinite(v) ? int(v) : NULL_GLYPH);

/** A refusal, in the place a number would have gone. Never a dash. */
function refusalCell(reason, short = 'REFUSED') {
  return h('td.xs-refusal', { title: reason || '' },
    h('span.chip.st-failed', short));
}

function caveat(text, kind) {
  return h(`div.caveat${kind ? '.' + kind : ''}`, { title: text }, h('span', text));
}

function labelled(label, text) {
  return h('div.caveat.crit', { title: text }, h('span', h('b', label + ': '), text));
}

/* ------------------------------------------------------------------ view */

export function createCrossSectionPanel({ store } = {}) {
  let compounds = [];          // the SET this server has run, alphabetical
  let compound = null;         // the ONE compound on the axis
  let payload = null;
  let signature = null;
  let metric = 'kill_fraction';
  let busy = false;
  let failure = null;

  /* ------------------------------------------------------------- styles */
  // Injected here rather than added to app/web/css, so this feature owns every
  // file it writes. Colour literals are forbidden by the workspace contract
  // §7: every value below is a token from css/tokens.css.
  if (!document.getElementById('xs-style')) {
    const css = `
#xs-dlg { width: min(1360px, 96vw); height: min(880px, 92vh); overflow: hidden; }
.xs-body { height: calc(100% - var(--h-head)); overflow: auto; min-height: 0;
           padding: var(--sp-4) var(--sp-5) var(--sp-7); }
.xs-axis { display: flex; gap: var(--sp-3); align-items: baseline;
           flex-wrap: wrap; margin-bottom: var(--sp-3); }
.xs-axis .xs-fixed { font-family: var(--font-mono); color: var(--accent);
                     font-size: var(--fs-md); }
.xs-verdict { border: 1px solid var(--line-1); border-left-width: 3px;
              border-radius: var(--radius); padding: var(--sp-3) var(--sp-4);
              margin: var(--sp-3) 0; background: var(--bg-2); }
.xs-verdict h3 { margin: 0 0 var(--sp-2); font-size: var(--fs-md); }
.xs-verdict p { margin: 0; color: var(--fg-1); line-height: var(--lh-body); }
.xs-verdict.xs-crit { border-left-color: var(--crit); background: var(--crit-bg); }
.xs-verdict.xs-ok { border-left-color: var(--ok); }
.xs-verdict.xs-info { border-left-color: var(--info); }
.xs-refusal { text-align: center; }
.xs-denom { color: var(--fg-2); font-size: var(--fs-micro); font-family: var(--font-mono); }
.xs-flip { color: var(--crit); font-weight: 700; }
.xs-hold { color: var(--ok); }
.xs-rung { font-family: var(--font-mono); }
.xs-conv { color: var(--fg-2); }
.xs-sub { color: var(--fg-2); font-size: var(--fs-xs); line-height: var(--lh-body); }
.xs-tbl { margin-bottom: var(--sp-3); }
`;
    document.head.appendChild(h('style#xs-style', css));
  }

  /* -------------------------------------------------------------- shell */

  const bodyEl = h('div.xs-body');
  const compoundSel = h('select.sel.sm', {
    title: 'The ONE compound held fixed on this axis. This list is a SET, '
      + 'alphabetical — never ordered by any metric, and it carries no "best" '
      + 'mark, because ordering compounds by predicted effect is withdrawn.',
    onchange: (e) => { compound = e.target.value || null; signature = null; reload(); },
  });
  const metricSel = h('select.sel.sm', {
    title: 'Which ABSOLUTE column is shown below the ladder. It selects a '
      + 'column, not a conclusion, and cannot change the ladder or the verdict.',
    onchange: (e) => { metric = e.target.value; reload(); },
  }, ...['kill_fraction', 'mean_tumor_auc', 'therapeutic_coverage_1x',
    'penetration_length_um'].map((m) => h('option', { value: m }, m)));

  const dlg = h('dialog#xs-dlg',
    h('div.panel-head',
      h('h2', 'One compound, many sections'),
      h('span.hint', 'where the answer holds across tissue, and where it flips'),
      h('span.vrule'),
      h('span.lbl', 'compound'), compoundSel,
      h('span.grow'),
      h('span.lbl', 'absolute column'), metricSel,
      h('button.btn.ghost.sm', {
        onclick: () => reload(), title: 'GET /api/cross-section',
      }, 'refresh'),
      h('button.btn.ghost.sm', { id: 'xs-close', onclick: () => dlg.close() }, 'close')),
    bodyEl);
  document.body.appendChild(dlg);

  const btn = h('button.btn.ghost.sm', {
    id: 'btn-crosssection',
    title: 'One compound across every section that ran it: where its reach '
      + 'verdict holds across tissue and where it flips',
    onclick: open,
  }, 'sections');
  const slot = $('.topbar-group.last');
  if (slot) slot.insertBefore(btn, slot.firstChild);

  /* ---------------------------------------------------------- data flow */

  const offline = () => !!(api && api.stub);

  async function open() {
    if (!dlg.open) dlg.showModal();
    await loadCompounds();
    await reload();
  }

  async function loadCompounds() {
    if (offline()) return;
    try {
      const out = await getJSON('/api/cross-section/compounds');
      compounds = out.compounds || [];
      if (!compound) {
        const first = compounds.find((c) => c.comparable_sections);
        compound = first ? first.compound : (compounds[0] || {}).compound || null;
      }
    } catch (err) {
      failure = `GET /api/cross-section/compounds failed: ${err.message}`;
      compounds = [];
    }
    fill(compoundSel, ...compounds.map((c) => h('option', {
      value: c.compound,
      selected: c.compound === compound ? true : null,
      title: c.comparable_sections ? `${c.n_sections} sections` : (c.reason || ''),
    }, `${c.compound} — ${c.n_sections} section${c.n_sections === 1 ? '' : 's'}`)));
  }

  async function reload() {
    failure = null;
    if (offline()) { payload = null; render(); return; }
    if (!compound) { payload = null; render(); return; }
    busy = true; render();
    try {
      const q = new URLSearchParams({ compound, metric });
      if (signature) q.set('signature', signature);
      payload = await getJSON(`/api/cross-section?${q.toString()}`);
    } catch (err) {
      payload = null;
      failure = `GET /api/cross-section failed: ${err.message}`;
    } finally {
      busy = false;
      render();
    }
  }

  /* ------------------------------------------------------------ render */

  function render() {
    clear(bodyEl);

    if (offline()) {
      bodyEl.appendChild(labelled('NO SERVER',
        'The backend is unreachable and this view has no stub. It will not '
        + 'draw a cross-section axis out of synthetic data: the entire output '
        + 'is whether two REAL tissues agree, and seeded noise would answer '
        + 'that question with a fabrication.'));
      return;
    }
    if (failure) { bodyEl.appendChild(labelled('REQUEST FAILED', failure)); return; }
    if (busy && !payload) { bodyEl.appendChild(h('div.hint', 'loading…')); return; }
    if (!payload) {
      bodyEl.appendChild(caveat(
        'No compound has been run on a section on this server yet, so there is '
        + 'nothing to hold fixed. Run one compound on two sections at the same '
        + 'RunSpec and this view becomes a comparison.', 'info'));
      return;
    }

    // THE GUARD, before a single row is built. A payload naming two compounds
    // is not rendered in any reduced form: the refusal replaces the view.
    let axisCompound;
    try {
      axisCompound = assertOneCompoundAxis(payload);
    } catch (err) {
      bodyEl.appendChild(labelled('AXIS REFUSED', err.message));
      return;
    }

    const sections = payload.sections || [];
    const ids = sections.map((s) => s.section);
    const labelOf = (id) => (sections.find((s) => s.section === id) || {}).label || id;

    /* ---- axis statement ------------------------------------------------ */
    bodyEl.appendChild(h('div.xs-axis',
      h('span.lbl', 'compound FIXED'),
      h('span.xs-fixed', axisCompound),
      h('span.lbl', 'sections VARY'),
      h('span.xs-fixed', `${sections.length}`),
      h('span.lbl', 'signature'),
      h('span.mono.xs-sub', { title: payload.signature_note || '' },
        payload.signature || NULL_GLYPH)));
    bodyEl.appendChild(h('div.xs-sub', (payload.axis || {}).statement || ''));
    bodyEl.appendChild(caveat(payload.signature_note || '', 'info'));

    /* ---- 1. sections --------------------------------------------------- */
    bodyEl.appendChild(h('div.hr-lbl',
      '1 · Sections — what each panel and platform can support'));
    bodyEl.appendChild(h('div.tblwrap.xs-tbl', h('table.tbl.zebra',
      h('thead', h('tr',
        h('th', 'section'), h('th', 'panel'), h('th', 'platform'),
        h('th.right', 'viable tumour cells'), h('th.right', 'of tumour cells'),
        h('th', 'reach map'), h('th', 'per-cell statistics'), h('th', 'run'))),
      h('tbody', sections.map((s) => h('tr',
        h('td', h('b', s.label || s.section), h('div.xs-denom', s.section)),
        h('td.mono.xs-sub', s.panel || NULL_GLYPH),
        h('td.mono.xs-sub', (s.spot || {}).platform || 'unknown'),
        h('td.right.num', countText(s.n_viable)),
        h('td.right.num', countText(s.n_tumor)),
        h('td.mono.xs-sub', s.window_mode || NULL_GLYPH),
        h('td', (s.spot || {}).per_cell_statistics_refused
          ? h('span.chip.st-failed', { title: (s.spot || {}).refusal || '' }, 'REFUSED')
          : (s.spot || {}).spatial_spread_valid === null
            ? h('span.chip.unknown', {
              title: 'This section’s platform could not be named, so '
                + 'whether per-cell statistics are valid on it was NOT '
                + 'established. That is a different claim from "valid".',
            }, 'not established')
            : h('span.chip.measured', 'permitted')),
        h('td.mono.xs-sub', { title: `${s.run_id} · ${ts(s.created_at)}` },
          s.run_id || NULL_GLYPH)))))));

    for (const s of sections) {
      for (const lab of (s.labels || [])) {
        bodyEl.appendChild(labelled(s.label || s.section, lab));
      }
      if (!s.absolute_kill_supported && s.absolute_kill_reason) {
        bodyEl.appendChild(caveat(`${s.label || s.section}: ${s.absolute_kill_reason}`));
      }
      if ((s.spot || {}).per_cell_statistics_refused) {
        bodyEl.appendChild(labelled(
          `${s.label || s.section} — PER-CELL STATISTICS REFUSED`,
          (s.spot || {}).refusal || ''));
        for (const q of ((s.spot || {}).refused_quantities || [])) {
          bodyEl.appendChild(caveat(`withheld on this section — ${q.quantity}: ${q.why}`));
        }
      }
    }

    /* ---- 2. the ladder ------------------------------------------------- */
    const reach = payload.reach || {};
    bodyEl.appendChild(h('div.hr-lbl',
      `2 · ${axisCompound} on ${ids.length} section${ids.length === 1 ? '' : 's'} — `
      + 'under-reached fraction F(x), the ladder the runs reported'));
    bodyEl.appendChild(h('div.xs-sub', reach.ladder_is || ''));

    const refusedBy = {};
    for (const r of (reach.sections_refused || [])) refusedBy[r.section] = r.reason;

    if (!(reach.rows || []).length) {
      bodyEl.appendChild(caveat(
        'No section on this axis can supply an F(x) ladder, so there is no '
        + 'curve to draw. The blockers below are the output.', 'crit'));
    } else {
      bodyEl.appendChild(h('div.tblwrap.xs-tbl', h('table.tbl.zebra',
        h('thead', h('tr',
          h('th.right', 'x (reach ratio)'),
          ...ids.map((id) => h('th.right', { title: labelOf(id) }, id)))),
        h('tbody', (reach.rows || []).map((row) => h('tr',
          h('td.right.xs-rung', num(row.x, 2),
            row.is_convention_point
              ? h('span.xs-conv', {
                title: 'x = 0.1 is a REPORTING CONVENTION, not a physical '
                  + 'constant, and it is NOT the headline here. The curve is '
                  + 'the evidence; read it, not one point on it.',
              }, ' · convention')
              : null),
          ...ids.map((id) => (refusedBy[id]
            ? refusalCell(refusedBy[id])
            : h('td.right.num',
              fracText((row.under_reached_fraction || {})[id]),
              h('div.xs-denom',
                `${countText((row.under_reached_cells || {})[id])} of `
                + `${countText((row.of_viable_tumor_cells || {})[id])} viable`))))))))));
      bodyEl.appendChild(h('div.xs-sub',
        ((reach.rows || [])[0] || {}).cells_are_derived || ''));
    }

    for (const r of (reach.sections_refused || [])) {
      bodyEl.appendChild(labelled(
        `${r.label || r.section} — NO NUMBER ON THIS LADDER`, r.reason || ''));
    }
    for (const b of (reach.blockers || [])) {
      bodyEl.appendChild(labelled('NOT COMPARABLE', b));
    }

    /* ---- 3. the disagreement — the primary result ---------------------- */
    const stab = payload.stability || {};
    bodyEl.appendChild(h('div.hr-lbl',
      '3 · Where the answer is NOT stable across tissue'));
    bodyEl.appendChild(h(`div.xs-verdict.xs-${stab.level || 'info'}`,
      h('h3', stab.verdict || NULL_GLYPH),
      h('p', stab.text || '')));

    if ((stab.pairs || []).length) {
      bodyEl.appendChild(h('div.tblwrap.xs-tbl', h('table.tbl.zebra',
        h('thead', h('tr',
          h('th', 'section pair'), h('th.right', 'rungs compared'),
          h('th', 'ordering'), h('th', 'higher throughout'),
          h('th', 'fraction vs count'))),
        h('tbody', stab.pairs.map((p) => h('tr',
          h('td.mono.xs-sub', (p.sections || []).join('  vs  ')),
          h('td.right.num', int(p.n_rungs_compared)),
          h('td', h('span', { class: p.order_flips ? 'xs-flip' : 'xs-hold' },
            p.verdict || '')),
          h('td.mono.xs-sub', p.order_flips
            ? h('span.xs-flip', 'no ordering exists')
            : (p.higher_throughout || 'identical')),
          h('td', (p.fraction_count_divergences || []).length
            ? h('span.chip.st-failed', {
              title: 'The FRACTION and the COUNT disagree about which section '
                + 'is more under-reached, at the rungs listed below. '
                + 'results/heldout_patient.md §1.1 is the archetype: a '
                + 'fraction fell 0.67x while the areal density rose 4.19x, '
                + 'because the denominator changed.',
            }, `${p.fraction_count_divergences.length} rung(s) DISAGREE`)
            : h('span.chip.measured', 'agree'))))))));

      for (const p of stab.pairs) {
        if (p.order_flips) {
          bodyEl.appendChild(labelled(
            `ORDER FLIPS — ${(p.sections || []).join(' vs ')}`,
            'Which of these two tissues is more under-reached by '
            + `${axisCompound} depends on which rung of the ladder is read `
            + `(sign changes across x = ${(p.flip_rungs || []).map((x) => num(x, 2)).join(', ')}). `
            + 'The curves cross, so no ordering of these two tissues exists '
            + 'for this compound. This is a finding, not a failure.'));
        }
        for (const d of (p.fraction_count_divergences || [])) {
          const [a, b] = p.sections || [];
          bodyEl.appendChild(labelled(
            `FRACTION AND COUNT POINT IN OPPOSITE DIRECTIONS at x = ${num(d.x, 2)}`,
            `${a} ${fracText((d.fractions || {})[a])} of `
            + `${countText((row_denom(reach, d.x, a)))} viable vs `
            + `${b} ${fracText((d.fractions || {})[b])} of `
            + `${countText((row_denom(reach, d.x, b)))} viable — yet the CELL `
            + `COUNTS are ${countText((d.counts || {})[a])} and `
            + `${countText((d.counts || {})[b])}. The fraction says one tissue `
            + 'is worse and the count says the other is. Neither is wrong: the '
            + 'denominators differ. results/heldout_patient.md §1.1 measured '
            + 'exactly this shape on two real GBM sections and its conclusion '
            + 'was "these two rows point in opposite directions, and that is '
            + 'the finding".'));
        }
      }
    }
    if (stab.how_read) bodyEl.appendChild(caveat(stab.how_read, 'info'));
    if (stab.noise_floor) bodyEl.appendChild(caveat(stab.noise_floor, 'crit'));

    /* ---- 4. the absolute column ---------------------------------------- */
    const abs = payload.absolute || {};
    bodyEl.appendChild(h('div.hr-lbl',
      `4 · ${abs.metric || metric} — this compound's own level per section. `
      + 'NEVER a headline'));
    if (abs.withheld) {
      bodyEl.appendChild(labelled('WITHHELD',
        abs.withheld_reason
        || 'the backend verdict on this pair is comparable=false, so no value '
           + 'is rendered.'));
    } else {
      bodyEl.appendChild(h('div.tblwrap.xs-tbl', h('table.tbl',
        h('thead', h('tr', h('th', 'section'), h('th.right', abs.metric || metric),
          h('th', 'cells.state'))),
        h('tbody', ids.map((id) => h('tr',
          h('td', labelOf(id), h('div.xs-denom', id)),
          h('td.right.num', fmtMetric(abs.metric || metric, (abs.values || {})[id])),
          h('td', (abs.state_measured || {})[id]
            ? h('span.chip.measured', 'MEASURED')
            : h('span.chip.assumed', 'ASSUMED'))))))));
    }
    for (const b of (abs.blockers || [])) bodyEl.appendChild(labelled('BLOCKER', b));
    for (const c of (abs.caveats || [])) bodyEl.appendChild(caveat(c));

    /* ---- caveats that travel with the ladder --------------------------- */
    bodyEl.appendChild(h('div.hr-lbl', 'Caveats — each travels with the numbers above'));
    for (const c of (reach.caveats || [])) bodyEl.appendChild(caveat(c));
    for (const p of (payload.problems || [])) bodyEl.appendChild(caveat(p, 'crit'));
    bodyEl.appendChild(caveat(payload.withdrawn || '', 'crit'));

    /* ---- other signatures, last: provenance for the whole view --------- */
    if ((payload.matched_sets || []).length > 1) {
      bodyEl.appendChild(h('div.hr-lbl',
        'Other RunSpec signatures this compound was run at — not merged in'));
      bodyEl.appendChild(h('div.tblwrap', h('table.tbl',
        h('thead', h('tr', h('th', 'signature'), h('th.right', 'sections'),
          h('th', 'on this axis'))),
        h('tbody', payload.matched_sets.map((s) => h('tr',
          h('td.mono.xs-sub', s.signature),
          h('td.right.num', int(s.n_sections)),
          h('td', s.is_primary
            ? h('span.chip.measured', 'yes')
            : h('button.btn.ghost.sm', {
              onclick: () => { signature = s.signature; reload(); },
            }, 'switch to it'))))))));
    }
  }

  /** The denominator a rung reported for one section — looked up rather than
   *  recomputed, so the count under a fraction and the count in this sentence
   *  are the same number. */
  function row_denom(reach, x, section) {
    const row = (reach.rows || []).find((r) => r.x === x) || {};
    return (row.of_viable_tumor_cells || {})[section];
  }

  function fmtMetric(name, v) {
    if (typeof v !== 'number' || !Number.isFinite(v)) return NULL_GLYPH;
    if (name === 'kill_fraction' || name === 'therapeutic_coverage_1x') return pct(v, 2);
    if (name === 'mean_tumor_auc') return sci(v, 3);
    if (name === 'penetration_length_um') return num(v, 1) + ' um';
    return auto(v);
  }

  return {
    el: dlg,
    button: btn,
    open,
    reload,
    get compound() { return compound; },
    get payload() { return payload; },
  };
}
