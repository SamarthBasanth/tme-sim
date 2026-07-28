// history.js — the bottom dock: run history and A/B comparison.
//
// Two explicit radio columns pick the comparison pair. Selecting a run for
// comparison is deliberately NOT the same gesture as opening it in the
// viewer — conflating "look at this" with "compare against this" is how
// people end up reading the wrong pair of numbers.

import { h, clear } from '../util/dom.js';
import { METRICS, metricProvenance, chip } from './provenance.js';
import { auto, num, int, pct, sci, dur, ts, NULL_GLYPH, relPct } from '../util/fmt.js';

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
  ['mgmt', (s) => (s.mgmt === null || s.mgmt === undefined ? 'null' : num(s.mgmt, 2) + ' log2CPM')],
  ['pk.plasma_peak', (s) => num(s.pk && s.pk.plasma_peak, 2) + ' uM'],
  ['pk.half_life_s', (s) => int(s.pk && s.pk.half_life_s) + ' s'],
  ['pk.n_doses', (s) => int(s.pk && s.pk.n_doses)],
  ['pk.dose_interval_s', (s) => int(s.pk && s.pk.dose_interval_s) + ' s'],
  ['total_time_s', (s) => int(s.total_time_s) + ' s'],
];

export function createHistoryPanel({ onOpen, onRefresh }) {
  let runs = [];
  let activeId = null;
  let aId = null;
  let bId = null;

  const tbody = h('tbody');
  const countEl = h('span.hint', '0 runs');
  const cmpEl = h('div.panel-body.pad');

  const table = h('table.tbl.zebra',
    h('thead', h('tr',
      h('th', { style: { width: '22px' }, title: 'comparison slot A' }, 'A'),
      h('th', { style: { width: '22px' }, title: 'comparison slot B' }, 'B'),
      h('th', 'run'),
      h('th', 'status'),
      h('th', 'section'),
      h('th', 'drug'),
      h('th.right', 'dx / um'),
      h('th.right', 'kill / %'),
      h('th.right', 'AUC / uM.s'),
      h('th.right', 'pen / um'),
      h('th.right', 'runtime'),
      h('th', 'submitted'),
    )),
    tbody,
  );

  function render() {
    clear(tbody);
    countEl.textContent = `${runs.length} run${runs.length === 1 ? '' : 's'}`;
    if (!runs.length) {
      tbody.appendChild(h('tr', h('td', { colspan: '12', class: 'dim' },
        'No runs yet. Submit a RunSpec from the left rail.')));
      return;
    }
    for (const r of runs) {
      const m = r.metrics || {};
      const tr = h('tr.sel-able', {
        tabindex: '0',
        'aria-selected': String(r.run_id === activeId),
        onclick: (e) => { if (e.target.tagName !== 'INPUT') open(r.run_id); },
        onkeydown: (e) => {
          if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(r.run_id); }
          if (e.key === 'ArrowDown' && tr.nextElementSibling) { e.preventDefault(); tr.nextElementSibling.focus(); }
          if (e.key === 'ArrowUp' && tr.previousElementSibling) { e.preventDefault(); tr.previousElementSibling.focus(); }
        },
      },
        h('td', h('input', {
          type: 'radio', name: 'cmp-a', checked: aId === r.run_id,
          'aria-label': `set ${r.run_id} as comparison A`,
          onchange: () => { aId = r.run_id; renderCompare(); },
        })),
        h('td', h('input', {
          type: 'radio', name: 'cmp-b', checked: bId === r.run_id,
          'aria-label': `set ${r.run_id} as comparison B`,
          onchange: () => { bId = r.run_id; renderCompare(); },
        })),
        h('td.mono', r.run_id),
        h('td', h('span.chip.st-' + (r.status || 'queued'), r.status || '?')),
        h('td', { title: r.spec && r.spec.section }, (r.spec && r.spec.section) || NULL_GLYPH),
        h('td', (r.spec && (r.spec.drug ?? 'none')) || NULL_GLYPH),
        h('td.right.num', num(r.spec && r.spec.dx_um, 1)),
        h('td.right.num', fmtMetric('kill_fraction', m.kill_fraction)),
        h('td.right.num', fmtMetric('mean_tumor_auc', m.mean_tumor_auc)),
        h('td.right.num', fmtMetric('penetration_length_um', m.penetration_length_um)),
        h('td.right.num', dur(r.runtime_s)),
        h('td.mono', { style: { fontSize: 'var(--fs-micro)' } }, ts(r.created_at)),
      );
      tbody.appendChild(tr);
    }
  }

  function open(id) {
    activeId = id;
    render();
    onOpen && onOpen(id);
  }

  function renderCompare() {
    const A = runs.find((r) => r.run_id === aId) || null;
    const B = runs.find((r) => r.run_id === bId) || null;
    clear(cmpEl);

    if (!A || !B) {
      cmpEl.appendChild(h('div.empty', { style: { maxWidth: 'none', textAlign: 'left', padding: '4px 0' } },
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
    cmpEl.appendChild(h('div.tblwrap', specTbl));
    if (!changed.length) {
      cmpEl.appendChild(h('div.caveat.info', h('span',
        'The two specs are identical. Any metric difference is run-to-run nondeterminism, which for a seeded engine should be exactly zero — if it is not, that is a finding.')));
    }

    /* --- metric comparison --------------------------------------------- */
    cmpEl.appendChild(h('div.hr-lbl', 'Metrics'));
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
    cmpEl.appendChild(h('div.tblwrap', mTbl));

    /* --- caveats union -------------------------------------------------- */
    const ca = (A.provenance && A.provenance.caveats) || [];
    const cb = (B.provenance && B.provenance.caveats) || [];
    const only = (x, other) => x.filter((c) => !other.includes(c));
    const onlyA = only(ca, cb), onlyB = only(cb, ca);
    cmpEl.appendChild(h('div.hr-lbl', 'Caveats that differ'));
    if (!onlyA.length && !onlyB.length) {
      cmpEl.appendChild(h('div.caveat.info', h('span',
        ca.length ? 'Both runs carry the same caveats; the comparison inherits all of them.'
          : 'Neither run reported caveats.')));
    } else {
      for (const c of onlyA) cmpEl.appendChild(h('div.caveat', h('span', h('b', 'A only: '), c)));
      for (const c of onlyB) cmpEl.appendChild(h('div.caveat', h('span', h('b', 'B only: '), c)));
    }

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
        h('button.btn.ghost.sm', { onclick: () => onRefresh && onRefresh(), title: 'GET /api/runs' }, 'refresh')),
      h('div.panel-body.tblwrap', table)),
    h('div.dock-div'),
    h('div.panel',
      h('div.panel-head', h('h2', 'Comparison'), h('span.grow'),
        h('span.hint', 'A vs B')),
      cmpEl),
  );

  renderCompare();

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
