// reports.js — the therapy-report panel.
//
// A REPORT is the deliverable: a ranked, mechanistically attributed drug panel
// for ONE specimen, with the provenance of every input and the limits of the
// output attached. This panel does four things and nothing else:
//
//   1. trigger a report for the section currently in the RunSpec editor,
//   2. watch it build (it is 2 runs per drug on the shared job queue),
//   3. view the rendered document inline,
//   4. download that document as one self-contained HTML file.
//
// It is deliberately ADDITIVE: it appends itself as a third dock column and
// brings its own stylesheet (css/reports.css). It does not modify the viewer,
// the run history, or the spec editor.
//
// Routes consumed:  POST /api/reports · GET /api/reports · GET /api/reports/{id}
//                   GET /api/reports/{id}/html · DELETE /api/reports/{id}

import { h, clear } from '../util/dom.js';
import * as api from '../api.js';
import { toast, notifyError } from './toast.js';
import { num, pct, dur, NULL_GLYPH } from '../util/fmt.js';

/** The five-drug GBM panel experiments/patient_report.py ranks. Pre-selected. */
export const DEFAULT_PANEL = [
  'temozolomide', 'lomustine', 'cisplatin', 'methotrexate', 'doxorubicin',
];

const STATUS_CHIP = {
  queued: 'st-queued', running: 'st-running', done: 'st-done',
  failed: 'st-failed', cancelled: 'st-cancelled',
};

export function createReportsPanel({ store, getSection }) {
  let reports = [];
  let selected = null;
  let picked = new Set(DEFAULT_PANEL);
  let submitting = false;

  /* ---------------------------------------------------- drug picker ---- */
  const drugBox = h('div.rp-drugs');

  function buildDrugs(drugs) {
    clear(drugBox);
    const names = (drugs || []).map((d) => d.name).filter(Boolean);
    const list = names.length ? names : DEFAULT_PANEL;
    for (const n of list) {
      const cb = h('input', {
        type: 'checkbox', checked: picked.has(n),
        'aria-label': `include ${n} in the panel`,
        onchange: (e) => {
          if (e.target.checked) picked.add(n); else picked.delete(n);
          paintHead();
        },
      });
      drugBox.appendChild(h('label.chk.rp-drug', { title: n }, cb, h('span', n)));
    }
    paintHead();
  }

  const countEl = h('span.hint');
  function paintHead() {
    const sec = getSection() || NULL_GLYPH;
    countEl.textContent = `${picked.size} drug${picked.size === 1 ? '' : 's'} · ${sec}`;
    btnNew.disabled = submitting || picked.size === 0;
  }

  /* ------------------------------------------------------- submission -- */
  const btnNew = h('button.btn.primary.sm', {
    title: 'POST /api/reports — enqueues 2 runs per drug (the real arm and the resistance-off counterfactual) on the same reconstruction',
    onclick: submit,
  }, 'generate report');

  async function submit() {
    const section = getSection();
    if (!section) {
      toast('warn', 'No section selected',
        'A report is a property of one specimen. Pick a section in the RunSpec editor first.');
      return;
    }
    const drugs = Array.from(picked);
    if (!drugs.length) return;
    submitting = true; paintHead();
    const spec = (store && store.get().spec) || {};
    // Carry the operator's own grid / window / PK across; the backend refuses
    // and reports back any field a report owns (drug, sensitivity_mode).
    const spec_overrides = {};
    for (const k of ['dx_um', 'total_time_s', 'vascular_mode', 'oxygen_mode', 'state_mode', 'mgmt', 'pk']) {
      if (spec[k] !== undefined) spec_overrides[k] = spec[k];
    }
    try {
      const res = await api.postReport({ section, drugs, spec_overrides });
      toast('ok', `Report ${res.report_id} queued`,
        `${(res.run_ids || []).length} runs: the real arm and the resistance-off counterfactual for each of ${drugs.length} drugs. The reconstruction is built once and shared.`,
        6000);
      for (const w of (res.warnings || [])) toast('warn', 'Report spec warning', w, 0);
      selected = res.report_id;
      await refresh();
    } catch (err) {
      notifyError('POST /api/reports failed', err);
    } finally {
      submitting = false; paintHead();
    }
  }

  /* ------------------------------------------------------------ list --- */
  const tbody = h('tbody');
  const table = h('table.tbl',
    h('thead', h('tr',
      h('th', 'report'),
      h('th', 'status'),
      h('th', 'section'),
      h('th.right', 'arms'),
      h('th', 'top-ranked'),
      h('th.right', 'kill'),
      h('th.right', 'runtime'),
      h('th', ''))),
    tbody);

  function renderList() {
    clear(tbody);
    if (!reports.length) {
      tbody.appendChild(h('tr', h('td', { colspan: '8', class: 'dim' },
        'No reports yet. Pick a drug panel and press ',
        h('b', 'generate report'), '.')));
      return;
    }
    for (const r of reports) {
      const p = r.progress || {};
      const top = r.top_ranked;
      const isDone = r.status === 'done';
      const tr = h('tr.sel-able', {
        tabindex: '0',
        'aria-selected': String(r.report_id === selected),
        onclick: (e) => { if (e.target.tagName !== 'BUTTON') open(r.report_id); },
        onkeydown: (e) => {
          if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(r.report_id); }
        },
      },
        h('td.mono', r.report_id),
        h('td',
          h('span.chip.' + (STATUS_CHIP[r.status] || 'st-queued'), r.status),
          r.status === 'running'
            ? h('span.rp-prog', { title: `${p.runs_done || 0} of ${p.runs_total || 0} runs done${p.stage ? ' · ' + p.stage : ''}` },
              h('i', { style: { width: pct(p.fraction || 0, 0) + '%' } }))
            : null),
        h('td', { title: r.section }, r.section || NULL_GLYPH),
        h('td.right.num', `${p.runs_done ?? 0}/${r.n_runs ?? 0}`),
        h('td', top ? h('b', top.drug) : h('span.dim', NULL_GLYPH)),
        h('td.right.num', top && top.kill !== null && top.kill !== undefined
          ? num(top.kill, 4) : NULL_GLYPH),
        h('td.right.num', dur(r.runtime_s)),
        h('td.right',
          isDone
            ? h('button.btn.ghost.sm', {
              title: 'Download the self-contained HTML document',
              onclick: (e) => { e.stopPropagation(); download(r.report_id); },
            }, 'download')
            : (r.status === 'running' || r.status === 'queued')
              ? h('button.btn.ghost.sm', {
                title: 'DELETE /api/reports/{id} — cancels every queued arm',
                onclick: (e) => { e.stopPropagation(); cancel(r.report_id); },
              }, 'cancel')
              : null),
      );
      if (r.status === 'failed' && r.error) tr.title = r.error;
      tbody.appendChild(tr);
    }
  }

  /* --------------------------------------------------------- viewer ---- */
  const frame = h('iframe.rp-frame', {
    title: 'Rendered therapy report',
    // The document is served from this same origin and loads nothing external;
    // sandboxing it anyway means a report can never script the workstation.
    sandbox: '',
  });
  const dlgHead = h('div.panel-head');
  const dlgNote = h('div.rp-note');
  const dlg = h('dialog#rp-dlg', dlgHead, dlgNote, frame);

  function paintViewer(r) {
    clear(dlgHead);
    dlgHead.appendChild(h('h2', 'Therapy report'));
    dlgHead.appendChild(h('span.mono.strong', r.report_id));
    dlgHead.appendChild(h('span.chip.' + (STATUS_CHIP[r.status] || 'st-queued'), r.status));
    dlgHead.appendChild(h('span.chip.crit', 'not for clinical use'));
    dlgHead.appendChild(h('span.grow'));
    dlgHead.appendChild(h('button.btn.ghost.sm', {
      onclick: () => download(r.report_id),
      title: 'GET /api/reports/{id}/html — one self-contained file, no external requests',
    }, 'download html'));
    dlgHead.appendChild(h('button.btn.ghost.sm', { onclick: () => dlg.close() }, 'close'));

    clear(dlgNote);
    if (r.status === 'done') {
      frame.src = api.reportHtmlUrl(r.report_id);
      frame.style.display = '';
    } else {
      frame.removeAttribute('src');
      frame.style.display = 'none';
      const p = r.progress || {};
      dlgNote.appendChild(h('div.caveat' + (r.status === 'failed' ? '.crit' : '.info'),
        h('span', r.status === 'failed'
          ? `This report failed: ${r.error || 'no message returned'}`
          : `Building — ${p.runs_done || 0} of ${p.runs_total || 0} arms complete${p.stage ? ` (${p.stage})` : ''}. The document renders when every arm is done; a partial panel is not a ranking.`)));
    }
    // The limits travel with the summary too, so they are visible before the
    // document has finished rendering.
    for (const lim of (r.limits || [])) {
      if (!lim.applies || lim.severity === 'info') continue;
      dlgNote.appendChild(h('div.caveat' + (lim.severity === 'critical' ? '.crit' : ''),
        h('span', h('b', lim.title + ' '), lim.body)));
    }
  }

  async function open(id) {
    selected = id;
    renderList();
    let r;
    try {
      r = await api.getReport(id);
    } catch (err) {
      notifyError(`GET /api/reports/${id} failed`, err);
      return;
    }
    paintViewer(r);
    if (!dlg.open) dlg.showModal();
  }

  async function download(id) {
    try {
      const res = await fetch(api.reportHtmlUrl(id), { cache: 'no-store' });
      if (!res.ok) throw new Error(`${res.status} ${res.statusText}`);
      const blob = new Blob([await res.text()], { type: 'text/html' });
      const url = URL.createObjectURL(blob);
      const a = h('a', { href: url, download: `tmesim_therapy_report_${id}.html` });
      document.body.appendChild(a);
      a.click();
      a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 4000);
    } catch (err) {
      notifyError(`GET /api/reports/${id}/html failed`, err);
    }
  }

  async function cancel(id) {
    try {
      await api.cancelReport(id);
      await refresh();
    } catch (err) {
      notifyError(`DELETE /api/reports/${id} failed`, err);
    }
  }

  /* ----------------------------------------------------------- poll ---- */
  const errEl = h('div.rp-err');

  async function refresh() {
    try {
      reports = await api.listReports();
      clear(errEl);
      renderList();
      if (dlg.open && selected) {
        const cur = reports.find((x) => x.report_id === selected);
        if (cur && cur.status !== 'done') {
          try { paintViewer(await api.getReport(selected)); } catch { /* keep the last view */ }
        }
      }
    } catch (err) {
      clear(errEl);
      errEl.appendChild(h('div.caveat', h('span',
        err && err.status === 404
          ? 'Therapy reports need the live backend — the stub in js/mock.js does not synthesise them, and a synthesised report would be a fabricated clinical document. Start the server with ./.venv/bin/python -m app.server.'
          : `GET /api/reports failed: ${(err && err.message) || err}`)));
    }
  }

  /** True while any report still has arms in flight — drives the poll rate. */
  function busy() {
    return reports.some((r) => r.status === 'queued' || r.status === 'running');
  }

  const root = h('div.panel.rp',
    h('div.panel-head',
      h('h2', 'Therapy reports'),
      countEl,
      h('span.grow'),
      h('button.btn.ghost.sm', { onclick: refresh, title: 'GET /api/reports' }, 'refresh'),
      btnNew),
    h('div.panel-body.rp-body',
      h('div.rp-lead',
        h('b', 'One specimen in, one ranked drug panel out.'),
        ' Each drug is simulated twice on the same cached reconstruction — the real run, and a resistance-off counterfactual — so the O',
        h('sub', '2'),
        '-resistance term of the attribution is measured, not asserted. The document carries its own limits: ',
        h('b', 'the ranking is the output; the absolute kill is not.')),
      drugBox,
      errEl,
      h('div.tblwrap', table)),
  );

  document.body.appendChild(dlg);
  renderList();
  paintHead();

  return {
    el: root,
    dialog: dlg,
    setDrugs: buildDrugs,
    refresh,
    busy,
    setSection: paintHead,
    get selected() { return selected; },
  };
}
