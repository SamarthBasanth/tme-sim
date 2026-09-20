// reports.js — the therapy-report panel.
//
// A REPORT is the deliverable: ONE compound on ONE specimen, mechanistically
// attributed, with the provenance of every input and the limits of the output
// attached. (This header said "a ranked, mechanistically attributed drug
// panel" until 2026-09-20 — two words that the rest of this file, the backend
// and `results/` all stopped meaning on 2026-08-21, when compound ordering was
// measured and withdrawn. The code was right and its own first paragraph was
// wrong, which is this repository's most durable defect: the claim surface
// drifting away from the engine.)
//
// This panel does five things and nothing else:
//
//   1. trigger a report for the compound and section in the RunSpec editor,
//   2. watch it build (it is 2 runs of that compound on the shared job queue),
//   3. say what a report that will NEVER produce a document says, and why,
//   4. view the rendered document inline,
//   5. download that document as one self-contained HTML file.
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

/* DEFAULT_PANEL is GONE. It was the five-drug GBM panel, pre-selected, and
 * exporting it made "which five drugs" a decision this module owned. A report
 * is now ONE compound — the one already in the RunSpec editor — so there is no
 * panel to default and nothing to pre-select. The backend rejects a `drugs[]`
 * body outright rather than silently taking its first element, so a client that
 * still sends one gets an error naming the withdrawal instead of a report about
 * a compound nobody chose. Withdrawn 2026-08-21;
 * results/pk_ranking_robustness.md, results/capacity_ranking_robustness.md. */

const STATUS_CHIP = {
  queued: 'st-queued', running: 'st-running', done: 'st-done',
  failed: 'st-failed', cancelled: 'st-cancelled',
};

/** What each state MEANS for the document, which is the thing this panel is
 *  about. "queued" and "failed" are both "no document", and before this they
 *  looked identical on the row: one chip, and a dash in every other column. */
const STATUS_TITLE = {
  queued: 'on the queue — no arm has started, and no document exists yet',
  running: 'building — the document is written when every arm is done',
  done: 'the document exists and can be downloaded',
  failed: 'no document was produced. The server’s reason is on the row and in the viewer',
  cancelled: 'stopped before it finished — no document was produced',
};

/** The same five-state stripe the run history uses (js/ui/history.js). Five
 *  lines, deliberately copied rather than imported: this module must not take
 *  a dependency on the run-history panel to know what colour "failed" is, and
 *  there is no shared style module to put it in that is not a stylesheet. */
const STATUS_INK = {
  queued: 'var(--line-2)',
  running: 'var(--accent)',
  done: 'var(--ok)',
  failed: 'var(--crit)',
  cancelled: 'var(--fg-3)',
};

/** Row-sized form of a long server sentence. The full text is never only
 *  truncated: it is on the row's tooltip and, in full, in the viewer. */
function firstSentence(text, max = 90) {
  const s = String(text || '').trim().replace(/\s+/g, ' ');
  if (!s) return '';
  const stop = s.indexOf('. ');
  if (stop > 0 && stop < max) return s.slice(0, stop + 1);
  if (s.length <= max) return s;
  // On a word boundary: "refused by t …" reads as a typo, not as a cut.
  const sp = s.lastIndexOf(' ', max);
  return s.slice(0, sp > max * 0.6 ? sp : max).trim() + ' …';
}

/** Load a run into the viewer through the shell's published hook — the same
 *  guarded call js/ui/protocol.js makes, and for the same reason: a report's
 *  arms ARE runs, and the only way to look at one is the run the shell owns. */
function openRunInViewer(runId) {
  const t = window.__tmesim;
  if (t && typeof t.openRun === 'function') { t.openRun(runId); return; }
  toast('warn', 'Viewer not reachable',
    `The shell did not expose an openRun hook, so run ${runId} could not be loaded `
    + 'into the viewer from here. It is in the run history, in the dock to the left.',
    7000);
}

/**
 * `", (N running, M queued)"`, or `''` — the per-arm breakdown the backend's
 * `progress.by_status` already carries (report.py `ReportManager.progress`)
 * but that, before this, nothing on screen ever read: the aggregate
 * `runs_done / runs_total` fraction cannot tell an operator "which arms",
 * only "how many". Shared by the list row's tooltip and the viewer's
 * in-progress note so the two don't drift into two different phrasings of
 * the same number, the way the dock's two `refresh` buttons once did.
 */
function armBreakdown(progress) {
  const by = (progress && progress.by_status) || {};
  const parts = [];
  if (by.running) parts.push(`${by.running} running`);
  if (by.queued) parts.push(`${by.queued} queued`);
  return parts.length ? ` (${parts.join(', ')})` : '';
}

export function createReportsPanel({ store, getSection }) {
  let reports = [];
  let selected = null;
  let submitting = false;
  /** What the open viewer was last painted FROM, so the poll can tell a real
   *  change from a no-op — including the done transition. */
  let viewerStatus = null;
  let viewerFraction = null;

  /* ---------------------------------------------------- the compound --- *
   * There is no picker here any more, and that is the point of the change.
   * The multi-select checkbox box asked the operator to assemble a panel, and
   * a panel exists to be ordered. The report's compound is now simply the one
   * in the RunSpec editor — the same compound the viewport, the metrics, the
   * reach map and the mechanism breakdown are all already about — so the
   * document and the workstation cannot disagree about their subject. */
  const drugBox = h('div.rp-drugs');

  function currentDrug() {
    const spec = (store && store.get().spec) || {};
    return spec.drug || null;
  }

  /** Kept in the returned interface under its old name so main.js's catalogue
   *  refresh keeps working; it no longer builds a picker, it repaints the
   *  subject line. */
  function buildDrugs() { paintHead(); }

  const countEl = h('span.hint');
  function paintHead() {
    const sec = getSection() || NULL_GLYPH;
    const drug = currentDrug();
    const label = `${drug || 'no compound'} · ${sec}`;
    countEl.textContent = label;
    // This column is the narrowest of the dock's three (reports.css:
    // minmax(300px, 1.05fr)), and the head row also carries `refresh
    // reports` and `generate report`. `.panel-head .hint` truncates with an
    // ellipsis at that width (layout.css) — confirmed on screen shrinking to
    // "5..." for a real section id, same as the project title pass 3 found
    // wrapping/clipping and fixed with a hover title; this one just never
    // got the same fix. One hover away instead of gone.
    countEl.title = label;
    clear(drugBox);
    drugBox.appendChild(h('div', {
      style: {
        fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)',
        color: 'var(--fg-1)', padding: 'var(--sp-2) 0',
      },
    }, drug
      ? ['The report is about ', h('b', drug),
        ', the compound currently in the RunSpec editor. Change it there and '
        + 'the report follows — a document and the workstation that produced '
        + 'it should never disagree about their subject.']
      : ['No compound is selected. A report is one compound on one specimen, '
        + 'so pick a drug in the RunSpec editor first — this used to accept a '
        + 'panel of several, and that was withdrawn on 2026-08-21 after the '
        + 'ordering it produced was measured and did not hold.']));
    btnNew.disabled = submitting || !drug;
  }

  /* ------------------------------------------------------- submission -- */
  const btnNew = h('button.btn.primary.sm', {
    title: 'POST /api/reports — enqueues 2 runs of ONE compound (the real arm and the resistance-off counterfactual) on the same reconstruction',
    onclick: submit,
  }, 'generate report');

  async function submit() {
    const section = getSection();
    if (!section) {
      toast('warn', 'No section selected',
        'A report is a property of one specimen. Pick a section in the RunSpec editor first.');
      return;
    }
    const drug = currentDrug();
    if (!drug) {
      toast('warn', 'No compound selected',
        'A report is one compound on one specimen. Pick a drug in the RunSpec '
        + 'editor first.');
      return;
    }
    submitting = true; paintHead();
    const spec = (store && store.get().spec) || {};
    // Carry the operator's own grid / window / PK across; the backend refuses
    // and reports back any field a report owns (drug, sensitivity_mode).
    const spec_overrides = {};
    for (const k of ['dx_um', 'total_time_s', 'vascular_mode', 'oxygen_mode', 'state_mode', 'mgmt', 'pk']) {
      if (spec[k] !== undefined) spec_overrides[k] = spec[k];
    }
    try {
      const res = await api.postReport({ section, drug, spec_overrides });
      toast('ok', `Report ${res.report_id} queued`,
        `${(res.run_ids || []).length} runs of ${drug}: the real arm and the `
        + 'resistance-off counterfactual, on one shared reconstruction. The '
        + 'counterfactual is what MEASURES the resistance term rather than '
        + 'asserting it.',
        6000);
      // A non-default spec can trigger several of these at once (grid
      // coarseness, PK window, sensitivity mode, ...). They used to toast one
      // per warning with ttl=0 (sticky, never auto-dismissed) so the text
      // would not be lost before anyone opened the report — but `.toasts` is
      // `position: fixed; right; bottom` (components.css), which is exactly
      // where this dock column lives, so 3-4 real warnings buried the entire
      // Therapy Reports panel — status chips, the progress bar, the table,
      // `generate report` itself — under an unremovable stack the operator
      // had to dismiss one at a time before they could see the panel they
      // just used. They now have a permanent home in the report viewer
      // (paintViewer, below), so the toast can be transient like every other
      // 'warn' toast in the app instead of standing in for storage.
      const warns = res.warnings || [];
      if (warns.length === 1) {
        toast('warn', 'Report spec warning', warns[0]);
      } else if (warns.length > 1) {
        toast('warn', `${warns.length} report spec warnings`, 'Open the report to read them.');
      }
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
  const COL_COUNT = 8;
  const table = h('table.tbl',
    h('thead', h('tr',
      h('th', 'report'),
      h('th', 'status'),
      h('th', 'section'),
      h('th.right', { title: 'arms done of arms in this report. One report is 2 runs: the real arm and the resistance-off counterfactual.' }, 'arms'),
      h('th', 'compound'),
      h('th.right', { title: 'fraction of tumour cells the compound never reached in a pharmacologically meaningful amount. It exists only once every arm is done.' }, 'unreached'),
      h('th.right', 'runtime'),
      h('th', { 'aria-label': 'report action', title: 'download, cancel or open' }, ''))),
    tbody);

  /** "3 reports · 1 building · 1 failed" — the panel's own state, above the
   *  table, so "is anything happening" is answered without reading rows. */
  const tallyEl = h('span.hint');

  function paintTally() {
    const by = {};
    for (const r of reports) by[r.status] = (by[r.status] || 0) + 1;
    const bits = [`${reports.length} report${reports.length === 1 ? '' : 's'}`];
    const building = (by.queued || 0) + (by.running || 0);
    if (building) bits.push(`${building} building`);
    if (by.done) bits.push(`${by.done} with a document`);
    if (by.failed) bits.push(`${by.failed} failed`);
    if (by.cancelled) bits.push(`${by.cancelled} cancelled`);
    tallyEl.textContent = bits.join(' · ');
    tallyEl.title = reports.length
      ? 'Every report this server holds, by state. A report with a document can be '
        + 'downloaded from its row; one that is building has no document yet.'
      : 'No report has been generated on this server yet.';
  }

  // This column was `top-ranked` + `kill`: a RANKING and its winning number,
  // both readable off this list without ever opening the document. Both are
  // gone. What replaces them is the compound the report is about and the
  // per-cell unreached fraction — the one number on this row that is a property
  // of the tissue and the compound rather than of a comparison.
  //
  // The delivery-law disclosure STAYS, and it is re-pointed rather than
  // removed. It used to qualify the top pick; it now qualifies the map, which
  // is the stronger version of the same warning — the source term sets the
  // SHAPE of the exposure field, and a shape is the one thing a within-map
  // ratio does not divide out. Populated from the report summaries, which carry
  // the backend's runtime read of PDEConfig (app/server/delivery_law.py), so
  // nothing about the law is written here and this line cannot go stale if the
  // default flips.
  const lawEl = h('div');

  function renderLaw() {
    clear(lawEl);
    const law = (reports.find((r) => r.delivery_law) || {}).delivery_law;
    if (!law || !law.headline) return;
    lawEl.appendChild(h('div.caveat.crit', { title: law.caveat },
      h('span', h('b', 'THE MAP IS LAW-DEPENDENT. '), law.headline)));
  }

  function focusReportRow(fromId, step) {
    const rows = Array.from(tbody.querySelectorAll('tr[data-report]'));
    const i = rows.findIndex((el) => el.dataset.report === fromId);
    const next = step === 'first' ? rows[0]
      : step === 'last' ? rows[rows.length - 1]
        : rows[i + step];
    if (next) next.focus();
  }

  /** The build state, as a bar AND as words. A report is 2 runs on a shared
   *  queue, so "queued" and "running" are genuinely different answers to "is
   *  anything happening": the first has nothing in flight at all. */
  function progressBit(r) {
    const p = r.progress || {};
    if (r.status === 'running') {
      return h('span.rp-prog', {
        title: `${p.runs_done || 0} of ${p.runs_total || 0} arms done${armBreakdown(p)}`
          + `${p.stage ? ' · ' + p.stage : ''}`,
      }, h('i', { style: { width: pct(p.fraction || 0, 0) + '%' } }));
    }
    if (r.status === 'queued') {
      // A stub, not an animation: prefers-reduced-motion must not turn an
      // indeterminate bar into a strobe (docs/UI_ARCHITECTURE.md §4.2), and
      // `.progress.indet` already carries that contract.
      return h('span.rp-prog', {
        title: `0 of ${p.runs_total || r.n_runs || 2} arms have started${armBreakdown(p)}`,
      }, h('i', { style: { width: '8%' } }));
    }
    return null;
  }

  function renderList() {
    clear(tbody);
    renderLaw();
    paintTally();
    if (!reports.length) {
      tbody.appendChild(h('tr', h('td', { colspan: String(COL_COUNT), class: 'dim' },
        'No reports yet. Pick a compound in the RunSpec editor and press ',
        h('b', 'generate report'), '.')));
      return;
    }
    for (const r of reports) {
      const p = r.progress || {};
      const head = r.headline || {};
      const isDone = r.status === 'done';
      const isSel = r.report_id === selected;
      // The stripe is skipped on the selected row: components.css already
      // draws the accent inset there, and an inline style would win and take
      // the selection cue away.
      const firstStyle = isSel ? {}
        : { boxShadow: `inset 2px 0 0 ${STATUS_INK[r.status] || 'var(--line-2)'}` };
      const tr = h('tr.sel-able', {
        tabindex: '0',
        dataset: { report: r.report_id, status: r.status || 'queued' },
        'aria-selected': String(isSel),
        'aria-label': `report ${r.report_id}, ${r.status}`,
        onclick: (e) => { if (e.target.tagName !== 'BUTTON') open(r.report_id); },
        onkeydown: (e) => {
          if (e.target !== tr) return;
          if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(r.report_id); }
          else if (e.key === 'ArrowDown') { e.preventDefault(); focusReportRow(r.report_id, 1); }
          else if (e.key === 'ArrowUp') { e.preventDefault(); focusReportRow(r.report_id, -1); }
          else if (e.key === 'Home') { e.preventDefault(); focusReportRow(r.report_id, 'first'); }
          else if (e.key === 'End') { e.preventDefault(); focusReportRow(r.report_id, 'last'); }
        },
      },
      h('td.mono', { style: firstStyle }, r.report_id),
      h('td',
        h('span.chip.' + (STATUS_CHIP[r.status] || 'st-queued'),
          { title: STATUS_TITLE[r.status] || r.status }, r.status),
        progressBit(r)),
      h('td', { title: r.section }, r.section || NULL_GLYPH),
      h('td.right.num', {
        title: `${p.runs_done ?? 0} of ${r.n_runs ?? 0} arms done${armBreakdown(p)}`,
      }, `${p.runs_done ?? 0}/${r.n_runs ?? 0}`),
      h('td', r.drug ? h('b', r.drug) : h('span.dim', NULL_GLYPH)),
      // The unreached fraction, or the backend's REASON it is undefined —
      // never a dash standing in for a refusal nobody read. And for a report
      // that is still building, "pending" rather than a refusal it never made.
      h('td.right.num',
        Number.isFinite(head.unreached_fraction_of_tumor)
          ? h('span', { title: 'fraction of tumour cells the drug never reached in a pharmacologically meaningful amount (per-cell mechanism attribution)' },
            pct(head.unreached_fraction_of_tumor, 2) + '%')
          : (r.status === 'queued' || r.status === 'running')
            ? h('span.dim', { title: 'This number is read off the finished attribution; it exists when every arm is done.' }, 'pending')
            : h('span.dim', {
              title: head.unreached_null_reason
                || (r.status === 'failed' || r.status === 'cancelled'
                  ? `This report is ${r.status}, so no attribution was assembled and `
                    + 'there is no number to read. That is not the backend declining '
                    + 'to report one.'
                  : 'not reported by the backend'),
            }, NULL_GLYPH)),
      h('td.right.num', {
        title: isDone ? 'wall-clock seconds from submission to document'
          : 'elapsed so far',
      }, dur(r.runtime_s)),
      h('td.right',
        isDone
          ? h('button.btn.ghost.sm', {
            dataset: { act: 'download' },
            title: 'GET /api/reports/{id}/html — download the self-contained HTML document',
            onclick: (e) => { e.stopPropagation(); download(r.report_id); },
          }, 'download')
          : (r.status === 'running' || r.status === 'queued')
            ? h('button.btn.ghost.sm', {
              dataset: { act: 'cancel' },
              title: 'DELETE /api/reports/{id} — cancels every queued arm',
              onclick: (e) => { e.stopPropagation(); cancel(r.report_id); },
            }, 'cancel')
            : h('button.btn.ghost.sm', {
              dataset: { act: r.status === 'failed' ? 'why' : 'open' },
              title: r.status === 'failed'
                ? 'Open this report: the failure in full, the spec warnings and the limits'
                : 'Open this report',
              onclick: (e) => { e.stopPropagation(); open(r.report_id); },
            }, r.status === 'failed' ? 'why' : 'open')),
      );
      tbody.appendChild(tr);

      // A FAILED report states its failure on the list, not only behind a
      // click: it is the one row whose numbers will never arrive, and it used
      // to be a row of dashes with the reason hidden in a `title` attribute
      // that no keyboard and no touch device can reach.
      // FAILED only. A cancelled report is something the operator did on
      // purpose; a second row explaining it to them on every repaint is noise,
      // and its chip's tooltip plus the viewer already carry the sentence. A
      // failure is the one that arrives unasked.
      if (r.status === 'failed') {
        tbody.appendChild(h('tr', { dataset: { detail: r.report_id } },
          h('td', {
            colspan: String(COL_COUNT),
            style: {
              height: 'auto', whiteSpace: 'normal', overflow: 'visible',
              background: 'var(--bg-1)', padding: 'var(--sp-3) var(--sp-4)',
              boxShadow: `inset 2px 0 0 ${STATUS_INK[r.status]}`,
            },
          },
          h('div', { style: { width: '420px', maxWidth: '100%' } },
            h('div.caveat.crit', h('span',
              h('b', 'No document was produced. '),
              r.error
                ? [h('span', { title: r.error }, firstSentence(r.error, 160)),
                  ' Open the report for this sentence in full, the spec '
                  + 'warnings it was submitted with, and the arms that ran.']
                : 'The server reported this report as failed and sent no reason '
                  + 'with it. Read the console of the process running Vivum.'))))));
      }
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
    const isDone = r.status === 'done';
    viewerStatus = r.status;
    viewerFraction = (r.progress || {}).fraction ?? null;
    clear(dlgHead);
    dlgHead.appendChild(h('h2', 'Therapy report'));
    dlgHead.appendChild(h('span.mono.strong', r.report_id));
    dlgHead.appendChild(h('span.chip.' + (STATUS_CHIP[r.status] || 'st-queued'),
      { title: STATUS_TITLE[r.status] || r.status }, r.status));
    dlgHead.appendChild(h('span.chip.crit', 'not for clinical use'));
    dlgHead.appendChild(h('span.grow'));
    // Disabled rather than absent, and it says why: a button that vanishes
    // while a report builds reads as "this report has no download", which is
    // a different claim from "not yet".
    dlgHead.appendChild(h('button.btn.ghost.sm', {
      disabled: !isDone,
      onclick: () => download(r.report_id),
      title: isDone
        ? 'GET /api/reports/{id}/html — one self-contained file, no external requests'
        : `There is no document to download: this report is ${r.status}. The HTML is `
          + 'written once every arm is done.',
    }, 'download'));
    dlgHead.appendChild(h('button.btn.ghost.sm', { onclick: () => dlg.close() }, 'close'));

    clear(dlgNote);
    if (isDone) {
      frame.src = api.reportHtmlUrl(r.report_id);
      frame.style.display = '';
    } else {
      frame.removeAttribute('src');
      frame.style.display = 'none';
      const p = r.progress || {};
      if (r.status === 'failed') {
        // BUILDING and FAILED were one branch and one caveat colour apart.
        // They are opposite states — one ends in a document and one never
        // will — so the failed one leads with that, quotes the server in
        // full, and offers the only action that can change it.
        dlgNote.appendChild(h('div.caveat.crit', h('span',
          h('b', 'No document was produced, and none will be. '),
          r.error
            ? ['The server said: ', h('em', r.error)]
            : 'The server sent no reason with the failure. It is recorded '
              + 'server-side: read the console of the process running Vivum.')));
        dlgNote.appendChild(h('div.btn-row',
          h('button.btn.sm', {
            title: 'POST /api/reports with this report’s own section, compound and '
              + 'spec overrides — the same request that produced this one.',
            onclick: () => regenerate(r),
          }, 'generate again with the same spec')));
      } else if (r.status === 'cancelled') {
        dlgNote.appendChild(h('div.caveat', h('span',
          h('b', 'Cancelled. '),
          'Its queued arms were stopped before the document could be assembled, '
          + 'so there is nothing to render. Any arm that had already finished is '
          + 'still a run in the run history.')));
      } else {
        const frac = Number.isFinite(p.fraction) ? p.fraction : 0;
        const prog = h('div.progress' + (r.status === 'queued' ? '.indet' : ''),
          { style: { maxWidth: '260px', margin: 'var(--sp-4) 0' } },
          h('i', { style: { width: pct(frac, 0) + '%' } }));
        dlgNote.appendChild(h('div.caveat.info', h('span',
          h('b', r.status === 'queued' ? 'Queued. ' : 'Building. '),
          `${p.runs_done || 0} of ${p.runs_total || r.n_runs || 0} arms complete`
          + `${armBreakdown(p)}${p.stage ? ` · ${p.stage}` : ''}. `
          + 'The document is written when every arm is done — a half-finished '
          + 'attribution is not a result, so nothing partial is rendered here. '
          + 'This panel re-reads the report while it builds; the download button '
          + 'above turns on by itself.',
          prog)));
      }
      // The arms ARE runs, and while a report builds they are the only thing
      // with any detail to show. Naming them turns "watch it build" into
      // something you can actually open.
      const arms = r.arms || {};
      const armRows = Object.keys(arms).length
        ? Object.keys(arms).map((k) => [k, arms[k]])
        : (r.run_ids || []).map((id, i) => [`arm ${i + 1} of ${(r.run_ids || []).length}`, id]);
      if (armRows.length) {
        dlgNote.appendChild(h('div.hr-lbl', 'Arms'));
        dlgNote.appendChild(h('div', {
          style: {
            display: 'flex', flexWrap: 'wrap', gap: 'var(--sp-3) var(--sp-4)',
            padding: '0 var(--sp-5) var(--sp-4)',
          },
        }, armRows.map(([name, id]) => h('button.btn.ghost.sm', {
          title: `Open run ${id} — the ${name} arm — in the viewer and the run history.`,
          onclick: () => openRunInViewer(id),
        }, `${name}: ${id}`))));
      }
    }
    // Spec warnings from submission time (POST /api/reports) — see the
    // comment in submit() for why they no longer live in a stack of sticky
    // toasts. They apply for the life of the report, same as the limits
    // below, so they render the same way every time this paints.
    for (const w of (r.warnings || [])) {
      dlgNote.appendChild(h('div.caveat', h('span', h('b', 'spec warning: '), w)));
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
      const res = await fetch(api.reportHtmlUrl(id), {
        cache: 'no-store',
        headers: api.api.apiKey ? { 'x-api-key': api.api.apiKey } : {},
      });
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

  /**
   * Re-submit a failed report EXACTLY as it was submitted.
   *
   * `r` here is the FULL view (`GET /api/reports/{id}`), which is the only
   * place `spec_overrides` exists — the list summary does not carry it. That
   * matters: re-generating from the summary would silently substitute the
   * RunSpec editor's CURRENT grid and PK for the ones this report used, and
   * call the result the same request.
   */
  async function regenerate(r) {
    if (!r.section || !r.drug) {
      toast('warn', 'Cannot re-generate this report',
        'This report record does not carry both a section and a compound, so the '
        + 'request that produced it cannot be reproduced from it.');
      return;
    }
    try {
      const res = await api.postReport({
        section: r.section, drug: r.drug, spec_overrides: r.spec_overrides || {},
      });
      toast('ok', `Report ${res.report_id} queued`,
        `Same section (${r.section}), same compound (${r.drug}) and the same spec `
        + `overrides as ${r.report_id}, which failed. If it fails the same way, the `
        + 'failure is in the request or the specimen, not in the queue.', 7000);
      selected = res.report_id;
      await refresh();
      try { paintViewer(await api.getReport(res.report_id)); } catch { /* the list still updated */ }
    } catch (err) {
      notifyError('POST /api/reports failed', err);
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
        // REPAINT ON ANY CHANGE, not on "not done yet". The old condition was
        // `cur.status !== 'done'`, which is false at exactly the moment the
        // document appears — so a viewer opened on a building report kept
        // saying "Building — 1 of 2 arms" for as long as it was left open,
        // with the finished document one poll away and never rendered. The
        // one transition this panel exists to show was the one it skipped.
        const changed = cur && (cur.status !== viewerStatus
          || ((cur.progress || {}).fraction ?? null) !== viewerFraction);
        if (changed) {
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
      h('span.grow'),
      h('button.btn.ghost.sm', {
        onclick: refresh,
        title: 'Re-read the report list. GET /api/reports',
      }, 'refresh reports'),
      btnNew),
    h('div.panel-body.rp-body',
      h('div.rp-lead',
        // The headline is read before the mitigation that follows it, so it
        // states the deliverable as what it is: a hypothesis to test at the
        // bench. It used to read "One specimen in, one ranked hypothesis" and
        // to end "the ranking is the output". Compound ordering was withdrawn
        // on 2026-08-21 after it was measured (results/pk_ranking_robustness.md,
        // results/capacity_ranking_robustness.md), and the document is now
        // about one compound on one specimen.
        h('b', 'One specimen and one compound in, one discrimination out: '
          + 'was it delivery, was it dose, or was it the cell’s own state?'),
        ' The compound is simulated twice on the same cached reconstruction — the real run, and a resistance-off counterfactual — so the O',
        h('sub', '2'),
        '-resistance term of the attribution is measured, not asserted. The document carries its own limits: ',
        h('b', 'the map is the output; the absolute level is not, and neither '
          + 'is a comparison against another compound.')),
      // This line used to be a SELECTION COUNT ("5 drugs · section"), living in
      // the panel head, where it competed with two unshrinkable buttons for a
      // 300 px column: it ellipsed to "5..." at the real dock width, and the
      // overflow it caused pushed `generate report` off the right edge of the
      // viewport entirely. There is nothing to count now, so it names the
      // subject instead — but it stays out of the head for the same reason.
      h('div.rp-count', h('span.lbl', 'subject'), countEl),
      drugBox,
      errEl,
      // The list's own state, in the same place every time, so "is anything
      // building" never requires reading rows. Same line shape as `subject`.
      h('div.rp-count', h('span.lbl', 'reports'), tallyEl),
      h('div.tblwrap', table),
      lawEl),
  );

  document.body.appendChild(dlg);
  // The report's subject is the RunSpec's compound, so this panel has to
  // repaint when that changes. Before, the subject was a checkbox set this
  // module owned and nothing outside it could move; now it is owned elsewhere,
  // and a stale subject line would let an operator press `generate report`
  // believing it is about a compound it is not.
  if (store && store.on) store.on('spec', paintHead);
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
