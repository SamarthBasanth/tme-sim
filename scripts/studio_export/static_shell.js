// static_shell.js — the honesty layer of the static export.
//
// It runs AFTER js/main.js has built the workstation, and does four things:
//
//   1. Puts a banner above the top bar saying, in plain words, that this is a
//      snapshot of a local application, what was precomputed, and the two
//      commands that run the real thing.
//   2. Removes every control that cannot work without a server (submit,
//      validate, cancel) rather than leaving it on screen looking live.
//   3. States the section identity and the sampling rate in the top bar, where
//      it is visible at all times and cannot be scrolled away from.
//   4. Re-purposes the "Active job" panel, which has nothing to poll, into a
//      read-out of the RECORDED engine log for whichever run is open.
//
// It adds no numbers of its own. Everything it prints comes from
// window.__TMESIM_BAKED, which came from the runs.

import { h, clear, fill } from '../util/dom.js';

const B = window.__TMESIM_BAKED;
const M = B.export_meta;

/* ------------------------------------------------------------- banner --- */

function buildBanner() {
  const cmd = (t) => h('code.badge-code', t);

  const runLines = B.runs.map((r) => h('div.xb-run',
    h('span.mono.strong', r.run_id),
    h('span.xb-sep', '·'),
    h('span.mono', (r.spec && r.spec.drug) || 'no drug'),
    h('span.xb-sep', '·'),
    h('span.dim.mono', `${Math.round(r.runtime_s)} s engine time`),
  ));

  return h('div.xbanner', { role: 'note' },
    h('div.xb-mark', 'STATIC EXPORT'),
    h('div.xb-body',
      h('p',
        h('b', 'A single-file snapshot of a local application, not a live one. '),
        'tmesim Studio is a FastAPI + WebGL2 workstation that runs on your machine; '
        + 'this page has no server behind it. Pan, zoom, cell picking, channels, field '
        + 'overlays, run switching and the A/B comparison are live and operate on real '
        + 'data. Submitting, cancelling and validating a run are not — nor are the '
        + 'protocol workspace, therapy reports and the project/specimen browser, which '
        + 'are backend jobs, sweeps and mutable server state. Those surfaces were '
        + 'removed rather than left on screen looking live.'),
      h('p',
        h('b', 'The three runs were precomputed by the production engine '),
        '(', cmd('reconstruct_real'), ' → ', cmd('apply_drug_sensitivity'), ' → ',
        cmd('simulate_response'), ' → ', cmd('delivery_metrics.compute_all'),
        ') on ', h('b', M.section_label), ' at ', cmd(`dx = ${M.dx_um} um`), ' ',
        cmd(`sensitivity_mode = ${M.sensitivity_mode}`), ' ',
        cmd(`plasma_peak = ${M.plasma_peak_uM} uM`), ' ',
        cmd(`total_time = ${M.total_time_h} h`),
        '. Every RunSpec field except ', cmd('drug'), ' is identical, so the comparison '
        + 'isolates transport and potency.'),
      h('div.xb-runs', runLines),
      h('p.xb-warn',
        h('b', `${M.plasma_peak_uM} uM is a temozolomide-scale plasma level and is NOT `
          + 'clinically attainable for cisplatin or doxorubicin. '),
        'Holding PK fixed is what makes the three comparable; read the deltas as model '
        + 'behaviour, not as dosing advice. ',
        h('b', 'NOT FOR CLINICAL USE'),
        ' — this is an untreated resection with no outcome label, and nothing here has '
        + 'been scored against a measured treatment response. The per-run caveats in the '
        + 'right rail are the most important thing on screen and are shown in full.'),
      h('p',
        h('b', 'Sampling: '), M.sampling_note,
        ' Field grids are the full ', h('span.mono', `${M.field_rows}×${M.field_cols}`),
        ' reconstruction mesh at ', h('span.mono', `${M.dx_um} um`), '. ',
        h('span.dim', `Exported ${M.built_at} from ${M.repo}; ${M.payload_mib} MiB of `
          + 'TMSG/TMSF binary inlined, no network request is made by this page.')),
    ),
    h('div.xb-side',
      h('div.hr-lbl', 'run the real thing'),
      h('pre.xb-cmd', M.commands.join('\n'))),
    h('button.btn.ghost.sm.xb-x', {
      title: 'Collapse this notice (the STATIC EXPORT marker in the top bar stays)',
      onclick: (e) => {
        const b = e.target.closest('.xbanner');
        b.classList.toggle('collapsed');
        e.target.textContent = b.classList.contains('collapsed') ? 'what is this?' : 'collapse';
      },
    }, 'collapse'),
  );
}

/* ------------------------------------------------------------ top bar --- */

function patchTopbar() {
  const bar = document.querySelector('.topbar');
  if (!bar) return;

  const spacer = bar.querySelector('.topbar-spacer');
  const group = h('div.topbar-group.keep',
    h('span.lbl', 'sampling'),
    h('span.mono.strong', { title: M.sampling_note }, M.sampling_label),
  );
  bar.insertBefore(group, spacer);

  // The backend indicator: main.js paints it after its probe resolves, which
  // may be after this module runs, so hold it rather than set it once.
  const txt = document.getElementById('hdr-backend');
  const led = document.getElementById('hdr-led');
  const LABEL = 'static export · no backend';
  const TITLE = 'There is no FastAPI backend behind this page. Every response '
    + 'comes from data baked into the file by scripts/build_studio_export.py. '
    + 'Nothing is synthesised: js/mock.js is not bundled.';
  const hold = () => {
    if (txt.textContent !== LABEL) {
      txt.textContent = LABEL;
      txt.style.color = 'var(--fg-1)';
      txt.title = TITLE;
    }
    if (led.className !== 'led') led.className = 'led';
  };
  hold();
  new MutationObserver(hold).observe(bar, {
    subtree: true, childList: true, characterData: true, attributes: true,
  });

  // The status bar's "api ms" cell measures a network round trip that does not
  // happen here. Relabel it for what it now measures.
  const cell = document.getElementById('sb-latency');
  if (cell && cell.parentElement) {
    const lbl = cell.parentElement.querySelector('.lbl');
    if (lbl) {
      lbl.textContent = 'resolve';
      lbl.title = 'Local time to look up and decode a baked payload. There is no '
        + 'network request; this is not an API latency.';
    }
  }
}

/* ------------------------------------------------------- dead controls -- */

function removeDeadControls() {
  // --- submit / validate ------------------------------------------------ //
  const submitBtn = Array.from(document.querySelectorAll('.rail-left .btn'))
    .find((b) => /submit run/i.test(b.textContent));
  const row = submitBtn && submitBtn.parentElement;
  if (row) {
    fill(row,
      h('div.xnote',
        h('b', 'Run submission is disabled.'),
        ' This export has no engine. The RunSpec form and source editor below are '
        + 'live so you can read exactly what produced each run, and edit it to copy '
        + 'out — but ', h('code.badge-code', 'POST /api/runs'), ' needs the local server.'),
    );
  }

  // --- the active-job panel: nothing to poll, so show what was recorded -- //
  const sect = document.getElementById('sect-job');
  if (!sect) return;
  const summary = sect.querySelector('summary');
  if (summary) {
    clear(summary);
    summary.appendChild(document.createTextNode('Engine log — recorded'));
    summary.appendChild(h('span.hint', 'RunResult.log_tail'));
  }
  const body = sect.querySelector('.sect-body');
  if (!body) return;
  const pre = h('pre.xlog', { 'aria-label': 'Recorded engine log for the open run' }, '');
  const head = h('div.frow.xlog-head',
    h('span.lbl', 'run'), h('span.mono.strong.xlog-id', '—'),
    h('span.grow'),
    h('span.chip.st-done', 'terminal'));
  fill(body,
    h('div.xnote',
      'Live ', h('code.badge-code', 'GET /api/runs/{id}/status'), ' polling and ',
      h('code.badge-code', 'DELETE /api/runs/{id}'),
      ' are gone with the server. What is left is the real stage-by-stage log the '
      + 'worker wrote while this run executed.'),
    head, pre);

  // main.js owns the active run; read it rather than duplicating the routing.
  let shown = null;
  setInterval(() => {
    const st = window.__tmesim && window.__tmesim.store && window.__tmesim.store.get();
    const run = st && st.activeRun;
    if (!run || run.run_id === shown) return;
    shown = run.run_id;
    head.querySelector('.xlog-id').textContent = run.run_id;
    const chip = head.querySelector('.chip');
    chip.className = 'chip st-' + (run.status || 'done');
    chip.textContent = run.status || '—';
    pre.textContent = (run.log_tail || ['(no log recorded for this run)']).join('\n');
    pre.scrollTop = 0;
  }, 250);
}

/* -------------------------------------------- mechanism: unhidden notes -- */

/**
 * The mechanism legend prints the limiter roll-up as "the largest share is
 * DELIVERY-limited", but UNREACHED and UNDERDOSED are two bands on ONE axis
 * (E/IC50), so that sentence can read as a spatial-delivery claim when the
 * cause is potency. tmesim.sim.mechanism ships the disclaimer itself, in the
 * TMSM aggregate, under `exposure_axis_note` — the legend just does not render
 * it, and it ships `exposure_shortfall_fold` with it.
 *
 * On temozolomide here that matters: every tumour cell sits at E/IC50 ~ 0.06
 * against a 606 uM IC50, so 97% of survivors land in the "unreached" band while
 * the cell-to-cell SPREAD is only ~2.4x. The shortfall is potency, not
 * transport. This adds no interpretation of its own — it prints the backend's
 * own text and its own numbers, which is the whole point of an export that
 * claims not to hide caveats.
 */
function surfaceMechanismNotes() {
  const note = () => document.querySelector('.mechleg .m-note');
  let shownFor = null;
  let agg = null;

  const paint = () => {
    const el = note();
    if (!el || !agg || el.querySelector('.x-axisnote')) return;
    const p = agg.exposure_ratio_percentiles_tumor || {};
    const fold = agg.exposure_shortfall_fold;
    const box = h('div.x-axisnote');
    if (Number.isFinite(fold) && Number.isFinite(p.p05) && Number.isFinite(p.p95)) {
      box.appendChild(h('div',
        h('b', 'exposure vs potency: '),
        `the median tumour cell needs ${fmtx(fold)} its exposure to reach IC50 `
        + `(${agg.drug && agg.drug.ic50_uM ? nsig(agg.drug.ic50_uM) + ' uM' : 'IC50'}). `
        + `Cell-to-cell spread in E/IC50 is ${fmtx(p.p95 / p.p05)} `
        + `(p05 ${nsig(p.p05)} → p95 ${nsig(p.p95)}).`));
    }
    if (agg.exposure_axis_note) {
      box.appendChild(h('div', { style: null }, agg.exposure_axis_note));
    }
    if (box.childNodes.length) el.appendChild(box);
  };

  const fmtx = (v) => (v >= 10 ? Math.round(v) : v.toFixed(2)) + '\u00d7';
  const nsig = (v) => (v >= 100 ? String(Math.round(v))
    : v >= 1 ? v.toFixed(2) : v.toPrecision(2));

  // The legend re-fills its note element on every render, which drops whatever
  // we appended; put it back rather than racing it.
  const host = document.querySelector('.mechleg');
  if (host) {
    new MutationObserver(paint).observe(host, { childList: true, subtree: true });
  }

  setInterval(async () => {
    const st = window.__tmesim && window.__tmesim.store && window.__tmesim.store.get();
    const id = st && st.activeRunId;
    if (!id || id === shownFor) return;
    shownFor = id;
    agg = null;
    try {
      agg = (await window.__tmesim.api.getMechanism(id)).aggregate;
    } catch {
      return;   // this run has no attribution; main.js already says so
    }
    paint();
  }, 400);
}

/* ----------------------------------------------------------- keyboard --- */

function patchHelpDialog() {
  const dlg = document.getElementById('help-dlg');
  if (!dlg) return;
  for (const tr of dlg.querySelectorAll('tbody tr')) {
    if (/submit the current RunSpec/i.test(tr.textContent)
        || /protocol workspace/i.test(tr.textContent)) tr.remove();
  }
}

/**
 * The protocol workspace ships its button in index.html's markup and wires it
 * from a SECOND module entry point, js/ui/protocol.js, which this export does
 * not bundle: every one of its actions is a POST to /api/protocol/*. A button
 * whose module is absent is inert, so it goes.
 */
function removeProtocolEntryPoint() {
  const btn = document.getElementById('btn-protocol');
  if (btn) btn.remove();
}

/* --------------------------------------------------------------- boot --- */

function install() {
  const shell = document.getElementById('shell');
  if (!shell) return false;
  shell.classList.add('has-xbanner');
  shell.insertBefore(buildBanner(), shell.firstChild);
  patchTopbar();
  removeDeadControls();
  removeProtocolEntryPoint();
  patchHelpDialog();
  surfaceMechanismNotes();
  // Set from the baked value so document.title and the <title> element the
  // Artifact publisher reads can never disagree.
  if (M.document_title) document.title = M.document_title;
  return true;
}

if (!install()) {
  // main.js builds synchronously before its first await, so this should never
  // be needed — but a silent blank export would be worse than a retry.
  const t = setInterval(() => { if (install()) clearInterval(t); }, 50);
  setTimeout(() => clearInterval(t), 5000);
}

export const installed = true;
