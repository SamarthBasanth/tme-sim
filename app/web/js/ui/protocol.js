// protocol.js — the TPL workspace: script pane, drug upload, plan-before-run,
// and the sweep result grid.
//
// It mounts itself over the shell's body region (the top bar and the status
// strip stay visible and keep working) and is toggled from the top bar or with
// P. It talks to app/server only through js/protocol/api.js and reads the
// existing shell through `window.__tmesim` — it does not reach into the viewer,
// the spec panel, or the history dock.
//
// THE FOUR THINGS THIS SCREEN EXISTS TO ENFORCE
//   1. A diagnostic lands on its LINE. A refusal from the executor is a
//      diagnostic in the gutter with a remedy you can apply, not a toast that
//      disappears while you are reading it.
//   2. An uploaded compound shows its resolved DrugParams with PER-FIELD
//      provenance BEFORE anything runs. A flagged default is visible at upload
//      time or the upload was pointless.
//   3. Nothing launches without a plan. A two-axis sweep is dozens of PDE
//      solves; the run count and the cost are on screen with an explicit
//      confirm before a single job is queued.
//   4. Results are a grid, every cell carries its provenance, and the standing
//      caveats sit next to the numbers rather than behind a disclosure.

import { h, fill, clear, debounce } from '../util/dom.js';
import { createEditor } from '../editor/editor.js';
import { TEMPLATES } from '../editor/templates.js';
import * as papi from '../protocol/api.js';
import { METRICS as TPL_METRICS } from '../protocol/stub.js';
import { chip, metricProvenance } from './provenance.js';
import { toast, notifyError } from './toast.js';
import { makeSplitter } from './splitter.js';
import { num, int, pct, sci, dur, bytes, NULL_GLYPH, nice125 } from '../util/fmt.js';
import { sample, gradientCss } from '../gl/ramps.js';

const RAMP = 'viridis';
const SVG_NS = 'http://www.w3.org/2000/svg';

/* ------------------------------------------------------------ helpers --- */

function svg(tag, attrs, ...kids) {
  const el = document.createElementNS(SVG_NS, tag);
  for (const k in attrs || {}) if (attrs[k] !== null && attrs[k] !== undefined) el.setAttribute(k, String(attrs[k]));
  for (const c of kids) {
    if (c === null || c === undefined || c === false) continue;
    el.appendChild(c.nodeType ? c : document.createTextNode(String(c)));
  }
  return el;
}

const metricMeta = (key, fromPlan) => {
  const fromServer = (fromPlan || []).find((m) => m.key === key);
  return fromServer || TPL_METRICS[key] || { key, label: key, unit: '', fmt: 'auto' };
};

function fmtMetric(meta, v) {
  if (v === null || v === undefined || Number.isNaN(v)) return NULL_GLYPH;
  switch (meta.fmt) {
    case 'pct': return pct(v, 1);
    case 'sci': return sci(v, 2);
    case 'num1': return num(v, 1);
    default: return num(v, 3);
  }
}

/**
 * The value in the units TPL itself uses. `expect kill > 0.05` compares against
 * a FRACTION, so an assertion row that printed "39.6" beside "> 0.05" would
 * read as a pass by a factor of a hundred. Assertions are shown raw.
 */
function fmtRaw(meta, v) {
  if (v === null || v === undefined || Number.isNaN(v)) return NULL_GLYPH;
  return meta.fmt === 'pct' ? num(v, 4) : fmtMetric(meta, v);
}

function fmtQty(q) {
  if (!q) return NULL_GLYPH;
  if (typeof q === 'number') return num(q, q % 1 === 0 ? 0 : 2);
  const v = Math.abs(q.value) >= 100 ? num(q.value, 0) : num(q.value, Math.abs(q.value) < 1 ? 3 : 1);
  return q.unit ? `${v} ${q.unit}` : v;
}

/** Text colour that survives on a given ramp sample. */
function inkOn(rgb) {
  const m = /rgb\((\d+),(\d+),(\d+)\)/.exec(rgb);
  if (!m) return 'var(--fg-0)';
  const l = (0.2126 * +m[1] + 0.7152 * +m[2] + 0.0722 * +m[3]) / 255;
  return l > 0.55 ? 'hsl(213 30% 8%)' : 'hsl(0 0% 100%)';
}

/** Provenance of one metric on one run, in the shell's result taxonomy. */
function cellProvenance(tplKey, run) {
  if (run && run.provenance && run.provenance.kind) {
    return { kind: run.provenance.kind, why: run.provenance.why || run.provenance.note || '' };
  }
  const specKey = (TPL_METRICS[tplKey] || {}).spec;
  if (specKey && run && run.result) return metricProvenance(specKey, run.result);
  return { kind: 'unknown', why: 'The executor did not state a provenance for this cell.' };
}

/* ================================================================ panel == */

export function createProtocolWorkspace() {
  let planResult = null;
  let planSource = '';
  let sweepId = null;
  let sweepData = null;
  let pollTimer = 0;
  let selectedIndex = -1;
  let metricKey = null;
  let validateToken = 0;
  let uploadData = null;

  /* --------------------------------------------------------- editor ---- */

  const editor = createEditor({
    value: TEMPLATES[0].source,
    title: 'protocol.tpl',
    ariaLabel: 'TPL protocol source',
    onChange: (text) => runValidate(text),
    onInput: () => { if (planResult) markPlanStale(); },
    onRun: () => doPlan(),
    onValidateNow: () => runValidate(editor.value, true),
  });

  const diagList = h('div.pdiags', { role: 'list', 'aria-label': 'Protocol diagnostics' });

  /* ------------------------------------------------------ templates ---- */

  const tplList = h('div.ptpl', TEMPLATES.map((t) => h('button.ptpl-item', {
    type: 'button',
    title: `Replace the buffer with the "${t.label}" template`,
    onclick: () => loadTemplate(t),
  },
    h('span.t', t.label),
    h('span.r', `${t.runs} run${t.runs === 1 ? '' : 's'}`),
    h('span.n', t.note),
  )));

  function loadTemplate(t) {
    if (editor.value.trim() && !window.confirm(`Replace the editor buffer with the "${t.label}" template? Unsaved edits are lost.`)) return;
    editor.setValue(t.source);
    editor.flush();
    clearPlan();
    editor.focus();
  }

  /* --------------------------------------------------------- upload ---- */

  const fileInput = h('input', {
    type: 'file',
    accept: '.smi,.smiles,.mol,.sdf,.json,.yaml,.yml',
    style: { display: 'none' },
    onchange: (e) => { if (e.target.files && e.target.files[0]) doUpload(e.target.files[0]); e.target.value = ''; },
  });

  const drop = h('div.pdrop', {
    tabindex: '0',
    role: 'button',
    'aria-label': 'Upload a drug definition file',
    onclick: () => fileInput.click(),
    onkeydown: (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); fileInput.click(); } },
    ondragover: (e) => { e.preventDefault(); drop.classList.add('over'); },
    ondragleave: () => drop.classList.remove('over'),
    ondrop: (e) => {
      e.preventDefault();
      drop.classList.remove('over');
      const f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
      if (f) doUpload(f);
    },
  },
    h('span', 'Drop a drug file here, or click to choose'),
    h('span.fmt', '.smi  .mol  .sdf  .json  .yaml'),
    h('span.fmt', 'POST /api/protocol/upload → resolved DrugParams with per-field provenance'),
  );

  const uploadOut = h('div');

  /* ----------------------------------------------------------- plan ---- */

  const planOut = h('div');
  const confirmBar = h('div.pconfirm.hidden');

  /* -------------------------------------------------------- results ---- */

  const resultsOut = h('div');

  /* ------------------------------------------------------------ shell -- */

  const stubNote = h('div.pws-stub.hidden');

  const btnValidate = h('button.btn.sm', {
    onclick: () => { editor.flush(); runValidate(editor.value, true); },
    title: 'POST /api/protocol/validate',
  }, 'validate');
  const btnPlan = h('button.btn.sm', {
    onclick: () => doPlan(),
    title: 'POST /api/protocol/plan — expands the sweep and states the cost. Nothing runs.',
  }, 'plan');
  const btnRun = h('button.btn.primary.sm', {
    onclick: () => doPlan(true),
    title: 'Plan, then ask for confirmation before queueing anything (Ctrl+Enter)',
  }, 'Run…');
  const btnClose = h('button.btn.ghost.sm', { onclick: () => hide(), title: 'Close the protocol workspace (Esc)' }, 'close');

  const leftCol = h('div.pws-col',
    editor.el,
    h('details.sect', { open: true },
      h('summary', 'Diagnostics', h('span.hint', '/api/protocol/validate')),
      diagList),
    h('details.sect', { open: false },
      h('summary', 'Templates', h('span.hint', `${TEMPLATES.length} programs`)),
      tplList),
  );

  const splitEl = h('div.pws-split');

  const rightCol = h('div.pws-col.right',
    h('div.pws-scroll',
      h('details.sect', { open: true },
        h('summary', 'Drug file', h('span.hint', '/api/protocol/upload')),
        h('div.sect-body', drop, fileInput, uploadOut)),
      h('details.sect', { open: true, id: 'pws-plan' },
        h('summary', 'Plan', h('span.hint', 'before anything runs')),
        h('div.sect-body', planOut)),
      h('details.sect', { open: true, id: 'pws-results' },
        h('summary', 'Sweep results', h('span.hint', 'grid')),
        resultsOut),
    ),
    confirmBar,
  );

  const root = h('section.pws', {
    id: 'protocol-ws',
    hidden: true,
    'aria-label': 'Protocol editor and sweep workspace',
  },
    h('div.pws-top',
      h('div.pws-head',
        h('h2', 'Protocol'),
        h('span.hint.mono', 'TPL v1'),
        h('span.vrule', { 'aria-hidden': 'true' }),
        btnValidate, btnPlan, btnRun,
        h('span.grow'),
        h('span.mono.dim', { id: 'pws-mode' }, ''),
        btnClose),
      stubNote),
    h('div.pws-body', leftCol, splitEl, rightCol),
  );

  makeSplitter(splitEl, {
    axis: 'x', varName: '--pws-l', min: 320, max: 1100,
    onChange: () => { editor.repaint(); scheduleCurveRedraw(); },
  });

  /* ================================================== validation ======= */

  function normDiag(d, fallbackSeverity) {
    if (typeof d === 'string') {
      const m = /^(?:line\s+)?(\d+)\s*[:,]\s*(.+)$/i.exec(d);
      return m
        ? { line: Number(m[1]), message: m[2], severity: fallbackSeverity }
        : { line: 0, message: d, severity: fallbackSeverity };
    }
    if (!d || typeof d !== 'object') return { line: 0, message: String(d), severity: fallbackSeverity };
    const line = Number(d.line ?? (Array.isArray(d.loc) ? d.loc[0] : 0)) || 0;
    const col = Number.isFinite(d.col) ? d.col : (Array.isArray(d.loc) ? d.loc[1] : undefined);
    return {
      line,
      col: Number.isFinite(col) ? col : undefined,
      end_col: Number.isFinite(d.end_col) ? d.end_col : undefined,
      severity: d.severity || fallbackSeverity,
      message: d.message || d.msg || d.detail || JSON.stringify(d),
      code: d.code || null,
      remedy: d.remedy || null,
      remedy_text: d.remedy_text || null,
      remedy_kind: d.remedy_kind || (d.remedy_text ? 'insert_line' : null),
      remedy_line: Number.isFinite(d.remedy_line) ? d.remedy_line : line,
      source: d.source || null,
    };
  }

  async function runValidate(text, explicit) {
    const token = ++validateToken;
    editor.setStatus('validating…', 'busy');
    let res;
    try {
      res = await papi.validate(text);
    } catch (err) {
      if (token !== validateToken) return;
      editor.setStatus('validate failed', 'err');
      renderDiags([{
        line: 1, severity: 'error', message: `POST /api/protocol/validate failed: ${err.message}`,
        code: err.status ? `HTTP ${err.status}` : 'network',
      }]);
      if (explicit) notifyError('POST /api/protocol/validate failed', err);
      return;
    }
    if (token !== validateToken) return;
    paintMode();
    const errors = (res.errors || []).map((d) => normDiag(d, 'error'));
    const warnings = (res.warnings || []).map((d) => normDiag(d, 'warning'));
    const all = errors.concat(warnings);
    renderDiags(all);
    const ok = res.ok !== undefined ? !!res.ok : errors.length === 0;
    btnRun.disabled = !ok;
    btnPlan.disabled = !ok;
    editor.setStatus(
      ok ? (warnings.length ? `valid · ${warnings.length} warning${warnings.length === 1 ? '' : 's'}` : 'valid')
        : `${errors.length} error${errors.length === 1 ? '' : 's'}`,
      ok ? (warnings.length ? 'warn' : 'ok') : 'err');
  }

  function renderDiags(list) {
    editor.setDiagnostics(list);
    const frag = document.createDocumentFragment();
    if (!list.length) {
      frag.appendChild(h('div.pdiag-empty',
        'No diagnostics. The executor accepts this program — which is not the same as the program being a good idea.'));
    }
    const order = { error: 0, warning: 1, info: 2 };
    const sorted = list.slice().sort((a, b) =>
      (order[a.severity] ?? 3) - (order[b.severity] ?? 3) || (a.line - b.line));
    for (const d of sorted) {
      const kind = d.severity === 'error' ? 'err' : d.severity === 'warning' ? 'warn' : 'info';
      const btn = h('button.pdiag.' + kind, {
        type: 'button',
        role: 'listitem',
        title: d.line ? `jump to line ${d.line}` : 'no source line',
        onclick: () => editor.gotoLine(d.line || 1, d.col, d.end_col),
      },
        h('span.g', kind === 'err' ? '×' : kind === 'warn' ? '!' : 'i'),
        h('span.ln', d.line ? `ln ${d.line}` : '—'),
        h('span.msg',
          d.message,
          d.code ? h('span.code', ` [${d.code}]`) : null,
          d.source === 'stub' ? h('span.code', ' [stub]') : null),
      );
      if (d.remedy_text) {
        const fix = h('button.fix', {
          type: 'button',
          title: `Apply: ${d.remedy_text}`,
          onclick: (e) => { e.stopPropagation(); applyRemedy(d); },
        }, `apply  ${d.remedy_text}`);
        btn.appendChild(h('span', { style: { gridColumn: '3' } }, fix));
      }
      frag.appendChild(btn);
    }
    fill(diagList, frag);
  }

  function applyRemedy(d) {
    const lines = editor.value.split('\n');
    const at = Math.max(1, Math.min(lines.length, d.remedy_line || d.line || 1));
    if (d.remedy_kind === 'replace_line') {
      const ind = (/^\s*/.exec(lines[at - 1]) || [''])[0];
      lines[at - 1] = ind + d.remedy_text;
      editor.setValue(lines.join('\n'));
      editor.gotoLine(at);
    } else {
      lines.splice(at, 0, d.remedy_text);
      editor.setValue(lines.join('\n'));
      editor.gotoLine(at + 1);
    }
    editor.flush();
    toast('info', 'Applied to the program', d.remedy_text, 3500);
  }

  /* ====================================================== upload ======== */

  async function doUpload(file) {
    fill(uploadOut, h('div.caveat.info', h('span', `Uploading ${file.name} (${bytes(file.size)}) …`)));
    let res;
    try {
      res = await papi.upload(file);
    } catch (err) {
      fill(uploadOut, h('div.caveat.crit', h('span',
        h('b', 'Upload refused. '),
        err.message || String(err))));
      notifyError('POST /api/protocol/upload failed', err);
      return;
    }
    uploadData = res;
    paintMode();
    renderUpload(res);
  }

  const PROV_CLASS = {
    measured: 'measured', computed: 'computed',
    'user-supplied': 'supplied', 'flagged-default': 'flagged',
  };
  const PROV_WHY = {
    measured: 'Measured: an instrument or a published measured dataset produced this value directly.',
    computed: 'Computed: derived by the molecular layer from the structure you supplied. Reproducible; only as good as its inputs.',
    'user-supplied': 'User-supplied: you gave this number in the file or in the drug block. It overrides derivation and is not checked against a measurement.',
    'flagged-default': 'FLAGGED DEFAULT: nothing in your file determines this field, so a default stands in. It is shown here, before the run, precisely so it is not discovered inside a result.',
  };

  function provChip(kind) {
    const cls = PROV_CLASS[kind] || 'unknown';
    return h('span.chip.' + cls, { title: PROV_WHY[kind] || 'Provenance not stated.' }, kind);
  }

  function renderUpload(res) {
    const d = res.drug || {};
    const fields = d.fields || {};
    const flagged = Object.entries(fields).filter(([, f]) => f.provenance === 'flagged-default');

    const rows = Object.entries(fields).map(([k, f]) => h('tr' + (f.provenance === 'flagged-default' ? '.flagged' : ''),
      h('td', k),
      h('td.v' + (f.value === null || f.value === undefined ? '.null' : ''),
        f.value === null || f.value === undefined ? 'not resolved' : num(f.value, Math.abs(f.value) >= 100 ? 2 : 4),
        f.unit ? h('span.unit', f.unit) : null),
      h('td.p', provChip(f.provenance)),
    ));
    const notes = Object.entries(fields).map(([k, f]) => h('tr',
      h('td', { colspan: '3', class: 'note' }, h('b', k + ': '), f.source ? `${f.source} — ` : '', f.note || '')));

    fill(uploadOut,
      h('div.hr-lbl', 'Parsed preview — nothing has run'),
      h('dl.kv',
        h('dt', 'file'), h('dd', { title: res.file.name }, res.file.name),
        h('dt', 'size'), h('dd', bytes(res.file.size)),
        h('dt', 'format'), h('dd', res.file.format),
        h('dt', 'drug'), h('dd', d.name || NULL_GLYPH),
        h('dt', 'source'), h('dd', { title: d.source || '' }, d.source || NULL_GLYPH),
      ),
      flagged.length
        ? h('div.caveat', h('span',
          h('b', `${flagged.length} field${flagged.length === 1 ? '' : 's'} resolved to a FLAGGED DEFAULT: `),
          flagged.map(([k]) => k).join(', '),
          '. That is visible now, before the run, rather than inside a result. Override it in the drug block or accept it knowingly.'))
        : h('div.caveat.info', h('span', 'Every DrugParams field was either measured, computed from the structure, or supplied by you. No flagged defaults.')),
      h('div.tblwrap',
        h('table.tbl.pfields',
          h('thead', h('tr', h('th', 'field'), h('th.right', 'value'), h('th', 'provenance'))),
          h('tbody', rows, notes))),
      (res.warnings || []).map((w) => h('div.caveat', h('span', w.message || String(w)))),
      (res.caveats || []).map((c) => h('div.caveat.info', h('span', c))),
      h('div.btn-row', { style: { marginTop: 'var(--sp-4)' } },
        h('button.btn.sm', { onclick: () => insertDrugBlock(res), title: 'Write a `drug` block for this file into the editor' },
          'insert drug block')),
    );
  }

  function insertDrugBlock(res) {
    const d = res.drug || {};
    const name = (d.name || 'UPLOADED').toUpperCase().replace(/[^A-Z0-9_]/g, '_');
    const src = d.smiles ? `smiles "${d.smiles}"` : `file "uploads/${res.file.name}"`;
    const lines = [`drug ${name} from ${src}`];
    for (const [k, f] of Object.entries(d.fields || {})) {
      if (f.provenance === 'user-supplied' && f.value !== null) {
        lines.push(`  ${k} ${f.value}${f.unit && f.unit !== 'g/mol' && f.unit !== '1/s' ? ' ' + f.unit : ''}`);
      } else if (f.provenance === 'flagged-default') {
        lines.push(`  # ${k}: FLAGGED DEFAULT — give a value here or accept the executor's`);
      }
    }
    lines.push('end', '');
    editor.setValue(lines.join('\n') + editor.value);
    editor.flush();
    editor.gotoLine(1);
    toast('info', 'Drug block inserted', `${name} — ${lines.length - 2} field line(s)`, 4000);
  }

  /* ======================================================== plan ======== */

  function clearPlan() {
    planResult = null;
    planSource = '';
    fill(planOut, h('div.pdiag-empty', 'No plan yet. `plan` expands the sweep and states the run count and cost; it launches nothing.'));
    confirmBar.classList.add('hidden');
    clear(confirmBar);
  }

  function markPlanStale() {
    if (!planResult || planSource === editor.value) return;
    planResult = null;
    confirmBar.classList.add('hidden');
    const note = h('div.caveat', h('span',
      'The program changed after this plan was made. Re-plan before running — the run count and the cost below are stale.'));
    planOut.insertBefore(note, planOut.firstChild);
  }

  async function doPlan(thenConfirm) {
    editor.flush();
    fill(planOut, h('div.caveat.info', h('span', 'POST /api/protocol/plan …')));
    confirmBar.classList.add('hidden');
    let p;
    try {
      p = await papi.plan(editor.value);
    } catch (err) {
      fill(planOut, h('div.caveat.crit', h('span', h('b', 'Plan failed. '), err.message || String(err))));
      notifyError('POST /api/protocol/plan failed', err);
      return;
    }
    paintMode();
    planResult = p;
    planSource = editor.value;

    // Plan errors and refusals are diagnostics on their line, not a toast.
    const diags = []
      .concat((p.errors || []).map((d) => normDiag(d, 'error')))
      .concat((p.warnings || []).map((d) => normDiag(d, 'warning')));
    if (p.refusal) {
      diags.unshift(normDiag({
        line: p.refusal.line || 1,
        severity: 'error',
        message: p.refusal.message,
        code: p.refusal.code,
        remedy_text: p.refusal.remedy_text || p.refusal.remedy,
        remedy_kind: 'insert_line',
        remedy_line: p.refusal.remedy_line || p.refusal.line,
      }, 'error'));
    }
    if (diags.length) renderDiags(diags);

    renderPlan(p);
    if (p.ok && !p.refusal && thenConfirm) showConfirm(p);
  }

  function renderPlan(p) {
    const overCap = p.refusal && p.refusal.code === 'max_runs_exceeded';
    const nodes = [];

    nodes.push(h('div.pplan-nums',
      h('div.pplan-num',
        h('span.k', 'runs'),
        h('span.v' + (overCap ? '.over' : ''), int(p.run_count)),
        h('span.u', `of ${int(p.max_runs)} max_runs`)),
      h('div.pplan-num',
        h('span.k', 'est. cost'),
        h('span.v', dur(p.estimated_cost_s)),
        h('span.u', 'wall, serial')),
      h('div.pplan-num',
        h('span.k', 'axes'),
        h('span.v', int((p.axes || []).length)),
        h('span.u', (p.axes || []).map((a) => a.name).join(' × ') || 'none')),
    ));
    if (p.cost_note) nodes.push(h('div.caveat.info', h('span', p.cost_note)));

    if (!p.ok && !p.refusal) {
      nodes.push(h('div.caveat.crit', h('span',
        'The program does not validate, so there is nothing to plan. The diagnostics beside the script carry the reasons and the lines.')));
    }

    if (p.refusal) {
      nodes.push(h('div.prefusal',
        h('h4', `refused — ${p.refusal.code}`),
        h('div', p.refusal.message),
        p.refusal.remedy
          ? h('div', { style: { marginTop: '4px' } },
            'Remedy: add ', h('code', p.refusal.remedy), ' to the program, then plan again.',
            h('button.btn.sm', {
              style: { marginLeft: '6px' },
              onclick: () => applyRemedy({
                remedy_text: p.refusal.remedy_text || p.refusal.remedy,
                remedy_kind: 'insert_line',
                remedy_line: p.refusal.remedy_line || p.refusal.line,
              }),
            }, 'insert it'))
          : null,
        h('div', { style: { marginTop: '4px', color: 'var(--fg-2)' } },
          `It is also a diagnostic on line ${p.refusal.line} — it does not vanish while you read it.`),
      ));
    }

    for (const a of p.axes || []) {
      nodes.push(h('div.hr-lbl', `axis · ${a.name}${a.unit ? ' / ' + a.unit : ''} — ${a.values.length} values`));
      nodes.push(h('div.mono', {
        style: { fontSize: 'var(--fs-micro)', color: 'var(--fg-1)', overflowWrap: 'anywhere', padding: '0 0 3px' },
      }, a.values.map((v) => num(v, Math.abs(v) >= 100 ? 0 : 2)).join('  ')));
    }

    if ((p.assumptions || []).length) {
      nodes.push(h('div.hr-lbl', 'dose → PK conversion · stated assumptions'));
      nodes.push(h('ul.passump', p.assumptions.map((a) => h('li', a))));
    }

    if ((p.runs || []).length) {
      nodes.push(h('div.hr-lbl', `expanded run list — ${p.runs.length} runs`));
      const shown = p.runs.slice(0, 200);
      nodes.push(h('div.tblwrap', { style: { maxHeight: '210px' } },
        h('table.tbl.zebra',
          h('thead', h('tr',
            h('th', '#'), h('th', 'drug'), h('th', 'combination'),
            h('th.right', 'dose'), h('th.right', 'every'), h('th.right', 'dx / um'),
            h('th.right', 'est / s'))),
          h('tbody', shown.map((r) => h('tr',
            h('td.mono', String(r.index)),
            h('td', r.drug || NULL_GLYPH),
            h('td', { title: r.label }, r.label),
            h('td.right.num', fmtQty(r.params && r.params.dose)),
            h('td.right.num', fmtQty(r.params && r.params.every)),
            h('td.right.num', num(r.spec && r.spec.dx_um, 1)),
            h('td.right.num', num(r.est_s, 1)),
          ))))));
      if (p.runs.length > shown.length) {
        nodes.push(h('div.pdiag-empty', `${p.runs.length - shown.length} further runs not listed.`));
      }
    }

    for (const c of p.caveats || []) nodes.push(h('div.caveat', h('span', c)));

    fill(planOut, nodes);
  }

  function showConfirm(p) {
    clear(confirmBar);
    confirmBar.classList.remove('hidden');
    const go = h('button.btn.primary', {
      onclick: () => { confirmBar.classList.add('hidden'); doRun(p); },
    }, `Confirm and run ${p.run_count} run${p.run_count === 1 ? '' : 's'}`);
    confirmBar.appendChild(h('span.lbl', 'confirm'));
    confirmBar.appendChild(h('span.mono.strong',
      `${int(p.run_count)} runs · est ${dur(p.estimated_cost_s)}`));
    confirmBar.appendChild(h('span.dim',
      p.run_count > 1 ? 'each run is a full PDE solve' : 'one full PDE solve'));
    confirmBar.appendChild(h('span.grow'));
    confirmBar.appendChild(h('button.btn.sm', { onclick: () => confirmBar.classList.add('hidden') }, 'Cancel'));
    confirmBar.appendChild(go);
    go.focus();
  }

  /* ======================================================= execute ====== */

  async function doRun(p) {
    let res;
    try {
      res = await papi.start(editor.value, p.run_count);
    } catch (err) {
      const line = err.payload && err.payload.refusal ? err.payload.refusal.line : 1;
      renderDiags([normDiag({ line, severity: 'error', message: err.message, code: err.status ? `HTTP ${err.status}` : 'refused' }, 'error')]);
      fill(resultsOut, h('div.caveat.crit', h('span', h('b', 'Refused. '), err.message)));
      return;
    }
    sweepId = res.sweep_id;
    selectedIndex = -1;
    toast('ok', `Sweep ${sweepId} queued`, `${res.run_count} runs`, 4000);
    document.getElementById('pws-results').open = true;
    poll();
  }

  // Repainting a live grid every 900 ms would steal focus from a cell the
  // operator is reading, so a poll that changes nothing changes nothing.
  let lastSig = '';
  const sweepSignature = (sw) => sw.status + '|' + sw.runs
    .map((r) => `${r.status}:${(r.progress || 0).toFixed(2)}:${r.metrics ? 1 : 0}`).join(',');

  async function poll() {
    if (!sweepId) return;
    clearTimeout(pollTimer);
    try {
      sweepData = await papi.sweep(sweepId);
      const sig = sweepSignature(sweepData);
      if (sig !== lastSig) { lastSig = sig; renderResults(sweepData); }
    } catch (err) {
      fill(resultsOut, h('div.caveat.crit', h('span',
        h('b', `GET /api/protocol/sweeps/${sweepId} failed. `), err.message)));
    }
    const terminal = sweepData && ['done', 'failed', 'cancelled', 'partial'].includes(sweepData.status);
    if (!terminal) pollTimer = setTimeout(poll, 900);
  }

  async function cancelSweep() {
    if (!sweepId) return;
    try {
      sweepData = await papi.cancel(sweepId);
      renderResults(sweepData);
      toast('warn', `Sweep ${sweepId} cancelled`, 'Runs already in flight stop at their next stage boundary.', 8000);
    } catch (err) {
      notifyError(`DELETE /api/protocol/sweeps/${sweepId} failed`, err);
    }
  }

  /* ======================================================== results ===== */

  /**
   * Dimensions of the result grid. A swept parameter is a dimension; so is the
   * drug, when a `compare` produced more than one. Everything downstream — the
   * curve, the heatmap, the axis labels — is built from this one list.
   */
  function dimsOf(sw) {
    const axes = sw.axes || [];
    const drugs = [...new Set(sw.runs.map((r) => r.drug).filter(Boolean))];
    const dims = axes.map((a, i) => ({
      name: a.name, unit: a.unit, values: a.values, numeric: true,
      indexOf: (r) => (r.coords ? r.coords[i] : 0),
      label: (v) => `${num(v, Math.abs(v) >= 100 ? 0 : 2)}`,
    }));
    if (drugs.length > 1) {
      dims.push({
        name: 'drug', unit: null, values: drugs, numeric: false,
        indexOf: (r) => drugs.indexOf(r.drug),
        label: (v) => String(v),
      });
    }
    return dims;
  }

  function renderResults(sw) {
    const nodes = [];
    const dims = dimsOf(sw);
    const metrics = sw.metrics && sw.metrics.length
      ? sw.metrics
      : Object.keys(TPL_METRICS).map((k) => ({ key: k, ...TPL_METRICS[k] }));
    if (!metricKey || !metrics.some((m) => m.key === metricKey)) metricKey = metrics[0].key;
    const meta = metricMeta(metricKey, metrics);

    /* --- header: status, metric selector, cancel ---------------------- */
    const c = sw.counts || {};
    nodes.push(h('div.pws-head', { style: { height: '24px', borderBottom: '1px solid var(--line-1)' } },
      h('span.lbl', 'sweep'),
      h('span.mono.strong', sw.sweep_id),
      h('span.chip.st-' + (sw.status === 'partial' ? 'failed' : sw.status), sw.status),
      h('span.vrule', { 'aria-hidden': 'true' }),
      h('span.mono', `${int(c.done || 0)} done`),
      c.running ? h('span.mono.accent', `${int(c.running)} running`) : null,
      c.queued ? h('span.mono.dim', `${int(c.queued)} queued`) : null,
      c.failed ? h('span.mono.crit', `${int(c.failed)} failed`) : null,
      h('span.dim', `/ ${int(sw.run_count)}`),
      h('span.vrule', { 'aria-hidden': 'true' }),
      h('span.lbl', 'elapsed'),
      h('span.mono', dur(sw.elapsed_s)),
      h('span.grow'),
      h('div.seg', { role: 'group', 'aria-label': 'metric shown in the grid' },
        metrics.map((m) => h('button', {
          type: 'button',
          'aria-pressed': String(m.key === metricKey),
          title: `${m.label}${m.unit ? ' / ' + m.unit : ''}`,
          onclick: () => { metricKey = m.key; renderResults(sweepData); },
        }, m.key))),
      ['done', 'failed', 'cancelled', 'partial'].includes(sw.status)
        ? null
        : h('button.btn.danger.sm', { onclick: cancelSweep, title: 'Cancel every queued and running job in this sweep' }, 'Cancel sweep'),
    ));

    /* --- standing caveats sit WITH the numbers ------------------------ */
    for (const cav of sw.caveats || []) nodes.push(h('div.caveat', h('span', cav)));

    /* --- the grid ------------------------------------------------------ */
    const vals = sw.runs.map((r) => (r.metrics ? r.metrics[metricKey] : null)).filter((v) => Number.isFinite(v));
    const lo = vals.length ? Math.min(...vals) : 0;
    const hi = vals.length ? Math.max(...vals) : 1;

    if (!dims.length) {
      nodes.push(h('div.sect-body', renderStrip(sw, dims, meta, lo, hi)));
    } else if (dims.length === 1) {
      if (dims[0].numeric) nodes.push(h('div.sect-body', { style: { paddingBottom: '0' } }, curveHost(sw, dims[0], meta)));
      nodes.push(h('div.sect-body', { style: { paddingTop: '2px' } }, renderMatrix(sw, dims, meta, lo, hi)));
    } else {
      nodes.push(h('div.sect-body', renderMatrix(sw, dims, meta, lo, hi)));
    }

    nodes.push(legendFor(meta, lo, hi, vals.length));

    /* --- assumptions + selected run ----------------------------------- */
    if ((sw.assumptions || []).length) {
      nodes.push(h('details.sect',
        h('summary', 'Dose → PK conversion', h('span.hint', `${sw.assumptions.length} assumptions`)),
        h('ul.passump', sw.assumptions.map((a) => h('li', a)))));
    }
    nodes.push(renderDetail(sw, metrics));

    fill(resultsOut, nodes);
    scheduleCurveRedraw();
  }

  function cellFor(sw, run, meta, lo, hi) {
    const v = run && run.metrics ? run.metrics[meta.key] : null;
    const status = run ? run.status : 'queued';
    const prov = run && status === 'done' ? cellProvenance(meta.key, run) : { kind: 'unknown', why: 'not finished' };
    const el = h('button.pcell.st-' + status + '.pv-' + prov.kind, {
      type: 'button',
      'aria-selected': String(run && run.index === selectedIndex),
      onclick: () => { if (run) selectRun(sw, run.index); },
      title: run
        ? [
          run.label,
          `run_id ${run.run_id || 'not submitted'}`,
          `status ${status}`,
          v !== null && v !== undefined ? `${meta.label} = ${fmtMetric(meta, v)} ${meta.unit}` : 'no value yet',
          `provenance: ${prov.kind}`,
          prov.why,
          run.error ? 'ERROR: ' + run.error : '',
        ].filter(Boolean).join('\n')
        : 'no run at this combination',
    });
    if (!run) { el.disabled = true; el.textContent = NULL_GLYPH; return el; }
    if (status === 'done' && Number.isFinite(v)) {
      const t = hi > lo ? (v - lo) / (hi - lo) : 0.5;
      const bg = sample(RAMP, t);
      el.style.background = bg;
      el.style.color = inkOn(bg);
      el.appendChild(h('span', fmtMetric(meta, v)));
    } else if (status === 'failed') {
      el.appendChild(h('span.st', 'failed'));
    } else if (status === 'running') {
      el.appendChild(h('span.st', 'running'));
      el.appendChild(h('i.bar', { style: { width: Math.round((run.progress || 0) * 100) + '%' } }));
    } else if (status === 'cancelled') {
      el.appendChild(h('span.st', 'cancel'));
    } else {
      el.appendChild(h('span.st', 'queued'));
    }
    return el;
  }

  function renderStrip(sw, dims, meta, lo, hi) {
    return h('div.btn-row', sw.runs.map((r) => cellFor(sw, r, meta, lo, hi)));
  }

  function renderMatrix(sw, dims, meta, lo, hi) {
    const col = dims[0];
    const rest = dims.slice(1);
    const rowCombos = rest.reduce((acc, d) => {
      const out = [];
      for (const base of acc) for (let i = 0; i < d.values.length; i++) out.push(base.concat(i));
      return out;
    }, [[]]);

    const byKey = new Map();
    for (const r of sw.runs) byKey.set(dims.map((d) => d.indexOf(r)).join('|'), r);

    const head = h('tr',
      h('th.corner', rest.length ? rest.map((d) => d.name).join(' × ') : ''),
      col.values.map((v, i) => h('th', { title: `${col.name} = ${col.label(v)}${col.unit ? ' ' + col.unit : ''}` },
        col.label(v))));

    const body = rowCombos.map((combo) => h('tr',
      h('td.rowhead', {
        title: rest.map((d, i) => `${d.name} = ${d.label(d.values[combo[i]])}${d.unit ? ' ' + d.unit : ''}`).join('\n'),
      }, rest.length
        ? rest.map((d, i) => d.label(d.values[combo[i]])).join(' · ')
        : meta.key),
      col.values.map((_, ci) => h('td', cellFor(sw, byKey.get([ci].concat(combo).join('|')), meta, lo, hi))),
    ));

    return h('div',
      h('div.plegend',
        h('span.lbl', 'columns'),
        h('span.mono.strong', col.name),
        col.unit ? h('span.unit', col.unit) : null,
        rest.length ? h('span.lbl', { style: { marginLeft: '8px' } }, 'rows') : null,
        rest.length ? h('span.mono.strong', rest.map((d) => d.name).join(' × ')) : null,
        rest.length && rest.some((d) => d.unit) ? h('span.unit', rest.map((d) => d.unit || '—').join(' × ')) : null),
      h('div.pgrid-wrap',
        h('table.pgrid', h('thead', head), h('tbody', body))));
  }

  function legendFor(meta, lo, hi, n) {
    const ramp = h('div.ramp');
    ramp.style.backgroundImage = gradientCss(RAMP, 32);
    return h('div.plegend',
      h('span.lbl', meta.label),
      h('span.unit', meta.unit || ''),
      h('span.tick', n ? fmtMetric(meta, lo) : NULL_GLYPH),
      ramp,
      h('span.tick', n ? fmtMetric(meta, hi) : NULL_GLYPH),
      h('span.dim', `· ${RAMP}, scaled across the ${n} finished run${n === 1 ? '' : 's'}`),
      h('span.grow'),
      h('span.lbl', 'cell edge = provenance'),
      h('span', h('i.sw.pv-measured'), ' measured'),
      h('span', h('i.sw.pv-computed'), ' computed'),
      h('span', h('i.sw.pv-assumed'), ' assumed'),
    );
  }

  /* --- one-axis response curve -------------------------------------- */

  let curveState = null;
  function curveHost(sw, dim, meta) {
    const host = h('div', { style: { minHeight: '186px' } });
    curveState = { host, sw, dim, meta };
    drawCurve();
    return host;
  }

  const scheduleCurveRedraw = debounce(() => drawCurve(), 120);
  window.addEventListener('resize', scheduleCurveRedraw);

  function drawCurve() {
    if (!curveState || !curveState.host.isConnected) return;
    const { host, sw, dim, meta } = curveState;
    const W = Math.max(260, host.clientWidth || 420);
    const H = 186;
    const M = { l: 52, r: 10, t: 10, b: 32 };
    const iw = W - M.l - M.r;
    const ih = H - M.t - M.b;

    const pts = sw.runs
      .map((r) => ({ r, x: dim.values[dim.indexOf(r)], y: r.metrics ? r.metrics[meta.key] : null }))
      .filter((p) => Number.isFinite(p.x))
      .sort((a, b) => a.x - b.x);
    const xs = dim.values;
    const x0 = Math.min(...xs);
    const x1 = Math.max(...xs);
    const ys = pts.map((p) => p.y).filter(Number.isFinite);
    const y0 = ys.length ? Math.min(...ys) : 0;
    const y1 = ys.length ? Math.max(...ys) : 1;
    const yPad = (y1 - y0) * 0.12 || Math.abs(y1) * 0.12 || 1;
    const ylo = y0 - yPad;
    const yhi = y1 + yPad;

    const px = (v) => M.l + (x1 > x0 ? (v - x0) / (x1 - x0) : 0.5) * iw;
    const py = (v) => M.t + ih - (yhi > ylo ? (v - ylo) / (yhi - ylo) : 0.5) * ih;

    const kids = [];
    // y grid + ticks
    const yStep = nice125((yhi - ylo) / 4) || 1;
    for (let v = Math.ceil(ylo / yStep) * yStep; v <= yhi; v += yStep) {
      kids.push(svg('line', { class: 'gridline', x1: M.l, x2: M.l + iw, y1: py(v), y2: py(v) }));
      kids.push(svg('text', { x: M.l - 5, y: py(v) + 3, 'text-anchor': 'end' }, fmtMetric(meta, v)));
    }
    // x ticks: one per swept value, thinned to fit
    const stride = Math.max(1, Math.ceil(xs.length / Math.max(2, Math.floor(iw / 46))));
    xs.forEach((v, i) => {
      if (i % stride) return;
      kids.push(svg('text', { x: px(v), y: H - M.b + 12, 'text-anchor': 'middle' },
        num(v, Math.abs(v) >= 100 ? 0 : 2)));
    });
    kids.push(svg('line', { class: 'ax', x1: M.l, x2: M.l + iw, y1: M.t + ih, y2: M.t + ih }));
    kids.push(svg('line', { class: 'ax', x1: M.l, x2: M.l, y1: M.t, y2: M.t + ih }));
    kids.push(svg('text', { class: 'axlabel', x: M.l + iw / 2, y: H - 3, 'text-anchor': 'middle' },
      `${dim.name}${dim.unit ? ' / ' + dim.unit : ''}`));
    kids.push(svg('text', {
      class: 'axlabel', x: 10, y: M.t + ih / 2, 'text-anchor': 'middle',
      transform: `rotate(-90 10 ${M.t + ih / 2})`,
    }, `${meta.label}${meta.unit ? ' / ' + meta.unit : ''}`));

    const done = pts.filter((p) => Number.isFinite(p.y));
    if (done.length > 1) {
      kids.push(svg('path', { class: 'path', d: done.map((p, i) => `${i ? 'L' : 'M'}${px(p.x).toFixed(1)},${py(p.y).toFixed(1)}`).join(' ') }));
    }
    for (const p of pts) {
      const sel = p.r.index === selectedIndex;
      const cls = 'pt' + (sel ? ' sel' : '') + (p.r.status === 'failed' ? ' failed' : '')
        + (Number.isFinite(p.y) ? '' : ' pending');
      const c = svg('circle', {
        class: cls,
        cx: px(p.x).toFixed(1),
        cy: Number.isFinite(p.y) ? py(p.y).toFixed(1) : (M.t + ih).toFixed(1),
        r: sel ? 4.5 : 3.2,
        tabindex: '0',
        role: 'button',
      });
      c.appendChild(svg('title', null,
        `${dim.name} = ${num(p.x, 2)}${dim.unit ? ' ' + dim.unit : ''}\n`
        + `${meta.label} = ${fmtMetric(meta, p.y)} ${meta.unit}\nstatus ${p.r.status}`));
      c.addEventListener('click', () => selectRun(sw, p.r.index));
      c.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); selectRun(sw, p.r.index); } });
      kids.push(c);
    }

    const el = svg('svg', {
      class: 'pcurve', viewBox: `0 0 ${W} ${H}`, width: '100%', height: H,
      role: 'img',
      'aria-label': `${meta.label} against ${dim.name}, ${done.length} of ${pts.length} runs finished`,
    }, ...kids);
    fill(host, el);
  }

  /* --- selected run -------------------------------------------------- */

  function selectRun(sw, index) {
    selectedIndex = index;
    renderResults(sw);
    const run = sw.runs.find((r) => r.index === index);
    if (run && run.run_id && run.status === 'done') openInViewer(run.run_id);
  }

  function renderDetail(sw, metrics) {
    const run = sw.runs.find((r) => r.index === selectedIndex);
    if (!run) {
      return h('div.pdetail', h('div.bd', h('span.dim',
        'Select a cell to load that run into the viewer behind this workspace and see its parameters, metrics, provenance and RunSpec here.')));
    }
    const expects = (sw.expects || []).map((x) => {
      const v = run.metrics ? run.metrics[x.metric] : null;
      if (!Number.isFinite(v)) return { ...x, pass: null };
      const pass = x.op === '>' ? v > x.value : x.op === '<' ? v < x.value
        : x.op === '>=' ? v >= x.value : x.op === '<=' ? v <= x.value
          : x.op === '==' ? v === x.value : v !== x.value;
      return { ...x, pass, actual: v };
    });

    return h('div.pdetail',
      h('div.hd',
        h('span.lbl', 'run'),
        h('span.mono.strong', run.label),
        h('span.chip.st-' + run.status, run.status),
        h('span.grow'),
        h('span.mono.dim', run.run_id || 'not submitted'),
        run.run_id && run.status === 'done'
          ? h('button.btn.sm', {
            onclick: () => { openInViewer(run.run_id); hide(); },
            title: 'Close this workspace and show the tissue for this run',
          }, 'show in viewer')
          : null),
      h('div.bd',
        run.error ? h('div.caveat.crit', h('span', h('b', 'Run failed: '), run.error)) : null,
        h('div.hr-lbl', 'metrics'),
        h('table.tbl',
          h('thead', h('tr', h('th', 'metric'), h('th.right', 'value'), h('th', 'provenance'))),
          h('tbody', metrics.map((m) => {
            const v = run.metrics ? run.metrics[m.key] : null;
            const p = run.status === 'done' ? cellProvenance(m.key, run) : { kind: 'unknown', why: 'the run has not finished' };
            return h('tr',
              h('td', { title: m.label }, m.key),
              h('td.right.num', fmtMetric(m, v), h('span.unit', m.unit || '')),
              h('td', chip(p.kind === 'user-supplied' || p.kind === 'flagged-default' ? 'assumed' : p.kind, p.why)));
          }))),
        expects.length ? h('div.hr-lbl', 'expectations — recorded, never enforced') : null,
        expects.length ? h('div', expects.map((x) => h('div.frow', { style: { fontSize: 'var(--fs-xs)' } },
          // Status vocabulary, NOT provenance vocabulary. A passing assertion
          // is not a measurement, and spending the green `measured` chip on it
          // devalues the one chip that is supposed to mean "an instrument
          // produced this".
          h('span.chip.' + (x.pass === null ? 'st-queued' : x.pass ? 'st-done' : 'st-failed'),
            x.pass === null ? 'pending' : x.pass ? 'pass' : 'fail'),
          h('span.mono', `expect ${x.metric} ${x.op} ${x.value}`),
          h('span.dim', x.actual !== undefined
            ? `actual ${fmtRaw(metricMeta(x.metric, metrics), x.actual)}`
            : 'the run has not produced this metric yet')))) : null,
        h('div.hr-lbl', 'parameters'),
        h('dl.kv',
          h('dt', 'drug'), h('dd', run.drug || NULL_GLYPH),
          h('dt', 'dose'), h('dd', fmtQty(run.params && run.params.dose)),
          h('dt', 'every'), h('dd', fmtQty(run.params && run.params.every)),
          h('dt', 'repeat'), h('dd', fmtQty(run.params && run.params.repeat)),
          h('dt', 'cycles'), h('dd', fmtQty(run.params && run.params.cycles)),
          h('dt', 'dx'), h('dd', fmtQty(run.params && run.params.dx)),
          h('dt', 'route'), h('dd', (run.params && run.params.route) || NULL_GLYPH),
          h('dt', 'runtime'), h('dd', dur(run.runtime_s))),
        h('div.hr-lbl', 'submitted RunSpec'),
        h('pre.mono', {
          style: {
            margin: '0', padding: '4px 6px', maxHeight: '150px', overflow: 'auto',
            fontSize: 'var(--fs-micro)', lineHeight: '1.4', color: 'var(--fg-2)',
            background: 'var(--bg-void)', border: '1px solid var(--line-1)',
          },
        }, JSON.stringify(run.spec, null, 2)),
      ));
  }

  /* --------------------------------------------------------- viewer ---- */

  function openInViewer(runId) {
    window.dispatchEvent(new CustomEvent('tmesim:open-run', { detail: { run_id: runId } }));
    const t = window.__tmesim;
    if (t && typeof t.openRun === 'function') { t.openRun(runId); return; }
    toast('warn', 'Viewer not reachable',
      `The shell did not expose an openRun hook, so run ${runId} could not be loaded into the viewer from here. It is in the run history.`, 7000);
  }

  /* ----------------------------------------------------------- mode ---- */

  function paintMode() {
    const el = document.getElementById('pws-mode');
    if (el) el.textContent = papi.protoApi.stub ? 'STUB' : 'live';
    if (papi.protoApi.stub) {
      stubNote.classList.remove('hidden');
      fill(stubNote,
        h('b', 'STUB EXECUTOR ·'),
        h('span',
          `${papi.protoApi.reason}. The program is parsed and expanded in this browser by js/protocol/stub.js against PROTOCOL_LANGUAGE.md, and each planned run is submitted as a v1 RunSpec through POST /api/runs. `,
          'The dose→PK arithmetic and the cost estimate are the stub\'s, not the executor\'s — read them as a rehearsal of the interface, not as a result.'));
    } else {
      stubNote.classList.add('hidden');
    }
  }

  /* ------------------------------------------------------ visibility --- */

  let visible = false;
  function show() {
    if (visible) return;
    visible = true;
    root.hidden = false;
    paintMode();
    editor.repaint();
    editor.focus();
    scheduleCurveRedraw();
    const b = document.getElementById('btn-protocol');
    if (b) b.setAttribute('aria-pressed', 'true');
  }
  function hide() {
    if (!visible) return;
    visible = false;
    root.hidden = true;
    const b = document.getElementById('btn-protocol');
    if (b) b.setAttribute('aria-pressed', 'false');
  }
  const toggle = () => (visible ? hide() : show());

  // Keys pressed inside the workspace must not also drive the viewer behind it.
  root.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && !e.defaultPrevented) { hide(); }
    e.stopPropagation();
  });

  clearPlan();
  fill(resultsOut, h('div.pdiag-empty',
    'No sweep yet. `Run…` plans first, shows the run count and the estimated cost, and waits for an explicit confirmation before queueing anything.'));
  renderDiags([]);

  return {
    el: root, show, hide, toggle, editor,
    get visible() { return visible; },
    validateNow: () => runValidate(editor.value, true),
    // Exposed for the same reason window.__tmesim is: a headless driver has to
    // be able to reach the plan step without synthesising a click.
    plan: () => doPlan(),
    planAndConfirm: () => doPlan(true),
    get lastPlan() { return planResult; },
    get lastSweep() { return sweepData; },
  };
}

/* ============================================================== mount ==== */

function mount() {
  if (document.getElementById('protocol-ws')) return;
  const ws = createProtocolWorkspace();
  document.body.appendChild(ws.el);

  let btn = document.getElementById('btn-protocol');
  if (!btn) {
    btn = h('button#btn-protocol.btn.ghost.sm', {
      title: 'Protocol editor and parameter sweeps (P)',
      'aria-pressed': 'false',
    }, 'protocol');
    const group = document.querySelector('.topbar-group.last');
    if (group) group.insertBefore(btn, group.firstChild);
  }
  btn.addEventListener('click', () => ws.toggle());

  window.addEventListener('keydown', (e) => {
    const t = e.target;
    const typing = t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA'
      || t.tagName === 'SELECT' || t.isContentEditable);
    if (typing || e.ctrlKey || e.metaKey || e.altKey) return;
    if (e.key === 'p' || e.key === 'P') { e.preventDefault(); ws.toggle(); }
  });

  window.__tmesimProtocol = ws;
  // Kick off the first validation so the pane is never blank-but-silent.
  ws.validateNow();
}

if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', mount);
else mount();
