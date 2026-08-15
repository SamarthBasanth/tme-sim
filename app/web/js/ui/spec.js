// spec.js — RUN CONTROL: the RunSpec form and the RunSpec source editor,
// kept in two-way sync, both validated live against POST /api/validate.
//
// Sync contract between the two views:
//   form edit  -> spec object -> re-serialise editor buffer (unless the
//                 operator is mid-keystroke in the editor)
//   editor edit-> parse; on success adopt the object and repaint the form;
//                 on parse failure the form is LEFT ALONE and the parse error
//                 is shown as a diagnostic. We never silently discard text.
//
// Validation is remote, not local: the backend owns the rules. The local
// mock implements the same rules only so the panel is exercisable offline.

import { h, clear, fill, append, debounce } from '../util/dom.js';
import * as api from '../api.js';
import { toast, notifyError } from './toast.js';
import { railGroup, railNav } from './rail.js';
import { createCompoundPanel } from './compound.js';
import { createTissuePanel, NOT_VALIDATED, NOT_RUNNABLE } from './tissue.js';
import { dur, num } from '../util/fmt.js';

const MODES = {
  sensitivity_mode: ['uniform', 'determinant'],
  vascular_mode: ['discrete', 'kde'],
  oxygen_mode: ['pde', 'vessel_distance', 'markers'],
  state_mode: ['oxygen', 'vessel_distance', 'markers'],
};

// Fallback for a RunSpec that predates state_mode (a restored run, or a hand-typed
// buffer). It is the pipeline default and is bit-identical to the pre-option
// behaviour, so assuming it here changes nothing.
const STATE_MODE_DEFAULT = 'oxygen';

const NO_DRUG = '— none (reconstruction only)';

/** Normalise whatever the backend puts in errors[]/warnings[]. */
function normDiag(d) {
  if (typeof d === 'string') {
    const m = /^([A-Za-z_][\w.]*)\s*:\s*(.+)$/.exec(d);
    return m ? { path: m[1], message: m[2] } : { path: '', message: d };
  }
  if (!d || typeof d !== 'object') return { path: '', message: String(d) };
  const path = d.path || (Array.isArray(d.loc) ? d.loc.filter((x) => x !== 'body').join('.') : '') || d.field || '';
  const message = d.message || d.msg || d.detail || JSON.stringify(d);
  return { path, message };
}

export function createSpecPanel({ store, onSubmit, onDrugsChanged, onSectionsChanged }) {
  let editorFocused = false;
  let syncing = false;
  const fieldEls = new Map();

  // Bring your own compound. Owned by ./compound.js; mounted in the Define
  // group below. When the operator adds or removes one, the drug catalogue has
  // changed for every panel that shows drugs, so the shell is asked to refetch
  // rather than this panel patching its own <select> and leaving the reports
  // panel stale.
  const compoundPanel = createCompoundPanel({
    onChanged: () => { if (onDrugsChanged) onDrugsChanged(); },
  });

  // Bring your own tissue. Owned by ./tissue.js, mounted beside the compound
  // compiler. Ingesting or removing a section changes the section CATALOGUE,
  // not just this picker (the projects panel and the protocol validator both
  // hold their own copy), so the shell is asked to refetch rather than this
  // panel patching its own <select> and leaving the others stale.
  const tissuePanel = createTissuePanel({
    onRegistered: (sectionId) => {
      if (onSectionsChanged) onSectionsChanged(sectionId || null);
    },
  });

  /* ------------------------------------------------------------- form -- */

  // Last list from GET /api/sections, kept so the caveat line can be repainted
  // on selection change without another round trip.
  let sections_ = [];
  const sectionCaveats = h('div.section-caveats', { style: 'display:none' });
  const sectionSel = h('select.inp', {
    id: 'f-section',
    onchange: () => { fromForm(); paintSectionCaveats(); },
  });
  const drugSel = h('select.inp', { id: 'f-drug', onchange: fromForm });
  const dxInp = h('input.inp.num', { id: 'f-dx', type: 'number', step: '0.5', min: '2', max: '100', oninput: fromForm });
  const totalInp = h('input.inp.num', { id: 'f-total', type: 'number', step: '10', min: '1', oninput: fromForm });
  const mgmtChk = h('input', { type: 'checkbox', id: 'f-mgmt-on', onchange: fromForm });
  const mgmtInp = h('input.inp.num', { id: 'f-mgmt', type: 'number', step: '0.1', min: '-2', max: '12', oninput: fromForm });
  const pkPeak = h('input.inp.num', { id: 'f-pk-peak', type: 'number', step: '0.1', min: '0', oninput: fromForm });
  const pkHalf = h('input.inp.num', { id: 'f-pk-half', type: 'number', step: '600', min: '1', oninput: fromForm });
  const pkN = h('input.inp.num', { id: 'f-pk-n', type: 'number', step: '1', min: '1', max: '30', oninput: fromForm });
  const pkInt = h('input.inp.num', { id: 'f-pk-int', type: 'number', step: '3600', min: '1', oninput: fromForm });

  const segs = {};
  for (const [key, opts] of Object.entries(MODES)) {
    const buttons = opts.map((o) => h('button', {
      type: 'button',
      'aria-pressed': 'false',
      dataset: { val: o },
      onclick: () => {
        for (const b of segs[key].children) b.setAttribute('aria-pressed', String(b.dataset.val === o));
        fromForm();
      },
    }, o));
    segs[key] = h('div.seg', { role: 'group', 'aria-label': key.replace(/_/g, ' ') }, buttons);
  }

  const unit = (t) => h('span.unit', t);
  const lbl = (text, forId, title) => h('label.lbl', { for: forId, title: title || '' }, text);

  function row(label, forId, control, unitText, title) {
    const holder = unitText
      ? h('div.with-unit', control, unit(unitText))
      : control;
    // id + role="alert" so the error is BOTH reachable via aria-describedby
    // (a screen reader visiting the field reads it alongside aria-invalid)
    // AND announced on its own the moment validation fills it in, for a
    // field the operator isn't currently on.
    const errId = forId + '-err';
    const errEl = h('div.field-err.hidden', { id: errId, role: 'alert' });
    // The control we mark aria-invalid on must be the focusable field itself,
    // not a wrapper: a div with aria-invalid tells assistive tech nothing.
    // The target must be the field the error PATH names, not merely the first
    // focusable descendant. In the MGMT row the first descendant is the "set"
    // checkbox and the field is the number input beside it, so a backend error
    // on path `mgmt` marked the checkbox invalid and left the offending value
    // unmarked. Prefer an exact id match; fall back to first-focusable.
    const target = /^(INPUT|SELECT|TEXTAREA)$/.test(control.tagName)
      ? control
      : (control.querySelector('#' + CSS.escape(forId))
        || control.querySelector('input, select, textarea')
        || control);
    target.setAttribute('aria-describedby', errId);
    fieldEls.set(forId, { control: target, errEl });
    return [lbl(label, forId, title), holder, errEl];
  }

  const formEl = h('div.fgrid',
    row('Section', 'f-section', h('div', sectionSel, sectionCaveats), null,
      'GET /api/sections'),
    row('Grid dx', 'f-dx', dxInp, 'um', 'reconstruction grid pitch'),
    row('Drug', 'f-drug', drugSel, null, 'GET /api/drugs; null = reconstruction only'),
    [lbl('Sensitivity', null), segs.sensitivity_mode, h('div.field-err.hidden')],
    [lbl('Vascular', null), segs.vascular_mode, h('div.field-err.hidden')],
    [lbl('Oxygen', null), segs.oxygen_mode, h('div.field-err.hidden')],
    [lbl('Cell state', null), segs.state_mode, h('div.field-err.hidden')],
    row('MGMT',
      'f-mgmt',
      h('div.frow',
        h('label.chk', { title: 'null when unset — the model then assumes no MGMT information' },
          mgmtChk, h('span.lbl', 'set')),
        h('div.with-unit.grow', mgmtInp, unit('log2 CPM'))),
      null, 'per-patient scalar; spatial MGMT was measured and found to be noise'),
    h('div.span2.hr-lbl', 'Pharmacokinetics'),
    row('Plasma peak', 'f-pk-peak', pkPeak, 'uM'),
    row('Half life', 'f-pk-half', pkHalf, 's'),
    row('Doses', 'f-pk-n', pkN, 'n'),
    row('Interval', 'f-pk-int', pkInt, 's'),
    h('div.span2.hr-lbl', 'Integration'),
    row('Total time', 'f-total', totalInp, 's', 'simulated wall time of the treatment course'),
  );

  /* ----------------------------------------------------------- editor -- */

  const gutter = h('div.editor-gutter', { 'aria-hidden': 'true' });
  const ta = h('textarea.editor-ta', {
    id: 'spec-editor',
    spellcheck: 'false',
    autocomplete: 'off',
    autocapitalize: 'off',
    autocorrect: 'off',
    wrap: 'off',
    'aria-label': 'RunSpec JSON source',
    oninput: () => { fromEditor(); },
    onfocus: () => { editorFocused = true; },
    onblur: () => { editorFocused = false; },
    onscroll: () => { gutter.scrollTop = ta.scrollTop; },
    onkeydown: (e) => {
      if (e.key === 'Tab') {
        e.preventDefault();
        const s = ta.selectionStart, en = ta.selectionEnd;
        ta.value = ta.value.slice(0, s) + '  ' + ta.value.slice(en);
        ta.selectionStart = ta.selectionEnd = s + 2;
        fromEditor();
        return;
      }
      // Tab is claimed above for indentation (same convention as the TPL
      // editor, js/editor/editor.js), which would otherwise be a hard
      // keyboard trap: nothing else here restores normal Tab behaviour.
      // Escape blurs the field, mirroring the TPL editor's escape hatch, so
      // the operator always has a keyboard-only way out.
      if (e.key === 'Escape') {
        e.stopPropagation();
        ta.blur();
      }
    },
  });

  const parseState = h('span', 'json ok');
  const editorEl = h('div.editor',
    h('div.editor-head',
      h('span', 'RunSpec'),
      h('span.grow'),
      parseState,
    ),
    h('div.editor-scroll', gutter, ta),
    h('div.editor-foot',
      h('button.btn.ghost.sm', { onclick: format, title: 'Re-indent the buffer' }, 'format'),
      h('button.btn.ghost.sm', { onclick: copySpec, title: 'Copy JSON to clipboard' }, 'copy'),
      h('button.btn.ghost.sm', { onclick: resetDefaults, title: 'GET /api/defaults' }, 'defaults'),
      h('span.grow'),
      h('span.dim', { id: 'spec-bytes' }),
    ),
  );

  const diagList = h('div.diag-list');

  /* ------------------------------------------------------------ submit -- */

  const submitBtn = h('button.btn.primary', {
    id: 'btn-submit-run',
    onclick: () => onSubmit(store.get().spec),
    title: 'POST /api/runs',
  }, 'Submit run');
  const validateBtn = h('button.btn', { onclick: () => runValidate(true), title: 'POST /api/validate' }, 'Validate');
  // Validation state is a STATUS, so it wears the st-* vocabulary. It used to
  // borrow measured/assumed/unknown, which put a green "measured" chip on the
  // fact that a JSON object passed a schema check.
  const okBadge = h('span.chip.st-queued', 'unvalidated');

  const submitRow = h('div.frow', { style: { padding: 'var(--sp-4) var(--sp-5)' } },
    submitBtn, validateBtn, h('span.grow'), okBadge);

  /* ------------------------------------------------------- job monitor -- */

  const jobBar = h('i', { style: { width: '0%' } });
  const jobProgress = h('div.progress', jobBar);
  const jobStage = h('span.mono', '—');
  const jobPct = h('span.mono', '—');
  const jobElapsed = h('span.mono', '—');
  const jobId = h('span.mono.dim', '—');
  const jobStatus = h('span.chip.st-queued', 'idle');
  const cancelBtn = h('button.btn.danger.sm', {
    disabled: true,
    onclick: () => cancelTracked(),
    title: 'DELETE /api/runs/{id}',
  }, 'Cancel');
  const logEl = h('pre.mono', {
    style: {
      margin: '0', padding: '4px 6px', maxHeight: '92px', overflow: 'auto',
      fontSize: 'var(--fs-micro)', lineHeight: '1.4', color: 'var(--fg-2)',
      background: 'var(--bg-void)', border: '1px solid var(--line-1)',
    },
    'aria-label': 'Job log tail',
  }, '');

  const jobEl = h('div.stack.tight',
    h('div.frow', jobId, jobStatus, h('span.grow'), cancelBtn),
    jobProgress,
    h('div.frow',
      h('span.lbl', 'stage'), jobStage, h('span.grow'),
      jobPct, h('span.unit', '%'),
      h('span.dim', '·'), jobElapsed),
    logEl,
  );

  /* -------------------------------------------------------------- root -- */

  // A scrolling body with a fixed commit row pinned to the bottom.
  //
  // Two groups, matching the READ rail opposite: DEFINE is everything that
  // mutates the RunSpec, COMMIT is what happens to it afterwards. The commit
  // row itself is outside the scroll on purpose — it is the one control whose
  // reachability must not depend on scroll position.
  // See docs/UI_ARCHITECTURE.md section 2.2.
  const body = h('div.rail-body',
    railGroup('define', 'Define', 'RunSpec v1',
      h('details.sect', { open: true },
        h('summary', 'Specimen, therapy, model', h('span.hint', 'form')),
        h('div.sect-body', formEl)),
      // Both BYO panels sit directly under the form that picks a section and a
      // drug, because that is the moment an operator discovers the list does
      // not contain their specimen or their molecule. Tissue first: it is the
      // outer of the two choices — a compound is compared ON a section. Both
      // stay collapsed by default; each is a deliberate act, not part of the
      // ordinary path.
      tissuePanel.el,
      compoundPanel.el,
      h('details.sect', { open: true },
        h('summary', 'Source', h('span.hint', 'JSON · two-way sync')),
        h('div.sect-body', { style: { padding: '0' } }, editorEl))),

    railGroup('job', 'Commit', 'POST /api/runs',
      h('details.sect', { open: true, id: 'sect-job' },
        h('summary', 'Active job', h('span.hint', '/status poll')),
        h('div.sect-body', jobEl))),
  );

  const nav = railNav([
    { id: 'define', label: 'Define' },
    { id: 'job', label: 'Commit' },
  ], body);

  const root = h('div.rail.rail-left',
    h('div.panel-head', h('h2', 'Setup'), h('span.grow'),
      h('span.hint', { title: 'This rail defines what will be simulated. It never interprets a result — that is the Read rail on the right.' },
        'what will run')),
    nav.el,
    body,
    // The validator read-out sits WITH the commit control, outside the scroll,
    // not at the tail of the Source section inside it.
    //
    // Inside the scroll it was the last thing in a 1364 px column shown
    // through an 815 px viewport, so at rest the operator saw one amber box
    // sawn off mid-sentence by the bottom edge of the rail — text that is cut,
    // not text that is scrolled, which reads as a rendering bug rather than as
    // "there is more below". It is also the one thing you are supposed to have
    // read BEFORE pressing the button directly beneath it, and a control whose
    // reachability must not depend on scroll position (the reason `submitRow`
    // is already out here) should not sit above a warning whose visibility
    // does. It costs nothing on a clean spec: `.diag-list` is empty, and an
    // empty one has no border and no height.
    h('div', { style: { borderTop: '1px solid var(--line-2)', background: 'var(--bg-2)' } },
      diagList, submitRow),
  );

  /* --------------------------------------------------------- behaviour -- */

  // What a section CANNOT support, shown where the section is chosen rather
  // than discovered afterwards in a run's caveat panel. `caveats` comes from
  // GET /api/sections (see api.py::_with_caveats): platform-level refusals for
  // every section, plus the stored capability verdicts of an ingested one.
  function sectionTitle(s) {
    const lines = [];
    if (s.panel) lines.push(s.panel);
    if (s.ingested) lines.push(NOT_VALIDATED);
    if (s.runnable === false) lines.push(NOT_RUNNABLE);
    if (!s.available && s.reason) lines.push(s.reason);
    for (const c of s.caveats || []) {
      lines.push(`[${c.status}] ${c.capability}: ${c.reason}`);
    }
    return lines.join('\n\n');
  }

  function setSections(sections) {
    clear(sectionSel);
    for (const s of sections) {
      // A blocking caveat is not the same as unavailable: the section still
      // runs, but something it could otherwise report is refused on it.
      const blocking = s.n_caveats_blocking || 0;
      const mark = s.available && blocking ? `  [${blocking} refused]` : '';
      // An ingested section is marked IN THE OPTION, not only in the panel that
      // ingested it: this picker is where it gets chosen, and by then that panel
      // may be closed — or the section may have been ingested by the CLI in
      // another session entirely. The distinction the mark carries is the one
      // that matters: a built-in specimen is what this repo's evidence was
      // measured from, and this one is not.
      const own = s.ingested ? '  [ingested — not validated]' : '';
      // `runnable: false` is a MEX bundle: it ingested, it is inspectable, and
      // the reconstruction pipeline has no route for it (ingest_registry.py).
      // This picker chooses what will be SUBMITTED, so such a section is not
      // offered here at all — the alternative is a run that is guaranteed to
      // fail after the operator has waited for it.
      const noRun = s.runnable === false ? '  [no run route — inspect only]' : '';
      sectionSel.appendChild(h('option', {
        value: s.id,
        disabled: !s.available || s.runnable === false,
        title: sectionTitle(s),
      }, `${s.label}${s.available ? '' : '  [unavailable]'}${own}${noRun}${mark}`));
    }
    sections_ = sections;
    paintSectionCaveats();
  }

  // The picker's own caveat line, repainted whenever the selection changes.
  // Silent on a built-in section with nothing refused, so it stays meaningful.
  function paintSectionCaveats() {
    if (!sectionCaveats) return;
    clear(sectionCaveats);
    const s = (sections_ || []).find((x) => x.id === sectionSel.value);
    const blocking = (s && s.caveats ? s.caveats : []).filter(
      (c) => c.status === 'NO');
    // An ingested section says so HERE, at the moment of selection, in full and
    // in the same type as every other refusal in this block. It leads, because
    // it qualifies everything under it: the capability caveats below say what
    // this specimen cannot support, and this line says that nothing it CAN
    // support has been checked against a measurement.
    const ingested = !!(s && s.ingested);
    if (!blocking.length && !ingested) { sectionCaveats.style.display = 'none'; return; }
    sectionCaveats.style.display = '';
    if (ingested) {
      sectionCaveats.appendChild(h('div.hint', { title: NOT_VALIDATED },
        h('b', 'This specimen was ingested here, not shipped with Codon. '),
        NOT_VALIDATED.replace(/^INGESTED — NOT VALIDATED\.\s*/, '')));
    }
    if (!blocking.length) return;
    sectionCaveats.appendChild(h('div.hint',
      `This section cannot support ${blocking.length} readout` +
      `${blocking.length === 1 ? '' : 's'}:`));
    for (const c of blocking) {
      sectionCaveats.appendChild(
        h('div.hint', { title: c.reason }, `• ${c.capability} — REFUSED`));
    }
  }

  function setDrugs(drugs) {
    clear(drugSel);
    drugSel.appendChild(h('option', { value: '' }, NO_DRUG));
    for (const d of drugs) {
      // A self-compiled compound is marked IN THE OPTION, not only in the panel
      // that produced it: the picker is where it gets chosen, and by then the
      // compile screen may be closed. The suffix says how many of its fields are
      // predictions, because that is the fact that distinguishes it from a
      // curated row — a curated row's IC50 is a GDSC measurement, this one's is
      // the operator's own number or a flagged placeholder.
      const label = `${d.name}  ·  IC50 ${num(d.ic50_uM, 1)} uM  ·  BBB ${num(d.bbb_penetrance, 2)}`;
      const suffix = d.is_user_compound
        ? (d.reason
          ? '   [YOURS — does not compile]'
          : `   [YOURS · ${d.predicted_fields ?? '?'}/${d.total_fields ?? '?'} fields predicted]`)
        : '';
      drugSel.appendChild(h('option', {
        value: d.name,
        disabled: !!(d.is_user_compound && d.reason),
        title: d.is_user_compound
          ? (d.reason
            || 'Compiled from a SMILES you supplied. Every transport parameter '
             + 'except molecular weight and TPSA is a prediction. Its absolute '
             + 'kill, AUC and penetration depth are NOT defensible. Compare it '
             + 'against other compounds on this section and read the ordering — '
             + 'noting that diffusivity comes from a two-route model scoring '
             + 'R^2 = 0.718 blind against 71 measured compounds, which still '
             + 'recovers only about a quarter of the measured spread from one '
             + 'cell line, so similar-mass analogues are separated weakly.')
          : '',
      }, label + suffix));
    }
  }

  function toForm(spec) {
    if (!spec) return;
    syncing = true;
    if (sectionSel.value !== spec.section) sectionSel.value = spec.section ?? '';
    drugSel.value = spec.drug ?? '';
    dxInp.value = spec.dx_um ?? '';
    totalInp.value = spec.total_time_s ?? '';
    const hasMg = spec.mgmt !== null && spec.mgmt !== undefined;
    mgmtChk.checked = hasMg;
    mgmtInp.disabled = !hasMg;
    mgmtInp.value = hasMg ? spec.mgmt : '';
    const pk = spec.pk || {};
    pkPeak.value = pk.plasma_peak ?? '';
    pkHalf.value = pk.half_life_s ?? '';
    pkN.value = pk.n_doses ?? '';
    pkInt.value = pk.dose_interval_s ?? '';
    for (const [key, opts] of Object.entries(MODES)) {
      for (const b of segs[key].children) {
        b.setAttribute('aria-pressed', String(b.dataset.val === spec[key]));
      }
      if (!opts.includes(spec[key])) {
        for (const b of segs[key].children) b.setAttribute('aria-pressed', 'false');
      }
    }
    syncing = false;
    paintSectionCaveats();
  }

  function readForm() {
    const prev = store.get().spec || {};
    const pick = (key) => {
      const btn = Array.from(segs[key].children).find((b) => b.getAttribute('aria-pressed') === 'true');
      return btn ? btn.dataset.val : prev[key];
    };
    const numOr = (el, fallback) => {
      const v = parseFloat(el.value);
      return Number.isFinite(v) ? v : fallback;
    };
    return {
      section: sectionSel.value || prev.section,
      dx_um: numOr(dxInp, prev.dx_um),
      drug: drugSel.value === '' ? null : drugSel.value,
      sensitivity_mode: pick('sensitivity_mode'),
      vascular_mode: pick('vascular_mode'),
      oxygen_mode: pick('oxygen_mode'),
      state_mode: pick('state_mode') ?? STATE_MODE_DEFAULT,
      mgmt: mgmtChk.checked ? numOr(mgmtInp, 0) : null,
      pk: {
        plasma_peak: numOr(pkPeak, 1),
        half_life_s: numOr(pkHalf, 14400),
        n_doses: Math.round(numOr(pkN, 1)),
        dose_interval_s: numOr(pkInt, 86400),
      },
      total_time_s: numOr(totalInp, 600),
    };
  }

  function serialise(spec) {
    // Stable key order so the editor buffer does not shuffle under the cursor.
    const o = {
      section: spec.section,
      dx_um: spec.dx_um,
      drug: spec.drug ?? null,
      sensitivity_mode: spec.sensitivity_mode,
      vascular_mode: spec.vascular_mode,
      oxygen_mode: spec.oxygen_mode,
      state_mode: spec.state_mode ?? STATE_MODE_DEFAULT,
      mgmt: spec.mgmt ?? null,
      pk: spec.pk,
      total_time_s: spec.total_time_s,
    };
    return JSON.stringify(o, null, 2);
  }

  function fromForm() {
    if (syncing) return;
    mgmtInp.disabled = !mgmtChk.checked;
    const spec = readForm();
    const text = serialise(spec);
    store.set({ spec, specText: text, specTextValid: true });
    if (!editorFocused) setEditorText(text);
    paintGutter();
    scheduleValidate();
  }

  function fromEditor() {
    const text = ta.value;
    store.set({ specText: text });
    let parsed = null, err = null;
    try { parsed = JSON.parse(text); } catch (e) { err = e; }
    if (err) {
      store.set({ specTextValid: false });
      parseState.textContent = 'json invalid';
      parseState.className = 'crit';
      renderDiags([{ path: '', message: 'JSON parse error: ' + err.message, line: guessErrLine(text, err) }], []);
      okBadge.className = 'chip st-failed';
      okBadge.textContent = 'unparsed';
      submitBtn.disabled = true;
      paintGutter();
      return;
    }
    parseState.textContent = 'json ok';
    parseState.className = 'ok';
    store.set({ spec: parsed, specTextValid: true });
    toForm(parsed);
    paintGutter();
    scheduleValidate();
  }

  function setEditorText(text) {
    if (ta.value === text) return;
    const s = ta.selectionStart, e = ta.selectionEnd;
    ta.value = text;
    try { ta.setSelectionRange(Math.min(s, text.length), Math.min(e, text.length)); } catch { /* detached */ }
  }

  function guessErrLine(text, err) {
    const m = /position\s+(\d+)/.exec(err.message);
    if (!m) return 0;
    return text.slice(0, Number(m[1])).split('\n').length;
  }

  let errLines = new Set();
  function paintGutter() {
    const lines = ta.value.split('\n').length;
    const frag = document.createDocumentFragment();
    for (let i = 1; i <= lines; i++) {
      const s = h('span', String(i));
      if (errLines.has(i)) s.className = 'gl-err';
      frag.appendChild(s);
      frag.appendChild(document.createTextNode('\n'));
    }
    fill(gutter, frag);
    gutter.scrollTop = ta.scrollTop;
    const b = document.getElementById('spec-bytes');
    if (b) b.textContent = `${lines} ln · ${new Blob([ta.value]).size} B`;
  }

  function lineOfPath(path) {
    if (!path) return 0;
    const key = path.split('.').pop();
    const lines = ta.value.split('\n');
    for (let i = 0; i < lines.length; i++) {
      if (new RegExp(`"${key}"\\s*:`).test(lines[i])) return i + 1;
    }
    return 0;
  }

  function focusLine(n) {
    if (!n) { ta.focus(); return; }
    const lines = ta.value.split('\n');
    let pos = 0;
    for (let i = 0; i < n - 1 && i < lines.length; i++) pos += lines[i].length + 1;
    ta.focus();
    ta.setSelectionRange(pos, pos + (lines[n - 1] || '').length);
    const lh = 15;
    ta.scrollTop = Math.max(0, (n - 4) * lh);
    gutter.scrollTop = ta.scrollTop;
  }

  function renderDiags(errors, warnings) {
    errLines = new Set();
    // The diagnostics get their own bounded scroll (css .diag-scroll) rather
    // than running on to the end of the rail. Four multi-sentence validator
    // warnings are 475 px of prose at the tail of a 1364 px column in an
    // 815 px viewport: the last one was sawn in half by the bottom of the
    // scroller with nothing to say it continued, which reads as broken text
    // rather than as scrolled text. The count header says how many there are
    // before you scroll, so bounding the box hides no fact — same contract the
    // protocol panel's .pdiags already runs under.
    //
    // The box also gets an explicit "more below" footer. It would not need one
    // if a scrollbar were visible, but base.css asks for an always-visible
    // track and the platform overrides it with an overlay bar that occupies no
    // width and paints nothing at rest (measured: offsetWidth === clientWidth
    // on an element with scrollHeight 483 in 116). Without a cue, bounded and
    // clipped look identical — which is the exact complaint this box exists to
    // answer, so the cue is drawn rather than delegated.
    const scroll = h('div.diag-scroll');
    const frag = document.createDocumentFragment();
    const add = (d, kind) => {
      const line = d.line !== undefined ? d.line : lineOfPath(d.path);
      if (kind === 'err' && line) errLines.add(line);
      // Three fragments, not one string. The spec path and the line number are
      // machine tokens and keep the mono face; the message is a sentence of
      // English and takes the sans face, the same way the identical text does
      // in the Caveats rail. Nothing is added, removed or reworded — the same
      // characters in the same order, told apart by what they are.
      frag.appendChild(h('button.diag.' + kind, {
        type: 'button',
        onclick: () => focusLine(line),
        title: line ? `jump to line ${line}` : 'no source line',
      }, h('span',
        d.path ? h('span.d-path', d.path + ': ') : null,
        d.message,
        line ? h('span.d-ln', `  (ln ${line})`) : null)));
    };
    for (const e of errors) add(e, 'err');
    for (const w of warnings) add(w, 'warn');
    scroll.appendChild(frag);
    clear(diagList);
    if (errors.length || warnings.length) {
      diagList.appendChild(h('div.diag-head',
        errors.length ? h('span.n-err', `${errors.length} error${errors.length === 1 ? '' : 's'}`) : null,
        warnings.length ? h('span.n-warn', `${warnings.length} warning${warnings.length === 1 ? '' : 's'}`) : null,
        h('span.grow'),
        h('span', 'click one to jump to its line')));
      diagList.appendChild(scroll);

      const more = h('button.diag-more', {
        type: 'button', hidden: true,
        title: 'Scroll the diagnostics. Nothing is hidden — every warning is in the list above.',
        onclick: () => scroll.scrollBy({ top: scroll.clientHeight - 20, behavior: 'smooth' }),
      }, h('span.g', '▾'), h('span.t', ''));
      diagList.appendChild(more);

      const syncMore = () => {
        const over = scroll.scrollHeight - scroll.clientHeight;
        const left = over - scroll.scrollTop;
        more.hidden = over <= 2;
        more.querySelector('.t').textContent = left > 2 ? 'more below — scroll' : 'end of list';
        more.classList.toggle('at-end', left <= 2);
      };
      scroll.addEventListener('scroll', syncMore, { passive: true });
      requestAnimationFrame(syncMore);
    }

    // mirror onto the form controls
    for (const [, { control, errEl }] of fieldEls) {
      control.removeAttribute('aria-invalid');
      if (control.classList) control.classList.remove('invalid');
      errEl.classList.add('hidden');
      errEl.textContent = '';
    }
    const idFor = {
      section: 'f-section', dx_um: 'f-dx', drug: 'f-drug', mgmt: 'f-mgmt',
      total_time_s: 'f-total',
      'pk.plasma_peak': 'f-pk-peak', 'pk.half_life_s': 'f-pk-half',
      'pk.n_doses': 'f-pk-n', 'pk.dose_interval_s': 'f-pk-int',
    };
    for (const e of errors) {
      const id = idFor[e.path];
      if (!id) continue;
      const t = fieldEls.get(id);
      if (!t) continue;
      t.control.setAttribute('aria-invalid', 'true');
      t.control.classList.add('invalid');
      t.errEl.textContent = e.message;
      t.errEl.classList.remove('hidden');
    }
    paintGutter();
  }

  const scheduleValidate = debounce(() => runValidate(false), 340);

  async function runValidate(explicit) {
    const st = store.get();
    if (!st.spec || !st.specTextValid) return;
    store.set({ validating: true });
    okBadge.className = 'chip st-running';
    okBadge.textContent = 'validating';
    try {
      const res = await api.postValidate(st.spec);
      const errors = (res.errors || []).map(normDiag);
      const warnings = (res.warnings || []).map(normDiag);
      const ok = res.ok !== undefined ? !!res.ok : errors.length === 0;
      store.set({ validation: { ok, errors, warnings }, validating: false });
      renderDiags(errors, warnings);
      okBadge.className = 'chip st-' + (ok ? 'done' : 'failed');
      // `valid · 5 warn` used to be one string inside one green chip, so the
      // ok hue — the loudest "nothing to see here" signal the palette has —
      // was carrying a warning count that on this section includes a hard
      // stop. Same words, same chip, same order; the half that is a warning
      // count is now in the warning colour. The verdict is unchanged: the
      // spec IS valid, and it DOES carry five warnings.
      clear(okBadge);
      if (!ok) {
        okBadge.textContent = `${errors.length} error${errors.length === 1 ? '' : 's'}`;
      } else if (warnings.length) {
        append(okBadge, ['valid', h('span.sep', '·'),
          h('span.n-warn', `${warnings.length} warn`)]);
      } else {
        okBadge.textContent = 'valid';
      }
      submitBtn.disabled = !ok;
    } catch (err) {
      store.set({ validating: false });
      okBadge.className = 'chip st-failed';
      okBadge.textContent = 'validate failed';
      submitBtn.disabled = false;   // do not block the operator on our own outage
      if (explicit) notifyError('POST /api/validate failed', err);
    }
  }

  function format() {
    try {
      const o = JSON.parse(ta.value);
      setEditorText(serialise(o));
      fromEditor();
    } catch (e) { notifyError('Cannot format: JSON does not parse', e); }
  }

  async function copySpec() {
    try { await navigator.clipboard.writeText(ta.value); } catch { /* clipboard blocked */ }
  }

  async function resetDefaults() {
    try {
      const d = await api.getDefaults();
      loadSpec(d);
    } catch (err) { notifyError('GET /api/defaults failed', err); }
  }

  function loadSpec(spec) {
    store.set({ spec, specText: serialise(spec), specTextValid: true });
    setEditorText(serialise(spec));
    toForm(spec);
    paintGutter();
    runValidate(false);
  }

  /* ------------------------------------------------------ job monitor -- */

  let pollTimer = 0;
  let tracked = null;
  let onJobDone = null;
  // Cancellation is COOPERATIVE on the backend: DELETE sets a flag and a
  // running job only stops at the next stage boundary. The operator must be
  // told that, otherwise a click that visibly does nothing reads as a bug and
  // they hit it again (or worse, assume the run stopped when it has not).
  let cancelPending = false;
  let cancelNote = '';

  function track(runId, doneCb) {
    tracked = runId;
    onJobDone = doneCb;
    cancelPending = false;
    cancelNote = '';
    jobId.textContent = runId;
    cancelBtn.disabled = false;
    cancelBtn.textContent = 'Cancel';
    document.getElementById('sect-job').open = true;
    clearTimeout(pollTimer);
    poll();
  }

  async function poll() {
    if (!tracked) return;
    try {
      const s = await api.getStatus(tracked);
      store.set({ status: s, trackedRunId: tracked });
      const pct = Math.round((s.progress || 0) * 100);
      jobBar.style.width = pct + '%';
      jobPct.textContent = String(pct).padStart(3, ' ');
      jobStage.textContent = s.stage || '—';
      jobElapsed.textContent = dur(s.elapsed_s);
      const terminal = s.status === 'done' || s.status === 'failed' || s.status === 'cancelled';
      if (cancelPending && !terminal) {
        jobStatus.className = 'chip st-cancelled';
        jobStatus.textContent = `cancelling · still ${s.status}`;
        jobStatus.title = cancelNote;
      } else {
        jobStatus.className = 'chip st-' + (s.status || 'queued');
        jobStatus.textContent = s.status || '—';
        jobStatus.title = '';
      }
      jobProgress.classList.toggle('done', s.status === 'done');
      jobProgress.classList.toggle('failed', s.status === 'failed');
      jobProgress.classList.toggle('indet', s.status === 'queued' && !cancelPending);
      logEl.textContent = [
        ...(s.log_tail || []),
        cancelPending && !terminal ? '[ui] cancel requested — ' + (cancelNote || 'waiting for the next stage boundary') : null,
      ].filter(Boolean).join('\n');
      logEl.scrollTop = logEl.scrollHeight;

      if (s.status === 'done' || s.status === 'failed' || s.status === 'cancelled') {
        cancelBtn.disabled = true;
        const id = tracked;
        tracked = null;
        if (onJobDone) onJobDone(id, s.status);
        return;
      }
    } catch (err) {
      jobStatus.className = 'chip st-failed';
      jobStatus.textContent = 'poll error';
      logEl.textContent = String(err && err.message);
    }
    pollTimer = setTimeout(poll, 900);
  }

  async function cancelTracked() {
    if (!tracked) return;
    const id = tracked;
    cancelBtn.disabled = true;
    cancelBtn.textContent = 'Cancelling…';
    try {
      const res = await api.cancelRun(id);
      cancelPending = true;
      // The backend explains WHY a running job does not stop instantly. Show
      // its own words rather than inventing a reassurance of our own.
      cancelNote = (res && (res.reason || res.detail || res.message)) || '';
      const stillRunning = res && res.status && res.status !== 'cancelled';
      toast(stillRunning ? 'warn' : 'info',
        stillRunning ? `Cancel requested for ${id} — not stopped yet` : `Run ${id} cancelled`,
        cancelNote, stillRunning ? 12000 : 4000);
    } catch (err) {
      notifyError(`DELETE /api/runs/${id} failed`, err);
      cancelBtn.disabled = false;
      cancelBtn.textContent = 'Cancel';
    }
  }

  return {
    el: root,
    setSections, setDrugs, loadSpec, track, validate: () => runValidate(true),
    /** Repaint the operator's own compound list (GET /api/compounds). */
    refreshCompounds: () => compoundPanel.refresh(),
    /** Repaint the ingested-section list (GET /api/sections, `ingested` rows). */
    refreshTissue: () => tissuePanel.refresh(),
    get trackedId() { return tracked; },
    focusEditor: () => ta.focus(),
    setSubmitting(on) {
      submitBtn.disabled = on;
      submitBtn.textContent = on ? 'Submitting…' : 'Submit run';
    },
  };
}
