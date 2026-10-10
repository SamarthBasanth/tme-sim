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
import { railGroup, railNav, sectTitle } from './rail.js';
import { createCompoundPanel } from './compound.js';
import {
  createTissuePanel, NOT_VALIDATED, NOT_RUNNABLE, sectionScaleRow,
} from './tissue.js';
import { plasmaPeakNote } from './doseinput.js';
import { createPreviewPanel } from './preview.js';
import { dur, num, int } from '../util/fmt.js';

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

/* ===================================================== the three questions ===
 *
 * WHAT THIS FORM ASKED FOR, BEFORE 2026-09-20. Fourteen controls — twenty tab
 * stops — in one flat grid, every one of them arriving pre-filled from
 * GET /api/defaults and every one of them presented as a decision:
 *
 *   section · dx · drug · sensitivity_mode · vascular_mode · oxygen_mode ·
 *   state_mode · mgmt set · mgmt · plasma_peak · half_life_s · n_doses ·
 *   dose_interval_s · total_time_s
 *
 * Nine of the fourteen have a defensible default that a first-time operator has
 * no basis to move and no way to evaluate: the four mode segments select
 * pipeline routes, MGMT is a per-patient scalar that is null unless you have
 * one, and three of the four PK fields are the shape of the plasma curve rather
 * than its height. Presenting nine defaults as nine blanks is what makes this
 * screen read as a form to be filled rather than a question to be answered.
 *
 * WHAT FOLDS, AND WHAT THAT COSTS. Nothing is removed and nothing moves to
 * another surface; the nine go behind ONE `<details>` whose summary PRINTS
 * EVERY ONE OF THEIR VALUES, so the fold hides the widgets and never the
 * state. It opens itself when a validator error lands on a field inside it, and
 * its summary counts how many of them differ from the stock defaults the
 * session booted with. docs/UI_ARCHITECTURE.md §2.2 argues against splitting
 * this form into four `<details>` — "four header strips and buys nothing" —
 * and that argument is about four; this is one, and it buys back six rows and
 * ten tab stops on the first screen a stranger sees.
 *
 * WHAT STAYS VISIBLE is the set that changes the answer: which specimen, which
 * compound, how much of it, for how long, and on what mesh — the last of which
 * stays because it is the one field whose cost is quadratic and the only one
 * the server refuses a run over.
 */
const QUESTIONS = [
  ['1', 'Which specimen?'],
  ['2', 'Which therapy?'],
  ['3', 'At what resolution?'],
];

/** The ids inside the folded block, for the auto-open rule in renderDiags(). */
const FOLDED_FIELDS = new Set(['f-mgmt', 'f-pk-half', 'f-pk-n', 'f-pk-int']);

/** A validator diagnostic's `path` -> the form control that field is edited in.
 *
 * Module scope since 2026-09-20: it was a literal inside renderDiags(), and the
 * submit gate needs the same mapping to answer "press here to go to the thing
 * that is wrong". Two copies of this object would disagree the first time a
 * field moved. */
const FIELD_OF_PATH = {
  section: 'f-section', dx_um: 'f-dx', drug: 'f-drug', mgmt: 'f-mgmt',
  total_time_s: 'f-total',
  'pk.plasma_peak': 'f-pk-peak', 'pk.half_life_s': 'f-pk-half',
  'pk.n_doses': 'f-pk-n', 'pk.dose_interval_s': 'f-pk-int',
};

/* ========================================================== grid budget ===
 *
 * WHAT THIS EXISTS FOR. `dx_um` is the only field in this form whose cost is
 * quadratic, it is the field a first-time operator reaches for first, and
 * until now nothing in the editor said what a given value would ask the
 * server to allocate. The server refuses above a ceiling; the operator found
 * that out by pressing Submit.
 *
 * THE FORMULA IS THE SERVER'S, COPIED. app/server/api.py::validate_spec:
 *
 *     if (dx is not None and section_row
 *             and section_row.get("extent_um") and dx > 0):
 *         ex, ey = section_row["extent_um"]
 *         n_cells_grid = (ex / dx) * (ey / dx)
 *         if n_cells_grid > _MAX_GRID_CELLS:   -> an ERROR, the run is refused
 *
 * Deliberately NOT rounded to whole columns and rows. `js/mock.js` line 461
 * rounds — `Math.round(EX/dx) * Math.round(EY/dx)` — and the two disagree by a
 * fraction of a percent. This estimate has to predict a REFUSAL, so it mirrors
 * the code that refuses, not the code that reconstructs.
 *
 * THE DRIFT, STATED. This is a hand-copy of a Python expression into
 * JavaScript. Nothing enforces that the two stay equal: there is no JS test
 * runner in this repo, no endpoint that returns the count, and no shared
 * source for the expression. If api.py's expression changes — a different
 * rounding, a third dimension, a per-source correction — this function keeps
 * predicting the old number until a human reads both files. The blast radius
 * is bounded: the estimate is advisory, POST /api/validate remains the only
 * thing that decides, and its verdict is already rendered below the form.
 */
function gridCells(extentUm, dxUm) {
  if (!Array.isArray(extentUm) || extentUm.length < 2) return null;
  const ex = Number(extentUm[0]);
  const ey = Number(extentUm[1]);
  const dx = Number(dxUm);
  if (!Number.isFinite(ex) || !Number.isFinite(ey)) return null;
  if (!Number.isFinite(dx) || dx <= 0) return null;
  return (ex / dx) * (ey / dx);
}

/** Python's `{x:g}` closely enough to PREVIEW a sentence the server will write.
 *
 * Six significant figures, trailing zeros dropped — so 12.0 prints `12`, the
 * way the refusal does. It decides nothing: the refusal itself is the
 * server's, and this only predicts how it will read. Python switches to
 * exponent form outside roughly [1e-5, 1e6) and this does not; `dx_um` is a
 * micrometre pitch and never lands there in practice, and if it ever does the
 * preview reads oddly while the real refusal is unaffected.
 */
function pyG(x) {
  if (!Number.isFinite(x)) return String(x);
  return String(Number(x.toPrecision(6)));
}

/** The refusal the server would write, reproduced from api.py's own f-string.
 *
 * Reproduced rather than paraphrased so that the sentence the operator reads
 * before submitting is the sentence they would read after. Source, verbatim:
 *
 *     field_b = n_cells_grid * 8
 *     field = (f"{field_b / 1e9:,.1f} GB" if field_b >= 1e9
 *              else f"{field_b / 1e6:,.1f} MB")
 *     tuned = _MAX_GRID_CELLS != _DEFAULT_MAX_GRID_CELLS
 *     errors.append(
 *         f"dx_um = {dx:g} um on this section gives a grid of "
 *         f"~{n_cells_grid:,.0f} voxels, above the "
 *         f"{_MAX_GRID_CELLS:,} ceiling"
 *         + (" this server is configured with "
 *            "(TMESIM_MAX_GRID_CELLS)" if tuned else "")
 *         + f" — one float64 field alone would be {field}, and a "
 *         "solve allocates several plus the cell population. "
 *         "Use a coarser dx_um.")
 *
 * `tuned` is a fact only the server holds (it compares its ceiling against its
 * own default), so this client can only reproduce that clause when a server
 * payload has told it so. When it has not, the clause is omitted and the
 * preview is marked as a preview.
 */
function refusalText(dx, n, ceiling, tuned) {
  const bytes = n * 8;
  const field = bytes >= 1e9
    ? `${num(bytes / 1e9, 1)} GB`
    : `${num(bytes / 1e6, 1)} MB`;
  return `dx_um = ${pyG(dx)} um on this section gives a grid of ~${int(n)} `
    + `voxels, above the ${int(ceiling)} ceiling`
    + (tuned ? ' this server is configured with (TMESIM_MAX_GRID_CELLS)' : '')
    + ` — one float64 field alone would be ${field}, and a solve allocates `
    + 'several plus the cell population. Use a coarser dx_um.';
}

/** The field name the server side is expected to grow. See `ceiling_` below. */
const CEILING_KEY = 'max_grid_cells';

/** Pull a grid-cell ceiling out of ANY server payload, or return null.
 *
 * Three shapes are accepted because the field does not exist yet and the
 * session adding it has not committed to one: a top-level `max_grid_cells`,
 * or the same key under a `limits` / `server_limits` object. Anything absent,
 * non-numeric or non-positive reads as "this server did not say", which is a
 * state the panel renders rather than papers over.
 */
function readCeiling(o) {
  if (!o || typeof o !== 'object') return null;
  const candidates = [
    o[CEILING_KEY],
    o.limits && o.limits[CEILING_KEY],
    o.server_limits && o.server_limits[CEILING_KEY],
  ];
  for (const c of candidates) {
    const v = Number(c);
    if (Number.isFinite(v) && v > 0) return v;
  }
  return null;
}

/** Recover the ceiling from a refusal the server has already written.
 *
 * Last resort, and strictly an improvement on knowing nothing: once the server
 * has refused, its own sentence carries the number. Coupled to the wording in
 * api.py — if that is reworded the match simply fails and the panel goes back
 * to saying the ceiling is unknown, which is what it was saying anyway. Never
 * preferred over a payload field.
 */
const CEILING_FROM_REFUSAL = /above the ([\d,]+) ceiling/;
const CEILING_TUNED_MARK = 'TMESIM_MAX_GRID_CELLS';

/** The ONLY runtime this repository records for a Studio run, with its spec.
 *
 * docs/DEPLOYMENT.md, on the container acceptance run: "A real end-to-end run
 * — reconstruction plus transport plus metrics, submitted over the published
 * port and polled to `done` — completes in-container on 40,887 cells. On the
 * stock `GET /api/defaults` spec (`xenium_319`, temozolomide, `dx 12 um`, 6 h
 * window, `plasma_peak 50 uM`) that is **`runtime_s` 47.13 cold**, on a 4-CPU
 * VM-backed daemon."
 *
 * The same paragraph attaches two conditions and both are honoured here: it
 * is not a benchmark, and "a runtime was quoted with no RunSpec attached" is
 * named as the defect — "the window, the mesh and the compound each move it.
 * Quote the spec or don't quote the time." So the spec travels with the
 * number, and the number is shown as NOT APPLYING the moment the form moves
 * off that spec. NOTHING here scales it with dx: this repo records no measured
 * dx-to-runtime relation for a Studio run, and inventing one is forbidden.
 */
const RECORDED_RUN = {
  seconds: 47.13,
  where: 'docs/DEPLOYMENT.md',
  host: 'a 4-CPU VM-backed daemon, cold',
  spec: {
    section: 'xenium_319',
    drug: 'temozolomide',
    dx_um: 12,
    total_time_s: 21600,
    plasma_peak: 50,
  },
};

/** Does the form currently hold the exact spec RECORDED_RUN was measured on? */
function matchesRecordedRun(spec) {
  if (!spec) return false;
  const r = RECORDED_RUN.spec;
  const pk = spec.pk || {};
  return spec.section === r.section
    && spec.drug === r.drug
    && Number(spec.dx_um) === r.dx_um
    && Number(spec.total_time_s) === r.total_time_s
    && Number(pk.plasma_peak) === r.plasma_peak;
}

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

  // The coarse-solve preview. Owned by ./preview.js, mounted in the Commit
  // group below. It is handed the SPEC AS THE FORM HOLDS IT, unedited, so what
  // is previewed is what Submit would send — the server substitutes its own
  // coarse dx_um and refuses any spec whose pair it has not measured, and no
  // copy of those rules lives in the browser.
  const previewPanel = createPreviewPanel({ getSpec: () => store.get().spec });

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

  /* ------------------------------------------------------ grid budget -- */

  // The server-reported grid-cell ceiling, in grid cells. DELIBERATELY NOT A
  // CONSTANT IN THIS FILE.
  //
  // app/server/api.py holds it as `_MAX_GRID_CELLS`, which is
  // `_DEFAULT_MAX_GRID_CELLS = 40_000_000` unless the operator sets the
  // environment variable `TMESIM_MAX_GRID_CELLS` — and a constrained host is
  // supposed to set it, because the 40 M default is the "no request should
  // ever ask for this" number and not the "this machine can afford it" one.
  // A number typed here would therefore be wrong on exactly the deployments
  // that need it most, and wrong silently.
  //
  // `GET /api/limits` exposes it, added 2026-08-30 alongside the Account
  // panel. It is a route of its OWN rather than a field on /api/defaults, for
  // the reason the previous version of this comment gave: the /api/defaults
  // body is adopted verbatim as the RunSpec and posted straight back to
  // /api/validate, which warns on unrecognised fields, so a ceiling published
  // there would make every boot report an unrecognised-field warning.
  //
  // `api.getLimits()` is `noFallback`: it never answers from js/mock.js. A
  // fabricated ceiling is worse than no ceiling, because it would make this
  // panel confidently green on a spec the real server refuses. When the call
  // fails, this stays null and the panel says the ceiling is unknown.
  //
  // adoptCeilingFromDiags() remains the second source: once the server has
  // actually refused a spec, its own refusal text carries the real number, and
  // that is the one source that cannot be out of date.
  let ceiling_ = null;
  let ceilingSrc_ = '';
  let ceilingTuned_ = false;

  const gbChip = h('span.chip.st-queued', {
    title: 'How this spec stands against the server\'s grid-cell ceiling.',
  }, 'no estimate');
  const gbArith = h('div.gb-arith', '—');
  const gbCount = h('div.gb-count', '—');
  const gbCeil = h('div.gb-line');
  const gbRefusal = h('div.gb-refusal.hidden');
  const gbFacts = h('div.gb-line');

  /* TWO STANDING FACTS, EACH BEHIND ITS OWN CLAIM.
   *
   * The block above this — the arithmetic, the count, the ceiling, and the
   * refusal the server would write — is the part that stops a run being wasted,
   * and it is never folded, never conditional and never smaller than anything
   * around it.
   *
   * These two are different in kind. They do not move with dx: one is a
   * measurement of what a 512 MB host does, the other is the single runtime
   * this repository records. Both were rendered in full on every repaint, which
   * put ~190 words of unchanging prose between GRID DX and the bottom of the
   * form and cost the refusal preview above them the attention they were
   * spending.
   *
   * The FOLD HIDES THE ELABORATION, NOT THE CLAIM. Each summary carries the
   * sentence that matters — "a grid under the ceiling can still kill the
   * server", and whether the recorded runtime is this spec's — at the same type
   * size as the body, with the numbers that qualify inside. A reader who never
   * opens either one has still read both claims. */
  const gbMemBody = h('div');
  const gbMem = h('details.gb-line.gb-fold',
    h('summary', {
      style: { cursor: 'pointer', lineHeight: 'var(--lh-body)' },
      title: 'The measurement is recorded in app/server/api.py.',
    },
    h('b', 'A grid under the ceiling can still kill the server.'),
    ' Measured: dx = 12 um OOM-kills a 512 MB host at under 1 % of the '
    + 'default ceiling.'),
    gbMemBody);

  const gbRunTail = h('span');
  const gbRunBody = h('div');
  const gbRun = h('details.gb-line.gb-fold',
    h('summary', {
      style: { cursor: 'pointer', lineHeight: 'var(--lh-body)' },
      title: 'docs/DEPLOYMENT.md, the container acceptance run.',
    },
    h('b', 'Runtime: '),
    'this repository records exactly one end-to-end Studio run time. ',
    gbRunTail),
    gbRunBody);

  // Sits directly under GRID DX, as the whole of question 3, because it is
  // about this section crossed with this dx — the control immediately above it
  // and the one chosen at the top of the form. (It used to sit between GRID DX
  // and DRUG, when those two were adjacent; the questions reordered around it
  // and it stayed attached to the field it is about.)
  //
  // IT IS STILL NOT A SECTION OF ITS OWN, and the whole block is still
  // uncollapsed: it has no state to configure, and collapsing the estimate or
  // the refusal preview would hide a limit (UI_ARCHITECTURE.md §2.4). The two
  // <details> inside it hold standing facts that do not move with dx, and each
  // one states its claim in its own summary — see the comment on gbMem.
  const gridBudget = h('div.span2.gridbudget',
    h('div.gb-head', h('span.gb-t', 'Grid budget'), h('span.grow'), gbChip),
    gbArith, gbCount, gbCeil, gbRefusal, gbFacts, gbMem, gbRun);

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

  /* -------------------------------------------- what a choice will mean --
   *
   * Four read-outs, painted under the control they belong to, every one of them
   * from a payload the server already sent and not one of them a prediction:
   *
   *   sectionScale  what this specimen IS — platform, count, extent, panel
   *                 size — rendered by tissue.js::sectionScaleRow, so the
   *                 picker and the ingest panel cannot word it differently.
   *                 The picker used to say nothing but the label until the run
   *                 came back.
   *   drugNote      the two numbers that decide whether a compound can reach
   *                 anything: its IC50 with provenance, and its BBB
   *                 penetrance. Both were crammed into the <option> text,
   *                 where they are unreadable and gone the moment the list
   *                 closes.
   *   doseNote      plasma_peak over that IC50 — js/ui/doseinput.js's own
   *                 wording, moved to the field that fixes the ratio. It
   *                 predicts no verdict; see plasmaPeakNote().
   *   windowNote    total_time_s in hours, because 21600 is not a duration
   *                 anyone reads.
   */
  const noteStyle = {
    // No font-size. These qualify the value in the control directly above
    // them, `.inp` is --fs-sm, and a qualifier set smaller than the number it
    // qualifies is the one thing this rail may not do. Inheriting body prose
    // size makes them exactly equal.
    lineHeight: 'var(--lh-body)', color: 'var(--fg-2)',
    margin: 'var(--sp-1) 0 var(--sp-2)',
  };
  const sectionScale = h('div.sec-scale');
  const drugNote = h('div.span2.q-note', { style: noteStyle });
  const doseNote = h('div.span2.q-note', { style: noteStyle });
  const windowNote = h('div.span2.q-note', { style: noteStyle });

  /* ------------------------------------------------------ the fold ------
   * The nine defaults. See the QUESTIONS comment at the top of this file for
   * why they fold and what the fold is forbidden from hiding. */
  const advValues = h('span.mono', {
    style: {
      color: 'var(--fg-1)', fontWeight: '400', textTransform: 'none',
      letterSpacing: '0',
    },
  }, '—');
  const advCount = h('span.chip.st-queued', {
    title: 'How many of the nine fields in this block differ from the spec '
      + 'this session booted with (GET /api/defaults).',
  }, 'all default');

  const advBody = h('div.fgrid', { style: { padding: 'var(--sp-3) 0 var(--sp-2)' } },
    h('div.span2.hr-lbl', 'Model routes'),
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
    h('div.span2.hr-lbl', 'Dose schedule'),
    row('Half life', 'f-pk-half', pkHalf, 's', 'pk.half_life_s — plasma elimination half life'),
    row('Doses', 'f-pk-n', pkN, 'n', 'pk.n_doses — how many doses in the course'),
    row('Interval', 'f-pk-int', pkInt, 's', 'pk.dose_interval_s — time between doses'),
  );

  const advFold = h('details.span2.advfold', { id: 'sect-spec-defaults' },
    h('summary', {
      style: {
        cursor: 'pointer', padding: 'var(--sp-3) 0',
        lineHeight: 'var(--lh-body)',
      },
      title: 'Nine fields that arrive from GET /api/defaults and change the '
        + 'model route, the MGMT scalar and the shape of the plasma curve. '
        + 'Their values are printed here whether this is open or closed.',
    },
    h('span.lbl', 'Defaults in use'), ' ', advCount, ' ', advValues),
    advBody);

  const qHead = (n) => {
    const [num_, text] = QUESTIONS[n];
    return h('div.span2.hr-lbl.fq',
      // A heading, not a styled div: this rail had no heading structure at all
      // below the group name, so a screen-reader user had no route from the top
      // of the form to "at what resolution". `font: inherit` keeps the UA
      // heading styles out of a row whose type is already set by .hr-lbl —
      // the same reset rail.js::sectTitle documents for .sect-h.
      h('h4', {
        style: {
          font: 'inherit', margin: '0', letterSpacing: 'inherit',
          textTransform: 'inherit', color: 'inherit',
        },
      }, `${num_} · ${text}`));
  };

  const formEl = h('div.fgrid',
    qHead(0),
    row('Section', 'f-section', h('div', sectionSel, sectionScale, sectionCaveats), null,
      'GET /api/sections'),

    qHead(1),
    row('Drug', 'f-drug', drugSel, null, 'GET /api/drugs; null = reconstruction only'),
    drugNote,
    row('Plasma peak', 'f-pk-peak', pkPeak, 'uM',
      'pk.plasma_peak — the peak concentration the plasma curve reaches. It is '
      + 'a concentration, not an administered dose.'),
    doseNote,
    row('Total time', 'f-total', totalInp, 's', 'simulated wall time of the treatment course'),
    windowNote,

    qHead(2),
    row('Grid dx', 'f-dx', dxInp, 'um', 'reconstruction grid pitch'),
    gridBudget,

    advFold,
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

  // The shortcut is main.js's, bound on window and live even while a field has
  // focus (`wireKeys`: the Ctrl/Cmd+Enter branch is tested BEFORE the "is the
  // operator typing" early return). It is printed here because a keyboard path
  // nobody is told about is not a keyboard path.
  const SUBMIT_KEYS = /Mac|iPhone|iPad/.test(
    (typeof navigator !== 'undefined' && navigator.userAgent) || '') ? '⌘⏎' : 'Ctrl+Enter';

  const submitBtn = h('button.btn.primary', {
    id: 'btn-submit-run',
    // THE PREVIEW IS DROPPED THE INSTANT THE REAL SOLVE STARTS, before the
    // POST, not after it answers. A preview number kept beside a real one is
    // one readout with two live values, and the gap between them is not a
    // finding — it is the measured grid difference, which is already stated
    // next to the preview in its own words.
    onclick: () => { previewPanel.supersede(); onSubmit(store.get().spec); },
    title: `POST /api/runs — ${SUBMIT_KEYS} submits from anywhere in the Studio`,
  }, 'Submit run');
  const validateBtn = h('button.btn', { onclick: () => runValidate(true), title: 'POST /api/validate' }, 'Validate');
  // Validation state is a STATUS, so it wears the st-* vocabulary. It used to
  // borrow measured/assumed/unknown, which put a green "measured" chip on the
  // fact that a JSON object passed a schema check.
  const okBadge = h('span.chip.st-queued', 'unvalidated');

  /* WHY SUBMIT IS OFF, IN WORDS, AND AS A WAY TO THE THING THAT IS WRONG.
   *
   * The gate is not touched. `submitBtn.disabled` is still set by the server's
   * own verdict on POST /api/validate, and the one branch that deliberately
   * leaves it ENABLED — our own validator outage, where blocking the operator
   * on our failure would be the wrong refusal — still does.
   *
   * What was missing is the sentence. A greyed primary button is the only
   * control in this rail that can be off, and it said why in nothing at all:
   * not a title, not a label, not a line. A first-time operator whose
   * catalogue has not loaded sees a dead Submit and no way to learn that the
   * reason is that no section is chosen.
   *
   * It is a BUTTON, not a caption, for two reasons: it is focusable, so the
   * reason is reachable by keyboard while the disabled button it explains is
   * not, and pressing it moves focus to the field that is blocking — which is
   * the actual next action in every case. */
  const submitWhy = h('button.btn.ghost.sm.submit-why', {
    type: 'button', hidden: true,
    onclick: () => focusBlocker(),
    // `.btn` is one --h-ctl row of `white-space: nowrap`, which is right for a
    // verb and wrong for a sentence: unstyled, "Waiting for the section
    // catalogue (GET /api/sections)." is a 320 px word that pushes the rail
    // into a horizontal scroll — the one thing docs/UI_ARCHITECTURE.md §4.3
    // forbids outright. It takes its own full-width line under the buttons and
    // wraps there. `.submit-why` is named for the stylesheet owner; until it
    // exists these five properties are the whole of it.
    style: {
      flex: '1 1 100%', maxWidth: '100%',
      height: 'auto', minHeight: 'var(--h-ctl-sm)',
      whiteSpace: 'normal', textAlign: 'left', justifyContent: 'flex-start',
      lineHeight: 'var(--lh-body)', padding: 'var(--sp-2) var(--sp-3)',
    },
  }, '');

  const submitKeys = h('span.hint.mono', {
    title: `${SUBMIT_KEYS} submits the current RunSpec from anywhere in the `
      + 'Studio, including from inside a field.',
  }, SUBMIT_KEYS);

  const submitRow = h('div.frow', {
    style: { padding: 'var(--sp-4) var(--sp-5)', flexWrap: 'wrap' },
  },
  submitBtn, validateBtn, submitKeys, h('span.grow'), okBadge, submitWhy);

  /* ------------------------------------------------- the submit state --
   *
   * Four states, one of which is "nothing is wrong". Each one owns a sentence
   * and the control that answers it; `setBlock()` is the only writer of
   * `submitBtn.disabled` in this file, so the button and its explanation
   * cannot disagree.
   *
   * The SECTION gate is evaluated here rather than pushed in by a caller,
   * because it is a property of the form and not of any response: a RunSpec
   * with no `section` names no specimen, and at boot — before the first
   * validation lands — the button was enabled over exactly that spec.
   */
  let submitting_ = false;
  let blockWhy_ = '';
  let blockTarget_ = null;

  function sectionMissing() {
    const spec = store.get().spec || {};
    if (spec.section) return '';
    // Three states, and they are not the same sentence. "The catalogue has not
    // arrived" is a wait; "the catalogue arrived with nothing runnable in it"
    // is a dead end with a next action; "pick one" is neither. Saying the
    // second while the first is true would be this panel making a claim about
    // a request that has not come back.
    if (!sectionsSeen_) return 'Waiting for the section catalogue (GET /api/sections).';
    const selectable = Array.from(sectionSel.options).some((o) => !o.disabled);
    return selectable
      ? 'Choose a section first — this RunSpec names no specimen.'
      : 'GET /api/sections returned no section this Studio can run. Ingest one '
        + 'under Your tissue.';
  }

  /** Set (or clear) the reason Submit is off, and the control that fixes it. */
  function setBlock(why, target) {
    blockWhy_ = why || '';
    blockTarget_ = why ? (target || null) : null;
    applySubmit();
  }

  function applySubmit() {
    const missing = sectionMissing();
    const why = missing || blockWhy_;
    const target = missing ? sectionSel : blockTarget_;
    submitBtn.disabled = submitting_ || !!why;
    submitBtn.textContent = submitting_ ? 'Submitting…' : 'Submit run';
    // While the run is being posted the button says so itself, so a second
    // sentence beside it would be noise about a state that lasts one request.
    submitWhy.hidden = submitting_ || !why;
    submitWhy.textContent = why;
    submitWhy.title = target
      ? 'Press to go to the control that is blocking the run.'
      : '';
  }

  /** The control the current block is about. Derived, never stored: the
   *  section gate can take over from a validator error and back again, and a
   *  target cached at the moment one of them fired outlives it. */
  function currentTarget() {
    return sectionMissing() ? sectionSel : blockTarget_;
  }

  /**
   * The control the FIRST validator error is about, or the editor.
   *
   * A path the form has no control for is not a dead end: every field of the
   * RunSpec is editable in the source buffer, so an error on one the form does
   * not render (an unrecognised key, a nested PK field the form flattens) still
   * has somewhere to send the operator.
   */
  function blockerFor(errors) {
    for (const e of errors || []) {
      const id = FIELD_OF_PATH[e.path];
      const t = id && fieldEls.get(id);
      if (t && t.control) return t.control;
    }
    return ta;
  }

  /** Put the cursor on whatever is blocking the run. */
  function focusBlocker() {
    const t = currentTarget();
    if (!t) return;
    // A control inside a closed <details> cannot take focus, so anything that
    // folds is opened before it is asked to.
    if (t === ta) srcSect.open = true;
    if (FOLDED_FIELDS.has(t.id)) advFold.open = true;
    try { t.scrollIntoView({ block: 'nearest' }); } catch { /* older engines */ }
    t.focus();
  }

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

  // Where the job actually lives. A reconstruction is minutes, not seconds,
  // which makes "can I go and do something else" the first question an
  // operator has and the one this panel never answered: it showed a progress
  // bar and said nothing about whose process was making it move. It is the
  // SERVER's — POST /api/runs returns an id and the work continues with no
  // browser attached — and main.js::reattachLiveRun now picks a run back up on
  // load, so this states a behaviour the product has rather than a
  // reassurance. The cancel clause is the same cooperative-cancellation fact
  // cancelTracked() already explains when you press the button; it is worth
  // knowing BEFORE you press it.
  const jobNote = h('div.dim', {
    style: { fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)' },
  },
  'The run executes in the Vivum server process, not in this tab. Closing the '
  + 'tab or reloading does not stop it, and this panel picks up a run that is '
  + 'still in flight when the page loads. Cancel is what stops one, and it '
  + 'takes effect at the next stage boundary rather than instantly.');

  const jobEl = h('div.stack.tight',
    h('div.frow', jobId, jobStatus, h('span.grow'), cancelBtn),
    jobProgress,
    h('div.frow',
      h('span.lbl', 'stage'), jobStage, h('span.grow'),
      jobPct, h('span.unit', '%'),
      h('span.dim', '·'), jobElapsed),
    logEl,
    jobNote,
  );

  /* -------------------------------------------------------------- root -- */

  // A scrolling body with a fixed commit row pinned to the bottom.
  //
  // Two groups, matching the READ rail opposite: DEFINE is everything that
  // mutates the RunSpec, COMMIT is what happens to it afterwards. The commit
  // row itself is outside the scroll on purpose — it is the one control whose
  // reachability must not depend on scroll position.
  // See docs/UI_ARCHITECTURE.md section 2.2.

  // THE SOURCE EDITOR IS CLOSED ON ARRIVAL, since 2026-09-20.
  //
  // It is the same RunSpec the form above it holds, in JSON, kept in two-way
  // sync — so on the first screen a stranger sees, the Setup rail asked them to
  // read the same twelve fields twice, the second time as source code with a
  // line gutter. It is roughly a third of the rail's height and it is not a
  // decision; it is a second view of decisions taken above.
  //
  // Nothing about it is hidden or moved. It keeps its section, its place in the
  // Define group and its summary; it opens on a click, on Ctrl/Cmd+L
  // (`focusEditor`, wired in main.js::wireKeys), and BY ITSELF whenever the
  // buffer stops parsing or the validator points at a line in it — the two
  // states where the form cannot show you what is wrong.
  const srcSect = h('details.sect', { id: 'sect-source' },
    h('summary', sectTitle('Source'), h('span.hint', 'JSON · two-way sync')),
    h('div.sect-body', { style: { padding: '0' } }, editorEl));

  const body = h('div.rail-body',
    railGroup('define', 'Define', 'RunSpec v1',
      h('details.sect', { open: true },
        h('summary', sectTitle('Specimen, therapy, model'), h('span.hint', '3 questions')),
        h('div.sect-body', formEl)),
      // Both BYO panels sit directly under the form that picks a section and a
      // drug, because that is the moment an operator discovers the list does
      // not contain their specimen or their molecule. Tissue first: it is the
      // outer of the two choices — a compound is compared ON a section. Both
      // stay collapsed by default; each is a deliberate act, not part of the
      // ordinary path.
      tissuePanel.el,
      compoundPanel.el,
      srcSect),

    railGroup('job', 'Commit', 'POST /api/runs',
      // The preview sits in COMMIT, above the active job, because the question
      // it answers — "roughly what will this say" — is a question about a
      // decision that has not been taken yet. It is deliberately NOT in the
      // Read rail: that rail reports RESULTS, with provenance chips, exports
      // and copy-a-quotable-sentence buttons, and an unquotable number does
      // not belong on the surface whose job is to produce quotable ones.
      // js/ui/preview.js carries the full argument.
      previewPanel.el,
      h('details.sect', { open: true, id: 'sect-job' },
        h('summary', sectTitle('Active job'), h('span.hint', '/status poll')),
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
    sectionsSeen_ = true;
    // A section row is the most natural place for the server to publish the
    // ceiling — the refusal is per-section arithmetic — so the rows are swept
    // for it before the budget is painted. No row carries it today.
    for (const s of sections || []) if (adoptCeiling(s, 'GET /api/sections')) break;
    paintSectionCaveats();
    paintGridBudget();
    paintNotes();
    // The catalogue landing can only ever RELAX the section block, never add
    // one; recompute so the commit bar stops saying "waiting" the moment it is
    // no longer waiting.
    applySubmit();
  }

  /* ------------------------------------------- the consequence read-outs --
   *
   * Everything below paints from a payload the server already sent. Nothing is
   * computed, nothing is predicted, and every branch that has no number to
   * print says so rather than printing a shape with nothing in it.
   */

  // Last list from GET /api/drugs. Kept for the same reason `sections_` is: the
  // facts that decide whether a compound can reach anything — its IC50 and its
  // penetrance — arrive with the catalogue and were being spent on one line of
  // <option> text that is unreadable and gone the moment the list closes.
  let drugs_ = [];

  /** The spec this session booted with, for the fold's "N changed" count. */
  let stock_ = null;
  /** How many folded fields were off-default at the last paint. */
  let changed_ = 0;
  /** Has GET /api/sections come back at all? See sectionMissing(). */
  let sectionsSeen_ = false;

  function paintNotes() {
    /* --- 1. the specimen: what it IS. -------------------------------- */
    clear(sectionScale);
    const sec = (sections_ || []).find((x) => x.id === sectionSel.value) || null;
    // appendChild(null) is a TypeError that takes the whole render with it, so
    // the nullable branch is guarded at the call site, here, not inside.
    if (sec) sectionScale.appendChild(sectionScaleRow(sec));

    /* --- 2. the compound, and the dose against its potency. ---------- */
    clear(drugNote);
    clear(doseNote);
    const chosen = drugSel.value;
    const drug = chosen ? (drugs_ || []).find((d) => d.name === chosen) || null : null;
    if (!chosen) {
      drugNote.appendChild(h('span',
        h('b', 'No compound — reconstruction only. '),
        'The section is reconstructed and nothing is transported: no delivery '
        + 'map, no mechanism attribution, no kill. The dose and the window '
        + 'below still size the run, and every panel on the Read rail that '
        + 'describes a therapy stays empty.'));
    } else if (!drug) {
      drugNote.appendChild(h('span',
        'This spec names ', h('span.mono.strong', chosen),
        ', which is not in the catalogue this Studio loaded (GET /api/drugs). '
        + 'Its IC50 and penetrance are unknown here, so the ratio below is not '
        + 'shown. The server decides whether the name resolves.'));
    } else {
      const bbb = Number(drug.bbb_penetrance);
      const bits = [];
      bits.push(Number.isFinite(Number(drug.ic50_uM))
        ? h('span', h('span.lbl', 'IC50 '),
          h('span.mono.strong', num(drug.ic50_uM, 1)), h('span.unit', 'uM'),
          drug.ic50_source ? h('span', ` — ${drug.ic50_source}`) : null)
        : h('span', 'This server states no IC50 for this compound.'));
      bits.push(h('div',
        h('span.lbl', 'BBB penetrance '),
        Number.isFinite(bbb)
          ? h('span', h('span.mono.strong', num(bbb, 2)),
            h('span.unit', 'of plasma'))
          : h('span.dim', 'not stated'),
        // The one value of this field that decides the run before it starts.
        // js/ui/reach.js states the same arithmetic on the result: a zero
        // penetrance gives an identically zero source term, so the exposure
        // field is zero everywhere and the reach map is refused. Saying it
        // here costs nothing; discovering it after a reconstruction costs the
        // run.
        Number.isFinite(bbb) && bbb === 0
          ? h('span', ' — ', h('b', 'zero'),
            '. The source term is then identically zero at every vessel, so '
            + 'the exposure field is zero everywhere and the reach map is '
            + 'refused. That is the compound\'s parameterisation, not a '
            + 'finding this run will make.')
          : null));
      if (drug.is_user_compound) {
        bits.push(h('div',
          h('b', 'This is a compound you compiled. '),
          `Every transport parameter except molecular weight and TPSA is a `
          + `prediction (${drug.predicted_fields ?? '?'}/`
          + `${drug.total_fields ?? '?'} fields), and setting it beside another `
          + 'compound is not defensible — compound ordering was measured and '
          + 'withdrawn.'));
      }
      append(drugNote, bits);

      // The ratio that is fixed the moment the dose is typed. One wording,
      // authored in js/ui/doseinput.js, which is also where the result rail
      // reads it from.
      const note = plasmaPeakNote({
        peak_uM: parseFloat(pkPeak.value),
        ic50_uM: drug.ic50_uM,
        ic50Source: drug.ic50_source,
        userCompound: !!drug.is_user_compound,
      });
      if (note) doseNote.appendChild(note);
    }

    /* --- 2b. the window, in hours rather than in seconds. ------------ */
    clear(windowNote);
    const tot = parseFloat(totalInp.value);
    windowNote.appendChild(Number.isFinite(tot) && tot > 0
      ? h('span',
        h('span.mono.strong', `${int(tot)} s = ${dur(tot)}`),
        ' of simulated treatment time. The solve runs over this window, and '
        + 'every exposure figure the run reports is an average or a total '
        + 'across it, so changing it moves them.')
      : h('span', 'Enter a positive number of seconds to set the window.'));

    /* --- 3. the nine folded defaults, printed whether open or shut. -- */
    paintFoldSummary();

    // An empty note is still a grid row and still carries its margins, so a
    // spec with no drug left two ghost gaps between DRUG and PLASMA PEAK.
    // `display:none` takes it out of the grid entirely.
    for (const el of [sectionScale, drugNote, doseNote, windowNote]) {
      el.style.display = el.firstChild ? '' : 'none';
    }
  }

  /** Which option of a segmented control is pressed, or '' if none is. */
  function segValue(key) {
    const btn = Array.from(segs[key].children)
      .find((b) => b.getAttribute('aria-pressed') === 'true');
    return btn ? btn.dataset.val : '';
  }

  /**
   * The folded block's summary line: every value it holds, and how many of them
   * have been moved off the spec this session booted with.
   *
   * This is the whole licence for folding them. A `<details>` that hid nine
   * values would be a menu; one that prints all nine in its summary hides the
   * widgets and nothing else.
   */
  function paintFoldSummary() {
    const mg = mgmtChk.checked && mgmtInp.value !== ''
      ? `MGMT ${num(parseFloat(mgmtInp.value), 2)} log2 CPM`
      : 'MGMT unset';
    const halfS = parseFloat(pkHalf.value);
    const intS = parseFloat(pkInt.value);
    const nDose = parseFloat(pkN.value);
    const parts = [
      segValue('sensitivity_mode') || 'sensitivity —',
      segValue('vascular_mode') || 'vascular —',
      segValue('oxygen_mode') || 'oxygen —',
      segValue('state_mode') || 'state —',
      mg,
      Number.isFinite(halfS) ? `half life ${dur(halfS)}` : 'half life —',
      Number.isFinite(nDose) ? `${int(nDose)} dose${nDose === 1 ? '' : 's'}` : 'doses —',
      Number.isFinite(intS) ? `every ${dur(intS)}` : 'interval —',
    ];
    advValues.textContent = parts.join(' · ');

    // "Changed" is measured against the spec this session booted with, and only
    // when there is one. With no snapshot the chip says so rather than
    // claiming everything is default — the two are not the same statement.
    if (!stock_) {
      advCount.className = 'chip st-queued';
      advCount.textContent = '9 fields';
      return;
    }
    const cur = readForm();
    // null and undefined are the same statement here (MGMT unset), and a
    // number that arrived as a string from an <input> is the same value as the
    // one in the snapshot.
    const same = (a, b) => (a === null || a === undefined
      ? (b === null || b === undefined)
      : (a === b || Number(a) === Number(b)));
    const sp = stock_.pk || {};
    const cp = cur.pk || {};
    let n = 0;
    for (const k of ['sensitivity_mode', 'vascular_mode', 'oxygen_mode', 'state_mode']) {
      if (!same(cur[k], stock_[k])) n++;
    }
    if (!same(cur.mgmt, stock_.mgmt)) n++;
    for (const k of ['half_life_s', 'n_doses', 'dose_interval_s']) {
      if (!same(cp[k], sp[k])) n++;
    }
    advCount.className = 'chip ' + (n ? 'assumed' : 'st-queued');
    advCount.textContent = n
      ? `${n} of 9 changed`
      : '9 fields, all default';
    // A field that has just been moved off its default is not a detail any
    // more, so the block opens ONCE on that transition — on the transition and
    // not on the state, or a fold the operator deliberately closed would spring
    // back open on their next keystroke anywhere in the form.
    if (n && !changed_) advFold.open = true;
    changed_ = n;
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
        h('b', 'This specimen was ingested here, not shipped with Vivum. '),
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

  /* ------------------------------------------------------ grid budget -- */

  /** Adopt a ceiling from a server payload. Higher-quality sources win. */
  function adoptCeiling(payload, where) {
    const v = readCeiling(payload);
    if (v === null) return false;
    ceiling_ = v;
    ceilingSrc_ = where;
    // "Tuned" means the operator lowered it with TMESIM_MAX_GRID_CELLS, which
    // is the clause the server adds to its own refusal. The server reports the
    // PROVENANCE (`max_grid_cells_source`), not a boolean, so this reads that
    // rather than a `..._tuned` flag that was guessed at before the endpoint
    // existed and that nothing has ever sent. If the field is absent we do not
    // infer it: the preview then simply omits the clause.
    const src = payload && (payload.max_grid_cells_source
      || (payload.limits && payload.limits.max_grid_cells_source));
    ceilingTuned_ = src === 'TMESIM_MAX_GRID_CELLS';
    return true;
  }

  /** Adopt a ceiling out of a refusal the server has already written. */
  function adoptCeilingFromDiags(diags) {
    for (const d of diags || []) {
      const msg = (d && d.message) || '';
      const m = CEILING_FROM_REFUSAL.exec(msg);
      if (!m) continue;
      const v = Number(m[1].replace(/,/g, ''));
      if (!Number.isFinite(v) || v <= 0) continue;
      ceiling_ = v;
      ceilingSrc_ = 'the server\'s own refusal for this spec';
      ceilingTuned_ = msg.includes(CEILING_TUNED_MARK);
      return true;
    }
    return false;
  }

  function paintGridBudget() {
    const sec = (sections_ || []).find((x) => x.id === sectionSel.value) || null;
    const dx = parseFloat(dxInp.value);
    const dxOk = Number.isFinite(dx) && dx > 0;
    const extent = sec && Array.isArray(sec.extent_um) && sec.extent_um.length >= 2
      && Number.isFinite(Number(sec.extent_um[0]))
      && Number.isFinite(Number(sec.extent_um[1]))
      ? [Number(sec.extent_um[0]), Number(sec.extent_um[1])]
      : null;
    const n = gridCells(extent, dx);

    // gbMem and gbRun are <details>: clearing the element itself would take
    // the summary — the claim — with it, and reset an open fold to shut on
    // every keystroke in dx. Only their bodies are repainted.
    for (const el of [gbArith, gbCount, gbCeil, gbRefusal, gbFacts,
      gbMemBody, gbRunBody, gbRunTail]) clear(el);
    gbRefusal.classList.add('hidden');

    const chip = (klass, text, title) => {
      gbChip.className = 'chip ' + klass;
      gbChip.textContent = text;
      gbChip.title = title;
    };

    /* --- what dx MEANS, stated whatever else is or is not known. --------- */
    gbFacts.appendChild(h('span',
      dxOk
        ? h('span.mono.strong', `dx = ${pyG(dx)} um`)
        : h('span.mono.strong', 'dx'),
      ' is the edge length of one square grid cell — the mesh the '
      + 'reconstruction lays over the section and every transport solve steps '
      + 'on. Halving it quadruples the grid.'));

    /* --- which section, and its REAL cell count. ------------------------ */
    if (!sec) {
      gbFacts.appendChild(h('div',
        'No section is selected in this form, or the section list has not '
        + 'loaded yet (GET /api/sections).'));
    } else {
      const ident = h('div',
        'Section ', h('span.mono.strong', sec.id), ' — ');
      if (sec.n_cells === null || sec.n_cells === undefined) {
        ident.appendChild(h('span', 'this section reports no cell count'
          + (sec.n_cells_reason ? ': ' + sec.n_cells_reason : '.')));
      } else {
        ident.appendChild(h('span.mono.strong', int(sec.n_cells)));
        ident.appendChild(h('span', ' measured cells'));
        if (n !== null && n > 0) {
          // Arithmetic on two numbers the server supplied, labelled as such.
          // It is the honest answer to "is this mesh fine or coarse for this
          // specimen" and it is not a claim about the physics.
          ident.appendChild(h('span',
            `, i.e. ${num(sec.n_cells / n, 2)} measured cells per grid cell.`));
        } else {
          ident.appendChild(h('span', '.'));
        }
      }
      gbFacts.appendChild(ident);
    }

    /* --- the count, the ceiling, the refusal. --------------------------- */
    if (!sec) {
      gbArith.textContent = '—';
      gbCount.textContent = '— grid cells';
      chip('st-queued', 'no section', 'No section row to take an extent from.');
      gbCeil.textContent = 'The grid size cannot be computed until a section is '
        + 'selected: it is that section\'s measured extent divided by dx.';
    } else if (!extent) {
      gbArith.textContent = '—';
      gbCount.textContent = '— grid cells';
      chip('st-queued', 'extent unknown',
        'GET /api/sections returned no extent_um for this section.');
      gbCeil.appendChild(h('span',
        'This section reports no ', h('span.mono', 'extent_um'),
        ', so the grid size cannot be computed here'
        + (sec.extent_um_reason ? ': ' + sec.extent_um_reason : '.')
        + ' The server evaluates its ceiling from the same field, so it will '
        + 'not refuse on grid size either — the run proceeds and whatever the '
        + 'reconstruction allocates is unbounded by that check.'));
    } else if (!dxOk) {
      gbArith.textContent = '—';
      gbCount.textContent = '— grid cells';
      chip('st-queued', 'dx not set', 'dx_um must be a positive number.');
      gbCeil.textContent = 'Enter a positive dx to size the grid.';
    } else {
      const [ex, ey] = extent;
      gbArith.textContent = `(${num(ex, 0)} / ${pyG(dx)}) x (${num(ey, 0)} / ${pyG(dx)}) um`;
      clear(gbCount);
      append(gbCount, [h('span.mono.strong', int(n)), h('span.gb-u', ' grid cells')]);

      if (ceiling_ === null) {
        chip('st-queued', 'ceiling unknown',
          'No server payload has reported ' + CEILING_KEY + '.');
        gbCeil.appendChild(h('span',
          'This server does not report its grid-cell ceiling, so how close '
          + 'this spec is to a refusal is ', h('b', 'unknown here'),
          '. The ceiling is real — ', h('span.mono', '_MAX_GRID_CELLS'),
          ' in ', h('span.mono', 'app/server/api.py'),
          ', which a constrained host is expected to lower through ',
          h('span.mono', 'TMESIM_MAX_GRID_CELLS'),
          ' — and no number is printed here because a hardcoded one would be '
          + 'wrong on exactly those hosts. Press Validate: the server answers '
          + 'with the ceiling in its refusal, and this line then shows it.'));
      } else {
        const frac = n / ceiling_;
        // The dx at which this section reaches the ceiling, from the inverse
        // of the server's own expression: n = ex*ey/dx^2 > C  <=>  dx <
        // sqrt(ex*ey/C). It answers "how close am I" in the units of the
        // control the operator is actually turning.
        const dxAtCeiling = Math.sqrt((ex * ey) / ceiling_);
        const provenance = ceilingTuned_
          ? 'configured on this server through TMESIM_MAX_GRID_CELLS'
          : `reported by ${ceilingSrc_}`;
        gbCeil.appendChild(h('span',
          'Ceiling ', h('span.mono.strong', int(ceiling_)),
          ` grid cells (${provenance}), reached on this section at `,
          h('span.mono.strong', `dx = ${num(dxAtCeiling, 2)} um`),
          '. This spec is at ',
          h('span.mono.strong', `${num(frac * 100, frac < 0.1 ? 3 : 1)} %`),
          ' of it.'));

        if (n > ceiling_) {
          chip('st-failed', 'will be refused',
            'The server refuses this spec on grid size.');
          gbRefusal.classList.remove('hidden');
          append(gbRefusal, [
            h('div.gb-r-lbl', 'the server will refuse this spec'),
            h('div.gb-r-msg', refusalText(dx, n, ceiling_, ceilingTuned_)),
            h('div.gb-r-note',
              'Reproduced from the refusal in app/server/api.py so that what '
              + 'you read before submitting is what you would read after. The '
              + 'server writes the real one; if the two ever differ, the '
              + 'server\'s is the true one.'),
          ]);
        } else {
          // NOT green. Under the ceiling is not the same as affordable — see
          // the memory note directly below, where a spec at 0.7 % of the
          // default ceiling kills a 512 MB container. A green "ok" chip here
          // would be the loudest reassurance the palette has, attached to a
          // fact that does not support it.
          chip('st-queued', 'under ceiling',
            'Below the reported ceiling. That is not the same as affordable '
            + 'on this host — read the memory note.');
        }
      }
    }

    /* --- grid size is not the same as memory. Measured, and recorded. ---- */
    gbMemBody.appendChild(h('span',
      'app/server/api.py records the measurement: on a 512 MB instance '
      + '(Render\'s free tier), flagship Xenium GBM, in a real container, ',
      h('span.mono.strong', 'dx = 30 um'), ' peaks at 220 MiB, and ',
      h('span.mono.strong', 'dx = 12 um'),
      ' OOM-KILLS the container — the process dies, and /healthz stops '
      + 'answering, so one ordinary request takes the service down for '
      + 'everyone. That spec is under 1 % of the 40,000,000 default ceiling. A '
      + 'constrained host has to lower TMESIM_MAX_GRID_CELLS to what it can '
      + 'actually afford.'));

    /* --- runtime, only where the repo records one. ----------------------- */
    const st = store.get();
    const matches = matchesRecordedRun(st && st.spec);
    const r = RECORDED_RUN;
    // The half of this claim that MOVES WITH THE FORM goes in the summary, so
    // "that is not this spec" is read whether or not the fold is opened. It is
    // the half that decides whether the number means anything here.
    append(gbRunTail, [
      h('span.mono.strong', `${num(r.seconds, 2)} s`),
      matches
        ? h('span', ' — ', h('b', 'that is this spec.'))
        : h('span', ' — ', h('b', 'that is NOT this spec, and it does not '
          + 'transfer.')),
    ]);
    gbRunBody.appendChild(h('span',
      'Measured ',
      h('span.mono.strong', `${num(r.seconds, 2)} s`),
      ` on ${r.host} (${r.where}), for the stock defaults spec: `,
      h('span.mono', `${r.spec.section} · ${r.spec.drug} · dx ${r.spec.dx_um} um`
        + ` · ${r.spec.total_time_s / 3600} h · plasma peak ${r.spec.plasma_peak} uM`),
      '. ',
      matches
        ? h('span', h('b', 'That is this spec.'),
          ' It is one machine and one cold start, not a benchmark.')
        : h('span', h('b', 'That is not this spec, and the number does not '
          + 'transfer. '),
        'The window, the mesh and the compound each move it, and this repo '
        + 'records no measured relation between dx and run time for a Studio '
        + 'run, so none is estimated here.')));
  }

  function setDrugs(drugs) {
    // Kept, not just rendered: the IC50 and the penetrance in these rows are
    // what the note under the picker and the ratio under the dose field are
    // made of, and a <select> remembers none of it.
    drugs_ = Array.isArray(drugs) ? drugs : [];
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
             + 'kill, AUC and penetration depth are NOT defensible, and neither '
             + 'is setting it beside another compound — compound ordering was '
             + 'measured and withdrawn. What it supports is a map: run it on '
             + 'this section and read where it reaches relative to the '
             + 'best-perfused tissue in that same run. Note that the shape of '
             + 'that map rests on a diffusivity from a two-route model scoring '
             + 'R^2 = 0.718 blind against 71 measured compounds, which still '
             + 'recovers only about a quarter of the measured spread from one '
             + 'cell line, so a similar-mass analogue moves the gradient less '
             + 'than it should.')
          : '',
      }, label + suffix));
    }
    // The picker is repainted by main.js on every compound change, and the
    // selection survives it — so the notes have to be repainted with it or
    // they describe the previous catalogue.
    paintNotes();
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
    paintGridBudget();
    paintNotes();
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
    // Live, on every keystroke in `dx` and every change of section — the whole
    // point is that the operator sees the cost while turning the dial, not
    // 340 ms later when the debounced validator lands, and not after Submit.
    paintGridBudget();
    // Same reason the budget is painted live: the ratio under the dose field
    // and the window in hours are both answers to what the operator is typing
    // RIGHT NOW, not to what the validator will say about it in 340 ms.
    paintNotes();
    // The section block is a property of the form, so it answers on the
    // keystroke rather than 340 ms later when the debounced validator lands.
    applySubmit();
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
      // The form is deliberately LEFT ALONE on a parse failure (see the sync
      // contract at the top of this file), which means the form cannot show
      // what is wrong — so the editor that can is opened rather than left shut
      // with a dead button above it.
      srcSect.open = true;
      setBlock('The RunSpec source does not parse — fix the line marked in '
        + 'Source.', ta);
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
    // The Source section is closed by default and focus cannot enter a closed
    // <details>, so every route INTO the editor opens it first. Without this a
    // click on a diagnostic silently did nothing.
    srcSect.open = true;
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
        // `scroll-behavior: auto !important` in base.css's reduced-motion
        // block does NOT govern a scriptive smooth scroll — a `behavior`
        // passed to scrollBy wins over the CSS property, so this was the one
        // animation in the Studio that ignored the setting entirely. It is
        // also the only JS-driven motion in app/web, which is why the check is
        // made here rather than wrapped in a helper.
        onclick: () => scroll.scrollBy({
          top: scroll.clientHeight - 20,
          behavior: window.matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth',
        }),
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
    for (const e of errors) {
      const id = FIELD_OF_PATH[e.path];
      if (!id) continue;
      // An error on a folded field is an error nobody can see. The fold is a
      // convenience for fields that are behaving.
      if (FOLDED_FIELDS.has(id)) advFold.open = true;
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
      // The verdict body is the third place the ceiling could be published,
      // and the refusal itself is the fourth — a payload field wins, and the
      // refusal text is only read when nothing has published one.
      if (!adoptCeiling(res, 'POST /api/validate') && ceiling_ === null) {
        adoptCeilingFromDiags(errors);
      }
      store.set({ validation: { ok, errors, warnings }, validating: false });
      renderDiags(errors, warnings);
      paintGridBudget();
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
      setBlock(ok ? '' : `${errors.length} error`
        + `${errors.length === 1 ? '' : 's'} from the validator — press to go `
        + 'to the first one.', ok ? null : blockerFor(errors));
    } catch (err) {
      store.set({ validating: false });
      okBadge.className = 'chip st-failed';
      okBadge.textContent = 'validate failed';
      // Do not block the operator on our own outage: the server is the thing
      // that refuses a spec, and it has not been asked yet.
      setBlock('');
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
    // `loadSpec` receives the GET /api/defaults body verbatim (main.js), and
    // that body is still read for a ceiling in case a future server publishes
    // one there. It is NOT where the ceiling comes from: /api/defaults is
    // adopted verbatim as the RunSpec and posted straight back to
    // /api/validate, which warns on unrecognised fields, so a ceiling on it
    // would make every boot report an unrecognised-field warning. That is why
    // the server put it on /api/limits, fetched just below.
    adoptCeiling(spec, 'GET /api/defaults');
    // Fire-and-forget: the ceiling is an ANNOTATION on the form, so waiting for
    // it would delay the editor for a number that only qualifies what is in it.
    // A failure is swallowed HERE and nowhere else -- the panel's "ceiling not
    // known" state is the honest rendering of it, and a toast on every boot of
    // an older server would be noise about a field that server never had.
    api.getLimits().then((lim) => {
      if (adoptCeiling(lim, 'GET /api/limits')) paintGridBudget();
    }).catch(() => { /* stays unknown; the panel says so */ });
    // THE FIRST SPEC THIS SESSION SEES IS THE STOCK ONE, and the folded block
    // counts its nine fields against it. Only the first, and only when it
    // parses as an object: main.js calls loadSpec() again to restore a run's
    // own spec (and on a catalogue refresh, with the CURRENT spec), and
    // adopting one of those as "the defaults" would silently redefine what
    // "changed" means halfway through a session.
    if (!stock_ && spec && typeof spec === 'object') {
      try { stock_ = JSON.parse(JSON.stringify(spec)); } catch { stock_ = null; }
    }
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
      // A failed poll is a failure to ASK, not a failure of the run — the job
      // is in the server process and keeps going. Two things this used to get
      // wrong, both in the alarming direction:
      //
      //  * the chip read `poll error`, which names our machinery rather than
      //    what happened to the operator;
      //  * the raw message REPLACED the log tail, so a momentary blip wiped
      //    the only record of what the run had done and left `no such run
      //    'r_00079'` sitting where the progress log had been — a sentence
      //    that reads as "your run is gone" about a run that is fine.
      //
      // The log tail is now left alone and the state is annotated instead.
      // Polling continues on the same interval: the next tick repairs it.
      jobStatus.className = 'chip st-failed';
      jobStatus.textContent = 'status unavailable';
      jobStatus.title = 'The last status request did not get an answer, so the '
        + 'numbers above are from the previous one. This says nothing about the '
        + 'run itself — it is executing in the server process and is unaffected. '
        + 'The next poll retries.\n\n'
        + String((err && err.message) || err);
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

  // THE PANEL IS BUILT BEFORE ANY CATALOGUE HAS LANDED, and for the first
  // seconds — or for the whole session, if /api/sections and /api/defaults both
  // fail, which is the case main.js catches and toasts — the spec names no
  // section. Submit was ENABLED over that spec, because nothing had set its
  // disabled property yet. It is evaluated once here so the button's first
  // state is a computed one rather than an unset one.
  applySubmit();
  paintNotes();

  return {
    el: root,
    setSections, setDrugs, loadSpec, track, validate: () => runValidate(true),
    /** Repaint the operator's own compound list (GET /api/compounds). */
    refreshCompounds: () => compoundPanel.refresh(),
    /** Repaint the ingested-section list (GET /api/sections, `ingested` rows). */
    refreshTissue: () => tissuePanel.refresh(),
    get trackedId() { return tracked; },
    /** Ctrl/Cmd+L, from main.js::wireKeys. The section it lives in is closed
     *  by default now, and focus cannot enter a closed <details>. */
    focusEditor: () => { srcSect.open = true; ta.focus(); },
    setSubmitting(on) {
      submitting_ = on;
      applySubmit();
    },
  };
}
