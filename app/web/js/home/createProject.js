// createProject.js — the four dialogs that make the workspace a place you can
// keep work in: create, rename, duplicate, delete.
//
// Injects its own DOM, the way js/ui/projects.js does for its own dialog: a
// `<dialog>` built on demand, appended to `<body>`, and REMOVED on close. It
// adds nothing to home.html and owns no layout.
//
// THREE THINGS THIS FILE IS ACCOUNTABLE FOR
// -----------------------------------------
//
// 1. DELETE TELLS THE TRUTH ABOUT RUNS. docs/WORKSPACE_CONTRACT.md §2: "a
//    project owns nothing but ids", and "the consequence, which the UI must not
//    hide: deleting a project does not delete its runs. The home page says so at
//    the point of deletion, not in a help file." A destructive dialog that
//    overstates what it destroys teaches the user to distrust every other
//    number on the page; one that understates it loses their work. So the
//    delete dialog states both halves — what goes, and what stays — reads the
//    counts off the project's own document rather than guessing, and prints the
//    server's `runs_retained` list length afterwards.
//
// 2. DUPLICATE TELLS THE TRUTH ABOUT WHAT IS COPIED. Contract §1:
//    `POST /api/projects/{id}/duplicate` "copies specimens and their run ids.
//    Runs are NOT re-executed and NOT copied — a run belongs to the store, and
//    two specimens may reference the same run id." A user who believes
//    duplicate re-runs the panel will wait for results that are not coming.
//
// 3. A FIELD ERROR ATTACHES TO ITS FIELD. Never a toast — and here that is not
//    only a style rule, it is a correctness one: a modal `<dialog>` opened with
//    `showModal()` is in the TOP LAYER, which paints above every z-index on the
//    page, and css/components.css puts the toast stack at `z-index: 120`. A
//    toast raised while one of these dialogs is open is behind the backdrop and
//    is never seen. So everything that can go wrong while a dialog is open is
//    reported INSIDE it — on the field at fault where there is one, in the
//    form-level region where there is not. `toast()` and `notifyError()` are
//    used for exactly the two cases where no dialog is open: the outcome after
//    one closes, and a failure to gather what a dialog needs before it opens.
//
// WHAT IS DELIBERATELY ABSENT
// ---------------------------
// Any notion of an owner. Contract §9 and `GET /api/account` agree: one shared
// key, no identities, no attribution. There is no "your project", no creator
// field and no avatar in any of these four dialogs, because a dialog that
// implied one would be the UI asserting something the server denies.

import { h, clear } from '../util/dom.js';
import { int, ts, NULL_GLYPH } from '../util/fmt.js';
import * as api from '../api.js';
import { toast, notifyError } from '../ui/toast.js';
import { confirmable } from '../ui/confirm.js';

/* ======================================================== focus plumbing == */

const FOCUSABLE = [
  'a[href]',
  'button:not([disabled])',
  'input:not([disabled])',
  'select:not([disabled])',
  'textarea:not([disabled])',
  '[tabindex]:not([tabindex="-1"])',
].join(', ');

function focusables(root) {
  return Array.from(root.querySelectorAll(FOCUSABLE))
    .filter((el) => !el.hasAttribute('hidden') && el.getClientRects().length > 0);
}

/* ============================================================== the sheet = */

/**
 * One `<dialog>`, built, shown modal, and resolved when it closes.
 *
 * `showModal()` already traps focus and already answers Escape, and this adds
 * the three things it does not guarantee across browsers or across a dialog
 * whose contents change while it is open:
 *
 *   * focus MOVES IN, to a named control, after the first paint;
 *   * Tab and Shift+Tab wrap inside the dialog even after a control has been
 *     enabled, disabled or replaced;
 *   * focus RETURNS to the trigger, and the trigger is remembered explicitly
 *     rather than left to the UA — a trigger that was re-rendered away while
 *     the dialog was open is checked with `isConnected` before it is focused.
 *
 * While a request is in flight the dialog refuses to close. Escaping out of a
 * half-finished POST would leave the workspace and the screen disagreeing about
 * whether a project exists.
 */
function sheet({ cls, width, title, hint, trigger }) {
  const prevFocus = (trigger && trigger.isConnected) ? trigger
    : (document.activeElement instanceof HTMLElement ? document.activeElement : null);

  const closeBtn = h('button.btn.ghost.sm.cp-close', {
    type: 'button',
    title: 'Close without changing anything (Escape)',
    'aria-label': 'Close without changing anything',
  }, 'close');

  const bodyEl = h('div.cp-body', {
    style: {
      padding: 'var(--sp-6) var(--sp-7) var(--sp-7)',
      overflowY: 'auto',
      // The shared `dialog` rule in css/components.css already caps the whole
      // sheet; this caps the SCROLLING part of it so the head stays put.
      maxHeight: 'calc(100vh - var(--h-bar) - var(--h-status) - var(--h-head) - 2 * var(--sp-5))',
      // Every sentence in these dialogs is set at --fs-sm, one step above the
      // --fs-disclosure floor named in css/tokens.css. It matters below: the
      // caveats are lifted to match, so that no number in here is ever set
      // larger than the sentence that qualifies it.
      fontSize: 'var(--fs-sm)',
      lineHeight: 'var(--lh-body)',
      color: 'var(--fg-1)',
    },
  });

  const dlg = h('dialog.cp-dlg' + (cls ? '.' + cls : ''), {
    'aria-label': title,
  },
    h('div.panel-head.cp-head',
      h('h2', title),
      hint ? h('span.hint', hint) : null,
      h('span.grow'),
      closeBtn),
    bodyEl);
  // Never set `display` on the <dialog> itself: an inline display would beat
  // the UA's `dialog:not([open]) { display: none }` and leave a closed dialog
  // painted on the page.
  dlg.style.width = width || 'min(520px, 94vw)';

  let busy = false;
  let done = false;
  let result = null;
  let resolve = null;
  const closed = new Promise((r) => { resolve = r; });

  function cleanup() {
    if (done) return;
    done = true;
    dlg.remove();
    if (prevFocus && prevFocus.isConnected) {
      try { prevFocus.focus(); } catch { /* a detached or unfocusable trigger */ }
    }
    resolve(result);
  }

  dlg.addEventListener('close', cleanup);
  dlg.addEventListener('cancel', (e) => {
    // Escape. Refused mid-request, for the reason in the doc comment.
    if (busy) { e.preventDefault(); return; }
    result = null;
  });
  dlg.addEventListener('keydown', (e) => {
    if (e.key !== 'Tab') return;
    const items = focusables(dlg);
    if (!items.length) return;
    const first = items[0];
    const last = items[items.length - 1];
    const active = document.activeElement;
    if (e.shiftKey && (active === first || !dlg.contains(active))) {
      e.preventDefault(); last.focus();
    } else if (!e.shiftKey && (active === last || !dlg.contains(active))) {
      e.preventDefault(); first.focus();
    }
  });
  closeBtn.addEventListener('click', () => { if (!busy) finish(null); });

  function finish(r) {
    result = r;
    if (dlg.open) dlg.close();
    else cleanup();
  }

  function setBusy(on) {
    busy = !!on;
    closeBtn.disabled = busy;
    if (busy) dlg.setAttribute('aria-busy', 'true');
    else dlg.removeAttribute('aria-busy');
  }

  function open(initialFocus) {
    document.body.appendChild(dlg);
    dlg.showModal();
    // After the first paint: a control inside a dialog that has only just
    // entered the top layer is not reliably focusable in the same task.
    requestAnimationFrame(() => {
      const target = (initialFocus && initialFocus.isConnected)
        ? initialFocus : focusables(dlg)[0];
      if (target) { try { target.focus(); } catch { /* ignore */ } }
      if (target && target.select && target.dataset.cpSelect === '1') target.select();
    });
    return closed;
  }

  return { dlg, body: bodyEl, open, finish, setBusy, isBusy: () => busy };
}

/* ============================================================== fields ==== */

let fieldSeq = 0;

/**
 * A label, a control, an always-present explanation and an error slot, laid out
 * on the existing `.fgrid` (css/components.css: a 96px label column, with
 * `.field-err` already spanning `1 / -1`).
 *
 * The error carries `role="alert"` and is wired to the control by
 * `aria-describedby` + `aria-invalid`, the same pairing js/ui/compound.js and
 * js/ui/tissue.js already use. Copied on purpose: a screen-reader user who
 * tabs back to a rejected field hears why it was rejected.
 */
function field(label, control, { note, title } = {}) {
  const id = control.id || ('cp-f' + (++fieldSeq));
  control.id = id;
  const errId = id + '-err';
  const errEl = h('div.field-err.cp-err.hidden', { id: errId, role: 'alert' });
  control.setAttribute('aria-describedby', errId);
  // `.field-err` and `.field-note` are both --fs-micro in css/components.css,
  // and css/tokens.css names the rule they break here: "ANY ELEMENT CARRYING A
  // SENTENCE is set at --fs-disclosure or larger. […] --fs-micro keeps
  // everything it is actually for: […] labels and machine tokens of a few
  // characters, read as marks rather than sentences." In the Studio's rails
  // those two carry marks; in these dialogs they carry sentences — the error
  // is a refusal that has to say what to do next — so they are lifted. The
  // error goes to --fs-sm, the size of the input text it is correcting; the
  // note, which qualifies nothing, sits on the floor.
  errEl.style.fontSize = 'var(--fs-sm)';

  const labelEl = h('label.lbl.cp-label', { for: id, title: title || '' }, label);
  const noteEl = note ? h('div.field-note.cp-note', note) : null;
  if (noteEl) noteEl.style.fontSize = 'var(--fs-disclosure)';
  const holder = h('div.cp-control', control, noteEl);

  const api_ = {
    id,
    control,
    nodes: [labelEl, holder, errEl],
    fail(message) {
      errEl.textContent = message;
      errEl.classList.remove('hidden');
      control.setAttribute('aria-invalid', 'true');
      try { control.focus(); } catch { /* ignore */ }
      if (control.select) control.select();
    },
    clear() {
      errEl.textContent = '';
      errEl.classList.add('hidden');
      control.removeAttribute('aria-invalid');
    },
  };
  control.addEventListener('input', () => api_.clear());
  control.addEventListener('change', () => api_.clear());
  return api_;
}

/* ============================================================== errors ==== */

/**
 * The form-level error region. See the header: a toast raised over a modal
 * dialog is painted under the backdrop and is never seen, so a failure that
 * happens while one of these is open has to be reported inside it.
 *
 * It says the same three things `js/ui/toast.js` says — what failed, what to do
 * about it, and the route plus status to quote — because those are the three
 * things a failure has to say. What it does not do is guess: an unclassified
 * failure is captioned neutrally rather than blamed on anyone.
 */
function errorRegion() {
  const el = h('div.caveat.crit.cp-form-err.hidden', { role: 'alert' });
  // Lifted from --fs-disclosure to the size of the dialog's own prose. The
  // rule is that a caveat is never SMALLER than what it qualifies; it may be
  // larger, and in a sheet whose body runs at --fs-sm, --fs-xs would be
  // smaller than the sentence it is correcting.
  el.style.fontSize = 'var(--fs-sm)';
  const span = h('span');
  el.appendChild(span);

  return {
    el,
    clear() { clear(span); el.classList.add('hidden'); },
    show(route, err) {
      const status = (err && typeof err.status === 'number') ? err.status : null;
      const said = (err && (err.message || String(err))) || 'no detail given';
      let remedy;
      if (status === 0 || status === null) {
        remedy = 'The server did not answer. Check that it is running — '
          + './.venv/bin/python -m app.server — and try again. Nothing was changed.';
      } else if (status === 401 || status === 403) {
        remedy = 'This server rejected the key this browser is sending. Open '
          + 'account in the top bar to set it, then try again. Nothing was changed.';
      } else if (status === 404) {
        // Two different things produce a 404 here and the user cannot tell
        // them apart from the outside, so both are named rather than one
        // guessed at: the project is gone, or the route was never built.
        remedy = 'Either this project no longer exists — reload the workspace — '
          + 'or this server predates the route (docs/WORKSPACE_CONTRACT.md §1 '
          + 'adds it) and needs updating. Nothing was changed.';
      } else if (status === 409) {
        remedy = 'The server refused this as a conflict. Reload the workspace — '
          + 'someone or something else may have changed this project already.';
      } else if (status >= 500) {
        remedy = 'The server failed while handling this. Whether anything '
          + 'changed is not knowable from here: reload the workspace and look '
          + 'before retrying.';
      } else {
        remedy = 'The server refused this request. The detail below is what it said.';
      }
      clear(span);
      span.appendChild(h('b', 'FAILED: '));
      span.appendChild(document.createTextNode(remedy + ' '));
      span.appendChild(h('span.mono.dim.cp-tech',
        `${route}${status === null ? '' : ' → ' + status} · ${said}`));
      el.classList.remove('hidden');
    },
  };
}

/**
 * The server states a field error as `"<field>: <what is wrong>"` —
 * app/server/projects.py raises `ProjectError("name: required, a non-empty
 * project name")` and `("name: must not be empty")`, and api.js's `detailOf`
 * carries FastAPI's own `loc`-prefixed messages in the same shape. So a 400
 * whose message names a field we have on screen is attached to that field
 * instead of dumped into the form-level region.
 *
 * @returns {boolean} true if it was attached to a field.
 */
function attachFieldError(err, fields) {
  const status = (err && typeof err.status === 'number') ? err.status : null;
  if (status !== 400 && status !== 422) return false;
  const said = (err && err.message) || '';
  const m = /^\s*([a-z_]+)\s*:\s*(.+)$/i.exec(said);
  if (!m) return false;
  const target = fields[m[1].toLowerCase()];
  if (!target) return false;
  target.fail(m[2].trim() + ' — the server refused it, so nothing was changed.');
  return true;
}

/* ============================================================== counting == */

/**
 * What a project holds, from whichever of the two shapes we were handed.
 *
 * `GET /api/projects` returns a SUMMARY (`n_specimens`, `n_runs`); `GET
 * /api/projects/{id}` returns the DOCUMENT (`specimens[].run_ids`). They do not
 * answer the same question: `n_runs` in the summary is the number of run
 * REFERENCES — `sum(len(s["run_ids"]))` in app/server/projects.py — and two
 * specimens are allowed to reference the same run, so it is not the number of
 * distinct runs. `runs_retained` in the delete response IS distinct
 * (`sorted(set(run_ids))`). Conflating them would have the dialog promise to
 * retain more runs than exist.
 */
function counts(p) {
  const specimens = Array.isArray(p && p.specimens) ? p.specimens : null;
  if (specimens) {
    const refs = specimens.flatMap((s) => (s && s.run_ids) || []);
    return {
      exact: true,
      nSpecimens: specimens.length,
      nRefs: refs.length,
      nDistinct: new Set(refs).size,
      sections: specimens.map((s) => s && s.section).filter(Boolean),
    };
  }
  const nRefs = Number.isFinite(p && p.n_runs) ? p.n_runs : null;
  const nSpecimens = Number.isFinite(p && p.n_specimens) ? p.n_specimens : null;
  return {
    exact: false,
    nSpecimens,
    nRefs,
    nDistinct: null,
    sections: Array.isArray(p && p.sections) ? p.sections.filter(Boolean) : [],
  };
}

const plural = (n, one, many) => `${int(n)} ${n === 1 ? one : (many || one + 's')}`;

/* ============================================================= refusals === */

function offline() { return !!(api.api && api.api.stub); }

/**
 * Every project route rejects in stub mode — js/mock.js has no route for any of
 * them (`mock: no route for POST /api/projects`), deliberately, because a
 * fabricated study on screen is worse than no study. Refusing here, before the
 * form, turns that into an answer instead of a 404 after the user has typed.
 */
function refusal() {
  const el = h('div.caveat.crit.cp-refusal', h('span',
    h('b', 'THE BACKEND IS NOT ANSWERING. '),
    'Projects are real state on the server — specimens, and the run ids that '
    + 'belong to each. js/mock.js does not synthesise one, so there is nothing '
    + 'here to create, rename, duplicate or delete until the server answers. '
    + 'Start it with ',
    h('span.mono', './.venv/bin/python -m app.server'),
    ' and reload this page.'));
  el.style.fontSize = 'var(--fs-sm)';
  return el;
}

/* ============================================================== resolve === */

/**
 * The three editing dialogs take a project. Accept an id as well, because a
 * caller holding only a row's id should not have to learn which shape each
 * dialog wants — and fetch the DOCUMENT when given one, since the document is
 * what carries the exact per-specimen run ids that the delete dialog states.
 *
 * This is the one failure path with no dialog on screen, so it is the one that
 * belongs in a toast: `notifyError` classifies it and says what to do next.
 */
async function resolveProject(p) {
  if (p && typeof p === 'object' && p.id) return p;
  const id = typeof p === 'string' ? p : null;
  if (!id) throw new Error('createProject.js: a project or a project id is required');
  try {
    return await api.getProject(id);
  } catch (err) {
    notifyError(`GET /api/projects/${id} failed`, err);
    return null;
  }
}

/* ================================================================ create == */

function suggestName(projects) {
  const d = new Date();
  const p2 = (n) => String(n).padStart(2, '0');
  const base = `Study ${d.getFullYear()}-${p2(d.getMonth() + 1)}-${p2(d.getDate())}`;
  const taken = new Set((projects || [])
    .map((x) => String((x && x.name) || '').trim().toLowerCase()));
  if (!taken.has(base.toLowerCase())) return base;
  for (let i = 2; i < 200; i++) {
    const n = `${base} (${i})`;
    if (!taken.has(n.toLowerCase())) return n;
  }
  return base;
}

/**
 * Create a project.
 *
 * THIS IS THE FLOW THAT DECIDES WHETHER ANY OF THE REST HAPPENS. A first-time
 * visitor's first act is creating a project, and a form that reads as an empty
 * database row — id, name, notes, Save — teaches them that this is filing, not
 * work. So the dialog arrives already answering itself:
 *
 *   * the name is PRE-FILLED with a suggestion that does not collide with an
 *     existing one, and selected, so typing replaces it and not typing is a
 *     valid choice;
 *   * a first specimen is PRE-SELECTED — the first section this server can
 *     actually run — so the project is born with something in it;
 *   * the footer states, in a sentence that updates as the form changes,
 *     exactly what the button will do, including what it will NOT do: no run
 *     is started, and nothing is simulated by creating a project.
 *
 * One press of `Create project` therefore produces a project with a specimen
 * bound to a real section, ready to open in the Studio.
 *
 * @param {object}  [o]
 * @param {Element} [o.trigger]   the control that opened this; focus returns.
 * @param {Array}   [o.projects]  the workspace's current projects, for the
 *   name suggestion and the duplicate-name note. Optional.
 * @param {Array}   [o.sections]  `GET /api/sections`. Fetched if absent.
 * @returns {Promise<{action:'create', project:object, specimen:object|null}|null>}
 */
export async function openCreateProject(o = {}) {
  const s = sheet({
    cls: 'cp-dlg-create',
    width: 'min(560px, 94vw)',
    title: 'New project',
    hint: 'a study that groups runs by specimen',
    trigger: o.trigger,
  });

  if (offline()) {
    s.body.appendChild(refusal());
    s.body.appendChild(h('div.btn-row.cp-foot', { style: { marginTop: 'var(--sp-6)' } },
      h('button.btn.cp-cancel', { type: 'button', onclick: () => s.finish(null) }, 'Close')));
    return s.open();
  }

  let sections = Array.isArray(o.sections) ? o.sections : null;
  if (!sections) {
    // `GET /api/sections` is a catalogue, not project state, so a failure here
    // is not a reason to refuse the create — it costs the default specimen and
    // nothing else, and the form says so where the control would have been.
    try { sections = await api.getSections(); } catch { sections = []; }
  }
  const available = (sections || []).filter((x) => x && x.available);

  /* --- name ------------------------------------------------------------- */

  const nameIn = h('input.inp.cp-name', {
    id: 'cp-create-name',
    type: 'text',
    maxlength: '200',
    autocomplete: 'off',
    spellcheck: 'false',
    required: true,
    value: suggestName(o.projects),
    dataset: { cpSelect: '1' },
  });
  const nameField = field('Name', nameIn, {
    title: 'How you will recognise this study in the workspace. '
      + 'The id (p_00001, …) is assigned by the server and is the identity.',
    note: 'Suggested — type over it. Names do not have to be unique; the id is '
      + 'the identity.',
  });
  const nameNote = nameField.nodes[1].querySelector('.cp-note');

  function checkNameCollision() {
    const v = nameIn.value.trim().toLowerCase();
    const clash = (o.projects || []).find(
      (p) => p && String(p.name || '').trim().toLowerCase() === v && v);
    if (!nameNote) return;
    nameNote.textContent = clash
      ? `${clash.id} already uses this name. That is allowed — the id is the `
        + 'identity — but two rows reading the same is a workspace you have to '
        + 'hover to navigate.'
      : 'Suggested — type over it. Names do not have to be unique; the id is '
        + 'the identity.';
  }
  nameIn.addEventListener('input', checkNameCollision);

  /* --- notes ------------------------------------------------------------ */

  const notesIn = h('textarea.inp.cp-notes', {
    id: 'cp-create-notes',
    rows: '3',
    maxlength: '2000',
    spellcheck: 'true',
    placeholder: 'optional — what question this study is meant to answer',
  });
  notesIn.style.height = 'auto';
  notesIn.style.minHeight = '46px';
  notesIn.style.padding = 'var(--sp-3)';
  const notesField = field('Notes', notesIn, {
    title: 'Free text kept with the project. Nothing reads it but you.',
    note: 'Optional. Kept with the project; no part of the model reads it.',
  });

  /* --- the useful default: a first specimen ----------------------------- */

  const specChk = h('input', {
    type: 'checkbox',
    id: 'cp-create-spec',
    checked: available.length > 0,
    disabled: available.length === 0,
  });
  const sectionSel = h('select.inp.cp-section', {
    id: 'cp-create-section',
    disabled: available.length === 0,
    'aria-label': 'Section to bind as the first specimen',
  }, available.length
    ? available.map((x, i) => h('option', {
        value: x.id, selected: i === 0, title: x.panel || '',
      }, `${x.id} — ${x.label}`))
    : [h('option', { value: '' }, 'no section on this server is ready')]);

  const specField = field('First specimen', sectionSel, {
    title: 'A specimen binds one tissue section to this project. '
      + 'Sections are registered on the server by filesystem path.',
    note: available.length
      ? 'A specimen binds one section this server already has to the project. '
        + 'Nothing is copied and nothing is run; you pick the therapy in the Studio.'
      : 'No section on this server is ready to run, so a specimen cannot be '
        + 'bound yet. The project is still created — add a specimen from the '
        + 'Studio once a section is registered.',
  });
  // The label points at the select, but the ROW is a checkbox above it, so the
  // row goes in first and then adopts the select out of the plain holder that
  // `field()` built. Built in this order on purpose: `field()` moves the
  // control into its own holder, so a row assembled BEFORE the call would have
  // the select taken back out of it again.
  const specRow = h('div.cp-spec',
    h('label.chk.cp-spec-chk', { for: 'cp-create-spec' }, specChk,
      h('span', 'Start with a specimen')));
  specRow.style.display = 'grid';
  specRow.style.gap = 'var(--sp-3)';
  specField.nodes[1].insertBefore(specRow, specField.nodes[1].firstChild);
  specRow.appendChild(sectionSel);
  sectionSel.disabled = available.length === 0 || !specChk.checked;
  specChk.addEventListener('change', () => {
    sectionSel.disabled = !specChk.checked || available.length === 0;
    paintSummary();
  });
  sectionSel.addEventListener('change', paintSummary);

  /* --- what the button will do ------------------------------------------ */

  const summaryEl = h('div.cp-summary');
  summaryEl.style.marginTop = 'var(--sp-6)';
  summaryEl.style.color = 'var(--fg-2)';

  function paintSummary() {
    const withSpec = specChk.checked && !!sectionSel.value;
    const sid = sectionSel.value;
    clear(summaryEl);
    summaryEl.appendChild(document.createTextNode(
      withSpec
        ? 'Creates a project holding 1 specimen — '
        : 'Creates an empty project. '));
    if (withSpec) {
      summaryEl.appendChild(h('span.mono.strong', sid));
      summaryEl.appendChild(document.createTextNode('. '));
    }
    summaryEl.appendChild(document.createTextNode(
      'No run is started and nothing is simulated: the server assigns an id '
      + '(p_…) and you choose the therapy in the Studio.'));
  }
  paintSummary();

  /* --- assembly ---------------------------------------------------------- */

  const errors = errorRegion();
  const createBtn = h('button.btn.primary.cp-submit', { type: 'submit' }, 'Create project');
  const cancelBtn = h('button.btn.cp-cancel', {
    type: 'button', onclick: () => s.finish(null),
  }, 'Cancel');
  const footEl = h('div.btn-row.cp-foot', createBtn, cancelBtn);
  footEl.style.marginTop = 'var(--sp-6)';
  errors.el.style.marginTop = 'var(--sp-5)';

  const form = h('form.cp-form', { novalidate: true },
    h('div.fgrid.cp-grid',
      nameField.nodes, notesField.nodes, specField.nodes),
    summaryEl,
    errors.el,
    footEl);

  form.addEventListener('submit', async (e) => {
    e.preventDefault();
    if (s.isBusy()) return;
    errors.clear();
    nameField.clear();

    const name = nameIn.value.trim();
    if (!name) {
      nameField.fail('Type a name for the study — it is how you will find it in '
        + 'the workspace. The id is assigned by the server, not by you.');
      return;
    }

    s.setBusy(true);
    createBtn.disabled = true;
    cancelBtn.disabled = true;
    const restore = () => {
      s.setBusy(false);
      createBtn.disabled = false;
      cancelBtn.disabled = false;
      createBtn.textContent = 'Create project';
    };

    createBtn.textContent = 'Creating…';
    let project;
    try {
      project = await api.createProject(name, notesIn.value.trim());
    } catch (err) {
      restore();
      if (!attachFieldError(err, { name: nameField, notes: notesField })) {
        errors.show('POST /api/projects', err);
      }
      return;
    }

    let specimen = null;
    if (specChk.checked && sectionSel.value) {
      createBtn.textContent = 'Adding specimen…';
      try {
        specimen = await api.addSpecimen(project.id, sectionSel.value, null);
      } catch (err) {
        // The project EXISTS. Reporting this as a failed create would send the
        // user back to make a second one.
        restore();
        s.finish({ action: 'create', project, specimen: null });
        toast('warn', `${project.id} created, without a specimen`,
          `The project exists. Binding ${sectionSel.value} to it failed: `
          + `${(err && err.message) || 'no detail'}. Add a specimen from the `
          + 'Studio, or from the projects dialog inside it.', 12000);
        return;
      }
    }

    restore();
    s.finish({ action: 'create', project, specimen });
    toast('ok', `${project.id} created`,
      specimen
        ? `${project.name} · 1 specimen (${specimen.id} — ${sectionSel.value}). `
          + 'No run has been started.'
        : `${project.name} · no specimen yet. No run has been started.`,
      5000);
  });

  s.body.appendChild(form);
  return s.open(nameIn);
}

/* ================================================================ rename == */

/**
 * Edit a project's name and its notes.
 *
 * Exported under both names on purpose. The action a workspace row offers is
 * `rename`, which is what the task is called; the dialog also carries the notes
 * field, because the notes are the only thing on screen that says what the name
 * MEANS and renaming a study without them is renaming it blind.
 * UI_ARCHITECTURE.md §5 — "a label must not describe behaviour the code does
 * not have" — is why the dialog's own title is `Edit project`, not `Rename`,
 * and why the recommended trigger label is `edit`.
 *
 * @param {object|string} project  a Project, a ProjectSummary, or a project id.
 * @param {object} [o]
 * @param {Element} [o.trigger]
 * @returns {Promise<{action:'rename', project:object}|null>}
 */
export async function openEditProject(project, o = {}) {
  const p = await resolveProject(project);
  if (!p) return null;

  const s = sheet({
    cls: 'cp-dlg-edit',
    width: 'min(520px, 94vw)',
    title: 'Edit project',
    hint: 'name and notes',
    trigger: o.trigger,
  });

  if (offline()) {
    s.body.appendChild(refusal());
    s.body.appendChild(h('div.btn-row.cp-foot', { style: { marginTop: 'var(--sp-6)' } },
      h('button.btn.cp-cancel', { type: 'button', onclick: () => s.finish(null) }, 'Close')));
    return s.open();
  }

  const nameIn = h('input.inp.cp-name', {
    id: 'cp-edit-name',
    type: 'text',
    maxlength: '200',
    autocomplete: 'off',
    spellcheck: 'false',
    value: p.name || '',
    dataset: { cpSelect: '1' },
  });
  const nameField = field('Name', nameIn, {
    title: 'The name only. The id never changes and nothing that references '
      + 'this project is affected.',
    note: 'The id stays ' + p.id + '. Renaming touches nothing else: not the '
      + 'specimens, not the run associations, not a single run.',
  });

  const notesIn = h('textarea.inp.cp-notes', {
    id: 'cp-edit-notes', rows: '4', maxlength: '2000', spellcheck: 'true',
    placeholder: 'optional — what question this study is meant to answer',
  });
  notesIn.value = p.notes || '';
  notesIn.style.height = 'auto';
  notesIn.style.minHeight = '60px';
  notesIn.style.padding = 'var(--sp-3)';
  const notesField = field('Notes', notesIn, {
    note: 'Optional. Kept with the project; no part of the model reads it.',
  });

  const c = counts(p);
  const identEl = h('div.cp-ident',
    h('span.mono.strong', p.id),
    h('span.dim', ' · '),
    h('span.mono.dim',
      `${plural(c.nSpecimens, 'specimen')} · `
      + `${plural(c.nRefs, 'run reference')} · created ${ts(p.created_at)}`));
  identEl.style.marginBottom = 'var(--sp-5)';

  const errors = errorRegion();
  const saveBtn = h('button.btn.primary.cp-submit', { type: 'submit' }, 'Save');
  const cancelBtn = h('button.btn.cp-cancel', {
    type: 'button', onclick: () => s.finish(null),
  }, 'Cancel');

  const footEl = h('div.btn-row.cp-foot', saveBtn, cancelBtn);
  footEl.style.marginTop = 'var(--sp-6)';
  errors.el.style.marginTop = 'var(--sp-5)';

  const form = h('form.cp-form', { novalidate: true },
    identEl,
    h('div.fgrid.cp-grid', nameField.nodes, notesField.nodes),
    errors.el,
    footEl);

  form.addEventListener('submit', async (e) => {
    e.preventDefault();
    if (s.isBusy()) return;
    errors.clear();
    nameField.clear();

    const name = nameIn.value.trim();
    if (!name) {
      nameField.fail('A project cannot have an empty name. Type one, or press '
        + 'Cancel to leave it as ' + JSON.stringify(p.name || '') + '.');
      return;
    }
    const notes = notesIn.value;
    if (name === (p.name || '') && notes === (p.notes || '')) {
      // Not an error, and not worth a round trip that would stamp updated_at
      // and reorder the workspace for a change nobody made.
      s.finish(null);
      return;
    }

    s.setBusy(true);
    saveBtn.disabled = true; cancelBtn.disabled = true;
    saveBtn.textContent = 'Saving…';
    let updated;
    try {
      updated = await api.patchProject(p.id, { name, notes });
    } catch (err) {
      s.setBusy(false);
      saveBtn.disabled = false; cancelBtn.disabled = false;
      saveBtn.textContent = 'Save';
      if (!attachFieldError(err, { name: nameField, notes: notesField })) {
        errors.show(`PATCH /api/projects/${p.id}`, err);
      }
      return;
    }
    s.setBusy(false);
    s.finish({ action: 'rename', project: updated });
    toast('ok', `${p.id} updated`,
      `${updated.name} · specimens, run associations and runs are unchanged.`, 4000);
  });

  s.body.appendChild(form);
  return s.open(nameIn);
}

/** The workspace row's verb. Same dialog; see `openEditProject`. */
export const openRenameProject = openEditProject;

/* ============================================================= duplicate == */

/**
 * Duplicate a project.
 *
 * The one thing this dialog exists to prevent is the belief that duplicating a
 * study re-runs it. Contract §1: the route "copies specimens and their run ids.
 * Runs are NOT re-executed and NOT copied — a run belongs to the store, and two
 * specimens may reference the same run id." A user who thinks otherwise will
 * sit waiting for results that were never queued, and will later believe two
 * independent runs agree when they are one run counted twice.
 *
 * @param {object|string} project
 * @param {object} [o]
 * @param {Element} [o.trigger]
 * @returns {Promise<{action:'duplicate', project:object, source:object}|null>}
 */
export async function openDuplicateProject(project, o = {}) {
  const p = await resolveProject(project);
  if (!p) return null;

  const s = sheet({
    cls: 'cp-dlg-dup',
    width: 'min(560px, 94vw)',
    title: 'Duplicate project',
    hint: 'copies the specimen list, not the runs',
    trigger: o.trigger,
  });

  if (offline()) {
    s.body.appendChild(refusal());
    s.body.appendChild(h('div.btn-row.cp-foot', { style: { marginTop: 'var(--sp-6)' } },
      h('button.btn.cp-cancel', { type: 'button', onclick: () => s.finish(null) }, 'Close')));
    return s.open();
  }

  const c = counts(p);

  const nameIn = h('input.inp.cp-name', {
    id: 'cp-dup-name',
    type: 'text',
    maxlength: '200',
    autocomplete: 'off',
    spellcheck: 'false',
    value: `${p.name || p.id} (copy)`,
    dataset: { cpSelect: '1' },
  });
  const nameField = field('New name', nameIn, {
    note: 'The copy gets its own id from the server. ' + p.id + ' is not changed '
      + 'in any way by this.',
  });

  const identEl = h('div.cp-ident',
    h('span.lbl', 'copying'), ' ',
    h('span.mono.strong', p.id), ' ',
    h('span', p.name || NULL_GLYPH));
  identEl.style.marginBottom = 'var(--sp-5)';

  /* What is copied, and what is emphatically not. Both at the dialog's prose
     size, so neither sentence is set smaller than the counts inside it. */
  const copies = h('div.cp-copies',
    h('div.cp-copies-head', 'The copy gets'),
    h('ul.cp-list',
      h('li', `${plural(c.nSpecimens, 'specimen')} — the same sections, bound again.`),
      h('li', `${plural(c.nRefs, 'run reference')}`
        + (c.exact && c.nDistinct !== c.nRefs
          ? ` to ${plural(c.nDistinct, 'distinct run')}.`
          : '.')
        + ' The references are copied; the runs they point at are not.')));
  copies.style.marginTop = 'var(--sp-5)';
  const head = copies.querySelector('.cp-copies-head');
  head.style.color = 'var(--fg-2)';
  head.style.marginBottom = 'var(--sp-2)';
  const ul = copies.querySelector('.cp-list');
  ul.style.margin = '0';
  ul.style.paddingLeft = 'var(--sp-7)';

  const notCopied = h('div.caveat.crit.cp-keep', h('span',
    h('b', 'RUNS ARE NEITHER RE-EXECUTED NOR DUPLICATED. '),
    'A run lives in app/runs/ and a project holds only its id, so after this '
    + 'both projects point at the SAME runs. Nothing is queued, nothing is '
    + 'simulated, and the copy shows finished results immediately because they '
    + 'are the original’s results. Two projects agreeing about a run is '
    + 'one run counted twice, not two runs that agree.'));
  notCopied.style.fontSize = 'var(--fs-sm)';
  notCopied.style.marginTop = 'var(--sp-5)';

  const errors = errorRegion();
  const dupBtn = h('button.btn.primary.cp-submit', { type: 'submit' }, 'Duplicate');
  const cancelBtn = h('button.btn.cp-cancel', {
    type: 'button', onclick: () => s.finish(null),
  }, 'Cancel');

  const footEl = h('div.btn-row.cp-foot', dupBtn, cancelBtn);
  footEl.style.marginTop = 'var(--sp-6)';
  errors.el.style.marginTop = 'var(--sp-5)';

  const form = h('form.cp-form', { novalidate: true },
    identEl,
    h('div.fgrid.cp-grid', nameField.nodes),
    copies,
    notCopied,
    errors.el,
    footEl);

  form.addEventListener('submit', async (e) => {
    e.preventDefault();
    if (s.isBusy()) return;
    errors.clear();
    nameField.clear();

    const name = nameIn.value.trim();
    if (!name) {
      nameField.fail('Name the copy, so the workspace does not show you two '
        + 'rows you cannot tell apart.');
      return;
    }

    s.setBusy(true);
    dupBtn.disabled = true; cancelBtn.disabled = true;
    dupBtn.textContent = 'Duplicating…';
    let copy;
    try {
      copy = await api.duplicateProject(p.id, name);
    } catch (err) {
      s.setBusy(false);
      dupBtn.disabled = false; cancelBtn.disabled = false;
      dupBtn.textContent = 'Duplicate';
      if (!attachFieldError(err, { name: nameField })) {
        errors.show(`POST /api/projects/${p.id}/duplicate`, err);
      }
      return;
    }
    s.setBusy(false);
    s.finish({ action: 'duplicate', project: copy, source: p });
    toast('ok', `${copy.id} is a copy of ${p.id}`,
      `${copy.name} · ${plural(c.nSpecimens, 'specimen')} and their run `
      + 'references were copied. No run was re-executed and no run was '
      + 'duplicated — both projects point at the same run ids.', 7000);
  });

  s.body.appendChild(form);
  return s.open(nameIn);
}

/* ================================================================ delete == */

/**
 * Delete a project.
 *
 * WORKSPACE_CONTRACT.md §2, in full, because this dialog is the only place the
 * product says it out loud: "A project owns nothing but ids. […] The
 * consequence, which the UI must not hide: deleting a project does not delete
 * its runs. The home page says so at the point of deletion, not in a help
 * file."
 *
 * Three decisions follow from that:
 *
 *   * The dialog reads the counts off the project DOCUMENT (`specimens[]
 *     .run_ids`), not off a summary row, because the summary's `n_runs` counts
 *     REFERENCES and the server's `runs_retained` is a SET. Promising to retain
 *     seven runs and then reporting five retained would read as data loss.
 *   * If the count cannot be established, the dialog says so and still offers
 *     the delete — the guarantee does not depend on the number, and a delete
 *     button that disappears when a GET fails is a worse answer than an honest
 *     "not known".
 *   * Nothing in here is set larger than the sentence that qualifies it. The
 *     counts sit inside the same prose as the guarantee, at the same size, and
 *     the kept-runs caveat is lifted from --fs-disclosure to --fs-sm to match.
 *
 * Focus opens on `Cancel`, not on the destructive control, and the destructive
 * control is `confirmable` (js/ui/confirm.js) — armed once, confirmed on a
 * second press, disarmed by Escape, a click elsewhere, or six seconds.
 *
 * @param {object|string} project
 * @param {object} [o]
 * @param {Element} [o.trigger]
 * @param {boolean} [o.armed]  false makes `Delete project` a SINGLE press.
 *   Default true. This exists because js/home/projectCard.js raises its own
 *   inline confirmation, with its own §2 disclosure, BEFORE it calls
 *   `onDelete` — so on the wiring js/home/main.js ships today, deleting a
 *   project is three presses and the guarantee is stated twice. That is one
 *   confirmation too many and one owner too many; whichever module keeps the
 *   disclosure, the other should stop asking. Passing `{armed: false}` here is
 *   the one-word way to resolve it in favour of the card, and dropping the
 *   card's inline confirm is the way to resolve it in favour of this dialog.
 *   Reported rather than decided unilaterally: it spans two territories.
 * @returns {Promise<{action:'delete', id:string, runs_retained:string[]}|null>}
 */
export async function openDeleteProject(project, o = {}) {
  // Deliberately NOT `resolveProject`: the dialog must open promptly even if
  // the document read is slow or fails, so a summary is enough to open with and
  // the exact counts are filled in underneath.
  const seed = (project && typeof project === 'object' && project.id)
    ? project
    : (typeof project === 'string' ? { id: project } : null);
  if (!seed) throw new Error('createProject.js: a project or a project id is required');

  const s = sheet({
    cls: 'cp-dlg-del',
    width: 'min(560px, 94vw)',
    title: 'Delete project',
    hint: 'the project document only',
    trigger: o.trigger,
  });

  if (offline()) {
    s.body.appendChild(refusal());
    s.body.appendChild(h('div.btn-row.cp-foot', { style: { marginTop: 'var(--sp-6)' } },
      h('button.btn.cp-cancel', { type: 'button', onclick: () => s.finish(null) }, 'Close')));
    return s.open();
  }

  const identIdEl = h('span.mono.strong', seed.id);
  const identNameEl = h('span.cp-ident-name', seed.name || '');
  const identEl = h('div.cp-ident', identIdEl, ' ', identNameEl);
  identEl.style.marginBottom = 'var(--sp-5)';

  const losesEl = h('div.cp-lose');
  const keepsEl = h('div.caveat.crit.cp-keep');
  // Lifted from --fs-disclosure to the dialog's own prose size, so that the
  // guarantee is never set smaller than the counts it qualifies. A caveat may
  // be larger than what it corrects; it may never be smaller.
  keepsEl.style.fontSize = 'var(--fs-sm)';
  keepsEl.style.marginTop = 'var(--sp-5)';

  function paint(doc, unknownWhy) {
    const c = counts(doc || seed);
    if (doc && doc.name) identNameEl.textContent = doc.name;
    clear(losesEl);
    losesEl.appendChild(h('div.cp-lose-head', 'Deleted'));
    const ul = h('ul.cp-list',
      h('li', 'The project document — its name, its notes, and its id '
        + `${seed.id}, which is not reused.`),
      h('li', c.nSpecimens === null
        ? 'Its specimen list.'
        : `Its ${plural(c.nSpecimens, 'specimen')} — the bindings between this `
          + 'project and those sections.'),
      h('li', c.nRefs === null
        ? 'The run references those specimens hold.'
        : `The ${plural(c.nRefs, 'run reference')} those specimens hold — the `
          + 'associations, not the runs.'));
    ul.style.margin = '0';
    ul.style.paddingLeft = 'var(--sp-7)';
    losesEl.appendChild(ul);
    const head = losesEl.querySelector('.cp-lose-head');
    head.style.color = 'var(--fg-2)';
    head.style.marginBottom = 'var(--sp-2)';
    losesEl.style.marginTop = 'var(--sp-5)';

    clear(keepsEl);
    const kept = (c.nDistinct === null) ? null : c.nDistinct;
    keepsEl.appendChild(h('span',
      h('b', 'RUNS ARE NOT DELETED. '),
      kept === null
        ? 'Every run this project references stays in app/runs/, which is the '
          + 'single store for run results. '
        : (kept === 0
            ? 'This project references no run, so there is nothing to retain. '
            : `The ${plural(kept, 'run')} this project references stay in `
              + 'app/runs/, which is the single store for run results. '),
      'A project holds run IDS, never results, so deleting one detaches its '
      + 'runs and never deletes them. They stay in the Studio’s run '
      + 'history, and in any other project that references them.',
      unknownWhy
        ? h('span.dim', ' The exact count could not be read from this server '
          + `(${unknownWhy}); the guarantee above does not depend on it.`)
        : null));
  }

  paint(seed, null);

  const errors = errorRegion();
  const cancelBtn = h('button.btn.cp-cancel', {
    type: 'button', onclick: () => s.finish(null),
  }, 'Cancel');
  const deleteBtn = h('button.btn.danger.cp-delete', {
    type: 'button',
    title: `Delete the project document ${seed.id}. Its runs are NOT deleted.`,
    'aria-label': `Delete project ${seed.id}`,
  }, 'Delete project');

  const doDelete = async () => {
    if (s.isBusy()) return;
    errors.clear();
    s.setBusy(true);
    deleteBtn.disabled = true; cancelBtn.disabled = true;
    let out;
    try {
      out = await api.deleteProject(seed.id);
    } catch (err) {
      s.setBusy(false);
      deleteBtn.disabled = false; cancelBtn.disabled = false;
      errors.show(`DELETE /api/projects/${seed.id}`, err);
      return;
    }
    s.setBusy(false);
    const retained = (out && out.runs_retained) || [];
    s.finish({ action: 'delete', id: seed.id, runs_retained: retained });
    toast('ok', `${seed.id} deleted`,
      retained.length
        ? `${plural(retained.length, 'run')} retained in app/runs/ — a project `
          + 'owns ids, not results.'
        : 'It referenced no run, so nothing was retained and nothing was lost '
          + 'from app/runs/.',
      6000);
  };

  if (o.armed === false) {
    deleteBtn.addEventListener('click', doDelete);
  } else {
    confirmable(deleteBtn, {
      confirmLabel: `delete ${seed.id}?`,
      armedTitle: `Deletes the project document ${seed.id} and its specimen list. `
        + 'Its runs stay in app/runs/. Press again to confirm; Escape cancels.',
      onConfirm: doDelete,
    });
  }

  // Cancel first in the DOM as well as in the tab order: the safe choice is the
  // one the keyboard reaches first, and the one Enter lands on when focus opens
  // here.
  const foot = h('div.btn-row.cp-foot', cancelBtn, deleteBtn);
  foot.style.marginTop = 'var(--sp-6)';
  errors.el.style.marginTop = 'var(--sp-5)';

  s.body.appendChild(h('div.cp-del-body',
    identEl, losesEl, keepsEl, errors.el, foot));

  const opened = s.open(cancelBtn);

  // Fill in the exact counts behind the already-open dialog. A summary row
  // carries `n_runs`, which is the number of REFERENCES; only the document
  // gives the distinct runs the guarantee is about.
  if (!Array.isArray(seed.specimens)) {
    api.getProject(seed.id)
      .then((doc) => { if (!s.isBusy()) paint(doc, null); })
      .catch((err) => paint(seed, (err && err.message) || 'the read failed'));
  }

  return opened;
}

/* ============================================================== grouped === */

/**
 * The four dialogs with shared defaults bound once, for a caller that has the
 * workspace in hand and does not want to pass `projects` and `sections` on
 * every call.
 *
 * Each returns the same result object its standalone twin does, or null when
 * the user backed out. The change callback is called only on a REAL change, so
 * a caller can make it the single place it refetches `GET /api/workspace`. Both
 * names are accepted: js/home/main.js calls
 * `createProjectDialogs({mode:'create', onDone: () => location.reload()})` and
 * then `.create()`, and a callback that is named differently by the two sides
 * of a handshake is a callback that silently never fires.
 */
export function createProjectDialogs(defaults = {}) {
  const shared = () => ({
    projects: typeof defaults.projects === 'function' ? defaults.projects() : defaults.projects,
    sections: typeof defaults.sections === 'function' ? defaults.sections() : defaults.sections,
  });
  const changed = defaults.onChanged || defaults.onDone || defaults.onChange || null;
  const announce = async (p) => {
    const r = await p;
    if (r && changed) await changed(r);
    return r;
  };
  return {
    create: (o = {}) => announce(openCreateProject({ ...shared(), ...o })),
    rename: (project, o = {}) => announce(openEditProject(project, o)),
    edit: (project, o = {}) => announce(openEditProject(project, o)),
    duplicate: (project, o = {}) => announce(openDuplicateProject(project, o)),
    remove: (project, o = {}) => announce(openDeleteProject(project, o)),
  };
}
