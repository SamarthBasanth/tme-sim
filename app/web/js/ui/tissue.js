// tissue.js — BRING YOUR OWN TISSUE: point Vivum at a spatial-transcriptomics
// export on the server's disk, get a verdict on what that specimen can and
// cannot support, and — only if it survives every check — a section the picker
// can select.
//
// Why this panel exists at all
// ----------------------------
// The capability existed in `POST /api/sections/scan` and in `vivum ingest`
// before this panel did, and appeared NOWHERE in the Studio: the section picker
// was a read-only list from GET /api/sections. A field tester who will not
// hand-POST JSON therefore could not load their own sample, which is the entire
// value proposition. api.py::ingest_page put a self-contained stopgap form at
// /ingest and says in its own docstring that the destination is "the Studio
// picker itself should learn to ingest". This is that.
//
// Three rules this panel is built to, in order:
//
//  1. AN INGESTED SECTION IS NOT A VALIDATED ONE, AND THAT IS SAID BEFORE
//     ANYTHING IS INGESTED, not after. The built-in specimens are the ones this
//     repo ships measured evidence for; a section ingested here has been scored
//     against nothing. That statement is a permanent block in this panel, it is
//     repeated on the registration result, and js/ui/spec.js repeats it a third
//     time in the picker at the moment such a section is SELECTED — because
//     that is the moment it starts affecting a number. None of the three is
//     collapsible and none is behind a hover.
//
//  2. THE REFUSAL IS THE PRODUCT. A directory that is not ingestable comes back
//     with a reason, and the reason is rendered verbatim — it is the same
//     sentence `vivum ingest` prints, so this surface cannot drift from the
//     CLI. What this panel ADDS is one line of "what to do about it" chosen
//     from the failure class, because "extent 2391947.3 x 1795825.1 um exceeds
//     100000 um" tells a first-time user what is wrong and not what to type
//     next. The remedy line never restates a number the server did not send.
//
//  3. NOTHING IS UPLOADED FROM THIS BROWSER. There is no multipart endpoint and
//     no file picker: the operator names a directory and the SERVER PROCESS
//     reads it off its own disk. That is a security-relevant fact about who can
//     read what, and it is stated on the form rather than left to be inferred
//     from a placeholder that happens to look like a path.
//
// Stub mode: ingesting is REFUSED outright when js/api.js has fallen back to
// the mock, for the same reason js/ui/compound.js refuses to compile there — a
// synthesised capability report would look exactly like a measured one.

import { h, clear } from '../util/dom.js';
import * as api from '../api.js';
import { toast, notifyError } from './toast.js';
import { confirmable } from './confirm.js';
import { caveatNode } from './inspect.js';
import { STANDING_LIMITS } from './provenance.js';
import { int, num } from '../util/fmt.js';

// A full ingest reads and decompresses the whole count matrix, so it is minutes
// on a large section, not seconds — api.js's ordinary 30 s budget aborts it
// mid-load and leaves the operator with a timeout for a request the server is
// still happily working on. This is the deliberate exception.
const INGEST_TIMEOUT_MS = 15 * 60 * 1000;

// The refusal classes the backend actually produces, each matched on a phrase
// that is part of the sentence's meaning rather than its punctuation, so a
// re-wording upstream degrades to "no remedy line" and never to a WRONG one.
// Every remedy is procedural — what to do next — because the facts belong to
// the server's own sentence, which is rendered underneath it in full.
const REMEDIES = [
  [/is not a directory and does not exist|is not a directory/i,
    'Nothing is at that path. It has to be a path on the machine running this '
    + 'Vivum server, not on the computer your browser is running on, and it has '
    + 'to be the directory itself rather than a file inside it.'],
  [/is a CONTAINER, not a sample/i,
    'This directory holds your samples rather than being one. Ingest one of the '
    + 'child directories named above, one at a time.'],
  [/no recognised layout/i,
    'Vivum could not tell what this directory is. Check you are pointing at the '
    + 'sample directory itself — for a Xenium run that is the folder holding '
    + 'cells.csv.gz, for Space Ranger the one holding spatial/ and the matrix, '
    + 'and for anything else a MEX triple beside a coordinates table. If a '
    + 'component is still in its .tar.gz, extract it first.'],
  [/exceeds .* um|not in micrometres|are not in micrometres|almost certainly PIXELS/i,
    'The coordinate units are wrong. This pipeline works in micrometres '
    + 'throughout; convert the centroid columns before ingesting rather than '
    + 'letting the reconstruction run at the wrong scale, where nothing '
    + 'downstream would notice.'],
  [/below the .* floor|cells \(floor/i,
    'There is not enough tissue here to reconstruct. If this is a crop, ingest a '
    + 'larger region of the same section.'],
  [/is a BUILT-IN section/i,
    'That id belongs to a built-in specimen, and Vivum will not let an ingested '
    + 'section shadow one — a published number must keep meaning what it was '
    + 'measured from. Give this sample a different section id below.'],
  [/is already ingested from/i,
    'That id is taken by a section ingested from a different directory. Give '
    + 'this one a different section id, or remove the existing one from the list '
    + 'at the bottom of this panel first.'],
  [/expression matrix is all zeros|negative values/i,
    'The count matrix is not raw counts of the kind this pipeline needs. Ingest '
    + 'the unnormalised matrix your instrument produced.'],
  [/panel gene list unreadable/i,
    'The panel could not be read without unpacking the matrix. Tick "measure '
    + 'detection rates" (which unpacks it) or ingest, and the report will be '
    + 'complete.'],
];

function remedyFor(message) {
  const m = String(message || '');
  for (const [re, text] of REMEDIES) if (re.test(m)) return text;
  return '';
}

/* The four states the ingest gate can be closed in, each as the sentence that
 * is shown beside the button. The gate is unchanged — nothing may be registered
 * until a scan of the CURRENT path has come back without a refusal — and what
 * is new is only that the closed state says which of the four it is in.
 * "Scan first" is a different instruction from "scan AGAIN", and the button
 * could not tell them apart. */
const WHY_UNSCANNED = 'Ingest opens once a scan of this directory comes back '
  + 'without a refusal. Scan first: it reads filenames and the gene panel and '
  + 'registers nothing.';
const WHY_STALE = 'Ingest is closed: the report below describes the directory '
  + 'that was scanned, not the one now in the field. Scan again.';
const WHY_REFUSED = 'Ingest is closed: the server refused this directory. Its '
  + 'own reason is below, in full.';
const WHY_REGISTERED = 'Ingest is closed: this directory is registered. A '
  + 'second press would only collide with the id it just took — point the '
  + 'field at another directory to load one more.';

/** The fields inside the collapsed "Name and depth" block. */
const FOLDED_FIELDS = new Set(['t-id', 't-label', 't-deep']);

/**
 * Is this record spot transcriptomics?
 *
 * The one platform test in this file, because it changes what a count IS
 * (spots, not cells) and it is what app/server/api.py keys its per-cell
 * refusal on — `_SPOT_SOURCES = ("visium",)`, matched against `platform` for an
 * ingested record and `source` for a built-in one. Same key, same word, so the
 * panel and the refusal cannot disagree about which sections are affected.
 */
function isSpot(rec) {
  return String((rec && (rec.platform || rec.source)) || '').toLowerCase() === 'visium';
}

/**
 * The one sentence that distinguishes a section ingested here from a built-in
 * one. Exported because js/ui/spec.js says it again at the point of SELECTION —
 * one string, one claim, three places it cannot be missed.
 */
export const NOT_VALIDATED = (
  'INGESTED — NOT VALIDATED. The built-in specimens are the ones this repo '
  + 'ships measured evidence for, and every published number was measured from '
  + 'them. A section you ingest has been scored against nothing: no benchmark, '
  + 'no scorecard, no comparison against a measurement. The capability report is '
  + 'the only check it has passed, and that is a check on the FILES and the '
  + 'gene panel, not on the biology.'
);

/**
 * The one constraint this panel places on the OPERATOR rather than on the
 * output. Every claim in it is from docs/DEPLOYMENT.md §0, which is the
 * authority for this wording; the point of stating it here is that §0's four
 * checks are about strings this form is asking for, and the moment to read them
 * is before they are typed rather than after they are registered.
 */
export const DEIDENTIFIED_ONLY = (
  'DE-IDENTIFIED TISSUE ONLY. Vivum has one shared secret and no accounts, no '
  + 'audit log, no encryption at rest and no per-record access control, so it '
  + 'must not hold data that identifies a person. The tissue itself never '
  + 'leaves this machine — but three of the strings you are about to type do: '
  + 'GET /api/sections returns the DIRECTORY PATH to every client that can '
  + 'reach this server and writes it to app/runs/_ingested_sections.json, and '
  + 'the SECTION ID and LABEL both default to the directory name, appear in the '
  + 'picker, and put the label in the heading of every exported report. Copy or '
  + 'rename to a study code first. The full checklist, including the '
  + 'identifiers that ride inside a 10x export without Vivum ever reading them, '
  + 'is docs/DEPLOYMENT.md §0.'
);

/**
 * A section's measured scale, with the right noun on the count.
 *
 * `n_cells` on a VISIUM record is a count of SPOTS, not of cells: cache.py's
 * `_ingested_metadata` copies `n_cells_ingested` into `n_spots` for exactly
 * that platform, and tmesim/data/ingest.py states that a Visium row "is a
 * ~55 um spot pooling several cells, NOT a cell". This panel printed all of
 * them as "cells".
 *
 * MODULE SCOPE, AND EXPORTED, since 2026-09-20. It was a closure inside
 * `createTissuePanel`, which meant the one line in the product that says what a
 * specimen actually IS — platform, count, extent, panel size — existed only on
 * the panel an operator opens to ingest their own tissue. The picker on the
 * other rail, where a section is CHOSEN, said nothing but the label until after
 * the run. `js/ui/spec.js` now renders this under the section select, so the
 * consequence of the choice is legible before it is made and the two surfaces
 * cannot word it differently.
 *
 * Takes any GET /api/sections row, built-in or ingested.
 */
export function sectionScaleRow(rec) {
  const spot = isSpot(rec);
  const n = Number.isFinite(rec.n_spots) ? rec.n_spots
    : (Number.isFinite(rec.n_cells) ? rec.n_cells : null);
  const ext = Array.isArray(rec.extent_um) && rec.extent_um.length >= 2
    ? rec.extent_um : null;
  const bits = [];
  if (rec.platform || rec.source) {
    bits.push(h('span.mono', String(rec.platform || rec.source)));
    bits.push(h('span.dim', '·'));
  }
  bits.push(n === null
    ? h('span.dim', 'count not stated')
    : h('span', h('span.mono', int(n)),
      h('span.unit', spot ? 'spots' : 'cells')));
  bits.push(h('span.dim', '·'));
  // The extent used to be `(array || []).map(...).join(' × ') + ' mm'`, which
  // on a record without one printed the string " mm" — a unit with no number
  // in front of it, which is worse than a dash.
  bits.push(ext === null
    ? h('span.dim', 'extent not stated')
    : h('span', h('span.mono', `${num(ext[0] / 1000, 2)} × ${num(ext[1] / 1000, 2)}`),
      h('span.unit', 'mm')));
  if (Number.isFinite(rec.n_genes)) {
    bits.push(h('span.dim', '·'));
    bits.push(h('span', h('span.mono', int(rec.n_genes)),
      h('span.unit', 'genes')));
  }
  return h('div.frow', {
    style: {
      gap: 'var(--sp-3)', flexWrap: 'wrap', padding: 'var(--sp-1) 0',
      fontSize: 'var(--fs-xs)', color: 'var(--fg-2)',
    },
    title: spot
      ? 'A Visium row is a ~55 um spot pooling several cells, so this is a '
        + 'count of spots. The cells a run works on are PLACED from the '
        + 'per-spot deconvolved fractions and do not exist until the section '
        + 'is reconstructed.'
      : 'Measured at ingest time by tmesim.data.ingest and stored with the '
        + 'section; not re-derived on this screen.',
  }, bits);
}

/**
 * What `runnable: false` means, in the operator's terms. Set by
 * ingest_registry.register() on a `mex` bundle: it ingested and it is real, and
 * reconstruct_real accepts 'xenium' or 'visium' only, so there is no route from
 * it to a simulation. Exported for the same reason NOT_VALIDATED is — the
 * picker has to say it too.
 */
export const NOT_RUNNABLE = (
  'INSPECTABLE, NOT RUNNABLE. This section loaded and its capability report is '
  + 'real, but the reconstruction pipeline has no route for its platform: a '
  + 'generic MEX bundle can be listed and inspected and cannot be simulated. It '
  + 'is not offered in the run picker rather than failing after you wait for it.'
);

export function createTissuePanel({ onRegistered }) {
  /* ------------------------------------------------------------- form -- */

  const fieldEls = new Map();

  const pathIn = h('input.inp.mono', {
    id: 't-path', type: 'text', autocomplete: 'off', spellcheck: 'false',
    placeholder: '/data/xenium_run_2026_08/output-XETG00001',
    oninput: () => { clearFieldErrors(); markStale(); paintNameEcho(); },
    onkeydown: (e) => { if (e.key === 'Escape') { e.stopPropagation(); pathIn.blur(); } },
  });
  const idIn = h('input.inp.mono', {
    id: 't-id', type: 'text', autocomplete: 'off', spellcheck: 'false',
    placeholder: 'defaults to the directory name',
    oninput: () => { clearFieldErrors(); markStale(); paintNameEcho(); },
  });
  const labelIn = h('input.inp', {
    id: 't-label', type: 'text', autocomplete: 'off',
    placeholder: 'optional — shown in the section picker',
    oninput: () => { paintNameEcho(); },
  });
  const deepChk = h('input', {
    type: 'checkbox', id: 't-deep', onchange: () => paintNameEcho(),
  });

  const lbl = (text, forId, title) => h('label.lbl', { for: forId, title: title || '' }, text);

  function row(label, forId, control, unitText, title) {
    const holder = unitText ? h('div.with-unit', control, h('span.unit', unitText)) : control;
    const errId = forId + '-err';
    const errEl = h('div.field-err.hidden', { id: errId, role: 'alert' });
    control.setAttribute('aria-describedby', errId);
    fieldEls.set(forId, { control, errEl });
    return [lbl(label, forId, title), holder, errEl];
  }

  /* --------------------------------------------------- naming, folded ---
   *
   * ONE QUESTION ON THIS PANEL, not four. A first-time operator has exactly one
   * thing to supply — WHERE the sample is — and the other three fields all have
   * a real default, set by app/server/ingest_registry.py::register: the id is
   * the directory name slugified (`_make_id`), the label is
   * `f"{platform.title()} · {sid} (ingested)"`, and detection rates are
   * measured anyway on the ingest that matters (`register` implies a full
   * ingest). Three fields with defaults, presented as blanks, read as three
   * more decisions before anything can happen.
   *
   * THE FOLD HIDES THE WIDGETS, NEVER THE FACT. The summary prints the id and
   * the label that WILL be used, recomputed from the path on every keystroke —
   * which is strictly more than the open form said, because the open form
   * showed two empty boxes and the word "defaults to the directory name"
   * without ever naming the string. That matters here more than anywhere else
   * on this panel: DEIDENTIFIED_ONLY, two blocks above, is about exactly these
   * two strings leaving the machine, and it is now possible to READ the string
   * that will leave before pressing anything.
   *
   * It opens itself whenever it stops being a detail: a field error inside it
   * (`fieldError`), or a server remedy that names one of its controls
   * (`showRefusal`).
   */
  const nameEcho = h('span.mono', {
    style: { color: 'var(--fg-1)', fontWeight: '400', textTransform: 'none', letterSpacing: '0' },
  }, 'the directory name');
  const nameFold = h('details.span2.tnames', { id: 'sect-tissue-names' },
    h('summary', {
      style: { cursor: 'pointer', padding: 'var(--sp-2) 0', lineHeight: 'var(--lh-body)' },
      title: 'Section id, label and scan depth. All three have defaults; the '
        + 'id and the label that will actually be used are printed here.',
    },
    h('span.lbl', 'Name and depth'), ' — ', nameEcho),
    h('div.fgrid', { style: { padding: 'var(--sp-3) 0 var(--sp-2)' } },
      row('Section id', 't-id', idIn, null,
        'Identifies the section in the picker, in every RunSpec and in every export. Defaults to the directory name, so set it to a study code if the directory name is identifying. It may not shadow a built-in specimen.'),
      row('Label', 't-label', labelIn, null,
        'Free text shown in the picker, and the heading of every report exported from this section. Not used by the model. Must not identify a person.'),
      [lbl('Detection', null),
        h('label.chk', { title: 'Reads the count matrix so the report can state the MEASURED firing rate of each marker instead of only whether the gene is on the panel. Slower. Ingesting always does this.' },
          deepChk, h('span.lbl', 'measure detection rates')),
        h('div.field-err.hidden')]));

  const formEl = h('div.fgrid',
    row('Directory', 't-path', pathIn, null,
      'An absolute path ON THE MACHINE RUNNING THIS SERVER. No tissue is uploaded from this browser; the server process reads the directory off its own disk. This path is returned by GET /api/sections to every client that can reach this server, so it must not identify a person.'),
    nameFold,
  );

  /** The basename of a directory path, in either slash convention. */
  function baseName(p) {
    const s = String(p || '').trim().replace(/[/\\]+$/, '');
    if (!s) return '';
    const parts = s.split(/[/\\]/);
    return parts[parts.length - 1] || '';
  }

  /**
   * THE SERVER'S OWN ID RULE, COPIED — and the copy is declared, the way
   * js/ui/spec.js::gridCells declares its copy of the grid ceiling.
   *
   * app/server/ingest_registry.py::_make_id, verbatim:
   *
   *     base = sample_id or os.path.basename(os.path.abspath(path).rstrip(sep))
   *     slug = _ID_RE.sub("_", str(base).strip().lower()).strip("_")
   *     return slug or "ingested_sample"
   *
   * with `_ID_RE = re.compile(r"[^a-z0-9_]+")`.
   *
   * It is NOT the bare directory name: `output-XETG00001` registers as
   * `output_xetg00001`. Printing the bare basename here would show the operator
   * a string the server never uses, on the one panel whose job is to let them
   * read the identifier before it leaves the machine. If that rule ever changes
   * upstream this echo goes stale and the registered id is still whatever the
   * server says — the result panel prints it on every ingest, in full.
   */
  function slugId(base) {
    const slug = String(base || '').trim().toLowerCase()
      .replace(/[^a-z0-9_]+/g, '_').replace(/^_+|_+$/g, '');
    return slug || 'ingested_sample';
  }

  /** What the id and the label will be, as strings, before anything is sent. */
  function paintNameEcho() {
    const base = baseName(pathIn.value);
    const id = idIn.value.trim();
    const label = labelIn.value.trim();
    const deep = deepChk.checked
      ? 'detection rates measured'
      : 'filenames and panel only';
    if (!base && !id && !label) {
      nameEcho.textContent = 'id from the directory name, label from the '
        + `platform and the id · ${deep}`;
      return;
    }
    const effId = base || id ? slugId(id || base) : '(no path yet)';
    // The label default is NOT the directory name: ingest_registry.register
    // writes `f"{record.platform.title()} · {sid} (ingested)"`, and the
    // platform is not known until the scan has read the directory. Say the
    // shape and leave the platform as the unknown it is.
    const effLabel = label || `… · ${effId} (ingested)`;
    nameEcho.textContent = `id ${effId} · label ${effLabel} · ${deep}`;
  }

  /* --------------------------------------------------------- controls -- */

  // Scan is the PRIMARY control and ingest is gated behind it, in the same
  // shape js/ui/compound.js uses for compile-then-add: the operator sees what
  // the directory actually is, and what it cannot support, BEFORE anything
  // enters the section list. Editing the path or the id retires that answer —
  // a capability report that describes a different directory is worse than no
  // report, because it is the one thing on this panel meant to be trusted.
  const scanBtn = h('button.btn.primary', {
    onclick: () => go(false),
    title: 'POST /api/sections/scan — reads filenames and the gene panel, registers nothing',
  }, 'Scan');
  const ingestBtn = h('button.btn', {
    disabled: true,
    onclick: () => go(true),
    title: 'POST /api/sections/scan with register=true — full ingest, then adds the section to the picker. Scan first.',
  }, 'Ingest + add to picker');
  const statusChip = h('span.chip.st-queued', 'not scanned');

  // WHY INGEST IS OFF, IN WORDS, BESIDE THE BUTTON THAT IS OFF.
  //
  // The gate itself is unchanged and stays: nothing may be registered until a
  // scan of THIS path has come back without a refusal. What was missing is the
  // sentence. A greyed button with the reason in a `title` is a reason no
  // keyboard user reaches and no hurried user hovers, and "Ingest + add to
  // picker" looked broken rather than gated.
  const ingestWhy = h('span.dim', {
    role: 'status',
    style: {
      fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)',
      flex: '1 1 100%', minWidth: '0',
    },
  }, '');

  const controlsEl = h('div.frow', {
    style: { padding: 'var(--sp-4) var(--sp-5)', flexWrap: 'wrap' },
  }, scanBtn, ingestBtn, h('span.grow'), statusChip, ingestWhy);

  /* ---------------------------------------------------------- output -- */

  const outEl = h('div.stack.tight');
  const storedEl = h('div.stack.tight');

  /* ------------------------------------------------------------ root -- */

  // The standing block. It is the first thing in the panel, it is outside every
  // <details> in it, and it is rendered with the same component and the same
  // severity vocabulary the result panel uses for a hard stop — a limit that is
  // quieter here than it is there would be a limit the operator learns to skip.
  const standingEl = h('div.caveats',
    caveatNode(NOT_VALIDATED, 1, 'standing'),
    // Verbatim from provenance.js rather than re-worded: the project-level limit
    // is the same claim wherever it is shown, and an ingested section does not
    // soften it.
    caveatNode(STANDING_LIMITS[STANDING_LIMITS.length - 1], 2, 'standing'),
  );

  const root = h('details.sect', { open: false, id: 'sect-tissue' },
    h('summary', 'Your tissue', h('span.hint', 'directory → section')),
    standingEl,
    h('div.sect-body',
      h('p.dim', {
        style: {
          margin: '0 0 var(--sp-4)', fontSize: 'var(--fs-xs)',
          lineHeight: 'var(--lh-body)',
        },
      },
      'Point Vivum at a 10x Xenium bundle, a Space Ranger (Visium) output, or a '
      + 'bare MEX triple beside a coordinates table, and it will say what the '
      + 'sample is and what it cannot support before anything is run. '),
      h('p.dim', {
        style: {
          // --fs-disclosure, not --fs-micro. This paragraph carries the
          // security-relevant fact about WHOSE DISK IS READ, and rule 3 at the
          // top of this file makes it a stated claim rather than something to
          // be inferred. A stated claim set in the unit-suffix size is a
          // claim that reads as a footnote. See tokens.css --fs-disclosure.
          margin: '0 0 var(--sp-4)', fontSize: 'var(--fs-disclosure)',
          lineHeight: 'var(--lh-body)',
        },
      },
      // Scoped to tissue on purpose. The protocol panel DOES upload files (a
      // real multipart POST for a drug file), so the unqualified version of
      // this sentence would be a false statement about the product.
      h('b', 'Nothing is uploaded from this browser for tissue.'),
      ' You name a directory and the server process reads it from its own disk, '
      + 'so the path must exist on the machine running Vivum and that process '
      + 'must be able to read it. (Drug files, in the protocol panel, are a real '
      + 'upload.)'),
      // The de-identification constraint sits IMMEDIATELY above the fields it
      // is about. The leak surface for tissue is not the file contents — those
      // never move — it is the three strings the operator is about to type:
      // GET /api/sections hands `cache_dir` to every client that can reach the
      // server, and the id and label both default to the directory name.
      // Wording follows docs/DEPLOYMENT.md §0 rather than being re-derived.
      h('div.caveats', caveatNode(DEIDENTIFIED_ONLY, 1)),
      formEl),
    controlsEl,
    // The ANSWER sits directly under the controls that asked for it. It used to
    // sit under the bounds list below, which put the server's own reason for a
    // refusal a scroll away from the button that caused it — verified in the
    // rail at 1512x950, where the refusal block was entirely below the fold.
    h('div.sect-body', { style: { padding: '0' } }, outEl),
    // The published refusal boundaries, stated up front. Every number here is
    // read from tmesim/data/ingest.py (MIN_CELLS, MIN_EXTENT_UM,
    // MAX_EXTENT_UM) — this panel does not get to invent a limit.
    h('div.sect-body', { style: { paddingTop: '0' } },
      h('div.hr-lbl', { style: { margin: '0 0 var(--sp-2)' } }, 'What is refused, and where the bounds are'),
      h('ul', {
        style: {
          margin: '0', paddingLeft: '1.1em', display: 'flex',
          flexDirection: 'column', gap: 'var(--sp-2)',
          // "What is refused, and where the bounds are" IS the published-limits
          // surface of this panel — the list a first-time user reads to find
          // out whether their sample can be loaded at all. It is the last
          // thing that should be set smaller than the app's body prose.
          fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)',
          color: 'var(--fg-2)',
        },
      },
      h('li', 'Xenium: a per-cell table (cells.csv.gz) carrying x_centroid and '
        + 'y_centroid, plus a cell_feature_matrix MEX triple, extracted or still '
        + 'in its .tar.gz. gene_panel.json is optional and only carries provenance.'),
      h('li', 'Visium: a Space Ranger output — spatial/ with tissue_positions '
        + 'and scalefactors_json.json, beside filtered_feature_bc_matrix/. The '
        + 'two must be siblings. ',
      // Said BEFORE the path is typed, not discovered after the ingest. A
      // Visium directory registers and runs; what it cannot do is carry a
      // per-cell number, and this is the whole class of refusal a first-time
      // operator is most likely to walk into. app/server/api.py keys it on the
      // platform alone (_SPOT_SOURCES), so it lands on every Visium section
      // regardless of panel content, and its full text appears against the
      // section in the list below once it is registered.
      h('b', 'A Visium ingest is SPOT-level: '),
      'a row is a ~55 um spot pooling several cells, so it has no measured '
        + 'cells at all — a run PLACES them by deconvolution — and per-cell '
        + 'spatial statistics, the perivascular reference and the delivery '
        + 'spread are all REFUSED on it. It still ingests, and the counts it '
        + 'reports are spots.'),
      h('li', 'MEX: matrix.mtx + features.tsv + barcodes.tsv (gzipped or not) '
        + 'beside any CSV/TSV with a barcode column and an x/y pair. Ingestable '
        + 'and inspectable, but NOT runnable: the reconstruction pipeline has no '
        + 'route for it.'),
      h('li', 'Coordinates must be micrometres. Under 100 um across, or over '
        + '100,000 um (10 cm), is refused as a unit error rather than run.'),
      h('li', 'Fewer than 50 cells is refused: below that every band fraction is '
        + 'dominated by single cells.'),
      h('li', 'An ingested section can never take the id of a built-in one, so '
        + 'nothing you load can change what a published number was measured from.'))),
    h('div.sect-body', { style: { paddingTop: '0' } }, storedEl),
  );

  /* -------------------------------------------------------- behaviour -- */

  // The rail is a scroll container taller than the viewport, so an answer
  // rendered off-screen reads as a button that did nothing. `block: 'nearest'`
  // moves the rail only as far as it has to and leaves the page alone.
  function reveal() {
    try { outEl.scrollIntoView({ block: 'nearest' }); } catch { /* older engines */ }
  }

  function clearFieldErrors() {
    for (const [, { control, errEl }] of fieldEls) {
      control.removeAttribute('aria-invalid');
      control.classList.remove('invalid');
      errEl.classList.add('hidden');
      errEl.textContent = '';
    }
  }

  /**
   * The rendered report describes the directory it was read from. The moment
   * the path or the id changes it stops describing what is on screen, so the
   * report is visibly retired and "Ingest" is withdrawn until a fresh scan has
   * run. Same rule, same reason, as the compound compiler's markStale().
   */
  function markStale() {
    // `scanned` is cleared here too, not only the button. The flag is what
    // setBusy() re-derives the gate from at the end of the NEXT call, so
    // leaving it true meant a path edited after a clean scan could be ingested
    // without a fresh one — the exact walk-around this function exists to
    // close, held shut only by the button's own disabled state.
    scanned = false;
    closeGate(outEl.firstChild ? WHY_STALE : WHY_UNSCANNED);
    if (!outEl.firstChild) return;
    if (statusChip.textContent !== 'stale — scan again') {
      statusChip.className = 'chip st-queued';
      statusChip.textContent = 'stale — scan again';
      outEl.style.opacity = '0.45';
    }
  }

  function fieldError(id, message) {
    const t = fieldEls.get(id);
    if (!t) return;
    t.control.setAttribute('aria-invalid', 'true');
    t.control.classList.add('invalid');
    t.errEl.textContent = message;
    t.errEl.classList.remove('hidden');
    // An error on a folded field is an error nobody can see or fix. The fold
    // is a convenience for fields that are behaving; the moment one of them is
    // the problem it stops being one.
    if (FOLDED_FIELDS.has(id)) nameFold.open = true;
  }

  /**
   * A refusal, rendered as one. `message` is the backend's own sentence and is
   * shown in full; `remedy` is this panel's one line of what to do next.
   *
   * The two are visually separated on purpose: the operator must be able to see
   * which words are the server's finding and which are the UI's advice, because
   * only the first is evidence.
   */
  function showRefusal(message, { fieldId = 't-path', heading = 'Not ingested' } = {}) {
    clearFieldErrors();
    const remedy = remedyFor(message);
    fieldError(fieldId, remedy || 'Refused — see the reason below.');
    // Two of the nine remedies tell the operator to change a control that is
    // inside the fold ("Give this sample a different section id", 'Tick
    // "measure detection rates"'). An instruction pointing at a control the
    // reader cannot see is not a remedy, so the fold opens on those.
    if (/section id|detection rates/i.test(remedy)) nameFold.open = true;
    statusChip.className = 'chip st-failed';
    statusChip.textContent = 'refused';
    // A refused directory is never ingestable, so the gate closes again.
    scanned = false;
    closeGate(WHY_REFUSED);

    clear(outEl);
    outEl.style.opacity = '';
    outEl.appendChild(h('div', {
      style: {
        padding: 'var(--sp-4) var(--sp-5)',
        background: 'var(--crit-bg)',
        borderTop: '1px solid var(--line-2)',
        borderBottom: '1px solid var(--line-2)',
        fontSize: 'var(--fs-xs)', lineHeight: 'var(--lh-body)',
        color: 'var(--fg-0)',
      },
    },
    h('b', heading + '. '),
    remedy || 'The server refused this directory. Its reason is below.'));
    outEl.appendChild(h('div.hr-lbl',
      { style: { padding: 'var(--sp-4) var(--sp-5) var(--sp-1)' } },
      'What the server found'));
    outEl.appendChild(h('pre.mono', {
      style: {
        margin: '0 var(--sp-5) var(--sp-5)', padding: 'var(--sp-3)',
        fontSize: 'var(--fs-disclosure)', lineHeight: '1.45', color: 'var(--fg-1)',
        background: 'var(--bg-void)', border: '1px solid var(--line-1)',
        whiteSpace: 'pre-wrap', overflowWrap: 'anywhere', overflowX: 'auto',
      },
    }, String(message || '').trim()));
    outEl.appendChild(h('div.dim', {
      style: {
        padding: '0 var(--sp-5) var(--sp-5)',
        fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)',
      },
    },
    'This is a refusal, not a crash. Vivum stops on a sample it cannot honestly '
    + 'load rather than reconstructing something that would look like tissue. '
    + 'Nothing was registered and nothing was changed.'));
    reveal();
  }

  /* ----------------------------------------------------------- render -- */

  // `confidence` here is detect()'s confidence that it RECOGNISED THE FILE
  // LAYOUT — 1.0 means every component of that layout was present and named as
  // expected. It is not a confidence in any result, and a bare "confidence
  // 0.85" beside a tissue section reads as one. This project has shipped a
  // "Model confidence: High" readout before, on a product that makes no
  // clinical claim; it is not shipping another one. The number therefore never
  // appears without the words "format" and "files" next to it.
  const FORMAT_CONF_TITLE = (
    'How completely the DIRECTORY LAYOUT was recognised: 1.0 means every file '
    + 'that layout expects was present and named as expected, and it drops for '
    + 'each component that had to be inferred or is missing. It says nothing '
    + 'about the tissue, the data quality, or any number produced from it.'
  );

  function detectionRow(det) {
    const conf = det.confidence ?? det.detect_confidence;
    const hasConf = Number.isFinite(conf);
    const cells = [];
    if (det.platform) cells.push(['format detected', det.platform]);
    if (det.layout) cells.push(['layout', det.layout]);
    return h('div', {
      style: { borderBottom: '1px solid var(--line-1)' },
    },
    h('div.frow', {
      style: {
        gap: 'var(--sp-4)', padding: 'var(--sp-3) var(--sp-5) var(--sp-1)',
        fontSize: 'var(--fs-xs)', flexWrap: 'wrap',
      },
    },
    cells.flatMap(([k, v]) => [h('span.lbl', k), h('span.mono.strong', String(v))]),
    !hasConf ? null : h('span.frow', {
      style: { gap: 'var(--sp-4)' }, title: FORMAT_CONF_TITLE,
    }, h('span.lbl', 'file-format match'),
    // OUT OF WHAT. It was a bare "0.85" whose ceiling existed only in a hover:
    // a lone number between 0 and 1 beside a tissue section reads as a
    // confidence in the sample, and this one is a count of which files of a
    // known layout were found. The denominator is now printed next to it.
    h('span.mono.strong', num(conf, 2), h('span.unit', '/ 1.00')))),
    // The same sentence that was hover-only. It is prose, so it takes the
    // disclosure floor, and it is on screen because "confidence" beside a
    // patient sample is the one word on this panel most likely to be
    // misread — this project does not ship another "Model confidence: High".
    !hasConf ? null : h('div', {
      style: {
        padding: '0 var(--sp-5) var(--sp-3)',
        fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)',
        color: 'var(--fg-2)',
      },
    }, 'That is how completely the FILE LAYOUT was recognised — 1.00 means '
      + 'every file this layout expects was present and named as expected. It '
      + 'says nothing about the tissue, the data quality, or any number '
      + 'produced from it.'));
  }

  /* ------------------------------------------- per-section caveats -------
   * WHAT A SECTION CANNOT SUPPORT, BEFORE A RUN IS LAUNCHED.
   *
   * GET /api/sections returns `caveats` on every record — a list of
   * `{capability, status, reason}` built by app/server/api.py::_with_caveats —
   * plus `spatial_spread_valid` and `n_caveats_blocking`. Until this pass this
   * panel read that response for one thing (which ids are ingested) and threw
   * the rest away, so the list of registered sections was an id, a chip and a
   * path, and the refusals attached to them were reachable only from a `title`
   * on an <option> in the picker on the other rail.
   *
   * WHICH SECTIONS CARRY WHICH, verified rather than assumed:
   *
   *  * The SPOT refusal is not in ingest_registry.py at all. It is added at
   *    the HTTP boundary, keyed on the platform — `_SPOT_SOURCES = ("visium",)`
   *    matched against `platform` (ingested) or `source` (built-in) — so it
   *    lands on every Visium section, the four built-in ones and any Visium
   *    directory ingested here alike. Its text is api.py::SPOT_SPATIAL_REFUSAL
   *    and it says per-cell spatial statistics AND the delivery map are invalid
   *    on such a section, because the endothelial positions its source term is
   *    built from are synthesised like every other cell's.
   *  * ingest_registry.py adds exactly one caveat of its own: on a `mex`
   *    bundle it appends an "INSPECTABLE, not runnable" detect_issue and sets
   *    `runnable: False`.
   *  * everything else comes from the stored CAPABILITY report, which only an
   *    ingested section carries, and from `detect_issues`.
   *
   * The reasons are rendered IN FULL and never behind a toggle, with the same
   * component the run's own caveats use, because this panel's second rule is
   * that the refusal is the product. */
  function sectionCaveats(rec) {
    const list = Array.isArray(rec.caveats) ? rec.caveats : null;
    const box = h('div');
    if (!list) {
      // Not silence. An older server, or one that could not describe this
      // section, returns the record without the field, and "no caveats shown"
      // must never be read as "no caveats".
      box.appendChild(h('div.caveats', caveatNode(
        'CAVEATS NOT STATED FOR THIS SECTION. This server did not return a '
        + 'caveats list on GET /api/sections, so what this specimen cannot '
        + 'support is unknown from this screen — that is not the same as '
        + 'nothing being refused. Scan the directory again, or read the '
        + 'capability report above.', 1)));
      return box;
    }
    const blocking = list.filter((c) => c.status === 'NO');
    const qualified = list.filter((c) => c.status !== 'NO');
    box.appendChild(h('div.cv-index', {
      title: 'Every caveat below is shown in full — none is collapsed, clamped '
        + 'or dismissible. "Refused" counts the readouts this section cannot '
        + 'support at all; "qualified" counts the ones that come with a '
        + 'condition on how they are read.',
    },
    h('span' + (blocking.length ? '.hard' : ''),
      h('b', String(blocking.length)), ' refused'),
    h('span.sep', '·'),
    h('span', h('b', String(qualified.length)), ' qualified'),
    h('span.grow'),
    h('span.dim', 'none hidden')));
    if (!list.length) {
      box.appendChild(h('div.caveat.info', h('span',
        'This server refused nothing and qualified nothing on this section. '
        + 'Its capability report, above, is still the only check it has '
        + 'passed.')));
      return box;
    }
    const cav = h('div.caveats');
    // Refusals first: a reader who stops after one line should have stopped
    // after the one that says a number cannot be produced at all.
    [...blocking, ...qualified].forEach((c, i) => {
      const text = `${c.capability || 'capability'} — ${c.status || 'CAVEAT'}: `
        + `${c.reason || 'no reason given by the server'}`;
      cav.appendChild(caveatNode(text, i + 1, c.status === 'NO' ? 'crit' : ''));
    });
    box.appendChild(cav);
    return box;
  }

  function render(res, registered) {
    clear(outEl);
    outEl.style.opacity = '';

    const sec = res.section || null;
    if (registered && res.ok) {
      // The banner on a SUCCESSFUL ingest is the not-validated one, in the
      // colour the rest of the Studio uses for a limit — not a green "added".
      // What just happened is that an unscored specimen entered the picker.
      outEl.appendChild(h('div', {
        style: {
          padding: 'var(--sp-4) var(--sp-5)',
          background: 'var(--warn-bg)',
          borderTop: '1px solid var(--line-2)',
          borderBottom: '1px solid var(--line-2)',
          fontSize: 'var(--fs-xs)', lineHeight: 'var(--lh-body)',
          color: 'var(--fg-0)',
        },
      },
      h('b', `Registered as ${res.section_id} — and NOT VALIDATED. `),
      'It is in the section picker now. It is not one of the built-in specimens '
      + 'this repo ships evidence for, and nothing has been scored against it. '
      + 'Read the capability report below before you read any number it produces.'));
      if (sec && sec.runnable === false) {
        outEl.appendChild(h('div.caveats', caveatNode(NOT_RUNNABLE, 1)));
      }
    }

    const det = res.detection || (sec ? {
      platform: sec.platform, layout: sec.layout,
      confidence: sec.detect_confidence, issues: sec.detect_issues,
    } : null);
    if (det) {
      outEl.appendChild(detectionRow(det));
      if (sec) {
        outEl.appendChild(h('div.frow', {
          style: {
            gap: 'var(--sp-4)', padding: 'var(--sp-3) var(--sp-5)',
            borderBottom: '1px solid var(--line-1)', fontSize: 'var(--fs-xs)',
            flexWrap: 'wrap',
          },
        },
        // "cells" was the label on every platform, including Visium — where
        // `n_cells_ingested` is a count of ~55 um SPOTS, each pooling several
        // cells (tmesim/data/ingest.py; cache.py copies it into `n_spots` for
        // exactly this platform). One word, and it was the wrong one on a
        // quarter of the built-in catalogue.
        h('span.lbl', isSpot(sec) ? 'spots' : 'cells'),
        Number.isFinite(sec.n_cells_ingested)
          ? h('span.mono.strong', int(sec.n_cells_ingested),
            h('span.unit', isSpot(sec) ? 'spots' : 'cells'))
          : h('span.dim', 'not stated'),
        h('span.lbl', 'genes'),
        Number.isFinite(sec.n_genes_ingested)
          ? h('span.mono.strong', int(sec.n_genes_ingested), h('span.unit', 'genes'))
          : h('span.dim', 'not stated'),
        h('span.lbl', 'extent'),
        // `(x || []).map(...).join(' × ') + ' mm'` printed the bare string
        // " mm" when the record carried no extent: a unit with nothing in
        // front of it, which reads as a measurement of zero.
        Array.isArray(sec.extent_um_ingested) && sec.extent_um_ingested.length >= 2
          ? h('span.mono.strong',
            sec.extent_um_ingested.map((v) => num(v / 1000, 2)).join(' × '),
            h('span.unit', 'mm'))
          : h('span.dim', 'not stated')));
        if (isSpot(sec)) {
          outEl.appendChild(h('div', {
            style: {
              padding: 'var(--sp-1) var(--sp-5) var(--sp-3)',
              fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)',
              color: 'var(--fg-2)', borderBottom: '1px solid var(--line-1)',
            },
          }, 'A Visium row is a ~55 um spot pooling several cells, so the count '
            + 'above is spots. This section has no measured cells at all: a run '
            + 'PLACES them from the per-spot deconvolved fractions, which is why '
            + 'per-cell statistics are refused on it — see the refusals under '
            + 'its entry in the list at the bottom of this panel.'));
        }
      }
      const issues = det.issues || det.detect_issues || [];
      if (issues.length) {
        outEl.appendChild(h('div.hr-lbl',
          { style: { padding: 'var(--sp-4) var(--sp-5) var(--sp-1)' } },
          `Detection caveats · ${issues.length}`));
        const box = h('div.caveats');
        issues.forEach((c, i) => box.appendChild(caveatNode(String(c), i + 1)));
        outEl.appendChild(box);
      }
    }

    if (res.capability_text) {
      outEl.appendChild(h('div.hr-lbl',
        { style: { padding: 'var(--sp-4) var(--sp-5) var(--sp-1)' } },
        'Capability report'));
      // Verbatim, monospaced, and IN FULL: this is the exact text
      // `vivum ingest` prints, and it is the only place the operator is told
      // what their own sample cannot support. It gets its full height — no
      // clamp, no inner scroller, no "show more" — because a limit behind a
      // scrollbar is a limit that does not get read. Only the horizontal axis
      // scrolls, so the fixed-width report keeps its columns without widening
      // the rail.
      outEl.appendChild(h('pre.mono', {
        style: {
          margin: '0 var(--sp-5) var(--sp-5)', padding: 'var(--sp-3)',
          overflowX: 'auto',
          fontSize: 'var(--fs-disclosure)', lineHeight: '1.45', color: 'var(--fg-1)',
          background: 'var(--bg-void)', border: '1px solid var(--line-1)',
          whiteSpace: 'pre',
        },
      }, res.capability_text));
    }

    if (!outEl.firstChild) {
      showRefusal('The server answered with neither a capability report nor a '
        + 'reason. That is a bug — do not work around it.');
      return;
    }
    reveal();
  }

  /* ------------------------------------------------------------ calls -- */

  let busy = false;
  let tick = 0;
  // True once a scan of the CURRENT path+id has come back without a refusal.
  // Cleared by markStale() on any edit, so the gate cannot be walked around by
  // scanning one directory and ingesting another.
  let scanned = false;
  // The sentence beside the disabled Ingest button. It is state, not decoration:
  // every path that closes the gate names its own reason, and setBusy() restores
  // that reason when the request it was covering finishes.
  let whyOff_ = WHY_UNSCANNED;

  /** Close the ingest gate, and say in words why it is closed. */
  function closeGate(reason) {
    whyOff_ = reason;
    ingestBtn.disabled = true;
    ingestWhy.textContent = busy ? '' : reason;
  }

  /** Open it. Only a scan of the current path that came back clean may. */
  function openGate() {
    whyOff_ = '';
    scanned = true;
    ingestBtn.disabled = busy;
    ingestWhy.textContent = '';
  }

  function setBusy(on, verb) {
    busy = on;
    scanBtn.disabled = on;
    ingestBtn.disabled = on || !scanned;
    ingestWhy.textContent = on ? '' : (scanned ? '' : whyOff_);
    clearInterval(tick);
    if (!on) return;
    const t0 = Date.now();
    statusChip.className = 'chip st-running';
    const paint = () => {
      statusChip.textContent = `${verb} · ${Math.round((Date.now() - t0) / 1000)} s`;
    };
    paint();
    tick = setInterval(paint, 1000);
  }

  async function go(registered) {
    if (busy) return;
    const path = pathIn.value.trim();
    clearFieldErrors();
    if (!path) {
      fieldError('t-path', 'Enter the sample directory as an absolute path on '
        + 'the machine running this server.');
      pathIn.focus();
      return;
    }
    if (api.api && api.api.stub) {
      showRefusal('The Studio is running on stub data because the backend is '
        + 'not answering, so no directory was read and no capability report was '
        + 'produced. A synthesised one would look exactly like a measured one. '
        + 'Start the backend and try again.',
      { heading: 'Refused — no backend' });
      return;
    }

    const body = { path, measure_detection_rates: !!deepChk.checked, register: !!registered };
    if (registered) body.measure_detection_rates = true;  // a register is always a full ingest
    if (idIn.value.trim()) body.section_id = idIn.value.trim();
    if (labelIn.value.trim()) body.label = labelIn.value.trim();

    setBusy(true, registered ? 'ingesting' : 'scanning');
    clear(outEl);
    outEl.style.opacity = '';
    outEl.appendChild(h('div.dim', {
      style: {
        padding: 'var(--sp-4) var(--sp-5)', fontSize: 'var(--fs-xs)',
        lineHeight: 'var(--lh-body)',
      },
    },
    registered || body.measure_detection_rates
      ? 'Reading the count matrix. This is one synchronous request and a large '
        + 'section takes minutes — leave this panel open.'
      : 'Reading filenames and the gene panel only.'));
    reveal();

    try {
      const res = await api.postSectionScan(body, { timeoutMs: INGEST_TIMEOUT_MS });
      if (res && res.error) {
        showRefusal(res.error);
        return;
      }
      render(res, registered);
      statusChip.className = 'chip st-done';
      statusChip.textContent = registered ? 'registered' : 'scanned';
      // A clean scan arms the ingest; a completed ingest disarms it again,
      // because a second press would only collide with the id it just took.
      if (registered) { scanned = false; closeGate(WHY_REGISTERED); } else openGate();
      if (registered && res.ok) {
        toast('warn', `Section ${res.section_id} ingested — not validated`,
          'It is in the picker. Nothing has been scored against it: read its '
          + 'capability report before you read any number it produces.', 12000);
        await refreshStored();
        if (onRegistered) onRegistered(res.section_id);
      }
    } catch (err) {
      const status = err && err.status;
      if (status === 401) {
        showRefusal('This server requires an access key and the browser did not '
          + 'send a valid one. Reload the Studio and enter the key.',
        { heading: 'Refused — not authenticated' });
      } else if (status === 400) {
        showRefusal(String((err && err.message) || 'the request was rejected'));
      } else if (err && err.name === 'AbortError') {
        showRefusal('The browser stopped waiting after 15 minutes. THE SERVER MAY '
          + 'STILL BE INGESTING — check the section picker before starting over, '
          + 'and do not assume nothing happened.',
        { heading: 'Timed out waiting' });
      } else {
        // Anything else is a transport failure rather than a verdict on the
        // sample, and must not be dressed up as one.
        statusChip.className = 'chip st-failed';
        statusChip.textContent = 'failed';
        notifyError('POST /api/sections/scan failed', err);
      }
    } finally {
      setBusy(false);
      if (statusChip.className === 'chip st-running') {
        statusChip.className = 'chip st-queued';
        statusChip.textContent = 'idle';
      }
    }
  }

  /* --------------------------------------------------- ingested so far -- */

  /** One sentence, at the disclosure floor, in the panel's own voice. */
  function storedNote(text) {
    return h('div', {
      style: {
        fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)',
        color: 'var(--fg-2)', padding: 'var(--sp-2) 0',
      },
    }, text);
  }

  async function refreshStored() {
    let list = [];
    try {
      list = await api.getSections();
    } catch (err) {
      // Was `clear(storedEl); return;` — a silent blank. The list of sections
      // this operator has registered vanishing without a word is exactly the
      // shape of "there are none", and the two must not look the same.
      clear(storedEl);
      storedEl.appendChild(h('div.hr-lbl', { style: { margin: '0 0 var(--sp-2)' } },
        'Ingested sections · unknown'));
      storedEl.appendChild(storedNote(
        'GET /api/sections did not answer, so this list could not be read. '
        + 'Nothing has been removed — what is registered on the server is '
        + 'unchanged and simply unknown from here. '
        + String((err && err.message) || '')));
      return;
    }
    // Array.isArray, not `list || []`: a server that answers 200 with an object
    // instead of a list threw a TypeError out of this function and left the
    // whole list — including every refusal on it — blank, with the failure
    // visible only in the console.
    if (!Array.isArray(list)) {
      clear(storedEl);
      storedEl.appendChild(h('div.hr-lbl', { style: { margin: '0 0 var(--sp-2)' } },
        'Ingested sections · unknown'));
      storedEl.appendChild(storedNote(
        'GET /api/sections answered with something that is not a list of '
        + 'sections, so this list could not be read. Nothing has been removed; '
        + 'what is registered on the server is unchanged and unknown from here.'));
      return;
    }
    const mine = list.filter((s) => s && s.ingested);
    clear(storedEl);
    if (!mine.length) {
      // An empty list is a state with two causes, and the operator's next
      // action differs between them.
      storedEl.appendChild(h('div.hr-lbl', { style: { margin: '0 0 var(--sp-2)' } },
        'Ingested sections · 0'));
      storedEl.appendChild(storedNote(api.api && api.api.stub
        ? 'The Studio is running on stub data because the backend is not '
          + 'answering, so this list is not a statement about what is '
          + 'registered on any server. Start the backend to read the real one.'
        : 'No section has been ingested on this server. The built-in specimens '
          + 'are still in the picker; nothing here removes or replaces them. '
          + 'Sections registered by a SEPARATE process are not picked up until '
          + 'this server restarts (app/server/ingest_registry.py::load_store).'));
      return;
    }
    const refusedTotal = mine.reduce(
      (a, s) => a + (Array.isArray(s.caveats)
        ? s.caveats.filter((c) => c.status === 'NO').length : 0), 0);
    storedEl.appendChild(h('div.hr-lbl', { style: { margin: '0 0 var(--sp-2)' } },
      `Ingested sections · ${mine.length} · none validated`
      + (refusedTotal ? ` · ${refusedTotal} refused readout`
        + `${refusedTotal === 1 ? '' : 's'}` : '')));
    for (const rec of mine) {
      storedEl.appendChild(h('div.frow', {
        style: {
          gap: 'var(--sp-3)', alignItems: 'baseline',
          padding: 'var(--sp-2) 0',
          fontSize: 'var(--fs-xs)',
        },
      },
      h('span.mono', { style: { color: 'var(--fg-0)' } }, rec.id),
      h('span.chip.assumed', {
        title: NOT_VALIDATED,
      }, rec.runnable === false ? 'inspect only' : 'not validated'),
      h('span.grow'),
      h('span.dim.mono', {
        title: rec.cache_dir || '',
        style: { maxWidth: '16ch', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap', direction: 'rtl' },
      }, rec.cache_dir || ''),
      // Two-step: removing a section silently changes what the picker offers,
      // and a run already made against it keeps referring to an id that is no
      // longer resolvable. The bytes on disk are untouched either way, which is
      // why this needs no undo — re-ingesting the same directory restores it.
      confirmable(h('button.btn.ghost.sm', {
        title: `Remove this ingested section from the picker. DELETE /api/sections/${rec.id}. The directory on disk is not touched.`,
        'aria-label': `Remove ingested section ${rec.id}`,
      }, 'remove'), {
        confirmLabel: `remove ${rec.id}?`,
        armedTitle: `This takes ${rec.id} out of the section picker. Nothing on `
          + 'disk is deleted and re-ingesting the same directory brings it back. '
          + 'Press again to confirm; Escape cancels.',
        onConfirm: async () => {
          try {
            await api.deleteSection(rec.id);
            await refreshStored();
            if (onRegistered) onRegistered(null);
            toast('info', `${rec.id} removed from the picker`,
              'The directory on disk was not touched. Runs already made against '
              + 'it are unchanged, but they now name a section the picker no '
              + 'longer offers.', 9000);
          } catch (err) { notifyError(`DELETE /api/sections/${rec.id} failed`, err); }
        },
      }),
      ));
      // What this specimen IS — scale, in the right noun.
      storedEl.appendChild(sectionScaleRow(rec));
      // Whether it can be reached at all. `available: false` means the bytes
      // the record points at are gone or no longer ingestable; the row used to
      // look identical to a healthy one.
      if (rec.available === false) {
        storedEl.appendChild(h('div.caveats', caveatNode(
          'SECTION UNAVAILABLE: ' + String(rec.reason
            || 'the server did not say why. The record is still registered; '
              + 'the data it points at could not be read.'), 1, 'crit')));
      }
      // What it CANNOT support — in full, before anything is run against it.
      storedEl.appendChild(sectionCaveats(rec));
      storedEl.appendChild(h('div', {
        style: { borderBottom: '1px solid var(--line-2)', margin: '0 0 var(--sp-5)' },
      }));
    }
  }

  // The two read-outs that describe a state rather than a result, painted once
  // at construction so the panel is never blank about its own gate: the id and
  // label that would be used, and why Ingest is not available yet.
  paintNameEcho();
  ingestWhy.textContent = WHY_UNSCANNED;

  return {
    el: root,
    refresh: refreshStored,
    open() { root.open = true; pathIn.focus(); },
  };
}
