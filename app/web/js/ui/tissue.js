// tissue.js — BRING YOUR OWN TISSUE: point Codon at a spatial-transcriptomics
// export on the server's disk, get a verdict on what that specimen can and
// cannot support, and — only if it survives every check — a section the picker
// can select.
//
// Why this panel exists at all
// ----------------------------
// The capability existed in `POST /api/sections/scan` and in `codon ingest`
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
//     sentence `codon ingest` prints, so this surface cannot drift from the
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
    + 'Codon server, not on the computer your browser is running on, and it has '
    + 'to be the directory itself rather than a file inside it.'],
  [/is a CONTAINER, not a sample/i,
    'This directory holds your samples rather than being one. Ingest one of the '
    + 'child directories named above, one at a time.'],
  [/no recognised layout/i,
    'Codon could not tell what this directory is. Check you are pointing at the '
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
    'That id belongs to a built-in specimen, and Codon will not let an ingested '
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
  'DE-IDENTIFIED TISSUE ONLY. Codon has one shared secret and no accounts, no '
  + 'audit log, no encryption at rest and no per-record access control, so it '
  + 'must not hold data that identifies a person. The tissue itself never '
  + 'leaves this machine — but three of the strings you are about to type do: '
  + 'GET /api/sections returns the DIRECTORY PATH to every client that can '
  + 'reach this server and writes it to app/runs/_ingested_sections.json, and '
  + 'the SECTION ID and LABEL both default to the directory name, appear in the '
  + 'picker, and put the label in the heading of every exported report. Copy or '
  + 'rename to a study code first. The full checklist, including the '
  + 'identifiers that ride inside a 10x export without Codon ever reading them, '
  + 'is docs/DEPLOYMENT.md §0.'
);

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
    oninput: () => { clearFieldErrors(); markStale(); },
    onkeydown: (e) => { if (e.key === 'Escape') { e.stopPropagation(); pathIn.blur(); } },
  });
  const idIn = h('input.inp.mono', {
    id: 't-id', type: 'text', autocomplete: 'off', spellcheck: 'false',
    placeholder: 'defaults to the directory name',
    oninput: () => { clearFieldErrors(); markStale(); },
  });
  const labelIn = h('input.inp', {
    id: 't-label', type: 'text', autocomplete: 'off',
    placeholder: 'optional — shown in the section picker',
  });
  const deepChk = h('input', { type: 'checkbox', id: 't-deep' });

  const lbl = (text, forId, title) => h('label.lbl', { for: forId, title: title || '' }, text);

  function row(label, forId, control, unitText, title) {
    const holder = unitText ? h('div.with-unit', control, h('span.unit', unitText)) : control;
    const errId = forId + '-err';
    const errEl = h('div.field-err.hidden', { id: errId, role: 'alert' });
    control.setAttribute('aria-describedby', errId);
    fieldEls.set(forId, { control, errEl });
    return [lbl(label, forId, title), holder, errEl];
  }

  const formEl = h('div.fgrid',
    row('Directory', 't-path', pathIn, null,
      'An absolute path ON THE MACHINE RUNNING THIS SERVER. No tissue is uploaded from this browser; the server process reads the directory off its own disk. This path is returned by GET /api/sections to every client that can reach this server, so it must not identify a person.'),
    row('Section id', 't-id', idIn, null,
      'Identifies the section in the picker, in every RunSpec and in every export. Defaults to the directory name, so set it to a study code if the directory name is identifying. It may not shadow a built-in specimen.'),
    row('Label', 't-label', labelIn, null,
      'Free text shown in the picker, and the heading of every report exported from this section. Not used by the model. Must not identify a person.'),
    [lbl('Detection', null),
      h('label.chk', { title: 'Reads the count matrix so the report can state the MEASURED firing rate of each marker instead of only whether the gene is on the panel. Slower. Ingesting always does this.' },
        deepChk, h('span.lbl', 'measure detection rates')),
      h('div.field-err.hidden')],
  );

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

  const controlsEl = h('div.frow', { style: { padding: 'var(--sp-4) var(--sp-5)' } },
    scanBtn, ingestBtn, h('span.grow'), statusChip);

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
      'Point Codon at a 10x Xenium bundle, a Space Ranger (Visium) output, or a '
      + 'bare MEX triple beside a coordinates table, and it will say what the '
      + 'sample is and what it cannot support before anything is run. '),
      h('p.dim', {
        style: {
          margin: '0 0 var(--sp-4)', fontSize: 'var(--fs-micro)',
          lineHeight: 'var(--lh-body)',
        },
      },
      // Scoped to tissue on purpose. The protocol panel DOES upload files (a
      // real multipart POST for a drug file), so the unqualified version of
      // this sentence would be a false statement about the product.
      h('b', 'Nothing is uploaded from this browser for tissue.'),
      ' You name a directory and the server process reads it from its own disk, '
      + 'so the path must exist on the machine running Codon and that process '
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
          fontSize: 'var(--fs-micro)', lineHeight: 'var(--lh-body)',
          color: 'var(--fg-2)',
        },
      },
      h('li', 'Xenium: a per-cell table (cells.csv.gz) carrying x_centroid and '
        + 'y_centroid, plus a cell_feature_matrix MEX triple, extracted or still '
        + 'in its .tar.gz. gene_panel.json is optional and only carries provenance.'),
      h('li', 'Visium: a Space Ranger output — spatial/ with tissue_positions '
        + 'and scalefactors_json.json, beside filtered_feature_bc_matrix/. The '
        + 'two must be siblings.'),
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
    ingestBtn.disabled = true;
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
    statusChip.className = 'chip st-failed';
    statusChip.textContent = 'refused';
    // A refused directory is never ingestable, so the gate closes again.
    scanned = false;
    ingestBtn.disabled = true;

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
        fontSize: 'var(--fs-micro)', lineHeight: '1.45', color: 'var(--fg-1)',
        background: 'var(--bg-void)', border: '1px solid var(--line-1)',
        whiteSpace: 'pre-wrap', overflowWrap: 'anywhere', overflowX: 'auto',
      },
    }, String(message || '').trim()));
    outEl.appendChild(h('div.dim', {
      style: {
        padding: '0 var(--sp-5) var(--sp-5)',
        fontSize: 'var(--fs-micro)', lineHeight: 'var(--lh-body)',
      },
    },
    'This is a refusal, not a crash. Codon stops on a sample it cannot honestly '
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
    const cells = [];
    if (det.platform) cells.push(['format detected', det.platform]);
    if (det.layout) cells.push(['layout', det.layout]);
    return h('div.frow', {
      style: {
        gap: 'var(--sp-4)', padding: 'var(--sp-3) var(--sp-5)',
        borderBottom: '1px solid var(--line-1)', fontSize: 'var(--fs-xs)',
        flexWrap: 'wrap',
      },
    },
    cells.flatMap(([k, v]) => [h('span.lbl', k), h('span.mono.strong', String(v))]),
    conf === undefined || conf === null ? null : h('span.frow', {
      style: { gap: 'var(--sp-4)' }, title: FORMAT_CONF_TITLE,
    }, h('span.lbl', 'file-format match'), h('span.mono.strong', num(conf, 2))));
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
        h('span.lbl', 'cells'), h('span.mono.strong', int(sec.n_cells_ingested)),
        h('span.lbl', 'genes'), h('span.mono.strong', int(sec.n_genes_ingested)),
        h('span.lbl', 'extent'),
        h('span.mono.strong', (sec.extent_um_ingested || [])
          .map((v) => num(v / 1000, 2)).join(' × ') + ' mm')));
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
      // `codon ingest` prints, and it is the only place the operator is told
      // what their own sample cannot support. It gets its full height — no
      // clamp, no inner scroller, no "show more" — because a limit behind a
      // scrollbar is a limit that does not get read. Only the horizontal axis
      // scrolls, so the fixed-width report keeps its columns without widening
      // the rail.
      outEl.appendChild(h('pre.mono', {
        style: {
          margin: '0 var(--sp-5) var(--sp-5)', padding: 'var(--sp-3)',
          overflowX: 'auto',
          fontSize: 'var(--fs-micro)', lineHeight: '1.45', color: 'var(--fg-1)',
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

  function setBusy(on, verb) {
    busy = on;
    scanBtn.disabled = on;
    ingestBtn.disabled = on || !scanned;
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
      scanned = !registered;
      ingestBtn.disabled = registered;
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

  async function refreshStored() {
    let list = [];
    try {
      list = await api.getSections();
    } catch {
      clear(storedEl);
      return;
    }
    const mine = (list || []).filter((s) => s.ingested);
    clear(storedEl);
    if (!mine.length) return;
    storedEl.appendChild(h('div.hr-lbl', { style: { margin: '0 0 var(--sp-2)' } },
      `Ingested sections · ${mine.length} · none validated`));
    for (const rec of mine) {
      storedEl.appendChild(h('div.frow', {
        style: {
          gap: 'var(--sp-3)', alignItems: 'baseline',
          padding: 'var(--sp-2) 0', borderBottom: '1px solid var(--line-1)',
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
    }
  }

  return {
    el: root,
    refresh: refreshStored,
    open() { root.open = true; pathIn.focus(); },
  };
}
