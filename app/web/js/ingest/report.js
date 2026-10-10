// report.js — what the selected directory IS, and the one action on it.
//
// THE DEFECT THIS FIXES, as much as the typing did
// The old form had two buttons and one outcome for everything that was not a
// perfect bundle: a bare failure. Meanwhile `tmesim/data/ingest.py::detect`
// had already worked out WHY — "cell_feature_matrix/ is INCOMPLETE: missing
// barcodes.tsv[.gz]", "spatial: no scalefactors_json.json - spot_diameter_
// fullres is the only way to convert pixels to micrometres" — and that
// sentence reached a console and no screen. So this pane leads with the
// detection verdict for the directory in view, before anyone presses anything.
//
// Every sentence here that judges a sample comes from the server verbatim.
// `capability_text` is the exact report `vivum ingest` prints, rendered as
// received, so this surface cannot come to disagree with the command line
// about what a specimen supports. The only words this file authors are the
// labels around them.

import { h, fill, clear } from '../util/dom.js';
import { inspect, register, ServerRefusal } from './net.js';

const PLATFORM_LABEL = {
  xenium: '10x Xenium',
  visium: 'Space Ranger / Visium',
  mex: 'MEX triple + coordinates',
  multi_sample: 'a container of samples',
  unknown: 'not recognised',
};

function pct(x) {
  return typeof x === 'number' ? Math.round(x * 100) + '%' : null;
}

/** The lines `detect` produced, as a list. Never flattened to one string:
 *  each one is a separate missing component and reads as one. */
function lines(cls_, title, items) {
  if (!items || !items.length) return null;
  return h('div.rp-lines',
    h('h3', title),
    h('ul', items.map((s) => h('li', { class: cls_ }, s))));
}

export function createReport({ onRegistered }) {
  const state = { path: '', diagnosis: null, busy: false };

  const headEl = h('div#rp-head.rp-head');
  const verdictEl = h('div#rp-verdict.rp-verdict');
  const actionsEl = h('div#rp-actions.rp-actions');
  const outEl = h('div#rp-out.rp-out', { 'aria-live': 'polite' });

  const idEl = h('input#rp-id.inp.sm', {
    type: 'text', spellcheck: 'false',
    placeholder: 'defaults to the directory name',
    'aria-label': 'Section id',
  });
  const labelEl = h('input#rp-label.inp.sm', {
    type: 'text', spellcheck: 'false',
    placeholder: 'shown in the section picker',
    'aria-label': 'Label',
  });
  const ratesEl = h('input#rp-rates', { type: 'checkbox' });

  const inspectBtn = h('button#rp-inspect.btn.sm', { type: 'button' },
    'Inspect — read only');
  const registerBtn = h('button#rp-register.btn.primary.sm', { type: 'button' },
    'Add as a section');

  const el = h('section#rp.rp', { 'aria-label': 'What this directory is' },
    headEl, verdictEl, actionsEl,
    h('details#rp-naming.rp-naming',
      h('summary', 'Name it something other than the directory name'),
      h('div.rp-field', h('label', { for: 'rp-id' }, 'Section id'), idEl),
      h('div.rp-field', h('label', { for: 'rp-label' }, 'Label'), labelEl),
      h('p.rp-hint',
        'The id is how runs and reports refer to this section, so it cannot '
        + 'be one a built-in section already uses — the server refuses that '
        + 'rather than letting an added section change what a published '
        + 'number was measured from.')),
    outEl);

  // --- rendering --------------------------------------------------------- //
  function renderEmpty() {
    clear(verdictEl);
    clear(actionsEl);
    clear(outEl);
    fill(headEl, h('p.rp-hint',
      'Pick a directory on the left. This panel says what it is, and what it '
      + 'cannot support, before anything is added.'));
  }

  function renderHead() {
    fill(headEl,
      h('h2.rp-title', 'This directory'),
      h('code.rp-path', state.path));
  }

  function renderVerdict() {
    const d = state.diagnosis;
    if (!d) { clear(verdictEl); return; }
    const bits = [];

    if (d.registered_as) {
      bits.push(h('p.rp-already',
        'Already a section on this server, as ',
        h('code', d.registered_as),
        '. Adding it again is refused; the server will not let one directory '
        + 'be two sections.'));
    }

    const chips = [
      h('span.chip', PLATFORM_LABEL[d.platform] || d.platform || 'unknown'),
      d.layout ? h('span.chip', d.layout) : null,
      pct(d.confidence) ? h('span.chip', 'detected ' + pct(d.confidence)) : null,
      h('span', {
        class: 'chip ' + (d.ingestable ? 'ok-chip' : 'no-chip'),
      }, d.ingestable ? 'can be read' : 'cannot be read as a section'),
    ];
    bits.push(h('div.rp-chips', chips));

    if (d.error) {
      bits.push(h('p.rp-refusal', d.error));
    }
    bits.push(lines('issue', d.ingestable
      ? 'Reads, with these caveats' : 'Why it cannot be read', d.issues));
    bits.push(lines('note', 'Worth knowing', d.notes));

    const names = Object.keys(d.samples || {});
    if (names.length) {
      bits.push(h('p.rp-contained',
        names.length === 1
          ? 'One directory inside this one is itself a recognised bundle. '
          : names.length + ' directories inside this one are themselves '
            + 'recognised bundles. ',
        'Each is a separate section: open the one you want on the left.'));
    }
    fill(verdictEl, bits);
  }

  function renderActions() {
    const d = state.diagnosis;
    if (!state.path) { clear(actionsEl); return; }
    registerBtn.disabled = !(d && d.ingestable) || !!(d && d.registered_as)
      || state.busy;
    inspectBtn.disabled = state.busy;
    fill(actionsEl,
      h('div.btn-row', inspectBtn, registerBtn),
      h('label.rp-check', ratesEl,
        'measure detection rates — reads the whole sample, so the report can '
        + 'say "MKI67 fires in 1.34% of cells" instead of "MKI67 is on the '
        + 'panel". A minute on a large section.'),
      d && !d.ingestable && !d.registered_as
        ? h('p.rp-hint', 'Adding is held back because the detection above '
          + 'says this directory cannot be read as a section. Fix what it '
          + 'names, or open a directory that can.')
        : null);
  }

  function render() {
    if (!state.path) { renderEmpty(); return; }
    renderHead();
    renderVerdict();
    renderActions();
  }

  // --- the two actions --------------------------------------------------- //
  function busy(on, what) {
    state.busy = on;
    renderActions();
    if (on) {
      fill(outEl, h('p.rp-working', what));
    }
  }

  function renderCapability(body) {
    const bits = [];
    if (body.error) {
      // A DELIBERATE refusal (IngestError, the capability report, a detection
      // issue) arrives here verbatim and is shown verbatim. An unexpected
      // internal failure has already been replaced server-side by a stable
      // sentence carrying a reference id — which is the part worth quoting,
      // so it must survive to the screen intact.
      bits.push(h('p.rp-refusal', body.error));
    }
    if (body.capability_text) {
      bits.push(h('h3.rp-sub', 'Capability report'));
      bits.push(h('pre.rp-pre', body.capability_text));
    }
    if (!bits.length) {
      bits.push(h('p.rp-hint',
        'The server returned no capability report for this directory.'));
    }
    fill(outEl, bits);
  }

  async function doInspect() {
    busy(true, ratesEl.checked
      ? 'Reading the whole sample. On a large section this takes about a minute.'
      : 'Reading the file names and the panel. This is milliseconds.');
    try {
      renderCapability(await inspect(state.path,
        { measureRates: ratesEl.checked }));
    } catch (err) {
      if (!(err instanceof ServerRefusal)) throw err;
      fill(outEl, h('p.rp-refusal', err.message));
    } finally {
      busy(false);
    }
  }

  async function doRegister() {
    busy(true, 'Reading and validating the sample. Registration is refused, '
      + 'not degraded, if it does not load — so this runs the whole ingest '
      + 'first.');
    try {
      const body = await register(state.path, {
        sectionId: idEl.value.trim(), label: labelEl.value.trim(),
      });
      // 200 with ok:false is the ordinary refusal shape on this endpoint: the
      // id collided, or the sample failed validation and nothing was
      // registered. Reading the status alone would report success.
      if (!body.ok) { renderCapability(body); return; }
      const sid = body.section_id;
      const runnable = !(body.section && body.section.runnable === false);
      fill(outEl,
        h('p.rp-ok',
          'Added as ', h('code', sid), '. It is in ',
          h('code', 'GET /api/sections'), ' now, and in the Studio\'s section '
          + 'picker after a reload.'),
        runnable ? null : h('p.rp-refusal',
          'Inspectable only: the reconstruction pipeline has no route for '
          + 'this platform, so the section can be examined but not run.'),
        h('div.btn-row',
          h('a.btn.sm', { href: '/studio' }, 'Open the Studio'),
          h('a.btn.ghost.sm', { href: '/api/sections' }, 'See the section list')),
        body.capability_text
          ? h('div', h('h3.rp-sub', 'Capability report'),
            h('pre.rp-pre', body.capability_text))
          : null);
      if (state.diagnosis) state.diagnosis.registered_as = sid;
      renderActions();
      if (onRegistered) onRegistered(sid);
    } catch (err) {
      if (!(err instanceof ServerRefusal)) throw err;
      fill(outEl, h('p.rp-refusal', err.message));
    } finally {
      busy(false);
    }
  }

  inspectBtn.addEventListener('click', doInspect);
  registerBtn.addEventListener('click', doRegister);

  /** Called by the browser pane on every listing, including a re-read.
   *
   *  THE OUTPUT IS CLEARED ONLY WHEN THE DIRECTORY CHANGES, and that is not a
   *  nicety. Registering ends by asking the browser pane to re-read the
   *  directory, so the newly-added section shows as registered in the list —
   *  which calls straight back into here. Clearing unconditionally therefore
   *  wiped the "Added as ..." confirmation and the capability report within a
   *  frame of them appearing: the action succeeded and the screen went blank.
   *  Caught by the end-to-end test in tests/test_browse_confinement.py, which
   *  is why that test clicks the real button instead of asserting on a
   *  payload. */
  function show({ path, diagnosis }) {
    const moved = (path || '') !== state.path;
    state.path = path || '';
    state.diagnosis = diagnosis || null;
    if (moved) clear(outEl);
    render();
  }

  renderEmpty();
  return { el, show, state };
}
