// inspect.js — the READ rail. Everything that interprets the active run, and
// nothing that changes what will be run (that is the SETUP rail, ui/spec.js).
//
// The twelve sections are organised into three named groups with a sticky
// header each and a navigator at the top of the rail:
//
//   RESULT   metrics, mechanism, caveats, provenance, reconstruction report
//   VIEW     channels, field overlay, window/level, render
//   MEASURE  line profile, regions, distance
//
// RESULT is first because it is the payoff of the primary verb. VIEW is second
// because it is used continuously. MEASURE is last because it is episodic. The
// navigator makes the order irrelevant for REACHING a section, so the order is
// free to express priority. Nothing is closed and nothing is hidden — see
// docs/UI_ARCHITECTURE.md section 2.4.

import { h, clear, fill } from '../util/dom.js';
import { CELL_TYPES } from '../gl/viewer.js';
import { MECH_META } from '../mechanism.js';
import { FIELD_NAMES, fieldMeta } from '../fields.js';
import { RAMP_NAMES, RAMP_META } from '../gl/ramps.js';
import { METRICS, metricProvenance, reportProvenance, fieldProvenance, chip, STANDING_LIMITS, claimSetNotice, claimSetEraLabel } from './provenance.js';
import { auto, num, int, pct, sci, dur, NULL_GLYPH } from '../util/fmt.js';
import { railGroup, railNav, sectTitle } from './rail.js';
import { createReachPanel } from './reach.js';
import { createFindingPanel } from './finding.js';

const CHAN_VARS = ['--ch-tumor', '--ch-stromal', '--ch-immune', '--ch-endothelial', '--ch-necrotic'];

function fmtMetric(key, v) {
  const m = METRICS[key];
  if (v === null || v === undefined) return NULL_GLYPH;
  switch (m && m.fmt) {
    case 'pct': return pct(v, 2);
    case 'sci': return sci(v, 3);
    case 'num1': return num(v, 1);
    default: return auto(v);
  }
}
function metricUnit(key) {
  const m = METRICS[key];
  if (!m) return '';
  return m.fmt === 'pct' ? '%' : m.unit;
}

/**
 * WHAT EACH METRIC IS A FRACTION *OF*, in cells, on the same line as the
 * percentage.
 *
 * Three of the four metrics are per-TUMOUR-CELL statistics and the panel used
 * to print them as bare percentages, with the only tumour-cell count on the
 * rail eight rows down in `Reconstruction report`. A reader who wants "11.2% of
 * what?" had to hold one number, scroll, and multiply — and this project has
 * already mis-stated one of its own findings by attaching a percentage to the
 * wrong denominator (see the header of tests/test_first_session.py, item 3).
 *
 * Nothing here is computed: the denominator is `report.n_tumor`, read off the
 * same RunResult, printed beside the percentage rather than derived from it.
 * The mapping is not a guess — each line names where the tumour-cell
 * denominator is asserted upstream:
 *
 *   kill_fraction            reconciles bit-exactly with the per-cell labels,
 *                            whose own killed_fraction_of_tumor divides by
 *                            n_tumor (tmesim/sim/mechanism.py)
 *   therapeutic_coverage_1x  "the fraction of tumor cells whose time-averaged
 *                            exposure exceeds mult * IC50"
 *                            (tmesim/validation/delivery_metrics.py)
 *   mean_tumor_auc           "the same statistics evaluated on the real tumor
 *                            CELLS" (same file)
 *
 * `penetration_length_um` is deliberately absent: it is a fitted decay distance
 * of a field, not a per-cell statistic, and inventing a denominator for it
 * would be worse than leaving the question open.
 */
const METRIC_DENOMINATOR = {
  kill_fraction: (n) => `of the ${n} tumour cells in this reconstruction`,
  therapeutic_coverage_1x: (n) => `of the ${n} tumour cells in this reconstruction`,
  mean_tumor_auc: (n) => `mean over the ${n} tumour cells in this reconstruction`,
};

/**
 * An absence, stated. The em dash on its own is the shape of a missing number,
 * a zero and a refusal all at once, and the reader cannot tell which — so the
 * glyph never appears in this rail without the sentence that says which one it
 * is.
 *
 * @param {string} reason  why there is no value. Never composed from a number.
 */
function why(reason) {
  return h('div.kv-why', reason);
}

/** The null placeholder, marked as one rather than left as loose punctuation. */
function nullValue() {
  return h('span.nullv', NULL_GLYPH);
}

/**
 * A titled band of quantities that share a provenance kind.
 *
 * The distinction between "an instrument measured this on this specimen" and
 * "the model produced this from something an instrument measured" is the
 * product's entire claim, and until this pass it was carried in the READ rail
 * by a chip's HUE and a chip's `title` — a colour and a hover. This is the
 * structural version: two boxes, each stating its kind in a word, a glyph and
 * a rule (css/base.css `.prov-band`).
 *
 * @param {'measured'|'derived'|'assumed'|'refused'} kind
 * @param {string} head  the band's own claim about itself, in words
 * @param {string} note  one sentence naming where that claim comes from
 */
function provBand(kind, head, note, ...body) {
  return h('div.prov-band.' + kind,
    h('h4.pb-head', head),
    // ABOVE the values, not under them. On a spot section this note is the
    // refusal that says the two coordinates below it were synthesised, and the
    // rule this rail already applies to a superseded claim set holds here too:
    // a reader who is going to be shown a number that is not what it looks like
    // has to be told before they read it, not after.
    note ? h('div.pb-note', note) : null,
    ...body);
}

/* ---------------------------------------------------------- caveats ------
 * The caveats panel must stay COMPLETE — every word of every caveat on screen,
 * nothing collapsed, nothing dismissible, nothing clamped. It does not have to
 * stay a paragraph dump, and as a paragraph dump it was losing: ten tinted
 * blocks of prose occupied more of the rail than the metrics did, so the panel
 * read as "this run is mostly warnings" long before anyone had read one.
 *
 * These strings are not arbitrary prose. They are written to a house shape:
 *
 *     <subject>: <VERDICT> (<citation>): <argument>
 *
 * with the subject usually an identifier (`oxygen_mode=markers`,
 * `cell_state (oxygen)`) or a shouted headline, and the citation always a
 * parenthesis of repo paths. Three fields, run together into one sentence.
 * splitCaveat() takes them back apart so they can be laid out as three lines
 * instead of parsed out of a paragraph by eye.
 *
 * Every branch is conservative: if a caveat does not match the shape it falls
 * through with `subject = null` and the whole string is rendered as the body.
 * Nothing is ever dropped — the citation is MOVED, not deleted, and the
 * original untouched string goes on the element's `title` so the verbatim
 * sentence is always one hover away from the restructured one.
 */

/** A caveat that says a number must not leave the screen, not merely that it
 *  is qualified. These are the only ones that keep a tinted field. */
const CV_HARD = /\bMUST NOT\b|\bNOT FOR CLINICAL USE\b|\bNOT DEFENSIBLE\b|\bNOT BE QUOTED\b|\bNEVER BE\b/;

function splitCaveat(text) {
  const raw = String(text).trim();
  let rest = raw;

  // 1. the citation. Always a parenthesis containing a repo path; inside a
  //    sentence it reads as noise, above one it reads as provenance.
  let ref = null;
  const m = rest.match(/\(([^()]*\.md[^()]*)\)/);
  if (m) {
    ref = m[1].trim();
    rest = (rest.slice(0, m.index) + rest.slice(m.index + m[0].length))
      .replace(/\s+([:;,.])/g, '$1')
      .replace(/\s{2,}/g, ' ')
      .trim();
  }

  // 2. the subject. Three shapes, tried in order of how sure we can be:
  //      a) an explicit "<label>: rest" prefix short enough to be a label
  //         rather than the first half of the sentence;
  //      b) a shouted headline — a run of at least TWO consecutive all-caps
  //         words ("ABSOLUTE KILL IS NOT DEFENSIBLE"). Two, because a single
  //         shouted word is emphasis inside a sentence, not a title: lifting
  //         "ABSOLUTE" out of "ABSOLUTE kill_fraction is an upper bound" gives
  //         a heading that says nothing and a body that starts mid-clause;
  //      c) a leading identifier — `key=value` or lower_snake_case, which in
  //         these strings is always the thing being qualified.
  let subject = null;
  const colon = rest.indexOf(':');
  if (colon > 0 && colon <= 88 && !/\.\s/.test(rest.slice(0, colon))) {
    subject = rest.slice(0, colon).trim();
    rest = rest.slice(colon + 1).trim();
  } else {
    const W = "[A-Z][A-Z0-9_'/-]*";
    const caps = rest.match(new RegExp(`^(${W}(?:[ ,]+${W}){1,18})`));
    const ident = rest.match(/^([a-z][a-z0-9]*(?:_[a-z0-9]+)+(?:=[A-Za-z0-9_.'-]+)?)(?=\s)/);
    if (caps && caps[1].length >= 8 && caps[1].length <= 96) {
      subject = caps[1].replace(/[,\s]+$/, '');
      rest = rest.slice(subject.length);
    } else if (ident) {
      subject = ident[1];
      rest = rest.slice(subject.length);
    }
  }
  // The split lands wherever the sentence happened to break, so the body can
  // start on the punctuation that used to join it to its heading.
  rest = rest.replace(/^[\s.,;:—–-]+/, '');

  // A subject that swallowed the entire caveat is not a subject.
  if (subject && !rest) { rest = subject; subject = null; }

  return { subject, ref, body: rest, raw, hard: CV_HARD.test(raw) };
}

/**
 * @param {string} text  the caveat, verbatim from the backend
 * @param {number} n     1-based index within its block
 * @param {string} extra extra class(es) — 'standing' for a project-level limit
 */
export function caveatNode(text, n, extra) {
  const c = splitCaveat(text);
  const cls = '.caveat' + (extra ? '.' + extra : '') + (c.hard ? '.sev-2' : '.sev-1');
  return h('div' + cls, { title: c.raw },
    c.subject
      ? h('div.cv-key',
        h('span.cv-n', String(n).padStart(2, '0')),
        h('span.cv-subj', c.subject))
      : null,
    c.ref ? h('div.cv-ref', c.ref) : null,
    h('div.cv-body',
      // With no subject line there is no numbered row to hang the index on, so
      // it leads the body instead. The list stays countable either way.
      c.subject ? null : h('span.cv-n', String(n).padStart(2, '0') + '  '),
      c.body),
  );
}

export function createInspectPanel({ viewer, onChange, onBenchmark, onResetView,
  windowLevelEl, profileEl, regionEl, distanceEl, mechanismEl, onIsolateMechanism }) {
  /* ------------------------------------------------------- channels ---- */
  const chanList = h('div');
  const chanDen = h('div.kv-why');
  const chanRows = [];

  function buildChannels(counts, total) {
    clear(chanList);
    chanRows.length = 0;
    const n = Number.isFinite(total) ? total : null;

    /* THE DENOMINATOR, ONCE, ABOVE THE COLUMN OF PERCENTAGES.
     * Each row used to read "12,345  58.2%" with no statement anywhere in the
     * panel of what the 58.2 was a percentage OF — the only cell count on this
     * screen is in the top bar, which is a different panel and can be showing a
     * different run. It is one sentence and it removes the arithmetic. */
    clear(chanDen);
    if (!counts) {
      chanDen.appendChild(h('span',
        'No geometry is loaded, so there are no cells to count. These five '
        + 'channels still toggle what the viewer would draw; the counts appear '
        + 'when a run\'s cells load over GET /api/runs/{id}/geometry.'));
    } else if (!n) {
      chanDen.appendChild(h('span',
        'The loaded geometry reports zero cells, so every share below would '
        + 'divide by zero and none is shown. The counts are the payload\'s own.'));
    } else {
      chanDen.appendChild(h('span',
        'Counts are cells in the loaded geometry. Every share below is of the '
        + 'same denominator: ', h('b', int(n)), ' cells.'));
    }

    CELL_TYPES.forEach((name, i) => {
      const cb = h('input', {
        type: 'checkbox', checked: true,
        'aria-label': `channel ${name}`,
        // Solo, from the keyboard. The gesture existed as alt-click only, which
        // is a mouse-only path to a control that has no other route: there is
        // no "show only this type" button anywhere in the app. Alt+Enter on the
        // focused channel does exactly what alt-click does, and both are named
        // on the row and in the panel hint rather than left to be discovered.
        onkeydown: (e) => {
          if (!e.altKey || (e.key !== 'Enter' && e.key !== ' ')) return;
          e.preventDefault();
          solo(i);
        },
        onchange: (e) => { viewer.setVisible(i, e.target.checked); syncRow(i); onChange && onChange(); },
      });
      const c = counts ? counts[i] : null;
      const share = (c !== null && c !== undefined && n) ? `${pct(c / n, 1)}%` : null;
      const row = h('label.chan', {
        title: c === null || c === undefined
          ? `${name} — click, or Space, to toggle this channel. `
            + 'Alt-click, or Alt+Enter, shows only this one. No cell counts '
            + 'are available until a geometry is loaded.'
          : `${name} — ${int(c)} of ${n ? int(n) : 'an unknown number of'} cells`
            + `${share ? ' (' + share + ')' : ''}. Click, or Space, to toggle `
            + 'this channel. Alt-click, or Alt+Enter, shows only this one.',
      },
        cb,
        h('span.sw', { style: { background: `var(${CHAN_VARS[i]})` } }),
        h('span.nm', name),
        // The unit is on the count and the share is only printed when there is
        // a denominator to divide by — "—%" was a null glyph wearing a unit,
        // which reads as a measured zero.
        c === null || c === undefined
          ? h('span.ct', nullValue())
          : h('span.ct', int(c), h('span.unit', 'cells'),
            share ? h('span', { style: { marginLeft: 'var(--sp-4)' } }, share) : null),
      );
      row.addEventListener('click', (e) => {
        if (!e.altKey) return;
        e.preventDefault();
        solo(i);
      });
      chanRows.push(row);
      chanList.appendChild(row);
    });
  }

  /** Show only channel `i`. The one gesture with no button anywhere else. */
  function solo(i) {
    CELL_TYPES.forEach((_, j) => {
      viewer.setVisible(j, j === i);
      const cb = chanRows[j].querySelector('input');
      if (cb) cb.checked = j === i;
      syncRow(j);
    });
    onChange && onChange();
  }
  function syncRow(i) {
    chanRows[i].classList.toggle('off', !viewer.visible[i]);
  }
  function setAll(on) {
    CELL_TYPES.forEach((_, i) => {
      viewer.setVisible(i, on);
      chanRows[i].querySelector('input').checked = on;
      syncRow(i);
    });
    onChange && onChange();
  }

  const chanHead = h('div.frow', { style: { marginBottom: 'var(--sp-2)' } },
    h('button.btn.ghost.sm', {
      title: 'Draw every cell type.', onclick: () => setAll(true),
    }, 'all'),
    h('button.btn.ghost.sm', {
      title: 'Draw no cell type. The pinned cell stays marked either way.',
      onclick: () => setAll(false),
    }, 'none'),
    h('span.grow'),
    // Both routes to the same gesture, named. It was "alt-click = solo", which
    // is a mouse-only instruction for a control with no other route.
    h('span.dim', {
      style: { fontSize: 'var(--fs-disclosure)' },
      title: 'Solo shows exactly one channel and hides the other four. On a '
        + 'focused channel checkbox, Alt+Enter (or Alt+Space) does the same '
        + 'thing as alt-clicking the row.',
    }, 'solo: alt-click · alt+enter'),
  );

  /* ---------------------------------------------------- field overlay -- */

  const fieldSel = h('select.inp', {
    id: 'f-overlay',
    onchange: () => onFieldChange(fieldSel.value || null),
  },
    h('option', { value: '' }, '— none'),
    FIELD_NAMES.map((n) => h('option', { value: n, title: fieldMeta(n).note },
      `${fieldMeta(n).short} · ${fieldMeta(n).unit}`)),
  );

  const rampSel = h('select.inp', {
    id: 'f-ramp',
    onchange: () => { viewer.setRamp(rampSel.value); onChange && onChange(); },
  }, RAMP_NAMES.map((n) => h('option', { value: n, title: RAMP_META[n].note }, RAMP_META[n].label)));

  const opacity = h('input', {
    type: 'range', min: '0', max: '100', value: '55', id: 'f-opacity',
    'aria-label': 'Field overlay opacity',
    oninput: () => { viewer.fieldOpacity = opacity.valueAsNumber / 100; opacityVal.textContent = opacity.value; viewer.dirty = true; },
  });
  const opacityVal = h('span.mono', '55');

  const gamma = h('input', {
    type: 'range', min: '30', max: '250', value: '100', id: 'f-gamma',
    'aria-label': 'Field ramp gamma',
    oninput: () => { viewer.fieldGamma = gamma.valueAsNumber / 100; gammaVal.textContent = (gamma.valueAsNumber / 100).toFixed(2); viewer.dirty = true; },
  });
  const gammaVal = h('span.mono', '1.00');

  // Colour source. This was a two-state checkbox labelled "Colour cells / by
  // field value, not type" while the stage toolbar carried a THREE-state
  // segmented control over the same viewer property. Two vocabularies for one
  // setting, and the checkbox could not express `mech` at all — with the viewer
  // in mechanism mode it read as unchecked, i.e. as "colouring by type", which
  // was simply false. Same control, same three states, same words as the
  // toolbar. The rail keeps it because this is where the field controls are.
  const COLOUR_MODES = [
    ['0', 'type', 'Colour each cell by its measured phenotype call (C)'],
    ['1', 'field', 'Colour each cell by the field value at its position, through the current window (C)'],
    ['2', 'mech', 'Colour each tumour cell by WHY it died or survived — the per-cell mechanism attribution (C)'],
  ];
  const colourBtns = COLOUR_MODES.map(([mode, label, title]) => h('button', {
    type: 'button', title, 'aria-pressed': mode === '0' ? 'true' : 'false',
    dataset: { mode },
    onclick: () => { viewer.setColorMode(Number(mode)); onChange && onChange(); },
  }, label));
  const colourSeg = h('div#f-colour.seg', { role: 'group', 'aria-label': 'Cell colour source' }, colourBtns);

  const fieldStats = h('dl.kv');
  const fieldNote = h('div.field-note');

  let onFieldChange = () => {};
  // The overlay currently loaded, kept so the SELECTED CELL panel can say which
  // field the value under the cell belongs to, and in what unit. `viewer.field`
  // is the raw grid and carries neither.
  let currentField = null;

  function setFieldInfo(name, f) {
    const meta = fieldMeta(name || '');
    if (!name || !f) {
      // Not "state: —". Nothing is wrong and nothing is missing: no overlay is
      // selected, which is a choice the operator has made and can unmake with
      // the control four rows above this line.
      fill(fieldStats, why(
        'No field overlay is selected, so there is no grid to describe. Pick '
        + 'one in Field, above, to load it over GET /api/runs/{id}/field/{name} '
        + 'and read its range, its pitch and its provenance here.'));
      clear(fieldNote);
      currentField = null;
      renderCell();
      return;
    }
    const fp = fieldProvenance(name);
    // payloads.py is explicit that a grid with NO finite entry reports nan for
    // both ends of its range, and that "the front end must treat a non-finite
    // range as 'no data', not as 0..1". `auto(NaN)` returns the null glyph, so
    // before this the panel printed "— vol frac" — a dash wearing a unit, which
    // is exactly the shape of a measured zero.
    const finiteRange = Number.isFinite(f.min) && Number.isFinite(f.max);
    fill(fieldStats,
      h('dt', 'min'), h('dd', finiteRange ? auto(f.min) : nullValue(),
        finiteRange ? h('span.unit', meta.unit) : null),
      h('dt', 'max'), h('dd', finiteRange ? auto(f.max) : nullValue(),
        finiteRange ? h('span.unit', meta.unit) : null),
      finiteRange ? null : why(
        'This grid carries no finite value, so it has no range. The payload '
        + 'reports nan for both ends rather than substituting 0..1 '
        + '(app/server/payloads.py), and nothing is drawn from it.'),
      h('dt', 'grid'), h('dd', `${int(f.rows)} × ${int(f.cols)}`,
        h('span.unit', 'nodes')),
      h('dt', 'pitch'), h('dd', num(f.dx, 2), h('span.unit', 'um/node')),
      h('dt', 'provenance'), h('dd', chip(fp.kind, fp.why)),
      // Where min and max come from, which is not the same question as where
      // the field comes from: they are the payload header's own finite range
      // over the whole grid, not a statistic of the tissue under it.
      why('min and max are the payload header\'s range over the FINITE entries '
        + 'of the whole grid — every node, tissue or not — so they are the '
        + 'colour-ramp bounds, not a statistic of this section.'));
    // --fs-disclosure, not the stylesheet's --fs-micro. `meta.note` is prose,
    // and for the oxygen field it is the sentence that says the whole overlay
    // is a SURROGATE. A limit set in the size reserved for unit suffixes is a
    // limit that reads as a footnote; tokens.css names the floor.
    fill(fieldNote, h('span', { style: { fontSize: 'var(--fs-disclosure)' } },
      `${meta.source}. ${meta.note}`));
    currentField = { name, f, meta };
    renderCell();
  }

  /* ---------------------------------------------------------- render -- */

  const sizeBoost = h('input', {
    type: 'range', min: '30', max: '400', value: '100', id: 'r-size',
    'aria-label': 'Point size',
    oninput: () => { viewer.sizeBoost = sizeBoost.valueAsNumber / 100; sizeVal.textContent = (sizeBoost.valueAsNumber / 100).toFixed(2); viewer.dirty = true; },
  });
  const sizeVal = h('span.mono', '1.00');
  const alphaR = h('input', {
    type: 'range', min: '10', max: '100', value: '100', id: 'r-alpha',
    'aria-label': 'Cell opacity',
    oninput: () => { viewer.pointAlpha = alphaR.valueAsNumber / 100; alphaVal.textContent = alphaR.value; viewer.dirty = true; },
  });
  const alphaVal = h('span.mono', '100');
  const lodChk = h('input', {
    type: 'checkbox', checked: true, id: 'r-lod',
    onchange: () => { viewer.lodAuto = lodChk.checked; if (!lodChk.checked) viewer.lod = 1; viewer.dirty = true; },
  });
  const roundChk = h('input', {
    type: 'checkbox', checked: true, id: 'r-round',
    onchange: () => { viewer.roundPoints = roundChk.checked; viewer.dirty = true; },
  });
  const benchOut = h('div', { style: { fontSize: 'var(--fs-micro)', fontFamily: 'var(--font-mono)', color: 'var(--fg-2)', lineHeight: '1.45' } });

  /* --------------------------------------------------------- metrics -- */

  const metricsEl = h('div');
  // Clicking a survival row on THE FINDING isolates that bucket in the tissue
  // view. The callback is supplied by main.js rather than done here, because
  // switching colour mode also has to sync the top-bar control -- doing it
  // directly on the viewer would leave the picture and the control disagreeing
  // about what is on screen.
  const findingPanel = createFindingPanel({
    onIsolate: onIsolateMechanism || null,
  });
  const reachPanel = createReachPanel();
  const provEl = h('dl.kv');
  const caveatsEl = h('div.caveats');
  const tallyEl = h('span.cv-tally');
  const reportEl = h('dl.kv');
  const errorEl = h('div.hidden');

  /* --------------------------------------------------- selected cell ----
   * WHY A CELL IS IN THE STATE IT IS IN.
   *
   * Clicking a cell already did something: the stage's top-left card names its
   * phenotype, its coordinates, the field value under the pointer and the local
   * composition. What it never said is the thing the click is actually asked in
   * order to find out — WHY this cell is where the run put it — and it could
   * not say it, because that answer is three payloads wide: the mechanism label
   * (TMSM), the band thresholds that produced the label (the same payload's
   * JSON block), and the exposure field the bands were evaluated on (TMSF).
   *
   * It also could not make the distinction this product exists to make. On the
   * stage card a cell's POSITION and its FIELD VALUE are two rows of the same
   * `dl` in the same type at the same weight, and one of them came off a
   * segmentation while the other is a bilinear interpolation of a solved grid.
   * Here they are in two separately-titled bands — see provBand() — so the
   * measured/derived split is structural rather than a chip's hue.
   *
   * ON A SPOT SECTION THE FIRST BAND IS NOT A MEASUREMENT AT ALL. A Visium
   * reconstruction has no measured cells: it PLACES them from the per-spot
   * deconvolved fractions at spot centroids with a seeded jitter
   * (app/server/cache.py::_visium_metadata, tmesim/reconstruct/real_pipeline.py
   * ::reconstruct_real_visium), which is why app/server/api.py refuses per-cell
   * spatial statistics on every `visium` section outright. `report.source` is
   * the run's own record of which reconstructor produced it, so this panel
   * reads that rather than guessing from a section id.
   *
   * Nothing here is computed in the browser. Every number is read off the
   * viewer's own payload buffers and every sentence about the mechanism is the
   * backend's `label_descriptions` / `threshold_derivation` string. */
  const cellEl = h('div');
  let currentRun = null;

  /** The run's own attribution block, or null. */
  function mechAggregate() {
    return (currentRun && currentRun.mechanism) || null;
  }

  /**
   * Is this run's per-cell geometry MEASURED, or synthesised from spots?
   * @returns {{spot: boolean, source: string|null}}
   */
  function geometryOrigin() {
    const src = currentRun && currentRun.report && currentRun.report.source
      ? String(currentRun.report.source).toLowerCase() : null;
    return { spot: src === 'visium', source: src };
  }

  function cellEmpty(msg) {
    fill(cellEl, h('div', {
      style: {
        fontSize: 'var(--fs-disclosure)', lineHeight: 'var(--lh-body)',
        color: 'var(--fg-2)',
      },
    }, msg));
  }

  function renderCell() {
    // The stub viewer (no WebGL2) has no picking at all, and neither does a
    // viewer with no geometry loaded. Both are ordinary states, and each says
    // which one it is rather than sharing one blank panel.
    if (typeof viewer.cellInfo !== 'function') {
      cellEmpty('The cell viewer is unavailable on this workstation (no '
        + 'WebGL2), so no cell can be picked. Every other panel in this rail '
        + 'still reads the run.');
      return;
    }
    if (!viewer.geom) {
      cellEmpty('No cell geometry is loaded, so there is nothing to pick. Open '
        + 'a finished run to load its cells.');
      return;
    }
    const sel = Number.isInteger(viewer.selected) ? viewer.selected : -1;
    if (sel < 0) {
      cellEmpty('No cell pinned. Click one: this panel then shows what was '
        + 'measured on it, what the model derived, and its mechanism band. '
        + 'Click empty tissue, or unpin on the stage toolbar, to clear.');
      return;
    }
    let info = null;
    try {
      info = viewer.cellInfo(sel);
    } catch {
      info = null;   // a payload swapped out from under the selection
    }
    if (!info) {
      cellEmpty('The pinned cell is no longer in the loaded geometry — the '
        + 'stage has been given a different run\'s cells since it was pinned. '
        + 'Nothing is shown for it, because the index would name a different '
        + 'cell.');
      return;
    }

    clear(cellEl);
    const origin = geometryOrigin();

    /* ---- band 1: what an instrument put here, or what did not ---------- */
    const posRows = h('dl.kv',
      h('dt', 'phenotype'), h('dd', info.typeName,
        chip(origin.spot ? 'computed' : 'measured', origin.spot
          ? 'On a spot section a cell TYPE is a deconvolution product, not a segmentation call.'
          : 'Cell phenotype calls come from the segmented spatial assay.')),
      h('dt', 'x'), h('dd', num(info.x, 1), h('span.unit', 'um')),
      h('dt', 'y'), h('dd', num(info.y, 1), h('span.unit', 'um')),
      why('Coordinates are micrometres in the reconstruction\'s own domain '
        + 'frame, with the origin at the corner of the grid — not slide '
        + 'coordinates and not a distance from anything.'));

    if (origin.spot) {
      cellEl.appendChild(provBand('refused',
        'not measured — synthesised from a spot',
        'This section is spot-level, so this cell was PLACED by deconvolution '
        + 'at a spot centroid with a seeded jitter; nothing segmented it and '
        + 'nothing measured where it sits. Per-cell spatial statistics are '
        + 'refused outright on such a section — the full refusal is on the '
        + 'section itself, in the picker on the Setup rail '
        + '(app/server/api.py::SPOT_SPATIAL_REFUSAL). Read the coordinates '
        + 'below as the synthesis, not as this tissue.',
        posRows));
    } else {
      cellEl.appendChild(provBand('measured',
        'measured on this specimen',
        origin.source
          ? `Reconstructed by the ${origin.source} route: positions and `
            + 'phenotype calls come from the segmented spatial assay, and '
            + 'nothing in the model chain moved them.'
          : 'This run does not record which reconstructor produced it, so '
            + 'whether these two rows are segmented or synthesised is not '
            + 'established from this screen.',
        posRows));
    }

    /* ---- band 2: what the model made of it ----------------------------- */
    const derived = [];
    const agg = mechAggregate();
    const code = (info.label === null || info.label === undefined)
      ? null : Number(info.label);
    const meta = (code !== null && MECH_META[code]) ? MECH_META[code] : null;

    if (meta) {
      const desc = (agg && agg.label_descriptions
        && agg.label_descriptions[meta.name]) || null;
      derived.push(h('div', {
        style: {
          display: 'flex', alignItems: 'baseline', gap: 'var(--sp-4)',
          marginBottom: 'var(--sp-2)',
        },
      },
      h('span', {
        style: {
          width: '9px', height: '9px', flex: '0 0 auto',
          background: `var(${meta.var})`,
          border: '1px solid color-mix(in srgb, var(--fg-0) 30%, transparent)',
        },
      }),
      h('span.mono.strong', meta.name),
      h('span.grow'),
      chip('computed', 'The per-cell mechanism attribution, decoded from the '
        + 'TMSM payload for this run. It is the engine\'s own label, and the '
        + 'whole label array reconciles bit-exactly with kill_fraction.')));

      /* FROM THIS ONE CELL TO EVERY CELL LIKE IT.
       *
       * The third leg of the round trip. THE FINDING goes from a category to
       * the cells (click a survival row, the specimen shows only that band);
       * the Mechanism section goes back to the claim (the button above it).
       * This goes from ONE cell to its whole category — which is the question
       * a reader has the moment they read a label off a cell they picked:
       * where are the others?
       *
       * It routes through the SAME callback the finding rows use
       * (main.js onIsolateMechanism), so the colour mode, the top-bar control,
       * the mechanism legend and the stage cannot disagree about what is on
       * screen — and the legend's filter band then states the filter and
       * offers the one action that clears it. Same words as the finding row's
       * own title, because it is the same action (docs/UI_ARCHITECTURE.md 5).
       *
       * Rendered only when a callback exists, so the stub viewer and any
       * caller that did not wire one get no dead control. */
      if (onIsolateMechanism) {
        derived.push(h('div.frow', { style: { marginBottom: 'var(--sp-4)' } },
          h('button.btn.ghost.sm', {
            type: 'button',
            title: 'Show only these cells in the tissue view — every cell the '
              + `run labelled ${meta.name}, and no others. It changes the `
              + 'picture only; no number on this rail moves. The Mechanism '
              + 'section, which this scrolls to, undoes it with "show all '
              + 'categories".',
            onclick: () => onIsolateMechanism(code),
          }, 'show only these cells'),
          h('span.grow')));
      }
      derived.push(h('div', {
        style: {
          fontSize: 'var(--fs-sm)', lineHeight: 'var(--lh-body)',
          color: 'var(--fg-1)', marginBottom: 'var(--sp-4)',
        },
      }, desc || 'This run\'s attribution payload carries no description for '
        + 'this label, so only the label\'s name is shown. Nothing is '
        + 'substituted for the sentence.'));

      // WHICH BAND, AND WHERE ITS EDGES ARE. The label is a bucket, and a
      // bucket without its boundaries is a verdict without its test. Both the
      // thresholds and the sentence explaining how each was arrived at are the
      // backend's (tmesim/sim/mechanism.py::MechanismThresholds.resolve).
      const th = (agg && agg.thresholds) || null;
      if (th) {
        /* A DASH WEARING A UNIT IS THE SHAPE OF A MEASURED ZERO.
         *
         * These three rows printed `num(v, 2)` — which returns the em dash for
         * an absent value — and then appended the unit span unconditionally,
         * so a threshold the payload does not carry rendered as `— × kill
         * prob.`: a null glyph with a unit on it, in a list of band edges,
         * with no sentence anywhere saying which one was missing. This panel
         * names that exact failure twice elsewhere in this file ("`—%` was a
         * null glyph wearing a unit, which reads as a measured zero") and was
         * still committing it here.
         *
         * Found 2026-09-20 by rendering the panel rather than by reading it.
         * The unit now travels only with a value, the glyph is marked as one,
         * and every absent edge is named in a sentence under the list — which
         * matters more here than in most of the rail: the cell's LABEL was
         * decided by these edges, so an edge that is not on screen is a
         * verdict whose test is not on screen. */
        const edge = (label, v, unitTxt, field) => [
          h('dt', label),
          h('dd', Number.isFinite(v) ? num(v, 2) : nullValue(),
            Number.isFinite(v) ? h('span.unit', unitTxt) : null),
          Number.isFinite(v) ? null
            : why(`No ${label} edge: this run's attribution payload carries no `
              + `thresholds.${field}. The band above was resolved with one, so `
              + 'the edge exists — it is this record that does not state it. '
              + 'Read the dash as absent, never as zero.'),
        ];
        derived.push(h('dl.kv',
          edge('unreached below', th.unreached_ratio, 'E/IC50', 'unreached_ratio'),
          edge('underdosed below', th.adequate_ratio, 'E/IC50', 'adequate_ratio'),
          edge('resistant at or below', th.resistant_modulation,
            '× kill prob.', 'resistant_modulation'),
          why('The bands this run sorted every surviving tumour cell into. '
            + 'E is the cell\'s own time-averaged exposure and IC50 is the '
            + 'compound\'s, so E/IC50 is a ratio and carries no unit of its '
            + 'own. These are this run\'s resolved values, not defaults.')));
      }
    } else if (code !== null) {
      derived.push(why(`This cell carries mechanism code ${code}, which this `
        + 'client has no name for. The attribution payload and this build '
        + 'disagree about the label set; nothing is guessed.'));
    } else {
      derived.push(why(currentRun && currentRun.mechanism_reason
        ? 'No mechanism label for this cell: ' + currentRun.mechanism_reason
        : 'No per-cell mechanism attribution is loaded for this run, so there '
          + 'is no answer to "why is this cell in this state". A '
          + 'reconstruction-only run has none by construction; otherwise the '
          + 'status bar carries the backend\'s reason.'));
    }

    /* the local field value — the other half of "why" */
    if (currentField && currentField.f) {
      const fp = fieldProvenance(currentField.name, currentRun);
      const v = info.fieldValue;
      const hasV = v !== null && v !== undefined && Number.isFinite(v);
      derived.push(h('dl.kv',
        h('dt', currentField.meta.short), h('dd',
          hasV ? auto(v) : nullValue(),
          hasV ? h('span.unit', currentField.meta.unit) : null,
          chip(fp.kind, fp.why)),
        hasV
          ? why(`${currentField.meta.label} at this cell, bilinearly `
            + `interpolated between the four nearest nodes of the ${int(currentField.f.rows)}`
            + ` × ${int(currentField.f.cols)} solved grid at `
            + `${num(currentField.f.dx, 2)} um per node. It is a value OF THE `
            + 'FIELD read at the cell\'s position, never a per-cell '
            + 'measurement of the cell.')
          : why('This cell lies outside the loaded field grid, so there is no '
            + 'value to sample. Nothing is extrapolated to the edge.')));
    } else {
      derived.push(why('No field overlay is loaded, so the local field value '
        + 'this cell\'s band was decided on is not on screen. Load one in '
        + 'Field overlay, under VIEW, to read it here.'));
    }

    /* The band head and the chips inside it use two words for one idea —
     * "derived" on the band, "computed" on every chip — because they come from
     * two vocabularies: `.prov-band` states the KIND OF EVIDENCE, and
     * js/ui/provenance.js::KINDS states HOW A SINGLE VALUE WAS ARRIVED AT.
     * Both words are load-bearing where they are and neither can be renamed
     * from this file, so the note reconciles them once, here, rather than
     * leaving a reader to decide whether a "computed" row inside a "derived"
     * band is a third thing. */
    cellEl.appendChild(provBand('derived',
      'derived by the model',
      'Everything in this band is a model output computed from the band above '
      + 'plus the RunSpec. None of it was measured on this cell, and none of '
      + 'it has been scored against a measured response. The chips below read '
      + '"computed" or "assumed" — that is how each single value was arrived '
      + 'at; "derived" is what the whole band is. Nothing in it is measured, '
      + 'whichever chip it wears.',
      ...derived));
  }

  const headRunEl = h('span.mono.strong', {
    style: { marginLeft: 'var(--sp-3)', flex: '0 0 auto' },
  }, 'no run');

  function renderRun(run) {
    currentRun = run || null;
    findingPanel.render(run || null);
    reachPanel.render(run && run.status === 'done' ? run : null);
    // Two words, not a dash. This slot answers "whose numbers am I reading",
    // and `—` is not an answer to it — it is the same glyph the rail uses for
    // a missing measurement, in the one place where nothing is missing.
    headRunEl.textContent = (run && run.run_id) || 'no run';
    headRunEl.title = run
      ? `Every number in this rail was produced by ${run.run_id}.`
      : 'No run is loaded, so every panel in this rail is showing its empty '
        + 'state rather than a stale result.';

    /* metrics */
    clear(metricsEl);
    if (!run) {
      metricsEl.appendChild(h('div.dim', { style: { fontSize: 'var(--fs-xs)', padding: '4px 0' } },
        'No run selected. Metrics appear once a run reaches status=done.'));
    } else if (run.status === 'failed') {
      // NOT the error text again. `errorEl` sits immediately above this block
      // (see railBody) and already carries `run.error` verbatim, so printing it
      // here rendered the same red sentence twice, stacked, in one panel — an
      // operator reading two identical failures had to compare them character
      // by character to find out they were one. This slot answers the question
      // the panel is asking, which is where the metrics went.
      metricsEl.appendChild(h('div.caveat.crit', h('span', run.error
        ? 'This run failed, so it produced no metrics. The server’s reason is '
          + 'directly above.'
        : 'This run failed, so it produced no metrics — and the server sent no '
          + 'reason with it. The failure is recorded on the console of the '
          + 'process running Vivum.')));
    } else if (!run.metrics) {
      metricsEl.appendChild(h('div.dim', { style: { fontSize: 'var(--fs-xs)', padding: '4px 0' } },
        `Run ${run.run_id} is ${run.status}; no metrics yet.`));
    } else {
      /* NO METRIC IN THIS TABLE GETS `.lead` ANY MORE. 2026-08-22.
       *
       * It used to: "the first entry is the one the run was submitted to
       * answer", set at --fs-xl while the rest stayed at --fs-lg, so the panel
       * had "a place for the eye to land instead of four equal numbers". The
       * reasoning was sound and the result was the worst defect found in a
       * first session with this product.
       *
       * The first entry is `kill_fraction`. So the first and largest number on
       * the READ rail — the thing a reader who has never seen this Studio
       * lands on, at the largest type size the panel has — was an absolute
       * kill percentage, on a section whose own caveat four sections below
       * reads "ABSOLUTE KILL MUST NOT BE QUOTED FROM THIS SECTION UNDER ANY
       * STATE ROUTE". A stranger reads top-down and quotes the one number this
       * project forbids. The product exists to prevent exactly that failure
       * and was committing it in its own UI, in the largest available type.
       *
       * The fix is NOT to remove the number: it is a real output, it is the
       * denominator of the survivor arithmetic, and hiding it would be its own
       * dishonesty. The fix is to stop it being where the eye lands. The eye
       * now lands on the FINDING block at the top of this group, which states
       * the defensible claim — the delivery/potency split among survivors, in
       * counts of cells — and which cannot print an absolute kill rate as a
       * headline by construction.
       *
       * Why no OTHER metric was promoted into the vacancy: there is no
       * defensible candidate. All four read an absolute level. `kill_fraction`
       * and `therapeutic_coverage_1x` are no-resistance upper bounds that read
       * the IC50; `mean_tumor_auc` scales with the plasma level; and
       * `penetration_length_um` is a fitted parameter of a gradient whose
       * transport length this repo measures as 2.45x too long -- about
       * 2.2-3.0x on this product's own map, NOT the 14x an earlier reading
       * gave (results/kyle2004_penetration.md sec. 6, which corrects it: the
       * 14x was a ratio of two dimensionless ratios run through an
       * exponential estimator, and that amplification does not transfer
       * here).
       * Promoting any of them would move the violation rather than fix it.
       * Four equal numbers is the honest layout for four quantities of which
       * none is the answer. */
      const nTumor = (run.report && Number.isFinite(run.report.n_tumor))
        ? run.report.n_tumor : null;
      const reasons = run.metric_reasons || {};
      for (const key of Object.keys(METRICS)) {
        const v = run.metrics[key];
        const m = METRICS[key];
        const p = metricProvenance(key, run);
        const isNull = v === null || v === undefined;
        // THE DENOMINATOR, BESIDE THE PERCENTAGE. See METRIC_DENOMINATOR above
        // for which count each metric divides by and where that is asserted.
        // Printed only when the run actually carries the count — this composes
        // no number and substitutes none.
        const den = (!isNull && nTumor !== null && METRIC_DENOMINATOR[key])
          ? METRIC_DENOMINATOR[key](int(nTumor)) : null;
        // A null metric is a REFUSAL with a recorded reason, not an unknown.
        // jobs.py writes `metric_reasons[name]` every time it declines to
        // report one ("requires a drug; this RunSpec has drug=null", "no TUMOR
        // cells in this reconstruction", "the exp+floor fit did not converge").
        // The rail dropped all of it and printed a dash.
        const reason = isNull ? (reasons[key] || null) : null;
        const metaNote = (text) => h('span', {
          // Not --fs-micro: these are the sentences that stop the percentage
          // above them being attached to the wrong population, or a dash being
          // read as a zero.
          style: {
            gridColumn: '1 / -1', fontFamily: 'var(--font-ui)',
            fontSize: 'var(--fs-disclosure)', color: 'var(--fg-1)',
            paddingTop: 'var(--sp-2)',
          },
        }, text);
        metricsEl.appendChild(h('div.metric',
          h('div.m-name', m.label),
          h('div.m-val' + (isNull ? '.null' : ''),
            isNull ? nullValue() : fmtMetric(key, v),
            isNull ? null : h('span.unit', metricUnit(key))),
          h('div.m-meta',
            chip(p.kind, p.why),
            h('span.m-src', m.src),
            den ? metaNote(den) : null,
            isNull ? metaNote(reason
              ? 'No value: ' + reason + '. Nothing is substituted for it.'
              : 'No value, and this run recorded no reason for it in '
                + 'metric_reasons. That is a backend omission; do not read the '
                + 'dash as a zero.') : null),
        ));
        /* THE HARD STOP, ON THE NUMBER, NOT A THOUSAND WORDS BELOW IT.
         *
         * `m.hard_stop` marks a metric whose own caveat says it must not leave
         * this screen at all, as distinct from one that is merely qualified.
         * The full statement is still in Caveats, in full, unedited — this
         * removes nothing and adds a pointer that cannot be scrolled past,
         * because it is rendered between this metric and the next one. A
         * caveat nobody reads protects nobody, and a reader who quotes a
         * number never reaches the block that forbids it. */
        if (m.hard_stop && !isNull) {
          metricsEl.appendChild(h('div.caveat.crit', { title: p.why },
            h('span', h('b', 'DO NOT QUOTE THIS NUMBER. '), m.hard_stop)));
        }
      }
      const extra = Object.keys(run.metrics).filter((k) => !(k in METRICS));
      for (const key of extra) {
        const v = run.metrics[key];
        // `=== null` missed `undefined`, which rendered as the null glyph in a
        // cell not marked null. Both are "the backend did not give a value".
        const isNull = v === null || v === undefined;
        metricsEl.appendChild(h('div.metric',
          h('div.m-name', key),
          h('div.m-val' + (isNull ? '.null' : ''), isNull ? nullValue() : auto(v)),
          h('div.m-meta',
            chip('computed', 'Metric not described by the front-end table; shown verbatim from the backend.'),
            // No unit, because this client has no table entry for this key and
            // guessing one would be worse than saying so.
            h('span.m-src', isNull
              ? (reasons[key] || 'No value, and no reason recorded for it.')
              : 'Not described by this client\'s metric table, so it is shown '
                + 'verbatim from the backend WITHOUT a unit — the unit is not '
                + 'carried on the wire and is not guessed here.')),
        ));
      }
      // The delivery law, on the same screen as the numbers it produced. Every
      // metric above is an integral of one transport solve, and that solve's
      // vascular source term is a chosen default that is measured to change
      // which drug comes first. Putting it only in Caveats — one section down —
      // would let a delivery number be read without it. The full statement,
      // with its provenance, is still down there; this is the pointer that
      // cannot be missed. Text and law come from the backend's runtime read of
      // PDEConfig (app/server/delivery_law.py), never from a constant here.
      const dlaw0 = (run.provenance || {}).delivery_law;
      if (dlaw0 && dlaw0.headline) {
        metricsEl.appendChild(h('div.caveat.crit', { title: dlaw0.caveat },
          h('span', dlaw0.headline,
            ' Full statement in Caveats below.')));
      }
    }

    /* error surface */
    clear(errorEl);
    if (run && run.error) {
      errorEl.classList.remove('hidden');
      errorEl.appendChild(h('div.caveat.crit', h('span', String(run.error))));
    } else {
      errorEl.classList.add('hidden');
    }

    /* provenance block */
    const pv = (run && run.provenance) || null;
    if (!pv) {
      // "state —". Three different situations wore that one dash: no run open,
      // a run that has not finished, and a finished run whose record carries no
      // provenance block at all — and only the third is a defect.
      fill(provEl, why(!run
        ? 'No run is loaded, so there is nothing to state the provenance of. '
          + 'Submit the RunSpec on the left, or open a run from the history '
          + 'dock below.'
        : run.status !== 'done'
          ? `Run ${run.run_id} is ${run.status}. The provenance block is written `
            + 'when the job finishes, so there is none yet.'
          : `Run ${run.run_id} finished without a provenance block. Every number `
            + 'in this rail is therefore unattributed, which is a backend '
            + 'omission rather than a result — do not read the absence as '
            + '"nothing to declare".'));
    } else {
      // IC50 is LITERATURE, not MEASURED. It is a real measurement — on GDSC
      // cell lines, not on this section. Wearing the same green chip as the
      // cell count implies both came off this specimen; only one did.
      const ic50Measured = pv.ic50_source && /measured/i.test(pv.ic50_source);
      const hasIc50 = pv.ic50_uM !== null && pv.ic50_uM !== undefined;
      // The delivery law leads the block. It is not a detail of the run — it is
      // the source term of the transport solve every number below is integrated
      // out of, and the only input measured to reorder the drug panel. The
      // string is the backend's, read off PDEConfig at run time
      // (app/server/delivery_law.py), so this row cannot describe a law the
      // engine has stopped using.
      const dlaw = pv.delivery_law || null;
      fill(provEl,
        ...(dlaw ? [
          h('dt', 'delivery law'),
          h('dd', h('span.mono', dlaw.law),
            chip('assumed', dlaw.caveat)),
        ] : []),
        // Dimensionless, and it never said so: a bare "0.4312" beside four
        // quantities that all carry units reads as a fifth measured value.
        h('dt', 'alpha'), h('dd', auto(pv.alpha), h('span.unit', 'dimensionless'),
          chip('assumed', pv.alpha_source
            ? `${pv.alpha_source}\n\nFitted once against a published penetration distance, then held fixed; it is not re-measured per section.`
            : 'Tortuosity calibration constant. Fitted once, then held fixed; it is not re-measured per section.')),
        h('dt', 'IC50'), h('dd',
          hasIc50 ? auto(pv.ic50_uM) : nullValue(),
          hasIc50 ? h('span.unit', 'uM') : null,
          chip(ic50Measured ? 'literature' : 'assumed',
            `${pv.ic50_source || 'no source given'}${pv.ic50_note ? '\n\n' + pv.ic50_note : ''}`
            + '\n\nMeasured in cell culture, NOT on this section. This tissue has no measured drug response.')),
        // The backend records WHY there is no potency (jobs.py writes
        // `ic50_reason` on a reconstruction-only run). The dash used to be the
        // whole message, and a missing IC50 is the difference between "this run
        // had no drug" and "we could not find one for this drug".
        hasIc50 ? null : why(pv.ic50_reason
          ? 'No IC50: ' + pv.ic50_reason
          : 'No IC50 was resolved for this run and no reason was recorded. '
            + 'Every potency-dependent number above is therefore undefined '
            + 'rather than zero.'),
        h('dt', 'runtime'), h('dd', dur(run.runtime_s), chip('measured', 'Wall-clock time the backend spent on this job.')),
      );
      if (pv.ic50_source) {
        // Full width and wrapping, at the disclosure floor. It was a --fs-micro
        // `dd`, and a `.kv > dd` is nowrap with an ellipsis — so the one line
        // saying WHERE the potency came from was both the smallest text in the
        // panel and truncated, with the rest of the sentence on a hover.
        provEl.appendChild(why('IC50 source: ' + pv.ic50_source));
      }
      if (pv.state_mode) {
        provEl.appendChild(h('dt', 'state route'));
        provEl.appendChild(h('dd', pv.state_mode,
          chip(pv.state_is_measured && !(pv.state_unmeasured || []).length ? 'measured' : 'assumed',
            pv.state_drives_kill === false
              ? String(pv.state_drives_kill_reason || '')
              : String(pv.state_source || ''))));
      }
      if (pv.sensitivity_note) {
        provEl.appendChild(h('dt', 'sensitivity'));
        // Sans: this is the one value in the block that is a paragraph rather
        // than a quantity, and in the mono face it was three lines of code-
        // coloured text hanging off a right-aligned key/value list.
        provEl.appendChild(h('dd', {
          title: pv.sensitivity_note,
          style: {
            // The sensitivity note is where the upper-bound reading of every
            // kill number is stated — the product's own thesis about which of
            // its outputs is defensible. It used to say "COMPARABLE across
            // drugs"; that comparison was measured and withdrawn on
            // 2026-08-21, and the note now says what survives it. Either way
            // it is prose stating a limit, so it takes the disclosure floor.
            fontFamily: 'var(--font-ui)', fontSize: 'var(--fs-disclosure)',
            whiteSpace: 'normal', textAlign: 'left', justifyContent: 'flex-start',
            lineHeight: 'var(--lh-body)', padding: '3px 0',
          },
        }, pv.sensitivity_note));
      }
    }

    /* caveats — always fully visible, never truncated, never behind a toggle.
       Restructured, not reduced: see splitCaveat() above. */
    clear(caveatsEl);
    const cav = (pv && pv.caveats) || [];
    const hard = [...cav, ...STANDING_LIMITS].filter((c) => CV_HARD.test(String(c))).length;
    // Two ways a record can put a withdrawn claim on this screen: it predates
    // the current claim set (`superseded`), or it is stamped with a claim set
    // its own caveats disagree with (`contradicted` — written by a build that
    // had not caught up). Both mark; the backend decides which one leads.
    const cs = (pv && pv.claim_set) || null;
    const stale = !!(cs && (cs.superseded || cs.contradicted));
    // The marker names WHICH finding. Calling a contradicted record
    // "superseded" would be a one-word summary that is false: it is not old,
    // it is stamped with a claim set its own text disagrees with, and the two
    // call for different actions (read it as history vs. re-run it).
    const staleWord = cs && cs.contradicted ? 'contradicted' : 'superseded';

    // The claim-set banner goes ABOVE the index strip, i.e. first in the panel
    // and before any stored caveat. A reader who is going to be shown a
    // withdrawn sentence has to be told before they read it, not after. It adds
    // to the panel and removes nothing: every caveat below is still the run's
    // own words, complete and in order — see claimSetNotice() for why marking
    // rather than hiding is the only honest option here.
    if (stale) caveatsEl.appendChild(claimSetNotice(pv));

    // An index strip, not a summary — it counts what is below it and hides
    // none of it. This is the one line that answers "how qualified is this
    // run" before three hundred words have been read, which is exactly the
    // question a wall of prose was refusing to answer at a glance.
    caveatsEl.appendChild(h('div.cv-index', {
      title: 'Every caveat below is shown in full — none is collapsed, clamped or dismissible. '
        + '"Hard stop" counts the ones that say a number must not be quoted off this screen at '
        + 'all, as distinct from the ones that qualify how it should be read.'
        + (stale ? '\n\nThis run\'s caveats were written against a claim set that has since been '
          + 'superseded. They are shown unedited; the banner above says what has been withdrawn.' : ''),
    },
      h('span', h('b', String(cav.length + STANDING_LIMITS.length)), ' stated'),
      h('span.sep', '·'),
      h('span' + (hard ? '.hard' : ''), h('b', String(hard)), ' hard stop', hard === 1 ? '' : 's'),
      // One word. "claim set superseded" pushed this fixed-height, nowrap
      // strip 21 px past the rail and clipped "none hidden" — the one promise
      // in the panel that must not be the thing that falls off the edge.
      stale ? h('span.sep', '·') : null,
      stale ? h('span.stale', staleWord) : null,
      h('span.grow'),
      h('span.dim', 'none hidden')));

    // The run's own caveats get an era label when — and only when — there is
    // more than one era in play. It is what turns the block below from advice
    // into a record: same words, stated as what this run said.
    if (stale) caveatsEl.appendChild(claimSetEraLabel(pv));

    if (!cav.length) {
      caveatsEl.appendChild(h('div.caveat.info', h('span',
        run ? 'The backend returned no caveats for this run. That is itself worth checking — a run with zero stated caveats is unusual.'
          : 'Caveats appear here for the selected run.')));
    } else {
      cav.forEach((c, i) => caveatsEl.appendChild(caveatNode(String(c), i + 1)));
    }
    // Standing limits are true of the PROJECT, not of this run, so they cannot
    // arrive in RunResult.provenance.caveats — that list is per-run and may
    // legitimately be empty. They are rendered unconditionally, below the run's
    // own caveats, so there is no run whose panel omits them.
    caveatsEl.appendChild(h('div.hr-lbl', { style: { margin: 'var(--sp-6) 0 var(--sp-1)', padding: '0 var(--sp-5)' } },
      'standing limits · independent of this run'));
    STANDING_LIMITS.forEach((s, i) => caveatsEl.appendChild(caveatNode(s, i + 1, 'standing')));

    // The section header gets the severity count only. `--h-head` is one fixed
    // 24 px row (tokens.css 3.3) and it already carries a title and a hint —
    // anything longer than two tokens wraps out of its own box.
    clear(tallyEl);
    if (hard) tallyEl.appendChild(h('span.hard', `${hard} hard`));
    // "This run's caveats are historical" is exactly the kind of thing the
    // header exists to say before the panel is read, so it says it here too.
    if (stale) tallyEl.appendChild(h('span.stale', staleWord));

    /* reconstruction report */
    const rep = (run && run.report) || null;
    clear(reportEl);
    if (!rep) {
      fill(reportEl, why(!run
        ? 'No run is loaded, so no specimen has been reconstructed. This block '
          + 'describes the grid and the tissue a run was built on.'
        : run.status !== 'done'
          ? `Run ${run.run_id} is ${run.status}. The reconstruction report is `
            + 'written when the job finishes.'
          : `Run ${run.run_id} finished without a reconstruction report, so the `
            + 'grid, the cell counts and the tissue means it was solved on are '
            + 'unstated. That is a backend omission, not an empty specimen.'));
    } else {
      /**
       * @param {string} k        the report key, for reportProvenance()
       * @param {string} label    the row label
       * @param {*} value         the rendered value, or null when absent
       * @param {string} unitTxt  the unit — every number in this block has one
       * @param {string} absent   what to say INSTEAD of a dash when value is null
       */
      const kv = (k, label, value, unitTxt, absent) => {
        const p = reportProvenance(k, run);
        reportEl.appendChild(h('dt', label));
        reportEl.appendChild(h('dd',
          value === null || value === undefined ? nullValue() : value,
          (unitTxt && value !== null && value !== undefined) ? h('span.unit', unitTxt) : null,
          chip(p.kind, p.why)));
        if (value === null || value === undefined) {
          reportEl.appendChild(why(absent
            || `${label} is not carried in this run's reconstruction report. `
              + 'Read the dash as absent, not as zero.'));
        }
      };
      const fin = (v) => (Number.isFinite(v) ? v : null);
      const pair = (a) => (Array.isArray(a) && a.length >= 2 ? a : null);

      // Units and denominators, on every row. Before this pass the block read
      // "vasc 0.021 frac / ECM 0.310 frac / O2 0.640 norm / coverage 74.2 %" —
      // four different quantities sharing three unexplained abbreviations, none
      // of which said what the fraction was OF or over what the mean was taken.
      // real_pipeline.reconstruction_report()._mean averages over the TISSUE
      // MASK, not the whole grid; tissue_coverage_fraction is that mask's own
      // share of the grid nodes. Both are said here, once, under the rows.
      kv('n_cells', 'reconstructed', fin(rep.n_cells) === null ? null : int(rep.n_cells), 'cells');
      // The two-number problem, removed: "tumour 21,140" next to "cells 89,436"
      // made the share an exercise for the reader. `celltype_fractions` is the
      // backend's own division and is printed when the run carries it; nothing
      // is divided here.
      const tumorFrac = rep.celltype_fractions
        && Number.isFinite(rep.celltype_fractions.TUMOR)
        ? rep.celltype_fractions.TUMOR : null;
      kv('n_tumor', 'of which tumour', fin(rep.n_tumor) === null ? null
        : h('span', int(rep.n_tumor), h('span.unit', 'cells'),
          tumorFrac === null ? null
            : h('span', { style: { marginLeft: 'var(--sp-4)' } },
              pct(tumorFrac, 1), h('span.unit', '%'))),
      null,
      'This reconstruction does not state a tumour-cell count, so every '
      + 'per-tumour-cell metric above has no denominator on this screen.');
      kv('grid', 'grid', pair(rep.grid) ? `${int(rep.grid[0])} × ${int(rep.grid[1])}` : null, 'nodes');
      kv('dx_um', 'dx', fin(rep.dx_um) === null ? null : num(rep.dx_um, 1), 'um/node');
      kv('extent_um', 'extent', pair(rep.extent_um)
        ? `${num(rep.extent_um[0], 0)} × ${num(rep.extent_um[1], 0)}` : null, 'um');
      kv('mean_vascular_density', 'vasc', fin(rep.mean_vascular_density) === null
        ? null : auto(rep.mean_vascular_density), 'vol frac');
      kv('mean_ecm_density', 'ECM', fin(rep.mean_ecm_density) === null
        ? null : auto(rep.mean_ecm_density), 'vol frac');
      kv('mean_oxygen', 'O2', fin(rep.mean_oxygen) === null
        ? null : auto(rep.mean_oxygen), 'norm. pO2');
      kv('tissue_coverage_fraction', 'coverage', fin(rep.tissue_coverage_fraction) === null
        ? null : pct(rep.tissue_coverage_fraction, 1), '% of grid nodes');
      reportEl.appendChild(why(
        'vasc, ECM and O2 are means over the TISSUE MASK — the nodes coverage '
        + 'counts — not over the whole grid. Volume fractions and normalised '
        + 'pO2 both run 0 to 1; pO2 is normalised to the arterial supply, not '
        + 'to a measured tension on this tissue.'));
      if (rep.oxygen) {
        reportEl.appendChild(h('dt', 'O2 mode'));
        reportEl.appendChild(h('dd', rep.oxygen.mode || nullValue(),
          chip(rep.oxygen.is_surrogate ? 'assumed' : 'computed', rep.oxygen.caveat || rep.oxygen.source || '')));
        if (!rep.oxygen.mode) {
          reportEl.appendChild(why('This run does not name the route that '
            + 'produced its oxygen field, so the chip beside it is the only '
            + 'thing standing between that field and being read as measured.'));
        }
      }
      if (rep.cell_state) {
        const cs = rep.cell_state;
        const unmeasured = Array.isArray(cs.unmeasured_states) ? cs.unmeasured_states : [];
        // is_measured is the pipeline's own verdict; an UNMEASURED state on the
        // panel demotes the chip regardless, because the route then assigns that
        // state to NO cell and a reader must not read that as a measured zero.
        const kind = unmeasured.length ? 'assumed' : (cs.is_measured ? 'measured' : 'assumed');
        reportEl.appendChild(h('dt', 'state mode'));
        reportEl.appendChild(h('dd', cs.mode || nullValue(),
          chip(kind, cs.caveat || cs.source || '')));
        if (!cs.mode) {
          reportEl.appendChild(why('This run does not name the route that '
            + 'labelled cells.state, and that label is what the resistance '
            + 'multipliers act on. Read the dash as unstated.'));
        }
        if (unmeasured.length) {
          reportEl.appendChild(h('dt', 'UNMEASURED'));
          reportEl.appendChild(h('dd', unmeasured.join(', '),
            chip('assumed', 'Declared UNMEASURED by the panel-coverage check.')));
          // Was a --fs-micro `dd` with the whole explanation on a hover, in a
          // cell that is nowrap-with-ellipsis. This is the single most
          // misreadable row in the block: a state that is assigned to NO cell
          // looks identical to a state that was computed and came out zero.
          reportEl.appendChild(why(
            'These states are assigned to NO cell because this panel cannot '
            + 'measure them. That is NOT the same as computed and found to be '
            + 'zero, and their resistance multipliers therefore act on nothing.'));
        }
        const fr = cs.fractions_tumor;
        if (fr) {
          // "tumour V/Q/H  12.0 / 74.1 / 13.9 %" asked the reader to expand
          // three initials and then to remember, from a row eleven above, what
          // the percentages were of. The words are written out and the
          // denominator is named on the same screen as the numbers.
          reportEl.appendChild(h('dt', 'tumour bands'));
          reportEl.appendChild(h('dd',
            `${pct(fr.viable, 1)} / ${pct(fr.quiescent, 1)} / ${pct(fr.hypoxic, 1)}`,
            h('span.unit', '%'),
            chip(kind, 'Tumour-cell band occupancy under this state route. These occupancies are what the PDConfig resistance multipliers act on.')));
          reportEl.appendChild(why(
            'viable / quiescent / hypoxic, in that order, as percentages of '
            + (fin(rep.n_tumor) === null
              ? 'the tumour cells in this reconstruction (whose count this run '
                + 'does not state)'
              : `the ${int(rep.n_tumor)} tumour cells above`)
            + ' — not of all cells.'));
        }
      }
      if (rep.sensitivity) {
        const s = rep.sensitivity;
        reportEl.appendChild(h('dt', 'sens mode'));
        reportEl.appendChild(h('dd', s.mode || nullValue(),
          chip(s.is_uniform ? 'assumed' : 'computed', s.caveat || '')));
        if (!s.mode) {
          reportEl.appendChild(why('This run does not name the per-cell '
            + 'sensitivity model it applied. Every kill number above is '
            + 'conditional on it.'));
        }
        if (s.mgmt_applied) {
          reportEl.appendChild(h('dt', 'MGMT scalar'));
          // Dimensionless, and specifically a MULTIPLIER on the per-cell
          // sensitivity array — a bare "0.62" beside a column of unit-carrying
          // values read as a fifth measured quantity.
          reportEl.appendChild(h('dd',
            Number.isFinite(s.mgmt_scalar) ? auto(s.mgmt_scalar) : nullValue(),
            Number.isFinite(s.mgmt_scalar) ? h('span.unit', '× sensitivity') : null,
            chip('computed', 'Derived from the supplied per-patient MGMT log2 CPM through the calibrated logistic.')));
          if (!Number.isFinite(s.mgmt_scalar)) {
            reportEl.appendChild(why('This run reports that an MGMT scalar was '
              + 'applied but does not carry its value, so what it multiplied '
              + 'the sensitivity by is unknown from this screen.'));
          }
        }
      }
    }

    // The selected cell's explanation is made of THIS run's label descriptions,
    // THIS run's band thresholds and THIS run's reconstructor, so it is repainted
    // whenever the run under it changes — otherwise a pinned cell keeps the
    // previous run's reasons under the new run's picture.
    renderCell();
  }

  /* ------------------------------------------------------------ root -- */

  const na = (what) => h('div.dim', { style: { padding: 'var(--sp-4) var(--sp-5)' } },
    `${what} is unavailable without WebGL2.`);

  /* THE RETURN LEG OF THE ROUND TRIP.
   *
   * Clicking a survival row on THE FINDING isolates that mechanism in the
   * specimen, and main.js routes it through setColourMode, which calls
   * revealMechanism() — so the rail SCROLLS from the top of RESULT down to
   * this section. That is the right destination (this is the key for the
   * picture the click just produced) and it is also several screens from what
   * the reader clicked, with nothing on screen offering the way back.
   *
   * One control, in the section the trip lands in, pointing at the section it
   * started from. It is a real <button>, so it is in the tab order and takes
   * the global focus ring; it is always present rather than appearing only
   * after an isolate, because a navigation aid that exists only in one state
   * is one a reader cannot learn. The isolate itself is undone by "show all
   * categories" in the legend below — see js/ui/mechanism.js. */
  const mechSect = h('details.sect', { open: true, id: 'sect-mech' },
    h('summary', sectTitle('Mechanism'), h('span.hint', 'reconciles with kill_fraction')),
    h('div.sect-body',
      h('div.frow', { style: { marginBottom: 'var(--sp-3)' } },
        h('span.dim', {
          style: { fontSize: 'var(--fs-disclosure)' },
          title: 'This breakdown decomposes the claim stated at the top of '
            + 'this rail. Clicking a survival row there isolates one of these '
            + 'categories in the specimen and brings you here.',
        }, 'decomposes the claim under'),
        h('button.btn.ghost.sm', {
          type: 'button',
          title: 'Scroll back up to The finding, the panel these categories '
            + 'decompose. Nothing about the picture changes.',
          onclick: () => revealFinding(),
        }, 'The finding'),
        h('span.grow')),
      mechanismEl || na('Mechanism')));

  const railBody = h('div.rail-body',
    /* ---- RESULT: what did the run find, and how much of it is measured -- */
    railGroup('result', 'Result', 'RunResult',
      // THE FINDING LEADS THE RAIL. Everything below it is supporting detail,
      // and until 2026-08-22 the rail opened on the supporting detail instead:
      // four absolute metrics, led by the one number the project forbids
      // quoting. This section states what the run found, in counts of cells,
      // in words nobody has to learn first — and names the measurement that
      // would refute it, which is the whole reason this claim replaced a
      // compound ranking. js/ui/finding.js carries the full argument.
      h('details.sect', { open: true, id: 'sect-finding' },
        h('summary', sectTitle('The finding'), h('span.hint', 'what this run says, in cells')),
        h('div.sect-body', findingPanel.el)),

      // SECOND, because it is the one panel in this rail that answers a question
      // about ONE cell rather than about the run, and the gesture that fills it
      // — clicking the specimen — is the first thing anybody does. It sat
      // nowhere at all: the stage's hover card names a cell's phenotype and
      // coordinates and stops there, so "why is this cell alive" had no answer
      // anywhere in the product. See the block comment on renderCell().
      h('details.sect', { open: true, id: 'sect-cell' },
        h('summary', sectTitle('Selected cell'),
          h('span.hint', 'click a cell · measured vs derived')),
        h('div.sect-body', cellEl)),

      h('details.sect', { open: true },
        h('summary', sectTitle('Metrics'), h('span.hint', 'provenance on every value')),
        h('div.sect-body', errorEl, metricsEl)),

      // REACH sits between Metrics and Mechanism because that is the order the
      // questions are asked in. Metrics says WHAT the run produced; Reach says
      // WHERE the drug got to and how much of the tumour it did not get to;
      // Mechanism says WHY each surviving cell survived. Reach is the section
      // that replaced the ranked drug panel as the product's answer, so it sits
      // high, is open by default, and is never collapsed to make room.
      h('details.sect', { open: true, id: 'sect-reach' },
        h('summary', sectTitle('Delivery or dose'), h('span.hint', 'the discrimination')),
        h('div.sect-body', reachPanel.el)),

      // Was a floating card over the stage's top-right quarter. It lives here
      // now, DIRECTLY under Metrics, because it is the decomposition of the
      // `Kill fraction` value printed immediately above it and reconciles with
      // it bit-exactly — the two belong on the same screen, in that order, and
      // neither belongs on top of the tissue. js/ui/mechanism.js carries the
      // full argument.
      mechSect,

      h('details.sect', { open: true },
        h('summary', sectTitle('Caveats'), tallyEl, h('span.hint', 'not collapsible for a reason')),
        h('div.sect-body', { style: { padding: '0' } }, caveatsEl)),

      h('details.sect', { open: true },
        h('summary', sectTitle('Provenance'), h('span.hint', 'RunResult.provenance')),
        h('div.sect-body', provEl)),

      h('details.sect', { open: true },
        h('summary', sectTitle('Reconstruction report'), h('span.hint', 'RunResult.report')),
        h('div.sect-body', reportEl))),

    /* ---- VIEW: how is the picture being drawn -------------------------- */
    railGroup('view', 'View', 'display only — never the run',
      h('details.sect', { open: true },
        h('summary', sectTitle('Channels'), h('span.hint', 'cell type')),
        h('div.sect-body', chanHead, chanList, chanDen)),

      h('details.sect', { open: true },
        h('summary', sectTitle('Field overlay'), h('span.hint', '/field/{name}')),
        h('div.sect-body',
          h('div.fgrid',
            h('label.lbl', { for: 'f-overlay' }, 'Field'), fieldSel,
            h('label.lbl', { for: 'f-ramp' }, 'Ramp'), rampSel,
            h('label.lbl', { for: 'f-opacity' }, 'Opacity'),
            h('div.with-unit', opacity, h('span.unit', opacityVal, '%')),
            h('label.lbl', { for: 'f-gamma' }, 'Gamma'),
            h('div.with-unit', gamma, h('span.unit', gammaVal)),
            h('span.lbl', 'Colour by'), colourSeg,
          ),
          fieldStats,
          fieldNote)),

      h('details.sect', { open: true },
        h('summary', sectTitle('Window / level'), h('span.hint', 'transfer function')),
        h('div.sect-body', windowLevelEl || na('Window / level'))),

      h('details.sect', { open: false },
        h('summary', sectTitle('Render'), h('span.hint', 'WebGL2')),
        h('div.sect-body',
          h('div.fgrid',
            h('label.lbl', { for: 'r-size' }, 'Point size'),
            h('div.with-unit', sizeBoost, h('span.unit', sizeVal, '×')),
            h('label.lbl', { for: 'r-alpha' }, 'Cell alpha'),
            h('div.with-unit', alphaR, h('span.unit', alphaVal, '%')),
            h('label.lbl', { for: 'r-round' }, 'Round'),
            h('label.chk', roundChk, h('span.lbl', 'circular sprites')),
            h('label.lbl', { for: 'r-lod' }, 'Adaptive LOD'),
            h('label.chk', lodChk, h('span.lbl', 'subsample while panning')),
          ),
          h('div.frow', { style: { marginTop: 'var(--sp-3)' } },
            // Same action as the toolbar `fit` and the F key, so: same word.
            h('button.btn.sm', {
              title: 'Fit section to view. Same action as the toolbar fit and the F key.',
              onclick: () => onResetView && onResetView(),
            }, 'fit'),
            h('button.btn.sm', {
              title: 'Render 180 frames over a synthetic 850k-cell section and report '
                + 'median fps, frame time and upload time. A diagnostic for this machine, '
                + 'not a property of any run.',
              onclick: () => onBenchmark && onBenchmark(benchOut),
            }, 'benchmark'),
          ),
          benchOut))),

    /* ---- MEASURE: what did I measure on it ----------------------------- */
    railGroup('measure', 'Measure', 'drawn on the specimen',
      h('details.sect', { open: true },
        h('summary', sectTitle('Line profile'), h('span.hint', 'L · drag')),
        h('div.sect-body.tool-body', profileEl || na('Line profile'))),

      h('details.sect', { open: true },
        h('summary', sectTitle('Regions'), h('span.hint', 'R · drag')),
        h('div.sect-body.tool-body', regionEl || na('Regions'))),

      // Distance had a toolbar button and a shortcut but no panel, no read-out
      // and no clear. Every tool with a stage gesture has a home here now.
      h('details.sect', { open: true },
        h('summary', sectTitle('Distance'), h('span.hint', 'M · drag')),
        h('div.sect-body.tool-body', distanceEl || na('Distance')))),
  );

  const nav = railNav([
    { id: 'result', label: 'Result' },
    { id: 'view', label: 'View' },
    { id: 'measure', label: 'Measure' },
  ], railBody);

  // The rail head names the run it is interpreting. It used to say only
  // "active run", which is true and useless the moment a second run is in
  // flight: the top bar names the run you just submitted, this rail shows the
  // one before it, and nothing on screen said which numbers belonged to which.
  const root = h('div.rail.rail-right',
    h('div.panel-head', h('h2', 'Read'), h('span.grow'),
      h('span.hint', { title: 'This rail interprets the active run. It never changes what will be run — that is the Setup rail on the left.' },
        'showing'),
      headRunEl),
    nav.el,
    railBody,
  );

  buildChannels(null, 0);
  renderRun(null);
  setFieldInfo(null, null);

  /* --------------------------------------------- selection watcher ------
   * The pin lives on the viewer (`viewer.selected`, set by the stage's own
   * pointer handler in js/main.js) and there is no seam through which that
   * module tells a panel about it — the read-out it does drive is an overlay
   * it owns. Rather than add a second pin owner that could disagree with the
   * first, this reads the one the stage already writes.
   *
   * The cost is two integer comparisons and three identity comparisons per
   * animation frame, and `renderCell()` runs only when one of them changes.
   * The labels and the field are watched alongside the index because both
   * arrive AFTER the geometry — the mechanism payload is a separate request,
   * and a time-lapse frame swaps the field without any panel being told — so a
   * cell pinned before either lands would otherwise keep saying they are
   * missing. None of the four is reassigned per frame (js/gl/viewer.js writes
   * them only in setGeometry / setField / setLabels), so a still stage costs
   * four comparisons and no work. Started only when the viewer can actually
   * pick, so the stub viewer (no WebGL2) never enters the loop at all. */
  if (typeof viewer.cellInfo === 'function'
      && typeof requestAnimationFrame === 'function') {
    let lastSel = -1;
    let lastLabels = viewer.labels || null;
    let lastGeom = viewer.geom || null;
    let lastFieldRef = viewer.field || null;
    const watch = () => {
      const sel = Number.isInteger(viewer.selected) ? viewer.selected : -1;
      const labels = viewer.labels || null;
      const geom = viewer.geom || null;
      const fieldRef = viewer.field || null;
      if (sel !== lastSel || labels !== lastLabels || geom !== lastGeom
          || fieldRef !== lastFieldRef) {
        lastSel = sel; lastLabels = labels; lastGeom = geom;
        lastFieldRef = fieldRef;
        renderCell();
      }
      requestAnimationFrame(watch);
    };
    requestAnimationFrame(watch);
  }

  return {
    el: root,
    buildChannels,
    renderRun,
    setFieldInfo,
    set onFieldSelect(fn) { onFieldChange = fn; },
    get fieldName() { return fieldSel.value || null; },
    selectField(name) { fieldSel.value = name || ''; },
    setRampFor(name) {
      const meta = fieldMeta(name || '');
      if (meta.ramp && RAMP_NAMES.includes(meta.ramp)) {
        rampSel.value = meta.ramp;
        viewer.setRamp(meta.ramp);
      }
    },
    /** Keep the rail checkbox in step with the toolbar's colour-mode control. */
    /**
     * Mirror the stage toolbar. `disabled` is the same array main.js computes
     * for #tb-colour, so the two controls cannot disagree about whether a mode
     * is even available.
     */
    setColorMode(m, disabled) {
      colourBtns.forEach((b, i) => {
        b.setAttribute('aria-pressed', String(i === m));
        b.disabled = !!(disabled && disabled[i]);
      });
    },
    get rampName() { return rampSel.value; },
    /**
     * Bring the mechanism breakdown into view.
     *
     * Called when the operator selects the `mech` colour source: that turns the
     * breakdown into the KEY for what is on the stage, and a key you have to
     * go looking for is not one.
     */
    revealMechanism() { revealSection(mechSect); },
    /**
     * Bring THE FINDING back into view — the other direction of the same trip.
     *
     * Exposed as well as used internally so a caller that sends the reader to
     * the mechanism breakdown can send them back by the same route, rather
     * than a second module inventing its own scroll arithmetic against a rail
     * it does not own.
     */
    revealFinding,
    benchOut,
  };

  /**
   * Scroll the rail so `sect`'s own summary sits below the sticky group
   * header, opening it first if the operator had closed it.
   *
   * ONE implementation, because there are now two callers going in opposite
   * directions and a rail that scrolls differently depending on which way you
   * are travelling is a rail you cannot learn.
   *
   * scrollTop arithmetic on bounding rects, not scrollIntoView, for the reason
   * spelled out in ui/rail.js: the nearest scrollable ancestor is the shell,
   * and the shell must never scroll. The sticky-header offset is not optional
   * — without it "Mechanism · reconciles with kill_fraction" is scrolled
   * underneath RESULT and the panel appears to have lost its title.
   */
  function revealSection(sect) {
    if (!sect) return;
    if (!sect.open) sect.open = true;
    if (!railBody.isConnected || !railBody.clientHeight) return;
    const group = sect.closest ? sect.closest('.rail-group') : null;
    const head = group ? group.querySelector('.rg-head') : null;
    railBody.scrollTop += sect.getBoundingClientRect().top
      - railBody.getBoundingClientRect().top
      - (head ? head.getBoundingClientRect().height : 0);
    nav.sync();
  }

  /* Looked up rather than held as a binding: `#sect-finding` is built inline
   * in railBody and hoisting it out would move the string `'The finding'`
   * above `const railBody =`, where tests/test_first_session.py reads the
   * section order from. The lookup is once per click. */
  function revealFinding() {
    revealSection(railBody.querySelector('#sect-finding'));
  }
}
