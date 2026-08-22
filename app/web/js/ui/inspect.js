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
import { FIELD_NAMES, fieldMeta } from '../fields.js';
import { RAMP_NAMES, RAMP_META } from '../gl/ramps.js';
import { METRICS, metricProvenance, reportProvenance, fieldProvenance, chip, STANDING_LIMITS, claimSetNotice, claimSetEraLabel } from './provenance.js';
import { auto, num, int, pct, sci, dur, NULL_GLYPH } from '../util/fmt.js';
import { railGroup, railNav } from './rail.js';
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
  windowLevelEl, profileEl, regionEl, distanceEl, mechanismEl }) {
  /* ------------------------------------------------------- channels ---- */
  const chanList = h('div');
  const chanRows = [];

  function buildChannels(counts, total) {
    clear(chanList);
    chanRows.length = 0;
    CELL_TYPES.forEach((name, i) => {
      const cb = h('input', {
        type: 'checkbox', checked: true,
        'aria-label': `channel ${name}`,
        onchange: (e) => { viewer.setVisible(i, e.target.checked); syncRow(i); onChange && onChange(); },
      });
      const c = counts ? counts[i] : null;
      const row = h('label.chan', { title: `${name} — click to toggle, alt-click to solo` },
        cb,
        h('span.sw', { style: { background: `var(${CHAN_VARS[i]})` } }),
        h('span.nm', name),
        h('span.ct', c === null ? NULL_GLYPH : `${int(c)}  ${total ? pct(c / total, 1) : NULL_GLYPH}%`),
      );
      row.addEventListener('click', (e) => {
        if (!e.altKey) return;
        e.preventDefault();
        CELL_TYPES.forEach((_, j) => {
          viewer.setVisible(j, j === i);
          chanRows[j].querySelector('input').checked = j === i;
          syncRow(j);
        });
        onChange && onChange();
      });
      chanRows.push(row);
      chanList.appendChild(row);
    });
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

  const chanHead = h('div.frow', { style: { marginBottom: '3px' } },
    h('button.btn.ghost.sm', { onclick: () => setAll(true) }, 'all'),
    h('button.btn.ghost.sm', { onclick: () => setAll(false) }, 'none'),
    h('span.grow'),
    h('span.dim', { style: { fontSize: 'var(--fs-micro)' } }, 'alt-click = solo'),
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

  function setFieldInfo(name, f) {
    const meta = fieldMeta(name || '');
    if (!name || !f) {
      fill(fieldStats,
        h('dt', 'state'), h('dd.dim', 'no overlay'));
      fieldNote.textContent = '';
      return;
    }
    const fp = fieldProvenance(name);
    fill(fieldStats,
      h('dt', 'min'), h('dd', auto(f.min), h('span.unit', meta.unit)),
      h('dt', 'max'), h('dd', auto(f.max), h('span.unit', meta.unit)),
      h('dt', 'grid'), h('dd', `${f.rows} × ${f.cols}`),
      h('dt', 'dx'), h('dd', num(f.dx, 2), h('span.unit', 'um')),
      h('dt', 'provenance'), h('dd', chip(fp.kind, fp.why)),
    );
    fieldNote.textContent = `${meta.source}. ${meta.note}`;
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
  const findingPanel = createFindingPanel();
  const reachPanel = createReachPanel();
  const provEl = h('dl.kv');
  const caveatsEl = h('div.caveats');
  const tallyEl = h('span.cv-tally');
  const reportEl = h('dl.kv');
  const errorEl = h('div.hidden');
  const headRunEl = h('span.mono.strong', {
    style: { marginLeft: 'var(--sp-3)', flex: '0 0 auto' },
  }, NULL_GLYPH);

  function renderRun(run) {
    findingPanel.render(run || null);
    reachPanel.render(run && run.status === 'done' ? run : null);
    headRunEl.textContent = (run && run.run_id) || NULL_GLYPH;
    headRunEl.title = run
      ? `Every number in this rail was produced by ${run.run_id}.`
      : 'No run is loaded.';

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
          + 'process running Codon.')));
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
       * `penetration_length_um` is a fitted parameter of a gradient this repo
       * measures as ~14x too shallow (results/kyle2004_penetration.md sec. 6).
       * Promoting any of them would move the violation rather than fix it.
       * Four equal numbers is the honest layout for four quantities of which
       * none is the answer. */
      for (const key of Object.keys(METRICS)) {
        const v = run.metrics[key];
        const m = METRICS[key];
        const p = metricProvenance(key, run);
        const isNull = v === null || v === undefined;
        metricsEl.appendChild(h('div.metric',
          h('div.m-name', m.label),
          h('div.m-val' + (isNull ? '.null' : ''),
            fmtMetric(key, v),
            isNull ? null : h('span.unit', metricUnit(key))),
          h('div.m-meta',
            chip(p.kind, p.why),
            h('span.m-src', m.src)),
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
        metricsEl.appendChild(h('div.metric',
          h('div.m-name', key),
          h('div.m-val' + (v === null ? '.null' : ''), auto(v)),
          h('div.m-meta', chip('computed', 'Metric not described by the front-end table; shown verbatim from the backend.')),
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
      fill(provEl, h('dt', 'state'), h('dd.dim', NULL_GLYPH));
    } else {
      // IC50 is LITERATURE, not MEASURED. It is a real measurement — on GDSC
      // cell lines, not on this section. Wearing the same green chip as the
      // cell count implies both came off this specimen; only one did.
      const ic50Measured = pv.ic50_source && /measured/i.test(pv.ic50_source);
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
        h('dt', 'alpha'), h('dd', auto(pv.alpha),
          chip('assumed', pv.alpha_source
            ? `${pv.alpha_source}\n\nFitted once against a published penetration distance, then held fixed; it is not re-measured per section.`
            : 'Tortuosity calibration constant. Fitted once, then held fixed; it is not re-measured per section.')),
        h('dt', 'IC50'), h('dd', pv.ic50_uM === null || pv.ic50_uM === undefined ? NULL_GLYPH : auto(pv.ic50_uM), h('span.unit', 'uM'),
          chip(ic50Measured ? 'literature' : 'assumed',
            `${pv.ic50_source || 'no source given'}${pv.ic50_note ? '\n\n' + pv.ic50_note : ''}`
            + '\n\nMeasured in cell culture, NOT on this section. This tissue has no measured drug response.')),
        h('dt', 'runtime'), h('dd', dur(run.runtime_s), chip('measured', 'Wall-clock time the backend spent on this job.')),
      );
      if (pv.ic50_source) {
        provEl.appendChild(h('dt', 'IC50 src'));
        provEl.appendChild(h('dd', { title: pv.ic50_source, style: { fontSize: 'var(--fs-micro)' } }, pv.ic50_source));
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
      fill(reportEl, h('dt', 'state'), h('dd.dim', NULL_GLYPH));
    } else {
      const kv = (k, label, value, unitTxt) => {
        const p = reportProvenance(k, run);
        reportEl.appendChild(h('dt', label));
        reportEl.appendChild(h('dd', value, unitTxt ? h('span.unit', unitTxt) : null, chip(p.kind, p.why)));
      };
      kv('n_cells', 'cells', int(rep.n_cells));
      kv('n_tumor', 'tumour', int(rep.n_tumor));
      kv('grid', 'grid', Array.isArray(rep.grid) ? `${rep.grid[0]} × ${rep.grid[1]}` : NULL_GLYPH);
      kv('dx_um', 'dx', num(rep.dx_um, 1), 'um');
      kv('extent_um', 'extent', Array.isArray(rep.extent_um) ? `${num(rep.extent_um[0], 0)} × ${num(rep.extent_um[1], 0)}` : NULL_GLYPH, 'um');
      kv('mean_vascular_density', 'vasc', auto(rep.mean_vascular_density), 'frac');
      kv('mean_ecm_density', 'ECM', auto(rep.mean_ecm_density), 'frac');
      kv('mean_oxygen', 'O2', auto(rep.mean_oxygen), 'norm');
      kv('tissue_coverage_fraction', 'coverage', pct(rep.tissue_coverage_fraction, 1), '%');
      if (rep.oxygen) {
        reportEl.appendChild(h('dt', 'O2 mode'));
        reportEl.appendChild(h('dd', rep.oxygen.mode || NULL_GLYPH,
          chip(rep.oxygen.is_surrogate ? 'assumed' : 'computed', rep.oxygen.caveat || rep.oxygen.source || '')));
      }
      if (rep.cell_state) {
        const cs = rep.cell_state;
        const unmeasured = Array.isArray(cs.unmeasured_states) ? cs.unmeasured_states : [];
        // is_measured is the pipeline's own verdict; an UNMEASURED state on the
        // panel demotes the chip regardless, because the route then assigns that
        // state to NO cell and a reader must not read that as a measured zero.
        const kind = unmeasured.length ? 'assumed' : (cs.is_measured ? 'measured' : 'assumed');
        reportEl.appendChild(h('dt', 'state mode'));
        reportEl.appendChild(h('dd', cs.mode || NULL_GLYPH,
          chip(kind, cs.caveat || cs.source || '')));
        if (unmeasured.length) {
          reportEl.appendChild(h('dt', 'UNMEASURED'));
          reportEl.appendChild(h('dd', {
            title: 'These states are assigned to NO cell because this panel cannot measure them. That is NOT the same as computed and found to be zero.',
            style: { fontSize: 'var(--fs-micro)' },
          }, unmeasured.join(', '), chip('assumed', 'Declared UNMEASURED by the panel-coverage check.')));
        }
        const fr = cs.fractions_tumor;
        if (fr) {
          reportEl.appendChild(h('dt', 'tumour V/Q/H'));
          reportEl.appendChild(h('dd',
            `${pct(fr.viable, 1)} / ${pct(fr.quiescent, 1)} / ${pct(fr.hypoxic, 1)}`,
            h('span.unit', '%'),
            chip(kind, 'Tumour-cell band occupancy under this state route. These occupancies are what the PDConfig resistance multipliers act on.')));
        }
      }
      if (rep.sensitivity) {
        const s = rep.sensitivity;
        reportEl.appendChild(h('dt', 'sens mode'));
        reportEl.appendChild(h('dd', s.mode || NULL_GLYPH,
          chip(s.is_uniform ? 'assumed' : 'computed', s.caveat || '')));
        if (s.mgmt_applied) {
          reportEl.appendChild(h('dt', 'MGMT scalar'));
          reportEl.appendChild(h('dd', auto(s.mgmt_scalar), chip('computed', 'Derived from the supplied per-patient MGMT log2 CPM through the calibrated logistic.')));
        }
      }
    }
  }

  /* ------------------------------------------------------------ root -- */

  const na = (what) => h('div.dim', { style: { padding: 'var(--sp-4) var(--sp-5)' } },
    `${what} is unavailable without WebGL2.`);

  const mechSect = h('details.sect', { open: true, id: 'sect-mech' },
    h('summary', 'Mechanism', h('span.hint', 'reconciles with kill_fraction')),
    h('div.sect-body', mechanismEl || na('Mechanism')));

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
        h('summary', 'The finding', h('span.hint', 'what this run says, in cells')),
        h('div.sect-body', findingPanel.el)),

      h('details.sect', { open: true },
        h('summary', 'Metrics', h('span.hint', 'provenance on every value')),
        h('div.sect-body', errorEl, metricsEl)),

      // REACH sits between Metrics and Mechanism because that is the order the
      // questions are asked in. Metrics says WHAT the run produced; Reach says
      // WHERE the drug got to and how much of the tumour it did not get to;
      // Mechanism says WHY each surviving cell survived. Reach is the section
      // that replaced the ranked drug panel as the product's answer, so it sits
      // high, is open by default, and is never collapsed to make room.
      h('details.sect', { open: true, id: 'sect-reach' },
        h('summary', 'Delivery or dose', h('span.hint', 'the discrimination')),
        h('div.sect-body', reachPanel.el)),

      // Was a floating card over the stage's top-right quarter. It lives here
      // now, DIRECTLY under Metrics, because it is the decomposition of the
      // `Kill fraction` value printed immediately above it and reconciles with
      // it bit-exactly — the two belong on the same screen, in that order, and
      // neither belongs on top of the tissue. js/ui/mechanism.js carries the
      // full argument.
      mechSect,

      h('details.sect', { open: true },
        h('summary', 'Caveats', tallyEl, h('span.hint', 'not collapsible for a reason')),
        h('div.sect-body', { style: { padding: '0' } }, caveatsEl)),

      h('details.sect', { open: true },
        h('summary', 'Provenance', h('span.hint', 'RunResult.provenance')),
        h('div.sect-body', provEl)),

      h('details.sect', { open: true },
        h('summary', 'Reconstruction report', h('span.hint', 'RunResult.report')),
        h('div.sect-body', reportEl))),

    /* ---- VIEW: how is the picture being drawn -------------------------- */
    railGroup('view', 'View', 'display only — never the run',
      h('details.sect', { open: true },
        h('summary', 'Channels', h('span.hint', 'cell type')),
        h('div.sect-body', chanHead, chanList)),

      h('details.sect', { open: true },
        h('summary', 'Field overlay', h('span.hint', '/field/{name}')),
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
        h('summary', 'Window / level', h('span.hint', 'transfer function')),
        h('div.sect-body', windowLevelEl || na('Window / level'))),

      h('details.sect', { open: false },
        h('summary', 'Render', h('span.hint', 'WebGL2')),
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
        h('summary', 'Line profile', h('span.hint', 'L · drag')),
        h('div.sect-body.tool-body', profileEl || na('Line profile'))),

      h('details.sect', { open: true },
        h('summary', 'Regions', h('span.hint', 'R · drag')),
        h('div.sect-body.tool-body', regionEl || na('Regions'))),

      // Distance had a toolbar button and a shortcut but no panel, no read-out
      // and no clear. Every tool with a stage gesture has a home here now.
      h('details.sect', { open: true },
        h('summary', 'Distance', h('span.hint', 'M · drag')),
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
     * go looking for is not one. Opens the section if the operator closed it,
     * then scrolls the rail body — scrollTop arithmetic on bounding rects, not
     * scrollIntoView, for the reason spelled out in ui/rail.js: the nearest
     * scrollable ancestor is the shell, and the shell must never scroll.
     */
    revealMechanism() {
      if (!mechSect.open) mechSect.open = true;
      if (!railBody.isConnected || !railBody.clientHeight) return;
      // Land the section's own SUMMARY below the group's sticky header, not
      // under it. railNav gets away with a bare alignment because the thing it
      // scrolls to IS the sticky header; a section inside the group has one
      // sitting on top of it, and without this offset "Mechanism · reconciles
      // with kill_fraction" was scrolled underneath RESULT and the panel
      // appeared to have lost its title.
      const head = mechSect.closest('.rail-group')
        ? mechSect.closest('.rail-group').querySelector('.rg-head')
        : null;
      railBody.scrollTop += mechSect.getBoundingClientRect().top
        - railBody.getBoundingClientRect().top
        - (head ? head.getBoundingClientRect().height : 0);
      nav.sync();
    },
    benchOut,
  };
}
