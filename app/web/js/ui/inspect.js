// inspect.js — the READ rail. Everything that interprets the active run, and
// nothing that changes what will be run (that is the SETUP rail, ui/spec.js).
//
// The eleven sections are organised into three named groups with a sticky
// header each and a navigator at the top of the rail:
//
//   RESULT   metrics, caveats, provenance, reconstruction report
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
import { METRICS, metricProvenance, reportProvenance, fieldProvenance, chip, STANDING_LIMITS } from './provenance.js';
import { auto, num, int, pct, sci, dur, NULL_GLYPH } from '../util/fmt.js';
import { railGroup, railNav } from './rail.js';

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

export function createInspectPanel({ viewer, onChange, onBenchmark, onResetView,
  windowLevelEl, profileEl, regionEl, distanceEl }) {
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
  const provEl = h('dl.kv');
  const caveatsEl = h('div.caveats');
  const reportEl = h('dl.kv');
  const errorEl = h('div.hidden');

  function renderRun(run) {
    /* metrics */
    clear(metricsEl);
    if (!run) {
      metricsEl.appendChild(h('div.dim', { style: { fontSize: 'var(--fs-xs)', padding: '4px 0' } },
        'No run selected. Metrics appear once a run reaches status=done.'));
    } else if (run.status === 'failed') {
      metricsEl.appendChild(h('div.caveat.crit', h('span', run.error || 'Run failed; the backend returned no error message.')));
    } else if (!run.metrics) {
      metricsEl.appendChild(h('div.dim', { style: { fontSize: 'var(--fs-xs)', padding: '4px 0' } },
        `Run ${run.run_id} is ${run.status}; no metrics yet.`));
    } else {
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
      fill(provEl,
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
        provEl.appendChild(h('dd', {
          title: pv.sensitivity_note,
          style: { fontSize: 'var(--fs-micro)', whiteSpace: 'normal', textAlign: 'left', justifyContent: 'flex-start', lineHeight: '1.35', padding: '3px 0' },
        }, pv.sensitivity_note));
      }
    }

    /* caveats — always fully visible, never truncated, never behind a toggle */
    clear(caveatsEl);
    const cav = (pv && pv.caveats) || [];
    if (!cav.length) {
      caveatsEl.appendChild(h('div.caveat.info', h('span',
        run ? 'The backend returned no caveats for this run. That is itself worth checking — a run with zero stated caveats is unusual.'
          : 'Caveats appear here for the selected run.')));
    } else {
      for (const c of cav) caveatsEl.appendChild(h('div.caveat', h('span', String(c))));
    }
    // Standing limits are true of the PROJECT, not of this run, so they cannot
    // arrive in RunResult.provenance.caveats — that list is per-run and may
    // legitimately be empty. They are rendered unconditionally, below the run's
    // own caveats, so there is no run whose panel omits them.
    caveatsEl.appendChild(h('div.hr-lbl', { style: { margin: '6px 0 1px', padding: '0 var(--sp-5)' } },
      'standing limits · independent of this run'));
    for (const s of STANDING_LIMITS) {
      caveatsEl.appendChild(h('div.caveat.standing', h('span', s)));
    }

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

  const railBody = h('div.rail-body',
    /* ---- RESULT: what did the run find, and how much of it is measured -- */
    railGroup('result', 'Result', 'RunResult',
      h('details.sect', { open: true },
        h('summary', 'Metrics', h('span.hint', 'provenance on every value')),
        h('div.sect-body', errorEl, metricsEl)),

      h('details.sect', { open: true },
        h('summary', 'Caveats', h('span.hint', 'not collapsible for a reason')),
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

  const root = h('div.rail.rail-right',
    h('div.panel-head', h('h2', 'Read'), h('span.grow'),
      h('span.hint', { title: 'This rail interprets the active run. It never changes what will be run — that is the Setup rail on the left.' },
        'active run')),
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
    benchOut,
  };
}
