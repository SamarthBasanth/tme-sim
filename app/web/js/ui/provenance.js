// provenance.js — the measured / computed / assumed decision, in ONE place.
//
// The product claim is that a number never appears without saying where it
// came from. That is only credible if the classification is a single readable
// rule set rather than a chip hard-coded next to each render site.
//
//   measured  — comes out of an instrument or a published measured dataset.
//               Nothing in the model chain touched its value.
//   computed  — a deterministic function of measured inputs through the
//               model. Reproducible; only as good as its inputs.
//   assumed   — a value we chose. Includes calibration constants, defaults,
//               surrogates, and any number that is an upper bound because a
//               resistance term was set to unity rather than observed.
//
// The escalation rule matters more than the taxonomy: a computed quantity
// that depends on an assumed input is reported as ASSUMED, not computed.

import { h } from '../util/dom.js';

export const KINDS = ['measured', 'computed', 'assumed', 'unknown'];

export function chip(kind, why) {
  const k = KINDS.includes(kind) ? kind : 'unknown';
  return h('span.chip.' + k, { title: why || describe(k) }, k);
}

export function describe(kind) {
  switch (kind) {
    case 'measured': return 'Measured: an instrument or published measured dataset produced this value directly.';
    case 'computed': return 'Computed: a deterministic model function of measured inputs.';
    case 'assumed': return 'Assumed: a value we chose — default, calibration constant, surrogate, or an upper bound standing in for an unmeasured term.';
    default: return 'Provenance not stated by the backend for this quantity.';
  }
}

/** Metric display metadata: label, unit, formatter hint, base provenance. */
export const METRICS = {
  kill_fraction: {
    label: 'Kill fraction', unit: 'frac', fmt: 'pct', better: 'higher',
    base: 'computed',
    src: 'Hill kill on per-cell drug AUC (cells/pkpd)',
  },
  mean_tumor_auc: {
    label: 'Mean tumour AUC', unit: 'uM.s', fmt: 'sci', better: 'higher',
    base: 'computed',
    src: 'transport integral over tumour-labelled cells',
  },
  therapeutic_coverage_1x: {
    label: 'Coverage @ 1x IC50', unit: 'frac', fmt: 'pct', better: 'higher',
    base: 'computed',
    src: 'fraction of tumour cells whose exposure reaches 1x IC50',
  },
  penetration_length_um: {
    label: 'Penetration length', unit: 'um', fmt: 'num1', better: 'higher',
    base: 'computed',
    src: '1/e decay distance of drug AUC from the vessel wall',
  },
};

/**
 * Provenance for a metric on a specific run, applying the escalation rule.
 * @returns {{kind:string, why:string}}
 */
export function metricProvenance(key, run) {
  const m = METRICS[key];
  if (!m) return { kind: 'unknown', why: describe('unknown') };
  const spec = (run && run.spec) || {};
  const prov = (run && run.provenance) || {};
  const rep = (run && run.report) || {};
  const reasons = [];
  let kind = m.base;

  if (spec.sensitivity_mode === 'uniform' && (key === 'kill_fraction' || key === 'therapeutic_coverage_1x')) {
    kind = 'assumed';
    reasons.push('sensitivity_mode=uniform sets every cell to 1.0, so this is a no-resistance UPPER BOUND, not an expected response');
  }
  if (rep.oxygen && rep.oxygen.is_surrogate) {
    kind = 'assumed';
    reasons.push('the oxygen field is a transcriptional surrogate under oxygen_mode=markers, and this metric depends on it');
  }
  if (spec.vascular_mode === 'kde') {
    kind = 'assumed';
    reasons.push('vascular_mode=kde did not converge under grid refinement');
  }
  // The cell-state label gates the kill through the PDConfig resistance
  // multipliers and is the largest single term in any absolute kill number here.
  if (key === 'kill_fraction' && rep.cell_state && !rep.cell_state.is_measured) {
    kind = 'assumed';
    reasons.push(`cells.state is labelled by state_mode=${rep.cell_state.mode}, which is not a measurement; across the three state routes this kill moves by up to 0.3252 absolute (1.77x) with tissue, drug, transport and RNG held identical`);
  }
  if (prov.state_unmeasured && prov.state_unmeasured.length && key === 'kill_fraction') {
    kind = 'assumed';
    reasons.push(`${prov.state_unmeasured.join(', ')} is UNMEASURED on this panel and is assigned to NO cell — that is not "computed and found to be zero", and its resistance multiplier therefore acts on nothing`);
  }
  if (prov.ic50_source && !/measured/i.test(prov.ic50_source) && key !== 'penetration_length_um') {
    kind = 'assumed';
    reasons.push(`IC50 is not a measured value (${prov.ic50_source})`);
  }
  const why = [describe(kind), m.src, ...reasons].join('\n\n');
  return { kind, why };
}

/** Provenance of a reconstruction-report field. */
export function reportProvenance(key, run) {
  const rep = (run && run.report) || {};
  switch (key) {
    case 'n_cells':
    case 'n_tumor':
    case 'celltype_counts':
    case 'extent_um':
      return { kind: 'measured', why: 'Cell positions and phenotype calls come from the segmented spatial assay.' };
    case 'mean_oxygen':
      return rep.oxygen && rep.oxygen.is_surrogate
        ? { kind: 'assumed', why: rep.oxygen.caveat || 'Marker-derived oxygen surrogate.' }
        : { kind: 'computed', why: 'Steady-state reaction-diffusion solve on the reconstructed vasculature.' };
    case 'mean_vascular_density':
    case 'mean_ecm_density':
    case 'tissue_coverage_fraction':
      return { kind: 'computed', why: 'Kernel-density field derived from measured marker positions.' };
    case 'dx_um':
    case 'grid':
      return { kind: 'assumed', why: 'Grid pitch is an operator choice in the RunSpec, not a property of the tissue.' };
    default:
      return { kind: 'computed', why: describe('computed') };
  }
}
