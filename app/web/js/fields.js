// fields.js — units and provenance for the four overlay fields.
//
// The TMSF payload deliberately carries no unit string (it is a raw numeric
// grid), so the units live here, on the client, and are stated on every
// legend. If a unit is wrong it is wrong HERE — that is the point of keeping
// it in one table instead of sprinkling strings through the render code.

export const FIELDS = {
  vascular_density: {
    label: 'Vascular density',
    unit: 'vol frac',
    short: 'vasc',
    ramp: 'viridis',
    provenance: 'computed',
    source: 'endothelial marker positions -> discrete vessel reconstruction',
    note: 'Dimensionless volume fraction in [0,1]. Vessel radius is assumed, not imaged.',
  },
  ecm_density: {
    label: 'ECM density',
    unit: 'vol frac',
    short: 'ecm',
    ramp: 'cividis',
    provenance: 'computed',
    source: 'COL1A1 / FN1 / SPARC transcript KDE',
    note: 'Dimensionless. Drives the tortuosity reduction of effective diffusivity.',
  },
  oxygen: {
    label: 'Oxygen',
    unit: 'norm. pO2',
    short: 'O2',
    ramp: 'inferno',
    provenance: 'computed',
    source: 'steady-state reaction-diffusion solve (mode-dependent)',
    note: 'Normalised to arterial supply, [0,1]. Under oxygen_mode=markers this is a SURROGATE and the run caveats say so.',
  },
  drug_auc: {
    label: 'Drug AUC',
    unit: 'uM.s',
    short: 'AUC',
    ramp: 'magma',
    provenance: 'computed',
    source: 'convection-diffusion-reaction transport integral',
    note: 'Time-integrated interstitial concentration. Conditional on the PK profile in the RunSpec.',
  },
};

export const FIELD_NAMES = Object.keys(FIELDS);

/**
 * The time-lapse quantity. Deliberately NOT in FIELDS, because FIELD_NAMES is
 * the list of overlays `/api/runs/{id}/field/{name}` will serve and this one is
 * not among them — it arrives in the TMSL payload, one frame at a time. Putting
 * it in the table would add a field-selector entry that 404s.
 *
 * It is also a DIFFERENT QUANTITY from drug_auc: concentration is instantaneous
 * and can fall, AUC is its running integral and cannot. Sharing a legend label
 * between them would be the kind of small lie that ends up in a figure caption.
 */
export const TIMELAPSE_FIELD = {
  label: 'Drug concentration',
  unit: 'uM',
  short: 'C',
  ramp: 'magma',
  provenance: 'computed',
  source: 'transport solve snapshot at t (TMSL frame)',
  note: 'Instantaneous interstitial concentration, NOT the AUC integral. With a static plasma source it rises toward a quasi-steady state and does not clear; time_resolved=true is what lets it fall.',
};

export function fieldMeta(name) {
  if (name === 'drug_concentration') return TIMELAPSE_FIELD;
  return FIELDS[name] || {
    label: name, unit: '', short: name, ramp: 'viridis',
    provenance: 'unknown', source: 'not described by the contract', note: '',
  };
}
