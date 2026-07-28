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

export function fieldMeta(name) {
  return FIELDS[name] || {
    label: name, unit: '', short: name, ramp: 'viridis',
    provenance: 'unknown', source: 'not described by the contract', note: '',
  };
}
