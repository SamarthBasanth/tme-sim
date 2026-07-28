// mechanism.js — the vocabulary of the per-cell mechanism attribution.
//
// This file is metadata only: it must stay in lockstep with
// tmesim/sim/mechanism.py, whose LABELS tuple defines the wire codes. The
// backend also ships `label_names` and `label_descriptions` inside every TMSM
// payload, so the authoritative text always travels with the data; what lives
// here is the presentation layer — the CSS variable each code paints with, the
// short name a dense legend can fit, and which limiter each survival rolls up
// into.
//
// Append only. A code's numeric value is a wire value.

/** Wire codes, in the order tmesim.sim.mechanism.LABELS defines them. */
export const MECH = {
  KILLED: 0,
  SURVIVED_UNREACHED: 1,
  SURVIVED_UNDERDOSED: 2,
  SURVIVED_RESISTANT: 3,
  SURVIVED_STOCHASTIC: 4,
  UNATTRIBUTED: 5,
  ALREADY_DEAD: 6,
  NOT_TUMOR: 7,
};

/**
 * One entry per code, index === code.
 *  var      CSS custom property carrying the colour (see css/mechanism.css)
 *  short    what the legend row says
 *  limiter  which bucket of the delivery/dose/resistance/chance roll-up, or null
 *  alarm    true if a non-zero count is itself the finding
 */
export const MECH_META = [
  { name: 'KILLED', short: 'KILLED', var: '--mech-killed', limiter: null, alarm: false },
  { name: 'SURVIVED_UNREACHED', short: 'UNREACHED', var: '--mech-unreached', limiter: 'delivery', alarm: false },
  { name: 'SURVIVED_UNDERDOSED', short: 'UNDERDOSED', var: '--mech-underdosed', limiter: 'dose', alarm: false },
  { name: 'SURVIVED_RESISTANT', short: 'RESISTANT', var: '--mech-resistant', limiter: 'resistance', alarm: false },
  { name: 'SURVIVED_STOCHASTIC', short: 'STOCHASTIC', var: '--mech-stochastic', limiter: 'chance', alarm: false },
  { name: 'UNATTRIBUTED', short: 'UNATTRIBUTED', var: '--mech-unattributed', limiter: 'unattributed', alarm: true },
  { name: 'ALREADY_DEAD', short: 'ALREADY DEAD', var: '--mech-alreadydead', limiter: null, alarm: false },
  { name: 'NOT_TUMOR', short: 'non-tumour', var: '--mech-nottumor', limiter: null, alarm: false },
];

export const MECH_VARS = MECH_META.map((m) => m.var);
export const MECH_NAMES = MECH_META.map((m) => m.name);

/** The four survival categories, in legend order. */
export const SURVIVAL_CODES = [
  MECH.SURVIVED_UNREACHED,
  MECH.SURVIVED_UNDERDOSED,
  MECH.SURVIVED_RESISTANT,
  MECH.SURVIVED_STOCHASTIC,
  MECH.UNATTRIBUTED,
];

/** Human sentence for the limiter roll-up, used under the stacked bar. */
export const LIMITER_LABEL = {
  delivery: 'DELIVERY',
  dose: 'DOSE',
  resistance: 'RESISTANCE',
  chance: 'CHANCE',
  unattributed: 'UNATTRIBUTED',
};

/**
 * One sentence naming the dominant reason a drug failed on this specimen,
 * built ONLY from the counts the backend reported. Returns null when there is
 * nothing to say rather than inventing a story.
 */
export function dominantReason(aggregate) {
  if (!aggregate || !aggregate.survival_breakdown) return null;
  const b = aggregate.survival_breakdown;
  const entries = Object.entries(b).filter(([, v]) => v > 0);
  if (!entries.length) return null;
  entries.sort((a, c) => c[1] - a[1]);
  const [key, frac] = entries[0];
  return { limiter: key, label: LIMITER_LABEL[key] || key.toUpperCase(), fraction: frac };
}
