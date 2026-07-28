// templates.js — starting points for the protocol editor.
//
// A blank editor for a language nobody has written before is useless, so every
// template here is a REAL, runnable TPL program against a section this
// deployment actually carries (`xenium_319`), not placeholder prose. Each one
// is annotated with the thing it exists to teach: units are mandatory, the
// dose -> PK conversion is stated, a sweep is a cross product, `sensitivity
// uniform` is an upper bound.

export const TEMPLATES = [
  {
    id: 'single',
    label: 'Single drug · one run',
    note: 'The smallest complete program: define a drug, an administration, run it.',
    runs: 1,
    source: `# Single run — temozolomide on the 319-panel Xenium section.
# Every physical quantity carries a unit. \`dose 150\` is an error.

drug TMZ from library "temozolomide"
  ic50 from gdsc                  # measured; errors if TMZ is absent from the screen
end

protocol single_dose
  route oral
  dose 150 mg/m2
  every 24 h
  repeat 1
end

run TMZ with single_dose on xenium_319
  dx 12 um
  sensitivity uniform             # NO-RESISTANCE UPPER BOUND, not an expected response
  report coverage, kill, mean_auc, penetration_length
end

# An assertion is RECORDED, never enforced. State the hypothesis first.
expect coverage > 0.2
`,
  },

  {
    id: 'dose-sweep',
    label: 'Dose sweep · one axis',
    note: 'Eight doses, one axis. Produces a response curve against dose.',
    runs: 8,
    source: `# One-axis sweep: response against administered dose.
# 50..400 mg/m2 in steps of 50 -> 8 runs. The planner states the count and the
# cost before anything is launched.

drug TMZ from library "temozolomide"
  ic50 from gdsc
end

protocol adjuvant
  route oral
  dose 150 mg/m2                  # overridden by the sweep below
  every 24 h
  repeat 5
end

sweep dose from 50 to 400 mg/m2 step 50
sweep dx 12 um                    # a fixed override, NOT a sweep axis

run TMZ with adjuvant on xenium_319
  sensitivity uniform
  report coverage, kill, mean_auc, penetration_length
end

expect kill > 0.05
`,
  },

  {
    id: 'stupp',
    label: 'Stupp-like · multi-cycle',
    note: '5 days on, 23 days rest, six cycles — the adjuvant schedule shape.',
    runs: 1,
    source: `# Stupp-like adjuvant schedule: 5 daily doses, 23 days rest, six cycles.
# The executor converts this to PKConfig(plasma_peak, half_life_s, n_doses,
# dose_interval_s) and MUST state the conversion and its assumptions — body
# surface area, bioavailability, volume of distribution — in provenance.

drug TMZ from library "temozolomide"
  ic50 from gdsc
end

protocol stupp_adjuvant
  route oral
  dose 150 mg/m2
  every 24 h
  repeat 5
  then rest 23 d
  cycles 6
end

run TMZ with stupp_adjuvant on xenium_319
  dx 12 um
  sensitivity uniform
  report coverage, kill, mean_auc, penetration_length
end

expect kill > 0.05
`,
  },

  {
    id: 'compare',
    label: 'Two-drug comparison',
    note: 'Ranking two alkylators on one specimen. The ranking is the defensible output.',
    runs: 2,
    source: `# Two drugs, one schedule, one specimen.
# Cross-drug RANKING is the defensible output on this section; the absolute
# kill is not, because the panel cannot support cell state (see the caveats
# attached to every result).

drug TMZ from library "temozolomide"
  ic50 from gdsc
end

drug CCNU from library "lomustine"
  ic50 from gdsc
end

protocol adjuvant
  route oral
  dose 150 mg/m2
  every 24 h
  repeat 5
end

compare TMZ, CCNU with adjuvant on xenium_319
  dx 12 um
  sensitivity uniform
  rank by coverage
  report coverage, kill, mean_auc, penetration_length
end
`,
  },

  {
    id: 'smiles',
    label: 'Drug from SMILES',
    note: 'A compound the library does not know, compiled through the molecular layer.',
    runs: 1,
    source: `# A compound that is not in the built-in library: compile it from structure.
# Anything not given is DERIVED by the molecular layer and labelled \`computed\`;
# anything given here is labelled \`user-supplied\`. Upload a .smi/.mol/.sdf and
# the panel resolves the same fields with the same provenance before you run.

drug MYCPD from smiles "CN1C(=O)N2C=NC(=C2N=N1)C(=O)N"
  ic50 12 uM                      # explicit override -> user-supplied
  bbb_penetrance 0.31             # any DrugParams field may be overridden
end

protocol probe_iv
  route iv
  dose 100 mg/m2 over 24 h        # \`over\` = infusion duration, not a bolus
  repeat 1
end

run MYCPD with probe_iv on xenium_319
  dx 12 um
  sensitivity uniform
  report coverage, kill, mean_auc, penetration_length
end
`,
  },

  {
    id: 'grid',
    label: 'Dose × schedule · two axes',
    note: '7 doses × 4 intervals = 28 PDE solves. Rendered as a heatmap.',
    runs: 28,
    source: `# Two-axis sweep. Sweeps are a CROSS PRODUCT: 7 doses x 4 intervals = 28
# runs, each one a full PDE solve. The planner refuses more than max_runs (64)
# without an explicit \`allow runs N\` — widen the sweep and you will see it.

drug TMZ from library "temozolomide"
  ic50 from gdsc
end

protocol adjuvant
  route oral
  dose 150 mg/m2
  every 24 h
  repeat 5
end

sweep dose from 100 to 400 mg/m2 step 50
sweep every from 12 to 48 h in 4 steps
sweep dx 14 um

run TMZ with adjuvant on xenium_319
  sensitivity uniform
  report coverage, kill, mean_auc, penetration_length
end
`,
  },
];

export const templateById = (id) => TEMPLATES.find((t) => t.id === id) || null;
