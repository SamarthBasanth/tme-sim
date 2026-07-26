export const meta = {
  name: 'drug-selection-rigor',
  description: 'Rigor sweep: does spatial-transport modeling beat a cheap ML baseline at per-patient drug selection, across many INDEPENDENT ground-truth models?',
  phases: [
    { title: 'GroundTruths', detail: 'many agents each author a distinct biologically-motivated ground-truth model and score sim vs baselines' },
    { title: 'Adversarial', detail: 'agents hunt for plausible ground truths where the sim wins/loses' },
    { title: 'Synthesis', detail: 'map when transport helps; assess circularity; honest verdict' },
  ],
}

const ROOT = '/Users/samarthbasanth/tme-sim'
const VENV = `${ROOT}/.venv/bin/python`

// A spectrum of DISTINCT biological ground-truth motifs, deliberately spanning
// cases where the sim SHOULD help (delivery-limited) and SHOULD NOT (potency-dominant,
// adversarial). Coverage of both is what makes the test honest.
const MOTIFS = [
  { key: 'multiplicative',   brief: 'kill = sens * delivery(vasc,ecm,bbb,mw) * Hill(dose, ic50). Generic multiplicative delivery x potency.' },
  { key: 'liebig_min',       brief: "Liebig's law of the minimum: kill = sens * min(delivery_term, potency_term). Limited by whichever is worse." },
  { key: 'threshold',        brief: 'Sharp delivery threshold: a drug only works if its perivascular delivery exceeds a cutoff; below it, kill ~0 regardless of potency.' },
  { key: 'hypoxia_resist',   brief: 'Hypoxic tumors resist most drugs; delivery + a hypoxia penalty gate kill; some drugs less hypoxia-sensitive.' },
  { key: 'ecm_sieving',      brief: 'Size-dependent ECM sieving dominates: large drugs (high MW) are excluded from dense-ECM tumors almost independent of potency.' },
  { key: 'bbb_gated',        brief: 'BBB-gated: only BBB-penetrant drugs reach tumor unless vascularization is high (leaky); low-BBB potent drugs fail in well-differentiated tumors.' },
  { key: 'subtype_pref',     brief: 'Each GBM subtype has an intrinsic drug preference (a subtype x drug affinity matrix), MODULATED by delivery.' },
  { key: 'saturating',       brief: 'Saturating dose-response with tissue binding: delivery sets local AUC, kill saturates; strongly-bound drugs plateau early.' },
  { key: 'potency_dominant', brief: 'Potency dominates, delivery is a minor multiplier (~0.8-1.0). The most potent drug is almost always best. Sim SHOULD NOT help here (null/honesty case).' },
  { key: 'delivery_dominant',brief: 'Delivery dominates, potency compressed (ic50 spread flattened). Best drug is set by which reaches the tumor. Sim SHOULD help here.' },
  { key: 'noisy_realistic',  brief: 'Blend of delivery + potency + subtype + moderate noise; a messy realistic mixture with no single dominant factor.' },
  { key: 'anti_delivery',    brief: 'ADVERSARIAL: the true-best drug is the LEAST deliverable one (e.g. an efflux/priming quirk). The sim SHOULD FAIL here; report it honestly.' },
]

const METRIC_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['motif', 'description', 'acc_sim', 'acc_potency', 'acc_feature_ml', 'acc_random',
             'regret_sim', 'regret_potency', 'regret_feature_ml',
             'frac_best_not_most_potent', 'sim_vs_gt_rank_corr', 'sim_beats_ml', 'sim_beats_potency', 'notes'],
  properties: {
    motif: { type: 'string' },
    description: { type: 'string', description: 'the biological assumption + key params of the ground-truth model you built' },
    acc_sim: { type: 'number' }, acc_potency: { type: 'number' }, acc_feature_ml: { type: 'number' }, acc_random: { type: 'number' },
    regret_sim: { type: 'number' }, regret_potency: { type: 'number' }, regret_feature_ml: { type: 'number' },
    frac_best_not_most_potent: { type: 'number', description: 'how delivery-limited the GT is (0=potency decides, 1=best drug never the most potent)' },
    sim_vs_gt_rank_corr: { type: 'number', description: 'mean per-patient rank corr of sim scores vs this GT (circularity proxy)' },
    sim_beats_ml: { type: 'boolean' }, sim_beats_potency: { type: 'boolean' },
    notes: { type: 'string', description: 'honest caveats; did the harness run cleanly; anything degenerate' },
  },
}

function gtPrompt(motif, idx) {
  return `You are one agent in a rigorous study of whether a mechanistic tumor-microenvironment PDE simulator beats cheap baselines at PER-PATIENT DRUG SELECTION. Your job: independently author ONE biologically-motivated GROUND-TRUTH drug-response model, then score the fixed policies against it with the SHARED harness. Be honest — some motifs are cases where the sim SHOULD lose.

Working dir: ${ROOT} . Run python as: ${VENV} . Run commands FROM ${ROOT}.

The precomputed cache results/selection_cache.npz (load via numpy, allow_pickle=True) has, for a fixed 30-patient cohort x 5-drug panel:
  panel[D] (drug names), K_sim[P,D] (sim predicted kill -- DO NOT use inside your ground truth), potency[P,D],
  feat[P,F], feat_names[F], vasc[P], ecm[P], sens[P], drug_mw[D], drug_ic50[D], drug_bbb[D], dose (scalar), alpha.
Observable patient properties you MAY use to build the ground truth: vasc, ecm, sens, and the features in feat (vasc, ecm, sens, hypoxic, subtype one-hots). Drug properties: drug_mw, drug_ic50, drug_bbb. You must NOT use K_sim to build the ground truth (that would be circular).

The shared harness experiments/selection_lib.py exposes:
  score_policies(gt_kill: np.ndarray[P,D], seed=0) -> dict with acc{sim,potency,feature_ml,random,oracle}, regret{...}, diag{frac_best_not_most_potent, sim_vs_gt_rank_corr, best_distribution}, sim_beats_ml, sim_beats_potency.

YOUR MOTIF (#${idx}): "${motif.key}" -- ${motif.brief}
Design a DISTINCT ground-truth kill matrix gt_kill[P,D] in [0,1] that faithfully encodes THIS biological assumption (do not just copy the multiplicative form unless that IS your motif). Use your own reasoning about GBM drug delivery/response. Add mild noise with a fixed seed. Build a delivery term from vasc/ecm/drug_bbb/drug_mw as your motif dictates, and a potency term from drug_ic50/dose/sens.

STEPS:
1. Write experiments/gt_${motif.key}.py that: loads the cache, builds gt_kill[P,D] for your motif (deterministic, seeded), calls selection_lib.score_policies(gt_kill), and prints it via selection_lib.print_json(result). Put the file in experiments/ and run it as \`${VENV} experiments/gt_${motif.key}.py\` from ${ROOT} (so \`import selection_lib\` resolves and the cache path is correct).
2. Read the RESULT_JSON line it prints.
3. Return the metrics in the required schema, flattening from the RESULT_JSON (acc_sim=result.acc.sim, etc.; frac_best_not_most_potent and sim_vs_gt_rank_corr from result.diag). In "description" state your biological assumption + key parameters; in "notes" give honest caveats (e.g. if your motif made one drug trivially dominant, say so).

Be rigorous and honest. If your motif is a case where the sim should lose (potency_dominant, anti_delivery), REPORT the loss plainly — that is the point.`
}

phase('GroundTruths')
const gtResults = await parallel(MOTIFS.map((m, i) => () =>
  agent(gtPrompt(m, i), { schema: METRIC_SCHEMA, label: `gt:${m.key}`, phase: 'GroundTruths' })
))
const clean = gtResults.filter(Boolean)
log(`ground-truth models scored: ${clean.length}/${MOTIFS.length}`)

phase('Adversarial')
const advPrompt = (goal, i) => `You are an adversarial agent probing the drug-selection study at ${ROOT} (run python: ${VENV} from ${ROOT}). Cache + harness as described: results/selection_cache.npz and experiments/selection_lib.py (score_policies(gt_kill)). You must NOT use K_sim inside a ground truth.

GOAL: ${goal}
Construct a BIOLOGICALLY PLAUSIBLE ground-truth kill matrix gt_kill[P,D] (from vasc/ecm/sens/hypoxic/subtype and drug mw/ic50/bbb) that achieves this goal, iterating on parameters. Write experiments/adv_${i}.py, run it, read RESULT_JSON. Return schema metrics. In "description" explain the mechanism and whether it is biologically defensible or contrived; in "notes" state clearly whether you SUCCEEDED at the goal and what it implies about the sim's robustness.`

const adv = await parallel([
  () => agent(advPrompt('Make the mechanistic SIM WIN as decisively as possible over BOTH the potency shortcut AND the feature-ML baseline, while staying biologically plausible (delivery-limited regime).', 'sim_wins'), { schema: METRIC_SCHEMA, label: 'adv:sim_wins', phase: 'Adversarial' }),
  () => agent(advPrompt('Make the mechanistic SIM LOSE to the feature-ML baseline using a plausible ground truth (find a blind spot of the sim -- e.g. response driven by an observable patient feature the sim ignores).', 'sim_loses_ml'), { schema: METRIC_SCHEMA, label: 'adv:sim_loses_ml', phase: 'Adversarial' }),
  () => agent(advPrompt('Make the potency shortcut WIN (a regime where delivery barely matters), to map the boundary where the sim stops adding value.', 'potency_wins'), { schema: METRIC_SCHEMA, label: 'adv:potency_wins', phase: 'Adversarial' }),
]).then(r => r.filter(Boolean))

phase('Synthesis')
const all = [...clean, ...adv]
const synthesis = await agent(
  `You are synthesizing a rigorous drug-selection study for the tmesim glioblastoma simulator. Below are metrics from ${all.length} INDEPENDENTLY-authored ground-truth models (12 biological motifs + 3 adversarial), each scoring three policies at per-patient drug selection: the mechanistic PDE sim, a potency-only shortcut, and a cross-validated feature-ML baseline (top-1 accuracy + regret; higher acc / lower regret is better).

DATA (JSON):
${JSON.stringify(all, null, 1)}

Write a rigorous, HONEST synthesis to ${ROOT}/results/DRUG_SELECTION_RIGOR.md covering:
1. Headline: across how many of the ${all.length} independent ground truths does the sim beat (a) the potency shortcut and (b) the feature-ML baseline? Give the counts and the mean/median accuracy gap.
2. THE CIRCULARITY TEST: does the sim's advantage correlate with 'frac_best_not_most_potent' (how delivery-limited the GT is) and with 'sim_vs_gt_rank_corr' (how much the GT happens to agree with the sim's own physics)? If the sim only wins when sim_vs_gt_rank_corr is high, that is CIRCULAR and you must say so. If it wins even at moderate correlation in delivery-limited regimes, that is more genuine. Report the pattern with numbers.
3. Where the sim LOSES (potency_dominant, anti_delivery, any feature-ML wins) and why -- its blind spots.
4. An honest bottom line: under what conditions does modeling spatial transport add value to drug selection, is that answerable synthetically at all, and what real data would settle it. Do NOT oversell; a "only helps in delivery-limited regimes, and whether real GBM is there needs clinical data" conclusion is the expected honest answer.
Also return a 250-word summary of the same, with the key counts.`,
  { label: 'synthesis', phase: 'Synthesis' })

return { n_ground_truths: all.length, synthesis_head: synthesis.slice(0, 1200), all_metrics: all }
