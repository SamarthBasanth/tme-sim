# Synthetic-era experiments (RETIRED)

These scripts belong to the **synthetic-validation phase** of tmesim — Stage 0
(the honest kill-test on a synthetic GBM cohort with a *known* generative ground
truth) and the drug-selection rigor sweep (many independent synthetic
ground-truth models, plus adversarial hunts for cases where spatial-transport
modeling wins or loses vs. a cheap ML baseline).

They are **retired from the pipeline** per the real-only build
(`docs/REAL_BUILD_PLAN.md`, hard rule #1: *no synthetic ground truth, no invented
patients*). They all depended on the deleted synthetic cohort generator
(`tmesim/data/synthetic.py`) either directly (`generate_synthetic_patient` /
`generate_synthetic_cohort`) or transitively via `selection_lib` /
`results/selection_cache.npz` (which was precomputed from a synthetic cohort).

## What lives here

- `dose_diagnostic.py`, `dose_regime_contrast.py`, `transport_limited.py` —
  Stage-0 dose/transport diagnostics on the synthetic cohort.
- `drug_selection.py`, `precompute_selection.py`, `selection_lib.py`,
  `plot_selection_rigor.py`, `selection_rigor_workflow.js` — the drug-selection
  rigor workflow and its shared scoring library / cache builder / plotter.
- `gt_*.py` — the distinct biologically-motivated synthetic ground-truth models
  (multiplicative, Liebig-minimum, threshold, BBB-gated, ECM-sieving,
  hypoxia-resistant, potency/delivery-dominant, saturating, subtype-preference,
  noisy-realistic, anti-delivery).
- `adv_*.py` — adversarial ground truths sought to make the sim win / lose vs. ML.
- `finalize_report.py` — assembled the Stage-0 markdown report from a synthetic run.
- `molecular_demo.py` — the molecular-layer headline demo, run on a single
  synthetic patient (`generate_synthetic_patient`).

## Why they are still worth keeping

Their **honest findings remain valid as conclusions about the *limits* of
synthetic validation** and are preserved in:

- `results/STAGE0_FINDINGS.md`
- `results/DRUG_SELECTION_RIGOR.md`

Those write-ups stand as the record of what synthetic validation could and could
not establish (notably its circularity risk when the sim and the ground-truth
generator share transport assumptions). The real pipeline supersedes them:
`tmesim.reconstruct.real_pipeline.reconstruct_real` +
`tmesim.sim.engine.simulate_response`, driven from the `reconstruct` / `simulate`
CLI commands, on real open spatial data only.
