# Rigor Audit — what the verification pass found, and what still stands

This session put the real-data build through an adversarial verification pass. It found
**eleven distinct defects** (nine at first writing, two more since), several of which were mutually cancelling, and produced the
project's **first validation against measured drug concentrations in human tissue**.

Every headline number reported before this audit was wrong in at least one way. The
engine and the reconstruction were sound; the *claims* were not.

---

## 1. The first real external validation — human glioma drug exposure

Scored against **125 human glioma tissue samples across 42 patients** plus 27
PMID-anchored regional PK measurements (`results/glioma_regional_pk.md`).

| Axis | Result |
|---|---|
| **Partition magnitude** | 5/6 drugs within 10x, 4/6 within 3x. Aggregate bias 2.02x — **but 1.01x (essentially unbiased) excluding one outlier** |
| **Cross-drug ordering** | Measured penetrant-vs-excluded separation 4.66x; **model predicts 3.65x from structure alone** |
| **Regional ordering** | **FAIL.** Measured CE-vs-non-enhancing 1.45x-4.26x (p<0.0001, within-patient permutation); model predicts exactly 1.00x |

The regional failure has a precise cause — there is no barrier-integrity term anywhere in
the model — and the validation **derived the fix from the data**: a multiplicative barrier
factor cannot work (3.2x spread across drugs), but an **additive leak is drug-independent
to within 1.3x** (L = 0.211 / 0.230 / 0.272 across three independent studies), giving
`source = vascular_density * (bbb_penetrance + L * bbb_disruption) * plasma`, L ~ 0.24.
That is the signature of size-nonselective paracellular opening. **NOW IMPLEMENTED** in
`tmesim/reconstruct/barrier.py` (L=0.238): regional ratios go from a flat 1.00x to 1.19x and
5.87x against measured 1.45x and 4.26x; anchor error 2.49x -> 1.30x. Unfitted necrosis
ratios land at 1.02x/0.98x.

It also traced a 66x outlier (eribulin) to a genuine code defect: `bbb_penetrance` takes
the **mean** of four desirability scores, so a molecule with two disqualifying properties
is rescued by two favourable ones. *A mean cannot express a veto.* **NOW FIXED** (geometric
mean): eribulin 65.6x -> 6.9x, within-10x 5/6 -> 6/6, GM fold error 3.66x -> 2.00x.

---

## 2. Defects found and fixed

| # | Defect | Status |
|---|---|---|
| 1 | Endothelial fraction 2x too high (single-marker argmax bias): 20.6% -> 9.7% | fixed |
| 2 | `penetration_index` not grid-converged **and definitionally unable to discriminate drugs** | retired, replaced |
| 3 | Hypoxic fraction a grid artifact (94% at dx=30 -> 13% at dx=12 -> 7% at dx=8, not converged) | retracted |
| 4 | alpha off by exactly 1/ln^2(2) — half-distance vs e-fold convention (309.6 -> **148.6**) | fixed + guarded by test |
| 5 | Per-cell sensitivity a placeholder suppressing all kills ~5x | fixed (drug-specific, real markers) |
| 6 | Grid-dependent oxygen sink (raw count, scales as dx^2) | fixed |
| 7 | Grid-dependent tissue mask — **was partially cancelling #6** | **FIXED** (geometric 30um-disc mask; coverage spread 3.64x -> 1.007x; legacy mask was corrupting production kill 2.25x) |
| 8 | Unfair cross-drug comparison: TMZ advantaged because its resistance determinant is *unmeasured* | **FIXED** (uniform/determinant modes; plus a measured per-patient MGMT scalar that returns TMZ to #5) |
| 9 | `bbb_penetrance` mean cannot express a veto | **FIXED** (geometric mean; eribulin 65.6x -> 6.9x; within-10x 5/6 -> 6/6) |
| 10 | `membrane_permeability` returned 4.2e7 um/s for imatinib (tens of m/s) — unbounded logP | **FIXED** (logP clamp + unstirred aqueous boundary layer in series) |
| 11 | ~~Oxygen PDE **falsified as value-add** — beaten by distance-to-nearest-vessel (rho -0.150 vs -0.113), and by its own undiffused source (-0.127)~~ **WITHDRAWN — this "defect" was itself a defect: an inverted inequality** | **DEMOTION REVERSED** (`results/oxygen_demotion.md`, now the record of the correction). Every predictor is sign-aligned to oxygen, so the stronger one has the larger \|ρ\| — and it is the solve: re-derived on 816,769 cells at dx = 12 um, **-0.1523** against the distance baseline's **-0.1204** and its own undiffused source's **-0.1111**, i.e. the solve wins by **1.26x**. The quoted numbers above were backwards by their own magnitudes. The `vessel_distance` surrogate and the measured-marker route remain available via `reconstruct_real(oxygen_mode=...)`, and the default stays `'pde'` — now **on** the evidence rather than despite it. Deliberately modest: \|ρ\| ≈ 0.15 is weak, within-section, carried by VEGFA alone (HIF1A +0.0140), so this is *the best of four weak predictors*, not a validated oxygen field. Still true and untouched by the reversal: `OxygenConfig.consumption` is unfitted, so the solve's diffusion length is 2-5x the published ~150 um — a Spearman rho is invariant to monotone rescaling and cannot see that |

**#7 is the most instructive.** Two grid-dependence bugs pulled in opposite directions, so
the output looked *more* stable than it was. Compensating errors that produce plausible
results are invisible to a passing test suite — 55 tests were green throughout. Only
writing the math down (`docs/MATHEMATICAL_MODEL.md`) and checking invariants found them.

---

## 3. CORRECTED — the "vasculature blur" diagnosis was largely WRONG

**This section is superseded. Read `results/vessel_network_impact.md`.** The original
diagnosis below attributed three symptoms to the smooth endothelial KDE. Discrete
vessel segmentation tested it directly and refuted two of the three:

1. **The flat exposure floor is NOT a KDE artifact — it is real diffusion physics.**
   With a discrete source the field is exactly zero between vessels, yet far-field
   exposure moves only 1.12x (25.99% -> 23.18%). Reaction-diffusion length is 73-296 um
   and diffusion length ~4,000 um against a ~40 um intervessel gap, so drug genuinely
   reaches everywhere in a densely vascularised section. No tumour cell in any arm
   receives zero drug. The floor is a property of the tissue, not of our smoothing.

2. **"Penetration length cannot discriminate drugs" was STALE, and the cause was
   defect #7, not the vasculature.** Under the current geometric mask the fitted lengths
   are 178/109/58 um — a 3.07x spread. The old 26-28 um collapse was reproduced without
   running a solve by refitting the source null under the RETIRED occupancy mask: it
   returns 23.68 um, matching the published 23.7 exactly. Fixing the tissue mask cured it.

3. **The perivascular circularity (`niche_response.md`) DOES stand** and is eliminated by
   the discrete field.

A dose-matched control (absent from the original brief) was decisive: the KDE arm delivers
3.36x more total drug, so the naive KDE-vs-discrete coverage change for doxorubicin is
-0.586 while the dose-controlled change is -0.0166 — the naive comparison overstates the
structural effect 35x.

First vascular quantity to score against a MEASUREMENT rather than a threshold: vascular
area fraction 2.50% vs measured 2.1+-0.3% (PMID 36611806), 1.19x, untuned.

Status: **discrete vessels ARE now the global default** (`results/vessel_promotion.md`).
Both blockers were tested and neither survived. **alpha was never calibrated against the
KDE** — `calibrate_alpha` is a closed-form solve on the compiled drug plus the published
anchor, contains no vascular field, and returns 148.608... under either source;
`pk_multidrug` is byte-identical and still passes 2/2 out-of-sample. Forcing an *empirical*
section-scale re-fit instead BREAKS that check under both fields (mitoxantrone falls out of
band under the KDE-empirical alpha 455.2, topotecan under the discrete-empirical 77.6), so
148.6 stands. **Grid convergence came out inverted**: the discrete source converges on
coverage, AUC and kill by dx = 12 um and on fitted L by dx = 8 um, while the
**max-normalised KDE's total delivered dose does NOT converge** (-6.4% from dx = 8 to
dx = 6, outside the 5% tolerance) because it is divided by a grid-dependent maximum.

Two real costs, both quantified: distance-resolved claims need **dx <= 8 um** (the fitted L
of the shortest-L drug sits 24% low at dx = 12), and a concentrated source steepens the IFP
gradient so the advection CFL needs **dt <= 5 s** at dx <= 8 (the historical dt = 10 s
violates it). Promotion also changes the **unit** of every absolute exposure by 3.36x, so
KDE-era absolute numbers across `results/` are stale — see `vessel_promotion.md` §7.3 for
the list. §8.3's "the discrete field is the right substrate for the BBB leak" is WITHDRAWN:
the additive leak is a scalar on the source and both substrates transmit it to within 0.19%.

---

### Original (superseded) diagnosis, retained for the record

Three independent analyses converged on one root cause. The vascular source field is a
**smooth Gaussian KDE of endothelial cell positions** (bandwidth ~40 um), so:

- fitted penetration length is **26-28 um for ALL five drugs** against a **source-null of
  23.7 um** — the decay measures KDE bandwidth, not drug transport (`delivery_metrics.md`);
- the exposure profile has a **flat floor at ~12% of peak**, so no tissue is truly
  drug-starved and *delivery is never the binding constraint* (`dose_schedule.md`);
- "perivascular regions get more drug" is **circular** — the region label and the drug
  source share the endothelial signal (`niche_response.md`).

Real tumour vasculature is **discrete**: distinct vessels with genuinely avascular tissue
between them. Replacing the KDE with a segmented vessel network is the single
highest-value structural fix outstanding. Literature basis is established:
**~20 endothelial cells per mm of vessel, 1-2 EC per capillary cross-section**
(Shapson-Coe 2024 *Science* volume EM; Sargent 2023) — i.e. one EC ~ 50 um of vessel.

---

## 4. What survives the audit

- **Reconstruction fidelity**: cross-platform Pearson 0.997 — **RETRACTED, see results/revalidation.md**.
  Re-run under the fixed classifier it is 0.980 on the 4-vector and **0.016** on the non-tumour
  subvector (rank rho 1.00 -> 0.40). Attribution proven: reverting ONLY the classifier
  aggregation to legacy mean_z reproduces 0.9974/0.9450 to four decimals, so the arity fix
  accounts for 100% of the move. The 0.997 was agreement between a BIASED estimator and
  another estimator; its tightness was coincidental. Cross-platform concordance is an OPEN
  QUESTION, not a restored pass. Do not cite 0.997. Original text follows:
  cross-platform Pearson 0.997; endothelial agrees to 0.03%
  between two independent platforms on the SAME section.
  **RETRACTED (held-out test, results/heldout_patient.md): the external agreement does NOT
  replicate.** On a second real patient the endothelial fraction is 4.18% vs the 8.0%
  volume-EM value (0.52x), against section 1's 9.74% — and the areal endothelial density
  actually rose 2.67x, so the fraction fell only because section 2's tissue is 6.2x more
  cellular. The vascular area fraction fails harder: 5.79% vs measured 2.1+-0.3% (2.76x),
  against section 1's 2.50%. Cause is exact — area fraction = vessel count x a FIXED 7.5 um
  painted disc / tissue area, and 7.5 um is GRID-IMPOSED (a dx=12 um cell already has a
  6.8 um equivalent radius). Section 1's 1.19x agreement was partly a coincidence of its
  endothelial density matching that grid constant. MVD is high on BOTH sections, i.e.
  systematic rather than patient-specific. The two quantities this project cited as its
  external validations are exactly the two that fail, and both fail through constants
  entangled with the grid and with panel sparsity in a way one section could never expose.
  **UPDATE (this pass this text is stale on one point):** `results/heldout_patient.md` was
  regenerated after later fixes (discrete-vessel promotion `ceb96c0`, classifier-arity fix
  `4ad221f`, plus this pass's own vessel-distance/oxygen fixes). Re-run, the endothelial-
  fraction-vs-volume-EM claim now **REPLICATES** (7.07% vs 8.0%, against section 1's
  10.50%) — a genuine, incidental improvement, not yet earned by design. The vascular
  area fraction claim (now 10.98% vs 2.1%, 5.2x) and MVD (not reproduced on either
  section) are unchanged and still fail for the reasons described above. Do not cite
  4.18%/0.52x for the endothelial-fraction check; cite the current 7.07%/0.88x instead.
- **Numerics**: mass conserved to ~1e-13 on the real irregular mask; dt-converged; kill
  fraction grid-converged beyond dx<=20 um. Production resolution is **dx = 12 um** (also
  the finest grid satisfying the advection CFL at dt = 10 s).
- **Physics vs real PK**: one constant fit to doxorubicin predicts mitoxantrone (48 um) and
  chemically-distinct topotecan (84 um) inside their real cited bands, out-of-sample.
- **Potency**: GDSC vs HGCC agree at Spearman +0.83 (p=2.4e-5), real-vs-real.
- **The delivery x potency result**: ranking by **AUC/IC50** reproduces simulated kill
  ordering exactly. Delivery spans 120x, potency 1212x; the product decides. Doxorubicin
  is 49x more potent than cisplatin, 120x less delivered, and loses.
- **Recommended headline metric**: `therapeutic_coverage @ 1x IC50` — the fraction of tumour
  cells whose exposure exceeds the drug's measured IC50. Dimensionally sound,
  grid-converged (+-0.03), discriminating (0.311-0.925), and it recovers the kill ordering
  without touching the Hill curve, sensitivities, hypoxia factors or the Bernoulli draw.
  It is computed over the *fixed set of real cells*, which is why it does not inherit the
  KDE geometry confound.

---

## 5. The data landscape, settled

A 39-agent sweep (164 datasets, 143 open) established (`results/DATA_HUNT.md`):

- **No open dataset anywhere pairs measured spatial drug distribution with measured
  treatment response in the same tumour.** Not one.
- **No glioma + temozolomide spatial drug map exists**, open or otherwise (METASPACE: zero).
- **MGMT is absent from every Xenium panel** — so temozolomide sensitivity is not merely
  unmodelled but *unmodellable* from this modality.
- The experiment that would validate this project end-to-end already exists — the
  Agar/White adavosertib GBM platform, drug concentration per pixel in phase-1 patient
  tissue with molecular response — and has **no accession**. It is a months-long
  collaboration request, not a download.

What *is* open and downloaded: GIST imatinib MSI (6 matched human tumour/normal pairs, with
a quantitative calibration series), Moncayo HGSOC explants (114 ROIs pairing measured drug
with 11,224 genes), the human glioma regional PK compendium, and GSE232173 (glioma + TMZ +
measured ex-vivo response).

---

## 6. Outstanding work — status of the seven items listed when this section was written

All seven were listed as open after the first spend-limit interruption. Six have since been
closed, several with answers opposite to the ones expected:

| # | item | status |
|---|---|---|
| 1 | veto-shaped `bbb_penetrance` | **CLOSED** — defect #9; a mean cannot veto. Eribulin's error fell 66× → 6.9× |
| 2 | geometric tissue mask (fix #7) | **CLOSED** — and it turned out to be *partially cancelling* defect #6, which is why both had to be fixed together |
| 3 | discrete vessel network | **CLOSED** — [`vessel_promotion.md`](vessel_promotion.md); now the default. But the premise was wrong: the flat floor is **real diffusion physics** (L_rd 73–296 µm vs a 40 µm intervessel gap), not KDE blur, and drugs *did* separate — `penetration_length_L` spread went 1.07× → 6.91× |
| 4 | GIST imatinib, Moncayo, GSE232173 TMZ | **CLOSED** — [`gist_imatinib.md`](gist_imatinib.md), [`kill_calibration.md`](kill_calibration.md) |
| 5 | cell-state vs MKI67; oxygen vs VEGFA/HIF1A | **CLOSED, negative** — [`cell_state_impact.md`](cell_state_impact.md); [`oxygen_demotion.md`](oxygen_demotion.md) falsified the oxygen PDE as a value-add |
| 6 | cross-drug comparability modes | **CLOSED** — defect #8; `sensitivity_mode` uniform vs determinant |
| 7 | ranking robustness (cisplatin vs doxorubicin ~4% apart) | **STILL OPEN in spirit.** On the discrete-vessel default the gap is now wide (0.436 vs 0.259), but that is a *different* configuration, not a resolution of the original tie. Re-state before any top-1 call is made |

### Added since, and now open

- **Four mechanisms tested and rejected**, each by a different kind of argument, none by a
  failure to fit: efflux (wrong sign; provably inert on ratios), uptake (excluded by a 7.8×
  literature ceiling; OCT1 falsified), ion-trapping (excluded by a 3.09× measured bound, and a
  uniform R gives zero dispersion at equilibrium), immune sequestration (premise unsupported —
  the best-quantified case is a *depot* whose depletion **lowered** intratumoral drug 24% —
  and it worsens the anchor monotonically, 1.995× → 2.237×). The residual glioma error is
  **not** explained by any transport mechanism we have been able to name.
- **Therapeutic coverage carries a second, larger error bar** — up to 0.137 absolute from the
  2-D-section assumption, ~7× its grid bar ([`dimensionality.md`](dimensionality.md)). Mean
  tumour AUC is unbiased analytically; the transport physics needs no fix.
- **Schedule is structurally inert in this model — do not run regimen-optimisation on it.**
  Across ten real temozolomide regimens at matched dose ([`regimen.md`](regimen.md)), the
  total spread in mean tumour exposure is **0.000359%**, of which only **0.000013%** is
  transport; the remainder is the observation window clipping the last plasma tail, an
  artefact **27:1 larger than the physics**. All ten arms return *the same* drawn kill
  fraction. Staggering cannot change combination overlap either (largest Jaccard change
  0.014%), and the reason is structural rather than empirical: **the transport operator has
  no drug–drug coupling and Bliss is timing-blind**, so a stagger acts only through
  amplitude, never geometry. What would have to exist for schedule to matter — MGMT/AGAT
  depletion, inter-fraction repair, redistribution, regrowth, immune effects — is absent by
  construction. As-prescribed the ranking is monotone in cumulative dose and nothing else,
  which is precisely the prediction **RTOG 0525 falsified in humans** (2.1× dose, no OS
  benefit, double grade ≥3 toxicity).
- **The mechanism labels conflate two axes.** `UNREACHED` and `UNDERDOSED` are bands on one
  variable (E/IC50), so a drug that is merely *weak* is reported as *undelivered* — which is
  exactly what temozolomide does ([`mechanism_attribution.md`](mechanism_attribution.md)).
  The attribution cannot separate "too little arrived" from "the drug is too weak", and the
  UI now says so on screen. Separating them needs a second axis, not a re-labelling.

---

## 7. The honest position

The physics engine is numerically sound. The delivery/partition axis is validated in human
glioma with a quantified error budget. The spatial axis is **no longer confounded by our own
smoothing** — discrete vessels are the default and the flat floor turned out to be real
physics. The response axis still cannot be validated with open data at all, and that has not
moved: no open dataset pairs spatial drug distribution with measured response.

Two things must be said that were not in the original version of this section.

**The reconstruction is *not* externally validated — that claim was retracted.** Cross-platform
concordance fell 0.997 → 0.980 on the 4-vector and **0.016** on the non-tumour subvector, and
the drop is 100% attributable to fixing the classifier. All three reconstruction-level external
validations have now been retracted; every drug- and physics-level validation has held. That
asymmetry is the most informative pattern in the audit and it is not flattering to the part of
the model that is easiest to make look impressive.

**The residual glioma error survived four independent attempts to explain it mechanistically.**
That is a real limit, not a to-do item. It should be quoted as a standing ~2× geometric-mean
fold error rather than treated as something the next mechanism will absorb.

That is a real model with a real error budget — a fundamentally different object from a
simulator that produces confident numbers. The defects above were found because the model was
attacked, not because it failed; a version of this project that skipped the audit would still
be reporting nine wrong numbers with beautiful figures attached.
