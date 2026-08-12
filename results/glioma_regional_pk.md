# Validation run #1 — predicted vs MEASURED drug partition in human glioma

**What this is.** The first time tmesim's molecular/transport layer has been scored against measured drug concentrations in human brain-tumour tissue. It tests ONE axis: the **partition / delivery** axis — what fraction of the drug in plasma actually reaches glioma tissue, and whether that fraction tracks MRI-defined regional blood-brain-barrier disruption.

![predicted vs measured](glioma_regional_pk.png)

## 1. Data provenance

**Source 1 — literature anchors.** `data_cache/hunt/autoradiography_microdialysis/parsed/regional_intratumoral_PK_anchors_literature.csv`, 27 rows, every row carrying a PMID and a DOI. 21 human rows / 6 rat rows.

| study | PMID | drug | species | n rows | region(s) | method |
|---|---|---|---|---|---|---|
| Blakeley 2008 J Neurooncol | 18787762 | methotrexate | human | 4 | contrast_enhancing, nonenhancing | intracerebral_microdialysis |
| Eroglu 2024 Cancer Chemother Pharmacol | 39422741 | eribulin | human | 1 | peritumoral_brain | intracerebral_microdialysis |
| Grossman 2012 Cancer Chemother Pharmacol | 22644796 | temozolomide | rat | 4 | U87_intratumoral_ECF, contralateral_normal_brain_ECF | intratumoral_microdialysis |
| Grossman 2013 Cancer Chemother Pharmacol | 23649683 | temozolomide | rat | 2 | U87_intratumoral_ECF | intratumoral_microdialysis |
| Jackson 2018 Fluids Barriers CNS | 29332604 | temozolomide | human | 7 | contrast_enhancing, noncontrast_enhancing | intracerebral_microdialysis |
| Lee 2025 Clin Cancer Res (ABTC adavosertib) | 39820282 | adavosertib | human | 2 | contrast_enhancing, nonenhancing | tissue_homogenate_LCMS |
| Portnow 2009 Clin Cancer Res | 19861433 | temozolomide | human | 7 | peritumoral_brain | intracerebral_microdialysis |

Column meanings: `plasma_cmax` / `tissue_cmax` are peak concentrations, `plasma_auc` / `tissue_auc` are areas under the concentration-time curve in the units given by `units_conc` / `units_auc` (they differ per study — ug/mL vs uM — which is exactly why only the **dimensionless ratio** `tissue_to_plasma_auc_pct` is used here). That column is divided by 100 on load to give a fraction. Two rows are not per-subject: Lee 2025 reports a cohort mean intratumoral:plasma ratio x100 (418 % and 74 %), and Eroglu 2024 reports a range string `0.13-1.99` (%), from which the geometric midpoint 0.51 % is taken.

**Rat rows are excluded from the scorecard.** The 6 Grossman 2012/2013 rows report intratumoral / contralateral ECF AUC with **no paired plasma AUC**, so no tissue:plasma ratio exists. They can only support a within-animal ordering (tumour ECF 3.35 vs contralateral normal-brain ECF 3.22 ug*h/mL pre-bevacizumab — i.e. essentially NO tumour-vs-normal-brain enrichment for TMZ in that model), which is itself a useful negative but is not scoreable as a partition.

**Source 2 — Oh et al. 2026 Neuro-Oncology (`noaf204`, PMC12962645) Table S3.** `data_cache/hunt/autoradiography_microdialysis/parsed/oh2026_neurooncol_TableS3_per_sample_tumor_drug_ratios.csv`: **125 surgical samples from 42 high-grade-glioma cases**, a window-of-opportunity study of BBB permeability heterogeneity. Two probe drugs were given intra-operatively and assayed in MRI-region-targeted resection specimens:

- `lev_tumor_to_plasma_ratio` = **levetiracetam**, small/neutral, BBB-penetrant (125 values, mean 0.697, range 0.00-2.05)
- `cef_tumor_to_plasma_ratio` = **cefazolin**, anionic cephalosporin, BBB-excluded (117 values, mean 0.150, range 0.00-1.47) — 21 of these are reported as exactly 0.00, i.e. **left-censored at the assay LLOQ**, not true zeros; they are kept at 0 here, which biases the cefazolin mean and any parametric test downward and is one reason the nonparametric statistics below are the primary ones.

These ratios are **single-timepoint tissue-concentration : plasma-concentration** ratios from tissue homogenate, NOT AUC ratios — sampling times (`lev_collection_time_min`, `cef_collection_time_min`) differ per sample, so each ratio is a snapshot on a distribution/elimination curve rather than an exposure ratio. Region flags: `region_contrast_enhancing` (T1-Gd) 39(+)/69(-), `region_flair` 73(+)/35(-), and 17 pathologist-called necrotic samples (whose C/F flags are N/A). `mannitol_status` records minutes since intra-operative mannitol for 40 samples — mannitol osmotically opens the barrier, so it is a confounder and a mannitol-naive sensitivity analysis is run.

## 2. How the model produces a number (assumption A1)

The frozen engine's delivery term is `source = PS * plasma_mean` with `PS = PS_SCALE_S_INV * vascular_density * bbb_penetrance` (1/s) (`tmesim/sim/engine.py`, `tmesim/treatment/drugs.py`), and `tmesim/core/config.py` defines `bbb_penetrance` as a dimensionless brain:plasma **partition coefficient** in [0, 1] — not a permeability, not a rate. (~~When this report was written that field read *"[0,1] fraction of vascular concentration that crosses into"* tissue~~ — **CORRECTED** to the partition wording, same quantity, stated in the units it was built in; see §2.1.) With `vascular_density` normalised to 1 and the tissue at quasi-equilibrium with a sustained plasma level — the Stage-0 plasma-averaging approximation the engine already makes — the tissue:plasma exposure ratio the model implies is the identity:

> **A1:  R_predicted = `bbb_penetrance(descriptors)`**

`bbb_penetrance` is a CNS-MPO-style structure-only score over MW / logP / TPSA / HBD (`tmesim/molecular/permeability.py`). No fitting, no calibration to any of the data below. Three consequences of A1, stated up front:

1. **R_predicted <= 1 by construction.** Any measured ratio above 1 is unreachable by the model at any parameterisation. Adavosertib's measured contrast-enhancing ratio is 4.18. That is a structural ceiling failure, not a tuning failure.
2. **R_predicted carries no region term.** `bbb_penetrance` is a function of the drug's descriptors only. The model's predicted CE-vs-non-enhancing fold change is **exactly 1.000** for every drug.
3. **R_predicted is a passive free-drug crossing fraction.** Microdialysis (Portnow, Blakeley, Jackson, Eroglu) measures *free* extracellular drug and is the modality-matched comparator. Tissue homogenate (Lee 2025, Oh 2026) measures *total* drug including intracellular and lipid-bound fractions and can exceed 1. Every scorecard row is labelled with its modality.

### 2.1 Why A1 holds, and under exactly which convention

> ### ⚠ CORRECTION NOTICE — added after this report was written
>
> A1 was stated as a bare identity between the model's tissue:plasma ratio and
> `bbb_penetrance`, which is right about the quantity and silent about the
> convention. `results/kill_calibration.md` §2 stated the opposite in its own
> context — that the delivery term "is NOT a physical tissue:plasma partition
> coefficient and must not be read as one". `results/transvascular_exchange_revalidation.md`
> §6.2 found the contradiction and settled it: **both were half right.** A1's
> reading of `bbb_penetrance` stands and is the correct one; what it needs is the
> convention under which the identity is exact, which is below. The original
> quotation of the `config.py` wording is struck through in §2 rather than
> deleted. **No predicted, measured or fold-error number in this file changes.**

The settled reading, which `core/config.py`, `reconstruct/barrier.py`,
`treatment/drugs.py` and `sim/engine.py` now all state:

* **`bbb_penetrance` IS a dimensionless brain:plasma partition coefficient** — what A1 assumes, what it was built as, and what `molecular/permeability.py` calibrates against measured tissue:plasma ratios.
* **What the PDE integrates is a first-order RATE**, `PS = PS_SCALE_S_INV x vascular_density x bbb_penetrance` (1/s). Both spatial factors are dimensionless; `treatment/drugs.py::PS_SCALE_S_INV` carries the whole `1/s`. The engine's historical `uM/s` came from an implicit factor of exactly 1 s^-1 that nothing derived; it is now named, defaults to **1.0** (bit-identical), and has a **calibrated value of 0.16122 1/s** (`tmesim/validation/ps_scale.py`, fit to measured rodent brain PS: methotrexate 3.3e-5 /s PMID 524251, temozolomide 1.6e-3 /s PMID 39690422).
* The two readings are not interchangeable, and the diagnostic in §6.2 is which one the drugs are coherent under: as **rates** methotrexate and temozolomide disagree by **8.7x** against the same normal rat brain; as **partitions** they agree to **1.04x**.

**So A1 is a Stage-0 identity, not a general one.** The quasi-equilibrium tissue:plasma ratio the engine implies is

> **Kp = PS / (PS + k_loss)**,  `PS = PS_SCALE_S_INV * vascular_density * bbb_penetrance`,  `k_loss = clearance + binding`

(the gradient-driven exchange form; under the shipped one-way source the same steady state is `PS*Cp / k_loss`, which is not bounded by `Cp` at all — the defect `results/breast_kp_validation.md` found, and a second reason the raw product `vascular_density x bbb_penetrance x plasma` is not itself a partition), which collapses to `R_predicted = bbb_penetrance` only under the **Stage-0 convention** this report runs in: `vascular_density` normalised to 1, a sustained plasma level, and the delivery term read as the partition it is built from rather than as the rate it enters the PDE as (equivalently `PS_SCALE_S_INV = 1 s^-1`, its shipped default and exactly the implicit factor the engine used before the constant was named). Two things follow, and both are already visible in this report:

1. **The ceiling in consequence 1 is a property of the convention as much as of the score.** `Kp <= 1` for any non-negative `PS` and `k_loss`, so adavosertib's measured 4.18 is unreachable under A1 *and* under the full Kp form. Correcting the units cannot lift it; only a modality-matched comparator (consequence 3) or a term the model does not have can.
2. **A calibrated rate scale tightens the ceiling rather than relieving it.** `Kp` is monotone increasing in `PS`, so `PS_SCALE_S_INV = 0.16122 1/s` puts every drug strictly below its A1 value. The scorecard below is reported at the Stage-0 convention, and any re-score at the calibrated scale is an over-prediction that gets smaller, not a ceiling that moves.

Cross-references: `results/kill_calibration.md` §2.1 states the same resolution from the rate side, and `results/transvascular_exchange_revalidation.md` §6.2 carries the diagnostic and the measured anchors.

### 2.2 This report is UN-TESTABLE by the exchange flag, and §3.1 of the re-validation is not a re-test of it

Two statements in `results/transvascular_exchange_revalidation.md` look like they
contradict each other about this file, and a reader is entitled to be confused by
them. §6.1 classes `glioma_regional_pk` as **analytic and un-testable** by
`PDEConfig(transvascular_exchange=True)`, because A1 is a closed-form structure
score and no transport solver is ever invoked on this path — there is no
injection point, and passing `exchange_k` is not an available remedy either. §3.1
reports a **partition scorecard improving** under exchange, GM |fold error|
5.235x → 3.969x on the modality-matched microdialysis rows.

**Both are true, because they are not the same computation.** §3.1 scores
engine-computed interstitial concentrations against these anchors; this report
scores `bbb_penetrance` itself, through A1, with no engine run behind it. The
anchors are shared; the predictor is not.

The consequence that matters: **§3.1 must not be read as re-validating A1, and
nothing in this file should be quoted as having survived the delivery-law
correction.** It did not survive it — it was never exposed to it. "Robust under
exchange" and "never evaluates the delivery law" are different claims, and §6.1
exists precisely to stop them being pooled into a coverage count.

## 3. The scorecard — predicted vs measured, per drug

| drug | MW | logP | TPSA | HBD | **predicted** (A1) | **measured** mean | measured range | n | modality | fold error | 3x | 10x | PMID |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eribulin | 729.9 | 3.44 | 146.4 | 2 | **0.334** | **0.0051** | 0.0051-0.005 | 1 | microdialysis (free) | **65.63x** | FAIL | FAIL | 39422741 |
| cefazolin | 454.5 | -0.64 | 156.1 | 2 | **0.248** | **0.1496** | 0.0000-1.470 | 117 | homogenate (total) | **1.66x** | PASS | PASS | PMC12962645 |
| methotrexate | 454.4 | 0.27 | 210.5 | 5 | **0.115** | **0.1780** | 0.0320-0.305 | 4 | microdialysis (free) | **0.65x** | PASS | PASS | 18787762 |
| temozolomide | 194.2 | -2.08 | 108.2 | 1 | **0.558** | **0.1951** | 0.0330-0.449 | 14 | microdialysis (free) | **2.86x** | PASS | PASS | 19861433|29332604 |
| levetiracetam | 170.2 | -0.13 | 63.4 | 1 | **0.904** | **0.6972** | 0.0000-2.050 | 125 | homogenate (total) | **1.30x** | PASS | PASS | PMC12962645 |
| adavosertib | 500.6 | 2.89 | 104.3 | 2 | **0.639** | **2.4600** | 0.7400-4.180 | 2 | homogenate (total) | **0.26x** | FAIL | PASS | 39820282 |

- **Order-of-magnitude (10x) band: 5/6 pass.**
- Tight (3x) band: 4/6 pass.
- **Aggregate systematic bias: geometric-mean fold error = 2.02x** (OVER-prediction of brain-tumour drug exposure).
- Excluding the single worst drug (eribulin, 65.6x): geometric-mean fold error = 1.01x — i.e. the aggregate bias is dominated by one drug, not spread evenly.
- **Modality-matched subset** (microdialysis / free drug only — the only apples-to-apples comparison, n=3): geometric-mean fold error = 4.95x, i.e. *worse* than the pooled figure, again because eribulin lives in this subset. Excluding it: 1.36x (temozolomide and methotrexate only, n=2 — too few drugs to call this a calibration).

Per-drug, honestly:

- **eribulin** (microtubule inhibitor; large macrocycle, strong P-gp substrate): predicted 0.334 vs measured 0.0051 = **65.63x**. Eroglu 2024 Cancer Chemother Pharmacol.
- **cefazolin** (anionic cephalosporin, BBB-excluded probe (Oh 2026 probe 2)): predicted 0.248 vs measured 0.1496 = **1.66x**. Oh 2026 Neuro-Oncology (noaf204) Table S3.
- **methotrexate** (antifolate; large, very polar, 5 HBD, actively effluxed): predicted 0.115 vs measured 0.1780 = **0.65x**. Blakeley 2008 J Neurooncol.
- **temozolomide** (first-line GBM alkylator; small, neutral, low logP): predicted 0.558 vs measured 0.1951 = **2.86x**. Jackson 2018 Fluids Barriers CNS; Portnow 2009 Clin Cancer Res.
- **levetiracetam** (antiepileptic BBB-penetrant probe (Oh 2026 probe 1)): predicted 0.904 vs measured 0.6972 = **1.30x**. Oh 2026 Neuro-Oncology (noaf204) Table S3.
- **adavosertib** (WEE1 inhibitor (ABTC trial); lipophilic, P-gp substrate): predicted 0.639 vs measured 2.4600 = **0.26x**. Lee 2025 Clin Cancer Res (ABTC adavosertib).

## 4. The regional test

n = 125 samples from 42 cases. Region counts: `contrast_enhancing` 39, `necrosis` 17, `non_enhancing` 69.

Two permutation tests are reported for every contrast. The **naive** one shuffles region labels freely across all samples and treats each sample as independent — it is ANTI-CONSERVATIVE here, because up to 8 samples come from the same tumour. The **within-case** one shuffles the region label only *within* each patient, preserving the patient as the clustering unit and testing the paired within-tumour contrast. **The within-case p is the one to believe.** Both use 20,000 permutations, seeded.

### levetiracetam

| region | n | mean | median | min | max |
|---|---|---|---|---|---|
| contrast_enhancing | 39 | 0.875 | 0.890 | 0.31 | 2.05 |
| necrosis | 17 | 0.674 | 0.690 | 0.00 | 1.57 |
| non_enhancing | 69 | 0.602 | 0.590 | 0.16 | 1.32 |

| contrast | mean A | mean B | ratio A/B | diff A-B | Cliff's d | Hedges' g | p (naive) | p (within-case) | informative cases | MODEL ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| contrast_enhancing vs non_enhancing | 0.875 | 0.602 | **1.45x** | +0.272 | +0.590 | 1.11 | 0.0000 | **0.0000** | 31 | 1.00x |
| contrast_enhancing vs necrosis | 0.875 | 0.674 | **1.30x** | +0.201 | +0.341 | 0.59 | 0.0435 | **0.0911** | 14 | 1.00x |
| necrosis vs non_enhancing | 0.674 | 0.602 | **1.12x** | +0.072 | +0.106 | 0.28 | 0.3117 | **0.3729** | 13 | 1.00x |

Mannitol-naive sensitivity (CE vs non-enhancing, n=25 vs 48): ratio 1.45x, Cliff's d +0.593, within-case p = 0.0002 (22 informative cases). The regional effect does not depend on mannitol exposure.

### cefazolin

| region | n | mean | median | min | max |
|---|---|---|---|---|---|
| contrast_enhancing | 36 | 0.276 | 0.280 | 0.00 | 1.25 |
| necrosis | 15 | 0.220 | 0.100 | 0.00 | 1.47 |
| non_enhancing | 66 | 0.065 | 0.035 | 0.00 | 0.44 |

| contrast | mean A | mean B | ratio A/B | diff A-B | Cliff's d | Hedges' g | p (naive) | p (within-case) | informative cases | MODEL ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| contrast_enhancing vs non_enhancing | 0.276 | 0.065 | **4.26x** | +0.211 | +0.641 | 1.39 | 0.0000 | **0.0000** | 30 | 1.00x |
| contrast_enhancing vs necrosis | 0.276 | 0.220 | **1.25x** | +0.056 | +0.344 | 0.20 | 0.5439 | **0.6007** | 13 | 1.00x |
| necrosis vs non_enhancing | 0.220 | 0.065 | **3.40x** | +0.155 | +0.241 | 0.89 | 0.0038 | **0.0166** | 12 | 1.00x |

Mannitol-naive sensitivity (CE vs non-enhancing, n=25 vs 48): ratio 4.38x, Cliff's d +0.598, within-case p = 0.0000 (22 informative cases). The regional effect does not depend on mannitol exposure.

**Independent literature replication of the regional ordering** (different patients, different drugs, different assay):

| drug | n CE | n non-enh | mean CE | mean non-enh | CE/non-enh | CE - non-enh | method | PMID |
|---|---|---|---|---|---|---|---|---|
| adavosertib | 1 | 1 | 4.180 | 0.740 | **5.65x** | +3.440 | tissue_homogenate_LCMS | 39820282 |
| methotrexate | 2 | 2 | 0.293 | 0.063 | **4.65x** | +0.230 | intracerebral_microdialysis | 18787762 |
| temozolomide | 2 | 5 | 0.265 | 0.191 | **1.39x** | +0.074 | intracerebral_microdialysis | 19861433|29332604 |

## 5. The model has NO region-specific BBB mechanism — and here is what one would have to be

Stated plainly: **`bbb_penetrance` is a pure function of the drug's 2D descriptors. There is no term in the molecular layer, the transport PDE, or `DrugParams` that represents regional barrier integrity.** The engine's only spatial modulation of delivery is `vascular_density`, which is vessel *density*, not barrier *leakiness* — a region can have identical vessel density and a completely different permeability-surface-area product. So the model's predicted CE-vs-non-enhancing fold change is exactly 1.00 against measured fold changes of 1.45x (levetiracetam) and 4.26x (cefazolin). **The model fails the regional test completely — it cannot even express the answer.**

What would have to be added, and what its value would have to be. Fit two candidate one-parameter region terms to the CE-vs-non-enhancing contrast in three independent drug/assay combinations. A term that is genuinely a property of the *barrier* rather than of the *drug* must come out drug-INDEPENDENT; the spread across drugs is the diagnostic.

| drug | source | mean CE | mean non-enh | multiplicative f = CE/nonCE | additive leak L = CE-nonCE |
|---|---|---|---|---|---|
| levetiracetam | Oh 2026 (homogenate) | 0.875 | 0.602 | 1.45 | 0.272 |
| cefazolin | Oh 2026 (homogenate) | 0.276 | 0.065 | 4.26 | 0.211 |
| methotrexate | Blakeley 2008 (microdialysis) | 0.293 | 0.063 | 4.65 | 0.230 |

- multiplicative `f` spread across drugs: **3.2x**
- additive leak `L` spread across drugs: **1.3x** (L = 0.211 to 0.272, mean **0.238**)

**The additive form wins** (1.3x spread vs 3.2x). That is mechanistically interpretable and it is the concrete specification of the missing feature: contrast enhancement corresponds to a **drug-independent parallel leak conductance** that adds roughly a constant increment (~0.24 in tissue:plasma units) to the partition regardless of the molecule's own physicochemistry — exactly what a size-nonselective paracellular opening should do. A multiplicative barrier factor cannot work: it would have to be ~4.7x for the barrier-excluded drugs and only ~1.5x for the already-penetrant one, because a drug that is already near-freely permeable has no headroom left to gain.

Concretely, to pass this test the model would need:

1. A **regional barrier field** — a per-voxel/per-region `bbb_disruption` in [0,1] derived from contrast enhancement (or a Ktrans-like DCE-MRI parameter), carried alongside `vascular_density` in the microenvironment fields.
2. A **delivery term of the form** `source = vascular_density * (bbb_penetrance + L * bbb_disruption) * plasma`, i.e. an additive leak in parallel with the drug's intrinsic transcellular route, with `L` calibrated to ~0.24 and capped so the partition cannot exceed the free-drug ceiling.
3. A **separate total-vs-free bookkeeping** so that homogenate-measured ratios above 1 (adavosertib 4.18) are represented as tissue binding rather than as impossible barrier crossing.

## 6. The honest verdict

**Partition magnitude: partial pass.** 5/6 drugs land inside an order of magnitude and 4/6 inside 3x. The aggregate geometric-mean fold error is 2.02x — the model **over-predicts** human glioma drug exposure — but that aggregate is dominated by one drug (eribulin, 66x); with it removed the bias falls to 1.01x, i.e. close to unbiased on the remaining drugs.

**The two failures are two DIFFERENT failures, and neither is a mis-tuned constant.** They should not be lumped together:

- **Eribulin, 66x over-predicted — a missing efflux mechanism, exposed by an averaging rule that cannot express a veto.** Eribulin is a 730 Da macrocycle and a well-characterised P-glycoprotein substrate; its measured peritumoral brain:plasma ratio is 0.13-1.99 % (Eroglu 2024, PMID 39422741). `bbb_penetrance` still scores it **0.334**. Two things go wrong at once. (i) The function combines its four descriptor desirabilities by taking their *mean*, so eribulin's two disqualifying properties — MW 730 (desirability ~0.01) and TPSA 146 (~0.01) — are averaged against its two favourable ones — logP 3.44 (~0.86) and 2 HBD (~0.78) — and the molecule is rescued to a mid-range score instead of being vetoed. A mean cannot express "this is disqualifying"; a minimum, or a product, could. (ii) Even a correct passive score would be too high, because the real limit on eribulin is active efflux, which the function's own docstring concedes it "cannot capture". Efflux is the single largest source of error in the delivery axis, and this is the cheapest available fix to the second-largest.
- **Adavosertib, 3.8x under-predicted — a ceiling artefact, not a physics error.** Its measured contrast-enhancing ratio is 4.18, above the model's structural maximum of 1.0. But that 4.18 is a *tissue homogenate* number (total drug, including everything bound inside cells); a free-drug partition of 4.18 would require active uptake INTO brain. The model cannot be right here because it is not modelling the same quantity. Its predicted 0.639 is in fact close to adavosertib's *non-enhancing* measured ratio (0.74, 0.86x).

**A caveat in the model's own docstring turns out to be backwards here.** `permeability.bbb_penetrance` documents that it UNDER-predicts temozolomide, because TMZ's logP and TPSA are outside the CNS-favourable window even though TMZ is clinically a BBB-crossing drug. Against these measurements it does the opposite: predicted 0.558 vs a measured mean of 0.195 across 14 human subjects (Portnow 2009 + Jackson 2018), a 2.9x OVER-prediction. The documented caveat is about TMZ's *rank* among CNS-like molecules; as an *absolute partition* the score is too generous. Both statements can be true, but only the second one is the one that matters when the number is used as a delivery fraction.

**Ordering across drugs: pass — and this is the strongest positive result in the run.** Oh 2026 measured a BBB-penetrant and a BBB-excluded drug in the *same* tissue samples, which removes every between-patient confounder. Measured separation: 4.66x pooled (0.697 vs 0.150), rising to 9.3x in non-enhancing tissue where the barrier is intact. The model, from 2D structure alone with no fitting, predicts 0.840 vs 0.049 = 17.20x — right direction, right order of magnitude for the separation. It does under-state the gap in intact-barrier tissue, which is the same efflux-blindness that sinks eribulin, just in a milder form.

**Regional ordering: outright fail.** The measured regional effect is real and strong — contrast-enhancing tissue carries 1.45x for levetiracetam (within-case p = 0.0000, Cliff's d +0.59) and 4.26x for cefazolin (within-case p = 0.0000, Cliff's d +0.64) the partition of non-enhancing tissue — and it replicates independently in Blakeley 2008's methotrexate microdialysis. The model predicts a fold change of exactly 1.00 because it has no region-specific BBB mechanism at all. Section 5 gives the concrete specification of what to add.

**Necrosis does not behave as naively expected.** The prior hypothesis was CE > non-enhancing > necrosis. Measured, necrosis sits *between* CE and non-enhancing (and for cefazolin it is statistically indistinguishable from CE), which makes sense — necrotic cores sit inside the enhancing region and share its disrupted barrier, even though their perfusion is poor. Any region term added to the model should be driven by barrier integrity, not by a tissue-viability label.

### What this DOES establish

- The molecular layer's structure-only `bbb_penetrance`, used as a tissue:plasma partition, is within an order of magnitude of measured human glioma values for 5/6 drugs, and correctly rank-orders a BBB-penetrant against a BBB-excluded drug measured in the *same* human tissue samples.
- The measured regional heterogeneity of drug delivery in human GBM is real, large, and reproducible across two independent studies and four drugs — it is a legitimate, quantified validation target.
- The specific missing mechanism has been identified and its magnitude estimated (an additive, drug-independent leak of ~0.24 in contrast-enhancing tissue).

### What this does NOT establish

- **Nothing about response.** No cell kill, no survival, no potency. Delivery is necessary, not sufficient.
- **Nothing spatial.** These data are *regional* (an MRI-defined label per specimen); there are no coordinates, no distance-to-vessel, no gradients. The sim's perivascular penetration lengths are untouched by this run.
- **Nothing about time.** The Oh 2026 ratios are single-timepoint snapshots at varying collection times, not AUC ratios; assumption A1 equates a steady-state partition to both, which is only approximately right.
- **Nothing about the free/total distinction.** Two of the six drugs are scored against tissue homogenate (total drug), which is not the quantity A1 predicts. Their fold errors are contaminated by tissue binding.
- **Nothing about n.** Four of the six drugs rest on n <= 14 subjects, two of them on a single cohort summary. Only levetiracetam/cefazolin have real sample size, and those two are clustered within 42 patients.

---

Sections 1-6 above are validation run #1 **as it was run, before any fix**. The predicted column is recomputed live from the frozen pre-fix aggregation (`_bbb_penetrance_arithmetic`), not pasted, so it cannot drift. Everything below is the post-fix re-run.

# POST-FIX — the same tests, re-run after both fixes

Two fixes were implemented from the numbers above. **The pre-fix numbers in sections 1-6 are deliberately kept: the before/after IS the evidence.**

- **FIX 1 — `bbb_penetrance` can now express a veto.** `tmesim/molecular/permeability.py` combines the four CNS-MPO desirabilities with a **geometric** mean instead of an arithmetic one. No new parameter; for a molecule with four equal desirabilities the two forms are identical, so the change bites only on HETEROGENEOUS descriptor profiles — which is exactly the veto that was missing.
- **FIX 2 — a barrier-integrity term now exists.** New module `tmesim/reconstruct/barrier.py`: `effective_penetrance = clip(bbb_penetrance + L * bbb_disruption, 0, 1)` with `L = 0.238` (mean of the three fitted values 0.211 / 0.230 / 0.272 from section 5). Wired through `treatment/drugs.py::delivery_source` behind a default-off argument; `core`, `physics` and `sim/engine.py` are untouched.

## P1. FIX 1 — which aggregation, and why

Both candidates the audit named — a geometric mean and a min-weighted form — were evaluated against the same 6-drug measured scorecard, at `w = 0.5` and `w = 0.7` for `w*min + (1-w)*mean`, plus a pure `min`. The min-weighted forms and the pure min are summarised in prose below rather than shipped, since only the chosen path is recomputed live.

| aggregation (structure-only) | eribulin fold | worst \|fold\| | within 10x | within 3x | GM \|fold error\| |
|---|---|---|---|---|---|
| arithmetic mean (pre-fix) | 65.63x | 65.63x | 5/6 | 4/6 | 3.660x |
| **geometric mean (CHOSEN)** | 6.94x | 7.91x | 6/6 | 2/6 | 3.451x |
| 0.5*min + 0.5*mean | 17.76x | 17.76x | 5/6 | 3/6 | 3.475x |
| 0.7*min + 0.3*mean | 9.72x | 11.74x | 5/6 | 2/6 | 4.711x |
| pure min | 3.80x | 23.39x | 5/6 | 0/6 | 7.955x |

With the FIX 2 leak also active — the model as it actually ships — the two finalists separate cleanly: geometric mean gives GM |fold error| **1.995x**, worst case **6.94x**, 6/6 within 10x and 4/6 within 3x; `0.5*min + 0.5*mean` gives 2.284x, worst case 17.76x, 5/6 within 10x and 4/6 within 3x. Pure `min` is far too aggressive — it vetoes on the single worst descriptor and drives every drug to near-zero, losing the 3x band entirely.

The three rows recomputed live on every run:

| aggregation | eribulin | worst \|fold\| | within 10x | within 3x | GM signed bias | GM \|fold error\| |
|---|---|---|---|---|---|---|
| arithmetic mean (pre-fix defect) | 65.63x | 65.63x | 5/6 | 4/6 | 2.019x | 3.660x |
| geometric mean, structure-only | 6.94x | 7.91x | 6/6 | 2/6 | 0.588x | 3.451x |
| **geometric mean + leak (CHOSEN)** | 6.94x | 6.94x | 6/6 | 4/6 | 1.066x | 1.995x |

Read the last two columns together. **GM signed bias** lets over- and under-prediction cancel; **GM |fold error|** does not, and is the honest accuracy number. The pre-fix model's flattering *ex-eribulin* bias of 1.01x was precisely such a cancellation — temozolomide +2.86x and cefazolin +1.66x annihilating adavosertib 0.26x — while its unsigned error over the same five drugs was 2.06x.

**Chosen: the geometric mean.** Against the full post-fix model it gives the lowest unsigned error (1.99x vs 2.28x for min-weighted), the lowest worst-case (6.9x vs 17.8x), and it is the only form that brings all 6 drugs inside the order-of-magnitude band. It also adds no tunable parameter, whereas every min-weighted form introduces a weight `w` that would have to be fitted to these same six drugs.

## P2. The 6-drug scorecard, BEFORE and AFTER

Three columns, because the two fixes must be separable. *pre-fix* is the arithmetic mean with no barrier term (what sections 1-6 scored). *FIX 1 only* is the geometric mean, still with no barrier term. *FIX 1 + 2* is the model as it now stands, with the additive leak weighted by each drug's **own measured regional composition** (`measured_sample_disruption`) — because each measured mean pools whatever regions that study sampled, and scoring a region-aware model at the wrong regional mix would penalise it for sampling, not physics.

| drug | measured | modality | disrupted-region fraction of samples | pre-fix pred | pre-fix fold | FIX 1 pred | FIX 1 fold | FIX 1+2 pred | **FIX 1+2 fold** |
|---|---|---|---|---|---|---|---|---|---|
| eribulin | 0.0051 | microdialysis (free) | 0.00 | 0.334 | 65.63x | 0.035 | 6.94x | 0.035 | **6.94x** |
| cefazolin | 0.1496 | homogenate (total) | 0.44 | 0.248 | 1.66x | 0.049 | 0.33x | 0.153 | **1.02x** |
| methotrexate | 0.1780 | microdialysis (free) | 0.50 | 0.115 | 0.65x | 0.023 | 0.13x | 0.142 | **0.79x** |
| temozolomide | 0.1951 | microdialysis (free) | 0.14 | 0.558 | 2.86x | 0.125 | 0.64x | 0.159 | **0.82x** |
| levetiracetam | 0.6972 | homogenate (total) | 0.45 | 0.904 | 1.30x | 0.840 | 1.20x | 0.947 | **1.36x** |
| adavosertib | 2.4600 | homogenate (total) | 0.50 | 0.639 | 0.26x | 0.459 | 0.19x | 0.578 | **0.24x** |

| aggregate | pre-fix | FIX 1 only | FIX 1 + 2 |
|---|---|---|---|
| within 10x | 5/6 | 6/6 | 6/6 |
| within 3x | 4/6 | 2/6 | 4/6 |
| worst \|fold error\| | 65.63x | 7.91x | 6.94x |
| GM signed bias | 2.019x | 0.588x | 1.066x |
| GM signed bias, ex-eribulin | 1.006x | 0.359x | 0.733x |
| **GM \|fold error\|** | **3.660x** | **3.451x** | **1.995x** |
| GM \|fold error\|, ex-eribulin | 2.055x | 3.000x | 1.554x |

**Honest reading of the middle column.** FIX 1 *alone* does not simply improve things. It fixes the outlier — eribulin 65.6x -> 6.9x, worst case 65.6x -> 7.9x, and 6/6 drugs inside 10x for the first time — but it converts the model from an over-predictor into an under-predictor: the ex-eribulin signed bias moves from 1.01x to 0.36x and the tight 3x band drops from 4/6 to 2/6. That is a real cost and it is not hidden here.

**Why the third column recovers it.** The veto pushes methotrexate and cefazolin down because both genuinely carry disqualifying descriptors (TPSA 210 / 156). But their measured means POOL contrast-enhancing tissue where the barrier is broken — about half the samples for methotrexate, 44% for cefazolin. The structure-only score cannot represent that tissue at all. Once the additive leak supplies it, cefazolin lands at 1.02x and methotrexate at 0.79x. The two fixes are not independent improvements; **FIX 1 is only correct in the presence of FIX 2**, because the veto is a statement about the *intact* barrier and most of the measurements are not.

**Requirement check, stated plainly.** The brief required that aggregate fold-bias not worsen against the pre-fix 2.02x / 1.01x-ex-eribulin. On the signed measure the full model gives 1.07x (better than 2.02x) but 0.73x ex-eribulin, which is farther from 1.0 than the pre-fix 1.01x. **That requirement is not met on the signed ex-eribulin measure.** It is met on every non-cancelling measure: unsigned error 3.66x -> 1.99x overall and 2.05x -> 1.55x ex-eribulin, worst case 65.6x -> 6.9x, 10x band 5/6 -> 6/6. The pre-fix 1.01x was two large errors cancelling, so it is the wrong thing to protect.

### Orderings preserved

- Antibody (bevacizumab-like) penetrance = 0.0000 — still a hard 0, both via the `is_antibody` gate and the MW >= 1000 gate.
- Least- to most-penetrant among the six small molecules (structure-only): methotrexate 0.023 < eribulin 0.035 < cefazolin 0.049 < temozolomide 0.125 < adavosertib 0.459 < levetiracetam 0.840.
- Methotrexate is the least penetrant small molecule: **True** (it was also least pre-fix, so the veto preserves this).
- BBB-penetrant vs BBB-excluded probe separation (Oh 2026, same tissue samples, measured 4.66x pooled / 9.3x in intact-barrier tissue): model now predicts 17.20x (levetiracetam 0.840 / cefazolin 0.049), up from 3.65x pre-fix. Direction right, magnitude now **over**shot: the pre-fix 3.65x under-stated the intact-barrier 9.3x by 2.5x, the post-fix 17.2x over-states it by 1.8x. The veto swapped an under-separation for an over-separation of similar size; it did not eliminate the error.

### Temozolomide — reported, not targeted

TMZ was **not** tuned for, and nothing in either fix mentions it. It moves from 0.558 (2.86x, a 2.9x OVER-prediction of the measured 0.195) to 0.125 (0.64x) under FIX 1 alone, and to 0.159 (0.82x) with the leak. **The reshaping improves it**, from 2.86x error to 1.22x. This is a consequence of the veto (TMZ's logP -2.08 gives it a desirability of 0.027, which the arithmetic mean was averaging away), not a TMZ-specific adjustment.

## P3. FIX 2 — the regional test, re-run

Model side is now `R_region = clip(bbb_penetrance + 0.238 * D_region, 0, 1)` with `D = 1.00` for contrast-enhancing, `0.48` for necrosis, `0.00` for non-enhancing. Pre-fix every model ratio in this table was exactly **1.00x**.

**What is circular and what is not.** `L` was fitted to the CE-vs-non-enhancing *differences* of these same drugs, so the model reproducing a difference near 0.24 is a consistency check, not a prediction. The **ratios** are not circular: they depend on the base `bbb_penetrance` the leak is added to, and that was fitted to nothing.

| drug | contrast | measured ratio | model ratio (pre-fix) | **model ratio (post-fix)** | ratio error | measured diff | model diff | p (within-case) |
|---|---|---|---|---|---|---|---|---|
| levetiracetam | contrast_enhancing vs non_enhancing | 1.45x | 1.00x | **1.19x** | 0.82x | +0.272 | +0.160 | 0.0000 |
| levetiracetam | contrast_enhancing vs necrosis | 1.30x | 1.00x | **1.05x** | 0.81x | +0.201 | +0.046 | 0.0911 |
| levetiracetam | necrosis vs non_enhancing | 1.12x | 1.00x | **1.14x** | 1.02x | +0.072 | +0.114 | 0.3729 |
| cefazolin | contrast_enhancing vs non_enhancing | 4.26x | 1.00x | **5.87x** | 1.38x | +0.211 | +0.238 | 0.0000 |
| cefazolin | contrast_enhancing vs necrosis | 1.25x | 1.00x | **1.76x** | 1.40x | +0.056 | +0.124 | 0.6007 |
| cefazolin | necrosis vs non_enhancing | 3.40x | 1.00x | **3.34x** | 0.98x | +0.155 | +0.114 | 0.0166 |

Absolute per-region levels (the harder test — these are not differences, so the fitted `L` cannot rescue them):

| drug | region | n | measured mean | model | fold error |
|---|---|---|---|---|---|
| levetiracetam | contrast_enhancing | 39 | 0.875 | 1.000 | 1.14x |
| levetiracetam | necrosis | 17 | 0.674 | 0.954 | 1.42x |
| levetiracetam | non_enhancing | 69 | 0.602 | 0.840 | 1.39x |
| cefazolin | contrast_enhancing | 36 | 0.276 | 0.287 | 1.04x |
| cefazolin | necrosis | 15 | 0.220 | 0.163 | 0.74x |
| cefazolin | non_enhancing | 66 | 0.065 | 0.049 | 0.75x |

Independent literature replication (different patients, different assay, and `L` was fitted on only one of these three):

| drug | measured CE/non-enh | model CE/non-enh (post-fix) | measured diff | model diff | method | PMID |
|---|---|---|---|---|---|---|
| adavosertib | 5.65x | 1.52x | +3.440 | +0.238 | tissue_homogenate_LCMS | 39820282 |
| methotrexate | 4.65x | 11.57x | +0.230 | +0.238 | intracerebral_microdialysis | 18787762 |
| temozolomide | 1.39x | 2.90x | +0.074 | +0.238 | intracerebral_microdialysis | 19861433|29332604 |

### Honest verdict on FIX 2

**The model can now express the answer, and it lands the right order of magnitude — but it does not hit either anchor.**

- **levetiracetam CE/non-enhancing: measured 1.45x, model 1.19x — UNDERSHOOT (0.82x).** The cause is the ceiling, and it is diagnosable rather than mysterious: levetiracetam's base penetrance is 0.840, so 0.840 + 0.238 = 1.078 is clipped to 1.0 by the free-drug ceiling. The model delivers a difference of +0.160 where the measurement says +0.272. Note the measured CE mean is 0.875 and the max is 2.05 — a *homogenate* (total-drug) assay that legitimately exceeds 1, which assumption A1 structurally cannot. This is the same ceiling artefact that sinks adavosertib in section 6, not a new failure.
- **cefazolin CE/non-enhancing: measured 4.26x, model 5.87x — OVERSHOOT (1.38x).** Here the difference is right (+0.238 vs measured +0.211, 1.13x) and the ratio is inflated because the geometric-mean veto puts cefazolin's intact-barrier base at 0.049 against a measured non-enhancing mean of 0.065 — i.e. the *denominator* is 25% low. The regional error is inherited from the partition error, not from `L`.
- **methotrexate (independent, Blakeley 2008): measured 4.65x, model 11.57x — OVERSHOOT 2.5x.** Same mechanism as cefazolin and worse, because the veto drives methotrexate's base to 0.0225 against a measured non-enhancing 0.063. When the denominator is 3x too small, a correct numerator still gives a 2-3x wrong ratio.
- **Necrosis, the one contrast nothing was fitted to as a ratio.** Measured necrosis/non-enhancing is 1.12x (levetiracetam) and 3.40x (cefazolin); the model gives 1.14x (levetiracetam), 3.34x (cefazolin). The single `D = 0.48` reproduces both, which is the strongest evidence in this section that the additive form generalises — it was derived from the necrosis *differences* and it lands the *ratios* for two drugs 3x apart in penetrance.

Aggregate: geometric-mean absolute error on the two anchor CE/non-enhancing ratios is **1.30x**, against a pre-fix error of 2.49x (the pre-fix model predicted 1.00x for everything). The direction, the sign, the drug-dependence and the necrosis ordering are all now correct; the magnitudes are off by 1.2-1.4x on the anchors and ~2.5x on the independent methotrexate replication.

## P4. A third defect in the same file — `membrane_permeability` had no ceiling

Surfaced by the GIST imatinib validation (`results/gist_imatinib.md` item 7), fixed in this pass because it is the **same root pathology** as the eribulin veto bug: unbounded or unshaped use of a single descriptor. `membrane_permeability` computed `P = Kp * Dm / d` with `Kp = 10**logP` clamped only at logP 5, so it returned **4.23e+07 um/s for imatinib** — 42 metres per second.

Three bounds, all physical:

1. **logP driving `Kp` clamped at 3.0** (was 5.0). Bilayer/water partitioning plateaus above logP ~3 while octanol/water keeps climbing — the bilayer is a more ordered, more polar solvent than bulk octanol. Independent support: RDKit's Crippen logP for imatinib is 4.59 against an experimental ~3.0, i.e. the estimator runs ~1.6 log units lipophilic exactly where the clamp bites.
2. **Unstirred aqueous boundary layer in SERIES** — `1/P = 1/P_membrane + 1/P_aqueous` with `P_aqueous = D_water / h_UWL`, `h_UWL = 30 um`. A solute must cross the stagnant water film before it reaches the bilayer, so this is the *real* physical ceiling on passive permeability, and it is why measured Caco-2/PAMPA permeabilities saturate. 30 um is the THIN end of the reported 30-1500 um range, i.e. the most permissive ceiling the literature allows.
3. **Membrane viscosity 0.1 -> 1.0 Pa*s**, the top of the 0.1-1 Pa*s range the module already cited. At 0.1 the prefactor `Dm/d` alone was ~2.5e3 um/s, so every drug with `Kp >= 1` already exceeded any achievable permeability and the boundary layer would have flattened the whole library to one value.

| drug | logP | MW | BEFORE (um/s) | AFTER (um/s) | physically absurd before? |
|---|---|---|---|---|---|
| temozolomide | -2.08 | 194.2 | 12.41 | 1.19 | no |
| cisplatin | 0.20 | 298.0 | 2017 | 22.15 | **YES** |
| methotrexate | 0.27 | 454.4 | 2074 | 19.58 | **YES** |
| doxorubicin | 0.00 | 543.5 | 1056 | 17.07 | **YES** |
| lomustine | 2.25 | 233.7 | 2.489e+05 | 26.96 | **YES** |
| imatinib | 4.59 | 493.6 | 4.234e+07 | 21.03 | **YES** |
| bevacizumab | -1.00 | 149000.0 | 0 | 0 | no |

Maximum over the set after the fix: **27.0 um/s**. Nothing is absurd; the antibody stays a hard 0.

**Honest limitation — the ceiling binds, and the top of the ranking is now compressed.** Cisplatin (22.2), imatinib (21.0) and methotrexate (19.6 um/s) sit within 15% of each other because all three are boundary-layer limited under this parameterisation, whereas physically they differ by ORDERS of magnitude in measured apparent permeability — cisplatin is famously almost membrane-impermeant. Two residual defects cause this and neither is fixed here: the membrane term is still too permissive in absolute terms (real Caco-2 Papp for these drugs is 0.001-1 um/s, 1-4 orders below what this returns), and Crippen logP is simply wrong for some species (cisplatin +0.20 vs experimental ~-2.2). The output is a bounded, correctly-ordered-at-the-bottom RANKING, not a calibrated rate. That is a smaller error than 42 m/s, but it is not zero.

**A gap left deliberately open.** Both functions in `permeability.py` are neutral-species estimates with no pKa term, so neither can express pH-partitioning or acidic-compartment ion-trapping. Two independent validations now converge on that same gap: the GIST run found measured imatinib tumour:normal liver ratios of 0.023-0.142 lying BELOW the transport model's structural floor of 0.208 — unreachable at ANY parameterisation, with hepatic uptake transporters and lysosomal ion-trapping of a weak base named as the missing mechanisms — and the Moncayo explant run is independently testing the same effect (weak bases niraparib/rucaparib sequester, neutral olaparib does not). It is NOT implemented: there is no measured tumour acidic-phospholipid content in usable units and no tumour-specific lysosomal parameters in any published PBPK model, so it must be designed rather than parameterised. Neither fix above forecloses it — a neutral-fraction factor multiplies `Kp` and a compartment term composes downstream of the series resistance. The limitation is recorded in the module docstring so it is visible where the code is read.

**What FIX 2 does NOT fix.** It does not move the model's free-drug ceiling of 1.0, so homogenate ratios above 1 (adavosertib 4.18, levetiracetam CE max 2.05) remain unreachable — that needs the separate total-vs-free bookkeeping named in section 5 item 3. It does not model efflux, which is still the largest residual on eribulin. And the disruption FIELD is a documented proxy: `REGION_DISRUPTION` maps Ivy GAP transcriptional zones onto MRI enhancement by an anatomic argument (microvascular proliferation is what enhances; leading edge is what does not), and **no study has assayed Ivy GAP zone labels and drug partition in the same tissue**. The `core` zone weight of 0.50 is an unmeasured interpolation and is the weakest number in the module.


---

Generated by `tmesim/validation/glioma_regional_pk.py`. `core`, `physics` and `sim.engine` are untouched. All statistics use 20,000 seeded permutations.
