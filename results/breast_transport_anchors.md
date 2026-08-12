# Breast-tissue transport anchors — a literature survey, and one honest null

**Purpose.** GBM's transport physics has an external, measured anchor chain
(doxorubicin penetration fit → mitoxantrone predicted out-of-sample, see
`tmesim/validation/calibration.py` and `results/pk_reality.md`). Breast has
had no equivalent: every breast validation in this repo to date
(`results/kill_calibration_breast.md`) tests the 3-D:2-D **potency-ratio**,
never the **transport** side. This document surveys what measured
breast-tissue transport data actually exists, so that gap can be closed with
real numbers rather than borrowed ones.

Nothing here is wired into the engine yet. This is the evidence review that
has to come first.

**Provenance convention used below.** Literature retrieved via PubMed.
Entries marked **[verified here]** were read directly from the PubMed record
during this survey. Entries marked **[survey-reported]** came from a
delegated literature sweep and have NOT been independently re-read — they are
recorded as leads at a lower confidence tier, exactly like
`molecular/potency.py`'s literature-vs-measured tiering.

---

## 1. The headline result: a breast-specific tissue:plasma anchor DOES exist

Two independent human cohorts, two chemically different drugs, one number.
Both measured by **in vivo microdialysis directly in the primary breast
tumour interstitium** — i.e. the actual interstitial compartment this
simulator solves for, not a whole-tissue homogenate.

| drug | interstitial AUC as % of plasma AUC | n | PMID |
|---|---:|---:|---|
| 5-fluorouracil | **61 ± 11 %** (mean ± SE) | 10 | **9205062** [verified here] |
| methotrexate | **60 ± 20 %** (mean ± SE, CV = 100 %) | 9 | **9679959** [verified here] |

* Müller M *et al.* "5-fluorouracil kinetics in the interstitial tumor space:
  clinical response in breast cancer patients." *Cancer Res* 1997;57(13):2598-601.
  PMID 9205062.
* Müller M *et al.* "Interstitial methotrexate kinetics in primary breast
  cancer lesions." *Cancer Res* 1998;58(14):2982-5. PMID 9679959.

**Why this is a better anchor than a penetration depth would have been.** It
is a *ratio of two measured AUCs in humans*, unit-free, and it constrains
exactly the quantity this engine computes (`mean_tumor_auc` against the
plasma source term). Two unrelated drugs landing at 0.60 independently is a
convergence, not a coincidence.

**Two caveats that must travel with the number, both from the papers' own
abstracts:**

1. **The spread is enormous and is itself the finding.** MTX carries a
   coefficient of variation of 100 %, and the paper states plainly there is
   *no* correlation between plasma AUC and interstitial AUC (P = 0.93),
   concluding "plasma levels of MTX are not predictive of intratumor levels."
   A model that predicts tight plasma→interstitium coupling is falsified by
   this. The 0.60 is a population central value, not a per-patient constant.
2. **The two papers disagree about what the number MEANS clinically.** The
   5-FU paper found high interstitial AUC associated with increased tumour
   response; the MTX paper found interstitial load NOT associated with
   response, and concluded interstitial access "is not a rate-limiting step
   for clinical response" for MTX. Both can be true — delivery limits
   response for some drugs and not others — and that is precisely this
   project's own recurring "ranking is defensible, absolute response is not"
   boundary, arrived at independently by clinical pharmacologists.

## 2. Interstitial fluid pressure — the Jain effect, in human breast

The engine already models interstitial pressure. These are the human breast
numbers to test it against.

| tissue | IFP (mmHg) | n | PMID |
|---|---:|---:|---|
| primary breast carcinoma | **15 ± 9** (mean ± SD) | 17 patients / 21 tumours | **1423283** [verified here] |
| invasive ductal carcinoma | **29 ± 3** (mean ± SE) | 25 patients | **7850532** [verified here] |
| normal breast parenchyma | **−0.3 ± 0.1** | — | 7850532 [verified here] |
| benign breast tumours | 3.6 ± 0.8 | — | 7850532 [verified here] |
| **non-invasive (DCIS) carcinoma** | **−0.3 ± 0.2** | — | 7850532 [verified here] |
| normal tissue (all sites, same study) | 2 ± 4 | 11 | 1423283 [verified here] |

* Less JR, Posner MC, Boucher Y, Borochovitz D, Wolmark N, Jain RK.
  "Interstitial hypertension in human breast and colorectal tumors."
  *Cancer Res* 1992;52(22):6371-4. PMID 1423283.
* Nathanson SD, Nelson L. "Interstitial fluid pressure in breast cancer,
  benign breast conditions, and breast parenchyma." *Ann Surg Oncol*
  1994;1(4):333-8. PMID 7850532.

**The DCIS row is the most valuable line in this table.** Non-invasive
carcinoma shows *no* interstitial hypertension (−0.3 mmHg, statistically
indistinguishable from normal parenchyma), while invasive carcinoma in the
same study, same operating room, same technique reads 29 mmHg. That is a
built-in negative control: any mechanism this engine claims produces elevated
IFP must NOT produce it from tumour cellularity alone, since DCIS is a
tumour and has none. Few validation targets come with their own null case.

Note the two studies disagree by ~2x (15 vs 29 mmHg) on ostensibly the same
quantity — different technique (wick-in-needle vs intraoperative) and
different case mix. Use the range, not either point.

## 3. Penetration depth — the honest null, and a warning about the literature

**There is no breast-specific drug penetration half-distance of comparable
quality to GBM's doxorubicin 40-50 µm anchor (Primeau 2005, PMID 16361566).**
This was the original target of this survey and it was not found. Recorded as
a null rather than filled with a substitute.

The reason is not absence of data — it is that **three mutually
non-convertible quantities are all published under the word "penetration"**:

1. **Exponential half-distance / decay length** — Primeau 2005 (PMID 16361566),
   the format this repo's calibration actually uses.
2. **Thresholded mean penetration depth by image analysis** — the Jain-lab
   orthotopic breast measurements, e.g. Liu J *et al.* *PNAS*
   2012;109(41):16618-23, PMID 22996328 [survey-reported], giving doxorubicin
   ~9.6-45.1 µm across orthotopic MDA-MB-231/4T1.
3. **Intensity ratio at a fixed distance** — Tannock-lab measurements
   reporting doxorubicin still at 70-75 % of perivascular peak at 60 µm in
   MDA-MB-231/MCF-7 [survey-reported], which implies a decay length well over
   100 µm.

Taken at face value these say doxorubicin in breast penetrates 7 µm, 45 µm,
and >100 µm. They are not in conflict about biology; they are measuring
different things. **Do not blend them, and do not convert between them.**

Two further notes for whoever picks this up:

* Primeau 2005's own models (16C, EMT6) are **mammary** tumours, so the
  existing anchor is less GBM-specific than this repo has been treating it —
  it is better described as a cross-tissue small-molecule anchor that happens
  to have been validated here on glioma transport.
* **A number to actively avoid.** The figure "mammary carcinoma
  intercapillary distance 80-135 µm / normal breast 215 µm" circulates in the
  modelling literature. The survey traced its cited sources and found every
  one to be an MVD or vessel-area study, not an intercapillary-distance
  measurement — the value appears back-derived from ~1/sqrt(MVD) with no
  primary measurement behind it. It must not be adopted as measured data.
  Genuine human-breast intercapillary distance appears not to be published at
  all.

## 4. Secondary quantities (leads, not yet verified here)

All **[survey-reported]** — recorded so the trail is not lost, at the lower
confidence tier, and NOT to be used as anchors until read directly.

| quantity | value | source |
|---|---|---|
| extracellular water diffusivity, human breast tumour, in vivo | 1810 ± 460 µm²/s (Dex, IMPULSED MRI, n=7) | PMID 31765494 |
| whole-tissue ADC | malignant 1020 ± 170 vs normal breast 1780 ± 130 µm²/s | PMID 27242965 |
| microvessel density, human breast (CD34, 10-field mean) | 87 ± 46 /mm² (median 74, range 28-301) | PMID 9376265 |
| interstitial volume fraction | ~0.6 (from v_in = 38.9 ± 8.0 %) | PMID 31765494 |
| collagen area fraction, invasive ductal carcinoma | ~25-40 % across three studies | PMIDs 37729354, 30018725, 33791282 |
| tumour pO2, human breast (median) | 30 mmHg vs normal breast 65 mmHg | PMID 2040005 |
| paclitaxel tumour:plasma, human breast | 4-70x (total tissue, n=5) | PMID 24670687 |

Caveats carried from the survey: the MRI diffusivities measure **water, not
drug**; MVD numbers are hot-spot biased upward by design; the paclitaxel
ratio is **total** tissue concentration including tubulin-bound intracellular
drug, so it is not the same quantity as the §1 microdialysis free-interstitial
numbers and must not be pooled with them; and the pO2 figure has an
unresolved discrepancy (a later compilation from the same group is quoted at
10 mmHg rather than 30 for what appears to be the same cohort).

**Genuine nulls found** (searched, not present in accessible literature):
hydraulic conductivity of human breast tumour tissue; interstitial fluid
velocity in breast; intercapillary/intervessel distance in human breast; and
tissue:plasma partition for docetaxel, epirubicin or cisplatin in human
breast.

## 5. What this changes

The breast flagship's validation gap is narrower than it looked, and it is
a *different* gap than assumed. The transport side is not unanchorable — §1
gives a real, human, breast-specific, two-drug convergent target
(interstitial AUC ≈ 0.60 x plasma AUC), and §2 gives a pressure target with
its own negative control. What does NOT exist is the specific *penetration
half-distance* format GBM was validated in.

The natural next step is therefore **not** to keep hunting for a breast
penetration depth, but to test the engine's predicted mean tumour interstitial
AUC against the 0.60 ratio, with the CV = 100 % spread as the error bar and
the plasma-does-not-predict-interstitium finding as an explicit falsification
condition.

## 6. That test has now been run, and the model FAILED it

**See `results/breast_kp_validation.md`** (driver:
`experiments/breast_kp_anchor.py`). Summary of the outcome:

| drug | predicted | measured | verdict |
|---|---:|---|---|
| 5-fluorouracil | 7.7956 | 0.61 ± 0.11 | **12.8x too HIGH** |
| methotrexate | 0.0571 | 0.60 ± 0.20 | **10.5x too LOW** |

The two errors run in **opposite** directions, so no single constant fixes
them. The model spreads the panel over 970x where the measurement finds one
number. Two structural causes were diagnosed: a blood-**brain**-barrier
permeability applied to breast tissue, and a one-way delivery source with no
back-flux, which leaves the interstitial:plasma ratio unconstrained.

The anchor did its job on first contact. This is the value of §1 — not that
it confirmed the model, but that it did not.

**Both diagnosed causes were subsequently tested, and the residual is now
localised to a third term: the alpha-scaled binding sink
`k = clearance + binding`.** `results/transvascular_exchange_revalidation.md`
fixed the one-way source; `results/breast_tissue_permeability.md` replaced the
blood-brain-barrier permeability with a measured breast one and it made things
*worse*, which is what localised the defect — the PS each drug would need to
hit its own anchor differs by 4.58x while measured breast permeability differs
by 1.02x, and that 4.58x is inherited from `k`, which differs by 4.78x.
`results/binding_sink_refuted.md` then refuted that term in magnitude and
direction.

**This survey therefore has a sequel.** `k` has never been constrained by a
delivery-magnitude measurement in any tissue — `alpha` is fit on penetration
*depth* via `L = sqrt(D_eff/k)`, in which PS does not appear. The equivalent
literature search for `k`, in the same format as this document and with the
same honest nulls, is **`results/breast_binding_anchors.md`**. It reports four
new candidate anchors (three of them breast-specific), the direction the
evidence points, and the null that survives: no measured interstitial
loss-rate ratio exists for these two drugs in any tissue.
