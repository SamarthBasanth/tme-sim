# Does the compound ORDERING survive realistic per-compound plasma dosing?

**A sensitivity study on the shipping panels. No model change.** Nothing under
`tmesim/` was modified, no default moved, `calibrate_alpha()` is still
`22.223477711172563` (asserted against live code at every stage), and every arm
below is the unmodified engine driven with a per-run `PKConfig`. At the shipping
PK the arm IS the shipping arm, bit-for-bit.

* script: `experiments/pk_ranking_robustness.py`
  (`--stage linearity | gbm | breast | confirm | analyse | report`)
* frozen execution record: `results/sweep/pk_ranking_robustness.json` +
  `results/sweep/pk_ranking_robustness_cells.npz`
* the finding this answers: `results/axis_inventory.md` sec. 7 item 1 and sec. 9
  item 1, which ranked this the highest untested axis in the repo and flagged its
  own cost estimate as an unverified inference.
* ranking-stability conventions are `results/ranking_robustness.md`'s; sec. 4
  introduces one addition (the RESTORE FACTOR) and argues for it.

---

## 0. The answer

> ## NO.
>
> **A swap in the published ordering sits INSIDE A SINGLE COMPOUND'S OWN LICENSED
> DOSE RANGE — say this first.** Irinotecan's FDA label publishes Cmax at two
> doses in one table: 1,660 ng/mL at 125 mg/m2 (N=64) and 3,392 ng/mL at
> 340 mg/m2 (N=6), a **2.04x** span. The published GBM order's
> `irinotecan > temozolomide` relation needs only **1.961x** to reverse. One
> regulatory table contains both sides of the swap. No comparison across
> compounds, no assumption about regimens, and no argument from us is required.
>
> **The full published GBM ordering does not survive real dosing.**
> Shipping: `cisplatin` > `doxorubicin` > `etoposide` > `irinotecan` > `procarbazine` > `temozolomide`
> Real PK:  `cisplatin` > `etoposide` > `doxorubicin` > `temozolomide` > `irinotecan` > `procarbazine`
>
> Per-compound plasma exposure across the retained compounds spans
> **113.1x**. The order first breaks at a spread of
> **1.53x**.
>
> > **The published GBM ordering survives per-compound plasma exposure
> > differences of up to about 1.5-fold. Real human Cmax across the
> > compounds it ranks differs by 113-fold.**
>
> **THE GBM TOP-1 PICK HOLDS, and it is the first thing a per-compound axis
> has ever left standing in this repo.**
> `cisplatin` keeps rank 1 on the GBM panel in every arm, at
> every point of the spread family, with no breakpoint anywhere — including
> the unbound-drug arm and the adversarial arm. **It does NOT generalise: the
> breast top-1 moves** (below), so the defence is specific to this panel and
> must be stated as such.
>
> **The effect is ENTIRELY per-compound, and the common-mode control is null.**
> The geometric-mean exposure factor is **0.1120**
> — the panel is driven at roughly **8.9x**
> the exposure these compounds are actually given at. Applying that factor
> UNIFORMLY to every compound reproduces the published order **exactly**.
> Applying only the per-compound part, with the geometric mean forced to
> exactly 1, breaks the order **identically to the full real arm**.
> `results/axis_inventory.md` sec. 1.1 argued that the right test is whether a
> knob's ACTION factors through a per-compound quantity rather than whether
> the knob is per-compound. **This is that test executed with the common-mode
> control run and shown to be null.**
>
> **The adversarial control does not rescue it.** Re-running with every
> published dose choice that FAVOURS the published order applied at once —
> irinotecan: 340 mg/m2, 90-min IV (the label's OTHER published dose);
> doxorubicin: 75 mg/m2, EXTRAPOLATED from 20 mg/m2 assuming the label's dose-proportionality -- and 75 is ABOVE the 30-70 mg/m2 range the label validates, so this is extrapolation, flagged;
> — the order still moves, with 2 swaps surviving:
> `cisplatin` > `etoposide` > `doxorubicin` > `irinotecan` > `temozolomide` > `procarbazine`
>
> **Breast panel.** Top-1 **MOVES** (`vinorelbine` -> `cisplatin`); the full order **moves**.
> Shipping: `vinorelbine` > `cisplatin` > `paclitaxel` > `doxorubicin` > `docetaxel` > `epirubicin` > `5-fluorouracil` > `palbociclib` > `tamoxifen`
> Real PK:  `cisplatin` > `5-fluorouracil` > `epirubicin` > `paclitaxel` > `docetaxel` > `vinorelbine` > `doxorubicin` > `palbociclib` > `tamoxifen`
> 13 of 36
> pairs swap. **The breast result differs from the GBM one in a way that must
> be said rather than smoothed over: on breast the common-mode control is NOT
> null** — the geometric-mean factor alone moves the order. So the clean
> "entirely per-compound" decomposition is a GBM finding, not a general one.
> The reason is visible in the numbers: seven of the nine breast compounds sit
> within a shipping-`E[kill]` band of 0.041,
> so the published breast ordering is close to a tie before any perturbation.
>
> **The headline GBM five-drug ordering cannot be evaluated on this axis at all.**
> `results/ranking_robustness.md` scores `cisplatin > doxorubicin > lomustine >
> methotrexate > temozolomide`. **Two of those five compounds cannot be assigned
> a human plasma exposure** (sec. 6). An ordering that cannot in principle be
> dose-corrected is a different object from one that survives dose correction,
> and that is a property of the claim, not a limitation of this test.

**What changes on the marketing surface.**

* **C2 — *"the cross-drug ranking is invariant to common-mode perturbations"*
  — is FALSIFIED on breast, and this is the most consequential line in this
  document.** `results/axis_inventory.md` sec. 2 establishes that the surviving
  robustness claims C1, C2 and C4 are drawn **entirely** from the common-mode
  and numerical block. On the breast panel the geometric-mean exposure factor
  **alone** — the same number applied to every compound, which is exactly what
  "common-mode" means — moves the order (sec. 8.1). **This is therefore not
  one more per-compound axis failing. It is the first failure INSIDE the class
  the surviving claims rest on.** The inventory's own conclusion changes: *"the
  per-compound layer fails, the common-mode layer holds"* is no longer true as
  stated.
  **Scope, precisely, and do not generalise it in either direction:** *not
  falsified on GBM* (sec. 7.1, where the common-mode control is exactly null),
  *falsified on breast*. And the marketing clause enumerates four specific
  axes — cell-state route, oxygen route, alpha, grid — **none of which this
  result touches**. Those may still hold. What has gone is the GENERAL
  invariance, and with it the right to state C2 without naming a panel.
* **"The compound ranking is robust" cannot be said.** This is the SECOND
  independent axis in two days to break it, after
  `results/capacity_ranking_robustness.md`. Counting with
  `results/axis_inventory.md` sec. 2, that makes **seven per-compound axes ever
  scored against an ordering, and seven that moved it. There are still no
  exceptions.**
* **"Our GBM top pick is robust" survives this axis, and that is real.** It is
  the first per-compound axis to leave a claim standing anywhere in this corpus,
  and it should be said with the same emphasis as the failure. **It must be said
  panel-specifically**: the breast top-1 moves on the same axis (`vinorelbine` -> `cisplatin`), so "our top pick" is not a
  property of the method, it is a property of the GBM panel.
* **Every published GBM number in this repo was computed at roughly 9x supratherapeutic exposure.** It provably does not change the ordering (the
  common-mode control above), but it bears on every absolute quantity, and
  nothing anywhere in `results/` says it.

---

## 1. Step 1: the free-lunch claim, verified before anything was concluded

`results/axis_inventory.md` sec. 7 #1 stated, and correctly flagged as its own
unverified inference, that under the two shipping defaults a per-compound PK is
*exactly a per-compound multiplicative scale on each drug's AUC*, so a re-rank
costs no new PDE solves. **That inference is load-bearing and it is correct.**
Getting it wrong in the convenient direction would have invalidated everything
downstream, so it was established first, by reading the code AND numerically.

### 1.1 The algebra, from the code

* `PDEConfig.time_resolved` defaults to `False` (`core/config.py:87`), so
  `sim/engine.py:304-318` collapses the whole plasma profile to its window mean
  and builds ONE static source `S = vascular_density * bbb_penetrance *
  plasma_mean`.
* `treatment.drugs.delivery_source` is exactly linear in `plasma_conc` — one
  multiplication, no clamp, no saturation.
* The transport PDE is linear in `C`, and its implicit operator
  `A = I - dt*L + dt*K` does not involve the source at all, so `C` and therefore
  `AUC` are exactly proportional to `S`.
* The ONLY nonlinearity in the loop is the non-negativity clamp
  `np.maximum(c_new, 0.0)` at `physics/diffusion.py:455`.

**The repo already asserts this in two places and had never measured it.**
`validation/kill_calibration.py:841-843` says outright *"The PDE is linear in the
source and starts from C = 0, so `AUC(P) = P * AUC(1)` exactly"* and OPERATES on
that basis; `decompose`'s docstring repeats it. So this was a shipped, load-bearing
assumption with no measurement attached to it. The measurement below is the first.

### 1.2 The measurement, on real tissue at the SHIPPING grid

Real Xenium GBM section, dx = 12 um, dt = 10 s,
convection ON, 6 h window — the literal shipping configuration, not a cheaper
one. Each drug is solved at the shipping plasma level and again at a scaled
level; the scaled AUC field is compared with `f x` the baseline field.

| drug | plasma factor | max rel. dev. of AUC field | support identical? | E[kill] solved | E[kill] rescaled |
|---|---:|---:|:--:|---:|---:|
| cisplatin | 0.005 | 1.21e-14 | yes | 0.0349606075 | 0.0349606075 |
| cisplatin | 0.033 | 1.41e-14 | yes | 0.2086645097 | 0.2086645097 |
| cisplatin | 0.66 | 1.20e-14 | yes | 0.5165689406 | 0.5165689406 |
| cisplatin | 26 | 1.28e-14 | yes | 0.5411709037 | 0.5411709037 |
| temozolomide | 0.005 | 8.98e-15 | yes | 0.0002486128 | 0.0002486128 |
| temozolomide | 0.033 | 8.00e-15 | yes | 0.0023825806 | 0.0023825806 |
| temozolomide | 0.66 | 8.38e-15 | yes | 0.0740399947 | 0.0740399947 |
| temozolomide | 26 | 7.92e-15 | yes | 0.4996659123 | 0.4996659123 |

**Worst relative deviation anywhere: 1.41e-14.** That is
machine precision. The support of the AUC field is identical in every arm, so
the zero-clamp never fires. **Linearity holds; the re-rank is exact and costs
no new solves.**

### 1.3 Confirmation at the ACTUAL factors used (gbm panel)

The strongest available form of the check: every retained compound re-solved
at its OWN real PK, and the solved `E[kill]` compared with the rescaled one
that every number in this document is built from.

| drug | f | E[kill] solved | E[kill] rescaled | abs. diff |
|---|---:|---:|---:|---:|
| temozolomide | 0.48459 | 0.0535778492 | 0.0535778492 | 0.0e+00 |
| procarbazine | 0.0058527 | 0.0005016144 | 0.0005016144 | 1.1e-19 |
| cisplatin | 0.4877 | 0.5064839542 | 0.5064839542 | 0.0e+00 |
| methotrexate | 25.563 | 0.5005302416 | 0.5005302416 | 1.1e-16 |
| etoposide | 0.66178 | 0.2841335637 | 0.2841335637 | 5.6e-17 |
| irinotecan | 0.064959 | 0.0255796112 | 0.0255796112 | 6.9e-18 |
| doxorubicin | 0.033154 | 0.0788649422 | 0.0788649422 | 0.0e+00 |

Worst absolute difference **1.1e-16**. Had any row
disagreed, the corresponding result would have been withdrawn rather than
explained.

### 1.4 Confirmation at the ACTUAL factors used (breast panel)

The strongest available form of the check: every retained compound re-solved
at its OWN real PK, and the solved `E[kill]` compared with the rescaled one
that every number in this document is built from.

| drug | f | E[kill] solved | E[kill] rescaled | abs. diff |
|---|---:|---:|---:|---:|
| doxorubicin | 0.033154 | 0.1604203476 | 0.1604203476 | 2.8e-17 |
| paclitaxel | 0.12424 | 0.3608205566 | 0.3608205566 | 5.6e-17 |
| tamoxifen | 0.010407 | 0.0036559575 | 0.0036559575 | 2.2e-18 |
| docetaxel | 0.12735 | 0.3588562431 | 0.3588562431 | 5.6e-17 |
| epirubicin | 0.29414 | 0.3698064755 | 0.3698064755 | 1.1e-16 |
| 5-fluorouracil | 0.71176 | 0.4344334901 | 0.4344334901 | 0.0e+00 |
| vinorelbine | 0.056279 | 0.3501327393 | 0.3501327393 | 0.0e+00 |
| cisplatin | 0.4877 | 0.4570135911 | 0.4570135911 | 0.0e+00 |
| palbociclib | 0.0077657 | 0.0057457846 | 0.0057457846 | 0.0e+00 |

Worst absolute difference **1.1e-16**. Had any row
disagreed, the corresponding result would have been withdrawn rather than
explained.

### 1.5 What the cheapness of this test actually means

**The free lunch is a consequence of `time_resolved = False` specifically, and
that is worth more than the result it enables.** Under `time_resolved = True` a
per-compound half-life changes the SHAPE of the driver and not only its
amplitude, and the rescale stops being exact. A test that is cheap *because of a
shipped simplification* should say so, because the cheapness evaporates the
moment that simplification is revisited — and `tmesim/protocol/dosing.py:50`
already sets `time_resolved = True` for every TPL protocol, so the two halves of
this repo do not agree about which regime they are in.

---

## 2. FINDING: every published interstitial:plasma ratio in Codon is in an undeclared unit on its denominator

This has to be settled before a Cmax can be chosen, because published Cmax is
almost always TOTAL plasma drug while only UNBOUND drug crosses into tissue, and
the unbound fraction across an oncology panel spans a factor of ~90 (tamoxifen
<0.01 to temozolomide 0.85). That is a per-compound factor of the same order as
the axis under test, so driving the engine in the wrong unit would mean measuring
our own unit error rather than the model's PK sensitivity.

**Nothing in the codebase declares the convention.** `PKConfig.plasma_peak`'s
field doc (`core/config.py:59`) says only *"peak plasma concentration after a
dose, uM"*; `treatment.drugs.plasma_profile` says *"the vascular (blood) drug
concentration"*; `delivery_source` says *"the current blood concentration"*.
Total or free is stated in none of them.

**The one operational definition in the repo is TOTAL.**
`protocol/dosing.py:20` derives `plasma_peak_uM = 1000 x amount_mg x F / (Vd_L x
MW_g_per_mol)` — the one-compartment bolus relation `C0 = F*D/Vd`, which is a
total plasma concentration, with no free-fraction term anywhere in it.

**What the calibration anchors on, checked rather than assumed.**

| anchor | quantity fitted | carries a total/free convention? |
|---|---|---|
| `calibrate_alpha()` | doxorubicin penetration LENGTH (Primeau 2005 IF) | **no** — a length is invariant to the plasma level |
| `validation/ps_scale.py` | rodent brain permeability-surface-area products | **no** — a clearance |
| `results/breast_kp_validation.md` | interstitial:plasma AUC ratio vs microdialysis | **YES, and it is inherited and undeclared** |

`breast_kp_validation.md` sec. 1 is explicit that the engine's state variable is
the *free interstitial* concentration — *"exactly what a microdialysis probe
samples — not total tissue homogenate, which would include intracellular and
protein-bound drug"* — and then divides by `plasma_mean`. The denominator's
convention comes from Müller 1997 (PMID 9205062) and Müller 1998 (PMID 9679959)
and is stated nowhere in this repo.

> **Every published interstitial:plasma ratio in Codon is in an undeclared unit
> on its denominator.**

This document therefore reports **two clearly-labelled arms and never folds them
together**: **Arm T** at published TOTAL Cmax, matching `dosing.py`'s own
derivation and what the shipped 50 uM nominal is; and **Arm F** at unbound Cmax
= total Cmax x a separately sourced measured `fu`, so a reader can see whether a
swap is a DOSING result or a BINDING result.

---

## 3. FINDING: `free_fraction` is computed for every compound and dropped on the floor

`tmesim/molecular/binding.py`'s own module docstring, line 6, lists
`plasma_free_fraction` as **"the unbound fraction driving delivery"**.
`compile_drug.py:217` computes it for every compiled molecule, `types.py:77`
stores it on `MolecularProfile.free_fraction`, and `gist_imatinib.py:1049` prints
it in a published report.

**`compile_drug.molecule_to_drugparams` (lines 362-371) does not pass it into
`DrugParams`, and `DrugParams` has no such field.** Asserted at runtime: the
dataclass fields are `name, diffusivity, molecular_weight, ic50, hill, binding,
clearance, bbb_penetrance, hypoxia_resistance`. So the quantity the repo itself
names as driving delivery never reaches the engine. This is the same shape as
`formal_charge` in `results/axis_inventory.md` sec. 5.2: computed, documented as
load-bearing, consumed nowhere.

**The useful half: there is therefore no double-count risk.** Nothing anywhere in
the delivery path corrects for protein binding, so Arm F is the FIRST such
correction rather than a second one, and it is interpretable as stated.

### 3.1 Scoring the heuristic against measurement — the first time it has been scored

`plasma_free_fraction` is a logistic in logP, `ff = 1/(1 + 10^((logP-3)/1.5))`,
tagged **CONFIDENCE: low** and never scored against data. Having sourced
measured protein binding for this panel, it costs nothing to score it (n = 12).

| drug | logP | heuristic ff | measured fu | fold error | fu source tier |
|---|---:|---:|---:|---:|---|
| etoposide | 1.339 | 0.9276 | 0.0300 | **30.92x** | FDA label |
| docetaxel | 3.260 | 0.4015 | 0.0300 | **13.38x** | US PI (in three cancer patients) |
| epirubicin | 0.001 | 0.9901 | 0.2300 | **4.30x** | FDA label |
| doxorubicin | 0.001 | 0.9901 | 0.2500 | **3.96x** | FDA label |
| paclitaxel | 3.736 | 0.2442 | 0.0650 | **3.76x** | FDA label (midpoint of a published 89-98% range) |
| palbociclib | 2.966 | 0.5130 | 0.1500 | **3.42x** | FDA label |
| irinotecan | 4.091 | 0.1578 | 0.5100 | **3.23x** | FDA label (midpoint of a published range 30-68% bound) |
| methotrexate | 0.268 | 0.9851 | 0.5000 | **1.97x** | FDA label |
| vinorelbine | 4.754 | 0.0634 | 0.1100 | **1.73x** | FDA label (a directly measured free fraction) |
| temozolomide | -2.078 | 0.9996 | 0.8500 | **1.18x** | FDA label |
| 5-fluorouracil | -0.798 | 0.9971 | 0.9000 | **1.11x** | TERTIARY -- BC Cancer Drug Manual only |
| tamoxifen | 5.996 | 0.0100 | 0.0100 | **1.00x** | MHRA SmPC (a '>99%' bound statement, so fu is an upper bound) |

**Geometric-mean fold error 3.23x, worst 30.9x**, Spearman rank correlation against the measured
values **rho = 0.650**.

The module's docstring claims *"What we CAN defend is the direction and rough
magnitude of the trends"*. On this panel the direction is **defensible** and the magnitude is not: the heuristic is not merely imprecise, it puts
etoposide (measured fu 0.03, 97 % bound) at ff = 0.93. It is fortunate that
nothing consumes it (sec. 3) — but a value that is both wrong and unused is a
worse artefact than one that is only unused, because it is printed in reports.

---

## 4. The RESTORE FACTOR — how a swap is reported here

A swap reported as a bare fact ("A and B exchange places") tells a reader
nothing about how close the call was, and quantifying it the usual way would mean
inventing a PK uncertainty distribution, which this study refuses to do.

**Construction.** For a pair where the published order says `A > B` and the real
PK says `B > A`, the restore factor is the multiplicative scale on `A`'s
PK-adjusted per-cell exposure that makes `E[kill](A) = E[kill](B)` again, found
by bisection. `E[kill]` is strictly increasing in the exposure scale, so it is
unique. Multiplying it by `A`'s published Cmax converts it into a statement in
units a pharmacologist can check:

> **the published order requires this compound's Cmax to be at least X uM.**

It is distribution-free, needs no PK uncertainty model, and — the part that
matters — **it is falsifiable by a reader who knows more pharmacology than we
do.** A clinical pharmacologist reads "doxorubicin would need 5.90 uM" and knows
immediately whether that is reachable. No probability this project could quote
has that property. **It supersedes reporting a swap as a bare fact and should be
how ranking sensitivity is reported here from now on.**

---

## 5. Where the numbers come from, and the rules declared before the data

Four literature agents sourced human plasma Cmax, half-life and protein binding
under a strict rule: **verbatim quote plus PMID / label document, or NOT
SOURCED.** No value was estimated, dose-scaled, or recalled. The full quotes and
sources for every number are in `PK_TABLE` and `PK_EXCLUDED` in the script.

* **R1.** Retain a compound only if a peak plasma concentration OF THE COMPOUND
  THE ENGINE SIMULATES has been measured in adult humans at a stated clinical
  dose, in a retrievable primary source or regulatory label, with a verbatim quote.
* **R2.** Where a source publishes a RANGE, both endpoints are carried and the
  headline uses the geometric midpoint, labelled as an aggregation of a published
  range and never as a measured value.
* **R3.** Nothing is estimated, extrapolated or dose-scaled to fill a gap. A
  compound with no sourceable Cmax is EXCLUDED and the exclusion is a result.
* **R4.** Cmax is converted to uM with the molecular weight THE ENGINE ITSELF
  compiled, so the conversion cannot drift from the object being dosed. Platinum
  is converted on a platinum-atom basis (MW 195.084).

### 5.1 The rule earned its keep: a fabricated citation was caught and rejected

A web-search layer returned, attributed to the FDA lomustine label, the sentence
*"the parent drug lomustine has not been detected in blood, but the peak
concentrations of the ring hydroxylated metabolites are approximately 3 uM after
doses of 130 mg/m2"*. The sourcing agent downloaded the 2013, 2014 and 2016 FDA
lomustine labels and searched all three: **the text is not present in any of
them.** The underlying facts happen to be close to Lee 1985's, but the citation
is fabricated. It was rejected. This is the strongest available evidence that the
verbatim-quote-plus-PMID rule is doing work rather than adding ceremony, and it is
recorded here rather than quietly dropped.

---

## 6. The compounds that CANNOT be assigned a human exposure

These are results, not gaps.

### 6.1 `lomustine` — EXCLUDED under R1

The parent compound has never been detected in human plasma after oral dosing. This is a MEASURABILITY exclusion, not a mechanism claim -- see sec. 6.6, which states exactly what the evidence supports and no more.

> Parent CCNU was not detected in the plasma of any of the patients, probably due to complete conversion to monohydroxylated metabolites during the 'first pass' through liver and gut.

*Lee FY, Workman P, Roberts JT, Bleehen NM, Cancer Chemother Pharmacol 1985;14(2):125-31, PMID 3971475 (n=4, 130 mg/m2 oral)*

### 6.2 `carmustine` — EXCLUDED under R1

No Cmax at a standard 150-200 mg/m2 dose exists, and no plasma elimination half-life has ever been published. The only PK study is at 300-750 mg/m2 dose-normalised to 1 g/m2, and reports a rate constant rather than a half-life. The current BiCNU label's sec. 12.3 contains a Distribution subsection and no elimination data at all.

> Intravenously administered carmustine is rapidly degraded, with no intact drug detectable after 15 minutes. [2011 label; FDA DELETED this sentence in the 2017 revision, which makes no measurability claim at all.]

*Henner WD et al., Cancer Treat Rep 1986;70(7):877-80, PMID 3719578; FDA BiCNU label 2011 and 2017 (017422s038, 017422s055)*

### 6.3 `cyclophosphamide` — EXCLUDED under R1

No observed Cmax for the parent compound at a standard 600 mg/m2 dose exists in any label or retrievable primary source. The best available number is a BACK-EXTRAPOLATED ZERO-TIME INTERCEPT C0 = 254.4 uM at 1000 mg/m2, which is not an observed peak. Five studies dose at exactly 600 mg/m2 and none reports a Cmax in its abstract. Separately, cyclophosphamide is a PRODRUG, so the parent's Cmax may not be the exposure that matters even if it were measured.

> Table 1 Relevant pharmacokinetic parameters of CP ... C0 (uM) | 254.4 (44.98) | 232.0 (133.4) | 6

*FDA Cyclophosphamide Injection label (210852s000) sec. 12.3 -- states no Cmax; Chan KK et al., Cancer Res 1994;54(24):6421-9, PMID 7987837, Table 1*

### 6.4 `vincristine` — EXCLUDED under R1

No numeric adult peak plasma/serum concentration is stated in the US label or in any primary PK study retrieved. Sethi 1981 reports alpha/beta/gamma half-lives, Vc, Vdss, clearance and AUC, and no Cmax. A '100-400 nM peak' figure circulating on secondary aggregator pages was checked and is UNCONFIRMED at primary source.

> The initial, middle and terminal half-lives are 5 minutes, 2.3 hours, and 85 hours respectively [FDA label -- no Cmax is given]

*FDA Vincristine Sulfate Injection label rev. 07/2013; Sethi VS et al., Cancer Res 1981;41(9 Pt 1):3551-5, PMID 7260915; Nelson 1982 PMID 7070351; Jackson 1981 PMID 7296550*

### 6.5 What this does to the headline five-drug GBM ordering

`results/ranking_robustness.md` scores `cisplatin > doxorubicin > lomustine >
methotrexate > temozolomide` through 10,240 re-rankings. **`lomustine` cannot be
assigned a human plasma exposure at all**, and `methotrexate`'s clinical exposure
is a regimen rather than a number (sec. 9). So:

> **The headline-5 ordering cannot be evaluated on this axis, because 40 % of its
> members have no measurable parent plasma concentration or no single one.**

The model currently drives, at a nominal 50 uM plasma level, a compound that has
never been detected in human plasma.

### 6.6 What the nitrosourea evidence does and does not support

It is tempting to conclude that simulating the transport of parent lomustine and
carmustine is category-wrong. **The evidence does not support that, and it is not
claimed here.**

* **Lomustine: the label is silent, leaning parent.** Current GLEOSTINE and all
  three generic PLR labels carry byte-identical sec. 12.1 — *"Lomustine alkylates
  DNA and RNA."* Every occurrence of `metabolit*` in the full SPL was checked: it
  appears three times (lactation, geriatric excretion, half-life) and **never in
  an activity-attribution sentence.**
* **Carmustine: metabolite-attributed, but WEAKENED over time.** The 2011 BiCNU
  label reads *"the antineoplastic and toxic activities of carmustine may be DUE
  TO metabolites"*; the 2017 and all six current SPLs read *"The metabolites may
  CONTRIBUTE TO antitumor activity."* The revision moved AWAY from sole-metabolite
  attribution, and FDA also deleted the *"no intact drug detectable after 15
  minutes"* sentence, so the current label makes no measurability claim at all.
* **"Parent vs metabolite" is a false dichotomy.** The 4-hydroxy-CCNUs are
  themselves intact nitrosoureas — hydroxylation is on the cyclohexyl ring and the
  N-nitroso group is untouched — and both parent and metabolites decompose to the
  same class of alkylating species. Kastrissios 1996 (PMID 8765435) distinguishes
  alkylating moieties formed FROM THE PARENT from those formed FROM THE METABOLITE.
  The parent is not inert; it is one of two precursors feeding the same pool.
* **"Parent unmeasurable" is route- and species-specific, not intrinsic.** In mice
  given CCNU i.p. the parent is readily measurable. Human non-detectability is
  attributed to first pass after ORAL dosing.

**The claim the evidence does support is a calibration claim, not a mechanism
claim:**

> **The object the model transports is not the object measured in human plasma,
> and the two differ materially in the property transport is most sensitive to.**

Parent CCNU is ~90-96 % protein bound; the major circulating metabolite is
~49-53 % bound (Lee & Workman 1983, PMID 6849803, Table IV, in vitro and in vivo)
— **roughly a ten-fold difference in free fraction between what is modelled and
what circulates**, landing directly on Arm F. Note this also contradicts the
tertiary "lomustine 50 % bound" figure in the BC Cancer monograph, which appears
to describe the metabolite rather than the parent.

**Kept as a gap:** Wheeler 1977 (PMID 145865), the definitive head-to-head on the
metabolites' own antitumour properties, is paywalled with no abstract. Zero
content was retrieved. It is recorded as unsourced, not as support.

**The exclusions stand on measurability grounds regardless**, which is a
validation problem independent of mechanism, and always the stronger argument.

---

## 7. GBM panel

Real Xenium GBM section, 40,887 measured cells
(23,812 tumour), grid (424, 675),
dx = 12 um, dt = 10 s, convection ON,
6 h window, sensitivity mode `uniform`. **This is the shipping configuration of
`experiments/real_drug_panel.py`, not a cheaper one.** One solve per compound; the
kill-kernel identity against the frozen `cells.pkpd.apply_drug_kill` is
**0.0e+00**.

*Cross-check.* The shipping-arm `E[kill]` values below reproduce
`results/capacity_ranking_robustness.md`'s `R = 1` column to 4 decimal places,
although that study ran at dx = 15 um with convection OFF. Two grids and both
convection settings agree on the baseline.

### 7.1 Arm T (published TOTAL plasma Cmax) — the base arm

Retained: `temozolomide`, `procarbazine`, `cisplatin`, `etoposide`, `irinotecan`, `doxorubicin`

| drug | Cmax used (uM) | t1/2 (h) | dose the Cmax is at | `f` | E[kill] shipping | E[kill] real PK |
|---|---:|---:|---|---:|---:|---:|
| **cisplatin** | 15.43 | 120 | 100 mg/m2, 6 h IV infusion | 0.4877 | 0.526025 | 0.506484 |
| **doxorubicin** | 1.102 | 30.8 | 20 mg/m2, 30-min IV infusion (AIDS-related Kaposi sarcoma regimen) | 0.033154 | 0.478390 | 0.078865 |
| **etoposide** | 27.69 | 6.63 | 100 mg/m2 IV, day 1 of CHOEP | 0.66178 | 0.344598 | 0.284134 |
| **irinotecan** | 2.829 | 5.8 | 125 mg/m2, 90-min IV infusion | 0.064959 | 0.293498 | 0.025580 |
| **procarbazine** | 1.098 | 0.69 | 185 +/- 4 mg/m2/day oral, day 1 (+EIASD cohort) | 0.0058527 | 0.157135 | 0.000502 |
| **temozolomide** | 38.63 | 1.8 | 150 mg/m2, single oral dose | 0.48459 | 0.111307 | 0.053578 |

* shipping order: `cisplatin` > `doxorubicin` > `etoposide` > `irinotecan` > `procarbazine` > `temozolomide`
* real-PK order:  `cisplatin` > `etoposide` > `doxorubicin` > `temozolomide` > `irinotecan` > `procarbazine`
* top-1 **held**; full order **MOVED**; per-compound exposure spread **113.1x**
* geometric-mean exposure factor **0.1120** (the panel is driven at ~8.9x real exposure)
* COMMON-MODE control (every compound at the geometric mean): order **reproduced EXACTLY**
* PER-COMPOUND-only arm (geometric mean pinned to 1): order **MOVED**, and it is the SAME order as the full real arm
* **the order first breaks at a per-compound spread of 1.53x** (family exponent lambda = 0.09)
* **top-1 has no breakpoint anywhere on the spread family**
* a PURE COMMON multiplier does not break the order anywhere over 1e-4x to 1e3x — the 'a constant error cancels out of a ranking' defence is VALID for this axis, and is simply not what this axis is

**Every swap, with the restore factor (sec. 4).**

| pair | shipping winner | real-PK winner | shipping margin | new margin | restore factor on the loser | Cmax the published order requires |
|---|---|---|---:|---:|---:|---:|
| `temozolomide` vs `procarbazine` | `procarbazine` | `temozolomide` | 0.0458 | 0.0531 | 54.802x | 60.17 uM (published 1.098 uM) |
| `temozolomide` vs `irinotecan` | `irinotecan` | `temozolomide` | 0.1822 | 0.0280 | 1.961x | 5.549 uM (published 2.829 uM) |
| `etoposide` vs `doxorubicin` | `doxorubicin` | `etoposide` | 0.1338 | 0.2053 | 5.356x | 5.903 uM (published 1.102 uM) |

### 7.2 Arm T with high-dose methotrexate added

Methotrexate's only sourceable adult-scale Cmax is at the high-dose regime, in a children-and-young-adults cohort on target-concentration-guided dosing (Nagamine 2025, PMID 41372843). It is reported as its own arm and never merged with the others. **The conclusion holds either way** — the point of showing it is that a result which survives its own most contested row is stronger than one that quietly drops it.

Retained: `temozolomide`, `procarbazine`, `cisplatin`, `methotrexate`, `etoposide`, `irinotecan`, `doxorubicin`

| drug | Cmax used (uM) | t1/2 (h) | dose the Cmax is at | `f` | E[kill] shipping | E[kill] real PK |
|---|---:|---:|---|---:|---:|---:|
| **cisplatin** | 15.43 | 120 | 100 mg/m2, 6 h IV infusion | 0.4877 | 0.526025 | 0.506484 |
| **doxorubicin** | 1.102 | 30.8 | 20 mg/m2, 30-min IV infusion (AIDS-related Kaposi sarcoma regimen) | 0.033154 | 0.478390 | 0.078865 |
| **etoposide** | 27.69 | 6.63 | 100 mg/m2 IV, day 1 of CHOEP | 0.66178 | 0.344598 | 0.284134 |
| **irinotecan** | 2.829 | 5.8 | 125 mg/m2, 90-min IV infusion | 0.064959 | 0.293498 | 0.025580 |
| **procarbazine** | 1.098 | 0.69 | 185 +/- 4 mg/m2/day oral, day 1 (+EIASD cohort) | 0.0058527 | 0.157135 | 0.000502 |
| **methotrexate** | 955.2 | 11 | median total 15.32 g/m2, 1 h loading + 5 h maintenance infusion | 25.563 | 0.129749 | 0.500530 |
| **temozolomide** | 38.63 | 1.8 | 150 mg/m2, single oral dose | 0.48459 | 0.111307 | 0.053578 |

* shipping order: `cisplatin` > `doxorubicin` > `etoposide` > `irinotecan` > `procarbazine` > `methotrexate` > `temozolomide`
* real-PK order:  `cisplatin` > `methotrexate` > `etoposide` > `doxorubicin` > `temozolomide` > `irinotecan` > `procarbazine`
* top-1 **held**; full order **MOVED**; per-compound exposure spread **4367.7x**
* geometric-mean exposure factor **0.2433** (the panel is driven at ~4.1x real exposure)
* COMMON-MODE control (every compound at the geometric mean): order **reproduced EXACTLY**
* PER-COMPOUND-only arm (geometric mean pinned to 1): order **MOVED**, and it is the SAME order as the full real arm
* **the order first breaks at a per-compound spread of 1.29x** (family exponent lambda = 0.03)
* top-1 first breaks at a per-compound spread of 7222.14x
* a PURE COMMON multiplier does not break the order anywhere over 1e-4x to 1e3x — the 'a constant error cancels out of a ranking' defence is VALID for this axis, and is simply not what this axis is

**Every swap, with the restore factor (sec. 4).**

| pair | shipping winner | real-PK winner | shipping margin | new margin | restore factor on the loser | Cmax the published order requires |
|---|---|---|---:|---:|---:|---:|
| `temozolomide` vs `procarbazine` | `procarbazine` | `temozolomide` | 0.0458 | 0.0531 | 54.802x | 60.17 uM (published 1.098 uM) |
| `temozolomide` vs `irinotecan` | `irinotecan` | `temozolomide` | 0.1822 | 0.0280 | 1.961x | 5.549 uM (published 2.829 uM) |
| `procarbazine` vs `methotrexate` | `procarbazine` | `methotrexate` | 0.0274 | 0.5000 | 3297.155x | 3620 uM (published 1.098 uM) |
| `methotrexate` vs `etoposide` | `etoposide` | `methotrexate` | 0.2148 | 0.2164 | 8.298x | 229.8 uM (published 27.69 uM) |
| `methotrexate` vs `irinotecan` | `irinotecan` | `methotrexate` | 0.1637 | 0.4750 | 115.175x | 325.9 uM (published 2.829 uM) |
| `methotrexate` vs `doxorubicin` | `doxorubicin` | `methotrexate` | 0.3486 | 0.4217 | 45.709x | 50.37 uM (published 1.102 uM) |
| `etoposide` vs `doxorubicin` | `doxorubicin` | `etoposide` | 0.1338 | 0.2053 | 5.356x | 5.903 uM (published 1.102 uM) |

### 7.3 Arm F (unbound Cmax = total x measured fu)

A SECOND per-compound axis, reported separately so a reader can tell a dosing result from a binding result. Cisplatin uses MEASURED ultrafilterable platinum rather than `Cmax x fu`, because platinum-albumin adducts are irreversible and a reversible unbound fraction is the wrong object for it (the label says so outright). Procarbazine is omitted from this arm: its plasma protein binding is NOT SOURCED — the BC Cancer monograph's PK table literally reads "no information found".

Retained: `temozolomide`, `cisplatin`, `etoposide`, `irinotecan`, `doxorubicin`  ·  omitted from this arm: `procarbazine`

| drug | Cmax used (uM) | t1/2 (h) | dose the Cmax is at | `f` | E[kill] shipping | E[kill] real PK |
|---|---:|---:|---|---:|---:|---:|
| **cisplatin** | 2.102 | 0.663 | 100 mg/m2, 6 h IV infusion | 0.010775 | 0.526025 | 0.078957 |
| **doxorubicin** | 0.2755 | 30.8 | 20 mg/m2, 30-min IV infusion (AIDS-related Kaposi sarcoma regimen) | 0.0082886 | 0.478390 | 0.017835 |
| **etoposide** | 0.8308 | 6.63 | 100 mg/m2 IV, day 1 of CHOEP | 0.019853 | 0.344598 | 0.010348 |
| **irinotecan** | 1.443 | 5.8 | 125 mg/m2, 90-min IV infusion | 0.033129 | 0.293498 | 0.011787 |
| **temozolomide** | 32.83 | 1.8 | 150 mg/m2, single oral dose | 0.4119 | 0.111307 | 0.044949 |

* shipping order: `cisplatin` > `doxorubicin` > `etoposide` > `irinotecan` > `temozolomide`
* real-PK order:  `cisplatin` > `temozolomide` > `doxorubicin` > `irinotecan` > `etoposide`
* top-1 **held**; full order **MOVED**; per-compound exposure spread **49.7x**
* geometric-mean exposure factor **0.0300** (the panel is driven at ~33.4x real exposure)
* COMMON-MODE control (every compound at the geometric mean): order **reproduced EXACTLY**
* PER-COMPOUND-only arm (geometric mean pinned to 1): order **MOVED**, and it is the SAME order as the full real arm
* **the order first breaks at a per-compound spread of 7.47x** (family exponent lambda = 0.515)
* top-1 first breaks at a per-compound spread of 85.86x
* a PURE COMMON multiplier does not break the order anywhere over 1e-4x to 1e3x — the 'a constant error cancels out of a ranking' defence is VALID for this axis, and is simply not what this axis is

**Every swap, with the restore factor (sec. 4).**

| pair | shipping winner | real-PK winner | shipping margin | new margin | restore factor on the loser | Cmax the published order requires |
|---|---|---|---:|---:|---:|---:|
| `temozolomide` vs `etoposide` | `etoposide` | `temozolomide` | 0.2333 | 0.0346 | 3.677x | 3.055 uM (published 0.8308 uM) |
| `temozolomide` vs `irinotecan` | `irinotecan` | `temozolomide` | 0.1822 | 0.0332 | 3.263x | 4.709 uM (published 1.443 uM) |
| `temozolomide` vs `doxorubicin` | `doxorubicin` | `temozolomide` | 0.3671 | 0.0271 | 2.306x | 0.6353 uM (published 0.2755 uM) |
| `etoposide` vs `irinotecan` | `etoposide` | `irinotecan` | 0.0511 | 0.0014 | 1.118x | 0.9291 uM (published 0.8308 uM) |

### 7.4 Dose audit — which factors may be rescaled, and which may not

Compounds are being compared at whatever dose their own source published,
which is not necessarily a coherent clinical regimen across the panel. That is
a real threat to the comparison's meaning, so it is made visible rather than
buried: the dose each factor is derived at is printed in every table above,
and every compound whose sources publish MORE THAN ONE dose is re-run at the
alternative here.

| compound | alternative published exposure | Cmax (uM) | resulting order | swaps |
|---|---|---:|---|---:|
| `irinotecan` | 340 mg/m2, 90-min IV (the label's OTHER published dose) | 5.782 | `cisplatin` > `etoposide` > `doxorubicin` > `irinotecan` > `temozolomide` > `procarbazine` | 2 |
| `doxorubicin` | 75 mg/m2, EXTRAPOLATED from 20 mg/m2 assuming the label's dose-proportionality -- and 75 is ABOVE the 30-70 mg/m2 range the label validates, so this is extrapolation, flagged | 4.133 | `cisplatin` > `etoposide` > `doxorubicin` > `temozolomide` > `irinotecan` > `procarbazine` | 3 |

**Dose-proportionality, per each compound's own source — this table is the
most reusable thing here, because it says which factors a future study may
rescale and which it may not.**

| compound | may a Cmax be linearly rescaled? | evidence |
|---|---|---|
| etoposide | **yes**, over 100-600 mg/m2 | label: *"the areas under the plasma concentration vs. time curves (AUC) and the maximum plasma concentration (Cmax) values increase linearly with dose"* |
| epirubicin | **yes**, over 60-150 mg/m2 (but the 60 and 75 mg/m2 Cmax are non-monotonic at n=6) | ELLENCE label: *"Epirubicin pharmacokinetics are linear over the dose range of 60 to 150 mg/m2"* |
| palbociclib | **yes**, over 25-225 mg | ClinPharm review: *"the AUC and Cmax increased proportionally with dose"* |
| doxorubicin | **only within 30-70 mg/m2**, on n=4 — and the sourced value is at 20 mg/m2, BELOW that range, so scaling it is extrapolation | label: *"In four patients, doxorubicin demonstrated dose-independent pharmacokinetics across a dose range of 30 mg/m2 to 70 mg/m2"* |
| irinotecan | **no need to** — the label publishes two measured doses | CAMPTOSAR Table 9, 125 and 340 mg/m2 |
| paclitaxel | **no** — Cmax is NOT dose-proportional and infusion duration dominates | label: *"with a 3-hour infusion, for a 30% increase in dose, the Cmax and AUC were increased by 68% and 89%"*; 175 mg/m2 gives 3650 ng/mL over 3 h and 365 ng/mL over 24 h |
| vincristine | **no** — explicitly non-proportional, clearance rises with dose | (and it is excluded anyway, sec. 6.3) |
| 5-fluorouracil | **no** — bolus and infusional regimens differ ~120x in circulating level | Casale 2004 vs Joel 2004 |
| cisplatin | **no** — free-platinum Cmax is dominated by infusion duration at the same 100 mg/m2 | Vermorken 1986: *"Peak levels were both dose and schedule dependent"* |

### 7.5 The ADVERSARIAL control — the arm built to save the published order

This is **not** the default arm. Every published dose choice that FAVOURS
the published ordering is applied at once, including one that is explicit
extrapolation past the top of a validated range:

* irinotecan: 340 mg/m2, 90-min IV (the label's OTHER published dose)
* doxorubicin: 75 mg/m2, EXTRAPOLATED from 20 mg/m2 assuming the label's dose-proportionality -- and 75 is ABOVE the 30-70 mg/m2 range the label validates, so this is extrapolation, flagged

* shipping order: `cisplatin` > `doxorubicin` > `etoposide` > `irinotecan` > `procarbazine` > `temozolomide`
* **best case for the published claim:** `cisplatin` > `etoposide` > `doxorubicin` > `irinotecan` > `temozolomide` > `procarbazine`
* order **MOVED**, 2 swap(s) surviving; top-1 **held**

**The ordering breaks in every arm that can be constructed, including the
one constructed to save it.** That removes every objection about dose
choice at once, and it is a stronger result than the base case.

---

## 8. Breast panel

### 8.0 The breast ordering is not separated enough to support a ranking claim in the first place

**Read this before the swaps, because the swaps are a symptom and this is the
diagnosis.**

> Seven of the nine retained breast compounds sit inside a shipping-`E[kill]`
> band of **0.041** (0.4967 down to 0.4555).

The restore factors (sec. 4) make that concrete. These published pairwise
calls are decided by less than 15 % in exposure:

| published call | restore factor on the loser |
|---|---:|
| `vinorelbine` > `docetaxel` | **1.079x** |
| `vinorelbine` > `paclitaxel` | **1.099x** |
| `paclitaxel` > `epirubicin` | **1.099x** |
| `docetaxel` > `epirubicin` | **1.131x** |

A reader who accepts every physics assumption in this model should still not
act on those calls. **Almost anything moves the breast order because there is
almost nothing to move.** That is a different and more fundamental problem
than "a perturbation reorders the panel", and it means the 13
swaps below are evidence of the separation problem rather than of a PK effect
specifically. The PK effect is real and is quantified below; the separation
problem would remain even if it were not.

Real Janesick Xenium breast section (`rep1`),
167,780 measured cells
(89,785 tumour), grid (458, 628),
dx = 12 um, dt = 10 s,
convection ON, 6 h window — the shipping configuration of
`experiments/breast_drug_panel.py`. Kill-kernel identity
**0.0e+00**.

### 8.1 Arm T (published TOTAL plasma Cmax)

Retained: `doxorubicin`, `paclitaxel`, `tamoxifen`, `docetaxel`, `epirubicin`, `5-fluorouracil`, `vinorelbine`, `cisplatin`, `palbociclib`

| drug | Cmax used (uM) | t1/2 (h) | dose the Cmax is at | `f` | E[kill] shipping | E[kill] real PK |
|---|---:|---:|---|---:|---:|---:|
| **vinorelbine** | 1.857 | 34.8 | 30 mg/m2 as a 10-min infusion (Navelbine reference formulation) | 0.056279 | 0.496740 | 0.350133 |
| **cisplatin** | 15.43 | 120 | 100 mg/m2, 6 h IV infusion | 0.4877 | 0.483454 | 0.457014 |
| **paclitaxel** | 4.274 | 20.2 | 175 mg/m2 IV over 3 h | 0.12424 | 0.479046 | 0.360821 |
| **doxorubicin** | 1.102 | 30.8 | 20 mg/m2, 30-min IV infusion (AIDS-related Kaposi sarcoma regimen) | 0.033154 | 0.478724 | 0.160420 |
| **docetaxel** | 4.58 | 13.9 | 100 mg/m2 as a one-hour infusion | 0.12735 | 0.473344 | 0.358856 |
| **epirubicin** | 9.751 | 32.1 | 75 mg/m2 IV over 6-10 min | 0.29414 | 0.459880 | 0.369806 |
| **5-fluorouracil** | 426.2 | 0.215 | 400 mg/m2 by intravenous BOLUS injection | 0.71176 | 0.455494 | 0.434433 |
| **palbociclib** | 0.2592 | 29 | 125 mg orally once daily, 21-days-on / 7-days-off | 0.0077657 | 0.386580 | 0.005746 |
| **tamoxifen** | 0.3284 | 142 | 20 mg orally once daily for three months | 0.010407 | 0.297272 | 0.003656 |

* shipping order: `vinorelbine` > `cisplatin` > `paclitaxel` > `doxorubicin` > `docetaxel` > `epirubicin` > `5-fluorouracil` > `palbociclib` > `tamoxifen`
* real-PK order:  `cisplatin` > `5-fluorouracil` > `epirubicin` > `paclitaxel` > `docetaxel` > `vinorelbine` > `doxorubicin` > `palbociclib` > `tamoxifen`
* top-1 **MOVED**; full order **MOVED**; per-compound exposure spread **91.7x**
* geometric-mean exposure factor **0.0855** (the panel is driven at ~11.7x real exposure)
* COMMON-MODE control (every compound at the geometric mean): order **MOVED**
* PER-COMPOUND-only arm (geometric mean pinned to 1): order **MOVED**
* **the order first breaks at a per-compound spread of 1.50x** (family exponent lambda = 0.09)
* top-1 first breaks at a per-compound spread of 5.32x
* a pure common multiplier breaks the order at 0.0001x

**Every swap, with the restore factor (sec. 4).**

| pair | shipping winner | real-PK winner | shipping margin | new margin | restore factor on the loser | Cmax the published order requires |
|---|---|---|---:|---:|---:|---:|
| `doxorubicin` vs `docetaxel` | `doxorubicin` | `docetaxel` | 0.0054 | 0.1984 | 4.761x | 5.246 uM (published 1.102 uM) |
| `doxorubicin` vs `epirubicin` | `doxorubicin` | `epirubicin` | 0.0188 | 0.2094 | 5.273x | 5.812 uM (published 1.102 uM) |
| `doxorubicin` vs `5-fluorouracil` | `doxorubicin` | `5-fluorouracil` | 0.0232 | 0.2740 | 11.289x | 12.44 uM (published 1.102 uM) |
| `paclitaxel` vs `epirubicin` | `paclitaxel` | `epirubicin` | 0.0192 | 0.0090 | 1.099x | 4.697 uM (published 4.274 uM) |
| `paclitaxel` vs `5-fluorouracil` | `paclitaxel` | `5-fluorouracil` | 0.0236 | 0.0736 | 2.593x | 11.08 uM (published 4.274 uM) |
| `paclitaxel` vs `vinorelbine` | `vinorelbine` | `paclitaxel` | 0.0177 | 0.0107 | 1.099x | 2.04 uM (published 1.857 uM) |
| `docetaxel` vs `epirubicin` | `docetaxel` | `epirubicin` | 0.0135 | 0.0110 | 1.131x | 5.179 uM (published 4.58 uM) |
| `docetaxel` vs `5-fluorouracil` | `docetaxel` | `5-fluorouracil` | 0.0179 | 0.0756 | 2.865x | 13.12 uM (published 4.58 uM) |
| `docetaxel` vs `vinorelbine` | `vinorelbine` | `docetaxel` | 0.0234 | 0.0087 | 1.079x | 2.004 uM (published 1.857 uM) |
| `epirubicin` vs `5-fluorouracil` | `epirubicin` | `5-fluorouracil` | 0.0044 | 0.0646 | 2.141x | 20.88 uM (published 9.751 uM) |
| `epirubicin` vs `vinorelbine` | `vinorelbine` | `epirubicin` | 0.0369 | 0.0197 | 1.193x | 2.215 uM (published 1.857 uM) |
| `5-fluorouracil` vs `vinorelbine` | `vinorelbine` | `5-fluorouracil` | 0.0412 | 0.0843 | 2.506x | 4.653 uM (published 1.857 uM) |
| `vinorelbine` vs `cisplatin` | `vinorelbine` | `cisplatin` | 0.0133 | 0.1069 | 3.705x | 6.88 uM (published 1.857 uM) |

### 8.2 Dose / schedule audit

| compound | alternative published exposure | Cmax (uM) | resulting order | swaps |
|---|---|---:|---|---:|
| `doxorubicin` | 75 mg/m2, EXTRAPOLATED from 20 mg/m2 assuming the label's dose-proportionality -- and 75 is ABOVE the 30-70 mg/m2 range the label validates, so this is extrapolation, flagged | 4.133 | `cisplatin` > `5-fluorouracil` > `epirubicin` > `paclitaxel` > `docetaxel` > `vinorelbine` > `doxorubicin` > `palbociclib` > `tamoxifen` | 13 |
| `paclitaxel` | 175 mg/m2 over 24 h instead of 3 h (SAME dose, same table) | 0.4274 | `cisplatin` > `5-fluorouracil` > `epirubicin` > `docetaxel` > `vinorelbine` > `doxorubicin` > `paclitaxel` > `palbociclib` > `tamoxifen` | 14 |
| `5-fluorouracil` | infusional de Gramont: 2400 mg/m2 over 46 h, steady-state 0.46 +/- 0.08 mg/L -- NOT a Cmax, and ~120x below the bolus | 3.536 | `cisplatin` > `epirubicin` > `paclitaxel` > `docetaxel` > `vinorelbine` > `doxorubicin` > `5-fluorouracil` > `palbociclib` > `tamoxifen` | 8 |

### 8.3 Two breast compounds are dosed in a way the engine cannot represent

`tamoxifen` and `palbociclib` are ORAL and chronically dosed, and reach steady
state. A single-bolus `PKConfig` cannot express that. The closest available
objects are their published STEADY-STATE concentrations (tamoxifen 122 ng/mL
after 20 mg daily for three months; palbociclib 116 ng/mL at steady state on
125 mg daily), used with their long half-lives so the 6 h window sees a
near-constant level — which is the right shape, but it is a substitution and
is labelled as one. This is the same problem as an oral lomustine Cmax, in a
milder form. Tamoxifen additionally has an active metabolite,
N-desmethyl tamoxifen, circulating ~3x higher than the parent and not modelled.

### 8.4 Arm F (unbound Cmax)

Tamoxifen's `fu` is derived from a *">99 % bound"* statement, so 0.01 is an UPPER BOUND on its unbound fraction, not a measurement. 5-fluorouracil's `fu` is TERTIARY-source-only (BC Cancer Drug Manual; no FDA or MHRA label states a percentage). Both are flagged in the table rather than silently used.

Retained: `doxorubicin`, `paclitaxel`, `tamoxifen`, `docetaxel`, `epirubicin`, `5-fluorouracil`, `vinorelbine`, `cisplatin`, `palbociclib`

| drug | Cmax used (uM) | t1/2 (h) | dose the Cmax is at | `f` | E[kill] shipping | E[kill] real PK |
|---|---:|---:|---|---:|---:|---:|
| **vinorelbine** | 0.2042 | 34.8 | 30 mg/m2 as a 10-min infusion (Navelbine reference formulation) | 0.0061907 | 0.496740 | 0.083944 |
| **cisplatin** | 2.102 | 0.663 | 100 mg/m2, 6 h IV infusion | 0.010775 | 0.483454 | 0.053694 |
| **paclitaxel** | 0.2778 | 20.2 | 175 mg/m2 IV over 3 h | 0.0080757 | 0.479046 | 0.069376 |
| **doxorubicin** | 0.2755 | 30.8 | 20 mg/m2, 30-min IV infusion (AIDS-related Kaposi sarcoma regimen) | 0.0082886 | 0.478724 | 0.044803 |
| **docetaxel** | 0.1374 | 13.9 | 100 mg/m2 as a one-hour infusion | 0.0038205 | 0.473344 | 0.035852 |
| **epirubicin** | 2.243 | 32.1 | 75 mg/m2 IV over 6-10 min | 0.067651 | 0.459880 | 0.184243 |
| **5-fluorouracil** | 383.6 | 0.215 | 400 mg/m2 by intravenous BOLUS injection | 0.64058 | 0.455494 | 0.426585 |
| **palbociclib** | 0.03888 | 29 | 125 mg orally once daily, 21-days-on / 7-days-off | 0.0011648 | 0.386580 | 0.000597 |
| **tamoxifen** | 0.003284 | 142 | 20 mg orally once daily for three months | 0.00010407 | 0.297272 | 0.000015 |

* shipping order: `vinorelbine` > `cisplatin` > `paclitaxel` > `doxorubicin` > `docetaxel` > `epirubicin` > `5-fluorouracil` > `palbociclib` > `tamoxifen`
* real-PK order:  `5-fluorouracil` > `epirubicin` > `vinorelbine` > `paclitaxel` > `cisplatin` > `doxorubicin` > `docetaxel` > `palbociclib` > `tamoxifen`
* top-1 **MOVED**; full order **MOVED**; per-compound exposure spread **6155.1x**
* geometric-mean exposure factor **0.0076** (the panel is driven at ~130.7x real exposure)
* COMMON-MODE control (every compound at the geometric mean): order **MOVED**
* PER-COMPOUND-only arm (geometric mean pinned to 1): order **MOVED**
* **the order first breaks at a per-compound spread of 1.42x** (family exponent lambda = 0.04)
* top-1 first breaks at a per-compound spread of 18.59x
* a pure common multiplier breaks the order at 0.0001x

**Every swap, with the restore factor (sec. 4).**

| pair | shipping winner | real-PK winner | shipping margin | new margin | restore factor on the loser | Cmax the published order requires |
|---|---|---|---:|---:|---:|---:|
| `doxorubicin` vs `epirubicin` | `doxorubicin` | `epirubicin` | 0.0188 | 0.1394 | 4.851x | 1.337 uM (published 0.2755 uM) |
| `doxorubicin` vs `5-fluorouracil` | `doxorubicin` | `5-fluorouracil` | 0.0232 | 0.3818 | 40.256x | 11.09 uM (published 0.2755 uM) |
| `paclitaxel` vs `epirubicin` | `paclitaxel` | `epirubicin` | 0.0192 | 0.1149 | 3.370x | 0.9364 uM (published 0.2778 uM) |
| `paclitaxel` vs `5-fluorouracil` | `paclitaxel` | `5-fluorouracil` | 0.0236 | 0.3572 | 35.021x | 9.73 uM (published 0.2778 uM) |
| `paclitaxel` vs `cisplatin` | `cisplatin` | `paclitaxel` | 0.0044 | 0.0157 | 1.282x | 2.694 uM (published 2.102 uM) |
| `docetaxel` vs `epirubicin` | `docetaxel` | `epirubicin` | 0.0135 | 0.1484 | 6.830x | 0.9385 uM (published 0.1374 uM) |
| `docetaxel` vs `5-fluorouracil` | `docetaxel` | `5-fluorouracil` | 0.0179 | 0.3907 | 82.871x | 11.39 uM (published 0.1374 uM) |
| `epirubicin` vs `5-fluorouracil` | `epirubicin` | `5-fluorouracil` | 0.0044 | 0.2423 | 8.298x | 18.61 uM (published 2.243 uM) |
| `epirubicin` vs `vinorelbine` | `vinorelbine` | `epirubicin` | 0.0369 | 0.1003 | 2.569x | 0.5247 uM (published 0.2042 uM) |
| `epirubicin` vs `cisplatin` | `cisplatin` | `epirubicin` | 0.0236 | 0.1305 | 3.889x | 8.172 uM (published 2.102 uM) |
| `5-fluorouracil` vs `vinorelbine` | `vinorelbine` | `5-fluorouracil` | 0.0412 | 0.3426 | 20.368x | 4.16 uM (published 0.2042 uM) |
| `5-fluorouracil` vs `cisplatin` | `cisplatin` | `5-fluorouracil` | 0.0280 | 0.3729 | 28.010x | 58.87 uM (published 2.102 uM) |

---

## 9. FINDING: the engine cannot express what a plasma curve is

Two independent expressivity limits surfaced while doing this, and neither is
visible in any current output.

**9.1 The window-mean collapse.** With `time_resolved = False` the plasma profile
is replaced by its arithmetic mean over the exposure window before the solver
sees it. Therefore **the engine cannot distinguish a short-half-life, high-Cmax
compound from a long-half-life, low-Cmax one at matched window mean.** Two
compounds with completely different plasma kinetics become bit-identical inputs.
`results/dose_schedule.md` records the SCHEDULE half of this consequence; the
HALF-LIFE half — that Cmax and t1/2 are not separately identifiable at all — is
recorded nowhere. On this panel half-life enters only through the bounded factor
`(1 - exp(-lambda T))/(lambda T)`, which compresses a 5-day terminal half-life
and a 1.8 h one into a factor of about 2.5.

**9.2 A compound's real exposure is a regimen, not a number.** Methotrexate spans
roughly three orders of magnitude between low-dose weekly and 12 g/m2 with
leucovorin rescue. 5-fluorouracil spans ~120x between bolus and infusional
de Gramont. Paclitaxel spans 10x between a 3 h and a 24 h infusion of the SAME
dose. Methotrexate's Cmax is comparable across oral and subcutaneous routes while
the route difference lands in AUC — which is exactly the distinction a
window-mean driver cannot express. **`PKConfig` has one `plasma_peak` per
compound, so the model asks a question — "which drug is better?" — that has no
regimen-independent answer for several of the compounds on its own panels.**

---

## 10. The asymmetry: which published claims this touches, and which it cannot

`results/axis_inventory.md` sec. 7 #1 asserted that the 3-D:2-D kill verdicts are
immune to this axis. **Verified rather than inherited.**
`validation/kill_calibration.decompose` (lines 585-625) takes `e_unit`, the
per-cell exposure **per uM of plasma**, and sets
`plasma = c_bath / percentile(e_unit, 99)`. The plasma level is chosen so the
best-perfused tumour cells sit exactly at the bath concentration. **The scored
quantity therefore never reads `PKConfig` at all** — it is not merely robust to a
per-drug plasma scale, it is constructed independently of the plasma level.

| claim | exposed to per-compound PK? | why |
|---|:--:|---|
| K1/K2/K3, the 3-D:2-D kill verdicts | **NO** | per-drug p99 bath-matching normalises the plasma scale away by construction |
| penetration length `L`, interstitial:plasma ratio, converged gradient shape | **NO** | steady-state / equilibrium-referenced; a source scale cancels |
| C2, *"cross-drug ranking is invariant to common-mode perturbations"* | **not falsified on GBM; FALSIFIED on breast** | the common-mode control is exactly null on GBM (sec. 7.1) — C2 is true there and simply does not cover this axis. On breast the geometric-mean factor ALONE moves the order (sec. 8.1), so C2 does not hold panel-agnostically |
| the full compound ORDERING on both panels | **YES, and it moves** | sec. 7, sec. 8 |
| the top-1 pick, GBM | **YES, and it holds** | sec. 7.1 |
| the top-1 pick, breast | **YES** | sec. 8.1 |
| the dose-limited-vs-resistance-limited reading in `dose_schedule.md` sec. 6.3 | **YES** | it is stated at a nominal 1x that is ~9x supratherapeutic |

**A correction to the record, stated as a correction.** An earlier reading of the
GBM result inside this same investigation treated *"top-1 holds"* as a defence of
the METHOD. It is not; it is a property of the GBM panel, and the breast panel
contradicts it on the same axis with the same code. The same over-reading is
available for the common-mode decomposition in sec. 7.1, which is why sec. 8.1's
contradiction of it is reported at the top rather than in a limitations list.

**That asymmetry is the point.** The same input is harmless to one published claim
and fatal to another, and until this document a reader could not tell which. The
scored 3-D:2-D panel has a defence the ranking panel does not, and nothing on any
surface said so.

---

## 11. What is NOT claimed here

* **This is not evidence that the engine's PK model is wrong in form.** A
  one-compartment bolus with a window mean is a defensible Stage-0 choice. The
  finding is narrower and harder to escape: **the ranking claim is conditional on
  every compound sharing a plasma curve, which no two of them do.**
* **No PK uncertainty distribution is invented.** Every statement about how close
  a swap is takes the restore-factor form (sec. 4), which needs none.
* **Arm F is not a correction the engine should adopt.** It is a second,
  separately reported axis. Whether unbound drug is the right driver depends on a
  convention the repo has never declared (sec. 2), and declaring it is a
  re-validation exercise, not a flag flip.
* **Convection is ON in every arm** (unlike
  `results/capacity_ranking_robustness.md`, which had to run it off). Nothing
  about the PDE operator changes between arms here, so no compromise was needed.
* **One tissue section per panel.** `results/heldout_patient.md` shows absolutes
  move a lot between sections; this study does not re-run the second section.
* **`n_doses = 1` throughout**, matching the shipping default. Multi-dose courses
  are a separate axis, swept in `results/dose_schedule.md`.

---

## 12. Reproduction

```
./.venv/bin/python experiments/pk_ranking_robustness.py --stage solves
./.venv/bin/python experiments/pk_ranking_robustness.py --stage linearity
./.venv/bin/python experiments/pk_ranking_robustness.py --stage confirm
./.venv/bin/python experiments/pk_ranking_robustness.py --stage report
```

`calibrate_alpha()` asserted live at every stage = **22.223477711172563**. `PKConfig` defaults
asserted: `plasma_peak = 50.0` uM, `half_life = 4` h, `n_doses = 1`.
The shipping window-mean plasma level is
**31.094993 uM**, which is the denominator of every `f`
in this document.

`experiments/real_drug_panel.py`, `experiments/breast_drug_panel.py` and
`tmesim/validation/ranking_robustness.py` were re-checked at run time and none
passes a `pk=` argument (`_assert_panels_pass_no_pk`), so the finding this study
exists for is still true of the live sources.

*No sweep changed a default and nothing under `tmesim/` was modified. Test suite
not executed by this pass and not expected to move: 1171 passed / 0 failed.*
