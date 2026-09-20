# Packing density as an external anchor on `φ` — Grantab & Tannock 2012

**Date: 2026-09-20. Route R4 of `results/verification_routes_2026-09-20.md`.**

`calibrate_alpha()` at START = **22.223477711172563** (asserted equal to
22.223477711172563 before any of the work below).
`calibrate_alpha()` at END = **22.223477711172563** (asserted again after it).
`two_route.DEFAULT_CELL_VOLUME_FRACTION` = 0.7892, unchanged, asserted at both ends.

---

## HEADLINE — **FAIL**, on the pre-registered dispersion criterion P3

> *This block was written after the run. Section 0 below is the pre-registration and
> has not been edited since it was written to disk before anything was computed.*

The shipped cell-volume-fraction term `φ = 0.7892` was scored, for the first time,
against an external measurement: two independent perturbations of cellular packing
density in multicellular cultures, with matched within-drug penetration changes
(Grantab & Tannock 2012, PMC3407510). Eight within-drug ratios, two mechanistically
independent perturbations, two sublines, a laboratory and cell line the law was never
fit to.

**Three pre-registered criteria. Two pass, one fails, and the pre-registered rule is
that all three must pass.**

| criterion | measured | bar | result |
|---|---|---|---|
| **P1** accuracy — `GMFE(engine, measured)` | **1.2678** | ≤ 2.0 | **PASS** *(demoted; see §0.1)* |
| **P2** anti-vacuity — engine vs the φ-blind null | **1.2678** vs null **2.1878** | ≤ 1.7503 | **PASS** |
| **P3** dispersion — `ln(S_pred)/ln(S_meas)` | **0.1962** | ≥ 0.50 | **FAIL** |

**What failed, precisely.** The engine gets the *central tendency* of this experiment
right and is not a vacuous pass — it beats "packing density does nothing" by 42% in
geometric mean fold error, which is a real result and is reported as one in §4.2. What
it cannot do is produce the *spread*. The measured ratios range over **2.2572×**; the
engine's predictions range over **1.1732×**. It accounts for **19.6%** of the log
dynamic range of the phenomenon it is being credited with explaining.

**The sharpest form of the failure** is in §4.4: within a *single* cell line under a
*single* perturbation — one packing-density change, for which the engine necessarily
emits one number — the measured responses differ by up to **1.9053×** where the engine
predicts a spread of **1.0037×**. The law is a function of `φ` alone, so this is a
structural limit, not a calibration error. Packing density is not a sufficient
explanation of these measurements.

**Three escape routes were checked and all three are closed.** (i) The unreported
chamber geometry cannot rescue it: the transient arm reaches the largest measured ratio
only at a slab thickness of **3200 µm**, a 3.2 mm multicellular culture (§4.5).
(ii) No choice of packing mapping rescues it — the alternative mapping makes the spread
*smaller*, 1.047× (§5.2). (iii) Dropping the weakest rows does not change it: P3 fails
at 0.1918 with all doxorubicin excluded, and at 0.0775 and 0.0133 with the two
perturbations scored separately (§5).

**What survives, and it is not nothing.** Inverting the law per row gives the
`φ_control` that would reproduce each measurement. All eight lie inside the physically
admissible band, spanning **0.7242–0.9293**, geometric mean **0.8318**. **The shipped
0.7892 sits inside that range.** This is the first external constraint of any kind on a
number derived from this repository's own tissue, and it does not refute it (§4.3). It
is also a loose constraint — a 1.28× span — and §6 says why it is weaker than it looks.

**The single most important caveat, stated in the headline rather than buried.** The
mapping from the paper's "% nuclear area" to the engine's cell volume fraction is
**not a verified identity of quantities**, and it is the weakest link in this document.
It is argued in §0.4 before any number was computed and re-examined in §6.1. The P3
failure is robust to it; the P1/P2 passes are not.

---

## 0. Pre-registration

**Status at the time this section was written: PRE-REGISTRATION ONLY. Nothing had been computed.**

**This section was written to disk before the experiment script was run, and has not
been edited since.** Everything below §1 was appended afterwards. The rule this obeys
is `results/transport_anchor_prereg.md` §0 and the repository's standing instruction
that a pre-registered threshold may never be moved after the residual is visible — a
rule that exists because a previous session declared a ≥2× bar, measured 1.56×, and
substituted an unrelated 1.74 floor.

### 0.1 Full disclosure of what was already known when the thresholds were set

Pre-registration is worthless if the analyst has already seen the answer. I had, in
part, and I am recording exactly how much rather than claiming a purity I do not have:

* I had fetched the paper and extracted every number in §0.3 below.
* I had computed, by hand and then in a scratch script, the **impermeable-limit
  hindrance ratios** `h(φ_after)/h(φ_before)` for the four perturbations under the
  mapping in §0.4. I knew they were in the neighbourhood of 1.7–2.0 before I fixed the
  thresholds.
* I had run the three compounds through `two_route_diffusivity` and established that
  all three sit at the impermeable limit (§0.5), so the predicted ratio is
  compound-independent.
* I had **not** computed any GMFE, any spread statistic, any implied `φ`, or any
  transient solve.

Because criterion **P1** below could therefore have been chosen to fit, P1 is set to
the field-standard value (the 2-fold rule) rather than to anything derived from these
numbers, and **P1 is explicitly demoted: it is not the criterion this document turns
on.** The criteria that decide the verdict are **P2** and **P3**, which are comparisons
against a null and against the measurement's own dynamic range, and I had computed
neither the null nor the spread when I set them.

### 0.2 What is under test

`DEFAULT_CELL_VOLUME_FRACTION = 0.7892` (`tmesim/molecular/two_route.py:125`). It enters
every compiled compound's effective diffusivity through the two-route law

```
hindrance = (1 + 2·β·φ) / (1 − β·φ)          (Maxwell / Hasselman–Johnson)
β         = (σ_c* − σ_m) / (σ_c* + 2·σ_m)
```

and in the impermeable-cell limit (`σ_c* → 0`, `β → −1/2`) reduces exactly to

```
hindrance_impermeable = (1 − φ) / (1 + φ/2)
```

The value 0.7892 was obtained by Delesse's principle on **this repository's own Xenium
breast section**. That makes it self-consistent, not verified. Its own source comment
concedes it is *"An UPPER estimate — Xenium's expansion segmentation tiles the section.
Literature tumour band implies 0.40–0.72"*, i.e. the shipped value sits **above** the
band the repo itself cites. It has never been scored against an external measurement.

### 0.3 The measurement, extracted by me from the primary text

Grantab RH, Tannock IF. *Penetration of anticancer drugs through tumour tissue as a
function of cellular packing density and interstitial fluid pressure and its
modification by bortezomib.* BMC Cancer 2012;12:214. PMID 22672469, PMC3407510,
DOI 10.1186/1471-2407-12-214. Open access.

Fetched 2026-09-20 from `https://pmc.ncbi.nlm.nih.gov/articles/PMC3407510/`
(HTTP 200, 195,531 bytes). Every figure below was read out of that document by me, with
the sentence it came from recorded. Nothing here is taken from the scouting survey's
transcription.

**(a) Packing density under EGTA, with the sentence:**

> "As determined by clonogenic assays, exposure of MCC to 10 mM EGTA for 90 min caused
> no toxicity to either cell line; at this dose mean packing density (+/-SD) was reduced
> from 50 ± 7.3% to 40 ± 1.6% and from 40 ± 6.8% to 31 ± 4.5%, in MCC derived from
> HCT-8Ea and HCT-8E11 cells respectively (p = 2.2 × 10-9 and p = 8.1 × 10-6 for HCT-8Ea
> and HCT-8E11 cells respectively)."

**(b) That the change is a bulk change, not a surface artefact** — load-bearing, because
a surface-only change would not be a change in bulk `φ` at all:

> "The changes in packing density were uniform throughout the MCC and not limited to the
> MCC surface (p = 0.78 and p = 0.82 for HCT-8Ea and HCT-8E11 cells respectively)."

**(c) Packing density under bortezomib** — Table 1, "MCC Packing Density (% nuclear
area)", 24 h exposure, control vs 1.0 μM:

| bortezomib (μM) | HCT-8Ea | HCT-8E11 |
|---|---|---|
| 0 | 53.2 ± 5.5 | 50.8 ± 3.8 |
| 0.25 | 50.9 ± 6.2 | 48.6 ± 3.2 |
| 0.50 | 50.2 ± 5.6 | 46.6 ± 4.1 |
| **1.0** | **44.6 ± 5.7\*** | **42.5 ± 3.8\*** |
| 2.5 | 45.9 ± 5.6\* | MCC disaggregated |

with "\* p≤0.05 for comparison with control conditions" and

> "Packing density changes were shown to be significant after pre-treatment with 1 μM of
> bortezomib (p = 1.3 × 10-7 and p = 1.8 × 10-6 for HCT-8Ea and HCT-8E11 derived MCC
> respectively as assessed using a t-test)."

**(d) Penetration change under EGTA** — same drug, same dose, same MCC:

> "EGTA pre-treatment improved penetration of 5-fluorouracil and gemcitabine (measured at
> 6 h) by ~3.3-4.4 fold in MCC derived from HCT-8Ea cells and by ~1.5-1.9 fold in MCC
> derived from HCT-8E11 cells (Figure 1). EGTA also led to ~2-fold increase in
> penetration of doxorubicin (p = 0.02 for HCT-8E11, p = 0.002 for HCT-8Ea; data not
> shown)."

**(e) Penetration change under bortezomib:**

> "Bortezomib pre-treatment improved penetration of 5-fluorouracil by ~1.7-fold and
> gemcitabine by ~3-fold in MCC derived from HCT-8Ea cells and increased penetration of
> both drugs by ~2-fold in MCC derived from HCT-8E11 cells at 6 h (Figure 1).
> Penetration of doxorubicin was also increased by about 1.8-fold at 6 h (p = 0.02 for
> HCT-8Ea, p = 0.005 for HCT-8E11, data not shown)."

**(f) What the reported fold-change actually is** — this is the single most important
methods sentence for the comparison, because it fixes the observable:

> "Drug penetration at 6 h was calculated as the ratio of the concentration of drug in
> compartment 2 in treated MCL to the concentration of the drug in compartment 2 in
> non-treated MCL."

So the reported fold is `C₂(6 h)_treated / C₂(6 h)_control` — a dimensionless,
within-drug, same-dose ratio. Dose asymmetry, detectability and absolute calibration
all cancel out of it identically. This is what makes the route safe from the
cross-compound ranking withdrawal of 2026-08-21.

**(g) The chamber, and the internal standard:**

> "A volume of 0.5 mL was added to one side of the MCC (compartment 1), and initial drug
> concentrations (using a combination of radio-labelled and unlabelled drug) were: 10 μM
> doxorubicin, 77 μM 5-fluorouracil and 100 μM gemcitabine."

> "…permit sensitive detection of the drug in compartment 2, which contained 18 mL of
> stirred culture media."

> "[14C]-sucrose was included at a concentration of 3 μM in all experiments as an
> internal standard, with the exception of those conducted with [14C]-doxorubicin. Only
> MCC with a maximum variation of ±20% in sucrose penetration for a given experimental
> condition were used in data analysis."

> "To study penetration of radio-labelled anticancer drugs through MCC, they were
> dissolved in 2 × α-MEM and mixed in a 1:1 ratio with 1% agar solution to prevent
> convection."

**(h) What the paper does NOT report, and I will not invent:**

* **No MCC thickness.** The only size datum is "Experiments were conducted 3-6 times
  with MCC ranging in cell number from 3-5 × 106". Grantab 2006, the sibling that does
  report thicknesses, is closed access (`verification_routes_2026-09-20.md` §R4 records
  Unpaywall `"oa_status":"closed"`, no repository copy) and I did not obtain it.
* **No membrane area.** "Semi-porous Teflon membrane culture inserts (Millipore, Bedford
  MA)" with no diameter given.
* **No per-drug SD/SEM on the fold-changes.** They are stated with a tilde. Figure 1
  carries "mean of 3 or more experiments ± SEM" but the folds themselves are
  approximate, and the doxorubicin folds are "data not shown" in both experiments.

**(i) A correction to the scouting survey, recorded here because it changes the
interpretation.** `verification_routes_2026-09-20.md` R4 describes "two mechanistically
independent **non-toxic** treatments". Only EGTA is non-toxic. Bortezomib is not:

> "We did not observe significant alterations in cellular packing density in MCC derived
> from HCT-8 colon carcinoma cell lines using sub-cytotoxic doses of bortezomib, but
> doses that induced less than 50% cell kill decreased cellular packing density
> (Table 1)."

> "For studies of drug penetration, MCC were exposed for 24 h to 1 μM bortezomib or its
> diluent; this dose induced modest cell kill and led to the greatest change in packing
> density."

The paper is explicit that **sub-cytotoxic bortezomib did not move packing density at
all**; the packing change is inseparable from cell kill. The bortezomib arm therefore
carries a confound the EGTA arm does not, and this is pre-registered as a reason to
report the two perturbations separately as well as together, and never to describe the
bortezomib arm as "non-toxic".

### 0.4 The mapping from "% nuclear area" to `φ` — the weakest link, argued up front

**They are not the same quantity.** The paper measures the areal fraction of DAPI-stained
**nuclei** on 4 μm sections:

> "Four-micron sections were stained with haematoxylin and eosin, or with DAPI to
> quantify packing density. Packing density is presented as percentage nuclear area for
> all cell lines and treatment conditions without taking into account differences in
> nuclear size."

The engine's `φ` is the **cell** volume fraction. A nucleus is a proper subset of a cell,
so nuclear volume fraction `ν` and cell volume fraction `φ` satisfy `ν = f_N · φ` with
`f_N` the nuclear volume fraction *of a cell*, `0 < f_N < 1`. By Delesse's principle an
areal fraction on a section estimates the corresponding volume fraction, so the measured
"% nuclear area" estimates `ν`, not `φ`.

Three consequences, all of which are carried as limits rather than buried:

1. **The identity mapping `φ = nuclear area fraction` is physically inadmissible.** It
   asserts `f_N = 1`, i.e. that the cells are entirely nucleus. It is nevertheless
   computed below as the extreme low-sensitivity bound, and labelled as inadmissible.

2. **The ratio mapping is forced, given one stated assumption.** If neither nuclear size
   nor cell size changes under the perturbation — only the number of cells per unit
   volume — then `f_N` is unchanged and

   ```
   φ_treated / φ_control  =  ν_treated / ν_control  =  NA_treated / NA_control
   ```

   i.e. the *relative* change in `φ` equals the measured relative change in nuclear area
   fraction, with `f_N` cancelling. **The paper explicitly did not verify this
   assumption** — "without taking into account differences in nuclear size" — and both
   perturbations could plausibly violate it: EGTA chelates calcium, and bortezomib
   induces apoptosis, which is accompanied by nuclear condensation. Nuclear condensation
   would *shrink* nuclei, so part of the measured drop in nuclear area could be nuclear
   shrinkage rather than cell unpacking. **That failure mode makes the true change in
   `φ` smaller than the measured change in `NA`, which means the mapping used here
   over-states the engine's predicted effect.** It is therefore conservative in the
   direction that matters if the engine passes, and not conservative if it fails.

3. **The absolute anchor `φ_control` does not cancel and is not measured.** The ratio
   mapping fixes `φ_treated/φ_control` but not `φ_control` itself. Because the hindrance
   function is strongly non-linear in `φ`, the predicted penetration ratio depends on
   where on the curve the control sits. The **scored** arm uses the shipped
   `φ_control = 0.7892` — because that is what the product actually ships and what the
   route exists to score. The **physically admissible band** is
   `NA_control < φ_control < 1`, and a sweep across it is reported as the primary
   diagnostic.

**Pre-registered mapping (scored arm):** `φ_control = 0.7892` (shipped, every row);
`φ_treated = 0.7892 × NA_treated / NA_control`.

Which gives, from §0.3(a) and §0.3(c):

| perturbation | line | NA_ctrl | NA_trt | NA ratio | φ_ctrl | φ_trt |
|---|---|---|---|---|---|---|
| EGTA 10 mM/90 min | HCT-8Ea | 50.0 | 40.0 | 0.80000 | 0.7892 | 0.631360 |
| EGTA 10 mM/90 min | HCT-8E11 | 40.0 | 31.0 | 0.77500 | 0.7892 | 0.611630 |
| bortezomib 1 μM/24 h | HCT-8Ea | 53.2 | 44.6 | 0.83835 | 0.7892 | 0.661624 |
| bortezomib 1 μM/24 h | HCT-8E11 | 50.8 | 42.5 | 0.83661 | 0.7892 | 0.660256 |

Note for the record that the **four NA ratios span only 0.775–0.838**, a range of 1.08×.
Whatever the mapping, the engine is being asked to reproduce the measured penetration
ratios from four nearly identical perturbations.

Note also an inconsistency internal to the paper that no mapping can remove: the control
packing density of HCT-8E11 MCC is **40%** in the EGTA experiment (§0.3a) and **50.8%**
in Table 1 (§0.3c). That is a 27% relative disagreement between two measurements of the
same quantity in the same cell line in the same paper. The ratio mapping is insulated
from it — only the within-experiment ratio is used — but it is a direct measurement of
the reproducibility of the "% nuclear area" assay and is carried into §0.6 as such.

### 0.5 What the engine predicts, and why it is compound-independent

Compiled through the shipped path
(`compile_drug.build_profile` → `two_route.two_route_diffusivity`), at `φ = 0.7892`:

| compound | MW | logP | interstitial D (μm²/s) | P (μm/s) | membrane resistance ratio | transcellular gain |
|---|---|---|---|---|---|---|
| 5-fluorouracil | 130.08 | −0.89 | 323.800 | 4.57e−2 | 7.41e+2 | 1.00668 |
| doxorubicin | 543.52 | +1.27 | 184.740 | 1.79e−5 | 1.17e+6 | 1.0000046 |
| gemcitabine | 263.20 | (see below) | 247.958 | 2.33e−5 | 1.15e+6 | 1.0000045 |

All three sit at the **impermeable-cell limit**: the transcellular route contributes
between 0.0005% and 0.67% of the hindrance. Therefore the engine's predicted penetration
ratio is, to four decimal places, **the same number for all three compounds**, and is a
function of `φ` alone.

Two things follow, and both are pre-registered:

* **The gemcitabine logP is not load-bearing and no citable value is needed.** Gemcitabine
  is not in the engine's compound library and has no entry in
  `lipophilicity.EXPERIMENTAL_LOGP`. A sensitivity run at logP = +0.5 (1.9 log units
  above the literature neighbourhood, chosen to be implausibly lipophilic) moves the
  transcellular gain only from 1.0000045 to 1.00035. The gemcitabine rows are therefore
  scored at the impermeable limit and the descriptor gap is declared irrelevant, with the
  sensitivity reported.
* **The engine structurally cannot predict a compound-dependent response to a packing
  perturbation.** This is a falsifiable consequence of the law's form and is tested as
  P3.

### 0.6 The observable, and the one assumption that cannot be avoided

The measurement is `C₂(6 h)_treated / C₂(6 h)_control`. The engine emits a **diffusivity**
ratio `r_D = hindrance(φ_treated)/hindrance(φ_control)`. Mapping one to the other requires
a transport solve, and the transport solve requires the slab thickness and the membrane
area — **neither of which the paper reports** (§0.3h).

The regime dependence is not a detail. In a donor/receiver cell:

* at times long compared with the diffusive lag but short compared with equilibration,
  receiver concentration is linear in `D` and the concentration ratio **equals** `r_D`;
* at times comparable to or below the lag `L²/6D`, the ratio is **larger** than `r_D`,
  without bound as `t → 0`;
* as either arm approaches equilibrium, the ratio is **smaller** than `r_D`, tending to 1.

Since the unreported geometry can place 6 h anywhere on that curve, a transient solve
with free thickness and area could be made to produce almost any predicted ratio. **Using
it as the scored arm would be indistinguishable from tuning a parameter to hit the
anchor, and it is therefore refused.**

**Pre-registered scored observable:** `r_D`, the engine's diffusivity ratio, compared
directly against the measured concentration ratio, under the explicitly stated and
**unverified** assumption that the 6 h measurement lies in the linear regime. A transient
sweep over a wide geometry range is computed and reported as a **band**, labelled
unconstrained, and is not scored.

**Declared measurement noise floor.** The paper's own acceptance criterion is "a maximum
variation of ±20% in sucrose penetration for a given experimental condition" (§0.3g).
Two such measurements enter every ratio, so the ratios carry a floor of roughly
`√2 × 20% ≈ 28%` before any biological effect. Differences between measured ratios
smaller than that are not resolved by this experiment.

### 0.7 The scored rows (n = 8 distinct reported measurements)

The paper reports eight distinct fold-changes. Two of them are stated once for a pair of
drugs and two are stated once for a pair of cell lines; that lumping is the paper's, not
mine, and is recorded per row.

| # | perturbation | line(s) | drug(s) | measured fold, as printed |
|---|---|---|---|---|
| 1 | EGTA | HCT-8Ea | 5-FU and gemcitabine (lumped by the paper) | ~3.3–4.4 |
| 2 | EGTA | HCT-8E11 | 5-FU and gemcitabine (lumped by the paper) | ~1.5–1.9 |
| 3 | EGTA | HCT-8Ea | doxorubicin | ~2 (data not shown) |
| 4 | EGTA | HCT-8E11 | doxorubicin | ~2 (data not shown) |
| 5 | bortezomib | HCT-8Ea | 5-FU | ~1.7 |
| 6 | bortezomib | HCT-8Ea | gemcitabine | ~3 |
| 7 | bortezomib | HCT-8E11 | 5-FU | ~2 |
| 8 | bortezomib | HCT-8E11 | gemcitabine | ~2 |

Rows 3 and 4 share one printed value ("~2-fold") reported for both lines with separate
p-values; the same is true of the bortezomib doxorubicin value ("~1.8-fold"), which is
reported as a supplementary ninth and tenth row and is **excluded from the primary n = 8**
to avoid crediting the same printed number twice within one perturbation. Sensitivity
S4 re-runs the aggregate with all doxorubicin rows dropped entirely, because the
doxorubicin experiments alone carried no sucrose internal standard.

Ranges are scored at their **geometric midpoint**; the endpoints are reported and
re-scored as sensitivity S5.

### 0.8 Pre-registered criteria — the thresholds, fixed now

Let `R_i` be the measured fold and `P_i` the engine's predicted `r_D` for row `i`, and
`GMFE(P, R) = exp( mean_i | ln(P_i / R_i) | )`.

**P1 — accuracy.** FAIL if `GMFE(engine) > 2.0`.
*Justification, independent of these data:* the 2-fold rule is the standard acceptance
criterion for predicted-versus-observed in DMPK, and it is stricter than the 3.0 bar
this repository pre-set for `results/brain_slice_anchor.md`. **Demoted per §0.1** — P1
could have been chosen to fit and is reported for completeness, not as the verdict.

**P2 — anti-vacuity. This is a decisive criterion.** Define the φ-blind null: a model
in which the packing term carries no information, so `P_i ≡ 1.0` for every row. FAIL if
`GMFE(engine) > 0.8 × GMFE(null)`. A φ term that does not beat "packing density does
nothing" by at least 20% has not been tested by this experiment, and a pass on P1 alone
would be one of the vacuous green checks this repository has shipped seven of in a week.

**P3 — dispersion. This is a decisive criterion.** The measurement has a dynamic range:
`S_meas = max R_i / min R_i`. So does the prediction: `S_pred = max P_i / min P_i`. FAIL
if `ln(S_pred) / ln(S_meas) < 0.5`, i.e. if the φ term accounts for less than half the
log dynamic range of the phenomenon it is being credited with explaining. The 0.5 bar is
a judgement call and is declared as one: below it, packing density is not the dominant
explanation of these measurements, whatever the central tendency does.

**P4 — direction.** Every `P_i > 1`. **This criterion is VACUOUS BY CONSTRUCTION and is
not counted.** `φ` decreases in every row and the hindrance is monotone decreasing in
`φ`, so P4 can only pass. It is listed so that it cannot later be quoted as evidence.

**P5 — admissibility (diagnostic, not scored).** For each row, invert the law for the
`φ_control` that would reproduce the measured ratio exactly. Report whether it lies in
the physically admissible band `(NA_control, 1)`, and report the spread of the implied
values across rows. A `φ_control` outside the band means no packing fraction whatsoever
reproduces that measurement under this law.

**Overall verdict rule, fixed now: the headline is PASS only if P1, P2 and P3 all pass.
Any one of them failing makes the headline FAIL, and the failure leads the document.**

### 0.9 Pre-registered sensitivities

* **S1** — `φ_control` swept over the admissible band; implied `φ_control` per row.
* **S2** — the inadmissible identity mapping `φ_control = NA_control`, as the
  extreme low-sensitivity bound.
* **S3** — transient two-compartment solve over a wide thickness sweep; reported as an
  unconstrained band, not scored.
* **S4** — aggregate with all doxorubicin rows dropped (no internal standard).
* **S5** — aggregate at the printed range endpoints rather than geometric midpoints.
* **S6** — gemcitabine logP perturbed to an implausibly lipophilic +0.5.
* **S7** — EGTA arm and bortezomib arm scored separately, since only EGTA is non-toxic.

### 0.10 What a PASS would and would not license

A PASS licenses exactly this sentence and no other: *"The two-route hindrance law
reproduces, within a pre-declared band, the measured change in multicellular-layer drug
penetration produced by a measured change in cellular packing density, in a laboratory
and cell line the law was not built on, under a stated and unverified mapping from
nuclear area fraction to cell volume fraction."*

It would **not** license any absolute diffusivity, any cross-compound ordering, any kill
fraction, any clinical or diagnostic statement, or the word "validated".

---

*(End of pre-registration. Everything below was produced after this section was written
to disk.)*

## 1. Method

`experiments/packing_density_anchor.py`, run once, output in
`results/packing_density_anchor_data.json`. No parameter was tuned. The only engine
value used is the shipped `φ = 0.7892`, and the only thing done to it is to scale it by
a ratio the paper measured.

1. **Compile each compound through the shipped path** —
   `compile_drug.build_profile(..., two_route=False)` for the interstitial (matrix)
   diffusivity, then `two_route.two_route_diffusivity`. `two_route=False` on the first
   call is required, not stylistic: since `results/two_route_promotion.md`,
   `build_profile` applies the two-route term itself, so the default would hand the law
   a `D` that already carries it and apply it twice. This is the same call
   `experiments/two_route_transport.py::_engine_interstitial_D` makes.
2. **Map the measured packing change onto `φ`** by the pre-registered ratio mapping of
   §0.4: `φ_treated = 0.7892 × NA_treated/NA_control`.
3. **Predict** `r_D = hindrance(φ_treated)/hindrance(φ_control)`.
4. **Score** against the paper's `C₂(6 h)_treated / C₂(6 h)_control`, on the criteria
   fixed in §0.8.
5. **Sensitivities S1–S7** of §0.9, plus two post-hoc diagnostics (D1, D2) that are
   labelled as post-hoc and are not criteria.

### 1.1 Why the predicted ratio is compound-independent

Compiled at `φ = 0.7892`:

| compound | MW | logP | interstitial D (µm²/s) | P (µm/s) | membrane resistance ratio | transcellular gain |
|---|---|---|---|---|---|---|
| 5-fluorouracil | 130.08 | −0.89 | 323.8 | 4.57e−02 | 7.41e+02 | 1.0066822 |
| doxorubicin | 543.52 | +1.27 | 184.7 | 1.79e−05 | 1.17e+06 | 1.0000046 |
| gemcitabine | 263.20 | (n/a) | 248.0 | 2.33e−05 | 1.15e+06 | 1.0000045 |

All three sit at the **impermeable-cell limit**: the transcellular route contributes at
most **0.67%** of the hindrance. The predicted ratio is therefore a function of `φ`
alone and is the same number for all three compounds to three decimal places.

Two consequences, both pre-registered in §0.5:

* **The gemcitabine descriptor gap is not load-bearing.** Gemcitabine has no entry in
  `lipophilicity.EXPERIMENTAL_LOGP`. Sensitivity **S6** moves its logP from −1.4 to
  +0.5 — implausibly lipophilic — and the predicted ratio moves from 1.85366069 to
  1.85332677, a relative change of **1.8e−04**. No citable gemcitabine logP was needed
  and none was invented.
* **The engine structurally cannot emit a compound-dependent response to a packing
  perturbation.** That is a falsifiable consequence of the law's form, and it is what
  P3 and §4.4 test.

### 1.2 Transport-harness verification (for the unscored S3 arm only)

The two-compartment transient solver used in S3 was verified before use, against the
quasi-steady analytic limit `C₂/C∞ = 1 − exp(−kt)` with `k = (A·D/L)(1/V₁ + 1/V₂)`:

* analytic **0.84276283** vs numeric **0.84227468**, relative error **5.79e−04**;
* grid halving (dx and dt both halved): 0.76756348 → 0.76743125, relative difference
  **1.72e−04**.

This verifies the *harness*, not the engine, and licenses nothing about the physics.

---

## 2. The raw table — predicted versus measured

`meas` is as printed in the paper. Ranges are scored at their geometric midpoint per
§0.7; both endpoints are re-scored as sensitivity S5. `pred` is the engine's `r_D`.

| # | perturbation | line | drug | measured | predicted | pred/meas | in primary n=8 |
|---|---|---|---|---|---|---|---|
| 1 | EGTA | HCT-8Ea | 5-fluorouracil | ~3.3–4.4 | 1.8473 | 0.485 | **yes** |
| 2 | EGTA | HCT-8Ea | gemcitabine | ~3.3–4.4 | 1.8537 | 0.486 | no (paper prints one value for rows 1–2) |
| 3 | EGTA | HCT-8E11 | 5-fluorouracil | ~1.5–1.9 | 1.9604 | 1.161 | **yes** |
| 4 | EGTA | HCT-8E11 | gemcitabine | ~1.5–1.9 | 1.9676 | 1.166 | no (paper prints one value for rows 3–4) |
| 5 | EGTA | HCT-8Ea | doxorubicin | ~2 | 1.8537 | 0.927 | **yes** |
| 6 | EGTA | HCT-8E11 | doxorubicin | ~2 | 1.9676 | 0.984 | **yes** |
| 7 | bortezomib | HCT-8Ea | 5-fluorouracil | ~1.7 | 1.6771 | 0.987 | **yes** |
| 8 | bortezomib | HCT-8Ea | gemcitabine | ~3 | 1.6821 | 0.561 | **yes** |
| 9 | bortezomib | HCT-8E11 | 5-fluorouracil | ~2 | 1.6847 | 0.842 | **yes** |
| 10 | bortezomib | HCT-8E11 | gemcitabine | ~2 | 1.6898 | 0.845 | **yes** |
| 11 | bortezomib | HCT-8Ea | doxorubicin | ~1.8 | 1.6821 | 0.935 | no (supplementary; see §0.7) |
| 12 | bortezomib | HCT-8E11 | doxorubicin | ~1.8 | 1.6898 | 0.939 | no (supplementary; see §0.7) |

**Read this table as within-drug response ratios only.** Every row is one drug at one
dose through one culture, before versus after a perturbation. Nothing in this table
orders compounds by penetration, by deliverability or by anything else, and no two
compounds' concentrations appear anywhere in this document. Cross-compound ranking is
withdrawn (2026-08-21) and this route does not touch it.

### 2.1 The mapping that produced the `pred` column

| perturbation | line | NA_ctrl (%) | NA_trt (%) | NA ratio | φ_ctrl | φ_trt |
|---|---|---|---|---|---|---|
| EGTA 10 mM / 90 min | HCT-8Ea | 50.0 ± 7.3 | 40.0 ± 1.6 | 0.80000 | 0.7892 | 0.631360 |
| EGTA 10 mM / 90 min | HCT-8E11 | 40.0 ± 6.8 | 31.0 ± 4.5 | 0.77500 | 0.7892 | 0.611630 |
| bortezomib 1 µM / 24 h | HCT-8Ea | 53.2 ± 5.5 | 44.6 ± 5.7 | 0.83835 | 0.7892 | 0.661623 |
| bortezomib 1 µM / 24 h | HCT-8E11 | 50.8 ± 3.8 | 42.5 ± 3.8 | 0.83661 | 0.7892 | 0.660256 |

Every `φ_trt` is inside [0, 1): admissible. **The four packing perturbations span only
1.082×** (NA ratios 0.775–0.838) — the engine is being asked to reproduce measured
penetration ratios from 1.5 to 4.4 out of four nearly identical inputs. That framing,
noted in the pre-registration before anything was computed, is the whole story.

---

## 3. Sources — every number, with its URL and its sentence

Fetched 2026-09-20 from **https://pmc.ncbi.nlm.nih.gov/articles/PMC3407510/**
(HTTP 200, 195,531 bytes). Grantab RH, Tannock IF, *BMC Cancer* 2012;12:214,
PMID 22672469, PMC3407510, DOI 10.1186/1471-2407-12-214, open access. Every sentence is
quoted in full in §0.3 and reproduced in the `sentence` field of every row of
`results/packing_density_anchor_data.json`, so the code and this document cannot drift
apart. Nothing was taken from the scouting survey's transcription; §0.3(i) records one
substantive correction to it.

The doxorubicin folds in rows 5, 6, 11 and 12 are stated in the paper as **"data not
shown"**, and the doxorubicin experiments alone carried **no `[14C]`-sucrose internal
standard**. Sensitivity **S4** re-runs the aggregate without them.

---

## 4. Verdict

### 4.1 P3 — FAIL. The engine cannot produce the spread.

| | value |
|---|---|
| measured spread across the primary 8, `max/min` | **2.2572×** |
| predicted spread | **1.1732×** |
| `ln(S_pred) / ln(S_meas)` | **0.1962** |
| pre-registered bar | ≥ 0.50 |

The engine accounts for **19.6%** of the log dynamic range. The bar was set at 0.5 in
§0.8 on the stated ground that a term explaining less than half the log spread is not
the dominant explanation of the phenomenon; the measured 0.1962 is not near the bar,
and no plausible re-setting of a judgement-call threshold changes the conclusion.

### 4.2 P1 and P2 — PASS. This is a real result and is not vacuous.

| | value |
|---|---|
| `GMFE(engine, measured)` | **1.2678** |
| `GMFE(φ-blind null, measured)` — the term contributes nothing, ratio ≡ 1 | **2.1878** |
| pre-registered P2 bar, `0.8 × GMFE(null)` | 1.7503 |

The engine's geometric mean fold error is **42% lower** than the null in which packing
density has no effect at all. On central tendency, the two-route law's `φ` sensitivity
is approximately the right size: it predicts ~1.68–1.97× where the measurements
average near 2×.

**P2 is the criterion that stops this being one of the seven vacuous green checks this
repository shipped in a week.** It was possible for P1 to pass on a term carrying no
information; it did not. **P4 (direction) is recorded as VACUOUS BY CONSTRUCTION and is
not counted** — `φ` falls in every row and the hindrance is monotone decreasing in `φ`,
so P4 could only ever pass. It must never be quoted as evidence.

### 4.3 P5 — the implied `φ`, and the first external constraint on 0.7892

Inverting the law per row for the `φ_control` that reproduces each measurement exactly:

| # | perturbation | line | drug | measured | implied `φ_control` | admissible band | admissible |
|---|---|---|---|---|---|---|---|
| 1 | EGTA | HCT-8Ea | 5-fluorouracil | 3.811 | 0.929313 | (0.500, 1) | yes |
| 3 | EGTA | HCT-8E11 | 5-fluorouracil | 1.688 | 0.724223 | (0.400, 1) | yes |
| 5 | EGTA | HCT-8Ea | doxorubicin | 2.000 | 0.815523 | (0.500, 1) | yes |
| 6 | EGTA | HCT-8E11 | doxorubicin | 2.000 | 0.794899 | (0.400, 1) | yes |
| 7 | bortezomib | HCT-8Ea | 5-fluorouracil | 1.700 | 0.795016 | (0.532, 1) | yes |
| 8 | bortezomib | HCT-8Ea | gemcitabine | 3.000 | 0.919532 | (0.532, 1) | yes |
| 9 | bortezomib | HCT-8E11 | 5-fluorouracil | 2.000 | 0.848164 | (0.508, 1) | yes |
| 10 | bortezomib | HCT-8E11 | gemcitabine | 2.000 | 0.846758 | (0.508, 1) | yes |

**8 of 8 admissible.** Range **0.7242–0.9293**, geometric mean **0.8318**, median
0.8311, spread 1.2832×. **The shipped 0.7892 lies inside the range.**

This is the first external constraint of any kind on a value derived from this
repository's own Xenium section, and the constraint does not refute it. Two things keep
it honest, though:

* The band is **1.28× wide** and every point in it is admissible, so this is a weak
  constraint. It excludes `φ ≲ 0.72` and `φ ≳ 0.93` *under this mapping*, nothing more.
* The `φ` values are **spread precisely because P3 failed**: each row needs its own
  `φ_control` because one `φ` cannot fit them all. The implied-`φ` table and the
  dispersion failure are the same finding stated forwards and backwards. It would be
  dishonest to quote the "0.7892 is inside the range" result without the sentence that
  the range exists only because the law could not fit the rows with one value.

Note separately that `φ = 0.7892` sits above the 0.40–0.72 literature tumour band that
`two_route.py:125` cites against itself; the implied values here sit above it too.

### 4.4 The sharpest statement of the failure — D2 (post-hoc diagnostic)

The engine is a function of `φ` alone, so for a given (perturbation, cell line) it must
emit one number whatever the compound. What the measurement does:

| perturbation | line | NA ratio | engine predicts | measured rows | measured spread | predicted spread |
|---|---|---|---|---|---|---|
| EGTA | HCT-8Ea | 0.80000 | 1.8505 | 3.81 and 2.00 | **1.9053×** | 1.0037× |
| EGTA | HCT-8E11 | 0.77500 | 1.9640 | 1.69 and 2.00 | 1.1847× | 1.0037× |
| bortezomib | HCT-8Ea | 0.83835 | 1.6796 | 1.70 and 3.00 | **1.7647×** | 1.0030× |
| bortezomib | HCT-8E11 | 0.83661 | 1.6873 | 2.00 and 2.00 | 1.0000× | 1.0030× |

**One packing-density change, one predicted number, measured responses differing by up
to 1.9053×.** The paper's own acceptance gate implies a measurement noise floor of
roughly `√2 × 20% ≈ 28%` (§0.6), so the 1.76× and 1.91× spreads are well outside noise,
while the 1.18× spread is not resolved by this experiment and the 1.00× is a single
printed value covering two rows.

This is a *dispersion* statement, not an ordering: these are within-drug response
ratios, each dimensionless and internally standardised, and nothing here says which
compound penetrates further than which.

### 4.5 S3 — the unreported geometry is not an escape route

The measured observable is a concentration ratio at 6 h; the engine emits a diffusivity
ratio. Converting one to the other needs the slab thickness and the membrane area, and
the paper reports **neither** (§0.3h). §0.6 pre-registered the refusal to score that arm
for exactly the reason that a free geometry could be made to produce almost any ratio —
which would be indistinguishable from tuning to the anchor. It was nevertheless computed
as a band, to see whether it could rescue the result. It cannot.

| slab thickness (µm) | C₂/C∞ control | C₂/C∞ treated | transient ratio | `r_D` | amplification |
|---|---|---|---|---|---|
| 40 | 0.995522 | 0.999744 | 1.0042 | 1.8473 | 0.5436 |
| 100 | 0.885480 | 0.981422 | 1.1084 | 1.8473 | 0.6000 |
| 150 | 0.762462 | 0.929641 | 1.2193 | 1.8473 | 0.6600 |
| 200 | 0.657544 | 0.862289 | 1.3114 | 1.8473 | 0.7099 |
| 400 | 0.404784 | 0.620505 | 1.5329 | 1.8473 | 0.8298 |
| 800 | 0.208280 | 0.364172 | 1.7485 | 1.8473 | 0.9465 |
| 1600 | 0.073281 | 0.164293 | 2.2420 | 1.8473 | 1.2136 |
| 2400 | 0.022539 | 0.076278 | 3.3843 | 1.8473 | 1.8320 |
| 3200 | 0.005459 | 0.031792 | 5.8238 | 1.8473 | 3.1526 |

Across every thickness in the physically plausible range for a multicellular culture
(roughly 100–200 µm), the transient **compresses** the predicted ratio rather than
amplifying it — amplification 0.60–0.71 — because at the engine's own `D_eff` of
49.27 µm²/s the 6 h point is near equilibrium. Including this arm would make the engine
predict **less**, not more, and would worsen P1 and P2 while leaving P3 failed.

The transient arm first reaches the largest measured ratio (3.811×) at a slab thickness
of **3200 µm**. A 3.2 mm multicellular culture is not a thing; the cultures here carry
3–5 × 10⁶ cells on a Millipore insert. **The unreported geometry therefore cannot
explain the gap, in either direction.**

---

## 5. Sensitivities — the verdict does not move

| sensitivity | n | GMFE | null | P1 | P2 | `ln`-fraction | P3 |
|---|---|---|---|---|---|---|---|
| **primary** (§0.7, geometric midpoints) | 8 | 1.2678 | 2.1878 | PASS | PASS | 0.1962 | **FAIL** |
| **S4** all doxorubicin rows dropped | 6 | 1.3512 | 2.2543 | PASS | PASS | 0.1918 | **FAIL** |
| **S5** measured at range low endpoints | 8 | 1.2637 | 2.1173 | PASS | PASS | 0.2026 | **FAIL** |
| **S5** measured at range high endpoints | 8 | 1.2718 | 2.2607 | PASS | PASS | 0.1680 | **FAIL** |
| **S7** EGTA arm alone (the only non-toxic one) | 4 | 1.2731 | 2.2522 | PASS | PASS | 0.0775 | **FAIL** |
| **S7** bortezomib arm alone | 4 | 1.2624 | 2.1252 | PASS | PASS | 0.0133 | **FAIL** |

**S6** — gemcitabine logP −1.4 → +0.5: predicted ratio 1.85366069 → 1.85332677,
relative change 1.8e−04. Not load-bearing.

### 5.1 S7 — and why the two perturbations had to be separated

§0.3(i) corrected the scouting survey: only EGTA is non-toxic. Bortezomib at 1 µM
"induced modest cell kill", and the paper is explicit that **sub-cytotoxic bortezomib
did not move packing density at all**. The bortezomib packing change is therefore
inseparable from cell killing, and its penetration gain could reflect loss of binding
sites or loss of cells rather than a change in `φ`. Scored alone, the clean EGTA arm
fails P3 at **0.0775** — harder than the pooled result. The confound is not what is
driving the failure.

### 5.2 S2 — no mapping choice rescues P3

The physically inadmissible identity mapping `φ_control = NA_control` (which asserts the
cells are entirely nucleus) gives **GMFE 1.7944**, which **fails P2** at its 1.7503 bar,
and a predicted spread of only **1.047×** — *smaller* than the shipped mapping's 1.1732×.
Both admissible and inadmissible mappings fail P3, and the inadmissible one fails P2 as
well. **The dispersion failure is mapping-robust.**

### 5.3 The paper's own packing-density SD, propagated

Propagating ±1 SD on both control and treated nuclear area through the mapping:

| perturbation | line | `r_D` nominal | `r_D` low | `r_D` high |
|---|---|---|---|---|
| EGTA | HCT-8Ea | 1.8473 | 1.1037 | 2.4540 |
| EGTA | HCT-8E11 | 1.9604 | **0.7284** | 2.9761 |
| bortezomib | HCT-8Ea | 1.6771 | **0.7854** | 2.4903 |
| bortezomib | HCT-8E11 | 1.6847 | 1.0598 | 2.2684 |

**Two of the four low ends fall below 1.0** — i.e. at the extremes of the paper's own
per-section SD, the engine does not even determine the *direction* of the effect. This
must be read correctly and is easy to overstate in either direction: these are SDs, not
SEMs, and the paper reports the mean differences at p = 2.2e−09 and p = 8.1e−06, so the
*means* are well determined. What the table measures is **section-to-section
heterogeneity in packing density**, and what it shows is that a single section's packing
density is not a sufficient statistic for its own penetration behaviour. It is a
statement about the anchor's resolving power, and it cuts against P1 and P2 — the two
criteria that passed — not against P3.

---

## 6. What this does not show

**6.1 It does not show that nuclear area fraction is cell volume fraction.** This is the
weakest link and it was named as such in §0.4 before anything was computed. They are
different quantities: nuclear volume fraction `ν = f_N·φ` with `f_N` the nuclear fraction
of a cell. The ratio mapping used here is forced *only* if neither nuclear nor cell size
changes under the perturbation, and **the paper explicitly did not check** — "without
taking into account differences in nuclear size". Both perturbations could violate it:
EGTA chelates calcium; bortezomib induces apoptosis, which brings nuclear condensation.
If nuclei shrank, part of the measured drop in nuclear area is not unpacking at all, the
true change in `φ` is smaller than assumed, and the engine's predicted effect here is
**over-stated**. That direction would worsen P1 and P2 and leave P3 failed. The mapping
is therefore a threat to the two criteria that passed, not to the one that failed, and
every claim in §4.2 and §4.3 must be read as conditional on it.

**6.2 It does not test the functional form of the law.** Post-hoc diagnostic D1 re-ran
the primary 8 under three laws with the same monotonicity:

| law | GMFE | predicted spread | `ln`-fraction |
|---|---|---|---|
| Hasselman–Johnson `(1−φ)/(1+φ/2)` — the shipped law | 1.2668 | 1.1697 | 0.1926 |
| naive volume exclusion `(1−φ)` — no tortuosity at all | 1.3160 | 1.1477 | 0.1693 |
| 2-D Maxwell–Garnett `(1−φ)/(1+φ)` | 1.2512 | 1.1833 | 0.2068 |

They are indistinguishable — 1.25 to 1.32, and the *best* of the three is not the shipped
one. **This experiment constrains the magnitude of the `φ` sensitivity and says nothing
about which effective-medium law generates it.** A P1/P2 pass is not evidence for
Maxwell/Hasselman–Johnson specifically. (These use the closed-form impermeable limit, so
they differ from §4.2's 1.2678 by the ≤0.67% transcellular gain.)

**6.3 It does not show the engine's absolute diffusivity is right — and the repo already
knows it is not.** Every quantity scored here is a ratio, by design, so absolute
calibration cancels identically. `results/two_route_promotion.md` records that the MCL
absolute `D_eff` **FAILS** (49.33 µm²/s against a measured 74–130 band) while the MCL
hindrance ratio **PASSES**. This route adds a third reading and is consistent with that
pattern: **the ratios are approximately right and the absolutes are not tested here.**
It does not adjudicate that contradiction, which was R4's stated hope in
`results/verification_routes_2026-09-20.md`; R5's mannitol arm is the route that can,
because a pure paracellular marker isolates the extracellular leg.

**6.4 It does not license the word "validated", nor anything downstream of transport.**
No absolute kill fraction, no clinical-response prediction, no diagnostic or
treatment-planning output, no "digital twin", no personalised medicine, no AI/ML claim.
A `φ` term that reproduces a central tendency in a colon-carcinoma multicellular culture
says nothing about a patient.

**6.5 It does not produce or support a cross-compound ordering.** Every comparison is
one drug, one dose, one culture, before versus after. §4.4 compares the *dispersion* of
within-drug response ratios against the engine's predicted dispersion; it does not place
two compounds' penetration levels side by side, and no such comparison appears anywhere
in this document. Cross-compound ranking remains withdrawn (2026-08-21).

**6.6 It is one paper, one laboratory, two sublines of one cell line.** HCT-8Ea and
HCT-8E11 are variants of the same parental colon carcinoma line, so "two cell lines"
overstates the independence. The two perturbations are genuinely mechanistically
independent (calcium chelation versus proteasome inhibition), and that is the strongest
structural feature of the dataset — but only one of them is non-toxic.

**6.7 Three numbers the paper does not report, which are limits rather than gaps in this
analysis**: the MCC thickness, the membrane area, and any SD or SEM on the fold-changes
themselves (they are printed with a tilde). Grantab 2006, the sibling that does report
thicknesses, is closed access and was not obtained. Nothing was substituted for any of
them.

**6.8 One internal inconsistency in the source, carried rather than smoothed.** The
control packing density of HCT-8E11 MCC is **40%** in the EGTA experiment and **50.8%**
in Table 1 — a 27% relative disagreement between two measurements of the same quantity
in the same cell line in the same paper. The ratio mapping is insulated from it, since
only within-experiment ratios are used, but it is a direct measurement of the
reproducibility of the "% nuclear area" assay that this entire route rests on.

---

## 7. What this changes

* **`φ = 0.7892` now has an external constraint where it had none.** It is loose
  (0.7242–0.9293 under this mapping) and it does not refute the shipped value. The
  sentence "never scored against any external measurement" is no longer true of it; the
  sentence "verified" is not true of it either.
* **A new, registered finding about the law's reach:** the two-route hindrance is a
  function of `φ` alone, and measured penetration responses to a packing perturbation
  are not. A second mechanism is required. The paper's own framing points at one — EGTA
  works by disrupting E-cadherin-mediated cell–cell adhesion, i.e. by changing the
  *geometry of the paracellular path*, which is not the same variable as the volume
  fraction of cells.
* **Nothing in the engine was changed by this document.** No constant was tuned, no
  default moved, no file outside the three created for this route was touched.
  `calibrate_alpha()` and `DEFAULT_CELL_VOLUME_FRACTION` were asserted unchanged at the
  start and at the end of the run.

---

## 8. Reproduction

```
.venv/bin/python experiments/packing_density_anchor.py
```

Writes `results/packing_density_anchor_data.json`. Asserts
`calibrate_alpha() == 22.223477711172563` and
`two_route.DEFAULT_CELL_VOLUME_FRACTION == 0.7892` at both the start and the end of the
run, and fails loudly if either has moved.

**`calibrate_alpha()` START = 22.223477711172563 — asserted.**
**`calibrate_alpha()` END   = 22.223477711172563 — asserted.**
