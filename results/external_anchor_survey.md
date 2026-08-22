# External anchor survey — can the transport physics be validated against public data?

> # ⚠️ CORRECTED 2026-08-14 — THE SURVEY'S STRUCTURAL DIAGNOSIS OF THE DIFFUSIVITY MODEL NO LONGER DESCRIBES THE SHIPPING ENGINE, AND ITS HIGHEST-VALUE BLOCKED ANCHOR IS NO LONGER BLOCKED
>
> **This banner is the only edit made to this file.** No number, quotation or
> verdict below it has been changed or re-run. Everything below was true on
> **2026-08-09**.
>
> **What changed, in two steps.**
> 1. `d3b9955` (2026-08-13, *"Pruijn 2008 retrieved"*) obtained the
>    compound-by-compound table this survey records as unobtainable, and scored
>    the anchor.
> 2. `fd2a9e7` (2026-08-14, *"Flip two-route transport ON"*) changed the
>    diffusivity model this survey diagnoses. `compile_drug.build_profile` now
>    routes `logp` and `hbd` into `effective_D`. Record:
>    `results/two_route_promotion.md`.
>
> **Which statements are now false.**
>
> | where | statement | status |
> |---|---|---|
> | §2 (Hicks/TPZ row) | *"The engine derives `D_eff` from a Stokes–Einstein radius that is a function of **MW alone**… times a Brownian hindrance factor that depends on the **ECM fraction, not on the molecule**."* | **FALSE of the shipping engine.** The hindrance now depends on the molecule, through `logp` and `hbd`. |
> | §2 (same row) | *"Scoring the spread… would be a cheap and decisive test… It is blocked solely on the compound-by-compound table."* | **NO LONGER BLOCKED, AND NO LONGER PENDING.** The table was retrieved and the test was run — see `results/diffusivity_falsification.md` (the MW-only model, falsified) and `results/two_route_promotion.md` (the replacement, scored). This survey's call was correct: it was the highest-value anchor identified, and it decided the physics. |
> | §3.2 | *"the engine's only obstacle population is 'ECM'… so **cellular volume exclusion — which is most of a tumour's volume — is not in the hindrance law**."* | **FALSE.** This defect was the subject of `results/cellular_exclusion.md` and is closed. |
> | §3.2 | *"**The engine hinders diffusion about 7× too weakly**"*, engine 1.77-fold vs measured 12-fold, ratio **0.148×** | **SUPERSEDED.** At the promotion the MCL hindrance factor is **10.3274-fold** against the same measured 12-fold, and that row went **FAIL → PASS** (`results/two_route_promotion.md`). |
>
> **What still stands, and is the more durable half of this file.** The survey's
> premise — that the project sells transport physics and externally validates
> very little of it — is what motivated everything above, and the retrieval
> discipline (every value quoted with its DOI, nothing inferred from a title) is
> unchanged. Note also that the **absolute** MCL `D_eff` row still FAILS after the
> promotion; it merely fails in the opposite direction. The hindrance *ratio*
> passing is not the absolute value passing.

**Question this file answers.** `vivum validate` scores 13 anchor groups and prints
about itself that **11 of its 11 scored quick-tier rows never invoke the transport
PDE**. `results/heldout_patient_validation.md` found that the one anchor hoped to
be the first out-of-sample PDE test measures reconstruction/vessel geometry
computed *before* the solver runs. So the project sells transport physics and
externally validates almost none of it. This file is a search for measurements
that would close that gap, and an honest verdict on how far they can close it.

**Method and provenance rules.** Literature retrieved via PubMed and the open web
during this survey. According to PubMed, every value below was read from the
retrieved record (abstract, PMC full text, or the publisher's rendered tables)
and is quoted with its DOI. Nothing is inferred from a title. Where the primary
source could not be retrieved the row is marked **UNVERIFIED** with the specific
thing that blocked it. Nothing here is registered in
`tmesim/validation/scorecard.py`; this is a recommendation.

**No clinical-outcome anchors were sought or are proposed.** Vivum does not
predict patient response and an anchor requiring it would be out of scope.

---

## 0. The rule this survey is written under

**A measured value is not an anchor. An anchor is a value PLUS its tissue, its
marker or modality, its counting/sampling convention, and its uncertainty.**
Any candidate missing one of those is not usable yet, and marking it
NOT COMMENSURABLE or UNVERIFIED is the correct outcome rather than a failure.

This project has manufactured a false regression twice by pooling incompatible
modalities. The four compartments that keep getting conflated, and which the
engine distinguishes:

| compartment | what it is | does the PDE compute it? |
|---|---|:--:|
| **free interstitial** | unbound drug in extracellular extravascular fluid | **YES — this is `C`** |
| total tissue / homogenate | interstitial + intracellular + bound, per gram | no |
| intracellular | inside cells, incl. DNA-bound / tubulin-bound | no |
| plasma total vs plasma unbound | denominator of every partition ratio | source term is one of these — **which one must be stated** |

Microdialysis measures the first. Biopsy homogenate measures the second. They
differ by an order of magnitude for a bound drug and **must never be pooled**.

### 0.1 Cautionary example — the same failure mode is live, not historical

Three anchors already cited by this repo as external validation, checked against
their own primary records. **These are catalogued, not fixed** — other agents own
the affected files.

| PMID | used as | the commensurability problem |
|---|---|---|
| **38723085** — Shapson-Coe et al., *Science* 2024;384(6696):eadk4858, [DOI](https://doi.org/10.1126/science.adk4858) | 8.05 % endothelial-fraction reference for a **glioblastoma** section in `heldout_patient.py` | According to PubMed the specimen is *"a cubic millimeter of human temporal cortex that was surgically removed to gain access to an underlying epileptic focus"* — **non-tumour cortex from an epilepsy resection**. Two independent mismatches: (a) wrong tissue, GBM vs healthy cortex; (b) a **3-D volume fraction from nanoscale EM** compared against a 2-D spatial-transcriptomics area quantity. That row is currently NOT TESTABLE for an unrelated reason (circularity in the typing rule), so its verdict is right by accident. |
| **27583442** — Marien et al., *PLoS One* 2016;11(9):e0161496, [DOI](https://doi.org/10.1371/journal.pone.0161496) | the MVD reference | The abstract fixes a **counting convention** the model may not share: *"Every CD31-positive object was counted, with exclusion of suspected CD31-positive monocytes, macrophages and tumor cells. Furthermore, if interconnected, stained objects were considered a single vessel."* The same tissue yields a different MVD depending on whether interconnected objects merge. |
| **36611806** — Tamma et al., *Cells* 2022;12(1):11, [DOI](https://doi.org/10.3390/cells12010011) | the vascular-area-fraction reference | 30 IDH1-wild-type GBM biopsies, **CD34** microvessels by IHC — while the MVD anchor above counts **CD31**. The repo's two vascular anchors are **not the same assay**, even though both are "vasculature". |

The generalisable point: two of the three were adopted without their convention
travelling with them. That is the same error class as pooling free-interstitial
with total-tissue concentrations, and it got through the adoption process twice.

---

## 1. What already exists — searched so this survey does not re-find it

Read before searching: `results/breast_transport_anchors.md`,
`results/breast_kp_validation.md`, `results/pk_reality.md`,
`results/calibration.md`, `results/logp_provenance.md`,
`results/alpha_identifiability.md`, `tmesim/validation/pk_anchors.py`, and the
24-entry `OUT_OF_SCOPE` table in `tmesim/validation/scorecard.py`.

Already held: doxorubicin perivascular half-distance (PMID 16361566, **in-sample**
— `alpha` is solved to it in closed form); the joint dox/mitoxantrone/topotecan
25–75 µm range (PMID 23680920); TMZ brain-ECF:plasma AUC 0.178 (PMID 19861433);
MTX and cisplatin CSF:plasma (PMIDs 3855310, 8548876); breast interstitial:plasma
AUC for 5-FU 0.61 ± 0.11 and MTX 0.60 ± 0.20 (PMIDs 9205062, 9679959 — the engine
**FAILED** both, 12.8× high and 10.5× low); breast IFP 15 ± 9 and 29 ± 3 mmHg
(PMIDs 1423283, 7850532); human glioma ECF MTX elimination (PMID 18787762).

`results/alpha_identifiability.md` §4 names exactly what is missing: **(a)** a
free-interstitial first-order loss rate `k` for the anchor drug from microdialysis
with a resolvable washout, and **(b)** a perivascular profile for a *weakly-bound,
deeply-penetrating* drug (temozolomide) from a different chemical class. This
survey looked hard for both.

---

## 2. Ranked candidates

Ranking principle, per brief: **an anchor on an ABSOLUTE quantity is worth more
than one on a ratio**, because the engine's ratios are protected by p99
normalisation while its absolutes are known to be poorly constrained
(`results/pk_reality.md`: `k` was ~300× too small; `alpha_identifiability.md`:
one identified number in the whole calibration).

| # | candidate | measured quantity | value | modality | absolute? | commensurable? | actionable now? |
|---|---|---|---|---|:--:|:--:|:--:|
| **1** | Hicks 1998 / Hicks 2007 MCL | tissue diffusion coefficient `D` **and** first-order consumption `k` | D 74–130 µm²/s; hindrance 12-fold; k 0.31–0.52 /min | MCL flux fitted to a diffusion-reaction model | **YES (both)** | **YES for D and hindrance**; PARTIAL for k | **YES — prototyped, §3** |
| **2** | Boucher 1997 | human intracranial-tumour IFP | **2.0 ± 2.5 mmHg**, n=11 | wick-in-needle / manometry | **YES** | **YES** (engine solves IFP in mmHg) | **YES**, with a stated confound |
| **3** | Mader 2003 | capecitabine tumour:plasma AUC, human breast lesion | **1.76 (0.54–7.07)**, n=10 | microdialysis vs **ultrafiltered** plasma | ratio | **YES for the parent drug only** | **YES** — reuses `breast_kp_anchor.py` |
| 4 | Tunggal 1999 | relative MCL flux, 4 panel drugs | dox & mitox **<10 %**; MTX & 5-FU **~30–50 %** of Teflon-alone | MCL flux vs blank membrane | ratio | PARTIAL — needs 1-D slab geometry | NO — new geometry |
| 5 | Cowan 2001 | MTX MCL penetration, **concentration-invariant** 10 µM→3 mM | ~25 % of Teflon-alone at 6 h | MCL flux | ratio | PARTIAL — but a clean **linearity** test | NO — new geometry |
| 6 | Pruijn 2008 | `D_mcl` for **67 + 4 compounds** | range **2.7–187 µm²/s** (from PMID 16912264) | MCL flux | **YES** | **YES in principle** | **NO — publisher-blocked** |
| 7 | Jørgensen 2024 | doxorubicin unbound tissue AUC vs plasma | ratios 0.03–0.22 | microdialysis | ratio | **NOT COMMENSURABLE as published** | NO |
| 8 | Ekstrøm 1997 | human tumour MTX interstitial : unbound plasma | **9.3–14 %**, n=**1** | microdialysis, **no recovery correction** | ratio | lower bound only | NO |
| 9 | Grossman 2013 | TMZ tumour ECF Cmax / AUC / t½ | 0.59 ± 0.22 µg/mL; 1.82 ± 0.64 µg·h/mL; 1.7 ± 0.3 h | microdialysis, recovery-corrected | absolute | **NO plasma measured** | **NO** |
| 10 | Brana 2014 | human doxorubicin perivascular gradient | *no number retrieved* | IF vs vessel distance, 1 biopsy | — | **UNVERIFIED** | NO |

### 2.1 — #1 Multicellular-layer diffusion coefficient and consumption rate **[HIGHEST VALUE]**

> Hicks KO, Fleming Y, Siim BG, Koch CJ, Wilson WR. "Extravascular diffusion of
> tirapazamine: effect of metabolic consumption assessed using the multicellular
> layer model." *Int J Radiat Oncol Biol Phys* 1998;42(3):641-9. **PMID 9806526**,
> [DOI](https://doi.org/10.1016/s0360-3016(98)00268-5)
>
> Hicks KO, Myint H, Patterson AV, Pruijn FB, Siim BG, Patel K, Wilson WR.
> *Int J Radiat Oncol Biol Phys* 2007;69(2):560-71. **PMID 17869669**,
> [DOI](https://doi.org/10.1016/j.ijrobp.2007.05.049)

| quantity | value | units | uncertainty | system |
|---|---:|---|---|---|
| `D_mcl` tirapazamine | 7.4e-7 → **74** | cm²/s → µm²/s | not stated | V79-171b MCL, ~120 µm |
| `D_mcl` tirapazamine | 1.3e-6 → **130** | cm²/s → µm²/s | not stated | MGH-U1 MCL, ~360 µm |
| `D_mcl` tirapazamine | 1.30e-6 → **130** | cm²/s → µm²/s | **± 0.05e-6 (± 5)** | SiHa MCL (independent replicate) |
| `D_mcl` PR-104A | 4.42e-7 → **44.2** | cm²/s → µm²/s | ± 0.15e-7 (± 1.5) | SiHa MCL |
| hindrance vs culture medium | **12** | fold | not stated | V79-171b MCL |
| anoxic consumption `k` | 0.52 → **8.67e-3** | /min → /s | not stated | V79-171b MCL |
| anoxic consumption `k` | 0.31 → **5.17e-3** | /min → /s | not stated | MGH-U1 MCL |

**Why this is the single most commensurable measurement found.** The engine's
transport PDE has exactly two per-drug coefficients — `D_eff` and
`k = clearance + alpha·binding` — and Hicks fitted flux through a multicellular
layer **to a diffusion-reaction model of the same form** to obtain both. Same
units, same functional role, same equation. There is no unit conversion, no
compartment translation, and no modality bridge to get wrong. It is also
**absolute on both axes**, which is precisely where the project is weakest.

**COMMENSURABILITY VERDICT.**
* `D_mcl` vs `MolecularProfile.effective_diffusivity` — **COMMENSURABLE.** Both
  are the tissue-scale effective diffusion coefficient of a small molecule in the
  same continuum PDE. One caveat must travel with it: an MCL's obstacles are
  **packed cells with little stroma**, whereas the engine's only obstacle
  population is "ECM" at a fraction of 0.4. That is not a reason to reject the
  anchor — it is the reason the anchor is informative, because it asks whether
  the hindrance law accounts for cellular volume exclusion at all.
* hindrance factor (12-fold) — **COMMENSURABLE and dimensionless.** It is
  `D_medium / D_tissue`, exactly `1 / MolecularProfile.hindrance`. This is the
  cleanest single row in the whole survey: no units, no compartment, no cohort.
* `k` (anoxic consumption) — **PARTIALLY COMMENSURABLE, stated explicitly.** It is
  a first-order loss of free drug from the interstitium, in 1/s, in the same
  position in the same PDE. But the *process* is anoxic bioreductive metabolism of
  a hypoxia-activated prodrug, not the engine's generic binding+clearance sink.
  Read it as a **measured order-of-magnitude bound** on how large a real
  interstitial consumption rate is, **not** as a like-for-like score.
* Tirapazamine is not on the engine's panel, but the engine compiles arbitrary
  molecules from SMILES (`molecular/compile_drug.py`), so this is fully
  **out-of-sample**: nothing in the repo was ever fitted to a TPZ measurement.

### 2.2 — #2 Human intracranial-tumour interstitial fluid pressure

> Boucher Y, Salehi H, Witwer B, Harsh GR, Jain RK. "Interstitial fluid pressure
> in intracranial tumours in patients and in rodents." *Br J Cancer*
> 1997;75(6):829-36. **PMID 9062403**, [DOI](https://doi.org/10.1038/bjc.1997.148)

**Measured: mean IFP 2.0 ± 2.5 mmHg in 11 patients with intracranial tumours.**

**COMMENSURABLE — YES.** `tmesim/physics/pressure.py` solves a screened-Poisson
Darcy model and returns IFP **in mmHg**, unlike the oxygen field whose
mmHg↔dimensionless mapping was never derived (which is why
`results/oxygen_calibration_*.md` are honest nulls). The engine's IFP is a real
dimensioned field and nothing external currently scores it.

**Why it is valuable.** It is *low*. The engine's GBM flagship is expected to
produce elevated core IFP — the Jain effect is the module's stated purpose — and
the human intracranial measurement says 2.0 mmHg, an order of magnitude below the
breast numbers already surveyed (15 ± 9 and 29 ± 3 mmHg). A tissue-specific
anchor that runs *against* the module's headline mechanism is worth more than one
that agrees with it.

**Two caveats that must travel with it, both from the paper itself.**
1. The authors state the low value is attributable in part to **anti-oedema
   therapy**: *"the results also suggest that antioedema therapy in patients with
   brain tumours is responsible in part for the low tumour IFP."* Patients were
   treated with dexamethasone/mannitol/furosemide. The engine models no such
   intervention. This is a real confound and makes 2.0 mmHg a **treated-state**
   measurement, not an untreated-tumour baseline.
2. **Near-circularity on the engine side.** `pressure.py` hardcodes
   `_PV = 25.0` mmHg (effective microvascular driving pressure) and its own
   docstring says the parameters are *"order-of-magnitude biological values, not
   patient-calibrated."* Predicted core IFP is bounded by `_PV`. So a score
   against the *breast* 15–29 mmHg range would be close to circular — `_PV` was
   chosen at 25, sitting inside it. The **brain** value at 2.0 mmHg is the
   non-circular direction, because nothing in the model was tuned toward it.

### 2.3 — #3 Capecitabine in human breast lesions, four compounds, one modality

> Mader RM, Schrolnberger C, Rizovski B, Brunner M, Wenzel C, Locker G,
> Eichler HG, Mueller M, Steger GG. "Penetration of capecitabine and its
> metabolites into malignant and healthy tissues of patients with advanced breast
> cancer." *Br J Cancer* 2003;88(5):782-7. **PMID 12618890**,
> [DOI](https://doi.org/10.1038/sj.bjc.6600809) (PMC2376334, full text retrieved)

Ratio of AUC(0–5 h), **median (range), n = 10 patients**, read from the paper's
Table 2 (*"Transcapillary transfer of capecitabine and its metabolites to
malignant and healthy tissues"*):

| ratio | capecitabine | DFCR | DFUR | 5-FU |
|---|---|---|---|---|
| **tumour / plasma** | **1.76 (0.54–7.07)** | 0.67 (0.45–2.97) | 0.53 (0.28–2.56) | 0.66 (0–5.56) |
| s.c. tissue / plasma | 1.31 (0.46–16.7) | 0.67 (0.41–2.60) | 0.44 (0.25–1.36) | 1.53 (0.32–3.38) |
| **tumour / s.c. tissue** | **1.11 (0.36–1.99)** | 1.12 (0.55–1.93) | 1.20 (0.60–3.07) | 0.66 (0–17.5) |

**The denominator is right, and this is the rare case where that is documented.**
The paper states: *"100 µl plasma were centrifuged through Centrisart-C4 filter
units … to separate free drug from protein-bound drug. After ultrafiltration, the
analytical results obtained from blood samples were directly compared with those
obtained from microdialysates, which contained only free drug."* So this is
**free interstitial : free plasma** — exactly the engine's state variable over
exactly the right source-term units. Probes were recovery-calibrated by
retrodialysis (in-tumour recovery ranges given per compound).

**COMMENSURABILITY VERDICT — split, and the split is the whole point.**
* **capecitabine (parent): COMMENSURABLE.** Absorbed, delivered from plasma,
  measured free in tumour interstitium against free plasma. This is a delivery
  measurement.
* **DFCR, DFUR, 5-FU: NOT COMMENSURABLE.** These are metabolites, and the paper's
  own discussion establishes they are generated *en route and intracellularly*
  (thymidine phosphorylase inside tumour cells) and then effluxed:
  *"once converted to FU intracellularly, this agent is hardly redistributed to
  the systemic circulation."* Their interstitial concentration therefore contains
  a **local generation term the engine does not model**. Scoring the engine's
  delivery-only prediction against them would be a modality conflation of exactly
  the kind §0 warns about. **Do not score them.**
* Tissue is a **cutaneous chest-wall metastasis**, not primary breast. State it.

**Two things make it worth adopting anyway.** First, the `tumour / s.c. tissue`
column is a **within-patient negative control**: same probe technology, same
patient, same 5 h window, tumour vs healthy connective tissue. The engine, given
tumour vs normal geometry, should predict ≈1.1 for capecitabine, and a model that
manufactures a large tumour-vs-normal contrast is falsified without any absolute
calibration being involved. Second, the 5-FU tumour/plasma median of **0.66**
lands beside the repo's existing primary-breast 5-FU anchor of **0.61 ± 0.11**
(PMID 9205062) from a different cohort, different route and different lesion type
— a genuine corroboration of the anchor the engine already fails by **12.8×**.

### 2.4 — #4/#5 Multicellular-layer flux on the engine's own panel drugs

> Tunggal JK, Cowan DS, Shaikh H, Tannock IF. *Clin Cancer Res* 1999;5(6):1583-6.
> **PMID 10389947** — no DOI in the PubMed record.
>
> Cowan DS, Tannock IF. *Int J Cancer* 2001;91(1):120-5. **PMID 11149410**,
> [DOI](https://doi.org/10.1002/1097-0215(20010101)91:1<120::aid-ijc1021>3.0.co;2-y)

Tunggal 1999, EMT6 and MCF7 MCLs ~200 µm: *"The penetration of doxorubicin and
mitoxantrone was limited and very slow (<10% of the rate of penetration through
the Teflon membrane alone). The penetration of methotrexate and 5-FU was more
rapid (approximately 30–50% of the rate…)."*

Cowan 2001, EMT6 and MCF7 MCLs: MTX penetration *"only about 25% of that through
the Teflon membrane alone after 6 hr. This was true for all methotrexate
concentrations examined from 10 µM to 3 mM."*

**Why these matter more than their coarseness suggests.** All four compounds —
doxorubicin, mitoxantrone, methotrexate, 5-FU — are **on the engine's panel**, and
the measurement says doxorubicin and mitoxantrone are **3–5× slower** than
methotrexate. The engine's calibrated panel puts doxorubicin at `L` = 64.9 µm and
methotrexate at 75.2 µm (`results/alpha_identifiability.md` §2) — essentially
equal. That is a directional disagreement on the engine's own drugs.

Cowan's **concentration invariance from 10 µM to 3 mM** is separately valuable: it
is a direct experimental test of the engine's **linearity assumption** (a
first-order sink `−kC` with `k` independent of `C`). The measurement supports
linearity over 300-fold in concentration. That is one of the very few places the
literature confirms rather than refutes a structural choice here.

**COMMENSURABLE — PARTIAL. NOT ACTIONABLE without new code.** "% of the rate
through the blank membrane" is a ratio of steady-state fluxes through a 1-D slab
with reaction. The engine solves a 2-D reconstructed tissue with vessel sources
and has **no 1-D slab geometry**. Computing this would require a new solve
configuration, which is out of scope for this survey and touches `tmesim/`.

### 2.5 — #6 The Pruijn 71-compound dataset — **still publisher-blocked**

> Pruijn FB, Patel K, Hay MP, Wilson WR, Hicks KO. "Prediction of Tumour Tissue
> Diffusion Coefficients of Hypoxia-Activated Prodrugs from Physicochemical
> Parameters." *Aust J Chem* 2008;61(9):687-693.
> [DOI](https://doi.org/10.1071/CH08240)

Retried via a different route, as instructed. **Still not obtained.** The
publisher's article page returns **HTTP 403 Forbidden** to automated retrieval,
and the journal is not indexed in PubMed for this record, so there is no PMC or
PubMed abstract path either. What is verifiable from the publisher's and
Monash University's public listings: 67 structurally diverse benzotriazine
di-N-oxides plus 4 miscellaneous compounds, `D_mcl` measured in HT29 MCLs, and a
QSAR from MW, logP(7.4), HBD and HBA with R² = 0.907 and Q² = 0.879.

**The 71 individual `D_mcl` values were NOT retrieved and are recorded as
UNVERIFIED.** The identity of the "four miscellaneous compounds" was not
retrieved either — that matters, because if any of them are conventional
cytotoxics the dataset would reach the engine's actual panel.

**However, a bounded version of the same information IS verified**, from a paper
that is retrievable:

> Hicks KO, Pruijn FB, Secomb TW, Hay MP, Hsu R, Brown JM, Denny WA, Dewhirst MW,
> Wilson WR. *J Natl Cancer Inst* 2006;98(16):1118-28. **PMID 16912264**,
> [DOI](https://doi.org/10.1093/jnci/djj306)

*"A large range of tissue diffusion coefficients (0.027 × 10⁻⁶ – 1.87 × 10⁻⁶
cm²/s) was observed for the TPZ analogs"* — i.e. **2.7 to 187 µm²/s across 16
compounds measured by one method in one cell line (HT29)**.

**This is the most falsifying single number in the survey and it is verified.**
The engine derives `D_eff` from a Stokes–Einstein radius that is a function of
**MW alone** (`descriptors.hydrodynamic_radius_nm`: `r = 0.066·MW^(1/3)` for small
molecules) times a Brownian hindrance factor that depends on the **ECM fraction,
not on the molecule**. Across a congeneric series of similar MW, the engine must
predict a near-constant `D_eff`, while the measurement spans **69-fold**.
Scoring the spread — not the absolute values — needs only the compound list, and
would be a cheap and decisive test of the diffusivity model. It is blocked solely
on the compound-by-compound table.

### 2.6 — #7 Doxorubicin microdialysis — the right drug, the wrong denominator

> Jørgensen AR, Bue M, Hanberg P, Petersen EK, Harlev C, Hansen J,
> Baad-Hansen T, Safwat A, Stilling M. "Doxorubicin concentrations in bone
> tumour-relevant tissues after bolus and continuous infusion: a randomized
> porcine microdialysis study." *Cancer Chemother Pharmacol* 2024;93(6):555-564.
> **PMID 38332155**, [DOI](https://doi.org/10.1007/s00280-023-04637-1)
> (PMC11130026, tables retrieved)

Published `AUC_tissue / AUC_plasma`, mean (95 % CI), 16 pigs, 150 mg doxorubicin:

| compartment | bolus | continuous |
|---|---|---|
| cancellous bone | 0.03 (−0.01; 0.07) | 0.04 (0.01; 0.08) |
| subcutaneous | 0.15 (0.06; 0.24) | 0.12 (0.03; 0.20) |
| muscle | 0.09 (0.04; 0.13) | 0.06 (0.02; 0.11) |
| synovial fluid | 0.22 (0.16; 0.28) | 0.10 (0.05; 0.16) |

**COMMENSURABILITY VERDICT — NOT COMMENSURABLE AS PUBLISHED.** Three independent
blockers, any one of which is disqualifying:

1. **Mixed denominators.** The paper states: *"All concentrations quantified in
   the dialysate are of unbound doxorubicin, while the plasma samples represent
   the total concentrations (bound + unbound)."* The published ratio is therefore
   **unbound-tissue over total-plasma**. That is not `Kp,uu` and it is not the
   engine's quantity. This is textbook modality conflation and it is sitting
   inside the paper's headline table.
2. **No tumour.** All four compartments are normal porcine tissues (bone,
   subcutis, muscle, synovium). The engine reconstructs tumour tissue from spatial
   transcriptomics; it cannot represent pig cancellous bone.
3. **No washout kinetics.** No terminal elimination half-life is reported for any
   compartment, so this does **not** deliver the free-interstitial loss rate `k`
   that `results/alpha_identifiability.md` §4(a) asks for.

**What is salvageable, flagged as DERIVED not measured.** An intravenous
microdialysis catheter in the same animals gives **unbound plasma** exposure
(AUC 14.9 bolus / 12.5 continuous, in the paper's units). Dividing tissue AUC by
*that* rather than by total plasma would give an unbound:unbound ratio. **This
was not done here and must not be quoted as a measured value** — a ratio of group
summary statistics is not the paper's per-animal ratio, and the paper's own CIs
cross zero in two compartments.

### 2.7 — #9/#10 Two near-misses, recorded so they are not re-found

**Temozolomide interstitial ECF in glioma — the plasma leg is missing.**
> Grossman R, Tyler B, Rudek MA, Kim E, Zadnik P, Khan U, Blakeley JO,
> Pathak AP, Brem H. *Cancer Chemother Pharmacol* 2013;72(1):93-100.
> **PMID 23649683**, [DOI](https://doi.org/10.1007/s00280-013-2172-3) (PMC4596243)

TMZ-alone arm, U87 intracerebral glioma in nude rats, n = 7, recovery-corrected
(recovery 70.1 ± 10.8 %): **ECF Cmax 0.59 ± 0.22 µg/mL**, **AUC(0–∞)
1.82 ± 0.64 µg·h/mL**, **tumour-ECF t½ 1.7 ± 0.3 h**.

This is a genuine free-interstitial measurement of a panel drug in a glioma — and
it is **NOT ACTIONABLE**, because the study measured **no plasma pharmacokinetics
at all**. There is no denominator. Manufacturing one from a different TMZ cohort
is precisely the cross-study pooling that produced this project's earlier false
regressions. Recorded as a near-miss; the missing piece is a paired plasma AUC in
the same animals.

**A human doxorubicin perivascular gradient exists but could not be read.**
> Brana I, Ocana A, Chen EX, Razak ARA, Haines C, Lee C, Douglas S, Wang L,
> Siu LL, Tannock IF, Bedard PL. *Invest New Drugs* 2014;32(6):1269-77.
> **PMID 25213162**, [DOI](https://doi.org/10.1007/s10637-014-0159-5)

The abstract states: *"A single on-treatment tumor biopsy showed a sharply
decreasing gradient in doxorubicin concentration and associated activity markers
with increasing distance from tumor blood vessels."* This would be a **human**
counterpart to the murine Primeau anchor the whole calibration rests on.

**UNVERIFIED. Blocked by:** the Springer full text redirects to an authentication
endpoint and is not retrievable; there is no PMC deposit. Per the standing rule, a
paywalled paper is a stated gap and **no number is inferred from the abstract** —
note the abstract gives no distance, no decay length and n = 1 biopsy, so even if
retrieved it may not yield a scoreable value.

### 2.8 — Genuine nulls: searched, not found

* **A perivascular penetration profile for temozolomide, or for any weakly-bound
  deeply-penetrating drug, by the Tannock immunofluorescence-vs-CD31 method.**
  This is request (b) of `alpha_identifiability.md` §4 and it **does not exist in
  the accessible literature.** The method depends on drug autofluorescence or a
  specific antibody; TMZ has neither. The nine unconstrained binding ratios in the
  calibration remain unconstrained, and no search strategy found here changes that.
* **A free-interstitial washout-derived `k` for doxorubicin in tumour tissue.**
  Request (a) of the same section. Not found. The nearest approaches are the MCL
  consumption rate (right units and right equation, wrong process, different drug)
  and the porcine study (right drug, no washout, no tumour).
* **A public MALDI-MSI dataset of drug distribution in tumour sections with
  paired vessel staining and released raw data.** The modality is real and
  growing, and qualitative co-registration with vessel IHC is published, but no
  openly downloadable dataset with a quantitative drug-vs-vessel-distance profile
  was located. Recorded as the highest-upside *future* class.
* **A tumour IFP measurement not confounded by anti-oedema or cytotoxic therapy in
  the brain.** PMID 9062403's own authors attribute their low value partly to
  therapy; PMID 15774788 (Taghian 2005) measures a *treatment-induced change*
  (paclitaxel −36 % IFP, n = 54) which the engine cannot represent because it does
  not model a drug's effect on vessel or stromal mechanics.

---

## 3. Prototype — scoring the engine against candidate #1

Driver: `experiments/anchor_prototype.py`. **Nothing is registered.** It compiles
tirapazamine through the unmodified molecular pipeline
(`descriptors.from_smiles` → `compile_drug.build_profile`, ECM reference 0.4, the
repo default) and reads off the two PDE coefficients. SMILES validated by
round-trip: RDKit returns MW **178.151** against the published tirapazamine MW
178.15.

Compiled: `free_diffusivity` 509.47 µm²/s, `effective_diffusivity` **287.32
µm²/s**, hindrance 0.5640, `clearance` 6.946e-05 /s, `tissue_binding_rate`
2.800e-05 /s, `alpha` 148.6082.

### 3.1 Effective diffusion coefficient — ABSOLUTE

| MCL cell line | measured | engine | ratio | PMID |
|---|---:|---:|---:|---|
| V79-171b | 74 µm²/s | 287.3 | **3.88× over** | 9806526 |
| MGH-U1 | 130 µm²/s | 287.3 | **2.21× over** | 9806526 |
| SiHa | 130 ± 5 µm²/s | 287.3 | **2.21× over** | 17869669 |

Measured range across three independent cell lines: **74–130 µm²/s**. The engine
predicts **287.3**, i.e. **OUTSIDE** the measured range, though inside the loose
3× band this repo elsewhere calls an order-of-magnitude pass. Note the three
measurements agree with each other to 1.8× while the engine misses all three in
the same direction — a systematic bias, not scatter.

### 3.2 Hindrance factor — ABSOLUTE, dimensionless, **the sharpest result**

| | value |
|---|---:|
| measured, V79 MCL vs culture medium (PMID 9806526) | **12.0-fold** |
| engine, Brownian obstacles at ECM = 0.4 | **1.77-fold** |
| ratio engine/measured | **0.148×** |

**The engine hinders diffusion about 7× too weakly.** This row needs no units, no
compartment translation and no cohort matching — it is a pure dimensionless ratio
of the same drug's diffusivity in tissue versus in free medium. The diagnosis is
structural and visible in the code: the engine's only obstacle population is
"ECM", and `hydrodynamic_radius_nm` makes the tracer a function of MW alone, so
**cellular volume exclusion — which is most of a tumour's volume — is not in the
hindrance law**. An MCL is almost pure packed cells, which is exactly why it
exposes this.

### 3.3 First-order loss rate — ABSOLUTE, and an unexpected result

| MCL cell line | measured `k` | engine `k_raw` | ratio | engine `k` with `alpha` | ratio |
|---|---:|---:|---:|---:|---:|
| V79-171b | 8.667e-03 /s | 9.747e-05 | **0.0112×** | 4.231e-03 | **0.488×** |
| MGH-U1 | 5.167e-03 /s | 9.747e-05 | **0.0189×** | 4.231e-03 | **0.819×** |

Two findings, and the second must be read carefully.

1. **The uncalibrated sink is ~53–89× too small.** This independently reproduces
   `results/pk_reality.md`'s ~300× diagnosis on a *different drug*, in a
   *different system*, against a *directly measured rate constant* rather than
   against a penetration depth. The Stage-0 binding/clearance priors are too small
   by one to two orders of magnitude, and that is now supported by two
   unrelated measurements instead of one.
2. **With `alpha` applied, the engine lands within a factor of 1.2–2.0 of a
   genuinely measured, fully out-of-sample tissue consumption rate.** `alpha` was
   fitted in closed form to a doxorubicin *penetration depth* in mouse tumours;
   nothing about tirapazamine, MCLs or a rate constant entered that fit.

**Do not oversell finding 2, and this file will not.** The measured constant is
anoxic bioreductive metabolism of a hypoxia-activated prodrug — a real,
drug-specific consumption process that doxorubicin does not have — while the
engine's `k` is a generic binding+clearance sink. Same units, same slot in the
same equation, **different physical process.** The honest statement is: *the
magnitude `alpha` implies for a small-molecule interstitial loss rate is not
absurd, and is within 2× of one independently measured rate constant of the same
class.* That is weaker than validation and stronger than nothing, and it is the
first out-of-sample evidence about `alpha`'s **magnitude** rather than its
**ratios** that this project has.

### 3.4 What the two coefficients imply together

| | `L = sqrt(D/k)` |
|---|---:|
| V79-171b, **measured D and measured k** | **92.4 µm** |
| MGH-U1, **measured D and measured k** | **158.6 µm** |
| engine, predicted D and `alpha`-calibrated k | 260.6 µm |
| engine, predicted D and unamplified k | 1716.9 µm |

Computed coefficient-for-coefficient so that no measured value is ever mixed with
a predicted one. The engine over-predicts penetration length by **1.6–2.8×**
against a fully measured `L` — much better than the 17.7× in
`results/pk_reality.md`, and the improvement is attributable to `alpha`. The
residual is the diffusivity error of §3.1–3.2, which `alpha` cannot absorb because
`alpha` multiplies only the sink.

---

## 4. Recommendation

**Adopt (in this order).**

1. **The MCL hindrance factor, 12-fold, PMID 9806526** — as a **FAIL** row on
   `MolecularProfile.hindrance`. Dimensionless, absolute, out-of-sample, one
   number, no convention to get wrong, and it isolates a structural gap
   (no cellular volume exclusion) rather than a fitted constant. Cheapest and
   most informative row available.
2. **The MCL diffusion coefficient, 74–130 µm²/s, PMIDs 9806526 + 17869669** — as
   a **FAIL** row on `effective_diffusivity`, carrying the MCL-vs-stroma caveat.
   Three independent cell lines make the target band real rather than nominal.
3. **Human intracranial-tumour IFP, 2.0 ± 2.5 mmHg, PMID 9062403** — the first
   external score of the pressure solver, in real mmHg, on the GBM flagship's own
   tissue. Must carry both the anti-oedema-therapy confound and the `_PV = 25.0`
   near-circularity note.
4. **Capecitabine tumour:plasma 1.76 and tumour:s.c. 1.11, PMID 12618890** — a
   free:free human breast ratio with a within-patient normal-tissue control.
   Score **capecitabine only**; DFCR/DFUR/5-FU are explicitly excluded as
   metabolites with an unmodelled local generation term.

**Do not adopt.** The porcine doxorubicin ratios (mixed denominators, no tumour);
the TMZ tumour-ECF numbers (no plasma leg); the Ekstrøm human case report
(n = 1, no recovery correction); anything from Brana 2014 (unretrievable).

**Obtain next, in priority order.**

1. **The Pruijn 2008 compound-by-compound `D_mcl` table** (doi:10.1071/CH08240).
   It is 71 measured diffusion coefficients by one method in one cell line
   spanning ~69-fold, against an engine whose `D_eff` is a function of MW and ECM
   fraction only. This is a **decisive, cheap test of the diffusivity model** and
   it is blocked purely on institutional access to one 2008 *Aust J Chem* paper.
   **This is the single highest-value obtainable measurement identified.**
2. A plasma AUC measured **in the same animals** as an intratumoural TMZ
   microdialysis experiment — which would convert PMID 23649683 from a
   near-miss into the free-interstitial glioma anchor the project lacks.
3. Institutional access to PMID 25213162 (Brana 2014) to see whether the single
   human doxorubicin biopsy carries a quantitative distance.

---

## 5. Verdict — can the transport physics be externally validated with public data?

**Partially, and more than this project currently believes — but not in the way
it has been trying.**

**Yes for the solver's inputs.** The two coefficients of the transport PDE —
`D_eff` and `k` — have both been directly measured, in the same units, by fitting
the same diffusion-reaction equation to tissue flux. Those measurements are
public, retrievable, absolute and out-of-sample, and §3 scores against them today
without one line of `tmesim/` changing. On first contact the engine **misses**:
hindrance 7× too weak, `D_eff` 2.2–3.9× too high, uncalibrated `k` 53–89× too
small. The pressure field can also be scored, in real mmHg, against a human
brain-tumour measurement.

**No for the solver's outputs.** Nothing public was found that measures a
**spatially resolved free-interstitial concentration field** in tumour tissue that
this engine could reconstruct. Microdialysis gives one point per probe with no
coordinates. MCLs give flux through a slab with no vasculature. Perivascular
immunofluorescence exists only for autofluorescent or antibody-detectable drugs —
which is why the anchor set is stuck on doxorubicin and mitoxantrone, two
chemically adjacent intercalators, one of which the calibration is fitted to. The
project's own stated need — a penetration profile for a weakly-bound,
deep-penetrating drug — **has no public answer.** MALDI-MSI is the modality that
could change this, and no open dataset with quantitative drug-versus-vessel
distance was found.

**So the honest framing.** Vivum can defensibly claim that its transport
*coefficients* are externally testable and currently fail against measurement in
identifiable, structural ways. It cannot yet claim that its predicted
concentration *fields* have been validated against anything, and nothing found in
this survey would let it. Closing the first gap is worth doing because it is real,
cheap, and points at a fixable modelling omission. Claiming it closes the second
would be the same category error that §0.1 catalogues three live instances of.
