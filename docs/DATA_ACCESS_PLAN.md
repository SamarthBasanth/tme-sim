# Data Access Plan — unlocking the validation frontier

*Written 2026-07-27. Fact-checked and updated 2026-08-12 (§0, §4, §8, §11). Companion to
`results/DATA_HUNT.md` (the 39-agent, 164-dataset sweep), `results/SESSION_AUDIT.md` (the
rigour audit), and `results/ranking_robustness.md` §8 (which names the decisive measurement).*

**This document adds nothing to the code and changes no result.** It is the access plan for
the four or five measurements that would move `tmesim` from "a transport model with a
quantified error budget" to "a model whose response prediction has been tested." Every one of
them is a partnership, purchase or application — not a download.

Outreach drafts live in [`docs/outreach/`](outreach/). **They are drafts for human review and
must not be sent unedited.**

> ## Status, in one block
>
> | | |
> |---|---|
> | **Axis that cannot be validated from open data** | **Response** — the delivery→kill link, in space. |
> | **Why** | No dataset anywhere, open or controlled, pairs a measured spatial drug distribution with a measured treatment response in the same tumour. 164 datasets swept 2026-07-27; landscape re-checked 2026-08-12 (§11) with the same result. For glioma + temozolomide there is no spatial drug map at all — METASPACE returns zero. |
> | **Is this a to-do item?** | **No.** It is a structural limit of the data that exists. No further code, search or compute in this repository removes it. |
> | **The single experiment that would close it** | The Agar/White adavosertib GBM platform — drug µM per pixel at 30 µm raster in ABTC phase-1 patient tissue, co-registered with phosphoproteomics and t-CyCIF response markers. **No accession.** Both published data-availability statements (2022, 2025) say "available on request". §4. |
> | **What it would take** | A collaboration request on human-subject tissue from a closed trial: IRB pathway + DUA + a counterparty institution. 6–18 months, real probability of no reply. |
> | **State of the request** | **Written, fact-checked, and UNSENT.** Four drafts in [`docs/outreach/`](outreach/). No person or group named in them has been contacted; none knows this project exists. |
> | **What is blocking each draft** | 01 (Agar/Santagata) and 03 (dbGaP): institutional affiliation, gate G4. 02 (Erasmus MC) and 04 (GSE324268 depositors): nothing but a human decision to send, plus writing the pre-registration they promise, which does not yet exist. |

---

## 0. The honest starting position

State this, unchanged, in every approach. It is the reason anyone would say yes.

| | |
|---|---|
| **What is built** | A patient-specific GBM microenvironment reconstruction (10x Xenium 319-gene and 480-gene panels, 4 Visium sections) feeding a mechanistic reaction–diffusion transport solver and a Hill/PD kill layer. All inputs measured; the synthetic-data era is deleted from the package. |
| **What was found by attacking it** | **11 distinct defects**, several mutually cancelling and invisible to a green 55-test suite (`results/SESSION_AUDIT.md`). Every headline number published before the audit was wrong in at least one way. |
| **What is externally validated** | **Delivery/partition, in human glioma**: scored against 125 tissue samples across 42 patients plus 27 PMID-anchored regional PK measurements. **~2× geometric-mean fold error** (1.995×), 6/6 drugs within 10×, worst case 6.94×. Regional contrast-enhancing-vs-non-enhancing ratios 1.19× / 5.87× predicted against 1.45× / 4.26× measured, after the additive barrier leak was implemented (L = 0.238). *Corrected 2026-08-12: this row previously read "separation 3.65× predicted against 4.66× measured" and "Reconstruction cross-platform Pearson 0.997". **Both are wrong as current claims** — 3.65× is the pre-fix separation (post-fix the model predicts 17.20×, an over-separation), and the 0.997 was **retracted** in `results/revalidation.md`. See §8.* |
| **What is NOT validated at the reconstruction level** | Cross-platform composition concordance is an **open question**, not a pass. Every reconstruction-level external validation this project has made has been retracted or fails on a held-out section; every transport-level one has held. That asymmetry is the most informative pattern in the audit and it is not flattering to the part of the model that is easiest to make look impressive. |
| **What was falsified** | The **oxygen PDE** as a value-add — beaten by plain distance-to-nearest-vessel (ρ −0.150 vs −0.113) and by its own undiffused source. The **`penetration_index`** metric — definitionally unable to discriminate drugs; retired. Both retractions are published in the repo. |
| **What is NOT validated** | **Response.** The delivery→kill link has never been tested in space, because no open dataset pairs a measured spatial drug distribution with a measured treatment response in the same tumour. Not one, anywhere. |
| **Stage** | Research prototype. Not clinically validated, not a device, no regulatory pathway, no patient has been or will be treated on its output. |

**The single sentence that defines this document:** the physics is scored and the response is
not, and the reason is access, not code.

### 0.1 The gate that applies to half this plan: institutional affiliation

Three of the routes below — dbGaP (which now also governs controlled TCIA imaging), the NCI
PDMR MTA, and HTAN Level 1/2 — independently require:

- an **eRA Commons account**, which an institution must create ("a principal investigator
  cannot create their own eRA Commons account");
- a **PI at tenure-track professor / senior scientist level** (dbGaP explicitly bars graduate
  students and postdocs from submitting);
- an **institutional Signing Official** who co-signs before NIH sees the request, plus a named
  IT Director attesting to NIH security practice;
- for PDMR, an **institutional Authorized Official** and, for PDX, an institutional ACUC
  protocol.

**An unaffiliated individual cannot pass these gates.** NIH publishes an
"unaffiliated/independent applicant" registration route, but it is written for *grant
submission* and says nothing about data access requests; a request in which the applicant is
their own Signing Official and own IT Director is exactly the conflict the requirement exists
to prevent. Treat it as unavailable.

**Consequence for sequencing:** the two targets that do *not* require an institution — the
Erasmus MC collaboration and the Agar/Santagata collaboration — are also the two highest-value
targets. Pursue them first. If either succeeds, the collaborating institution becomes the
route to everything else. That is the actual unlock.

### 0.2 What we can offer, honestly

Not money, not tissue, not a clinical claim. What exists and is transferable:

1. **The reconstruction pipeline** — Xenium/Visium → typed cells → vascular field → ECM/oxygen
   → transport solve, on real coordinates, mass-conserving to ~1e-13 on an irregular mask,
   grid-converged at dx = 12 µm.
2. **The validation harness** — `tmesim/validation/` (24 modules) is built to *score a model
   against a measurement and report where it fails*, including permutation tests, LOO
   cross-validation, and an explicit uncertainty/ranking-robustness Monte Carlo. It is
   reusable against any collaborator's data and it is designed to publish negative results.
3. **A pre-registered, falsifiable prediction.** Before seeing any new measurement we can
   commit, in writing and with a hash, to the predicted per-region exposure, the predicted
   drug ordering, and the predicted MGMT effect size. That is worth more to a collaborator
   than a post-hoc fit, and it costs us the option of quietly tuning.
4. **Co-authorship on their terms, or none.** Including the option that the modelling appears
   as a supplementary analysis, or not at all if it fails.

---

## 1. Prioritised target list

Ranked by (what it unlocks) × (probability of getting it) ÷ (time and cost). Full detail in
§2–§6.

| # | Target | Unlocks | Mechanism | Timeline | Cost |
|---|---|---|---|---|---|
| **1** | **Erasmus MC / GLIOscreen — per-patient ex-vivo IC50 for cisplatin + doxorubicin on a specimen that is also spatially profiled** | The **decisive** measurement named in `ranking_robustness.md` §8. Collapses the dominant remaining uncertainty in the top-1 drug call, and creates the world's first same-specimen spatial + ex-vivo-potency glioma dataset. | Collaboration email → MTA + Erasmus non-WMO/NWBTC review | 3–9 months to data | Low for us (compute + analysis). Assay is a marginal add-on to a screen they already run. Spatial run €10–20k if new. |
| **2** | **A paired MGMT readout on a spatially profiled section** — largely **already solved**, see §3 | Converts TMZ sensitivity from *unmodellable* to *measured*. Removes the "no intrinsic resistance modelled" flag on every TMZ number. | **GSE324268 is a free download today.** Residual ask: a methylation call on the same block. | Days (download) / weeks (methylation ask) | $0 / low |
| **3** | **Agar & Santagata (BWH) — adavosertib MALDI-MSI + t-CyCIF from the ABTC phase-1 GBM specimens** | The **only** existing measurement of drug concentration per pixel in human glioma with a co-registered molecular response. Would validate the delivery→kill link in space, in the right disease, for the first time. | Cold approach to the lab; then IRB + DUA on phase-1 patient tissue | 6–18 months, may be never | Low cash, high patience |
| **4** | **TCIA controlled imaging via dbGaP `phs004225` — Burdenko-GBM-Progression (RTDOSE + 651 response labels), Brain-Tumor-Progression (pre-computed perfusion)** | The structural `dose(x,y,z) → local control` test, at n = 180. Wrong modality (radiation, not drug) but the right shape. | **dbGaP DAR — requires institutional SO.** Blocked without affiliation. | 1–2 months with an SO; blocked otherwise | $0 |
| **5** | **HTAN + NCI PDMR** | *Substantially less than the brief assumed — see §6.* HTAN has **no GBM atlas**; PDMR has **no glioma models and no public drug-response data at all**. | dbGaP `phs002371` / PDMR MTA | Months | $0–$250/vial |

**Two of the five targets shrank on inspection, and one dissolved into a free download.** That
is the plan working: §7 is the register of what changed.

---

## 2. Target 1 — the decisive cheap experiment (Erasmus MC / GLIOscreen)

### 2.1 What it would let us validate that we currently cannot

`results/ranking_robustness.md` §8 states the problem precisely. The panel ranking is decided
by **delivery × potency**, and the top two candidates sit at opposite ends of both axes:
cisplatin has ~157× the mean AUC and ~49× the IC50 of doxorubicin. Cisplatin (0.241) and
doxorubicin (0.231) differ by ~4%, well inside plausible parameter uncertainty; under the
stress scenario P(cisplatin = rank 1) falls to 0.582.

> **Until a per-patient potency readout exists, the reportable output of this pipeline for the
> top of the panel is a SET, not a winner.**

A per-patient IC50 pair for cisplatin and doxorubicin collapses the widest prior in the Monte
Carlo — "the median glioma cell line" — from a population spread to a measurement error. It
also removes a tier asymmetry: doxorubicin's 1/3×–3× literature prior exists *purely* because
doxorubicin is absent from the GDSC glioma screen.

Second, and larger than the ranking question: **if the potency assay is run on a specimen that
is also spatially profiled, it is the first dataset of its kind in glioma.** Our sweep found
no published paper anywhere in which a glioma resection received both an ex-vivo drug
sensitivity screen and a spatial molecular assay (Visium / Xenium / CosMx / GeoMx / MERFISH /
MALDI-MSI). The closest is Shekarian et al., *Sci Adv* 2022 (PMID 35776791) — seven GBM
patients, 5-ALA-guided core and periphery into perfusion bioreactors, treated with
anti-CD47/nivolumab, then CODEX on the same explants. Spatial proteomics, immunotherapy arms,
not a dose-response screen. Two independent 2025 reviews (PMID 39953052, PMID 39799339) cite
Shekarian as the sole example of the design, which corroborates the gap.

The gap has been named by the field's own leaders. Yu KKH et al., *Nat Commun* 2025
(PMID 40295505, Agar lab) took stereotactic GBM cores through spatial transcriptomics,
t-CyCIF, MALDI-MSI and scRNA-seq **and derived PDX/PDCL from the same biopsy — then never
drug-screened them.** The paper says, verbatim: *"Having models that are matched to
multi-omics data from the same needle biopsy would allow for functional integration with
multi-omics results."* That is this experiment, written down by the group best placed to run
it, unexecuted.

### 2.2 Who to ask, and why them

**Dr Martine L. M. Lamfers**, PI, **GLIOscreen**, Department of Neurosurgery / Brain Tumor
Center, Erasmus MC Cancer Institute, Rotterdam.
Group page: <https://www.erasmusmc.nl/en/research/groups/glioscreen> · <https://glioscreen.com>
Corresponding-author email on PMID 37620410: **m.lamfers@erasmusmc.nl** (verified from three
independent sources — the PubMed affiliation string, her Erasmus MC group page, and PMID
40425813).

Why this group:

- They already run the exact assay. Ntafoulis I, Kleijn A, Ju J, …, **Lamfers MLM**. *Ex vivo
  drug sensitivity screening predicts response to temozolomide in glioblastoma patients and
  identifies candidate biomarkers.* **Br J Cancer** 2023;129(8):1327–1338. **PMID 37620410**,
  GEO **GSE232173**. 6-day CellTiter-Glo ATP viability on patient-derived GSC cultures; TMZ
  across 66 cultures; the % viability @100 µM readout out-predicted both AUC and IC50 against
  patient survival.
- We have already **used and audited that dataset** (`results/gbm_tmz_exvivo.md`,
  `results/mgmt_scalar.md`) and independently reproduced their MGMT finding from the raw
  counts (our ρ = +0.516, p = 0.034, n = 17; theirs ρ = +0.71). Approaching them with a
  completed, published-in-repo reanalysis of their own data is the credibility that makes this
  email answerable.
- They have the **Erasmus Brain Tumor Cell Culture Repository** — 350+ patient-derived brain
  tumour samples.
- Their follow-up (Posthoorn-Verheul C, …, Lamfers MLM, *npj Precis Oncol* 2025;9:156,
  **PMID 40425813**) reports GSC cultures from **114 consecutive high-grade glioma specimens**,
  **96% establishment success**, and a 20-agent screen returned inside a clinically relevant
  window for 16/18 samples. This is an operating, high-throughput pipeline, not a one-off.
- **Spatial capability is already inside the building, and already inside their author lists.**
  This is the finding that makes the ask realistic rather than aspirational:
  - **Johan M. Kros** (Pathology & Clinical Bioinformatics) is a **co-author on PMID 37620410
    itself**, and his lab with **Dr Dana A. M. Mustafa** runs **NanoString GeoMx DSP on brain
    tissue** in-house (PMID 37061716, PMID 40621453, PMID 36685537).
  - **Dr Miao-Ping Chien** (Molecular Genetics, Erasmus MC / Oncode) is a **co-author on their
    own 2025 culture paper (PMID 40425813)** and builds spatially resolved single-cell
    technologies (PMID 37626852, PMID 35273957).
  - **van Leenders' group** (Pathology) has published **Visium** on radical prostatectomy
    specimens (PMID 42236150) — in-house spatial transcriptomics on surgical tumour material.
  - **PARTS** (Pathology Research and Trial Service,
    <https://www.erasmusmc.nl/en/research/groups/pathology-parts>) explicitly serves
    researchers "within and outside the Erasmus MC." Secondary sources report a 2024
    acquisition of 10x spatial transcriptomics; *whether that is Visium, Xenium or both is not
    verified* and is worth one question.

  So the ask is **not** "acquire a capability you don't have." It is "pair two things you
  already do, on one block, for the first time."

### 2.3 The exact ask

One specimen, minimum. Ideally 3–5.

1. From a GBM resection already going into the GSC pipeline, retain an **adjacent FFPE or OCT
   block** from the same tumour region.
2. Run their **standard 6-day CellTiter-Glo dose-response** on the derived culture, adding
   **cisplatin** and **doxorubicin** to whatever panel they run. A full curve (not the
   single-dose format of the 2023 paper) is what the decisive test needs — an IC50 per drug
   per patient.
   *Open question for them:* are cisplatin and doxorubicin already in the 20-agent panel of
   PMID 40425813? If yes, the incremental cost is zero.
3. Run **any** spatial assay on the retained block — Xenium, Visium, GeoMx or CosMx. The
   pipeline consumes all four. Xenium/Visium are preferred (single-cell or spot coordinates);
   GeoMx ROIs are usable but coarser.
4. Provide the routine clinical covariates already recorded: IDH, **MGMT methylation status**,
   age, extent of resection, and — if available and permissible — PFS/OS.

**Cheapest viable variant, if a new spatial run is not fundable:** they run the two extra
drugs on cultures derived from specimens for which a spatial section *already exists*
somewhere in Erasmus MC's archive. The retrospective version costs the assay only.

### 2.4 Access mechanism and timeline

Dutch route, factual:

- Research using **only coded, already-collected leftover tissue**, with no procedure imposed
  on the patient, is **not WMO-liable by definition**, so no accredited-METC review under the
  WMO is required (CCMO: <https://www.ccmo.nl/onderzoekers/aanvullende-informatie-over-bepaalde-soorten-onderzoek/niet-wmo-onderzoek>;
  Health-RI ELSI: *"Een onderzoek met (alleen) nader of secundair gebruik van beschikbare
  gegevens of lichaamsmateriaal valt per definitie niet onder de WMO."*).
- It must instead be covered by the **"geen bezwaar" (no-objection/opt-out) system** for
  residual tissue and supplied **coded**, per the COREON **Gedragscode Gezondheidsonderzoek
  2022** (<https://www.coreon.org/gedragscode-gezondheidsonderzoek/>), §5.5.1; statutory anchor
  art. 7:467 BW (WGBO).
- In practice you still need (a) review by Erasmus MC's non-WMO biobank committee, the
  **NWBTC**, which the Central Biobank regulation states "is part of the METC," plus
  scientific/collection-owner sign-off; and (b) a signed **MTA/MDTA** before anything leaves.
  Joint Dutch UMC MDTA template:
  <https://elsi.health-ri.nl/sites/elsi/files/2020-07/20.30817%20OenO%20Template%20-%20Material%20and%20associated%20Data%20Transfer%20Agreement.pdf>
- External collaborators submit via **Research Manager**
  (<https://erasmusmc.myresearchmanager.com/>); an account is requested from
  **metc@erasmusmc.nl**.
- Note **§8.3.3** of the Central Biobank regulation: **any commercial involvement routes
  through the Erasmus MC TTO** rather than the Central Biobank's non-commercial MTA. If this
  project ever becomes commercial, say so at the outset — discovering it later is the fastest
  way to lose the collaboration.
- The **Wet zeggenschap lichaamsmateriaal (WZL, bill 35844) is not in force** (still with the
  Tweede Kamer as of early 2026) and would move body material from opt-out to statutory
  consent. Re-check before relying on the current regime.

**Realistic timeline:** 2–6 weeks to a reply or silence; 1–3 months to agree a design;
2–4 months for NWBTC + MTA; then the assay. **3–9 months to first data**, and that assumes
they are interested. Data transfer is more likely to be *us receiving derived tables* (IC50s,
a spatial matrix) than raw tissue, which keeps the MTA simple.

### 2.5 What we provide, and cost

We provide: a written pre-registered prediction before unblinding; the reconstruction pipeline
run on their section at no cost; the validation harness; and full analysis code. We ask for
co-authorship only if the modelling contributes; we will publish a negative result if the
prediction fails, and we will say so in the initial email so it is not a surprise later.

Cost to us: compute and time only. Cost to them: the two extra dose-response curves are
marginal on an existing screen. A *new* spatial run is the real expense — as a public
reference point, core-facility list prices for Xenium are roughly $11k (v1 base) to $17k
(Prime 5K) per 2-slide run (e.g. UAB's published rate card), and a 100-gene custom add-on adds
roughly $3.5k–$7k. **These are third-party list prices, not a quote, and Erasmus MC's internal
rates will differ.**

Draft: [`docs/outreach/02_erasmus_lamfers.md`](outreach/02_erasmus_lamfers.md).

---

## 3. Target 2 — MGMT on a spatially profiled section (mostly already solved)

### 3.1 What it unlocks

From `results/gbm_tmz_exvivo.md` and `results/mgmt_scalar.md`: MGMT is absent from both Xenium
panels we use (319-gene brain, 480-gene IO; 739-gene union, MGMT present = **False**), **no
panel gene surrogates for it** (0/396 pass both FDR bars; the panel's MGMT-correlation base
rate is 14.1% against 12.3% genome-wide, i.e. no enrichment at all), and **GDSC cannot resolve
it** (ρ = −0.068, p = 0.71, n = 33 — GDSC dosed TMZ at 0.03–30 µM against an effective range of
25–300 µM, so its TMZ arm has no dynamic range). So `SENSITIVITY_DETERMINANTS['temozolomide']
is None` and the model's between-patient variance in TMZ sensitivity is **exactly zero by
construction**.

Measured cost of that forfeit: MGMT explains **27%** of the rank variance in measured ex-vivo
response (50% by the paper's own coefficient) across an **81-percentage-point** viability
spread and a **700×** spread in implied IC50, separating cohorts whose median OS is 27.8 vs
12.3 months. Supplying a measured MGMT moves predicted TMZ kill from 0.1151 to 0.023–0.040 —
a 2.9×–5.0× change, i.e. the dominant uncertainty in TMZ's absolute number.

### 3.2 The mechanism exists — and so, already, does the data

The brief anticipated a custom-probe purchase. Confirmed mechanism, for the record:

- The tool is the **10x Xenium Panel Designer** (<https://www.10xgenomics.com/support/software/xenium-panel-designer/latest>),
  "add-on custom" configuration, **up to 100 genes** on a pre-designed base panel. Advertised
  turnaround **4–6 weeks** (<https://www.10xgenomics.com/platforms/xenium/panels>).
- **MGMT is designable, verified against 10x's own reference files** — `ENSG00000170430`
  appears in both `human_and_mouse_2020-A-ref-yesprobe-genes_v1assay.csv` and the Prime
  equivalent, and in neither exclusion list.
- Two premises to drop: the canonical MGMT transcript is **4,739 nt** across 5 exons, so
  transcript length is a non-issue against a 40-base probe requirement; and the **"MGMT
  pseudogene" does not appear to exist** — HGNC returns one MGMT-family symbol, Ensembl reports
  zero paralogues, NCBI has a single human gene record. Paralog homology is 10x's main
  sequence-level exclusion reason and it does not apply.
- **MGMT is already on the off-the-shelf Xenium Prime 5K Human Pan Tissue panel** (5,001 genes,
  3 probe sets) — no custom design needed at all if Prime 5K is the base.

**But the ask is smaller than that, because the data already exists and is free:**

| Dataset | What it is | MGMT | Access |
|---|---|---|---|
| **GEO GSE324268** (Hölzel lab, Bonn; public Mar 2026) | **Xenium, human primary GBM (`GBM_tumor07`, 851,509 cells)** + a GBM-organoid-on-brain-slice sample. Panel = Xenium Human Brain v1.1 (266) **+ a custom 100-gene add-on literally named "Glioblastoma Panel"** (design_id ECQNGH) = 366 targets. | **YES — 8 probe sets. Verified in the data: 15,191 counts across 14,675 cells (1.7%).** | **Free direct download, no login, no DUA** |
| **GEO GSE305393** | CosMx 6K Discovery, **human glioma TMA, 18 patients**, 39 FOVs, 32,953 cells | **YES**, verified in the expression matrix (~2.9% of cells) | Free direct, per-file |
| 10x Visium HD Human Brain Cancer (6.5 mm FFPE) | Whole-transcriptome probe set v2.0 at 2/8/16 µm bins — but the tissue is **anaplastic oligodendroglioma WHO III, not GBM** | YES, 3 probes | Free direct |
| GSE318562 / Zenodo 10.5281/zenodo.15014838 (Petti lab) | Visium, 22 samples, 9 IDH-WT GBM, core-to-edge | Present but **near the floor** (46 counts in 46/1,787 spots) | Free, CC BY |

**Therefore the residual ask is small and specific:**

1. **Do this week, no permission needed:** re-run the MGMT scalar analysis against **GSE324268**
   in place of Visium pseudobulk. This is a *targeted 8-probe* measurement on 851k single cells,
   not a whole-transcriptome detection-floor artefact — it is the first honest test of whether
   MGMT is a resolvable spatial field or only a per-section scalar. Our current position
   (`mgmt_scalar.md` §2) is that Visium can only support a scalar; this dataset can falsify or
   confirm that.
2. **Ask the Hölzel lab** (GSE324268) for the MGMT **promoter-methylation** status of that
   specimen, and for the add-on panel's design rationale. A short, specific, low-burden email.
3. **Only if 1 and 2 fail** should a custom panel be purchased — and if it is purchased,
   **use Prime 5K and skip the custom design entirely.**

### 3.3 The caveat that must travel with every MGMT claim

Every one of these platforms measures MGMT **transcript abundance**. The clinically actionable
biomarker in GBM is **promoter methylation**, which no probe- or imaging-based spatial platform
measures. Expression is a correlate, not the assay. Our own analysis actually argues for
expression on the merits — in the Erasmus cohort MGMT *expression* was prognostic for OS
(HR 1.85, p = 6.1e-05, externally replicated at HR 1.21) while *promoter methylation* was not
significant in the same cohort (p = 0.16) — but that argument must be made explicitly every
time, not assumed. And at 1.7–2.9% per-cell detection, single-cell MGMT calls are noisy;
anything defensible aggregates over regions or populations.

Draft: [`docs/outreach/04_holzel_gse324268.md`](outreach/04_holzel_gse324268.md).

---

## 4. Target 3 — Agar & Santagata, adavosertib MALDI-MSI in phase-1 GBM

### 4.1 What it unlocks — the whole thing, at once

This is the only existing measurement of **drug concentration per pixel in human glioma tissue
with a co-registered molecular response.** It is the experiment `tmesim` would design if it
could design one.

Full citation, confirmed: Lopez BGC, Kohale IN, Du Z, …, Ligon KL, Wen PY, Alexander B, Lee EQ,
**Santagata S**, Sarkaria J, **White FM**, **Agar NYR**. *Multimodal platform for assessing drug
distribution and response in clinical trials.* **Neuro-Oncology** 2022 Jan 5;24(1):64–77.
doi:10.1093/neuonc/noab197. **PMID 34383057, PMCID PMC8730776.** Funded by NIH **U54 CA210180**
(MIT/Mayo Physical Sciences Oncology Center for Drug Distribution and Drug Efficacy in Brain
Tumors).

It reports measured adavosertib at **1–19 µM per pixel** in GBM84 with intratumoral **CV
20.9–85.6%**, in GBM PDX *and* in specimens from **6 GBM patients enrolled in NCT01849146** —
the Adult Brain Tumor Consortium phase-1 study of adavosertib + radiation + temozolomide, which
included a surgical window-of-opportunity arm in recurrent GBM specifically to measure
intratumoral drug levels. (Clinical companion: *Clin Cancer Res* 2025;31(6):983, PMC11913579.)

With it, and nothing else, we could test:

- **predicted vs measured intratumoral drug CV** in human GBM — the first direct test of the
  spatial delivery axis in the right disease;
- **predicted local exposure vs measured local phosphoproteomic response** — the first test of
  the delivery→kill link in space, in any tumour;
- whether our BBB additive-leak term (L ≈ 0.238, derived in `glioma_regional_pk.md` from three
  independent studies) predicts the *within-section* variation, not just the region means.

### 4.2 The access reality — no accession, and no institution to route through

The data-availability statement, verbatim from PMC:

> "The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium
> via the PRIDE partner repository with the dataset identifier PXD018782 and
> 10.6019/PXD018782. **The MALDI-MSI data are available upon request.** The processed t-CyCIF
> montage images were uploaded to Omero."

**And the trial report says the same thing again, three years later.** Lee EQ, Alexander BM,
Romo CG, Supko JG, **Agar NYR**, …, **Santagata S**, …, **Wen PY**. *Phase I Study of
Adavosertib with Radiotherapy and Temozolomide in Newly Diagnosed Glioblastoma and Intratumoral
Drug Levels in Recurrent Glioblastoma.* **Clin Cancer Res** 2025;31(6):983–992.
**PMID 39820282, PMC11913579**, doi:10.1158/1078-0432.CCR-24-2311. Verified from PMC
2026-08-12; its data-availability statement, verbatim:

> "The data that support the findings of this study are available on request from the
> corresponding author. The data are not publicly available due to privacy restrictions."

This paper is the better primary citation for the outreach draft, and it adds three facts the
2022 paper does not:

- **The MSI quantitation is absolute, not arbitrary ion counts.** MALDI-FT-ICR MSI on a 9.4 T
  instrument at a **30 µm laser raster**, quantified against a gelatin/brain-homogenate mimetic
  calibration series (AZD1775 at 1, 2, 5, 10, 15, 20, 100 µM plus a blank). Measured range
  **970 nM–7,270 nM** across **10 specimens from 6 patients** (3 specimens from 2 patients at
  200 mg — 2 CE, 1 NCE; 7 from 4 patients at 425 mg — 3 CE, 4 NCE).
- **Their two modalities disagree on region, and that disagreement is exactly what this model
  outputs.** LC-MS/MS on homogenates gives mean intratumoral:plasma 4.18 ± 3.36 in
  contrast-enhancing vs 0.74 ± 0.63 in non-enhancing tissue (200 mg cohort; 7.30 ± 3.13 vs
  2.62 ± 3.00 at 425 mg), while the MSI shows "no apparent correlation between drug
  concentration and contrast-enhancing versus non-contrast enhancing tumour". The paper
  attributes this to specimen allocation and MSI sensitivity.
- **The local drug↔response coupling is reported but not deposited.** "Regions of high
  adavosertib levels showed increased levels of DNA damage response as measured by Nibrin (NBN)
  pS343." That sentence is the axis this project cannot validate, measured, in the right
  disease, and unavailable.
- **The trial is closed.** NCT01849146 enrolled August 2013 – May 2020, terminated 2020-05-05,
  71 patients, 12 through the intratumoral-drug-distribution arm, 2 through microdialysis (with
  catheter failures). Nothing further will be collected.

One of three modalities is self-serve:

- **PRIDE PXD018782 is public, CC0, no auth** — <https://www.ebi.ac.uk/pride/archive/projects/PXD018782>,
  FTP <https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/10/PXD018782/>. **But our
  verification pass established it contains no quantified phosphosite response** — all 80 files
  are raw `.raw` spectra or peptide identifications, zero reporter-ion/abundance/ratio fields.
  To make it interpretable you additionally need the **TMT channel → experimental arm design
  map**, which is not deposited and must come from the White lab at MIT.
- **MALDI-MSI: not deposited anywhere.** Request only.
- **t-CyCIF: "uploaded to Omero" with no accession, no URL.** Not found on any public Minerva
  story, cycif.org/data, or HTAN. Assume closed.

Contacts, as publicly listed (do not treat as an invitation):

| | |
|---|---|
| Agar lab | **Surgical Molecular Imaging Laboratory (SMIL)**, BWH Neurosurgery — <https://agarlab.bwh.harvard.edu/>. **Web contact form only**; no email or phone printed on the lab site. Corresponding address in the paper: nathalie_agar@dfci.harvard.edu |
| Santagata | <https://sysbio.med.harvard.edu/sandros> (HMS Systems Biology) prints ssantagata@partners.org and assistant cbuzzell@bwh.harvard.edu. The paper prints ssantagata@bics.bwh.harvard.edu. **santagata.bwh.harvard.edu is currently down.** |
| Laboratory of Systems Pharmacology | <https://labsyspharm.org/about/contact/> |

**Two structural problems to plan around.** First, **the PSOC is defunct** — `psoc.mit.edu`
returns NXDOMAIN, U54 CA210180 ran 2016-08-29 to 2022-07-31, and NCI lists it under "Past
Projects." There is no PSOC intermediary to route through. (PS-ON itself is still active and
has an **Associate Membership** mechanism for outside investigators —
nci.physics@mail.nih.gov — but that is a network-participation route, not a data route.)
Second, SMIL is **not a fee-for-service core**; there is no rate card and no service menu, only
a stated openness to collaboration and a web form.

### 4.3 Realistic assessment

Phase-1 patient tissue means **IRB and DUA, not a link**. A US IRB reliance or a data-use
agreement requires a counterparty institution, which returns us to §0.1. The honest framing:
this is a **6–18 month collaboration ask with a real probability of no reply**, and it should
be sent anyway, because the downside is one email and the upside is the only dataset in the
world that closes the loop.

**Sequence it after Target 1.** If Erasmus MC (or any institution) is engaged, the Agar
approach can be made jointly and has an institution behind it. Sent cold and unaffiliated, it
is a long shot; sent with a collaborating PI, it is a normal academic request.

**A softer, higher-probability first ask** — and the one the draft leads with — is *not* the
patient MSI. It is:

1. the **TMT channel → arm design map** for PXD018782, which is low-sensitivity, involves no
   patient tissue, and would let the community actually use an already-public deposit; and
2. the **PDX arm** of the MALDI-MSI (GBM22/GBM84), which is not human-subject data and
   therefore not IRB-encumbered.

Getting either establishes contact and demonstrates good faith before anything touching
NCT01849146 is discussed.

Draft: [`docs/outreach/01_agar_santagata.md`](outreach/01_agar_santagata.md).

---

## 5. Target 4 — TCIA controlled collections (now dbGaP)

### 5.1 The mechanism changed and the old instructions are wrong

**The "TCIA Restricted License" no longer exists.** TCIA ceased hosting controlled-access
datasets on **2025-07-08** and moved access administration to dbGaP and hosting to the NCI
Cancer Research Data Commons, to comply with EO 14117 / 28 CFR Part 202 / NIH Guide Notice
NOT-OD-25-083. Authoritative page:
<https://www.cancerimagingarchive.net/nih-controlled-data-access-policy/>. The old
`TCIA-Restricted-License-*.pdf` email-the-signed-form workflow is superseded.
`help@cancerimagingarchive.net` is still valid but is now the wrong first stop — it handles
Data Retriever and manifest problems, not access.

**Current route:** a single dbGaP DAR against **`phs004225` — "Radiology Imaging Data
Containing Faces from The Cancer Imaging Archive"** (v1.p1, released 2026-01-16). **One request
covers all affected TCIA face datasets**; you do not apply per collection. On approval: log in
at <https://nci-crdc.datacommons.io/login> → Profile → Create API Key, download the collection
manifest, and point the TCIA Data Retriever at the JSON key.

Requirements are the standard NIH DAR — eRA Commons, PI at professor/senior-scientist level,
institutional SO co-signature, Research Use Statement, DUC, 1-year renewable term. **See §0.1:
this is blocked without an affiliation.** Hard exclusion: applicants at institutions in China
(incl. HK/Macau), Russia, Iran, North Korea, Cuba and Venezuela are prohibited outright.
DAC review itself averages ~2 weeks; end-to-end 2–6 weeks with a cooperative SO.

Scale of the change: of 236 TCIA collections, 180 are Public, **50 are "Limited"** and 6 are in
transition. Essentially every brain collection with intact face data is now Limited — including
TCGA-GBM, TCGA-LGG, CPTAC-GBM, REMBRANDT, IvyGAP, RHUH-GBM, GLIS-RT, ACRIN-FMISO-Brain and
RIDER NEURO MRI.

### 5.2 Burdenko-GBM-Progression — and the part that is already free

<https://www.cancerimagingarchive.net/collection/burdenko-gbm-progression/> · DOI 10.7937/E1QP-D183

- **180 patients**, 645 studies / 4,956 series / 344,003 images, 131 GB.
- Modalities **CT, MR, RTPLAN, RTDOSE, RTSTRUCT** — RTDOSE DICOM is present, giving a true 3D
  delivered-dose field. Every RTSTRUCT has GTV/CTV/PTV plus 10 OAR contours; a subset has GTV
  annotated on **follow-up** studies. MRI→CT registration files are **not** provided.
- **The response labels are already open.** The clinical CSV is separately **CC BY 4.0** and
  downloadable with no auth:
  <https://www.cancerimagingarchive.net/wp-content/uploads/Burdenko-GBM-clinical-data.csv> —
  columns `Response_1st_fup` … `Response_10th_fup` with values
  `response`/`stable`/`progression`/`pseudoprogression`, plus IDH1/2 (110 pts), MGMT (92 pts),
  and surgery/RT/follow-up/death dates.

**What the gated half unlocks:** the structural `dose(x,y,z) → local control` test at n = 180.
Radiation dose is not drug concentration — the mechanism of cell kill, the time course and the
delivery physics all differ — so this validates the *shape* of the model's spatial
dose→response coupling and not its pharmacology. That distinction must be in the caption of any
figure that uses it, or the claim is dishonest.

### 5.3 Brain-Tumor-Progression

<https://www.cancerimagingarchive.net/collection/brain-tumor-progression/> · DOI 10.7937/K9/TCIA.2018.15quzvnb

**20 patients**, 40 studies, 8,798 images, 3.16 GB — the only TCIA GBM collection shipping
**pre-computed** DICOM perfusion maps (ADC, normalised CBF, normalised rCBV, standardised rCBV
from GRE-EPI DSC with contrast preload), plus T1w pre/post, FLAIR, T2w, and binary tumour
masks, all co-registered to T1+C. Two timepoints per patient: within 90 days after CRT, and at
progression. Small enough to actually work with.

Caveat: 40 studies = exactly 2 per patient, so **every patient is a progressor by
construction** — no control arm and no tabular outcome file. And the standing caveat from
`DATA_HUNT.md` §6.4 applies: **contrast-agent perfusion is not drug delivery.** Gd is ~550 Da,
hydrophilic, BBB-excluded, bolus-administered. Nothing in the open imaging world measures a
therapeutic in tissue.

### 5.4 QIN-GBM-Treatment-Response — confirmed trap

<https://www.cancerimagingarchive.net/collection/qin-gbm-treatment-response/>

Confirmed verbatim from the collection page: it contains *"'double baseline' multi-parametric
MRI images"* and *"the patients were scanned after surgery but prior to the start of therapy,
typically 2–5 days apart."* **Both scans are pre-therapy. There is no response outcome data of
any kind.** The name refers to the QIN programme's goal, not the data's content. Its only
legitimate use is quantifying the Ktrans/rCBV measurement noise floor. 54 patients, 589,314
images. **Do not request it expecting response labels.**

### 5.5 The open substitutes — do these first

Both are CC BY 4.0, no application, and together cover most of what Burdenko offers:

| Collection | n | DOI | Why |
|---|---|---|---|
| **MU-Glioma-Post** | 203 | 10.7937/7K9K-3C83 | Longitudinal post-treatment glioma MR + segmentations; clinical CSV carries date and characterisation of first and subsequent progression — **Burdenko's labels, open** |
| **CFB-GBM** | 264 | 10.7937/v9pn-2f72 | Pre/post-treatment MRI **plus radiotherapy plans including RTDOSE** — the closest open analogue to Burdenko's RT content |
| UCSD-PTGBM | 178 | 10.7937/fwv2-dt74 | Already characterised in `DATA_HUNT.md` T1.6 |
| UPENN-GBM | 630 | 10.7937/TCIA.709X-DN49 | Largest open GBM cohort; `DATA_HUNT.md` Tier 2 #10 |

**Recommendation: run the dose→response experiment on CFB-GBM + MU-Glioma-Post first.** If it
produces a signal, the Burdenko DAR becomes a well-motivated replication request with a
concrete Research Use Statement — which is also what makes a DAR likely to be approved. If it
produces nothing, the DAR was never worth the affiliation problem.

---

## 6. Target 5 — HTAN and NCI PDMR (substantially smaller than assumed)

Both were carried forward from `docs/real_data/drug_response.md`, which proposed:
*reconstruct HTAN GBM tissue → predict → score against paired PDMR/organoid measured kill.*
**Checking the two repositories directly, that plan is not currently executable, for reasons
that have nothing to do with access.**

### 6.1 HTAN — there is no GBM atlas

- dbGaP accession is **`phs002371`** (currently v9.p2, 1,894 consented subjects), a single
  umbrella study for the whole network — no per-centre accessions.
  <https://docs.humantumoratlas.org/data_access/db_gap/>
- **The only primary brain-tumour atlas is HTA202 (Caltech + UCLA Neurosurgery), and it is
  low-grade glioma progression profiled by seqFISH** — <https://humantumoratlas.org/center/hta202>.
  HTA8 (MSKCC) covers CNS **metastasis**, not glioma. WUSTL, Vanderbilt, Duke, HTAPP and OHSU
  cover prostate, colorectal, breast DCIS and others — **none covers brain.**
- **No brain/glioma/GBM data appears in any HTAN release through Release 7 (26 Nov 2025)**
  (<https://humantumoratlas.org/data-updates>). HTA202 is a Phase 2 award early in its
  submission cycle.
- Separately worth knowing: **HTAN's imaging and Level 3/4 data are largely open CC BY 4.0**
  and need no dbGaP at all, via the portal (<https://humantumoratlas.org/explore>), Synapse,
  and the NCI Imaging Data Commons. The dbGaP gate is on Level 1/2 raw sequence. One
  inconsistency to check per-file: the HTAN data-sharing paper describes controlled access as
  covering "level 1 and 2 genomic and transcriptomic data, as well as level 2 imaging data,"
  while the portal describes imaging as open.

**Action: monitor, do not apply.** Set a reminder to re-check HTA202's release status in
6 months. There is nothing to request today.

### 6.2 NCI PDMR — no glioma models, and no public drug-response data

Portal moved: `pdmr.cancer.gov` now redirects to
<https://dctd.cancer.gov/drug-discovery-development/reagents-materials/pdmr>.

- **The public data is genuinely frictionless** — the database at
  <https://pdmdb.cancer.gov/web/apex/f?p=101:1> needs no login (1,318 patients, 3,242
  distribution lots, CSV export on every report), and bulk omics sit on **anonymous FTP** at
  `ftp://dctdftp.nci.nih.gov/pub/pdm/` (~65,000 files, WES and RNA-seq FASTQ, RSEM TPM, VCF).
  *Correction to a premise in our earlier notes: PDMR raw FASTQ is **not** dbGaP-controlled.*
- **There is no public drug-response data.** "Preclinical Drug Study Results" is still listed
  as *In Development*. Mechanically the data flows inward — MTA Article 9 obliges recipients to
  report drug sensitivity and Article 10 lets NCI republish after a year — but nothing has
  surfaced.
- **CNS content across all 1,318 patients is three cases**: one anaplastic ependymoma
  (PDX + organoid), one paediatric pilocytic astrocytoma (PDC only), one neuroblastoma.
  **Zero glioblastoma. Zero high-grade glioma.**
- MTA route, for completeness: PI signs the Request Form, the institution's Authorized Official
  signs the MTA, "modifications to NCI's MTA Terms of Agreement will not be considered."
  Pricing is $0 NCI / $250 academic / $2,500 commercial per cryopreserved PDX or organoid vial,
  recipient-paid shipping, 5 years of annual reporting. Turnaround **not published**. Contact
  NCI_PDM_Repository@mail.nih.gov.

**Action: do not apply.** Retained here only so that the next person does not re-spend the
time, and as the worked example in the generic application skeleton.

### 6.3 What replaces this target

The Erasmus MC ask (§2) *is* the replacement, and it is strictly better: a real institution,
a real GBM biobank of 350+ models, a validated screen, and spatial capability in the same
building. The HTAN+PDMR plan was two half-datasets from two repositories that would have had to
be joined across different patients; Target 1 is one specimen.

Draft skeleton for whenever a controlled-access application *is* the right tool:
[`docs/outreach/03_controlled_access_skeleton.md`](outreach/03_controlled_access_skeleton.md).

---

## 7. What we can do with NO new access

This plan must not be purely blocked on other people. Ranked by value, all executable now with
data already on disk or freely downloadable.

1. **Re-run the MGMT analysis against GSE324268.** Free download, human primary GBM, 851,509
   cells, MGMT measured with **8 dedicated probes** rather than inferred from detection-floor
   whole-transcriptome counts. This is a direct test of the claim in `mgmt_scalar.md` §2 that
   MGMT can only be a per-section scalar. **Highest-value item in this section** — it either
   upgrades MGMT to a genuine spatial field or confirms the scalar, and either answer is
   publishable in-repo. Also gives a second GBM section with a 366-gene panel to reconstruct.
2. **Finish the three interrupted external validations** already listed in `SESSION_AUDIT.md`
   §6.4 — GIST imatinib (6 matched human tumour/normal pairs, per-patient plasma boundary
   condition, and a measured tumour:plasma range of 0.02–0.72 that systemic exposure explains
   none of: r = 0.03, p = 0.93); the Moncayo HGSOC explants (376,736 MSI pixels at 20 µm with a
   heavy-labelled internal standard on every one, three compounds, plus the pre-specified
   compound-specificity falsification test — weak bases should show accumulation heterogeneity
   and olaparib should not); and GSE232173. All three are on disk.
3. **Run the dose→response shape test on the OPEN TCIA collections** — CFB-GBM (RTDOSE, n=264)
   and MU-Glioma-Post (progression labels, n=203). This is §5.5: do the experiment on open data
   first, and let the result decide whether the Burdenko DAR is worth the affiliation problem.
4. **Score against Burdenko's clinical CSV alone.** The 651 per-follow-up response labels,
   IDH, MGMT and the treatment timeline are CC BY and downloadable today. Anything that does
   not need the pixels can be done now.
5. **Finish the two interrupted internal validations** — cell state vs measured MKI67, and
   simulated oxygen vs measured VEGFA/HIF1A. There is a partial result on the table (VEGFA is
   the single most anti-correlated gene of all 480 against simulated oxygen; HIF1A is exactly
   null) and it needs completing before either is claimed. Note this bears directly on the
   falsified oxygen PDE (defect #11): if hypoxia is better taken from measured markers, that
   should be implemented, not just observed.
6. **Pre-register the predictions.** Before any collaboration data arrives, commit — file,
   timestamp, git hash — to the predicted cisplatin/doxorubicin ordering, the predicted
   per-region exposure ratios, and the predicted MGMT effect size. This costs nothing, takes an
   afternoon, and is the difference between "we predicted it" and "we fitted it." It is also
   the single most persuasive thing that can be put in front of a prospective collaborator.
7. **Second penetration anchor** (`ranking_robustness.md` §8 item 2). α is currently one number
   from one study. Mining a second anchor for any panel drug other than doxorubicin from the
   existing PK compendium would let α be fitted rather than assumed-transferable, shrinking the
   widest transport prior. Affects levels more than order, so it ranks below the potency ask —
   but it needs no one's permission.
8. **Publish the negative results properly.** Two components falsified, eleven defects, and a
   documented `QIN-GBM-Treatment-Response` trap. A short, honest write-up of what does not work
   is itself a credential when approaching Targets 1 and 3, and costs no access.

---

## 8. Register of corrections

Assumptions carried into this plan that direct checking overturned. Recorded so they are not
re-propagated.

| Assumption | Finding |
|---|---|
| "MGMT is absent from every Xenium panel, so a custom add-on probe is needed." | True of the v1 panels we use, but **MGMT is on the off-the-shelf Prime 5K panel (3 probe sets)** and — more importantly — **GSE324268 is a free, already-public human GBM Xenium dataset with MGMT on an 8-probe custom add-on, 851,509 cells.** The purchase is not needed. |
| "HTAN GBM spatial, gated behind dbGaP + HTAN DUA." | **HTAN has no GBM atlas.** The only primary brain-tumour atlas is HTA202, low-grade glioma by seqFISH, with no brain data in any release through Nov 2025. |
| "NCI PDMR organoid/PDX drug screens, via MTA." | **PDMR has zero glioblastoma models and no public drug-response data at all** (3 CNS cases total across 1,318 patients; "Preclinical Drug Study Results" still *In Development*). Its raw FASTQ is on anonymous FTP, not dbGaP. |
| "TCIA restricted-data request form, weeks." | **The TCIA Restricted License was retired 2025-07-08.** Access is now a dbGaP DAR against `phs004225` covering all face-containing collections at once — which imports the institutional-SO requirement. |
| "Burdenko is gated." | **Half of it is not.** The 651 response labels, IDH, MGMT and the treatment timeline are CC BY 4.0 and downloadable now; only the pixels are controlled. |
| "The Agar/White platform sits inside an active PSOC we could route through." | **The PSOC is defunct** — U54 CA210180 ended 2022-07-31, `psoc.mit.edu` is NXDOMAIN, NCI lists it under Past Projects. There is no intermediary; SMIL is not a fee-for-service core and its only advertised contact is a web form. |
| "Erasmus MC does the ex-vivo screen; someone else would have to do the spatial." | **Spatial capability is inside Erasmus MC and already inside their author lists** — Kros (GeoMx, and a co-author on PMID 37620410 itself), Chien (spatial single-cell, a co-author on their 2025 paper), van Leenders (Visium on surgical specimens), and PARTS, which serves external researchers. |
| "MGMT probe design may be blocked by a pseudogene or transcript length." | Neither. Canonical MGMT is 4,739 nt across 5 exons; HGNC/Ensembl/NCBI show **no MGMT pseudogene and zero paralogues**; `ENSG00000170430` is on 10x's designable list for both v1 and Prime. |
| "QIN-GBM-Treatment-Response contains response data." | Confirmed **trap** — a 2–5 day double-baseline test-retest study, both scans pre-therapy. |

---

## 9. Decision gates

| Gate | Test | If it fails |
|---|---|---|
| **G1** — 6 weeks after the Erasmus email | Any substantive reply? | Widen: Zurich (Weller/Snijder pharmacoscopy + Tumor Profiler, PMID 39304781, PMID 40425842 — the group with the exact same-specimen capability, just never pointed at glioma), Golebiewska/LIH, Hutter/Basel (PMID 35776791). |
| **G2** — 8 weeks after the Agar email | Any reply, even to the PDX/TMT-map ask? | Stop pursuing. Record it. Reconsider only with an institution behind us. |
| **G3** — after §7 items 1–3 | Did the open-data dose→response test produce a signal? | If no: do not file the Burdenko DAR. The affiliation cost is not justified by a null. |
| **G4** — any point | Has an institutional affiliation or a collaborating PI been secured? | Until then, treat every dbGaP/MTA route in this document as closed, and say so rather than filing an application that will be rejected at the Signing Official step. |

---

## 10. Referenced source documents

- `results/DATA_HUNT.md` — the 164-dataset sweep, Tier 3 gated table, confirmed dead ends
- `results/SESSION_AUDIT.md` — the 11 defects, what survives, the honest position
- `results/ranking_robustness.md` §8 — the decisive-measurement analysis
- `results/gbm_tmz_exvivo.md` — the GSE232173 reanalysis and the MGMT forfeit
- `results/mgmt_scalar.md` — the MGMT→sensitivity mapping and its CI
- `results/glioma_regional_pk.md` — the human glioma delivery validation and the additive-leak fix
- `docs/real_data/drug_response.md` — the earlier (now partly superseded) HTAN/PDMR plan
