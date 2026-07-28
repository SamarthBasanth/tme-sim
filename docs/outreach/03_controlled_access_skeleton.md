# DRAFT — NOT SENT. FOR HUMAN REVIEW AND EDITING.

> **Status: unsent skeleton.** No application has been filed with dbGaP, HTAN, NCI PDMR or
> TCIA. Fill every `[BRACKET]`.
>
> ## ⚠️ THIS CANNOT BE FILED AS THINGS STAND
>
> Every route below requires an **institutional Signing / Authorized Official**, an
> **eRA Commons account that an institution must create**, and a **PI at tenure-track
> professor or senior-scientist level**. dbGaP explicitly bars graduate students and
> postdoctoral fellows from submitting. An unaffiliated individual cannot pass these gates.
>
> NIH does publish an "unaffiliated/independent applicant" registration route, but it is
> written for **grant submission**, says nothing about data access requests, and would have
> the applicant act as their own Signing Official and own IT Director — which is precisely the
> conflict the requirement exists to prevent. **Treat it as unavailable.**
>
> **Do not file any of these until an affiliation or a collaborating PI exists** (decision
> gate G4 in `../DATA_ACCESS_PLAN.md` §9). Filing an application that will be rejected at the
> Signing Official step wastes a reviewer's time and burns the sender's name.
>
> **Priority 4.** This file exists so the text is ready the day the gate opens, and so nobody
> re-researches the mechanics.

---

## 1. Read this before deciding to apply at all

`../DATA_ACCESS_PLAN.md` §5 and §6 record what direct checking of these four repositories
found. Summarised, because it changes whether an application is worth filing:

| Route | What we assumed | What is actually there |
|---|---|---|
| **HTAN** | GBM spatial multi-omics behind dbGaP + HTAN DUA | **No GBM atlas exists.** The only primary brain-tumour atlas is HTA202 (Caltech/UCLA), **low-grade glioma by seqFISH**, and no brain data appears in any release through Nov 2025. HTAN's imaging and Level 3/4 data are largely open CC BY 4.0 anyway. **Nothing to request today — monitor instead.** |
| **NCI PDMR** | Organoid / PDX drug screens via MTA | **Zero glioblastoma models** (3 CNS cases across 1,318 patients) and **no public drug-response data at all** — "Preclinical Drug Study Results" is still *In Development*. Omics are on **anonymous FTP**, not dbGaP. **Do not apply.** |
| **TCIA restricted** | A TCIA restricted-license form, weeks | **Retired 2025-07-08.** Now a dbGaP DAR against **`phs004225`**, which imports the institutional-SO requirement. |
| **Burdenko-GBM-Progression** | Gated | **Half of it is already open.** The 651 per-follow-up response labels, IDH, MGMT and the treatment timeline are CC BY 4.0. Only the pixels are controlled. |

**So of the four, only one is worth an application, and only conditionally**: the TCIA/dbGaP
route for Burdenko's RTDOSE volumes and Brain-Tumor-Progression's perfusion maps — and
`../DATA_ACCESS_PLAN.md` §5.5 argues for running the same experiment on the **open** CFB-GBM
(n = 264, includes RTDOSE) and MU-Glioma-Post (n = 203, progression labels) collections first.
A null there means the application was never worth filing; a signal there gives the Research
Use Statement below a concrete, defensible motivation, which is also what makes a DAR likely to
be approved.

---

## 2. Mechanics, per route

### 2.1 dbGaP DAR — covers controlled TCIA imaging and HTAN Level 1/2

- **Accessions.** TCIA face-containing imaging: **`phs004225`** ("Radiology Imaging Data
  Containing Faces from The Cancer Imaging Archive", v1.p1, released 2026-01-16). **One request
  covers all affected TCIA collections** — Burdenko, Brain-Tumor-Progression, QIN-GBM,
  RHUH-GBM, TCGA-GBM/LGG, CPTAC-GBM, IvyGAP, ACRIN-FMISO-Brain and the rest of the 50 now
  marked "Limited". HTAN: **`phs002371`** (single umbrella study for the whole network,
  currently v9.p2).
- **Prerequisites**: eRA Commons account (institution-created); PI at professor/senior-scientist
  level; institutional **Signing Official** who co-signs before NIH sees the request; a named
  **IT Director** attesting to NIH Security Best Practices; a study-specific **Data Use
  Certification**; and, for datasets so flagged, IRB approval (*whether `phs002371` is
  IRB-flagged is not verified — check on the study page*).
- **Written content**: a **Research Use Statement (≤2,000 characters)** and a **non-technical
  summary (≤1,100 characters)**. Drafts in §3.
- **Turnaround**: DAC review averages ~2 weeks. End to end, realistically 2–6 weeks with a
  cooperative SO, 3–6 months from a cold start. Approval term 1 year, renewable; closeout
  requires SO-certified data destruction.
- **Hard exclusion**: applicants at institutions in China (incl. Hong Kong and Macau), Russia,
  Iran, North Korea, Cuba and Venezuela are prohibited outright.
- **After approval** (TCIA specifically): log in at <https://nci-crdc.datacommons.io/login> →
  Profile → Create API Key (JSON) → download the collection manifest from the collection page →
  point the **TCIA Data Retriever** at the JSON key. `help@cancerimagingarchive.net` is still
  valid but no longer grants access; use it only for Data Retriever and manifest problems.
- References: <https://www.cancerimagingarchive.net/nih-controlled-data-access-policy/> ·
  <https://grants.nih.gov/policy-and-compliance/policy-topics/sharing-policies/accessing-data/dbgap> ·
  <https://docs.humantumoratlas.org/data_access/db_gap/>

### 2.2 NCI PDMR MTA — documented for completeness; do not file

Portal: <https://dctd.cancer.gov/drug-discovery-development/reagents-materials/pdmr>.
Four documents, all under `.../using-pdmr/request-material/` (note: the MTA moved from
`/sops/` — older links 404): the **Request Form – Domestic**, the **Material Transfer
Agreement**, the **Material Transfer Terms of Agreement** (v. 2026-04-01, incorporated by
reference), and the **Material Request Procedures**. **The PI signs the Request Form; the
institution's Authorized Official signs the MTA.** "Modifications to NCI's MTA Terms of
Agreement will not be considered." Pricing $0 NCI / $250 academic / $2,500 commercial per
cryopreserved PDX or organoid vial, recipient-paid shipping (FedEx account number mandatory),
no human-subjects use, no commercial use, no onward transfer, and an **annual report for 5
years** covering STR profiles, NGS, drug sensitivity and metabolomics. Turnaround not
published. Contact `NCI_PDM_Repository@mail.nih.gov`.

**Again: there are no glioblastoma models. This route is a dead end for this project.**

### 2.3 HTAN — monitor, do not apply

Set a calendar reminder to re-check **HTA202** (<https://humantumoratlas.org/center/hta202>)
and the release log (<https://humantumoratlas.org/data-updates>) in six months. If low-grade
glioma seqFISH data lands, revisit. One thing to check per-file when it does: the HTAN
data-sharing paper describes controlled access as covering "level 1 and 2 genomic and
transcriptomic data, as well as level 2 imaging data," while the portal describes imaging as
open CC BY 4.0. Confirm on the file's portal record rather than assuming either.

---

## 3. Draft application text

### 3.1 Research Use Statement — dbGaP `phs004225` (TCIA imaging)

*Target ≤2,000 characters. Current draft is close to the limit; trim the last paragraph first.*

> We are validating an open-source mechanistic model that predicts the spatial distribution of
> therapeutic dose within glioblastoma tissue and the resulting local response. The model is
> reconstructed from spatial transcriptomics and solves reaction–diffusion transport through
> patient-specific tissue architecture. Its delivery predictions have been scored against 125
> human glioma tissue samples across 42 patients and 27 published regional pharmacokinetic
> measurements, achieving approximately 2-fold geometric-mean error on tumour:plasma partition.
> Its **response** predictions have never been tested, because no open dataset pairs a measured
> spatial dose distribution with a measured local outcome in the same tumour.
>
> Burdenko-GBM-Progression is, to our knowledge, the only accessible collection that supplies
> both. Its DICOM RTDOSE objects give a true three-dimensional delivered-dose field, and its
> openly licensed clinical table supplies up to ten per-follow-up response assessments
> (`response` / `stable` / `progression` / `pseudoprogression`) across 180 patients, with
> IDH1/2 and MGMT status. We will test one pre-registered hypothesis: that voxels receiving low
> delivered dose are enriched for subsequent local progression, and that our model's spatial
> coupling function predicts the observed dose–local-control relationship better than a
> dose-independent baseline. Brain-Tumor-Progression adds pre-computed perfusion maps (nCBF,
> nrCBV, srCBV, ADC) in the same patients post-chemoradiation and at progression, letting us
> test whether measured perfusion heterogeneity improves the transport model's source term.
>
> We will state explicitly and prominently that radiation dose is not drug concentration: the
> mechanism of cell kill, the time course, and the delivery physics all differ. This request
> validates the **shape** of the model's spatial dose→response coupling, not its pharmacology.
> Analysis is retrospective and computational. No attempt will be made to re-identify any
> individual; face-containing volumes will be used only for registration and never rendered,
> published, or redistributed. Findings will be published openly, including if the model fails,
> and all analysis code will be released.

### 3.2 Non-technical summary

*Target ≤1,100 characters.*

> Brain tumours are hard to treat partly because treatment does not reach every part of the
> tumour equally. We have built a computer model that predicts where treatment reaches inside
> an individual patient's tumour, using detailed measurements of the tumour's own tissue
> structure. We have already checked that the model gets the amount of drug reaching tumour
> tissue roughly right, by comparing it with published measurements from 42 patients. What we
> have not been able to check is whether the places the model says receive less treatment are
> the places where the tumour actually comes back. Answering that requires medical images that
> record both the dose delivered and what happened afterwards. We are requesting two such
> imaging collections. The work is entirely retrospective and computational; we will not
> contact any patient and will make no attempt to identify anyone. Results will be published
> openly whether or not the model turns out to be correct, and the software will be free for
> anyone to use and check.

### 3.3 The specific data elements needed

Include as an appendix or in correspondence. Being specific shortens review and shows the
request is not a fishing expedition.

**From Burdenko-GBM-Progression** (DOI 10.7937/E1QP-D183):

| Element | Why |
|---|---|
| **DICOM RTDOSE** volumes | The delivered 3D dose field — the independent variable. Nothing else in the request works without it. |
| **RTSTRUCT** — GTV/CTV/PTV, and the subset with **GTV contoured on follow-up studies** | Defines where the tumour was and where it recurred. The follow-up GTV subset is the actual outcome geometry. |
| Planning CT | The frame RTDOSE and RTSTRUCT are defined in. |
| Follow-up MR series (T1ce, FLAIR) | Localising progression; also our reconstruction's input modality. |
| *Already open, no request needed* | The clinical CSV — `Response_1st_fup` … `Response_10th_fup`, IDH1/2 (110 pts), MGMT (92 pts), surgery/RT/follow-up/death dates. Download it now and cite it in the application. |

**Known limitation to declare, not discover later:** MRI→CT registration files are **not**
provided in this collection. We will register them ourselves and report the registration error.

**From Brain-Tumor-Progression** (DOI 10.7937/K9/TCIA.2018.15quzvnb, 20 patients, 3.16 GB):

| Element | Why |
|---|---|
| Pre-computed DICOM **nCBF, nrCBV, srCBV, ADC** maps | The only TCIA GBM collection shipping perfusion maps ready-made. Tests whether measured perfusion improves our vascular source term. |
| T1w pre/post, FLAIR, T2w | Co-registered to T1+C already. |
| Binary tumour masks | Provided; saves a segmentation step. |
| Both timepoints (≤90 days post-CRT, and at progression) | The within-patient contrast is the whole point. |

**Limitation to declare:** 40 studies across 20 patients is exactly two per patient, so **every
patient is a progressor by construction** — there is no control arm and no tabular outcome
file. This collection can support a within-patient spatial comparison, not a
progressor-vs-non-progressor classification. Say so in the application.

**Explicitly NOT requested: QIN-GBM-Treatment-Response.** Despite its name it contains no
response data — it is a double-baseline test–retest repeatability study in which both scans are
acquired 2–5 days apart and **before** therapy begins. Stating this in the application is
worthwhile: it shows the collections were read rather than keyword-matched.

---

## 4. Pre-filing checklist

- [ ] An institutional affiliation or a collaborating PI exists (gate **G4**). **If not, stop.**
- [ ] The open-data version of the experiment has been run on **CFB-GBM** (DOI 10.7937/v9pn-2f72,
      n = 264, includes RTDOSE) and **MU-Glioma-Post** (DOI 10.7937/7K9K-3C83, n = 203,
      progression labels) — gate **G3**. A null means do not file.
- [ ] Burdenko's **open** clinical CSV has been downloaded and analysed, and anything answerable
      without the pixels has been answered.
- [ ] The PI meets the professor / senior-scientist bar and holds an eRA Commons account.
- [ ] The Signing Official has agreed **in advance** to co-sign, and an IT Director is named.
- [ ] The institution is not in an excluded jurisdiction.
- [ ] The hypothesis is **pre-registered** with a timestamp and git hash before the data arrives.
- [ ] The "radiation dose is not drug concentration" caveat appears in the RUS, in the
      pre-registration, and in the caption of every resulting figure.
