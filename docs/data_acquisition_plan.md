# Data Acquisition Plan — real glioblastoma data for `tmesim`

> # ⚠ SUPERSEDED — DO NOT FOLLOW THIS DOCUMENT
>
> **Checked 2026-09-03. This page describes a pipeline that no longer exists, and it
> was the only stale file in `docs/` carrying no banner saying so.**
>
> * **`tmesim/data/loaders.py` does not exist**, and neither does the
>   `load_real_cohort()` it documents in §B. `grep -rn "load_real_cohort"` over the
>   tree returns no definition — only this file and `ARCHITECTURE_CONTRACT.md`
>   referring to one.
> * Their **absence is now a pass condition**, not a gap.
>   `scripts/reproduce_real.sh:176-180` fails the harness outright if
>   `tmesim/data/loaders.py` or `tmesim/data/synthetic.py` reappears on the import
>   path, and prints `synthetic modules absent … OK` when they do not.
> * The **synthetic 2-D layout** in §B is part of what that check exists to keep out.
> * The controlled-access section (§C) is superseded by
>   [`DATA_ACCESS_PLAN.md`](DATA_ACCESS_PLAN.md) §6–§7, which re-checked these routes
>   and found several of them do not hold the data assumed here.
>
> **Read [`REAL_BUILD_PLAN.md`](REAL_BUILD_PLAN.md) and
> [`DATA_ACCESS_PLAN.md`](DATA_ACCESS_PLAN.md) instead.** This file is kept as the
> record of what was planned in July 2026 and what replaced it — the §A source table
> is still a useful account of which GDC/GDSC endpoints were reachable on that date.

*Last verified live: 2026-07-25 (GDC Data Release 45.0; GDSC release 8.5).*

This document records (A) what the open-tier live acquisition in
`tmesim/data/acquire.py` actually pulls, (B) how `tmesim/data/loaders.py` turns
it into `PatientRecord`s, its approximations, and (C) the **controlled-access**
datasets a funded team would use for genuine spatial validation.

---

## A. What downloads today (open tier, live)

`acquire_all()` writes to `./data_cache/` and returns a per-source status dict.
On the last verified run, **all four sources succeeded**:

| source | status | what it is | records | cache path |
|---|---|---|---|---|
| `tcga_gbm_clinical` | ok | TCGA-GBM clinical via GDC `cases` API (age, sex, vital status, survival, per-patient therapeutic agents) | 80 patients | `data_cache/tcga_gbm_clinical.json` |
| `tcga_gbm_expression` | ok | OPEN STAR gene-count TSVs via GDC `files` + `data`, filtered to the same clinical patients | 8 files | `data_cache/tcga_expression/` |
| `ivy_gap` | ok | Ivy GAP (Allen) RNA-seq gene-expression matrix zip (anatomic-structure atlas) | 270 samples | `data_cache/ivy_gap/gene_expression_matrix.zip` |
| `gdsc_glioma` | ok | GDSC2 fitted dose-response (Sanger), streamed + filtered to glioma cell lines | 12,001 rows / 51 cell lines | `data_cache/gdsc_glioma_dose_response.csv` |

### Endpoints used (all verified reachable)
- GDC cases: `https://api.gdc.cancer.gov/cases` (filter `project.project_id = TCGA-GBM`)
- GDC files: `https://api.gdc.cancer.gov/files` (open `Gene Expression Quantification`, `STAR - Counts`)
- GDC data: `https://api.gdc.cancer.gov/data/{file_id}`
- Ivy GAP zip: `https://glioblastoma.alleninstitute.org/api/v2/well_known_file_download/305873915`
- GDSC2 CSV: `https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.csv`

### Sources considered and rejected
- **`ftp.sanger.ac.uk` / `www.cancerrxgene.org/api`** — not reachable / returns HTTP 410.
  The `cog.sanger.ac.uk` object store (used above) is the working mirror.
- **DepMap/PRISM via figshare** — reachable but the direct `ndownloader` links are
  short-lived signed-S3 redirects and the glioma subset needs a separate cell-line
  annotation join; GDSC2 already carries `TCGA_DESC=GBM/LGG` inline, so it is the
  cleaner open drug-response ground truth and DepMap was dropped to avoid fragility.

Every network call has a `requests` timeout, a single retry, and is wrapped in
`try/except`; a failed source is recorded as `failed`/`skipped` and never crashes
the pipeline.

---

## B. From cache to `PatientRecord`s (and the approximations)

`load_real_cohort()` produced **131 real records** on the verified run:
**51 GDSC glioma cell lines + 80 TCGA-GBM patients** (63/131 labelled responders).

### Response labels (real, but proxied)
- **GDSC cell lines** — `observed_response` from the real Temozolomide LN_IC50
  drug screen (median split; lower IC50 = sensitive). This is a genuine
  drug-response ground truth, just in monoculture, not tissue.
- **TCGA-GBM patients** — the open tier has **no per-drug RECIST response**, so
  `observed_response` is an **overall-survival median-split proxy** (longer
  survival ≈ better therapy response — imperfect, documented in-code).

### The microenvironment is a COARSE APPROXIMATION
Real open GBM data has **no spatial imaging attached**. So each record's
`Microenvironment` is a **synthetic 2-D layout parameterised by the real
features**, not a reconstruction of that patient's tissue:
- tumour-core size ← proliferation signal (MKI67);
- vascular density ← angiogenesis signal (VEGFA), sparser in the hypoxic core;
- immune fraction ← macrophage/T-cell markers (CD68, CD3E);
- per-cell `sensitivity` ← anchored to the real response (GDSC IC50; TCGA
  MGMT-methylation proxy = low MGMT expression).

Marker genes (MGMT, EGFR, PTEN, NF1, VEGFA, CD68, MKI67, OLIG2, CHI3L1, CD3E)
are parsed from the STAR-count TSVs (`fpkm_unstranded`) for the 7–8 patients who
have both clinical and expression cached. Everything is deterministic (seed from
the patient/cell-line id) and never raises.

**Bottom line:** the *labels and molecular features are real*; the *spatial
arrangement is synthetic*. This suffices to exercise the Stage-0 pipeline on
genuine cohorts but is **not** a spatial validation.

---

## C. Controlled-access datasets for REAL spatial validation

A funded team with data-use agreements should validate against tissue with true
spatial coordinates + measured drug response. Priority targets:

### 1. HTAN — Human Tumor Atlas Network (spatial multi-omics atlases)
- Portal: <https://data.humantumoratlas.org> (CDS/Synapse-gated Level-3/4 access).
- Contains GBM atlases with CODEX/MIBI/Visium and matched single-cell — the
  right substrate for real `reconstruct/` segmentation + phenotyping.
- Access: dbGaP + HTAN DUA; imaging levels require an approved data-access request.

### 2. GEO spatial-transcriptomics GBM series (accessions verified live)
- **GSE194329** — *"Spatial multi-omics … radial glia-like cancer-initiating
  cells … malignant gliomas"*: 10x **Visium** on GBM/DMG surgical resections
  (5 GBM incl. 2 secondary + peritumor). Spatially resolved expression.
- **GSE174554** — *"A single-cell atlas of human glioma under therapy"*: 80
  IDH-wt GBM (40 primary + 40 patient-matched recurrent), all standard-of-care
  (temozolomide + RT), with single-cell **digital-spatial** and
  **spatial-transcriptomic** assays — pairs spatial data with a defined therapy
  and recurrence outcome.
- **10x Genomics public GBM datasets** — Visium and **Xenium** human glioblastoma
  demonstration panels: <https://www.10xgenomics.com/datasets> (filter: brain /
  glioblastoma). Open but no drug-response labels; useful for spatial-model
  calibration of `physics/` and `reconstruct/`.

### 3. Patient-derived organoid / explant drug screens (response ground truth)
- **GBM organoid biobanks** (e.g. Jacob et al., *Cell* 2020) — patient-derived
  GBM organoids preserving microenvironment, with drug-response readouts;
  associated data via the authors' supplements / EGA.
- **HCMI** (Human Cancer Models Initiative) GBM models via GDC
  (<https://portal.gdc.cancer.gov> → HCMI-CMDC) — models with matched molecular
  data; some sensitivity data controlled.
- **PDMR** (NCI Patient-Derived Models Repository, <https://pdmr.cancer.gov>) —
  GBM PDX/organoid models with drug-treatment data.

### 4. Longitudinal / treatment-response clinical cohorts
- **GLASS consortium** (Glioma Longitudinal Analysis) — Synapse `syn17038081`
  / <https://www.synapse.org> — primary→recurrent matched GBM with treatment.
- **CGGA** (Chinese Glioma Genome Atlas, <http://www.cgga.org.cn>) — glioma
  expression + **clinical treatment/response and survival** fields; largely open,
  a stronger real response label than TCGA's survival proxy.

### What "real validation" would replace
Swap the synthetic-spatial `_build_coarse_micro` step for tissue reconstructed
from HTAN/Visium/Xenium images via `reconstruct/segmentation.py`, and swap the
survival proxy for measured organoid/PDX drug response (PDMR/HCMI/organoid
biobanks). Only then does the Stage-0 "does the mechanistic TME beat a non-spatial
baseline?" test become a genuine claim rather than a plumbing check.
