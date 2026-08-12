# Real spatial-data status — how far open data closes the loop

*Live run: 2026-07-26. Author artifact: `tmesim/data/spatial_ingest.py`.
Sections C and D rewritten 2026-08-12 after the 164-dataset sweep
(`results/DATA_HUNT.md`) falsified three of the routes they originally named.*

> ## The response axis cannot be validated from open data. This is a property of the data that exists, not a task on a list.
>
> A 39-agent sweep over 164 datasets (143 of them open), plus a re-check of the landscape on
> **2026-08-12**, establishes that **no open dataset anywhere pairs a measured spatial drug
> distribution with a measured treatment response in the same tumour.** Not one, in any tumour
> type. For glioma + temozolomide specifically, no spatial drug map exists at all, open or
> otherwise — METASPACE returns zero.
>
> The single experiment that would close it already exists and is not downloadable: the
> **Agar/White adavosertib GBM platform** — drug concentration per pixel in phase-1 patient
> tissue with co-registered molecular response — which has **no accession** and whose two
> published data-availability statements (2022 and 2025) both say "available on request".
> Obtaining it is a months-long collaboration request on human-subject tissue from a closed
> trial, requiring an IRB pathway and a data-use agreement, and therefore an institution.
>
> Drafts of that request, and of three others, are written, fact-checked and **UNSENT**:
> [`docs/outreach/`](outreach/). Nobody named in them has been contacted.
>
> Everything in sections A and B below is a **geometry** result. None of it is a response
> validation, and no amount of further engineering makes it one.

This document is the HONEST record of running the tmesim
reconstruct -> sim loop on a **real spatial-transcriptomics** GBM sample instead
of synthetic geometry. It states exactly what downloaded, what ran, and — the
crux — what is still blocked from a genuine drug-response validation.

TL;DR
-----
- **Reached a genuinely open, no-login, directly-downloadable spatial GBM
  dataset** (10x Visium Human Glioblastoma) and pulled the real bytes here.
- **Closed the loop on real-derived geometry:** measured spot coordinates +
  measured marker expression -> `reconstruct` (real vasculature/ECM/assemble) ->
  `sim.engine.simulate_response` with a compiled temozolomide. The engine accepts
  the real geometry and returns a differentiated prediction.
- **The loop is closed for GEOMETRY, not for VALIDATION.** Open spatial GBM data
  has **no paired drug-response outcome label**. We can predict a kill fraction
  on real tissue; we cannot check it against a measured response, because no
  dataset — open *or* controlled — pairs a measured spatial drug distribution
  with a measured response in the same tumour. **Controlled access does not fix
  this**; that was the original assumption in section C and it was wrong. The gap
  requires an experiment that has been run once and not shared (section C3).

---

## A. What actually downloaded and ran here

### Dataset reached (OPEN, direct download, verified HTTP 200)
**10x Genomics — "Human Glioblastoma, Whole Transcriptome Analysis"** (Visium v1,
Space Ranger 1.2.0). A single human GBM surgical section.
- Landing: <https://www.10xgenomics.com/datasets/human-glioblastoma-whole-transcriptome-analysis-1-standard-1-2-0>
- Direct files (no login, no DUA), under
  `https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/`:

| file | bytes (downloaded) | contents used |
|---|---|---|
| `..._spatial.tar.gz` | 9,134,896 | `tissue_positions_list.csv` (spot px coords), `scalefactors_json.json` (spot diameter -> um/px) |
| `..._filtered_feature_bc_matrix.tar.gz` | 54,134,945 | MEX sparse matrix: `matrix.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz` |

Both tarballs downloaded in full and are cached under
`data_cache/visium_gbm/`. Parsed with **stdlib + scipy/pandas only** — no
`scanpy`/`anndata`/`h5py` (not installed in the venv), which is why the MEX
tarball route (scipy.io.mmread) was chosen over the `.h5`.

### What the run produced (live, 2026-07-26)
- **3,468 in-tissue Visium spots** read; all 10 panel marker genes present
  (SOX2, OLIG2, EGFR, MKI67, GFAP, CD68, AIF1, CD3E, PECAM1, VEGFA).
- Tissue extent ~5.7 x 4.9 mm -> Domain `(152, 132)` at `dx=37.8 um`.
- Phenotype mixture (relative-enrichment argmax): **TUMOR 825, STROMAL 820,
  IMMUNE 871, ENDOTHELIAL 952**.
- `reconstruct` built the vascular field (`build_vascular_network` with a
  tumor-density mask) and ECM (`build_tissue_fields`), assembled via
  `assemble_microenvironment`.
- Compiled **temozolomide** from structure (`molecular.compile_drug`) ->
  `DrugParams(eff_D=278 um^2/s, ic50~606 uM, bbb=0.558)`.
- `simulate_response` ran on the real-derived geometry:
  **predicted kill fraction 0.354**, penetration index 0.298, mean oxygen 0.59.

Re-running `python -m tmesim.data.spatial_ingest` reproduces this end to end.

---

## B. What is REAL vs APPROXIMATED in the reconstruction (no overclaiming)

**Real (measured bytes drive these):**
- Spot **(x, y) coordinates** — from the Visium tissue-position file, converted to
  micrometers with the real scalefactor. This is the whole point: the geometry is
  no longer synthetic.
- Per-spot **marker-gene expression** — from the real count matrix.

**Approximated (documented in code + metadata):**
- **One agent per Visium spot.** A Visium spot is a ~55 um well pooling ~1-10
  cells; it is not a single cell. We model each spot as one representative agent
  at its centroid. This is a resolution approximation, not single-cell truth.
- **Phenotype = relative-enrichment marker argmax, NOT deconvolution.** Each spot
  is labelled by the marker module it is most enriched for (z-scored across spots
  before argmax). A raw argmax collapses ~all spots to "stroma" because **GFAP is
  broadly high in GBM** — reactive astrocytes *and* astrocytic tumor cells express
  it — so GFAP is not a clean tumor/stroma discriminator. Relative enrichment is a
  standard coarse dominant-module call; genuine per-spot composition needs
  reference-based deconvolution (cell2location / RCTD), which needs
  `scanpy`/`anndata` (not available here) and a matched single-cell reference.
- **Per-agent `sensitivity` is a heuristic**, tied to the tumor-marker signal
  (mirrors `reconstruct.segmentation._sensitivity_from_tumor`). It is **not** a
  measured drug sensitivity — no such label exists in this data (see C).
- **ECM + vasculature are model estimates**, not measured vessel/matrix stains.

So: the *spatial arrangement and expression are real*; the *phenotype call is
coarse and the sensitivity + fields are modelled*. This is a real
**geometry reconstruction**, and a real **mechanistic prediction** on it — it is
**not** a validated drug-response result.

---

## C. The crux: no paired drug-and-response measurement exists in the same tumour

*Rewritten 2026-08-12. The original version of this section named HTAN and NCI PDMR as the
routes that would close the loop. Direct checking falsified both. They are kept below, struck
through, because a plan that quietly deletes its wrong answers cannot be audited.*

### C1. Why the sample above cannot score the model

- The 10x Visium GBM sample (and the 10x **Xenium** human-brain/GBM preview, and the GEO
  Visium series **GSE194329**) are **untreated resections with no drug administered and no
  response measured.** They give tissue structure, full stop.
- **GSE174554** ("single-cell atlas of human glioma under therapy") pairs spatial assays with a
  *defined standard-of-care* (TMZ + RT) and a *recurrence* outcome — but that is a
  cohort-level survival/recurrence signal, **not** a per-sample measured drug kill, and the
  matched spatial + outcome levels are access-gated.
- So `tmesim` can predict a kill fraction **on** real tissue and cannot score it against a
  measured response **from** that tissue.

### C2. The general result, from a 164-dataset sweep

`results/DATA_HUNT.md` (2026-07-27; 39 agents, 164 datasets, 143 open, every URL HTTP-checked)
and a landscape re-check on **2026-08-12** establish:

1. **No open dataset anywhere pairs a measured spatial drug distribution with a measured
   per-specimen treatment response in the same tumour.** Every candidate breaks at that link.
2. **No glioma + temozolomide spatial drug map exists**, open or otherwise. METASPACE returns
   zero. Only three papers in all of PubMed cover drug-distribution MSI in glioma, and none
   deposits its MSI data.
3. **MGMT is absent from every standard Xenium panel**, so TMZ sensitivity was not merely
   unmodelled but unmodellable from that modality. One custom-panel deposit (GSE324268) closed
   the *detection* half of this — and `results/mgmt_spatial.md` then showed MGMT has no usable
   spatial structure at that depth, below the assay's own negative controls. The per-specimen
   scalar is the empirically correct object, not a compromise.
4. **Ion intensity is not concentration**, and **contrast-agent perfusion is not drug
   delivery.** Nothing in the open imaging world measures a therapeutic in tissue.
5. **Ex vivo systems delete the source term.** The one open dataset that comes closest
   (Moncayo et al., HGSOC explants, 376,736 MSI pixels at 20 µm with a heavy-labelled internal
   standard) is bath-dosed, has no vasculature, is ovarian rather than glioma, n = 3 patients,
   and its "response" is a transcriptional apoptosis signature. It validates diffusion, uptake
   and sequestration. It cannot validate perfusion-driven delivery.

**This is a structural limit of the data that exists.** It is not a backlog item, and no
further code, search or compute in this repository will remove it.

### C3. The single experiment that would close it

The **Agar/White adavosertib GBM platform**: measured adavosertib µM per pixel by MALDI-FT-ICR
MSI at a 30 µm raster, in GBM PDX *and* in specimens from patients on the ABTC phase-1 trial
NCT01849146, co-registered with phosphoproteomics and t-CyCIF DNA-damage-response markers.
Right disease, right modality, drug and response in the same tissue.

- Method paper: Lopez BGC et al., *Neuro-Oncology* 2022;24(1):64–77, PMID 34383057, PMC8730776.
- Trial report: Lee EQ et al., *Clin Cancer Res* 2025;31(6):983–992, PMID 39820282,
  PMC11913579 — measured 970 nM–7,270 nM across 10 specimens from 6 patients.
- **Neither deposits the MSI.** 2022: "The MALDI-MSI data are available upon request." 2025:
  "available on request from the corresponding author … not publicly available due to privacy
  restrictions." The only self-serve component, PRIDE PXD018782, was verified to contain no
  quantified phosphosite response — 80 files of raw spectra and peptide identifications, zero
  reporter-ion or ratio fields.

It is a **collaboration request on human-subject tissue from a closed trial**, which means an
IRB pathway and a data-use agreement, which means a counterparty institution. Timeline
6–18 months, with a real probability of no reply.

**A second glioma-specific asset has the same shape**: the Jonas Lab / Kibur Medical
implantable-microdevice study (Peruzzi & Jonas, BWH; PMID 37672566), which pairs per-reservoir
drug concentration with per-site caspase-3 / pH2AX in human glioma — and publishes it only as
plotted points, with the table available on request
(`results/zenodo_kibur_deposit_check.md`).

### C4. Status of the request: written, checked, UNSENT

Four outreach drafts exist in [`docs/outreach/`](outreach/), fact-checked against the primary
record on 2026-08-12. **None has been sent. No person or group named in them has been contacted
and none knows this project exists.** Two of the four (01 Agar/Santagata, 03 dbGaP) are further
blocked on institutional affiliation; the two that are not (02 Erasmus MC, 04 GSE324268
depositors) are limited only by a human deciding to send them and by writing the
pre-registration they promise, which does not yet exist.

### C5. Superseded — kept for audit

> ~~1. **HTAN (Human Tumor Atlas Network)** — GBM spatial multi-omics (CODEX/MIBI/Visium) with
> matched single-cell; access via dbGaP + HTAN DUA.~~
> **FALSIFIED: HTAN has no GBM atlas.** The only primary brain-tumour atlas is HTA202
> (low-grade glioma, seqFISH), with no brain data in any release through Nov 2025. Most HTAN
> imaging and Level 3/4 data is open CC BY 4.0 anyway. Monitor; do not apply.
>
> ~~2. **NCI PDMR** — patient-derived organoid / PDX drug screens for the measured response.~~
> **FALSIFIED: PDMR has zero glioblastoma models** (3 CNS cases across 1,318 patients) **and no
> public drug-response data at all** — "Preclinical Drug Study Results" is still *In
> Development*. Its omics are on anonymous FTP, not dbGaP. Dead end for this project.
>
> ~~3. The genuine validation is then: reconstruct HTAN GBM tissue, predict, and score against
> the paired organoid/PDX kill from the same model.~~
> Not available by that route. See C3 for the route that would work.

---

## D. Exact residual gaps (so nothing is hidden)

| gap | why it exists here | what removes it | status |
|---|---|---|---|
| **No paired spatial-drug + response measurement to score against** | **No such open dataset exists, in any tumour type** (164-dataset sweep; re-checked 2026-08-12) | The Agar/White phase-1 MSI, or an equivalent purpose-run experiment. C3. | **BLOCKED — structural. Requires a collaboration, an IRB pathway and an institution.** |
| No per-patient ex-vivo potency on a spatially profiled specimen | Never been run in glioma, as far as can be established | The Erasmus MC / GLIOscreen ask (`outreach/02`) — an add-on to a screen they already run | **BLOCKED on sending an email.** Draft ready, unsent. |
| Coarse phenotype (argmax, not deconvolution) on the Visium section | no `scanpy`/`anndata`; no matched sc reference | install spatial stack + cell2location/RCTD with a GBM sc atlas | Engineering. Not blocked. |
| One agent per spot (~55 µm), not per cell | Visium resolution | single-cell spatial — done elsewhere in the repo on three Xenium sections, incl. GSE324268 at 851,509 cells | Largely addressed. |
| ECM/vasculature modelled, not measured | Visium has no vessel/matrix channel | multiplexed imaging with real vessel/ECM masks | Open. |
| Compiled TMZ potency is low-confidence | structure alone cannot give IC50, and GDSC's TMZ arm has no dynamic range (all 51 glioma IC50s above the 30 µM ceiling) | a measured per-patient IC50 — the same ask as row 2 | **BLOCKED on the same email.** |

**Bottom line:** open data lets `tmesim` reconstruct a real GBM microenvironment geometry and
run the full physics pipeline on it, reproducibly, here. The delivery axis is separately
validated against human tissue measurements with a quantified error budget. The **response
axis is not validated and cannot be validated from open data** — that is a statement about the
world's datasets, not about this repository's engineering backlog, and it will stay true until
someone runs an experiment that does not currently exist or shares one that is not currently
shared.
