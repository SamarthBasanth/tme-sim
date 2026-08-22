# Real spatial-data status — how far open data closes the loop

*Last run live: 2026-07-26. Author artifact: `tmesim/data/spatial_ingest.py`.*

> ## ⚠ ARCHIVED RECORD — DO NOT QUOTE ANY NUMBER FROM THIS FILE
>
> This is a historical snapshot of a **superseded** pipeline, kept for the record.
> Every quantitative claim below predates corrections that have since landed, and
> at least four are now known-wrong:
>
> * `bbb=0.558` for temozolomide is the **retracted arithmetic-mean** value; the
>   shipped geometric-mean veto (`molecular/permeability.py`) returns **0.1254**.
> * `penetration index 0.298` quotes a **RETIRED** metric — `penetration_index` is
>   a scale-invariant shape statistic that reports the section's vascular geometry,
>   not the drug's delivery (`results/delivery_metrics.md` §4). Use
>   `therapeutic_coverage` at 1× IC50.
> * `predicted kill fraction 0.354` is an **absolute kill**, which is not a
>   defensible output under any configuration, and is quoted here with no interval
>   (the propagated interval spans ~30×).
> * The reconstruction path described (`build_vascular_network`,
>   `build_tissue_fields`, `assemble_microenvironment`,
>   `segmentation._sensitivity_from_tumor`) is classified **legacy, not on the real
>   path** by `docs/architecture.md`, and the `ENDOTHELIAL 952` call comes from the
>   relative-enrichment argmax that was the project's defect #1.
>
> For the current state read [`../README.md`](../README.md),
> [`../results/CLAIM_INVENTORY.md`](../results/CLAIM_INVENTORY.md) and
> `vivum validate`.

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
  on real tissue; we cannot yet check it against a measured response, because no
  open dataset pairs spatial tissue with a per-sample drug readout. That gap
  requires controlled-access data (section C).

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

## C. The crux: NO paired outcome label in open spatial data

This is the honest blocker, and it is a data problem, not a code problem.

- The 10x Visium GBM sample (and the 10x **Xenium** human-brain/GBM preview, and
  the GEO Visium series **GSE194329**) are **untreated resections with no drug
  administered and no response measured.** They give tissue structure, full stop.
- **GSE174554** ("single-cell atlas of human glioma under therapy") pairs spatial
  assays with a *defined standard-of-care* (TMZ + RT) and a *recurrence* outcome —
  but that is a cohort-level survival/recurrence signal, **not** a per-sample
  measured drug-kill, and the matched spatial + outcome levels are access-gated.
- Consequently, tmesim can today predict a kill fraction **on** real tissue, but
  **cannot score it** against a measured response from the same tissue. Every open
  spatial GBM dataset we could reach lacks the label.

### What actually closes the validation loop (controlled access)
To turn "predicts on real geometry" into "predicts the *right* answer" you need
spatial tissue **paired with a measured drug response**, which lives behind DUAs:

1. **HTAN (Human Tumor Atlas Network)** — GBM spatial multi-omics (CODEX/MIBI/
   Visium) with matched single-cell; the correct substrate for real
   `reconstruct/` segmentation + deconvolution. Access via dbGaP + HTAN DUA;
   imaging levels need an approved data-access request.
   Portal: <https://data.humantumoratlas.org>.
2. **Patient-derived organoid / PDX drug screens for the MEASURED response** —
   NCI **PDMR** (<https://pdmr.cancer.gov>), **HCMI** GBM models via GDC, and GBM
   organoid biobanks (e.g. Jacob et al., *Cell* 2020). These provide the
   drug-response ground truth the spatial atlases lack.
3. The genuine validation is then: reconstruct HTAN GBM tissue with
   deconvolution-grade phenotyping, predict response with tmesim, and score
   against the paired organoid/PDX kill from the same model. Only that makes the
   Stage-0 "does mechanistic TME beat a non-spatial baseline?" a real claim.

---

## D. Exact residual gaps (so nothing is hidden)

| gap | why it exists here | what removes it |
|---|---|---|
| No response label to score against | open spatial GBM data has none | HTAN tissue + PDMR/HCMI/organoid drug screen (DUAs) |
| Coarse phenotype (argmax, not deconvolution) | no `scanpy`/`anndata`; no matched sc reference | install spatial stack + run cell2location/RCTD with a GBM sc atlas |
| One agent per spot (~55 um), not per cell | Visium resolution | use **Xenium**/CODEX single-cell spatial (10x Xenium GBM is open but panel-limited & unlabelled) |
| ECM/vasculature modelled, not measured | Visium has no vessel/matrix channel | multiplexed imaging (CODEX/MIBI, in HTAN) for real vessel/ECM masks |
| Compiled TMZ potency is low-confidence | structure alone can't give IC50 | anchor IC50 to a measured GDSC/organoid assay (already wired in `molecular.potency`) |

**Bottom line:** open data lets tmesim reconstruct a *real* GBM microenvironment
geometry and run the full physics/sim pipeline on it — that half of the loop is
genuinely closed and reproducible here. The *validation* half stays open until a
DUA unlocks spatial tissue paired with a measured drug response; that is a data
access step, not remaining engineering.
