# Real single-cell-resolution spatial GBM — Xenium Human Brain (Glioblastoma)

*Acquired live 2026-07-26. Loader: `tmesim/data/real/singlecell_spatial.py`
(`load(path) -> {cell_xy_um, expression, genes, ...}`).*

This is the **single-cell-resolution** spatial layer the Visium ingest could
only approximate. Every row is one segmented cell with a measured `(x, y)`
centroid in micrometers and a targeted-panel count vector — true single-cell
positions, the ideal input for a cell-resolved microenvironment.

## Dataset

- **Name:** 10x Genomics **Xenium Human Brain Preview Data** — the
  **glioblastoma** section.
- **Assay:** Xenium In Situ (imaging-based single-cell spatial transcriptomics),
  Xenium Analyzer **v1.3.0**.
- **Sample id / accession:** `Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon`.
  One adult human GBM FFPE section (tissue from Avaden Biosciences). The same
  release also ships an Alzheimer's and a healthy brain section (not downloaded).
- **Panel:** Xenium Human Brain base panel **`hBrain_254g`** (254 genes, design
  `PD_326`) + a custom **65-gene add-on** (`PD_327`) = **319 gene-expression
  targets**. The matrix additionally carries 41 negative-control codewords, 20
  negative-control probes, and 161 unassigned codewords (541 features total);
  the loader drops these by default so `genes` is the real 319-gene panel.
- **Open status:** **OPEN.** Licensed **CC BY 4.0**, no login, no DUA. Individual
  files are directly downloadable from `cf.10xgenomics.com` (verified HTTP 200,
  content-length matched on download).
- **Landing page:**
  <https://www.10xgenomics.com/datasets/xenium-human-brain-preview-data-1-standard>

## What downloaded (into `data_cache/xenium_gbm/`)

Base URL:
`https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon`

| local file | url suffix | bytes | contents |
|---|---|---|---|
| `Xenium_GBM_cells.csv.gz` | `_cells.csv.gz` | 2,167,253 | per-cell centroids + QC metadata (40,887 rows) |
| `Xenium_GBM_cell_feature_matrix.tar.gz` | `_cell_feature_matrix.tar.gz` | 14,279,458 | MEX triple: `matrix.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz` |
| `Xenium_GBM_gene_panel.json` | `_gene_panel.json` | 152,134 | panel identity + 339 target definitions |

**Total ~16.6 MB — this is the SMALLEST COMPLETE real unit** (a full GBM section:
all cells + full panel counts + panel spec). The complete downloads are the real
bytes; nothing here is untested-on-full-bytes.

Not downloaded (large, not needed for cell-resolved geometry): the full
`_outs.zip` bundle (**6.82 GB**, contains transcripts.parquet, morphology OME-TIFF,
cell/nucleus boundary polygons, Xenium Explorer files) and the H&E image. The
richer sibling GBM dataset (below) was also left as documented-only.

### Sibling GBM dataset (documented, verified URLs, NOT downloaded)
10x **Xenium_V1_Human_Brain_GBM_FFPE** (Xenium v2.0.0), GBM FFPE with the
**Human Immuno-Oncology Profiling Panel** (380 genes) — richer for immune niche
work. Individual files verified reachable (HTTP 200):
`cells.csv.gz` 44,907,408 B; `cell_feature_matrix.tar.gz` 242,459,483 B;
`gene_panel.json` 256,404 B; `xe_outs.zip` 13.15 GB. Same on-disk format as
below, so the loader ingests it unchanged by pointing `path` at its cache.
Landing:
<https://www.10xgenomics.com/datasets/ffpe-human-brain-cancer-data-with-human-immuno-oncology-profiling-panel-and-custom-add-on-1-standard>

## Format (what ingestion targets)

**`cells.csv.gz`** — one row per cell, header:
```
cell_id, x_centroid, y_centroid, transcript_counts, control_probe_counts,
control_codeword_counts, unassigned_codeword_counts, total_counts,
cell_area, nucleus_area
```
- `x_centroid`, `y_centroid` are in **micrometers** (Xenium global coordinate
  frame). Measured extent here: **8,087.8 x 5,065.4 um** (~8.1 x 5.1 mm),
  x in [2.3, 8090.2], y in [5.1, 5070.4].
- `cell_area`, `nucleus_area` are in **um^2**.

**`cell_feature_matrix/`** — standard 10x MEX sparse matrix (`scipy.io.mmread`):
- `matrix.mtx.gz`: **features x cells** = **541 x 40,887** raw integer counts.
- `features.tsv.gz`: 3 columns `ensembl_id, symbol, feature_type`. Feature-type
  counts: 319 `Gene Expression`, 20 `Negative Control Probe`, 41 `Negative
  Control Codeword`, 161 `Unassigned Codeword`.
- `barcodes.tsv.gz`: 40,887 `cell_id`s. **Column i of the matrix == row i of
  `cells.csv`** (verified: identical `cell_id` order). The loader still reindexes
  defensively by `cell_id`.

**`gene_panel.json`** — `payload.panel` gives identity/description/species/tissue;
`payload.targets` lists 339 target definitions.

## Parse result (venv, verified)

- **n_cells = 40,887**; **319 panel genes**; expression returned as a
  `scipy` CSR matrix, cells x genes = `(40887, 319)` (or dense float32 on request).
- **Coordinate units: micrometer** (native; no conversion needed — contrast with
  Visium, which needs `scalefactors_json` to go px -> um).
- GBM-relevant panel genes confirmed present: EGFR, SOX2, OLIG2, MKI67, CD68,
  AIF1 (IBA1), PTPRC (CD45), PECAM1 (CD31), CD163, MOG. (Absent from this brain
  panel: GFAP, VEGFA, CD3E — use SOX2/OLIG2, PTPRC/CD68, PECAM1 as substitutes.)

Reproduce: `python -m tmesim.data.real.singlecell_spatial`.

## Honest limitations

- **No paired drug-response outcome.** Like the open Visium GBM section, this is
  a single untreated resection with no per-sample RECIST / drug readout. It
  supports true single-cell geometry + a mechanistic prediction, **not**
  supervised response validation — the same access gap tracked in
  `docs/real_data_status.md`.
- **Targeted panel, not whole transcriptome** (319 genes). Rich for cell typing
  and niche markers; not a substitute for the Visium WTA when a gene is off-panel.
- One section, one patient. No replicates in this download.
