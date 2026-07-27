# Real Visium (spot-level spatial transcriptomics) glioblastoma sections

All data here is **REAL** Space Ranger output (measured bytes on disk): a MEX
sparse expression matrix (`matrix.mtx.gz` + `features.tsv.gz` +
`barcodes.tsv.gz`) plus a `spatial/` folder with per-spot pixel coordinates
(`tissue_positions*.csv`) and `scalefactors_json.json`. No synthetic stand-ins.

Loader: `tmesim/data/real/visium.py`
  - `list_samples()` -> ids present locally.
  - `load_sample(sample_id)` -> `{coords_um, expression, genes, barcodes, ...}`
    (`expression` is a scipy CSR of shape `(n_spots, n_genes)`, integer UMI
    counts; `coords_um` is `(n_spots, 2)` micrometers aligned to `barcodes`).

Cache root: `data_cache/visium_gbm/`.

## Honesty note on what is / isn't here
- These are real GBM **tissue architectures + expression**. They carry **no
  paired drug-response / outcome label** — open spatial GBM data never does.
  They support geometry reconstruction and mechanistic prediction, not
  supervised drug-response validation.
- One Visium spot is a ~55 um well pooling ~1-10 cells (Visium v1), NOT a single
  cell. The loader returns raw per-spot counts; it does not deconvolve.

## Resolution (Visium v1, all samples)
- Spot diameter: 55 um. Center-to-center spot pitch: **100 um** (fixed hex grid).
- `um-per-pixel (umpp) = 55 / spot_diameter_fullres` (from `scalefactors_json.json`).
  This is the image scale used to convert spot pixel coords to micrometers; it
  differs per sample because each section was imaged at a different magnification.

## Samples

| sample_id        | accession                              | source / open status | n_spots (in-tissue) | n_genes | umpp (um/px) | tissue extent (um) | median UMI/spot |
|------------------|----------------------------------------|-----------------------|--------------------:|--------:|-------------:|--------------------|----------------:|
| `10x_parent_gbm` | 10x `Parent_Visium_Human_Glioblastoma` | 10x Genomics public, direct download, no login (OPEN) | 3468 | 36601 | 0.614 | 4908 x 5668 | 11596 |
| `GBM1`           | GEO `GSM5833533` (series `GSE194329`)  | NCBI GEO, direct download (OPEN) | 4764 | 36601 | 0.256 | 5708 x 5408 | 6614 |
| `GBM2`           | GEO `GSM5833534` (series `GSE194329`)  | NCBI GEO, direct download (OPEN) | 2891 | 36601 | 0.258 | 5447 x 5376 | 7269 |
| `GBM3`           | GEO `GSM5833535` (series `GSE194329`)  | NCBI GEO, direct download (OPEN) | 1775 | 36601 | 1.011 | 4491 x 3637 | 21580 |

All four samples contain the full GBM marker panel used by the pipeline
(tumor: SOX2/OLIG2/EGFR/MKI67/GFAP; immune: CD68/AIF1/CD3E/PTPRC; endothelial:
PECAM1/VEGFA; oligodendrocyte: MBP/PLP1) and those markers carry real,
biologically-sensible signal (e.g. GFAP detected in ~98-100% of spots as
expected in glioma; EGFR and MKI67 vary by tumor).

## Sources / provenance

### 10x_parent_gbm (pre-existing; fetched by `tmesim/data/spatial_ingest.py`)
- Dataset: 10x Genomics "Human Glioblastoma, Whole Transcriptome Analysis"
  (Visium v1, Space Ranger 1.2.0), a single fresh-frozen human GBM resection.
- Landing: https://www.10xgenomics.com/datasets/human-glioblastoma-whole-transcriptome-analysis-1-standard-1-2-0
- Files (direct, no login):
  - `https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz`
  - `https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_spatial.tar.gz`
- On disk: `data_cache/visium_gbm/{spatial, filtered_feature_bc_matrix}` (root level)
  plus the two source tarballs.

### GBM1 / GBM2 / GBM3 (added; GEO GSE194329)
- Series: **GSE194329** — Ren et al., a spatial transcriptomic atlas of human
  H3K27M-DMG and glioblastoma (10x Visium, Space Ranger 1.1.0, genome hg38).
  Samples are snap-frozen surgical GBM sections; each is an independent primary
  human GBM tumor (three different tissues, distinct from the 10x sample above).
- Series landing: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194329
- Per-sample direct downloads (NCBI GEO FTP, no login) — full Space Ranger output
  tarballs; only the `spatial/` folder and `filtered_feature_bc_matrix/` MEX dir
  were extracted into the cache (the large `molecule_info.h5` /
  `filtered_feature_bc_matrix.h5` were not retained):
  - GBM1: `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833533/suppl/GSM5833533_GBM1_spaceranger_out.tar.gz` (~259 MB)
  - GBM2: `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833534/suppl/GSM5833534_GBM2_spaceranger_out.tar.gz` (~189 MB)
  - GBM3: `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833535/suppl/GSM5833535_GBM3_spaceranger_out.tar.gz` (~252 MB)
- On disk: `data_cache/visium_gbm/<sample_id>/{spatial, filtered_feature_bc_matrix}`.
- Download status: **all three downloaded fully** (byte sizes matched the GEO
  `content-length` exactly) and parse cleanly.

Additional GBM sections available in the same series if more are wanted:
`GBM4` (`GSM5833536`) and `GBM5_1` / `GBM5_2` (`GSM5833537` / `GSM5833538`,
peri-tumoral pair), plus five H3K27M-DMG sections (`DMG1..DMG5`).

## Other open Visium GBM sources noted (not downloaded)
- 10x `Targeted_Visium_Human_Glioblastoma_Pan_Cancer` (Space Ranger 1.2.0,
  direct download, HTTP 200) — a **targeted ~1000-gene panel of the SAME tissue
  block** as `10x_parent_gbm`, so it is a different assay of one tissue rather
  than a new tissue; skipped to keep "multiple real tissues" honest.
  (`Targeted_..._Neuroscience` at the analogous path returned HTTP 403.)

## Reproduce
```bash
./.venv/bin/python -m tmesim.data.real.visium   # lists + summarizes every cached sample
```
