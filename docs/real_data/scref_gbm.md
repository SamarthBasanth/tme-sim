# Real GBM single-cell RNA-seq reference atlas

Deconvolution reference for resolving spatial (Visium) spots into real
cell-type fractions and malignant vs. non-malignant content.

## Source / accession

- **Study:** Darmanis S, Sloan SA, Croote D, et al. (2017), *Cell Reports*
  21(5):1399-1410. "Single-Cell RNAseq analysis of diffuse neoplastic
  infiltrating cells at the migrating front of human glioblastoma."
- **Accession:** GEO **GSE84465** (NCBI Gene Expression Omnibus).
- **Platform:** SmartSeq2, plate-based full-length scRNA-seq (GPL16791,
  Illumina HiSeq 2500), HTSeq raw counts.
- **Cohort:** 4 primary human GBM patients (BT_S1, BT_S2, BT_S4, BT_S6),
  sampling both the tumour core (`Tumor`) and peritumoral brain (`Periphery`).

## License / open status

Public / open. NCBI GEO is a public repository and this series is openly
released (no controlled-access / dbGaP tier). Files are downloaded directly
over the public GEO FTP/HTTPS endpoint with no authentication. Verified
downloadable with `curl` before caching.

## Files cached (`data_cache/scref_gbm/`)

| File | Size | Format | Contents |
|------|------|--------|----------|
| `GSE84465_GBM_All_data.csv.gz` | ~20.5 MB | gzip, **whitespace-delimited** text (despite `.csv`) | genes x cells raw HTSeq counts; first column = gene symbol, column headers = `<plate_id>.<well>` cell barcodes |
| `GSE84465_series_matrix.txt.gz` | ~146 KB | GEO series matrix | per-cell annotations (cell type, tissue, patient, neoplastic flag, tSNE cluster) |

Direct URLs (verified working):
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_GBM_All_data.csv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/matrix/GSE84465_series_matrix.txt.gz`

The two files are joined on the `<plate_id>.<well>` barcode (the counts-matrix
column names, e.g. `1001000173.G8`), which the series matrix reconstructs from
its `plate id` + `well` characteristics.

## Shape

- **n_genes = 23,460** (after dropping 5 trailing HTSeq summary rows:
  `no_feature`, `ambiguous`, `too_low_aQual`, `not_aligned`,
  `alignment_not_unique`).
- **n_cells = 3,589** (100% of cells join cleanly to annotations).
- Values are raw integer UMI/read counts (min 0, observed max in the thousands).

## No subsampling

The full matrix is small (~20 MB gzipped, 3,589 cells) so the **entire dataset
is cached and loaded** — no subsample. `load_reference(max_cells=N)` offers an
optional documented uniform subsample for quick experiments, off by default.

## Cell-type label vocabulary (real, author-provided)

Field `cell_type` in the returned labels (verbatim from GEO metadata):

| cell type | n cells | role in the niche |
|-----------|--------:|-------------------|
| `Immune cell` | 1847 | tumour-associated macrophages / microglia + lymphocytes |
| `Neoplastic` | 1091 | **malignant** GBM cells |
| `OPC` | 406 | oligodendrocyte precursor cells |
| `Astocyte` | 88 | astrocytes (note: `Astocyte` is the original GEO spelling / author typo) |
| `Oligodendrocyte` | 85 | oligodendrocytes |
| `Vascular` | 51 | endothelium / vasculature |
| `Neuron` | 21 | neurons |

Additional real per-cell fields also returned in `cell_labels`:
- `neoplastic` — binary malignant state: `Neoplastic` (1091) vs `Regular` (2498).
  This is the dataset's **malignant-state** annotation. (This atlas does not
  carry Neftel-style AC/MES/OPC/NPC meta-modules; those live in the
  login-gated Broad SCP503 / GSE131928 per-cell metadata. The honest scope
  here is malignant-vs-normal plus the 7 cell types above.)
- `tissue` — `Tumor` (2343) vs `Periphery` (1246): sampling location.
- `patient_id` — BT_S1 / BT_S2 / BT_S4 / BT_S6.
- `tsne_cluster`, `plate_id`, `well`, `gsm`.

## Loader

`tmesim/data/real/scref_gbm.py`:

```python
from tmesim.data.real import load_reference
expression, cell_labels, gene_names = load_reference()
#   expression  : pandas DataFrame, genes x cells raw counts (23460 x 3589)
#   cell_labels : pandas DataFrame aligned to expression.columns
#   gene_names  : np.ndarray of gene symbols
```

Dependency-light (pandas + scipy/stdlib only; no scanpy/anndata).

## How it serves as a deconvolution reference

For each spatial Visium spot (bulk mixture of a few cells), the reference gives
a per-cell-type expression signature: group these 3,589 annotated single cells
by `cell_type` (and by the `neoplastic` malignant flag) to build reference
signature profiles, then estimate the cell-type / malignant fractions that best
explain each spot's expression (NNLS / regression-based deconvolution). The
`Neoplastic` vs `Regular` flag additionally lets us separate malignant tumour
content from the normal-brain + immune + vascular stroma of the niche, which is
exactly the L2 "spot -> real cell-type composition" step in the build plan.
