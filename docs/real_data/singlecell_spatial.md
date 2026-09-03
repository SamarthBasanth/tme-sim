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

### Sibling GBM dataset — NOW DOWNLOADED (2026-07-27)
10x **Xenium_V1_Human_Brain_GBM_FFPE** (Xenium v2.0.0), GBM FFPE with the
**Human Immuno-Oncology Profiling Panel**. Cached in `data_cache/xenium_gbm_io/`.
Full details in the section **"Second section: Xenium GBM FFPE, Immuno-Oncology
panel"** at the end of this document.

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

---

# Second section: Xenium GBM FFPE, Immuno-Oncology panel

*Acquired live 2026-07-27 into `data_cache/xenium_gbm_io/`. Same loader,
different cache dir: `load(path="data_cache/xenium_gbm_io")`.*

This is the **richer** of the two real single-cell spatial GBM sections: **20x
the cells** and a **deeper immune/vascular panel**. It does **not** replace the
319-gene brain-panel section — the two panels share only 60 genes and each
carries markers the other lacks. Both are cached; both load.

## Dataset

- **Name:** 10x Genomics **FFPE Human Brain Cancer Data with Human
  Immuno-Oncology Profiling Panel and Custom Add-on**.
- **Assay:** Xenium In Situ, Xenium Analyzer **v2.0.0** (vs v1.3.0 for the
  brain-panel section).
- **Sample id:** `Xenium_V1_Human_Brain_GBM_FFPE`. One adult human GBM FFPE
  section.
- **Panel (measured from the shipped `gene_panel.json` and `features.tsv.gz`):**
  base panel **Xenium Human Immuno-Oncology Profiling** (`hImmune_v1`, **380
  gene targets**) + custom **Add-on B** (**100 gene targets**) = **480
  gene-expression targets**. Note the background figure of "380 genes" is the
  *base* panel only; the shipped matrix carries **480**. Plus 41 negative-control
  codewords and 20 negative-control probes = 541 features total. (Unlike the
  v1.3.0 section, this v2.0.0 matrix has **no** `Unassigned Codeword` features.)
- **Open status:** **OPEN.** CC BY 4.0, no login, no DUA. Direct download from
  `cf.10xgenomics.com` (verified HTTP 200; content-length matched byte-exactly
  on download).
- **Landing page:**
  <https://www.10xgenomics.com/datasets/ffpe-human-brain-cancer-data-with-human-immuno-oncology-profiling-panel-and-custom-add-on-1-standard>

## What downloaded (into `data_cache/xenium_gbm_io/`)

Base URL:
`https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Brain_GBM_FFPE/Xenium_V1_Human_Brain_GBM_FFPE`

| local file | url suffix | bytes | sha256 (first 16) | contents |
|---|---|---|---|---|
| `Xenium_GBM_cells.csv.gz` | `_cells.csv.gz` | 44,907,408 | `f789533ce5913fed` | per-cell centroids + QC metadata (816,769 rows) |
| `Xenium_GBM_cell_feature_matrix.tar.gz` | `_cell_feature_matrix.tar.gz` | 242,459,483 | `b1dbf0305fc3dc04` | MEX triple: `matrix.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz` |
| `Xenium_GBM_gene_panel.json` | `_gene_panel.json` | 256,404 | `d59f40ab391b1f16` | panel identity + 500 target definitions |

**Total 287,623,295 B (~288 MB) downloaded; all three byte-exact against the
advertised `content-length`.** Local basenames deliberately match the v1.3.0
cache so the loader needs no filename special-casing.

Not downloaded: `_xe_outs.zip` (**13.15 GB** — transcripts.parquet, morphology
OME-TIFF, boundary polygons, Xenium Explorer bundle). Not needed for
cell-resolved geometry.

Reproduce: the URLs and verified byte counts are in
`tmesim/data/real/singlecell_spatial.py` as `XENIUM_GBM_IO_BASE` /
`XENIUM_GBM_IO_FILES`.

## Parse result (venv, verified)

- **n_cells = 816,769** (vs 40,887 — **20.0x** more cells).
- **480 panel genes**; expression is a `scipy` CSR matrix, cells x genes =
  `(816769, 480)`, **52,190,559 non-zeros**.
- **Coordinate units: micrometer** (native, same Xenium global frame).
  Measured extent **11,455.4 x 10,836.9 um** (~11.5 x 10.8 mm), x in
  [4.4, 11459.7], y in [3.8, 10840.7]. Bounding box ~124 mm^2, ~6,580 cells/mm^2
  (vs ~41 mm^2 for the v1.3.0 section — **~3x the tissue area as well**).
- Depth per cell: transcript_counts median **157**, mean 174.4, range 0–1630.
  Median `cell_area` 75.1 um^2, median `nucleus_area` 35.3 um^2.
- Negative-control burden (a real QC number): control probes 4,614 and control
  codewords 8,693 out of 1.42e8 total counts = **0.0032% / 0.0061%** — clean.
- `barcodes.tsv.gz` order is **identical** to `cells.csv.gz` order (verified);
  the loader still reindexes by `cell_id` defensively.
- `cells.csv.gz` adds one column vs v1.3.0: `deprecated_codeword_counts` (all
  zero here). The loader reads the table generically, so this is a non-event.
- Load cost: **~8 s, ~1.7 GB peak RSS** on the venv (scipy 1.13.1). Sparse by
  default; `dense=True` would be 816,769 x 480 float32 = **1.57 GB** — prefer
  the CSR.

Reproduce: `python -m tmesim.data.real.singlecell_spatial data_cache/xenium_gbm_io`.

## Panel comparison vs the 319-gene brain panel

The two panels overlap **far less than expected**: only **60 shared genes**.
The IO section **gains 420** genes and **loses 259**.

| | brain panel (v1.3.0) | IO panel (v2.0.0) |
|---|---|---|
| gene targets | 319 | **480** |
| cells | 40,887 | **816,769** |
| tissue bbox | 8.1 x 5.1 mm | **11.5 x 10.8 mm** |
| shared genes | 60 | 60 |
| unique to it | 259 | 420 |

### Target genes explicitly checked

| gene | brain 319 | **IO 480** | note |
|---|---|---|---|
| **MGMT** | absent | **ABSENT** | TMZ methylation determinant — on **neither** panel |
| **ABCB1** (P-gp) | absent | **ABSENT** | efflux — on **neither** panel |
| ABCG2 (BCRP) | absent | **ABSENT** | efflux — on neither panel |
| **GFAP** | absent | **ABSENT** | astro marker — on neither panel |
| PECAM1 (CD31) | present | **LOST** | endothelial — IO uses VWF/FLT1/PLVAP instead |
| CD68 | present | present | |
| CD163 | present | present | |
| AIF1 (IBA1) | present | present | |
| PTPRC (CD45) | present | present | |
| CD8A | present | present | |
| MKI67 | present | present | |
| EGFR | present | present | |
| **CD3E** | absent | **GAINED** | pan-T cell |
| **FOXP3** | absent | **GAINED** | Treg |
| **PDCD1** (PD-1) | absent | **GAINED** | checkpoint |
| **CD274** (PD-L1) | absent | **GAINED** | checkpoint |
| **VEGFA** | absent | **GAINED** | angiogenesis |
| **HIF1A** | absent | **GAINED** | hypoxia TF |
| CA9 | absent | **ABSENT** | hypoxia marker — on neither panel |
| COL1A1 | absent | **ABSENT** | (brain panel has COL1A2/COL12A1; IO has none) |
| TOP2A | present | **LOST** | topo-II — see "what this costs" below |
| SOX2 | present | **LOST** | glioma stem/lineage |
| OLIG2 | present | **LOST** | glioma lineage |

### What the IO panel adds (420 genes) — by axis

- **T / NK cells and checkpoints (new):** CD3D, CD3E, CD247, CD4, CD8B, CD2,
  CD27, CD28, CD38, CD69, LCK, TRAT1, TCF7, TBX21, EOMES, RORC, RUNX3, TOX,
  ZNF683, ITGAE, ENTPD1, **FOXP3, CTLA4, PDCD1, PDCD1LG2, CD274, LAG3, TIGIT,
  HAVCR2, VSIR, CD276, IDO1**, GZMB/GZMH/GZMK, PRF1, KLRC1/D1/F1/K1, FGFBP2,
  SPON2, SAMD3, CTSW, IL2RA/2RB/2RG, IL7R, TNFRSF9, ICOSLG, FASLG, FAS.
- **B / plasma cells (new, entirely absent before):** MS4A1, CD19, CD79A, CD79B,
  PAX5, BANK1, TCL1A, FCRL3, FCER2, TNFRSF13C, MZB1, JCHAIN, IGKC, IGLC3,
  IGHM, IGHG1-4, TNFRSF17, SDC1, DERL3, XBP1, FKBP11, PRDX4, SSR4.
- **Myeloid / TAM depth (new):** CD14, MRC1 (CD206), MARCO, CSF1, CSF1R, MPEG1,
  MNDA, FCN1, CLEC10A, CLEC12A, LILRA4, LILRB2, LILRB4, MCEMP1, CHIT1, MMP9,
  MMP12, S100A9, S100A12, ARG1, PLA2G7, ANPEP; dendritic-cell programs BATF3,
  CD1A, CD1C, LAMP3, SPIB, TCF4, PLD4; mast CPA3, MS4A2, HPGDS.
- **Vascular / perivascular niche (new):** VWF, FLT1, PLVAP, RAMP2, ACKR1, SELE,
  ICAM1, HIGD1B, MMRN1, TMEM100, GNG11, CAV1, CAVIN2 — plus mural RGS5, PDGFRB,
  NOTCH3, ACTA2, MYH11, DES.
- **Hypoxia / angiogenesis (new):** **VEGFA, HIF1A**, EDN1, EDNRB, GDF15,
  CDKN1A, EPO, EPOR, IGFBP7, FGF2, A2M.
- **Chemokine / cytokine axes (new, ~70 genes):** the CCL2/3/4/5/7/8/11/13-21/26/28
  and CXCL1/2/3/5/6/9/10/11/12/13/14/16 ligand families with their receptors
  CCR1-7, CXCR1-6; IL1A, IL1B, IL2, IL4, IL5, IL6, IL10, IL12B, IL13, IL15,
  IL18, IL23A, IL13RA2, TNF, TGFB1; and the signalling layer JAK1/2/3,
  STAT1/2/3/4/5B/6, IRF1-8, SOCS1-3, ISG15, MX1, IFIT2/3, IFITM3, MYD88.
- **Oncogenic drivers (new or newly complete):** PTEN, TP53, RB1, NF1,
  CDKN2A/2B/2C/2D, CDKN1A/1B/1C, MYC, KRAS, BRAF, PIK3CA, MET, ERBB2/3/4,
  NOTCH1/2/3, CTNNB1, ATM, ARID1A, KMT2C, CDK1/2/4/6/12, CCND1, JAK/STAT above;
  proliferation STMN1, UBE2C, CENPF, ORC6, TUBA1B, TUBB.
- **Glial context retained:** S100B, SLC1A3, PTN, PTGDS, SPARCL1, APOE, SOX9,
  VCAN, ID4, TREM2, CX3CR1, GPR34, C1QB.

### What this costs (259 genes lost) — do NOT discard the brain panel

The IO panel is immune-designed, so it drops most of the neuro/glioma-lineage
program that the brain panel carries:

- **Glioma lineage / stemness:** SOX2, SOX10, SOX4, SOX11, OLIG1, OLIG2,
  PTPRZ1, BCAN, CSPG4 (NG2), NES, PAX6, TTYH1.
- **IDH1 and IDH2** — the WHO-defining glioma classifier genes. Present on the
  brain panel, **absent** on the IO panel.
- **Astrocyte:** AQP4, GJA1, CHI3L1, SERPINA3.
- **Oligodendrocyte program (complete loss):** MOG, MAG, MOBP, MYRF, CLDN11,
  OPALIN, ERMN, ST18, UGT8, KLK6.
- **Homeostatic microglia:** P2RY12, P2RY13.
- **Neuronal:** GAD1, GAD2, SLC17A6, SLC17A7, RELN, SNCA, NRGN, VIP, NPTX1.
- **Proliferation / drug-target genes used by the current sensitivity model:**
  **TOP2A** and **PCNA**.
- **ECM:** COL1A2, COL12A1, COL25A1, POSTN, THBS1, LOX, SULF1, TGFBI.
- **Other hypoxia/stress:** HILPDA, HMOX1.

**Practical read:** the brain-panel section stays the reference for tumour
lineage, ECM and oligodendrocyte biology; the IO section becomes the reference
for immune, vascular and hypoxia biology. They are complementary sections of
**different patients**, so they cannot be merged cell-for-cell — only compared
as two independent real microenvironments.

## Loader compatibility

`tmesim.data.real.singlecell_spatial.load(path="data_cache/xenium_gbm_io")`
reads this dataset with **no change to the parsing path** — the MEX read, the
`cell_id` alignment, the `Gene Expression` feature filter and the um coordinates
all work as written. Verified: both sections load, and the full test suite is
green (57 passed).

One **provenance** fix was required, because `tmesim/reconstruct/real_xenium.py`
stamps `data["source"]["dataset"]` and `["panel"]` into downstream result
metadata. That dict was hardcoded to the v1.3.0 brain-panel section, so loading
this cache would have labelled IO-panel results as "hBrain_254g ... 319 gene
targets" — a silent provenance falsification. `source` is now resolved from the
**base-panel name inside the cache's own `gene_panel.json`**:

- `hBrain_254g` -> the v1.3.0 brain-panel entry (string-identical to before, so
  existing results/metadata are unchanged);
- `Xenium Human Immuno-Oncology Profiling` -> the IO entry;
- anything else -> an explicitly `UNREGISTERED ...` source rather than
  inheriting someone else's provenance.

`source` additionally now carries `cache_dir` and `n_panel_genes`.

### Known downstream gap — CLOSED since this was written (2026-09-03)

> **This subsection is superseded and is kept as the record of the gap, not of the
> state.** Every clause below has been overtaken: `reconstruct/marker_panels.py:137`
> now defines `IO_PANEL_MARKER_MAP` using **exactly the substitutes named here as
> available-but-out-of-scope** (endothelial VWF/FLT1/PLVAP; astro-lineage
> SLC1A3/S100B/PTN; IL13RA2 for tumour), it is registered as `io_480` at `:233`,
> `select_marker_map()` picks it for this section, and `results/classifier_fix.md`
> scores the IO section with **four** classes — TUMOR, STROMAL, IMMUNE,
> ENDOTHELIAL — not two. Making the marker map panel-aware was not "left open"; it
> was done.

`real_xenium.MARKER_MAP` was written against the brain panel. On the IO panel
its ENDOTHELIAL marker (PECAM1) and all three STROMAL markers (MOG, AQP4, GJA1)
are absent, as are TUMOR's SOX2/OLIG2. Running `phenotype_cells` on the IO
section therefore yields only two classes — **TUMOR 760,223 (93.1%) / IMMUNE
56,546 (6.9%)** — with `missing = [AQP4, GJA1, MOG, OLIG2, PECAM1, SOX2]`. This
is a real limitation of the *existing* marker map, not of the data: the IO panel
has strictly better substitutes available (endothelial VWF/FLT1/PLVAP/ACKR1;
mural/stromal RGS5/PDGFRB/NOTCH3/ACTA2; astro-lineage tumour
SLC1A3/S100B/PTN/SPARCL1; plus IL13RA2 as a glioma-specific tumour marker).
Making the marker map panel-aware is out of scope for this ingest and is left
open.

## What this unlocks

Measured per-cell prevalence on the 816,769 cells (fraction of cells with >=1
transcript) — these are real detections, not estimates:

| gene | % cells positive | mean counts/cell |
|---|---|---|
| S100B | 96.95% | 11.196 |
| EGFR | 95.81% | 10.321 |
| SLC1A3 | 90.83% | 5.320 |
| VEGFA | 63.36% | 3.154 |
| HIF1A | 54.93% | 0.988 |
| PDGFRA | 51.77% | 1.170 |
| MKI67 | 33.73% | 0.792 |
| CD163 | 18.20% | 0.485 |
| RGS5 | 12.45% | 0.467 |
| FLT1 | 9.42% | 0.408 |
| TREM2 | 9.11% | 0.130 |
| IL13RA2 | 7.12% | 0.081 |
| CD68 | 6.16% | 0.076 |
| PTPRC | 4.97% | 0.062 |
| VWF | 4.15% | 0.123 |
| CD274 (PD-L1) | 2.37% | 0.027 |
| IDO1 | 1.48% | 0.026 |
| CD8A | 0.99% | 0.014 |
| TIGIT | 0.88% | 0.009 |
| LAG3 | 0.68% | 0.007 |
| CD3E | 0.61% | 0.009 |
| CTLA4 | 0.59% | 0.006 |
| FOXP3 | 0.51% | 0.005 |
| PDCD1 (PD-1) | 0.08% | 0.001 |

Newly possible, on measured bytes:

1. **A real hypoxia -> angiogenesis field.** VEGFA (63% of cells) and HIF1A
   (55%) are both measured per cell at um resolution. The vascular/perivascular
   readout (VWF, FLT1, PLVAP, RGS5, PDGFRB) is co-measured, so a
   *hypoxia-gradient -> VEGFA -> vessel-proximity* relationship can be tested as
   a spatial correlation on real data rather than assumed. The current
   brain-panel section measures **none** of VEGFA/HIF1A.
2. **A genuine immune niche.** For the first time: T cells (CD3E/CD3D/CD8A/CD4),
   Tregs (FOXP3/CTLA4/IL2RA), exhaustion (PDCD1/LAG3/TIGIT/HAVCR2/TOX), and
   ligand side (CD274/PDCD1LG2/IDO1) in the same coordinate frame. A real
   "immune-cold vs immune-infiltrated region" call becomes measurable — and the
   numbers above already tell the honest story: this tumour is **immune-cold**
   (T cells 0.6% of cells; PD-1+ 0.08%) with a **dominant myeloid/TAM**
   compartment (CD163 18%, TREM2 9%). Any checkpoint-axis modelling must respect
   how rare those cells are; at 625 PDCD1+ cells, per-region statistics will be
   noise-dominated.
3. **Cell-cell interaction modelling with real ligand-receptor pairs.** ~70
   chemokine/cytokine genes with matched receptors (CXCL12/CXCR4, CCL2/CCR2,
   CXCL9-10-11/CXCR3, CSF1/CSF1R) at single-cell um resolution — the input for
   spatial ligand-receptor inference the 319 panel cannot support.
4. **20x the cells, 3x the area.** 816,769 cells over ~124 mm^2 gives enough
   statistical power for rare-population and niche-level spatial statistics that
   40,887 cells cannot support.
5. **Driver-gene heterogeneity:** PTEN, TP53, NF1, CDKN2A, EGFR, PDGFRA, MET
   co-measured per cell — real intratumoural subclonal-marker structure.

### What it does NOT unlock (be blunt)

- **TMZ sensitivity is STILL UNMODELLED *on this section*.** **MGMT is absent from
  this panel**, exactly as it is from the 319-gene panel, so temozolomide
  sensitivity remains **unmeasured here** and any per-cell TMZ sensitivity on *this*
  dataset would be invented. **Corrected 2026-09-03:** this bullet generalised that
  to "no Xenium panel in this project measures MGMT", which is no longer true —
  `data_cache/gse324268_xenium/` is a 366-gene Xenium GBM panel that **does** carry
  MGMT (`ENSG00000170430`, present in `GSE324268_gene_panel.json.gz`), it is a
  registered section (`gse324268` in `SECTION_RECIPES`), and
  `docs/outreach/04_holzel_gse324268.md` is built entirely on it. The gap is
  panel-specific, not project-wide.
- **Efflux is STILL UNMODELLED.** **ABCB1 (P-gp) and ABCG2 (BCRP) are both
  absent.** No real per-cell efflux term is available.
- **GFAP is absent** (as on the brain panel). Astrocytic identity must still be
  proxied — here by SLC1A3/S100B/PTN/SPARCL1 (the brain panel's AQP4/GJA1 are
  also gone).
- **CA9 and COL1A1 are absent**, so the canonical hypoxia-marker and
  fibrillar-collagen readouts are proxies (HIF1A/VEGFA; and no fibrillar
  collagen at all on this panel — the brain panel's COL1A2/POSTN is the better
  ECM source).
- **TOP2A is absent**, so the existing TOP2A-based doxorubicin/etoposide
  sensitivity determinant **cannot be evaluated on this section**. MKI67, EGFR,
  TP53 and CDK1 remain, so the cell-cycle-dependent drugs still have their
  measured determinant.
- **No paired drug-response outcome.** Same gap as every other section here: one
  untreated resection, no RECIST / survival label. Mechanistic prediction only,
  **not** supervised response validation.
- **One section, one patient**, and a *different* patient from the 319-gene
  section — no cross-panel per-cell join is possible.
