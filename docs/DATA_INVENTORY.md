# DATA_INVENTORY — the backup record for `data_cache/`

`data_cache/` is **gitignored with nine deliberate exceptions**, and holds **3.9 GB**
of real open data that nothing in this repository can regenerate from source code.
**This file is the backup record**: if `data_cache/` is deleted, everything below must
be re-downloadable from what is written here.

> **Re-measured 2026-09-03.** This line said 6.1 GB and said the directory was
> gitignored without qualification. Both were wrong. `du -sh data_cache` → **3.9G**.
> `git ls-files data_cache` → **9 tracked files**: the `brain_slice/` and
> `pruijn2008/` provenance notes and CSVs, which are small, hand-curated and
> deliberately versioned. Every size in this file is a measurement with a date, not a
> constant — re-measure before quoting one.

Scope: the datasets recorded below, with path, source/accession, licence, size, what
it enables, and **whether it is reproducible via `python -m tmesim.cli acquire` or
needs a manual/scripted fetch**.

> **This is NOT "every dataset on disk", which is what this line claimed until
> 2026-09-03.** Checked by listing `data_cache/` and grepping each entry against this
> file: **9 top-level entries are not named here at all** — `gbmspace_visium`,
> `greenwald_gbm_visium`, `gse324268_xenium`, `xenium_breast_janesick`,
> `gdsc_breast_dose_response.csv` and four `*_runs.json` caches. The most serious is
> **`xenium_breast_janesick`**, which `tmesim/validation/scorecard.py:2686` declares a
> hard `needs=`, so an anchor cannot run without it — and this backup record carries no
> URL, no size and no recipe for it. **If `data_cache/` were deleted today, this file
> would not be enough to rebuild it.** Re-run that listing before trusting the coverage.

| Tier | Path | Size | Needed for |
|---|---|---:|---|
| **A — core pipeline** | `data_cache/*` (excluding `hunt/`) | **3.6 GiB** (3.87 GB decimal), i.e. nearly all of it | Everything. Without it no reconstruction, no simulation, no internal validation runs. |
| **B — external validation corpus** | `data_cache/hunt/` | **360 MB**, 24 top-level dirs, 1,217 files | The audit-era external validations (`glioma_regional_pk`, `gist_imatinib`, `moncayo_explant`, `gbm_tmz_exvivo`) and the landscape survey in `results/DATA_HUNT.md`. |

> **Both rows re-measured 2026-09-03** (`du -sh`, `find … -type d/-type f`). Tier A was
> given as ~0.9 GB and Tier B as 5.3 GB / 38 dirs / ~1,950 files. The two were close to
> **inverted**: the bulk is in Tier A, not in the hunt corpus.

Three acquisition classes are used throughout:

| Class | Meaning |
|---|---|
| **AUTO** | `python -m tmesim.cli acquire` fetches it. Includes the `xenium_319` tissue section since 2026-08-15 (`--sections`); the four tabular sources alone leave every section unavailable. Measured on a clean clone 2026-08-15: **1 m 50 s and 3 m 52 s on two runs**, 5/5 sources OK — the earlier "17 s" in this row did not survive re-measurement and timing is dominated by the GDC endpoints, so treat it as minutes, not seconds. |
| **MANUAL** | Open, no login, no DUA — but not a copy-paste-and-go: too large, or needing extraction, a format conversion, or a companion dataset. A `curl` recipe is given below, and `python -m tmesim.cli sections` prints the same commands per section. |
| **SCRIPTED** | Needs a harvest script or an API walk. The script, endpoint or accession list is given. |

---

# TIER A — core pipeline data

Everything here is **open**: no login, no DUA, no MTA. Per-dataset narrative
provenance (what was deliberately *not* downloaded, byte counts, sha256 prefixes)
is in [`real_data/`](real_data/).

## A1. `data_cache/xenium_gbm/` — 30 MB — **AUTO** — *the primary tissue*

> Re-checked 2026-09-03: `SECTION_RECIPES['xenium_319']['mode']` is `auto`, so
> `vivum acquire` fetches this section. This heading said **MANUAL**.

| | |
|---|---|
| **Dataset** | 10x Genomics **Xenium In Situ**, FFPE Human Brain Glioblastoma with Add-on (Xenium Analyzer v1.3.0) |
| **Sample id** | `Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon` |
| **Licence** | **CC BY 4.0**, open direct download |
| **Contents** | **40,887** segmented cells, measured µm centroids, **319** gene targets (hBrain_254g + 65 add-on), 8.1 × 5.1 mm |
| **Enables** | Every headline result. The flagship reconstruction, the drug panel, `delivery_metrics`, `vessel_profile`, `vascular_validation`, `numerical_rigor`, `uncertainty`, and the Xenium half of `reconstruction_fidelity`. |
| **Loader** | `tmesim.data.real.singlecell_spatial.load()` (this is the default cache) |
| **Provenance doc** | [`real_data/singlecell_spatial.md`](real_data/singlecell_spatial.md) |

```bash
mkdir -p data_cache/xenium_gbm
BASE=https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon
curl -L -o data_cache/xenium_gbm/Xenium_GBM_cells.csv.gz               "${BASE}_cells.csv.gz"                # 2,167,253 B
curl -L -o data_cache/xenium_gbm/Xenium_GBM_cell_feature_matrix.tar.gz "${BASE}_cell_feature_matrix.tar.gz"  # 14,279,458 B
curl -L -o data_cache/xenium_gbm/Xenium_GBM_gene_panel.json            "${BASE}_gene_panel.json"             # 152,134 B
```

The loader extracts `cell_feature_matrix.tar.gz` on first use into
`cell_feature_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}`.

## A2. `data_cache/xenium_gbm_io/` — 536 MB — **AUTO** — *the deep panel*

> Re-checked 2026-09-03: `SECTION_RECIPES['xenium_io']['mode']` is `auto`
> (reachable via `vivum acquire --sections all`). This heading said **MANUAL**.

| | |
|---|---|
| **Dataset** | 10x Genomics **Xenium In Situ**, FFPE Human Brain Cancer with Human Immuno-Oncology Profiling Panel and Custom Add-on (Xenium Analyzer v2.0.0) |
| **Sample id** | `Xenium_V1_Human_Brain_GBM_FFPE` |
| **Licence** | **CC BY 4.0**, open direct download |
| **Contents** | **816,769** cells (20× the 319-panel section), **480** gene targets (hImmune_v1 380 + Add-on B 100), 11.5 × 10.8 mm (~3× the area), 52,190,559 non-zeros |
| **Enables** | The deep immune/vascular readout. Gains **CD3E, FOXP3, PDCD1, CD274, VEGFA, HIF1A**; loses PECAM1 (uses VWF/FLT1/PLVAP instead). Only **60 genes** overlap with the 319 panel. It is the substrate for the not-yet-written oxygen-vs-VEGFA/HIF1A validation. |
| **Loader** | `singlecell_spatial.load(path="data_cache/xenium_gbm_io")` — **the directory is the only selector**; both caches use identical local basenames and provenance is read out of the shipped `gene_panel.json` |
| **Landing page** | <https://www.10xgenomics.com/datasets/ffpe-human-brain-cancer-data-with-human-immuno-oncology-profiling-panel-and-custom-add-on-1-standard> |

```bash
mkdir -p data_cache/xenium_gbm_io
BASE=https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Brain_GBM_FFPE/Xenium_V1_Human_Brain_GBM_FFPE
curl -L -o data_cache/xenium_gbm_io/Xenium_GBM_cells.csv.gz               "${BASE}_cells.csv.gz"                # 44,907,408 B  sha256 f789533ce5913fed…
curl -L -o data_cache/xenium_gbm_io/Xenium_GBM_cell_feature_matrix.tar.gz "${BASE}_cell_feature_matrix.tar.gz"  # 242,459,483 B sha256 b1dbf0305fc3dc04…
curl -L -o data_cache/xenium_gbm_io/Xenium_GBM_gene_panel.json            "${BASE}_gene_panel.json"             # 256,404 B     sha256 d59f40ab391b1f16…
```

**Deliberately not downloaded:** `_xe_outs.zip` (**13.15 GB** — transcripts.parquet,
morphology OME-TIFF, boundary polygons, Xenium Explorer bundle). Not needed for
cell-resolved geometry. The canonical URL list also lives in code as
`singlecell_spatial.XENIUM_GBM_IO_BASE` / `XENIUM_GBM_IO_FILES`.

## A3. `data_cache/visium_gbm/` — 243 MB — **MANUAL** — *4 spot-level sections*

| | |
|---|---|
| **Datasets** | 10x public "parent" Visium GBM section + **GEO GSE194329** samples GBM1/GBM2/GBM3 (three further primary human GBM sections) |
| **Licence** | 10x: CC BY 4.0. GEO GSE194329: open, no DUA |
| **Contents** | 1,775–4,764 spots each, whole transcriptome, `spatial/` positions + scalefactors |
| **Enables** | The Visium NNLS-deconvolution half of `reconstruction_fidelity` — **the cross-platform composition check, whose Pearson 0.997 is RETRACTED**: after the classifier-arity fix the 4-vector reads 0.980 and the non-tumour subvector collapses **0.945 → 0.016**, indistinguishable from no agreement (`results/revalidation.md` §4). Concordance here is an open question, not evidence. Also `experiments/cross_tissue.py` (5 tumours at matched dx), `experiments/niche_response.py`, `patient_report --source visium` |
| **Loader** | `tmesim.data.real.visium.list_samples() / load_sample(id)` — ids `10x_parent_gbm`, `GBM1`, `GBM2`, `GBM3` |
| **Provenance doc** | [`real_data/visium_gbm.md`](real_data/visium_gbm.md) |

```bash
mkdir -p data_cache/visium_gbm && cd data_cache/visium_gbm
# 10x parent section (~63 MB)
BASE=https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma
curl -L -O "${BASE}_filtered_feature_bc_matrix.tar.gz"    # 54,134,945 B
curl -L -O "${BASE}_spatial.tar.gz"                       #  9,134,896 B
tar xzf Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz
tar xzf Parent_Visium_Human_Glioblastoma_spatial.tar.gz

# GSE194329 GBM1/2/3 (~700 MB of tarballs; keep only spatial/ + filtered_feature_bc_matrix/)
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833533/suppl/GSM5833533_GBM1_spaceranger_out.tar.gz  # 258,814,702 B
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833534/suppl/GSM5833534_GBM2_spaceranger_out.tar.gz  # ~189 MB
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833535/suppl/GSM5833535_GBM3_spaceranger_out.tar.gz  # ~252 MB
for s in GBM1 GBM2 GBM3; do mkdir -p "$s"; tar xzf GSM5833*_${s}_spaceranger_out.tar.gz -C "$s" --strip-components=1; done
cd ../..
```

## A4. `data_cache/scref_gbm/` — 21 MB — **MANUAL** — *the deconvolution reference*

| | |
|---|---|
| **Dataset** | Darmanis et al. 2017, **GEO GSE84465** — human GBM scRNA-seq atlas |
| **Licence** | Open (GEO), no DUA |
| **Contents** | 23,460 genes × 3,589 annotated cells, 7 cell types |
| **Enables** | Every Visium reconstruction. `reconstruct.deconvolve` builds its NNLS signatures from this; without it the Visium path cannot run at all. |
| **Loader** | `tmesim.data.real.scref_gbm.load_reference()` |
| **Provenance doc** | [`real_data/scref_gbm.md`](real_data/scref_gbm.md) |

```bash
mkdir -p data_cache/scref_gbm
curl -L -o data_cache/scref_gbm/GSE84465_GBM_All_data.csv.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_GBM_All_data.csv.gz   # 20,519,629 B
curl -L -o data_cache/scref_gbm/GSE84465_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/matrix/GSE84465_series_matrix.txt.gz #    145,533 B
```

## A5. `data_cache/ivy_gap/` — 36 MB on disk; `data_cache/ivygap/` **ABSENT** — **AUTO (+ unzip) / MANUAL**

> Re-measured 2026-09-03. `data_cache/ivy_gap/` is 36 MB and holds the **unextracted**
> `gene_expression_matrix.zip`. `data_cache/ivygap/` — the ~89 MB extracted tree the
> loader actually reads — **does not exist**, so `load_structures()` raises
> `FileNotFoundError` today. The unzip step in this section has not been run on this
> machine. See [`real_data/ivygap.md`](real_data/ivygap.md), which now carries the same
> status and the list of which files the zip does and does not contain.

| | |
|---|---|
| **Dataset** | **Ivy Glioblastoma Atlas Project** (Allen Institute; Puchalski et al. 2018) |
| **Licence** | Open, Allen Institute terms of use |
| **Contents** | 270 RNA-seq samples across 7 anatomic structures, 37 tumours |
| **Enables** | `reconstruct.region_label` — the anatomic niche zones (leading edge / cellular tumour / microvascular proliferation / perinecrotic), used by `experiments/niche_response.py` and the region permutation tests in `reconstruction_fidelity` |
| **Loader** | `tmesim.data.real.ivygap.load_structures()` / `structure_signatures()` |
| **Provenance doc** | [`real_data/ivygap.md`](real_data/ivygap.md) |

`acquire` fetches the zip into `data_cache/ivy_gap/`; the loader reads the
**extracted** files from `data_cache/ivygap/`:

```bash
python -m tmesim.cli acquire                      # -> data_cache/ivy_gap/gene_expression_matrix.zip
mkdir -p data_cache/ivygap
unzip -o data_cache/ivy_gap/gene_expression_matrix.zip -d data_cache/ivygap
# -> fpkm_table.csv, rows-genes.csv, columns-samples.csv, README.txt
```

Four summary tables are **not in the zip** and need a separate fetch:

```bash
for n in rna_seq_samples_details sub_block_details gene_expression_details tumor_details; do
  curl -L -o "data_cache/ivygap/$n.csv" "https://glioblastoma.alleninstitute.org/api/v2/gbm/$n.csv"
done
```

## A6. `data_cache/gdsc_glioma_dose_response.csv` — 2.0 MB — **AUTO**

| | |
|---|---|
| **Dataset** | **GDSC2** fitted dose-response, glioma (GBM/LGG) subset, release 8.5 |
| **Licence** | Open (Sanger COG bucket) |
| **Contents** | 12,001 glioma dose-response rows |
| **Enables** | The measured-potency axis: every IC50 the molecular layer uses. `molecular.potency`, `pharmacology_real`, `ranking_robustness`, `patient_report`. |
| **Source** | `https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/` (walked by `tmesim/data/acquire.py`) |

> **⚠ `cancerrxgene.org` returned HTTP 410 on 2026-07-27.** Do not cite the GDSC
> web portal as a live URL. The Sanger COG bucket above still serves the file and
> is what `acquire` uses.

## A7. `data_cache/drugresponse/` — 1.5 MB — **MANUAL**

| | |
|---|---|
| **Datasets** | **HGCC** portal JSON (Johansson et al. 2020, PMID 32668248) + the Johansson `mmc2–mmc7` and Jacob 2020 `TableS1–S6` supplements |
| **Licence** | Open; HGCC portal is public |
| **Contents** | 262 drugs × 85 patient-derived GBM cultures |
| **Enables** | `pharmacology_real` — the real-vs-real potency corroboration (GDSC vs HGCC, Spearman +0.83, p = 2.4e-5), the only validation in the repo with no model in the loop |
| **Loader** | `tmesim.data.real.drug_response.load_hgcc()` |
| **Provenance doc** | [`real_data/drug_response.md`](real_data/drug_response.md) |

```bash
mkdir -p data_cache/drugresponse
curl -L -o data_cache/drugresponse/hgcc_drugresponse.json http://hgcc.se/portal/json/drugsummary.json
curl -L -o data_cache/drugresponse/hgcc_bootstrap.json    http://hgcc.se/portal/json/bootstrap.json
# also available: http://hgcc.se/portal/json/{expr,cna,mutation}.json
```

The `johansson2020_mmc*.xlsx` and `jacob2020_TableS*.xlsx` files are journal
supplements; see `real_data/drug_response.md` for the article DOIs.

## A8. `data_cache/tcga_gbm_clinical.json` + `data_cache/tcga_expression/` — 68 KB + 30 MB — **AUTO**

| | |
|---|---|
| **Dataset** | **TCGA-GBM** via the GDC REST API |
| **Licence** | Open tier, no dbGaP needed for these endpoints |
| **Contents** | 80 clinical cases; **7** open STAR gene-count files + `manifest.json` (counted 2026-09-03; this said 8, which was the directory entry count including the manifest) |
| **Enables** | Cohort context only. **No current validation module depends on it** — it is background, not load-bearing. |
| **Source** | `https://api.gdc.cancer.gov/{cases,files,data}` |

## A9. PK anchors — **NO DOWNLOAD** (hard-coded, PMID-cited)

The tissue-PK ground truth is **not a cached file**. It lives in
`tmesim/validation/pk_anchors.py` as a registry of PMID-verified literature
values (penetration lengths and brain:plasma ratios for 6 drugs), each carrying
its PMID. Inspect with `python -m tmesim.validation.pk_anchors`. The anchor that
sets the single calibration constant α is **Primeau 2005, PMID 16361566**
(doxorubicin perivascular half-distance). Narrative provenance:
[`real_data/pk_anchors.md`](real_data/pk_anchors.md).

---

# TIER B — `data_cache/hunt/` — the external-validation corpus

**360 MB, 24 top-level directories, 1,217 files** — re-measured 2026-09-03; this read
5.3 GB / 38 directories / ~1,950 files, and the corpus on disk is roughly **15× smaller**
than that. Treat the per-dataset entries below as a record of what the hunt *retrieved*,
not as a manifest of what is on disk today: files listed here may have been pruned since,
and at least one named below has. Produced by
`experiments/data_hunt_workflow.js`: **28 hunter angles** (one directory each) +
**10 independent verification agents** (`verified_0` … `verified_9`) + a synthesis
pass → [`../results/DATA_HUNT.md`](../results/DATA_HUNT.md), which surveyed
**164 datasets, 143 of them open**.

> **Read this before trusting the section below.** Only **8 of 38** directories
> contain a machine-readable manifest, README or provenance JSON. For the other
> 30, re-download depends on accessions embedded in filenames plus
> `results/DATA_HUNT.md`. **Nothing records a checksum or an exact `curl`
> invocation for any hunt file.** URLs below are marked:
> **[ON DISK]** recorded verbatim in a cache file · **[HUNT]** recorded only in
> `results/DATA_HUNT.md` (all `curl`-checked HTTP 200 on 2026-07-27) ·
> **[DERIVED]** reconstructed here from an accession — verify before relying on it.

## B0. The four directories a validation module actually reads

Only these are load-bearing for code in `tmesim/validation/`. Everything else in
Tier B is survey material.

| Directory | Read by | Status |
|---|---|---|
| `autoradiography_microdialysis/` | `validation.glioma_regional_pk` | **runs** → `results/glioma_regional_pk.md` |
| `maldi_msi_drug/` + `verified_0/` | `validation.gist_imatinib` | **runs** → `results/gist_imatinib.md` |
| `geo_organoid_screens/` | `validation.gbm_tmz_exvivo` | **runs** → `results/gbm_tmz_exvivo.md` |
| `geomx_dsp/` | `validation.moncayo_explant` | **runs** → `results/moncayo_explant.md` |

---

## B1. `autoradiography_microdialysis/` — 5.2 MB — ★ **the best glioma delivery asset**

**SCRIPTED** (literature curation + EuropePMC supplement fetch).

Powers the project's first validation against measured drug concentrations in
human brain tumour tissue, and the derivation of the additive-leak barrier fix
(L ≈ 0.24).

| File | Contents |
|---|---|
| `parsed/regional_intratumoral_PK_anchors_literature.csv` | **27 rows**, every row PMID- and DOI-anchored. Columns: `source_pmid, source_doi, study, species, drug, subject_id, region, plasma_cmax, plasma_auc, tissue_cmax, tissue_auc, tissue_to_plasma_auc_pct, units_conc, units_auc, method, notes`. 21 human / 6 rat. |
| `parsed/oh2026_neurooncol_TableS3_per_sample_tumor_drug_ratios.csv` | **125 surgical samples across 42 human high-grade-glioma cases**. Columns: `case, mannitol_status, lev_collection_time_min, cef_collection_time_min, lev_tumor_to_plasma_ratio, cef_tumor_to_plasma_ratio, region_contrast_enhancing, region_flair, region_necrosis`. Probes: levetiracetam (BBB-penetrant) and cefazolin (BBB-excluded). |
| `parsed/jackson2018_regadenoson_TMZ_microdialysis_TableS1_raw.csv` | Jackson 2018 raw table |
| `epmc_PMC12962645_supp.zip` + `PMC12962645_supp/` | Oh 2026 *Neuro-Oncology* `noaf204` supplement |
| `epmc_PMC5767971_supp.zip` + `PMC5767971_supp/` | Jackson 2018 *Fluids Barriers CNS* `12987_2017_88` supplement |

**Underlying studies** (each row carries its own PMID): Portnow 2009 (PMID 19861433,
TMZ peritumoral microdialysis), Blakeley 2008 (18787762, methotrexate CE vs
non-enhancing), Jackson 2018 (29332604, TMZ CE vs non-CE), Grossman 2012/2013
(22644796 / 23649683, TMZ rat ECF), Lee 2025 ABTC (39820282, adavosertib), Eroglu
2024 (39422741, eribulin).

**Re-obtain** — **[HUNT]** `https://pmc.ncbi.nlm.nih.gov/articles/PMC12962645/`
and `https://pmc.ncbi.nlm.nih.gov/articles/PMC5767971/`;
**[DERIVED]** supplements via
`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12962645/supplementaryFiles`.
The two `parsed/*.csv` files were **hand-curated from those papers** — they are not
a download. **If lost, they must be re-extracted from the supplements above.**

## B2. `geomx_dsp/` — 42 MB — ★ **highest-value single directory**

**MANUAL + SCRIPTED.** The Moncayo HGSOC PARPi explant set — the only open asset
in existence carrying a per-pixel measured drug field *and* a molecular state on
the same tissue (`results/DATA_HUNT.md` T1.1).

| File | Contents |
|---|---|
| `NatCommun2026_SourceData_SupFig2_MSI_perpixel.csv.gz` (13 MB) | **376,736 MALDI-MSI pixels × 19 cols**, 20 µm pitch, X/Y coordinates, 3 drugs (olaparib 152,146 px / niraparib 149,237 px / rucaparib 75,353 px), raw + heavy-labelled internal standard + normalised, 18 slices, 3 patients |
| `zenodo17610220_ROI_drug_measurements.csv` | 114 ROIs × 13 cols |
| `zenodo17610220_ROI_drug_plus_expression.csv.gz` | the same 114 ROIs × 11,224 genes of GeoMx WTA |
| `zenodo17610220_drugMeasurement_MainTable.rds`, `zenodo17610220_modelling.Rmd`, `zenodo17610220_README.txt` | the authors' own `lme4` spec |
| `GSE281519_RAW.tar` + `GSE281519_raw/` | **132 `.dcc.gz`** ROI-level DSP counts |
| `GSE232050_*`, `GSE254145_*`, `GSE271255_*` | three further GeoMx accessions (GBM progression vs pseudoprogression; primary/recurrent GBM; GeoMx counts) |

**Re-obtain:**

```bash
# [HUNT] Zenodo deposit (CC-BY-4.0)
#   https://zenodo.org/records/17610220
# [HUNT] GEO
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281519
# [DERIVED] GEO raw tar
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281519/suppl/GSE281519_RAW.tar

# [HUNT] THE PER-PIXEL MAP — 130,354,714 B, 70 sheets. Springer ESM ONLY.
curl -L -o 41467_2026_70558_MOESM6_ESM.xlsx \
  'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-70558-1/MediaObjects/41467_2026_70558_MOESM6_ESM.xlsx'
```

> **The single most important file here is NOT re-downloadable from any
> accession.** `NatCommun2026_SourceData_SupFig2_MSI_perpixel.csv.gz` was
> extracted from sheet **`SupF. 2 (MSI)`** of that 130 MB workbook. The Zenodo and
> GEO deposits contain only code and ROI-level scalars. **The PMC `bin/` mirror
> serves a ~1.8 KB HTML stub — you must use the Springer endpoint.**

Paper: Moncayo et al., *Nat Commun* 2026, PMID **41844641**,
DOI `10.1038/s41467-026-70558-1`, <https://pmc.ncbi.nlm.nih.gov/articles/PMC13144616/>.
Open, CC-BY-4.0, no DUA. **Ovarian, not glioma; bath-dosed ex vivo (20 µM, 24 h),
so it validates diffusion + uptake + sequestration, never vascular delivery.**

## B3. `maldi_msi_drug/` — 20 MB — ★ **best in-vivo human transport asset**

**SCRIPTED — and the harvest scripts are on disk, which makes this the
best-documented hunt directory.**

| File | Contents |
|---|---|
| `harvest_gist_imatinib.py`, `harvest_gist_quantitative.py` | **the re-download recipe** |
| `gist_imatinib_msi_index.csv` | 55 rows × 20 cols, incl. a per-row CDN `image_url` |
| `gist_imatinib_quantitative_index.csv` | the quantitative dilution series index |
| `ion_images/` | 110 `.npy` (55 intensity + 55 `.mask.npy`) |
| `ion_images_quant/` | 13 `.npy` |

**Provenance [ON DISK]:** GraphQL endpoint `https://metaspace2020.eu/graphql`;
METASPACE project UUID **`2b505b2a-ca3a-11e9-9d77-9b72ac66ec3c`** ("GIST");
quantitative dataset ids `2022-12-22_11h50m30s`, `2022-07-11_17h11m07s`,
`2026-07-16_18h53m52s`; annotation filter `databaseId:22, fdrLevel:0.5`; compounds
`imatinib` (m/z 494.266285047, C29H31N7O) and `n-desmethylimatinib` (480.2506);
per-row CDN images of the form
`https://d3s59uoti91q90.cloudfront.net/iso/<dataset_id>/<hash>`.
Paper: Abu Sammour et al., *Sci Rep* 2019;9:10698,
<https://pmc.ncbi.nlm.nih.gov/articles/PMC6650609/>.

**Licence/access:** the derived 16-bit ion-image PNGs are anonymously fetchable;
**raw imzML/ibd is registration-gated** (free account at metaspace2020.eu).

Pairs with **`verified_0/`** (4.2 MB): `gist_S1_clinical.csv` (27×10),
`gist_spatial_plus_clinical.csv` (55×31 — the spatial↔clinical join, verified
**55/55 with zero orphans**), `sampleA_tumor_imatinib.npy` + mask, and the Sci Rep
supplement `41598_2019_47089_MOESM1_ESM.docx`.

## B4. `geo_organoid_screens/` — 72 MB — **MANUAL**

Read by `validation.gbm_tmz_exvivo` → `results/gbm_tmz_exvivo.{md,png}`. **GSE232173
is the only cached asset that is the right disease, the right drug and a measured
response at once**: patient-derived GSC cultures under temozolomide with ATP
viability (Ntafoulis et al., *Br J Cancer* 2023;129:1327-1338, **PMID 37620410**,
doi `10.1038/s41416-023-02402-y`; Erasmus MC / GLIOTRAIN). It scores the potency
axis only — no drug concentration was measured, and dissociated cultures have no
microenvironment.

| File | Accession |
|---|---|
| `GSE232173_count_data_24_cells.txt.gz`, `JOINED_GSE232173_TMZviability.csv` | **GSE232173** — glioma + TMZ + measured ex-vivo viability |
| `GSE148842_family.soft.gz`, `GSE148842_sample_metadata.csv`, `GSE148842/` (5-file subset) | **GSE148842** |
| `BJC2023_MOESM1–5` (.docx/.xlsx) | *British Journal of Cancer* 2023 supplement |
| `Golebiewska2020_ActaNeuropathol_MOESM1.pdf` | Golebiewska et al. 2020 supplement |

**[DERIVED]** `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232173`;
supplements under `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232173/suppl/`.

## B5. The remaining hunter angles

Survey material. Not read by any validation module today; catalogued so it can be
re-obtained if a future validation needs it.

| Directory | Size | What it is | Accession / source | Class |
|---|---:|---|---|---|
| `sc_drug_perturb/` | 774 MB | **GSE148842** — 43 acute GBM slice cultures, 10 patients, scRNA-seq under vehicle / etoposide / panobinostat / 4 single arms; plus **GSE226202**, **GSE239651** | **[DERIVED]** `ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148842/suppl/GSE148842_RAW.tar` | MANUAL |
| `hgcc_expression/` | 755 MB | HGCC full resource: `HGCC_DNA_methylation.txt` (638 MB), RNA expression, CNV, `HGCC_drug_response_AUC.txt` (262 drugs), clinical; plus **GSE72217/72218/72209** | **[HUNT]** `https://hgcc.se/`, `http://hgcc.se/portal/json/{bootstrap,expr,cna,mutation,drugsummary}.json` | MANUAL |
| `pdx_encyclopedia/` | 631 MB | **GSE132650** Mayo GBM PDX expression (539 MB); Vaubel 2020 Tables S2–S7 incl. **`TableS7_PDX_response.csv`** (37 PDX lines, median survival under Placebo/TMZ/RT/RT+TMZ/BEV, MGMT status, patient OS); Novartis **PDX Encyclopedia** Gao 2015 `pdxe_gao2015_MOESM10.xlsx` (122 MB) | **[HUNT]** `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132650` | MANUAL |
| `depmap_gdsc_ctrp/` | 505 MB | Potency prior at scale: PRISM 19Q4 secondary (264 MB) + 24Q2 extended primary, DepMap **24Q4** expression/model, **CTRPv2** PharmacoSet (**Zenodo 3905470**), NCI-60 / DTP **ALMANAC** combination scores | **[HUNT]** `https://depmap.org/portal/`, `https://zenodo.org/records/3905470`; **[DERIVED]** NCI-60 via `https://discover.nci.nih.gov/cellminer/` | MANUAL |
| `zenodo_figshare/` | 412 MB | 4 deposits: GBM anti-PD-1 GeoMx (**Zenodo 16839828**); pamiparib multimodal MSI (**Zenodo 13375209**, CC-BY-4.0); PARPi MSI 10-workbook set (**Zenodo 5158808**); PDAC avatar — **`urls.txt` only, data never downloaded** (Dryad DOI `10.5061/dryad.qnk98sfr2`, file endpoints `https://datadryad.org/api/v2/files/{3269529,3269535,3269536,3269537,3269707}/download`) | **[ON DISK]** for the Dryad URLs; `https://zenodo.org/api/records/<id>/files/<name>/content` for Zenodo | MANUAL |
| `slice_organotypic/` | 361 MB | **duplicate** `GSE148842_RAW.tar`; Mann 2023 *Cell Rep Med* OBSC dose-response (`Mann2023_OBSC_doseresponse_extracted.csv`, 572 values, 7 models × 11 drugs); Horowitz 2020 *npj Precision Oncology* microfluidic slice | supplement URLs **not recorded** | MANUAL |
| `synapse_open/` | 227 MB | HTAN CHOP pediatric HGG myeloid scRNA `.h5ad` (233 MB, **no accession recorded**); **GLASS** Varn 2022 *Cell* `mmc1–4.xlsx`; HTAN OHSU SMMART `mmc2/3.xlsx` | **[HUNT]** `https://www.synapse.org/glass`, `https://humantumoratlas.org/` | MANUAL |
| `biorxiv_recent/` | 224 MB | **S-BIAD2027** (DOI `10.6019/S-BIAD2027`, **CC BY 4.0**) — IMC of human TNBC, NACT responders vs non-responders, 35-plex; `quiche_nt_preprocessed.h5ad` (216 MB, **no accession recorded**); Moncayo PARPi bioRxiv supplement | **[ON DISK]** accession + licence; **[DERIVED]** `https://www.ebi.ac.uk/biostudies/files/S-BIAD2027/<file>` | MANUAL |
| `pdmr_hcmi/` | 217 MB | HCMI glioma models: 44 STAR gene-count TSVs + clinical/treatment tables; NCI **PDMR** patient/specimen/treatment-history-and-best-response exports; `pdmr_ftp_file_listing.txt` | **[ON DISK]** every row of `hcmi_glioma_paired_expression_manifest.tsv` carries a GDC `file_id` → `https://api.gdc.cancer.gov/data/{file_id}`. Catalog `https://hcmi-searchable-catalog.nci.nih.gov/`; PDMR `https://pdmr.cancer.gov` | SCRIPTED |
| `imc_codex_treated/` | 206 MB | GBM imaging mass cytometry `SpatialExperiment` `.rds` (185 MB, 34-plex, **real µm stage coordinates** — **no accession recorded**); **IDR `idr0076`** METABRIC breast IMC annotation; MIBI melanoma with ICB response; CODEX melanoma | **[ON DISK]** CODEX: `https://doi.org/10.5061/dryad.k0p2ngfcc`. **[DERIVED]** IDR: `https://ftp.ebi.ac.uk/pub/databases/IDR/idr0076-ali-metabric/` | MANUAL |
| `cgga_clinical/` | 177 MB | **CGGA** mRNAseq_693 / mRNAseq_325 / mRNA_array_301 / scRNA_6148 expression + clinical (release 20200506 / 20220104); **GSE72951** = the **BELOB** trial (bevacizumab/lomustine, recurrent GBM) | **[HUNT]** `http://www.cgga.org.cn` | MANUAL |
| `upenn_lumiere/` | 85 MB | **UPENN-GBM** CaPTk radiomics — 67 CSVs of {ap-rCBV, PH, PSR} and {AD, FA, RD, TR} × {ET, ED, NC}, 467 subjects; **LUMIERE** (CC0) 616 expert RANO timepoints / 91 patients | **[ON DISK]** `UPENN-GBM_DownloadManifest20221129.tcia` carries `downloadServerUrl=https://nbia.cancerimagingarchive.net/nbia-download/servlet/DownloadServlet` + an explicit SeriesInstanceUID list. **[HUNT]** LUMIERE `https://doi.org/10.6084/m9.figshare.c.5904905.v1` | SCRIPTED |
| `msi_repositories/` | 75 MB | The Zenodo link records themselves: `links_5158808.txt`, `links_5776125.txt`, `links_5795371.txt` — **filename \| bytes \| URL per line, including items deliberately not downloaded** (e.g. `Rawdata_UbezioIJBS.zip`, 426,319,381 B). Plus METASPACE annotation dumps. | **[ON DISK]** `https://zenodo.org/api/records/<id>/files/<name>/content`; landing pages `https://zenodo.org/records/{5158808,5776125,5795371}`, all **CC-BY-4.0** | MANUAL |
| `geo_organoid_screens/` | 72 MB | see B4 | | MANUAL |
| `htan_open_tier/` | 31 MB | HTAN **Data Release v7.0** entity table (20 MB) + metadata, participants clinical (incl. `TherapeuticAgents`, `BestOverallResponse`, `TreatmentDose`), atlas index, access-tier summary | **[ON DISK]** `syn_meta.json` → Synapse **`syn26156529`**. **[HUNT]** `https://humantumoratlas.org/`, `https://data.humantumoratlas.org` | MANUAL |
| `idr_bioimage/` | 26 MB | **S-BSST808** intravital + ultramicroscopy vessel data with a derived `vessel_geometry_summary.json` (`vessel_vol_frac`, `dist_to_vessel_um_p50/p90/p99`, 5 nm voxels); **S-BIAD2130**; TONIC nivolumab breast MIBI-TOF file list (**no accession recorded**) | **[DERIVED]** `https://www.ebi.ac.uk/biostudies/studies/S-BSST808`, `/S-BIAD2130` | MANUAL |
| `organ_on_chip/` | 26 MB | BBB-chip **Zenodo 12804228** (incl. `HPLC__Drug_crossing_tests.xlsx`); lung-fibrosis MPS **Zenodo 15475937**; RSC *Lab on a Chip* `d5lc01105a`; Yi 2019 GBM-on-chip source data | **[ON DISK]** `https://zenodo.org/records/12804228`, `/15475937`; `https://doi.org/10.1039/d5lc01105a` | MANUAL |
| `ebi_arrayexpress/` | 24 MB, **719 files** | **E-GEOD-84010 = the AVAglio trial** (bevacizumab, newly diagnosed GBM): 349 NanoString `.RCC` + 349 sample tables + `clinical_response.csv`. Plus E-GEOD-7696/36782/58400, E-MTAB-1852/9282 (BIOMEDE, **NCT02233049**), E-TABM-579; an EGA dump (all `access_type: controlled`) | **[DERIVED]** `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/<ACC>`, files at `https://www.ebi.ac.uk/biostudies/files/<ACC>/<file>` | MANUAL |
| `dce_mri_response/` | 13 MB | Clinical tables only (the imaging is Aspera/gated): **UCSD-PTGBM** v3, Burdenko-GBM-Progression, RHUH-GBM, plus LUMIERE/UPENN duplicates | **[HUNT]** `https://www.cancerimagingarchive.net/collection/{ucsd-ptgbm,upenn-gbm,rhuh-gbm,burdenko-gbm-progression}/` — clinical CSVs are separately **CC BY 4.0** | MANUAL |
| `geo_treated_spatial/` | 12 MB | **GSE280951** (AMG-232 / navtemadlin, MDM2i) normalized counts; **GSE287631** GeoMx WTA + protein | **[DERIVED]** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE280nnn/GSE280951/suppl/` | MANUAL |
| `pharmacoscopy_exvivo/` | 11 MB | **DARTT-1** ex-vivo pharmacoscopy trial (33 patients × ~166 drug/combination columns with clinical follow-up); *Nature Medicine* `s41591-024-03224` ESM `MOESM3–17` | **[DERIVED]** `https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-024-03224-y/MediaObjects/41591_2024_3224_MOESM3_ESM.xlsx` | MANUAL |
| `tcga_treatment_response/` | 4.1 MB | `tcga_gbm_lgg_drug_response_RECIST.csv` (340 rows: patient / drug / best response / start-end day / dose) + TCGA legacy biotab clinical for GBM and LGG | **[HUNT]** `https://portal.gdc.cancer.gov/`, `https://api.gdc.cancer.gov/cases` | SCRIPTED |
| `spatial_perturbation/` | 3.3 MB | **GSE249198** — *"Single-cell spatial pharmacobiology identifies conserved stromal barriers to therapeutic antibody delivery in human solid tumors"* (Stanford) + the *Nat Biotech* source-data `MOESM3–13.xlsx`. **A measured antibody-delivery dataset; not assigned a tier in DATA_HUNT.md — arguably under-ranked.** | **[DERIVED]** `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249198` | MANUAL |
| `clinical_trial_open/` | 2.9 MB | **Catalogues, not data**: Vivli glioblastoma DataCite DOIs, CSDR sponsor study index, NCI Project Data Sphere CNS metadata (all `"downloadable": false`), PDS access PDF | **[ON DISK]** e.g. DOI `10.25934/pr00011517.0`; requests via `https://vivli.org` | gated |
| `radiopharm_pet/` | 1.7 MB | **ACRIN-6684** ¹⁸F-FMISO hypoxia PET in GBM: 31 SDTM clinical CSVs + `acrin6684_fmiso_pet_vs_outcome.csv` (38 patients, `SUVmean_tumor … days_to_death`) | **[ON DISK]** dbGaP **`phs004225`**; GDC/CRDC GUIDs `dg.4DFC/<uuid>` → `https://nci-crdc.datacommons.io/ga4gh/drs/v1/objects/dg.4DFC/<uuid>`. **[HUNT]** `https://www.cancerimagingarchive.net/collection/acrin-fmiso-brain/` | SCRIPTED |

## B6. `verified_0` … `verified_9` — 69 MB total

> Re-measured 2026-09-03 (`du -sch data_cache/hunt/verified_*`). This heading said
> 371 MB.

Independent re-download caches from the 10 verification agents. Two of them
**falsified** hunter claims, and those falsifications are load-bearing.

| Dir | Size | Contents | Accession | Outcome |
|---|---:|---|---|---|
| `verified_0` | 4.2 MB | GIST clinical + spatial join, sample-A tumour arrays, Sci Rep supplement | METASPACE "GIST"; Sci Rep 2019;9:10698 | Join verified **55/55, zero orphans**. **FALSIFIED: the GIST study has no treatment-response endpoint of any kind** — S1 carries only organ, met/primary, sex, age, dose regimen, plasma imatinib ng/mL, KIT/PDGFRA mutations, mitotic count. |
| `verified_1` | 2.6 MB | IS-normalised imatinib maps + FTICR calibration ladder | METASPACE `2022-07-11_17h11m07s`, `2022-12-22_11h50m30s` | supports the quantitative arm |
| `verified_2` | 23 MB | `.mzid.gz` + `.mzML` for GBM22/GBM84 AZD1775 (adavosertib) TMT phospho-Y | **PRIDE PXD018782** | **FALSIFIED: no quantified phosphosite response** — all 80 files are raw spectra or peptide identifications, zero reporter-ion/abundance/ratio fields. |
| `verified_3` | 612 KB | METASPACE annotation dump + one ion image | dataset `2018-11-20_15h25m20s` | negative check |
| `verified_4` | 15 MB | `Efficacy study.xlsx`, `IALU157bis_2.xlsb`, PK workbooks | **Zenodo 5776125** (PEGPH20 + paclitaxel) | **Critical caveat: the MSI-imaged animals are NOT the efficacy-followed animals** — zero ID overlap in SKOV3-HAS3 and BxPC3, two coincidental matches in SKOV3. The drug-distribution → response link is arm-level, never per-animal. |
| `verified_5` | 18 MB | `NIR_MDR734_1.xlsb` | **Zenodo 5158808** | file-present check |
| `verified_6` | 6.3 MB | **Full METASPACE metadata verbatim**: TOF `2022-12-22_11h50m30s` (Liver, metastatic GIST, DHB matrix, **Pixel_Size 50×50**, grid 830×337, `isPublic: true`); FTICR `2022-07-11_17h11m07s` (small intestine, RP 160,859 @ m/z 200) | METASPACE | **confirms the 50 µm pixel figure** that every micrometre number in `gist_imatinib.md` inherits |
| `verified_7` | 180 KB | pioglitazone MSI ↔ HPLC bridge + `parsed_MSI_tidy.csv`, `derived_tumor_liver_partition.csv` | **Zenodo 5795371** | the only open MSI→absolute-µg/g calibration (MSI 12.9 / 28.5 µg/g vs HPLC 8.7 / 6.2 / 9.2 µg/g) |
| `verified_8` | 304 MB | 46 sections × {`.imzML`, `.ibd`, `.mask.nrrd`, `.RMS.nrrd`} + a 3D stack | **Zenodo 17374882** | 3D spheroid MALDI fixture for exercising a 3D solver — **no drug, no response** |
| `verified_9` | 108 KB | UCSD-PTGBM clinical v3 + negative cases + bval/bvec | TCIA **UCSD-PTGBM**, **CC BY 4.0, no DUA** | clinical only; the 33.1 GB NIfTI package is Aspera faspex |

---

# Known gaps in this record

Read this before relying on Tier B.

1. **No manifest exists for 30 of 38 hunt directories.** Re-download depends on
   accessions in filenames plus `results/DATA_HUNT.md`. **No checksum and no exact
   `curl` invocation is recorded for any hunt file.**
2. **Provenance permanently lost** — these cannot be re-obtained from anything
   written down anywhere in the repo:
   - `organ_on_chip/dryad_*.zip` — three failed/empty Dryad downloads; the README
     is literally `{"error":"Unauthorized, must have current bearer token"}` and no
     DOIs were recorded
   - `imc_codex_treated/gbm_imc_all_samples_spe.rds` (185 MB) and
     `mibi_melanoma_*.csv` — no accession
   - `biorxiv_recent/quiche_nt_preprocessed.h5ad` (216 MB) — no accession
   - `synapse_open/htan_chop_phgg_myeloid.h5ad` (233 MB) — no Synapse id
   - `idr_bioimage/TONIC_*` — no accession
3. **Broken artefacts — do not re-fetch as-is:**
   `msi_repositories/imatinib_annotation_TOF-GIST.json` (contains a GraphQL error,
   not data), `organ_on_chip/dryad_README.md` (auth-error blob),
   `organ_on_chip/dryad_data.zip` and `dryad_organoid_array.zip` (56 bytes each),
   `cgga_clinical/__MACOSX/` (junk).
4. **Duplication was the largest single waste when this was written, and it is gone.**
   Re-checked 2026-09-03: `GSE148842_RAW.tar` is **no longer anywhere under
   `data_cache/`**, so the ~716 MB figure below and the total it was a fraction of are
   both historical. The list is kept as the record of what to avoid re-downloading.
   `GSE148842_RAW.tar` (358 MB) existed
   twice (`sc_drug_perturb/`, `slice_organotypic/`); `GSE148842/*.cts.txt.gz` three
   times; `NIR_wt943_2.xlsb` three times; the Zenodo 5776125 and 5795371 sets twice;
   LUMIERE/UPENN clinical CSVs twice; `zenodo_figshare/gbm_pd1_geomx_z16839828/`
   stores both extracted and `.tar.gz` copies. A restore can safely skip the
   duplicates.
5. **The hand-curated files are not downloads.**
   `autoradiography_microdialysis/parsed/*.csv` and
   `geomx_dsp/NatCommun2026_SourceData_SupFig2_MSI_perpixel.csv.gz` were extracted
   from papers and a journal workbook. Losing them means redoing the extraction,
   not re-running a `curl`. **These are the two highest-value files in Tier B and
   the two least recoverable — back them up separately.**
6. **Two falsified claims must not be re-propagated** (see B6): the GIST/imatinib
   study has **no** treatment-response endpoint; PRIDE PXD018782 has **no**
   quantified phosphosite response.

# What cannot be downloaded at all

Nothing in the open tier pairs **spatial GBM tissue** with that patient's
**measured drug response**. The experiment that would close the loop — the
Agar/White adavosertib GBM MALDI-MSI platform, drug concentration per pixel in
phase-1 patient tissue with molecular response — has **no accession**; the data
statement says verbatim *"available upon request"*. That is a months-long
collaboration request to the MIT–Mayo PSOC, not a download. Ranked access paths
for the gated tier (HTAN WUSTL Visium Level 1/2 via dbGaP; TCIA restricted
collections; Vivli/CSDR clinical-trial IPD; METASPACE raw imzML behind free
registration) are tabulated in
[`../results/DATA_HUNT.md`](../results/DATA_HUNT.md) §4, together with the
confirmed dead ends so nobody re-spends the time.

# Deleted on purpose — 2026-08-30 disk sweep

The repository reached 8.4 GB on a working laptop and was swept. Two rules
governed it: nothing a TEST depends on was removed, and nothing was removed
without a recorded way back.

**`data_cache/greenwald_gbm_codex/ZH1007_INF.ome.tif` — 1,548,785,153 bytes,
DELETED.** The second of the two CODEX pairs. Verified unreferenced by any
code, test, experiment or CI workflow before removal — the only mentions in the
repository are in `results/gbm_3d_registration.md`, which is a dated record of a
pass that, in its own words, added no test because no new code existed to pin.
Its derived output, `results/gbm_3d_registration_zh1007inf_overlay.png`, is
still present. **To get it back:** Zenodo `10.5281/zenodo.12624860`, served
under exactly that filename, exactly that byte count.

**`data_cache/greenwald_gbm_codex/ZH1007_NEC_B.ome.tif` — 920,859,118 bytes,
KEPT.** `tests/test_codex_adapter.py` skips its whole module without it, so
deleting it would have converted eleven live tests into eleven silent skips —
a green suite testing nothing, which is the defect class this project is most
prone to. Confirmed still running after the sweep: 11 passed.

**Tarballs in `data_cache/*` — KEPT**, though their extracted directories sit
beside them. They look redundant and are not: `tmesim/data/ingest.py` looks for
`*cell_feature_matrix.tar.gz` by name, and `tmesim/data/real/singlecell_spatial.py`
records each one's exact expected byte count. Deleting them triggers a
re-download, not a saving.

**`data_cache/hunt/` — KEPT.** It reads like scratch from a finished search and
is not: `tmesim/validation/gbm_tmz_exvivo.py`, `scorecard.py`,
`moncayo_explant.py` and `kill_calibration.py` all read directories under it by
path.

**`app/runs/*/timelapse.bin` — 54 runs pruned, 782,579,352 bytes**, by
`scripts/prune_run_payloads.py`. Runs cited in `results/`, `archive/`, `docs/`
or `README.md` are protected, as are the ten most recent. A pruned run records
the fact in its own `result.json`, and the API answers **410** naming the date
and the byte count rather than the pre-existing 404 reason, which would have
claimed the run never captured one. Re-running the spec in the record
reproduces the payload.

**`.export-drive/` (72 MB), `landing-3d/.next`, `__pycache__`,
`.pytest_cache`, `tmesim.egg-info`, two zero-byte `curl_codex*.log` — DELETED.**
All regenerated on demand; `.export-drive/` is rebuilt by
`scripts/studio_export/drive.mjs`.

Result: 8.4 GB → 6.2 GB.

## One casualty of the 2026-08-30 prune, recorded rather than quietly fixed

`app/runs/r_00028/timelapse.bin` **was deleted despite being cited.**
`docs/RUNBOOK_DEPLOY.md:368` cites three runs as a range —
``r_00027`-`r_00029`` — and the prune script's citation scan matched run ids
one at a time. It found the two endpoints, protected them, and never saw the
middle: `r_00028` is named nowhere on that line. The scan found *something*
there, so it looked like it had worked.

Found by the claim-surface audit of the same day, not by the script. The script
now expands ranges (ASCII hyphen, `--`, en dash, em dash and the word "to"),
`r_00028` is protected going forward, and
`tests/test_app_server.py::TestPruneCitationProtection` pins it against a
temporary tree — the first version of that test passed with the bug still in
place, because the audit report itself names `r_00028` in prose and `results/`
is a citation root.

The file is not recoverable as itself: re-running the RunSpec produces a new
run id, not `r_00028`. The run record is intact and the API answers `410` with
the date, the byte count and how to regenerate, so a reader following that
citation gets a straight answer rather than a missing file. That is the whole
reason the prune records itself; it is not a substitute for not having deleted
it.
