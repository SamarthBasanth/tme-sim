# Ivy GAP — anatomic-structure RNA-seq (real GBM regional architecture)

**Ivy Glioblastoma Atlas Project (Ivy GAP)**, Allen Institute for Brain Science.
The real regional architecture of human glioblastoma: tissue laser-microdissected
into named **anatomic structures** and each profiled by RNA-seq. This is our
source for **real niche/region composition** — margin, core, perivascular, and
necrotic zones — so the reconstructed microenvironment's regions are labelled
against genuine measured transcriptional signatures (not just the tumour bulk).

- Citation: Puchalski et al. 2018, *Science* 360(6389):660-663,
  "An anatomic transcriptional atlas of human glioblastoma."
- Portal: http://glioblastoma.alleninstitute.org (Downloads tab:
  https://glioblastoma.alleninstitute.org/static/download.html)
- GEO mirror (RNA-seq): GSE107559.

## Open status — FULLY OPEN

All files below are freely downloadable from the Allen Institute portal with **no
registration and no controlled-access tier** (Allen Institute Terms of Use;
freely redistributable for research). No dbGaP / HTAN gate. Verified live and
re-downloaded 2026-07-26.

## Files in `data_cache/ivygap/`

> ### ⚠ STATUS ON THIS MACHINE, MEASURED 2026-09-03: none of these files is on disk
>
> `data_cache/ivygap/` **does not exist**. What exists is
> `data_cache/ivy_gap/gene_expression_matrix.zip` (~37 MB), **unextracted**. The
> loader falls back from `ivygap/` to `ivy_gap/` for the *directory*
> (`tmesim/data/real/ivygap.py:61,64`) but does **not** extract the archive, so
> `load_structures()` raises today:
>
> ```
> FileNotFoundError: Missing cached Ivy GAP file:
>   …/data_cache/ivygap/fpkm_table.csv. Extract gene_expression_matrix.zip
>   into data_cache/ivygap/ first (see docs/real_data/ivygap.md).
> ```
>
> **Four of the eight files below are recoverable by extracting that zip; four are
> not in it at all.** The archive holds exactly `fpkm_table.csv`,
> `rows-genes.csv`, `columns-samples.csv` and `README.txt`. The four
> `/api/v2/gbm/<name>.csv` summary tables — `rna_seq_samples_details.csv`,
> `sub_block_details.csv`, `gene_expression_details.csv` and `tumor_details.csv` —
> are neither on disk nor in the zip, and must be re-fetched from the portal.
>
> **What that means for the tables further down this page.** The marker-signature
> table is derived from the FPKM matrix and *is* reproducible once the zip is
> extracted. The **fractional-area** table and the **subtype / MGMT** counts are
> derived from `sub_block_details.csv` and `tumor_details.csv` respectively, and are
> **not obtainable from anything this repository currently holds** — they are
> reported here on the authority of the 2026-07-26 run, not on anything a reader can
> re-derive today. Until 2026-09-03 this page presented all of it as cached and
> present.

Extracted from the open `gene_expression_matrix.zip`
(`/api/v2/well_known_file_download/305873915`, ~37 MB zipped) plus four open
summary CSVs (`/api/v2/gbm/<name>.csv`).

| File | Size | Contents |
|---|---|---|
| `fpkm_table.csv` | 82 MB | genes x samples normalised FPKM matrix. Row index = `gene_id`, header row = `rna_well_id`. 25,873 genes x 270 samples. |
| `rows-genes.csv` | 1.6 MB | `gene_id` -> `gene_symbol`, `gene_entrez_id`, `chromosome`, `gene_name`. |
| `columns-samples.csv` | 52 KB | `rna_well_id` -> `structure_id`/`structure_abbreviation`/`structure_name`/`structure_color`, plus `tumor_id`/`block_id`/`specimen_id`. 270 rows. |
| `rna_seq_samples_details.csv` | 100 KB | Per-sample structure ontology (`structure_acronym`, `structure_name`, `structure_parent_acronym`) joined to donor clinical fields. |
| `sub_block_details.csv` | 232 KB | Tumour-feature (H&E histology) **fractional areas** per tumour sub-block — the 11-feature region composition. 946 sub-blocks x 42 donors. |
| `gene_expression_details.csv` | 4.8 MB | ISH expression "energy" per anatomic feature per sub-block (image-based, separate from RNA-seq FPKM). |
| `tumor_details.csv` | 4 KB | Per-donor `molecular_subtype`, `mgmt_methylation`, `survival_days`, `egfr_amplification`, `age_in_years`, KPS. 42 donors. |
| `README.txt` | 2 KB | Original Allen Institute matrix README. |

Note: acquisition (`tmesim/data/acquire.py`) also caches the raw zip under
`data_cache/ivy_gap/`; the loader accepts either spelling.

## Format / how the tables join

- FPKM matrix: `gene_id` (int) x `rna_well_id` (int). The loader maps genes to
  symbol via `rows-genes.csv` and keeps sample columns as string `rna_well_id`.
- Samples join `columns-samples.csv` (`rna_well_id`) <-> FPKM columns 1:1 — all
  270 samples present on both sides, zero mismatch.
- Clinical/molecular fields join on `rna_seq_samples_details.csv.sample_id` ==
  `rna_well_id` (122/122 reference-histology samples annotated).

## Anatomic-structure vocabulary (RNA-seq) + sample counts

270 RNA samples, **37 tumours**, two study designs sharing one FPKM matrix
(README quotes 122 + 148 = 270; the 148 span 34 tumours by ISH):

- **Reference histology (H&E), 122 samples / 10 tumours** — the 5 canonical
  anatomic structures. This is the clean, well-powered set.
- **CSC-cluster (ISH), 148 samples** — putative cancer-stem-cell / niche clusters
  sampled by HIGH (`-<GENE>`, 89 samples) or LOW (`-control-...`, 59 samples)
  expression of an 18-gene reference set. Adds two vascular/necrotic niches.

By **base anatomic structure** (all 270 samples), with region mapping:

| Base | Structure | Region | n (all) | n (H&E ref) |
|---|---|---|---|---|
| LE | Leading Edge | margin | 19 | 19 |
| IT | Infiltrating Tumor | margin | 24 | 24 |
| CT | Cellular Tumor | core | 111 | 30 |
| CTmvp | Microvascular Proliferation | perivascular | 28 | 25 |
| CThbv | Hyperplastic Blood Vessels (in CT) | perivascular | 22 | 0 (CSC only) |
| CTpan | Pseudopalisading Cells Around Necrosis | necrotic | 40 | 24 |
| CTpnz | Perinecrotic Zone | necrotic | 26 | 0 (CSC only) |

Coarse **microenvironment-region** counts (all 270): core 111, necrotic 66,
perivascular 50, margin 43. (H&E-ref only: core 30, margin 43, perivascular 25,
necrotic 24.) `CThbv`/`CTpnz` exist only among CSC-cluster samples, so set
`reference_histology_only=False` to populate the perivascular/necrotic groups
when grouping by base/region.

The full fine-grained vocabulary is 80 `structure_abbreviation` labels (e.g.
`CTpnz-PI3`, `CThbv-TGFBR2`, `CT-control-...`); the base token before the first
`-` gives the 7 structures above.

## Tumour-feature (H&E histology) annotations

`sub_block_details.csv` gives machine-classified **fractional area** of each of
11 tumour features per sub-block (columns `normalized_area_*`), i.e. the
whole-slide region composition:

`le` Leading Edge · `lehbv` LE hyperplastic blood vessels · `it` Infiltrating
Tumor · `ithbv` IT hyperplastic blood vessels · `ct` Cellular Tumor ·
`ctpnz` perinecrotic zone · `ctpnn` perinecrotic (narrow) · `ctpan`
pseudopalisading cells around necrosis · `cthbv` CT hyperplastic blood vessels ·
`ctmvp` microvascular proliferation · `ctne` necrosis.

Mean fractional areas (946 sub-blocks): ct 0.46, ctne 0.29, le 0.29, it 0.27,
ctpnz 0.07, ctpan 0.06, cthbv 0.05, ctmvp 0.04 — a realistic whole-tumour region
mix usable as a prior for region abundances in reconstruction.

## Per-donor molecular context (`tumor_details.csv`, 42 donors)

Molecular subtype: Proneural 8, Classical 7, Mesenchymal 5, mixed (Classical/
Mesenchymal, Classical/Neural, Mesenchymal/Neural, Neural/Proneural) 13, Neural 4,
NA 5. MGMT methylation: No 23 / Yes 17 / NA 2. Also survival_days, EGFR
amplification, KPS, age.

## How per-structure signatures label the microenvironment

`structure_signatures()` averages FPKM over the samples of each structure/region
to give a **real measured transcriptional signature per zone**. Marker behaviour
(mean FPKM) is consistent with the region mapping being real rather than assumed —
each marker is most abundant in the zone its biology predicts:

| Marker (biology) | margin | core | perivascular | necrotic |
|---|---|---|---|---|
| PECAM1 / CD31 (endothelium) | 6.8 | 7.6 | **90.5** | 3.3 |
| CA9 (hypoxia) | 0.3 | 4.5 | 1.0 | **29.6** |
| VEGFA (angiogenesis / hypoxia) | 6.8 | 38.9 | 6.3 | **268.7** |

(region-grouped, all samples). PECAM1 dominates the perivascular zone
(CTmvp/CThbv), CA9 and VEGFA dominate the necrotic zone (CTpan/CTpnz) — exactly
the expected GBM biology. So reconstruction can score any reconstructed region
against these signatures to assign it a real label: **margin** (LE/IT — the
infiltrative brain-tumour interface), **core** (CT — dense proliferative bulk),
**perivascular** (CTmvp/CThbv — angiogenic/vascular niche), **necrotic**
(CTpan/CTpnz — hypoxic perinecrotic/pseudopalisading border).

## Loader — `tmesim/data/real/ivygap.py`

```python
from tmesim.data.real import load_structures, structure_signatures

d = load_structures()                      # 270 samples; .samples, .expression, .genes
d = load_structures(reference_histology_only=True)   # clean 122-sample H&E set
d.samples[["structure_abbreviation", "base_structure", "region", "molecular_subtype"]]
d.expression        # genes x samples FPKM (index = gene symbol)

sig = structure_signatures(by="base_structure", reference_histology_only=False)  # genes x 7
sig = structure_signatures(by="region", reference_histology_only=False)          # genes x 4
```

- `load_structures(cache_dir=None, reference_histology_only=False) -> IvyGapStructures`
  with `.samples` (sample x structure + region + clinical), `.expression`
  (genes x samples FPKM), `.genes`, and `.region_for(sample_id)`.
- `structure_signatures(by="base_structure"|"region"|"structure_abbreviation")`
  -> genes x groups mean FPKM.
- `REGION_MAP` (base structure -> margin/core/perivascular/necrotic) and
  `TUMOR_FEATURES` (11-feature histology vocabulary) are exported constants.
- Missing cache raises `FileNotFoundError` with acquisition guidance — no
  synthetic fallback.
