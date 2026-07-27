# DATA_HUNT — Open data that could validate a patient-specific TME drug-transport/response simulator

Date: 2026-07-27
Scope: 28 search angles; 3 candidates put through an independent adversarial verification pass; remaining
angles' download caches inspected directly on disk for this synthesis.

**Provenance labels used below**

| Label | Meaning |
|---|---|
| `[VERIFIED]` | Put through the independent verification pass, which re-downloaded from source and re-parsed. Claims here are the *post-correction* claims. |
| `[INSPECTED]` | Not in the verification pass. I opened the local files myself for this document and report only what I directly observed (shapes, column names, value ranges), plus HTTP-200-checked URLs. |
| `[REPORTED]` | Hunter-reported only; not independently re-checked. Treat as a lead. |

Every URL in this document was `curl`-checked and returned HTTP 200 (or 202 for the figshare DOI) on
2026-07-27.

---

## 1. Headline

**Yes — but exactly one, and it is not in glioma.**

One open dataset genuinely pairs a measured spatial drug distribution with spatially matched molecular
tissue state and a co-located molecular response signal, in human tumour tissue:

> **Moncayo et al., *Nat Commun* 2026 — PARP inhibitors in human high-grade serous ovarian carcinoma
> patient-derived explants.** MALDI-MSI drug maps (section 0) → GeoMx whole-transcriptome ROIs placed on
> high-drug and low-drug regions (section +1) → IHC for cleaved caspase-3 / γH2AX on flanking sections.
> Zenodo 17610220 (CC-BY-4.0) + GEO GSE281519 + the journal's open Source Data file.
> PubMed PMID 41844641, DOI [10.1038/s41467-026-70558-1](https://doi.org/10.1038/s41467-026-70558-1).

That is the good news, and it is real. Now the four things that must be said in the same breath:

1. **It is ovarian, not glioma.** Every other asset in `tme-sim` is glioma.
2. **The explants were dosed ex vivo in a bath** (20 µM, 24 h, per the GEO treatment protocol). There is
   no vasculature and no perfusion source term. It validates diffusion + cellular uptake + sequestration.
   It does **not** validate vascular delivery.
3. **The pixel-level drug map is real and open, but it is in the journal's Source Data, not in any
   repository.** The Zenodo/GEO deposits carry only ROI-level scalars (114 ROIs, 3 patients, 12 slices).
   The **376,736-pixel, X/Y-coordinate, 20 µm MALDI-MSI table** is sheet `SupF. 2 (MSI)` of the 130 MB
   Source Data workbook — which is why every repository-oriented search missed it. n is still 3 patients.
4. **The "response" in the open files is a transcriptional apoptosis signature**, not viability, not
   RECIST, not survival.

And the two harder headline findings, which matter more than the win:

- **No open dataset anywhere pairs a measured spatial drug distribution with a measured per-specimen
  treatment response in the same tumour.** Not one. Every candidate breaks at that link.
- **Nothing exists for glioma + temozolomide.** METASPACE (≈19–50k datasets) contains **zero** datasets
  matching temozolomide `[REPORTED, consistent with two independent angles]`. Only three papers in all of
  PubMed cover drug-distribution MSI in glioma, and none deposits its MSI data. The one study that did the
  exact experiment in GBM (Agar/White adavosertib platform) deposited only raw phosphoproteomics spectra —
  see §4.

Two claims from the raw hunt were **falsified** by verification and must not be re-propagated:

- The GIST/imatinib study has **no treatment-response endpoint of any kind** — no RECIST, no Choi, no PFS,
  no OS. `[VERIFIED]` Its Supplementary Table S1 carries organ, met/primary, sex, age, dose regimen,
  plasma imatinib ng/mL, KIT/PDGFRA mutations, and mitotic count. Nothing else.
- PRIDE PXD018782 contains **no quantified phosphosite response** — all 80 files are raw spectra or
  peptide identifications, with zero reporter-ion/abundance/ratio fields. `[VERIFIED]`

---

## 2. Tier 1 — can close the loop

Tier 1 splits cleanly into two halves that, with one exception, **do not overlap**: datasets with a
measured spatial *drug* field, and datasets with spatial tissue + a measured *response*. Only T1.1 and
T1.2 have anything on both sides.

### 2A — Measured spatial drug distribution

#### T1.1 — Moncayo et al. 2026, PARP inhibitors in HGSOC patient-derived explants ★ best overall
- **URL**: https://zenodo.org/records/17610220 (CC-BY-4.0) ·
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281519 ·
  paper https://pmc.ncbi.nlm.nih.gov/articles/PMC13144616/ ·
  Source Data (130 MB, 70 sheets) `41467_2026_70558_MOESM6_ESM.xlsx` from the Springer ESM endpoint
- **Open status**: open_direct_download, no registration, no DUA. Zenodo and GEO are direct; the Source
  Data workbook downloads cleanly (HTTP 200, 130,354,714 bytes) from the Springer ESM endpoint
  `https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-70558-1/MediaObjects/41467_2026_70558_MOESM6_ESM.xlsx`
  — note the PMC `bin/` mirror returns a ~1.8 KB HTML stub, so use Springer. `[INSPECTED]`
- **What it contains** (all directly observed):
  - **Per-pixel MALDI-MSI drug map, with coordinates.** Source Data sheet `SupF. 2 (MSI)` (identical
    content in `Fig. 1c`) is **376,736 rows × 19 columns**:
    `file | Patient | Drug | Tissue | Slice | ROI | Nir[M+Na] | Ola[M+Na] | Ruc[M+H] |
    Heavy labelled_Nir[M+Na] | Heavy labelled_Ola[M+Na] | Heavy labelled_Ruc[M+H] |
    norm_Nir[M+Na] | norm_Ola[M+Na] | norm_Ruc[M+H] | Scan | X | Y | ID`.
    X spans 0.18–10.90 mm, Y spans 0.04–7.66 mm, and the coordinate step is **0.02 mm — a 20 µm pixel
    pitch**. 18 slices; 2,531 to 44,512 on-tissue pixels each (median 19,197). Example slice
    `P23-015_Nir4` is a 157 × 203 grid with 18,556 tissue pixels. All three drugs are present
    (Olaparib 152,146 px, Niraparib 149,237 px, Rucaparib 75,353 px) across L.Ovary / Omentum / Bladder.
    Both the raw drug ion and its **heavy-labelled internal standard** are given per pixel, so the
    IS-normalisation is auditable rather than taken on trust.
    *Independent specificity check I ran:* on niraparib-dosed slices mean `norm_Nir` = 0.421 while the
    off-target channels sit at 1.7e-4 (Ola) and 6.8e-5 (Ruc) — a 2,500–6,000× channel separation. Same
    pattern for olaparib (0.106 vs ~1e-5) and rucaparib (0.0147 vs ~1e-6/5e-5). The signal is real and
    channel-specific, not matrix background.
    Extracted for convenience to
    `data_cache/hunt/geomx_dsp/NatCommun2026_SourceData_SupFig2_MSI_perpixel.csv.gz` (13 MB).
  - `zenodo17610220_ROI_drug_measurements.csv`: **114 ROIs × 13 columns**. Per-ROI `drug_measure`
    (range 0.0010–1.961), `drug_content` High/Low (57/57), `measure_norm_factor`, `drug_measure_norm`,
    `Patient` (3: P23_015, P23_037, P23_071), `Slice` (12), `Tissue` (L.Ovary 73, Omentum 28, Bladder 13),
    `Compound` (Niraparib 62, Rucaparib 52).
  - `zenodo17610220_ROI_drug_plus_expression.csv.gz`: the **same 114 ROIs × 11,224 genes** of GeoMx WTA
    expression, i.e. drug and transcriptome on one row.
  - `drug_measure` is confirmed to be the MALDI-MSI drug ion intensity **relative to a drug-matched
    deuterated internal standard**, averaged over the ROI footprint: the Source Data sheet `Fig. 2b`
    reproduces the identical values (0.728396, 0.925717, …) under the header
    `PDE ID | ROI | drug intnsity | drug_content | Patient | Q3 | drug | q3norm`.
  - Method chain per the paper: section 0 = MALDI-MSI (niraparib [M+Na]⁺ 343.1535, olaparib [M+Na]⁺
    457.1652, rucaparib [M+H]⁺ 324.1512), then H&E; section +1 = GeoMx; sections −2/−1/+2/+3 = IHC for
    cleaved caspase-3 and γH2AX plus WT1/Pax8.
  - The deposit ships the authors' own `lme4` model:
    `gene ~ drug_measure_norm + (1|Patient) + (1|Patient:Slice)`.
  - Source Data also carries FACS rucaparib-fluorescence-sorted high vs low cells (`Fig. 3f`) and a
    21-day reculture outgrowth assay (`Fig. 3g`) — a genuine drug-level → clonogenic-response measurement,
    but in cell lines (OVCAR8/PEO4/PEO1/ES2), not in the tissue.
- **Validation it enables**:
  1. Score predicted per-ROI intracellular drug accumulation against measured per-ROI drug intensity,
     inside the authors' pre-specified mixed-effects frame (patient and slice as random effects).
     n = 114 observations — the largest per-observation count of any Tier 1 asset.
  2. A hard, falsifiable, compound-specific test: the weak bases **niraparib and rucaparib** show
     lysosomal-sequestration-driven accumulation heterogeneity; **olaparib does not**. A simulator with
     pH-partitioning / lysosomal trapping should reproduce that split; one without it structurally cannot.
     (Note: olaparib ROIs were dropped from the GeoMx arm because its drug hotspots were too small to hit
     the nuclei minimum — so the olaparib contrast lives in the MSI/Source Data, not the ROI table.)
  3. Drug-level → molecular-response coupling: regress your predicted local exposure against the measured
     apoptotic and lysosomal transcriptional signatures on the same ROIs.
  4. **Pixel-level penetration morphology.** With 20 µm X/Y coordinates on 18 slices you can score
     predicted vs measured intratumoral gradients, hotspot size distributions, and coefficient of
     variation per slice — the same class of test as the GIST maps in T1.3, but with an internal standard
     per pixel and three chemically distinct compounds in the same tissue.
- **What it does NOT enable**: absolute concentrations (values are IS-normalised ratios, unitless);
  joining the pixel map to the transcriptome at pixel resolution (the ROI table has no X/Y — only
  `Slice`, `ROILabel` and High/Low, so the drug↔expression link is at ROI granularity and across adjacent
  sections 0 and +1, not co-registered pixels); anything about vascular delivery; anything about patient
  outcome.

#### T1.2 — PEGPH20 + paclitaxel, Istituto Mario Negri
- **URL**: https://zenodo.org/records/5776125 (CC-BY-4.0), open_direct_download. `[REPORTED]`, partially
  `[VERIFIED]` (verification cache at `data_cache/hunt/verified_4/` holds `Efficacy study.xlsx`,
  `IALU157bis_2.xlsb`, `PK BxPC3.xls`, `PK_HAS3.xls`, `PK_SKOV3.xls`).
- **Contains**: 26 per-tumour MALDI-MSI workbooks (`.xlsb`, not imzML — which is why keyword searches for
  imzML miss it entirely) with a full 2D per-pixel paclitaxel/internal-standard grid at 100 µm plus
  tumour/stroma/ROI masks; bulk HPLC paclitaxel in plasma/tumour/liver (µg/g, µg/mL) at 1/3/4/24 h;
  per-mouse tumour growth curves for 112 mice across 3 xenograft models × 4 arms
  (CTRL / PEGPH20 / PTX / PEGPH20+PTX).
- **Why it is nearly the best thing here**: it is the only open dataset that *manipulates an explicit ECM
  transport barrier* — hyaluronan, varied genetically (SKOV3 vs SKOV3-HAS3) and enzymatically (PEGPH20) —
  and then measures penetration, bulk PK, and tumour shrinkage.
- **The caveat that must never be dropped**: the MSI-imaged animals are **not** the efficacy-followed
  animals. Cross-referencing animal IDs gives zero overlap in SKOV3-HAS3 and BxPC3 and only two
  coincidental matches in SKOV3. The drug-distribution → response link is **arm-level and model-level,
  not per-animal**. "Predicted penetration vs observed shrinkage per tumour" would be an overstatement.
- **Validation it enables**: does the simulator get the *direction and magnitude* of the arm-level
  penetration shift right when you degrade the HA barrier, and does the predicted delivered dose rank the
  four efficacy arms correctly, in all three models?

#### T1.3 — METASPACE "GIST": imatinib MSI in human GIST, tumour + matched normal ★ best in-vivo human
- **URL**: https://metaspace2020.eu/project/2b505b2a-ca3a-11e9-9d77-9b72ac66ec3c ·
  paper https://pmc.ncbi.nlm.nih.gov/articles/PMC6650609/ (Sci Rep 2019;9:10698)
- **Open status**: derived 16-bit ion-image PNGs are anonymously fetchable, no auth. Raw imzML/ibd is
  **registration-gated** ("Please Sign in to download"). `[VERIFIED]`
- **Contains** (post-verification, corrected): 46 public human datasets. **11 patients** have a tumour
  imatinib map; **6 patients (A, AA, J, K, L, U) have matched tumour + normal imatinib maps**. All ~20 TOF
  datasets have zero drug annotations — the drug maps are FTICR. Imatinib (m/z 494.2663) and
  N-desmethylimatinib (m/z 480.2506), per-pixel float arrays with a bitwise-consistent tissue mask.
  Pixel size is **not** in the deposit (`MS_Analysis.Pixel_Size` is null); the 50 µm figure comes from the
  paper text only. Values are **arbitrary ion counts**, not ng/g.
- **The pairing is real and was proved, not assumed**: the sample letter in each dataset name
  (`GIST_sample<X>_tumor/normal_FTICR`) is the same identifier as `Organ(ID)` in Supplementary Table S1.
  The join matched **55/55 spatial rows with zero orphans**.
- **Validation it enables** — this is a *transport* loop, not a response loop:
  - Per-lesion tumour:normal partition coefficient, with a within-patient control (6 patients).
  - Intratumoral gradient / coefficient of variation within each section.
  - Boundary condition is given: per-patient plasma imatinib ng/mL + mg/day dose. And there is real
    unexplained variance for a model to earn — **tumour-section mean imatinib vs plasma concentration is
    Spearman r = 0.03, p = 0.93 (n = 9)**. Systemic exposure predicts essentially nothing about tissue
    exposure. That gap *is* the thing a microenvironment transport model exists to fill.
  - Independent reproduction of the paper's headline from raw pixels: tumour:normal imatinib ratios
    0.02–0.14 for liver metastases vs 0.72 for the one primary small-intestine GIST.
- **Local**: `data_cache/hunt/maldi_msi_drug/` (55 index rows, 110 `.npy`, ~20 MB) and
  `data_cache/hunt/verified_0/` (`gist_S1_clinical.csv` 27×10, `gist_spatial_plus_clinical.csv` 55×31).
- **Hard limit**: n = 6 matched pairs, 9 tumour sections with a plasma value. Enough to falsify a badly
  wrong transport model; **not** enough to rank two plausible ones.

#### T1.4 — PARP-inhibitor MSI with an efflux-transporter contrast
- **URL**: https://zenodo.org/records/5158808 (CC-BY-4.0), open_direct_download.
- **Contains**: same Mario Negri `.xlsb` format — 10 per-tumour workbooks, niraparib in wild-type vs
  **MDR734 / MDR740 P-glycoprotein-overexpressing** xenografts, plus an olaparib series, plus a 426 MB raw
  archive. Verified parse of `NIR_wt943_2.xlsb`: `distr` grid 90×95, 8,512 finite values, drug/IS ratio
  0.000–3.141, mean 0.276, with `mask`/`maskT`/`maskS`/`ROI` sheets. `[REPORTED]`; verification cache holds
  `verified_5/NIR_MDR734_1.xlsb` (18 MB). `[INSPECTED — file present]`
- **Validation it enables**: a clean test of whether the transport model handles **active efflux**,
  isolated from the ECM axis in T1.2, plus a two-drug lipophilicity/permeability sweep. **No response data
  in the deposit.**

#### T1.5 — Pamiparib multimodal drug distribution
- **URL**: https://zenodo.org/records/13375209 (CC-BY-4.0). `[INSPECTED — local files listed]`
- **Contains**: `annotated_positive.xlsx`, `annotated_negative.xlsx`, `LC_MSMS_data.xls`, and two
  whole-slide TIFFs (`retta_pam.tif`, `retta_ola.tif`). Pamiparib is explicitly framed as *not* an
  ABC-transporter substrate, so it is the natural negative control for the T1.4 efflux test.
- Local: `data_cache/hunt/zenodo_figshare/pamiparib_msi_z13375209/`. Not yet parsed in depth.

### 2B — Spatial tissue + measured response (but **no drug measurement anywhere in them**)

#### T1.6 — UCSD-PTGBM
- **URL**: https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/ — CC BY 4.0, **no DUA**.
- **Contains**: 178 patients / 243 post-treatment timepoints. Voxelwise leakage-corrected DSC-CBV, ASL
  perfusion, RSI multishell-diffusion cellularity maps, T1/T1ce/FLAIR, all in 256³ 1 mm space; BraTS
  4-compartment plus neuroradiologist-corrected "cellular tumour" masks. Per-timepoint hard outcome labels:
  192 residual/recurrent (45 pathology-proven) vs 51 post-treatment-change (4 pathology-proven
  pseudoprogression, 9 clinically-confirmed PsP, 10 radiation necrosis, 28 non-specific). OS 243/243, PFS
  242/243, chemo agent (Temodar n=139, CCNU n=3), bevacizumab flag and dates. IDH/MGMT/1p19q/ATRX.
- **Validation it enables**: whether predicted low-delivery voxels coincide with residual/recurrent voxels
  and predicted high-delivery voxels with treatment-effect voxels, in the right disease, at voxel level.
- **The disqualifying caveat, stated plainly**: CBV is a **gadolinium** measurement. Gd (~550 Da,
  hydrophilic, BBB-excluded, bolus kinetics) is not temozolomide (194 Da, lipophilic, BBB-penetrant, oral).
  Using CBV to validate a drug-penetration map requires an unvalidated leap. Say so, or the claim is
  dishonest.
- **Local**: clinical tables only — `data_cache/hunt/dce_mri_response/` and `verified_9/`. The 33.1 GB
  NIfTI package is Aspera faspex (needs Connect/`ascli`, no subsetting).

#### T1.7 — LUMIERE
- **URL**: https://doi.org/10.6084/m9.figshare.c.5904905.v1 — CC0, no registration.
- **Contains**: **616 expert RANO-rated timepoints across 91 patients** (PD 253, SD 97, CR 27, PR 20,
  Pre-Op 92, Post-Op 124), each with a written rationale carrying target-lesion millimetre measurements;
  longitudinal T1/T1ce/T2/FLAIR with DeepBraTumIA and HD-GLIO-AUTO segmentations; 4,792 × 152 pyradiomics
  rows; IDH and MGMT (qualitative + % methylation); survival in weeks. The study-level join was verified:
  563 rows match on (Patient, Timepoint), of which **373 rows across 78 patients carry a real PD/SD/PR/CR
  label**. `[REPORTED, join independently reproduced by the hunter]`
- **Validation it enables**: the response/tumour-burden trajectory axis only.
- **What it cannot do**: LUMIERE is **structural MRI only** — confirmed from the paper and from
  `LUMIERE-datacompleteness.csv`, which has exactly four sequence columns (CT1, T1, T2, FLAIR). No
  perfusion, no diffusion, no transport content of any kind. Also: no per-patient treatment record beyond a
  note that two patients received Avastin. The "treatment" is an assumed uniform Stupp protocol.
- **Local**: `data_cache/hunt/upenn_lumiere/`, `data_cache/hunt/dce_mri_response/`.

---

## 3. Tier 2 — partial (potency or delivery axis only)

Ranked by how much a simulator can actually be scored against them.

**1. Human regional intratumoral PK anchors — glioma, real drugs, region-labelled.** `[INSPECTED]`
The strongest *delivery-axis* asset in the correct disease, and the only Tier-2 entry that is
region-resolved.
- `data_cache/hunt/autoradiography_microdialysis/parsed/regional_intratumoral_PK_anchors_literature.csv`
  — 27 rows, PMID- and DOI-anchored, columns `plasma_cmax | plasma_auc | tissue_cmax | tissue_auc |
  tissue_to_plasma_auc_pct | method`. Sources: Portnow 2009 (TMZ, peritumoral brain, intracerebral
  microdialysis, e.g. subject 1 tissue:plasma AUC = 16.4%), Blakeley 2008 (methotrexate,
  contrast-enhancing vs nonenhancing), Jackson 2018 (TMZ, CE vs non-CE), Grossman 2012/2013 (TMZ,
  intratumoral vs contralateral normal-brain ECF), Lee 2025 ABTC (adavosertib, CE vs nonenhancing, tissue
  homogenate LC-MS), Eroglu 2024 (eribulin, peritumoral).
- `parsed/oh2026_neurooncol_TableS3_per_sample_tumor_drug_ratios.csv` — **125 samples across 42 human
  cases**, per-sample tumour:plasma ratios for two probe drugs with **region flags**
  (`region_contrast_enhancing` +39/−69, `region_flair` +73/−35, `region_necrosis`). Ratios: probe 1
  mean 0.697 (range 0–2.05), probe 2 mean 0.150 (range 0–1.47) — a ~4.7× separation between a
  BBB-penetrant and a BBB-excluded compound in the *same* tissue samples. Source
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12962645/ (Neuro-Oncology, `noaf204`).
- **What it validates**: the predicted tumour:plasma partition coefficient, per region, in human glioma.
  This is the numeric target for the delivery axis. It is regional, not spatial — no coordinates.

**2. Mayo GBM PDX panel (Vaubel 2020) — molecular + measured in-vivo response, in glioma.** `[INSPECTED]`
`data_cache/hunt/pdx_encyclopedia/mayo_vaubel2020_TableS7_PDX_response.csv` — **37 PDX lines**, median
survival (days) under Placebo / TMZ / RT / RT+TMZ / BEV with per-arm ratios, MGMT methylation status,
tumour type, **patient OS in months** and the patient's actual treatment. Paired expression at
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132650. This is the best molecular→response
validation set in the right disease with a real in-vivo drug effect. No spatial, no drug concentration.

**3. HGCC — GBM patient-derived cultures, 262 drugs.** `[INSPECTED]` https://hgcc.se/
`HGCC_drug_response_AUC.txt` (262 drug rows × ~40 cultures + control lines) with matched
`HGCC_RNA_expression.txt`, `HGCC_DNA_methylation.txt`, `HGCC_DNA_copy_number_gene_level.txt`,
`HGCC_clinical.txt`. Validates the potency axis (IC50/AUC given free drug concentration) with a molecular
covariate, in glioma. Local: `data_cache/hunt/hgcc_expression/`.

**4. GSE148842 — drug response measured *in intact human GBM tissue*.** `[INSPECTED]`
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148842 — 43 acute slice cultures from **10 GBM
patients**, scRNA-seq under vehicle (n=20), etoposide 2.5 µM (n=8), panobinostat 0.2 µM (n=8), and single
arms of RO4929097 / tazemetostat / ispinesib / Ana-12. The response is transcriptome-wide and
cell-type-deconvolvable, in tissue that still contains its own microenvironment. No drug concentration
measured; slice culture removes the vasculature. Local:
`data_cache/hunt/sc_drug_perturb/`, `data_cache/hunt/slice_organotypic/`.

**5. Zenodo 5795371 — pioglitazone MSI ↔ HPLC absolute-concentration bridge.** `[REPORTED]`
https://zenodo.org/records/5795371 (CC-BY-4.0). Not spatial (aggregated per sample), but it is the only
open artefact that converts an MSI ion abundance into **ug/g** and then cross-validates it against
orthogonal HPLC on the *same* tissue (MSI 12.9 / 28.5 µg/g at 2 h vs HPLC 8.7 / 6.2 / 9.2 µg/g). The
untreated control quantifies to ≈0, bounding the assay noise floor. Local: `verified_7/`.

**6. Large-scale potency reference: DepMap / PRISM / CTRPv2 / NCI-60 ALMANAC.** `[INSPECTED]`
https://depmap.org/portal/ — local cache holds `PRISM_19Q4_secondary_GLIOMA_subset.csv`,
`PRISM_19Q4_secondary-screen-dose-response-curve-parameters.csv`, `CTRPv2_PharmacoSet_zenodo3905470.rds`,
`DepMap_24Q4_Expression_ProteinCoding_TPMlogp1_CNSBRAIN.csv`, `DTP_NCI60_ALMANAC_COMBO_SCORE.zip`.
This is the potency-axis prior, at scale, with molecular covariates — and ALMANAC is the open combination
reference. Cell lines only. (Note: `cancerrxgene.org` returned HTTP 410 on 2026-07-27; do not cite GDSC's
web portal as a live URL.)

**7. Mann 2023 organotypic brain-slice culture dose-response.** `[INSPECTED]`
`data_cache/hunt/slice_organotypic/Mann2023_OBSC_doseresponse_extracted.csv` — 572 values, 7 models ×
11 drugs × {MaxKill, EC10/25/50/75/90, Slope, AUC}. **Caveat I observed directly**: the companion
`Mann2023_TableS3_patient_tumors_UNALIGNED.csv` is explicitly flagged `UNRESOLVED — only N of 11 drug
columns populated; blank cells lost in PDF text layer`. The cell-line-on-OBSC table (TableS1) is usable;
the **patient-tumour table is not**, because it was scraped from a PDF text layer that dropped blanks.
Do not use the patient rows.

**8. TCGA GBM/LGG drug-level RECIST.** `[INSPECTED]` https://portal.gdc.cancer.gov/ —
`tcga_treatment_response/tcga_gbm_lgg_drug_response_RECIST.csv`, 340 rows of
`patient | drug name | best response | start/end day | dose`, i.e. per-regimen RECIST in a cohort with
full molecular profiling. Bulk molecular → clinical response only; no spatial, no drug concentration, and
`prescribed_dose` is `[Not Available]` on the rows I inspected.

**9. ACRIN-6684 FMISO PET.** `[INSPECTED]`
https://www.cancerimagingarchive.net/collection/acrin-fmiso-brain/ —
`radiopharm_pet/acrin6684_fmiso_pet_vs_outcome.csv`, 38 patients with `SUVmean_tumor | SUVmax_tumor |
TBmax | hypoxic_vol_cc | T_cereb_max | hypoxic_vol_cereb_cc | days_to_death | last_assess_day`. A measured
*microenvironment* variable (hypoxia) against survival, in GBM. It is a tracer, not a drug; but if the
simulator models hypoxia as a driver of both transport and response, this is the one open cohort that puts
a number on hypoxic volume and an outcome next to it.

**10. UPENN-GBM CaPTk regional DSC + DTI.** `[REPORTED]`
https://www.cancerimagingarchive.net/collection/upenn-gbm/ — 67 CSVs of {ap-rCBV, PH, PSR} and
{AD, FA, RD, TR} × {enhancing tumour, necrotic core, edema}, 467 subjects, joining **467/467** to the
clinical table. Use it for **transport-coefficient priors**: DTI trace/AD/RD stratified by tumour
sub-region is a measured per-compartment tissue-diffusivity dataset. The endpoint side is thin — OS under
*assumed* standard of care (the open clinical CSV has 12 columns and **not one is a drug, dose, or
regimen**), and the PsP-vs-true-progression label is only n = 60. That powers a correlation, never a
validation. TCIA's NBIA REST API needs no auth, so individual `ep2d_perf` (290 series) and DTI-30dir
(381 series) can be pulled without the 139 GB whole collection.

**Also present but weaker / unranked**: GBM PDX + GeoMx (`zenodo_figshare/gbm_pd1_geomx_z16839828`),
CGGA 693/325 expression + clinical (`cgga_clinical/`), GLASS longitudinal (https://www.synapse.org/glass),
HTAN open tier (https://humantumoratlas.org/), the DARTT-1 pharmacoscopy ex-vivo screen (33 patients ×
~166 drug/combination columns, `pharmacoscopy_exvivo/DARTT-1_Extended_Data_Table_1.xlsx` — image-based
ex-vivo drug response with clinical follow-up, but pan-cancer and blood/effusion-weighted), and the
3D spheroid MALDI imzML test fixture (https://zenodo.org/records/17374882, `verified_8/`, useful for
exercising a 3D diffusion solver against real compartment masks — **no drug, no response**).

---

## 4. Tier 3 — gated

| Asset | Gate | Precise request path | Realistic timeline |
|---|---|---|---|
| **Agar/White adavosertib GBM MALDI-MSI** — measured adavosertib µM per pixel (reported 1–19 µM in GBM84, intratumoral CV 20.9–85.6%) co-registered with phosphoproteomics and t-CyCIF, in GBM PDX **and phase-1 patient specimens**. The exact experiment, in the exact disease. | **No accession exists.** Data statement says verbatim "The MALDI-MSI data are available upon request"; t-CyCIF "uploaded to Omero" with no identifier. | Email the corresponding authors — Nathalie Y. R. Agar (Brigham & Women's) and Sandro Santagata — at the MIT-Mayo PSOC. Phase-1 patient tissue ⇒ expect IRB + DUA, not a link. Separately, the TMT channel→arm design map needed to make PXD018782 interpretable must come from the White lab (MIT). | Months. This is a collaboration request, not a download. Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC8730776/ |
| **TCIA imaging: RHUH-GBM, Burdenko-GBM-Progression, Brain-Tumor-Progression, QIN-GBM-Treatment-Response** | NIH Controlled Data Access (face-reconstruction risk). Clinical CSVs are separately CC BY 4.0 and already downloaded. | TCIA restricted-data request form linked from each collection page: https://www.cancerimagingarchive.net/collection/rhuh-gbm/ · /burdenko-gbm-progression/ · /brain-tumor-progression/ · /qin-gbm-treatment-response/ | Weeks. |
| ↳ **Burdenko** specifically | as above (131 GB) | same | Worth it: DICOM **RTDOSE** gives a true 3D delivered-dose field paired with **651 per-follow-up response labels** across 180 patients. Wrong modality (radiation, not drug) but structurally the exact `dose(x,y,z) → local control` test. |
| ↳ **Brain-Tumor-Progression** specifically | as above (only 3.16 GB) | same | The only TCIA GBM collection shipping **pre-computed** nCBF/nrCBV/srCBV maps, in the same patient post-CRT and at progression. Small enough to actually work with. |
| ↳ **QIN-GBM-Treatment-Response** | as above | — | **A trap.** Despite the name it contains **no response data** — it is a 2–5 day double-baseline test-retest repeatability study, both scans before therapy. Its only use is quantifying the Ktrans/rCBV measurement noise floor. Do not request it expecting response labels. |
| **HTAN WUSTL 10x Visium Level 1/2** | dbGaP | dbGaP controlled-access application via the HTAN portal (https://humantumoratlas.org/); Level 3 is on Synapse and open. 25 participants. | Weeks–months. |
| **METASPACE raw imzML/ibd** | `open_registration` — free account; the API also returns a `Download_Limit_Reached.txt` stub under load. | Register at https://metaspace2020.eu/ . The derived 16-bit ion-image PNGs need no auth and are what has been harvested. | Same day, plus patience. |
| **Clinical-trial IPD (Vivli / CSDR / NCI PDS)** | Data-request application + statistical analysis plan; usually a secure enclave, no export. | Local cache holds `clinical_trial_open/Vivli_glioblastoma_DataCite_DOIs.json`, `CSDR_All-Sponsor-Studies.xls`, `PDS_CNS_datasets_metadata.json`, `PDS_how_to_access_and_download.pdf`. | 1–6 months, and the enclave usually blocks the workflow a simulator needs. |
| **UPENN-GBM / UCSD-PTGBM analysis-ready NIfTI** | Not gated — but IBM Aspera faspex, 69 GB and 33.1 GB, all-or-nothing, no subsetting. | Needs Aspera Connect or `ascli`. For subsets use the NBIA REST API instead (no auth): `getSeries?Collection=UPENN-GBM`, then `getImage?SeriesInstanceUID=<uid>`. | Hours of transfer. |

**Confirmed dead ends** (so nobody re-spends the time): ProteomeXchange contains exactly **one** imzML
dataset in the entire federation (PXD001283 — mouse bladder phospholipids, no drug, no tumour). MassIVE
holds no MSI-with-drug. Zenodo full-text `imzML` returns 36 records; `imzML AND (drug OR dosed OR
pharmacokinetic)` returns **zero**. The famous 3D paclitaxel-penetration paper (Sci Rep srep37027)
deposited three `.avi` rotation videos and a 13 KB PDF — its concentrations exist only as five binned
colour ranges in figures. A figure is not a dataset.

**Lesson from the one success**, worth encoding into future hunts: the best asset found (T1.1) is invisible
to every repository search. Its repository deposits contain code and ROI summaries; the actual 376,736-pixel
drug map is sheet 34 of a 130 MB journal Source Data workbook, in CSV-shaped Excel, not imzML. Likewise the
second-best (T1.2, T1.4) are `.xlsb` workbooks on Zenodo that no `imzML` query returns. **Search the
Source Data, not just the accession.**

---

## 5. Recommended immediate action

**Run on the Moncayo HGSOC PDE set now** — `Zenodo 17610220` + `GSE281519` + the Nat Commun Source Data.
It is already on disk, it needs no gate, and it is the only open asset in existence that carries a
per-pixel measured drug field **and** a molecular state on the same tissue.

Local paths (all present and parsed):
- `/Users/samarthbasanth/tme-sim/data_cache/hunt/geomx_dsp/NatCommun2026_SourceData_SupFig2_MSI_perpixel.csv.gz`
  — **376,736 MSI pixels × 19 cols, 20 µm pitch, X/Y, 3 drugs, raw + heavy-labelled IS + normalised** (13 MB,
  extracted here from Source Data `41467_2026_70558_MOESM6_ESM.xlsx`)
- `/Users/samarthbasanth/tme-sim/data_cache/hunt/geomx_dsp/zenodo17610220_ROI_drug_measurements.csv` (114×13)
- `/Users/samarthbasanth/tme-sim/data_cache/hunt/geomx_dsp/zenodo17610220_ROI_drug_plus_expression.csv.gz` (114×11,237)
- `/Users/samarthbasanth/tme-sim/data_cache/hunt/geomx_dsp/zenodo17610220_modelling.Rmd` (the authors' own `lme4` spec)
- `/Users/samarthbasanth/tme-sim/data_cache/hunt/geomx_dsp/GSE281519_raw/` (133 `.dcc.gz`, ROI-level DSP counts)

**Exactly what analysis it supports, and nothing more:**

0. **Pixel-level penetration fit (do this first).** Reconstruct each of the 18 slices onto its 20 µm grid
   from `X`/`Y`, solve diffusion-into-explant from a 20 µM Dirichlet bath boundary for 24 h, and score the
   predicted field against `norm_<drug>` per pixel. Report per-slice Spearman ρ, the CV of the measured vs
   predicted field, and hotspot-size distributions. This is the single most falsifiable open test
   available: 376,736 measured points, an internal standard on every one of them, and three compounds.
1. **Accumulation-heterogeneity fit.** Simulate niraparib and rucaparib uptake into a diffusion-limited
   explant geometry at 20 µM bath concentration for 24 h, with pH-partitioning / lysosomal sequestration
   switched on and off. Score predicted per-ROI accumulation against `drug_measure_norm` inside the
   authors' own frame: `measure ~ prediction + (1|Patient) + (1|Patient:Slice)`. Report the variance
   explained beyond the patient/slice random effects, because that is the only part your model can claim.
2. **The compound-specificity falsification test.** The model must produce *more* accumulation
   heterogeneity for the weak bases (niraparib, rucaparib) than for olaparib. If it produces the same
   heterogeneity for all three, the lysosomal-reservoir mechanism is missing and the result is negative —
   which is a publishable outcome, not a failure.
3. **Exposure → molecular response coupling.** Regress predicted local exposure against the apoptotic and
   lysosomal signature scores computed from the 11,224-gene ROI matrix. The paper's own GSEA
   (Source Data `Fig. 2c`: APOPTOSIS, size 154, NES 2.278, FDR q = 0) is the benchmark to beat or match.

**Mandatory companion run, same week:** the METASPACE GIST imatinib set
(`data_cache/hunt/maldi_msi_drug/`, `verified_0/`). Moncayo has no vasculature; GIST is in vivo, in human,
with a per-patient plasma boundary condition and a within-patient matched-normal control. Predict the
tumour:normal partition coefficient for the 6 matched patients from plasma + dose, and report it against
the measured 0.02–0.72 range. Because plasma explains none of the tissue variance (r = 0.03, p = 0.93),
any predictive skill here is attributable to the transport model and nothing else. State the n (6 pairs)
in the same sentence as the result.

**If the project must stay in glioma**, the immediate run is instead
`autoradiography_microdialysis/parsed/regional_intratumoral_PK_anchors_literature.csv` +
`oh2026_neurooncol_TableS3_per_sample_tumor_drug_ratios.csv`: predict tissue:plasma AUC% per region
(contrast-enhancing vs non-enhancing vs peritumoral) for TMZ, methotrexate and adavosertib, and score
against the PMID-anchored measured values. It is regional rather than spatial, but it is the only measured
human glioma drug-exposure data that exists openly.

---

## 6. Honest limits — what still cannot be validated with open data

1. **The delivery→kill link cannot be validated in space, at all.** You can validate spatial delivery
   (non-glioma, small n, ex vivo or xenograft) or you can validate response (glioma, no drug concentration,
   assumed uniform Stupp protocol). No open dataset does both in the same tissue. Any figure that puts a
   predicted penetration map next to an observed response must be captioned to say the two came from
   different specimens, or different studies, or both.
2. **There is no glioma + temozolomide spatial drug map, open or otherwise.** METASPACE has zero TMZ
   datasets. Every spatial-drug asset above requires a cross-indication and usually cross-drug transfer
   argument. That argument is defensible for the *physics* — diffusion, partition, efflux, sequestration —
   and is not defensible for the *biology*. It must be made explicitly and defended, not assumed.
3. **Ion intensity is not concentration.** Only the Mario Negri pioglitazone deposit provides an
   MSI↔HPLC absolute calibration, and it is a different drug in a different tissue. The GIST imatinib maps
   are arbitrary ion counts; the Moncayo values are IS-normalised ratios. Reporting absolute µM or ng/g
   from any of them without a stated calibration path would be fabrication.
4. **Contrast-agent perfusion is not drug delivery.** Every "perfusion" map in TCIA is gadolinium. Ktrans
   and CBV describe the vasculature's transport capacity for a hydrophilic, BBB-excluded, ~550 Da bolus
   tracer. Nothing in the open imaging world measures a therapeutic in tissue.
5. **Ex vivo systems delete the source term.** Moncayo (bath-dosed explants), GSE148842 (slice cultures)
   and Mann 2023 (organotypic slices) all remove the vasculature. They validate diffusion, uptake and
   sequestration. They cannot validate perfusion-driven delivery, interstitial pressure, or vessel
   normalisation — which are usually the whole point of a TME transport model.
6. **n is small everywhere it matters.** 3 patients / 114 ROIs (Moncayo); 6 matched tumour-normal pairs
   and 9 plasma-linked sections (GIST); 42 cases (regional PK); 60 PsP labels (UPENN-GBM); 37 PDX lines
   (Mayo). These support falsification of a clearly-wrong model. None of them supports a claim that one
   plausible model outperforms another.
7. **Pairing is arm-level, not specimen-level, in the one dataset with efficacy.** Zenodo 5776125's imaged
   mice are not its followed mice (zero animal-ID overlap in two of three models).
8. **The strongest matching study is closed.** The Agar/White adavosertib GBM platform is the experiment
   this project would design if it could design one, and its MSI arm has no accession and its
   phosphoproteomics deposit has no quantitation. Getting it is a months-long collaboration ask on phase-1
   patient tissue.
9. **Coverage caveat on this document.** It synthesises 4 fully-reported hunt angles, 3 independent
   verifications, and my own direct inspection of the remaining on-disk caches. Entries marked
   `[INSPECTED]` were checked by opening the files, not by re-running the original searches; entries marked
   `[REPORTED]` have not been independently re-checked and should be treated as leads until they are.
