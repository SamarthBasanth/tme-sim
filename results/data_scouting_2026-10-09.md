# Data scouting, 2026-10-09 — what is newly available, what is already held, and what would actually move a claim

**Status: SCOUTING PASS. Nothing here is a result, a test, or a pre-registration.**
No code was changed. No existing file in `results/` or `archive/` was edited.
This file is the only thing created. Nothing below is registered in
`tmesim/validation/scorecard.py`.

**Brief, in the owner's words:** *"see if you can scour the internet for better
resources / databases / tests to train on."*

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

*(Closing guard in §8.)*

---

## 0. The rules this pass was run under, and the method

**Scope rule.** Vivum is a transport simulator validated on the molecular-dynamics
force-field analogy: verify the transport, never claim clinical response. So a
dataset's value here is *not* how much signal it carries — it is whether it bears
on a transport claim. **No clinical-outcome or treatment-response dataset was
sought, and §6 says plainly why the obvious "train on this" answer is the wrong
direction for this product.**

**Provenance rule.** Every accession, byte count, sample count, gene-panel
membership and PMID below was retrieved live on **2026-10-09** from a primary
record, not from a title or a search snippet. The specific instruments used:

| what | how it was checked |
|---|---|
| dataset existence, sample counts, design | NCBI E-utilities `esearch`/`esummary` on `db=gds`, and the GEO `targ=self&form=text&view=brief` record for individual GSE/GSM |
| exact file sizes and per-sample file lists | the series' own `suppl/filelist.txt` on `ftp.ncbi.nlm.nih.gov`, which carries a byte count per file |
| **gene-panel membership** | downloaded the actual feature table and grepped it — `features.tsv.gz` for MEX deposits, and for Xenium zarr deposits the `cell_features/.zattrs` `feature_keys` array unzipped out of `cell_feature_matrix.zarr.zip` |
| Visium HD bin pitch | downloaded `barcodes.tsv.gz` and read the barcode string format |
| imaging-MS corpus | METASPACE GraphQL endpoint, `allDatasets`/`countDatasets`/`dataset` queries |
| literature | PubMed `efetch` abstracts and EuropePMC `fullTextXML` |

Anything that could not be retrieved is marked **UNVERIFIED** with the specific
thing that blocked it. **No accession, PMID or number in this file was inferred
or reconstructed.** Where a claim rests on a pattern rather than a statement
(§2.3, the matched-pair claim), that is said explicitly.

**This pass ran three strands in parallel**, and not all of it is first-hand.
I took the deposited-tissue strand myself; two parallel strands covered the
transport literature and the single-cell reference atlases. **§7 states exactly
which statements I verified against a primary record and which were reported to me
and are not independently checked.** Read that before quoting any number onward.
As a rule of thumb: accessions, byte counts, panel membership, PMIDs, journals,
trial identifiers and principal findings are verified; **per-patient values and
table contents are not.**

**Value-per-effort rule, as instructed.** A 2 GB open dataset beats an 800 GB
gated one. This turned out to matter more than expected: three of the best finds
have *nominal* sizes of 40–89 GB and *actual* costs of 190–460 MB, because the
bulk is morphology OME-TIFFs that this project does not need and GEO allows
per-file fetches. The nominal size is a bad proxy and §2 gives the real number.

---

## 1. What the repo already holds — read this first so nothing is double-bought

Measured from `data_cache/` and cross-read against `docs/DATA_INVENTORY.md`,
`docs/real_data/*.md`, `results/DATA_HUNT.md` and
`results/external_anchor_survey.md`. **This project has already been scouted
hard** — a 28-angle hunt (`results/DATA_HUNT.md`, 2026-07-27, 164 datasets
surveyed) and an anchor survey (`results/external_anchor_survey.md`, 2026-08-09).
The opening this pass exploits is simply **time**: the finds in §2.1–§2.3 were
deposited or released *after* those passes.

### 1.1 Tissue already on disk

| path | what it is | size |
|---|---|---|
| `xenium_gbm/` | 10x Xenium GBM, 319-gene panel, 40,887 cells — *the flagship section* | 30 MB |
| `xenium_gbm_io/` | 10x Xenium GBM FFPE, 480-gene IO panel, 816,769 cells | 517 MB |
| `visium_gbm/` | 4 Visium sections: 10x parent + GSE194329 GBM1/2/3 | 242 MB |
| `gbmspace_visium/` | GBM-Space (Sanger), 31 section dirs materialised | 1.1 GB |
| `greenwald_gbm_visium/` + `greenwald_gbm_codex/` | Greenwald GBM: 2 Visium (ZH1007 inf + nec), 1 CODEX OME-TIFF | 151 MB + 878 MB |
| `gse324268_xenium/` | Xenium GBM, hypoxia-oriented panel (primary GBM + organoid-on-slice) | 320 MB |
| `xenium_breast_janesick/` | matched Visium+Xenium breast reference | 183 MB |
| `ivy_gap/` | Ivy GAP anatomic-structure RNA-seq (**zip unextracted, loader raises**) | 36 MB |

### 1.2 Reference / tabular already on disk

`scref_gbm/` (Darmanis GSE84465, 3,589 cells — the deconvolution reference),
`drugresponse/` (HGCC + Jacob supplements), `gdsc_{glioma,breast}_dose_response.csv`,
`tcga_expression/` + `tcga_gbm_clinical.json`, `brain_slice/`, `pruijn2008/`.

### 1.3 Anchors already held — do not re-find these

Primeau 2005 (16361566, **in-sample**: α is solved to it), Patel/Trédan/Tannock
2013 (23680920), Portnow 2009 (19861433), Morse 1985 (3855310), Nakagawa 1996
(8548876), Hicks 1998/2007 (9806526, 17869669), Pruijn 2008, Boucher 1997,
Mader 2003, Blakeley 2008 (18787762), Jackson 2018 (29332604), Oh 2026
(125 samples / 42 HGG cases), Lee 2025 ABTC (39820282), Eroglu 2024 (39422741),
Abu Sammour 2019 GIST imatinib MSI, Moncayo 2026 PARPi MSI (41844641), plus the
breast set and the 24-entry `OUT_OF_SCOPE` table in `scorecard.py`.

### 1.4 The limitation this pass was pointed at

`results/CLAIM_INVENTORY.md` §6 ranks **"a second section, same tumour type"** as
item 1 of what would most change the project's self-assessment, calling it *"the
cheapest, largest improvement available."* §2.1 below is the direct answer to
that item, and it is larger and cheaper than expected.

---

## 2. Ranked candidates

Ranked by value-per-effort, not by size or novelty.

**Read §2.0 first.** The two highest-value items in this file turned out to come
from **the same clinical trial**, which neither strand of this pass was looking
for and which changes what the top find actually is.

| # | candidate | what it adds | real cost | licence | bears on |
|---|---|---|---:|---|---|
| **0** | **NCT04528680, both halves** — Habashy 2025 (PMID **40495421**) free-interstitial carboplatin/fluorescein time-courses **+** GSE339006 Xenium 5K tissue, same trial, same compartment | **a measured free-interstitial loss rate `k` AND reconstructable tissue of the compartment it was measured in** | ~382 MiB + free | GEO open / PMC | `alpha_identifiability.md` §4(a), §2.0 |
| **1** | **GSE318562** — 22 Visium GBM sections, 9 donors, core→margin | **7× the Visium tissue, drop-in format** | **674 MB** | GEO open | the *n*=1-in-tissue limitation, §2.1 |
| **2** | **GSE335275 + GSE316917** — 15 Xenium GBM sections, 362-gene panel **with CA9/VEGFA/ADM/LOX + PECAM1** | measured hypoxia instead of imputed; 15 sections | **187 MiB + 452 MiB** | GEO open | hypoxia imputation, §2.2 |
| **3** | **the 6 matched Visium↔Xenium pairs** inside #1 and #2 | a **6-fold replicated held-out deconvolution test** | included in #1+#2 | GEO open | the RETRACTED Pearson 0.997, §2.3 |
| **4** | **GSE339006** — Xenium **Prime 5K** (5,101 targets), 10 sections, 5 recurrent GBM patients, paired sonicated/non-sonicated peritumoral brain | the **only** panel found carrying PECAM1+CLDN5+PLVAP **and** HIF1A/CA9/SLC2A1 **and** ABCB1/ABCG2 **and** MGMT at once | **381.6 MiB** | GEO open | panel gap + the TMZ anchor's own tissue, §2.4 |
| **5** | **Evans 2008**, PMID **19138031** — EF5 vs radial distance from individual vessels, 6 human brain tumours | a **falsification-grade** result on the distance-from-vessel paradigm, in human GBM | free | CC (PMC) | §2.5 — read this one even if nothing is downloaded |
| 6 | **GSE290724** — **Visium HD**, 2 µm bins, human GBM, peritumoral/edge/core | a 4th platform at 2 µm | 252 MB | GEO open | §2.6, with a real caveat |
| **7** | **GSE302502** — Xenium 468-gene, 12 TLS+ gliomas | **CD31 *and* CD34 on one panel**, + PLVAP/CLDN5/VWF, full hypoxia, mural markers — but **cropped to TLS regions** | **29.3 MiB** | GEO open | the CD31-vs-CD34 convention gap, §2.7 |
| 8 | GSE305393 — CosMx SMI **6K** panel, 18 glioma + 9 brain-met, TMA | a 5th platform; TMA cores limit extent | 4.4 GB | GEO open | §2.8 |
| **9** | **Wälchli 2024 GSE256493** — human brain vasculature atlas **including GBM vessels** | **49,999 GBM endothelial cells** in place of the current **51** | **118 MB** | GEO + UCSC open | the deconvolution reference, §2.9 |
| **10** | **Johnson/Sanai 2026** (PMID **41206763**), **Schiff 2026** (PMID **41099363**) | **absolute unbound drug concentration in human glioma**, two compounds 3 orders apart | free | PMC | §2.10 |
| 11 | Salavati 2024 (PMID **38435818**) + Elmghirbi 2018 (PMID **29524243**) | measured hydraulic conductivity `K` and an IFP boundary-layer length — **25–40× apart** | free | PMC | §2.10 |
| 12 | Siletti 2023 CELLxGENE — `Supercluster: Vascular` | 9,932 normal-brain vascular cells, and repairs the 88/85/21-cell columns too | **47 MB** | open | §2.9 |
| — | Cribaro 2021, PMID 33579378 — 3D confocal human GBM vasculature | **checked and rejected**, §4 | — | CC BY | nothing extractable |

---

### 2.0 — The convergence: one clinical trial supplies both halves **[READ THIS FIRST]**

Two independent strands of this pass — one scouting deposited tissue, one scouting
published transport measurements — landed on **the same trial, NCT04528680**,
without being pointed at each other. Both were verified from primary records:

| half | record | what it is |
|---|---|---|
| **tissue** | GEO **GSE339006**, released 2026-10-01 | Xenium **Prime 5K** (5,101 gene targets), **10 sections**, paired **sonicated / non-sonicated peritumoral brain** from **5 recurrent GBM patients**, biopsied 45 min after sonication. §2.4. |
| **drug** | **Habashy KJ et al., Clin Cancer Res 2025 Aug 14;31(16):3562-3570. PMID 40495421**, DOI [10.1158/1078-0432.CCR-25-0080](https://doi.org/10.1158/1078-0432.CCR-25-0080), PMC12310118 | **microdialysis** free-interstitial carboplatin and fluorescein, **4 patients**, catheters in **sonicated and non-sonicated brain regions**, **drug levels over 24 hours**. |

The Habashy abstract names the trial verbatim: *"We analyzed four patients who
underwent intraoperative LIPU/MB with intravenous administration of carboplatin
and fluorescein in the NCT04528680 clinical trial. Microdialysis catheters were
implanted into sonicated and nonsonicated brain regions, and drug levels were
measured over 24 hours."* The GSE339006 series summary names the same trial and
the same perturbation. Both are **human peritumoral brain**; both carry the
**paired sonicated/non-sonicated** design.

**Why this is the most valuable thing in this file.**
`results/alpha_identifiability.md` §4 names exactly what the calibration is
missing, and item **(a)** is *a free-interstitial first-order loss rate `k` for
the anchor drug from microdialysis with a resolvable washout*. Habashy is that
measurement — free interstitial (microdialysate, not homogenate, not CSF),
absolute, human brain, with a 24-hour time course. And GSE339006 is
**reconstructable tissue of the exact compartment the `k` was measured in, under
the same perturbation**, which no anchor in this repo has ever had.

`external_anchor_survey.md` §5 concluded that nothing public measures what this
engine needs and that *"this may require generating data, not finding it."*
**This does not overturn that** — see the limit below — but it is the closest
anyone has come, and it was generated by someone else.

**The limit, stated before anyone gets excited.** Microdialysis measures a
**time course at a point** — wherever the catheter sits. It is **not a spatially
resolved field**, so it still cannot score the solver's *spatial* prediction. What
it can do is pin the **loss term**, which `results/pk_reality.md` records as having
been ~300× wrong, and which `alpha_identifiability.md` says is one of the
unidentified numbers in the whole calibration. That is a parameter the engine
currently has no measured value for, in the right compartment, in the right
species. It is also a genuinely *pre-registrable* test: the reported sonication
invariance of the loss rate (influx up 3.1-fold, `k` statistically unchanged)
is a structural prediction a solver either reproduces or does not.

**Two things I have NOT verified and which matter:**
1. **Whether the 4 microdialysis patients and the 5 Xenium patients overlap.** The
   GSE339006 sample titles name patients `01-01-110`, `01-01-121`, `01-01-135`;
   Habashy reports 4 patients without, in the abstract, giving matching
   identifiers. Same trial is confirmed; **same patients is UNVERIFIED.** If they
   overlap, this becomes a per-patient tissue↔PK pairing, which would be
   unprecedented for this project. If not, it is still the same protocol and
   compartment.
2. **The per-patient numbers.** The t½, AUC and clearance values reported to me by
   the scouting strand are **not reproduced in this file**, because I verified the
   paper's identity, trial, design, compartment and the direction of its findings
   from the primary abstract, but **did not open Tables 1–3 myself**. The abstract
   confirms the 3.1-fold brain-to-plasma AUC increase and that *"TMZ concentrations
   remained below their plasma levels, as parenchymal drug clearance mirrored
   plasma clearance"*. **Anyone adopting a number from this paper must read the
   tables; do not quote a value via this file.**
3. **The TMZ arm is not independent.** The abstract says *"Published microdialysis
   data of TMZ without LIPU/MB were used for comparison"* — that is Portnow 2009
   (PMID 19861433), which this repo **already holds and which is in-sample for the
   partition prior**. Only carboplatin and fluorescein are new.

---

### 2.1 — #1 GSE318562: 22 Visium sections, 9 GBM donors, core to margin **[BEST SINGLE FIND]**

- **Accession:** GEO **GSE318562**, part of the study *"Mapping the spatial
  architecture of glioblastoma from core to edge delineates niche-specific tumor
  cell states and intercellular interactions"*. BioProject **PRJNA1405549**.
  Preprint: bioRxiv **2025.04.04.647096**. **No PubMed ID is linked in the GEO
  record as of 2026-10-09** — treat it as a preprint, not a peer-reviewed source.
- **Species + disease:** *Homo sapiens*, **IDH-wildtype glioblastoma**, stated per
  sample as `Brain tumor section (IDH wildtype Glioblastoma)`.
- **Assay:** 10x Visium Spatial Gene Expression, platform GPL24676. Per the
  series' own extraction protocol, verbatim: *"Using image-guided surgical
  biopsies and 5-aminolevulinic acid (5-ALA) fluorescence, we collected 22 samples
  from tumor regions ranging from the core to the margin in nine patients with
  IDHWT GBM"*. **The first 13 samples are fresh-frozen/OCT (Visium 3' Rev. D); the
  next 9 are FFPE (Visium FFPE Rev. E)** — a mixed-chemistry cohort, which matters
  for any pooled comparison.
- **Size:** `GSE318562_RAW.tar` = **706,416,640 B (674 MB)** for all 22.
  Per sample: `filtered_feature_bc_matrix.tar.gz` 4.9–43 MB +
  `spatial.tar.gz` 5.9–12 MB.
- **Licence / access:** NCBI GEO, open, direct HTTPS download, no login, no DUA.
- **Why this is the top find — format, and this was tested, not assumed.** I
  downloaded the smallest sample, **GSM9497620 (`TWBK-B186-01`)**, both files
  (4,859,662 B + 7,036,399 B) and listed and extracted them:

```
TWBK-B186-01_spatial/{scalefactors_json.json, tissue_positions_list.csv,
                      tissue_hires_image.png, tissue_lowres_image.png,
                      aligned_fiducials.jpg, detected_tissue_image.jpg}
TWBK-B186-01_filtered_feature_bc_matrix/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}
```

  That is **exactly the stock Space Ranger shape `tmesim/data/real/visium.py`
  reads** — a `spatial/` folder with `tissue_positions*.csv` + `scalefactors_json.json`
  and a MEX triple. Measured contents of that one sample:
  `scalefactors_json.json` gives `spot_diameter_fullres` **74.78797327502073**
  → `umpp` = 55/74.788 = **0.7354 µm/px** on the loader's own convention;
  `tissue_positions_list.csv` has 4,992 rows of which **1,420 are in-tissue**;
  the matrix is **36,601 genes × 1,420 spots, 1,214,948 non-zeros**. The only
  wrinkle is that the inner directories are name-prefixed
  (`TWBK-B186-01_spatial`, not `spatial`), so extraction needs the same
  `--strip-components=1` the existing GSE194329 recipe in
  `docs/DATA_INVENTORY.md` A3 already uses.
  Contrast `results/gbmspace_ingestion.md`, where GBM-Space was found, inspected
  and *blocked by format* because it ships `.h5ad`. **This one needs no new
  ingester**, and that is now a measurement rather than an inference from a
  filename.
- **Whole transcriptome, so the hypoxia and vessel genes are all there.** Verified
  present in that sample's `features.tsv.gz` (36,601 features): **PECAM1, VEGFA,
  CA9, HIF1A**, GFAP, SOX2, MKI67, CD68, MBP. Unlike the Xenium panels, a Visium
  deposit has no panel gap at all — see §2.2 for what that does and does not mean.
- **Donor structure** (parsed from the 22 sample titles — all verified present in
  `filelist.txt`):

| donor | sections |
|---|---|
| B172 | 01, 02 |
| B176 | 01-02, 02-01, 02-04 |
| B178 | 02-01 |
| B185 | 01, 02, 03 |
| B186 | 01, 02, 03, 04 |
| ST073021 | 02, 04 |
| WU1220 | 01, 02, 03 |
| WU1221 | **Core2, Edge1** |
| WU1227 | **Core5, Edge2** |

  Two donors carry **explicitly labelled matched core/edge pairs**, and five carry
  3–4 sections of the same tumour. That is within-tumour replication, which is a
  different and stronger thing than more tumours: it is what would let a verdict
  be separated from section-to-section noise. `results/sectioning_bias_2026-09-20.md`
  and `results/multisection_generalisation_2026-09-20.md` are exactly the files
  that needed it — the latter found the per-cell verdict takes three different
  values across sections under byte-identical thresholds, on 7 sections of mixed
  provenance. Here are 22 from one protocol.

**One sentence on the claim:** this is the direct, cheap answer to
`CLAIM_INVENTORY.md` §6 item 1, and it would let every "*n* = 1 in tissue" caveat
in the inventory be either removed or confirmed as a real spread.

---

### 2.2 — #2 GSE335275 + GSE316917: 15 Xenium GBM sections on a panel that *contains the hypoxia axis*

Same study as #1, two sibling series.

- **Accessions:** GEO **GSE335275** (6 sections) and **GSE316917** (9 sections).
- **Assay:** Xenium In Situ, FFPE, 5 µm sections. Panel, verbatim from the
  GSE316917 design: *"a panel of oligonucleotides targeting transcripts from 362
  genes, comprising 266 genes from the predesigned 10X Genomics human
  brain-and-glioblastoma panel and 96 custom genes indicated by our scRNA-seq
  analyses"*.
- **Panel membership — measured, not assumed.** Unzipped
  `cell_features/.zattrs` from `GSM9809297_WU1221_Edge_cell_feature_matrix.zarr.zip`
  (3.5 MB): 542 feature rows, of which **362 of type `gene`**, 20,859 cells in
  that section. Grepping `feature_keys`:

| present | absent |
|---|---|
| **PECAM1**, FLT1 | VWF, PLVAP, CLDN5, CD34 |
| **VEGFA, CA9, ADM, LOX** | HIF1A, SLC2A1, NDRG1 |
| MKI67, AQP4, SOX2, EGFR, OLIG2 | GFAP, MBP, MGMT |
| CD68, AIF1, PTPRC, PLP1, POSTN | CD3E, ABCB1, ABCG2, COL1A1 |

  **This is the point, stated precisely, because the loose version of it is
  wrong.** `results/landscape_cancer_census.md` records, repo-verified against the
  panel JSONs, that **VEGFA / CA9 / HIF1A are absent from both stock 10x GBM
  Xenium panels** this project uses — which is why
  `results/hypoxia_imputation_model.md` and
  `results/panel_gap_imputation_sweep.md` exist. On this panel **CA9, VEGFA, ADM
  and LOX are measured**, alongside PECAM1 and FLT1. HIF1A is still absent, so the
  axis is partial.

  **What this is not:** it is *not* "the first measured hypoxia transcripts in this
  project." Visium is whole-transcriptome, so the four Visium sections already in
  `data_cache/visium_gbm/` — and the 22 in §2.1 — carry CA9/HIF1A/VEGFA already
  (verified above). **The panel gap is a Xenium problem specifically.** The real
  gain here is therefore narrower and should be described as what it is:
  **hypoxia and vessel transcripts at single-cell resolution rather than pooled
  over a 55 µm spot**, on 15 sections. That matters because the flagship
  reconstruction and every cell-resolved result run on Xenium, not Visium — so the
  imputation is in the path that produces the headline numbers, and this is the
  tissue that could retire it there.
- **Size — nominal vs real:**

| series | full tar | **what you actually need** |
|---|---:|---:|
| GSE335275 | 59,217,838,080 B (**55 GB**) | **196,199,849 B = 187.1 MiB** |
| GSE316917 | 96,007,546,880 B (**89 GB**) | **451.8 MiB** |

  The "what you need" column is `cells.parquet.gz` + `cell_feature_matrix.zarr.zip`
  + `cell_boundaries.parquet.gz` + `analysis.zarr.zip`, summed from `filelist.txt`.
  **96–99% of both deposits is `morphology.ome.tif.gz` at 4–12 GB per section**,
  which this project does not need for cell-resolved geometry — the same judgement
  already recorded for `xenium_gbm_io`'s 13.15 GB `_xe_outs.zip` in
  `docs/DATA_INVENTORY.md` A2. GEO serves per-file, so **15 Xenium GBM sections
  cost ~640 MiB, not 144 GB.**
- **Donors:** GSE335275 = 3 donors (B189/ST073021, WU1220, WU1221), with per-sample
  clinical fields present in the GEO record (sex, age, tumour location, tumour
  size, MGMT methylation status, EGFR). GSE316917 = sections GBM024D, GBM024E,
  GBM024H, GBM024I, GBM025A, GBM026I, GBM026M, GBM030, GBM034 → **6 donors, and
  GBM024 alone contributes 4 sections of one tumour** (GSM9462339 confirms
  `patient id: GBM024`, male, 50, right frontal, MGMT unmethylated).
- **GSE316917 has no matched Visium partner, and this is checked, not assumed.**
  Its donor labels (GBM024/025/026/030/034) do not intersect the GSE318562 Visium
  donor labels (B172/B176/B178/B185/B186/ST073021/WU1220/WU1221/WU1227) at all —
  these are separate cohorts (GSE335275 is deposited as `xenium_cohort3`). So the
  matched-pair test in §2.3 is confined to GSE335275's 6 sections; GSE316917's 9
  are additional tissue only.
- **Licence / access:** GEO, open, direct download, no DUA.
- **Format cost, stated honestly:** `cell_feature_matrix.zarr.zip` is **zarr, not
  the MEX triple** `singlecell_spatial.load()` reads, and `cells.parquet.gz` is
  parquet, not `cells.csv.gz`. This needs an adapter. It is a small adapter — the
  zarr holds a plain CSR (`indices`/`indptr`/`data`) plus `feature_keys`, which is
  exactly what the loader wants — but it is not zero, and `gbmspace_ingestion.md`
  is this repo's own cautionary tale about assuming a format fits.

**One sentence on the claim:** this would replace an imputed hypoxia field with a
measured one on 15 sections, which converts `hypoxia_imputation_model.md` from a
modelling assumption into something scoreable.

---

### 2.3 — #3 The matched pairs, and one donor matched THREE ways: the held-out test

This is the **benchmark/held-out-test** item in the brief, and it falls out of #1
and #2 at almost no extra download cost.

**There is a third sibling series, and it changes the design.** **GSE318564** —
*"...[scRNA-seq]"*, 13 samples, same study — is scRNA-seq from the same cohort.
`GSE318564_RAW.tar` is **1.0 GB**. Its sample titles, retrieved via `esummary`:
`B178_1`, `B178_2`, `B178_3`, `B183_1`, `B183_2`, `B183_3`, **`B189`**,
`WU-1225-core`, `WU-1225-edge`, `WU-1226-coreCD45neg`, `WU-1226-edgeCD45neg`,
`GBM034_Core`, `GBM034_Edge`.

Cross-referencing all three series by donor label gives:

| donor | scRNA (GSE318564) | Visium (GSE318562) | Xenium (GSE335275 / GSE316917) |
|---|---|---|---|
| **B189** | ✅ `B189` | ✅ `ST073021-02`, `ST073021-04` | ✅ GSE335275 — **`patient id: B189`** |
| B178 | ✅ `B178_1/2/3` | ✅ `B178-02-01` | — |
| GBM034 | ✅ `GBM034_Core/Edge` | — | ✅ GSE316917 `GBM034` |
| B183, WU-1225, WU-1226 | ✅ | — | — |

**Donor B189 is matched three ways — its own scRNA reference, its own Visium
sections, and its own Xenium sections.** The B189↔ST073021 identity is not
inferred from the label: the GSE335275 GEO record for `HumanBrain-Tumor-ST073021-02`
and `-04` carries the characteristic **`patient id: B189`** verbatim, and
GSE318564 has a scRNA sample titled `B189`.

**Why three ways matters more than six pairs.** The deconvolution step has two
error sources that a two-way comparison cannot separate: the **reference** may be
wrong for this tumour (Darmanis, 4 other donors), and the **unmixing** may be
wrong. With B189 you can build the signature matrix from *that patient's own*
scRNA, deconvolve *that patient's* Visium spots, and score against *that
patient's* Xenium single-cell fractions. That removes the reference-mismatch
confound entirely, and any residual disagreement is the unmixing itself. There is
no cheaper way to isolate it, and ~1.7 GB (GSE318564 + GSE318562) buys the whole
design.

All six GSE335275 Xenium sections carry a section label that **also appears, on the
same patient id, in the GSE318562 Visium set**:

| Xenium (GSE335275) | patient id | Visium (GSE318562) |
|---|---|---|
| `HumanBrain-Tumor-ST073021-02` | B189 | `ST073021-02` (GSM9497624) |
| `HumanBrain-Tumor-ST073021-04` | B189 | `ST073021-04` (GSM9497625) |
| `HumanBrain-Tumor-WU1220-01` | WU1220 | `WU1220-01` (GSM9497626) |
| `HumanBrain-Tumor-WU1220-3` | WU1220 | `WU1220-03` (GSM9497628) |
| `HumanBrain-Tumor-WU1221-Edge1` | WU1221 | `WU1221-Edge1` (GSM9497630) |
| `HumanBrain-Tumor-WU1221-Core2` | WU1221 | `WU1221-Core2` (GSM9497629) |

**Why it matters.** `docs/DATA_INVENTORY.md` A3 records that the cross-platform
composition check — Visium NNLS deconvolution vs Xenium single-cell fractions —
has a **RETRACTED** Pearson 0.997: after the classifier-arity fix the 4-vector
reads 0.980 and *the non-tumour subvector collapses 0.945 → 0.016,
indistinguishable from no agreement* (`results/revalidation.md` §4). That verdict
rests on **one** tissue pair. Here are **six**, on three donors, from one
protocol — enough to say whether the collapse is the truth or an artefact of a
single pair, and the only honest way to run it is to **declare the threshold
before downloading the Xenium half**, which is possible precisely because the two
halves are separate accessions.

**What is NOT verified, and it is the load-bearing thing.** I confirmed the
section labels and patient ids match. I did **not** confirm whether a matched pair
is the *same physical section* or *adjacent sections of the same block* — neither
GEO record I read states it, and the two cases are not equivalent for a
spot-to-cell comparison (adjacent sections differ by one cell diameter of
tissue). The bioRxiv preprint's methods would settle it. **Marked UNVERIFIED;
check before this is written up as a same-section comparison.** Note also that the
GSE318562 cohort is mixed fresh-frozen/FFPE while the Xenium sections are all
FFPE, so the chemistry of each specific pair needs checking too.

---

### 2.4 — #4 GSE339006: Xenium Prime 5K on the exact tissue the TMZ anchor came from

- **Accession:** GEO **GSE339006**, released **2026-10-01** — eight days before
  this pass, and two months after `DATA_HUNT.md`. Siblings: **GSE339010**
  (scRNA-seq), **GSE338996** (Multiome scATAC+RNA), **GSE338999** (not inspected).
- **Design, verbatim:** *"Paired sonicated and non-sonicated FFPE peritumoral brain
  regions from 5 recurrent glioblastoma patients were profiled by Xenium to assess
  spatial transcriptomic changes after LIPU/MB-mediated BBB opening."* Intraoperative
  low-intensity pulsed ultrasound + microbubbles, biopsy **45 minutes after
  sonication**, phase I/II trial **NCT04528680**. 10 sections.
- **Assay / panel — measured.** `features.tsv.gz` (64,216 B, identical across
  samples) has **9,681 rows: 5,101 `Gene Expression`**, 3,291 Deprecated Codeword,
  619 Unassigned Codeword, 609 Negative Control Codeword, 40 Negative Control
  Probe, 21 Genomic Control. That is a **Xenium Prime 5K-class panel**, and it is
  the only panel found in this pass that carries all four axes at once:

| axis | genes verified present |
|---|---|
| vessel / BBB endothelium | **PECAM1, FLT1, PLVAP, CLDN5** |
| hypoxia | **HIF1A, CA9, VEGFA, SLC2A1, NDRG1, ADM** |
| efflux | **ABCB1, ABCG2** |
| GBM / repair | **MGMT**, MKI67, SOX2, EGFR, GFAP, AQP4, POSTN, LOX |

  Verified absent: VWF, COL1A1. **MGMT is notable**: `landscape_cancer_census.md`
  records MGMT as absent from both stock 10x GBM Xenium panels and found only on
  one bespoke 366-gene GEO panel. Here it is on a 5K panel with the hypoxia axis.
- **Size:** `GSE339006_RAW.tar` = 42,638,120,960 B (**40 GB**) nominal. Real cost
  per sample: `matrix.mtx.gz` 4.5–25 MB + `features.tsv.gz` 64 KB +
  `barcodes.tsv.gz` 35–360 KB + `cells.parquet.gz` 0.13–0.67 MB +
  `cell_feature_matrix.h5` 1.4–7.3 MB. **Measured from `filelist.txt`: all 70
  non-morphology, non-transcript files across the 10 sections total
  400,100,370 B = 381.6 MiB.** (Adding the `transcripts.parquet.gz` files takes it
  to 4,071,360,416 B = 3.8 GiB; they are not needed for cell-resolved geometry.)
  The 40 GB nominal is `morphology.ome.tif.gz` (1.2–3.8 GB each).
- **Format advantage over #2:** this deposit ships the **MEX triple**
  (`matrix.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz`) **plus**
  `cells.parquet.gz`, so it is closer to what `singlecell_spatial.load()` already
  reads than the zarr deposits are.
- **Licence / access:** GEO, open, direct download, no DUA. Human trial tissue, but
  deposited in the open tier — no dbGaP, no DAC.

**The commensurability point, which is the real reason this ranks here.** The
tissue is **peritumoral brain, not tumour core** — a different compartment from
everything else in `data_cache/`, and that must travel with any use of it. But it
is *the same compartment as the project's highest-confidence partition anchor*:
Portnow 2009 (PMID 19861433), the TMZ brain-ECF:plasma AUC **0.178** that
`pk_anchors.py` carries at `high` confidence, was measured by microdialysis in
**human peritumoral brain**. Vivum has no tissue of that compartment at all. Five
patients' worth, at 5K plex, with CLDN5/PLVAP/ABCB1 in the panel, is the first
chance to reconstruct the tissue the anchor was measured in rather than applying a
peritumoral anchor to a tumour-core reconstruction.

**What it is NOT, stated plainly.** It measures **transcripts, not drug
concentration**. The sonicated/non-sonicated pairing is a BBB-permeability
perturbation — tempting, because `CLAIM_INVENTORY.md` §6 item 3 is
*"Flipping `transvascular_exchange` to True and re-validating"* — but **no drug
was measured and this dataset cannot score the transport PDE.** It can inform
the structural `bbb_penetrance` prior and the endothelial/vessel geometry term.
Treating a transcriptional response to BBB opening as a transvascular-permeability
measurement would be precisely the modality conflation
`external_anchor_survey.md` §0.1 catalogues three live instances of.

---

### 2.5 — #5 Evans 2008: the one literature find that could falsify something, and it is free

> Evans SM, Jenkins KW, Jenkins WT, Dilling T, Judy KD, Schrlau A, Judkins A,
> Hahn SM, Koch CJ. "Imaging and analytical methods as applied to the evaluation
> of vasculature and hypoxia in human brain tumors." *Radiat Res*
> 2008;170(6):677-90. **PMID 19138031**, PMCID **PMC2629387**,
> DOI [10.1667/RR1207.1](https://doi.org/10.1667/RR1207.1).

- **What was measured:** EF5 (a 2-nitroimidazole hypoxia marker) given IV at
  **21 mg/kg ~24 h before surgery**, then *"EF5 binding as a function of radial
  distance from individual blood vessels"*, converted to tissue pO₂.
- **Tissue:** **six human brain tumours** — one hemangiopericytoma, one WHO Grade 2
  and one WHO Grade 3 glial tumour, plus glioblastomas. **Patient tissue, not
  xenograft.** This is the right tissue and the right species, which almost no
  transport-adjacent anchor in this repo is.
- **The findings, quoted:** in the three low-grade/non-glial tumours *"all tissue
  pO(2) values calculated by EF5 binding were >20 mmHg"* and *"EF5 binding
  gradients (measured as a function of distance from each observed vessel) were
  low, with small positive and negative values averaging close to zero."* In the
  glioblastomas, *"much lower tissue oxygen levels were found, including near some
  vessels"*, and gradients were larger *"but positive and negative values again
  averaged to near zero."* The conclusion: ***"neither distance to the nearest
  blood vessel nor distance from each observed blood vessel provide reliable
  estimates of tissue pO(2)."***
- **Why it ranks above three downloadable datasets.** `external_anchor_survey.md`
  §6 item 2 says the project needs *"any external anchor that scores the solver on
  the strength of a solve"*, and that there is currently one and it fails. This is
  not that — it does not score the drug solver. It is something rarer and cheaper:
  **a measured, human, in-tumour result that is prima facie adverse to the
  distance-from-vessel field paradigm itself**, published, free, and sitting in the
  compartment this engine computes. The project has an oxygen field
  (`results/oxygen_calibration_landscape.md`,
  `results/oxygen_sink_form_data.json`) already demoted once
  (`results/oxygen_demotion.md`); this says the demotion was, if anything, not
  cautious enough, and offers a mechanism (in-flowing vs out-flowing blood) for
  *why* radial distance underdetermines pO₂.
- **The limit, honestly.** It is about **oxygen, not drug**. O₂ is consumed orders
  of magnitude faster than a cytotoxic is cleared and is carried by flowing blood
  in both directions, so "radial distance fails for pO₂" does **not** transfer to
  "radial distance fails for C". Primeau 2005 measured a clean exponential drug
  gradient in the same geometry. The honest reading is that the two coexist and the
  engine should not quietly borrow the drug result's credibility for its oxygen
  field.
- **Extractability: UNVERIFIED.** The EuropePMC `fullTextXML` endpoint for
  PMC2629387 returned **HTTP 500** on 2026-10-09, so whether per-distance numbers
  sit in a named table or only in figures was not established. PMC2629387 is an
  open PMC record and the figures are reachable; this needs one more fetch attempt
  by a different route before anyone claims the numbers are or are not extractable.

---

### 2.6 — #6 GSE290724: a Visium HD human GBM section, with a caveat that caps its value

- **Accession:** GEO **GSE290724**, *"Neurodevelopmental hijacking of
  oligodendrocyte lineage programs drives glioblastoma infiltration [Spatial
  transcriptomics]"*. 2 samples: GSM8820300 (`6425_A`), GSM8820301 (`6425_B`).
  Released 2025-04-01.
- **Platform — verified by inspection, because GEO's text index does not say it.**
  A `db=gds` query for `(glioblastoma OR glioma) AND "Visium HD"` returns **0
  hits**; this series was found only via the unspaced token `VisiumHD`. I
  downloaded `GSM8820300_6425_A_barcodes.tsv.gz` and read it: barcodes are
  **`s_002um_02058_00301-1`** — the 10x **Visium HD 2 µm square-bin** format.
  `features.tsv.gz` has **18,085 rows** (filtered whole transcriptome).
- **Species + disease:** *Homo sapiens*, glioblastoma. The study's design is
  **regional**: *"For each patient, we profile the peritumoral region, tumor edge,
  and tumor core"*. Only 2 GSMs are in this series, both from donor `6425`.
- **Size:** `GSE290724_RAW.tar` = 264,427,520 B (**252 MB**).
- **Licence / access:** GEO, open, direct download, no DUA.
- **The caveat, and it is a real one.** The deposit is **MEX only — there is no
  `spatial/` folder, no `tissue_positions`, and no `scalefactors_json.json`.**
  Micrometre coordinates are recoverable, because a Visium HD barcode encodes its
  own bin row and column on a fixed 2 µm grid, so `(col × 2 µm, row × 2 µm)` is an
  exact position with no image scale needed — arguably cleaner than Visium v1,
  where `umpp` differs per sample. But there is **no H&E registration and no
  in-tissue mask**, so tissue extent would have to be inferred from the occupied
  bins. 2 µm bins are also far below the ~55 µm scale the reconstruction pipeline
  is built around and would need binning up.
- **One sentence on the claim:** a fourth platform at a resolution that makes
  `results/dx_sensitivity_steady_data.json` and `CLAIM_INVENTORY.md` §6 item 4
  (*"Resolving dx"* — `L`/dx = 2.4 for the shortest-`L` drugs) answerable from
  measured tissue rather than from grid refinement alone.

### 2.7 — #7 GSE302502: the best-instrumented panel found, on the wrong shape of tissue

- **Accession:** GEO **GSE302502**, *"Spatial immune profiling defines a subset of
  human gliomas with functional tertiary lymphoid structures [Xenium]"*. Released
  2025-09-01. 12 TLS+ patients drawn from a 642-glioma Frankfurt University Cancer
  Center cohort, plus whole-slide multiplex IF and high-plex IF **on consecutive
  sections**, plus bulk RNA-seq on 64 gliomas.
- **Panel — measured.** Downloaded
  `GSM9105216_PAT-377_TLS_only_cell_feature_matrix_features.tsv.gz` (4,251 B):
  541 rows, **468 `Gene Expression`** + 41 Negative Control Codeword + 20 Negative
  Control Probe + 12 Unassigned Codeword. Verified present:

| axis | genes |
|---|---|
| vessel / endothelium | **PECAM1 (CD31), CD34, VWF, FLT1, PLVAP, CLDN5** |
| **mural / pericyte** | **ACTA2, PDGFRB, RGS5** |
| hypoxia | **VEGFA, HIF1A, CA9, SLC2A1** |
| ECM / efflux | **COL1A1, ABCG2** |
| GBM / glia | MKI67, GFAP, AQP4, SOX2, EGFR |

  Verified absent: ABCB1, MGMT, OLIG2. **This is the most completely instrumented
  panel found in this pass for transport purposes** — six endothelial markers, a
  mural-cell layer, the full hypoxia axis, and COL1A1 for the ECM hindrance term.
- **The one genuinely distinctive thing it offers.**
  `external_anchor_survey.md` §0.1 records that this repo's two vascular anchors
  **are not the same assay** — Marien (PMID 27583442) counts **CD31** microvessels,
  Tamma (PMID 36611806) counts **CD34** — and that the mismatch travelled
  unnoticed into adoption. **This panel carries PECAM1 and CD34 on the same cells
  in the same tissue.** It is the only dataset found here that could put a number
  on how far apart those two conventions actually are, which is a prerequisite for
  using either anchor honestly. That is a small, cheap, falsifiable question, and
  it is worth more than the dataset's size suggests.
- **Size: 29.3 MiB** for all 84 non-morphology, non-transcript files (measured from
  `filelist.txt`), against **89,612,175,360 B (83 GB)** nominal. Format is the
  **MEX triple** plus `cells.parquet.gz` — close to existing loaders.
- **The disqualifying caveat, stated plainly.** Every deposited matrix is named
  `*_TLS_only_*`. **The cell matrices appear to be cropped to the tertiary-lymphoid-
  structure regions, not the whole section** — a `cells.parquet.gz` of 317 KB is
  consistent with a few tens of thousands of cells, not a full section. If so there
  is **no full-section geometry here at all**, and a transport field cannot be
  solved on a crop with invented boundaries. **Marked UNVERIFIED — inferred from
  the filename convention and file size, not from a statement or from reading the
  parquet.** Verify before spending anything beyond the 4 KB features file.
- **And the selection bias on top of it:** the 12 sequenced patients were chosen
  *because* they are TLS+, i.e. the 15% tail of a 642-tumour cohort. Fine for a
  perivascular-ECM question, wrong for anything that must be representative.

**One sentence on the claim:** worth 29 MiB for the CD31-vs-CD34 convention
measurement alone, even if the TLS crop makes it useless for a transport solve.

### 2.8 — #8 GSE305393: a fifth platform, on the wrong geometry

*"Microbial signals in primary and metastatic brain tumors [CosMx]"*. **CosMx SMI
6K Discovery panel** + custom bacterial probes, two **tissue microarrays**: glioma
n=18, brain metastasis n=9, **39 fields of view per TMA**. 4.4 GB, and it ships a
`README.txt` — rare among these deposits. Adds a fifth platform and a 6K panel.
**The limit is geometry:** TMA cores are ~1 mm, and a transport simulator needs a
field spanning many inter-vessel distances with intact boundary conditions.
FOV-sized tissue is the wrong shape for this engine regardless of plex. Panel
membership **UNVERIFIED** (not inspected).

### 2.9 — #9, #12 The deconvolution reference: the 51-cell problem, and a 118 MB fix

This is **brief item 3** — the component that genuinely *is* fitted to data, and
so the only place "training on public databases" has an honest meaning here.

**The defect, restated from the repo's own record.** `docs/real_data/scref_gbm.md`
gives the reference's composition: Immune 1847, Neoplastic 1091, OPC 406,
Astrocyte 88, Oligodendrocyte 85, **Vascular 51**, Neuron 21. **The endothelial
signature column that sets the vessel source term in the transport PDE is built
from fifty-one cells**, from 4 donors, on SmartSeq2 — a different chemistry from
the 3'-biased UMI Visium data it is unmixing.

**#9 — Wälchli et al. 2024. Verified myself, every number below.**

> Wälchli T, Ghobrial M, Schwab M, et al. "Single-cell atlas of the human brain
> vasculature across development, adulthood and disease." *Nature* 2024
> Aug;632(8025):603-613. **PMID 38987604**, DOI
> [10.1038/s41586-024-07493-y](https://doi.org/10.1038/s41586-024-07493-y).

- **Accession:** GEO **GSE256493** (SuperSeries; GEO record confirms
  `!Series_pubmed_id = 38987604`).
- **Species + disease:** human. Fetal brain and periphery, adult temporal-lobe
  controls, and five vascular-dependent pathologies **including glioblastoma** —
  confirmed by the supplementary file names on the GEO FTP listing, which include
  `GSE256493_Glioblastoma_sorted_endothelial_cells_seurat_object.rds.gz` (2.8 G)
  and `GSE256493_Glioblastoma_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz`
  (900 M).
- **Assay:** 10x Chromium scRNA-seq, samples either FACS-sorted CD31+/CD45−
  (endothelial) or unsorted (endothelial + perivascular).
- **The cell counts, verified by me against the UCSC Cell Browser
  `dataset.json`:** `brain-vasc/facs/gbm` **`"sampleCount": 49999`**;
  `brain-vasc/unsorted/gbm` **`"sampleCount": 22297`**, shortLabel
  *"GBM Endothelial and Perivascular Cells"*. The metadata fields include
  `EC_clusters` and `Patient`.
- **The cheap route, verified by HTTP HEAD:**
  `https://cells.ucsc.edu/brain-vasc/facs/gbm/exprMatrix.tsv.gz` → HTTP 200,
  **Content-Length 123,811,242 (118 MB)**; `meta.tsv` → HTTP 200,
  **3,735,924 B**. The `dataset.json` reports the same size with md5 prefix
  `7f1d6317cf`.
- **Licence / access:** open direct anonymous HTTP from both GEO and UCSC. No
  registration, no DUA. (The scouting strand reports the authors' CELLxGENE
  collection currently returns `PRIVATE` with 0 datasets — **UNVERIFIED by me**,
  but irrelevant, because GEO and UCSC both serve it.)
- **Malignant cell states: no.** It is a vasculature atlas.

**For 118 MB, the endothelial column goes from 51 cells to 49,999 GBM-specific
endothelial cells.** That is a ~980× increase in the worst-powered column in the
signature matrix, on the quantity that sets the vessel source term, and it is
**independent of the current reference** (unlike GBmap — see below).

**#12 — Siletti et al. 2023 Human Brain Cell Atlas**, CELLxGENE collection
`283d65eb-dd53-496d-adb7-7570c7caa443`, *Science* 2023, DOI
[10.1126/science.add7046](https://doi.org/10.1126/science.add7046). snRNA-seq,
human **normal** brain. The `Supercluster: Vascular` object is reported as
**9,932 cells in a 47 MB H5AD**, with `endothelial cell`, `pericyte`,
`vascular associated smooth muscle cell` and `fibroblast` labels; a `Fibroblast`
supercluster of 9,156 cells; and `All non-neuronal cells` at 888,263 cells / 4.3
GB. **These counts and sizes are from the scouting strand's reading of the
CELLxGENE curation API and are UNVERIFIED by me** (I verified Wälchli's numbers
and spent the remaining budget there). Its distinct value: it is the only
candidate that **also repairs the Astrocyte-88, Oligodendrocyte-85 and Neuron-21
columns**, which matter because NNLS is a joint fit and noise in those columns
leaks into the vascular estimate through the simplex constraint. Normal brain, so
it will not carry GBM-specific endothelial phenotype.

**GBmap (CELLxGENE `999f2a15-3d7e-440b-96ae-2c806799c08c`, Ruiz-Moreno et al.,
*Neuro-Oncology* 2025, DOI 10.1093/neuonc/noaf113) — take it for the malignant
states, not the vessels, and never call it independent.** The scouting strand
reports core GBmap at 338,564 cells / 110 donors / 16 studies (7.6 GB H5AD) with
`endothelial cell`, `mural cell` and `malignant cell` types, Neftel-style
AC/MES/OPC/NPC labels, and endothelial cells split into tip / capillary /
arterial-like. **Two caveats, both important and both UNVERIFIED-by-me but
structurally sound:** (i) the vascular share is ~1% of cells, i.e. roughly 3,400
endothelial + mural *combined* — far below Wälchli, for 65× the download;
(ii) **core GBmap's constituent studies include GSE84465 — the reference this
project already uses — and GSE131928, and extended GBmap adds GSE182109.** So a
"GBmap vs Darmanis" comparison is a **superset against its own subset, not two
independent references**, and must never be presented as corroboration. Only
Wälchli and Siletti are genuinely independent of the current reference.

**The honest read on whether this changes a deconvolution result.** It will move
the numbers, and **the most likely direction is that the vascular fraction goes
DOWN, which is why the direction must be written down before the run.** With a
single `Vascular` column built from 51 cells, pericyte/vSMC/fibroblast signal in a
55 µm Visium spot has nowhere to go *except* the vascular column. Introduce a
separate mural column — which both Wälchli-unsorted and GBmap have — and NNLS will
reallocate it. If the vessel source term is currently calibrated against a
vascular fraction that was silently absorbing pericytes, **this change inverts
that calibration.** Given `results/accuracy_gap_decomposition.md` and the
2026-08-22 finding that a headline verdict inverts under a term shipped OFF, this
is the familiar shape of defect here, and it argues for running the swap as a
declared sensitivity arm rather than shipping it as an improvement.

Two further cautions that should travel with any adoption:
- **Do not change reference size and sequencing chemistry in the same step.**
  Darmanis is SmartSeq2 (full-length, no UMIs); Wälchli, Siletti and GBmap are all
  10x (3'-biased, UMI). A 10x reference is better *matched* to Visium independent
  of *n*, so the platform swap and the *n* increase are confounded unless run as
  separate arms. Any single-step swap produces a number nobody can attribute.
- **FACS sorting makes a clean basis vector and a biased prior.** CD31+/CD45−
  sorting is ideal for a signature column, but by construction excludes the
  endothelial–immune doublets that genuinely occupy Visium spots. If the goal is
  in-situ-faithful composition rather than a clean basis, the 22,297-cell
  *unsorted* GBM object is the right file.
- **And the boundary that does not move:** a better reference reduces bias and
  variance in the fitted fractions. **It does not make 55 µm spot-level unmixing
  identifiable**, and it does not license absolute per-spot claims. That is the
  same line as the 2026-08-13 settlement — relative comparison defensible,
  absolute numbers not — and a 980× bigger reference does not cross it.

### 2.10 — #10, #11 Transport anchors: absolute free concentrations in human glioma, and a 25–40× contradiction

These came from the literature strand. **I verified the identity, journal, design
and principal finding of #10's two papers from their primary PubMed records; I did
not open their supplementary tables.** #11's numbers are the strand's reading and
are **UNVERIFIED by me** — flagged inline.

**#10a — absolute unbound drug in human high-grade glioma, two compounds, three
orders apart.**

> Johnson KC, Tien AC, Jiang J, et al. "Single-nucleus transcriptomics,
> pharmacokinetics, and pharmacodynamics of CDK4/6 and mTOR inhibition in a Phase
> 0/1 trial of recurrent high-grade glioma." *Neuro Oncol* 2026
> Mar 1;28(3):659-671. **PMID 41206763**, DOI
> [10.1093/neuonc/noaf257](https://doi.org/10.1093/neuonc/noaf257), PMC13070491.
> NCT03834740.

Verified from the abstract: 24 patients with recurrent HGG; *"Median unbound
ribociclib concentrations in gadolinium non-enhancing tumor regions were
significantly above the biochemical IC50 for CDK4/6 inhibition at 400 and 600 mg
QD doses. **Unbound everolimus concentrations were undetectable (<0.1 nM) in tumor
regions across all dose levels.**"* The compartment is **unbound (free)** —
the quantity the PDE computes — in **human glioma**, separated into
**non-enhancing vs contrast-enhancing** regions.

**What it could falsify:** two drugs, the same patients, the same vasculature, and
free interstitial concentrations separated by **≥3 orders of magnitude**. A solver
whose free-C prediction is dominated by vascular delivery geometry rather than by
per-compound binding and efflux cannot produce both numbers. That is a sharp,
cheap, per-compound magnitude test — the surviving check form named in
`docs/real_data/pk_anchors.md` — and it is in the right tissue and species.
Per-patient values are reported by the strand to be in Supplementary Table S3;
**read it rather than quoting via this file.**

**#10b — a free partition ratio of ~1 in human glioma, behind both barrier states.**

> Schiff D, Ye X, Li J, et al. "Phase I/II and window-of-opportunity study of
> pamiparib and metronomic temozolomide for recurrent isocitrate dehydrogenase
> mutant gliomas." *Neuro Oncol* 2026 Feb 1;28(2):485-494. **PMID 41099363**, DOI
> [10.1093/neuonc/noaf246](https://doi.org/10.1093/neuonc/noaf246), PMC12979031.

Verified from the abstract: *"In non-enhancing and enhancing tumor, pamiparib
exhibited an **unbound tumor/plasma ratio of 0.92 and 0.98**, respectively."*
66 subjects enrolled; the intratumoral PK is an exploratory cohort.

**This is a partition ratio, not a gradient — it belongs to the `bbb_penetrance`
check and must never be pooled with the µm anchors**, per the standing rule in
`docs/real_data/pk_anchors.md`. Its value: a free partition of **~1** in
**non-enhancing** human glioma is at the very top of the 0.013–0.178 range the
repo's three existing partition anchors span, and it is IDH-mutant glioma rather
than GBM. Also note the trial co-dosed metronomic TMZ at 20 mg/day and, per the
strand, **reports no TMZ tissue concentration** — so it adds no TMZ anchor. The
strand flags that the paper writes *"µmol/L (or µmol/kg)"* without committing to a
volume-vs-mass basis; **resolve that before use**, because it is exactly the
convention-not-travelling-with-the-number failure catalogued in
`external_anchor_survey.md` §0.1.

**#11 — a measured hydraulic conductivity, an IFP boundary-layer length, and a
contradiction worth pre-registering. ALL NUMBERS IN THIS ITEM ARE UNVERIFIED BY
ME.**

- Salavati H et al., *Bioeng Transl Med* 2024;9(2):e10617, **PMID 38435818**,
  PMC10905546 — modified Ussing chambers on fresh human surgical specimens.
  Reported `K` for human peritoneal/pancreatic tumour of order
  **4–7 × 10⁻¹⁵ m²/(Pa·s)**, with normal peritoneum ~**2.8 × 10⁻¹³** (~40× higher),
  up to 4× heterogeneity within one tumour, and an inverse correlation with
  collagen fibre density. **Hard limitation: no brain or glioma tissue.**
- Elmghirbi R et al., *Magn Reson Med* 2018;80(5):2040-2052, **PMID 29524243**,
  PMC6107371 — wick-in-needle in orthotopic U251 rat cerebral tumour. Reported
  TIFP **6.0 ± 3.7 mmHg** (range 1.6–14.9), and — rarest and most useful — an
  **IFP boundary-layer length of 1.0–2.5 mm, mean 1.66 ± 0.48 mm**, over which IFP
  rises from atmospheric to plateau.

**The two things to take from this pair:**
1. **A pre-registrable contradiction.** The strand's unit conversion puts
   Elmghirbi's rat brain-tumour `K` ~25–40× above Salavati's human
   peritoneal/pancreatic values. **A single `K` in Vivum cannot satisfy both.** The
   honest procedure is to declare which one the solver claims *before* looking and
   let the other break it. (The conversion is the strand's, not the papers' —
   re-derive it.)
2. **A shape constraint independent of absolute value.** The 1.66 ± 0.48 mm IFP
   transition length constrains the *shape* of the pressure field without
   committing to its magnitude, which makes it unusually robust for a solver whose
   absolutes are known to be poorly constrained.
   **Caveat: it is rat, and the repo's one human intracranial IFP anchor (Boucher
   1997, 2.0 ± 2.5 mmHg) is ~3× lower than Elmghirbi's 6.0 ± 3.7.** Those two are
   not reconcilable by species alone and should not be averaged.

---

## 3. Verified negatives — worth as much as the finds, because they stop a re-hunt

1. **METASPACE still contains no glioma drug-distribution data.** Re-checked live
   on 2026-10-09 against the GraphQL endpoint. `countDatasets` for
   `simpleQuery: "temozolomide"` → **0**; `"glioblastoma"` → **0**;
   `"brain tumor"` → **0**; `"glioma"` → **1**. `DATA_HUNT.md` §1 reported this as
   `[REPORTED]` from two independent angles; it is now directly verified, 14 months
   on, and **unchanged**.
2. **The one new human glioma MSI dataset is not a drug map.** `glioma_tma_c`,
   dataset id `2025-06-03_11h25m11s`, *Homo sapiens* / Brain / Diseased, uploaded
   2025-06-03 (i.e. new since the last hunt). Queried its full record: MALDI,
   Orbitrap, 140,000 resolving power at m/z 200, **50 µm pixel size**,
   paraformaldehyde-fixed, and annotated against **HMDB v4 only** — endogenous
   metabolites, no drug database. It is a TMA. Not an anchor.
3. **No Visium HD glioma dataset is discoverable by that name in GEO.** The
   `"Visium HD"` text query returns 0 glioma hits; #2.6 exists but is not indexed
   as Visium HD. Anyone re-searching this later should query the platform's barcode
   format or the unspaced token, not the product name.
4. **No deposited human GBM light-sheet / cleared-tissue vessel geometry dataset
   was found.** Searched; the human-GBM hits are confocal (§4) and mouse for the
   deposited light-sheet data.

### 3.1 Categories the literature strand found genuinely EMPTY

Reported by the literature strand, each with the search it ran.
**UNVERIFIED by me** — but an empty category is a finding, and these say where
*not* to spend the next pass:

- **Human brain-tumour IFP in mmHg: EMPTY beyond Boucher 1997.** 14 PubMed records
  across all years; the only primary human intracranial measurement is the one
  this repo already holds. Everything since is preclinical, non-brain, or a
  DCE-MRI inversion.
- **Drug concentration vs distance from vessels, published after mid-2025:
  effectively EMPTY.** Attacked from five directions (MALDI/DESI + glioma + drug;
  MSI + penetration/gradient + vessel/CD31; MSI + heme co-registration, all years;
  immunofluorescence + vessel distance + drug; EuropePMC incl. preprints) — the
  largest yield was 3 records and none is a drug-vs-vessel profile. **This is the
  same conclusion `external_anchor_survey.md` §5 reached, re-reached independently
  14 months later.**
- **A temozolomide perivascular penetration profile in µm: EMPTY at any date.**
  34 records, nothing spatial at the vessel scale. This is
  `alpha_identifiability.md` §4 item **(b)** — a perivascular profile for a
  weakly-bound, deeply-penetrating drug — and **it does not exist in the published
  literature.** The engine's single calibration constant α is solved in closed form
  to an anchor (Primeau 2005) that is both in-sample and from the opposite chemical
  class, and **no published measurement can currently break that circularity.**
- **CED concentration-vs-distance in human glioma: EMPTY for concentration.**
  Human CED trials since 2020 report volume of distribution, not
  concentration-vs-distance.
- **PET quantification of absolute drug concentration in human glioma: EMPTY.**
  58 records surfaced only tracer development and diagnostic imaging.
- **Ex vivo human brain slice / organotypic drug-penetration gradients: EMPTY.**
  6 records, none measures a gradient.
- **The one live lead that is not a dataset:** the strand reports that Rendo et al.
  2025 (*Sci Transl Med* 17(786):eadn6274, **PMID 39970230**) ran MALDI-MSI of
  navtemadlin **co-imaged with heme b as a vessel/blood proxy in 13 human GBM
  patients**, and published it only as a region-averaged bar chart, not
  distance-binned, and did not deposit it. If that is right, the pixel-level
  drug-vs-vessel data for human brain tumour **exists and is unpublished**, and the
  only route to it is asking the authors. **UNVERIFIED by me**; worth one email
  before anyone concludes the category is unreachable.

### 3.2 Two numbers that will show up in searches and must NOT be used

Both flagged by the literature strand; recorded here precisely so they do not
propagate into this repo, which is the documented failure mode
(`CLAIM_INVENTORY.md`, the claim-surface-drift thread).

1. **Huber et al. 2014, *Anal Chem* 86(21):10568-75, PMID 25263480** — a web
   search surfaces a claim that it reports erlotinib declining over ~100 µm from
   the vessel and afatinib plateauing ~20 µm then declining to ~140 µm. The strand
   traced those µm figures to a **ResearchGate summary page, not the paper**; the
   primary record reports only a **region-based** association between drug level
   and vascularisation. **Mark UNVERIFIED; do not use.** Even if real, the tissue
   is subcutaneous xenograft, not brain. The paper's real contribution is the
   *method* — co-registered MALDI + CD31 IHC + digital image analysis — which is
   the correct recipe for the anchor this project cannot find, never run on brain
   or against a distance axis.
2. **Geiger et al. 2021, *Front Neurosci* 15:678978, PMID 34276289** — reports
   *"a diffusion coefficient of approximately 0.12 cm²·s⁻¹"* for ketamine
   diffusing into brain tissue. **That is physically impossible** — roughly 10⁴×
   free-water self-diffusion, a gas-phase magnitude. Treat as a published unit
   error. It is flagged because it will surface in any search for `D` in brain
   tissue, and this repo's `D_eff` work would be the thing that picked it up.

---

## 4. Checked and rejected — so nobody re-chases it

> Cribaro GP et al. "Three-dimensional vascular microenvironment landscape in
> human glioblastoma." *Acta Neuropathol Commun* 2021. **PMID 33579378**,
> PMCID **PMC7879533**, DOI
> [10.1186/s40478-020-01115-0](https://doi.org/10.1186/s40478-020-01115-0),
> licence **CC BY**.

On its title this is the ideal geometry anchor: 3D vessel architecture in human
GBM, which is what `CLAIM_INVENTORY.md` §6 item 5 (*"Re-deriving `supply` for the
discrete-vessel referent"*) needs, and which would improve on the repo's two
existing vascular anchors that `external_anchor_survey.md` §0.1 shows **are not
even the same assay** (Marien CD31 vs Tamma CD34). I pulled the full text and read
it. **It does not work, for three separate reasons:**

1. **The numbers are not in the text.** Vessel diameter, vessel volume, volumetric
   surface area and vascular volume fraction are all described as measured, and
   reported **only in figures** (Fig. 2c–d, 2a–b, 8b). Diameters are given as
   qualitative bands — "small, medium and large" — with the only numeric threshold
   being *"Small vessels ≤ 5 μm diameter insufficient to accommodate blood cell
   passage."* Vessel density, inter-vessel distance and tortuosity are **not
   reported at all**.
2. **The dispersion is SEM, not SD.** Stated as *"Data are represented as mean ±
   SEM"*. An anchor needs a population spread; SEM divided by an *n* that varies
   per analysis (50 data points for vessel characteristics, 24 for CD31) is not
   one.
3. **It is not deposited.** There is no formal data availability statement; the
   abstract says only *"Images and data sets are available on demand as a resource
   for public access."* That is an email, not a download.

Also worth recording: it would have needed the **3D-vs-2D** caution anyway. It
measures 60 µm-thick serial sections by confocal with **COL-IV** (basement
membrane) and **CD31** (endothelium) — a 3D volume quantity on a third marker
convention, which is the exact error class that put a healthy-cortex **3D EM
volume fraction** up against a **2D spatial-transcriptomics area fraction** in
`heldout_patient.py` (`external_anchor_survey.md` §0.1, PMID 38723085).

**Verdict: not usable without contacting the authors. Do not re-chase on the
title.**

---

## 5. Honest read — would any of this change a claim Vivum currently makes?

**Yes — three of them, and the biggest one is a parameter, not a dataset.**

**It would change a claim, and it is the most important item here: §2.0, a
measured free-interstitial loss rate.** `alpha_identifiability.md` §4 names two
missing things; Habashy 2025 (PMID 40495421) supplies item **(a)** — a
free-interstitial `k` from human-brain microdialysis with a 24-hour washout — and
`pk_reality.md` records `k` as having been **~300× too small**. A parameter that
the calibration currently cannot identify, and that the engine has been wrong
about by two and a half orders of magnitude, now has a measured value in the right
compartment and the right species. Nothing else in this file touches a number of
that consequence. And §3.1 says plainly that item **(b)** — the weakly-bound
perivascular profile that would break α's circularity — **still does not exist in
the literature at all.**

**It would change a claim: §2.3, the six matched pairs.** The retracted Pearson
0.997 and the non-tumour subvector collapsing to **0.016** is currently a
one-tissue result, and `docs/DATA_INVENTORY.md` A3 already says of it
*"Concordance here is an open question, not evidence."* Six pairs on three donors
can close that question in either direction. That is the only item in this file
that can change a **verdict** rather than widen a caveat — and it is the one most
able to make the project's position *worse*, which is why it is the one worth
running. It also has the rare property that it can be genuinely pre-registered:
the threshold can be fixed before the Xenium half is downloaded, because the two
halves are separate accessions. Given this project withdrew a headline claim on
2026-08-21 over a post-hoc substitution, that property is worth more than the data.

**It would change a claim: §2.2, measured CA9/VEGFA/ADM at single-cell
resolution.** `hypoxia_imputation_model.md` and `panel_gap_imputation_sweep.md`
exist because those transcripts are off the Xenium panels — and Xenium is the path
the flagship results run on. On these 15 sections they are measured. An imputation
whose target becomes measurable stops being an assumption and becomes a scoreable
step that can fail, which is the point. (§2.2 spells out the narrower, accurate
version of this: Visium already had these genes, so the gain is resolution, not
first measurement.)

**It could repair a known defect for 29 MiB: §2.7, CD31 and CD34 on one panel.**
The repo's two vascular anchors are different assays and that mismatch reached
adoption unnoticed (`external_anchor_survey.md` §0.1). GSE302502 carries both
markers on the same cells in the same tissue, which is the only way found here to
measure how far apart the two conventions are. This does not change a verdict; it
makes two existing anchors quotable without a silent assumption. Contingent on the
TLS-crop question in §2.7.

**It would weaken a claim, and should: §2.5, Evans 2008.** A measured human
in-tumour result saying radial distance from vessels does not reliably predict
pO₂ is adverse to the oxygen field, free, and already published. The oxygen field
was demoted once; this argues the demotion did not go far enough. Nobody has to
download anything to act on it.

**It would NOT change a claim, and saying so is the useful part:**

- **§2.1, the 22 Visium sections.** This is the single best *find* and it changes
  **no verdict**. It changes an *n*. Almost every defensible claim in
  `CLAIM_INVENTORY.md` is qualified by "*n* = 1 in tissue"; 22 sections from 9
  donors on one protocol replace that caveat with a measured spread. That is worth
  674 MB, but it is breadth, not a new kind of evidence — and
  `multisection_generalisation_2026-09-20.md` is the standing warning that more
  sections can just as easily reveal that a verdict was never stable.
- **§2.4, GSE339006.** Genuinely the best panel found and the right tissue for the
  TMZ anchor, and still **not a transport measurement**. It has no drug in it.
- **Nothing here closes the real gap, and the gap got better characterised rather
  than smaller.** `external_anchor_survey.md` §5 concluded that *nothing public
  measures a spatially resolved free-interstitial drug concentration field in
  tumour tissue this engine could reconstruct*, and that closing it *"may require
  generating data, not finding it."* **This pass confirms that conclusion from two
  independent directions rather than overturning it:** §3 item 1 (METASPACE still
  has zero temozolomide and zero glioblastoma datasets, now directly verified) and
  §3.1 (a five-angle literature sweep for drug-vs-vessel-distance since mid-2025
  returning nothing usable). §2.0 gets closer than anything before it — but
  microdialysis measures a **time course at a point**, not a field, so the
  *spatial* prediction remains unscored by any external measurement.
- **The structural imbalance, stated as plainly as it deserves.** Almost
  everything found this pass enriches the **input** side of the simulator —
  geometry, cell types, hypoxia, vessels, a 980× better endothelial reference.
  Only §2.0 and §2.10 touch the **output** side, and only as point values and
  ratios, never as a field. **A pass that improves the inputs much faster than the
  output checks makes the engine more elaborate and no more falsifiable.** On this
  project's own thesis — verification, not simulation — that is the wrong
  direction, and the fact that the input finds are the cheap ones and the output
  finds are the scarce ones is precisely why it keeps happening. The §2.3 held-out
  test and the §2.10 per-compound magnitude checks are the only items here that
  produce a verification result rather than a better input; they are what this pass
  is actually worth.

---

## 6. The wrong direction, named explicitly — as the brief asked

**If the answer to "what should we train on" comes back as a patient-outcome
dataset, that is the wrong answer for Vivum, and the reason is not squeamishness.**

The obvious "train on this" targets exist, are open, and are already on this
machine — `data_cache/hunt/tcga_treatment_response/`,
`data_cache/hunt/dce_mri_response/` (LUMIERE, RHUH-GBM, Burdenko, UCSD-PTGBM
clinical tables), `data_cache/hunt/pdx_encyclopedia/` with a
`TableS7_PDX_response.csv` carrying median survival per arm, HGCC's 262×85
response matrix, GDSC2. **Fitting a response predictor to any of them is out of
scope by standing prohibition**, and it would be a mistake even if it were
permitted, for three reasons:

1. **It is what the product's credibility is made of.** Vivum's claims survive
   scrutiny *because* it refuses this one. Cross-compound ranking was **withdrawn**
   on 2026-08-21 rather than defended; absolute kill fraction is never a
   defensible output. A response predictor re-acquires, in one step, every claim
   the project spent months giving up.
2. **The loop does not close in the data anyway.** `DATA_HUNT.md`'s headline
   finding, from 164 datasets: *"No open dataset anywhere pairs a measured spatial
   drug distribution with a measured per-specimen treatment response in the same
   tumour. Not one."* A predictor trained across that gap is fitting tissue state
   to outcome with the transport physics as decoration — and the simulator becomes
   unfalsifiable, because any error can be absorbed by the fitted layer.
3. **It is the known null.** `~/tme-sim`'s own Stage-0 kill test found the
   simulator does not beat a cheap baseline. That honest null is *why* the thesis
   is verification, not prediction. Training on outcomes would convert a clean
   negative result into an unauditable positive one.

The honest meaning of "training on public databases" for this product is narrow
and real: **the deconvolution reference** (§2.9) and **the panel-imputation model**
(§2.2, which §2.2 would let you partly retire). Those are components that already
learn from data. Outcome labels are not.

---

## 7. Open threads and what was not covered

- **What is verified by me vs. reported to me.** This pass ran three strands: I
  took the deposited-tissue strand myself, and two parallel strands covered the
  transport literature (§2.0, §2.10, §3.1, §3.2) and the reference atlases (§2.9).
  **Verified by me against a primary record:** everything in §1, §2.1–§2.8, §3,
  §4; the identity, journal, trial, design, compartment and principal finding of
  PMIDs 40495421 / 41206763 / 41099363 (PubMed `efetch`); Wälchli GSE256493's GEO
  record, its PMID 38987604 linkage, its GBM supplementary file names, and the
  UCSC cell counts (49,999 / 22,297) and file size (123,811,242 B) by
  `dataset.json` and HTTP HEAD; Evans 2008 PMID 19138031 by abstract; Cribaro by
  full text. **Reported to me and NOT independently verified:** every per-patient
  or per-table number in §2.0 and §2.10, all of §2.9's Siletti and GBmap counts,
  §3.1's search tallies, and §3.2's two warnings. Those are marked inline. **No
  number in a `reported-to-me` category may be quoted onward via this file** — read
  the primary source, which is the rule this project learned the hard way.
- **Tier-2 literature leads not written up here**, reported by the literature
  strand and worth a look by whoever picks this up: Rendo 2025 navtemadlin MSI
  (PMID 39970230, §3.1); Xiong 2024 (PMID 38287297) reporting that the tumour
  **margin** has the *lowest* extracellular diffusion coefficient — lower than core
  *and* than normal brain, which if true falsifies a uniform `D` *and* the usual
  assumed direction, and which the strand could not retrieve numbers for (Nano Lett
  paywalled, EuropePMC 500); Mayol Del Valle 2026 (PMID 41790275), whose abstract
  the strand says **contradicts its own Table 1** on tumour-vs-peritumoral
  direction for one of two drugs; and a PNAS correction (PMID 41564141) that must
  be read before using PMID 40901880.
- **Not inspected, and named so the gap is visible:** the CosMx 6K panel membership
  for GSE305393; GSE338999; whether GSE302502's matrices really are TLS-cropped
  (§2.7, the one thing that decides whether it is 29 MiB well spent); the Greenwald
  deposit beyond the one CODEX + two Visium sections already cached; HTAN v7 open
  tier beyond `data_cache/hunt/htan_open_tier/`; `data_cache/ivy_gap/`'s
  unextracted zip, which is a known-broken loader rather than a scouting item but
  is cheaper to fix than anything in this file.
- **The one thing to verify before §2.3 is quoted:** same physical section vs
  adjacent section for the matched pairs, from the bioRxiv methods. §2.3 says this
  too; it is repeated here because it is the kind of detail that gets dropped
  between a scouting file and a claim, which is this project's documented recurring
  defect (`results/CLAIM_INVENTORY.md`; the claim-surface-drift thread).

---

## 8. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

Unchanged from §0. No code, data loader, scorecard entry or existing dated record
was touched by this pass. The only file this pass created or modified is this one.

**One observation, recorded because another agent may be working in this tree.**
A sweep for files modified during this pass
(`find results archive tmesim docs app -newermt "2026-10-09 19:00"`) returned two
hits: this file, and **`docs/AUTH_CONTRACT.md`, mtime 2026-10-09 19:26:21**. This
pass did not open, read or write `docs/AUTH_CONTRACT.md`. Flagging it rather than
investigating it, since another owner may hold that file and this project's
fan-out discipline is exclusive file ownership.

**Nothing in this file is registered, pre-registered, or scored.** Every item is a
recommendation, and the two items worth acting on first (§2.0 and §2.3) are
specifically the ones whose value depends on declaring the threshold *before* the
data is in hand.
