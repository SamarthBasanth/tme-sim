# Landscape survey — functional/drug-response measurement platforms

Date: 2026-07-28
Scope: research-only survey, no code/data changes. Six parallel research passes (WebSearch/WebFetch/PubMed)
across PDO biobanks, explant/slice culture, PDX, microfluidic/organ-on-chip, ex vivo functional precision
medicine + intratumoral microdosing, and MALDI/DESI mass-spec imaging, followed by manual verification of
the highest-value leads (PMID/DOI cross-checks via NCBI eutils, direct GEO/Zenodo record inspection).

**Relationship to `results/DATA_HUNT.md`**: that document (2026-07-27) is the authoritative answer to "what
open *dataset* already exists." Its headline: exactly one open dataset anywhere pairs a measured spatial
drug field with a measured molecular response on the same tissue (Moncayo et al. 2026, HGSOC explants,
Zenodo 17610220 / GSE281519), and `results/moncayo_explant.md` already ran tme-sim's transport+kill layers
against it and got an honest null (the transport-only model predicts a flat field at 24 h bath-equilibrium;
the lysosomal-trapping mechanism gets the compound ordering right and the magnitude 175x wrong; the
model's predicted exposure does not predict the measured apoptosis signature). This document answers a
different, complementary question: since tme-sim is now scoped as **downstream of measurement**, which
existing *platforms* — not just already-deposited datasets — generate the pairing we need, or could be
asked to, and who should we approach. Two genuinely new leads surfaced in this pass that DATA_HUNT.md does
not contain (CIVO's `GSE235475`, and the Jonas Lab/Kibur Medical glioma-microdevice Zenodo deposit
`10.5281/zenodo.8226867`) — both flagged prominently below.

**Access-tier vocabulary used throughout**: `open` (download now, no account) · `registered` (free
account/DUA, self-service) · `controlled` (dbGaP/EGA-style Data Access Committee, non-commercial-only,
weeks–months) · `commercial` (fee-for-service, no research data-sharing path found) · `none` (not deposited
anywhere, figures/text only, or "available upon request" with no accession).

---

## 1. Headline answer

**No platform surveyed has an open, downloadable dataset that pairs spatial tumour structure with measured
drug response on the same specimen, in glioma.** That reconfirms DATA_HUNT.md's finding one layer up, at
the level of platforms rather than datasets: the capability to generate this pairing exists and is actively
being used in a handful of labs/companies, but in every GBM-relevant case the data is gated ("available
upon request," no accession) rather than deposited.

**The closest thing to a full, open, GBM-specific answer is not a dataset but a platform we had not
previously identified**: the Jonas Lab (Brigham & Women's)/Kibur Medical intratumoral microdevice (IMD),
first-in-human in glioma (Peruzzi et al., *Sci Transl Med* 2023, PMID 37672566, PMC10754230, DOI
[10.1126/scitranslmed.adi0069](https://doi.org/10.1126/scitranslmed.adi0069)). It implants up to 20
spatially segregated drug microdoses directly into a living glioma at the time of standard-of-care
resection, then measures **both** drug concentration (autofluorescence + MALDI-TOF at non-fluorescing
drugs) **and** biological response (cleaved caspase-3, phospho-H2AX IHC, spatial transcriptomics) at
matched coordinates on the *same* resected tissue, exploiting the distance-dependent concentration
gradient around each reservoir to correlate dose with response. A companion data deposit — **Zenodo
[10.5281/zenodo.8226867](https://doi.org/10.5281/zenodo.8226867), open, CC-BY-4.0, no registration** —
contains Nanostring transcriptomic data and MALDI metabolomic (glutathione/PEG) maps from the trial. This
is the single most important finding of this survey and is detailed in §2.6 and §4.

The second new finding is CIVO's own GEO deposit, **GSE235475** (open, no registration, verified directly
against the GEO record) — the first genuinely public spatial-exposure + spatial-response pairing this
project has located anywhere, in head-and-neck carcinoma, not glioma (§2.6, §4).

---

## 2. Platform-by-platform survey

### 2.1 Patient-derived organoid (PDO) biobanks and screening companies

| Platform | Spatial structure recorded at readout? | Readout granularity | Public data | Access tier |
|---|---|---|---|---|
| HUB (Hubrecht Organoid Technology) | No — bulk viability/luminescence screen correlated with DNA/RNA-seq | Bulk well-level | Molecular data in paper supplements only | commercial (screening service) |
| HCMI (NCI/CRUK/Sanger/Broad/Hubrecht, via ATCC + GDC) | No | N/A — sequencing/clinical, not functional | WGS/WES/RNA-seq/clinical open (masked) / controlled (germline) | open (sequencing) / controlled (germline); **drug-response data type not confirmed to exist in GDC despite ATCC marketing language** |
| Sanger/CRUK Tumour Organoid Biobank | Characterization only; functional screening stated as future intent | N/A yet | Characterization data via EGA / Cell Model Passports | open (characterization) / commercial (screening access) |
| GBM organoids (GBOs), Jacob et al./Song lab (JHU→UPenn) | Heterogeneity captured compositionally (scRNA-seq, clonal), not spatially registered to drug readout | Per-organoid/well, follow-up scRNA-seq on separate organoids | None found | commercial / collaboration-only |
| Xilis (MicroOrganoSphere) | No | Per-droplet, aggregated | None | commercial |
| SEngine Precision Medicine (PARIS test) | No | Bulk luminescence | None | commercial (CLIA test) |
| Cellesce | Per-organoid confocal morphometric imaging (not registered spatial-omics) | Per-organoid | None | commercial |
| Crown Bioscience (OrganoidXplore, PDXO) | No | Bulk IC50 | Model metadata browsable; IC50 datasets not freely downloadable | commercial |

Sources: van de Wetering et al., *Cell* 2015, PMID 25957691, DOI [10.1016/j.cell.2015.03.053](https://doi.org/10.1016/j.cell.2015.03.053); Jacob et al., *Cell* 2020, PMID 31883794, DOI [10.1016/j.cell.2019.11.036](https://doi.org/10.1016/j.cell.2019.11.036); GBO protocol, PMID 33169003, DOI [10.1038/s41596-020-0402-9](https://doi.org/10.1038/s41596-020-0402-9); Gobits et al., *JCO Precis Oncol* 2026, PMID 41538757, DOI [10.1200/PO-25-00501](https://doi.org/10.1200/PO-25-00501); Al-Aloosi et al., *Front Oncol* 2024, PMID 38239650; ATCC HCMI page (atcc.org/cell-products/collections-and-projects/human-cancer-models-initiative); GDC HCMI-CMDC (portal.gdc.cancer.gov/projects/HCMI-CMDC); Sanger Tumour Organoid Biobank (sanger.ac.uk/technology/tumour-organoid-biobank); Crown Bioscience organoid database (crownbio.com/databases/organoid-models).

**Bottom line**: no PDO platform records spatial tissue architecture at the moment of drug-response
measurement in a form usable for spatial validation. Best for informing cell-state/heterogeneity priors
(GBOs), not for the response axis directly.

### 2.2 Explant / precision-cut tumour slice culture (PCTS) — highest-priority category

This is the only category where the model system itself preserves native vasculature, stroma, and immune
architecture while being dosed with real drug ex vivo. A 2026 systematic review of GBM organotypic
brain-slice models found **zero** deposited public datasets across 26 studies screened (2010–2025), and
explicitly cites "substantial methodological heterogeneity" as the field's core limitation (Petralia,
D'Amico et al., *Cancers* 18(3):372, 2026, PMID 41681845, PMC12896505, DOI
[10.3390/cancers18030372](https://doi.org/10.3390/cancers18030372)).

**GBM/glioma-specific groups:**
- **Merz et al. 2013**, *Neuro-Oncology* 15(6):670–681 — human GBM slices dosed with TMZ (200 µM), X-ray,
  and carbon-ion irradiation; caspase-3 readout stratified by MGMT status. Foundational GBM-slice
  drug-response paper. Data availability UNVERIFIED (full text not accessible in this pass).
- **UNC (Satterlee/Hingtgen/Sancar) — Organotypic Brain Slice Culture (OBSC)**: passage-zero patient GBM,
  TMZ ± EdU dosing. *Mol Cancer Ther* 24(8):1213–1225, 2025, PMID 40130554, DOI
  [10.1158/1535-7163.MCT-24-1098](https://doi.org/10.1158/1535-7163.MCT-24-1098); *PNAS* 2026, PMID
  41474743, DOI [10.1073/pnas.2532187123](https://doi.org/10.1073/pnas.2532187123). No data deposition
  found. **Note**: DATA_HUNT.md independently holds `Mann2023_OBSC_doseresponse_extracted.csv` (572 values,
  7 models × 11 drugs) from this lineage, already extracted locally — the companion patient-tumour table
  is flagged unusable there (PDF-scrape artefact).
- **Ravi/Heiland lab (Freiburg)**, *Life Science Alliance* 2019, PMID 31249133, DOI
  [10.26508/lsa.201900305](https://doi.org/10.26508/lsa.201900305) — human organotypic brain slices +
  microinjected patient GBM cells, invasion/drug response + astrocyte RNA-seq. Closest architecture+molecular
  combination found; a GEO/SRA accession could not be confirmed in this pass (403 on full text) — worth a
  direct follow-up given Heiland's group is otherwise a major spatial-GBM-transcriptomics player.
- **GSE148842** (already in DATA_HUNT.md, re-confirmed relevant here): 43 acute GBM slice cultures, 10
  patients, scRNA-seq under vehicle/etoposide/panobinostat/others. Response is transcriptome-wide and
  cell-type-deconvolvable in tissue retaining its own microenvironment; no drug concentration measured;
  slice culture removes vasculature. **Open**, GEO.
- Comba/Lowenstein/Castro (Michigan) glioma explant + confocal time-lapse: analysis code public (GitHub
  smotsch/analysis_glioma, Zenodo DOI [10.5281/zenodo.7587235](https://doi.org/10.5281/zenodo.7587235));
  underlying raw imaging/spatial data availability UNVERIFIED.

**Non-GBM state-of-the-art technical template**: Dong/Aulitzky (Tübingen/AstraZeneca), Perfusion Air
Culture, *Cells* 12(5):807, 2023, DOI [10.3390/cells12050807](https://doi.org/10.3390/cells12050807),
PMC10001200 — cotton-mesh perfusion, >7-day viability, spatially resolved multiplex IHC (Ki67, HIF1α,
γH2AX, cleaved-caspase3, PD-L1, CD4/CD8/CD68/FOXP3) mapped across three tissue depths on the *same*
cisplatin-treated slice. Spatial structure and drug response genuinely co-recorded on the same specimen —
but "available upon request from the corresponding authors," not deposited.

**Drug distribution *within* the slice** (not just endpoint viability): this sub-question had almost no
literature. No paper was found applying MALDI-MSI (or any spatial drug-concentration method) to a
living, dosed PCTS/OBSC slice with a paired spatial-response readout on the same section — every MALDI-MSI
example found elsewhere in this survey (§2.7) uses fixed endpoint tissue, organoids, or xenografts, not
living slice culture.

**Public data check**: searched GEO, Zenodo, Image Data Resource, BioImage Archive. Only public artefact
found in the GBM-slice space beyond GSE148842 is analysis *code* (Comba/Michigan), not raw
imaging/spatial-omics data.

### 2.3 PDX models and encyclopedias

| Resource | Response readout | Spatial data on same specimens? | Data | Access |
|---|---|---|---|---|
| NCI PDMR | Tumour volume vs. time (bulk) | Baseline histology only, not paired with treated cohort | Web DB, FASTQ/VCF | open (browse); bulk-download details unclear |
| Novartis PDXE (Gao et al. 2015) | Bulk volumetric "1×1×1" response | No | GEO GSE78806 + paper supplements | open |
| Mayo GBM PDX (Sarkaria lab, Vaubel 2019) | Bulk survival-ratio (37/96 lines have RT/TMZ/RT+TMZ/BEV data) | Baseline TMA histology (58 lines), not paired with treated cohort | cBioPortal, SRA PRJNA543854/PRJNA548556, GEO GSE132650 | open |
| EurOPDX Data Portal | Bulk growth-curve (283/995 models) | No | dataportal.europdx.eu | registered (free) |

Sources: Vaubel et al., *Clin Cancer Res* 2019, PMID 31852831, PMC7056576; Gao et al., *Nat Med* 2015,
PMID 26479923, DOI [10.1038/nm.3954](https://doi.org/10.1038/nm.3954); Dudová et al., *BMC Genomics* 2022,
PMID 35193494, DOI [10.1186/s12864-022-08367-1](https://doi.org/10.1186/s12864-022-08367-1); NCI PDMR
(pdmr.cancer.gov).

**One PDX-adjacent finding matters more than the table above**: **Lopez, Kohale, Du et al., "Multimodal
platform for assessing drug distribution and response in clinical trials," *Neuro-Oncology* 24(1):64–77,
2021/2022, PMID 34383057, PMC8730776, DOI [10.1093/neuonc/noab197](https://doi.org/10.1093/neuonc/noab197)**
(Agar/White/Santagata labs, Brigham & Women's/MIT/Mayo PSOC). This is the Sarkaria-lineage GBM PDX platform
(GBM12/22/84) combined with MALDI-MSI drug-distribution mapping (30–100 µm) + t-CyCIF single-cell
DNA-damage/apoptosis phenotyping + phosphoproteomics on adjacent sections, applied to **both PDX flank
tumours and human phase-1 trial specimens (NCT01849146, adavosertib)**. This is genuinely per-cell +
per-region drug-response data co-registered with spatial structure, in GBM — but small (6 human patients,
one drug) and, critically, **the MALDI-MSI drug-distribution data itself is stated verbatim as "available
upon request"**, with no accession; only the phosphoproteomics reached PRIDE (PXD018782 — but per
DATA_HUNT.md's independent verification, PXD018782 contains only raw spectra/IDs, zero quantified
phosphosite response). t-CyCIF images are "uploaded to Omero" with no public identifier found. This is
already documented as a Tier-3 gated asset in `results/DATA_HUNT.md` §4; this survey independently
reconfirms every detail of that entry and adds nothing new on access, only confirmation.

### 2.4 Microfluidic tumour-on-chip / organ-on-chip

No platform surveyed — commercial or academic — publishes an open dataset pairing (a) tumour architecture,
(b) a directly imaged spatial drug-concentration gradient, and (c) spatially resolved cell-kill, on real
patient material.

- **Emulate Bio** — EAC-Chip uses real patient esophageal biopsies + matched CAFs (n=8), scores drug
  sensitivity via bulk viability/histology/secreted biomarkers over 2–3 months; no spatial drug-gradient
  imaging described. No public dataset found.
- **CN Bio** — core product is liver-on-chip; oncology use is indirect (metastasis-to-liver models using
  cell lines, not patient tissue).
- **Javelin Biotech** — scope correction: this is a liver ADME/PK chip company (Pfizer collaboration), not
  a tumour-on-chip platform. Does not apply to this survey category.
- **Mimetas OrganoPlate** (Graft format) — supports patient tumour grafts by design; no publication found
  that validates a spatial drug-gradient + position-resolved response pairing in an oncology context.
- **Horowitz et al. 2020** (Univ. Washington), *npj Precis Oncol*, PMID 32435696, PMC7237421, DOI
  [10.1038/s41698-020-0117-y](https://doi.org/10.1038/s41698-020-0117-y) — **strongest hit in this
  category**: microfluidic multiplexer perfusing intact human GBM-xenograft and patient tumour slices;
  directly imaged fluorescent doxorubicin distribution and separately quantified regional cell death
  (Sytox Green+), but as **separate experiments, not simultaneous on the same tissue region**. Data
  "available from corresponding authors on reasonable request" — not public.
- **Veith et al. 2024** (Institut Curie), *Cell Rep Med*, PMID 38703767, PMC11148770 — patient-derived NSCLC
  tumour-on-chip, single-cell live-imaging tracking (>1,500 trajectories) of anti-PD-1 response; uniform
  dosing, no spatial gradient measured; only analysis code deposited, not raw imaging.

### 2.5 Ex vivo functional precision medicine (non-spatial)

- **Dynamic BH3 profiling**: MOMP readout by flow cytometry on **dissociated** cells — tissue architecture
  destroyed before the assay. GBM application: Fernandez et al., *Nat Commun* 2024, PMID 39572533,
  PMC11582606, DOI [10.1038/s41467-024-54138-9](https://doi.org/10.1038/s41467-024-54138-9). No spatial
  variant found. Access: per-publication only, no repository.
- **FIMM DSRT**: 384-well viability against ~350–530 compounds, dissociated cells, bulk well-average.
  Analysis tool **Breeze**, open, no login: [breeze.fimm.fi/v2](https://breeze.fimm.fi/v2/); some datasets
  feed the open aggregator **DrugComb** ([drugcomb.org](https://drugcomb.org), DOI
  [10.1093/nar/gkab294](https://doi.org/10.1093/nar/gkab294)). No spatial component — structurally cannot
  supply the pairing tme-sim needs, but a legitimate open source of potency-axis priors.

Neither platform preserves or records spatial structure at readout; both discard it by design
(dissociation is required for flow cytometry / plate-based viability).

### 2.6 Intratumoral microdosing — CIVO and the GBM-specific alternative

**CIVO (Presage Biosciences)**. Handheld injector depositing up to 8 drugs as parallel micro-columns
(≤10 mm long, ~1.5 mm diffusion radius) directly into a living, surface-accessible tumour, each co-dosed
with fluorescent tracking microspheres for registration; tissue is resected during the patient's
already-planned surgery, sectioned transverse to the columns, and each site profiled by IHC/ISH and,
since ~2023, GeoMx DSP (bulk-per-ROI) and CosMx (single-cell). Six clinical trials, 68 subjects, three
indications: lymphosarcoma, soft-tissue sarcoma, head & neck cancer (Gundle et al., *Front Pharmacol* 2024,
PMID 38681192, PMC11048044, DOI [10.3389/fphar.2024.1367581](https://doi.org/10.3389/fphar.2024.1367581)).
**Explicitly not used in glioma** — the paper states the platform is "limited to surface-accessible
injections (~3 cm)" and has "only been used in LSA, STS, and HNC," which structurally excludes deep
intracranial tumours.

"Drug distribution" in CIVO is not literal LC-MS/MRI concentration mapping; it is inferred from fluorescent
microsphere position plus, where possible, a target-engagement IHC proxy (e.g., anti-SUMO adduct staining
for subasumstat). Response (SUMO pathway inhibition, IFN-stimulated gene expression, Ki67, cleaved
caspase-3) is scored on the same or serial FFPE sections, co-registered to the same ROI (Derry et al.,
*Clin Cancer Res* 2023, PMID 37389981, PMC10502463, DOI
[10.1158/1078-0432.CCR-23-0827](https://doi.org/10.1158/1078-0432.CCR-23-0827)).

**New finding — this dataset is genuinely public and was independently verified in this pass**: **GEO
`GSE235475`** ("Trackable Intratumor Microdosing and Spatial Profiling Provide Early Insights into
Activity of Investigational Agents in the Intact Tumor Microenvironment," public since 2023-06-28, PMID
37389981). Verified directly against the GEO record:
- 354 areas of interest (AOIs, GeoMx 400 µm² ROIs) from 9 head-and-neck-carcinoma patients, subasumstat
  (TAK-981) vs. vehicle.
- Per-sample GEO characteristics include `mil113_grade` (IHC-graded drug-exposure category: low/high) and
  `mil113_percent` (continuous exposure score) **alongside** a 1,812-gene (Cancer Transcriptome Atlas)
  expression profile for that same AOI — i.e., a spatial exposure grade and a spatial molecular-response
  readout, on the same specimen, at matched ROI coordinates.
- Processed data (`GSE235475_TMM_normalized-CPM_values.txt.gz`, `GSE235475_raw_data_file.txt.gz`) is
  **openly downloadable from the GEO FTP, no registration** (confirmed by direct fetch in this pass). Raw
  FASTQ is noted in the series record as "requested" (not yet deposited as of last update).
- Access tier: **open** (processed counts + exposure grade); raw sequencing reads not yet public.
- Caveats: head & neck, not glioma; exposure is a categorical IHC grade, not an absolute drug
  concentration; readout is AOI/regional (400 µm² footprints), not single-cell (CosMx single-cell data from
  the same trial does not appear to be in this GEO series).

**GBM-specific alternative — Jonas Lab (Brigham & Women's)/Kibur Medical intratumoral microdevice (IMD)**.
First-in-human pilot in high-grade glioma, n=6 patients, NCT04135807 (Peruzzi et al., *Sci Transl Med*
2023, PMID 37672566, PMC10754230, DOI
[10.1126/scitranslmed.adi0069](https://doi.org/10.1126/scitranslmed.adi0069)). Each patient received two
IMDs (10–15 mm apart) implanted intraoperatively during standard resection, each releasing up to ~20 drugs
into spatially distinct micro-reservoirs. At each reservoir, the authors measured: **drug concentration**
(autofluorescence imaging for doxorubicin/lapatinib; MALDI-TOF MS for non-fluorescing compounds), **and
response** (cleaved caspase-3 and phospho-H2AX IHC, spatial transcriptomics, MALDI-based metabolomics),
explicitly exploiting the distance-dependent concentration gradient around each reservoir to correlate
local dose with local response on serial sections of the **same** resected tissue. 11/12 specimens were
usable and results correlated with clinical/radiological TMZ response.

**Data verified directly in this pass**: **Zenodo [10.5281/zenodo.8226867](https://doi.org/10.5281/zenodo.8226867)**,
**open, CC-BY-4.0, no registration** (fetched the record and file list directly). Contents: `Transcriptomic
data - Nanostring.xlsx` (110 KB) and `MALDI data - glutathione and PEG.h5ad` (25 MB). This is the exact
experiment, in the exact disease, with an open deposit — but read the file list precisely: the deposited
MALDI layer covers glutathione and PEG (a metabolite and the delivery-vehicle polymer, used as
tissue-response/reservoir markers), **not** ion images of the therapeutic drugs themselves; the paper's
data-availability statement separately says "all data associated with this study are present in the paper
or the Supplementary Materials," implying the core per-reservoir drug-concentration and IHC-response
quantification (the actual pairing) may live in supplementary tables/figures of the *Sci Transl Med* paper
rather than as a machine-readable table in the Zenodo deposit. This needs a direct check of the paper's
Supplementary Materials (not completed in this pass) before treating the pairing as fully open — but it is
by a wide margin the most promising, most open, most GBM-relevant lead in this entire survey.

**Company status**: Presage Biosciences is active (AstraZeneca partnership Nov 2023 in head & neck,
FDA "Study May Proceed" for PBA-0405 2024, VIB-KU Leuven collaboration ~March 2025); external access is via
a commercial "Clinical Collaboration" contact, no open academic data-sharing program. Kibur Medical
(commercializing the Jonas Lab IMD) is an active company per its own press materials (2024 AACR
presentation on high-grade-glioma microdosing data).

### 2.7 MALDI-MSI / DESI-MSI of drug distribution in tumour tissue

See §5 for the itemized dataset list with URLs. Summary: the technique is well established in pharma for
exactly this measurement (drug ion images in tumour sections), but almost none of it is deposited as raw,
downloadable imaging data — nearly every paper found states "available upon request" or says nothing about
deposition. The one confirmed open MSI drug-distribution dataset in solid tumour tissue (Prasad et al.
2020, GigaScience) is non-CNS (ovarian/colon xenograft, paclitaxel). METASPACE's own description paper
confirms the platform hosts real, open pharma drug-distribution data (a cassette-dosed rat-liver
erlotinib/olanzapine dataset), but this pass could not confirm or deny a tumour/glioma-specific drug
dataset there — METASPACE's web app is JS-rendered and not scrapeable by the tools available in this pass;
this is a **negative-by-tool-limitation**, not a negative-by-search, and should be manually re-checked
(register at metaspace2020.org and browse/search by molecule).

---

## 3. Critical question 1 — ranking by how directly each platform could validate the RESPONSE axis

Ranked by directness of the pairing (spatial drug field ↔ spatial biological response, same specimen), with
the concrete validation experiment each would support.

1. **Jonas Lab/Kibur Medical GBM intratumoral microdevice** (§2.6). *Why #1*: real human GBM, spatially
   distinct dosing sites with matched drug-concentration and biological-response readouts on the same
   resected tissue, and a CC-BY-4.0 Zenodo deposit already exists. *Experiment*: request/inspect the Sci
   Transl Med Supplementary Materials for the per-reservoir concentration-vs-distance and
   caspase-3/pH2AX-vs-distance tables; if machine-readable, feed each patient's reservoir geometry and the
   measured local drug concentration into tme-sim's transport model as a boundary condition (point source,
   not a vascular network — closer to the CIVO/Moncayo bath-like local-source regime than to whole-tumour
   perfusion), and score predicted-vs-measured local response (caspase-3/pH2AX) at matched distances from
   each reservoir. This is the only route in this entire survey to a within-glioma, same-specimen
   delivery→kill validation.
2. **Agar/White/Santagata GBM PDX + human-trial multimodal platform** (Lopez et al. 2021/2022, §2.3). *Why
   #2*: methodologically the single closest match to what tme-sim predicts — a literal MALDI-MSI drug
   concentration map paired with single-cell (t-CyCIF) DNA-damage/apoptosis phenotyping, in GBM12/22/84 PDX
   and human phase-1 tissue. *Why not #1*: data is gated ("available upon request," no accession); this is
   a collaboration ask (contact Nathalie Agar, Brigham & Women's; per DATA_HUNT.md's Tier-3 entry, expect
   months and an IRB/DUA process for the phase-1 human tissue). *Experiment*: identical in structure to #1
   — feed PDX vascular/structural geometry (if obtainable) into tme-sim, predict the spatial adavosertib
   field, compare to the measured MALDI-MSI map; compare predicted local kill probability to the measured
   t-CyCIF apoptotic fraction at registered coordinates.
3. **CIVO / GSE235475** (§2.6). *Why #3*: genuinely open today, no request needed — but head & neck, not
   glioma, and the "drug distribution" side is a categorical IHC exposure grade, not a diffusion field with
   coordinates, which limits it to a transport-agnostic dose-response correlation rather than a true
   spatial-PDE validation. *Experiment*: as a **methodology proof-of-concept, not a glioma validation**,
   regress `mil113_percent` (exposure) against the IFN-response/cell-cycle gene signatures also present in
   the expression matrix, replicating the general "local exposure predicts local transcriptional response"
   logic tme-sim's kill model assumes, on a genuinely public dataset, before spending effort chasing gated
   GBM data.
4. **UNC OBSC / Merz-style / Freiburg GBM slice culture groups** (§2.2). *Why #4*: real GBM tissue, real
   TMZ dosing, architecture nominally preserved, but none of the found groups pairs a spatial *drug
   concentration* map with the spatial response — response is read out with the drug treated as a uniform
   bath, not a field with spatial structure. *Experiment*: a direct data-sharing ask to Merz/Freiburg/UNC
   for TMZ dose-response + IHC spatial-biomarker tables (Ki67/caspase-3/γH2AX by depth or region), useful
   for the kill-model's dose-response calibration even without a paired drug map — analogous to what
   `results/gbm_tmz_exvivo.md` in this repo likely already does with published aggregate TMZ ex vivo data
   (not reviewed in this pass; check for overlap before re-requesting).
5. **Horowitz (UW) tumour-slice microfluidic platform** (§2.4). *Why #5*: real GBM xenograft/patient slices,
   directly imaged drug (doxorubicin) distribution — but measured in a separate experiment from the cell-death
   readout, not simultaneously on the same tissue. *Experiment*: request the underlying imaging data; use
   the fluorescence-based penetration profile alone to validate transport (independent of response), the
   same category of test as the METASPACE GIST imatinib set in DATA_HUNT.md T1.3.
6. **Mayo GBM PDX panel (Vaubel/Sarkaria) and NCI PDMR** (§2.3). Bulk growth-curve/survival-ratio response
   only, no spatial pairing on the treated cohort. *Use*: potency-axis and molecular-response priors (as
   DATA_HUNT.md Tier 2 already uses the Vaubel table), not a spatial validation.
7. **GBM organoids (GBOs) / PDO biobanks** (§2.1). No spatial pairing at readout. *Use*: cell-state/clonal
   heterogeneity priors for the microenvironment reconstruction step, not the response axis.
8. **Dynamic BH3 profiling / FIMM DSRT** (§2.5). Structurally non-spatial (dissociated cells). *Use*:
   potency-axis priors only (as DepMap/PRISM already serve in DATA_HUNT.md Tier 2).

---

## 4. Critical question 2 — does any public dataset pair spatial tissue structure with measured drug response on the same specimen?

**Honest answer, searched hard**: **Not fully, and not in glioma — but the answer is closer to "barely
yes, twice, both outside glioma" than to a flat no**, and there is one credible near-miss inside glioma
that needs one more verification step before we can say no there too.

1. **Moncayo et al. 2026** (already established by `results/DATA_HUNT.md`, re-confirmed here, not
   re-searched from scratch): open, HGSOC, ex vivo bath-dosed explants, MALDI-MSI drug map + GeoMx
   transcriptome on matched ROIs. Already run through tme-sim (`results/moncayo_explant.md`) — honest null.
2. **CIVO GSE235475** (new in this pass, verified directly against GEO): open, head & neck carcinoma,
   in-human intratumoral microdosing, IHC-graded spatial exposure + spatial (AOI-level) transcriptome on
   matched ROIs. Not yet exploited by this project. Weaker than Moncayo on the "drug field" side (a
   categorical grade, not a concentration map with coordinates) but stronger on the "in-human, real
   perfused tumour" side (Moncayo is ex vivo bath-dosed, no vasculature; CIVO is in situ in a living,
   perfused tumour).
3. **Jonas Lab/Kibur GBM microdevice** (new in this pass): the platform demonstrably measures both sides
   of the pairing on matched coordinates in glioma, and a CC-BY-4.0 Zenodo deposit exists — but this pass
   could not confirm that the deposited files (Nanostring + glutathione/PEG MALDI) or the paper's
   Supplementary Materials actually carry the machine-readable per-reservoir drug-concentration ×
   response table, as opposed to that pairing existing only as summary figures. **This is the single most
   important open item from this survey**: read PMC10754230's Supplementary Materials in full before
   concluding either way. If the per-reservoir table is there, this becomes a third — and the first
   glioma-specific — public "yes."
4. **Agar/White GBM MALDI-MSI + t-CyCIF** (Lopez et al., confirmed gated, already known via DATA_HUNT.md):
   the methodologically closest match to what tme-sim needs, in the right disease — and confirmed **not**
   public (verbatim "available upon request," no accession for the imaging arm).

So: zero *fully open, machine-readable, glioma-specific* pairings confirmed today. Two non-glioma open
pairings exist (one already exploited with a null result, one newly identified and unexploited). One
glioma-specific platform almost certainly has the right measurement but its exact public-data completeness
is unresolved pending a Supplementary Materials read. One glioma-specific platform has the right
measurement but is fully gated. This is a materially more optimistic picture than "no open dataset pairs
these things at all," but it does not change DATA_HUNT.md's core verdict for glioma specifically: **there
is still no public, machine-readable, glioma-specific dataset pairing a measured spatial drug field with a
measured spatial response**, pending the one open item above.

---

## 5. Critical question 3 — MALDI/DESI datasets, itemized, with URLs

| Dataset / paper | Drug | Tissue | Downloadable? | Access tier | Source |
|---|---|---|---|---|---|
| Prasad et al. 2020, GigaScience | Paclitaxel ± bevacizumab | Mouse xenograft, ovarian (A2780) + colon (HCT116) cell lines — **not PDX, not GBM** | Yes — GigaDB + DNAS-KNAW; R package `CorrDrugTumorMSI` on GitHub | **open** | PMID 33241286, PMC7688471, DOI [10.1093/gigascience/giaa131](https://doi.org/10.1093/gigascience/giaa131) |
| Morosi et al. 2017, Sci Rep | Paclitaxel | Ovarian PDX (HOC84), mesothelioma PDX (MPM484), + cell-line xenografts | No data-availability statement found | **none** | PMC5175283 |
| Lopez et al. 2021/2022, Neuro-Oncology | Adavosertib | **GBM PDX (GBM12/22/84) + human phase-1 GBM trial specimens** | MSI data "available upon request"; only proteomics on PRIDE (PXD018782, itself non-quantitative per DATA_HUNT.md verification) | **none** (MSI) / open (raw spectra only, proteomics) | PMID 34383057, PMC8730776, DOI [10.1093/neuonc/noab197](https://doi.org/10.1093/neuonc/noab197) |
| Pokorny et al. 2015, Clin Cancer Res | MK-1775 (adavosertib precursor) | GBM22 PDX, orthotopic vs. heterotopic | Figures only | **none** | PMID 25609063, PMC4401631 |
| Sharifi et al. 2019, Clin Cancer Res | ZR2002 | Orthotopic EGFRvIII/GSC mouse brain tumour | Figures only | **none** | PMID 31540977 |
| METASPACE — general | multiple (incl. cassette-dosed rat-liver erlotinib/olanzapine confirmed via the platform's own description paper) | Liver (confirmed); tumour/glioma-specific holding **unconfirmed** in this pass | Ion images open; raw imzML/ibd registration-gated | **open** (images) / **registered** (raw) | Palmer et al., bioRxiv [10.1101/539478](https://www.biorxiv.org/content/10.1101/539478v2.full) → Nat Methods 2017; per DATA_HUNT.md, METASPACE holds **zero** temozolomide-annotated datasets |
| PRIDE / ProteomeXchange | — | — | Per DATA_HUNT.md's own verification: exactly **one** imzML dataset in the entire federation (PXD001283, mouse bladder phospholipids — no drug, no tumour) | **open** (that one, irrelevant, dataset) | Race et al., PMC4357645 |
| MassIVE (UCSD CCMS) | — | — | No imaging-MS drug-tumour dataset found | **none** | — |
| Mario Negri Zenodo deposits (T1.2/T1.4/T1.5 in DATA_HUNT.md) | Paclitaxel, niraparib, olaparib, pamiparib | Xenograft (ovarian, pancreatic, MDR-overexpressing lines) — **not GBM** | Yes | **open** (Zenodo 5776125, 5158808, 13375209, CC-BY-4.0) | see DATA_HUNT.md §2A for full column-level detail; already on disk in this repo's `data_cache/hunt/` |
| Moncayo et al. 2026 (already in DATA_HUNT.md) | Niraparib, olaparib, rucaparib | Human HGSOC explants | Yes — 376,736-pixel per-pixel table in journal Source Data + Zenodo/GEO ROI summaries | **open** | Zenodo 17610220, GSE281519, PMID 41844641 |

**Confirmed dead ends** (carried over from DATA_HUNT.md, independently re-confirmed in this pass): Zenodo
full-text search `imzML AND (drug OR dosed OR pharmacokinetic)` returns zero records; the well-known 3D
paclitaxel-penetration paper (Sci Rep srep37027) deposits only `.avi` rotation videos, no raw MSI data.

**Bottom line for Q3**: there is still no public, downloadable MALDI/DESI drug-distribution dataset in
glioma tissue (human or PDX), open or gated-but-accessioned. The Agar/White GBM PDX + human-trial dataset
(Lopez et al.) is the exact right experiment and is fully gated. METASPACE's tumour/glioma holdings could
not be exhaustively ruled out in this pass due to its JS-rendered interface — manually registering and
searching by molecule name is a same-day, zero-cost action worth doing before concluding nothing is there.

---

## 6. Critical question 4 — which platform to pursue, and what to ask for

**Primary recommendation: the Jonas Lab (Brigham & Women's, Oliver Jonas)/Kibur Medical GBM intratumoral
microdevice.** This is the only platform in the entire survey that (a) is glioma-specific, (b) has already
run in humans, (c) explicitly measures spatially resolved drug concentration and spatially resolved
biological response on matched coordinates in the same resected tumour, and (d) already has an open
CC-BY-4.0 data deposit — i.e., the relationship is not starting from zero. **Concrete ask**: contact Oliver
Jonas's lab (and/or Kibur Medical, which is commercializing the platform and presented new high-grade
glioma microdosing data at AACR 2024) and request (i) the per-reservoir, per-distance drug-concentration
table and the matched caspase-3/phospho-H2AX quantification underlying Fig. 2–4 of the *Sci Transl Med*
2023 paper, in machine-readable form, for all 6 patients/12 specimens; (ii) whether they would be willing
to run tme-sim's transport-model prediction against their next glioma cohort prospectively, using the
known reservoir geometry and drug PK as the model's boundary condition — this would be a genuine,
forward-looking validation rather than a retrospective data-mining ask, and gives them a novel
computational-interpretation offering in return (their published analysis does not appear to include a
mechanistic transport-PDE fit, only correlational dose-vs-distance-vs-response).

**Secondary recommendation, if the primary stalls**: Nathalie Agar (Brigham & Women's) and Sandro Santagata
(MIT-Mayo PSOC), for the adavosertib GBM PDX + phase-1 MALDI-MSI + t-CyCIF dataset (Lopez et al. 2021/2022).
This is already flagged as a Tier-3 gated asset in `results/DATA_HUNT.md` with the same contact path; this
survey adds no new access route, only independent confirmation that the ask is worth making. Expect
months, an IRB/DUA process for the phase-1 human tissue, and a narrower deliverable (6 patients, one drug)
than the Jonas Lab ask.

**Tertiary, lowest-effort, immediately actionable**: pull and reanalyze **CIVO's `GSE235475`** now — it
needs no outreach, no waiting, no DUA. It will not validate tme-sim's transport-PDE (no drug-concentration
field with coordinates, wrong disease), but it is a genuinely public dataset that pairs a spatial exposure
variable with a spatial molecular-response readout, and running the same class of analysis
`results/moncayo_explant.md` already ran on Moncayo (regress predicted/measured exposure against
response-signature scores, with the same honest permutation-null discipline) would either strengthen or
further stress-test the general "local exposure predicts local transcriptional response" assumption
underlying tme-sim's kill model, independent of glioma-specific transport physics.

**What not to pursue further, and why**: CIVO/Presage as a company partnership (commercial-only access
path, and the platform is structurally excluded from intracranial tumours by its needle-based, ≤3 cm-depth
design); PDO/PDX biobanks and companies (Xilis, SEngine, Cellesce, Crown Bioscience, HUB) — all
commercial-only and none records spatial structure at the response readout; dynamic BH3 profiling and FIMM
DSRT — structurally non-spatial by assay design, useful only as potency-axis priors already served better
by DepMap/PRISM per DATA_HUNT.md.

---

## 7. What this document adds to the project's existing knowledge, stated plainly

- Confirms, via independent re-verification, every DATA_HUNT.md claim touched in this pass (Moncayo
  pairing, Agar/White GBM MALDI-MSI gating, METASPACE's zero-TMZ finding, PRIDE's single irrelevant imzML
  record, Vaubel Mayo PDX panel).
- Adds two platforms/datasets not previously identified in this project: **CIVO's `GSE235475`** (open, HNC,
  spatial exposure + spatial transcriptome pairing) and the **Jonas Lab/Kibur Medical GBM intratumoral
  microdevice** with its Zenodo deposit `10.5281/zenodo.8226867` (open, glioma-specific, the platform most
  likely to close this project's core validation gap through direct collaboration rather than data mining).
- Leaves one specific, cheap, high-value follow-up unfinished for whoever picks this up next: read
  PMC10754230's Supplementary Materials in full to determine whether the per-reservoir
  drug-concentration × response table is already machine-readable in the paper or Zenodo deposit, which
  would upgrade §4's answer from "one glioma-specific near-miss" to "one glioma-specific yes."
