# Landscape survey: public serial/multi-section GBM/glioma tissue data, for a GBM-specific rerun of `registration_3d.py`

**Purpose.** `results/registration_3d_breast_poc.md` proved tme-sim's serial-section 3D
registration method (`tmesim/reconstruct/registration_3d.py`) works on real tissue —
but on breast cancer Xenium data, because `results/landscape_serial_sections.md` searched
broadly across cancer types and found the breast pair as the best public serial-section
spatial-transcriptomics data available. tme-sim's flagship cancer is glioblastoma (GBM),
chosen specifically because the blood-brain barrier makes drug *delivery* the dominant
clinical variable there (`results/landscape_cancer_census.md`). This document asks one
question, cast as wide as the task brief allows within the GBM/glioma constraint: **is
there any public, real, serial/multi-section GBM or glioma tissue dataset — spatial
transcriptomics, histology/IHC, vasculature, organoid/PDX imaging, or coarse multi-region
— that would let the registration pipeline be re-run on the project's actual flagship
cancer instead of breast?**

**Method.** Read `results/registration_3d_breast_poc.md`, `tmesim/reconstruct/registration_3d.py`,
`results/landscape_serial_sections.md`, `results/landscape_cancer_census.md`,
`results/landscape_pathology_3d.md`, `results/gbmspace_ingestion.md`, `results/heldout_patient.md`,
`tmesim/data/real/ivygap.py` + `docs/real_data/ivygap.md`, and `docs/outreach/04_holzel_gse324268.md`
first, to avoid re-deriving what this project already knows. Then PubMed (via the PubMed
MCP, full text fetched directly wherever possible — not abstracts alone) and WebSearch/WebFetch
were used to re-verify the two closest prior negative findings (Torres-Sánchez 2021, FineVess
2024) and to search fresh, specifically for GBM/glioma serial- or multi-section data in any
modality. Every claim is cited; anything not confirmed by direct primary-source fetch is marked
**UNVERIFIED**. No code was modified — this is a read-only research pass, one file written.

---

## Verdict — stated plainly, as the task asked

**There is no public GBM dataset of literal serial sections (same tissue block, consecutive
cuts, a stated or inferable micron z-spacing) in spatial transcriptomics, histology, or any
other modality.** Every candidate found — three independent GBM spatial-transcriptomics
atlases, the Ivy GAP anatomic atlas, a 3D-neuronavigated multi-omics atlas, two human-GBM
light-sheet vasculature studies — samples **multiple spatially distinct regions of the same
tumour**, not a continuous serial stack through one block. That specific thing
(`landscape_serial_sections.md`'s Tier 1/Tier 2 breast/DLPFC pairs) does not have a GBM
analogue in the public domain today, as far as this pass could determine.

**But this pass found something better than "nothing," and it is GBM-specific and real:**
Greenwald, Darnell, Hoefflin et al. (Tirosh/Suvà labs), *"Integrative spatial analysis reveals
a multi-layered organization of glioblastoma,"* *Cell* 187(10):2485–2501 (2024)
([PMID 38653236](https://pubmed.ncbi.nlm.nih.gov/38653236/), [DOI 10.1016/j.cell.2024.03.029](https://doi.org/10.1016/j.cell.2024.03.029))
profiled **13 IDH-wildtype human GBM tumours by 10x Visium and 12 by CODEX (40-marker spatial
proteomics) on *near-adjacent* tissue sections from the same blocks**, and — this is the part
that matters most for `registration_3d.py` specifically — **cross-registered the CODEX and
Visium sections of each tumour into one shared coordinate frame using STalign**, the exact
diffeomorphic registration tool `results/landscape_pathology_3d.md` and
`results/landscape_serial_sections.md` already vetted as a peer for this project's own method.
This is real, GBM-specific, near-adjacent-section, cross-modality registration, already done
and published by an independent group, with data on GEO (raw Visium,
[GSE237183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237183)) and Zenodo (processed
Visium + CODEX imaging, [10.5281/zenodo.8105466](https://zenodo.org/records/8105466)) — **no
DAC, no login required to at least inspect and attempt a download.** It is not the same shape as
the breast Xenium-to-Xenium pair (it's Visium spots to CODEX single-cell protein images, not two
identical-modality point clouds; z-spacing between the CODEX and Visium sections is not stated
in the paper, only that they are "near-adjacent"), and CODEX ingestion is genuinely new-code
scope for `ingest.py` (no CODEX loader exists today). But it is the single most concrete,
GBM-specific "next target" this survey found — see the recommendation at the end.

Everything else confirms and extends the negative findings already on record:
`results/landscape_pathology_3d.md`'s "no publicly downloadable 3D vascular network of human
glioblastoma anywhere" and "CODA has never been applied to brain/GBM" both **still hold** as of
this pass (re-checked directly, not assumed). No human GBM organoid/PDX light-sheet study was
found. Ivy GAP has **no** spatial x/y/z coordinates, confirmed by direct code/schema inspection
— it cannot support even a coarse 3D reconstruction on its own. One genuinely new, real-3D-coordinate
GBM resource turned up (Mathur et al. 2024, UCSF/Costello lab) but its raw data is EGA-gated and
it is bulk multi-omics on discrete stereotactic biopsies, not spatial imaging — a reference point,
not an ingestion target.

---

## 1. Spatial transcriptomics: GBM/glioma multi-section datasets, checked one by one

### 1.1 GSE324268 — checked directly: NOT serial sections, NOT usable for this purpose

This project's own data cache and outreach draft (`docs/outreach/04_holzel_gse324268.md`)
already uses GSE324268 for its MGMT probes. Fetched the GEO record directly in this pass to
answer the task's specific question — does it include multiple sections per patient — and the
answer is **no**. GSE324268 has exactly **2 samples**: `GSM9572843` (`GBM_tumor07`, a primary
resected GBM specimen, 851,509 cells) and `GSM9572842` (`GBM1_organoid_brain_slice`, a
patient-derived GBM organoid cultured on a brain slice). These are **two different specimens
from what the record indicates are different sources** (one primary resection, one organoid),
not two sections of the same block. This closes the question cleanly: GSE324268 has nothing to
offer the registration pipeline. Source: direct GEO fetch,
[GSE324268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324268), 2026-07-28.

### 1.2 GBM-Space (Sanger/EGA, already known to this project) — multi-region, not confirmed serial; format already blocked

`results/gbmspace_ingestion.md` already established that GBM-Space's processed deposit
(`gbmspace.org`, EGA [EGAD00001015527](https://ega-archive.org/datasets/EGAD00001015527)) is a
flat per-section `.h5ad` archive, not a Space Ranger `outs/` tree, and stopped there per its own
scope (no new ingester built). This pass did **not** re-attempt ingestion (per this task's
explicit instruction not to duplicate that work) and instead tried to answer the narrower
question: does GBM-Space's own methodology include literal serial/consecutive sections of the
same block at a known spacing? `results/landscape_cancer_census.md` already records GBM-Space as
"97 Visium sections across 12 GBM tumours, 4–15 sampling sites per tumour" — the phrase "sampling
sites" is the operative ambiguity. Direct fetch of GBM-Space's own bioRxiv preprints
([2025.05.13.653495](https://www.biorxiv.org/content/10.1101/2025.05.13.653495v1),
[2025.05.13.653733](https://www.biorxiv.org/content/10.1101/2025.05.13.653733v1)) was **blocked
by HTTP 403 in this pass** (bioRxiv's automated-fetch guard), so GBM-Space's own methods text
could not be directly quoted here. **Stated as UNVERIFIED, not assumed either way**: whether any
of GBM-Space's 4–15 sampling sites per tumour are literal serial sections of one block at a
stated µm spacing, or are always spatially distinct regions/blocks — the public-facing summaries
found (gbmspace.org, Sanger Institute collaboration page) describe multi-region sampling
("necrotic," "infiltrating," "T1 contrast-enhancing" by navigation-guided resection — language
that, on inspection, turned out in this pass to belong to a *different* paper, §1.3 below, not
GBM-Space; see the correction there) but do not, in what could be fetched, explicitly document a
serial-section design. Given §1.3's and §1.4's independent confirmation that the two other GBM
spatial atlases checked in this pass are explicitly multi-region and not serial, the balance of
evidence is that GBM-Space is the same, but this is not a primary-source-confirmed claim for
GBM-Space specifically. **Format remains blocked either way** — `.h5ad`, not `ingest.py`-native
— so this is not an actionable target regardless of the serial-section question.

### 1.3 Greenwald/Tirosh/Suvà et al. 2024 — the actual best finding, corrected attribution

**Important correction made during this pass.** An early WebSearch summary attributed the
sentence *"In four cases, [we] analyzed multiple frozen tissue blocks from the same tumor, each
isolated from a different region and annotated as necrotic, infiltrating, or T1
contrast-enhancing"* to GBM-Space. Direct full-text fetch of the primary source (via the PubMed
MCP, PMC11088502) shows this sentence is verbatim from **Greenwald, Darnell, Hoefflin et al.,
*Cell* 2024** (§ full citation above), a **separate, independent GBM spatial atlas** from the
Weizmann Institute / MGH / University Hospital Zurich groups (lead/corresponding authors Mario L.
Suvà and Itay Tirosh), unrelated to the Sanger/EGA GBM-Space project. This report treats the two
as distinct datasets throughout — conflating them would have been exactly the kind of error this
task's honesty standard is meant to catch.

**What this paper actually is, verified by direct full-text fetch:**

- **19 glioma samples by 10x Visium** (13 IDH-wildtype human GBM + 6 IDH-mutant lower-grade
  glioma) and **12 by CODEX** (Akoya PhenoCycler-Fusion, 40-antibody spatial-proteomics single-cell
  imaging), from patients at University Hospital Zurich, Massachusetts General Hospital, and
  Brigham and Women's Hospital.
- **In four cases, multiple frozen tissue blocks from the same tumour were profiled by Visium**,
  each block from a different anatomic region (necrotic / infiltrating / contrast-enhancing) —
  confirming this is the same **multi-region, not serial-section** pattern found everywhere else
  in this survey.
- **The part that matters most for `registration_3d.py`**: CODEX was performed on **"near-adjacent"
  sections to the Visium sections of the same tumours**, and the two were **cross-registered into
  one shared spatial frame using STalign** — quoted directly from the Methods: *"STalign uses
  Large Deformation Diffeomorphic Metric Mapping (LDDMM) to align CODEX to Visium datasets using
  image varifolds. A rasterized density image is created from the centroid positions of the CODEX
  data and serves as source image. In the first step, the source image was aligned to the target
  H&E Visium image based on manually assigned landmark points using affine transformations...
  The remaining [samples] were subjected to a second step in which the source (CODEX density
  image) and target (Visium H&E) images were aligned [by diffeomorphic warp]."* The registration
  was **validated quantitatively, not just visually**: *"In 80% of the samples, MP [metaprogram]
  assignments were significantly correlated between CODEX and Visium data from near-adjacent
  sections (p<0.05, t-test)"* and *"Correlation between relative abundance of structured spots vs.
  disorganized spots in near adjacent samples in Visium and CODEX, r=0.73, p=0.011 by Pearson
  correlation."*
- **Section thickness / spacing**: Visium sections were cut at **10 µm** on a cryostat (quoted
  directly: *"Fresh frozen samples were sectioned at 10uM thickness with a cryostat onto 10X
  Visium Spatial Transcriptomics slides"*). **The exact number of microns separating the CODEX
  section from its paired Visium section is not stated in the fetched text** — "near-adjacent" is
  the only spacing descriptor found, the same category of gap `registration_3d.py`'s own
  `Z_SPACING_CAVEAT` already documents honestly for the breast Xenium Rep1/Rep2 pair. **UNVERIFIED,
  flagged exactly as such — do not assume a number.**
- **Data availability, confirmed by direct fetch of the paper's own Resource Availability
  section**: *"Spatial transcriptomics data have been deposited at GEO (raw data) and Zenodo
  (processed data)... Imaging (CODEX) data and antibody validations have been deposited at
  Zenodo."* Confirmed accessions by direct WebFetch of the actual records:
  - **GEO [GSE237183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237183)** — 19 raw
    Visium samples (13 GBM + 6 IDH-mutant), one `GSE237183_RAW.tar` (500.9 MB) containing CSV, H5,
    JPG, JSON, MTX, PNG, TSV files — file-type mix consistent with (but not confirmed identical
    to) a Space Ranger `outs/` tree; **directory layout not verified against `ingest.py`'s
    `resolve_visium_dirs()` in this pass — plausible, not confirmed.**
  - **Zenodo [10.5281/zenodo.8105466](https://zenodo.org/records/8105466)** — confirmed by direct
    fetch to contain: processed Visium per-sample archives (13 GBM + 6 IDH-mutant, 12.6–49.9 MB
    each) plus H&E scans; **and** CODEX multiplex-protein OME-TIFF images (5 tumours' worth,
    718 MB–4.5 GB each file) plus two QuPath project archives (`CODEX_qupath_project_v2.zip`,
    2.2 GB; `CODEX_organoid_qupath_project.zip`, 1.8 GB) containing per-cell segmentation.
  - **Code**: [github.com/tiroshlab/Spatial_Glioma](https://github.com/tiroshlab/Spatial_Glioma),
    stated as publicly available.
- **Ingestible by `tme-sim/data/ingest.py` today?** **Partially, and unevenly.** The Visium layer
  is plausibly close to `ingest.py`'s native Space Ranger detector (file types match; exact
  directory shape not verified — this would need the same "download a small sample, inspect
  before committing" step `results/gbmspace_ingestion.md` used, not assumed). The CODEX layer is
  **not** ingestible without new code: it is OME-TIFF multiplex imagery plus a QuPath project, not
  a point-cloud/coordinate table `ingest.py` recognizes. This is exactly QuPath GeoJSON export —
  `results/landscape_pathology_3d.md`'s **#1-ranked ingestion target**, chosen there specifically
  because it is documented, scriptable, and has real public example files. **This dataset is a
  concrete reason to build that ingester**: it would not just be a generic pathology-tool
  adapter, it would directly unlock a real, GBM-specific, already-cross-registered second-modality
  section for `registration_3d.py`.

### 1.4 A third, independent GBM spatial atlas — same multi-region pattern, one more confirmation

Kwon et al. (Seoul National University group), *"A spatially resolved human glioblastoma atlas
reveals distinct cellular and molecular patterns of anatomical niches,"* *Nature Communications*
(Feb 2026), [PMID 41714633](https://pubmed.ncbi.nlm.nih.gov/41714633/),
[DOI 10.1038/s41467-026-69716-2](https://doi.org/10.1038/s41467-026-69716-2), PMC13031279 —
verified by direct full-text fetch. **28 patients** (18 GBM, 5 astrocytoma IDH-mutant, 1
oligodendroglioma IDH-mutant, 4 miscellaneous), **17 cases profiled by Visium CytAssist (13 GBM +
3 astrocytoma + 1 normal-brain control, 32 tissue slices total, 115,914 spots)**, plus **Xenium
on 4 GBM patients** with a custom 348-gene panel (729,307 cells). Quoted directly from the
Results: *"we analyzed patients across different age groups with diverse genomic alterations and
profiled ST on **at least two tissue sections prepared from different locations with distinct
anatomical features within the same tumor**."* This is, once again, **multi-region sampling, not
serial sectioning** — the same pattern as GBM-Space and the Greenwald/Tirosh paper, now confirmed
independently a third time across three unrelated research groups. **The paper itself explicitly
cites and cross-validates against Greenwald et al.'s CODEX data** (§1.3), which is a useful
secondary confirmation that the two datasets are compatible in spirit but does not add a new
serial-section resource. Data availability for this paper's own accessions was **not fully
confirmed in this pass** (the fetched full text did not surface a GEO/Zenodo accession number
before the token budget for that fetch was spent) — **UNVERIFIED**, flagged rather than guessed;
the paper's own interactive portal is referenced in-text but its URL was not captured.

### 1.5 Net finding for §1

Three independent, real, peer-reviewed GBM spatial-transcriptomics atlases were checked directly
against primary sources in this pass (GBM-Space via its already-known ingestion report plus a
fresh methods check; Greenwald/Tirosh 2024; Kwon et al. 2026). **All three use multi-region
sampling — different physical locations/blocks within one tumour — not continuous serial
sectioning of one block.** This is a structural fact about how GBM spatial atlases are built
(driven by neurosurgical sampling logistics: a piece from the core, a piece from the margin, a
piece from a contrast-enhancing region — not a bread-loafed serial series the way a resected
breast or pancreas specimen can be), not a gap this survey simply failed to find past. The one
genuine exception — near-adjacent CODEX/Visium pairs within the Greenwald/Tirosh dataset — is
real but is a same-block pair spanning **one adjacent step**, not a multi-slice stack, and is
cross-modality rather than same-modality.

---

## 2. Serial histology/IHC (non-CODA) — checked, and CODA's brain gap re-confirmed

**CODA + brain/GBM: re-confirmed, still no hits.** `results/landscape_pathology_3d.md` already
established, via a search independently verified by a second agent, that CODA (Kiemen et al.,
*Nat Methods* 2022, [PMID 36280719](https://pubmed.ncbi.nlm.nih.gov/36280719/)) has never been
applied to brain or GBM tissue. This pass re-ran that check fresh (not trusting the prior
finding blind, per the task's instruction) with a search specifically dated to 2026: no
publication surfaced applying CODA or a named derivative to brain/GBM tissue. The finding stands
unchanged.

**Other groups doing serial-section IHC/registration on GBM, without CODA specifically:** search
surfaced general serial-section 3D reconstruction methodology papers (e.g. a 3D primary-vessel
reconstruction framework from serial microscopy images; a fully-automatic serial-section
registration method, [PMC4598826](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4598826/)) but
**none of these were found applied to human GBM tissue specifically** — they are generic
methodology or applied to other organs (laryngeal cancer, mouse prostate). No group-specific,
non-CODA, serial-H&E-or-IHC 3D GBM reconstruction was found. This is consistent with — not a
contradiction of — the CODA finding: the gap is in the *tissue type*, not the *tool*, and it
holds across every tool checked, not just CODA.

---

## 3. 3D vascular network of human glioma — re-checked, still not public

Re-ran the two closest candidates from `results/landscape_pathology_3d.md` directly against
their own Data Availability Statements, fetched fresh in this pass rather than trusting the
prior survey's characterization unchecked:

- **Torres-Sánchez et al. 2021**, *"Three-dimensional vascular microenvironment landscape in
  human glioblastoma,"* *Acta Neuropathol Commun* 9:20,
  [PMID 33579378](https://pubmed.ncbi.nlm.nih.gov/33579378/),
  [PMC7879533](https://pmc.ncbi.nlm.nih.gov/articles/PMC7879533/). Direct fetch of the Data
  Availability Statement, verbatim: *"The datasets during and/or analysed during the current
  study are available from the corresponding author on reasonable request."* **Unchanged from
  the prior survey — still gated.** 3 human GBM specimens, confocal microscopy on 60 µm-thick
  sections.
- **FineVess (Chen et al.)**, deep-learning vessel extraction in 3D confocal images of cleared
  human glioma tissue, [PMID 38633068](https://pubmed.ncbi.nlm.nih.gov/38633068/),
  [PMC11019690](https://pmc.ncbi.nlm.nih.gov/articles/PMC11019690/). Direct fetch of the Data
  Availability Statement, verbatim: *"Data underlying the results presented in this paper are
  not publicly available at this time but may be obtained from the authors upon reasonable
  request."* **Unchanged — still gated, even as code (GitHub, PRBioimages/FineVess) is public.**
  12 human glioma specimens (right frontal cerebral region), CD31-immunolabeled, optically
  cleared (SBC delipidation + OPTIClear), confocal (Leica SP8) — note this is a **broader glioma
  cohort, grade not confirmed GBM-specific for all 12** in this pass; the task brief's own prior
  characterization ("including grade IV/GBM cases") is carried forward as previously established,
  not independently re-verified grade-by-grade here.

**No new public human-GBM 3D vascular network — graph, mask, or raw cleared-tissue volume —
surfaced in a fresh, independent search covering 2025–2026 literature.** The plain statement from
`results/landscape_pathology_3d.md` — *"no publicly downloadable 3D vascular network of human
glioblastoma anywhere"* — **holds, re-confirmed rather than merely repeated.**

---

## 4. Multi-region GBM data with real positional metadata — Ivy GAP has none; one dataset does

**Ivy GAP: confirmed, by direct inspection of the loader code and its own schema documentation
(`tmesim/data/real/ivygap.py`, `docs/real_data/ivygap.md`), to carry NO spatial x/y/z
coordinates.** Every join key in the dataset — `rna_well_id`, `tumor_id`, `block_id`,
`specimen_id`, `structure_abbreviation` — is categorical (which named anatomic structure, which
tumour, which sub-block) or a fractional-area composition (`sub_block_details.csv`'s 11
tumour-feature columns), never a coordinate. There is no MRI registration, no stereotactic
position, no within-slide pixel location published in the open Ivy GAP files this project
already uses. **Ivy GAP cannot support even a coarse 3D reconstruction on its own** — this
confirms and sharpens what `results/heldout_patient.md` and the census document already implied
(Ivy GAP is used in this project purely as an expression *signature* source per named region, per
`REGION_MAP` in `ivygap.py`) into an explicit, checked negative for the registration use case
specifically.

**One dataset does carry real 3D positional metadata, found fresh in this pass**: Mathur,
Wang, Schupp et al. (Costello lab, UCSF), *"Glioblastoma evolution and heterogeneity from a 3D
whole-tumor perspective,"* *Cell* 187(2):446-463 (2024),
[PMID 38242087](https://pubmed.ncbi.nlm.nih.gov/38242087/),
[DOI 10.1016/j.cell.2023.12.013](https://doi.org/10.1016/j.cell.2023.12.013), PMC10832360 —
verified by direct full-text fetch. **10 patients with newly diagnosed, treatment-naïve
IDH-wildtype GBM, 103 samples, each with a real, stated stereotactic (x, y, z) position**: quoted
directly, *"the location of each sample was recorded as a set of LPS coordinates registered to a
preoperative MR image using a Brainlab pointer and Brainlab Cranial Navigation software."* This
is genuinely the kind of positional metadata point 4 of the task asked about — real 3D
coordinates, not region labels — and it is explicitly designed to let a user visualize "360°"
spatial patterns across the whole tumour on **an interactive public platform** (referenced
in-text; exact URL not captured in the fetched excerpt).

**Caveats that matter for `registration_3d.py` specifically, stated plainly:**

1. **The underlying molecular data is bulk multi-omics per macrodissected biopsy** (WES, RNA-seq,
   ATAC-seq, Hi-C), **not** spatial transcriptomics/imaging — there is no per-cell or per-spot
   in-tissue coordinate the way Xenium/Visium provide, only one coordinate per whole biopsy
   sample. It cannot feed `register_serial_sections()`'s point-cloud-registration machinery
   directly; it is a different data shape entirely (sparse 3D point samples with rich molecular
   readouts per point, not two dense 2D tissue sections to align).
2. **Raw data is EGA-gated**, not open: *"All exome, RNA-Seq, ATAC-Seq, and Hi-C datasets have
   been deposited in the European Genome-Phenome Archive (EGA) database under accession number
   EGAS00001006785, EGAD00001005221/2, and EGAD00001009496/7."* Only the interactive
   visualization platform is public; the underlying coordinate + sequencing data requires EGA
   Data Access Committee approval, the same friction category `results/gbmspace_ingestion.md`
   already found blocks GBM-Space's raw CRAM deposit.
3. **Discrete, sparse biopsies (103 points across a whole tumour), not a continuous serial
   stack** — genuinely useful as an independent reference for what "real 3D GBM structure" looks
   like at a coarse scale (e.g., to sanity-check whether a future registration pipeline's implied
   tumour geometry is plausible), but not itself a registration target.

**Net for §4**: this is a real, valuable finding — the first dataset in this project's awareness
with genuine 3D positional GBM data beyond Ivy GAP's flat region labels — but it answers a
different question than the one `registration_3d.py` needs answered, and its raw form is
access-gated. Recorded here as a reference point and a candidate future outreach target (EGA DAC
application), not as an ingestible dataset today.

---

## 5. Glioma organoid/PDX + light-sheet/optical-clearing 3D imaging — no human hit found

Searched specifically for glioma organoid or PDX studies using light-sheet microscopy or optical
tissue clearing for 3D imaging, even without paired transcriptomics, per the task's explicit
fallback ("might at least provide a real 3D GBM tissue architecture reference to validate a
registration pipeline's geometry against, even without expression data").

**Found**: *"Integrated 3D Light-Sheet and 2D Multiplex Imaging for Deep Histological Profiling
of a Somatic Mouse Glioblastoma Model,"* bioRxiv
[2025.10.31.685751](https://www.biorxiv.org/content/10.1101/2025.10.31.685751v1.full) (2025) —
**mouse, not human**, and a preprint (not yet peer-reviewed as of this pass). This is the closest
match found in modality and recency, but fails the task's GBM-specificity requirement on species.

**Not found**: any human glioma organoid or PDX study using light-sheet microscopy or optical
clearing for 3D vasculature or tissue-architecture imaging. A general review,
*"3D Whole-Brain Imaging Approaches to Study Brain Tumors,"*
[PMC8071100](https://pmc.ncbi.nlm.nih.gov/articles/PMC8071100/), was fetched directly and
confirmed, in its own words, that clearing/light-sheet applications to human tumour tissue exist
for **other** cancers (breast, lung/lymph-node metastases) but that the brain-tumour-specific
examples in the review are **all mouse preclinical models**; the review states plainly that this
class of technology "is not yet widely used in the cancer field" for human application. The two
human-GBM cleared-tissue vasculature studies already found (Torres-Sánchez 2021, FineVess 2024,
§3) are the closest real matches to this category, and they remain the ones covered there — no
additional human GBM organoid/PDX light-sheet study was located.

**Net for §5**: confirmed negative, consistent with and extending
`results/landscape_pathology_3d.md`'s Part 5 finding that light-sheet/clearing pipelines for
human tumour vasculature are still boutique, single-lab efforts, and GBM specifically has even
less coverage than other cancers in this modality.

---

## Ranked summary table

| Candidate | GBM-specific? | Serial/consecutive same-block? | Public/downloadable today? | `ingest.py`-ready? |
|---|---|---|---|---|
| GSE324268 | one sample yes, one is an organoid | **No — 2 different specimens** | Yes, free | Yes (Xenium, already used) |
| GBM-Space (Sanger/EGA) | Yes (12 IDH-wt) | UNVERIFIED, likely multi-region not serial | Processed data yes; raw CRAM EGA-gated | **No — `.h5ad`, needs new reader** (`results/gbmspace_ingestion.md`) |
| Greenwald/Tirosh 2024 (Visium) | Yes (13 IDH-wt GBM of 19) | Multi-region (4 cases, multiple blocks/tumour) | Yes — GEO GSE237183, Zenodo 8105466 | Plausible, **not verified** against `resolve_visium_dirs()` |
| Greenwald/Tirosh 2024 (CODEX, near-adjacent to Visium) | Yes | **Near-adjacent, same block — closest real match found** | Yes — Zenodo 8105466 | **No — OME-TIFF + QuPath project, needs a QuPath/CODEX reader** |
| Kwon et al. 2026 (Nat Commun) | Yes (13 GBM Visium, 4 GBM Xenium) | Multi-region ("at least two... different locations") | UNVERIFIED accession | UNVERIFIED |
| Ivy GAP | Yes | No coordinates at all, categorical only | Yes, fully open | No spatial coordinate to ingest |
| Mathur et al. 2024 (Costello/UCSF) | Yes (10 IDH-wt) | Discrete stereotactic biopsies, real (x,y,z), not serial | Interactive viewer yes; raw data EGA-gated | No — bulk multi-omics, wrong data shape |
| Torres-Sánchez 2021 | Yes (3 specimens) | N/A (single-block confocal volumes) | **No — reasonable request only** | N/A |
| FineVess 2024 | Broader glioma (grade unconfirmed per-specimen) | N/A | **No — reasonable request only** | N/A |
| Mouse GBM light-sheet (bioRxiv 2025) | **No — mouse** | N/A | Preprint, species disqualifies | N/A |

---

## Recommendation: next concrete step for `registration_3d.py`

**There is no drop-in GBM equivalent of the breast Xenium Rep1/Rep2 pair.** The honest,
unhedged answer to the task's framing question is: **not yet, and the reason is structural, not
a search failure** — every GBM spatial atlas checked samples multiple *regions* of a tumour by
surgical necessity, not multiple *serial slices* of one block, because GBM tissue does not come
off the table as a bread-loafable specimen the way a breast lumpectomy or a pancreas resection
does.

**The concrete next target, in priority order:**

1. **Attempt the Greenwald/Tirosh 2024 Visium layer first** (GSE237183) — same modality
   `ingest.py` already reads, GBM-specific (13 IDH-wt sections), zero new ingestion code *if* the
   directory layout matches (unverified — the next step is exactly the "download a small sample,
   inspect before committing" move `results/gbmspace_ingestion.md` already modeled). Even a
   single pair of GBM Visium sections from different regions of one tumour, registered the way
   `registration_3d.py` already does for tissue-geometry alignment, would be a genuine
   GBM-specific rerun — with the honest caveat that it is a cross-region, not cross-serial-section,
   alignment, and that distinction should be stated as plainly in any future report as this one
   states it here.
2. **Build the QuPath ingester `results/landscape_pathology_3d.md` already recommends as its #1
   target, then use it to pull the Greenwald/Tirosh CODEX layer** (Zenodo 8105466) and register it
   against its near-adjacent Visium section using the same registration method that dataset's own
   authors already validated with STalign. This is more work (a genuinely new ingester) but is
   the single most on-thesis target found in this entire survey: real GBM tissue, near-adjacent
   sections of the same block, cross-modality registration already proven to work by an
   independent peer-reviewed group.
3. **Email FineVess and/or Torres-Sánchez's corresponding authors**, per
   `results/landscape_pathology_3d.md`'s existing recommendation — unchanged by this pass, still
   the only route to real human GBM 3D vasculature geometry, gated behind "reasonable request"
   rather than a technical barrier.
4. **Apply for EGA DAC approval** for either GBM-Space's raw CRAM deposit or the Mathur et al.
   2024 3D-coordinate multi-omics deposit, in parallel, given both are slow — not a near-term
   registration target but a longer-term one worth starting now if this project wants real GBM 3D
   structure beyond what near-adjacent sections alone can offer.

**What would change this verdict**: any lab willing to hand over two or more consecutive,
same-block GBM sections (any modality) at a stated or photographable spacing — the single
concrete ask this survey's negative result points toward, and exactly the kind of thing worth
raising the next time this project is in outreach mode (`docs/outreach/`).
