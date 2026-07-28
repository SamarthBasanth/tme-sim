# Landscape census: public spatial transcriptomics data across solid tumour types

**Purpose.** tme-sim currently reconstructs a tumour microenvironment (TME) from spatial
transcriptomics (10x Xenium/Visium) and runs drug-transport physics on it — glioblastoma (GBM)
only. The user asked directly whether it is worth pivoting away from GBM to a cancer type with
richer public spatial data. This document builds a comparative census across 8 major solid tumour
types on data availability, panel-marker coverage, multi-region/serial-section availability,
outcome-data pairing potential, and existing computational prior art — then makes an explicit,
unhedged recommendation.

**Method.** GEO (`ncbi.nlm.nih.gov/gds`) was queried with an identical pattern —
`"<cancer type> AND spatial transcriptomics"` — for all 8 cancer types on 2026-07-28, and the
**Series** count from the results-page facet is used as a cross-cancer, apples-to-apples,
order-of-magnitude proxy for dataset volume. This is a **noisy proxy, not an exact platform
count**: GEO's search matches any series whose text mentions both terms, so it includes some
false positives (e.g. bulk deconvolution papers that use the word "spatial" loosely) and excludes
some real spatial datasets that don't use the exact search terms. It is, however, the same
method applied identically to all 8 cancers, which is what a comparative census needs. Vendor
public-dataset portals (10x Genomics, Bruker/NanoString, Vizgen), landmark papers, and PubMed
(via the PubMed MCP tool) were used to sanity-check and add texture beyond the raw GEO count.
Every specific figure below is cited; where a figure could not be pinned to a primary source it
is marked **UNVERIFIED**. This document does not modify any code and does not touch the paired
drug-response-data question — that was answered separately and exhaustively in
`results/DATA_HUNT.md` (headline: exactly one open dataset anywhere pairs a measured spatial drug
distribution with spatially co-located molecular tissue state, and it is ovarian, not glioma —
see §7 below for how that result folds into this one).

Today's date for currency purposes: 2026-07-28.

---

## 1. The true glioma/GBM baseline (not just what's in this repo)

Before comparing cancers, it's worth being honest about what this repo has actually pulled
versus what exists. **`data_cache/` in this repo currently holds 4 distinct spatial glioma
datasets**: the 10x Xenium Human Brain Glioblastoma panel (319 genes), the 10x Xenium Brain
Cancer GBM IO panel (480 genes), one Visium GBM dataset, and one GEO deposit (GSE324268, 366
genes) — plus non-spatial supporting data (TCGA-GBM clinical/expression, Ivy GAP anatomic
RNA-seq). That is a handful, and it is not close to the ceiling of what's public.

**True baseline, GEO search** (`ncbi.nlm.nih.gov/gds/?term=glioma+AND+spatial+transcriptomics`,
2026-07-28): **94 Series / 182 total items** (2 DataSets, 94 Series, 86 Samples). Narrowing to
just "glioblastoma" as the search term: **64 Series / 139 total items** (2 DataSets, 64 Series,
73 Samples) — glioma is the broader, more accurate search since a meaningful share of the
spatial literature covers diffuse/low-grade glioma and pediatric glioma alongside GBM specifically.

That GEO count is a floor, not a ceiling, once landmark atlas papers are counted individually:

- **GBM-Space**: 97 Visium sections across 12 GBM tumours, 4–15 sampling sites per tumour, 1.35M
  nuclei/spot integration ([EGA EGAD00001015527](https://ega-archive.org/datasets/EGAD00001015527); [gbmspace.org](https://www.gbmspace.org/)).
- A late-2025 multi-cohort integration combined **44 GBM Visium slides / 128,176 spots** across
  three studies (the paper's own 14 slides from 7 patients + the Broad and Freiburg cohorts) —
  "one of the largest GBM ST datasets analyzed to date" per the authors
  ([bioRxiv 2025.10.17.683110](https://www.biorxiv.org/content/10.1101/2025.10.17.683110v2.full)).
- A February 2026 atlas paper reports Visium CytAssist + Xenium validation on an independent GBM
  cohort ([Nature Communications, PubMed 41714633](https://pubmed.ncbi.nlm.nih.gov/41714633/); [DOI 10.1038/s41467-026-69716-2](https://www.nature.com/articles/s41467-026-69716-2)).
- The 10x Genomics public Xenium dataset portal alone carries **16 public Xenium datasets** in
  the "human brain (cancer)" category (part of a larger 8M-cell, 12-tissue Xenium corpus)
  per a 2026 benchmarking preprint that inventoried the portal
  ([DeepSpaceDB 2.0, bioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.15.699623.full.pdf)) — **UNVERIFIED** exactly how many of the 16 are GBM specifically vs. other brain-cancer entities, but the category is GBM-dominated in practice.

**Conclusion for §1**: glioma/GBM's real public spatial-data baseline is on the order of
**several dozen to ~100 independent datasets/series**, several with 40–100+ sections each — this
repo has ingested roughly 4 of them. The gap between "true baseline" and "what tme-sim has
downloaded" is itself the single most actionable, lowest-risk finding in this document,
independent of any pivot decision (see §8).

---

## 2. Comparison table

| Cancer type | GEO series count (2026-07-28, `<type> AND spatial transcriptomics`) | Panel-native drug-transport markers | Cancer-specific resistance determinant on a stock commercial panel? | Multi-region / serial-section data | TCGA/cBioPortal-adjacent outcome data | Prior computational TME/drug-transport modelling |
|---|---|---|---|---|---|---|
| **Glioma/GBM** (current focus) | **94** (`glioma`) / 64 (`glioblastoma`) [1] | PECAM1 ✅, MKI67 ✅ present; VEGFA/CA9/HIF1A/ABCB1/ABCG2/COL1A1 **absent** — repo-verified against the actual panel JSON [2] | **NO** — MGMT absent from both stock 10x GBM Xenium panels (319-gene, 480-gene); found only on one bespoke 366-gene GEO panel [2] | Yes, extensively — GBM-Space (97 sections/12 tumours), core-to-edge (multi-site), combined 44-slide atlas [1][3] | TCGA-GBM + TCGA-LGG, order of several hundred patients each, treatment/outcome fields via cBioPortal [4] | **Deepest and oldest field of any cancer here** — Swanson/Alvord/Murray 2002 founded a lineage of BBB-specific reaction-diffusion drug-delivery models still active today [5][6] |
| **Breast** | **100** (capped display; 428 total items) [1] | ERBB2/PGR/ESR1/GATA3/MKI67 confirmed on the dedicated 280-gene Xenium Breast Panel [7]; PECAM1/hypoxia/efflux status **UNVERIFIED** gene-by-gene | **YES** — ER/PR/HER2 (ESR1/PGR/ERBB2), breast's actual clinical resistance/subtyping axis, ship natively on the standard panel [7] | Yes — Janesick et al. matched Visium+Xenium reference dataset [8]; Bassiouni et al. 28 sections/14 TNBC tumours (GSE210616) | TCGA-BRCA (~1,000+ patients) **plus** METABRIC (n≈2,509, rich treatment+survival, no matched spatial) [4] | QSP models exist (HER2+ mBC, TNBC-TAM immunotherapy) [9] but are PK/PD- and immune-focused, not transport-geometry-focused; a tissue-engineered 3D microvessel TME model exists [10] |
| **Lung / NSCLC** | 63 (`lung cancer`) / **23** (`non-small cell lung cancer`, narrower term) [1] | Xenium Human Lung Panel + Multi-Tissue and Cancer panel (377 genes) cover "select immune, proliferation, tumor markers" [11]; **EGFR/ALK/KRAS presence UNVERIFIED** — not confirmed gene-by-gene, plausibly need the +100 custom add-on | **UNVERIFIED** | Yes, indirectly — TRACERx/PEACE: 501 longitudinal multi-region primary+metastatic samples, 24 NSCLC patients (genomic, not itself spatial-transcriptomic) [12] | TCGA-LUAD + TCGA-LUSC (~500 patients each) [4] | QSP platforms for checkpoint inhibitors exist in NSCLC [9]; no transport-physics-specific literature surfaced in this pass |
| **Colorectal** | **66** [1] | Uses Multi-Tissue/Cancer panel (colon is one of the 9 covered tissues) [11]; KRAS status **UNVERIFIED** | **UNVERIFIED** | Cross-platform benchmark used colon adenocarcinoma as one of its 3 tissue types (Stereo-seq/Visium HD/CosMx 6K/Xenium 5K vs. CODEX/scRNA-seq, same serial sections) [13]; GSE303070 alone carries 352 samples on "stromal niches" | TCGA-COAD/READ (~450–630 combined) [4] | Liver-metastasis-specific perfusion/PK models exist (mechanistic diffusion model estimating tumour-site chemo concentration from patient data) [14]; the classic "tumour cord" transport model is generic, not colorectal-specific [15] |
| **Pancreatic (PDAC)** | **54** [1] | Multi-Tissue/Cancer panel used for the 10x public pancreatic Xenium dataset [16]; KRAS (PDAC's ~90%-penetrant driver) status **UNVERIFIED** | **UNVERIFIED** | **Yes, and growing fast** — a 2026 "mega-block" study profiled 39 treatment-naive patients with full tumour cross-sections: Visium (340K spots) + CosMx 969-gene (530K single cells) [17]; a separate integrated atlas adds scRNA n=88 + Visium n=20/67,933 spots + Xenium n=2/307,679 cells [18] | TCGA-PAAD (~150–185, smallest major TCGA cohort) [4] | **Second-deepest delivery-specific field after GBM** — dense desmoplastic stroma + hypovascularity + high interstitial pressure is PDAC's own well-documented physical delivery barrier, with active mechanistic modelling (ultrasound-cavitation-enhanced liposomal delivery [19], a Phase I focused-ultrasound drug-delivery trial [20], imaging-biomarker perfusion studies showing gemcitabine-induced hypo-vascularization [21]) |
| **Melanoma** | **73** [1] | Uses Multi-Tissue/Cancer panel for the 10x public skin dataset [11]; BRAF/NRAS status **UNVERIFIED** | **UNVERIFIED** | Yes — melanoma brain metastasis: 20 sequential 20 µm Visium sections across 4 regions plus multi-region bulk sequencing [22]; meningeal-stroma leptomeningeal-disease spatial paper [23] | TCGA-SKCM (~410–470) [4] | Thin and generic — nanoparticle transvascular-transport simulations exist but are not melanoma-specific [24]; the melanoma modelling literature is dominated by immune-checkpoint/immunotherapy QSP, not transport physics |
| **Ovarian (esp. HGSOC)** | **30** [1] | 10x public dataset uses the broad Xenium Prime 5K Pan Tissue and Pathways panel (5,001 genes) [25] — broad panels capture hypoxia/efflux genes "for free" per the general platform-landscape finding [26]; BRCA1/2 (HGSOC's PARPi-resistance determinant) status **UNVERIFIED** on this specific panel | **UNVERIFIED**, but plausible given the panel's breadth | Yes — Moncayo et al. 2026: MALDI-MSI drug maps + GeoMx WTA across 12 serial slices, 3 patients, explicit serial-section protocol (section 0/+1/−2/−1/+2/+3) [27] | TCGA-OV (~300–585 depending on data freeze) [4] | Thin computationally, but **uniquely strong empirically**: Moncayo is the one dataset in the entire open-data landscape (per `results/DATA_HUNT.md`) that pairs a measured per-pixel spatial drug distribution with spatially co-located transcriptomic response signal, anywhere, in any cancer [27] |
| **HNSCC** | **18** — smallest of the 8 [1] | No dedicated panel found; likely Multi-Tissue/Cancer or custom; EGFR status **UNVERIFIED** | **UNVERIFIED** | No landmark large multi-region public HNSCC spatial atlas surfaced in this pass — **UNVERIFIED**, likely genuinely thinner than the other 7 | TCGA-HNSC (~500–530) [4] | A computational TME model exists but is immune-subtype/checkpoint-response focused, not transport-physics [28]; a topical-drug-delivery tissue-engineered mucosa model exists but is a lab construct, not a computational transport simulator on patient spatial data [29] |

**Ranking by raw GEO series count** (order-of-magnitude, noisy proxy): Breast (100) ≈ Glioma (94)
> Melanoma (73) > Colorectal (66) > Lung-broad (63) ≈ Glioblastoma-narrow (64) > Pancreatic (54)
> Ovarian (30) > NSCLC-narrow (23) > HNSCC (18).

**This ranking, on its own, would say "pivot to breast or stay roughly where glioma already
sits, and pancreatic/ovarian/HNSCC are the weak options."** §3 explains why raw dataset count is
the wrong single axis for a transport-physics model specifically.

---

## 3. The central tension: data volume vs. scientific need for a transport model

This is the tension the task explicitly asked to be addressed head-on, not glossed over.

**Glioblastoma's transport physics is unusually well-motivated because the blood-brain barrier
makes drug delivery THE dominant clinical problem** — not one factor among several, but in a
real sense the rate-limiting one. This is not an assertion invented for this document; it is the
premise of an entire, decades-old sub-field of mathematical oncology that exists specifically
because of it. Swanson, Alvord & Murray's 2002 model of "quantifying efficacy of chemotherapy of
brain tumors with homogeneous and heterogeneous drug delivery" (PMID 12675529,
[DOI 10.1023/a:1022644031905](https://doi.org/10.1023/a:1022644031905)) is foundational, and the lineage is still active: Boujelben et al. 2016
explicitly model "blood flow, vascular permeability and diffusion within the tumour
microenvironment" for GBM because "vasculature within these tumours is typically abnormal... and
they typically exhibit a disrupted blood-brain barrier" (PMID 27708763, [DOI 10.1098/rsfs.2016.0039](https://doi.org/10.1098/rsfs.2016.0039)).
No other cancer in this census has a comparably deep, comparably delivery-focused modelling
literature — because in most of the other 7, the dominant driver of treatment response is
*intrinsic tumour-cell sensitivity* (a receptor status, a driver mutation, an immune
microenvironment), not whether the drug physically reaches the cell.

**Breast is the sharpest illustration of the opposite pole.** It has the largest or
second-largest raw spatial-data volume of the 8 (§2), and its stock commercial Xenium panel is
the one confirmed in this census to carry its own flagship resistance determinant natively
(ESR1/PGR/ERBB2 [7] — a genuine structural advantage over GBM, where MGMT is absent from both
stock panels and had to be sourced from one bespoke GEO deposit [2]). But breast tumours are, in
the median case, well-vascularized, and the dominant axis of "does this drug work" is receptor
status (HER2+/ER+/TNBC), not whether the drug penetrates the tissue. That is exactly why the
breast prior-art column in §2 is QSP/PK-PD work, not transport-geometry work [9][10] — the field
itself has already made the call that delivery physics is not where the leverage is for breast.
A transport-physics simulator built on breast data would be data-rich and scientifically
thin — precisely the failure mode this project's own honest-null result already flagged
(`results/tme-sim-stage0.md`: the synthetic kill-test never beats a cheap non-spatial baseline
when response is sensitivity-dominated rather than delivery-dominated).

**Pancreatic cancer (PDAC) is the one cancer type in this census that plausibly matches GBM's
delivery-dominance story on independent physical grounds, not by analogy.** PDAC's dense
desmoplastic stroma, near-uniform hypovascularity, and characteristically high interstitial
fluid pressure are themselves a well-documented, mechanistically distinct physical barrier
([Heid et al. 2022, gemcitabine-induced hypo-vascularization measured by imaging biomarkers,
PMID 36074156](https://doi.org/10.1007/s00259-022-05930-6); active Phase I ultrasound-mediated
delivery trials specifically because the stroma is the acknowledged bottleneck [20]). It is not
as data-rich as breast/melanoma/glioma by raw GEO count (54, §2), but it is the only cancer in
this census with **both** a real, independently-documented delivery-dominant biology **and** a
fast-growing, large, genuinely multi-region public dataset — the 2026 mega-block cohort alone
(39 patients, Visium 340K + CosMx 530K transcriptomes, full tumour cross-sections [17]) is
comparable in scale to anything glioma has, and it is recent enough that it likely postdates
whatever the last GBM-focused literature search assumed as the ceiling.

**Ovarian (HGSOC) does not have PDAC's or GBM's mechanistic delivery-dominance story as clearly**
(peritoneal spread and ascites give it real distribution complexity, but this is far less
established in the literature as *the* dominant response bottleneck than the BBB is for GBM or
stroma is for PDAC), **but it has something no other cancer in this census has: an actual
measured spatial drug-distribution dataset** (Moncayo et al. — MALDI-MSI drug maps
co-registered with GeoMx transcriptomic response across serial sections, 3 patients, 12 slices
[27], already the subject of exhaustive verification in `results/DATA_HUNT.md`). That is not a
reason to make ovarian the primary flagship, but it is a genuine, unique validation asset worth
keeping on the radar independent of the main recommendation.

**Lung, colorectal, melanoma, and HNSCC** each have real data (§2) but none of them combines (a)
a documented, mechanistically-established delivery-dominant biology and (b) a strong or growing
public spatial corpus, in the way GBM and PDAC do. Colorectal's most delivery-relevant story
(liver metastasis perfusion) is a real niche but is about metastatic-site perfusion, not primary
tumour transport, and already has its own dedicated modelling niche (§2) separate from a
Visium/Xenium-ingest use case. Melanoma and lung are dominated by immune-checkpoint and
driver-mutation biology respectively. HNSCC is simply the thinnest cancer type in this census on
every axis measured (§2) — smallest GEO count, no located multi-region atlas, thinnest prior art.

---

## 4. What this means for outcome-data validation potential (point 4, all 8 types)

Every cancer type in this census has a TCGA cohort accessible via cBioPortal/GDC with
treatment-history and survival-outcome fields — this is not a differentiator between the 8
([cBioPortal](https://www.cbioportal.org/); [GDC PanCanAtlas](https://gdc.cancer.gov/about-data/publications/pancanatlas) [4]). Cohort sizes are broadly comparable in
order of magnitude (roughly 150–1,100+ patients depending on cancer type and TCGA data-freeze
version — exact counts vary meaningfully by source and are **UNVERIFIED** beyond order of
magnitude in this pass): breast (largest, ~1,000+) and pancreatic (smallest major cohort,
~150–185) sit at the two ends, with lung, colorectal, melanoma, ovarian, HNSCC, and glioma/LGG
combined in the several-hundred-patient middle. Breast additionally has METABRIC (n≈2,509, rich
treatment+survival, no matched spatial data) as an outsized extra resource [4]. **None of this
changes the ranking in §3** — as `results/DATA_HUNT.md` established exhaustively, TCGA-style
outcome cohorts are never the *same specimens* as any spatial dataset for any cancer type
surveyed there or here; this is about future-validation potential for the field broadly, and on
that axis all 8 cancers are roughly equally well-positioned.

---

## 5. Sources

[1] GEO DataSets search, `ncbi.nlm.nih.gov/gds/?term=<query>`, accessed 2026-07-28; each URL
follows the pattern `https://www.ncbi.nlm.nih.gov/gds/?term=<cancer+type>+AND+spatial+transcriptomics`.

[2] `results/flagship_markers.md` and `results/landscape_spatial_platforms.md` (this repo) —
direct inspection of `data_cache/xenium_gbm/Xenium_GBM_gene_panel.json` (319 genes) and the
480-gene GBM IO panel; MGMT/hypoxia/efflux gene absence verified against the actual cached panel
files, not inferred.

[3] [EGA EGAD00001015527 / GBM-Space](https://ega-archive.org/datasets/EGAD00001015527); [gbmspace.org](https://www.gbmspace.org/); [bioRxiv 2025.04.04.647096](https://www.biorxiv.org/content/10.1101/2025.04.04.647096v2) (core-to-edge mapping); [bioRxiv 2025.10.17.683110](https://www.biorxiv.org/content/10.1101/2025.10.17.683110v2.full) (44-slide combined atlas); [PubMed 41714633](https://pubmed.ncbi.nlm.nih.gov/41714633/) / [Nat Commun DOI 10.1038/s41467-026-69716-2](https://www.nature.com/articles/s41467-026-69716-2); [DeepSpaceDB 2.0 bioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.15.699623.full.pdf) (16 public Xenium brain-cancer datasets).

[4] [cBioPortal](https://www.cbioportal.org/); [GDC PanCanAtlas publications](https://gdc.cancer.gov/about-data/publications/pancanatlas); cohort sizes per-cancer are approximate and vary by TCGA data-freeze/version across sources consulted — treated as order-of-magnitude only.

[5] Swanson KR, Alvord EC, Murray JD. "Quantifying efficacy of chemotherapy of brain tumors with
homogeneous and heterogeneous drug delivery." *Acta Biotheoretica* 2002;50(4):223-37. PMID
12675529, [DOI 10.1023/a:1022644031905](https://doi.org/10.1023/a:1022644031905). (Retrieved via PubMed MCP; per PubMed attribution requirement.)

[6] Boujelben A, Watson M, McDougall S, et al. "Multimodality imaging and mathematical modelling
of drug delivery to glioblastomas." *Interface Focus* 2016;6(5):20160039. PMID 27708763, [DOI 10.1098/rsfs.2016.0039](https://doi.org/10.1098/rsfs.2016.0039). (Retrieved via PubMed MCP.)

[7] [10x Genomics Xenium Human Breast Panel documentation](https://www.10xgenomics.com/products/xenium-panels) and product/dataset pages surfaced via web search (280-gene panel; ERBB2/PGR/ESR1/GATA3/MKI67 confirmed present per 10x's own panel description).

[8] Janesick et al., "High resolution mapping of the breast cancer tumor microenvironment using
integrated single cell, spatial and in situ analysis of FFPE tissue," [bioRxiv 2022.10.06.510405](https://www.biorxiv.org/content/10.1101/2022.10.06.510405.full.pdf); Bassiouni et al. 2023, GEO [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616) (28 sections/14 TNBC tumours, per search-surfaced summary — **UNVERIFIED** against the primary paper directly in this pass).

[9] Tumor-associated macrophage QSP model, TNBC immunotherapy: [PMC9287616](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9287616/) / PMID 35856032.

[10] 3D tissue-engineered microvessel TME drug-delivery model: PMID 33405827.

[11] [10x Genomics Xenium Human Multi-Tissue and Cancer Panel](https://www.10xgenomics.com/products/xenium-v1-panel) (377 genes, covers breast/lung/colon/pancreas/skin/etc. cell-typing; used for the 10x public pancreatic-cancer and skin/lung datasets per their respective dataset pages).

[12] "Evolutionary characterization of lung cancer metastasis," *Nature* 2026, [DOI 10.1038/s41586-026-10428-4](https://www.nature.com/articles/s41586-026-10428-4) (TRACERx/PEACE, 501 multi-region primary+metastatic samples, 24 NSCLC patients; genomic, not itself a spatial-transcriptomics dataset).

[13] Ren et al., "Systematic benchmarking of high-throughput subcellular spatial transcriptomics
platforms across human tumors," *Nat Commun* 2025, [DOI 10.1038/s41467-025-64292-3](https://doi.org/10.1038/s41467-025-64292-3), [PMC12534522](https://pmc.ncbi.nlm.nih.gov/articles/PMC12534522/) (also cited in `results/landscape_spatial_platforms.md` §6); GEO [GSE303070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303070) (352 samples, colorectal-inclusive "stromal niches" study, surfaced via GEO search).

[14] Mechanistic diffusion model of chemotherapy concentration at colorectal liver-metastasis
tumour sites from patient-derived data: [PMC7866038](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7866038/).

[15] Generic "tumour cord" drug-delivery computational simulation: [PMC5451806](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5451806/).

[16] [10x Genomics: Pancreatic Cancer with Xenium Human Multi-Tissue and Cancer Panel](https://www.10xgenomics.com/datasets/pancreatic-cancer-with-xenium-human-multi-tissue-and-cancer-panel-1-standard).

[17] "In situ multi-modal characterization of pancreatic cancer reveals tumor cell identity as a
defining factor of the surrounding microenvironment," [ScienceDirect / Cell family journal, S2211124725015992](https://www.sciencedirect.com/science/article/pii/S2211124725015992) — 39 treatment-naive patients, "mega-block" (~5×8 cm) full tumour cross-sections, Visium (~340K spots) + CosMx 969-gene panel (~530K single-cell transcriptomes); publication dated per search metadata to Jan 2026, journal identity **UNVERIFIED** beyond the ScienceDirect listing in this pass.

[18] Integrated PDAC atlas (scRNA n=88/187,520 cells; Visium n=20/67,933 spots; bulk RNA n=1,383;
Xenium n=2/307,679 cells): [bioRxiv 2025.02.11.637627](https://www.biorxiv.org/content/10.1101/2025.02.11.637627.full.pdf).

[19] "CavitoMod-UTMDNet: A Mechanistic Cavitation-Diffusion Framework for Ultrasound-Targeted
Microbubble Destruction-Enhanced Liposomal Doxorubicin Delivery in Pancreatic Cancer." PMID
41969145, [DOI 10.1177/10849785261431195](https://doi.org/10.1177/10849785261431195). (Retrieved via PubMed MCP.)

[20] "PanDox: a phase I study of targeted chemotherapy delivery to non-resectable primary
pancreatic tumours using thermosensitive liposomal doxorubicin and focused ultrasound." PMID
37741968, [DOI 10.1186/s12885-023-11228-z](https://doi.org/10.1186/s12885-023-11228-z). (Retrieved via PubMed MCP; ClinicalTrials.gov NCT04852367.)

[21] Heid I, Trajkovic-Arsic M, Lohöfer F, et al. "Functional biomarkers derived from CT and MRI
differentiate PDAC subgroups and reveal gemcitabine-induced hypo-vascularization." *Eur J Nucl
Med Mol Imaging* 2022. PMID 36074156, [DOI 10.1007/s00259-022-05930-6](https://doi.org/10.1007/s00259-022-05930-6). (Retrieved via PubMed MCP.)

[22] "In-depth patient-specific analysis of tumor heterogeneity in melanoma brain metastasis:
Insights from spatial transcriptomics and multi-region bulk sequencing," [ScienceDirect S1936523325001998](https://www.sciencedirect.com/science/article/pii/S1936523325001998) / [PMC12284558](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12284558/) — 20 sequential 20 µm sections across 4 regions.

[23] "Spatial transcriptomics analysis identifies a tumor-promoting function of the meningeal
stroma in melanoma leptomeningeal disease," [PMC11228800](https://pmc.ncbi.nlm.nih.gov/articles/PMC11228800/).

[24] "Simulation of transvascular transport of nanoparticles in tumor microenvironments for drug
delivery applications," *Sci Rep* 2024, [PMC10798967](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10798967/) / [Nature DOI 10.1038/s41598-024-52292-0](https://www.nature.com/articles/s41598-024-52292-0) — generic tumour model, not melanoma-specific.

[25] [10x Genomics: FFPE Human Ovarian Cancer with 5K Human Pan Tissue and Pathways Panel plus 100 Custom Genes](https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-ovarian-cancer).

[26] `results/landscape_spatial_platforms.md` §1 Q3 (this repo) — established that
whole-transcriptome/very-high-plex platforms (Visium, Stereo-seq, Xenium Prime 5K, CosMx WTx)
capture vasculature/hypoxia/efflux genes "for free," while targeted panels are hit-or-miss and
must be checked gene-by-gene; this general finding is applied here rather than re-verified
gene-by-gene for every cancer type in this pass, given time constraints.

[27] Moncayo et al., *Nature Communications* 2026 — PARP inhibitors in HGSOC patient-derived
explants, MALDI-MSI + GeoMx WTA. Zenodo [17610220](https://zenodo.org/records/17610220); GEO [GSE281519](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281519); PMID 41844641, [DOI 10.1038/s41467-026-70558-1](https://doi.org/10.1038/s41467-026-70558-1). Full independent verification (376,736-pixel MSI table, ROI-level drug×expression pairing, specificity checks) already done in `results/DATA_HUNT.md` §2A T1.1 — not re-verified here, cited by reference.

[28] "Tumor microenvironment governs the prognostic landscape of immunotherapy for head and neck
squamous cell carcinoma: A computational model-guided analysis," *PLOS Comput Biol*,
[PMC12162103](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12162103/).

[29] "A multilayered epithelial mucosa model of head neck squamous cell carcinoma for analysis of
tumor-microenvironment interactions and drug development," [ScienceDirect S0142961220305238](https://www.sciencedirect.com/science/article/abs/pii/S0142961220305238) — tissue-engineered construct, not a computational transport simulator on patient spatial data.

CosMx NSCLC FFPE public dataset (8 samples/5 patients, ~800K cells, 960-gene panel):
[Bruker Spatial Biology](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/nsclc-ffpe-dataset/), already documented in `results/landscape_spatial_platforms.md`.

---

## 6. Ranked recommendation

**STAY ON GLIOBLASTOMA — do not pivot — with PDAC flagged as the one strong candidate to ADD as
a second flagship if and when capacity exists, ranked strictly behind finishing the GBM data
build-out below.**

**Justification, both data availability and scientific fit, not data availability alone:**

1. **Scientific fit is not close.** GBM is the only cancer type in this census with a decades-old,
   still-active, delivery-focused mathematical-modelling literature that exists *because* the
   BBB makes transport the dominant clinical bottleneck (§3, [5][6]). That is the actual thesis
   tme-sim is built on, and no other cancer type in this census matches it as cleanly. PDAC is
   the only serious second candidate on this axis specifically, not breast/melanoma/lung despite
   their larger raw datasets (§3).

2. **Data availability was never actually the constraint.** §1 shows the true glioma/GBM public
   baseline is on the order of 60–100+ independent datasets/series, several with 40–100+ sections
   each — this repo has ingested roughly 4. Before concluding "GBM doesn't have enough data,"
   the cheap, low-risk, zero-pivot-cost move is to actually pull GBM-Space (97 sections, 12
   tumours) and/or the 44-slide combined Visium atlas [3] into `data_cache/` and re-run the
   existing pipeline against them. That is strictly higher-leverage than a cancer-type pivot: it
   tests the real bottleneck (does more/richer GBM data change the honest-null result already
   found in `results/tme-sim-stage0.md`?) without discarding the one asset (the BBB-delivery
   rationale) that makes the physics worth building at all.

3. **Breast (and to a lesser extent melanoma), despite the largest datasets, are the wrong
   pivot targets precisely because they have the least need for a transport model.** §3 argues
   this is not merely a data-volume tradeoff but a structural one: breast response is
   receptor-status-dominated, and the field's own prior art (§2) has already concentrated on
   PK/PD and QSP, not transport geometry, because that is where the real signal is. Building a
   transport-physics simulator there would very likely reproduce this project's own prior finding
   (`results/DRUG_SELECTION_RIGOR.md`) that a cheap feature-based classifier beats the mechanistic
   physics whenever the ground truth is sensitivity-dominated rather than delivery-dominated —
   except worse, because breast's ground truth (receptor status) is *known and reliably
   measured*, unlike GBM's synthetic ground truths, so the cheap baseline would be even harder to
   beat there.

4. **PDAC is the one add-a-second-flagship candidate that earns its place on both axes
   simultaneously**, not on data volume alone: a real, independently-documented physical delivery
   barrier (dense stroma, hypovascularity, high IFP — mechanistically distinct from the BBB, but
   analogous in being the acknowledged dominant bottleneck [19][20][21]), plus a large and rapidly
   growing multi-region public dataset that did not exist in this form even a year prior to this
   census (the 39-patient mega-block cohort [17]). It should be evaluated as a **second flagship
   alongside GBM**, not a replacement — the two diseases would stress-test the same transport
   physics against two structurally different barriers (vascular/BBB vs. stromal/ECM), which is
   scientifically more valuable than doubling down on GBM alone or abandoning it for a
   biomarker-dominated cancer.

5. **Ovarian (HGSOC) is worth a standing watch, not a pivot.** It is the only cancer in this
   entire census — and per `results/DATA_HUNT.md`, the only cancer type in the *entire open-data
   landscape searched to date* — with a real measured spatial drug-distribution dataset (Moncayo
   [27]). That is a unique validation asset. But its delivery-dominance story is weaker and less
   established than GBM's or PDAC's, and the existing dataset is ex vivo/bath-dosed with no
   vascular delivery component (per `results/DATA_HUNT.md`'s own caveats) — so it is a validation
   target worth revisiting if/when tme-sim needs a real ground-truth drug-distribution test, not a
   flagship-cancer candidate today.

6. **Lung, colorectal, melanoma, HNSCC are not recommended at any priority.** None combines a
   documented delivery-dominant biology with a standout data position (§3); HNSCC in particular
   is the weakest cancer type in this census on essentially every axis measured (§2) and should
   not be considered further without new evidence.
