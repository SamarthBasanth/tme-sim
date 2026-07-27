# The Real Build — Validation Summary

**What this is:** a patient-specific glioblastoma **microenvironment** simulator built
entirely from **real measured data**. Every synthetic element has been removed from the
pipeline (the synthetic cohort generator and all fabricated-geometry loaders are deleted;
the synthetic-validation era is archived under `archive/synthetic_era/` with its findings
preserved).

---

## 1. The real data foundation (all open, all downloaded)

| Layer | Real dataset | Scale |
|---|---|---|
| Cell reference | Darmanis 2017, GEO **GSE84465** | 23,460 genes x 3,589 annotated cells, 7 cell types |
| Single-cell spatial | **10x Xenium** GBM section (CC BY 4.0) | **40,887 cells**, true um coordinates, 319-gene panel |
| Spot spatial | 4 Visium GBM sections (10x public + **GSE194329**) | 1,775-4,764 spots each, whole-transcriptome |
| Niche architecture | **Ivy GAP** (Allen Institute) | 270 samples across 7 anatomic structures |
| Drug potency | **GDSC** glioma + **HGCC** (Johansson 2020, PMID 32668248) | 262 drugs x 85 patient-derived cultures |
| Physics ground truth | 7 PMID-verified tissue-PK anchors, 6 drugs | penetration + brain:plasma ratios |

## 2. Reconstruction: real sample -> real whole microenvironment

The Xenium section reconstructs into a real, cell-resolved niche — not just tumor:
**TUMOR 73.8% / STROMAL 7.9% / IMMUNE 8.5% / ENDOTHELIAL 9.7%** on a 8130 x 5100 um grid.
Vasculature is a KDE of the *real endothelial cell positions*; ECM is the expression-weighted
KDE of **11 real matrix transcripts** (COL1A2, POSTN, VCAN, DCN, LAMA2, THBS1, ...); oxygen
is solved by the physics engine from that real perfusion field (93.7% of tumor cells hypoxic).

### Fidelity — validated, with statistics
- **Cross-platform concordance:** Xenium (marker phenotyping) vs Visium (NNLS deconvolution),
  two independent platforms and methods: **Pearson 0.997**; endothelial agrees to **0.03%**
  (9.74% vs 9.77%). Mean absolute difference 4.44%.
- **Niche labels are marker-faithful** (Ivy GAP signatures, 10,000-permutation tests):
  PECAM1 in perivascular **d=1.29, p=1e-4**; VEGFA in necrotic **d=1.35, p=1e-4**;
  CA9 in necrotic **d=0.67, p=1e-4**.
- **Deconvolution is stable:** +-1.0% per-spot std under reference resampling; held-out
  pseudobulk recovery correlation **0.935** (RMSE 8.3%).

A methodological bug was found and fixed honestly: the initial marker-argmax classifier
inflated ENDOTHELIAL to 20.6% (single-marker PECAM1 bias). A threshold-gated classifier
corrects it to 9.7% — independently corroborated by the Visium deconvolution's 9.77%.

## 3. Physics: validated against real measured pharmacokinetics

One mechanistic constant (alpha) is fit **once** to doxorubicin's real perivascular anchor
(Primeau 2005, PMID 16361566), then held fixed:

| drug | role | predicted L | real cited band | verdict |
|---|---|---|---|---|
| doxorubicin | in-sample fit | 64.9 um | 58-72 um | fit by construction |
| mitoxantrone | **out-of-sample** | 48.2 um | 36-108 um | **PASS** |
| topotecan (distinct class) | **out-of-sample** | 83.6 um | 36-108 um | **PASS** |

BBB-penetrance ordering vs real brain:plasma ratios is only *partially* concordant
(Spearman 0.50): methotrexate correctly least-penetrant, but cisplatin/temozolomide are
inverted — a known limit of structure-only rules (no efflux/carrier transport).

## 4. Pharmacology: real-vs-real corroboration
- Structure -> potency prediction is **weak and confounded** (apparent signal is bulky
  natural-product cytotoxics vs inert small molecules, not true QSAR). This is *why* the
  layer uses measured IC50 rather than predicted.
- **GDSC vs HGCC agree at Spearman +0.83 (p=2.4e-5)** across overlapping drugs — two
  independent real screens corroborate the measured-potency axis, no model in the loop.

## 5. Simulation on real tissue

**Full drug panel, one real tumor** — delivery, not potency, decides the outcome:
doxorubicin is **50x more potent** than cisplatin (IC50 0.5 vs 24.6 uM) yet is out-performed,
because it is delivery-starved (BBB 0.051, largest MW, trapped ~45 um from vessels, lowest
tissue AUC by two orders of magnitude). Best penetration: temozolomide (0.396). Best
predicted kill: cisplatin (6.1%).

**Five real tumors, same drug** — real inter-tumor heterogeneity drives the physics:

| tissue | vascular density | predicted penetration | predicted kill |
|---|---|---|---|
| Xenium-GBM | 0.094 | 0.396 | 1.9% |
| GBM1 | 0.110 | 0.264 | 9.4% |
| 10x_parent | 0.132 | 0.283 | 11.6% |
| GBM2 | 0.302 | 0.509 | 24.1% |
| GBM3 | 0.325 | 0.566 | 23.5% |

Real vascularity spans **3.5x** across actual patients and tracks predicted penetration
almost monotonically; predicted kill spans **1.9% -> 24%** from microenvironment differences
alone. (Hypoxia is saturated at 94-100% across this cohort, so it does *not* discriminate
these tumors — vascularity does.)

## 6. What is validated, and what is NOT

**Validated on real data:** the reconstruction is faithful and cross-platform concordant;
the niche labels are marker-faithful; the deconvolution is stable; the transport physics
reproduces real measured penetration for 2 out-of-sample drugs; the potency axis is
corroborated by an independent patient-derived screen.

**NOT validated:** that predicted drug *response* matches real patient outcomes. Every open
spatial GBM dataset is an untreated resection with no paired measured response, and every
open drug-response dataset (GDSC, HGCC) is dissociated cells with no microenvironment. No
open dataset bridges the two.

**The remaining gap is data access, not engineering.** Closing it requires spatial tissue
paired with measured per-model drug kill: HTAN GBM spatial multi-omics (dbGaP + HTAN DUA)
with NCI PDMR / organoid-biobank drug screens (MTA). The pipeline is built to ingest exactly
that the day it is available.
