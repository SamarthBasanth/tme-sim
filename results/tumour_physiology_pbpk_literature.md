# Tumour physiology numbers for PBPK modelling — literature survey

**Status: reference material, not yet wired into any model.** This is a literature
survey (interstitial volume fraction, interstitial fluid pressure, tumour water
content, and imaging-derived extracellular volume fraction), gathered to support
future PBPK/IFP calibration work — most directly the Darcy interstitial-pressure
solver (`tmesim/physics/pressure.py`) and the oxygen-model calibration gap
identified in [`results/oxygen_calibration_landscape.md`](oxygen_calibration_landscape.md).
No code in this repository consumes these numbers yet. Every claim below carries a
PMID/DOI; where a primary source was paywalled or unreachable, that is stated
explicitly rather than guessed at.

Two citations were independently re-verified against PubMed before this file was
committed (Kim 2004, PMID 15334566, and Dabir 2014, PMID 25384607 — both match
their reported numbers exactly). The remainder carry the researching agent's own
stated confidence levels and are not independently re-checked here; treat them
with that in mind before using any number as a calibration target.

---

## A. Measured tumour interstitial / extracellular volume fraction

**Classic rat tumour interstitial space (Gullino/Jain lineage), % of tumour volume:**

| Tumour | Interstitial space |
|---|---|
| Walker 256 carcinoma | 36.3% |
| Hepatoma H5123 | 43.3% |
| Hepatoma H3683 | 50.6% |
| Fibrosarcoma 4956 | 52.6% |
| Novikoff hepatoma | 54.6% |
| Fibrosarcoma A-MC | 60% |
| Fibrosarcoma C-MC | 55% |

Source: Wilhelm S, Vishwas D et al., "Non-invasive imaging of interstitial fluid
transport parameters in solid tumors in vivo," *Sci Rep* 2023;13:6721, PMID
37130836, PMCID PMC10154396, [DOI](https://doi.org/10.1038/s41598-023-33651-9) —
full text read, citing this table to the earlier rat-tumour literature.

Primary Gullino papers behind these numbers (PubMed-confirmed, abstract-only —
no abstract text survives on PubMed for these pre-1975 records, titles/authors
verified, table values not independently retrieved — AACR returned 403 on every
attempt):
- Gullino PM, Grantham FH, Smith SH. "The Interstitial Water Space of Tumors."
  *Cancer Res* 1965;25:727-31. PMID 14347560.
- Gullino PM, Clark SH, Grantham FH. "The Interstitial Fluid of Solid Tumors."
  *Cancer Res* 1964;24:780-94. PMID 14190544.
- Gullino PM, Grantham FH. "The Vascular Space of Growing Tumors." *Cancer Res*
  1964;24:1727-32. PMID 14230919.

**Kim et al. 2004 (dual-contrast MRI, absolute quantification) — independently
re-verified:** tumour (MDA PCa 2b s.c. xenograft) interstitial volume fraction
**27.5 ± 9.1%** vs skeletal muscle **15.9 ± 0.7%** (P<0.05); absolute blood
volumes did not differ significantly between the two tissues. Kim YR, Savellano
MD, Savellano DH, Weissleder R, Bogdanov A. "Measurement of tumor interstitial
volume fraction: method and implication for drug delivery." *Magn Reson Med*
2004;52:485-94. PMID 15334566, [DOI](https://doi.org/10.1002/mrm.20182).

**Krol et al. 1999 — "available volume fraction" (K_AV) in rat fibrosarcoma, ex
vivo tracer-exclusion method:** the interstitial fluid space, quantified by K_AV
of inulin, was **50% of total tissue volume**. K_AV of dextran fell sharply from
**0.28 ± 0.14 (Mr 40,000) to 0.10 ± 0.06 (Mr 70,000)** — i.e. available volume
fraction is strongly size-dependent above ~40–70 kDa, and was heterogeneous
between tumours more than within them. Krol A, Maresca J, Dewhirst MW, Yuan F.
"Available volume fraction of macromolecules in the extravascular space of a
fibrosarcoma: implications for drug delivery." *Cancer Res* 1999;59:4136-41.
PMID 10463619.

**Thurber/Wittrup antibody-transport models** use a lumped "void fraction" ε =
interstitial space / total tumour volume as a model parameter (a parameterisation
convention, not a new measurement). Wittrup KD, Thurber GM, Schmidt MM, Rhoden
JJ. "Practical theoretic guidance for the design of tumor-targeting agents."
*Methods Enzymol* 2012;503:255-68. PMID 22230572, PMC3978464,
[DOI](https://doi.org/10.1016/B978-0-12-396962-0.00010-0) — full text confirms
the definition; the specific numeric parameter-table value used in their work was
not retrievable from the accessible text.

**El-Kareh, Braunstein, Secomb 1993 (theoretical):** models interstitial volume
fractions "for a range... from 10 to 65%" as physiologically relevant, and notes
cellular exclusion effects are large "for small interstitial volume fractions
(around 20%)." PMID 8324199, PMC1262492,
[DOI](https://doi.org/10.1016/S0006-3495(93)81532-7).

**Model-derived value (secondary, not a direct measurement):** Piranfar A et al.,
3-D image-based prostate tumour transport model, cite fractional interstitial
volume (FVi) = **39%**, fractional cellular volume (FVc) = 61%, from a cited
reference the researching agent could not itself trace to a primary source.
*NPJ Syst Biol Appl* 2024;10:39. PMID 38609421, PMC11015041,
[DOI](https://doi.org/10.1038/s41540-024-00362-4).

**Confidence assessment (as reported by the researching agent).** *Replicated /
well-established*: the core finding that tumour interstitial volume fraction is
markedly larger than most normal tissues (rat tumours ~36–60% vs muscle 15.9%
per Kim 2004) is consistent across three independent methods (Gullino rat-tumour
series, Krol's ex vivo tracer exclusion, Kim's in vivo MRI) — but the exact
per-tumour-type numbers rest on a small set of 1960s rat models plus one 2004
mouse prostate xenograft, not a broad modern multi-lab consensus. Working central
estimate for the bulk interstitial fluid space: **≈0.3–0.5**, dropping sharply as
an "available" fraction for molecules above ~40–70 kDa due to steric exclusion.

---

## B. Tumour interstitial fluid pressure (IFP)

**Human tumours — primary measurement studies (mmHg):**

| Tumour type | Value | n | Method | Source |
|---|---|---|---|---|
| Cervix carcinoma (most replicated) | median 19, range −3 to 48 | 77–107 (4 cohorts) | wick-in-needle (WIN) | Roh 1991 PMID 1742744; Milosevic 1998 PMID 9635535; Milosevic 2001 PMID 11522633; Fyles 2006 PMID 16920212 |
| Melanoma, s.c. mets | 14.3 ± 12.5 (small 5.8±2, large 22.8±13) | 12 lesions | WIN | Boucher 1991, PMID 1742743 |
| Melanoma nodules | 29.8 (baseline) | 22 | WIN | Curti 1993, PMID 8485703 |
| Non-Hodgkin lymphoma nodules | 4.7 | 7 | WIN | Curti 1993, PMID 8485703 |
| Head & neck SCC | 13.2 ± 8.8 (range 4–33) | 19 | WIN | Gutmann 1992, PMID 1551128 |
| Breast carcinoma | 15 ± 9 vs 29 ± 3 SE — **two studies disagree** | subsets | WIN / intra-op needle | Less 1992 PMID 1423283; Nathanson 1994 PMID 7850532 |
| Colorectal liver mets | 21 ± 12 | subset | WIN | Less 1992, PMID 1423283 |
| Renal cell carcinoma | 38 (single tumour) | 1 | WIN | Less 1992, PMID 1423283 |
| All neoplasms pooled | 20 ± 13, range 4–50 | 41/21 tumours | WIN | Less 1992, PMID 1423283 |
| Intracranial tumours | 2.0 ± 2.5 — **confounded**: measured after dexamethasone+mannitol+furosemide+hypocapnia | 11 | WIN | Boucher 1997, PMID 9062403, PMC2063404 |
| Lung cancer | 8.5 | 211 | intra-op needle sensor | Mori 2015, PMID 25641031 |
| Thymic epithelial tumours | 11.3 ± 6.0 (1.5–28.6), rising stage I→IV (8.9→19.2) | 44 | Millar piezoelectric catheter | Mori 2017, PMID 29050375 |

**Review-level synthesis:** Nia HT, Munn LL, Jain RK, "Physical traits of
cancer," *Science* 2020;370:eaaz0868, PMID 33122355, PMC8274378,
[DOI](https://doi.org/10.1126/science.aaz0868) — human/animal tumour IFP
"ranging from <1 kPa (7.5 mmHg) in brain tumors to 5 kPa (37 mmHg) in renal cell
carcinomas," fairly uniform within a tumour, dropping precipitously at the
margin where IFP≈0.

**Animal/xenograft tumours (mmHg):**

| Model | Value | Method | Source |
|---|---|---|---|
| Rat R3230AC mammary, tissue-isolated | centre 18.4±9.3; microvascular pressure 17.3±6.1 (IFP≈MVP) | servo-null micropuncture + WIN | Boucher & Jain 1992, PMID 1516068 |
| Human osteosarcoma xenograft | orthotopic 30±9 vs s.c. 17±6 (p<0.001) | WIN | Brekken 2000, PMID 11131654 |
| Autochthonous murine PDAC | 75–130, mean 99 (n=7); normal mouse pancreas 8–13, mean 10.4 (n=8); falls to 20–30 post-PEGPH20 | Millar piezoelectric catheter (captures gel-bound fluid too, unlike WIN) | Provenzano 2012, PMID 22439937, PMC3371414 |

**Normal tissue IFP (mmHg):** consistently near-zero/slightly negative: human
tissue adjacent to tumours 2±4 (Less 1992); normal breast −0.3±0.1 SE (Nathanson
1994); normal cervix 0–3 (Roh 1991); mouse skin/subcutis −0.9±0.8 (Leunig 1992,
PMID 1384965); Nia 2020 review states normal-organ IFP "close to 0 mmHg."

**Confidence assessment.** *Replicated, robust*: cervix carcinoma (4 independent
cohorts, n>100 total, consistent median ~19 mmHg) is the strongest human IFP
dataset available; the general elevation of tumour IFP over near-zero adjacent
normal tissue is consistent across many tumour types and methods. *Single-study
/ weak priors*: renal cell carcinoma (n=1), lymphoma (n=7, one study), breast
(two studies disagree ~2-fold), colorectal liver mets, head & neck, brain
(confounded by peri-operative anti-oedema drugs), lung, thymic, and the mouse
PDAC 99 mmHg figure (no verified human PDAC IFP measurement exists).

**Methodological caveat worth carrying forward**: wick-in-needle (WIN, used for
essentially all classic human numbers above) measures only the free-fluid
pressure phase and may systematically underestimate IFP in stroma/hyaluronan-rich
tumours relative to piezoelectric-catheter methods — DuFort CC et al., *Biophys
J* 2016;110:2106-19, PMID 27166818, measured both simultaneously in the same
tissue and demonstrated the discrepancy directly. IFP is also spatially
near-uniform through the tumour core and drops to normal only in the outer
0.2–1.1 mm — Boucher, Baxter, Jain, *Cancer Res* 1990;50:4478-84, PMID 2369726.

**Not retrievable**: Heldin CH, Rubin K, Pietras K, Östman A, "High interstitial
fluid pressure — an obstacle in cancer therapy," *Nat Rev Cancer* 2004;4:806-13,
PMID 15510161 — closed-access everywhere checked (Nature, ResearchGate,
Academia.edu, Semantic Scholar API). Not fabricated from secondary citations.

---

## C. Tumour water content (% wet weight)

| Tissue | Water content | n | Method | Source |
|---|---|---|---|---|
| Panc1 pancreatic tumour (mouse) | 83.88% | 2 | wet vs 72h-lyophilized dry weight | Penet MF et al., *Front Oncol* 2021;10:599204, PMID 33585215, PMC7873637 |
| BxPC3 pancreatic tumour (mouse) | 84.00% | 5 | same | same |
| Normal mouse pancreas | 70.53% | 4 | same | same |
| Mouse liver / lung / kidney / heart / muscle / spleen (normal) | ~68% / ~80% / ~82% / ~79% / ~76% / ~79% | 4 each | same | same |

This paper also reports (SHG microscopy, qualitative) substantially higher
collagen fibre density in both pancreatic tumour types than normal pancreas, as
a proposed mechanistic explanation for elevated water content.

**Confidence assessment.** The Penet 2021 pancreatic numbers are a single, direct,
well-described study (n=2–5 per line) — solid as a data point, not independently
replicated for that tumour type in this search. The broader qualitative finding
(tumours run wetter than most normal tissues) recurs across pancreatic and
melanoma mouse studies, but a large, cross-tumour-type, well-replicated dataset
akin to the IFP literature was not found. The classic Gullino papers almost
certainly contain water-content numbers directly, but their table values were
not retrievable (see §A).

---

## D. Extracellular/extravascular volume fraction by imaging (human tumours)

**DCE-MRI ve (extended Tofts model), by tumour type:**

| Tumour type | ve | vp | n | Source |
|---|---|---|---|---|
| Cervical SCC, pre-treatment | 0.64 ± 0.22 (non-residual) / 0.64 ± 0.20 (residual) | 0.14±0.08 / 0.13±0.09 | 48 | Liu B et al., *Front Oncol* 2020;10:585738, PMID 33194734 |
| Head & neck SCC | 0.548 ± 0.119 | not reported | 16 | Surov A et al., *Transl Oncol* 2017;10:17-21, PMID 27888709 |
| Soft-tissue tumours, malignant vs benign | Tofts: 0.365±0.142 vs 0.315±0.146; ExTofts: 0.345 (0.234–0.459) vs 0.248 (0.202–0.396) — no significant ve difference by any of 5 models (all P>0.05) | vp significantly higher in malignant, all models P<0.001 | 92 (61/31) | Gao A et al., *Cancer Imaging* 2024;24:64, PMID 38773660 |

**Reproducibility caveat, directly demonstrated**: in a glioblastoma test-retest
study, ve reproducibility coefficient was 81–100% and coefficient of variation
40–48% depending on model and arterial-input-function method — a repeat scan of
the same tumour a day later can give a ve estimate differing by close to half
its own value. Bane O et al., *Neurooncol Adv* 2021, PMID 34988454, PMC8715899.

**CT/MRI extracellular volume fraction (ECV), by tumour type:**

| Tumour type | ECV | n | Source |
|---|---|---|---|
| PDAC, resectable, post-op | median 32% (IQR 21–41%); strongly correlated with histological collagen ratio (Pearson R=0.860, p<0.001) | 166 | Ming S et al., *Eur J Radiol* 2025;191:112273, PMID 40618544 |

**Normal tissue ECV, for comparison — Dabir et al. 2014 independently
re-verified**: healthy myocardium **0.25 ± 0.04 (1.5T) / 0.26 ± 0.04 (3T)**,
multicentre reference-value study, n=102 healthy subjects. Dabir D et al., "Reference
values for healthy human myocardium using a T1 mapping methodology," *J
Cardiovasc Magn Reson* 2014;16:69. PMID 25384607, PMC4203908,
[DOI](https://doi.org/10.1186/s12968-014-0069-x). Also: liver/spleen/muscle
(equilibrium contrast-enhanced MRI, n=40 healthy volunteers) — spleen 0.34,
liver 0.29, muscle 0.09, significantly different between organs (χ²=31.0,
P<0.001). Bandula S et al., *Radiology* 2013;268:858-64, PMID 23674785.

**Confidence assessment.** Normal-tissue ECV reference values (myocardium ~25–26%
especially) are well-established multicentre consensus figures. Tumour ve values
are markedly less well-constrained: the same benign lesions produced ve ranging
0.248–0.476 across five different kinetic models in the Gao 2024 study, and
test-retest reproducibility is poor (~40–48% CoV in glioblastoma). Tumour ECV
(PDAC, CT/MRI equilibrium method) has a clear biological correlate (collagen
content) but rests on single-region cohort studies, not an independently
replicated multicentre standard.

---

## Summary of confidence levels

| Parameter | Confidence |
|---|---|
| Tumour interstitial volume fraction ~30–60%, well above normal tissue (~10–20%) | Well-established (3 independent methods converge); exact per-tumour numbers rest on old/sparse studies |
| Tumour IFP elevated (10–30 mmHg typical) vs normal (~0 mmHg) | Well-established; cervix data specifically is robust (n>100, 4 studies) |
| Tumour IFP for any specific tumour type other than cervix | Weak — mostly single studies, some contradictory (breast) |
| Tumour water content > many normal tissues | Directionally consistent but sparse (pancreatic + melanoma mouse data only) |
| Normal-tissue ECV (myocardium ~25%, liver ~29%, muscle ~9%) | Well-established, multicentre reference values |
| Tumour ve by DCE-MRI, absolute values | Poorly constrained — model-dependent, ~40–50% test-retest CoV; use with caution as a calibration target |
| Tumour ECV by CT/MRI (PDAC) | Moderately established, single-region literature, good biological validation (collagen correlation) |

## What this could feed, if wired up later

- **Interstitial volume fraction (§A)** is the natural calibration target for
  `tmesim/physics/pressure.py`'s porosity term, if that solver is ever
  parameterised against a real physiological range rather than a default.
- **IFP magnitudes (§B)**, especially the cervix dataset's robustness, could
  anchor a sanity-check range for whatever the Darcy solver currently predicts,
  the same way `results/oxygen_calibration_landscape.md` anchored real pO2
  ranges against the oxygen model.
- **Nothing here should be treated as GBM-specific** — the intracranial-tumour
  IFP row above is both a single study and confounded by anti-oedema medication,
  and no glioma-specific interstitial volume fraction was found. Any future use
  of this table for the GBM flagship should say so explicitly rather than borrow
  a number silently from a different tumour type.
