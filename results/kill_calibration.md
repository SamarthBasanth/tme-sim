# ABSOLUTE KILL: the 3-D:2-D ratio test

**The measured anchor we hold is a 2-D monolayer assay. The model is a 3-D microenvironment. That gap is not error — it is the phenomenon the model exists to reproduce, and this file tests whether it reproduces it.**

> ### Headline
>
> **The 3-D:2-D ratio test PASSES for all 3 drugs.** On the real Xenium section this model predicts a temozolomide 3-D:2-D kill ratio of **0.2626** at the matched ex-vivo bath concentration (100 uM) — equivalent to an IC50 shift of **3.26x**. The measured spheroid-vs-monolayer range for temozolomide in glioma is **2-5x** (median 3x, n = 6 line/oxygen conditions, PMID 30694423). The prediction lands essentially **on the literature median**, from a model that has never been fitted to any 3-D-vs-2-D measurement. It could have missed by an order of magnitude and did not.
>
> **This does NOT make absolute kill defensible yet, and the reason has moved.** The test constrains the *product* of delivery and cell-state resistance; it does not constrain either factor separately, and it says nothing about intrinsic potency. Anchoring intrinsic potency to the measured ex-vivo viability gives temozolomide a 3-D kill of **0.0326** on this section — but with an error bar of **0.0048-0.0831** (17.4x) propagated from the inter-culture spread of the anchor itself.
>
> **The surprise is which term blocks it — with one important qualification.** The cell-state label is *bounded by construction*: the multipliers are x0.35 and x0.5, so no revision of *which* cells are labelled quiescent can move absolute kill by more than **2.9x**, and over the full 0-100 % quiescent sweep it moves 2.00x. The measured anchor spread moves it 17.4x. **Absolute kill is blocked more by the potency anchor's inter-patient spread than by the cell-state label** — which inverts the expectation in `results/dose_schedule.md`. What the cell-state work must deliver is not primarily the *fraction* but a defence of the *multipliers themselves* (x0.35 / x0.5), which are unvalidated priors in `PDConfig` and are NOT bounded by anything above.
>
> **Independently corroborated.** The concurrent cell-state workstream (`results/cell_state_impact.md`, landed during this run) measures the same term three different ways on the same section and gets `F_state` = 0.5502-0.9550 (1.00x-1.74x the default) — inside the bound above. Its oxygen-route kills also reproduce this module's engine-path values bit-for-bit for all three drugs. Section 5.1 reconciles the two in full.
>
> **The qualification, and it matters.** Those two bars are not the same kind of quantity. The 17.4x is largely *biological* spread between patients — it does not shrink with a better model, it shrinks by measuring this patient's own sensitivity. The cell-state bar is *epistemic*: it is our ignorance, and it is the one a better model can remove. Asked instead how well the *cohort-typical* GBM's intrinsic potency is pinned down (bootstrap CI on the median rather than the inter-patient IQR) the potency bar narrows to 0.0205-0.0569 (2.8x) — comparable to the cell-state bar. So the correct statement is narrower than 'cell state no longer matters': **for a patient whose sensitivity is unmeasured the anchor dominates; for the cohort-typical patient the two terms are the same size and BOTH must be closed.**
>
> ### The bottom line
>
> **Absolute kill is defensible as an INTERVAL, not as a number, and the interval is wide.** Propagating both open terms — the measured potency spread and the measured cell-state route span — gives temozolomide **0.0048-0.1442** (30x), cisplatin 0.2667-0.8852, doxorubicin 0.1666-0.8065 (section 5.1). Every term in those intervals traces to a measurement or to a named unvalidated prior, which is new; a 30x interval is still not a headline result, which is not new. **The cross-drug ranking — invariant across all three cell-state routes and across the whole potency band — remains the honest primary output of this simulator, exactly as before.** What has changed is that the 3-D apparatus itself has now passed a real external test, and that the absolute number, when quoted, comes with an auditable bar instead of a disclaimer.

## 0. What was run

* **Tissue** — the real 10x Xenium GBM section, reconstructed ONCE at the grid-converged `dx = 12 um` (40,887 cells, 23,812 tumour), reused for every drug. `dt = 10 s`, window 6 h, convection on, seed 0.
* **Drugs** — `calibrated_drugparams(name, alpha = 148.61)` then `apply_barrier_to_drugparams`: the identical compilation pipeline the drug panel uses.
* **One transport solve per drug.** The PDE is linear in the source and starts from `C = 0`, so `AUC(P) = P x AUC(1)` exactly; every concentration below rescales that one field and re-runs the engine's own unmodified kill mapping (`cells.pkpd.apply_drug_kill`). Nothing in core / physics / engine / molecular / reconstruct is modified.
* **Cell-state bands on this section** (from the oxygen solve): VIABLE 10.0%, HYPOXIC 0.0% (x0.35), QUIESCENT 90.0% (x0.5).

## 1. The bridge parameter, from the literature

This is the **bridge parameter** and it comes from published measurements only — nothing here is fitted to our model, and nothing here is fitted to the ex-vivo anchor. `fold` is the IC50 ratio 3-D / 2-D for the *same cells* at the *same nominal exposure*; `fold > 1` means the 3-D format is more resistant. Rows marked 'full text Table 1' were read from the paper's own results table in the PubMed Central full text; rows marked '(figure)' were read from a plotted figure and are correspondingly softer. Every PMID was verified against PubMed metadata.

| PMID | study | system | drug | 2-D | 3-D | **fold 3D/2D** | read from | note |
|---|---|---|---|---|---|---:|---|---|
| 30694423 | Musah-Eroje A & Watson S. J Neurooncol 2019;142(2):231-240 | U87, Cultrex BME 3D-TGA vs 2D monolayer, normoxia | temozolomide | IC50 (2D) | IC50 (3D) | **3x** | Fig. 4c (figure) | 3-D TGA vs monolayer, same cells, 72 h AlamarBlue, TMZ 5-1500 uM |
| 30694423 | Musah-Eroje A & Watson S. J Neurooncol 2019;142(2):231-240 | U87, 3D vs 2D, hypoxia (1% O2) | temozolomide | IC50 (2D) | IC50 (3D) | **4x** | Fig. 4c (figure) | hypoxia had NO effect on TMZ IC50 in 2-D; it only acts in combination with 3-D |
| 30694423 | Musah-Eroje A & Watson S. J Neurooncol 2019;142(2):231-240 | U251, 3D vs 2D, normoxia | temozolomide | IC50 (2D) | IC50 (3D) | **2x** | Fig. 4c (figure) | smallest shift in the set |
| 30694423 | Musah-Eroje A & Watson S. J Neurooncol 2019;142(2):231-240 | U251, 3D vs 2D, hypoxia (1% O2) | temozolomide | IC50 (2D) | IC50 (3D) | **5x** | full text (stated: 'as high as fivefold in the U251 cells') | the LARGEST shift the paper reports, and the only per-line value stated in text |
| 30694423 | Musah-Eroje A & Watson S. J Neurooncol 2019;142(2):231-240 | SNB19, 3D vs 2D, normoxia | temozolomide | IC50 (2D) | IC50 (3D) | **3x** | Fig. 4c (figure) |  |
| 30694423 | Musah-Eroje A & Watson S. J Neurooncol 2019;142(2):231-240 | SNB19, 3D vs 2D, hypoxia (1% O2) | temozolomide | IC50 (2D) | IC50 not reached | **2x** | Fig. 4c (figure) | LOWER BOUND -- IC50 not achieved at the top dose |
| 27576873 | Gomez-Roman N, et al. Neuro Oncol 2017;19(2):229-241 | G7/E2/R10 patient-derived GBM stem-like cells, Alvetex polystyrene scaffold, clonogenic survival | temozolomide | significant in 1/2 lines | significant in 2/2 lines | — | abstract + full text | **COUNTER-EXAMPLE**: TMZ was MORE effective in 3-D than 2-D. Excluded from the band and reported separately -- the sign is not settled |
| 32365781 | Samiei E, et al. Int J Mol Sci 2020;21(9):3162 | U251/U87 3-D hydrogel GBM-on-chip vs 2-D monolayer | temozolomide | higher apoptosis | lower apoptosis | — | abstract + full text | 'cells were significantly less sensitive to drugs in the 3D model as compared to monolayer'; the gap is dose-dependent and closes at the top dose |
| 12569392 | Gunther W, et al. Br J Cancer 2003;88(3):463-469 | U87MG / GaMG multicellular spheroids | temozolomide | — | apoptosis highest in the CORE | — | abstract | argues AGAINST penetration limitation for TMZ specifically -- it is small, uncharged and 100% orally bioavailable |
| 34794398 | Poon MTC, et al. BMC Cancer 2021;21(1):1240 | systematic review, 212 studies, 248 glioma lines, **2-D only** | temozolomide | U87 72 h IC50 median 230 uM, IQR 34.1-650.0 uM | — | — | abstract | **the 2-D reference itself has a ~19x IQR** -- larger than the 2-5x bridge it is being used to correct |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | CaCo2 colorectal, Matrigel spheroid vs monolayer, 72 h RT-Glo | cisplatin | 25.0 uM | 30.0 uM | **1.2x** | full text Table 1 |  |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | Colo205 colorectal, Matrigel spheroid vs monolayer | cisplatin | 75.0 uM | 30.0 uM | **0.4x** | full text Table 1 | 3-D MORE sensitive |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | HCT116 colorectal, Matrigel spheroid vs monolayer | cisplatin | 7.5 uM | 40.0 uM | **5.33x** | full text Table 1 |  |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | SW480 colorectal, Matrigel spheroid vs monolayer | cisplatin | 5.0 uM | 30.0 uM | **6x** | full text Table 1 | averaged over the 4 lines the shift is NOT significant (p = 0.818) |
| 32152850 | Trigub MM, et al. Bull Exp Biol Med 2020;168(4):583-588 | HeLa spheroids vs monolayer | cisplatin | IC50 (monolayer) | IC50 (spheroid) | **5x** | abstract | 'resistance index' 5.0 +/- 1.5; Pt(IV) complexes in the same study 1.8-2.3x |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | CaCo2 colorectal, Matrigel spheroid vs monolayer | doxorubicin | 0.40 uM | 0.40 uM | **1x** | full text Table 1 |  |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | Colo205 colorectal | doxorubicin | 0.20 uM | 0.25 uM | **1.25x** | full text Table 1 |  |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | HCT116 colorectal | doxorubicin | 0.08 uM | 0.30 uM | **3.75x** | full text Table 1 |  |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | SW480 colorectal | doxorubicin | 0.10 uM | 0.30 uM | **3x** | full text Table 1 | averaged over the 4 lines NOT significant (p = 0.121) |
| 37869710 | Jubelin C, et al. Front Bioeng Biotechnol 2023;11:1260049 | MNNG/HOS osteosarcoma, liquid-overlay spheroids vs 2-D | doxorubicin | 0.8 +/- 0.4 uM | 15.07 +/- 0.3 uM | **18.8x** | abstract (values stated) | the largest cytotoxic shift found; doxorubicin is THE drug Durand 1986 identified as penetration-limited |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | CaCo2 / Colo205 / HCT116 / SW480 colorectal | 5-fluorouracil | 10 / 10 / 7.5 / 75 uM | 75 / 60 / 75 / 90 uM | **7.5x** | full text Table 1 | per-line folds 7.5 / 6.0 / 10.0 / 1.2; the only drug in this panel with a significant averaged shift (p = 0.027) |
| 33395433 | Koch J, et al. PLoS One 2021;16(1):e0244513 | CaCo2 / Colo205 / HCT116 / SW480 colorectal | mitomycin C | 7.5 / 1.0 / 0.5 / 1.0 uM | 5.0 / 1.0 / 0.5 / 0.5 uM | **0.67x** | full text Table 1 | per-line folds 0.67 / 1.0 / 1.0 / 0.5 -- 3-D EQUAL or MORE sensitive in all four lines (p = 0.297) |
| 41422245 | Mutahar AZI, et al. J Transl Med 2025;23(1):1428 | MDA-MB-231 / MDA-MB-468 / MCF-7 breast spheroids vs 2-D | paclitaxel | — | — | **25x** | abstract ('IC50 values elevated by roughly 25-fold') | upper tail of the observed distribution |
| 27512118 | Hoda MA, et al. Mol Cancer Ther 2016;15(10):2357-2369 | 6 mesothelioma lines + 13 primary cultures, adherent monolayer vs non-adherent spheroids | trabectedin | <= 2.6 nM | <= 2.6 nM | **1x** | abstract | NO shift at all |
| 25634491 | Imamura Y, et al. Oncol Rep 2015;33(4):1837-1843 | 6 breast lines, MCTS vs 2-D | paclitaxel, doxorubicin | — | — | — | abstract | **BIMODAL**: the 3 lines forming DENSE spheroids were resistant; the 3 forming LOOSE spheroids matched 2-D exactly. Compaction, not dimensionality, predicts the shift |
| 24194875 | Hongisto V, et al. PLoS One 2013;8(10):e77232 | JIMT1 breast, 102 compounds, Matrigel 3-D vs 2-D vs polyHEMA | 102-compound panel | — | — | — | abstract | **COUNTER-EXAMPLE at scale**: across the 63 most active compounds Matrigel-3-D was SIGNIFICANTLY MORE SENSITIVE than 2-D; polyHEMA 3-D matched 2-D. Authors: 'large variations in drug responses were observed between the models' |

* **temozolomide in glioma (the matched anchor)** (n = 6 comparisons): median fold **3x**, geometric mean 3x, IQR 2.2-3.8x, full range **2-5x**.
* **all cytotoxics, solid tumours generally** (n = 20 comparisons): median fold **3x**, geometric mean 2.9x, IQR 1.2-5.1x, full range **0.4-25x**.
* **cisplatin (all tumour types)** (n = 5 comparisons): median fold **5x**, geometric mean 2.4x, IQR 1.2-5.3x, full range **0.4-6x**.
* **doxorubicin (all tumour types)** (n = 5 comparisons): median fold **3x**, geometric mean 3.1x, IQR 1.2-3.8x, full range **1-19x**.

**Five things about this table that must not be glossed over.**

1. **Temozolomide-specific 3-D-vs-2-D data is thin.** Exactly ONE study (PMID 30694423) reports a per-line, matched 2-D-vs-3-D temozolomide IC50 fold table. Everything else is 3-D-only, cross-study, or qualitative. A band built from one paper is a weak band.
2. **The sign is not settled.** Gomez-Roman *et al.* (PMID 27576873) found temozolomide *more* effective in 3-D scaffold culture than in 2-D, and Hongisto *et al.* (PMID 24194875) found Matrigel-3-D systematically MORE sensitive than 2-D across the 63 most active of 102 compounds. In the best matched distribution source (PMID 33395433, 16 drug x line pairs read from Table 1) **5 of 16 pairs have fold <= 1**.
3. **The 2-D reference is noisier than the effect being corrected.** PMID 34794398 (systematic review, 212 studies) puts the between-study IQR of the U87 72 h temozolomide IC50 at **34-650 uM — a ~19x spread**. A 2-5x bridge parameter sits inside the measurement noise of the assay it bridges from. Passing this test is therefore a *weak* pass by construction: the target band is wide.
4. **Compaction, not dimensionality, is the predictive variable.** PMID 25634491 found the shift only in the 3 of 6 breast lines that formed DENSE spheroids; the loose ones matched 2-D exactly. PMID 16424039 isolates this causally via alpha-E-catenin.
5. **The gap closes at saturating dose** (PMID 32365781), so a 3-D:2-D kill ratio is only meaningful with its concentration attached — which is why section 2 reports a curve rather than a number.

*Literature retrieved from PubMed. Key primary sources: Musah-Eroje & Watson 2019 [DOI](https://doi.org/10.1007/s11060-019-03107-0), Koch *et al.* 2021 [DOI](https://doi.org/10.1371/journal.pone.0244513), Durand 1986 (no DOI; PMID 2425117), Kobayashi *et al.* 1993 [DOI](https://doi.org/10.1073/pnas.90.8.3294), Minchinton & Tannock 2006 [DOI](https://doi.org/10.1038/nrc1893).*

### 1.1 Which factor should dominate? The published mechanistic split

This is the qualitative prior against which section 2's decomposition has to be judged. The literature is unusually clear on the ordering, and it is **not** the ordering our model produces.

| PMID | study | factor | finding |
|---|---|---|---|
| 2425117 | Durand RE. J Natl Cancer Inst 1986;77(1):247-252 | penetration vs microenvironment | V79 spheroids sorted by DEPTH, plus pre-separated cells to separate inherent from locational sensitivity, over 7 drugs. Verbatim: 'penetration was a problem only for doxorubicin; in contrast, the microenvironment in the intact spheroid had a much greater influence on cell sensitivity' |
| 8475071 | Kobayashi H, et al. Proc Natl Acad Sci USA 1993;90(8):3294-3298 | contact / multicellular resistance | EMT-6 sublines made resistant in vivo showed NO resistance as monolayers, FULL resistance as spheroids, and LOST it again when the spheroids were trypsinised apart |
| 11033306 | Desoize B & Jardillier J. Crit Rev Oncol Hematol 2000;36(2-3):193-207 | taxonomy of multicellular resistance | splits MCR into contact resistance and resistance inherent in the spheroid structure; mechanisms listed as apoptosis inhibition, HIGH QUIESCENT FRACTION, altered topoisomerase/repair expression, permeability, hypoxic/necrotic centre |
| 16424039 | Grantab R, Sivananthan S, Tannock IF. Cancer Res 2006;66(2):1033-1039 | packing density -> penetration | same parental line, tight (epithelioid) vs loose (round) sublines differing in alpha-E-catenin: paclitaxel, doxorubicin, methotrexate and 5-FU all penetrated significantly better through the LOOSE multicellular layers, and doxorubicin killed better there |
| 16862189 | Minchinton AI & Tannock IF. Nat Rev Cancer 2006;6(8):583-592 | review | canonical review of incomplete drug distribution in tumour tissue |
| 30694423 | Musah-Eroje A & Watson S. J Neurooncol 2019;142(2):231-240 | what the TMZ shift is NOT | neurosphere (cancer-stem-cell-enriched) IC50 was SIMILAR to 2-D, so stemness alone does not produce the shift; hypoxia alone in 2-D had NO effect; and 3-D-grown cells re-plated as 2-D returned to baseline IC50 |

**The consensus ordering for a small, uncharged cytotoxic like temozolomide: quiescence / contact-mediated resistance FIRST, hypoxia a conditional multiplier SECOND, penetration LAST** (penetration is the leading term only for doxorubicin — PMID 2425117, PMID 16424039). Note that the mechanism the literature names first, a high quiescent fraction (PMID 11033306), is exactly the term our engine carries as the x0.5 QUIESCENT multiplier — so the model does contain the right mechanism; section 2 asks whether it is sized right relative to delivery.

## 2. Factor decomposition of the model's own kill

At a **matched nominal exposure** the plasma level is set so the best-perfused tumour cells (p99 of the tissue exposure distribution) see exactly the bath concentration — the analogue of a spheroid's outer surface sitting in the well. This deliberately puts the BBB and the engine's vascular-source convention OUTSIDE the 3-D:2-D comparison, because a spheroid in a dish has neither; that leg is reported separately below each table and is not part of the test.

```
kill_3D = mean_i [ p_hill(E_i) * sensitivity_i * resistance_i ]
kill_2D = p_hill(C_bath) * mean_i[ sensitivity_i ]      # delivery = 1, no cell-state resistance

F_penetration = mean_i[p_hill(E_i)] / p_hill(C_bath)    # (a) delivery / the transport solve
F_state       = mean_i[resistance_i]                    # (b) quiescent x0.5 / hypoxic x0.35
F_sensitivity = mean_i[sensitivity_i]                   # (c) intrinsic -- COMMON to 2-D and 3-D, cancels

predicted 3D:2D  =  F_penetration * F_state
```

### temozolomide — IC50 606 uM, hill 1.2, BBB 0.125

| bath concentration | plasma needed (uM) | median cell / bath | **F_penetration** (a) | **F_state** (b) | F_sens (c) | product a·b | exact ratio | interaction | implied IC50 fold |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1xIC50 (60.6 uM) | 26.77 | 0.478 | **0.4572** | **0.5502** | 1.000 | 0.2515 | **0.2568** | 1.021 | 3.23x |
| 0.3xIC50 (181.8 uM) | 80.31 | 0.478 | **0.4910** | **0.5502** | 1.000 | 0.2701 | **0.2754** | 1.020 | 3.34x |
| 1xIC50 (606 uM) | 267.7 | 0.478 | **0.5993** | **0.5502** | 1.000 | 0.3297 | **0.3348** | 1.015 | 3.81x |
| 3xIC50 (1818 uM) | 803.1 | 0.478 | **0.7693** | **0.5502** | 1.000 | 0.4233 | **0.4270** | 1.009 | 5.28x |
| 10xIC50 (6060 uM) | 2677 | 0.478 | **0.9179** | **0.5502** | 1.000 | 0.5050 | **0.5066** | 1.003 | 10.81x |
| 100 uM (ex-vivo assay) (100 uM) | 44.17 | 0.478 | **0.4678** | **0.5502** | 1.000 | 0.2574 | **0.2626** | 1.020 | 3.26x |

*Robustness of the p99 bath-reference convention:* using the **maximum** tumour-cell exposure instead of the p99 gives F_penetration 0.3340 (vs 0.4678) and an implied IC50 fold of 4.35x (vs 3.26x); using the p95 gives 2.53x. Both variants stay inside this drug's literature band, so the section 3 verdict does not turn on the convention.

*Plasma leg, reported but deliberately NOT part of the 3-D:2-D comparison:* the engine's source convention puts **2.264 uM** of time-averaged tissue exposure at the p99 reference per **1 uM** of window-mean plasma, so reaching the bath concentration above needs 44.17 uM of plasma. Note that ratio exceeds 1 for temozolomide and cisplatin: the engine's delivery term `S = vascular_density x bbb_penetrance x plasma` enters the PDE as a **rate**, not a concentration (the documented Stage-0 source convention in `treatment/drugs.py`), so it is NOT a physical tissue:plasma partition coefficient and must not be read as one. **This is precisely why the 3-D:2-D comparison is matched on tissue exposure rather than on plasma** — doing so removes this convention, and the BBB, from the test entirely.

### cisplatin — IC50 24.55 uM, hill 1.2, BBB 0.858

| bath concentration | plasma needed (uM) | median cell / bath | **F_penetration** (a) | **F_state** (b) | F_sens (c) | product a·b | exact ratio | interaction | implied IC50 fold |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1xIC50 (2.455 uM) | 0.372 | 0.380 | **0.3622** | **0.5502** | 1.000 | 0.1993 | **0.2041** | 1.024 | 3.92x |
| 0.3xIC50 (7.365 uM) | 1.116 | 0.380 | **0.3931** | **0.5502** | 1.000 | 0.2163 | **0.2212** | 1.023 | 4.05x |
| 1xIC50 (24.55 uM) | 3.72 | 0.380 | **0.4980** | **0.5502** | 1.000 | 0.2740 | **0.2791** | 1.019 | 4.55x |
| 3xIC50 (73.65 uM) | 11.16 | 0.380 | **0.6840** | **0.5502** | 1.000 | 0.3763 | **0.3807** | 1.012 | 6.07x |
| 10xIC50 (245.5 uM) | 37.2 | 0.380 | **0.8763** | **0.5502** | 1.000 | 0.4822 | **0.4843** | 1.005 | 11.60x |

*Robustness of the p99 bath-reference convention:* using the **maximum** tumour-cell exposure instead of the p99 gives F_penetration 0.3005 (vs 0.4980) and an implied IC50 fold of 7.29x (vs 4.55x); using the p95 gives 3.40x. **Not all variants stay inside this drug's literature band** — max-reference OUTSIDE, p95 INSIDE. The section 3 verdict for this drug is therefore convention-dependent and must be read as such.

*Plasma leg, reported but deliberately NOT part of the 3-D:2-D comparison:* the engine's source convention puts **6.6 uM** of time-averaged tissue exposure at the p99 reference per **1 uM** of window-mean plasma, so reaching the bath concentration above needs 3.72 uM of plasma. Note that ratio exceeds 1 for temozolomide and cisplatin: the engine's delivery term `S = vascular_density x bbb_penetrance x plasma` enters the PDE as a **rate**, not a concentration (the documented Stage-0 source convention in `treatment/drugs.py`), so it is NOT a physical tissue:plasma partition coefficient and must not be read as one. **This is precisely why the 3-D:2-D comparison is matched on tissue exposure rather than on plasma** — doing so removes this convention, and the BBB, from the test entirely.

### doxorubicin — IC50 0.5 uM, hill 1.2, BBB 0.021

| bath concentration | plasma needed (uM) | median cell / bath | **F_penetration** (a) | **F_state** (b) | F_sens (c) | product a·b | exact ratio | interaction | implied IC50 fold |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1xIC50 (0.05 uM) | 0.8875 | 0.263 | **0.2584** | **0.5502** | 1.000 | 0.1422 | **0.1462** | 1.029 | 5.19x |
| 0.3xIC50 (0.15 uM) | 2.662 | 0.263 | **0.2821** | **0.5502** | 1.000 | 0.1552 | **0.1595** | 1.027 | 5.37x |
| 1xIC50 (0.5 uM) | 8.875 | 0.263 | **0.3688** | **0.5502** | 1.000 | 0.2029 | **0.2077** | 1.023 | 6.03x |
| 3xIC50 (1.5 uM) | 26.62 | 0.263 | **0.5476** | **0.5502** | 1.000 | 0.3013 | **0.3063** | 1.016 | 7.78x |
| 10xIC50 (5 uM) | 88.75 | 0.263 | **0.7854** | **0.5502** | 1.000 | 0.4321 | **0.4355** | 1.008 | 13.56x |

*Robustness of the p99 bath-reference convention:* using the **maximum** tumour-cell exposure instead of the p99 gives F_penetration 0.1571 (vs 0.3688) and an implied IC50 fold of 12.91x (vs 6.03x); using the p95 gives 4.12x. Both variants stay inside this drug's literature band, so the section 3 verdict does not turn on the convention.

*Plasma leg, reported but deliberately NOT part of the 3-D:2-D comparison:* the engine's source convention puts **0.05634 uM** of time-averaged tissue exposure at the p99 reference per **1 uM** of window-mean plasma, so reaching the bath concentration above needs 8.875 uM of plasma. Note that ratio exceeds 1 for temozolomide and cisplatin: the engine's delivery term `S = vascular_density x bbb_penetrance x plasma` enters the PDE as a **rate**, not a concentration (the documented Stage-0 source convention in `treatment/drugs.py`), so it is NOT a physical tissue:plasma partition coefficient and must not be read as one. **This is precisely why the 3-D:2-D comparison is matched on tissue exposure rather than on plasma** — doing so removes this convention, and the BBB, from the test entirely.

**Reading the decomposition.**

* The two 3-D factors are **comparable in size**, not dominated by one. For temozolomide at the ex-vivo concentration, delivery contributes 56% of the log-penalty and cell state 44%. That matters, because the literature says the split should run the *other* way for temozolomide: Durand 1986 (PMID 2425117) found penetration limiting only for doxorubicin, and Gunther *et al.* (PMID 12569392) found temozolomide-induced apoptosis HIGHEST in the spheroid core. **Our model attributes more of the 3-D penalty to delivery than the spheroid literature does.** Three readings are possible and this analysis cannot separate them: (i) our gradient is too steep — plausible, since section 3 of `results/dose_schedule.md` shows the measured gradient on this section is largely inherited from the vascular-source field rather than produced by transport; (ii) our state penalty is too mild; or (iii) **both are right and the geometries genuinely differ** — a 200 um spheroid bathed at its surface has almost no gradient for a small molecule, whereas this section has vessel-to-vessel distances up to 786 um, so a larger delivery term here is what a correct model should produce. Reading (iii) is the charitable one and cannot be dismissed; a measured intratumoural drug gradient in human GBM tissue would decide it.
* The factorisation is nearly exact: the interaction term (exact ratio / product of marginals) is 1.003-1.029, i.e. exposure, sensitivity and cell state are close to independent across tumour cells on this section. The product is a fair summary.
* **The ratio rises steeply with dose.** For every drug F_penetration climbs towards 1 as the bath concentration saturates the Hill term, so the 3-D penalty shrinks — exactly the behaviour PMID 32365781 measured (the 2-D/3-D gap closed at the top temozolomide dose). Any single-number '3-D penalty' is meaningless without its concentration.
* **Intrinsic sensitivity cancels.** In `mode='uniform'` (the fair cross-drug mode) per-cell sensitivity is 1.0, so F_sensitivity = 1 and the whole intrinsic-potency question lives in the IC50 — which is where section 4 anchors it. Under `mode='determinant'` the marker-derived sensitivity would scale BOTH formats and would still largely cancel.

## 3. THE TEST — does the predicted 3D:2D ratio land inside the measured range?

| drug | reporting concentration | predicted 3D:2D kill ratio | implied IC50 fold | literature fold range | literature kill-ratio band | verdict |
|---|---|---:|---:|---|---|---|
| temozolomide | 100 uM (ex-vivo assay) | **0.2626** | 3.26x | 2-5x | 0.1590-0.4622 | **PASS** |
| cisplatin | 1xIC50 | **0.2791** | 4.55x | 0.4-6x | 0.2086-1.5004 | **PASS** |
| doxorubicin | 1xIC50 | **0.2077** | 6.03x | 1-19x | 0.0575-1.0000 | **PASS** |

**VERDICT: PASS, 3/3.** Every drug's predicted 3-D:2-D ratio lands inside its own drug-matched measured band.

**What this test is.** The model has never before been asked to predict a quantity that a published experiment measures independently of it. It predicts one here: the ratio of 3-D tissue kill to 2-D monolayer kill at matched nominal exposure, which is the product of its transport solve and its cell-state model, with intrinsic potency divided out. Nothing in that prediction was tuned; the transport constant `alpha` was fit once on doxorubicin penetration depth, the cell-state bands come from the oxygen PDE, and the literature band was assembled after the fact.

**What this test is NOT.**

1. **It is a weak pass, and the weakness is measurable.** The temozolomide band spans 2-5x (2.5x wide) and rests on a single paper, while the 2-D assay it bridges from has a ~19x between-study IQR (PMID 34794398). A model would have to be badly wrong to fail.
2. **It constrains only the product.** F_penetration x F_state = 0.2574 is what is tested. The same product is reachable with a much weaker gradient and a much harsher state penalty, or the reverse. **Passing does not validate either factor individually**, and section 2 gives independent literature reason to think our split is tilted too far towards delivery.
3. **The geometries are not the same.** A spheroid is 100-500 um across with a bath at its surface; this is a 1-mm-scale vascularised tissue section with vessel-to-vessel distances up to 786 um and a perivascular source. Our tissue could legitimately sit at the resistant end of the spheroid band, or outside it, without the model being wrong. Equally, this section has **0.0% hypoxic cells** while the hypoxic arm of PMID 30694423 (folds 4x, 5x, >=2x) had 1 % O2 throughout — so the normoxic subset of that study (2-3x) is arguably the closer comparator, and our 3.26x sits just above it.
4. **It is one section, one patient.** No inter-patient variance is tested.
5. **For some drugs the verdict is convention-dependent.** The bath-matching reference is the p99 of tumour-cell exposure. Re-running with the *maximum* exposure as the reference (a defensible alternative, and the one a spheroid's outer surface literally corresponds to) keeps 2 of 3 drugs inside their band; outside under that convention: cisplatin (7.29x vs band 0.4-6x). Temozolomide — the drug with the matched anchor — stays inside under every convention tested (p95 2.53x, p99 3.26x, max 4.35x, band 2-5x).

## 4. The anchor, and the absolute kills

Measured ex-vivo temozolomide response (BJC2023_MOESM2.docx Supplemental Table 2 (primary), n = 40): median viability **82.4%** at 100 uM (IQR 62.1-97.1%, range 16.2-130.9%, bootstrap 95% CI on the median 71.8-88.4%). Operator %CV median 8.3%, max 42.9%. That is a **kill of 0.1759**.

**What the ex-vivo number is allowed to set, and what it is not.** It is a dissociated 2-D monolayer measurement, so it anchors the **2-D-appropriate** quantity — intrinsic potency — and nothing else. The delivery and cell-state factors are left entirely to the model, which is the whole point: calibrating the 3-D prediction to hit the 2-D number would absorb the effect section 3 just tested.

Solving the model's own Hill mapping (hill = 1.2) for the IC50 that reproduces kill 0.1759 at 100 uM gives **362 uM**, against the library's 606 uM (GDSC2 median over 51 glioma lines) — the measured patient-derived cultures are **1.67x more sensitive** than the extrapolated cell-line potency. For scale, that anchored value sits comfortably inside GDSC's own between-line IQR (302-1442 uM), so this is a re-centring, not a contradiction.

**Cisplatin and doxorubicin get no new anchor, and this is stated rather than papered over.** Their library IC50s are *already* 2-D quantities (cisplatin: GDSC2 median over 34 glioma monolayer lines; doxorubicin: a quoted literature range in glioma lines, 0.1-1 uM, not a measurement in this disease). There is no measured GBM ex-vivo viability for either drug in the cache, so their error bars are the *spread of the anchor source*, not a re-anchoring.

| drug | 2-D anchor | anchored IC50 (uM) | band (uM) | 2-D kill at the reference conc. | **3-D kill on this section** | **error bar (between-patient)** | library IC50 (uM) |
|---|---|---:|---|---:|---:|---|---:|
| temozolomide | measured ex-vivo 2-D viability @100 uM (PMID 37620410, n=40) | 362.1 | 151.1-1882 | 0.1760 | **0.0326** | 0.0048-0.0831 | 606 |
| cisplatin | GDSC2 glioma monolayer median IC50 (already a 2-D anchor) | 24.55 | 8.929-84.7 | 0.5000 | **0.4382** | 0.2667-0.5100 | 24.55 |
| doxorubicin | literature 2-D range in glioma lines (~0.1-1 uM); NOT a measured anchor in this disease | 0.5 | 0.1-1 | 0.5000 | **0.2647** | 0.1666-0.4646 | 0.5 |

For temozolomide a **second, tighter bar answers a different question**: the bootstrap 95 % CI on the *cohort median* viability (10,000 resamples of n = 40) implies IC50 218-545 uM and 3-D kill **0.0205-0.0569** (2.8x). Use the wide band when quoting a kill for a patient whose own sensitivity has not been measured; use this narrow one when asking how well the *typical* GBM's intrinsic potency is pinned down.

All 3-D kills are at the repo-nominal plasma level (peak 50 uM, 4 h half-life, 6 h window -> window-mean plasma 31.09 uM), the same axis `results/dose_schedule.md` uses.

**Fast path verified.** The analytic expectation `mean_i p_kill` and the engine's own Bernoulli realisation through `apply_drug_kill` agree to 5.22e-03 absolute (max over 3 drugs), against a one-sigma Bernoulli sampling scale of up to 3.22e-03 on 23,812 tumour cells — i.e. within ~1.6 sigma, consistent with sampling noise and nothing else.
**Verified end to end against the FULL engine.** A complete `sim.engine.simulate_response` run for temozolomide (its own prologue, its own transport solve, its own kill) gives kill **0.0183520914** against this module's 0.0181278103 — absolute delta **2.24e-04**, against a one-sigma Bernoulli sampling scale of 8.65e-04 over 23,812 tumour cells. The engine draws kills stochastically and this module reports the expectation, so agreement at the sampling scale is the correct outcome. Nothing in the engine is re-implemented.

**How to read the error bar.** For temozolomide it is the inter-culture IQR of the measured ex-vivo viability (62.1-97.1 %, i.e. implied IC50 151-1882 uM) propagated through the model's own delivery and state factors. It is a **between-patient** bar, not a confidence interval on one patient: it says that if this Xenium tumour's intrinsic temozolomide sensitivity sits anywhere in the middle half of the measured GBM distribution, its predicted 3-D kill is somewhere in 0.0048-0.0831. On top of that the response variable is itself the operator-*reproducibility* table with %CV up to 43 %, so the bar is if anything optimistic.

**A second, independent bar exists and should be quoted with it.** Taking the anchored 2-D kill and multiplying by the *literature* 3-D:2-D ratio band rather than the model's own factors gives temozolomide 0.0300-0.0850. The model's own prediction (0.0326) lies inside that, which is section 3's verdict restated in absolute units — but note the literature band here is applied to a *bath* concentration and the model's number is at a *plasma* concentration behind a BBB, so they are close for a reason that is partly coincidental. Do not read the agreement as two independent confirmations.

**These are still predictions with nothing to score them against.** The Xenium section is an untreated resection.

## 5. What is still blocking absolute kill

The quiescent label is a **model assignment from the simulated oxygen field, not a measurement** (`results/dose_schedule.md` limitation 9; a concurrent workstream is testing it against measured MKI67). This is how far the absolute kill moves if that label is wrong — everything else (exposure field, sensitivity, anchored IC50) held fixed:

| quiescent fraction (forced) | F_state | 3-D kill (temozolomide, anchored IC50) | x vs the model's own labelling |
|---:|---:|---:|---:|
| 0.0% | 1.0000 | 0.0581 | 1.78x |
| 25.0% | 0.8750 | 0.0533 | 1.64x |
| 50.0% | 0.7500 | 0.0471 | 1.44x |
| 75.0% | 0.6250 | 0.0396 | 1.21x |
| 90.0% | 0.5502 | 0.0340 | 1.04x |
| 100.0% | 0.5000 | 0.0291 | 0.89x |

The sweep assigns the quiescent label to the LEAST-exposed cells (quiescence tracks the same oxygen gradient that tracks perfusion), which is the pessimistic ordering. The model's own labelling gives 0.0326; the swept row at the model's own 90.0% differs slightly because the real labels are not perfectly exposure-ordered.

**The bound.** Because the engine applies cell state as a bounded multiplier (x0.35 hypoxic, x0.5 quiescent, x1 viable), relabelling cells can move `F_state` only between 0.35 and 1.0. Over the full 0-100 % quiescent sweep above, temozolomide's absolute 3-D kill moves **2.00x** — against **17.4x** from the measured inter-patient potency spread in section 4. On this tissue, with 90.0 % of tumour cells already labelled quiescent and 0.0 % hypoxic, the model is near the resistant end of that bounded range, so the realistic remaining movement is smaller still.

**So the honest verdict on step 5 is a partial reversal of the prior expectation.** `results/dose_schedule.md` says the quiescent band is ~89 % of the temozolomide kill ceiling and that absolute kill therefore waits on the cell-state work. On a fresh run at dx = 12 um the ceiling is 0.5502 and the quiescent band is 81.8 % of it — the statement about the *ceiling* stands (the ceiling is definitionally nothing but `sensitivity x resistance`). But the *ceiling* is not the *kill*: at any achievable exposure the Hill term is far from saturated, and the kill is set jointly by potency, delivery and state. Measured against the kill, the cell-state label is the smaller of the two open uncertainties **when the comparison is against between-patient biological spread** (2.00x vs 17.4x) — and roughly EQUAL to it when the comparison is against how well the cohort-typical potency is pinned down (2.00x vs 2.8x, section 4's narrow bar). Neither reading makes cell state irrelevant; both make it non-dominant.

**What remains genuinely blocked, and on what.**

1. **The multipliers, not the labels.** x0.35 and x0.5 are unvalidated priors in `core/config.py:PDConfig` with no citation in this repo. The bound above is conditional on them. If the concurrent cell-state work revises the multiplier rather than the fraction, the bound moves with it. **This is the dependency to track.**
2. **The delivery/state split.** Section 3 constrains only the product. The literature (PMID 2425117, PMID 12569392) says penetration should be the *smaller* term for temozolomide; our model makes it the larger one (56% of the log-penalty). That could mean our gradient is too steep, or our state penalty too mild, or that the comparison is unfair because a spheroid is 200 um across and this section is 1 mm. The ratio test cannot distinguish those; a measured intratumoural drug gradient in human GBM could.
3. **Two of the three drugs have no measured 2-D anchor in this disease.** Doxorubicin's potency is a quoted literature range.
4. **No outcome label.** Untreated resection; zero response data.

**Therefore: comparisons remain the honest primary output.** Absolute kills in section 4 are now *quotable with a propagated error bar and a passed 3-D:2-D sanity test*, which they were not before — but they are predictions with a 17.4x bar and no ground truth, and the cross-drug ordering is still the more defensible claim.

### 5.1 Reconciliation with the concurrent cell-state workstream

`results/cell_state_impact.md` landed independently while this ran, on the same section at the same `dx = 12 um` with the same `alpha`. **The two agree numerically and they must be read together.**

First, the cross-check: that report's oxygen-route kills — temozolomide 0.0184, cisplatin 0.4362, doxorubicin 0.2595 — reproduce this module's engine-path values temozolomide 0.0184, cisplatin 0.4362, doxorubicin 0.2595 bit-for-bit. Two workstreams, one pipeline.

Second, it *measures* the quantity section 5 could only bound. It runs three different cell-state ROUTES and reports the band occupancies:

| cell-state route | what it is | VIABLE | QUIESCENT | HYPOXIC | **F_state** | x vs default |
|---|---|---:|---:|---:|---:|---:|
| `oxygen (shipping default)` | simulated oxygen field, banded | 10.04% | 89.96% | 0.00% | **0.5502** | 1.00x |
| `markers` | measured proliferation markers (MKI67 etc.) | 29.80% | 70.20% | 0.00% | **0.6490** | 1.18x |
| `vessel_distance` | distance-to-vessel banding | 91.08% | 8.67% | 0.26% | **0.9550** | 1.74x |

**Their measured span is inside the bound section 5 derived, which is the result both files needed.** `F_state` moves over 0.5502-0.9550, i.e. 1.00x-1.74x the shipping default — comfortably inside the [0.35, 1.0] range the multipliers permit, and close to the 2.00x this file's own 0-100 % sweep produced. The bound was derived analytically here and measured empirically there, and they agree.

**But their recommendation is stronger than a bound, and it stands.** That report concludes that on this 319-gene panel **no** cell-state route clears its own measurability bar — the oxygen route bands an unfitted simulated field, the marker route sits at a 0.35 % MKI67 detection rate, and the panel carries no hypoxia-response gene at all — and therefore that an absolute kill fraction should not be quoted as a *result* on this section. It also shows the cross-drug **ranking is invariant across all three routes** while the levels are not.

**Nothing in this file overturns that, and this file should not be read as doing so.** What sections 3 and 4 add is narrower and worth stating precisely: (i) the *product* of delivery and state now has an independent, passed, falsifiable check against measured spheroid-vs-monolayer data, which it did not have before; (ii) intrinsic potency now has a measured patient-derived 2-D anchor rather than an extrapolated cell-line one; and (iii) both remaining uncertainties are now **quantified and bounded** rather than open-ended. That converts 'do not quote an absolute kill' into 'quote it only as the band below, never as a point'.

**The combined error bar on absolute kill.** Multiplying the potency band by the measured state-route span:

| drug | 3-D kill (default route, point) | x potency band | x state-route span | **combined band** | total span |
|---|---:|---|---|---|---:|
| temozolomide | 0.0326 | 0.0048-0.0831 | 0.0326-0.0566 | **0.0048-0.1442** | 30.1x |
| cisplatin | 0.4382 | 0.2667-0.5100 | 0.4382-0.7606 | **0.2667-0.8852** | 3.3x |
| doxorubicin | 0.2647 | 0.1666-0.4646 | 0.2647-0.4595 | **0.1666-0.8065** | 4.8x |

The state-route span is applied as a linear scaling on kill (1.00x-1.74x), which is accurate because the Hill term carries no state dependence; checked against `results/cell_state_impact.md` section 3.2 it is ~1.3 % high, so the band is very slightly generous. The potency band is the between-patient one (section 4's wide bar).

**The bottom line for temozolomide: 0.0048-0.1442, a 30x span around a point estimate of 0.0326.** A number with a 30x error bar is not a number you dose a patient from, and it is not a headline result. It *is* now an honest interval with every term in it traceable to either a measurement or a named unvalidated prior — which is the difference between an uncertain result and an unquotable one. **The cross-drug ranking, invariant across all three cell-state routes and across the whole potency band, remains the defensible output of this section.**

## 6. Limitations

1. **No outcome validation.** The Xenium section is an untreated resection. Every absolute kill is a prediction with nothing to score it against.
2. **The 3-D:2-D test is a wide-band test.** The temozolomide literature band (2-5x) rests on a single paper (PMID 30694423) with per-line values read from a figure, and the 2-D assay it bridges from has a ~19x between-study IQR (PMID 34794398).
3. **The sign of the 2-D -> 3-D shift is not universally established** (PMID 27576873 found temozolomide MORE effective in 3-D; PMID 24194875 found Matrigel-3-D more sensitive across 63 active compounds; 5/16 matched pairs in PMID 33395433 had fold <= 1). The band used here is the resistant-direction consensus, which is a choice.
4. **Spheroid geometry != vascularised tissue geometry.** The literature measures 100-500 um spheroids bathed at the surface; this model is a 1-mm-scale section fed from a perivascular source. Matching them at the p99 tissue exposure is a defensible convention, not an equivalence.
5. **The perivascular gradient partly inherits the vascular-source field** (the confound documented in `tmesim/validation/vessel_profile.py` and `results/dose_schedule.md` section 3). F_penetration therefore carries that confound directly.
6. **Cell-state multipliers (x0.35, x0.5) are unvalidated priors** with no citation in this repo. The bound in section 5 is conditional on them.
7. **`mode='uniform'` throughout**, so per-cell intrinsic sensitivity is 1.0 and all intrinsic potency lives in the IC50. Under `mode='determinant'` the marker-derived sensitivity is near-binary at the panel's detection rate and would add noise rather than information to this comparison.
8. **MGMT is not applied.** `results/mgmt_scalar.md` shows a measured MGMT level would scale temozolomide kill down by a further 2.9-5.0x on real sections. The absolute temozolomide kill here is therefore an UPPER BOUND in the MGMT dimension, on top of its stated error bar.
9. **One calibration constant** (`alpha`, fit once on doxorubicin penetration depth) sets the perivascular sink for every drug, so F_penetration inherits its uncertainty for all three.
10. **Single section, single patient, 2-D grid.** No inter-patient variance, no out-of-plane transport.
11. **The bath-matching convention matters for some drugs.** Under a maximum-exposure reference rather than the p99, outside under that convention: cisplatin (7.29x vs band 0.4-6x). Section 3 reports this per drug rather than choosing the flattering convention.
12. **Bookkeeping note.** `results/dose_schedule.md` is internally inconsistent on the cell-state fractions: section 1.3 reports 90.0 % quiescent / 0.0 % hypoxic (ceiling 0.5502) while section 6 reports 86 % / 13.4 % (ceiling 0.4828). This run reproduces section 1.3 exactly — 90.0 % quiescent, 0.0 % hypoxic, ceiling 0.5502, and temozolomide kill 0.018352 at the nominal dose, bit-identical to that file's 1x row. Section 6 of that file is stale.

## References

*Literature retrieved from PubMed; DOI links given where the record carries one.*

* Musah-Eroje A & Watson S. J Neurooncol 2019;142(2):231-240. PMID **30694423** — [DOI](https://doi.org/10.1007/s11060-019-03107-0)
* Gomez-Roman N, et al. Neuro Oncol 2017;19(2):229-241. PMID **27576873** — [DOI](https://doi.org/10.1093/neuonc/now164)
* Samiei E, et al. Int J Mol Sci 2020;21(9):3162. PMID **32365781** — [DOI](https://doi.org/10.3390/ijms21093162)
* Gunther W, et al. Br J Cancer 2003;88(3):463-469. PMID **12569392** — [DOI](https://doi.org/10.1038/sj.bjc.6600711)
* Poon MTC, et al. BMC Cancer 2021;21(1):1240. PMID **34794398** — [DOI](https://doi.org/10.1186/s12885-021-08972-5)
* Koch J, et al. PLoS One 2021;16(1):e0244513. PMID **33395433** — [DOI](https://doi.org/10.1371/journal.pone.0244513)
* Trigub MM, et al. Bull Exp Biol Med 2020;168(4):583-588. PMID **32152850** — [DOI](https://doi.org/10.1007/s10517-020-04757-1)
* Jubelin C, et al. Front Bioeng Biotechnol 2023;11:1260049. PMID **37869710** — [DOI](https://doi.org/10.3389/fbioe.2023.1260049)
* Mutahar AZI, et al. J Transl Med 2025;23(1):1428. PMID **41422245** — [DOI](https://doi.org/10.1186/s12967-025-07585-x)
* Hoda MA, et al. Mol Cancer Ther 2016;15(10):2357-2369. PMID **27512118** — [DOI](https://doi.org/10.1158/1535-7163.MCT-15-0846)
* Imamura Y, et al. Oncol Rep 2015;33(4):1837-1843. PMID **25634491** — [DOI](https://doi.org/10.3892/or.2015.3767)
* Hongisto V, et al. PLoS One 2013;8(10):e77232. PMID **24194875** — [DOI](https://doi.org/10.1371/journal.pone.0077232)
* Durand RE. J Natl Cancer Inst 1986;77(1):247-252. PMID **2425117** — (no DOI -- 1986)
* Kobayashi H, et al. Proc Natl Acad Sci USA 1993;90(8):3294-3298. PMID **8475071** — [DOI](https://doi.org/10.1073/pnas.90.8.3294)
* Desoize B & Jardillier J. Crit Rev Oncol Hematol 2000;36(2-3):193-207. PMID **11033306** — [DOI](https://doi.org/10.1016/s1040-8428(00)00086-x)
* Grantab R, Sivananthan S, Tannock IF. Cancer Res 2006;66(2):1033-1039. PMID **16424039** — [DOI](https://doi.org/10.1158/0008-5472.CAN-05-3077)
* Minchinton AI & Tannock IF. Nat Rev Cancer 2006;6(8):583-592. PMID **16862189** — [DOI](https://doi.org/10.1038/nrc1893)
* Ntafoulis I, *et al.* Br J Cancer 2023;129(8):1327-1338. PMID **37620410**. GEO **GSE232173** (the 2-D ex-vivo anchor).
* Yang W, *et al.* Nucleic Acids Res 2013;41:D955-D961 (GDSC). PMID **23180760**.
* Prior work in this repo: `results/gbm_tmz_exvivo.md` (the measurement), `results/mgmt_scalar.md` (the intrinsic-resistance scalar), `results/dose_schedule.md` (the kill ceiling and its cell-state decomposition), `results/numerical_rigor.md` (why dx = 12 um).

*Generated by `tmesim/validation/kill_calibration.py` — one reconstruction, 3 transport solves, 1.5 min, seed 0, deterministic.*