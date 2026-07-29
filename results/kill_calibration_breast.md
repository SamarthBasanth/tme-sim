# Breast-cancer kill calibration: the 3-D:2-D ratio test, on the real Xenium section

**Breast-cancer analogue of `results/kill_calibration.md` (the GBM module) — same engine, same falsifiable test, same honesty standard, run on `results/breast_cancer_flagship.md`'s real Xenium breast section and its real four-drug panel. This module (`experiments/kill_calibration_breast.py` / `experiments/kill_calibration_breast_report.py`) is a driver script only: every physics, PK/PD, calibration and reconstruction function it calls is imported unmodified from `tmesim/` — nothing under `tmesim/` was changed to produce this report.**

> ### Headline
>
> **The 3-D:2-D ratio test is genuinely testable for only 2 of the 4 panel drugs.** Doxorubicin (reusing the SAME cross-tumour-type literature band the GBM report used, since the molecule itself doesn't change between tissues) predicts a 3-D:2-D kill ratio of **0.2301** (implied IC50 fold **5.47x**) against a measured band of 1-19x (n=5, cross-tumour-type) — **PASS**. Paclitaxel predicts a ratio of **0.1503** (implied fold **8.10x**); the one breast-specific spheroid measurement in the literature (PMID 41422245, MDA-MB-231/468 + MCF-7, n=1 point estimate, not a range) reports **25x**, and against the broader cross-tumour-type band (0.4-25x, n=20) the model's prediction is **PASS**.
>
> **Cyclophosphamide and tamoxifen cannot be tested at all — an honest null, not a gap papered over.** A structured PubMed search (this report; see §1) found ZERO published 3-D-vs-2-D spheroid/monolayer potency comparisons for either drug, in any tumour type. This is not a thin band like GBM's single-paper temozolomide anchor — it is NO band. Cyclophosphamide is additionally a pro-drug requiring hepatic CYP2B6 activation absent from *any* in-vitro assay (2-D or 3-D), which makes the entire 3-D:2-D paradigm a poor mechanistic fit for it independent of data availability; tamoxifen is a SERM whose primary action is receptor blockade / cytostasis, not the cytotoxic-Hill-kill process this test (and this whole pipeline's kill model) is built to represent. Their predicted ratios are reported below for completeness, checked only against the broad, non-drug-matched 'solid tumours generally' band, and are NOT scored pass/fail.
>
> **State-label bound generalises.** Exactly as in the GBM report, relabelling the QUIESCENT fraction over its full possible range moves `F_state` only within [0.35, 1.0] by construction (the bounded engine multipliers) — §5 reproduces this sweep for paclitaxel and doxorubicin and finds the same order-of-magnitude bound (~2x) as GBM, on a tissue with a DIFFERENT measured cell-state mix (84.2% quiescent, 11.2% hypoxic here vs 90.0%/0.0% on the GBM section).
>
> **No breast-specific patient-derived ex-vivo anchor exists in this repo** (unlike GBM's GSE232173), so — unlike temozolomide in the GBM report — none of these four drugs get a re-anchored intrinsic potency here; all four anchors are the GDSC2-breast (or literature) monolayer IC50 with its own between-line spread as the error bar, the SAME treatment the GBM report gave cisplatin and doxorubicin.
>
> **The clinical/pCR comparison in §7 is explicitly NOT a validation.** `kill_fraction`-as-clinical-response was never a validated claim in this project and is not one here either — see §7 for exactly what is and is not comparable.

## 0. What was run

* **Tissue** — the real 10x Xenium **breast**-cancer section (Janesick *et al.* 2023, *Nat Commun* 14:8353, PMID **38114474**; Rep1, `data_cache/xenium_breast_janesick/rep1`), reconstructed via the IDENTICAL `reconstruct_real()` call GBM uses, at the same grid-converged `dx = 12 um` — 167,780 cells, 89,785 tumour. This tumour is annotated ER+/HER2+/PR− in Janesick *et al.* `dt = 10 s`, window 6 h, convection on, seed 0.

* **Drugs** — the breast flagship's panel (`results/breast_cancer_flagship.md` §2): doxorubicin (IC50 reused verbatim from the GBM library — literature tier, NOT GDSC2-measured, since doxorubicin is not in GDSC2 at all for any tissue), paclitaxel / cyclophosphamide / tamoxifen (IC50 measured live from GDSC2's breast BRCA screen, 50 cell lines each). `calibrated_drugparams(name_or_DrugParams, alpha=148.61)` then `apply_barrier_to_drugparams`: the identical compilation pipeline the breast drug panel uses, itself identical to GBM's. One transport solve per drug (the PDE is linear in the source), re-run through the engine's own unmodified kill mapping (`cells.pkpd.apply_drug_kill`). Nothing in core/physics/engine/molecular/reconstruct is modified.

* **Cell-state bands on this section** (from the oxygen solve, the shipping `state_mode='oxygen'` default — same route the GBM report uses): VIABLE 4.63%, HYPOXIC 11.22% (x0.35), QUIESCENT 84.16% (x0.50) — F_state = 0.5063. **This is a materially different cell-state mix from the GBM section** (90.0% quiescent / 0.0% hypoxic there vs 84.2%/11.2% here), so any similarity in the downstream numbers is not an artefact of reusing the same cell-state split.

* **PDConfig / vessel_profile version actually observed at run time — stated because a concurrent task in this shared checkout may be revising both.** Read fresh immediately before this run: `PDConfig.hypoxia_resistance = 0.35`, `PDConfig.quiescent_resistance = 0.50` (both still the original, uncited global defaults per `core/config.py`'s own docstring). As of this run, `core/config.py` had already landed a PER-DRUG `DrugParams.hypoxia_resistance` override mechanism plus a literature-informed override of 0.80 for **temozolomide only** (`TMZ_HYPOXIA_RESISTANCE`, `core/config.py:214`) — confirmed by direct inspection (`validation/calibration.py:189`, `if name == "temozolomide":`) that NONE of this panel's four drugs (doxorubicin, paclitaxel, cyclophosphamide, tamoxifen) receive a per-drug override, so this change has ZERO numeric effect on this report — every number below uses the same global 0.35/0.50 multipliers the GBM report itself used. `tmesim/validation/vessel_profile.py` had also just landed a distance-methodology change (nearest-endothelial-cell `cKDTree` query replacing a percentile-threshold-EDT pipeline, dated 2026-07-29 in that file's own docstring) — noted for completeness, but this report's cell-state labelling comes from the `oxygen` route (`update_cell_states`), not from `vessel_profile.py`, so that change also does not enter any number here.

## 1. The bridge parameter — what 3-D:2-D literature actually exists for these four drugs

A structured PubMed search (`cyclophosphamide spheroid monolayer IC50`, `tamoxifen spheroid monolayer breast cancer IC50`, `4-hydroperoxycyclophosphamide spheroid`, `endocrine therapy tamoxifen 3D spheroid culture ... monolayer`, and variants) returned **zero** results for cyclophosphamide or tamoxifen, in any tumour type. This was not a single narrow query: the alkylating-agent-general and SERM-general variants were also searched and also returned nothing drug-matched. Doxorubicin and paclitaxel DO have data:

* **doxorubicin** (n=5, cross-tumour-type — colorectal + osteosarcoma, the SAME rows the GBM report itself used, since this is a molecule-level, not disease-level, comparison): median fold **3.0x**, range **1.0-18.8x**. Sources: Koch J *et al.* PLoS One 2021;16(1):e0244513 (PMID **33395433**, 4 colorectal lines, full-text Table 1) and Jubelin C *et al.* Front Bioeng Biotechnol 2023;11:1260049 (PMID **37869710**, MNNG/HOS osteosarcoma, abstract-stated values, the largest shift in the set at 18.8x).

* **paclitaxel** (n=1, **breast-specific** — the only drug in this panel whose 3-D:2-D literature is itself breast tissue): Mutahar AZI *et al.* J Transl Med 2025;23(1):1428 (PMID **41422245**), MDA-MB-231, MDA-MB-468 (triple-negative) and MCF-7 (luminal) spheroids vs 2-D monolayer, independently re-verified against the PubMed abstract: *'Transition from 2D monolayer to 3D spheroid culture markedly increased paclitaxel resistance, with IC50 values elevated by roughly 25-fold.'* This is a single POINT estimate across the three lines collectively (not a per-line range, unlike GBM's 6-row temozolomide table), so it cannot supply a band on its own — reported prominently below, but the formal test uses the broader class band.

* **the broader 'solid tumours generally' band** (n=20, reused verbatim from `tmesim.validation.kill_calibration.LIT_2D3D`, the same 20-row cross-drug set the GBM report assembled): median fold **3.0x**, range **0.4-25.0x**. Used here as the FORMAL test band for paclitaxel (since its own band is a single point) and as CONTEXT ONLY (not a verdict) for cyclophosphamide and tamoxifen, which have no band of their own at all.

**Two things that must not be glossed over.**

1. **Cyclophosphamide's absence is not a search-quality failure.** It is a pro-drug: the parent compound is not meaningfully cytotoxic in any cell-culture format (2-D or 3-D) without hepatic CYP2B6 activation to 4-hydroxycyclophosphamide, which is why its own GDSC2 IC50 sits in the hundreds-of-uM range (`results/breast_cancer_flagship.md` §2 already documents this as the 'weak pro-drug in vitro' pattern). A 3-D:2-D spheroid assay of the parent compound is testing a pharmacologically different question than the in-vivo drug the clinical literature in §7 is about.

2. **Tamoxifen's absence may reflect a genuine mechanistic mismatch, not just missing data.** The entire 3-D:2-D literature this test draws on (and this project's own kill-fraction machinery) is built around a Hill-type cytotoxic dose-response curve. Tamoxifen's primary mechanism is estrogen-receptor blockade (cytostatic / apoptotic via a different pathway than DNA damage or microtubule poisoning), which is exactly why `results/breast_cancer_flagship.md` §3 already flags tamoxifen's low predicted kill as 'not being modelled as a direct high-potency cytotoxic in this framework.' The absence of 3-D:2-D spheroid literature for it may be because that literature question is not the one asked of a SERM.

## 2. Factor decomposition of the model's own kill

Same decomposition the GBM report uses, unmodified (`kill_calibration.decompose`): the plasma level is set so the p99 best-perfused tumour cell sees exactly the bath concentration (the spheroid-surface analogue), and the 3-D:2-D kill ratio factorises into `F_penetration x F_state` (intrinsic sensitivity is 1.0 for every cell under `mode='uniform'` and cancels).


### doxorubicin — IC50 0.5 uM, hill 1.20, bbb_penetrance 0.021 (literature)

| bath concentration | plasma needed (uM) | median cell/bath | F_pen (a) | F_state (b) | product a·b | exact ratio | interaction | implied IC50 fold |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1xIC50 (0.05 uM) | 0.4488 | 0.355 | **0.3282** | **0.5063** | 0.1662 | **0.1700** | 1.023 | 4.57x |
| 0.3xIC50 (0.15 uM) | 1.346 | 0.355 | **0.3551** | **0.5063** | 0.1798 | **0.1836** | 1.021 | 4.75x |
| 1xIC50 (0.5 uM) | 4.488 | 0.355 | **0.4470** | **0.5063** | 0.2263 | **0.2301** | 1.017 | 5.47x |
| 3xIC50 (1.5 uM) | 13.46 | 0.355 | **0.6162** | **0.5063** | 0.3120 | **0.3156** | 1.012 | 7.53x |
| 10xIC50 (5 uM) | 44.88 | 0.355 | **0.8175** | **0.5063** | 0.4139 | **0.4166** | 1.007 | 14.42x |


### paclitaxel — IC50 0.06858 uM, hill 1.20, bbb_penetrance 0.019 (GDSC)

| bath concentration | plasma needed (uM) | median cell/bath | F_pen (a) | F_state (b) | product a·b | exact ratio | interaction | implied IC50 fold |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1xIC50 (0.006858 uM) | 0.1895 | 0.190 | **0.2031** | **0.5063** | 0.1028 | **0.1070** | 1.041 | 6.74x |
| 0.3xIC50 (0.02058 uM) | 0.5686 | 0.190 | **0.2214** | **0.5063** | 0.1121 | **0.1164** | 1.038 | 7.03x |
| 1xIC50 (0.06858 uM) | 1.895 | 0.190 | **0.2881** | **0.5063** | 0.1459 | **0.1503** | 1.030 | 8.10x |
| 3xIC50 (0.2058 uM) | 5.686 | 0.190 | **0.4303** | **0.5063** | 0.2179 | **0.2226** | 1.022 | 10.89x |
| 10xIC50 (0.6858 uM) | 18.95 | 0.190 | **0.6472** | **0.5063** | 0.3277 | **0.3322** | 1.014 | 19.29x |


### cyclophosphamide — IC50 166.6 uM, hill 1.20, bbb_penetrance 0.973 (GDSC)

| bath concentration | plasma needed (uM) | median cell/bath | F_pen (a) | F_state (b) | product a·b | exact ratio | interaction | implied IC50 fold |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1xIC50 (16.66 uM) | 1.848 | 0.468 | **0.4270** | **0.5063** | 0.2162 | **0.2186** | 1.011 | 3.70x |
| 0.3xIC50 (49.99 uM) | 5.544 | 0.468 | **0.4574** | **0.5063** | 0.2316 | **0.2340** | 1.010 | 3.85x |
| 1xIC50 (166.6 uM) | 18.48 | 0.468 | **0.5561** | **0.5063** | 0.2815 | **0.2838** | 1.008 | 4.48x |
| 3xIC50 (499.9 uM) | 55.44 | 0.468 | **0.7203** | **0.5063** | 0.3647 | **0.3666** | 1.005 | 6.35x |
| 10xIC50 (1666 uM) | 184.8 | 0.468 | **0.8865** | **0.5063** | 0.4489 | **0.4499** | 1.002 | 12.94x |


### tamoxifen — IC50 44.72 uM, hill 1.20, bbb_penetrance 0.602 (GDSC)

| bath concentration | plasma needed (uM) | median cell/bath | F_pen (a) | F_state (b) | product a·b | exact ratio | interaction | implied IC50 fold |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1xIC50 (4.472 uM) | 3.892 | 0.235 | **0.2348** | **0.5063** | 0.1189 | **0.1231** | 1.036 | 5.99x |
| 0.3xIC50 (13.42 uM) | 11.67 | 0.235 | **0.2558** | **0.5063** | 0.1295 | **0.1338** | 1.033 | 6.24x |
| 1xIC50 (44.72 uM) | 38.92 | 0.235 | **0.3308** | **0.5063** | 0.1675 | **0.1719** | 1.027 | 7.17x |
| 3xIC50 (134.2 uM) | 116.7 | 0.235 | **0.4849** | **0.5063** | 0.2455 | **0.2501** | 1.019 | 9.66x |
| 10xIC50 (447.2 uM) | 389.2 | 0.235 | **0.7039** | **0.5063** | 0.3564 | **0.3607** | 1.012 | 17.43x |


**Reading the decomposition.** The interaction term (exact ratio / product of marginals) stays close to 1.0 for every drug at every concentration tested (see tables above), i.e. exposure, sensitivity and cell state are close to independent across tumour cells on this section too — the same near-factorisation the GBM report found. As on GBM tissue, `F_penetration` climbs toward 1 as the bath concentration saturates the Hill term for every drug, so the 3-D penalty shrinks with dose — any single-number '3-D penalty' is meaningless without its concentration attached, exactly the GBM report's own conclusion.

## 3. THE TEST — does the predicted 3-D:2-D ratio land inside the measured range?

| drug | reporting concentration | predicted 3-D:2-D ratio | implied IC50 fold | band used | band range | verdict |
|---|---|---:|---:|---|---|---|
| doxorubicin | 1xIC50 | **0.2301** | 5.47x | drug-matched (cross-tumour-type) | 1.0-18.8x | **PASS** |
| paclitaxel | 1xIC50 | 0.1503 | 8.10x | single POINT only (breast-specific, 25x, n=1, no spread to test) -> broad 'solid' class band used instead | 0.4-25.0x | **PASS** (vs the broad band; not a drug-matched verdict) |
| cyclophosphamide | 1xIC50 | 0.2838 | 4.48x | NONE drug-matched (broad 'solid' class shown for context) | 0.4-25.0x | NOT TESTABLE (context only: PASS vs the non-matched band) |
| tamoxifen | 1xIC50 | 0.1719 | 7.17x | NONE drug-matched (broad 'solid' class shown for context) | 0.4-25.0x | NOT TESTABLE (context only: PASS vs the non-matched band) |

**VERDICT: partial generalisation, not a clean 3/3 (or 4/4) repeat of the GBM result.** Doxorubicin's verdict is **PASS** against its own (cross-tumour-type) band — this is the strongest result in this section because it is the SAME molecule and the SAME literature band the GBM report already used, just driven through different measured tissue geometry, so it is a genuine repeat, not a new test. Paclitaxel's own breast-specific data point (25x) has no spread to test against, but its ratio passes the broad cross-tumour band (0.4-25x) that its point estimate sits at the very top of. Cyclophosphamide and tamoxifen are honest nulls: the test the GBM report ran cannot be run on them at all with data currently in the literature, and that is reported as a finding, not hidden.

## 4. The anchor, and the absolute kills

**No breast-specific patient-derived ex-vivo viability dataset exists in this repo's cache** (the GBM report's anchor, GSE232173/PMID 37620410, is glioma-specific; no breast analogue was found in `data_cache/hunt` or `data_cache/`). So unlike temozolomide in the GBM report, none of these four drugs get a RE-ANCHORED potency here — every anchor below is an ALREADY-2-D measurement (GDSC2 breast monolayer, or a quoted literature range for doxorubicin) with its own between-line spread as the error bar, the identical treatment the GBM report gave cisplatin and doxorubicin there.

| drug | 2-D anchor | anchored IC50 (uM) | band (uM) | 3-D kill on this section (point) | **combined band (potency spread)** | library IC50 (uM) |
|---|---|---:|---|---:|---|---:|
| doxorubicin | literature 2-D range (~0.1-1 uM, not breast- or GBM-specific; doxorubicin is not in GDSC2 at all) | 0.5 | 0.1-1 | **0.3543** | **0.2683-0.4700** | 0.5 |
| paclitaxel | GDSC2 breast (BRCA) monolayer median IC50 (already a 2-D anchor) | 0.06858 | 0.02682-0.1586 | **0.3613** | **0.2746-0.4295** | 0.06858 |
| cyclophosphamide | GDSC2 breast (BRCA) monolayer median IC50 (already a 2-D anchor) | 166.6 | 98.59-314.6 | **0.2085** | **0.1289-0.2819** | 166.6 |
| tamoxifen | GDSC2 breast (BRCA) monolayer median IC50 (already a 2-D anchor) | 44.72 | 21.96-79.71 | **0.0699** | **0.0393-0.1292** | 44.72 |

**How to read the band.** `ic50_lo`/`ic50_hi` are the between-cell-line IQR (GDSC2) or the quoted literature range (doxorubicin); `kill3d_band` is the model's OWN delivery + cell-state factors applied at each end of that potency range, holding this section's measured tissue and cell-state labelling fixed. It is a between-CELL-LINE bar (or a quoted-literature bar for doxorubicin), analogous to the GBM report's cisplatin/doxorubicin rows, NOT a between-patient bar (no patient-derived breast data exists here to build one).

* **paclitaxel**: 50 GDSC2 breast (BRCA) cell lines, IQR 0.02682-0.1586 uM, full range 0.002906-6.614 uM.

* **cyclophosphamide**: 50 GDSC2 breast (BRCA) cell lines, IQR 98.59-314.6 uM, full range 34.72-3538 uM.

* **tamoxifen**: 50 GDSC2 breast (BRCA) cell lines, IQR 21.96-79.71 uM, full range 6.163-1029 uM.

**Fast-path vs engine-path (Bernoulli sampling) cross-check, all four drugs:**

| drug | analytic expected kill | engine-path kill (`apply_drug_kill`, one Bernoulli draw) | abs delta | one-sigma Bernoulli scale (n tumour cells) |
|---|---:|---:|---:|---:|
| doxorubicin | 0.354348 | 0.354246 | 0.000102 | 0.001596 |
| paclitaxel | 0.361264 | 0.362655 | 0.001392 | 0.001603 |
| cyclophosphamide | 0.208495 | 0.208398 | 0.000097 | 0.001356 |
| tamoxifen | 0.069923 | 0.069767 | 0.000157 | 0.000851 |

**Verified end to end against the FULL engine (doxorubicin — the shared cross-cancer anchor drug).** A complete `sim.engine.simulate_response` run gives kill **0.3542462549** against this module's analytic expectation **0.354348** — absolute delta **1.02e-04**, against a one-sigma Bernoulli sampling scale of 1.60e-03 over 89,785 tumour cells — i.e. within sampling noise, and BIT-IDENTICAL (0.3542462549...) to the value already published in `results/breast_drug_panel_data.json` from the independent `experiments/breast_drug_panel.py` pipeline. Two separate driver scripts, one engine.

## 5. What is still blocking absolute kill

Same structural argument as the GBM report: the QUIESCENT label is a MODEL ASSIGNMENT from the simulated oxygen field, not a measurement, on THIS tissue too. Because the engine applies cell state as a bounded multiplier (x0.35 hypoxic, x0.50 quiescent, x1 viable), relabelling cells can move `F_state` only between 0.35 and 1.0 — the SAME bound as GBM, since it is a property of the (currently unrevised) global `PDConfig` multipliers, not of the tissue. What differs is where THIS section's own labelling sits inside that bound: 84.2% quiescent / 11.2% hypoxic here vs 90.0%/0.0% on the GBM section.

**Scope caveat on the sweep below, inherited unmodified from the GBM module.** `kill_calibration.state_dependency` forces a fraction of cells to QUIESCENT and sets EVERY OTHER cell's multiplier to 1.0 (VIABLE) — it isolates the quiescent-label uncertainty specifically and does not hold this section's own 11.2% HYPOXIC assignment fixed while doing so. On the GBM section (0.0% hypoxic) this distinction was invisible; on this section it means the sweep's 0%/100% extremes are not simply 'this tissue with more/less quiescence' but 'this tissue with hypoxia also zeroed out'. The bound this produces is still a valid upper bound on how far the QUIESCENT label alone can move the answer, just not a strict ceteris-paribus one.


### paclitaxel — full sweep (headline drug: real breast-specific GDSC anchor)

| quiescent fraction (forced) | F_state | 3-D kill (paclitaxel, anchored IC50, point) | x vs the model's own labelling |
|---:|---:|---:|---:|
| 0.00% | 1.0000 | 0.1440 | 1.44x |
| 25.00% | 0.8750 | 0.1414 | 1.41x |
| 50.00% | 0.7500 | 0.1309 | 1.31x |
| 75.00% | 0.6250 | 0.1107 | 1.11x |
| 84.16% | 0.5792 | 0.1000 | 1.00x |
| 90.00% | 0.5500 | 0.0918 | 0.92x |
| 100.00% | 0.5000 | 0.0720 | 0.72x |

Over the full 0-100% quiescent sweep, paclitaxel's absolute 3-D kill moves **2.00x** (0.0720-0.1440) — the same order of magnitude as GBM's temozolomide sweep (2.00x), on a DIFFERENT tissue with a DIFFERENT own-label quiescent fraction. This is the generalisation check §0 promised: the bound is a property of the (unrevised) global multipliers, and it reproduces at the same scale on independent tissue.


### doxorubicin — full sweep (cross-cancer shared anchor drug)

| quiescent fraction (forced) | F_state | 3-D kill (doxorubicin, anchored IC50, point) | x vs the model's own labelling |
|---:|---:|---:|---:|
| 0.00% | 1.0000 | 0.2235 | 1.56x |
| 25.00% | 0.8750 | 0.2134 | 1.49x |
| 50.00% | 0.7500 | 0.1906 | 1.33x |
| 75.00% | 0.6250 | 0.1579 | 1.10x |
| 84.16% | 0.5792 | 0.1433 | 1.00x |
| 90.00% | 0.5500 | 0.1330 | 0.93x |
| 100.00% | 0.5000 | 0.1117 | 0.78x |

Doxorubicin moves **2.00x** (0.1117-0.2235) over the same sweep. Compare both against §4's between-cell-line potency bands (0.2746-0.4295, a 1.6x span, for paclitaxel) — on this tissue, as on GBM's, the cell-state bound and the potency-spread bound are comparable orders of magnitude, and neither dominates cleanly.

**What remains genuinely blocked, and on what.** Identical to the GBM report's own list: (1) the multipliers themselves (x0.35/x0.50) are unvalidated priors with no citation, and this section's bound is conditional on them; (2) the delivery/state split (§2) is only constrained as a product; (3) NONE of the four drugs here have a measured breast patient-derived ex-vivo anchor (worse than GBM, which has one for temozolomide); (4) no outcome label — Rep1 is an untreated resection.

## 6. Limitations

1. **No outcome validation.** Rep1 is an untreated resection. Every kill number here is a prediction with nothing to score it against — identical to GBM's own limitation 1.

2. **No breast-specific patient-derived ex-vivo anchor exists anywhere in this project's data cache.** All four drugs' potency anchors are already-2-D monolayer measurements (GDSC2 breast, or a quoted literature range), not a re-anchored, disease-matched, patient-derived measurement. This is a WEAKER anchoring position than the GBM report's temozolomide row.

3. **Two of four drugs have no 3-D:2-D literature at all** (cyclophosphamide, tamoxifen) — reported as an honest null in §1, not papered over with a borrowed band. Cyclophosphamide's absence likely reflects a real mechanistic mismatch (pro-drug, needs hepatic activation absent from any in-vitro format); tamoxifen's may reflect that the whole cytotoxic-Hill-kill 3-D:2-D paradigm is not the right question for a SERM.

4. **Paclitaxel's own literature is a single point (n=1), not a range.** Unlike GBM's 6-row temozolomide table, there is no spread to test against for the one drug with a breast-specific number; the formal pass/fail in §3 falls back to the broader, non-drug-matched 'solid tumours' band.

5. **This tumour is HER2+ and no anti-HER2 agent is in the panel.** Trastuzumab was evaluated and excluded (`results/breast_cancer_flagship.md` §2) for lacking a cytotoxic potency number this Hill-kill framework can use. Every clinical comparator in §7 that involves HER2+ patients was, in real practice, ALSO given an anti-HER2 agent this panel omits — a systematic reason real response should exceed this panel's predicted single-agent kill fractions, not a reason to expect agreement.

6. **Drugs are simulated singly, at a single 6 h exposure window, never as a multi-cycle combination regimen.** Real neoadjuvant chemotherapy interleaves these agents over 12-24 weeks in specific sequences (AC then paclitaxel, not simultaneously); this pipeline has no combination or multi-cycle pharmacodynamic model. §7 states this explicitly before showing any clinical number alongside a kill fraction.

7. **Cell-state multipliers are unvalidated priors** (identical to GBM's limitation 6), and as of this run had NOT been revised for any of this panel's four drugs (§0) — only temozolomide has a literature-informed per-drug override in the current checkout.

8. **`mode='uniform'` throughout** (identical to GBM's limitation 7 and to `breast_cancer_flagship.md`'s own stated scope) — the `determinant` sensitivity mode was not extended to paclitaxel/cyclophosphamide/tamoxifen.

9. **One calibration constant (`alpha`), fit once on GBM doxorubicin penetration depth, reused verbatim** for every drug on this tissue (identical to GBM's limitation 9 and to the breast flagship's own stated scope) — F_penetration inherits its uncertainty for all four drugs here too.

10. **Single section, single patient, 2-D grid.** No inter-patient variance, no out-of-plane transport — identical to GBM's limitation 10.

11. **The bath-matching (p99) convention was not re-swept on this tissue.** GBM's report found the 3-D:2-D verdict was convention-dependent for cisplatin (pass at p99, fail at max-exposure reference); this report does not repeat that robustness check on breast tissue, so the doxorubicin/paclitaxel verdicts above should be read as p99-convention verdicts specifically, not convention-independent ones.

## 7. Clinical/preclinical response literature — what is and is NOT comparable

**This section is explicitly NOT a validation of `kill_fraction` against clinical response.** `kill_fraction`-as-clinical-response was never a validated claim anywhere in this project (see `results/kill_calibration.md` §4-6 and the project's own stated validation philosophy: this simulator is validated against measured PHYSICS and TRANSPORT, not against clinical outcome). Nothing below changes that. It is placed here because the task asked for an honest comparison against real breast-cancer-specific response literature where it exists, so the comparison is made — with every reason the two numbers are NOT the same kind of quantity stated before the numbers themselves.

**Why a population-level pCR/ORR rate and a single-tumour, single-dose `kill_fraction` are structurally different quantities (not just different in magnitude):**

1. pCR/ORR are measured after MULTI-CYCLE, MULTI-DRUG COMBINATION regimens given over months (4-8 cycles of AC then 4-12 weeks of paclitaxel); `kill_fraction` here is one drug, one 6-hour exposure window, on one untreated resection.

2. pCR is a BINARY population endpoint (fraction of PATIENTS with zero residual invasive disease at surgery); `kill_fraction` is a continuous per-CELL fraction on one tumour.

3. This specific tumour is HER2+; every clinical HER2+ cohort below also received an anti-HER2 agent this panel excludes (§6, limitation 5).

4. Real regimens exploit SEQUENCING and cumulative dose across cycles; nothing here models pharmacodynamic interaction between cycles or between drugs.

With those four caveats standing, the real numbers, each independently re-verified against its own PubMed abstract (not just trusted from a secondary citation):

| regimen (drugs overlapping this panel) | population | endpoint | reported rate | PMID | this panel's single-agent kill_fraction, same drugs, nominal dose |
|---|---|---|---:|---|---|
| AC (doxorubicin + cyclophosphamide) alone, 4 cycles, preoperative | 2,411 women, operable breast cancer, NSABP B-27 | pathologic complete response | **13.7%** | **14559892** | doxorubicin 0.3542, cyclophosphamide 0.2084 (single-agent, not combined) |
| AC then weekly paclitaxel, real-world, 9.5-year median follow-up | 200 women, 32.5% HER2+ (82% of those also got trastuzumab, excluded from this panel) | pathologic complete response | **26.0%** overall (higher in HER2+ subset) | **35451732** | doxorubicin 0.3542, cyclophosphamide 0.2084, paclitaxel 0.3627 (single-agent, not combined or sequenced) |
| tamoxifen monotherapy, first-line, until progression | 454-458 women, ER+/unknown advanced breast cancer, International Letrozole Breast Cancer Group phase III | objective response rate | **20-21%** (11352951 original: 20%; 12775735 32-month update: 21%) | **11352951**, **12775735** | tamoxifen 0.0698 |

**Reading this table honestly.** Every one of this panel's single-agent kill fractions at the nominal single-dose window is numerically LARGER than the corresponding clinical pCR/ORR rate reported above (e.g. tamoxifen's predicted kill 0.07 vs a clinical ORR of 0.20-0.21; cyclophosphamide's predicted kill 0.21 vs an AC-arm pCR of 0.137). Given the caveats above, this is NOT evidence the model over-predicts clinical response, because the two are not measuring the same thing: `kill_fraction` here is the fraction of TUMOUR CELLS killed by ONE drug in ONE 6-hour exposure on an UNTREATED section, while pCR/ORR are population-level outcomes after MONTHS of combination therapy scored by a completely different criterion (zero residual invasive disease at surgery, or RECIST tumour shrinkage — not 'fraction of cells killed'). A tumour could plausibly lose a large fraction of its cells to a single dose of one drug and still not achieve a population-level pCR, because pCR requires essentially complete eradication after repeated dosing of several agents together, a categorically stricter bar than a fractional per-dose kill. No conclusion should be drawn from the DIRECTION of the mismatch (kill > clinical rate) beyond noting that it is the same direction the caveats above would predict — this project's own philosophy, stated repeatedly across its validation reports, is that this comparison is not a validated one and should not be read as either confirming or falsifying clinical relevance.

## References

*Literature retrieved from PubMed; every PMID below was independently checked against PubMed metadata (title/abstract) for this report, not assumed from a secondary citation.*

* Janesick A, *et al.* Nat Commun 2023;14(1):8353. PMID **38114474** — [DOI](https://doi.org/10.1038/s41467-023-43458-x) (the Xenium breast reconstruction)
* Koch J, *et al.* PLoS One 2021;16(1):e0244513. PMID **33395433** — [DOI](https://doi.org/10.1371/journal.pone.0244513) (doxorubicin 3-D:2-D, reused from the GBM report)
* Jubelin C, *et al.* Front Bioeng Biotechnol 2023;11:1260049. PMID **37869710** — [DOI](https://doi.org/10.3389/fbioe.2023.1260049) (doxorubicin 3-D:2-D, reused from the GBM report)
* Mutahar AZI, *et al.* J Transl Med 2025;23(1):1428. PMID **41422245** — [DOI](https://doi.org/10.1186/s12967-025-07585-x) (paclitaxel 3-D:2-D, breast-specific)
* Bear HD, *et al.* J Clin Oncol 2003;21(22):4165-4174. PMID **14559892** — [DOI](https://doi.org/10.1200/JCO.2003.12.005) (NSABP B-27, AC-alone pCR 13.7%)
* Dredze LM, *et al.* Breast Cancer Res Treat 2022;193(3):597-612. PMID **35451732** — [DOI](https://doi.org/10.1007/s10549-022-06598-0) (AC-paclitaxel real-world pCR 26.0%)
* Mouridsen H, *et al.* J Clin Oncol 2001;19(10):2596-2606. PMID **11352951** — [DOI](https://doi.org/10.1200/JCO.2001.19.10.2596) (tamoxifen first-line ORR 20%, original trial)
* Mouridsen H, *et al.* J Clin Oncol 2003;21(11):2101-2109. PMID **12775735** — [DOI](https://doi.org/10.1200/JCO.2003.04.194) (tamoxifen first-line ORR 21%, 32-month update)
* Yang W, *et al.* Nucleic Acids Res 2013;41:D955-D961 (GDSC). PMID **23180760** (GDSC2 breast BRCA cell-line panel, source of the paclitaxel/cyclophosphamide/tamoxifen IC50s)
* Prior work in this repo: `results/kill_calibration.md` (the GBM module this mirrors), `results/breast_cancer_flagship.md` (the reconstruction and drug panel this reuses), `results/breast_drug_panel_data.json` (the independent pipeline this report's doxorubicin full-engine number is bit-identical against), `tmesim/core/config.py` (the PDConfig/DrugParams state observed at run time, §0).


*Generated by `experiments/kill_calibration_breast.py` + `experiments/kill_calibration_breast_report.py` — one reconstruction, 4 transport solves, seed 0, deterministic.*
