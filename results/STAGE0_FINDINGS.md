# Stage 0 Kill-Test — Findings

**Cancer:** glioblastoma. **Question:** does a mechanistic patient-specific tumor-
microenvironment + drug-transport simulation predict drug response *better than a
cheap non-spatial baseline* (mean intrinsic drug sensitivity)? If not, the idea
is not ready to build on — and that's exactly what Stage 0 exists to find out
before any wet-lab or clinical spend.

## What was built
A complete, runnable pipeline (21 modules, 26 passing tests): reconstruct a 2-D
GBM microenvironment (cells, vasculature sparse in the core, ECM) → solve the
physics (convection-diffusion-reaction drug transport, steady-state oxygen +
hypoxia, interstitial fluid pressure with the Jain outward-flow effect) → kill
tumor cells via a Hill pharmacodynamic model modulated by per-cell sensitivity
and hypoxia/quiescence resistance → score predicted vs observed response against
three baselines. Data: a synthetic cohort with a **known, fair** ground truth
(response driven by BOTH transport-limitation and intrinsic sensitivity + noise),
plus a real open-data path (TCGA-GBM, Ivy GAP, GDSC glioma — 131 records).

## The result: the sim does NOT beat the baseline
Two experiments on the same fixed 40-patient synthetic cohort, sweeping dose from
saturating to transport-limited (IC50 = 25 µM).

**Saturating regime (E/IC50 = 10–370):**

| plasma_peak | mean kill | sim AUC | intrinsic-baseline AUC | Δ |
|---|---|---|---|---|
| 4–140 µM | ~0.52 (flat) | 0.883–0.887 | 0.883 | ≈ +0.005 (noise) |

**Transport-limited regime (E/IC50 = 0.05–1):**

| plasma_peak | E/IC50 | mean kill | sim AUC | intrinsic-baseline AUC | Δ |
|---|---|---|---|---|---|
| 0.10 | 0.05 | 0.15 | 0.728 | 0.883 | **−0.155** |
| 0.30 | 0.15 | 0.34 | 0.825 | 0.883 | −0.058 |
| 0.60 | 0.30 | 0.44 | 0.850 | 0.883 | −0.033 |
| 1.00 | 0.50 | 0.48 | 0.860 | 0.883 | −0.023 |
| 2.00 | 1.0 | 0.51 | 0.885 | 0.883 | +0.002 |

Baselines: intrinsic_sensitivity AUC 0.883, expression_logreg 0.850, mean (null) 0.500.

**Across every dose regime, the mechanistic transport simulation ties or trails
the cheap non-spatial intrinsic-sensitivity baseline.** The hypothesis that
transport modeling adds value in the transport-limited regime failed in the
*honest* direction: at low dose the sim gets **worse** (Δ = −0.155), because
delivery heterogeneity adds variance while the dominant signal (intrinsic
sensitivity) is exactly what the baseline already captures.

## Why (mechanistically, not a bug)
Dose propagates correctly — tissue AUC scales linearly with dose (verified). The
flat kill at high dose is genuine Hill **saturation**: above ~10× IC50, every
reached cell dies and kill is set by the resistant/hypoxic fraction, not delivery.
The transport signal only has room to matter near IC50 — and there it did not
help, because in this ground truth the intrinsic-sensitivity term dominates the
transport term.

## The honest caveat (this is not evidence about biology)
This is a **machinery** result, not a biological one. The synthetic ground truth
was constructed to depend strongly on intrinsic sensitivity (corr ≈ 0.77), so a
sensitivity baseline was always going to be hard to beat. A synthetic test can
only prove the pipeline runs, is stable, and that the honest harness will *report
a null when there is one* — which it did. It cannot tell you whether a real tumor
microenvironment carries transport-driven response signal beyond sensitivity.

**That question requires real paired data**: spatial tumor maps (HTAN / Visium /
Xenium GBM) reconstructed through the `reconstruct/` module, paired with measured
drug response (patient-derived organoid / PDX screens, PDMR/HCMI). The open data
we pulled has no spatial ground truth and only weak response labels, so it is an
integration smoke test, not validation. See `docs/data_acquisition_plan.md`.

## What this means for the idea
The most valuable thing Stage 0 produced is the **discipline**: a falsifiable
harness that a mechanistic-modeling story must clear before it earns wet-lab money
— and on the first honest test, the fancy physics did not clear it. That is the
whole thesis in miniature: the scarce, defensible asset is not the simulator, it
is the **verification that tells you when to believe it**. Before building the
transport engine into a product, the next step is to run this same harness on one
real paired dataset (organoid drug screen with spatial data) and see whether the
signal survives contact with reality.
