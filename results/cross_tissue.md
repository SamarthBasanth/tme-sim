# Cross-tissue predicted drug behaviour on REAL GBM microenvironments

The **same** real drug (temozolomide: doxorubicin-calibrated perivascular-consumption physics, alpha = 310, plus the binding-site barrier) and the **same** validated transport -> kill engine were run across every real glioblastoma tissue below. The drug and the physics are identical for all tumors; only the REAL microenvironment differs. Any spread in predicted penetration or kill is therefore **real inter-tumor heterogeneity**, not a tuned knob.

Drug physics (fixed for all tissues): diffusivity 278.0 um^2/s, ic50 606 uM, binding 6.169e-03/s, clearance 5.794e-05/s, bbb_penetrance 0.558.

## Per-tissue table

| Tissue | Platform | Diagnosis | Cells (tumor) | Grid @dx | Tumor | Endo | Vascular | Oxygen | Hypoxic tumor | Pred. penetration | Pred. kill |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Xenium-GBM** | Xenium (single-cell) | GBM (10x Xenium single-cell section) | 40,887 (30,174) | 170x271 @30um | 73.8% | 9.7% | 0.094 | 0.156 | 93.7% | 0.396 | 1.9% |
| **10x_parent_gbm** | Visium (spot-level) | Glioblastoma (primary resection, fresh-frozen) | 41,560 (34,434) | 57x50 @100um | 82.9% | 9.7% | 0.132 | 0.035 | 100.0% | 0.283 | 11.6% |
| **GBM1** | Visium (spot-level) | Glioblastoma, IDH-wildtype (primary) | 57,050 (44,079) | 55x58 @100um | 77.3% | 5.7% | 0.110 | 0.022 | 100.0% | 0.264 | 9.4% |
| **GBM2** | Visium (spot-level) | Glioblastoma, IDH-wildtype (recurrent) | 34,707 (18,042) | 54x55 @100um | 52.0% | 29.5% | 0.302 | 0.124 | 99.9% | 0.509 | 24.1% |
| **GBM3** | Visium (spot-level) | Glioblastoma, IDH-wildtype (primary) | 21,279 (12,635) | 37x45 @100um | 59.4% | 23.7% | 0.325 | 0.118 | 99.4% | 0.566 | 23.5% |

## Observed real spread across tumors

- **Vascular density**: 0.094 -> 0.325 (spread 0.231, 3.5x).
- **Mean oxygen**: 0.022 -> 0.156 (spread 0.133).
- **Hypoxic tumor fraction**: 93.7% -> 100.0% (spread 6.3 pts).
- **Predicted penetration_index**: 0.264 -> 0.566 (spread 0.301).
- **Predicted kill_fraction**: 1.9% -> 24.1% (spread 22.2 pts).

## Interpretation: how real microenvironment differences drive the predictions

The engine is a convection-diffusion-reaction transport solve followed by an oxygen-modulated Hill kill. Two real-tissue levers move the predictions:

1. **Vascularity -> delivery -> penetration.** The drug source is seeded at the real endothelial (vascular) field, so tumors with denser, more evenly distributed real vasculature get drug delivered to more of the tissue and show a higher predicted penetration_index. Panel (a) plots this real vascularity against predicted penetration across the tumors.

2. **Hypoxia -> oxygen-gated kill.** Hypoxic tumor regions (low real oxygen, inferred from the real-perfusion oxygen solve) are quiescent and less drug-sensitive, so tumors carrying a larger real hypoxic burden show a lower predicted kill for the identical exposure. Panel (b) plots real hypoxic fraction against predicted kill.

Because composition (tumor/stroma/immune/endothelial fractions), vascular density, ECM and oxygen all differ between these REAL resections, the same drug physics lands differently on each -- which is exactly the inter-tumor heterogeneity a per-patient simulator should surface.


**Which lever dominates in THIS cohort (honest read).** Hypoxic tumor fraction is near-saturated across all tissues (93.7%-100.0%), so it barely varies and does NOT discriminate between these tumors -- panel (b) is nearly flat on the x-axis. The real discriminator here is **vascularity -> penetration**: predicted penetration_index spans 0.264-0.566 and tracks the 0.094-0.325 spread in real vascular density, and predicted kill follows penetration far more than hypoxia. Read panel (a), not (b), as the heterogeneity story for this cohort; the hypoxia axis is kept for honesty, not because it separates these particular tumors.

## Honest caveats (read before trusting any number)

- **No outcome labels.** Every tissue here is an *untreated* open resection. There is NO paired measured temozolomide response for any of them. The kill and penetration numbers are what the validated physics predicts on real geometry -- they are **predictions, not validated outcomes**. Do not read the kill fraction as a clinical response.
- **Platform difference (Xenium vs Visium).** The Xenium section is single-cell resolution: real individual cells, per-cell threshold-gated marker typing, and an expression-weighted ECM field from real matrix transcripts. The four Visium sections are spot-level (~100 um spots): cell-type fractions are **deconvolved** (reference-based, coarser), each spot is a small neighbourhood rather than a cell, and the ECM field is a **stromal-fraction proxy** (no per-cell transcripts). Comparing the single-cell tumor to the spot-level tumors is therefore NOT apples-to-apples: differences between the Xenium point and the Visium points partly reflect platform + resolution + deconvolution, not just biology. In the figure Xenium is red and Visium is blue for exactly this reason.
- **Reconstruction, not measurement, of fields.** Vascular density, ECM and oxygen are inferred from real markers/perfusion, not directly measured concentrations. They are faithful to the data but are model fields.
- **One drug, one config.** Only temozolomide at one engine config is run here; the relative spread is the message, not any absolute value.


_Generated by `experiments/cross_tissue.py` on real data (5 tissues, runtime 24s). Figure: `results/cross_tissue.png`._
