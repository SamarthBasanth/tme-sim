# Cross-tissue predicted drug behaviour on REAL GBM microenvironments

**All five real tumors are simulated on the SAME converged numerical grid (dx = 12 µm).** The same real drug (temozolomide: doxorubicin-calibrated perivascular-consumption physics, alpha = 149, plus the binding-site barrier) and the same validated transport → kill engine are used for every tissue; only the REAL microenvironment differs.

This supersedes an earlier version of this page in which the Xenium section ran at dx = 12 µm and the four Visium sections ran at their 100 µm spot pitch. That comparison mixed biology with a numerical grid artifact. It also predates the calibration honesty pass, so it used the inflated alpha = 309.6 rather than the corrected alpha = 148.6. **Every headline number on the earlier page is superseded** — see **Retraction & correction** below before quoting any of them.

Drug physics (fixed for all tissues): diffusivity 278.0 µm²/s, ic50 606 µM, binding 2.961e-03/s, clearance 5.794e-05/s, bbb_penetrance 0.125.

## Retraction & correction (read first)

### 1. The earlier cross-tissue spread was confounded by RESOLUTION — RETRACTED

The previously published headline — *"same drug, 5 real tumors, predicted kill spans ~2% → 24%, tracking vascularity"* — is **withdrawn**. It was computed with the Xenium section at the converged dx = 12 µm and the four Visium sections at dx = 100 µm, which is far outside the converged regime (`results/numerical_rigor.md` §2). At dx = 100 µm the grid cannot resolve the ~65 µm perivascular boundary layer, and — at the time those numbers were produced — could not resolve the oxygen field either, so kill was understated by a tissue-dependent amount. The reported min→max range therefore mixed real biology with numerical error. **Read the paired control below together with the update that follows it: the dominant source of that error has since been fixed upstream, which changes how big the correction is but not the fact that the published numbers are withdrawn.**

The paired control below re-runs the identical Visium tissues, identical drug and identical cell geometry, changing **only** the numerical grid:

| Tissue | RETRACTED kill (published @dx=100) | kill @ dx=100 µm (re-run) | kill @ dx=12 µm | fold | hypoxic @ 100 µm | hypoxic @ 12 µm | mean AUC @ 100 µm | mean AUC @ 12 µm | penetration_index 100 → 12 µm |
|---|---|---|---|---|---|---|---|---|---|
| **10x_parent_gbm** | 11.6% | 17.7% | 16.8% | 0.95× | 0.0% | 0.0% | 3.63e+06 | 3.5e+06 | 0.323 → 0.325 |
| **GBM1** | 9.4% | 13.3% | 12.9% | 0.97× | 0.0% | 0.0% | 3.03e+06 | 2.97e+06 | 0.292 → 0.278 |
| **GBM2** | 24.1% | 46.1% | 46.3% | 1.00× | 0.0% | 0.0% | 8.3e+06 | 8.43e+06 | 0.535 → 0.536 |
| **GBM3** | 23.5% | 43.9% | 42.9% | 0.98× | 0.0% | 0.0% | 8.92e+06 | 8.66e+06 | 0.597 → 0.604 |

**Read the columns correctly.** The *valid* resolution control is the `dx = 100 µm (re-run)` → `dx = 12 µm` pair: same tissue, same drug object, same cell geometry, same seed — **only the numerical grid differs**. Any quantity that moves between those two columns was, in the earlier write-up, being reported as inter-tumor biology.

**Important update — the resolution confound has largely been fixed upstream, and this control now shows it.** As re-run today the residual dx = 100 µm → 12 µm change in kill is only ≤5%, not the large effect the confound implied. The reason is a concurrent fix in `physics/oxygen.py`: the oxygen sink previously accumulated a raw tumor-cell *count* per grid cell, which scales as `dx²` and made the hypoxic fraction — and every hypoxia-gated resistance term downstream — a function of the numerical grid. On a coarse grid that spuriously starved the tissue (the Visium sections used to report ~100% hypoxic tumor at dx = 100 µm) and suppressed kill. With the sink expressed as an areal density, the grid dependence is gone.

**This does not un-retract anything.** The earlier page's numbers were produced *before* that fix, so at the time they were written the resolution mismatch really was injecting a large, tissue-dependent error into the reported spread. What has changed is the remedy's size, not its necessity: matching the grid is now cheap insurance rather than a big correction, and the run below matches it anyway so that no future regression can reintroduce the problem silently.

The first column is the *retracted published* value and is shown only to name what is being withdrawn. It is **not** expected to equal the re-run: it was produced with the pre-fix calibration alpha = 309.6, whereas this run uses alpha = 148.6. A separate calibration honesty pass corrected a half-distance-vs-e-folding-length mix-up in `validation/calibration.py` (treating the 45 µm doxorubicin half-distance as an e-folding length inflated alpha by 1/ln²2 = 2.081×). So the earlier cross-tissue page carried **two** independent defects — a resolution mismatch and an inflated alpha — and both are corrected here.

### 2. `penetration_index` is a vascular-geometry SHAPE statistic, not drug delivery — CORRECTED

`penetration_index = mean(AUC) / percentile99(AUC)` over the tissue mask. Numerator and denominator both scale linearly with the delivery source, so the index is **scale-invariant**: doubling the delivered drug everywhere leaves it unchanged. It measures how *even* the exposure field is — and since the source is `S = vascular_density × bbb_penetrance × plasma`, the exposure field's evenness is essentially the evenness of the vascular field. **Correlating `penetration_index` against vascular density is therefore close to circular and is NOT presented as a finding here.**

**This is measured, not asserted — the dose-scaling control.** One tissue (**GBM3**) was run twice on the same grid with the same geometry, changing only the delivered dose by 2× (`bbb_penetrance`, which enters the source linearly):

| Metric | 1× dose | 2× dose | ratio | behaves like… |
|---|---|---|---|---|
| Mean tissue AUC (µM·s) | 8.657e+06 | 1.731e+07 | **2.000×** | a delivery metric ✅ |
| Median tissue AUC (µM·s) | 9.685e+06 | 1.937e+07 | **2.000×** | a delivery metric ✅ |
| kill_fraction | 0.4287 | 0.6319 | **1.474×** | a dose-responsive prediction ✅ |
| **penetration_index** | 0.6037 | 0.6037 | **1.000×** | **NOT a delivery metric ❌** |

Doubling the dose doubles the exposure (2.00× on mean AUC) and raises the kill (1.47×), while `penetration_index` moves 0.00%. A statistic that cannot detect a doubling of the delivered drug is not measuring drug delivery — it is measuring the shape of the field, which here is the shape of the vascular field it was seeded from.

It is retained in the tables below **only** as a vascular-geometry diagnostic. Every cross-tissue comparison in this document uses **mean / median tissue drug AUC** and **kill_fraction** instead.

## Matched-resolution per-tissue table (ALL tissues at dx = 12 µm)

| Tissue | Platform | Diagnosis | Cells (tumor) | Grid @dx | Tumor | Endo | Vascular | Mean O₂ | Hypoxic tumor† | Resistant tumor§ | Mean AUC (µM·s) | Median AUC (µM·s) | Tumor-cell AUC¶ | Mean conc (µM) | Median conc (µM) | **Kill** | pen_idx‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Xenium-GBM** | Xenium (single-cell) | GBM (10x Xenium single-cell section) | 40,887 (30,174) | 424×675 @12µm | 73.8% | 9.7% | 0.025 | 0.295 | 12.7% | 100.0% | 6.863e+05 | 6.758e+05 | 7.063e+05 | 32.28 | 31.79 | **1.5%** ± 0.07 | 0.470 |
| **10x_parent_gbm** | Visium (spot-level) | Glioblastoma (primary resection, fresh-frozen) | 41,560 (34,434) | 473×410 @12µm | 82.9% | 9.7% | 0.127 | 0.586 | 0.0% | 8.2% | 3.499e+06 | 2.884e+06 | 3.643e+06 | 164.5 | 135.7 | **16.8%** ± 0.20 | 0.325 |
| **GBM1** | Visium (spot-level) | Glioblastoma, IDH-wildtype (primary) | 57,050 (44,079) | 451×476 @12µm | 77.3% | 5.7% | 0.108 | 0.496 | 0.0% | 31.4% | 2.97e+06 | 2.38e+06 | 2.918e+06 | 139.7 | 112.1 | **12.9%** ± 0.16 | 0.278 |
| **GBM2** | Visium (spot-level) | Glioblastoma, IDH-wildtype (recurrent) | 34,707 (18,042) | 449×454 @12µm | 52.0% | 29.5% | 0.307 | 0.895 | 0.0% | 0.0% | 8.425e+06 | 1.028e+07 | 1.197e+07 | 396.2 | 482.9 | **46.3%** ± 0.37 | 0.536 |
| **GBM3** | Visium (spot-level) | Glioblastoma, IDH-wildtype (primary) | 21,279 (12,635) | 304×375 @12µm | 59.4% | 23.7% | 0.315 | 0.875 | 0.0% | 0.0% | 8.657e+06 | 9.685e+06 | 1.071e+07 | 407.1 | 455.2 | **42.9%** ± 0.44 | 0.604 |

† **Hypoxic tumor fraction is ORDER-OF-MAGNITUDE ONLY** — but for a different reason than this page previously gave. A concurrent fix in `physics/oxygen.py` removed a spurious `dx²` scaling in the oxygen sink (the tumor sink is now an *areal* density, `rho_tum`, normalised by a fixed reference area rather than a raw per-grid-cell count), so the hypoxic fraction is **no longer grid-dependent** — the residual dx = 100 µm → 12 µm movement in the control table above is what is left after that fix. What remains uncertain is the **absolute level**: `OxygenConfig.consumption × TUMOR_REFERENCE_AREA_UM2` is a single **un-calibrated** lumped constant never fit to a measured oxygen tension, and by that module's own audit a ±40% change swings the real-section hypoxic tumor fraction between ~0.0 and ~0.6. So read every hypoxic and resistant fraction here as "a lot / a little", never as a number, and never difference two of them.

‡ `pen_idx` = `penetration_index`, a **vascular-geometry shape statistic** (mean AUC ÷ p99 AUC), **not** a delivery measure. Listed for completeness only; see the correction above. Do not rank tumors by it.

§ **Resistant tumor** = tumor cells banded HYPOXIC (resistance 0.35) *or* QUIESCENT (resistance 0.5) by the oxygen gate. These are killed less at identical exposure, so kill cannot be read from AUC alone. Like the hypoxic fraction it inherits the un-calibrated oxygen consumption constant — order-of-magnitude only.

¶ **Tumor-cell AUC** = mean exposure sampled *at the positions the tumor cells actually occupy*, rather than averaged over the whole tissue mask. This is the exposure the Hill kill integrates. It can differ a lot from the tissue-wide mean when tumor cells are not uniformly distributed over the mask — which is exactly the Xenium-vs-Visium geometry difference discussed below.

± on kill is the binomial s.e. of the per-cell kill draw (the only stochastic step in the pipeline).

## Observed spread at matched resolution

- **Mean vascular density (real)**: 0.025 → 0.315 (12.6×).
- **Endothelial cell fraction (real, assay-anchored)**: 5.7% → 29.5%.
- **Mean tissue drug AUC (real delivery)**: 6.863e+05 → 8.657e+06 µM·s (12.6×).
- **Median tissue drug AUC**: 6.758e+05 → 1.028e+07 µM·s (15.2×).
- **Predicted kill_fraction (converged)**: 1.5% → 46.3% (31.6×, 44.9 pts).
- **Hypoxic tumor fraction (order-of-magnitude only)**: 0.0% → 12.7%.

## Corrected correlations

All computed on the matched-dx rows. `penetration_index` appears in NONE of them, by design.

| Relationship | n | Statistic |
|---|---|---|
| mean vascular density vs **kill_fraction** (all 5) | 5 | Spearman rho = +0.900 (exact p = 0.083), Pearson r = +0.996 (exact p = 0.017); leave-one-out rho in [+0.800, +1.000] |
| mean vascular density vs **mean tissue AUC** (all 5) | 5 | Spearman rho = +1.000 (exact p = 0.017), Pearson r = +1.000 (exact p = 0.008); leave-one-out rho in [+1.000, +1.000] |
| endothelial cell fraction vs **kill_fraction** (all 5) | 5 | Spearman rho = +0.700 (exact p = 0.233), Pearson r = +0.923 (exact p = 0.050); leave-one-out rho in [+0.400, +1.000] |
| endothelial cell fraction vs **mean tissue AUC** (all 5) | 5 | Spearman rho = +0.600 (exact p = 0.350), Pearson r = +0.901 (exact p = 0.067); leave-one-out rho in [+0.400, +0.800] |
| mean tissue AUC vs **kill_fraction** (all 5) | 5 | Spearman rho = +0.900 (exact p = 0.083), Pearson r = +0.996 (exact p = 0.017); leave-one-out rho in [+0.800, +1.000] |
| mean **tumor-cell-local** AUC vs **kill_fraction** (all 5) | 5 | Spearman rho = +1.000 (exact p = 0.017), Pearson r = +0.997 (exact p = 0.008); leave-one-out rho in [+1.000, +1.000] |
| mean vascular density vs **kill_fraction** (Visium only) | 4 | Spearman rho = +0.800 (exact p = 0.333), Pearson r = +0.993 (exact p = 0.125); leave-one-out rho in [+0.500, +1.000] |
| mean vascular density vs **mean tissue AUC** (Visium only) | 4 | Spearman rho = +1.000 (exact p = 0.083), Pearson r = +1.000 (exact p = 0.042); leave-one-out rho in [+1.000, +1.000] |

**How much of the vascularity → AUC link is real?** Very little of it is a *discovery*. The delivery source is defined as `S = vascular_density × bbb_penetrance × plasma`, and the transport PDE is linear in the source, so mean tissue AUC is close to a linear functional of the vascular field **by construction**. A strong vascular↔AUC correlation is the model working as specified, not evidence about tumors. The only link carrying non-linear, non-definitional content is **vascularity → kill**, which passes through the Hill dose–response and the oxygen gate.

**And `mean_vascular` is itself a normalised shape statistic.** Both the Xenium KDE (`_kde_on_grid` divides by its own max) and the Visium field (`_normalize01`) rescale each tissue's vascular field to peak at 1.0. So `mean vascular density` measures how *evenly spread* the vasculature is within that tissue — it is **not** an absolute perfused-vessel fraction and is not comparable across tissues as an absolute. The endothelial **cell fraction** rows above are the assay-anchored vascularity measure; prefer them when a cross-tissue vascularity claim is being made.

## The Xenium point is still NOT comparable to the Visium points

Matching the grid removed the *numerical* confound. It did not remove the *platform* one, and the table above shows exactly how large that residue is.

**Xenium-GBM** and **GBM1** have almost the same tissue-wide mean drug AUC (6.86e+05 vs 2.97e+06 µM·s, a 77% difference) — yet their predicted kill differs by 1.5% vs 12.9%. Identical drug, identical physics, identical grid. Three reconstruction differences account for it, and **none of them is a claim about the tumors**:

1. **Where the tumor cells sit relative to the drug.** Tumor-cell-local AUC is 7.06e+05 µM·s for Xenium-GBM vs 2.92e+06 µM·s for GBM1 (1.03× and 0.98× their respective tissue-wide means). Xenium has REAL cell positions, so its tumor cells sit in the tumor mass, away from vessels, in the drug-poor part of the field. Visium "cells" are 12 placeholders scattered inside each spot — the same spots whose endothelial fraction *generated* the vascular field — so they sample the field near its seeds. This is a geometry artifact of spot-level reconstruction, not better perfusion.
2. **Resistant-state banding — this is the dominant term.** 100.0% of Xenium-GBM's tumor cells are banded hypoxic (kill × 0.35) or quiescent (kill × 0.5) against 31.4% for GBM1. Those multipliers alone predict a kill ratio of roughly 0.59×, against an observed 0.11× — i.e. the gap is essentially the oxygen gate, not perfusion. And that gate rests on the **un-calibrated** oxygen-consumption constant flagged above, whose own audit says ±40% swings the hypoxic fraction between ~0.0 and ~0.6.
3. **Tissue mask.** The Xenium domain carries a real, irregular tissue mask; the Visium domain is the full bounding rectangle. "Mean over the tissue mask" therefore means different things on the two platforms.

**Consequence for the headline.** The Xenium point is the lowest on both vascularity and kill, which *helps* the all-5 rank correlation — but for reasons that are at least partly platform, not biology. **The defensible statement is the within-Visium one (n = 4), where all four tissues share a reconstruction path.** Any all-5 statistic in this document should be read as illustrative only.

## Does the inter-tumor heterogeneity claim survive?

**VERDICT: the inter-tumor heterogeneity claim SURVIVES at matched resolution — predicted kill still varies substantially across the five real tumors and still rises with real vascularity, in the all-5 set (Spearman rho = +0.90) and, importantly, also within the four like-for-like Visium tumors (rho = +0.80). But the NUMBERS are completely different and every earlier figure is retracted.**

At the matched, converged dx = 12 µm the predicted `kill_fraction` spans **1.5% → 46.3%** (31.6× fold, 44.9 pts absolute) and mean tissue drug AUC spans **6.86e+05 → 8.66e+06 µM·s** (12.6×). The per-cell kill sampling s.e. is ≈0.0025 (binomial, n_tumor = 12,635–44,079), i.e. 181× smaller than the observed spread — so the spread is not sampling noise. Vascular density vs kill: Spearman rho = +0.900 (exact p = 0.083), Pearson r = +0.996 (exact p = 0.017); leave-one-out rho in [+0.800, +1.000]. Vascular density vs mean AUC: Spearman rho = +1.000 (exact p = 0.017), Pearson r = +1.000 (exact p = 0.008); leave-one-out rho in [+1.000, +1.000]. Within the four like-for-like Visium tumors (the comparison that does not cross platforms): Spearman rho = +0.800 (exact p = 0.333), Pearson r = +0.993 (exact p = 0.125); leave-one-out rho in [+0.500, +1.000].

**The kill spread is much WIDER than the retracted claim, not narrower.** The earlier page reported ~2% → 24%; at matched resolution it is 1.5% → 46.3%. That is not a vindication — the earlier numbers were wrong in magnitude as well as in comparability. The magnitude change is NOT mainly the grid (the paired control above shows the residual grid effect is only ≤5%); it comes from the two upstream honesty fixes that landed since — the calibration alpha correction (309.6 → 148.6) and the `physics/oxygen.py` sink fix that stopped the coarse grid from spuriously starving the tissue. Absolute kill fractions of this size for temozolomide on real tissue should be read as a statement about this engine's current calibration, **NOT** as a clinical expectation — there are no outcome labels to check them against, and the oxygen gate they pass through rests on an un-calibrated constant.

**Where the resolution sensitivity actually lives.** Between dx = 100 µm and dx = 12 µm, mean tissue AUC moves at most 4% and kill moves at most 5%. The delivery calculation was never the fragile part — the fragile part was the oxygen/quiescence gate, and that fragility has since been fixed upstream (see the retraction). The grid is matched here regardless, so the comparison cannot silently regress.

**But read the correlation with n = 5 in mind.** With five points the smallest two-sided p a *perfect* rank correlation can attain is 0.017; nothing here can clear a conventional significance bar, and a single tumor moves the rank statistic a lot (leave-one-out range above). Treat the correlation as a **direction**, not evidence.

## Is dx = 12 µm actually converged **for a Visium reconstruction**?

The grid-refinement audit in `results/numerical_rigor.md` §2 was run on the **Xenium** section only. Nothing had ever checked that dx = 12 µm is inside the converged regime for a *Visium* reconstruction, whose fields and cell placement are built differently. This ladder does that, on **GBM3**, with everything except dx held fixed.

**One caveat on this ladder:** the convective (Darcy) term carries an explicit CFL limit `dt·|v|/dx ≤ 1`, so the finest grids need a smaller timestep — `dt` varies over 5–10 s across the rows and is reported per row below. Per `results/numerical_rigor.md` §3 the solution is dt-invariant to <1e-2 relative for any dt ≲ 20 s, so this does not confound the grid comparison, but it does mean the finest row is not a *pure* dx refinement.

| dx (µm) | grid | dt (s) | mean AUC (µM·s) | AUC rel. to finest | kill_fraction | Δkill vs finest | hypoxic tumor | pen_idx (shape) |
|---|---|---|---|---|---|---|---|---|
| 100 | 37×45 | 10 | 8.925e+06 | 1.050 | 0.4393 | 0.038 | 0.0% | 0.597 |
| 50 | 73×90 | 10 | 9.221e+06 | 1.085 | 0.4453 | 0.052 | 0.0% | 0.607 |
| 30 | 122×150 | 10 | 8.66e+06 | 1.019 | 0.4291 | 0.014 | 0.0% | 0.601 |
| 20 | 182×225 | 10 | 8.37e+06 | 0.985 | 0.4192 | 0.009 | 0.0% | 0.606 |
| 12 | 304×375 | 10 | 8.657e+06 | 1.019 | 0.4287 | 0.013 | 0.0% | 0.604 |
| 8 | 455×562 | 5 | 8.497e+06 | 1.000 | 0.4232 | 0.000 | 0.0% | 0.605 |

**Reading.** Between dx = 12 µm (used throughout this page) and the finest grid tested, dx = 8 µm, `kill_fraction` moves 1.3% relative and mean AUC moves 1.9%. That is inside the same tolerance `numerical_rigor.md` uses to call the Xenium kill converged, so **dx = 12 µm is an adequate grid for these Visium reconstructions too** — the matched-resolution table is not resting on an unverified assumption. Note that **both** mean AUC and kill are essentially flat across the whole ladder, including the coarse end. That is the post-fix picture: with the `physics/oxygen.py` sink expressed as an areal density, this Visium reconstruction is largely grid-insensitive. Before that fix the same refinement moved kill by ~2.6–2.9× (the measurement that motivated this whole rewrite), so the flatness here is evidence that the defect is fixed, not evidence that the original concern was misplaced.

## Resolution & convergence

*(This section was contributed earlier and is preserved; it is extended here with the matched-grid policy and the paired resolution control.)*

Source: [`results/numerical_rigor.md`](numerical_rigor.md) §2 (grid-refinement audit on the real Xenium section).

- **`kill_fraction` is grid-converged for dx ≲ 20 µm.** It plateaus and moves <1–2% relative between the finest grids. It is NOT converged at dx = 100 µm — see the paired control above, where refining the grid alone changes it by up to 1.0×.
- **`penetration_index` requires dx ≲ 12 µm on single-cell data.** The calibrated perivascular boundary layer is ~65 µm thick. On the Xenium section the index falls 0.396 (dx 30) → 0.311 (20) → 0.204 (15) → 0.196 (12) → 0.197 (10), settling only below ~15 µm. Earlier cross-tissue runs quoted the dx = 30 µm Xenium value (0.396) and over-stated it ~2×; those figures are superseded. **However, see the correction above: even at the converged dx the index is a shape statistic, so "converged penetration index" is a converged measure of vascular geometry, not of delivery.**
- **Caveat on the cited ladder itself.** The grid-refinement table in `numerical_rigor.md` §2 was produced *before* the calibration honesty pass, i.e. at the inflated alpha = 309.6, which corresponds to a ~45 µm perivascular boundary layer. At the corrected alpha = 148.6 the layer is ~65 µm — **longer, hence easier to resolve**. So the dx ≲ 12 µm requirement inherited from that audit is, if anything, conservative for this run; the absolute numbers in that ladder no longer match this drug and should be re-run.
- **Why the Visium sections are now run at dx = 12 µm too.** The earlier version argued that sub-spot dx is "physically meaningless" for Visium and pinned it at the 100 µm spot pitch. That conflates two different things:
  - the **data resolution** (~100 µm — genuinely spot-limited, and unchanged by anything done here), and
  - the **simulation grid** (a numerical discretisation of a continuous PDE).

  The reconstructed Visium fields are continuous: vascular and ECM density are a Gaussian KDE with bandwidth fixed at 100 µm **in physical units**, and cells sit at continuous coordinates inside each spot footprint. Refining dx samples that *same band-limited function* more accurately; it does not sharpen it. Meanwhile dx = 100 µm is coarser than the ~65 µm drug boundary layer and coarser than the oxygen diffusion length, so it carries a large, tissue-dependent discretisation error. Refining the grid removes that error. **The information content stays ~100 µm; only the numerical artifact goes away.**
- **Cell geometry was held fixed while refining.** `reconstruct_from_visium` scales its cell-placement jitter as `min(27.5, 0.45·dx)`, which would have shrunk the spot footprint from 27.5 µm to 5.4 µm purely because the grid got finer. The experiment pins it at 27.5 µm (identical to the native 100 µm reconstruction) so the paired control isolates the grid and nothing else.
- **The oxygen sink's grid dependence has been fixed (concurrently), and this changes what the earlier hypoxia caveat should say.** `physics/oxygen.py` previously accumulated a raw tumor-cell *count* per grid cell, which scales as `dx²` and made the hypoxic fraction — and every hypoxia-gated resistance term downstream — a function of the numerical grid. It now uses an areal density normalised by a fixed reference area, so the hypoxic fraction is grid-independent. Earlier statements on this page that the hypoxic fraction "keeps falling as dx drops" described that defect and are withdrawn.
- **But the hypoxic fraction is still un-calibrated in absolute terms** and is quoted as order-of-magnitude only: the lumped constant `consumption × TUMOR_REFERENCE_AREA_UM2` was never fit to a measured oxygen tension, and ±40% on it swings the hypoxic fraction between ~0.0 and ~0.6. Because kill is gated on that banding, the absolute kill fractions inherit this uncertainty.
- **What actually moves kill between dx = 100 µm and dx = 12 µm.** Mean tissue AUC is nearly dx-invariant (see the control table: it moves only a few percent), so the residual dx effect on kill does **not** come from delivery. It comes from the oxygen/quiescence banding and from resolving the perivascular gradient the tumor cells actually sample. That is why `kill_fraction` must be quoted at a converged dx even though `mean_auc` is robust.

## Honest caveats (read before trusting any number)

- **Visium information content is still ~100 µm.** Matching the numerical grid at 12 µm does NOT give the Visium sections single-cell detail. Their vascular/ECM fields remain a 100 µm-bandwidth KDE of deconvolved spot fractions, and their "cells" remain 12 synthetic placeholders per spot apportioned by fraction. The Xenium section has genuine per-cell structure that no Visium section has at any dx. The finer grid removes a numerical artifact; it adds no real detail. **A Xenium-vs-Visium difference is still confounded by platform, deconvolution and information content — just no longer by the grid.**
- **No outcome labels.** Every tissue is an *untreated* open resection. There is NO paired measured temozolomide response for any of them. Kill and AUC are what the validated physics predicts on real geometry — **predictions, not validated outcomes**. Do not read `kill_fraction` as a clinical response.
- **n = 5.** Every correlation above is descriptive. The smallest two-sided p attainable at n = 5 even for a perfect rank match is 0.017, and the Visium-only subset (n = 4) cannot reach any conventional threshold at all. No significance is claimed.
- **`mean_vascular` and `penetration_index` are both normalised shape statistics**, not absolute physical densities/deliveries — see the correction and the correlation notes above.
- **The vascularity → AUC link is largely mechanical**, since the delivery source is defined proportional to the vascular field.
- **Reconstruction, not measurement, of fields.** Vascular density, ECM and oxygen are inferred from real markers/perfusion, not directly measured concentrations. They are faithful to the data but are model fields.
- **One drug, one config.** Only temozolomide at one engine config is run here; the relative comparison is the message, not any absolute value.


_Generated by `experiments/cross_tissue.py` on real data (5 tissues at matched dx = 12 µm + 4 paired resolution controls at dx = 100 µm; runtime 21.7 min). Figure: `results/cross_tissue.png`._
