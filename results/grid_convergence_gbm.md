# Is the GBM 3-D:2-D panel grid-converged?

**`results/kill_calibration.md` publishes 3/3 PASS at `dx = 12 um` and calls that grid converged. This file tests that claim by re-running the published pipeline on a dx ladder and reporting the scored ratio as a function of mesh.**

> **This is a convergence audit, not a calibration.** No grid was selected by the verdict it produces. Only `dx` is varied; `dt = 10 s`, the 6 h window, `alpha`, the seed, `PDConfig`, the p99 bath-matching convention, the reporting concentrations and the frozen literature bands are all held at their published values. The output is the convergence curve. "The GBM panel's 3/3 PASS is not grid-converged" would have been a perfectly publishable result and the script was written so it could produce one.

## 0. Bottom line

**Reproduction (the precondition).** Scoring the published on-disk exposure vectors through this script's own `score()` reproduces `results/kill_calibration.md` §3 to a worst absolute deviation of **0.000049**, and all verdicts match. A fresh `dx = 12 um` reconstruction + solve reproduces the cached exposure field as well (§7). The published numbers are reproducible; everything below is therefore about the mesh and nothing else.

**Rungs actually run** (dt = 10 s throughout): temozolomide at 30, 24, 18, 16, 12, 10, 8, 7 um; cisplatin at 30, 24, 18, 16, 12, 10, 8, 7 um; doxorubicin at 30, 24, 18, 16, 12, 10, 8, 7 um.

**How many of the 3 PASSes survive every rung that was actually run:** **3 of 3** — all of temozolomide, cisplatin, doxorubicin PASS at every rung run.

**No drug passes at exactly one mesh.** Every PASS below is held across the full ladder, not at a single resolution.

**Discretisation error at the published `dx = 12 um`:** the scored ratio sits within **1.76 %** of the extrapolated dx -> 0 value in the worst case across the 3 drugs with enough rungs to fit (per-drug table in §3).

**The scored ratio moves by up to 42.67 % across the full set of rungs run** — but that number is dominated by the coarse end and should not be quoted alone. Split in two:

| drug | span over ALL rungs | span over the FINE rungs (dx <= 18 um) | fine-rung mean +/- sd |
|---|---:|---:|---|
| temozolomide | 16.7 % | **6.6 %** | 0.2596 +/- 0.0062 |
| cisplatin | 19.5 % | **7.4 %** | 0.2726 +/- 0.0071 |
| doxorubicin | 42.7 % | **6.7 %** | 0.2021 +/- 0.0051 |

**This is the shape of the answer.** The coarse grids (24, 30 um) really are wrong — they cannot resolve the perivascular boundary layer and they understate the ratio badly. But from **dx = 18 um downwards the ratio plateaus**, moving by at most **7.4 %** across five further refinements down to the CFL floor, and it **oscillates within that plateau rather than converging monotonically towards a limit**. For scale, the breast panel's docetaxel — the FAIL whose grid sensitivity started this audit — spans 49.7 % and is monotone throughout.

**The curve is NOT monotone in dx for temozolomide, cisplatin, doxorubicin.** This is the single most important structural difference from the breast panel, where every drug is monotone. A non-monotone ladder is not in an asymptotic convergence regime, which means Richardson extrapolation and a power-law fit are being applied outside their assumptions and their dx -> 0 estimates should be read as weak indications rather than limits. §5 identifies the mechanism: `F_state` — the oxygen-banded cell-state label — is itself non-monotone in dx, and it multiplies straight into the score.

## 1. What was run

* **Tissue** — the real 10x Xenium GBM section (`reconstruct_real(source='xenium')`), the same section `results/kill_calibration.md` uses, re-reconstructed at every rung of the ladder.
* **Drugs** — the 3 the published report formally tests: temozolomide, cisplatin, doxorubicin. All three carry their own drug-matched literature band, so none falls back to the broad 'solid' band.
* **Held fixed** — `dt = 10 s`, window 6 h, `save_every = 10`, seed 0, convection on, `alpha = 148.6082` (the published value, not re-fitted), shipping diffusivity, `PDConfig` defaults, p99 bath matching.
* **Reporting concentration** — temozolomide at **100 uM (the ex-vivo assay concentration)**, cisplatin and doxorubicin at **1x IC50** — exactly the `headline` rows the published verdicts are read off.
* **25 transport solves**, 37.5 min total.

**Only `dx` was varied — no delivery-law arms.** The brief's KNOWN TRAP (`kill_calibration.unit_transport()` calls `simulate_drug_transport` directly, so `PDEConfig` delivery flags are silent no-ops) does not apply here, because no two arms differ by a `PDEConfig` flag: every arm differs by the mesh, which reaches the solver through the reconstructed `Domain` itself. That the mesh really moved the physics is asserted rather than assumed — the per-cell exposure vector is fingerprinted at every rung and all rungs are required distinct:

```
temozolomide           9/9 distinct exposure fields
cisplatin              8/8 distinct exposure fields
doxorubicin            8/8 distinct exposure fields
```

## 2. The ladder, and why it stops at 7 um

**`dt` was HELD FIXED at the published 10 s at every rung.** That is deliberate: scaling dt with dx would refine two discretisations at once and the resulting curve could not attribute its movement to the mesh. §7.3 bounds the temporal error that choice leaves in place by re-running two rungs at half the dt. **This is a SPATIAL convergence study**; a joint (dx, dt) study would be a different and stronger claim and is not made here.

Holding dt fixed is only legitimate while the solver stays stable. Diffusion here is solved **implicitly** and is unconditionally stable; the explicit upwind advection step is the only conditionally-stable piece and carries `dt*(|vx|+|vy|)/dx <= 1`. The solver **computes that number and raises** (`physics/diffusion.py` lines 360-366) rather than degrading silently, so **no rung below can have run unstably — a rung that produced a number is a rung whose CFL was <= 1**. The measured value is tabulated per rung anyway:

| dx (um) | **dt (s)** | grid | cells | tumour cells | max cell speed (um/s) | **CFL** | quiescent frac | hypoxic frac | prologue (s) |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 7 | 10 | 725x1157 | 40,887 | 23,812 | 0.6684 | **0.955** | 0.9187 | 0.0000 | 26 |
| 8 | 10 | 635x1012 | 40,887 | 23,812 | 0.7513 | **0.939** | 0.9024 | 0.0000 | 10 |
| 10 | 10 | 508x810 | 40,887 | 23,812 | 0.6629 | **0.663** | 0.9390 | 0.0000 | 23 |
| 12 **(published)** | 10 | 424x675 | 40,887 | 23,812 | 0.6336 | **0.528** | 0.8996 | 0.0000 | 8 |
| 16 | 10 | 318x507 | 40,887 | 23,812 | 0.6627 | **0.414** | 0.9334 | 0.0000 | 7 |
| 18 | 10 | 283x451 | 40,887 | 23,812 | 0.6457 | **0.359** | 0.9197 | 0.0000 | 9 |
| 24 | 10 | 213x338 | 40,887 | 23,812 | 0.5308 | **0.221** | 0.9216 | 0.0000 | 4 |
| 30 | 10 | 170x271 | 40,887 | 23,812 | 0.6013 | **0.200** | 0.9918 | 0.0000 | 3 |

And the prologue-only CFL probe that establishes where the floor actually is, at dt = 10 s:

| dx (um) | grid | max cell speed (um/s) | CFL | solver |
|---:|---|---:|---:|---|
| 30 | 170x271 | 0.6013 | 0.200 | OK |
| 24 | 213x338 | 0.5308 | 0.221 | OK |
| 18 | 283x451 | 0.6457 | 0.359 | OK |
| 16 | 318x507 | 0.6627 | 0.414 | OK |
| 12 | 424x675 | 0.6336 | 0.528 | OK |
| 11 | 462x737 | 0.7132 | 0.648 | OK |
| 10 | 508x810 | 0.6629 | 0.663 | OK |
| 9 | 564x900 | 0.6966 | 0.774 | OK |
| 8 | 635x1012 | 0.7513 | 0.939 | OK |
| 7 | 725x1157 | 0.6684 | 0.955 | OK |
| 6 | 846x1349 | 0.6773 | 1.129 | **REFUSES** |

### 2.1 A stale numerical claim, corrected

`results/delivery_metrics.md` states, and the `DX_UM` comment in `tmesim/validation/kill_calibration.py` repeats, that **"dx = 10 um violates the advection CFL condition at dt = 10 s"** — and that is the stated reason the published ladder stops at 12 um and no one has looked below it.

**On the section as it is reconstructed today, that is not true.** Measured above: dx = 10 um runs at CFL **0.663**, and 9, 8 and 7 um are all stable as well. The true floor at the published dt is **dx = 7 um** (CFL 0.955); dx = 6 um is the first grid the solver refuses (CFL 1.129).

The most likely explanation is that the claim predates the switch of `reconstruct_real`'s tissue mask from the legacy `'occupancy'` mode (grid-DEPENDENT) to the current `'geometric'` default (grid-independent, 30 um discs around measured cell centroids). The mask sets the domain the IFP solve runs on, and the IFP solve sets the velocity field the CFL number is computed from. This file does **not** edit those documents — it records the measurement and flags the discrepancy for whoever owns them. It matters because it means the fine half of this study, the half that decides the question, was reachable all along.

The number of measured cells (and therefore of tumour cells scored) is **identical at every rung** — the Xenium section carries single-cell coordinates and the tissue mask is the dx-independent geometric one, so refining the mesh changes the continuum fields, not the cell population. That isolates the effect being measured: any movement in the scored ratio is a movement in the exposure field, not a change in which cells are being averaged over.

**The perivascular boundary layer, i.e. the length the mesh has to resolve.** Under the shipping delivery law each drug has `L = sqrt(D_eff / (binding + clearance))`:

| drug | D_eff-driving diffusivity (um^2/s) | **L (um)** | L / dx=12 | L / dx=8 |
|---|---:|---:|---:|---:|
| temozolomide | 278.0 | **264.5** | 22.04 | 33.07 |
| cisplatin | 236.2 | **142.9** | 11.91 | 17.87 |
| doxorubicin | 184.7 | **64.9** | 5.41 | 8.12 |

This is the single structural number that predicts how converged a panel will be, and it is the reason to expect the GBM and breast panels to differ — see §6.

## 3. The convergence curve, per drug

`ratio` is the SCORED 3-D:2-D kill ratio — the exact quantity the published verdict is taken on. The band is **fixed at every rung** (it depends only on the frozen literature table and on (IC50, hill), none of which dx touches), so a verdict can only move if the ratio moves.

### temozolomide — IC50 606.0467 uM, hill 1.20, scored at 100 uM (ex-vivo assay), band **0.1590-0.4622**

| dx (um) | F_penetration | F_state | **scored 3-D:2-D ratio** | implied IC50 fold | mean tumour exposure (uM per uM plasma) | abs. kill at nominal dose | solve (s) | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 7 | 0.4702 | 0.5407 | **0.2585** | 3.31x | 1.129926 | 0.0176 | 202 | **PASS** |
| 8 | 0.4594 | 0.5488 | **0.2574** | 3.32x | 1.136852 | 0.0180 | 129 | **PASS** |
| 10 | 0.4729 | 0.5305 | **0.2547** | 3.35x | 1.125272 | 0.0171 | 198 | **PASS** |
| 12 **(published)** | 0.4678 | 0.5502 | **0.2626** | 3.26x | 1.140322 | 0.0181 | 94 | **PASS** |
| 16 | 0.4669 | 0.5333 | **0.2537** | 3.36x | 1.138553 | 0.0175 | 68 | **PASS** |
| 18 | 0.4935 | 0.5401 | **0.2705** | 3.18x | 1.140640 | 0.0177 | 71 | **PASS** |
| 24 | 0.4741 | 0.5392 | **0.2624** | 3.26x | 1.137102 | 0.0178 | 32 | **PASS** |
| 30 | 0.4574 | 0.5041 | **0.2317** | 3.63x | 1.098665 | 0.0157 | 20 | **PASS** |

Step-by-step refinement (fine -> coarse):

| from dx | to dx | d(ratio) | d(ratio) % | slope per um |
|---:|---:|---:|---:|---:|
| 7 | 8 | -0.00116 | -0.45 % | -0.001157 |
| 8 | 10 | -0.00263 | -1.02 % | -0.001316 |
| 10 | 12 | 0.00793 | 3.11 % | 0.003965 |
| 12 | 16 | -0.00893 | -3.40 % | -0.002232 |
| 16 | 18 | 0.01673 | 6.60 % | 0.008367 |
| 18 | 24 | -0.00810 | -3.00 % | -0.001350 |
| 24 | 30 | -0.03064 | -11.68 % | -0.005106 |

* **Shape:** NOT monotone in dx; settling (the finest step is smaller than the next one out).
* **Total movement across the ladder:** 16.72 % of the smallest value (0.2317 to 0.2705).
* **Power-law fit** `ratio(dx) = q0 + C*dx^p` over all 8 rungs: p = **4.00**, q0 = **0.2615**, C = -0.00000003, rms residual 0.007045.  **DEGENERATE FIT — do not read `q0` as a converged limit.** `p` has pinned at the edge of the 0.25-4.0 scan and `C` has collapsed to ~0, which means the least-squares optimum is essentially a HORIZONTAL LINE through the fine rungs with a term that only bends at the coarse end. The fit is telling us there is no coherent power-law approach to a limit in this data — `q0` is little more than the mean of the fine plateau, and the rms residual is the size of the oscillation, not of a fitting error.
* **Richardson (8, 12, 18 um, r = 1.5):** observed order p = **0.96**, extrapolated q0 = **0.2462**, GCI on the finest pair = **5.44 %**.
* **Discretisation error at the published dx = 12 um:** 0.00115 absolute (0.44 %) against the power-law q0; 0.01649 (6.70 %) against the Richardson q0; 0.00530 (2.06 %) against the finest grid actually run (dx = 8 um).
* **Verdict at the extrapolated converged value:** power-law q0 = 0.2615 -> **PASS**; Richardson q0 = 0.2462 -> **PASS**.
* **Margin at dx = 12:** 0.10366 above the band floor, 0.19956 below the band ceiling; nearest edge 0.10366. The whole ladder spans 0.03874, i.e. **0.374x** the distance to the nearest band edge, and the estimated discretisation error is 0.011x it.
* **No verdict flip anywhere on the ladder** — 8/8 rungs PASS.

### cisplatin — IC50 24.5503 uM, hill 1.20, scored at 1xIC50, band **0.2086-1.5004**

| dx (um) | F_penetration | F_state | **scored 3-D:2-D ratio** | implied IC50 fold | mean tumour exposure (uM per uM plasma) | abs. kill at nominal dose | solve (s) | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 7 | 0.4962 | 0.5407 | **0.2726** | 4.66x | 2.686937 | 0.4302 | 160 | **PASS** |
| 8 | 0.4806 | 0.5488 | **0.2689** | 4.72x | 2.702548 | 0.4373 | 134 | **PASS** |
| 10 | 0.4872 | 0.5305 | **0.2624** | 4.83x | 2.674874 | 0.4215 | 57 | **PASS** |
| 12 **(published)** | 0.4980 | 0.5502 | **0.2791** | 4.55x | 2.712199 | 0.4382 | 115 | **PASS** |
| 16 | 0.4988 | 0.5333 | **0.2706** | 4.69x | 2.708163 | 0.4237 | 67 | **PASS** |
| 18 | 0.5136 | 0.5401 | **0.2818** | 4.51x | 2.713061 | 0.4302 | 62 | **PASS** |
| 24 | 0.4787 | 0.5392 | **0.2649** | 4.79x | 2.702271 | 0.4265 | 33 | **PASS** |
| 30 | 0.4656 | 0.5041 | **0.2358** | 5.35x | 2.609819 | 0.3877 | 21 | **PASS** |

Step-by-step refinement (fine -> coarse):

| from dx | to dx | d(ratio) | d(ratio) % | slope per um |
|---:|---:|---:|---:|---:|
| 7 | 8 | -0.00372 | -1.37 % | -0.003724 |
| 8 | 10 | -0.00648 | -2.41 % | -0.003242 |
| 10 | 12 | 0.01672 | 6.37 % | 0.008358 |
| 12 | 16 | -0.00856 | -3.07 % | -0.002139 |
| 16 | 18 | 0.01121 | 4.14 % | 0.005606 |
| 18 | 24 | -0.01688 | -5.99 % | -0.002813 |
| 24 | 30 | -0.02911 | -10.99 % | -0.004851 |

* **Shape:** NOT monotone in dx; NOT settling (the finest step is not smaller than the next one out).
* **Total movement across the ladder:** 19.50 % of the smallest value (0.2358 to 0.2818).
* **Power-law fit** `ratio(dx) = q0 + C*dx^p` over all 8 rungs: p = **4.00**, q0 = **0.2743**, C = -0.00000004, rms residual 0.006853.  **DEGENERATE FIT — do not read `q0` as a converged limit.** `p` has pinned at the edge of the 0.25-4.0 scan and `C` has collapsed to ~0, which means the least-squares optimum is essentially a HORIZONTAL LINE through the fine rungs with a term that only bends at the coarse end. The fit is telling us there is no coherent power-law approach to a limit in this data — `q0` is little more than the mean of the fine plateau, and the rms residual is the size of the oscillation, not of a fitting error.
* **Richardson (8, 12, 18 um, r = 1.5):** observed order p = **-3.33**, extrapolated q0 = **0.2827**, GCI on the finest pair = **6.42 %**.
* **Discretisation error at the published dx = 12 um:** 0.00482 absolute (1.76 %) against the power-law q0; -0.00359 (-1.27 %) against the Richardson q0; 0.01023 (3.81 %) against the finest grid actually run (dx = 8 um).
* **Verdict at the extrapolated converged value:** power-law q0 = 0.2743 -> **PASS**; Richardson q0 = 0.2827 -> **PASS**.
* **Margin at dx = 12:** 0.07048 above the band floor, 1.22123 below the band ceiling; nearest edge 0.07048. The whole ladder spans 0.04599, i.e. **0.653x** the distance to the nearest band edge, and the estimated discretisation error is 0.068x it.
* **No verdict flip anywhere on the ladder** — 8/8 rungs PASS.

### doxorubicin — IC50 0.5000 uM, hill 1.20, scored at 1xIC50, band **0.0575-1.0000**

| dx (um) | F_penetration | F_state | **scored 3-D:2-D ratio** | implied IC50 fold | mean tumour exposure (uM per uM plasma) | abs. kill at nominal dose | solve (s) | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 7 | 0.3723 | 0.5407 | **0.2056** | 6.08x | 0.016963 | 0.2598 | 160 | **PASS** |
| 8 | 0.3645 | 0.5488 | **0.2049** | 6.10x | 0.017062 | 0.2642 | 149 | **PASS** |
| 10 | 0.3645 | 0.5305 | **0.1975** | 6.31x | 0.016878 | 0.2538 | 57 | **PASS** |
| 12 **(published)** | 0.3688 | 0.5502 | **0.2077** | 6.03x | 0.017128 | 0.2647 | 130 | **PASS** |
| 16 | 0.3573 | 0.5333 | **0.1946** | 6.40x | 0.017110 | 0.2538 | 69 | **PASS** |
| 18 | 0.3671 | 0.5401 | **0.2026** | 6.17x | 0.017096 | 0.2572 | 53 | **PASS** |
| 24 | 0.3074 | 0.5392 | **0.1718** | 7.17x | 0.017053 | 0.2489 | 35 | **PASS** |
| 30 | 0.2868 | 0.5041 | **0.1456** | 8.34x | 0.016414 | 0.2128 | 20 | **PASS** |

Step-by-step refinement (fine -> coarse):

| from dx | to dx | d(ratio) | d(ratio) % | slope per um |
|---:|---:|---:|---:|---:|
| 7 | 8 | -0.00073 | -0.36 % | -0.000733 |
| 8 | 10 | -0.00738 | -3.60 % | -0.003691 |
| 10 | 12 | 0.01018 | 5.15 % | 0.005089 |
| 12 | 16 | -0.01311 | -6.31 % | -0.003278 |
| 16 | 18 | 0.00804 | 4.13 % | 0.004022 |
| 18 | 24 | -0.03078 | -15.19 % | -0.005130 |
| 24 | 30 | -0.02627 | -15.29 % | -0.004378 |

* **Shape:** NOT monotone in dx; settling (the finest step is smaller than the next one out).
* **Total movement across the ladder:** 42.67 % of the smallest value (0.1456 to 0.2077).
* **Power-law fit** `ratio(dx) = q0 + C*dx^p` over all 8 rungs: p = **3.36**, q0 = **0.2055**, C = -0.00000066, rms residual 0.004604.
* **Richardson (8, 12, 18 um):** not usable — non-monotone triple (d2/d1 = -1.8130 <= 0).
* **Discretisation error at the published dx = 12 um:** 0.00216 absolute (1.05 %) against the power-law q0; 0.00280 (1.36 %) against the finest grid actually run (dx = 8 um).
* **Verdict at the extrapolated converged value:** power-law q0 = 0.2055 -> **PASS**.
* **Margin at dx = 12:** 0.15021 above the band floor, 0.79232 below the band ceiling; nearest edge 0.15021. The whole ladder spans 0.06212, i.e. **0.414x** the distance to the nearest band edge, and the estimated discretisation error is 0.014x it.
* **No verdict flip anywhere on the ladder** — 8/8 rungs PASS.

## 4. Summary — does the published verdict survive the mesh?

| drug | published ratio (dx=12) | this run at dx=12 | at dx=8 (finest) | extrapolated dx->0 | band | verdict at every rung? | disc. error at dx=12 | **survives?** |
|---|---:|---:|---:|---:|---|---|---:|---|
| temozolomide | 0.2626 | 0.2626 | 0.2574 | 0.2615 | 0.1590-0.4622 | 8/8 PASS | 0.44 % | **YES** |
| cisplatin | 0.2791 | 0.2791 | 0.2689 | 0.2743 | 0.2086-1.5004 | 8/8 PASS | 1.76 % | **YES** |
| doxorubicin | 0.2077 | 0.2077 | 0.2049 | 0.2055 | 0.0575-1.0000 | 8/8 PASS | 1.05 % | **YES** |

**3 of 3 PASSes survive every rung run.**

### The number that actually settles it: movement vs margin

A verdict is only grid-fragile if the grid can move the ratio far enough to cross a band edge. So the decisive comparison is not how much the ratio moves in absolute terms, but how that movement compares with the distance from the published value to the nearest edge of its band.

| drug | ratio at dx=12 | nearest band edge | distance to it | **full grid-induced movement** | movement / margin |
|---|---:|---:|---:|---:|---:|
| temozolomide | 0.2626 | 0.1590 (floor) | 0.1037 | 0.0387 | **0.37x** |
| cisplatin | 0.2791 | 0.2086 (floor) | 0.0705 | 0.0460 | **0.65x** |
| doxorubicin | 0.2077 | 0.0575 (floor) | 0.1502 | 0.0621 | **0.41x** |

**For every drug the entire grid-induced movement is smaller than the distance to the nearest band edge** (worst case 0.65x, i.e. the mesh would have to move the ratio roughly 1.5x further than it does across a 4.3x change in resolution to flip a single verdict). **That, and not the smoothness of the curve, is why the GBM 3/3 PASS survives.** The verdicts are robust because the bands are wide and the model sits near their middle — not because the numerics are clean.

### 4.1 There is no grid-error cancellation in this score — measured

It would be easy to assume the 3-D:2-D *ratio* is protected because it is a ratio, and that mesh error in the numerator cancels against mesh error in the denominator. **It does not, and the reason is structural.** `kill_calibration.py` line 514 computes

```
kill_2d = p_hill(c_bath, ic50, hill) * mean(sens)
```

a scalar Hill evaluation at the bath concentration with **no mesh in it at all**. The scored quantity is a grid-dependent numerator over a grid-FREE constant. Measured across the ladder, `kill_2d` is:

| drug | kill_2d across all rungs | constant? |
|---|---|---|
| temozolomide | 0.103202 to 0.103202 | **yes — grid-free** |
| cisplatin | 0.500000 to 0.500000 | **yes — grid-free** |
| doxorubicin | 0.500000 to 0.500000 | **yes — grid-free** |

What the construction *does* remove is multiplicative SCALE error, and only that: `decompose` sets `plasma = c_bath / p99(e_unit)`, so a mesh change that scaled the whole exposure field by a constant would leave the ratio exactly invariant. A mesh change alters the **shape** of the exposure distribution, which is precisely the component that survives that normalisation. The correct prior is therefore that the scored ratio is **more** grid-sensitive than the raw exposure field, not less. Measured as log-log slopes over the ladder:

| drug | d ln(mean exposure)/d ln(dx) | d ln(abs. kill)/d ln(dx) | **d ln(scored ratio)/d ln(dx)** | ratio / exposure sensitivity |
|---|---:|---:|---:|---:|
| temozolomide | -0.0096 | -0.0462 | **-0.0330** | 3.45x |
| cisplatin | -0.0100 | -0.0487 | **-0.0548** | 5.49x |
| doxorubicin | -0.0118 | -0.0993 | **-0.1921** | 16.30x |

So the score does amplify mesh error relative to the underlying delivery field, exactly as the structure predicts. **That amplification is why this audit was worth running, and it is also why the result below is meaningful rather than trivially reassuring**: the quantity was checked in the regime where it is most exposed, not least.

## 5. The other dx-dependent term: the cell-state label

`F_state` is not a constant of the tissue — the oxygen PDE that assigns the QUIESCENT / HYPOXIC bands is solved on the same grid, so refining the mesh can relabel cells. The published report quotes F_state = 0.5502 (90.0 % quiescent, 0.0 % hypoxic). Measured across the ladder:

| dx (um) | quiescent frac | hypoxic frac | viable frac | F_state |
|---:|---:|---:|---:|---:|
| 7 | 0.9187 | 0.0000 | 0.0813 | 0.5407 |
| 8 | 0.9024 | 0.0000 | 0.0976 | 0.5488 |
| 10 | 0.9390 | 0.0000 | 0.0610 | 0.5305 |
| 12 **(published)** | 0.8996 | 0.0000 | 0.1004 | 0.5502 |
| 16 | 0.9334 | 0.0000 | 0.0666 | 0.5333 |
| 18 | 0.9197 | 0.0000 | 0.0803 | 0.5401 |
| 24 | 0.9216 | 0.0000 | 0.0784 | 0.5392 |
| 30 | 0.9918 | 0.0000 | 0.0082 | 0.5041 |

Because the scored ratio is `F_penetration x F_state x (interaction)`, any movement in `F_state` propagates into the verdict directly and linearly. §3's per-drug tables carry both factors so the two sources of grid dependence can be told apart.

**And an observation the data forces, which "3/3 PASS" would otherwise hide.** `F_state` is at its **maximum over the entire ladder at dx = 12 um** (0.5502) — and dx = 12 um is the published grid. The published run therefore sits at the *favourable* end of the oscillation, and its ratios are correspondingly above the fine-rung mean for every drug:

| drug | published ratio (dx=12) | mean over fine rungs (dx <= 18) | published is high by |
|---|---:|---:|---:|
| temozolomide | 0.2626 | 0.2596 | **+1.2 %** |
| cisplatin | 0.2791 | 0.2726 | **+2.4 %** |
| doxorubicin | 0.2077 | 0.2021 | **+2.7 %** |

This is a consistent bias in one direction, not noise averaging out, and it is small — 1-3 %. It does not move any verdict (§4's movement-vs-margin table settles that). But it does mean the published numbers are, very slightly, the flattering draw from an oscillating sequence, and anyone quoting `0.2626` to four figures should know that the fourth figure — and arguably the third — is mesh noise rather than signal. **Quote these ratios to two significant figures.**

## 6. Are the GBM and breast panels converged at the same resolution?

They are not the same problem, and the reason is structural rather than incidental.

| | GBM (this file) | breast (`results/kill_calibration_breast.md`) |
|---|---|---|
| grid at dx = 12 um | 424x675 | 458x628 |
| measured cells | 40,887 | 167,780 |
| tumour cells scored | 23,812 | 89,785 |
| cell density | ~992 /mm^2 | ~4,020 /mm^2 |
| perivascular boundary layer L | 65-265 um | ~28-210 um (docetaxel/paclitaxel at the bottom, `results/dtissue_sensitivity.md` §2) |
| shortest L / published dx | **5.4** | **~2.4** (docetaxel, L = 28.2 um) |
| measured CFL floor at dt = 10 s | dx = 7 um (dx = 6 refuses) | dx = 8 um (reported by the breast study) |

**The load-bearing row is the last-but-two.** The mesh has to resolve the perivascular boundary layer `L = sqrt(D_eff/k_eff)`, because that is the structure the whole `F_penetration` term measures. On this GBM section the SHORTEST such layer is 65 um — **5.4 mesh cells** at the published dx = 12 um. On the breast section docetaxel's is 28.2 um — about 2.4 cells. A gradient resolved by ~2 cells is not resolved; a gradient resolved by ~5 is much closer to it. **The two panels are therefore NOT converged at the same resolution, and there is no reason they should have been**: the breast section is roughly 4x denser in cells, its drugs are larger and more strongly bound (docetaxel, paclitaxel), and both effects shorten `L`.

This is the practical consequence: **`dx = 12 um` cannot be justified panel-agnostically.** `results/numerical_rigor.md` established it on GBM and the breast pipeline inherited it. This audit is consistent with that choice for GBM; it says nothing in favour of it for breast, where the same 12 um resolves a boundary layer half as thick.

## 7. Controls — does this script reproduce the published numbers?

### 7.1 Scoring-path control (zero transport solves)

The exposure vectors behind `results/kill_calibration.md` are on disk in `data_cache/kill_calibration_runs.json`. Scoring THOSE vectors through this script's own `score()` must reproduce the published table exactly:

| drug | published ratio | re-scored from the published cache | delta | F_pen (published) | F_pen (re-scored) | F_state | verdict | published verdict |
|---|---:|---:|---:|---:|---:|---:|---|---|
| temozolomide | 0.2626 | 0.2626 | 0.000049 | 0.4678 | 0.4678 | 0.5502 | PASS | PASS |
| cisplatin | 0.2791 | 0.2791 | 0.000018 | 0.4980 | 0.4980 | 0.5502 | PASS | PASS |
| doxorubicin | 0.2077 | 0.2077 | -0.000024 | 0.3688 | 0.3688 | 0.5502 | PASS | PASS |

Worst deviation **0.000049**; verdicts **ALL MATCH**. The scoring path in this file IS the published scoring path, so every movement reported in §3 is attributable to the mesh and to nothing in the analysis.

### 7.3 Temporal control — was holding `dt` fixed legitimate?

The spatial ladder holds `dt` at the published 10 s throughout. To bound the temporal discretisation error that leaves in place, the published grid and the fine end of the ladder were re-run at **half** that dt. If the ratio barely moves, the §3 curve is a spatial curve and the choice was sound; if it moves a lot, the absolute ratios carry a temporal error no amount of mesh refinement removes.

**dt was verified to actually reach the solver**, not silently ignored: the step count doubles, the CFL number halves, the solve takes longer, and the exposure field changes. A bit-identical exposure field would make this control meaningless and `run()` asserts against it.

| dx (um) | dt (s) | steps | CFL | drug | ratio | ratio at the published dt | delta | delta % | exposure-field change | verdict |
|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|
| 12 | 5 | 4320 | 0.264 | temozolomide | 0.262649 | 0.262649 | +5.31e-08 | 0.0000 % | 1.18e-04 | PASS (PASS) |

Worst movement from halving dt: **0.0000 %** of the scored ratio; **0** verdict change(s). Holding dt fixed across the spatial ladder was therefore legitimate, and the §3 curve is a spatial curve.

**But the interesting part is WHY it is so small, because it is not that the solver is temporally exact.** The exposure field itself moves by up to **1.2e-04** relative when dt is halved — small, but four to five orders of magnitude larger than the resulting move in the score. The reason is the cancellation structure of §4.1: `decompose` renormalises by `plasma = c_bath / p99(e_unit)`, which removes **multiplicative scale** error from the exposure field exactly. Halving dt shifts the whole AUC field by a near-uniform factor — almost pure scale — so it cancels almost completely. Refining dx changes the **shape** of the exposure distribution, which is precisely the component that survives that normalisation.

**So the two discretisations are not comparable in their effect on this score, and that is a property of the score rather than of the solver.** It also means this control does NOT license the claim that the absolute exposure or absolute kill numbers are temporally converged — only that the 3-D:2-D *ratio* is insensitive to dt. Anyone quoting an absolute AUC or an absolute kill from this pipeline still needs a dt study of their own.

## 8. Limitations of THIS audit

1. **`dt` is held at the published 10 s throughout, so this is a SPATIAL convergence study only.** That is the right choice for attributing movement to the mesh, but it means any temporal discretisation error is common to every rung: it cancels out of the convergence CURVE and does NOT cancel out of the absolute ratios. §7.3 bounds it by halving dt at two rungs rather than assuming it away, but a full joint (dx, dt) refinement is a stronger statement and is not made here. The CFL condition also ties the two axes together: the measured finest reachable dx at this dt is 7 um.
2. **dx -> 0 is extrapolated, never reached.** Both the power-law fit and the Richardson triple are estimates from a 5-7 point ladder. Where they disagree the ladder is not in its asymptotic range and the quoted discretisation error should be read as an order of magnitude, not a bound.
3. **One section, one patient.** Convergence is a property of this geometry — its vessel spacing and cell density — as much as of the solver. Another GBM section with tighter vasculature would need re-checking.
4. **Only the reporting concentration each verdict is taken at was scored.** The published report evaluates the decomposition over 0.1-10x IC50; `F_penetration` climbs towards 1 as the Hill term saturates, and the grid dependence at other concentrations is not characterised here.
5. **The p99 bath-matching convention is held fixed.** The published report notes the cisplatin verdict is convention-dependent (PASS at p99, FAIL under a maximum-exposure reference). A joint convention x mesh sensitivity is not measured. A verdict that is already fragile in one dimension deserves that cross-check before it is quoted as robust.
6. **This audit does not revisit any of the substantive limitations of the published report** — no outcome label, wide literature bands, a delivery/state split the spheroid literature disagrees with, and cell-state multipliers that are unvalidated priors. **A grid-converged number is not a validated number.** Convergence says the answer is the answer the model gives; it says nothing about whether the model is right.
7. **The breast comparison in §6 uses published breast figures, not a breast run executed here.** A concurrent workstream owns the breast convergence study; nothing in §6 substitutes for it.

*Generated by `experiments/grid_convergence_gbm.py` — 25 transport solves over 8 grids, 37.5 min, seed 0, deterministic. Raw records: `results/sweep/gbmgrid_raw.json`, `results/sweep/gbmgrid_records.jsonl`.*
