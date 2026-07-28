# Defensible delivery metrics — replacing `penetration_index`

> **Regenerated under the discrete vessel source.** `vascular_mode = discrete`
> (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in
> [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6**, the
> e-folding-length convention `L = d_half / ln 2` · per-cell sensitivity mode =
> **`uniform` (pipeline default — identical for every drug, so the cross-drug discrimination test is fair)** · **dx = 20 / 15 / 12 µm ladder (discrimination table quoted at dx = 12 µm)**.
>
> Every absolute exposure / AUC / coverage / kill below is on the **discrete-source
> unit** and is NOT comparable with the KDE-era numbers this file carried before this
> regeneration (the tissue-integrated source total is 3.36x smaller; that is a change
> of *unit*, not of physics — the dose-matched, geometry-only effect on coverage is
> <= 0.0166). To recover the old numbers, re-run with `mode='kde'`.
>
> **This regeneration corrects a claim that was wrong in the published file.** The KDE-era version listed `penetration_length_L` as *does NOT discriminate* on a 1.07x five-drug spread against a source null of 23.7 µm. On the discrete source the spread is **6.91x (39.6–273.3 µm)**, separation **12.2x** its own grid drift, and the source field **admits no null fit at all** (null L = nan) — so the metric now reads **YES, discriminates**. That reversal was predicted by `vessel_promotion.md` §5 (3.12x -> 6.91x) and flagged by `vessel_network_impact.md` §3; it is now actually fixed in this file. The headline metric is unchanged: `therapeutic_coverage` @1x IC50, range 0.0000–0.9936, separation **134x** grid drift.

**Module:** `tmesim/validation/delivery_metrics.py` (run: `./.venv/bin/python -m tmesim.validation.delivery_metrics`)

**Metric self-test:** `./.venv/bin/python -m tmesim.validation.delivery_metrics --selftest` (known-answer checks on the definitions, no simulation)

**Figure:** `results/delivery_metrics.png`

**Geometry:** the REAL open 10x Xenium GBM section — 40,887 measured cells (tumor 30,174), reconstructed by `reconstruct.real_pipeline.reconstruct_real('xenium', dx_um=…)`; real endothelial-KDE vascular field, real matrix-transcript ECM, physics-solved oxygen. **Drugs:** the 5-drug panel of [`results/real_drug_panel.md`](real_drug_panel.md), each compiled through the identical pipeline (`calibrate_alpha()` → α = 148.6, fit ONCE on doxorubicin → `calibrated_drugparams` → `apply_barrier_to_drugparams`). **Engine:** advection on, dt = 10 s, 6 h window, seed 0 — unchanged from the flagship run. Production resolution **dx = 12 µm**, grid (424, 675) (dx = 10 µm violates the advection CFL condition at dt = 10 s and is therefore not on the ladder).

---

## 0. The problem being fixed

`penetration_index = mean(AUC) / p99(AUC)` over the tissue mask is a **scale-invariant shape statistic**: multiply the entire exposure field by any constant and it does not move. At the grid-converged resolution it therefore reports the *vascular geometry of the section*, which is identical for all 5 drugs — and it duly collapses:

> at dx = 12 µm, `penetration_index` = **0.310–0.470** across the whole panel (spread 0.1606), while those same 5 drugs' mean delivered exposure differs by **156×**.

A delivery metric must be **absolute** (it scales when the delivered dose scales) and, ideally, **drug-specific** (it knows the drug's potency). The five metrics below are defined to be exactly that, and then audited on both axes that matter: grid convergence and drug discrimination.

---

## 1. Definitions, units, rationale

Let `AUC(x)` be the transport solve's exposure field (µM·s) and `T` the exposure-window length (s; here T = 21600 s = 6 h, taken from `transport.times[-1]`, the same `T` `cells.pkpd` uses). The **time-averaged effective concentration** is

```
E(x) = AUC(x) / T        [µM]
```

which is directly comparable to an IC50 because an IC50 *is* a concentration. Every per-cell statistic is evaluated at the **real measured cell positions** (`domain.to_index`, the same nearest-cell map the kill model uses), not over grid pixels — the tumor population is a FIXED set of 30,174 Xenium cells at every dx, so a per-cell statistic moves with dx only because the *field* moved, not because the rasterised mask changed.

| # | metric | units | what it answers | absolute or shape? |
|---|---|---|---|---|
| 1 | `mean_tissue_auc`, `median_tissue_auc`, `mean_tumor_auc`, `median_tumor_auc` | µM·s | how much drug arrived, in total | **absolute** |
| 2 | `therapeutic_coverage(auc, ic50, T)` — fraction of tumor cells with `E > m·IC50`, m ∈ {0.5, 1, 2} | dimensionless fraction ∈ [0,1] | what fraction of the tumor saw a **killing** dose | **absolute + drug-specific** |
| 3 | `exposure_percentiles` — P10/P50/P90 of per-tumor-cell `E` | µM | how *unevenly* the dose was delivered; P10 = the worst-served decile | **absolute** |
| 4 | `penetration_length_fit` — `L` from `AUC(d)=C₀·exp(−d/L)+floor` vs distance-from-nearest-vessel | µm | a real physical decay length | length (scale-invariant, but dimensional) |
| 5 | `hypoxic_vs_normoxic_exposure_ratio` — mean `E` in the worst- vs best-perfused quartile of tumor cells (plus the literal hypoxic/normoxic bands) | dimensionless | does the drug reach the poorly-perfused tissue? | *shape* (diagnostic only) |

**Rationale, metric by metric.**

1. **Absolute exposure (µM·s).** The integral the pharmacodynamics actually consumes. The median is reported next to the mean because the field is strongly right-skewed — a thin perivascular cuff carries much of the integral, so mean/median ≈ 1.1× on this section and the mean alone overstates what a typical tumor cell sees. Honest limitation: an exposure in µM·s is not interpretable without a potency, which is why it is not the headline.
2. **Therapeutic coverage.** The clinically meaningful reduction of a whole exposure *field* to one number: the fraction of tumor cells whose time-averaged concentration clears the drug's own IC50. Dimensionally sound (concentration vs concentration), bounded in [0,1], and — crucially — the same exposure field scores differently for a potent and a weak drug, so it **cannot** collapse to a geometry constant. The 0.5× and 2× variants say how hard the answer leans on the exact IC50 value.
3. **Exposure percentiles.** A mean hides the population that decides the outcome. P10 is the exposure the worst-served decile of the tumor receives — the reservoir that survives and regrows — and P90/P10 is the delivery heterogeneity.
4. **Penetration length.** A genuine physical length in µm, fitted to the binned AUC-vs-distance-from-vessel profile (the distance map, 10 µm binning and fits are **reused** from `validation/vessel_profile.py`, so the two reports cannot drift apart). It ships with its own **SOURCE NULL** — the identical fit on the raw `vascular_density` field — because the engine's delivery term is `vascular_density · bbb_penetrance · plasma`, i.e. a source spread over the whole continuous density field; with a fast perivascular sink the tissue field is close to locally reaction-dominated (`C ≈ S/k`) and largely **inherits the source's shape**.
5. **Hypoxic/normoxic ratio.** Two forms. The **primary, threshold-free** one ranks tumor cells by their own local oxygen and divides the mean `E` of the worst-perfused quartile by that of the best-perfused quartile — a partition that always exists and does not depend on where a hypoxia cutoff is drawn. The **banded** form uses the engine's own oxygen banding (`t = OxygenConfig.hypoxia_threshold` = 0.2; hypoxic `o < t` vs normoxic `o ≥ t`) and is **nan whenever a band is empty**. Stated at the definition: both are *ratios*, hence scale-invariant, hence structurally the same kind of statistic as `penetration_index`; the banded one additionally sits on a threshold crossing on an oxygen field whose hypoxic fraction is documented as **not** grid-converged. Included because the question is real, and reported as a diagnostic — not as a headline.

---

## 2. Grid convergence on the REAL section

Every drug re-simulated on the same section re-reconstructed at dx = 20, 15, 12 µm. The table reports the **panel-mean** value at each dx and, for the drift columns, the **worst drug** in the panel — the conservative reading.

**The convergence criterion, stated before the numbers.** Two scales are reported for the final refinement step, and which one decides is fixed by the metric's type, not by the answer:

- metrics **bounded in [0,1]** (the coverage and kill fractions) are judged on the **absolute** change, tolerance **0.03** (= 3 percentage points of tumor cells). Relative error is ill-conditioned for a bounded fraction — it diverges as the fraction goes to zero even when the answer moved by a fraction of a point.
- everything else (µM·s, µM, µm, dimensionless ratios) is judged on the **relative** change, tolerance **5%**.

Both columns are in the table for every metric, so any other rule can be applied to the same numbers.

| metric | units | dx 20 µm | dx 15 µm | dx 12 µm | Δrel 20→15 | Δrel 15→12 | Δabs 15→12 | converged? |
|---|---|---:|---:|---:|---:|---:|---:|:--:|
| `mean_tissue_auc` | uM*s | 6.359e+05 | 6.217e+05 | 6.351e+05 | 0.0222 | 0.0215 | 3.43e+04 | **YES** |
| `median_tissue_auc` | uM*s | 5.881e+05 | 5.885e+05 | 6.048e+05 | 0.0375 | 0.0364 | 3.94e+04 | **YES** |
| `mean_tumor_auc` | uM*s | 6.656e+05 | 6.475e+05 | 6.585e+05 | 0.0276 | 0.0176 | 2.81e+04 | **YES** |
| `median_tumor_auc` | uM*s | 5.995e+05 | 6.012e+05 | 6.144e+05 | 0.0284 | 0.0353 | 2.99e+04 | **YES** |
| `therapeutic_coverage_0.5x` | fraction | 0.5003 | 0.5218 | 0.5256 | 0.0845 | 0.0164 | 0.0133 | **YES** *(abs)* |
| `therapeutic_coverage_1x` | fraction | 0.3224 | 0.3152 | 0.319 | 0.132 | 0.0352 | 0.00739 | **YES** *(abs)* |
| `therapeutic_coverage_2x` | fraction | 0.1895 | 0.1868 | 0.1877 | 0.333 | 0.0786 | 0.00467 | **YES** *(abs)* |
| `median_E/IC50 (potency margin)` | - | 0.9058 | 0.9149 | 0.9313 | 0.0284 | 0.0353 | 0.0563 | **YES** |
| `P10_E/IC50 (worst-decile margin)` | - | 0.4804 | 0.5271 | 0.5252 | 0.24 | 0.0418 | 0.0292 | **YES** |
| `exposure_P10` | uM | 16.09 | 17.18 | 16.99 | 0.24 | 0.0418 | 0.716 | **YES** |
| `exposure_P50` | uM | 27.76 | 27.83 | 28.44 | 0.0284 | 0.0353 | 1.38 | **YES** |
| `exposure_P90` | uM | 48.73 | 43.95 | 44.55 | 0.124 | 0.0272 | 1.65 | **YES** |
| `exposure_P90/P10` | - | 4.157 | 3.212 | 3.228 | 0.294 | 0.0313 | 0.0796 | **YES** |
| `penetration_length_L` | um | 124.2 | 124.8 | 118.4 | 0.0285 | 0.162 | 12.9 | no |
| `lowO2/highO2_exposure (quartiles)` | - | 0.5901 | 0.6504 | 0.6546 | 0.132 | 0.015 | 0.00928 | **YES** |
| `hypoxic/normoxic_exposure (banded)` | - | 0.4668 | 0.5273 | 0.513 | 0.161 | 0.0286 | 0.0146 | **YES** |
| `penetration_index (RETIRED)` | - | 0.3629 | 0.3773 | 0.3803 | 0.112 | 0.0171 | 0.00706 | **YES** |
| `kill_fraction (reference)` | fraction | 0.1616 | 0.1616 | 0.1629 | 0.0578 | 0.0126 | 0.00245 | **YES** *(abs)* |

**Converged (worst drug, 15→12 µm):** `mean_tissue_auc`, `median_tissue_auc`, `mean_tumor_auc`, `median_tumor_auc`, `therapeutic_coverage_0.5x`, `therapeutic_coverage_1x`, `therapeutic_coverage_2x`, `median_E/IC50 (potency margin)`, `P10_E/IC50 (worst-decile margin)`, `exposure_P10`, `exposure_P50`, `exposure_P90`, `exposure_P90/P10`, `lowO2/highO2_exposure (quartiles)`, `hypoxic/normoxic_exposure (banded)`, `penetration_index (RETIRED)`, `kill_fraction (reference)`.

**NOT converged:** `penetration_length_L`.

**Read this honestly for the headline metric.** `therapeutic_coverage_1x` moves **0.0074 absolute** (0.7 percentage points of tumor cells) on its worst drug between 15 µm and 12 µm, which is inside the 0.03 tolerance — but that same move is **3.5% relative**, and it is the *least*-covered drug (temozolomide, coverage ≈ 0.31) that produces it, which is exactly the regime where relative error is misleading. The residual is real and monotone, not noise: carry it as a **±0.03 error bar on any quoted coverage**. The well-covered drugs are far tighter (cisplatin 0.0044, doxorubicin 0.0069 absolute).

**Amendment — the grid bar is not the largest bar this metric carries.** [`results/dimensionality.md`](dimensionality.md) quantified the 2-D-section assumption against a standalone 3-D solve and found `therapeutic_coverage` to be **the one published quantity that moves by materially more than the uncertainty it already carried**: up to **0.137 absolute at the operating points published in §5 below** (0.205 anywhere on the threshold ladder), i.e. **~7× the ±0.03 grid bar**. The bias is *structured*, not noise — 2-D under-counts at permissive thresholds and over-counts at strict ones, because the plane is blind to the out-of-plane vessels that fill the gaps and so sees a broader-tailed field than the tissue really has. **Quote coverage with both bars, the dimensionality one beside the grid one, not instead of it.** What does not change: the drug *ordering* (the two drugs at exactly 0.0000 are floored and cannot move; the gaps between the rest exceed the shifts), and `mean_tumor_auc`, which the same report proves **unbiased analytically** by a stereological identity. `penetration_length_L` picks up a further ~24% under-statement there — comparable to the 16% grid drift for which it is already listed as NOT converged above.

For contrast, `penetration_index` moves 1.7% relative over the same refinement step while its entire spread across the 5-drug panel at dx = 12 µm is only 42.2% — its grid drift is 0× larger than the entire signal it is supposed to carry.

---

## 3. Discrimination across the 5-drug panel (dx = 12 µm)

### 3a. The per-drug table

| drug | MW (Da) | real IC50 (µM) | mean tumor AUC (µM·s) | median tumor AUC (µM·s) | E P10 (µM) | E P50 (µM) | E P90 (µM) | **coverage @0.5× IC50** | **coverage @1× IC50** | **coverage @2× IC50** | median E/IC50 | P10 E/IC50 | L (µm) | lowO₂/highO₂ E | `penetration_index` | kill_fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194 | 606 | 7.06e+05 | 6.89e+05 | 20.8 | 31.9 | 43.6 | 0.0000 | **0.0000** | 0.0000 | 0.0526 | 0.0343 | 273.3 | 0.703 | 0.4703 | 0.0147 |
| lomustine | 234 | 50 | 8.7e+05 | 7.86e+05 | 20.1 | 36.4 | 62.8 | 0.8043 | **0.2175** | 0.0209 | 0.728 | 0.401 | 78.2 | 0.645 | 0.3698 | 0.1975 |
| cisplatin | 298 | 24.6 | 1.69e+06 | 1.58e+06 | 43.6 | 72.9 | 114 | 0.9978 | **0.9936** | 0.8535 | 2.97 | 1.78 | 152.5 | 0.662 | 0.4203 | 0.3699 |
| methotrexate | 454 | 12 | 1.48e+04 | 1.29e+04 | 0.293 | 0.597 | 1.13 | 0.0000 | **0.0000** | 0.0000 | 0.0498 | 0.0244 | 48.7 | 0.634 | 0.3313 | 0.0162 |
| doxorubicin | 544 | 0.5 | 1.08e+04 | 9.24e+03 | 0.194 | 0.428 | 0.86 | 0.8256 | **0.3836** | 0.0641 | 0.856 | 0.389 | 39.6 | 0.629 | 0.3097 | 0.2165 |

`L` is the decaying component of `C₀·exp(−d/L)+floor` fitted to the binned AUC-vs-distance profile at the p95 vessel threshold (2.5% of tissue area). **SOURCE NULL** (identical fit on the raw `vascular_density` field): L = nan µm.

### 3b. Spread vs numerical noise

A metric only ranks drugs credibly if its across-drug spread is much larger than its own grid drift. `separation = spread / worst-drug grid drift`, both taken on the **same scale** (absolute for the [0,1] fractions, relative for everything else — the rule stated in §2); we require ≥ 3. A metric that fails here is reporting less signal than it carries numerical noise.

| metric | min | max | fold range (max/min) | spread (abs) | spread (rel) | grid drift 15→12 µm (matched scale) | **separation** | discriminates? |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|
| `mean_tissue_auc` | 1.035e+04 | 1.63e+06 | 157.5 | 1.62e+06 | 2.55 | 0.0215 | **118** | **YES** |
| `median_tissue_auc` | 9109 | 1.551e+06 | 170.3 | 1.54e+06 | 2.55 | 0.0364 | **70** | **YES** |
| `mean_tumor_auc` | 1.083e+04 | 1.691e+06 | 156.1 | 1.68e+06 | 2.55 | 0.0176 | **145** | **YES** |
| `median_tumor_auc` | 9241 | 1.575e+06 | 170.5 | 1.57e+06 | 2.55 | 0.0353 | **72.1** | **YES** |
| `therapeutic_coverage_0.5x` | 0 | 0.9978 | n/a | 0.998 | 1.9 | 0.0133 | **75.1** | no |
| `therapeutic_coverage_1x` | 0 | 0.9936 | n/a | 0.994 | 3.12 | 0.00739 | **134** | no |
| `therapeutic_coverage_2x` | 0 | 0.8535 | n/a | 0.853 | 4.55 | 0.00467 | **183** | no |
| `median_E/IC50 (potency margin)` | 0.04979 | 2.971 | 59.66 | 2.92 | 3.14 | 0.0353 | **88.8** | **YES** |
| `P10_E/IC50 (worst-decile margin)` | 0.02439 | 1.777 | 72.88 | 1.75 | 3.34 | 0.0418 | **79.8** | **YES** |
| `exposure_P10` | 0.1943 | 43.63 | 224.6 | 43.4 | 2.56 | 0.0418 | **61.1** | **YES** |
| `exposure_P50` | 0.4278 | 72.93 | 170.5 | 72.5 | 2.55 | 0.0353 | **72.1** | **YES** |
| `exposure_P90` | 0.8596 | 114.4 | 133.1 | 114 | 2.55 | 0.0272 | **93.7** | **YES** |
| `exposure_P90/P10` | 2.097 | 4.424 | 2.11 | 2.33 | 0.721 | 0.0313 | **23** | **YES** |
| `penetration_length_L` | 39.58 | 273.3 | 6.906 | 234 | 1.97 | 0.162 | **12.2** | **YES** |
| `lowO2/highO2_exposure (quartiles)` | 0.6287 | 0.7029 | 1.118 | 0.0742 | 0.113 | 0.015 | **7.57** | **YES** |
| `hypoxic/normoxic_exposure (banded)` | 0.4946 | 0.5483 | 1.109 | 0.0537 | 0.105 | 0.0286 | **3.66** | **YES** |
| `penetration_index (RETIRED)` | 0.3097 | 0.4703 | 1.519 | 0.161 | 0.422 | 0.0171 | **24.7** | **YES** |
| `kill_fraction (reference)` | 0.01468 | 0.3699 | 25.19 | 0.355 | 2.18 | 0.00245 | **145** | **YES** |

**Discriminates:** `mean_tissue_auc`, `median_tissue_auc`, `mean_tumor_auc`, `median_tumor_auc`, `median_E/IC50 (potency margin)`, `P10_E/IC50 (worst-decile margin)`, `exposure_P10`, `exposure_P50`, `exposure_P90`, `exposure_P90/P10`, `penetration_length_L`, `lowO2/highO2_exposure (quartiles)`, `hypoxic/normoxic_exposure (banded)`, `penetration_index (RETIRED)`, `kill_fraction (reference)`.

**Does NOT discriminate:** `therapeutic_coverage_0.5x`, `therapeutic_coverage_1x`, `therapeutic_coverage_2x`. Three of these are informative failures, not bugs:

- **`penetration_length_L`** — the fitted drug `L` spans only 39.6–273.3 µm (6.91×) against a source-null `L` of nan µm, while its own grid drift is 16%. On real tissue the measured decay length is set by the endothelial-KDE geometry, not by drug transport — the confound `results/vessel_profile.md` quantifies at length. The metric is sound; this *section* cannot resolve drugs with it.
- **`lowO2/highO2_exposure`** — grid-stable (1.5% drift) but essentially identical for every drug (0.629–0.703). That is the correct answer, not a failure: the perfusion gradient is a property of the *tissue*, and all five drugs are delivered proportionally to the same vascular field, so they all under-serve the badly-perfused quartile by the same ~3×.
- **`penetration_index`** — the metric this module replaces. See §4.

---

## 4. The explicit contrast with `penetration_index`

| | `penetration_index` (retired) | `therapeutic_coverage` @1× IC50 (proposed headline) |
|---|---|---|
| definition | mean(AUC)/p99(AUC) over the tissue mask | fraction of tumor cells with AUC/T > IC50 |
| units | dimensionless *shape* | dimensionless *fraction of cells* |
| scales with delivered dose? | **NO** (scale-invariant by construction) | **YES** |
| knows the drug's potency? | **NO** | **YES** (IC50 enters explicitly) |
| range across the 5-drug panel | **0.3097 – 0.4703** (spread 0.1606, fold 1.519×) | **0.0000 – 0.9936** (spread 0.9936, fold nan×) |
| grid drift 15→12 µm (worst drug), absolute | 0.00706 | 0.00739 |
| grid drift 15→12 µm (worst drug), relative | 0.0171 | 0.0352 |
| **separation** (spread / drift, matched scale) | **24.7** | **134** |

Both are dimensionless numbers in [0,1] on the same tissue, so the absolute row is the like-for-like comparison: `penetration_index` drifts 0.00706 with grid refinement while its whole across-drug spread is 0.1606 — the numerics move it ~0× more than the drugs do. Coverage drifts 0.00739 against a spread of 0.994 — ~134× the other way round.

The scale-invariance claim is not rhetorical. `delivery_metrics.selftest()` (run it: `./.venv/bin/python -m tmesim.validation.delivery_metrics --selftest`) asserts arithmetically that multiplying an exposure field by 10⁻³ and by 10³ leaves `mean/p99` identical to 1e-12 while `therapeutic_coverage` moves from ~0.000 to 1.000 over the same 10⁶× dose change.

Per drug, side by side:

| drug | `penetration_index` | mean tumor AUC (µM·s) | coverage @1× IC50 |
|---|---:|---:|---:|
| temozolomide | 0.4703 | 7.06e+05 | 0.0000 |
| lomustine | 0.3698 | 8.7e+05 | 0.2175 |
| cisplatin | 0.4203 | 1.69e+06 | 0.9936 |
| methotrexate | 0.3313 | 1.48e+04 | 0.0000 |
| doxorubicin | 0.3097 | 1.08e+04 | 0.3836 |

The index moves by 0.1606 over a panel whose delivered exposure moves by 156×. It is not measuring delivery. Coverage moves by 0.994 on the same panel.

---

## 5. Recommendation — the headline delivery number

> **Headline: `therapeutic_coverage` at 1× IC50 — the fraction of tumor cells whose time-averaged exposure `E = AUC/T` exceeds the drug's real IC50.** Quote it at dx = 12 µm, always with the 0.5×/2× IC50 band beside it, with **`median E/IC50`** (the continuous potency margin) as its companion, and with `median_tumor_auc` (µM·s) as the supporting absolute number.

Justification, against the three requirements:

1. **Physically meaningful.** It compares a concentration to a concentration. `E = AUC/T` is exactly the quantity `cells.pkpd` feeds into the Hill dose-response, and `IC50` is by definition the concentration giving half-maximal effect — so the metric is the delivery half of the kill model, read out directly rather than through a stochastic Bernoulli draw. It is bounded in [0,1], monotone in delivered dose, and one sentence of plain English ("x% of the tumor saw at least a half-maximal-effect concentration").
2. **Grid-converged — with a stated caveat.** Worst-drug change from 15 µm to 12 µm is **0.0074 absolute** (0.7 percentage points of tumor cells), inside the 0.03 tolerance for a bounded fraction. That same step is 3.5% *relative*, driven entirely by the least-covered drug — so the honest statement is **converged to ±0.03 in coverage units, not to 5% relative**, and coverage should be quoted with that error bar. Unlike a mask-area average it is computed over a FIXED population of real cells, so grid refinement cannot change the denominator — which is why it is this stable at all.
3. **Discriminating.** It spans 0.0000–0.9936 across the panel (0.994 absolute, nan× fold) — a separation of **134×** its own grid drift, versus 24.7 for `penetration_index`. It also ranks the drugs the way the mechanism says it should: potency and delivery are combined, so a drug wins only when it both arrives AND is potent at the concentration it achieves.

**What to quote alongside it (and why not instead of it):**

- `median_tumor_auc` (µM·s) — the raw absolute delivery. Strongly discriminating (170× across the panel) though not itself grid-converged on this ladder (3.5% worst-drug drift; use `mean_tumor_auc`, 1.8%, if a converged absolute number is required), and in any case *not* interpretable without a potency: 1.6e+06 µM·s of cisplatin (IC50 24.6 µM) and 9.2e+03 µM·s of doxorubicin (IC50 0.5 µM) are not comparable doses. Report it as the input to coverage, not as the verdict.
- `median E/IC50` and `P10 E/IC50` (dimensionless) — the potency-normalised exposure margin at the median and at the worst-served decile. Same information content as coverage, continuous instead of thresholded, and immune to saturation; quote `median E/IC50` whenever coverage pins at 1.000. It inherits the median's grid sensitivity (3.5% worst-drug drift, versus 0.007 absolute for coverage), which is the price of being continuous — it is a ranking companion, not a replacement headline.
- `E_p10` (µM) — the worst-served decile, the honest floor of the delivery distribution and the number that should worry you about regrowth.
- `penetration_length_fit` `L` (µm) — quote only against the SOURCE NULL. On this section the fitted drug `L` spans 39.6–273.3 µm against a source-null `L` of nan µm (drug/null ratio nan–nan×), i.e. the measured length is dominated by vessel/KDE geometry rather than by drug transport — the same confound `results/vessel_profile.md` documents at length. The drug-intrinsic length separates cleanly only in the isolated single-vessel control, where doxorubicin's gradient is ~4× steeper than temozolomide's. `L` is a legitimate *physical* number and a legitimate *model* diagnostic; it is not a delivery ranking on this section.
- `lowO₂/highO₂ exposure ratio` — diagnostic only. On this section it is 0.629–0.703 across the panel. It is a ratio (hence scale-invariant, hence the same structural weakness as `penetration_index`) and the banded hypoxic/normoxic form sits on a partition that is itself not grid-converged. Read it as 'is delivery correlated with perfusion?', never as a drug ranking.

**What to retire.** `penetration_index` should no longer be quoted as a delivery number anywhere. It is not wrong — it is a faithful measurement of *exposure-field uniformity* — but on real tissue at converged resolution that uniformity is a property of the reconstructed vasculature, shared by every drug. If a uniformity statistic is wanted, quote it as `exposure_uniformity_index` and say what it is.

---

## 6. Honesty and caveats

- **No outcome label.** The Xenium section is an untreated resection. Every number here is a PREDICTION of the calibrated physics on real geometry, not a validated response. Coverage being grid-converged and drug-discriminating says it is a *well-posed* metric; it does not say the underlying delivery model is right.
- **Coverage inherits the IC50's uncertainty.** The IC50s are real GDSC / literature potencies measured in cell culture, not in this tissue. That is why the 0.5× and 2× variants are reported as a band and not buried: they are the metric's sensitivity to the one external number it depends on. Where the 0.5×–2× band is wide, the ranking at that position is soft.
- **Coverage can saturate — this is its one structural weakness, and it is why the continuous companion is reported next to it.** A drug delivering far above its IC50 *everywhere* scores 1.000 and stops discriminating. On this panel that does NOT bite — 0 of 5 drugs reach 0.999 and the maximum is 0.994 — so the weakness is latent here, not active. It will bite on a higher-dose or more potent panel. `median E/IC50` — the **potency-normalised exposure margin** — is the same comparison without the threshold: dimensionless, drug-specific, absolute, unbounded above, and it spans 60× across the panel. Read them together: coverage says *how much of the tumor* cleared the bar, the margin says *by how far*. Coverage is the right headline precisely in the regime that matters clinically — where part of the tumor is under-dosed — and the margin is the tiebreaker when it is not.
- **One anchor sets the absolute scale.** The single α (= 148.6 here, read live from `validation.calibration.calibrate_alpha()`) is fit to one doxorubicin penetration measurement; every absolute exposure — and therefore every coverage value — rides on it. Coverage is *relatively* robust because the same α applies to all 5 drugs, but the absolute level is not independently validated. **Older reports in `results/` may quote a different α** (the half-distance vs e-folding-length convention was corrected in `calibration.py`, a factor 1/ln²2 = 2.08×); the absolute µM·s figures here are therefore not directly comparable to those, while the *structural* conclusions — the index collapsing, coverage spreading — are about the relative behaviour of the metrics and are unaffected.
- **The tumor hypoxic fraction is not grid-converged** (`results/real_drug_panel.md`), which is exactly why metric 5 is a diagnostic and why its primary form is the threshold-free quartile contrast. On the reconstruction as it stands the hypoxic tumor fraction is: 20 µm → 12.2%, 15 µm → 12.7%, 12 µm → 12.7%. Where it is 0%, the banded hypoxic/normoxic ratio is legitimately undefined (`n/a`) — that is reported, not papered over.
- **2D section, straight-line distances, static plasma level.** All the standing limitations of the engine (`results/numerical_rigor.md`) apply unchanged; nothing in `core` / `physics` / `molecular` / `reconstruct` / `sim.engine` was modified for this audit.

*Total runtime 851 s (3 grids × 5 drugs).*
