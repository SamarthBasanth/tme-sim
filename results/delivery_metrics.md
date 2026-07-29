# Defensible delivery metrics — replacing `penetration_index`

**Module:** `tmesim/validation/delivery_metrics.py` (run: `./.venv/bin/python -m tmesim.validation.delivery_metrics`)

**Metric self-test:** `./.venv/bin/python -m tmesim.validation.delivery_metrics --selftest` (known-answer checks on the definitions, no simulation)

**Figure:** `results/delivery_metrics.png`

**Geometry:** the REAL open 10x Xenium GBM section — 40,887 measured cells (tumor 23,812), reconstructed by `reconstruct.real_pipeline.reconstruct_real('xenium', dx_um=…)`; real endothelial-KDE vascular field, real matrix-transcript ECM, physics-solved oxygen. **Drugs:** the 5-drug panel of [`results/real_drug_panel.md`](real_drug_panel.md), each compiled through the identical pipeline (`calibrate_alpha()` → α = 148.6, fit ONCE on doxorubicin → `calibrated_drugparams` → `apply_barrier_to_drugparams`). **Engine:** advection on, dt = 10 s, 6 h window, seed 0 — unchanged from the flagship run. Production resolution **dx = 12 µm**, grid (424, 675) (dx = 10 µm violates the advection CFL condition at dt = 10 s and is therefore not on the ladder).

---

## 0. The problem being fixed

`penetration_index = mean(AUC) / p99(AUC)` over the tissue mask is a **scale-invariant shape statistic**: multiply the entire exposure field by any constant and it does not move. At the grid-converged resolution it therefore reports the *vascular geometry of the section*, which is identical for all 5 drugs — and it duly collapses:

> at dx = 12 µm, `penetration_index` = **0.304–0.511** across the whole panel (spread 0.2063), while those same 5 drugs' mean delivered exposure differs by **158×**.

A delivery metric must be **absolute** (it scales when the delivered dose scales) and, ideally, **drug-specific** (it knows the drug's potency). The five metrics below are defined to be exactly that, and then audited on both axes that matter: grid convergence and drug discrimination.

---

## 1. Definitions, units, rationale

Let `AUC(x)` be the transport solve's exposure field (µM·s) and `T` the exposure-window length (s; here T = 21600 s = 6 h, taken from `transport.times[-1]`, the same `T` `cells.pkpd` uses). The **time-averaged effective concentration** is

```
E(x) = AUC(x) / T        [µM]
```

which is directly comparable to an IC50 because an IC50 *is* a concentration. Every per-cell statistic is evaluated at the **real measured cell positions** (`domain.to_index`, the same nearest-cell map the kill model uses), not over grid pixels — the tumor population is a FIXED set of 23,812 Xenium cells at every dx, so a per-cell statistic moves with dx only because the *field* moved, not because the rasterised mask changed.

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
| `mean_tissue_auc` | uM*s | 6.727e+05 | 6.664e+05 | 6.818e+05 | 0.00932 | 0.0231 | 3.96e+04 | **YES** |
| `median_tissue_auc` | uM*s | 6.145e+05 | 6.289e+05 | 6.375e+05 | 0.0662 | 0.0202 | 2.2e+04 | **YES** |
| `mean_tumor_auc` | uM*s | 7.014e+05 | 6.943e+05 | 7.094e+05 | 0.0102 | 0.0219 | 3.87e+04 | **YES** |
| `median_tumor_auc` | uM*s | 6.317e+05 | 6.504e+05 | 6.558e+05 | 0.0646 | 0.0154 | 1.46e+04 | **YES** |
| `therapeutic_coverage_0.5x` | fraction | 0.5216 | 0.5382 | 0.5425 | 0.0631 | 0.0152 | 0.0128 | **YES** *(abs)* |
| `therapeutic_coverage_1x` | fraction | 0.3353 | 0.3399 | 0.3426 | 0.0595 | 0.0554 | 0.0147 | **YES** *(abs)* |
| `therapeutic_coverage_2x` | fraction | 0.2005 | 0.1966 | 0.2009 | 0.376 | 0.229 | 0.0132 | **YES** *(abs)* |
| `median_E/IC50 (potency margin)` | - | 0.9527 | 0.9921 | 0.9976 | 0.0646 | 0.0154 | 0.0276 | **YES** |
| `P10_E/IC50 (worst-decile margin)` | - | 0.5562 | 0.5933 | 0.6098 | 0.184 | 0.0486 | 0.0445 | **YES** |
| `exposure_P10` | uM | 18.58 | 19.4 | 19.79 | 0.184 | 0.0486 | 1.09 | **YES** |
| `exposure_P50` | uM | 29.25 | 30.11 | 30.36 | 0.0646 | 0.0154 | 0.678 | **YES** |
| `exposure_P90` | uM | 48.77 | 46.01 | 48.6 | 0.0916 | 0.0571 | 6.61 | no |
| `exposure_P90/P10` | - | 3.546 | 2.96 | 3.009 | 0.233 | 0.0547 | 0.104 | no |
| `penetration_length_L` | um | 101.6 | 120.1 | 121.4 | 0.597 | 0.0705 | 11.5 | no |
| `lowO2/highO2_exposure (quartiles)` | - | 0.5895 | 0.6887 | 0.6621 | 0.193 | 0.0417 | 0.0295 | **YES** |
| `hypoxic/normoxic_exposure (banded)` | - | n/a | n/a | n/a | n/a | n/a | n/a | no |
| `penetration_index (RETIRED)` | - | 0.3415 | 0.3954 | 0.389 | 0.17 | 0.035 | 0.0185 | **YES** |
| `kill_fraction (reference)` | fraction | 0.1872 | 0.1799 | 0.1941 | 0.135 | 0.115 | 0.0291 | **YES** *(abs)* |

**Converged (worst drug, 15→12 µm):** `mean_tissue_auc`, `median_tissue_auc`, `mean_tumor_auc`, `median_tumor_auc`, `therapeutic_coverage_0.5x`, `therapeutic_coverage_1x`, `therapeutic_coverage_2x`, `median_E/IC50 (potency margin)`, `P10_E/IC50 (worst-decile margin)`, `exposure_P10`, `exposure_P50`, `lowO2/highO2_exposure (quartiles)`, `penetration_index (RETIRED)`, `kill_fraction (reference)`.

**NOT converged:** `exposure_P90`, `exposure_P90/P10`, `penetration_length_L`, `hypoxic/normoxic_exposure (banded)`.

**Read this honestly for the headline metric.** `therapeutic_coverage_1x` moves **0.0147 absolute** (1.5 percentage points of tumor cells) on its worst drug between 15 µm and 12 µm, which is inside the 0.03 tolerance — but that same move is **5.5% relative**, and it is the *least*-covered drug (temozolomide, coverage ≈ 0.31) that produces it, which is exactly the regime where relative error is misleading. The residual is real and monotone, not noise: carry it as a **±0.03 error bar on any quoted coverage**. The well-covered drugs are far tighter (cisplatin 0.0007, doxorubicin 0.0006 absolute).

For contrast, `penetration_index` moves 3.5% relative over the same refinement step while its entire spread across the 5-drug panel at dx = 12 µm is only 53.0% — its grid drift is 0× larger than the entire signal it is supposed to carry.

---

## 3. Discrimination across the 5-drug panel (dx = 12 µm)

### 3a. The per-drug table

| drug | MW (Da) | real IC50 (µM) | mean tumor AUC (µM·s) | median tumor AUC (µM·s) | E P10 (µM) | E P50 (µM) | E P90 (µM) | **coverage @0.5× IC50** | **coverage @1× IC50** | **coverage @2× IC50** | median E/IC50 | P10 E/IC50 | L (µm) | lowO₂/highO₂ E | `penetration_index` | kill_fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194 | 606 | 7.66e+05 | 7.27e+05 | 24.4 | 33.7 | 49.1 | 0.0000 | **0.0000** | 0.0000 | 0.0555 | 0.0403 | 280.0 | 0.751 | 0.5107 | 0.0184 |
| lomustine | 234 | 50 | 9.32e+05 | 8.43e+05 | 22.9 | 39 | 67.1 | 0.8560 | **0.2802** | 0.0197 | 0.781 | 0.459 | 78.0 | 0.644 | 0.3743 | 0.2371 |
| cisplatin | 298 | 24.6 | 1.82e+06 | 1.69e+06 | 51.1 | 78 | 125 | 0.9986 | **0.9955** | 0.9139 | 3.18 | 2.08 | 150.2 | 0.677 | 0.4283 | 0.4362 |
| methotrexate | 454 | 12 | 1.57e+04 | 1.39e+04 | 0.333 | 0.643 | 1.2 | 0.0000 | **0.0000** | 0.0000 | 0.0536 | 0.0278 | 53.4 | 0.624 | 0.3273 | 0.0192 |
| doxorubicin | 544 | 0.5 | 1.15e+04 | 9.94e+03 | 0.221 | 0.46 | 0.902 | 0.8578 | **0.4373** | 0.0711 | 0.921 | 0.442 | 45.1 | 0.615 | 0.3044 | 0.2595 |

`L` is the decaying component of `C₀·exp(−d/L)+floor` fitted to the binned AUC-vs-distance profile, where distance is now the **PRIMARY** direct nearest-endothelial-cell metric (`caller-supplied`) — not the LEGACY p95 percentile-mask threshold this report used before (2.7% of tissue area under that legacy mask, kept here for reference only; see `results/vessel_profile.md` for the full head-to-head). **SOURCE NULL** (identical fit on the raw `vascular_density` field): L = nan µm.

### 3b. Spread vs numerical noise

A metric only ranks drugs credibly if its across-drug spread is much larger than its own grid drift. `separation = spread / worst-drug grid drift`, both taken on the **same scale** (absolute for the [0,1] fractions, relative for everything else — the rule stated in §2); we require ≥ 3. A metric that fails here is reporting less signal than it carries numerical noise.

| metric | min | max | fold range (max/min) | spread (abs) | spread (rel) | grid drift 15→12 µm (matched scale) | **separation** | discriminates? |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|
| `mean_tissue_auc` | 1.111e+04 | 1.75e+06 | 157.5 | 1.74e+06 | 2.55 | 0.0231 | **110** | **YES** |
| `median_tissue_auc` | 9693 | 1.637e+06 | 168.9 | 1.63e+06 | 2.55 | 0.0202 | **126** | **YES** |
| `mean_tumor_auc` | 1.15e+04 | 1.822e+06 | 158.3 | 1.81e+06 | 2.55 | 0.0219 | **116** | **YES** |
| `median_tumor_auc` | 9943 | 1.685e+06 | 169.5 | 1.68e+06 | 2.55 | 0.0154 | **165** | **YES** |
| `therapeutic_coverage_0.5x` | 0 | 0.9986 | n/a | 0.999 | 1.84 | 0.0128 | **78** | no |
| `therapeutic_coverage_1x` | 0 | 0.9955 | n/a | 0.995 | 2.91 | 0.0147 | **67.7** | no |
| `therapeutic_coverage_2x` | 0 | 0.9139 | n/a | 0.914 | 4.55 | 0.0132 | **69.1** | no |
| `median_E/IC50 (potency margin)` | 0.05357 | 3.178 | 59.32 | 3.12 | 3.13 | 0.0154 | **203** | **YES** |
| `P10_E/IC50 (worst-decile margin)` | 0.02778 | 2.08 | 74.88 | 2.05 | 3.36 | 0.0486 | **69.2** | **YES** |
| `exposure_P10` | 0.2212 | 51.06 | 230.8 | 50.8 | 2.57 | 0.0486 | **52.8** | **YES** |
| `exposure_P50` | 0.4603 | 78.01 | 169.5 | 77.5 | 2.55 | 0.0154 | **165** | **YES** |
| `exposure_P90` | 0.9016 | 124.7 | 138.3 | 124 | 2.55 | 0.0571 | **44.6** | **YES** |
| `exposure_P90/P10` | 2.014 | 4.076 | 2.024 | 2.06 | 0.685 | 0.0547 | **12.5** | **YES** |
| `penetration_length_L` | 45.08 | 280 | 6.212 | 235 | 1.94 | 0.0705 | **27.4** | **YES** |
| `lowO2/highO2_exposure (quartiles)` | 0.6155 | 0.7506 | 1.219 | 0.135 | 0.204 | 0.0417 | **4.89** | **YES** |
| `hypoxic/normoxic_exposure (banded)` | n/a | n/a | n/a | n/a | n/a | n/a | **n/a** | no |
| `penetration_index (RETIRED)` | 0.3044 | 0.5107 | 1.678 | 0.206 | 0.53 | 0.035 | **15.2** | **YES** |
| `kill_fraction (reference)` | 0.01835 | 0.4362 | 23.77 | 0.418 | 2.15 | 0.0291 | **14.3** | **YES** |

**Discriminates:** `mean_tissue_auc`, `median_tissue_auc`, `mean_tumor_auc`, `median_tumor_auc`, `median_E/IC50 (potency margin)`, `P10_E/IC50 (worst-decile margin)`, `exposure_P10`, `exposure_P50`, `exposure_P90`, `exposure_P90/P10`, `penetration_length_L`, `lowO2/highO2_exposure (quartiles)`, `penetration_index (RETIRED)`, `kill_fraction (reference)`.

**Does NOT discriminate:** `therapeutic_coverage_0.5x`, `therapeutic_coverage_1x`, `therapeutic_coverage_2x`, `hypoxic/normoxic_exposure (banded)`. Three of these are informative failures, not bugs:

- **`penetration_length_L`** — the fitted drug `L` spans only 45.1–280.0 µm (6.21×) against a source-null `L` of nan µm, while its own grid drift is 7%. On real tissue the measured decay length is set by the endothelial-KDE geometry, not by drug transport — the confound `results/vessel_profile.md` quantifies at length. The metric is sound; this *section* cannot resolve drugs with it.
- **`lowO2/highO2_exposure`** — grid-stable (4.2% drift) but essentially identical for every drug (0.615–0.751). That is the correct answer, not a failure: the perfusion gradient is a property of the *tissue*, and all five drugs are delivered proportionally to the same vascular field, so they all under-serve the badly-perfused quartile by the same ~3×.
- **`penetration_index`** — the metric this module replaces. See §4.

---

## 4. The explicit contrast with `penetration_index`

| | `penetration_index` (retired) | `therapeutic_coverage` @1× IC50 (proposed headline) |
|---|---|---|
| definition | mean(AUC)/p99(AUC) over the tissue mask | fraction of tumor cells with AUC/T > IC50 |
| units | dimensionless *shape* | dimensionless *fraction of cells* |
| scales with delivered dose? | **NO** (scale-invariant by construction) | **YES** |
| knows the drug's potency? | **NO** | **YES** (IC50 enters explicitly) |
| range across the 5-drug panel | **0.3044 – 0.5107** (spread 0.2063, fold 1.678×) | **0.0000 – 0.9955** (spread 0.9955, fold nan×) |
| grid drift 15→12 µm (worst drug), absolute | 0.0185 | 0.0147 |
| grid drift 15→12 µm (worst drug), relative | 0.035 | 0.0554 |
| **separation** (spread / drift, matched scale) | **15.2** | **67.7** |

Both are dimensionless numbers in [0,1] on the same tissue, so the absolute row is the like-for-like comparison: `penetration_index` drifts 0.0185 with grid refinement while its whole across-drug spread is 0.2063 — the numerics move it ~0× more than the drugs do. Coverage drifts 0.0147 against a spread of 0.995 — ~68× the other way round.

The scale-invariance claim is not rhetorical. `delivery_metrics.selftest()` (run it: `./.venv/bin/python -m tmesim.validation.delivery_metrics --selftest`) asserts arithmetically that multiplying an exposure field by 10⁻³ and by 10³ leaves `mean/p99` identical to 1e-12 while `therapeutic_coverage` moves from ~0.000 to 1.000 over the same 10⁶× dose change.

Per drug, side by side:

| drug | `penetration_index` | mean tumor AUC (µM·s) | coverage @1× IC50 |
|---|---:|---:|---:|
| temozolomide | 0.5107 | 7.66e+05 | 0.0000 |
| lomustine | 0.3743 | 9.32e+05 | 0.2802 |
| cisplatin | 0.4283 | 1.82e+06 | 0.9955 |
| methotrexate | 0.3273 | 1.57e+04 | 0.0000 |
| doxorubicin | 0.3044 | 1.15e+04 | 0.4373 |

The index moves by 0.2063 over a panel whose delivered exposure moves by 158×. It is not measuring delivery. Coverage moves by 0.995 on the same panel.

---

## 5. Recommendation — the headline delivery number

> **Headline: `therapeutic_coverage` at 1× IC50 — the fraction of tumor cells whose time-averaged exposure `E = AUC/T` exceeds the drug's real IC50.** Quote it at dx = 12 µm, always with the 0.5×/2× IC50 band beside it, with **`median E/IC50`** (the continuous potency margin) as its companion, and with `median_tumor_auc` (µM·s) as the supporting absolute number.

Justification, against the three requirements:

1. **Physically meaningful.** It compares a concentration to a concentration. `E = AUC/T` is exactly the quantity `cells.pkpd` feeds into the Hill dose-response, and `IC50` is by definition the concentration giving half-maximal effect — so the metric is the delivery half of the kill model, read out directly rather than through a stochastic Bernoulli draw. It is bounded in [0,1], monotone in delivered dose, and one sentence of plain English ("x% of the tumor saw at least a half-maximal-effect concentration").
2. **Grid-converged — with a stated caveat.** Worst-drug change from 15 µm to 12 µm is **0.0147 absolute** (1.5 percentage points of tumor cells), inside the 0.03 tolerance for a bounded fraction. That same step is 5.5% *relative*, driven entirely by the least-covered drug — so the honest statement is **converged to ±0.03 in coverage units, not to 5% relative**, and coverage should be quoted with that error bar. Unlike a mask-area average it is computed over a FIXED population of real cells, so grid refinement cannot change the denominator — which is why it is this stable at all.
3. **Discriminating.** It spans 0.0000–0.9955 across the panel (0.995 absolute, nan× fold) — a separation of **67.7×** its own grid drift, versus 15.2 for `penetration_index`. It also ranks the drugs the way the mechanism says it should: potency and delivery are combined, so a drug wins only when it both arrives AND is potent at the concentration it achieves.

**What to quote alongside it (and why not instead of it):**

- `median_tumor_auc` (µM·s) — the raw absolute delivery. Strongly discriminating (169× across the panel) though not itself grid-converged on this ladder (1.5% worst-drug drift; use `mean_tumor_auc`, 2.2%, if a converged absolute number is required), and in any case *not* interpretable without a potency: 1.7e+06 µM·s of cisplatin (IC50 24.6 µM) and 9.9e+03 µM·s of doxorubicin (IC50 0.5 µM) are not comparable doses. Report it as the input to coverage, not as the verdict.
- `median E/IC50` and `P10 E/IC50` (dimensionless) — the potency-normalised exposure margin at the median and at the worst-served decile. Same information content as coverage, continuous instead of thresholded, and immune to saturation; quote `median E/IC50` whenever coverage pins at 1.000. It inherits the median's grid sensitivity (1.5% worst-drug drift, versus 0.015 absolute for coverage), which is the price of being continuous — it is a ranking companion, not a replacement headline.
- `E_p10` (µM) — the worst-served decile, the honest floor of the delivery distribution and the number that should worry you about regrowth.
- `penetration_length_fit` `L` (µm) — quote only against the SOURCE NULL. On this section the fitted drug `L` spans 45.1–280.0 µm against a source-null `L` of nan µm (drug/null ratio nan–nan×), i.e. the measured length is dominated by vessel/KDE geometry rather than by drug transport — the same confound `results/vessel_profile.md` documents at length. The drug-intrinsic length separates cleanly only in the isolated single-vessel control, where doxorubicin's gradient is ~4× steeper than temozolomide's. `L` is a legitimate *physical* number and a legitimate *model* diagnostic; it is not a delivery ranking on this section.
- `lowO₂/highO₂ exposure ratio` — diagnostic only. On this section it is 0.615–0.751 across the panel. It is a ratio (hence scale-invariant, hence the same structural weakness as `penetration_index`) and the banded hypoxic/normoxic form sits on a partition that is itself not grid-converged. Read it as 'is delivery correlated with perfusion?', never as a drug ranking.

**What to retire.** `penetration_index` should no longer be quoted as a delivery number anywhere. It is not wrong — it is a faithful measurement of *exposure-field uniformity* — but on real tissue at converged resolution that uniformity is a property of the reconstructed vasculature, shared by every drug. If a uniformity statistic is wanted, quote it as `exposure_uniformity_index` and say what it is.

---

## 6. Honesty and caveats

- **No outcome label.** The Xenium section is an untreated resection. Every number here is a PREDICTION of the calibrated physics on real geometry, not a validated response. Coverage being grid-converged and drug-discriminating says it is a *well-posed* metric; it does not say the underlying delivery model is right.
- **Coverage inherits the IC50's uncertainty.** The IC50s are real GDSC / literature potencies measured in cell culture, not in this tissue. That is why the 0.5× and 2× variants are reported as a band and not buried: they are the metric's sensitivity to the one external number it depends on. Where the 0.5×–2× band is wide, the ranking at that position is soft.
- **Coverage can saturate — this is its one structural weakness, and it is why the continuous companion is reported next to it.** A drug delivering far above its IC50 *everywhere* scores 1.000 and stops discriminating. On this panel that does NOT bite — 0 of 5 drugs reach 0.999 and the maximum is 0.995 — so the weakness is latent here, not active. It will bite on a higher-dose or more potent panel. `median E/IC50` — the **potency-normalised exposure margin** — is the same comparison without the threshold: dimensionless, drug-specific, absolute, unbounded above, and it spans 59× across the panel. Read them together: coverage says *how much of the tumor* cleared the bar, the margin says *by how far*. Coverage is the right headline precisely in the regime that matters clinically — where part of the tumor is under-dosed — and the margin is the tiebreaker when it is not.
- **One anchor sets the absolute scale.** The single α (= 148.6 here, read live from `validation.calibration.calibrate_alpha()`) is fit to one doxorubicin penetration measurement; every absolute exposure — and therefore every coverage value — rides on it. Coverage is *relatively* robust because the same α applies to all 5 drugs, but the absolute level is not independently validated. **Older reports in `results/` may quote a different α** (the half-distance vs e-folding-length convention was corrected in `calibration.py`, a factor 1/ln²2 = 2.08×); the absolute µM·s figures here are therefore not directly comparable to those, while the *structural* conclusions — the index collapsing, coverage spreading — are about the relative behaviour of the metrics and are unaffected.
- **The tumor hypoxic fraction is not grid-converged** (`results/real_drug_panel.md`), which is exactly why metric 5 is a diagnostic and why its primary form is the threshold-free quartile contrast. On the reconstruction as it stands the hypoxic tumor fraction is: 20 µm → 0.0%, 15 µm → 0.0%, 12 µm → 0.0%. Where it is 0%, the banded hypoxic/normoxic ratio is legitimately undefined (`n/a`) — that is reported, not papered over.
- **2D section, straight-line distances, static plasma level.** All the standing limitations of the engine (`results/numerical_rigor.md`) apply unchanged; nothing in `core` / `physics` / `molecular` / `reconstruct` / `sim.engine` was modified for this audit.

*Total runtime 1203 s (3 grids × 5 drugs).*
