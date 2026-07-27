# Defensible delivery metrics — replacing `penetration_index`

**Module:** `tmesim/validation/delivery_metrics.py` (run: `./.venv/bin/python -m tmesim.validation.delivery_metrics`)

**Metric self-test:** `./.venv/bin/python -m tmesim.validation.delivery_metrics --selftest` (known-answer checks on the definitions, no simulation)

**Figure:** `results/delivery_metrics.png`

**Geometry:** the REAL open 10x Xenium GBM section — 40,887 measured cells (tumor 30,174), reconstructed by `reconstruct.real_pipeline.reconstruct_real('xenium', dx_um=…)`; real endothelial-KDE vascular field, real matrix-transcript ECM, physics-solved oxygen. **Drugs:** the 5-drug panel of [`results/real_drug_panel.md`](real_drug_panel.md), each compiled through the identical pipeline (`calibrate_alpha()` → α = 148.6, fit ONCE on doxorubicin → `calibrated_drugparams` → `apply_barrier_to_drugparams`). **Engine:** advection on, dt = 10 s, 6 h window, seed 0 — unchanged from the flagship run. Production resolution **dx = 12 µm**, grid (424, 675) (dx = 10 µm violates the advection CFL condition at dt = 10 s and is therefore not on the ladder).

---

## 0. The problem being fixed

`penetration_index = mean(AUC) / p99(AUC)` over the tissue mask is a **scale-invariant shape statistic**: multiply the entire exposure field by any constant and it does not move. At the grid-converged resolution it therefore reports the *vascular geometry of the section*, which is identical for all 5 drugs — and it duly collapses:

> at dx = 12 µm, `penetration_index` = **0.195–0.196** across the whole panel (spread 0.0014), while those same 5 drugs' mean delivered exposure differs by **120×**.

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

1. **Absolute exposure (µM·s).** The integral the pharmacodynamics actually consumes. The median is reported next to the mean because the field is strongly right-skewed — a thin perivascular cuff carries much of the integral, so mean/median ≈ 1.4× on this section and the mean alone overstates what a typical tumor cell sees. Honest limitation: an exposure in µM·s is not interpretable without a potency, which is why it is not the headline.
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
| `mean_tissue_auc` | uM*s | 4.452e+06 | 4.672e+06 | 4.692e+06 | 0.0494 | 0.00412 | 5.18e+04 | **YES** |
| `median_tissue_auc` | uM*s | 3.953e+06 | 3.702e+06 | 3.477e+06 | 0.0749 | 0.067 | 6.7e+05 | no |
| `mean_tumor_auc` | uM*s | 4.401e+06 | 4.437e+06 | 4.317e+06 | 0.00931 | 0.0274 | 3.21e+05 | **YES** |
| `median_tumor_auc` | uM*s | 3.898e+06 | 3.449e+06 | 3.19e+06 | 0.128 | 0.0812 | 7.58e+05 | no |
| `therapeutic_coverage_0.5x` | fraction | 0.8816 | 0.8286 | 0.8113 | 0.146 | 0.0492 | 0.0318 | no *(abs)* |
| `therapeutic_coverage_1x` | fraction | 0.7127 | 0.6766 | 0.6567 | 0.0888 | 0.0828 | 0.0281 | **YES** *(abs)* |
| `therapeutic_coverage_2x` | fraction | 0.5251 | 0.5007 | 0.4876 | 0.628 | 0.049 | 0.0276 | **YES** *(abs)* |
| `median_E/IC50 (potency margin)` | - | 4.509 | 4.095 | 3.843 | 0.128 | 0.0812 | 0.684 | no |
| `P10_E/IC50 (worst-decile margin)` | - | 1.288 | 0.715 | 0.6131 | 0.515 | 0.153 | 0.258 | no |
| `exposure_P10` | uM | 54.93 | 27.8 | 23.59 | 0.515 | 0.153 | 11.4 | no |
| `exposure_P50` | uM | 180.5 | 159.7 | 147.7 | 0.128 | 0.0812 | 35.1 | no |
| `exposure_P90` | uM | 376.4 | 417.9 | 416.4 | 0.123 | 0.0125 | 3.47 | **YES** |
| `exposure_P90/P10` | - | 7.765 | 15.19 | 17.62 | 1.31 | 0.177 | 2.65 | no |
| `penetration_length_L` | um | 51.6 | 34.51 | 27.27 | 0.43 | 0.252 | 9.49 | no |
| `lowO2/highO2_exposure (quartiles)` | - | 0.3799 | 0.3319 | 0.3358 | 0.202 | 0.0129 | 0.00428 | **YES** |
| `hypoxic/normoxic_exposure (banded)` | - | n/a | 0.0001219 | 0.4673 | n/a | 3.83e+03 | 0.468 | no |
| `penetration_index (RETIRED)` | - | 0.2758 | 0.2022 | 0.1951 | 0.365 | 0.0859 | 0.0183 | no |
| `kill_fraction (reference)` | fraction | 0.5886 | 0.341 | 0.3019 | 0.429 | 0.139 | 0.0483 | no *(abs)* |

**Converged (worst drug, 15→12 µm):** `mean_tissue_auc`, `mean_tumor_auc`, `therapeutic_coverage_1x`, `therapeutic_coverage_2x`, `exposure_P90`, `lowO2/highO2_exposure (quartiles)`.

**NOT converged:** `median_tissue_auc`, `median_tumor_auc`, `therapeutic_coverage_0.5x`, `median_E/IC50 (potency margin)`, `P10_E/IC50 (worst-decile margin)`, `exposure_P10`, `exposure_P50`, `exposure_P90/P10`, `penetration_length_L`, `hypoxic/normoxic_exposure (banded)`, `penetration_index (RETIRED)`, `kill_fraction (reference)`.

**Read this honestly for the headline metric.** `therapeutic_coverage_1x` moves **0.0281 absolute** (2.8 percentage points of tumor cells) on its worst drug between 15 µm and 12 µm, which is inside the 0.03 tolerance — but that same move is **8.3% relative**, and it is the *least*-covered drug (temozolomide, coverage ≈ 0.31) that produces it, which is exactly the regime where relative error is misleading. The residual is real and monotone, not noise: carry it as a **±0.03 error bar on any quoted coverage**. The well-covered drugs are far tighter (cisplatin 0.0122, doxorubicin 0.0113 absolute).

For contrast, `penetration_index` moves 8.6% relative over the same refinement step while its entire spread across the 5-drug panel at dx = 12 µm is only 0.7% — its grid drift is 12× larger than the entire signal it is supposed to carry.

---

## 3. Discrimination across the 5-drug panel (dx = 12 µm)

### 3a. The per-drug table

| drug | MW (Da) | real IC50 (µM) | mean tumor AUC (µM·s) | median tumor AUC (µM·s) | E P10 (µM) | E P50 (µM) | E P90 (µM) | **coverage @0.5× IC50** | **coverage @1× IC50** | **coverage @2× IC50** | median E/IC50 | P10 E/IC50 | L (µm) | lowO₂/highO₂ E | `penetration_index` | kill_fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194 | 606 | 1.16e+07 | 8.57e+06 | 63.4 | 397 | 1.12e+03 | 0.6145 | **0.3108** | 0.0827 | 0.655 | 0.105 | 28.2 | 0.335 | 0.1949 | 0.1846 |
| lomustine | 234 | 50 | 3.17e+06 | 2.34e+06 | 17.3 | 108 | 305 | 0.8644 | **0.7557** | 0.5365 | 2.17 | 0.346 | 27.2 | 0.336 | 0.1953 | 0.3099 |
| cisplatin | 298 | 24.6 | 6.44e+06 | 4.76e+06 | 35.2 | 220 | 621 | 0.9571 | **0.9248** | 0.8684 | 8.98 | 1.43 | 27.7 | 0.335 | 0.1959 | 0.4142 |
| methotrexate | 454 | 12 | 2.73e+05 | 2.02e+05 | 1.49 | 9.36 | 26.3 | 0.6739 | **0.3876** | 0.1211 | 0.78 | 0.124 | 26.7 | 0.336 | 0.1948 | 0.2026 |
| doxorubicin | 544 | 0.5 | 9.68e+04 | 7.16e+04 | 0.528 | 3.32 | 9.3 | 0.9463 | **0.9047** | 0.8293 | 6.63 | 1.06 | 26.5 | 0.337 | 0.1946 | 0.3983 |

`L` is the decaying component of `C₀·exp(−d/L)+floor` fitted to the binned AUC-vs-distance profile at the p95 vessel threshold (5.0% of tissue area). **SOURCE NULL** (identical fit on the raw `vascular_density` field): L = 23.7 µm.

### 3b. Spread vs numerical noise

A metric only ranks drugs credibly if its across-drug spread is much larger than its own grid drift. `separation = spread / worst-drug grid drift`, both taken on the **same scale** (absolute for the [0,1] fractions, relative for everything else — the rule stated in §2); we require ≥ 3. A metric that fails here is reporting less signal than it carries numerical noise.

| metric | min | max | fold range (max/min) | spread (abs) | spread (rel) | grid drift 15→12 µm (matched scale) | **separation** | discriminates? |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|
| `mean_tissue_auc` | 1.055e+05 | 1.261e+07 | 119.5 | 1.25e+07 | 2.66 | 0.00412 | **646** | **YES** |
| `median_tissue_auc` | 7.823e+04 | 9.341e+06 | 119.4 | 9.26e+06 | 2.66 | 0.067 | **39.8** | **YES** |
| `mean_tumor_auc` | 9.682e+04 | 1.161e+07 | 119.9 | 1.15e+07 | 2.67 | 0.0274 | **97.4** | **YES** |
| `median_tumor_auc` | 7.163e+04 | 8.575e+06 | 119.7 | 8.5e+06 | 2.67 | 0.0812 | **32.8** | **YES** |
| `therapeutic_coverage_0.5x` | 0.6145 | 0.9571 | 1.558 | 0.343 | 0.422 | 0.0318 | **10.8** | **YES** |
| `therapeutic_coverage_1x` | 0.3108 | 0.9248 | 2.976 | 0.614 | 0.935 | 0.0281 | **21.9** | **YES** |
| `therapeutic_coverage_2x` | 0.08272 | 0.8684 | 10.5 | 0.786 | 1.61 | 0.0276 | **28.4** | **YES** |
| `median_E/IC50 (potency margin)` | 0.655 | 8.978 | 13.71 | 8.32 | 2.17 | 0.0812 | **26.7** | **YES** |
| `P10_E/IC50 (worst-decile margin)` | 0.1047 | 1.433 | 13.69 | 1.33 | 2.17 | 0.153 | **14.2** | **YES** |
| `exposure_P10` | 0.5283 | 63.44 | 120.1 | 62.9 | 2.67 | 0.153 | **17.5** | **YES** |
| `exposure_P50` | 3.316 | 397 | 119.7 | 394 | 2.67 | 0.0812 | **32.8** | **YES** |
| `exposure_P90` | 9.302 | 1121 | 120.5 | 1.11e+03 | 2.67 | 0.0125 | **214** | **YES** |
| `exposure_P90/P10` | 17.6 | 17.67 | 1.004 | 0.0704 | 0.00399 | 0.177 | **0.0226** | no |
| `penetration_length_L` | 26.5 | 28.25 | 1.066 | 1.74 | 0.064 | 0.252 | **0.254** | no |
| `lowO2/highO2_exposure (quartiles)` | 0.3345 | 0.3366 | 1.006 | 0.00212 | 0.00632 | 0.0129 | **0.49** | no |
| `hypoxic/normoxic_exposure (banded)` | 0.4663 | 0.468 | 1.004 | 0.00176 | 0.00377 | 3.83e+03 | **9.83e-07** | no |
| `penetration_index (RETIRED)` | 0.1946 | 0.1959 | 1.007 | 0.00138 | 0.00708 | 0.0859 | **0.0824** | no |
| `kill_fraction (reference)` | 0.1846 | 0.4142 | 2.244 | 0.23 | 0.76 | 0.0483 | **4.76** | **YES** |

**Discriminates:** `mean_tissue_auc`, `median_tissue_auc`, `mean_tumor_auc`, `median_tumor_auc`, `therapeutic_coverage_0.5x`, `therapeutic_coverage_1x`, `therapeutic_coverage_2x`, `median_E/IC50 (potency margin)`, `P10_E/IC50 (worst-decile margin)`, `exposure_P10`, `exposure_P50`, `exposure_P90`, `kill_fraction (reference)`.

**Does NOT discriminate:** `exposure_P90/P10`, `penetration_length_L`, `lowO2/highO2_exposure (quartiles)`, `hypoxic/normoxic_exposure (banded)`, `penetration_index (RETIRED)`. Three of these are informative failures, not bugs:

- **`penetration_length_L`** — the fitted drug `L` spans only 26.5–28.2 µm (1.07×) against a source-null `L` of 23.7 µm, while its own grid drift is 25%. On real tissue the measured decay length is set by the endothelial-KDE geometry, not by drug transport — the confound `results/vessel_profile.md` quantifies at length. The metric is sound; this *section* cannot resolve drugs with it.
- **`lowO2/highO2_exposure`** — grid-stable (1.3% drift) but essentially identical for every drug (0.335–0.337). That is the correct answer, not a failure: the perfusion gradient is a property of the *tissue*, and all five drugs are delivered proportionally to the same vascular field, so they all under-serve the badly-perfused quartile by the same ~3×.
- **`penetration_index`** — the metric this module replaces. See §4.

---

## 4. The explicit contrast with `penetration_index`

| | `penetration_index` (retired) | `therapeutic_coverage` @1× IC50 (proposed headline) |
|---|---|---|
| definition | mean(AUC)/p99(AUC) over the tissue mask | fraction of tumor cells with AUC/T > IC50 |
| units | dimensionless *shape* | dimensionless *fraction of cells* |
| scales with delivered dose? | **NO** (scale-invariant by construction) | **YES** |
| knows the drug's potency? | **NO** | **YES** (IC50 enters explicitly) |
| range across the 5-drug panel | **0.1946 – 0.1959** (spread 0.0014, fold 1.007×) | **0.3108 – 0.9248** (spread 0.6140, fold 2.98×) |
| grid drift 15→12 µm (worst drug), absolute | 0.0183 | 0.0281 |
| grid drift 15→12 µm (worst drug), relative | 0.0859 | 0.0828 |
| **separation** (spread / drift, matched scale) | **0.0824** | **21.9** |

Both are dimensionless numbers in [0,1] on the same tissue, so the absolute row is the like-for-like comparison: `penetration_index` drifts 0.0183 with grid refinement while its whole across-drug spread is 0.0014 — the numerics move it ~13× more than the drugs do. Coverage drifts 0.0281 against a spread of 0.614 — ~22× the other way round.

The scale-invariance claim is not rhetorical. `delivery_metrics.selftest()` (run it: `./.venv/bin/python -m tmesim.validation.delivery_metrics --selftest`) asserts arithmetically that multiplying an exposure field by 10⁻³ and by 10³ leaves `mean/p99` identical to 1e-12 while `therapeutic_coverage` moves from ~0.000 to 1.000 over the same 10⁶× dose change.

Per drug, side by side:

| drug | `penetration_index` | mean tumor AUC (µM·s) | coverage @1× IC50 |
|---|---:|---:|---:|
| temozolomide | 0.1949 | 1.16e+07 | 0.3108 |
| lomustine | 0.1953 | 3.17e+06 | 0.7557 |
| cisplatin | 0.1959 | 6.44e+06 | 0.9248 |
| methotrexate | 0.1948 | 2.73e+05 | 0.3876 |
| doxorubicin | 0.1946 | 9.68e+04 | 0.9047 |

The index moves by 0.0014 over a panel whose delivered exposure moves by 120×. It is not measuring delivery. Coverage moves by 0.614 on the same panel.

---

## 5. Recommendation — the headline delivery number

> **Headline: `therapeutic_coverage` at 1× IC50 — the fraction of tumor cells whose time-averaged exposure `E = AUC/T` exceeds the drug's real IC50.** Quote it at dx = 12 µm, always with the 0.5×/2× IC50 band beside it, with **`median E/IC50`** (the continuous potency margin) as its companion, and with `median_tumor_auc` (µM·s) as the supporting absolute number.

Justification, against the three requirements:

1. **Physically meaningful.** It compares a concentration to a concentration. `E = AUC/T` is exactly the quantity `cells.pkpd` feeds into the Hill dose-response, and `IC50` is by definition the concentration giving half-maximal effect — so the metric is the delivery half of the kill model, read out directly rather than through a stochastic Bernoulli draw. It is bounded in [0,1], monotone in delivered dose, and one sentence of plain English ("x% of the tumor saw at least a half-maximal-effect concentration").
2. **Grid-converged — with a stated caveat.** Worst-drug change from 15 µm to 12 µm is **0.0281 absolute** (2.8 percentage points of tumor cells), inside the 0.03 tolerance for a bounded fraction. That same step is 8.3% *relative*, driven entirely by the least-covered drug — so the honest statement is **converged to ±0.03 in coverage units, not to 5% relative**, and coverage should be quoted with that error bar. Unlike a mask-area average it is computed over a FIXED population of real cells, so grid refinement cannot change the denominator — which is why it is this stable at all.
3. **Discriminating.** It spans 0.3108–0.9248 across the panel (0.614 absolute, 2.98× fold) — a separation of **21.9×** its own grid drift, versus 0.0824 for `penetration_index`. It also ranks the drugs the way the mechanism says it should: potency and delivery are combined, so a drug wins only when it both arrives AND is potent at the concentration it achieves.

   **And it recovers the kill ordering exactly.** Ranked by coverage the panel is **cisplatin > doxorubicin > lomustine > methotrexate > temozolomide** — the *identical* order the simulated `kill_fraction` produces on the same runs (kill is the engine's own pharmacodynamic endpoint; on this ladder it is itself only grid-stable to 0.048 absolute, so the agreement is between two independently-computed orderings, not between a number and its own definition). Coverage is a pure *delivery* statistic: it never touches the Hill curve, the per-cell sensitivities, the hypoxia/quiescence resistance factors or the stochastic Bernoulli draw. Reproducing the kill ranking from delivery and potency alone is the strongest evidence on offer that it is measuring the right thing — and `penetration_index`, on the same panel, ranks the drugs **cisplatin > lomustine > temozolomide > methotrexate > doxorubicin**, which is noise.

**What to quote alongside it (and why not instead of it):**

- `median_tumor_auc` (µM·s) — the raw absolute delivery. Strongly discriminating (120× across the panel) though not itself grid-converged on this ladder (8.1% worst-drug drift; use `mean_tumor_auc`, 2.7%, if a converged absolute number is required), and in any case *not* interpretable without a potency: 8.6e+06 µM·s of temozolomide (IC50 606 µM) and 7.2e+04 µM·s of doxorubicin (IC50 0.5 µM) are not comparable doses. Report it as the input to coverage, not as the verdict.
- `median E/IC50` and `P10 E/IC50` (dimensionless) — the potency-normalised exposure margin at the median and at the worst-served decile. Same information content as coverage, continuous instead of thresholded, and immune to saturation; quote `median E/IC50` whenever coverage pins at 1.000. It inherits the median's grid sensitivity (8.1% worst-drug drift, versus 0.028 absolute for coverage), which is the price of being continuous — it is a ranking companion, not a replacement headline.
- `E_p10` (µM) — the worst-served decile, the honest floor of the delivery distribution and the number that should worry you about regrowth.
- `penetration_length_fit` `L` (µm) — quote only against the SOURCE NULL. On this section the fitted drug `L` spans 26.5–28.2 µm against a source-null `L` of 23.7 µm (drug/null ratio 1.12–1.19×), i.e. the measured length is dominated by vessel/KDE geometry rather than by drug transport — the same confound `results/vessel_profile.md` documents at length. The drug-intrinsic length separates cleanly only in the isolated single-vessel control, where doxorubicin's gradient is ~4× steeper than temozolomide's. `L` is a legitimate *physical* number and a legitimate *model* diagnostic; it is not a delivery ranking on this section.
- `lowO₂/highO₂ exposure ratio` — diagnostic only. On this section it is 0.335–0.337 across the panel. It is a ratio (hence scale-invariant, hence the same structural weakness as `penetration_index`) and the banded hypoxic/normoxic form sits on a partition that is itself not grid-converged. Read it as 'is delivery correlated with perfusion?', never as a drug ranking.

**What to retire.** `penetration_index` should no longer be quoted as a delivery number anywhere. It is not wrong — it is a faithful measurement of *exposure-field uniformity* — but on real tissue at converged resolution that uniformity is a property of the reconstructed vasculature, shared by every drug. If a uniformity statistic is wanted, quote it as `exposure_uniformity_index` and say what it is.

---

## 6. Honesty and caveats

- **No outcome label.** The Xenium section is an untreated resection. Every number here is a PREDICTION of the calibrated physics on real geometry, not a validated response. Coverage being grid-converged and drug-discriminating says it is a *well-posed* metric; it does not say the underlying delivery model is right.
- **Coverage inherits the IC50's uncertainty.** The IC50s are real GDSC / literature potencies measured in cell culture, not in this tissue. That is why the 0.5× and 2× variants are reported as a band and not buried: they are the metric's sensitivity to the one external number it depends on. Where the 0.5×–2× band is wide, the ranking at that position is soft.
- **Coverage can saturate — this is its one structural weakness, and it is why the continuous companion is reported next to it.** A drug delivering far above its IC50 *everywhere* scores 1.000 and stops discriminating. On this panel that does NOT bite — 0 of 5 drugs reach 0.999 and the maximum is 0.925 — so the weakness is latent here, not active. It will bite on a higher-dose or more potent panel. `median E/IC50` — the **potency-normalised exposure margin** — is the same comparison without the threshold: dimensionless, drug-specific, absolute, unbounded above, and it spans 14× across the panel. Read them together: coverage says *how much of the tumor* cleared the bar, the margin says *by how far*. Coverage is the right headline precisely in the regime that matters clinically — where part of the tumor is under-dosed — and the margin is the tiebreaker when it is not.
- **One anchor sets the absolute scale.** The single α (= 148.6 here, read live from `validation.calibration.calibrate_alpha()`) is fit to one doxorubicin penetration measurement; every absolute exposure — and therefore every coverage value — rides on it. Coverage is *relatively* robust because the same α applies to all 5 drugs, but the absolute level is not independently validated. **Older reports in `results/` may quote a different α** (the half-distance vs e-folding-length convention was corrected in `calibration.py`, a factor 1/ln²2 = 2.08×); the absolute µM·s figures here are therefore not directly comparable to those, while the *structural* conclusions — the index collapsing, coverage spreading — are about the relative behaviour of the metrics and are unaffected.
- **The tumor hypoxic fraction is not grid-converged** (`results/real_drug_panel.md`), which is exactly why metric 5 is a diagnostic and why its primary form is the threshold-free quartile contrast. On the reconstruction as it stands the hypoxic tumor fraction is: 20 µm → 0.0%, 15 µm → 0.0%, 12 µm → 13.4%. Where it is 0%, the banded hypoxic/normoxic ratio is legitimately undefined (`n/a`) — that is reported, not papered over.
- **2D section, straight-line distances, static plasma level.** All the standing limitations of the engine (`results/numerical_rigor.md`) apply unchanged; nothing in `core` / `physics` / `molecular` / `reconstruct` / `sim.engine` was modified for this audit.

*Total runtime 661 s (3 grids × 5 drugs).*
