# Does the P90/P10 delivery spread separate compounds above its own noise floor? — YES

**Status: RUN.** Driver: `experiments/delivery_spread_noise.py`. Data:
`results/delivery_spread_noise.json`.

**Nothing was changed.** No default, constant or shipped parameter.
`calibrate_alpha()` asserted live at `22.223477711172563`. Drugs compiled through
the flagship path the shipped real-tissue panel uses (`calibrated_drugparams` →
`apply_barrier_to_drugparams`). **Nothing was tuned to make the compounds
separate.**

---

## 0. Verdict up front

**The separation is real and it clears the noise floor by a wide margin.**

| | shipped config | claimed | **span over the full nuisance grid** |
|---|---:|---:|---|
| temozolomide | **2.111** | 1.99 | **1.969 – 2.171** |
| doxorubicin | **3.942** | 3.80 | **3.770 – 3.996** |

**DISJOINT. Gap = 1.599**, against compound spans of 0.202 and 0.226 — **the gap
between the two compounds is ~7× wider than either compound's own spread across
every nuisance setting tested.**

The sharpest way to state it: **in the worst-case pairing anywhere on the grid**
— temozolomide at its highest (2.171, `dx16_repair2`, 3 h) against doxorubicin at
its lowest (3.770, `dx20_native`, 12 h) — **doxorubicin's delivery spread still
exceeds temozolomide's by 1.74×.** There is no combination of grid, mask repair
or exposure window in the tested set that inverts or even closes the ordering.

**This is the opposite result from the penetration-length readout.**
`results/kyle2004_penetration.md` §8 found that the map's *half-distance* readout
could not separate doxorubicin from daunorubicin once the method was allowed to
move (39.4–63.1 vs 32.7–53.9 µm, overlapping). P90/P10 separates temozolomide
from doxorubicin with a 7× margin over the same class of perturbations. **The
delivery half of the panel is not resting on a difference the pipeline cannot
resolve.**

Both headline values also **reproduce**: the claimed 1.99 / 3.80 pair is
recovered at `dx12_native` (1.994 / 3.797), and the repaired-mask flagship
configuration gives 2.111 / 3.942.

---

## 1. Which nuisance axes actually enter this metric

This matters, and it is the reason the sweep is narrower than §8's rather than
weaker. `E_p90_over_p10` is a percentile ratio over per-tumour-cell exposure. Of
the four axes §8 varied:

| axis | enters P90/P10? | why |
|---|:--:|---|
| distance-bin width | **no** | there is no distance binning in this metric |
| distance metric (direct-cell vs percentile-mask) | **no** | no distance map is used |
| exposure vs AUC | **degenerate** | `exposure_uM = auc_uM_s / T` with `T` one scalar per run, and percentiles are homogeneous of degree one |
| grid resolution / tissue-mask repair | **yes** | through the field solve and through which grid cell each measured tumour cell samples |
| exposure window | **yes** | through the shape of the finite-window AUC field |

The degeneracy is **measured, not argued** —
`--stage degenerate` prints `E_p90/E_p10 = 1.9935720609` and
`auc_p90/p10 = 1.9935720609`, relative difference **0.000e+00**.

**So P90/P10 is structurally exposed to fewer reduction choices than the
half-distance readout was.** That is a real property of the metric and it is a
large part of why it survives where the other did not: it is a statistic of a
fixed cell population, not of a binned spatial reduction.

`delivery_metrics.tumor_cell_exposure` states the design intent — *"The tumor
cell population is a FIXED set of measured Xenium cells, identical at every dx —
which is what makes these statistics grid-comparable."* **That claim is now
tested rather than asserted:** `n_tumor = 23,812` at every one of the five grids,
and the metric moves by at most 10 % across them.

---

## 2. The grid

5 grid/mask configurations × 3 exposure windows × 2 compounds = **30 solves**,
all on the flagship 10x Xenium GBM section, all at the same dose.

| config | window | temozolomide | doxorubicin |
|---|---:|---:|---:|
| **dx12_repair2** *(flagship)* | **6 h** | **2.111** | **3.942** |
| dx12_native | 6 h | 1.994 | 3.797 |
| dx16_repair2 | 6 h | 2.134 | 3.977 |
| dx20_native | 6 h | 1.985 | 3.778 |
| dx30_native | 6 h | 1.994 | 3.826 |
| dx12_repair2 | 3 h | 2.153 | 3.961 |
| dx12_native | 3 h | 2.033 | 3.811 |
| dx16_repair2 | 3 h | 2.171 | 3.996 |
| dx20_native | 3 h | 2.022 | 3.788 |
| dx30_native | 3 h | 2.030 | 3.841 |
| dx12_repair2 | 12 h | 2.095 | 3.933 |
| dx12_native | 12 h | 1.980 | 3.788 |
| dx16_repair2 | 12 h | 2.117 | 3.969 |
| dx20_native | 12 h | 1.969 | 3.770 |
| dx30_native | 12 h | 1.980 | 3.817 |

**The dominant nuisance is mask repair, not grid resolution.** Repaired-mask
configurations sit systematically high (temozolomide 2.10–2.17, doxorubicin
3.93–4.00) and native-mask ones systematically low (1.97–2.03, 3.77–3.84), and
that split is larger than the effect of changing `dx` by a factor of 2.5 or the
exposure window by a factor of 4. It moves **both** compounds the same way, which
is why it does not threaten the separation — but it is worth knowing that the
published 1.99 / 3.80 pair is the **native**-mask reading and 2.11 / 3.94 is the
repaired one.

---

## 3. What this does and does not license

**Established:**

* The P90/P10 delivery spread separates temozolomide from doxorubicin on the
  flagship section by **1.87×** at the shipped configuration, and by at least
  **1.74×** anywhere on the tested nuisance grid.
* Both compounds' values move by **≤ 10 %** across a 2.5× change in grid
  resolution, mask repair on/off, and a 4× change in exposure window.
* The metric is **algebraically invariant** to the exposure-vs-AUC choice and
  **structurally independent** of the distance-binning and distance-metric
  choices that broke the half-distance readout.

**NOT established, and must not be claimed:**

* **That the MAGNITUDE of the spread is right.** This is a robustness test, not a
  validation. Nothing here compares 1.99 or 3.80 to a measurement.
  `results/kyle2004_penetration.md` §6 measured that the engine's perivascular
  gradient is **~14× too shallow** for the one compound where it is a genuine
  out-of-sample prediction, and **~3.7× too compressed** between compounds.
  A robust number can still be a robustly wrong number, and on the evidence of
  that file this one is probably compressed in the same direction.
* **That doxorubicin > temozolomide is a validated ordering.** It is a
  *reproducible* ordering. No measurement of these two compounds' delivery
  heterogeneity on this tissue exists to score it against.
* That this generalises to other compound pairs, tissues, or to pairs whose
  spreads are close together. The margin here is 1.87×; two compounds 1.1× apart
  would sit inside this same noise.

**The honest one-line statement.** *The delivery-spread separation is robust to
every analysis choice that can touch it — the ordering never inverts and never
closes — but robustness is not accuracy, and the same pipeline's absolute
gradient is measured ~14× too shallow.*

---

## 4. Reproducing this

```
./.venv/bin/python experiments/delivery_spread_noise.py --stage degenerate
./.venv/bin/python experiments/delivery_spread_noise.py --stage all
```

~16 min for the full grid (30 solves across 5 reconstructions).
