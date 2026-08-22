# Why is the reconstructed-section spatial readout noisy? — the axis that moves it is the BIN WIDTH, the mechanism is the estimator, and the published inversion warning is wrong

**Status: RUN.** Driver: `experiments/spatial_noise_decomposition.py`. Data:
`results/spatial_noise_decomposition.json`. Tests:
`tests/test_spatial_noise_estimators.py`.

**Nothing was tuned.** No shipped default, constant or parameter was changed.
`calibrate_alpha()` is asserted live at `22.223477711172563` in every stage and
holds. The two library additions (§7) are **additive only** — every pre-existing
key of `distance_profile` returns byte-identical values and no shipped number
moves. **No assertion was deleted or weakened.**

**Reproduction check, performed before any analysis:** replaying the cached
fields through this file's own binning reproduces
`results/kyle2004_penetration.md` §7's doxorubicin half-distance as
`45.18752771106239` — every digit of the published value. The replay path is the
shipped path.

---

## 0. Verdict up front

**Four stages were separated. One of them carries the spread, and it is not the
one the brief pointed at.**

| stage | axis | paired range of the doxorubicin d½, all other axes held fixed |
|---|---|---:|
| 3. BIN | **distance bin width** | **18.37 µm** |
| 2. DISTANCE | distance metric (primary vs legacy) | 5.19 µm |
| 1. SOLVE | tissue arm (dx 12→30 **and** mask repair) | 2.43 µm |
| 1. SOLVE | **mask repair alone, at fixed dx** | **0.47 µm** |
| 3. BIN | distance truncation (p99 / p99.5 / p100) — *an axis §8 never swept* | 0.00 µm † |
| — | AUC vs final concentration | 0.30 µm |

† **This axis is inert for one panel readout and not for the other, and that is a
property of the READOUT, not of the section.** 0.00 µm for the empirical
half-distance, which crosses at ~45 µm and never sees the tail; **2.80 µm for the
fitted `L`** (`exp+floor`, `dx12_repaired`/`dx12_native`/`dx30_native` × primary+
legacy × AUC+conc × bins 5/10/20, doxorubicin), which fits the whole profile out to
the truncation and therefore does. **Both are on the panel.** Anyone reading only
the half-distance row would conclude the truncation choice cannot matter; it
can, for the other number on the same figure. The truncation is currently an
unswept default (`max_distance_um = p99.5`) inside `distance_profile`.

**And the single most decision-relevant finding, which needed no new code at
all — it is in `results/kyle2004_penetration_data.json` as already saved:**

> **In all 36 registered nuisance cells, doxorubicin's half-distance exceeds
> daunorubicin's. Zero inversions. The paired difference is +6.14 to +13.89 µm
> and the paired ratio 1.135 to 1.287.**

`results/kyle2004_penetration.md` §8 states:

> *"A user who re-ran the same two compounds at `dx = 30 µm` on the native mask
> with 20 µm bins could get the opposite ordering."*

**That sentence is false, and the file's own data refutes it.** At
`dx30_native`, native mask, 20 µm bins, primary metric, AUC, the readouts are
doxorubicin **56.59** and daunorubicin **49.60** µm — correctly ordered, as they
are in every other cell of the grid. The C2 failure is real, but it is a failure
of the **marginal** bar (compare doxorubicin at one nuisance setting against
daunorubicin at a *different* one), not of the ordering. **Whether the marginal
bar or the paired bar is the right one for the product claim is a decision, and
it is §8's first decision.**

**The axis is the bin width; the mechanism is stage 4, not stage 3.** A matched
control — the same log-linear estimator run binned and unbinned, model held fixed
— shows that removing the binning changes the answer by 0.5 % and its span by 4 %
(§6.3). **Binning is not the disease.** What makes a readout bin-width sensitive
is what the *estimator* does with a bin average: taking its reference from the
first bin while measuring from zero (§2), plotting a bin at its centre rather than
where its value sits (§3), and carrying an additive background term through a bin
average (§6.3).

**Three defects were found and their mechanisms derived, not fitted.** Two are
fixed as additions (§7); the third needs a signature change and is left as a
decision (§8).

**One popular fix does not work, and the negative result matters:** making the
depth readout binning-free stabilises it by destroying the signal. The bin-free
exposure-weighted depth separates the pair by **1.03×** where the half-distance
separates it by **1.16×** (doxorubicin 38.14 µm against daunorubicin 37.01 µm,
`dx12_repaired`/`dx12_native`/`dx30_native` × primary+legacy × AUC+conc ×
p99/p99.5/p100, 108 cells per compound), and its median variant is exactly
degenerate — 54 of 108 paired cells tie to the last digit. **"Binning-free" is not
why `P90/P10` survives** (§6).

> **THIS TRAVELS BEYOND THIS FILE.** `results/delivery_spread_noise.md` §1
> justifies `P90/P10`'s robustness by its structural blindness to distance
> binning. **That justification does not survive this test.** A depth readout made
> binning-free is *not* thereby robust — it is robust the way a constant is. What
> actually protects `P90/P10` is that it is a statistic of the **value**
> distribution and never touches the distance axis at all (§6). **The metric may
> well still be the right one; the argument given for it is not.** This is the
> **second** justification for the same headline metric to fall this week — a
> sibling separately retracted the reach panel's claim that the viable-tumour
> denominator carries "~10× grid sensitivity" when the measured figure is 0.6 %.
> Two independent bad arguments for one metric is worth surfacing even though
> neither falsifies the metric.

---

## 1. What was measured, and how the attribution is exact

The reduction has four stages, and one moved number cannot say which one moved:

```
1. SOLVE     the PDE field on a tissue mask     (grid dx, mask repair)
2. DISTANCE  a distance-to-vessel per pixel     (primary vs legacy metric)
3. BIN       (d_i, C_i) pixels -> a profile     (bin width, truncation)
4. ESTIMATE  the profile -> one number          (empirical d½ / exp+floor fit)
```

Six real-section solves (3 tissue arms × 2 compounds, `dx = 12 µm` repaired,
`dx = 12 µm` native, `dx = 30 µm` native) were run **once** and cached. Every
downstream choice is then replayed on the *identical* fields, so stages 2–4 cost
no further solve and the attribution is measured rather than inferred from
marginal means. The registered arm-B grid of `kyle2004_penetration.md` §8 is
reproduced exactly, plus one axis that file left **implicit**: the distance
truncation inside `distance_profile` (`max_distance_um = p99.5` of the tissue
distance distribution, never swept). It turns out to matter not at all for the
half-distance — 0.00 µm — and modestly (2.80 µm) for the fitted `L`, which is
worth knowing precisely because it was never checked.

**Thirteen** estimators were run over all 108 cells per compound. Nothing was
selected on its score; §5 and §6 give the mechanical argument for each, and the
two promotions §8 asks for are argued by **convergence**, not by improvement.

---

## 2. Defect 1 — the half-distance is an interval measured from two different origins

`empirical_half_distance_um` normalises by `v[0]`, the **mean of the first bin**,
and then reports the crossing as a distance measured **from zero**. But `v[0]` is
not the profile's value at `d = 0`: it is its average over `[0, b)`, i.e. its
value at an effective distance of roughly `b/2`. **The reference is taken at
`b/2` and the answer is measured from 0.**

For a pure exponential the consequence is exact and closed-form. With bin width
`b`, the bin means of `C0·exp(−d/L)` are

```
v_j = C0 · (L/b) · (1 − exp(−b/L)) · exp(−j·b/L)
```

so `v_j / v_0 = exp(−j·b/L)` — the bin-averaging prefactor cancels — the profile
halves at fractional bin index `j* = L·ln2 / b`, and the estimator reports that
at the bin **centre** `(j*+0.5)·b`:

```
d_half_shipped  =  L·ln2  +  b/2          (exact; the truth is L·ln2)
```

`--stage bias` verifies this on re-binned synthetic profiles with a known answer,
and `tests/test_spatial_noise_estimators.py` pins it:

| bin `b` (µm) | shipped d½ | error | closed form `b/2` |
|---:|---:|---:|---:|
| 5 | 37.173 | **+2.516** | 2.5 |
| 10 | 39.907 | **+5.249** | 5.0 |
| 20 | 45.414 | **+10.757** | 10.0 |

(true d½ = `50·ln2` = 34.657 µm; the residual is second order, bounded by
`b²/(8L)` — the estimator chords a convex curve between bin centres.)

**This is not a modelling choice. It is an arithmetic error of `b/2`.** Over the
registered bin widths 5→20 µm it accounts for exactly **7.5 µm** — a third of the
23.7 µm arm-B span, and 7.5 of the bin axis's 18.37 µm. The fix is a change of
origin with no free parameter: measure the crossing from the anchor's own
location. `anchor_corrected_half_distance_um` (§7) does exactly that, recovers the
true value at every bin width on a profile with no floor, and on the real section
takes the bin axis from **18.369 → 10.869 µm** — a reduction of **exactly 7.500**,
the closed form to four significant figures, on real reconstructed tissue.
**The rest of the bin axis is §6.3's floor interaction, not this.**

---

## 3. Defect 2 — a bin's mean value does not sit at the bin's centre

`distance_profile` labels every bin by its geometric centre. On a 2-D section the
tissue area at distance `d` from a vessel **rises with `d`**, so the pixels inside
a bin are strongly skewed toward its far edge and the bin's mean value corresponds
to a distance systematically **larger** than the centre — by an amount that grows
with the bin width.

Measured on the flagship section (`dx12_repaired`, primary metric):

| bin `b` (µm) | bin centre | bin's **actual mean** distance |
|---:|---:|---:|
| 5 | 2.50 | **3.33** |
| 10 | 5.00 | **6.59** |
| 20 | 10.00 | **12.69** |

This distorts the **shape**, not just the anchor, so it hits the fitted
penetration length — the number on the panel — as well as the half-distance.

**The correctness test is convergence, and it is deliberately not a score.** As
`b → 0` the binning stops doing anything and the profile becomes the pixel cloud
itself, so the fine-bin value **is** the definition of the quantity. An estimator
is right if it still agrees with that limit at coarse bins.

`--stage converge`, doxorubicin, `dx12_repaired`, primary metric, AUC:

| bin (µm) | d½ shipped | d½ anchored | **`L` at bin centre** | **`L` at bin mean-distance** |
|---:|---:|---:|---:|---:|
| **2** | 39.31 | 37.95 | **46.38** | **46.27** |
| 3 | 42.16 | 40.15 | 46.67 | 46.41 |
| 4 | 39.38 | 36.69 | 46.81 | 46.38 |
| 5 | 40.03 | 36.69 | 47.14 | 46.55 |
| 7.5 | 41.94 | 36.92 | 47.93 | 46.57 |
| **10 (shipped)** | **45.19** | **38.60** | **49.30** | **46.89** |
| 15 | 51.97 | 42.35 | 52.94 | 48.12 |
| 20 | 58.54 | 45.85 | 58.75 | 50.88 |
| 30 | 68.94 | 50.27 | 73.07 | 57.38 |
| 40 | 79.20 | 54.81 | 82.51 | 65.30 |

Both labellings agree at `b = 2` (46.38 vs 46.27), as they must. **The
centre-labelled fit then walks away from its own converged value; the
mean-distance-labelled one does not.** At the shipped `BIN_UM = 10` the panel's
penetration length reads **49.30 µm** against a converged **46.27 µm** — a **+6.5 %
labelling bias** — while the mean-distance labelling reads 46.89 µm, **+1.3 %**.
At `b = 20` the centre labelling is **+27 %** off.

The same is true of the half-distance: **45.19 µm at the shipped `b = 10` against
39.31 µm at `b = 2`, i.e. +15 %.** (The shipped d½ column jitters by ±1.5 µm at
fine bins — 39.31, 42.16, 39.38, 40.03 — because a crossing-point estimator
interpolates between discrete bin centres. That jitter is irreducible for any
crossing statistic; it is ~4× smaller than the bias at the shipped bin width and
far smaller than it at coarse bins.)

**And the fix improves every other axis too**, which is the signature of a real
mechanism rather than a lucky one — the labelling error interacts with the grid
and the metric because those change the pixel-distance density inside each bin:

| estimator | arm | field | metric | **bin** | trunc |
|---|---:|---:|---:|---:|---:|
| `L`, bins at their centre (**shipped**) | 10.28 | 0.27 | 5.01 | **12.90** | 2.80 |
| `L`, bins at their mean distance | **7.31** | 0.27 | **3.81** | **5.05** | 2.74 |

---

## 4. Defect 3 — the legacy distance map is lattice-quantised, and the aliasing is silent

The PRIMARY metric (`direct_vessel_distance_um`, a cKDTree query against real
endothelial-cell coordinates) is continuous: on the flagship section its 201,674
tissue pixels take **201,475 distinct distance values**.

The LEGACY metric (`vessel_mask_from_density` → `distance_from_vessel_um`, an EDT
on the grid) takes **333**. Its values are an exact lattice comb —
`0, 12, 16.971 (=12√2), 24, 26.833 (=12√5), …` — because a Euclidean distance
transform on a lattice of pitch `dx` can only return lattice distances.

**Binning a comb finer than its pitch is meaningless, and `distance_profile`
hides it**: empty bins are skipped by `if n == 0: continue` with no record, so an
aliased profile is indistinguishable from a well-sampled one in the output. The
module constant is commented *"Distance-bin width (µm); comfortably above
dx = 12 µm"* — which stops being true the moment the `dx30_native` arm is reduced
at the same bin width.

**This is a defect of the LEGACY metric, not of coarse grids as such.** The
primary direct-cell distance is computed from continuous cell coordinates, so it
is unaffected: at `dx = 30 µm` and `b = 5 µm` it still yields **0 empty bins** and
a sensible first-bin mean distance of 3.26 µm. The quantisation enters only
through the EDT, which can return nothing but lattice distances.

Measured (doxorubicin, AUC, p99.5):

| arm | metric | bin | empty bins | first-bin mean distance | d½ |
|---|---|---:|---:|---:|---:|
| dx12_repaired | primary | 5 | 0 | 3.33 | 40.03 |
| dx12_repaired | legacy | 5 | **2** | **0.00** | 40.63 |
| dx12_repaired | legacy | 10 | 0 | **0.00** | 40.00 |
| dx30_native | primary | 5 | 0 | 3.26 | 39.36 |
| **dx30_native** | **legacy** | **5** | **24** | **0.00** | **61.71** |
| dx30_native | legacy | 10 | **6** | 0.00 | 57.69 |
| dx30_native | legacy | 20 | 0 | 0.00 | 62.69 |

`first-bin mean distance = 0.00` means the anchor bin contains **only the
zero-distance pixels of the vessel mask itself** — there is nothing else in
`[0, b)` to average with, because the next populated shell is at `dx`.

**These are the cells that produced arm B's extremes.** The largest paired
differences in the whole registered grid are all `dx30_native × legacy` cells, the
two extremes (13.72 and 13.89 µm) at bin 5. Excluding only the cells where
`bin_um < dx` — i.e. the provably aliased ones — the paired spread collapses:

| registered grid | n cells | paired diff (µm) | paired ratio |
|---|---:|---|---|
| all 36 | 36 | 6.14 – **13.89** | 1.1345 – **1.2870** |
| aliased legacy cells removed | 22 | 6.17 – **7.68** | 1.1345 – **1.2080** |

---

## 5. What is NOT the problem — four falsified leads

Reported because knowing where the noise is *not* is worth as much as knowing
where it is, and because two of these were the leads this investigation started
from.

* **Mask repair. 0.47 µm — essentially nothing.** The brief's strongest lead was
  `results/delivery_spread_noise.md` §2's finding that *"the dominant nuisance is
  mask repair, not grid resolution."* **That is true of `P90/P10` and does not
  transfer to the half-distance.** At fixed `dx = 12 µm`, repairing the mask moves
  the doxorubicin d½ by a median of −0.32 µm and never by more than 0.47 µm, on a
  45 µm number.

  `--stage repair` splits it further into its two mechanisms — repairing the mask
  changes **the solved field** (the domain connects, and `build_real_fields`
  re-derives the source) *and* **the pixel set** (7,922 added interstitial pixels
  enter the reduction):

  | | doxorubicin d½ | `P90/P10` |
  |---|---:|---:|
  | native solve, native pixels | 45.64 | 3.808 |
  | repaired solve, **native** pixels *(field change only)* | 45.80 | 3.805 |
  | repaired solve, repaired pixels *(+ pixel-set change)* | 45.19 | 3.838 |

  **The field change is +0.16 µm and −0.003; the pixel-set change is −0.61 µm and
  +0.033.** Whatever mask repair does to `P90/P10`, it does through *which samples
  enter the statistic*, not through the physics — and to the half-distance it does
  almost nothing at all. The `ARMS` docstring's concern that a fragmented mask
  "suppresses lateral transport" is not detectable in either readout at `dx = 12`.

* **`MIN_BIN_COUNT`. Inert where it ships, and a sensible guard where it is not.**
  The hypothesis was that a *fixed absolute* count threshold applied to bins whose
  population scales with bin width would silently shorten the analysed range at
  small `b`. **At the shipped p99.5 truncation it drops 0 bins on both `dx = 12 µm`
  arms and at most 1 on `dx30_native`**, and removing it moves the fitted `L` by at
  most 0.19 µm. Falsified as a cause.

  It is *not* inert at the p100 truncation this file added, where it drops 23–83
  bins and moves `L` by up to 2.0 µm — but that is the filter working, not failing:
  the section's distance distribution has a long thin tail (`max = 789 µm` against
  a p99.5 of `152 µm`) whose bins genuinely hold fewer than 30 pixels. **The
  interaction is real but it is bounded by the truncation, and the shipped
  truncation keeps it out.** The two guards are doing the same job twice.

* **The estimator refusing.** `fit_exponential_plus_floor` refuses below four
  usable bins. **It refused in 0 of 216 fits**, with `r² = 0.981–0.998` throughout.
  The refusal branch is not this readout's problem. (The refusal that *does* fire
  in practice is `empirical_half_distance_um` returning `None` for a profile that
  never halves — temozolomide — and that is a measured property of the profile,
  documented as such, not a defect.)

* **AUC vs final concentration (0.30 µm) and distance truncation (0.00 µm on the
  half-distance).** Both negligible for the readout that failed C2. The truncation
  result is new — it was an unswept implicit choice — and it is worth one caveat:
  it is 0.00 µm for the half-distance, which crosses long before the tail, but
  **2.80 µm for the fitted `L`**, which does see the tail. An implicit choice that
  does not matter for one readout on the panel matters modestly for the other.

---

## 6. Would a binning-free depth statistic be better? — NO, and the reason is instructive

`results/delivery_spread_noise.md` §1 attributes `P90/P10`'s survival to its
structural blindness: *"there is no distance binning in this metric."* The
obvious inference is that a binning-free **depth** statistic would inherit that
robustness. **It does not, and testing it is what shows why.**

Two binning-free depth statistics were computed straight from the pixel cloud —
the exposure-weighted **mean** and **median** of the distance distribution
(`Σ C_i d_i / Σ C_i` and its median), plus each divided by the same statistic
over the tissue itself so the geometry cancels. No bins, no bin width, no
truncation, no anchor bin, no fit, no refusal branch.

| estimator | DOX | DAU | worst span over the grid | sep/span | paired ratio | inversions |
|---|---:|---:|---:|---:|---|---:|
| empirical d½ (**shipped**) | 45.57 | 39.22 | 23.70 | 0.268 | 1.135 – 1.287 | **0/108** |
| `L`, bins at mean distance | 50.84 | 44.15 | **20.16** | **0.332** | 1.114 – 1.176 | **0/108** |
| bin-free mean depth | 38.14 | 37.01 | 26.24 | 0.043 | **1.024 – 1.035** | 0/108 |
| bin-free mean depth / tissue | 0.801 | 0.779 | **0.038** | **0.590** | 1.024 – 1.035 | 0/108 |
| bin-free **median** depth | 34.58 | 33.41 | 33.17 | 0.035 | **1.000** – 1.035 | **54/108** |
| **binned** log-linear `L` (no floor) | 69.76 | 61.40 | 56.45 | 0.148 | 1.084 – 1.139 | 0/108 |
| **unbinned** pixel log-linear `L` (same model) | 69.41 | 60.99 | 54.37 | 0.155 | 1.104 – 1.141 | 0/108 |
| `P90/P10` over the same pixels | 3.809 | 4.275 | 0.134 | 3.486 | — | — |

**Three things to read off this.**

1. **Binning-free stabilises by deleting the signal.** The bin-free mean depth
   separates the pair by **1.03×** where the half-distance separates it by
   **1.16×** and the fitted `L` by **1.15×**. A centroid of the whole distance
   distribution is dominated by the bulk of tissue at moderate distance, where
   both compounds' fields have nearly the same shape. It is stable the way a
   constant is stable.

2. **(§6.2) The bin-free median is degenerate on a quantised distance map.** 54 of 108
   paired cells give a difference of **exactly 0.000** — every legacy-metric cell,
   because a weighted median on a 333-value comb lands on the same atom for both
   compounds. A median cannot have more resolution than its axis.

3. **(§6.3) Deleting the bins does not delete the problem — and a matched control
   shows binning is barely the problem at all.** The log-linear estimator was run
   **binned and unbinned with the model held fixed** (pure exponential, log space,
   no floor term, count weights), which is the only way to attribute spread to the
   binning stage rather than to the model:

   | | central `L` | span | arm | **bin** | metric |
   |---|---:|---:|---:|---:|---:|
   | log-linear, **binned** | 69.76 | 56.45 | 24.77 | **2.76** | 0.83 |
   | log-linear, **unbinned pixels** | 69.41 | 54.37 | 25.87 | *blind* | 0.21 |
   | exp+floor, binned at centres (**shipped**) | 50.12 | 28.39 | 10.28 | **12.90** | 5.01 |
   | exp+floor, binned at mean distances | 50.84 | 20.16 | 7.31 | **5.05** | 3.81 |

   **Removing the binning entirely moves the central value by 0.5 % and the span
   by 4 %.** The 54–56 µm span of both log-linear rows is the *missing floor
   term* — a pure exponential over a profile with a plateau reads systematically
   longer, exactly as `vessel_profile.py` documents for `near_exp` vs `exp_floor`
   — not the binning.

   Read the other way: **the binning stage contributes only 2.76 µm to a
   log-linear estimator, but 12.90 µm to the exp+floor one.** The difference is
   that the bin-averaging prefactor `(L/b)(1−exp(−b/L))` cancels out of a
   log-space fit (it is a constant offset absorbed by `C0`) but *not* out of a fit
   with a floor, because it shrinks the exponential component while leaving the
   background alone — which is visible directly as the fitted floor fraction
   drifting 0.132 → 0.108 as `b` goes 5 → 20. Of that 12.90 µm, §3's labelling fix
   removes 7.85 and this prefactor–floor interaction is the ~5 µm that remains.

   **So the honest attribution inside stage 4 is: bin-width sensitivity is not a
   property of binning, it is a property of estimators that carry an additive
   background term through a bin average.**

**So `P90/P10` does not survive because it is binning-free. It survives because
it is a statistic of the VALUE distribution and never touches the distance axis
at all.** Every readout that reduces *along distance* inherits the distance map's
quantisation, its geometry and its anchor problem, bins or no bins. That is a
structural fact about the quantity, and it is the honest answer to "would a
different readout be inherently more stable": **yes — but only by no longer being
a depth.**

**What a defensible binning-free depth statistic would look like**, stated for the
record since it was asked: the **ratio** form is the only one worth anything —
`Σ C_i d_i / Σ C_i`, divided by the same quantity with `C ≡ 1`. It is
dimensionless, homogeneous of degree 0 in `C` (so the AUC-vs-concentration axis is
degenerate for it by the same algebra `delivery_spread_noise.md` §1 measures for
`P90/P10`), and the tissue geometry that the grid and mask change divides out —
which is why its span is 0.038 on a value of 0.80 and its sep/span of **0.590** is
the best of any distance readout here. **But its effect size is 2.8 %**
(0.801 vs 0.779). That is not a figure a customer can act on, and it should not be
promoted to the panel on the strength of a stability score.

---

## 7. What was fixed

Both changes are **additive**. Every pre-existing key of `distance_profile`
returns byte-identical values, nothing in the module calls the new estimator, and
**no shipped number moves.** `tests/test_spatial_noise_estimators.py` (14 tests)
pins the derivation rather than the outcome.

**`tmesim/validation/vessel_profile.py`**

1. **`anchor_corrected_half_distance_um`** — the shipped half-distance minus the
   distance at which its own reference sits, i.e. the half-distance as a *length*
   rather than as a *coordinate*. Removes the `b/2` term of §2 exactly. No free
   parameter; it is a change of origin. It propagates the shipped estimator's
   honest `None` for a profile that never halves rather than manufacturing a
   number.

2. **`distance_profile` now also returns `distance_mean_um` and `n_empty_bins`.**
   `distance_mean_um` is each bin's count-weighted mean distance — where the bin's
   value actually sits (§3) — and is what makes the §3 fix available to a caller.
   `n_empty_bins` counts the bins that the silent `if n == 0: continue` dropped,
   which is the direct signature of the lattice aliasing of §4 and was previously
   unobservable.

**Tests added** (`tests/test_spatial_noise_estimators.py`, 14 tests): the
closed-form `b/2` bias with a *derived* tolerance of `b²/(8L)` (the second-order
chord bound, not a number chosen to pass); the corrected estimator's recovery of
the true value at every bin width; that its span is >8× smaller; that
`distance_profile`'s original keys and values are unchanged; that `n_empty_bins`
fires on a 30 µm comb binned at 10 µm and is zero at 30 µm; and that the
correction propagates `None`.

### 7.1 The tests were mutation-checked, and one of them was failing to bite

A green check is not evidence until the defect it claims to catch has been planted
and seen to fail it. Each defect below was applied to a **copy** of
`vessel_profile.py` loaded under another module name (the working tree was never
mutated — sibling agents were running suites against it) and the assertions
re-run:

| planted defect | caught by |
|---|---|
| control — no mutation | *(all 8 assertions green)* |
| anchor correction does nothing (returns the shipped crossing) | `anchor_recovers_true_value`, `anchor_span_8x_smaller` |
| `distance_mean_um` silently returns the bin **centre** | `distance_mean_SKEWED_half` |
| `n_empty_bins` hardcoded to 0 | `n_empty_bins_detects_comb` |
| `empirical_half_distance_um` secretly de-biased (`crossing − d[0]`) | `shipped_bias_is_b_over_2` + 2 others |
| `distance_profile`'s `mean` becomes the bin **median** | `original_keys_unchanged` |
| `distance_profile`'s bin centres shifted by a half-bin | `original_keys_unchanged` + 4 others |
| anchor correction returns `0.0` instead of the honest `None` | `none_propagates` |

**Every mutation bites. One did not, until the test was rewritten.** The
`distance_mean_um` check as first written compared the bin's mean distance
against the bin centre on a *uniform-in-distance* cloud — where the two are
**equal by construction**. It passed against an implementation that simply
returned the centre, i.e. it certified exactly the defect it existed to catch.
The test now carries a second half on a cloud whose density rises with distance
(`p(d) ∝ d`, so the first bin's mean is analytically `40/3 = 13.33 µm` against a
centre of 10.00), which is the real 2-D situation and the one where the two must
differ. **Reported because it is the same failure mode this repository has been
finding all week, and it was in my own work.**

---

## 8. What is a DECISION, not a fix — brought rather than made

### 8.1 Is C2's marginal-disjointness bar the right bar? *(the important one)*

The registered C2 asks whether the two compounds' **marginal** intervals over the
whole nuisance grid are disjoint. That question compares doxorubicin at one
setting against daunorubicin at a *different* setting. The dominant nuisance —
bin width, 18.37 µm — is almost perfectly **common mode**: going from `b = 5` to
`b = 20` on the flagship arm moves doxorubicin by **+18.51 µm** and daunorubicin by
**+17.49 µm**, so **17.5 of the 18.4 µm cancels out of the comparison** any single
run actually makes, while inflating both marginals in full.

Measured: **0 inversions in 36/36 registered cells, and 0 in 108/108 of this
file's expanded grid, for every distance-binned estimator and for both bin-free
*mean* depth statistics.** (The two bin-free *median* statistics register 54/108,
but those are exact **ties**, not inversions — see §6.2. `P90/P10` registers
108/108 against a counter that assumes doxorubicin should be the larger; it is
legitimately the *smaller*, because a worse penetrator has the more heterogeneous
delivery. Read the JSON's `n_inversions` for those three rows with that in mind.)
Paired ratio
1.135–1.287 on the registered grid; 1.135–1.208 once the provably aliased cells
of §4 are removed.

* **The case for keeping the marginal bar:** the panel publishes an **absolute**
  half-distance. A customer who reads "45 µm" from one report and "39 µm" from
  another produced under different settings is genuinely misled, and the marginal
  bar is the one that protects them.
* **The case for a paired bar:** the surviving product claim is explicitly a
  *within-map comparison* — "this region got N× less than that region" — and
  `kyle2004_penetration.md` §1 justifies it precisely because *"within a single
  map the uncertain per-compound absolutes cancel."* A criterion that forbids the
  absolutes from cancelling is testing a claim the file says it is not making.

**My recommendation is to report both and retire neither.** The marginal-vs-paired
choice is left **open** — it is a product decision about what the panel publishes,
not something this file's measurements settle.

**DONE 2026-08-22:** `results/kyle2004_penetration.md` has been corrected in place
at three points — §0's "the compound difference disappears into the method", §8's
`dx = 30 µm` inversion scenario, and §10.1's summary bullet. Each retraction
**quotes the retracted sentence rather than deleting it** and is dated. **The C2
FAIL is not retracted anywhere** — only the claim about the *ordering*, which is a
different claim from the one C2 scores.

### 8.2 Should the bin-mean labelling be promoted into `penetration_length_fit`?

This is the one with a real payoff: it is a **defect** by the convergence argument
of §3, and fixing it takes the panel's penetration length from **+6.5 %** to
**+1.3 %** of its own converged value, cuts its bin sensitivity by **2.6×**, and
improves its grid and metric sensitivity as well.

**The cost is that it moves a published number**: doxorubicin's flagship `L` goes
from **49.30 → 46.89 µm** (−4.9 %), and every `L` in `results/vessel_profile.md`,
`results/delivery_metrics.md` and the panel moves with it. That is a
re-publication, not a patch, so I have not done it.

### 8.3 Should the anchor correction be promoted?

Same shape, larger blast radius: it moves **every** published half-distance down
by roughly half a bin width (45.19 → 38.60 µm at the shipped `BIN_UM = 10`). The
arithmetic in §2 is not in question. What is in question is whether "half-distance"
as this repository has published it means *the distance at which the profile
halves* (in which case the shipped number is wrong by `b/2`) or *the coordinate of
the half-crossing measured from the vessel* (in which case it is right and the
first bin's anchor is the thing that is mislabelled). **They are the same number
only in the limit `b → 0`.**

Note also that the correction does **not** remove all bin sensitivity — 18.37 →
10.87 µm with a centre anchor, → 9.02 with a mean-distance anchor. The residual is
the **floor**: with a non-zero background the bin-averaging prefactor stops
cancelling out of `v_j / v_0`, so the meaning of "half" itself becomes bin-width
dependent. `fit_exponential_plus_floor` already offers the alternative definition
(`d_half_total_um`, half of `C0 + floor`). **Which definition of "half" the panel
should publish under a non-zero background is a design question, not a bug.**

### 8.4 Should `distance_profile` refuse to bin below the grid pitch?

§4 shows that `bin_um < dx` produces a comb with up to **24 of 45 bins empty** and
a first bin containing nothing but the vessel mask itself. `distance_profile`
cannot currently detect this: it never receives `dx`. A guard needs either a new
argument or a caller-side check. `n_empty_bins` (§7) makes the condition
*observable* now; making it *fatal* changes a signature, and **14 of the 36 cells
of the registered arm-B grid are already in that condition** on the legacy metric.

There is a subtlety worth deciding deliberately rather than by default: a guard
keyed on `bin_um < dx` would also fire on the 14 **primary**-metric cells at the
same bin widths, where the distance map is continuous (201,475 distinct values)
and the binning is perfectly well posed. **The condition that actually matters is
whether the distance map is quantised, not what `dx` is** — and `n_empty_bins`,
which measures the symptom directly, is the better predicate for that.

---

## 9. What this licenses saying

**Established:**

* The arm-B instability is dominated by **one** reduction stage — the distance
  bin width (18.37 µm, against 5.19 for the distance metric, 2.43 for the tissue
  arm and 0.47 for mask repair alone).
* Its mechanism is **derived, not fitted**: an exact `+b/2` offset from measuring
  a half-distance from one origin and its reference from another (§2), plus a
  bin-labelling error that grows with bin width and biases the fitted penetration
  length by **+6.5 % at the shipped setting** and **+27 % at `b = 20`** (§3).
* The LEGACY distance map is **lattice-quantised to 333 distinct values** and
  binning it below the grid pitch aliases; the four most extreme cells in the
  registered grid are all provably aliased ones (§4).
* **The compound ordering never inverts** — 0 inversions in 36/36 registered cells
  and 108/108 expanded — so C2's failure is a property of the *marginal* bar, not
  of the ordering (§0, §8.1).
* Mask repair, `MIN_BIN_COUNT`, the fit's refusal branch, the AUC-vs-concentration
  choice and the distance truncation are all **falsified** as causes (§5).
* A binning-free depth statistic is stable but **signal-free** (1.03× against
  1.16×), and its median form is degenerate on a quantised axis (§6).
* **Binning itself is nearly innocent.** With the model held fixed, removing the
  binning stage entirely moves a log-linear `L` by 0.5 % and its span by 4 %. What
  makes an estimator bin-width sensitive is carrying an **additive background
  term** through a bin average — 2.76 µm of bin sensitivity without a floor
  parameter, 12.90 µm with one (§6.3).

**NOT established, and must not be claimed:**

* **That anything here makes the map more ACCURATE.** Every number in this file
  is a robustness measurement against the pipeline's own converged limit. It says
  nothing about the measurement. `results/kyle2004_penetration.md` §6 stands
  untouched: the engine's absolute perivascular gradient is **~14× too shallow**
  and its compound separation **~3.7× too compressed**. **A readout that is
  bin-width invariant and 14× too shallow is still 14× too shallow**, and the
  §3 fix moves the penetration length in the *shallower* direction (49.30 → 46.89),
  i.e. slightly further from the measurement.
* **That C2 should be re-scored as a pass.** It fails on the bar that was
  registered. Whether that bar is the right one is §8.1, and it is a decision for
  someone else.
* **That the ordering is validated.** It is *reproducible* across every analysis
  choice tested. `kyle2004_penetration.md` §5 remains the governing statement on
  why it came out right: an unidentified structure-only prior.
* That any of this generalises beyond one pair, one section, one tumour type.

---

## 10. Reproducing this

```
./.venv/bin/python experiments/spatial_noise_decomposition.py --stage bias       # seconds
./.venv/bin/python experiments/spatial_noise_decomposition.py --stage fields     # ~14 min, 6 solves
./.venv/bin/python experiments/spatial_noise_decomposition.py --stage decompose  # seconds
./.venv/bin/python experiments/spatial_noise_decomposition.py --stage converge   # seconds
./.venv/bin/python experiments/spatial_noise_decomposition.py --stage repair     # seconds
```

`--stage fields` caches the solved fields outside the repository (they are ~18 MB);
every other stage replays them, so the whole analysis re-runs in seconds once the
six solves are done. Measured solve times: `dx12_repaired` 378 s, `dx12_native`
382 s, `dx30_native` 64 s.

### 10.1 Suite state, stated precisely

**Files this work touched, and nothing else** (`git status`):

| file | status |
|---|---|
| `tmesim/validation/vessel_profile.py` | modified — **additive only** (§7) |
| `experiments/spatial_noise_decomposition.py` | new |
| `tests/test_spatial_noise_estimators.py` | new |
| `results/spatial_noise_decomposition.{md,json}` | new |
| `results/kyle2004_penetration.md` | modified — **three dated corrections, §8.1; no result retracted, no number changed** |

**No assertion was deleted or weakened anywhere.** No file under `app/`,
`landing-3d/`, `docs/`, `data_cache/` or any other `results/*.md` was touched.

| test selection | result |
|---|---|
| `tests/test_spatial_noise_estimators.py` | **14 passed** |
| + `test_vessel_network`, `test_vessel_source`, `test_claim_guard` | **84 passed** |
| `tests/test_propagation.py` | **54 passed** |

Plus the mutation check of §7.1: **8 planted defects, 8 caught, control green.**

**A full-suite number is deliberately not quoted here.** Concurrent agents are
editing `app/server/`, `app/web/`, `tmesim/molecular/descriptors.py`,
`tmesim/validation/propagation.py` and `tests/test_propagation.py` in the same
working tree, and the tree now carries ~60 tests that are not mine — so a total is
neither comparable to the `23e51f6` baseline of 1221 nor attributable. A full run
started at 03:03 reported 2 failures in `tests/test_propagation.py` that **pass in
isolation**; they were a file-changed-under-the-run artefact (my own edits to
`experiments/spatial_noise_decomposition.py`, which that guard globs). The clean
full run is sequenced with the coordinator.
