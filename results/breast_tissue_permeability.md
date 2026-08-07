# Measured breast permeability makes it WORSE — and that localises the defect

**A falsified fix, reported as one — and the failure is more informative than
a success would have been.** Replacing the blood-brain-barrier permeability
with an independently measured *breast* one does not close the breast anchor
gap. It widens it for one drug and barely moves the other. But the way it
fails identifies which parameter is actually wrong, which two sessions of
work had not.

Driver: `experiments/breast_tissue_permeability.py`. Data:
`results/breast_tissue_permeability_data.json`. Nothing under `tmesim/` was
modified — the substitution is applied by overriding `bbb_penetrance` on the
compiled `DrugParams`.

## 1. What was tried, and why it was the right thing to try

`PS = vascular_density * bbb_penetrance`, and `bbb_penetrance` is a
structure-derived **blood-BRAIN-barrier** term. Breast has no such barrier.
This is defect 2 of the three in
`results/transvascular_exchange_revalidation.md`.

The measured data says the replacement should be **drug-independent**: in
breast tumour the 5-FU:methotrexate permeability ratio is **1.02x**
(PMID 9205062 / 9679959) against **28.1x** from `bbb_penetrance`. The measured
ratio tracks barrier tightness monotonically — ~1.0x breast tumour, ~2x
skeletal muscle, 6-12x normal brain — so a single tissue-level permeability is
what the measurement supports, not a per-drug one.

Applied by setting `bbb_penetrance = PS_target / mean(vascular_density)` over
the tissue mask, which preserves vessel geometry exactly and replaces only the
level. Measured human breast tumour PS: **0.06-0.10 min^-1** by DCE-MRI
distributed-parameter fits, with 0.22 as an upper bound (n=18 dogs,
spontaneous mammary tumours; down-weighted because its extraction fraction
E = 0.56 breaches the E > 0.2 regime where K_trans underestimates PS).

## 2. Result — it gets worse

Real Janesick Xenium breast section, 167,780 cells, `dx = 12 um`, 6 h,
transvascular exchange ON, `alpha` unchanged.

| arm | 5-fluorouracil | methotrexate |
|---|---:|---:|
| baseline (BBB penetrance, as shipped) | **0.8141** (1.33x) | 0.0523 (0.09x) |
| measured PS 0.06 /min | 0.1725 (0.28x) | 0.0417 (0.07x) |
| measured PS 0.10 /min | 0.2533 (0.42x) | 0.0665 (0.11x) |
| measured PS 0.22 /min (upper bound) | 0.4127 (0.68x) | 0.1299 (0.22x) |
| **MEASURED ANCHOR** | **0.61 ± 0.11** | **0.60 ± 0.20** |

**5-fluorouracil was closest at baseline and every measured permeability makes
it worse.** Methotrexate improves at best to 0.22x of measured, still off by
4.6x. Neither lands inside its band at any point in the measured range.

The prediction was written into the driver before it was run, and it held.

## 3. The diagnostic — and it is decisive

At steady state `Kp = PS/(PS + clearance + binding)`, so the PS each drug
would NEED to reach its own measured anchor is `PS = Kp*k/(1 - Kp)`:

| drug | `clearance + binding` | PS needed | vs measured breast PS |
|---|---:|---:|---:|
| 5-fluorouracil | 5.586e-3 /s | 8.737e-3 /s | **2.4x** the top of the range |
| methotrexate | 2.670e-2 /s | 4.005e-2 /s | **10.9x** the top of the range |

Measured breast PS spans 1.0e-3 to 3.67e-3 /s. **Both drugs need more
permeability than has ever been measured in breast tumour.** With a measured
permeability AND a measured anchor, the system is **over-determined and
inconsistent** — so the error cannot be in permeability.

**Now the part that identifies the culprit.** The two drugs need PS values
differing by **4.58x**. But the measurement says their breast permeabilities
differ by **1.02x**. That 4.58x has to come from somewhere, and it does:
`PS_needed` is proportional to `k = clearance + binding`, and the model's
`k` differs between these two drugs by **4.78x**. The needed-PS spread is
simply the `k` spread, inherited.

So: measured permeability is flat across these drugs (1.02x). The measured
anchor is flat across these drugs (0.61 vs 0.60, also 1.02x). A model that
reproduces both **must** have a roughly flat `k` too. This one has `k` varying
by 4.78x.

> **The remaining error is in the alpha-scaled binding sink, not in the
> permeability — and specifically in its DRUG-TO-DRUG VARIATION.**

## 4. Why `alpha` was never constrained by this

`alpha` is fit once, on GBM doxorubicin **penetration depth**, through
`L = sqrt(D_eff / k)` with `k = clearance + binding`. That fit determines `k`
from a decay length and nothing else — **PS does not appear in the expression
at all** (see `results/transvascular_exchange_revalidation.md` §6.2). So `k`
has never been constrained by any delivery-magnitude measurement, in any
tissue. It could be badly wrong for breast and the GBM chain would still pass,
exactly as it does.

This also explains a result that looked puzzling earlier: the calibration
chain survived the transvascular-exchange fix untouched, for the structural
reason that PS vanishes from the analytic penetration length. The same
structural fact is why that chain cannot police `k`.

## 5. What this changes

**Not a fix, a localisation.** The breast anchor gap is now attributed to a
specific parameter, with the two alternatives ruled out by measurement:

* delivery-law FORM — fixed (transvascular exchange), and it helped: 5-FU went
  from 12.8x to 1.33x of measured.
* PERMEABILITY — tested here against measured values, and it is **not** the
  free parameter. Substituting it makes the fit worse.
* **BINDING/CLEARANCE (`alpha`)** — the remaining candidate, and the only one
  whose drug-to-drug spread (4.78x) matches the residual (4.58x).

**Do not respond by fitting `alpha` to the breast anchor.** That would be the
third parameter absorbing the error of the first two, on two data points, and
it would break the GBM out-of-sample mitoxantrone test that `alpha` currently
passes. The honest next step is an independent measurement of tissue binding
or clearance for these drugs in breast — the same category of external anchor
that made this diagnosis possible at all.

## 6. Status

Breast transport remains **UNVALIDATED**. The gap is narrower and much better
characterised than at the start of this work: 5-FU is 1.33x from its measured
anchor under the shipped configuration, against 12.8x before the delivery law
was corrected. Methotrexate remains 11x off and is now understood to be off
because its modelled loss rate is ~4.8x too high relative to 5-FU, not because
its permeability is wrong.
