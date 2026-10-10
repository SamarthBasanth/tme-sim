# Can a learned surrogate stand in for the transport solver? — NO on the two readouts the product defends

**Status: RUN AND MEASURED.** Drivers: `experiments/surrogate_feasibility/`
(`sweep.py`, `ladder_sweep.py`, `coarse_baseline.py`, `fit.py`, `analyse.py`,
`probe.py`, `probe_axes.py`). Data: `experiments/surrogate_feasibility/
{sweep_data,ladder_data,coarse_baseline,fit_results,analyse_results}.json` and
the `.log` beside each. Guards: `tests/test_surrogate_feasibility.py` (19
assertions; each class was confirmed to FAIL against a planted defect before
being kept — see §9.2).

**Nothing shipped was changed.** No file under `tmesim/`, `app/`, `docs/` or
`landing-3d/` was modified; no default, constant or threshold was moved.
`calibrate_alpha()` is asserted live at `22.223477711172563` at the start and
end of every stage and holds in all of them. Nothing was fitted to a
measurement of tissue. Every row of training data is a real
`app.server.jobs.execute_run` call — the same function the Studio queue calls.

**This is a feasibility study, not a feature.** No AI/ML claim reaches any
product surface from this work: nothing under `app/` was touched, and the
prohibition this repository enforces in `tests/test_home_page.py`
(`"machine learning"` on a user-visible surface is a CI failure) is untouched
and should stay.

---

## 0. Verdict up front

**A learned surrogate does not reproduce this solver's two defended readouts
well enough to show anyone, and the thing it would have to beat is not "no
preview" — it is a COARSER SOLVE, which beats it by a factor of ~5 and cannot
go stale.**

On the four sections where this repository holds per-cell spatial statistics
VALID, held out one donor at a time:

| | delivery spread `E_p90/E_p10` | fitted penetration length |
|---|---:|---:|
| across-donor range of the readout itself | **0.5357238964376265 x** (1.9936 – 2.5293) | **174.5058718458447 um** (103.18 – 277.69) |
| no model (train-set central value) | MAE **0.6899** | MAE **43.06 um** |
| ridge | MAE **0.4375**, p90 0.959 | MAE **32.89 um**, p90 88.24 |
| gradient boosting | MAE 0.5851, p90 1.194 | MAE 33.12 um, p90 145.6 |
| **the same solver on a 2.5x coarser grid** | **max abs err 0.0909**, median 0.0248 | **max abs err 17.32 um**, median 1.29 |

The best learned model lands **±0.4375 on a readout whose entire across-donor
range is 0.5357238964376265** — the error is 82 % of the full dynamic range.
The coarse re-solve lands inside **±0.0909**, which is 17 % of that range and
45 % of the width of the readout's own published nuisance band (1.969 – 2.171,
`results/delivery_spread_noise.md` §0). It is ~4.8x better than the best
model's MEAN error, it is the engine's own output so it carries the engine's
caveats, and it has no training set to go stale.

**Three secondary findings, each of which matters more than the model scores:**

1. **Three of the four PK knobs a user can turn are EXACTLY null for every
   shape readout** (§3). Doubling or quadrupling `plasma_peak` leaves
   `E_p90/E_p10` bit-identical; a 16x change in `half_life_s` moves it by
   `1.0000000000001659`; the whole reach ladder is bit-identical under a
   simultaneous 4x dose, 4x half-life and 6x window change (0 violations in 4
   checked pairs). A dose slider needs an exact scaling identity, not a model.
2. **The honest split is worth the entire apparent result** (§5). On
   `E_p90/E_p10` over all 921 rows a random 5-fold split reports gradient
   boosting at R² **0.940**, skill **+0.885**. Leave-one-donor-out on the same
   rows and the same model reports R² **−0.026**, skill **−0.022** — worse than
   predicting the training mean. Reporting the first number would have been a
   straightforward way to ship nothing.
3. **The one target that does transfer is the one this repository does not
   defend** (§5.3). Dose-normalised exposure magnitude is predicted across
   unseen donors to a median factor of **1.423x** (skill +0.745, R² 0.914),
   against 10.850x for the trivial baseline. That is a real result. It is also
   an ABSOLUTE exposure level, and this repository's position is that the
   delivery law's absolute levels are not defensible
   (`app/server/jobs.py`: tissue-over-plasma above 1 "IS THE LAW'S KNOWN
   UNBOUNDEDNESS MADE VISIBLE"). A surrogate that is accurate only about the
   quantity the product declines to stand behind is not a shippable preview.

**What a well-measured "no" is worth here:** the preview feature the product
owner wants is buildable today, without a model, by solving on a coarser grid —
and §6 measures what that costs and what it gets wrong.

---

## 1. The target, and why that one

**Scalar readouts, not the concentration field.** Three reasons, in order of
weight:

* The scalars are what the product prints and what its claims are made of.
  `app/server/reach.py::CONVENTION_NOTE` names the delivery spread
  (`P90/P10` of tumour-cell exposure) as the headline, in place of `F(0.1)`.
* A field surrogate's output is 424 x 675 = 286,200 numbers per run at the
  shipped `dx = 12 um`. Fitting that from 20 donors is not a sample-size
  problem one can caveat; it is the wrong shape of problem.
* A preview shows a number and a map. The number is the part a learned model
  could plausibly supply without the solve.

**The three targets, and why each is defensible to predict:**

| target | unit | why it is allowed |
|---|---|---|
| `E_p90_over_p10` | dimensionless | the product's own headline; a RATIO, so dose-invariant and scale-free |
| `penetration_length_um` | um | a transport length a user reads directly; the readout with the most across-donor headroom (2.691x) |
| `median_tumor_auc / plasma_peak` | s | dose-normalised, so it is a property of the tissue-and-compound rather than of the prescription |

**What was deliberately NOT made a target.** Absolute `kill_fraction` — never a
defensible output, and `tests/test_surrogate_feasibility.py` asserts no target
name contains `kill`. Nothing that implies an ordering between compounds; the
cross-compound ranking was withdrawn 2026-08-21 and compound descriptors enter
only as model INPUTS. No biological or clinical quantity: the ground truth here
is the solver, and the solver alone.

**Features are solve-free by construction.** The cached reconstruction's
summary statistics (`real_pipeline.py:1319-1329` builds
`mean_vascular_density`, `mean_ecm_density`, `mean_oxygen`,
`tissue_coverage_fraction`; `cache.py:1484` stores `max_darcy_speed_um_s` in
the cache entry itself), the compound's descriptors, the RunSpec's PK numbers,
and physics-motivated groups of those (`sqrt(D/k)`, `sqrt(D t)`, and the
closed-form `penetration_length_um(dp, ecm_density)` the product already
computes free on every run). 37 columns. Two guard tests assert that no name
from the solver's output blocks appears among them, in the declared lists and
in the built matrix.

---

## 2. The real n

**38 sections, 20 donors, 921 solver runs, 0 failures — and 4 sections where
the readout is defensible.**

| | count |
|---|---:|
| sections the backend can reconstruct | **38** (7 built-in + 31 ingested) |
| donor groups | **20** |
| solver runs executed, `status == "ok"` | **921 / 921** |
| rows on sections where `spatial_spread` is VALID | **99** |
| sections where `spatial_spread` is VALID | **4** — `xenium_319`, `xenium_io`, `gse324268`, `xenium_breast` |

**The 4 is the number that governs everything.** On the other 34 sections
`app.server.api.SPOT_SPATIAL_REFUSAL` declares the delivery map, its falloff,
its perivascular reference and the delivery spread INVALID as readings of that
tissue, and records that `penetration_length_um` moves by up to 19.4 % under a
pure change of the pseudo-cell placement seed. Those rows are kept and used,
deliberately: a surrogate's ground truth is the SOLVER, the solver produces a
number there, and reproducing it is a well-posed question. But no statement
about tissue is built on them, and every model score is reported twice — once
over all 921 rows and once over the 99 where the readout is defended.

**Grouping is by DONOR, not by section.** GBM-Space ships up to five sections
per patient (`AT10-BRA-5-FO-1_1`, `...-1_2`, `...-2_1`, ...). Leaving one out
while training on a sibling is a leak. 38 sections collapse to 20 donors; on
the valid subset each section is its own donor, so the section and donor splits
coincide there and are reported once.

**The effective n is smaller than 921.** Of the four PK shapes, two differ from
the shipped one only in `half_life_s`, which §3 measures to be an exactly null
axis for both shape targets. For `E_p90/E_p10` and `penetration_length_um` the
design is therefore 38 sections x 6 compounds x 2 windows = **456 distinct
cells**, and 465 of the 921 rows are duplicates of another row in those
targets. They were run anyway because `median_auc_over_peak` DOES move with
half-life (3.717x), and they are not counted twice in any claim.

**An incidental observation, offered as data rather than as a complaint.**
`tmesim/ml/panel_imputation.py` is this repository's existing learned component
and its guardrails are the right template (OPT-IN only, a provenance tier
`kind='ml_imputed'` that is "neither `measured` nor `computed` nor `assumed`",
`_ML_HYPOXIA_VALIDATED_PANELS = ("brain_319",)`, and "Every validation number
in this module is computed on a SECTION the model was not trained on"). Its
same-assay headline, `R² = 0.692` leave-one-section-out, is measured on
`AT10-BRA-5-FO-{2_1,3_1,4_2}` — three sections of **one patient**, AT10. The
file is straightforward about which number carries the weight: it names the
cross-patient cross-panel figure (`R² = 0.238`) as "the more relevant
cross-panel check" and says it "didn't clear a strict one". That is the
precedent this study followed in leading with the cross-donor number.

---

## 3. Three of the four user-facing knobs are EXACTLY null

Measured by grouping rows that are identical except along one axis and taking
the largest ratio the axis produces in each target.

| axis | range swept | comparison groups | `E_p90/E_p10` max fold change | `penetration_index` | reach `F(0.1)` |
|---|---|---:|---:|---:|---:|
| `plasma_peak` | 50 → 200 uM (4x) | 9 | **1.0** | **1.0** | **1.0** |
| `half_life_s` | 3600 → 57600 s (16x) | 228 | **1.0000000000001659** | 1.0000000000000817 | **1.0** |
| `total_time_s` | 3600 → 21600 s (6x) | 228 | 3.7169188850803536 | 1.3162582817734663 | **1.0** |

`E_p50_uM` moves by **exactly 4.0** under the 4x dose change and by
`3.7170896085598355` under the half-life change, so the field is responding —
it is the RATIOS that do not. The transport chain is exactly linear in the
plasma amplitude, so every relative readout is an exact algebraic invariant of
the dose.

`penetration_length_um` is the one near-exception, and only through its fit:
max fold change `1.0003492444001745` under the dose change and
`1.0012726288714509` under the half-life change — fourth-decimal wobble in an
`exp + floor` least-squares fit, not a physical response.

**What follows for the product.** The use case that motivated this study — a
user slides a dose or a schedule and wants the readout now — is answered for
three of the four knobs by an identity, at zero error and zero cost, and for
the fourth (`total_time_s`) by the solve. A model trained across those axes
would be learning to reproduce an invariance it is handed, and every duplicated
row would shrink its reported error without it having learned anything. That is
why the dose axis is nine confirmation rows in this design and not a sweep.

---

## 4. The reach ladder is a poor target, for two independent reasons

**It is invariant to every PK input.** `app/server/reach.py` computes the
shipped map at steady state (`window_mode = "steady_state"` on every row), and
a steady-state self-normalised ratio cannot depend on the schedule. Measured:
`reach_headline_F` max fold change **exactly 1.0** across all three axes above
(n = 149 groups for half-life and window, n = 2 for dose). `ladder_sweep.py`
re-asserts this on the full ladder rather than assuming it: 4 pairs, a
simultaneous 4x dose / 4x half-life / 6x window change, on one Xenium and one
Visium section — **0 violations**, headline and all seven rungs bit-identical
(`xenium_319` temozolomide stays
`[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42050520059435365]`).

**It is substantially degenerate, and refused outright where it matters most.**
Of 921 rows: **96** carry a `null_reason` and have no ladder at all; of the 825
finite headline values, **227 (27.5 %)** are exactly `0.0`; the median is
`0.017386155667994403`. The refusals are the engine's own, recorded verbatim —
18 of them "the perivascular reference exposure E_ref is nan", and 6 on
`xenium_io` because there is "only 1 VIABLE tumour cell in this reconstruction,
below the 30 this repository requires". So on the defended subset the ladder
exists for **3 of the 4** sections.

A model scored on a target that is 27.5 % a single constant, PK-invariant by
construction, and undefined on a quarter of the defended sections would be
reporting its own degeneracy. No ladder model was fitted.

---

## 5. Held-out error, in the units a user would read

### 5.1 The split is the whole experiment

`E_p90_over_p10`, all 921 rows, gradient boosting, same features, same
hyperparameters — only the split changes:

| split | folds | R² | skill vs the mean baseline |
|---|---:|---:|---:|
| random 5-fold over rows | 5 | **0.940** | **+0.885** |
| leave-one-section-out | 38 | 0.048 | −0.015 |
| **leave-one-donor-out** | 20 | **−0.026** | **−0.022** |

The random split is the wrong split and is computed only to show the size of
the flattery. Across sections and across donors, neither learner beats the
training mean.

### 5.2 On the sections where the readout is defended

Leave-one-donor-out over the 99 rows on 4 valid donors (the section and donor
splits coincide here). Errors are in the target's own unit, after inverting the
log fit the models were trained in.

**`E_p90_over_p10`** — across-donor range **0.5357238964376265** (1.9936 /
2.3211 / 2.3787 / 2.5293; fold spread **1.2687256241895875x**):

| | MAE | median AE | p90 AE | R² | skill |
|---|---:|---:|---:|---:|---:|
| no model | 0.6899 | 0.5933 | 1.659 | −0.116 | 0.000 |
| **ridge** | **0.4375** | 0.3375 | 0.959 | 0.590 | +0.366 |
| gradient boosting | 0.5851 | 0.512 | 1.194 | 0.286 | +0.152 |

**`penetration_length_um`** — across-donor range **174.5058718458447 um**
(277.69 / 133.16 / 125.56 / 103.18; fold spread 2.691214367086715x):

| | MAE | median AE | p90 AE | R² | skill |
|---|---:|---:|---:|---:|---:|
| no model | 43.06 | 29.51 | 73.74 | −0.122 | 0.000 |
| the shipped closed form | 54.2 | 28.59 | 143.1 | −0.411 | **−0.259** |
| **ridge** | **32.89** | 13.59 | 88.24 | 0.226 | +0.236 |
| gradient boosting | 33.12 | 13.36 | 145.6 | 0.097 | +0.231 |

Read the first row of each table against the third. **A ±0.44 error on a
readout whose entire across-donor range is 0.54 cannot tell a user which donor
they are looking at.** It is also larger than the whole published nuisance band
of that statistic on the flagship (0.202 wide). The penetration-length model is
better than nothing and better than the existing closed form — MAE 32.89 um
against 43.06 for the mean and 54.2 for the analytic formula — but its error is
19 % of the across-donor range and its p90 error (88 um) is half of it.

**Did a simple model suffice? Yes, and that is not good news.** Ridge beats
gradient boosting on BOTH defended targets on the honest split (0.4375 vs
0.5851; 32.89 vs 33.12 with a far better p90). The 719 kB boosted model buys
nothing a 3 kB linear model does not, and on the full 921-row donor split ridge
is numerically unstable: it overflows to MAE `8.094e+11` on `E_p90_over_p10`
and `633.2` on penetration length, reported here rather than quietly dropped.

### 5.3 The one target that transfers, and why it does not rescue the case

`median_auc_over_peak`, leave-one-donor-out, all 921 rows, 20 folds:

| | MAE (s) | median fold error | R² | skill |
|---|---:|---:|---:|---:|
| no model | 1.943e+05 | **10.850x** | −0.119 | 0.000 |
| ridge | 2.024e+05 | 1.536x | −3.408 | −0.042 |
| **gradient boosting** | **4.946e+04** | **1.423x** | **0.914** | **+0.745** |

On the 4 defended donors: median fold error 1.663x, skill +0.628. Predicting a
tissue's dose-normalised exposure magnitude to within ~1.4x on a donor never
seen is a genuine, reproducible result, and it is the only one in this study.

It does not make a shippable preview. The quantity is an ABSOLUTE exposure
level, and this repository's own position on absolute levels from the delivery
law is that they are not defensible: `app/server/jobs.py` states that
tissue-median-over-plasma-peak above 1 "IS THE LAW'S KNOWN UNBOUNDEDNESS MADE
VISIBLE, never the tissue concentrating the drug", and
`results/CLAIM_INVENTORY.md` §2 lists "Any absolute number" among what a
reviewer would not accept. A model that is accurate about the number the
product declines to stand behind, and no better than a constant about the two
numbers it does stand behind, is the wrong way round.

### 5.4 Does the model learn the tissue, or only the compound?

Feature-group ablation, leave-one-donor-out, skill vs the mean:

| features | n | `E_p90/E_p10` (hgb) | `penetration_length_um` (hgb) | `median_auc_over_peak` (hgb) |
|---|---:|---:|---:|---:|
| compound + PK only | 19 | −0.6434 | −0.3167 | +0.0047 |
| tissue + PK only | 26 | +0.0086 | −0.0706 | +0.7103 |
| PK only | 8 | −0.0008 | +0.0552 | +0.1809 |
| all | 37 | −0.0223 | −0.1073 | +0.7454 |

For the two defended targets no feature group reaches useful skill on a new
donor; the best cell in those columns is `+0.0086`. For the magnitude target
essentially all the skill is in the TISSUE columns (+0.7103 of +0.7454), which
is the right place for it to be and is consistent with that result being real.

### 5.5 A new compound on a seen tissue

Leave-one-compound-out on the 99 defended rows: `E_p90/E_p10` best MAE
**0.54** (hgb, skill +0.338) — equal to the entire across-donor range;
`penetration_length_um` best MAE **28.78 um** (hgb, skill +0.419), against a
within-split range of 77.44 um. Bring-your-own-SMILES is a shipped surface
(`results/bring_your_own_compound.md`), and a surrogate would serve it no
better than it serves a new donor.

---

## 6. The baseline that actually competes: the same solver, coarser

A surrogate is one way to show something now and solve properly behind it. The
other is to run the SAME solver on a coarser grid. That option has three
properties no learned model has: its output is the engine's, so it carries the
engine's caveats and refusals and cannot disagree with the real run about what
is defensible; it has no training set, so it cannot be stale and no new
section, compound or engine change invalidates it; and it makes no AI claim.

**What "agrees" means here.** NOT "approximates the true answer".
`results/CLAIM_INVENTORY.md` §1f RETRACTS claim N2: `dx = 12 um` is NOT grid
converged, the GBM ladder is non-monotone with observed orders 0.96 and −3.33,
and "neither panel is in an asymptotic regime". There is no converged answer
for either resolution to approximate. The question a preview has to answer is
narrower and needs no such thing: **does the number shown now match the number
the product will print when the real run finishes?**

Measured on the 2 sections where both resolutions are defensible and where `dx`
actually changes the grid — `xenium_319` (170x271 → 424x675) and
`xenium_breast` (184x252 → 458x628) — 6 compounds each, 12 pairs:

| readout | median rel err | p90 | max | median abs err | max abs err |
|---|---:|---:|---:|---:|---:|
| `E_p90_over_p10` | 0.7133 % | 1.2001 % | **1.762 %** | 0.0248 | **0.0909** |
| `penetration_index` | 0.2242 % | 0.4994 % | 1.0409 % | 0.000741 | 0.003268 |
| `median_auc_over_peak` | 0.1229 % | 0.3071 % | 0.3215 % | 22.671 s | 1039.66 s |
| `penetration_length_um` | 1.8053 % | 14.5047 % | **16.0405 %** | 1.293 um | **17.316 um** |

The headline delivery spread is reproduced to under 2 % everywhere and to
within **±0.0909 absolute** at worst — against a 0.5357238964376265
across-donor range and a 0.202-wide published nuisance band. This is consistent with, and extends to a
second section and six compounds, the figure `app/server/reach.py` already
publishes for one section and one compound (the spread "moves 1.0004x between
dx = 30 and dx = 12"); this study's `xenium_319` temozolomide pair reproduces
that section's shipped values bit-identically
(`1.9935720609199155` / `277.68964040999685` at dx = 30,
`1.9944090080526955` / `278.4919817884945` at dx = 12).

**The fitted penetration length is the exception and should be said plainly:**
up to 16 % off, 17.3 um absolute. That is the same order as the 19.4 % this
repository already records for that readout under a pure placement-seed change
on Visium sections. A coarse-grid preview of the delivery spread is sound; a
coarse-grid preview of the penetration length is not, and neither is a learned
one.

**On Visium sections a coarse preview is not even available, and this study
nearly mis-measured that.** See §9.1.

---

## 7. Cost

**A performance claim is a claim.** `results/performance_profile_2026-09-20.md`
§0.1–0.2 records one RunSpec on unchanged code measuring FIVE different
wall-clocks on this machine (47.1 s to 174.3 s at `dx = 12`, loadavg 4.6–48 on
10 CPUs) and concludes that what is stable is the SHAPE, not the absolutes.
Every figure below therefore carries its load, and the ratios are the part to
quote.

**Generating the training set.** 921 runs at `dx = 30 um` took **56.0 min**
(median run 2.19 s, p90 5.77 s, max 54.25 s; reconstruction median 1.79 s). The
reach-ladder collection added 228 runs in 17.1 min. At the resolution the
product actually ships, `dx = 12 um`, the measured cost ratio on the two Xenium
sections is **6.780x** on medians (58.77 s vs 8.67 s, loadavg 2.4–9.8 on 10
CPUs), so the same design at the shipped resolution is roughly **6 hours** —
and that understates it, because the two 850k-cell sections scale worst. For
calibration, `results/verdict_inversion_census.md` spent **25,502 s (7.1 h)**
on one section and 12 compounds.

**Training and serving.**

| | ridge | gradient boosting |
|---|---:|---:|
| train on all rows | 0.009 – 0.012 s | 0.878 – 1.241 s |
| serialised size | **3,012 B** | 719,165 – 719,693 B |
| single-row prediction, median | **0.419 – 0.425 ms** | 3.87 – 4.05 ms |

Inference is free on any scale that matters. The cost is entirely in the
training set and in keeping it true.

**What it would take to keep it current, and why this is the decisive cost.**
The surrogate is bound to one engine state: `alpha = 22.223477711172563`, the
one-way delivery law, `PDE_DT_S = 10.0`, convection ON, the steady-state reach
convention. Any change to those invalidates it while it continues to return
numbers. This repository already knows that failure mode and has built
machinery for it: `app/server/jobs.py` stamps `CLAIM_SET_VERSION` into every
run's provenance precisely so that a STORED record **banners itself** when a
claim is later withdrawn, and `tests/test_claim_set.py` parses every source
file to stop a withdrawn sentence being written into a fresh record.

A surrogate has no equivalent and is strictly worse than a stale stored run: a
stale run announces its age, whereas a surrogate prediction is generated fresh
and would carry the CURRENT claim-set version while having been trained on an
older engine. The withdrawal rate is not hypothetical — the grid-convergence
headline and the cross-compound ranking went in one week, and
`results/verdict_inversion_census.md` measures the dominant limiter INVERTING
for **7 of 12 compounds (58.3 %)** under a transvascular term the transport law
does not have, temozolomide among them. Shipping this would mean a ~6-hour
retrain gated on every engine change, plus a refusal path for the window
between a change and its retrain, plus a new claim surface of exactly the kind
`results/` already has no CI for.

---

## 8. What would have to be true to ship a surrogate

Written as falsifiable conditions, none of which this study meets:

1. **Beat the coarse solve on a defended readout, out of sample, on a new
   donor.** The bar is `±0.09` absolute on `E_p90/E_p10`. The best model
   measured here is `±0.44`, a factor of 4.8 away, and it is the simple model.
2. **Have enough donors to establish it.** 4 defended sections, 4 folds. This
   repository's own standard for a learned component, in
   `tmesim/ml/panel_imputation.py`, rejects a deep stack because "three
   independent tissue sections is not enough data to fit or honestly validate a
   network". Four is not materially more.
3. **Beat the cheap non-learned alternative.** The governing precedent is
   `results/subgrid_closure.md`: "a learned model that does not beat a
   two-parameter fit is not a result", sharpened there to a zero-parameter
   formula, with the conclusion "Nothing was trained, and training would not
   have been justified". The same conclusion holds here, against a
   zero-parameter alternative that is the solver itself.
4. **Target something the product defends.** The only target that transfers is
   an absolute exposure level (§5.3).
5. **Carry a staleness guard with teeth.** A version stamp on the training set
   plus a refusal when it does not match the live engine, in the shape of
   `CLAIM_SET_VERSION`. Buildable, but it converts the feature into a
   ~6-hour-per-engine-change maintenance obligation.

**The recommendation, which needs no model.** Build the preview as a coarse
solve, on Xenium sections only, for the delivery spread and the penetration
index but NOT the fitted penetration length, with the measured agreement in §6
shown beside it. Make the dose and half-life response exact algebraic scaling
(§3) rather than a re-solve or a prediction. And the unclaimed slot worth more
than any of this: there is no RunSpec-keyed RUN cache — only the reconstruction
cache (`API_CONTRACT.md:845`) — so a repeated spec re-solves from scratch
today. An exact memo on a frozen RunSpec plus the server-side solver constants
would make a repeat instant with zero error, and `jobs.py:86-97` already names
the caveat (the RunSpec alone does not pin the solver settings, so they must
enter the key).

---

## 9. Defects found in this study's own method

Both are recorded because a study that reports only its subject's problems is
not reporting its own.

### 9.1 A comparison that passed by comparing nothing

`coarse_baseline.py` first reported `median_rel_err = 0.000 %` on all four
readouts across 30 pairs — an apparently perfect result. It was vacuous. On the
three Visium sections in that set the reconstruction pins the grid to the spot
lattice, so `dx = 30` and `dx = 12` produce the **identical grid**
(`[54, 54]`, `[54, 54]`, `[54, 55]`) and 18 of the 30 pairs compared a run
against itself. The aggregate mixed 12 real comparisons with 18 tautologies and
the median landed on a tautology.

The fix was to restrict the comparison to pairs where the grid actually
differs, which is what §6 reports (12 pairs, 2 sections). The by-product is a
real finding: **on a Visium section `dx_um` is not a lever at all**, so a
coarse-grid preview does not exist there, and `dx = 30` and `dx = 12` give
bit-identical answers.

This is the `results/codon_vacuous_green_checks` pattern — a check that passes
without testing its subject — reproduced inside a study written to look for it.

### 9.2 A harvester that read a field that did not exist

The first pass recorded 822 empty reach ladders with no error anywhere.
`solver_runner.py` read `reach["ladder_x"]` / `reach["ladder_under_reached_
fraction"]`, but `app/server/reach.py` nests the ladder under
`reach["curve"]["x"]` / `["under_reached_fraction"]`. `.get()` on a wrong key
returns `None`, and a column of `None` is indistinguishable from a quantity the
solver declined to produce. Fixed, re-collected by `ladder_sweep.py`, and now
guarded by two tests — one asserting the SHAPE against `app.server.reach.
null_reach()` rather than against a copy of it, one asserting the source no
longer contains the flat keys.

**Every guard class was confirmed to bite.** Three defects were planted and the
suite re-run: `donor_of` returning the section id (2 failures), a solve output
added to `TISSUE_FEATURES` (2 failures), and a fold generator training on all
rows (7 failures). All reverted; 19 tests pass clean.

---

## 10. Reproduce

```sh
export TMESIM_STATE_DIR=experiments/surrogate_feasibility/_state
./.venv/bin/python experiments/surrogate_feasibility/probe.py          # dose invariance
./.venv/bin/python experiments/surrogate_feasibility/probe_axes.py xenium_319 12
./.venv/bin/python experiments/surrogate_feasibility/sweep.py          # 921 runs, 56 min
./.venv/bin/python experiments/surrogate_feasibility/ladder_sweep.py   # 228 runs, 17 min
./.venv/bin/python experiments/surrogate_feasibility/coarse_baseline.py
./.venv/bin/python experiments/surrogate_feasibility/analyse.py
./.venv/bin/python experiments/surrogate_feasibility/fit.py
./.venv/bin/python experiments/surrogate_feasibility/derived.py        # sec. 6 + the valid-only ranges
./.venv/bin/python experiments/surrogate_feasibility/digest.py         # every quotable number
./.venv/bin/python -m pytest tests/test_surrogate_feasibility.py
./.venv/bin/python experiments/surrogate_feasibility/check_report.py \
    results/surrogate_feasibility_2026-10-09.md
```

`_state` holds a copy of the ingested-section registry so the 31 ingested
sections resolve; the built-in sections need `data_cache/` populated as
`vivum sections` describes. `check_report.py` re-verifies this file: no
placeholder survives, and every decimal of 4+ significant digits in it appears
verbatim in one of the data files or in the repository record it cites.

**The claim this file makes, in one sentence:** a learned surrogate of this
solver is not worth shipping, because on the two readouts the product defends
it does not beat a coarser run of the solver itself — and on the one readout it
does predict well, the product does not defend the answer.
