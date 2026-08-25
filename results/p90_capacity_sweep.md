# Does `E_p90/E_p10` — the headline — survive a tissue-capacity term?

**Status: PRE-DECLARED, NOT YET RUN.** Driver: `experiments/p90_capacity_sweep.py`.
Data: `results/p90_capacity_sweep.json`.

This section was written and saved **before the sweep was started**. Nothing
below it existed when the predictions and the failure thresholds were fixed. If
the thresholds in §0.3 look convenient after the numbers land, that is what
pre-declaration is for: they cannot have been chosen to fit.

---

## 0. PRE-DECLARATION

### 0.1 Why this measurement exists

`E_p90/E_p10` — the P90/P10 delivery-spread ratio over per-tumour-cell exposure
(`tmesim/validation/delivery_metrics.exposure_percentiles`) — is the product's
headline delivery number. **Its capacity stability has never been measured.**
Two of its supporting arguments failed in the two days before this file was
written:

1. `app/web/js/ui/reach.js` claimed the viable-tumour denominator carried "~10x
   grid sensitivity". Measured at 0.6 % and retracted
   (`results/reach_invariances_measure.md` §7.1).
2. `results/delivery_spread_noise.md` argued P90/P10 is safe partly because it
   is structurally blind to distance binning.
   `results/spatial_noise_decomposition.md` showed binning was never the
   mechanism, so blindness to it is not the benefit claimed. Qualified in place.

And its nearest **measured** relatives are severely capacity-sensitive: on the
same section the normalised falloff moves **0.3772** absolute and the fitted
penetration length collapses **279.715 -> 76.740 µm** (3.6x) by R = 39
(`results/reach_invariances_measure.md` §4). That file says explicitly that
nothing in it licenses any claim about `E_p90/E_p10`, and that *"the presumption
should be that the headline is too, until someone sweeps it"*. This is that
sweep.

### 0.2 What I expect, written before running

The capacity factor `R` multiplies the whole AUC field and therefore **cancels
exactly** from a percentile ratio — the same cancellation the reach ratio enjoys.
That is not the question. The question is that capacity **rescales time**, so at
fixed wall-clock window `T` a capacity-R tissue has only reached the state the
R = 1 tissue reached at `T/R`. A younger field is less equilibrated and more
concentrated around vessels, and **P90/P10 reads exactly that shape**.

So my prediction, in the direction the mechanism points:

* **`E_p90/E_p10` will INCREASE monotonically with R, substantially.** For
  temozolomide I expect the R = 39 value to be **at least 2x** its R = 1 value
  (i.e. >= ~4.2 against a shipped ~2.11), on the reasoning that the fitted
  penetration length over the same range collapses 3.6x on the same section.
* **The two compounds will CONVERGE at large R.** At R in the hundreds to
  thousands the field is a thin perivascular shell for both compounds, and a
  shell's spread is set by the section's vascular geometry rather than by the
  compound. I expect the doxorubicin/temozolomide separation ratio to fall
  monotonically and to approach ~1 at R = 4184.
* I expect the **published parameter-division emulation arm to disagree
  visibly** with the correct arm at convection ON, because it is a capacity
  change plus an R-fold Péclet boost (see §0.5).

I would be surprised — and would report it as such — if the headline turned out
flat. Nothing here is arranged to make it fail; the thresholds below are
inherited, not invented.

### 0.3 What result makes me say the headline is UNSAFE

The thresholds are **not chosen by me**. They are the metric's own previously
measured nuisance envelope, taken verbatim from `results/delivery_spread_noise.md`
§0 and §2, which measured 30 solves over 5 grid/mask reconstructions and 3
exposure windows:

| inherited yardstick | value | source |
|---|---:|---|
| temozolomide, shipped flagship config | 2.111 | `delivery_spread_noise.md` §2 |
| temozolomide, full nuisance span | 1.969 – 2.171 (width **0.202**) | §0 |
| doxorubicin, shipped flagship config | 3.942 | §2 |
| doxorubicin, full nuisance span | 3.770 – 3.996 (width **0.226**) | §0 |
| worst-case separation anywhere on the nuisance grid | **1.74x** | §0 |
| separation at the shipped config | 1.87x | §3 |

**The three criteria, fixed now:**

* **C1 — VALUE STABILITY AT THE BEST-AVAILABLE CONSTANT.** Take the move of
  `E_p90/E_p10` from R = 1 to **R = 39** (`results/brain_slice_anchor.md`'s
  best-available constant, n = 33). **C1 FAILS if that move exceeds the
  compound's own full nuisance-grid width** — 0.202 for temozolomide, 0.226 for
  doxorubicin. In words: capacity fails C1 if it moves the headline further than
  every grid, mask-repair and exposure-window choice combined.
* **C2 — SEPARATION OVER THE MEASURED CAPACITY RANGE.** **C2 FAILS if the
  doxorubicin/temozolomide separation ratio drops below 1.74x at any R in the
  measured range 3.2 – 4184 mL/g.** 1.74x is the worst-case margin
  `delivery_spread_noise.md` established over its own nuisance grid; a capacity
  term that eats into it is doing something no analysis choice could.
* **C3 — IS THE SHIPPED NUMBER EVEN REPRESENTATIVE?** Report the full span of
  `E_p90/E_p10` over R = 3.2 – 4184 for each compound and state the shipped
  value's position in it. This one has no pass/fail bar attached and is not used
  to decide anything — it is descriptive, and it is declared as descriptive here
  so it cannot be promoted to a verdict afterwards.

**The verdict rule, fixed now:**

| C1 | C2 | verdict |
|---|---|---|
| pass | pass | **the headline survives capacity** |
| fail | pass | **the headline's VALUE is capacity-sensitive; the separation is not** |
| pass | fail | **the separation is capacity-sensitive** (would be surprising) |
| fail | fail | **the headline is capacity-sensitive and so is the claim built on it** |

### 0.4 A reproduction check that is not a threshold

R = 1 on the shipping arm **is** the shipping configuration, so it must
reproduce `delivery_spread_noise.md` §2's flagship row (2.111 / 3.942) to the
digits printed there. This is a check that I am measuring the same object, not a
criterion. **If it does not reproduce, this file reports that and stops** — it
does not adjust anything to make it reproduce. No parameter will be tuned to hit
any anchor, at any point, for any reason.

### 0.5 The instrument, and the trap it has to avoid — declared before running

Two routes to a capacity-R answer exist in this repo, and **only one of them is
a capacity change when convection is on.**

* **ARM A — the truncated-window route (CORRECT).** Run the shipping R = 1 model
  at `(dt/R, T/R)`. The solver's step is `c_new = A^-1[c + h*S - h*adv(c)]` with
  `A = I - h*L + h*k` (`tmesim/physics/diffusion.py:377`, `:447`, `:452`), so
  discretising `R dC/dt = div(D grad C) - div(v C) + S - k C` at timestep `dt`
  gives exactly the R = 1 solver at `h = dt/R` — **iterate for iterate, and with
  the advection term retarded by construction**, because the same `h` multiplies
  it. The step count is unchanged (2,160), so no time resolution is traded for
  the truncation. This is the route
  `results/reach_invariances_measure.md` §3.1 derives and §3.4 validates.
* **ARM B — the published parameter-division emulation (CONTAMINATED WITH
  CONVECTION ON).** `experiments/capacity_ranking_robustness.py::retard` divides
  `diffusivity`, `binding`, `clearance` and `bbb_penetrance` by R. It **cannot
  divide the velocity field**, which is built inside
  `tmesim.sim.engine.simulate_response`. Multiply the emulated equation through
  by R and the advection term keeps its R: with convection ON this is a capacity
  change **plus an R-fold boost to the Péclet number**, not a capacity change.
  `capacity_ranking_robustness.md` §2d runs convection OFF for exactly this
  reason. **The shipping config is convection ON.**

Both arms will be run, and **every number below names its arm in the same
sentence**. Arm B's numbers are reported as what they are — a capacity change
plus an R-fold Péclet boost — and are never quoted as the capacity answer.

The two routes are additionally cross-checked at **convection OFF**, where the
emulation is exact, so that Arm A's implementation is validated against an
independently derived one rather than trusted. **The convection-ON disagreement
is measured too, and it is the planted defect of this study** (§0.6).

### 0.6 The checks added, and the defect each one is made to catch

Every check below is exercised against a deliberately planted defect
(`--stage plant`), because a check that has never fired is not evidence.

| check | the defect it claims to catch | how it is planted |
|---|---|---|
| `_assert_alpha()` | the calibration moved under the study | `calibrate_alpha` monkeypatched to return a perturbed value |
| `_assert_population_invariant` | the tumour-cell population moved between R values, making percentiles incomparable | a row's `n_tumor` corrupted |
| `_assert_steps_invariant` | the truncation traded time resolution for window length | a row's `n_steps` corrupted |
| `_assert_ratio_defined` | `E_p10` underflowed to 0 at extreme R and the ratio silently became `inf`/`nan` | a row's `E_p10_uM` set to 0 |
| `_routes_agree` | **the trap: using Arm B as if it were a capacity change with convection ON** | the identity stage feeds it the convection-**ON** emulation and asserts the check FAILS, then the convection-OFF one and asserts it PASSES |

### 0.7 Configuration, stated once and repeated with every number

| | |
|---|---|
| section | flagship 10x Xenium GBM, `reconstruct_real(source='xenium')` via `vessel_profile.build_arm_microenvironment` |
| grid | dx = **12 µm**, mask **repair = 2** — the `delivery_spread_noise.md` **flagship** row |
| timestep | dt = **10 s** at R = 1 (`app/server/jobs.PDE_DT_S`); **dt/R** on Arm A |
| window | **21,600 s = 6 h** at R = 1 (`app/server/api.STUDIO_DEFAULT_TOTAL_TIME_S`); **T/R** on Arm A |
| steps | **2,160 in every run of both arms**, asserted |
| convection | **ON** (the shipping value), except in the identity cross-check which says so |
| compounds | temozolomide and doxorubicin, `calibrated_drugparams` -> `apply_barrier_to_drugparams` |
| PK / oxygen / seed | `PKConfig()` defaults; `OxygenConfig()` defaults; seed 0 |
| R ladder | 1, 3.2, 8, 20, 39, 60, 80, 120, 620, 4184 — the same ladder the sibling capacity files use, spanning the measured 3.2 – 4184 mL/g (n = 33, best-available constant 39, `results/brain_slice_anchor.md`) |

`calibrate_alpha()` is asserted live at `22.223477711172563` at the **top and
bottom of every stage**. Nothing under `tmesim/` or `app/` is written to.

---

*(Everything below this line is written after the sweep. The pre-declaration
above is frozen.)*
