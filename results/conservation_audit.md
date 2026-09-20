# The transport solver's five invariants, tested as properties — the nulls are exact, the solver's own CFL guard is not sufficient for positivity, and one monotonicity claim was never true as stated

**Status: RUN.** 2026-09-20. Driver: `experiments/conservation_audit.py`.
Property tests: `tests/test_conservation_invariants.py` (76 tests).

**Nothing was changed.** No default, constant or shipped parameter. No file in
`tmesim/` or `app/` was touched. `calibrate_alpha()` asserted live at
**22.223477711172563** at the start of the run and again at the end, and asserted
a third time inside the determinism worker's fresh subprocess. Three files were
added and nothing else: this document, the driver, and the test file.

**Deviation from the usual convention, stated rather than hidden.** Every other
`results/*.md` in this repo is paired with a `results/*.json` frozen record. This
one is not, because the work was scoped to exactly three new files and a fourth
would have been outside that scope. Every number below is reproduced from one
run of the driver, which prints all of them; the reproduction command is in
sec. 8.

---

## 0. The answer

Five invariant families, twenty-one pre-registered expectations, measured over a
class of inputs rather than at a point: a ladder of drugs / grids / timesteps /
source profiles / convection settings, a Peclet ladder spanning four decades,
all 12 compiled GBM compounds on all 5 loadable real sections, and four
independent `PYTHONHASHSEED` values in fresh subprocesses. A sixth stage plants
a defect against every threshold in the table and confirms each guard fires
(sec. 6): **15 of 15 plants reached their guard**, and the two that came back
green are reported as findings about the guards rather than as passes.

| # | Invariant | Pre-registered threshold | Measured | Verdict |
|---|---|---:|---:|:--|
| **A1** | discrete mass budget closes (right-endpoint loss) | < **1e-10** rel | **5.291e-13** worst of 26 | **PASS** |
| **A2** | the *trapezoidal* budget does NOT close (negative control) | must exceed **1e-5** | **1.429e-04 … 4.137e-03**, exactly O(dt) | **CONFIRMED** |
| **A3** | the solver's zero-clamp never fires on a physical input | exactly **0** | **0.000e+00** on all 26 | **PASS** |
| **B1** | no undershoot at any Peclet number, smooth flow | ≥ **-1e-9** rel | **0.000e+00** on all 16 rungs | **PASS** |
| **B2** | the shipped CFL guard is sufficient for positivity | < **-1e-3** ⇒ violated | **-1.000000** at the guard's own limit | **VIOLATED** |
| **B3** | the derived strict bound holds in production | ≤ **0.5** | **0.6506** at dx = 12 µm, dt = 10 s | **MARGIN LOST**, 1 of 14 |
| **C1** | more dose never less exposure anywhere; exposure exactly ∝ dose | < **1e-12** rel | **5.458e-15**; ordering exact on 55/55 | **PASS** |
| **C2** | a longer window never less AUC anywhere | ≥ **0** | **0.000e+00** (grid-wide and tissue-only) | **PASS** |
| **C3** | higher diffusivity never gives less anywhere | 0 cells decrease | **45.59 %** of tissue decreases | **FALSE AS STATED** (pre-registered as such) |
| **C3R** | …but monotone beyond a crossover r\* | 0 beyond r\* | 0 beyond r\*, **but only 0.12 % of tissue is beyond r\*** | **VACUOUS on 55/55** |
| **C3a** | in single-vessel geometry the crossing is single and clean | 1 crossing, > 50 % beyond | **1 crossing**, **69–94 %** beyond r\* | **PASS** |
| **C3b** | that repair transfers to real tissue | far band decreases → 0 | far band still **4.4 %** (dx 80) / **8.4 %** (dx 30) | **DOES NOT TRANSFER** |
| **C4** | kill fraction non-decreasing in dose at a fixed seed | ≥ **0** | smallest rung **+9.659e-04**, all 12 compounds | **PASS** |
| **D1** | bit-identical within a process | exact | identical | **PASS** |
| **D2** | bit-identical across restarts, independent of `PYTHONHASHSEED` | exact | **10/10** objects identical over 4 seeds | **PASS** |
| **D4** | `apply_drug_kill`'s default RNG is unseeded | non-zero spread | spread **4.074e-03** (**6.58 %** of the value), and **different every run** | **DEFECT CONFIRMED**, not reachable via the engine |
| **E1** | exposure exactly ∝ any scalar source factor | < **1e-12** | **3.800e-15** over a 1e6 span | **PASS** |
| **E2** | superposition | < **1e-12** | **2.367e-15** | **PASS** |
| **E3** | the capacity/time-rescale identity | < **1e-12** | **0.000e+00** (bit-identical) to **3.757e-15** | **PASS** |
| **E4** | the steady state is *exactly* capacity-invariant | < **1e-12** | **4.869e-13** worst over a 1e6 span | **PASS** |
| **E5** | the march's fixed point is the steady solve | < **1e-12** | **9.959e-15** at T ≥ 24 h | **PASS** |
| **P** | every guard above fires when its invariant is broken | all plants red | **15 / 15** plants reached their guard | **PASS** (sec. 6) |

> **THE FINDING is B2, and it is a real numerical defect.**
>
> `diffusion.simulate_drug_transport` guards its explicit advection step with
> `dt*max(|vx| + |vy|)/dx <= 1` — a maximum, over cells, of the **pointwise
> sum**. The upwind stencil's positivity condition is
> `dt*(max|vx| + max|vy|)/dx <= 1/2`, with the two maxima taken **independently
> over the grid**, because one cell can lose drug through all four of its faces
> at once while no single cell carries both maxima. The guard is therefore up to
> **2× too loose**, and on the minimal field that realises the factor exactly the
> scheme produces
>
> ```
> min(C_raw) / max(C) = -(2·CFL - 1)        for CFL > 1/2
> ```
>
> reproduced to **4.000e-07** over the three `D → 0` configurations — grids
> 9×9, 11×13 and 21×21, spacings 15, 7.5 and 25 µm, velocities 30, 2 and
> 100 µm/s. At `CFL = 1` — the guard's *own accepted
> limit* — that is **-100 % of the peak**, after which
> `np.maximum(c_new, 0.0)` silently puts the mass back: the clamp creates
> **exactly as much drug as was injected**. This is not a -1e-9 round-off
> overshoot. The minimal repro is 9 lines and is in sec. 2.2.
>
> **It is not reached by any run this repo makes.** Measured on all five real
> sections: nothing undershoots, at either timestep, on any grid, ever
> (**0.000e+00**). But at the app's own `PDE_DT_S = 10.0` s on the dx = 12 µm
> grid the strict bound reads **0.6506** — the *margin* is gone, and the only
> thing standing between the shipped configuration and the violation is that
> real Darcy velocity fields happen not to be adversarial. That is luck, not a
> guarantee, and it is the honest statement of B3.
>
> **THE SECOND FINDING is C3, and it is a claim that was never true.** "Higher
> diffusivity must not make a far-from-vessel cell see less" cannot hold
> everywhere, and that is a theorem, not a measurement: at the fixed point
> `sum(k·C*) = sum(S)`, so the total exposure is **exactly independent of the
> diffusivity** (measured: `sum(k·C_ss) = 1.0000000000` against `sum(S) = 1`, to
> **4.378e-13**, at every D on the ladder). If raising D raised exposure
> everywhere the total would have to rise. It cannot, so the curves must cross,
> and **45.59 %** of tissue does lose exposure. The repaired property
> ("monotone beyond r\*") is **true by construction and worthless**: only
> **0.12 %** of tissue lies beyond r\* on the median row, which is why it is
> reported VACUOUS rather than green. Asked in a form that can fail, the answer
> splits: in single-vessel geometry the crossing is a **single clean radius**
> with 69–94 % of the domain beyond it, and on a real section **there is no such
> radius at all** — even the outermost decile of distance-to-nearest-vessel
> still loses exposure, and refining the grid 2.7× moves the overall figure from
> **0.4225** to **0.4220**, so it is not a resolution artefact.
>
> **B2 also narrows a claim written the same day.**
> `results/mms_transport_verification.md` justifies its validity guard V1 with
> *"upwind is positivity-preserving at CFL ≤ 1, so … the solver's `max(C,0)`
> clamp provably never binds."* The M-matrix half of that is right and is
> re-derived here; the upwind half is false in general. **That study's own
> numbers are unaffected** — its velocity field's two component maxima are
> co-located, so its guard quantity and the strict one are *exactly equal*
> (measured ratio 1.0000, against 2.0000 for the counterexample), and it does
> not undershoot at any CFL up to 1. Sec. 2.4 states the narrowing precisely.
> Neither file is edited: `results/` entries are dated records.
>
> **THE THIRD is D4, minor and already contained.**
> `cells.pkpd.apply_drug_kill` defaults `rng` to an **unseeded**
> `np.random.default_rng()`. Six identical direct calls gave six different kill
> fractions spanning **6.58 %** of the value. `sim.engine.simulate_response`
> always passes a seeded generator, so no shipped run can reach it; every one of
> the 17 pre-existing call sites in `tmesim/` and `experiments/` passes `rng=`
> explicitly, and `app/` calls it at all only through the engine. The default is
> a trap for the next caller, not a live bug.
>
> **Everything else is a null, and the nulls are the point.** Mass balance,
> non-negativity under any smooth flow at any Peclet number, dose and window
> monotonicity across 12 compounds × 5 sections, determinism across process
> restarts and hash seeds, and all four scale invariances hold to between
> **0.000e+00** and **5e-13** — several of them *bit-identically*.

---

## 1. What was measured, and with what

Commit `f9d443d`. Python 3.9.6, numpy 2.0.2, scipy 1.13.1, macOS 15.7.3 arm64.
Total runtime 283.6 s. The driver was executed three times end to end; every
measured number was identical across all three except the one in sec. 4.3, which
is random by construction and is the finding.

**Sections** (all five the repo can currently load):
`xenium_319`, `xenium_io`, `gse324268` (Xenium, dx = 80 µm), `visium_10x_parent_gbm`,
`visium_gbm1` (Visium, dx = 100 µm), plus `xenium_319` re-reconstructed at
dx = 30 µm and dx = 12 µm for the resolution arms of B3 and C3b. The coarse grids
are deliberate: every invariant here is a statement about the **scheme**, and a
scheme that conserves mass at dx = 80 conserves it at dx = 12 for the same
reason. Where resolution could plausibly change the answer — B3's CFL margin and
C3b's crossover — it was measured at more than one dx and the finding is reported
per grid.

**Compounds**: all 12 in `molecular.library.GBM_DRUGS`, compiled through the same
flagship path the shipped panel uses (`calibrated_drugparams` →
`apply_barrier_to_drugparams`). One of them, **bevacizumab**, compiles to
`bbb_penetrance = 0.0`, so its delivery source is identically zero and every
field is exactly zero. Its monotonicity rows are **degenerate, not passing**, and
are recorded as degenerate rather than being counted as 12/12.

### 1.1 The instrument: a bit-identical mirror of the solver

Two of the five invariants cannot be measured from what
`simulate_drug_transport` returns:

* the **mass budget** needs the per-step loss `dt·Σ(k·c^{n+1})`, and
  `TransportResult` carries only the *trapezoidal* `auc` — which sec. 1.3 shows
  is the wrong quadrature for the budget by O(dt);
* the **non-negativity** question is about the field *before*
  `np.maximum(c_new, 0.0, out=c_new)` clamps it. The shipped result can never be
  negative. Asking whether the shipped result is negative answers a different
  question than "does this scheme overshoot".

So `march_mirror` re-runs the solver's loop using **that module's own private
operators** — `_effective_diffusivity`, `_build_diffusion_operator`,
`_advective_divergence` — with a `clamp` switch and per-step accounting added. No
discretisation is re-implemented.

`assert_mirror_is_production` checks the mirror reproduces the public solver
**bit-for-bit** (`np.array_equal` on both `final_concentration` and `auc`) on all
four code paths — static source / time-resolved profile × convection off / on —
and is called at the head of every stage that uses it. **4/4 paths bit-identical.**
If the mirror ever stops being the solver, every number it produced stops being
about the solver, and the run aborts rather than reporting.

The mirror's own guard was planted and watched (sec. 6): running it at
`dt + 1 ulp` breaks bit-identity, so the check is reached and not merely applied.

### 1.2 The budget, derived from the code rather than from the continuum equation

The shipped step is

```
(I - dt·L + dt·K) c^{n+1} = c^n + dt·P_n·S - dt·Adv(c^n)
```

`L` is assembled in finite-volume form with face coefficients set to zero at the
domain edge **and at every tissue / non-tissue face**, so its row sums are
identically zero. `Adv` is built from antisymmetric face fluxes zeroed on exactly
those same faces. Therefore `1ᵀL = 0` and `1ᵀAdv(c) = 0`, and summing the step
over every cell leaves

```
Σ c^{n+1} + dt·Σ(k·c^{n+1})  =  Σ c^n + dt·P_n·Σ S
```

and, cumulatively from `c⁰ = 0`,

```
FINAL MASS = INJECTED - REMOVED
INJECTED   = dt·Σ_n P_n·Σ S
REMOVED    = dt·Σ_n Σ(k·c^{n+1})          ← RIGHT-ENDPOINT, not trapezoid
```

`REMOVED` splits into *cleared* and *bound* in the exact ratio
`clearance : binding` (and, with `transvascular_exchange` on, a third *exchanged*
share), because `k` is the sum of those rates acting on the same cells.

**The load-bearing detail is the quadrature.** The sink is evaluated
**implicitly**, at the new time level, so the budget's loss term is a
right-endpoint sum. The `auc` the solver returns is a **trapezoidal** integral of
the same quantity. `k·Σ(auc)` is therefore a *different* rule and must miss by
O(dt). Writing the budget from the continuum equation — where any consistent
quadrature would do — produces a check that fails at 2e-3 and looks like a leak.
This is the concrete meaning of "derive it from the code".

**Tolerance, and why 1e-10.** The operator is LU-factorised once by SuperLU and
applied by triangular solves, so each step's residual is ~ε·cond(A).
`A = I - dt·L + dt·K` is a symmetric, strictly diagonally dominant M-matrix with
cond ≈ 1 + 8·dt·D_eff/dx², a few hundred at the ladder's worst point, and
N ≤ ~4 000 steps accumulate at worst ~N·ε·cond ≈ 1e-10. Anything at 1e-13 is
round-off; anything at 1e-6 is a leak.

### 1.3 A1 / A2 / A3 — measured

26 rows: a 10-row ladder over drug parameters, grid, dt, window, source profile,
ECM density, convection and an internal no-flux hole; one row with
`transvascular_exchange` on; and 15 rows of (5 real sections × 3 compounds).

| row | rel. residual, **right-endpoint** | rel. residual, **trapezoid** | clamp mass / injected |
|---|---:|---:|---:|
| inert tracer, no sink | 5.794e-14 | 5.794e-14 *(k = 0: the two rules coincide)* | 0.000e+00 |
| reactive | 1.704e-13 | **1.986e-03** | 0.000e+00 |
| fine dt (0.5 s) | 1.465e-14 | **1.429e-04** | 0.000e+00 |
| coarse dt (60 s) | 5.291e-13 | **4.137e-03** | 0.000e+00 |
| fast sink (k = 2e-2) | 1.088e-15 | **2.778e-03** | 0.000e+00 |
| slow diffusion (D = 3) | 5.263e-16 | **2.083e-03** | 0.000e+00 |
| convection on | 1.709e-13 | **1.986e-03** | 0.000e+00 |
| time-resolved source | 1.787e-13 | **1.918e-03** | 0.000e+00 |
| internal no-flux hole | 1.501e-13 | **1.974e-03** | 0.000e+00 |
| zero ECM | 1.219e-13 | **1.986e-03** | 0.000e+00 |
| transvascular exchange ON | 5.474e-14 | **2.777e-03** | 0.000e+00 |
| 15 × (real section, compound) | 1.080e-15 … 7.215e-15 | 1.647e-03 … 2.777e-03 | 0.000e+00 |

**A1 PASS** — worst **5.291e-13**, 189× inside the pre-registered 1e-10, over
every input class tried, including an interior non-tissue hole (so an internal
tissue boundary is a genuine no-flux wall, not a leak) and with the optional
back-flux term on (so the budget is not special-cased to the shipped default).

**A2 CONFIRMED** — the trapezoidal budget misses by 1.429e-04 to 4.137e-03, and the
error is first-order in dt: the test in
`test_conservation_invariants.py::test_the_trapezoidal_budget_does_not_close`
sweeps dt ∈ {2.5, 5, 10, 20} s and requires each 2× step change to move the error
by 1.7–2.3×. It does. A negative control that passed would have falsified the
derivation behind every other budget number here.

**A3 PASS** — the clamp created **exactly zero** mass on all 26 rows. The clamp
is not dead code, however: it fires, hard, on the B2 input (sec. 2.2). "Never
reached on a physical input" is the accurate statement, not "never reached".

**Known detection floor, measured rather than assumed.** A 1e-10 threshold does
not catch an arbitrarily small leak. Planting a per-step relative loss `ε` into a
copy of the march and measuring at 60 steps gives residual ≈ 26.7·ε:

| planted per-step leak ε | rel. residual | guard |
|---:|---:|:--|
| 0 | 5.4671e-14 | green (correct — this is the control) |
| 1e-13 | 2.6182e-12 | **green — missed** |
| 3e-12 | 8.0104e-11 | **green — missed** |
| 1e-11 | 2.6714e-10 | red |
| 1e-10 | 2.6719e-09 | red |
| 1e-9 | 2.6719e-08 | red |

**Measured bracket:** A1 fires at ε ≥ 1e-11 per step and is blind at ε ≤ 3e-12,
over a 60-step window. The response is linear in ε (residual ≈ 26.7·ε), so the
crossing extrapolates to ≈ 3.7e-12 per step — that figure is an interpolation
between the two measured rows, not a measurement, and is labelled as such. Either
way the floor is a property of the threshold, and it is stated here rather than
discovered later by whoever next relies on the guard. The full plant ladder is in
sec. 6.

---

## 2. NON-NEGATIVITY — a null everywhere it is physical, and a violation the guard permits

### 2.1 Why the implicit half cannot be the culprit, and B1

`A = I - dt·L + dt·K` has a strictly positive diagonal, non-positive off-diagonal
entries (−dt times `L`'s positive face coefficients) and strict diagonal
dominance, so it is an **M-matrix** and `A⁻¹ ≥ 0` elementwise. Since
`c^{n+1} = A⁻¹·rhs`, the field can only go negative if `rhs` does, and
`rhs = c^n + dt·P_n·S - dt·Adv(c^n)` has two non-negative terms. **The implicit
diffusion + reaction step cannot overshoot at any dt or any Peclet number. Only
the explicit upwind advection can.** That is a derivation, and B1 is its test.

Peclet ladder `Pe = v·dx/D_eff` over four decades, with dt chosen to hold
CFL = 0.9 at each rung (D_eff = 237 µm²/s):

| Pe | v (µm/s) | dt (s) | min raw C | relative undershoot |
|---:|---:|---:|---:|---:|
| 0 | 0 | 1 | 0.000e+00 | **0.000e+00** |
| 0.1 | 1.580 | 8.544 | 0.000e+00 | **0.000e+00** |
| 1 | 15.800 | 0.8544 | 0.000e+00 | **0.000e+00** |
| 10 | 158.00 | 0.08544 | 0.000e+00 | **0.000e+00** |
| 100 | 1 580.0 | 0.008544 | 0.000e+00 | **0.000e+00** |
| 1000 | 15 800 | 0.0008544 | 0.000e+00 | **0.000e+00** |

and the same ladder under a smoothly **diverging** radial field — the
configuration in which a conservative upwind scheme is most exposed — at
CFL = 0.9 and 0.99: ten more rungs, all **0.000e+00**.

**B1 PASS.** Note the value: it is not "small", it is **exactly zero**, sixteen
times over. The pre-registered threshold was -1e-9, so the measurement beat it by
an infinite margin rather than a numerical one. That is what the M-matrix
argument predicts, and it is why B2 had to be constructed rather than found.

### 2.2 B2 — the CFL guard admits inputs where positivity fails

The outflow fraction the explicit step removes from cell P is

```
(dt/dx)·[ max(vf_e,0) + max(-vf_w,0) + max(vf_n,0) + max(-vf_s,0) ]
```

with each face velocity the mean of its two cells. The x pair is bounded by
`(vx_E - vx_W)/2 ≤ max|vx|` and the y pair by `max|vy|`, so positivity requires

```
dt·(max|vx| + max|vy|)/dx  ≤  1          ← the two maxima taken INDEPENDENTLY
```

for stability, and — in the worst case where all four faces of one cell are
outflow faces — the stricter

```
dt·(max|vx| + max|vy|)/dx  ≤  1/2
```

The shipped guard computes `dt·max(|vx| + |vy|)/dx`: a maximum **over cells** of
the **pointwise sum**. Those two numbers are equal only when both maxima sit on
the same cell. When they do not, the guard is up to **2× smaller than the
quantity that matters**.

**The minimal reproducing input.** Nine lines, one grid cell, four non-zero
velocities:

```python
vel = np.zeros((2, ny, nx))
vel[0][i, j-1] = -v      # west neighbour pushes away in -x
vel[0][i, j+1] = +v      # east neighbour pushes away in +x
vel[1][i-1, j] = -v      # north neighbour pushes away in -y
vel[1][i+1, j] = +v      # south neighbour pushes away in +y
# every face of (i, j) is an outflow face:
#   east face velocity = (0 + v)/2 = +v/2, west = (-v + 0)/2 = -v/2, likewise y
# so the explicit step removes dt·2v/dx·c from (i, j),
# while max(|vx| + |vy|) over cells is only v.
```

driven by a one-step bolus (`source_profile = [1, 0, 0, …]`) so the source does
not keep refilling the cell it is being drained from, with `D → 0` so the
implicit smoothing does not mask the result. `divergent_cross_velocity()` in the
driver builds it.

| CFL | accepted by the shipped guard? | measured `min(C_raw)/max(C)` | derived law `-(2·CFL−1)` | clamp mass / injected |
|---:|:--:|---:|---:|---:|
| 0.40 | yes | **+0.000000** | +0.000000 | 0.0000 |
| 0.50 | yes | **+0.000000** | +0.000000 | 0.0000 |
| 0.60 | yes | **-0.200000** | -0.200000 | **0.2000** |
| 0.75 | yes | **-0.500000** | -0.500000 | **0.5000** |
| 0.90 | yes | **-0.800000** | -0.800000 | **0.8000** |
| 1.00 | yes | **-1.000000** | -1.000000 | **1.0000** |

The closed form is reproduced to **4.000e-07** over the 18 rows whose
diffusivity is negligible (`D = 1e-6`): grids 9×9, 11×13 and 21×21, spacings 15,
7.5 and 25 µm, velocities 30, 2 and 100 µm/s, each at all six CFL values. Two
further configurations at `D = 1e-3` and `D = 1e-1` were run and sit fractionally
inside the law, because the implicit step smooths part of the undershoot away;
they are counted in the total but not in the law's error, and the driver reports
both numbers separately. Across all five configurations, **20 of 30 rows are
accepted by the shipped CFL guard and undershoot past the -1e-3 threshold.**

Two things follow, and they are separate.

1. **The scheme is not positivity-preserving in the range its own guard
   accepts.** At CFL = 1 the field goes to −100 % of its peak. The task's own
   calibration — "a scheme that produces -1e-9 is arguably fine and a scheme that
   produces -0.3 is not" — puts this squarely in the second category.
2. **The zero-clamp converts that into a silent mass source.** At CFL = 1 it
   creates mass equal to **100 %** of everything injected: the domain ends the
   run with twice the drug that entered it, and nothing in the returned
   `TransportResult` records that it happened. A3's null holds only because no
   physical input reaches this regime.

The converse was measured too: on the *same* adversarial field, every CFL ≤ 1/2
gives undershoot exactly 0.000e+00. The derived bound is not merely necessary, it
is sufficient here.

**This audit does not fix it.** `tmesim/physics/diffusion.py` is not this work's
territory. The change that would close it is a one-line tightening of the guard's
threshold from `1` to `1/2` on the independent-maxima quantity — or a
positivity-preserving advection step, which is a larger job — and it is left to
whoever owns that file.

### 2.3 B3 — is it reached? Not yet, and the margin is gone

`strict_CFL = dt·(max|vx| + max|vy|)/dx`, on real Darcy velocity fields from
`physics.pressure.solve_ifp`:

| section | dx (µm) | dt = 0.5 s (the `PDEConfig()` default) | dt = 10 s (the app's `PDE_DT_S`) |
|---|---:|---:|---:|
| xenium_319 | 80 | 0.0015 | 0.0304 |
| xenium_io | 80 | 0.0005 | 0.0094 |
| gse324268 | 80 | 0.0010 | 0.0207 |
| visium_10x_parent_gbm | 100 | 0.0005 | 0.0096 |
| visium_gbm1 | 100 | 0.0005 | 0.0104 |
| xenium_319 | 30 | 0.0069 | 0.1374 |
| **xenium_319** | **12** (shipping) | 0.0325 | **0.6506** |

Thirteen of fourteen (section, dt) pairs sit inside the derived 1/2 bound with
orders of magnitude to spare. **One does not**: the app's own configuration —
`app/server/jobs.py`, `PDE_DT_S = 10.0` — on the dx = 12 µm grid the
steady-state and resolution work calls the shipping one. The interstitial
velocity grows as the grid is refined — back out of the table above,
`max|vx| + max|vy| = strict_CFL · dx / dt`, and `xenium_319` goes from
0.243 µm/s at dx = 80 to 0.781 µm/s at dx = 12 — so refining the grid **erodes
this margin**, and it is already spent.

**Nothing actually undershoots.** The raw, pre-clamp field was marched on every
section, at both timesteps, at dx = 80, 30 and 12, with a bolus profile chosen to
be the worst case, and the worst relative undershoot anywhere was
**0.000e+00**. Smooth Darcy fields from a pressure solve are simply not
adversarial in the way sec. 2.2's cross is: the local divergence never reaches
`2·max|v|/dx`.

The honest reading: **the invariant holds in production, and it holds by luck
rather than by the guard.** The guard would not have stopped a field that broke
it, and the shipped configuration is now inside the range where such a field
would be accepted.

### 2.4 A claim on the current claim surface that B2 narrows

`results/mms_transport_verification.md`, written the same day by separate work,
justifies its validity guard V1 with a general statement — twice, and the second
time with the word *provably*:

> *"With S ≥ 0, `A = I - dt L + dt K` is an M-matrix and upwind is
> positivity-preserving at CFL ≤ 1, so C ≥ 0 for every step and the solver's
> `max(C,0)` clamp provably never binds."*

**The M-matrix half is right and is re-derived independently in sec. 2.1. The
upwind half is not right as a general statement**, and sec. 2.2 is the
counterexample: at CFL = 0.6 on the divergent-cross field the clamp binds and
creates 20 % of the injected mass; at CFL = 1 it creates 100 %.

**That study's own measurements are unaffected, and this is not a correction to
any number in it.** Its E3 velocity field is `v = (V sin(πx/L), V sin(πy/L))`,
for which the guard's quantity and the strict one are *exactly equal* —
measured: `max(|vx| + |vy|) = 40.0000` and `max|vx| + max|vy| = 40.0000`, ratio
**1.0000** — because that field's two velocity maxima are co-located at the
domain centre. Marched at CFL = 0.6, 0.9 and 1.0 it undershoots by
**0.000e+00** at every step, consistent with the `min C` column of that
document, which never drops below 1.000. The divergent-cross field of sec. 2.2
has ratio **2.0000**, the worst case the stencil admits, and that is the entire
difference between the two.

So the narrowing is one clause: *for a velocity field whose two component maxima
are co-located, the guard's CFL ≤ 1 is the positivity bound; in general it is
CFL ≤ 1/2 on `max|vx| + max|vy|`, and the shipped guard does not measure that
quantity.* That document's "Limits" section already states it does not test the
CFL guard or the non-negativity clamp, so the two records are complementary
rather than in conflict — but the reason given for V1 is stated generally, and
generally it is false. Neither file is edited here; `results/` entries are dated
records, and this is the dated record that narrows it.

---

## 3. MONOTONICITY — two exact nulls and one claim that was never true

All three properties were tested the way the task specifies: for **every one of
the 12 compiled compounds on every one of the 5 sections**, not on a
representative pair. 55 (compound, section) rows; bevacizumab contributes 5
degenerate rows (zero source) recorded as degenerate.

### 3.1 C1 — more dose never less exposure, anywhere

Dose ladder λ ∈ {1e-3, 0.1, 1, 10, 1e3} — a 10⁶ span — on the delivery source.

* **Exact proportionality:** worst deviation of `auc(λS)` from `λ·auc(S)`
  anywhere on any grid: **5.458e-15**, against a 1e-12 threshold. The transport
  operator is linear in the source and the clamp never fires, so this is an
  *equality*, not an inequality.
* **Pointwise ordering:** `all(auc(λ_{i+1}·S) ≥ auc(λ_i·S))` held at **every grid
  point, on 55 of 55 rows**.

### 3.2 C2 — a longer window never less AUC, anywhere

Window ladder T ∈ {600, 1800, 3600} s at fixed dt. `auc` is a cumulative integral
of a non-negative field whose trajectory does not know where the window ends, so
extending it can only add. Worst increment over all rows and all grid points:
**0.000e+00** — i.e. never negative.

The grid-wide minimum is exactly 0 because non-tissue cells hold `auc = 0` in
both windows, which would make the inequality true for a vacuous reason. The
**tissue-restricted** minimum is reported separately and is also **0.000e+00**,
so the statement is not resting on the mask.

### 3.3 C3 — higher diffusivity: the claim is false, and the repair does not transfer

**Pre-registered as EXPECTED FALSE, from a derivation.** At the fixed point,
`1ᵀ(K - L + Adv)c* = 1ᵀS` reduces by the same zero-row-sum argument to

```
Σ (k · C*)  =  Σ S
```

so with uniform `k` the **total** steady exposure is `ΣS/k` — *exactly
independent of the diffusivity*. Measured directly in single-vessel geometry:
`Σ(k·C_ss) = 1.0000000000` against `ΣS = 1`, at every D from 25 to 400 µm²/s,
worst error **4.378e-13**. Measured through the finite-window solver on real
tissue: total AUC across the D ladder moves by **1.588e-14**.

If raising D raised exposure everywhere, the total would have to rise. It cannot.
**The curves must cross.**

They do: over 55 rows, the worst fraction of tissue that *loses* exposure when D
is raised is **45.59 %**.

**The obvious repair is vacuous, and is reported as vacuous.** Define r\* as the
outermost radius at which raising D still takes drug away; then "monotone beyond
r\*" is true by construction. It measures as 0 decreases beyond r\* — and the
region beyond r\* is **0.12 % of tissue on the median row, under 0.5 % on every
row, 55 of 55 rows below the 1 % vacuity guard**. A green check over a region
that does not exist is the exact defect this project has named
(`results/` has a standing record of seven of them), so the verdict is **VACUOUS**,
not PASS, and the question was re-asked in two forms that can fail.

**C3a — single-vessel geometry, at the steady state.** One vessel, 121×121 at
dx = 8 µm, `steady_state_concentration` (the exact fixed point of the same
stepper — E5 — so no finite-window confound):

| k (1/s) | D_lo → D_hi | r\* (µm) | r\* / L_lo | fraction of domain beyond r\* | sign crossings along a ray |
|---:|---|---:|---:|---:|---:|
| 7e-4 | 25 → 50 | 263.9 | 1.40 | 0.7674 | **1** |
| 7e-4 | 50 → 100 | 284.1 | 1.06 | 0.7278 | **1** |
| 7e-4 | 100 → 200 | 296.4 | 0.78 | 0.7051 | **1** |
| 7e-4 | 200 → 400 | 301.3 | 0.56 | 0.6948 | **1** |
| 5e-3 | 25 → 50 | 129.2 | 1.83 | 0.9434 | **1** |
| 5e-3 | 50 → 100 | 180.5 | 1.80 | 0.8904 | **1** |
| 5e-3 | 100 → 200 | 233.2 | 1.65 | 0.8161 | **1** |
| 5e-3 | 200 → 400 | 268.3 | 1.34 | 0.7579 | **1** |

**PASS.** Exactly one crossing in every case, and 69–94 % of the domain lies
beyond it. Here the repaired property is a real statement: *beyond r\*, higher
diffusivity never gives less.*

**C3b — it does not transfer to real tissue.** A *different and narrower*
measurement than the 45.59 % above: one section (`xenium_319`), one compound
(temozolomide), one D step (D → 4D), at the **steady state** rather than over a
finite window, binned by distance to the nearest high-delivery cell (top decile
of the source), at two grid spacings. The two figures are not the same
quantity and are not compared as one.

| distance band (µm) | n, dx = 80 | fraction decreasing | distance band (µm) | n, dx = 30 | fraction decreasing |
|---|---:|---:|---|---:|---:|
| 0.0 – 0.0 | 436 | 0.8968 | 0.0 – 0.0 | 3 112 | 0.7002 |
| 0.0 – 80.0 | 1 579 | 0.7258 | 0.0 – 30.0 | 11 419 | 0.6024 |
| 80.0 – 80.0 | 1 143 | 0.6605 | 30.0 – 30.0 | 8 307 | 0.5658 |
| 80.0 – 113.1 | 1 861 | 0.5895 | 30.0 – 42.4 | 13 618 | 0.5260 |
| 113.1 – 113.1 | 718 | 0.4763 | 42.4 – 42.4 | 5 311 | 0.4638 |
| 113.1 – 160.0 | 1 232 | 0.4091 | 42.4 – 60.0 | 9 261 | 0.4358 |
| 160.0 – 178.9 | 1 138 | 0.2566 | 60.0 – 67.1 | 8 495 | 0.3474 |
| 178.9 – 226.3 | 817 | 0.1934 | 67.1 – 67.1 | 4 545 | 0.3032 |
| 226.3 – 253.0 | 572 | 0.0769 | 67.1 – 94.9 | 8 335 | 0.2458 |
| **253.0 – 1088.1** | **549** | **0.0437** | **94.9 – 805.0** | **3 258** | **0.0844** |

The fraction falls monotonically with distance — the mechanism is real and
visible — **but it never reaches zero.** Even the outermost decile still loses
exposure, in 4.4 % of its cells at dx = 80 and 8.4 % at dx = 30. The tissue-wide
figure for *this* measurement is **0.4225** at dx = 80 and **0.4220** at dx = 30,
so a 2.7× refinement moves it by 0.0005: **this is not a discretisation
artefact.**

The reason is structural. With many vessels at many spacings, *distance to the
nearest vessel* does not order cells by their position relative to the
**aggregate** crossover: a cell 100 µm from its nearest vessel but ringed by
others sits in the near field, while a cell 100 µm from its only vessel sits in
the far field. A single radius cannot separate them.

**What may be said:** in single-vessel geometry, higher diffusivity never gives
less beyond a crossover radius of roughly 0.6–1.8 × the penetration length.
**What may not be said:** that this holds for a far-from-vessel cell on a real
section. It does not, it was measured at two resolutions, and
`test_the_far_field_repair_does_not_transfer_to_real_tissue` pins it so it cannot
be quietly claimed later.

### 3.4 C4 — kill fraction is monotone in dose at a fixed seed

With the seed fixed, `apply_drug_kill` draws the same per-cell uniform vector at
every rung, and the Hill map is monotone in exposure, so the killed **set** can
only grow. Measured on the engine, dose ladder over a 1000× span, **all 12
compounds**: every rung non-decreasing, smallest rung **+9.659e-04**.

*Only the sign of each rung is reported.* The invariant is the ordering **within**
one compound across its own dose ladder; it needs no absolute kill fraction, and
absolute kill fraction is not a defensible output of this project. Nothing here
compares one compound to another, and the compounds are listed alphabetically by
the driver for that reason.

---

## 4. DETERMINISM — exact, and one unseeded generator found

### 4.1 D1 / D2

Identical spec, `seed = 7`, `simulate_response` on `xenium_319` at dx = 80.
SHA-256 over the raw bytes of ten objects: the AUC field, the final concentration
field, the oxygen field, vascular density, ECM density, cell x-positions, cell
sensitivities, the viable-fraction trajectory, `repr(kill_fraction)` and
`repr(calibrate_alpha())`.

* **Same process, two calls:** bit-identical.
* **Four fresh subprocesses, `PYTHONHASHSEED` ∈ {0, 1, 12345, random}:**
  **10 of 10 objects identical** in all four. `auc = cc909893fd686078…`,
  `kill_fraction = 0.059339828657819615`, in every arm.

`PYTHONHASHSEED` is read once at interpreter start, so this had to be done in
fresh subprocesses; setting it inside the test process would prove nothing, and
the test does not do that.

### 4.2 The ordering search

`tmesim/` and `app/` were searched for iteration over an unordered `set`, for a
`set`/`dict` feeding a numeric output, and for an unseeded generator. Findings:

* Every `set(...)` that reaches an output is wrapped in `sorted(...)`.
* The three remaining bare `set(...)` uses
  (`reconstruct/region_label.py:314`, `reconstruct/marker_panels.py:318`,
  `validation/endothelial_density.py:1153-1155`) are **membership tests inside
  comprehensions that iterate an ordered list** — `[r for r in REGION_ORDER if r
  in set(labels)]`. A set cannot reorder anything when it is only on the
  right-hand side of `in`.
* `time.time()` appears in `physics/diffusion.py`, `physics/oxygen_proxy.py`,
  `sim/engine.py` and `reconstruct/vessel_network.py` — all four inside
  `__main__` self-test blocks, none on any library path.
* **One unseeded generator:** `cells/pkpd.py:157`, `apply_drug_kill`'s
  `rng = np.random.default_rng()` when the caller passes none.

### 4.3 D4 — the unseeded default, reached and measured

Six identical direct calls to `apply_drug_kill` with `rng` omitted, on the same
microenvironment, transport result and PD config:

> kill fraction ranged **0.057870 … 0.061944**, spread **4.074e-03**, which is
> **6.58 %** of the value.

The defect is real and is **reached**, not merely present in the source — and
the size of the spread is itself random, which is the cleanest possible
demonstration that it is not deterministic. Five executions of this same stage,
on identical inputs, gave spreads of **2.76 %**, **3.77 %**, **4.54 %**,
**6.58 %** and **9.31 %** of the value. Every other number in this document
reproduced bit-for-bit across those runs; this one is the only one that did not,
by construction.

It is **not reachable through the engine.** `sim.engine.simulate_response`
constructs `np.random.default_rng(cfg.seed)` and passes it, so three repeated
engine calls returned a single identical value. Enumerated mechanically: 17
pre-existing call sites of `apply_drug_kill` in `tmesim/` and `experiments/`,
four more in `tests/`, none in `app/` — and **every one of them passes `rng=`
explicitly**. The only two calls in the repo that omit it are this audit's own
(`experiments/conservation_audit.py` and `tests/test_conservation_invariants.py`),
which omit it deliberately in order to reach the defect. The default is a trap
for the next caller who does not, and the docstring already warns of it — but a
default that silently destroys reproducibility is a footgun regardless of its
docstring.

`tests/test_conservation_invariants.py::test_apply_drug_kill_default_rng_is_unseeded`
asserts the defect is **still present**, so that whoever fixes it is forced to
come back and supersede this section rather than leaving it stale.

---

## 5. SCALE INVARIANCE — re-derived here, not read off the record

The repo already claims two invariances, in `results/reach_invariances_measure.md`
(the barrier invariance and the capacity identity) and in
`results/steady_state_reach.md` (capacity is absent from the steady problem).
**Nothing below is transcribed from either.** Each was re-derived from the
discretisation and re-measured on independently constructed inputs, so that
agreement is corroboration rather than an echo.

### 5.1 E1 — exposure is exactly proportional to any scalar source factor

This is the invariance the reach map rests on: the plasma level and
`bbb_penetrance` both enter as a single multiplication on the source, so a ratio
between two points of the map cannot depend on either.

| λ | max relative deviation from an exact rescale |
|---:|---:|
| 1e-3 | 3.800e-15 |
| 0.1 | 3.730e-15 |
| 0.5 | **0.000e+00** (bit-identical) |
| 2 | **0.000e+00** (bit-identical) |
| 10 | 3.700e-15 |
| 1e3 | 3.049e-15 |

Worst **3.800e-15** over a 10⁶ span, against a 1e-12 threshold. Consistent with
the independent measurement in `results/pk_ranking_robustness.md` sec. 1.2
(1.41e-14 over a 5 200× plasma span at the shipping configuration) — different
grid, different compound, different ladder, same conclusion.

The exact-zero rows are worth noting: λ = 0.5 and λ = 2 are powers of two and
therefore exactly representable, so the rescale is bit-identical. λ = 0.1 is not,
which is precisely why `results/reach_invariances_measure.md` sec. 2.5 had to
narrow "bit-identical" for the barrier arm. This measurement reproduces that
structure from scratch.

### 5.2 E2 — superposition

`auc(S₁ + S₂) = auc(S₁) + auc(S₂)`, two disjoint vessel patterns: max relative
deviation **2.367e-15**. This is a stronger statement than E1 — a clamp that had
fired anywhere would break additivity while leaving proportionality intact — and
it is an independent confirmation of A3's null.

### 5.3 E3 — the capacity identity, tested without a capacity term

A tissue-capacity term `R` multiplies the time derivative and nothing else, so
the capacity-`R` problem **is** the ordinary problem with `(D, k, v, S)` and
`(T, dt)` rescaled:

```
auc(D, k, v, S; T, dt)  ==  s · auc(s·D, s·k, s·v, s·S; T/s, dt/s)      s = 1/R
```

No capacity term is implemented in this repo and none is needed to test the
identity one would satisfy.

| s | equivalent R | max relative deviation |
|---:|---:|---:|
| 0.25 | 4 | **0.000e+00** |
| 0.5 | 2 | **0.000e+00** |
| 2 | 0.5 | **0.000e+00** |
| 4 | 0.25 | **0.000e+00** |
| 10 | 0.1 | 3.757e-15 |

Bit-identical for every power-of-two `s`, and at 3.757e-15 for `s = 10` — again
the representability structure, not a physics difference. The identity holds.

The planted control (sec. 6) confirms the check is reached: rescaling
`(D, k, v, S)` but leaving `dt` alone moves the answer by **4.2324e-03**, nine
orders of magnitude above the tolerance.

### 5.4 E4 — the steady state is *exactly* capacity-invariant

`(K - L + Adv)c* = S` is homogeneous of degree one in a common factor on
`(D, k, v, S)`, so `c*` is exactly invariant to it. "Capacity is absent from the
steady problem" is therefore a statement that can fail, and here is the test:

| s | max relative deviation of `C_ss` |
|---:|---:|
| 1e-3 | 6.423e-14 |
| 0.1 | 2.753e-14 |
| 0.5 | **0.000e+00** |
| 2 | **0.000e+00** |
| 10 | 4.904e-14 |
| 1e3 | 4.869e-13 |

Worst **4.869e-13** over a 10⁶ span. The planted control moves it to **1.0244e+00**
when `k` is left unscaled — the guard is reached (sec. 6).

### 5.5 E5 — the march's fixed point *is* the steady solve

| window T | max relative deviation of the final concentration from `C_ss` |
|---:|---:|
| 1 h | 8.737e-02 |
| 24 h | **9.959e-15** |
| 96 h | **9.959e-15** |

The march converges to the steady solve and then stops improving at 9.959e-15 —
consistent with `app/server/reach.py`'s own note that the agreement floors at
backward-Euler error at `dt` and more simulated time does not remove it. The 1 h
row is the finite-window residual, not an error: at the shipping 6 h window the
field is nearer still.

---

## 6. The plants — every guard watched go red

A green check that cannot go red is not evidence. `stage_plants()` in the driver
re-runs every threshold above against a deliberately broken input **and against
the matching unplanted control**, and asserts which of the two it expects. A
guard that is red on everything is as useless as one that is green on everything,
so both halves are printed. Every plant is in the audit's own code or in the
input handed to the shipped code — **never in `tmesim/` or `app/`**, which are
not this work's territory. A guard that could only have been proved by editing
someone else's live file would have been reported as unproved instead; none had
to be.

**Result: 15 of 15 plants reached their guard, and every control stayed green.**

| guard | control (expected green) | plant (expected red) | measured |
|---|---|---|:--|
| A1 mass budget | ε = 0: 5.4671e-14 | per-step leak ε = 1e-11 | **RED** 2.6714e-10 |
| A1 mass budget | — | per-step leak ε = 1e-10 | **RED** 2.6719e-09 |
| A1 mass budget | — | per-step leak ε = 1e-9 | **RED** 2.6719e-08 |
| A3 zero-clamp | CFL 0.4: 0.0000 created | CFL 0.6 on the divergent cross | **RED** 0.2000 of injected created |
| A3 zero-clamp | — | CFL 1.0 | **RED** 1.0000 of injected created |
| B1 non-negativity | CFL 0.4: +0.000000 | CFL 0.6 | **RED** -0.200000 |
| B1 non-negativity | — | CFL 1.0 | **RED** -1.000000 |
| C1 dose ordering | 2× source: ordered | 0.5× source vs 1× | **RED** 900 grid points fell |
| C2 window | longer window: +1.2513e+04 | shorter window | **RED** -1.4597e+04 |
| C3 total invariance | D alone: 6.5967e-13 | move `k` with `D` | **RED** spread 8.2942e-01 |
| D determinism | two identical solves: equal | source at **+1 ulp** | **RED** hashes differ |
| E1 source scale | exact rescale: 7.9405e-16 | dose off by 1e-9 relative | **RED** 1.0000e-09 |
| E3 capacity identity | correct rescale: 0.0000e+00 | `dt` not rescaled | **RED** 4.2324e-03 |
| E4 steady invariance | correct scaling: 0.0000e+00 | `k` not scaled | **RED** 1.0244e+00 |
| mirror bit-identity | unplanted: reproduces the solver | mirror at `dt + 1 ulp` | **RED** arrays differ |

Two rows are recorded as GREEN *on purpose*, and they are the two most useful
rows in the table.

**The A1 detection floor — a guard that is reached, but only above a measured
threshold.** Planted per-step leaks of 1e-13 and 3e-12 produced residuals of
**2.6182e-12** and **8.0104e-11**, and the 1e-10 guard did **not** fire on
either. The response is linear (residual ≈ 26.7·ε over 60 steps), so the floor
sits between 3e-12 and 1e-11 per step. A1 therefore detects a systematic leak at
or above ≈ 3.7e-12 per step and is blind below it. That is the honest reach of
the 1e-10 threshold, it is printed by the driver as
`A1 detection floor at 60 steps: fires at eps >= 1e-11, blind at eps <= 3e-12`,
and it is stated in sec. 1.3 rather than left to be discovered by whoever next
relies on the guard.

**One plant was a NO-OP, and it is kept in the table rather than quietly
replaced.** The obvious determinism plant is to multiply the source by
`(1 + 1e-16)`. For values near 1.2 that rounds back to the original in IEEE-754:
the "perturbed" array is **bit-identical** to the source (`np.array_equal` →
`True`), so the guard staying green said nothing whatever about the guard. That
is a plant **applied but not reached** — precisely the failure mode the plant
exercise exists to catch, and it was caught in the exercise's own first row.
Replaced with `np.nextafter(src, np.inf)`, a genuine one-ulp change, the guard
fires immediately. Both rows are in the driver's output so the distinction cannot
be lost.

Beyond `stage_plants()`, two tests in `tests/test_conservation_invariants.py`
are **red-by-assertion**: `test_cfl_guard_is_not_sufficient_for_positivity` and
`test_apply_drug_kill_default_rng_is_unseeded` assert that the B2 and D4 defects
are **still present**. Each names the code change that would make it fail and
instructs that this document be superseded by a new dated record rather than
edited. They are not aspirational tests; they are a tripwire on a known defect,
so that fixing it cannot silently leave this file stale.

---

## 7. What this licenses, and what it does not

**Licensed.**

* The transport solve **conserves mass exactly**, to 5.291e-13 relative, over
  every drug / grid / timestep / window / source profile / convection setting
  tried and on every real section — including across internal no-flux boundaries
  and with the transvascular back-flux term on. That is a property of the
  scheme, not of a configuration.
* The exposure field is **exactly linear** in the delivery source: proportional
  to any scalar factor to 3.800e-15 over a 10⁶ span, and additive to 2.367e-15. Every
  ratio-between-two-points quantity this product reports is therefore free of the
  plasma level and of `bbb_penetrance` — re-derived here independently of
  `results/reach_invariances_measure.md`, and agreeing with it.
* The **capacity identity** and the **exact capacity-invariance of the steady
  state** are confirmed from the discretisation, to 3.757e-15 and 4.869e-13
  respectively, without implementing a capacity term.
* The pipeline is **bit-deterministic** across process restarts and independent of
  `PYTHONHASHSEED`; no object in the reconstruction or the solve depends on
  set or dict iteration order.
* **More dose never means less exposure**, anywhere, for any compound on any
  section; the same for a longer window; the same for kill fraction at a fixed
  seed.

**NOT licensed.**

* **That the advection step is positivity-preserving in the range the solver
  accepts.** It is not (sec. 2.2), by up to 100 % of the peak, and the clamp
  hides it.
* **That the shipped configuration is protected by its CFL guard.** It is not
  (sec. 2.3); it is protected by the smoothness of real Darcy fields, which the
  guard does not check and this audit did not prove is guaranteed.
* **That higher diffusivity means more exposure for a far-from-vessel cell on
  real tissue.** It does not (sec. 3.3), at either grid spacing tested. The
  statement survives only in single-vessel geometry.
* **Anything about accuracy, or about the discretisation's ORDER.** Every
  invariant here is *internal consistency*: a scheme can conserve mass perfectly,
  be exactly linear and be perfectly deterministic while solving the wrong
  equation at the wrong order. Nothing here touches whether `D_eff` or `k` are
  right — the descriptor path behind both remains unscored
  (`results/descriptor_falsification.md`) — and nothing here measures the order
  of the spatial or temporal discretisation. That second question is answered
  separately and independently by `results/mms_transport_verification.md`, whose
  method of manufactured solutions needs no tissue, no anchor and no
  calibration; the two studies are complementary, and sec. 2.4 records the one
  sentence of it this one narrows.
* **Any cross-compound statement.** The 12 compounds appear here only as 12
  independent instances of the same per-compound property. No ordering is
  computed, implied or reported, and no absolute kill fraction is published.

---

## 8. Reproducing this

```
cd ~/tme-sim
.venv/bin/python experiments/conservation_audit.py       # ~285 s, prints every number above
.venv/bin/python -m pytest tests/test_conservation_invariants.py -q
```

The driver asserts `calibrate_alpha() == 22.223477711172563` before the first
measurement and after the last, asserts its 21 pre-registered expectations as
data before any of them is measured, aborts if its solver mirror stops matching
`simulate_drug_transport` bit-for-bit, and finishes by planting a defect against
every threshold it just reported and checking each guard fires.

Sections are reconstructed from `data_cache/` and cached in-process; the five
loadable sections and the two extra grid spacings are listed in `SECTIONS` and
in `stage_nonnegativity`. The driver reports only its total runtime, so no
per-stage breakdown is quoted here.
