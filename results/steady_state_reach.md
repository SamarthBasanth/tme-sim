# The reach map as a steady-state object — and the previous disclosure was measuring the estimator

**Status: RUN.** Driver: `experiments/steady_state_reach.py`. Data:
`results/steady_state_reach.json`.

**Nothing was changed under `tmesim/`.** `calibrate_alpha()` asserted live at
`22.223477711172563`. Drugs compiled through the flagship path the shipped
real-tissue panel uses (`calibrated_drugparams` → `apply_barrier_to_drugparams`).
Nothing was tuned.

**Scope limit, stated up front because it bounds everything below.** Every arm
here — steady and finite-window alike — runs with
`PDEConfig.transvascular_exchange = False`, the shipped default. The
finite-window arm calls `simulate_drug_transport` bare, and the steady solve
takes `exchange_k=None`. So **both sides of every comparison are exchange-OFF**.
That makes the comparison internally consistent and licenses nothing about the
exchange term.

---

## 0. Verdict up front

| question | answer |
|---|---|
| numerical route | **one sparse linear solve** — the exact fixed point of the shipped stepper, no window, no `dt` |
| cost of the map's own field | **0.90 s** against **43.70 s** for a 6 h march of the same field — **48.6×** |
| cost added to a RUN | **+1.854 s**, i.e. **+4.2 %**. The run still marches; see §2 before quoting the 48.6× |
| does `P90/P10` still discriminate? | **YES.** Disjoint, and marginally *wider*: worst-case dox/tmz **1.789×** steady vs **1.771×** finite |
| capacity residual at `R = 39` | **0.3731 → 2.548e-14** (falloff); ladder **exactly 0.000e+00** at every `R` |
| was the shipping window unconverged? | **NO.** Field converged to **1.112e-05** at 6 h. `281/T` was measuring the estimator |
| counter-case hunt | **4 arms, all pass.** One silent-wrong regime found, located ~300 orders of magnitude below anything real |

**The one thing that got worse** is in §6: the reach map's `E_ref` no longer
comes from the same fit as the published penetration length, so the guarantee
that they cannot drift apart is gone. Measured divergence 4.8 % and 0.6 %.

---

## 1. What `281/T` was actually measuring

`app/server/reach.py`'s disclosure, as shipped 2026-08-14, said:

> *"THIS WINDOW IS NOT CONVERGED: the gap to the infinite-window limit is
> 0.01302 for temozolomide at the shipping 6 h and decays as 281/T with no
> plateau, so reaching 0.001 would need about 78 h."*

The arithmetic is right. **The subject is wrong.** Write

```
AUC(x,T) = C_ss(x)*T + J(x,T)
```

with `J` the integrated transient. If `C` converges at all then `J → const`, so
an AUC-based *ratio* carries an `O(1/T)` tail that **can never plateau no matter
how converged the field is**. `281/T` is a property of the estimator. The field
converges at the reaction rate — exponentially.

Flagship 10x Xenium GBM section, dx = 12 µm, repair 2, grid 424×675, convection
ON, temozolomide. `k = clearance + binding = 5.00761e-04 1/s`, so
`tau = 1/k = 1997.0 s` and the shipping 6 h window is **10.82 tau**.

| `T` | march | `‖C(T) − C_ss‖∞ / max C_ss` | AUC falloff gap | `gap × T` | mean `J` (µM·s) |
|---:|---:|---:|---:|---:|---:|
| 0.5 h | 21.2 s | 2.937e-01 | 0.202347 | 364.2 | −7530.9 |
| 1.0 h | 64.0 s | 1.070e-01 | 0.113440 | 408.4 | −10595 |
| **6.0 h** | 213.6 s | **1.112e-05** | 0.015649 | 338.0 | **−12698** |
| 12.0 h | 376.3 s | 2.293e-10 | 0.007479 | 323.1 | **−12698** |
| 24.0 h | 858.2 s | 1.992e-13 | 0.003659 | 316.1 | **−12698** |
| 48.0 h | 1397.5 s | 1.992e-13 | 0.001810 | 312.8 | **−12698** |

**Read the last column first. `mean J` is identical at 6, 12, 24 and 48 h.** The
transient integral has plateaued exactly. That is not the decomposition assumed,
it is the decomposition measured — and once `J` is constant, `gap ∝ 1/T` follows
by arithmetic, which the `gap × T` column duly shows settling to ~313.

**The field converges exponentially at `k`**, checked against `exp(−k·ΔT)`:

| interval | observed ratio | `exp(−k·ΔT)` |
|---|---:|---:|
| 0.5 → 1.0 h | 3.643e-01 | 4.060e-01 |
| 1.0 → 6.0 h | 1.039e-04 | 1.217e-04 |
| 6.0 → 12.0 h | 2.061e-05 | 2.007e-05 |
| 12 → 24 h | 8.688e-04 | 4.027e-10 |
| 24 → 48 h | 1.000e+00 | 1.622e-19 |

The last two rows are **not** a failure of the exponential law — they are the
march hitting its own floor, which §2 makes the argument out of.

**So: at the shipping window the field is converged to 1.112e-05 of its own
maximum, and the disclosure that said otherwise was describing a different
object.** This is the derivation the shipped `SUPERSESSION` clause points at.

---

## 2. Cost — and the march's accuracy ceiling

**Both timings are here, adjacent, because the obvious prior is wrong.** "Steady
state must be more expensive than a 6 h march" is what a reader will assume. It
is not.

The shipped stepper advances

```
(I − dt*L + dt*K) c_new  =  c + dt*S − dt*Adv(c)
```

so setting `c_new = c` cancels every `dt` and leaves the linear system
`(K − L + Adv) c = S`. The march pays **one factorisation + 2,160 triangular
solves**; the direct solve pays **one factorisation + one triangular solve**.

**Read §2.2 before quoting any ratio from this section.** The 48.6× is real and
it is about the FIELD. It is not a statement about how long a run takes, and the
difference matters.

Measured on an **idle machine** (best of 3; an earlier draft of this file quoted
3.37 s / 191.1 s, which were taken with two other heavy jobs on the box and were
inflated by roughly 4× on both sides — the ratio survived, the absolute numbers
did not, and they are corrected here). Same reconstruction, same compiled drug,
dx = 12, repair 2, convection ON:

| compound | **steady solve** | **6 h march of the same field** | ratio |
|---|---:|---:|---:|
| temozolomide | **0.90 s** | **43.70 s** | **48.6×** |
| doxorubicin | **0.94 s** | **45.20 s** | **48.1×** |

At every grid, with the linear system's own residual `‖Ac − S‖∞ / ‖S‖∞`:

| dx (µm) | unknowns | steady solve | rel. residual |
|---:|---:|---:|---:|
| 30 | 46,070 | 0.21 s | 1.835e-14 |
| 20 | 103,530 | 0.74 s | 3.549e-14 |
| 16 | 161,226 | 1.70 s | 7.626e-14 |
| **12** | **286,200** | **3.55 s** | **1.866e-13** |
| 8 | 642,620 | 6.10 s | 3.978e-13 |

Doubling the unknowns from dx = 12 to dx = 8 costs **1.7×**, not the
super-linear blow-up a sparse factorisation can show. There is no cost cliff in
the range this product uses.

### 2.1 The march has an accuracy ceiling the direct solve does not have

This was not planned; it fell out of a prototype left running.

The field error in §1 stops improving at **1.992e-13** and is *identical* at
24 h and 48 h. An independent convection-OFF prototype floors at **1.948e-13**
and is identical at 24 h and 72 h. That floor is the backward-Euler
time-discretisation error at `dt = 10 s`, which **does not vanish with `T`**.

So past about a day of simulated time **the march cannot get closer to its own
fixed point however long it runs.** The direct solve *is* that fixed point, to
the linear residual in the table above. This is stronger than "cheaper and
equivalent": integrating to a convergence criterion could never have reached
what one solve returns, at any cost.

**That also disposes of the "long fixed window" route.** It is dominated on
every axis — slower, and bounded away from the answer by a floor a longer window
cannot cross.

### 2.2 The steady map does NOT make a run faster. It makes it slightly slower.

This is the correction that matters most in this file, because the 48.6× above
invites exactly the wrong reading.

**The run still marches.** `simulate_response` is required for `kill_fraction`,
for `delivery_metrics`, and for the oxygen field that defines this map's own
viable band. The steady solve does not replace it — it is computed **in
addition** to it. So switching the reach map to steady state does not save a run
191 s, or 43 s, or anything. It **costs** a run:

| dx | unknowns | run's `simulate_response` | `compute_reach` finite | `compute_reach` steady | added |
|---:|---:|---:|---:|---:|---:|
| 30 | 46,070 | 4.30 s | 0.004 s | 0.127 s | **+0.123 s (+2.9 %)** |
| **12** | **286,200** | **44.10 s** | **0.014 s** | **1.867 s** | **+1.854 s (+4.2 %)** |

The steady figure includes the `solve_ifp` call `compute_reach` makes when it is
not handed a velocity, and the extra `penetration_length_fit` on the steady field
that §6 forces.

So the defensible sentences are:

* **TRUE:** obtaining this field by a direct solve is ~48× cheaper than
  obtaining it by a 6 h march, and *more accurate* (§2.1).
* **TRUE:** the steady reach map adds **4.2 %** to a run at the shipping grid.
* **FALSE:** "runs got 48× faster", or any framing in which steady state is a
  performance win for the product. It is a **correctness** win that happens to
  be cheap enough not to matter.

The reason the cost question was ever load-bearing is that it decided whether
the route was shippable at all — a route that made a run take twenty minutes
would have been ruled out regardless of its invariance properties. +4.2 % clears
that bar with room to spare. That is the whole claim.

---

## 3. The solve is the engine's, checked rather than argued

`steady_state_concentration` **imports** `_build_diffusion_operator` and
`_effective_diffusivity`, so the diffusion operator and the tortuosity law
cannot drift from the engine's. The upwind advection term is the one exception —
the engine applies it as a function of `c`, the steady solve needs it as a
matrix — and it is therefore the only place a silent disagreement could live.

Applied against `diffusion._advective_divergence` itself on random grids, masks
and velocity fields: **worst relative disagreement 2.140e-16** over 6 cases
(8 cases in the test, same tolerance).

And the fixed point is checked against the actual stepper, not just assembled:
the 24 h and 48 h rows of §1 land on the direct solve to 1.992e-13, i.e. to the
march's own floor.

Non-negativity matters here and is not cosmetic: the stepper clamps `C` at zero
every step and the direct solve clamps nothing, so if the discrete steady state
were not non-negative the two routes could not agree. Measured
`min C_ss = 0.000e+00` for both compounds on the flagship.

---

## 4. Capacity invariance, re-measured

A capacity term `R` enters as `R·dC/dt = div(D grad C) − div(v C) + S − k·C`.
**It multiplies the time derivative and nothing else, and the steady problem has
no time derivative.** `R` is absent by inspection — not small, not bounded,
absent.

`R` is emulated exactly as `experiments/capacity_ranking_robustness.py` does
(`D/R`, `binding/R`, `clearance/R`, `bbb_penetrance/R`), **convection OFF**, for
the reason that file's §2(d) gives. Flagship section, temozolomide.

| `R` | **steady falloff** | **steady ladder** | finite 6 h falloff | finite 6 h ladder |
|---:|---:|---:|---:|---:|
| 1.0 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.0000 |
| 3.2 | 3.975e-14 | 0.000e+00 | 0.0416 | 0.0461 |
| 8.0 | 0.000e+00 | 0.000e+00 | 0.1274 | 0.1221 |
| 20.0 | 2.637e-14 | 0.000e+00 | 0.2594 | 0.2497 |
| **39.0** | **2.548e-14** | **0.000e+00** | **0.3731** | 0.3569 |
| 120.0 | 4.452e-14 | 0.000e+00 | 0.5433 | 0.4749 |
| 4184.0 | 1.221e-14 | 0.000e+00 | 0.8476 | 0.9240 |

**The prediction was 0.3772 → ~0 at `R = 39`. Measured: 0.3731 → 2.548e-14.**

The finite-window arm reproduces `capacity_ranking_robustness.md`'s published
0.3772 to **1.1 %**. Exact agreement was not expected and its absence is the
point: that file ran dx = 15 / dt = 20, this runs dx = 12 / dt = 10, so the
closeness is evidence that the protocol being run is *theirs* rather than a
convenient variant.

**The ladder is exactly zero at every `R`**, not merely small. It counts viable
tumour cells below a threshold, and no cell crossed one. The falloff's ~1e-14 is
float64 round-off in forming `D/R`, `k/R` and `S/R` — worst 4.452e-14 over the
whole range.

> ### RECONCILED 2026-08-22 with `results/reach_invariances_measure.md`
>
> That file, written the same day, reports the ladder moving **0.4021** at
> `R = 39` where the finite column above reports **0.3569**. Both are right and
> the metric is not in dispute — **both compute `max |F(R) − F(1)|` over the
> seven rungs**, the identical reduction, verified in
> `steady_state_reach.py:485` against `reach_invariances_measure.py:585`.
>
> **What separates them is the ARM, and the grid is nearly irrelevant** — which
> is the opposite of the first guess, so it is worth stating with the numbers:
>
> | | this file | reach_invariances_measure.md |
> |---|---|---|
> | convection | **OFF** (§4 above, the published §2d protocol) | **ON** (the shipping configuration) |
> | grid | dx = 12 µm | dx = 15 µm |
> | value at `R = 39` | **0.3569** | **0.4021** |
>
> That file also measured its **convection-OFF** arm at dx = 15 and got
> **0.356253** — against this file's `0.3569` at dx = 12. So changing the grid
> by 25 % moves this quantity by ~0.0007, while turning convection on moves it
> by ~0.046. `0.3569` vs `0.4021` is **one arm apart, not one grid apart.**
>
> Note that this file's §4 runs convection OFF *deliberately*, because the
> published emulation is not a faithful capacity change with convection on
> (§5). That makes §4 the right comparison against the published claim and the
> *wrong* one against the shipping configuration — which is exactly why the
> other file's convection-ON number is the one the panel has to answer to.
>
> **The steady column is a different object again and is not comparable to
> either.** `0.000e+00` is exact because the steady problem has no `dC/dt` for a
> capacity term to multiply — it is not a smaller version of 0.3569, it is the
> absence of the mechanism that produces it.
>
> One consequence worth carrying: that file's per-rung breakdown shows the
> movement is **not spread across the ladder**. The bottom four rungs
> (`x ≤ 0.1`, which include the headline) are bit-identical out to `R = 80`;
> the top rung `F(1.0)` carries essentially all of it. A single max over rungs
> — which is what both files quote — therefore reports the worst rung, not the
> headline, and the two readings should never be summarised as one number
> without saying which rung moved.

---

## 5. The trap: that same sweep with convection ON

**This is the counter-case, and the shipping configuration walks into it.**

The published emulation cannot divide the velocity field. Multiply the emulated
steady equation through by `R`:

```
div(D grad C) − R*div(v C) + S − k*C = 0
```

The advection term **keeps its `R`**. So with convection on that sweep is not a
capacity change — it is a capacity change **plus an `R`-fold boost to the Péclet
number**. `capacity_ranking_robustness.md` ran convection OFF and never saw it.
`app.server.jobs.PDE_CONVECTION` is **True**.

Three arms, so the claim is identified rather than asserted:

| `R` | arm | falloff move | ladder move |
|---:|---|---:|---:|
| 39 | (a) convection OFF, published emulation | 2.548e-14 | 0.000e+00 |
| 39 | **(b) convection ON, published emulation** | **1.226** | **2.727e-01** |
| 39 | (c) convection ON, emulation **+ `v/R`** | 2.787e-14 | 0.000e+00 |
| 4184 | (a) convection OFF, published emulation | 1.221e-14 | 0.000e+00 |
| 4184 | **(b) convection ON, published emulation** | **4.380e+03** | 2.469e-01 |
| 4184 | (c) convection ON, emulation **+ `v/R`** | 1.177e-14 | 0.000e+00 |

Arm (c) is what identifies the residual as belonging to the **instrument**, not
the map: make the emulation faithful and the invariance returns to 1e-14.

**Two sentences, different subjects, both true:**

* *the map* is exactly capacity-invariant — a statement about the model;
* *the published sweep* is not a faithful capacity change when convection is
  on — a statement about an instrument.

The artefact is strongly grid-dependent (1.044e-01 on the dx = 80 test fixture
against 1.226 on the flagship), because a Péclet boost is. That is part of why
it is easy to walk into and hard to recognise, and it is why
`tests/test_reach.py::TestSteadyCapacityInvariance` keeps arm (b) as a **pinned
assertion** rather than a remark: the trap now has a test sitting in it.

---

## 6. What got WORSE: `E_ref` no longer shares a fit with the penetration length

In finite-window mode `E_ref` was **reused verbatim** from the caller's
`delivery_metrics.penetration_length_fit`, so the reach number and the published
penetration length could not drift apart *by construction*.

**That guarantee is gone.** `delivery_metrics` fits the finite-window AUC
profile; the steady map must be normalised by the first bin of the **steady**
profile. Different numbers off different fields. A caller-supplied fit is now
discarded, and the payload says so in `penetration_reused`.

The divergence is measured rather than waved at. Flagship, convection ON:

| compound | `L` from 6 h AUC | `L` from `C_ss` | ratio | analytic `sqrt(D_eff/k)` |
|---|---:|---:|---:|---:|
| temozolomide | 275.4364 µm | 288.6194 µm | 1.0479 | 291.07 µm |
| doxorubicin | 49.6505 µm | 49.9338 µm | 1.0057 | 74.47 µm |

So the loss is real but small — **4.8 %** and **0.6 %** — and it moves in the
defensible direction: the steady-fitted length is the one closer to the closed
form, 0.84 % from it for temozolomide against 5.4 % for the AUC fit. (Neither is
close for doxorubicin; the analytic form uses the bare `drug.diffusivity` and
ignores both the spatially varying `D_eff` and convection, so it is a sanity
bound, not a target.)

**This is a downgrade in a guarantee, traded for an upgrade in an invariance.**
It is worth naming as a trade rather than filing under "steady state is better".

---

## 7. Units, and the two consumer sites that hard-code them

`E_ref` is an **exposure (µM·s)** in finite-window mode and a **concentration
(µM)** in steady-state mode. The wire key stays `E_ref_uM_s` because
`app/web/js/ui/reach.js` and `app/server/report.py` read that name, and breaking
a wire contract to fix a unit label is the worse trade. The payload therefore
carries the true unit in **`E_ref_units`** and in **`units['E_ref']`**.

Two sites hard-coded the string instead of reading it. **Both were fixed on
2026-08-22** and now read `reach['E_ref_units']` with a fallback of `uM*s`:

* `app/server/report.py:2298` — the reach prose.
* `app/server/jobs.py:1007` — the `ctx.log` line.

The fallback is the right shape and is worth spelling out: a record with **no**
`E_ref_units` predates steady mode entirely, so `uM*s` is the correct reading
*for that record*. The fallback is era-correct, not a default.

### `window_s` needed no consumer fix, and the reason is a hazard

**Nothing renders the reach payload's `window_s`.** Checked rather than assumed,
because `report.py` has a field of *the same name* which is a **different
quantity** — the run's `total_time_s` (`report.py:1103`, used by
`attribute_compound` and `_compound_card`). That one still exists and is still
correct, because the run still marches (§2.2).

So there is no blank to render, and the absence is communicated through
`window_mode` in words instead. The hazard is worth recording: two fields with
one name, one of which is legitimately `None` in steady mode and one of which
must never be, is exactly the setup in which someone later "fixes" the wrong
one.

### Records written before the field existed

Runs stored before 2026-08-22 carry no `window_mode`. They **predate the field
and are silent about it**, which is not the same as saying `finite_window`.
Anything reading the key must distinguish *absent* from either value —
`WINDOW_MODES` deliberately contains no value meaning "unrecorded", so absence
cannot be round-tripped into a claim.

---

## 8. The delivery spread `P90/P10` — the stop-and-report

**The non-negotiable: if steady state destroys the discrimination, stop.** It
does not.

5 grid/mask configurations × 2 compounds, both objects measured in the same
process on the same reconstructions. The **window axis of
`delivery_spread_noise.md`'s grid is absent by construction** — there is no
window to vary — which removes a nuisance axis rather than ducking one.

| config | tmz steady | tmz finite 6 h | dox steady | dox finite 6 h |
|---|---:|---:|---:|---:|
| **dx12_repair2** *(flagship)* | **2.079** | **2.111** | **3.927** | **3.942** |
| dx12_native | 1.964 | 1.994 | 3.780 | 3.797 |
| dx16_repair2 | 2.106 | 2.134 | 3.962 | 3.977 |
| dx20_native | 1.956 | 1.985 | 3.767 | 3.778 |
| dx30_native | 1.965 | 1.994 | 3.812 | 3.826 |

| | temozolomide | doxorubicin | disjoint? | gap | worst-case dox/tmz | shipped dox/tmz |
|---|---|---|:--:|---:|---:|---:|
| **steady state** | 2.079 [1.956–2.106] | 3.927 [3.767–3.962] | **yes** | **1.662** | **1.789×** | **1.889×** |
| finite 6 h | 2.111 [1.985–2.134] | 3.942 [3.778–3.977] | yes | 1.644 | 1.771× | 1.868× |

**The separation survives and is marginally wider** — gap 1.662 against 1.644,
worst-case 1.789× against 1.771×. The steady values also sit within
`delivery_spread_noise.md`'s published nuisance spans (1.969–2.171 and
3.770–3.996), so this is not a new number, it is the same number with the window
taken out of it.

**What that does NOT license** is unchanged from `delivery_spread_noise.md` §3
and must not be quietly upgraded: robustness is not accuracy, no measurement of
these compounds' delivery heterogeneity on this tissue exists to score this
against, and the same pipeline's absolute gradient is measured ~14× too shallow
(`kyle2004_penetration.md` §6). A 1.789× margin also says nothing about two
compounds 1.1× apart.

---

## 9. The counter-case hunt

The evidence for steady state was one-sided, and one-sided evidence is when this
project has historically been most wrong. Four arms, run to look for the failure
rather than to confirm the result.

**(a) Fine grids — PASS.** §2's table. No cost cliff, residual grows from
1.8e-14 to 4.0e-13 across a 14× range of unknowns.

**(b) Clearance → 0 — PASS.** As `k → 0` the penetration length `L = sqrt(D/k)`
diverges and the operator approaches a **singular** pure-Neumann Laplacian. This
is the conditioning cliff; the question is how close the real panel gets to it.

| `k` factor | `k` (1/s) | `L = sqrt(D/k)` | rel. residual | max `C_ss` (µM) |
|---:|---:|---:|---:|---:|
| **1** *(shipping)* | 5.008e-04 | **291.1 µm** | **1.866e-13** | 18 |
| 1e-2 | 5.008e-06 | 2,910.7 µm | 1.381e-11 | 972.9 |
| 1e-4 | 5.008e-08 | 29,106.7 µm | 8.785e-10 | 9.702e+04 |
| 1e-6 | 5.008e-10 | 291,067 µm | 9.341e-08 | 9.703e+06 |
| 1e-8 | 5.008e-12 | 2,910,672 µm | 1.015e-05 | 9.703e+08 |

Degradation is **graceful over eight decades of `k`**. Even at `L = 2.9 m` —
three orders of magnitude larger than the section — the residual is 1e-05. The
shipping compound sits at `L = 291 µm` with residual 1.9e-13.

**(c) Sparse and absent vasculature — PASS, and one of these is a "no" that is
correct.**

*Absent* refuses, with no number:

```
vascular_density × 0  →  max C_ss = 0.000e+00, E_ref = None, headline = None,
null_reason = "the perivascular reference exposure E_ref is 0.0: ..."
```

With **zero endothelial cells** the refusal is earlier and harder —
`direct_vessel_distance_um` raises `ValueError("no endothelial cells
available")`, the typed error `jobs.py` already surfaces.

*Sparse* does **not** refuse, and should not. Thinning endothelial cells on the
dx = 30 section:

| endothelial cells | `E_ref` | usable bins | headline `F(0.1)` |
|---:|---:|---:|---:|
| 4,295 | 7.802 | 15 | 0.0 |
| 430 | 8.145 | 40 | 0.0 |
| 43 | 8.629 | 106 | 0.0 |
| 4 | 9.645 | 370 | 0.0 |

With 4 vessels the `MIN_BIN_COUNT = 30` rule is still satisfied per bin and
`E_ref` is a genuine perivascular mean. **A 4-vessel geometry is a real
geometry, and the map claims to describe the reach a geometry permits.**
Refusing here would have been wrong.

**(d) The one silent-wrong regime found — and deliberately not guarded.**

Scaling vascular density until `E_ref` goes denormal, against the unscaled
reference:

| `E_ref` | max &#124;Δfalloff&#124; | max &#124;Δladder&#124; |
|---:|---:|---:|
| 7.8e-08 | 2.220e-16 | 0.000e+00 |
| 7.8e-100 | 2.220e-16 | 0.000e+00 |
| 7.8e-250 | 2.220e-16 | 0.000e+00 |
| 7.8e-300 | 2.220e-16 | 0.000e+00 |
| 7.8e-310 | 3.804e-13 | 0.000e+00 |
| 7.8e-315 | 4.745e-08 | 0.000e+00 |
| **6.6e-320** | **5.050e-02** | **2.229e-02** |

At `E_ref = 6.6e-320` the map returns a plausible finite falloff whose last bin
reads 0.4866 against a true 0.5225, **with no refusal and no warning**. That is a
genuine silent-wrong regime.

**It is not guarded, on purpose.** It begins only in the IEEE denormal range —
exact to 2.2e-16 all the way down to `E_ref ~ 1e-300` — and the realised `E_ref`
on real tissue is O(1–10), because the source is `vascular_density ×
bbb_penetrance` and the plasma scalar is divided out by construction. The cliff
is ~300 orders of magnitude below anything a real section produces, and the one
axis that could in principle be tiny is exactly the one that never enters.

A positivity floor above `E_ref > 0` would be a threshold chosen to defend
against a case that cannot occur, and this repo has paid repeatedly for
thresholds chosen for reasons other than measurement. **Locating the cliff and
pinning it as a test is better than defending it**: the pin survives someone
later changing what feeds `E_ref`; a magic floor would be inherited and never
re-derived.

---

## 10. Honest read: is this better than the finite-window version?

**Yes, and the case is stronger than the three measurements that motivated it —
but it is not free, and one thing genuinely got worse.**

**Better:**

* The capacity residual is *removed*, not reduced: ladder exactly 0.000e+00,
  falloff 1e-14 round-off, against 0.3731 and 0.3569 at `R = 39` (§4).
* Obtaining the field costs **48.6×** less than marching to it, and the map
  adds only **+4.2 %** to a run (§2, §2.2) — cheap enough that the
  interactivity constraint never binds. It is NOT a performance win for the
  product and must not be sold as one.
* The march has an **accuracy ceiling** at 1.99e-13 that no window length
  crosses; the direct solve does not (§2.1). Neither alternative route —
  integrate-to-tolerance, or a long fixed window — could have reached this.
* One fewer nuisance axis: `delivery_spread_noise.md` swept exposure window as a
  real axis. It no longer exists (§8).
* No convergence criterion and no window convention to defend or to be
  quoted out of context later.
* The delivery spread, the headline, is **unharmed and marginally sharper** (§8).

**Worse, or newly risky:**

* **§6 is a real loss.** `E_ref` and the published penetration length no longer
  come from one fit. That guarantee was structural; its replacement is a
  measured 4.8 % / 0.6 % agreement, which is a weaker kind of assurance.
* **§5 is a trap the previous object did not have.** The shipped config is
  convection ON, and the standard capacity sweep produces a large residual
  against a map that is exactly invariant. That needs a disclosure and a test,
  and disclosures that need explaining are a cost.
* **A steady-state map is the wrong object for a question people will ask it.**
  It describes the field a *sustained* plasma level would produce. A compound
  cleared from plasma before the tissue equilibrates is still described by the
  steady field it would have reached. The finite-window map, for all its faults,
  answered "what did this schedule deliver" — and that is a question a clinician
  is more likely to be asking than "what does this geometry permit". **The new
  object is more defensible and less directly useful**, and if a customer
  eventually says so, they will be right and the answer will be to ship both
  rather than to argue.
* `E_ref` changing units mid-key is a wart with two live consumers (§7).

**Unchanged, and still the largest problem with this map:** none of the above is
validation. The transport law's descriptor path is unscored, the perivascular
gradient is measured ~14× too shallow in the one genuine out-of-sample test, and
an exactly capacity-invariant answer from a wrong `D` is exactly as wrong. Steady
state improved the *invariance* properties of a map whose *accuracy* is
unestablished, and it should not be quoted as having done anything else.

---

## 11. Reproducing this

```
./.venv/bin/python experiments/steady_state_reach.py --stage operator            #  ~5 s
./.venv/bin/python experiments/steady_state_reach.py --stage estimator           # ~49 min
./.venv/bin/python experiments/steady_state_reach.py --stage cost                # ~25 min
./.venv/bin/python experiments/steady_state_reach.py --stage capacity            # ~25 min
./.venv/bin/python experiments/steady_state_reach.py --stage capacity_convection # ~1 min
./.venv/bin/python experiments/steady_state_reach.py --stage spread              # ~17 min
./.venv/bin/python experiments/steady_state_reach.py --stage pathology           # ~4 min
```

Each stage merges into `results/steady_state_reach.json` rather than overwriting
it, so a single-stage re-run does not discard the sections it did not recompute;
`stages_run` records which are present.
