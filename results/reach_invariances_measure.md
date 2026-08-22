# The reach map's two invariances, measured — one is exact, one is not, and the shipped wording is right about one of them

**Status: RUN.** Driver: `experiments/reach_invariances_measure.py`. Data:
`results/reach_invariances_measure.json`.

**Nothing was changed.** No default, constant or shipped parameter.
`calibrate_alpha()` asserted live at `22.223477711172563`; `PDEConfig()`'s
defaults (`dt=0.5`, `total_time=43200.0`, `convection=True`, `save_every=200`,
`time_resolved=False`, `transvascular_exchange=False`), `REACH_LADDER` and
`HEADLINE_X` asserted at the start of every stage. Drugs compiled through the
flagship path the shipped real-tissue panel uses (`calibrated_drugparams` →
`apply_barrier_to_drugparams`). The two `bbb_penetrance` values in sec. 2 are
per-run `dataclasses.replace` copies; the library value is untouched. Two files
were added and nothing else: this document and the driver, plus the frozen
record beside them.

**Why this file exists.** Both measurements below were already being quoted —
in `app/server/reach.py`'s `SCALE_INVARIANCE_NOTE`, which ships inside every
run payload, and in `app/web/js/ui/reach.js`, which puts them on the Studio
panel — and neither had a `results/` record. By this repo's convention a number
the panel computed once and someone repeated is not a citation, so nothing else
could quote them. This is the record. **Nothing here is transcribed from those
disclosures; every figure was re-measured and is reproduced from the JSON.**

---

## 0. The answer

> **MEASUREMENT 1 — the barrier invariance — REPRODUCES.** The reach map is
> exactly invariant to `bbb_penetrance`, at all three grids tested. **One word
> in the shipped sentence has to be narrowed**: "bit-identical" is true of the
> under-reached ladder and the headline (max absolute difference **0.000e+00**,
> 0 of 7 rungs differing) and **false of the normalised falloff and the fitted
> penetration length**, which differ at 1–2 ulp (≤ **3.331e-16** absolute) and,
> once that rounding passes through the exponential fit, by **2.855e-05**
> relative at dx = 30 µm. The physics is exact; IEEE-754 is not.
>
> **MEASUREMENT 2 — the SHAPE half — REPRODUCES, to four significant figures.**
>
> | quantity | shipped | **re-measured** | verdict |
> |---|---:|---:|---|
> | normalised falloff, max absolute move at R = 39 | 0.3772 | **0.37723** (at 115 µm) | **confirmed** |
> | the range of the R = 1 falloff curve | 1.0 → 0.51 | **1.0 → 0.51284** | **confirmed** |
> | fitted penetration length, R = 1 → R = 39 | 279.7 → 76.7 µm | **279.715 → 76.740 µm** | **confirmed** |
>
> **MEASUREMENT 2 — the LADDER half — DOES NOT REPRODUCE, ON EITHER READING OF
> IT.** This is the finding, and it is reported before anything else because a
> shipped disclosure is resting on it.
>
> The shipped sentence is *"over R = 1 to 39 ... the whole ladder moves by at
> most **0.0011**, and it is exactly flat out to **R = 60**. ... across the full
> measured range of 3.2 to 4184 mL/g it moves by up to **0.80**."* The ladder is
> seven rungs. Read as written — the whole ladder — and read charitably — the
> headline rung `F(0.1)` alone — neither reading gives any of the three numbers:
>
> | | shipped | **the headline rung `F(0.1)`** | **the whole 7-rung ladder** |
> |---|---:|---:|---:|
> | max movement, R = 1 → 39 | 0.0011 | **0.000000** (exactly zero) | **0.402093** |
> | exactly flat out to R = | 60 | **80** (every R tested; first move at 100) | **1** (it moves at R = 2) |
> | max movement over R = 3.2 → 4184 | 0.80 | **0.740907** | **0.942701** |
>
> **Nor is the arm the explanation.** Re-run convection-OFF — the protocol
> `capacity_ranking_robustness.md` sec. 2d uses — the same three cells read
> 0.000000 / 0.356253, R = 80, and 0.741410 / 0.942701 (sec. 5). Neither arm
> produces any of the shipped figures.
>
> **The direction of the error depends on the reading, and both directions are
> stated rather than the worse one only.**
>
> * **Read as written — the whole ladder — all three understate the sensitivity,
>   which is the dangerous direction.** `0.0011` is not conservative-but-wrong:
>   it is **365x smaller** than the 0.402093 that actually occurs over that
>   range. "Flat to R = 60" describes a curve that is not flat past R = 1.
> * **Read as the headline rung alone, the errors are conservative and the claim
>   is safe to have relied on** — actual movement 0 against a stated 0.0011, and
>   0.741 against a stated 0.80, with the flat region understated (80, not 60).
>
> So the number was probably right about *something*, and the sentence built
> around it is wrong about the object it names. That is the defect: not a bad
> measurement, a claim that outgrew its measurement.
>
> **And the asymmetry the panel is built around is wrong as worded.** The panel
> says *"THE SHAPE IS NOT CAPACITY-INVARIANT — the ladder above nearly is, and
> this is not"*, and prints the two figures side by side on that basis. Measured:
> at R = 39 the falloff moves **0.3772** and the ladder moves **0.4021**. **They
> move by the same amount.** What is nearly immune is not the ladder — it is the
> bottom four rungs of it, `x ≤ 0.1`, which include the headline. The top rung
> `F(1.0)` moves 0.4021, and it is on the same chart.

> ### RECONCILED 2026-08-22 with `results/steady_state_reach.md`
>
> That file, written the same day, reports the finite-window ladder moving
> **0.3569** at `R = 39` where this one reports **0.4021**, and separately
> reports a **steady-state** ladder of exactly **0.000e+00**. Nothing here is in
> dispute:
>
> * **The metric is identical.** Both compute `max |F(R) − F(1)|` over the seven
>   rungs (`reach_invariances_measure.py:585` and `steady_state_reach.py:485`).
> * **The separator is the ARM, not the grid**, and the evidence is in this
>   file: the convection-**OFF** arm here measures **0.356253** at dx = 15,
>   against that file's **0.3569** at dx = 12. A 25 % change of grid moves the
>   quantity by ~0.0007; turning convection **ON** moves it by ~0.046 to
>   **0.4021**. That file's §4 runs convection OFF deliberately — it is
>   comparing against the published §2d protocol — so it is the right arm for
>   that comparison and the wrong one for the shipping configuration.
> * **The `0.000e+00` is a different object entirely** and must not be read as a
>   third value of this quantity. The steady problem has no `dC/dt` for a
>   capacity term to multiply, so the invariance there is structural rather than
>   numerical — an absence of the mechanism, not a small amount of it.
>
> The per-rung finding above is what survives all three readings and is the one
> to carry: the movement lives almost entirely in the top rung `F(1.0)`, while
> `F(0.1)` and below are bit-identical out to `R = 80`. **A single max over
> rungs reports the worst rung, not the headline**, so any summary that quotes
> one number without naming the rung is reporting something other than what the
> panel's headline does.

**What can honestly be said instead**, and it is a narrower claim that survives
measurement:

> *The under-reached fraction at thresholds `x ≤ 0.1` — which includes the
> reported headline — is EXACTLY unchanged by a tissue-capacity term at every
> capacity tested up to R = 80, and its first movement is at R = 100. The upper rungs are
> not: `F(0.5)` first moves at R = 13 and `F(1.0)` at R = 2, and by R = 39 the
> top rung has moved 0.40, which is as large as the falloff's move. The
> normalised falloff shape and the fitted penetration length are severely
> capacity-sensitive at every R (L collapses 279.7 → 76.7 µm by R = 39, a 3.6x
> change, well inside the measured 3.2–4184 range).*

**What is NOT claimed here.** Nothing in this file says a capacity term belongs
in the engine — `results/intracellular_retardation.md` refuted the only
available per-compound rule on measured Pruijn HT-29 rank agreement, and
`results/capacity_ranking_robustness.md` sec. 6.1 says so at length. This is a
sensitivity measurement on a shipped readout, one compound, one section.

---

## 1. The configuration these numbers describe, and the version they describe it at

**Read this section before quoting any number below.** The reach map is being
converted from a finite-window object to a steady-state one by a sibling agent
while this was measured. **Every number in this file is a FINITE-WINDOW number**
and describes the configuration in this table and no other.

| | |
|---|---|
| section | flagship 10x Xenium GBM, `reconstruct_real(source='xenium')` |
| cells | 40,887 measured; **23,812 tumour**; **2,007 viable** at dx = 15 µm |
| grid | dx = **15 µm**, (339, 541) |
| compound | **temozolomide**, `calibrated_drugparams` → `apply_barrier_to_drugparams` |
| timestep | dt = **10 s** (`app/server/jobs.PDE_DT_S`) |
| **window** | **21,600 s = 6 h** (`app/server/api.STUDIO_DEFAULT_TOTAL_TIME_S`) — **finite, and not converged** |
| convection | **ON** (`app/server/jobs.PDE_CONVECTION`), except where an arm says otherwise |
| steps | 2,160 in every arm, including every truncated one |
| PK | `PKConfig()` defaults: plasma_peak 50 µM, half_life 14,400 s, n_doses 1 |
| oxygen | `OxygenConfig().hypoxia_threshold` = 0.2; viable band `o >= 0.4` |
| seed | 0 |

**dx = 15 µm is the `results/ranking_robustness.md` Monte-Carlo grid, not the
studio's dx = 12 µm.** It is used here because it is the grid measurement 2 was
quoted at. Measurement 1 is additionally run at dx = 12 (the studio default) and
dx = 30.

### 1.1 The reach module is pinned to a commit, and that is not a formality

`app/server/reach.py` is owned by a sibling agent converting it to steady state.
Mid-run, its working-tree copy raised
`NameError: name 'FINITE_WINDOW_NOTE' is not defined` — a partially-applied
edit. Measuring the *current* disclosures against a file being rewritten
underneath the measurement would produce numbers that describe nothing in
particular, so the module is loaded from the **git blob** and verified by hash:

| | |
|---|---|
| pinned | `23e51f6:app/server/reach.py` |
| SHA-256 | `aaf5d5893b697780fd88af5de555dd11a78b4a241f8aeca07cf03fca7bbe7a66` |
| working tree at run time | `31002d866969eb1830ea3ae81cafb4ed397411b5de4553d3736ec035299f448c` — **differs** |

Both hashes are in the JSON. **If the steady-state rewrite lands, every number
in this file describes the PREVIOUS configuration**, and the pin is how a reader
establishes that after the fact rather than guessing. This repo has shipped
files that silently came to describe a configuration that no longer existed; the
pin exists so that this one cannot.

---

## 2. MEASUREMENT 1 — the reach map is exactly invariant to the barrier term

### 2.1 The claim, as shipped

> *"The reach ratio is EXACTLY invariant to transvascular penetrance — verified
> bit-identical at bbb_penetrance 0.1 and 1.0 — because penetrance multiplies
> every point of the map equally and therefore divides out of a ratio between
> points. ... So the delivery evidence here depends on the compound ONLY through
> its penetration length, L = sqrt(D_eff/k)."* — `app/web/js/ui/reach.js`

### 2.2 The premise is verified in the code, not assumed

`bbb_penetrance` reaching the transport solve through exactly one path is what
makes the claim true, so it was checked rather than trusted. In the simulation
path `drug.bbb_penetrance` is read in exactly one place —
`treatment.drugs.delivery_source`, `S = vascular_density * bbb_penetrance *
plasma_conc` — with no clipping on the default path (`bbb_disruption=None`,
`efflux=False`). It is a pure scalar multiplier on the source field. Every other
occurrence in `tmesim/` is a print, a serialisation, or the compile-time
`molecular.permeability.bbb_penetrance` that produces the value.

### 2.3 Method

One reach map per arm, through the same calls `app/server/jobs.py` makes:
`simulate_response` → `delivery_metrics.penetration_length_fit` →
`compute_reach`, with the distance transform computed once and shared. Two arms
per grid, identical in every respect except `bbb_penetrance` ∈ {0.1, 1.0}.

### 2.4 Result

| dx (µm) | n_viable / n_tumor | ladder bit-identical? | headline bit-identical? | falloff bit-identical? | max abs falloff diff | fitted L, rel. diff | `E_ref(0.1)/E_ref(1.0)` |
|---:|---:|:--:|:--:|:--:|---:|---:|---:|
| **15** | 2,007 / 23,812 | **YES** (0 of 7 rungs differ) | **YES** | **no** (10 of 15 bins differ) | **3.331e-16** | 3.778e-10 | 0.099999999999999992 |
| **12** | 2,008 / 23,812 | **YES** (0 of 7) | **YES** | **no** (10 of 15) | **3.331e-16** | 6.680e-10 | 0.10000000000000002 |
| **30** | 2,019 / 23,812 | **YES** (0 of 7) | **YES** | **no** (9 of 15) | **2.220e-16** | **2.855e-05** | 0.10000000000000001 |

`E_ref` itself scales by exactly the factor it should: relative error against
0.1 is **1.388e-16** at dx = 15 and dx = 12, and **0.000e+00** at dx = 30. At
dx = 15, `E_ref` is 3,774,275.051 µM·s at bbb = 0.1 and 37,742,750.51 µM·s at
bbb = 1.0. The headline `F(0.1)` is 0.0004982561036 in both arms at dx = 15,
0.0004980079681 in both at dx = 12, and 0.0 in both at dx = 30.

### 2.5 What has to be narrowed, and why it is worth narrowing

**The invariance is real and exact in exact arithmetic.** A ratio of two numbers
that are both linear in the same scalar does not depend on that scalar. The
measurement confirms it: the *counts* that produce `F(x)` are integers, and
scaling the source by 10 does not move a single cell across a single rung at any
of the three grids — hence a genuine 0.000e+00.

**But "bit-identical" is a claim about IEEE-754, and it is false for two of the
three published quantities.** The solve is a sequence of rounded floating-point
operations; multiplying the source by 0.1 (which is not exactly representable)
and then re-normalising does not return the same bits. The falloff differs by
1–2 ulp, and the exponential-plus-floor fit amplifies that into **2.855e-05**
relative on the fitted `L` at dx = 30 — small, but six orders of magnitude
larger than "identical", and large enough that anyone diffing two runs at that
grid would see it.

The correct shipped wording is the one that survives: *the ladder and the
headline are bit-identical; the falloff and the fitted length agree to
floating-point rounding.* That is not weaker as evidence — it is the same
evidence, stated at the strength it has.

### 2.6 What this licenses, and what it does not

**Licensed.**

* The reach ratio, the under-reached ladder and the reported headline carry
  **no information whatever** about transvascular penetrance. A compound with a
  10x weaker barrier produces the identical curve.
* Therefore the map depends on the compound **only** through the transport law
  — `L = sqrt(D_eff/k)` — and the same argument that makes it dose-free makes
  it barrier-blind. The two cannot be separated: they are one property.
* This is as much a **limit** as a feature, and the panel is right to say so.

**NOT licensed.**

* **Anything about brain penetration, barrier integrity, or how much compound
  crosses the vessel wall.** This map cannot see those quantities by
  construction. That is a stronger statement than "it was not measured".
* That the transport law itself is right. `bbb_penetrance` divides out; `D_eff`
  and `k` do not, and their descriptor path is unscored
  (`results/descriptor_falsification.md`).
* That the *absolute* exposure is invariant. `E_ref` moved by exactly 10x, as
  it should. Only ratios are protected.

---

## 3. The capacity identity, and the preconditions checked before relying on it

### 3.1 The identity, stated discretely rather than continuously

A tissue-capacity term `R` makes the transport equation
`R dC/dt = div(D grad C) - div(v C) + S - k C`, i.e. a pure rescale of time, so
`C_R(x, t) = C_1(x, t/R)` and `AUC_R(x, T) = R * AUC_1(x, T/R)`.

The **discrete** statement is stronger and is the one used. The solver's step is
`c_new = A^-1 [c + h*S - h*adv(c)]` with `A = I + h*(L + k)`. Discretising the
capacity equation at timestep `dt` gives `A = I + (dt/R)(L + k)` and
`rhs = c + (dt/R)S - (dt/R)adv(c)` — **the R = 1 solver at timestep
`h = dt/R`, iterate for iterate.** So

> the capacity-R reach map at window `T`, timestep `dt`
> **is** the R = 1 reach map at window `T/R`, timestep `dt/R`,

with the AUC field differing by the global factor `R`, which every reach
quantity divides out of exactly (sec. 2 is the measurement of that). The step
count is **2,160 in both**, so no time resolution is traded for the truncation.
No capacity term is implemented and no parameter is divided.

### 3.2 This route is not an approximation of the published emulation — it is the other way round

`results/capacity_ranking_robustness.md` emulates capacity by dividing
`diffusivity`, `binding`, `clearance` and `bbb_penetrance` by `R`. Its sec. 2d
states the limitation: the velocity field is built inside
`sim.engine.simulate_response` and **cannot be divided by R from outside**, so
that report runs convection **OFF**. Multiply the emulated steady equation
through by `R` and the advection term keeps its `R`: the emulation with
convection on is a capacity change *plus* an R-fold boost to the Péclet number.

The truncated-window route does not have that failure mode available to it,
because there is no separate velocity to forget: `h = dt/R` multiplies the
advection term too, giving `dt*div((v/R) C)` — exactly the retarded advection.
**Sec. 3.4 measures both claims rather than arguing them.**

### 3.3 Preconditions, measured

**(a) The non-negativity clamp never fires.** `physics/diffusion.py` applies
`np.maximum(c_new, 0.0, out=c_new)` every step. **A plasma-rescale test cannot
establish this** — `max(a*x, 0) = a*max(x, 0)` for `a > 0`, so the clamp is
positively homogeneous and an exact-rescale result survives it completely
intact. It was instrumented directly, with a counting shim on the diffusion
module's `np` that inspects the array *before* the clamp:

| arm | dt (s) | window (s) | clamp calls | steps with any negative | negative entries | worst undershoot |
|---|---:|---:|---:|---:|---:|---:|
| shipping, R = 1 | 10 | 21,600 | 2,160 | **0** | **0** | **0.000e+00** |
| most extreme truncation, R = 4184 | 0.00239006 | 5.16252 | 2,160 | **0** | **0** | **0.000e+00** |

So the one non-linearity in the solve is inert at both ends of the sweep, and
the field is linear in the source over the whole range this study uses.

**(b) The source is a single scalar.** `PDEConfig().time_resolved` is `False`,
asserted; the plasma level enters as one `plasma_mean` and `source_profile` is
`None`. A time-resolved source would break the identity, because `S(t)` would
need rescaling too.

**(c) `plasma_mean` is window-dependent, and it is harmless — measured, not
argued.** `simulate_response` averages the decaying plasma profile over
`arange(0, total_time, dt)`, so a truncated window genuinely changes the source
by a global factor (visible in the `E_ref` column of sec. 4, which falls from
4.73e+06 to 22.8 µM·s across the sweep). That factor is precisely what
measurement 1 shows every reach quantity divides out of, to 1.4e-16. **The
barrier invariance is not a side result here — it is the licence for the whole
method.**

**(d) The viable band does not move.** Oxygen is solved independently of
`PDEConfig`, so the population `F(x)` is a fraction *of* must be identical at
every R. Asserted at every row: **n_viable = 2,007 and n_tumor = 23,812 at all
19 capacities.**

**(e) The usable distance bins do not move.** They are set by per-bin
grid-point counts — pure geometry. Asserted identical at every R: 15 bins at
5…145 µm with counts 5,680 / 13,008 / 17,341 / 18,457 / 17,766 / 15,350 /
12,034 / 8,675 / 5,881 / 3,965 / 2,448 / 1,381 / 819 / 464 / 218.

### 3.4 The two routes, cross-checked

Truncated-window route vs the published parameter-division emulation, same
section, same R:

| R | convection | falloff, max abs diff | ladder, max abs diff | fitted L: truncated vs emulation |
|---:|:--:|---:|---:|---|
| 4 | **OFF** | **3.331e-16** | **0.000e+00** | 223.572093 vs 223.572094 µm (1.7e-09) |
| 39 | **OFF** | **3.331e-16** | **0.000e+00** | 67.320583 vs 67.320591 µm (1.1e-07) |
| 4 | **ON** | 7.326e-02 | 1.584e-01 | 235.100 vs **273.249** µm |
| 39 | **ON** | **7.211e-01** | **6.672e-01** | 76.740 vs **211.794** µm |

**Both halves land.** With convection off — the regime where the emulation is
exact — the two routes agree to solver precision, which validates this
implementation of the truncated route against an independently derived one.
With convection on they diverge enormously, and the divergence is the Péclet
artefact `capacity_ranking_robustness.md` sec. 2d warns about, now measured on
real tissue rather than reasoned about: at R = 39 the emulation reports a
penetration length of **211.8 µm** where the correct capacity-R answer is
**76.7 µm**.

**This also identifies the protocol the shipped numbers were measured under, and
it was the right one.** Only the convection-ON truncated route gives 76.7 µm
(convection-off gives 67.3; the emulation gives 211.8). Whoever made the
original measurement had convection on *and* the advection term correctly
retarded.

---

## 4. MEASUREMENT 2 — the capacity sweep

Shipping arm: **convection ON, truncated-window route**, dx = 15 µm, 6 h,
temozolomide. 19 capacities, 19 solves, 2,160 steps each. R = 1 **is** the
shipping map.

### 4.1 The sweep

| R | `F(0.1)` | ladder max move vs R=1 (all 7 rungs) | falloff max move vs R=1 | at (µm) | fitted L (µm) | `E_ref` (µM·s) |
|---:|---:|---:|---:|---:|---:|---:|
| **1** | 0.000498 | — | — | — | **279.715** | 4.732e+06 |
| 2 | 0.000498 | 0.019432 | 0.016054 | 125 | 264.669 | 2.696e+06 |
| 3.2 | 0.000498 | 0.044345 | 0.037404 | 125 | 246.441 | 1.639e+06 |
| 5 | 0.000498 | 0.101146 | 0.069554 | 125 | 222.194 | 9.389e+05 |
| 8 | 0.000498 | 0.156452 | 0.118137 | 125 | 191.326 | 4.834e+05 |
| 13 | 0.000498 | 0.244145 | 0.184248 | 125 | 157.472 | 2.289e+05 |
| 20 | 0.000498 | 0.309417 | 0.254758 | 115 | 128.746 | 1.142e+05 |
| 30 | 0.000498 | 0.366716 | 0.328359 | 115 | 92.589 | 5.850e+04 |
| **39** | **0.000498** | **0.402093** | **0.377235** | **115** | **76.740** | 3.788e+04 |
| 50 | 0.000498 | 0.433483 | 0.422789 | 115 | 65.172 | 2.512e+04 |
| 60 | 0.000498 | 0.459392 | 0.454949 | 115 | 58.219 | 1.860e+04 |
| 70 | 0.000498 | 0.473842 | 0.480862 | 115 | 53.129 | 1.445e+04 |
| **80** | **0.000498** | 0.483308 | 0.502432 | 105 | 49.205 | 1.162e+04 |
| 100 | **0.001495** | 0.497260 | 0.539130 | 95 | 43.481 | 8.091e+03 |
| 120 | 0.003986 | 0.506726 | 0.568269 | 95 | 39.449 | 6.033e+03 |
| 270 | 0.081714 | 0.731938 | 0.680451 | 75 | 26.406 | 1.674e+03 |
| 620 | 0.280518 | 0.848032 | 0.766973 | 65 | 18.202 | 4.626e+02 |
| 1420 | 0.538117 | 0.905830 | 0.824116 | 45 | 12.938 | 1.285e+02 |
| 4184 | 0.741405 | 0.942701 | 0.873002 | 35 | 8.590 | 2.278e+01 |

*Measured rat brain-slice `Vu,brain`: 3.2 – 4184 mL/g, n = 33, best available
constant 39 (`results/brain_slice_anchor.md`). The engine assumes 1.0.*

### 4.2 The falloff curve, R = 1 against R = 39

| distance (µm) | 5 | 15 | 25 | 35 | 45 | 55 | 65 | 75 | 85 | 95 | 105 | **115** | 125 | 135 | 145 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **R = 1** | 1.0000 | 0.9498 | 0.9227 | 0.8989 | 0.8718 | 0.8450 | 0.8147 | 0.7797 | 0.7461 | 0.7136 | 0.6789 | 0.6501 | 0.6117 | 0.5722 | 0.5128 |
| **R = 39** | 1.0000 | 0.8282 | 0.7331 | 0.6643 | 0.6023 | 0.5491 | 0.4957 | 0.4414 | 0.3932 | 0.3495 | 0.3064 | 0.2729 | 0.2392 | 0.2108 | 0.1719 |
| **move** | 0.0000 | −0.1215 | −0.1896 | −0.2346 | −0.2695 | −0.2959 | −0.3190 | −0.3383 | −0.3529 | −0.3641 | −0.3724 | **−0.3772** | −0.3725 | −0.3614 | −0.3409 |

The first bin is 1.0000 in both by construction — it is the reference. The R = 1
curve spans **1.0 → 0.51284**, and the largest move at R = 39 is **0.37723** at
115 µm. **Both shipped figures reproduce exactly.** The fitted length falls
**279.715 → 76.740 µm** (the near-field exponential fit, `near_exp_L_um`, falls
271.060 → 91.270 µm).

### 4.3 The ladder, rung by rung — this is where the disclosure fails

Movement of `F(x)` against R = 1, per rung. **The shipped claim is that the
whole of this table is ≤ 0.0011 down to the R = 39 row.**

| R | x=0.01 | x=0.02 | x=0.05 | **x=0.1 (headline)** | x=0.2 | x=0.5 | x=1.0 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0 | 0 | 0 | **0** | 0 | 0 | +0.01943 |
| 3.2 | 0 | 0 | 0 | **0** | 0 | 0 | +0.04434 |
| 5 | 0 | 0 | 0 | **0** | 0 | 0 | +0.10115 |
| 8 | 0 | 0 | 0 | **0** | 0 | 0 | +0.15645 |
| 13 | 0 | 0 | 0 | **0** | 0 | +0.00050 | +0.24415 |
| 20 | 0 | 0 | 0 | **0** | 0 | +0.02292 | +0.30942 |
| 30 | 0 | 0 | 0 | **0** | 0 | +0.06926 | +0.36672 |
| **39** | 0 | 0 | 0 | **0** | 0 | **+0.11510** | **+0.40209** |
| 50 | 0 | 0 | 0 | **0** | +0.00050 | +0.17588 | +0.43348 |
| 60 | 0 | 0 | 0 | **0** | +0.00199 | +0.23318 | +0.45939 |
| 70 | 0 | 0 | 0 | **0** | +0.00448 | +0.27354 | +0.47384 |
| **80** | 0 | 0 | 0 | **0** | +0.01445 | +0.32586 | +0.48331 |
| 100 | 0 | 0 | 0 | **+0.00100** | +0.03438 | +0.41903 | +0.49726 |
| 120 | 0 | 0 | 0 | +0.00349 | +0.05531 | +0.48779 | +0.50673 |
| 270 | +0.00050 | +0.00100 | +0.02143 | +0.08122 | +0.23667 | +0.73194 | +0.53612 |
| 620 | +0.02541 | +0.05331 | +0.14150 | +0.28002 | +0.51570 | +0.84803 | +0.55755 |
| 1420 | +0.14101 | +0.21126 | +0.37618 | +0.53762 | +0.70354 | +0.90583 | +0.56552 |
| 4184 | +0.42402 | +0.51968 | +0.64076 | +0.74091 | +0.85800 | +0.94270 | +0.57648 |

Reading it:

* **The headline rung is better than advertised.** `F(0.1)` is *bit-identical*
  at every sampled capacity from R = 1 to R = 80 — not "at most 0.0011", but
  exactly zero — and first moves at R = 100 (+0.00100). "Flat to R = 60" is
  true and understates the flat region. (Grid-limited: the boundary lies
  somewhere in (80, 100], which was not bisected — the claim is "flat at every
  R tested up to 80", not "flat at 99".)
* **The rung directly above it is not.** `F(0.5)` first moves at R = 13 and has
  moved **0.115** by R = 39.
* **The top rung is not remotely.** `F(1.0)` moves at **R = 2** — the second
  point on the grid — and has moved **0.402** by R = 39. It starts at 0.4056 at
  R = 1, so it nearly doubles.
* **`0.0011` matches nothing measured.** No rung, over R = 1…39, moves by that
  amount: the values available are 0 (five rungs), 0.115, and 0.402. The origin
  of the figure could not be reconstructed and is not guessed at here.
* **`0.80` over the measured range matches nothing either**: the headline rung
  moves 0.74091 and the whole ladder 0.94270 over R = 3.2…4184.

### 4.4 The consequence for the panel

The panel prints the ladder and the falloff adjacent, with a caveat block whose
whole purpose is to stop a reader treating them as equally solid: *"THE SHAPE IS
NOT CAPACITY-INVARIANT — the ladder above nearly is, and this is not."*

At R = 39 the falloff moves 0.3772 and the ladder moves 0.4021. **The asymmetry
is not there.** The true asymmetry is *within* the ladder — between the bottom
four rungs and the top three — and it is invisible in the current wording, which
is worse than no caveat, because it directs attention away from the rung that
actually moves while the chart plots that rung.

---

## 5. The convection-OFF arm

Run because the shipped ladder claim might have been made under
`capacity_ranking_robustness.md` sec. 2d's convection-off protocol, in which
case the shipping arm would not be the like-for-like comparison. Same section,
same route, same 2,160 steps, `convection=False`, on a 9-point **subset** of the
shipping grid — sharing R = 1 and R = 39 so the arms compare rather than merely
coexist.

### 5.1 The arm

| R | `F(0.1)` | ladder max move (all 7) | falloff max move | fitted L (µm) |
|---:|---:|---:|---:|---:|
| **1** | 0.000498 | — | — | **264.222** |
| 3.2 | 0.000498 | 0.052815 | 0.036346 | 233.959 |
| 20 | 0.000498 | 0.270553 | 0.247370 | 101.502 |
| **39** | **0.000498** | **0.356253** | **0.366271** | **67.321** |
| 60 | 0.000498 | 0.405580 | 0.441768 | 53.096 |
| **80** | **0.000498** | 0.426009 | 0.487822 | 45.702 |
| 120 | **0.003986** | 0.511211 | 0.551684 | 37.340 |
| 620 | 0.285501 | 0.850523 | 0.744462 | 17.857 |
| 4184 | 0.741903 | 0.942701 | 0.848185 | 8.534 |

Per-rung movement against R = 1 on this arm:

| R | x=0.01 | x=0.02 | x=0.05 | **x=0.1** | x=0.2 | x=0.5 | x=1.0 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3.2 | 0 | 0 | 0 | **0** | 0 | 0 | +0.05282 |
| 20 | 0 | 0 | 0 | **0** | 0 | +0.02940 | +0.27055 |
| **39** | 0 | 0 | 0 | **0** | 0 | **+0.13204** | **+0.35625** |
| 60 | 0 | 0 | 0 | **0** | +0.00299 | +0.24913 | +0.40558 |
| 80 | 0 | 0 | 0 | **0** | +0.01644 | +0.35227 | +0.42601 |
| 120 | 0 | 0 | 0 | +0.00349 | +0.05929 | +0.51121 | +0.44544 |
| 620 | +0.02541 | +0.05630 | +0.14250 | +0.28500 | +0.52466 | +0.85052 | +0.49377 |
| 4184 | +0.42601 | +0.52317 | +0.64175 | +0.74141 | +0.85899 | +0.94270 | +0.51271 |

### 5.2 It does not rescue the claim, and the two arms agree about why

**The same result, on both arms.** `F(0.1)` is bit-identical from R = 1 to
R = 80 and first moves at R = 120 on this arm (vs R = 100 on the shipping arm);
the whole ladder has moved **0.356253** by R = 39 (vs 0.402093); over
R = 3.2…4184 the whole ladder moves **0.942701** — the same figure as the
shipping arm to six decimals, because at that capacity both arms have collapsed
the field into the same near-vessel shell.

**So the shipped `0.0011` / flat-to-`R = 60` / `0.80` triple is not reproduced
under either arm**, and the choice of arm is not the explanation for the
discrepancy. The convection-off arm is also, independently, not the arm the
shipped numbers came from: it gives `L(R = 39) = 67.321 µm` against the shipped
76.7, and a R = 1 falloff spanning 1.0 → 0.49834 against the shipped 1.0 → 0.51.
The shipping arm matches both (76.740; 1.0 → 0.51284). **The shipped SHAPE
numbers were measured convection-ON, with the advection term correctly retarded
— and the ladder numbers, whatever they were measured on, do not match either
arm.**

**One genuine difference between the arms, stated because it cuts against
reading them as interchangeable:** at R = 1 the convection-off map is a
different map — `L` is 264.222 µm rather than 279.715, and `F(1.0)` is 0.4699
rather than 0.4056. Convection is not a small term on this section; it is the
capacity *emulation* that cannot represent it (sec. 3.4), not the physics that
lacks it.

---

## 6. What may and may not be claimed from this document

**Claimable.**

* The reach ratio, the under-reached ladder and the headline are **exactly
  invariant to `bbb_penetrance`**, measured at dx = 12, 15 and 30 µm; the ladder
  and the headline are bit-identical and the falloff and fitted `L` agree to
  floating-point rounding (sec. 2.4).
* The map therefore says **nothing** about transvascular penetrance or barrier
  integrity, and depends on the compound only through the transport law.
* The **normalised falloff moves 0.37723 absolute at R = 39** on a curve
  spanning 1.0 → 0.51284, and the **fitted penetration length collapses
  279.715 → 76.740 µm** — a 3.6x change at a capacity well inside the measured
  range. The shape must not be presented as capacity-invariant.
* The under-reached fraction at **x ≤ 0.1**, including the reported headline, is
  **exactly unchanged at every capacity tested from R = 1 to R = 80** and first
  moves at R = 100 (shipping arm) / R = 120 (convection-off arm).
* The published parameter-division capacity emulation is **invalid with
  convection on**, by a measured margin: 76.740 vs 211.794 µm at R = 39
  (sec. 3.4). Convection-off, it is exact to 3.331e-16.
* The solver's non-negativity clamp does not fire at any point of this study.

**NOT claimable.**

* **That the ladder is nearly capacity-invariant.** It is not: 0.402093 over
  R = 1…39. Only its bottom four rungs are. The shipped sentence must not be
  quoted.
* Any of the three shipped ladder figures (0.0011, flat-to-60, 0.80). None
  reproduced on either reading or either arm.
* **Anything about `E_p90/E_p10`, the current headline.** It was not measured
  here. It is a spread statistic of the same field, and while the capacity
  factor `R` cancels from it as a ratio exactly as it cancels from the reach
  ratio, cancelling the factor is not the issue — capacity rescales *time*, so
  the finite-window *shape* changes, and a spread statistic reads the shape.
  **The nearest measured relatives of that statistic (the falloff, the fitted
  length) are severely capacity-sensitive.** So the presumption should be that
  the headline is too, until someone sweeps it — which is one sweep, and cheap.
  Nothing here licenses adding "and it is capacity-robust" to its list of
  virtues.
* **The rest of the paragraph these two measurements live in.** The same
  `SCALE_INVARIANCE_NOTE` also states that the gap to the infinite-window limit
  is 0.01302 for temozolomide at 6 h, that it decays as 281/T with no plateau,
  and that about 78 h would be needed to reach 0.001. **None of those three
  were tested here**, and silence about them in this file is silence, not a
  pass. Two measurements in a paragraph were checked; the paragraph was not.
* That a capacity term belongs in the engine, or that any R here is the right
  one. `results/intracellular_retardation.md` refuted the only available
  per-compound rule.
* That these numbers describe the reach map that ships after the steady-state
  rewrite. They describe `23e51f6:app/server/reach.py` at a **finite 6 h
  window** (sec. 1.1), and a steady-state map is a different object — sec. 4's
  entire effect is a finite-window effect.
* Anything about a second compound, a second section, or a second tissue. One
  compound, one section.

---

## 7. Two things found on the way that contradict shipped copy

Neither is one of the two measurements; both were met while making them, and
both are recorded because a contradiction found and not written down is a
contradiction that gets rediscovered.

**7.1 The viable-band grid sensitivity was overstated by ~10.5x.**
`app/web/js/ui/reach.js` states the viable band *"is 89% of tumour cells at
dx = 30 µm and 2,008 of 23,812 at dx = 12 µm"* and infers *"~10x grid
sensitivity"* in the denominator of anything phrased as a fraction of viable
tumour. Measured here through the same reconstruction path:

| dx (µm) | n_viable | n_tumor | fraction |
|---:|---:|---:|---:|
| 12 | **2,008** | 23,812 | 8.43% |
| 15 | 2,007 | 23,812 | 8.43% |
| 30 | **2,019** | 23,812 | **8.48%** |

The dx = 12 figure is exactly right. **The dx = 30 figure is not 89% but 8.48%,
and across a 2.5x change in grid the viable band moves by 12 cells out of
~2,010 — 0.6%, not 10x.** The inference drawn from the pair does not survive.
(The observation that the *tumour* set is a fixed measured population and the
viable band is a computed one remains true and is a structural argument; the
10x was not.)

**7.2 `F(0.1)` is 0.0 at dx = 30 µm, not merely small.** At the coarse grid no
viable tumour cell falls below a tenth of the perivascular reference at all, so
the headline is exactly zero rather than 0.0005. This is consistent with
`reach.py`'s own `CONVENTION_NOTE` reasoning about why the threshold was
withdrawn as a headline, and strengthens it.

---

## 8. Reproducing this

```
./.venv/bin/python experiments/reach_invariances_measure.py --stage guards
./.venv/bin/python experiments/reach_invariances_measure.py --stage all
./.venv/bin/python experiments/reach_invariances_measure.py --stage capacity_off
```

`--stage all` is ~50 min (6 barrier solves across 3 grids, 10 precondition
solves, 19 sweep solves; 3,022 s measured). Stages accumulate into the one JSON
rather than clobbering it. `--stage guards` runs in seconds and is the cheap way
to check the frozen state before spending an hour.

**No assertion was deleted or weakened, and nothing outside `experiments/` and
`results/` was written to.** The suite baseline is 1,221 passed / 0 failed at
`23e51f6`; this change adds two files and touches no module the suite imports.

**Rider disclosure — the working tree was not clean.** Sibling agents were live
in this repository throughout. Two of their edits touch code this study reads,
and both were inspected rather than waved past:

1. `app/server/reach.py` is being rewritten from a finite-window object to a
   steady-state one, and briefly did not import. **This study does not read the
   working-tree copy at all** — it loads the `23e51f6` blob and verifies it by
   SHA-256 (sec. 1.1), which is recorded in the JSON alongside the working-tree
   hash so the divergence is visible.
2. `tmesim/validation/vessel_profile.py` is dirty, and this study reads
   `distance_profile`, `usable_bins` and `direct_vessel_distance_um` from it.
   The diff is **purely additive**: two new output keys (`distance_mean_um`,
   `n_empty_bins`) and one new function
   (`anchor_corrected_half_distance_um`) that nothing calls. No pre-existing key
   and no existing code path changes, so no number here moves. Saying the tree
   was dirty and why it does not matter is cheaper than asserting it was clean.

Every constant this study leans on (`calibrate_alpha()`, the six `PDEConfig`
defaults, `REACH_LADDER`, `HEADLINE_X`, `OxygenConfig().hypoxia_threshold`) is
asserted against **live** code at the start of every stage.
