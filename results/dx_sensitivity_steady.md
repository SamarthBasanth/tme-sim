# Spatial-resolution sensitivity of the shipped steady-state reach verdict

**Result: the pre-registered threshold PASSED — and the pass is WEAK, for a
reason stated in §4 rather than buried.** The shipped headline `F(0.1)` is a
**0-or-1-cell statistic** on this section: it is exactly one viable tumour cell
out of ~2,000 at every rung from `dx = 24 µm` down, and zero cells at 30 µm.
There was almost nothing on that rung for the mesh to move. The only rung of the
same ladder that carries a real population for temozolomide, `F(1.0)`, exceeds
the same 0.01 threshold at **all three** refinement steps below the shipped
default (−0.0103, −0.0108, **+0.0142**) and is non-monotone across the ladder.
Both facts are reported; the second is **not** substituted for the first.

Driver: `experiments/dx_sensitivity_steady.py`.
Raw output: `results/dx_sensitivity_steady_data.json` (16 rungs + 1 repeat, none
dropped).
Section 1 below was written and committed to this file **before** the sweep ran
and has not been edited since.

> This file exists because a **fabricated "10x grid sensitivity" claim was
> published from this repository and retracted**. It is replaced here by a
> measurement. Every number below section 1 comes from a run this repository
> executed; the raw payload of every rung is in the JSON beside this file.

---

## 1. Pre-registration

Written and committed to this file **before any swept quantity was measured**.

### 1.1 What is being measured

The object is the **shipped** reach map. Verified from source rather than
assumed:

* `app/server/jobs.py::execute_run` calls
  `compute_reach(micro, res, distance_um=dist, penetration=dm["penetration_length"], drug=dp, convection=PDE_CONVECTION)`.
  `compute_reach` itself is defined in **`app/server/reach.py`**, not in
  `jobs.py` — `jobs.py` imports it (line 584) and calls it (line 1062).
* Inside `compute_reach`, `steady = drug is not None`. The production call
  supplies `drug=dp`, so the shipped map is the **steady-state** one:
  `window_mode == "steady_state"`, built on
  `reach.steady_state_concentration`, and `E_ref` is a concentration (uM), not
  an exposure (uM·s).
* A tissue-capacity term `R` enters as `R·dC/dt`. The steady problem
  `(K − L + Adv)c* = S` has no time derivative, so the capacity term is
  **absent from everything measured here** by inspection.
* In steady mode a caller-supplied `penetration` is **discarded**
  (`if steady or penetration is None:` re-fits against the steady profile), so
  omitting that argument in this harness leaves the payload unchanged.
  `penetration_reused` is recorded per rung to show it was `False` in
  production too.

The steady field carries **no dt and no window**, so nothing in the reach map is
subject to the advection CFL. The finite-window march is still run at every rung,
because `compute_reach` reads `sim_result.fields['oxygen']` (the viable band) and
`delivery_metrics.tumor_cell_exposure` reads `fields['auc']` — and because the
*other* shipped verdict, the per-cell dominant limiter from `sim.mechanism`,
lives there. Both verdicts are therefore read off the **same** run at every rung.

### 1.2 Section

`xenium_319` — `tmesim.data.section_catalog.DEFAULT_SECTION`, the section every
`GET /api/defaults` spec names, reconstructed through
`validation.vessel_profile.build_arm_microenvironment(dx_um, repair=0)`, which is
literally `reconstruct_real(source="xenium", dx_um=dx)` at the shipped
`oxygen_mode='pde'` / `state_mode='oxygen'` — the same call
`app/server/cache.py::build_reconstruction` makes for this section (whose
registry `sample` is `None`, the loader default).

### 1.3 Compounds

| compound | why it is on the list |
|---|---|
| `temozolomide` | the compound the Studio **opens on**; one of the 7 of 12 whose delivery/dose verdict **inverts** under the missing transvascular exchange term (`jobs.VERDICT_STABLE_COMPOUNDS` excludes it) |
| `doxorubicin` | **is** in `VERDICT_STABLE_COMPOUNDS`, and is the second half of the headline pair `experiments/steady_state_reach.py` uses |

Chosen for what they are on the shipped panel, not for their numbers. If
resolution moved only the already-fragile compound, that would be a different
finding from moving both.

**No cross-compound comparison is made or implied.** Cross-compound ranking is
withdrawn (2026-08-21, absolute). The two compounds are two *independent*
single-compound sensitivity measurements that happen to be reported in one
table.

### 1.4 `dx` values

`dx_um ∈ {30, 24, 20, 16, 12, 10, 8, 7}` µm.

* **12 µm is the shipped default** — `app/server/api.py::defaults()` returns
  `"dx_um": 12.0`.
* **30 µm is the hosted-instance value** — `app/server/api.py` lines 161–165
  record the measurement on a 512 MB instance (Render's free tier), flagship
  Xenium GBM, in a real container: `dx = 30 µm` is 46,070 voxels and peaks at
  220 MiB (comfortable), while `dx = 12 µm` is 286,200 voxels and the container
  is **OOM-killed**. `app/web/js/ui/spec.js` line 833 renders the same number to
  the user.
* **20 µm** is on the ladder because a constrained host that sets
  `TMESIM_MAX_GRID_CELLS = 100_000` (the value `api.py` documents) refuses
  `dx = 12` and permits `dx ≳ 20` on this section.
* **7 µm is the fine end** because `results/grid_convergence_gbm.md` measured the
  advection CFL on this section at the production `dt = 10 s` and found
  `dx = 7 → 0.96` and `dx = 6 → 1.13 REFUSED`. 7 µm is the finest mesh the
  shipped pipeline can reach **without also changing `dt`**, which would refine
  two discretisations at once. The CFL is re-measured here at every rung rather
  than trusted.

### 1.5 What is held fixed

Only `dx` varies. `dt = 10 s` (`jobs.PDE_DT_S`), 6 h window
(`api.STUDIO_DEFAULT_TOTAL_TIME_S`), `save_every = 500`, `seed = 0`,
`convection = True`, `transvascular_exchange = False` (the shipped default),
`PKConfig()` defaults, the `calibrated_drugparams → apply_barrier_to_drugparams`
compile path, and `alpha = 22.223477711172563` are all at their shipped values
and are asserted at the start and the end of every stage.
`sensitivity_mode='uniform'` sets `cells.sensitivity` to 1.0 everywhere and is a
no-op against the reconstruction default, so it is not called — a documented
equivalence, not an omission.

Every solve goes through `sim.engine.simulate_response`, the one production
wiring point that turns `PDEConfig.transvascular_exchange` into an `exchange_k`
argument. This file calls `simulate_drug_transport` **nowhere**, so the known
defect — a bare call with a hand-built `PDEConfig` silently making both arms of a
comparison exchange-OFF — cannot arise.

### 1.6 Output quantities, per rung

grid `(ny, nx)` and voxel count · tissue voxels · reconstruction wall-clock ·
peak RSS (each rung in **its own subprocess**, because `ru_maxrss` is a
monotone high-water mark and is only honest per-process) · max Darcy speed and
the resulting CFL at the production `dt` · `kill_fraction` · the full reach
ladder `F(x)` at `x ∈ {0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0}` · `E_ref` and its
units · `n_viable` / `n_tumor` · `window_mode` · the normalised distance falloff
· `survival_breakdown` and the **dominant limiter** · the delivery-vs-dose call
among exposure-limited survivors and `cells_to_flip`.

### 1.7 THE PRE-REGISTERED FAILURE THRESHOLD

The primary quantity is **`F(0.1)`** — `headline_under_reached_fraction`, the
fraction of viable tumour cells whose steady-state concentration is below
`0.1 × E_ref`. It is the number the reach panel renders, formatted by
`app/web/js/ui/reach.js` as `pct(v, 1) + '%'`, i.e. to one decimal place of a
percentage. **The panel's own display resolution is therefore 0.001.**

> **The shipped steady-state reach verdict is declared RESOLUTION-DEPENDENT — and
> the shipped default `dx = 12 µm` declared NOT on the converged side — if
> EITHER:**
>
> **(1)** `|F(0.1; dx_{k+1}) − F(0.1; dx_k)| > 0.01` for **any single refinement
> step at `dx ≤ 12 µm`** (12→10, 10→8, 8→7), for either compound; **or**
>
> **(2)** the **categorical** verdict — the dominant limiter from
> `sim.mechanism.survival_breakdown`, or the delivery-vs-dose call among
> exposure-limited survivors — **differs between any two rungs** of the ladder,
> for either compound.
>
> **Secondary, reported either way and NOT a substitute for the above:** the
> full-ladder span `max F(0.1) − min F(0.1)` over 30→7 µm, called LARGE above
> 0.05.
>
> **If (1) or (2) fires, the headline of this document says FAILED and no other
> metric is substituted for `F(0.1)`.**

**Why 0.01.** It is ten times the panel's own display resolution, so a step that
clears it is a movement a user can *see* on screen. It is also far above the
ladder's arithmetic resolution `1/n_viable` (order 1e-4 on this section), so it
cannot be a discretisation artefact of the fraction itself.

### 1.8 What would make this study vacuous, and the checks that catch it

| check | the defect it catches |
|---|---|
| **A** `window_mode == "steady_state"` at every rung, no nulls | a dropped `drug=` argument, which would silently make every rung a **finite-window** map of a different field while every number still looked plausible |
| **B** per-rung `(E_ref, F-ladder)` fingerprints pairwise distinct | a sweep that reused one reconstruction — it would report perfect convergence with every rung green |
| **C** the `dx = 12` temozolomide rung reproduces the published container `kill_fraction` `0.11221232991768859` and grid `286,200` (`docs/DEPLOYMENT.md`) bit-for-bit | a harness that is not the shipped pipeline, whose convergence statement would be about something else |
| **D** `calibrate_alpha()` asserted `== 22.223477711172563` at the START and the END of every stage | a drifted calibration, which would make every number here incomparable with `results/` |

---

## 2. Method

### 2.1 The first attempt was contaminated and was DISCARDED

It has to be said before any number is quoted.

The first run of this ladder went against the live checkout. This repository is
worked on by several agents at once, and between the `dx = 20` and `dx = 16`
rungs two files in the compile/delivery chain — `tmesim/reconstruct/barrier.py`
(`apply_barrier_to_drugparams`, the drug compile) and `tmesim/treatment/drugs.py`
(the delivery source) — were edited, and later reverted, **underneath the running
sweep**. The rung that straddled the edit came back with

```
dx = 20  kill_fraction = 0.112044   dominant limiter = dose
dx = 16  kill_fraction = 0.317655   dominant limiter = resistance     <- NOT A MESH EFFECT
```

a 2.8× jump in `kill_fraction` and a categorical flip of the dominant limiter in
a single refinement step. Published as it stood, that would have been a
spectacular grid-sensitivity finding **and it would have been false** — which is
precisely the failure mode this document exists to replace. That attempt was
thrown away, not repaired, and the whole ladder was re-run. The re-run returns
`kill_fraction = 0.112044`, limiter `dose`, at `dx = 16`.

### 2.2 Where the published sweep was run

Against a **frozen, self-consistent snapshot** of the committed tree, extracted
read-only with `git archive HEAD` — no `checkout`, no `stash`, no `reset`; the
working tree is never touched — into a scratch directory, with `data_cache`
symlinked back to the repository.

```
HEAD = a679c0f3ebf52729761262c8ae080cb9d5030ccd
```

That snapshot is also the right object on the merits: uncommitted in-flight edits
by other agents are not what ships. And the choice is falsifiable rather than
asserted — **check C** requires the snapshot to reproduce the container-measured
shipped `kill_fraction` bit-for-bit, which it does.

Freezing was not paranoia. While this sweep was running, `app/server/reach.py` in
the live checkout had `steady_state_concentration` renamed to
`steady_state_concentration_v2`. Every function name quoted in this document is
the name in the frozen snapshot at the SHA above.

Two mechanisms make a repeat of the contamination visible rather than silent:

* `source_fingerprint()` — a SHA-256 over every `.py` under `tmesim/` and
  `app/server/` **of the tree actually being imported** — is recorded on every
  rung, and **check E** requires them all equal;
* the **repeat gate** re-runs the first rung *after* the whole ladder and
  requires it bit-identical in `kill_fraction`, `F(0.1)`, `E_ref` and source
  hash.

### 2.3 What each rung does

One subprocess per `(dx, compound)`, because `ru_maxrss` is a monotone
high-water mark and is only honest per process. Inside it:

1. `VP.build_arm_microenvironment(dx_um=dx, repair=0)` — the production
   `reconstruct_real(source="xenium", dx_um=dx)` at the shipped
   `oxygen_mode='pde'` / `state_mode='oxygen'`.
2. `solve_ifp(micro)` — the max Darcy speed and hence the advection CFL the
   production `dt` will actually be held to. The same number
   `app/server/api.py`'s CFL guard reports.
3. `simulate_response(micro, SimConfig(drug=dp, pk=PKConfig(), pde=PDEConfig(dt=10, total_time=6 h, convection=True, save_every=500), seed=0), attribute=True)`
   — the one production wiring point. `simulate_drug_transport` is called
   nowhere in this file, so the known defect (a bare call with a hand-built
   `PDEConfig` silently making both arms of a comparison exchange-OFF) cannot
   arise.
4. `direct_vessel_distance_um(micro.domain, endothelial_xy(micro.cells))` — the
   same distance map `jobs.py` builds and hands to both consumers.
5. `compute_reach(micro, res, distance_um=dist, drug=dp, convection=True)` — the
   shipped call. `penetration` is the only argument dropped, and `compute_reach`
   **discards** it in steady-state mode anyway; `penetration_reused` is recorded
   per rung and is `false` at every rung, which is what it is in production too.

The compound is compiled by `calibrated_drugparams(name, alpha)` then
`apply_barrier_to_drugparams(..., get_molecule(name))` — the path the Studio's
own `drug_rows` and the shipped real-tissue panel both use. Nothing is scaled and
nothing is overridden.

### 2.4 The two compounds' length scales, for reading the tables

Bulk `L = sqrt(D/(clearance + binding))` off the compiled `DrugParams` (the
solver additionally applies the ECM tortuosity law, so the realised length is
shorter; the *ratio* is what matters here):

| compound | D (µm²/s) | clearance (1/s) | binding (1/s) | `L = sqrt(D/k)` (µm) |
|---|---|---|---|---|
| temozolomide | 42.42 | 5.794e-05 | 4.428e-04 | **291.1** |
| doxorubicin | 27.93 | 6.308e-05 | 4.972e-03 | **74.5** |

`k = 5.008e-04 1/s` for temozolomide reproduces the `k = 5.00761e-04 1/s` quoted
in `app/server/reach.py`'s own steady-state note. Temozolomide's decay length is
~3.9× doxorubicin's, and `reach.py` states that the steady map depends on the
compound **only** through this length — so the two rows below are two different
points on one functional form, not two independent compounds.

### 2.5 Checks — all six pass

| check | result | evidence |
|---|---|---|
| **A** `window_mode == "steady_state"`, no nulls | **PASS** | `window_modes_seen = ["steady_state"]`, `null_rungs = []`, over all 16 rungs |
| **B** rungs pairwise distinct | **PASS** | 16 distinct `(compound, E_ref, F-ladder)` fingerprints from 16 rows, 0 duplicate pairs |
| **C** `dx = 12` reproduces the shipped container run | **PASS** | `kill_fraction = 0.11221232991768859` exactly, grid `286,200` exactly (`docs/DEPLOYMENT.md`) |
| **E** one source tree produced every rung | **PASS** | one SHA-256, `b13ef1960e4c259863a575027b0ebb34fa375bc019cb72c4e6c273d02af7f388`, over 151 `.py` files, at every rung |
| **F** repeat gate bit-identical | **PASS** | `dx = 30` temozolomide re-run after the whole ladder: `kill_fraction`, `F(0.1)` and `E_ref` all bit-identical |
| **D** `calibrate_alpha()` | **PASS** | `22.223477711172563` asserted at the start and the end of every stage and every rung subprocess |

`calibrate_alpha()` at the START of the sweep: **`22.223477711172563`**.
`calibrate_alpha()` at the END of the sweep: **`22.223477711172563`**.

---

## 3. Raw results

Frozen snapshot of `HEAD = a679c0f3ebf52729761262c8ae080cb9d5030ccd`, generated
`2026-08-30T06:32:43Z`. 16 rungs, **none dropped**; the ladder is complete as
pre-registered.

#### 3.1 Cost, grid and stability

The grid, voxel, tissue-voxel, Darcy-speed and CFL columns depend only on
the reconstruction and are identical for both compounds. **The timing and
memory columns are the temozolomide arm** — they are not compound-independent,
and the doxorubicin arm is given underneath.

| dx (µm) | grid `ny × nx` | voxels | tissue voxels | max Darcy \|v\| (µm/s) | CFL at `dt = 10 s` | reconstruct (s) | march (s) | steady reach (s) | rung wall (s) | peak RSS (MiB) |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 170 × 271 | 46,070 | 30,976 | 0.2762 | 0.092 | 2.2 | 12.8 | 0.28 | 15.3 | 505 |
| 24 | 213 × 338 | 71,994 | 48,445 | 0.4992 | 0.208 | 1.7 | 29.5 | 0.93 | 32.2 | 481 |
| 20 | 255 × 406 | 103,530 | 69,769 | 0.4591 | 0.230 | 2.9 | 24.5 | 0.34 | 27.7 | 511 |
| 16 | 318 × 507 | 161,226 | 108,973 | 0.4934 | 0.308 | 1.7 | 40.7 | 0.61 | 43.0 | 691 |
| **12** *(shipped default)* | 424 × 675 | **286,200** | 193,752 | 0.6252 | 0.521 | 2.8 | 64.6 | 1.23 | 68.8 | 1,036 |
| 10 | 508 × 810 | 411,480 | 278,968 | 0.6276 | 0.628 | 3.2 | 88.3 | 5.04 | 96.8 | 1,324 |
| 8 | 635 × 1012 | 642,620 | 435,924 | 0.7334 | 0.917 | 8.0 | 205.6 | 6.54 | 220.4 | 1,648 |
| 7 | 725 × 1157 | 838,825 | 569,385 | 0.6668 | 0.953 | 12.7 | 222.6 | 9.28 | 244.9 | 1,833 |

The `30 µm` row's `46,070` voxels and the `12 µm` row's `286,200` are the two
numbers `app/server/api.py` lines 162–163 quote from the container measurement,
reproduced here exactly. The CFL at `dx = 7` is **0.953**, against the 0.955
`results/CLAIM_INVENTORY.md` N3 and `results/grid_convergence_gbm.md` §2.1 report
— an independent reproduction of the fine-end floor, and confirmation that
`dx = 6` is genuinely out of reach at the production `dt`.

The doxorubicin arm's own cost, on the same reconstructions:

| dx (µm) | 30 | 24 | 20 | 16 | 12 | 10 | 8 | 7 |
|---|---|---|---|---|---|---|---|---|
| rung wall (s) | 9.5 | 29.7 | 24.6 | 55.9 | 117.7 | 157.0 | 242.6 | 289.9 |
| peak RSS (MiB) | 543 | 592 | 594 | 575 | 834 | 1,174 | 1,645 | 2,048 |

**Wall-clock is not a benchmark.** The machine was running other agents' `pytest`
suites throughout (load average 4.6–14.6), so the timing column shows how cost
scales with the grid and nothing else. Peak RSS is per-process and unaffected by
that, but it is a development checkout's full interpreter and import graph and is
**not** comparable to the 220 MiB the container measurement reports for `dx = 30`.

**A correction to a hosted-instance number.** §1.4 pre-registered `20 µm` as
"permitted" under a `TMESIM_MAX_GRID_CELLS = 100_000` ceiling, following the
`~20 um` in `app/server/api.py`. Measured: this section's cell-centroid extent is
`8087.84 × 5065.35 µm`, so `api.py`'s own `(ex/dx)·(ey/dx)` gives **102,419** at
`dx = 20` — **2.4 % over** a 100,000 ceiling, i.e. **refused**. The exact
break-even is `dx = 20.24 µm`. The real `Domain` grid is larger still (103,530).
The `~` in the source comment carries this, but a host copying "20" would get a
400.

#### 3.2 `temozolomide` — the steady-state reach ladder

| dx (µm) | `window_mode` | F(0.01) | F(0.02) | F(0.05) | **F(0.1)** | F(0.2) | F(0.5) | F(1.0) | **ΔF(0.1)** | `E_ref` (µM) | `n_viable` | 1st falloff bin (µm) | its count | falloff @ 105 µm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | `steady_state` | 0.000000 | 0.000000 | 0.000000 | **0.000000** | 0.000000 | 0.000000 | 0.420505 | — | 7.80242 | 2,019 | 5 | 1,376 | 0.7101 |
| 24 | `steady_state` | 0.000495 | 0.000495 | 0.000495 | **0.000495** | 0.000495 | 0.000495 | 0.379327 | +0.000495 | 7.68917 | 2,022 | 5 | 2,248 | 0.6957 |
| 20 | `steady_state` | 0.000498 | 0.000498 | 0.000498 | **0.000498** | 0.000498 | 0.000498 | 0.390547 | +0.000003 | 7.70668 | 2,010 | 5 | 3,184 | 0.6950 |
| 16 | `steady_state` | 0.000499 | 0.000499 | 0.000499 | **0.000499** | 0.000499 | 0.000499 | 0.382045 | +0.000001 | 7.69775 | 2,005 | 5 | 5,044 | 0.6936 |
| **12** | `steady_state` | 0.000498 | 0.000498 | 0.000498 | **0.000498** | 0.000498 | 0.000498 | 0.373008 | −0.000001 | 7.65849 | 2,008 | 5 | 8,875 | 0.6954 |
| 10 | `steady_state` | 0.000500 | 0.000500 | 0.000500 | **0.000500** | 0.000500 | 0.000500 | 0.362681 | +0.000002 | 7.65340 | 1,999 | 5 | 12,720 | 0.6964 |
| 8 | `steady_state` | 0.000501 | 0.000501 | 0.000501 | **0.000501** | 0.000501 | 0.000501 | 0.351852 | +0.000000 | 7.62429 | 1,998 | 5 | 19,919 | 0.6978 |
| 7 | `steady_state` | 0.000501 | 0.000501 | 0.000501 | **0.000501** | 0.000501 | 0.000501 | 0.366049 | +0.000000 | 7.64577 | 1,997 | 5 | 26,057 | 0.6946 |

#### 3.3 `temozolomide` — the dominant-limiter verdict (finite-window march, same run)

| dx (µm) | dominant limiter | its share | delivery | dose | resistance | chance | delivery : dose counts | delivery/dose verdict | cells to flip | `kill_fraction` |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | **dose** | 0.998534 | 0.001466 | 0.998534 | 0.000000 | 0.000000 | 31 : 21,111 | **dose** | 10,541 | 0.112128 |
| 24 | **dose** | 0.998345 | 0.001655 | 0.998345 | 0.000000 | 0.000000 | 35 : 21,110 | **dose** | 10,538 | 0.112002 |
| 20 | **dose** | 0.998439 | 0.001561 | 0.998439 | 0.000000 | 0.000000 | 33 : 21,111 | **dose** | 10,540 | 0.112044 |
| 16 | **dose** | 0.997966 | 0.002034 | 0.997966 | 0.000000 | 0.000000 | 43 : 21,101 | **dose** | 10,530 | 0.112044 |
| **12** | **dose** | 0.997824 | 0.002176 | 0.997824 | 0.000000 | 0.000000 | 46 : 21,094 | **dose** | 10,525 | **0.112212** |
| 10 | **dose** | 0.997777 | 0.002223 | 0.997777 | 0.000000 | 0.000000 | 47 : 21,093 | **dose** | 10,524 | 0.112212 |
| 8 | **dose** | 0.997730 | 0.002270 | 0.997730 | 0.000000 | 0.000000 | 48 : 21,093 | **dose** | 10,523 | 0.112170 |
| 7 | **dose** | 0.997682 | 0.002318 | 0.997682 | 0.000000 | 0.000000 | 49 : 21,091 | **dose** | 10,522 | 0.112212 |

#### 3.4 `doxorubicin` — the steady-state reach ladder

| dx (µm) | `window_mode` | F(0.01) | F(0.02) | F(0.05) | **F(0.1)** | F(0.2) | F(0.5) | F(1.0) | **ΔF(0.1)** | `E_ref` (µM) | `n_viable` | 1st falloff bin (µm) | its count | falloff @ 105 µm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | `steady_state` | 0.000000 | 0.000000 | 0.000000 | **0.000000** | 0.003962 | 0.280832 | 0.893512 | — | 0.213346 | 2,019 | 5 | 1,376 | 0.2207 |
| 24 | `steady_state` | 0.000495 | 0.000495 | 0.000495 | **0.000495** | 0.002967 | 0.274975 | 0.880811 | +0.000495 | 0.209931 | 2,022 | 5 | 2,248 | 0.2128 |
| 20 | `steady_state` | 0.000498 | 0.000498 | 0.000498 | **0.000498** | 0.003483 | 0.276617 | 0.884080 | +0.000003 | 0.210384 | 2,010 | 5 | 3,184 | 0.2123 |
| 16 | `steady_state` | 0.000499 | 0.000499 | 0.000499 | **0.000499** | 0.004489 | 0.264339 | 0.883791 | +0.000001 | 0.209550 | 2,005 | 5 | 5,044 | 0.2116 |
| **12** | `steady_state` | 0.000498 | 0.000498 | 0.000498 | **0.000498** | 0.002490 | 0.262948 | 0.875996 | −0.000001 | 0.206967 | 2,008 | 5 | 8,875 | 0.2128 |
| 10 | `steady_state` | 0.000500 | 0.000500 | 0.000500 | **0.000500** | 0.003002 | 0.259630 | 0.874937 | +0.000002 | 0.206484 | 1,999 | 5 | 12,720 | 0.2131 |
| 8 | `steady_state` | 0.000501 | 0.000501 | 0.000501 | **0.000501** | 0.003003 | 0.253754 | 0.870871 | +0.000000 | 0.205297 | 1,998 | 5 | 19,919 | 0.2144 |
| 7 | `steady_state` | 0.000501 | 0.000501 | 0.000501 | **0.000501** | 0.003005 | 0.254382 | 0.876314 | +0.000000 | 0.205728 | 1,997 | 5 | 26,057 | 0.2133 |

#### 3.5 `doxorubicin` — the dominant-limiter verdict (finite-window march, same run)

| dx (µm) | dominant limiter | its share | delivery | dose | resistance | chance | delivery : dose counts | delivery/dose verdict | cells to flip | `kill_fraction` |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | **resistance** | 0.980217 | 0.002413 | 0.003378 | 0.980217 | 0.013993 | 30 : 42 | **dose** | **7** | 0.477784 |
| 24 | **resistance** | 0.979431 | 0.002491 | 0.003857 | 0.979431 | 0.014221 | 31 : 48 | **dose** | **9** | 0.477322 |
| 20 | **resistance** | 0.979912 | 0.002411 | 0.003696 | 0.979912 | 0.013982 | 30 : 46 | **dose** | **9** | 0.477364 |
| 16 | **resistance** | 0.979357 | 0.002329 | 0.004337 | 0.979357 | 0.013976 | 29 : 54 | **dose** | **13** | 0.477154 |
| **12** | **resistance** | 0.979367 | 0.002569 | 0.004175 | 0.979367 | 0.013889 | 32 : 52 | **dose** | **11** | 0.476902 |
| 10 | **resistance** | 0.978804 | 0.002569 | 0.004496 | 0.978804 | 0.014131 | 32 : 56 | **dose** | **13** | 0.476944 |
| 8 | **resistance** | 0.978821 | 0.002647 | 0.004412 | 0.978821 | 0.014120 | 33 : 55 | **dose** | **12** | 0.476524 |
| 7 | **resistance** | 0.978737 | 0.002728 | 0.004493 | 0.978737 | 0.014042 | 34 : 56 | **dose** | **12** | 0.476608 |

---

## 4. Verdict against the pre-registration

### 4.1 The result

> ### **PRE-REGISTERED THRESHOLD: PASSED.**
>
> Neither criterion fired, for either compound.

**Criterion (1) — `|ΔF(0.1)| > 0.01` at any refinement step below `dx = 12 µm`:
DID NOT FIRE.** The three fine steps are, for *both* compounds identically:

| step | ΔF(0.1) | vs threshold 0.01 |
|---|---|---|
| 12 → 10 µm | +2.242e-06 | 4,460× below |
| 10 → 8 µm | +2.504e-07 | 39,900× below |
| 8 → 7 µm | +2.506e-07 | 39,900× below |

Largest step anywhere on the ladder: **+4.946e-04**, at 30 → 24 µm — still 20×
below the threshold, and above the shipped default in any case.

**Criterion (2) — categorical verdict differs between any two rungs: DID NOT
FIRE.** Over all eight rungs, from 30 µm to 7 µm:

* temozolomide — dominant limiter `dose` at every rung; delivery/dose verdict
  `dose` at every rung.
* doxorubicin — dominant limiter `resistance` at every rung; delivery/dose
  verdict `dose` at every rung.
* `window_mode` is `steady_state` at all 16 rungs, as the shipped call requires.

**Secondary (span, reported either way):** full-ladder span of `F(0.1)` is
**5.008e-04** for both compounds — 100× below the 0.05 "LARGE" band. NOT LARGE.

### 4.2 Why this PASS is weak, stated rather than buried

`F(x) = count / n_viable`. On this section `n_viable` is ~2,000, so one cell is
0.0005. Reading the counts back off the table:

| dx (µm) | 30 | 24 | 20 | 16 | 12 | 10 | 8 | 7 |
|---|---|---|---|---|---|---|---|---|
| viable cells below `0.1 × E_ref` | **0** | **1** | **1** | **1** | **1** | **1** | **1** | **1** |

**The shipped headline is a 0-or-1-cell statistic on this section.** The entire
movement of `F(0.1)` across a 4.3× mesh refinement is *one cell appearing between
30 and 24 µm*; every later "change" is the same single cell divided by an
`n_viable` that drifts from 2,022 to 1,997. For temozolomide the same is true of
rungs `F(0.01)` through `F(0.5)` — **six of the seven rungs are the same one
cell** — because its steady decay length (`L = sqrt(D/k) ≈ 291 µm` in bulk) is
long against this section's vessel spacing, so almost no viable tumour cell sits
below half the perivascular reference, let alone a tenth of it.

A test can only fail on a quantity that can move. This one barely can. The pass
is therefore **real but low-information**: it says the headline is not
mesh-sensitive, and it says so on a rung where a mesh effect would have had to
move a handful of cells to be visible at all.

### 4.3 The rung that DOES carry a population moves more than the threshold

Reported as a **separate, secondary observation about a different rung**. It is
**not** substituted for the pre-registered result above, which stands as PASSED.

`F(1.0)` — the fraction of viable tumour cells below the perivascular reference
itself — is the one rung with a real population for temozolomide (703–849 cells):

| dx (µm) | 30 | 24 | 20 | 16 | 12 | 10 | 8 | 7 |
|---|---|---|---|---|---|---|---|---|
| temozolomide `F(1.0)` | 0.420505 | 0.379327 | 0.390547 | 0.382045 | 0.373008 | 0.362681 | 0.351852 | 0.366049 |
| cells | 849 | 767 | 785 | 766 | 749 | 725 | 703 | 731 |

Span **0.0687**. Largest single step **0.0412** (30 → 24 µm). And **all three**
refinement steps below the shipped default exceed the 0.01 that criterion (1)
uses:

| step | ΔF(1.0) | vs 0.01 |
|---|---|---|
| 12 → 10 µm | **−0.010327** | exceeds |
| 10 → 8 µm | **−0.010829** | exceeds |
| 8 → 7 µm | **+0.014197** | exceeds |

The sequence is also **non-monotone** — it falls, rises, falls for four rungs,
then rises again at the finest step, reversing sign exactly where a converging
ladder should be settling. That is not an asymptotic regime, which is the same
thing `results/grid_convergence_gbm.md` reports for the finite-window quantities
on this section.

Had the pre-registration named `F(1.0)` instead of `F(0.1)`, criterion (1) would
have fired at every fine step and this document's headline would read **FAILED**.
It named `F(0.1)`, because `F(0.1)` is what the product prints; the threshold is
not moved and the result is not re-scored.

For doxorubicin the populated rungs move less — and stay **under** the threshold:
`F(0.5)` spans 0.0271 with a largest fine step of 0.0059, `F(1.0)` spans 0.0226
with a largest fine step of 0.0054. Both are still three to four orders of
magnitude above the fine-step movement of `F(0.1)`.

**Read together:** the mesh does move this map. It moves the part of the ladder
where cells actually are. It does not move the rung the product prints, because
almost nothing is on that rung.

### 4.4 Two other things the ladder moved

**`E_ref` drifts monotonically with refinement** — temozolomide 7.80242 →
7.64577 µM (**−2.01 %**), doxorubicin 0.213346 → 0.205728 µM (**−3.57 %**), over
30 → 7 µm. Every ratio in the payload is divided by this, so a −3.6 % drift in
the denominator is present in every published falloff number even where the
ratios themselves barely move. The normalised falloff at 105 µm moves only
−2.2 % (temozolomide 0.7101 → 0.6946) and −3.4 % (doxorubicin 0.2207 → 0.2133),
so numerator and denominator largely drift together — which is the scale-free
design working, and is worth saying explicitly rather than leaving as luck.

**`n_viable` itself is mesh-dependent**: 2,019 / 2,022 / 2,010 / 2,005 / 2,008 /
1,999 / 1,998 / 1,997 — a 1.2 % drift, non-monotone. The viable band is defined
by the oxygen field, which is solved on the grid, so the *denominator* of every
`F(x)` moves with `dx` independently of anything the drug does.

**`kill_fraction` barely moves at all**: 0.112002–0.112212 (span 0.000210,
**0.19 %**) for temozolomide and 0.476524–0.477784 (span 0.001260, **0.26 %**)
for doxorubicin, across the whole 4.3× refinement. That is a narrow spread and it
is reported as measured — it is **not** a retraction of the withdrawn
grid-convergence claim (`results/CLAIM_INVENTORY.md` N2), which was withdrawn on
*other* quantities, on a *breast* section, and for non-monotonicity rather than
for spread.

### 4.5 The doxorubicin verdict is stable here and one bad coin-flip from not being

The delivery/dose call is `dose` at every rung — but the counts are `30:42`,
`31:48`, `30:46`, `29:54`, `32:52`, `32:56`, `33:55`, `34:56`, so **`cells_to_flip`
is 7 to 13 at every rung**. Criterion (2) passed on a margin of roughly a dozen
cells out of 40,887. That is a fact about how thin the verdict is on this
compound, not about the mesh, and it is consistent with
`results/verdict_inversion_census.md`, which found doxorubicin sitting at
0.497 : 0.503 under the exchange-ON law. For temozolomide the same margin is
10,522–10,541 cells, i.e. not close.

---

## 5. Convergence: is there a `dx` beyond which the verdict stops moving, and is the shipped default on the converged side?

**For the categorical verdicts: yes, and the shipped default is on it.** Neither
the dominant limiter nor the delivery/dose call moves anywhere on the ladder —
`dose` for temozolomide and `resistance`/`dose` for doxorubicin from 30 µm all
the way to 7 µm. There is no rung at which the answer changes, so `dx = 12` is
not distinguishable from any other tested mesh on this axis. **Note what that
does *not* mean:** a verdict that never moves over the reachable range has not
been shown to be converged, only to be insensitive over that range — and for
doxorubicin it is insensitive while sitting a dozen cells from reversing.

**For the pre-registered headline `F(0.1)`: it stops moving at `dx = 24 µm`, and
the shipped default is well inside that.** After the single cell appears between
30 and 24 µm the count never changes again; the residual movement is `1/n_viable`
drift of order 1e-6 per step. But see §4.2 — this is convergence of a statistic
that holds one cell.

**For the populated part of the ladder: NO, and `dx = 12 µm` is NOT on a
converged side.** Temozolomide's `F(1.0)` moves more than 0.01 at *every*
refinement step below `dx = 12`, including 0.0142 at the finest step the shipped
`dt` permits, and it changes sign there. There is no `dx` in
`{30, 24, 20, 16, 12, 10, 8, 7}` beyond which it settles, and the ladder cannot
be extended — `dx = 6 µm` violates the advection CFL at the production
`dt = 10 s` (measured 0.953 at 7 µm here). **A converged answer for that rung is
not reachable without also changing `dt`, and this study deliberately did not.**

So the honest three-line summary of convergence:

1. the shipped *categorical* verdict does not move on any mesh this pipeline can
   run;
2. the shipped *headline number* does not move either, but only because it is one
   cell;
3. the part of the same object that has a population **is not converged at the
   shipped default and is not converged at the finest reachable mesh**.

---

## 6. What this does NOT show

**It is not a continuum-limit convergence proof.** The ladder stops at
`dx = 7 µm` because the advection CFL at the production `dt = 10 s` refuses
`dx = 6` (measured here at 0.953 for 7 µm, reproducing the 0.955 in
`results/CLAIM_INVENTORY.md` N3 and `results/grid_convergence_gbm.md` §2.1).
"Stops moving" below can therefore only ever mean *stops moving over 30 → 7 µm*.
No Richardson limit, no observed order and no GCI is claimed here — that is a
different study, and `results/grid_convergence_gbm.md` already reports that on
this section such fits come back **non-monotone and degenerate** for the
finite-window kill quantities.

**It is a SPATIAL study only.** `dt` is held at the production 10 s at every
rung, deliberately: scaling `dt` with `dx` would refine two discretisations at
once and the resulting curve could not be attributed to the mesh. The cost is
that the temporal discretisation error is common to every rung — it cancels out
of the movement *between* rungs and is inherited by every absolute number in the
tables.

**It is one section, one platform, one mask.** `xenium_319`, 10x Xenium,
`repair = 0` (the native reconstructed mask). Nothing here transfers to another
section, and this repository has already published the case where such a
transfer failed: on breast, docetaxel's verdict flips FAIL → PASS between
`dx = 12` and `dx = 16` (`results/grid_convergence_breast.md`), on a quantity
this study does not measure.

**It is two compounds, and they are not compared with each other.**
Cross-compound ranking is withdrawn (2026-08-21, absolute). The two rows are two
independent single-compound measurements printed in one table for convenience.
They are also not two independent probes of the model: `app/server/reach.py`
states that the steady map depends on the compound **only** through
`L = sqrt(D_eff/k)`, so this is one functional form sampled at two lengths, and
it cannot see a permeability barrier at all.

**The mesh moves the whole pipeline, not just the transport step.** The
geometric tissue mask is `dx`-independent by construction, but the ECM and
vascular KDE fields, the oxygen solve, the IFP/velocity solve and the cell-state
banding all live on the grid — and so does `n_viable`, the *denominator* of
`F(x)`. A movement in `F` is therefore never attributable to the transport
discretisation alone. `n_viable` is tabulated per rung for exactly this reason.

**`E_ref` is a moving reference.** `vessel_profile.BIN_UM = 10 µm` is FIXED and
independent of `dx`, and its own source comment describes it as "comfortably
above dx = 12 µm" — which it is not. `E_ref` is the mean of the first distance
bin holding at least `MIN_BIN_COUNT = 30` grid points, so both the number of
points in that bin and, in principle, which physical bin qualifies can move with
`dx`. Every rung's first usable bin centre and its count are tabulated so this
can be checked rather than assumed; on this section the usable bins were
contiguous from 5 µm at every rung, so the comb-aliasing failure mode
`distance_profile` warns about did not occur.

**Only the shipped delivery law.** `transvascular_exchange = False` at every
rung. `results/verdict_inversion_census.md` measured that switching that term ON
inverts the delivery/dose verdict for 7 of 12 shipped compounds, temozolomide
included. That axis is NOT crossed with `dx` here, so nothing below says how the
two interact.

**Nothing about tissue capacity.** A capacity term multiplies `dC/dt` and the
steady problem has no time derivative, so capacity is absent from the reach map
by construction. This study could not have detected a capacity effect and does
not claim to have ruled one out.

**Wall-clock is not a benchmark.** The sweep ran on a laptop with other agents'
`pytest` runs competing for CPU (load average 4.6–14.6 during the sweep), so the
wall-clock column indicates how cost *scales* with the grid and nothing more.
Peak RSS is per-process and unaffected by that contention, but it includes the
whole interpreter and import graph of a development checkout — it is **not**
comparable to the 220 MiB the container measurement in `app/server/api.py`
reports for `dx = 30`.

**No clinical reading.** `kill_fraction` appears only as a harness gate (check C)
and as a covariate. Absolute kill fraction is not a defensible output of this
project, no cross-compound ordering is stated, and nothing here is a treatment,
diagnostic or response prediction.

**And specifically: it does not replace a convergence study.** The retracted
"10× grid sensitivity" claim this file exists to replace was a claim about
magnitude. This file does not supply a corrected magnitude for it, because the
quantity that would carry one (`F(1.0)`) is non-monotone and unconverged over the
reachable range, and quoting a span from an unconverged non-monotone ladder as
"the grid sensitivity" would be the same error with a different number.

---

## 7. Reproducing this

```
git archive HEAD | tar -x -C <snap>          # a679c0f3ebf52729761262c8ae080cb9d5030ccd
ln -s <repo>/data_cache <snap>/data_cache
cd <snap> && PYTHONPATH=<snap> <repo>/.venv/bin/python \
    experiments/dx_sensitivity_steady.py --stage probe
cd <snap> && PYTHONPATH=<snap> <repo>/.venv/bin/python \
    experiments/dx_sensitivity_steady.py --stage sweep
cd <snap> && PYTHONPATH=<snap> <repo>/.venv/bin/python \
    experiments/dx_sensitivity_steady.py --stage report
```

`--stage sweep` takes roughly 45 minutes on an unloaded machine and about 75 on
the contended one this was measured on. It writes
`results/dx_sensitivity_steady_data.json` and nothing else. Every rung is a
separate subprocess and each prints its own payload after a `===RUNG_JSON===`
marker, so a single rung can be re-run and diffed on its own:

```
... --stage rung --dx 12 --compound temozolomide
```

**This document is not covered by `tests/test_claim_guard.py`.** Registering it
would mean adding entries to that file's `CLAIMS` tuple — a `Const` on
`calibrate_alpha`, a `Default` on `PDEConfig.transvascular_exchange`, and a
`CallEdge` for `app.server.jobs → app.server.reach` — which is outside the file
set this study was allowed to write. It is flagged here rather than left
unnoticed.
