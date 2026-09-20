# Where a Studio run's time and memory actually go

**2026-09-20.** A profile of the exact path `app/server/jobs.py::execute_run`
takes, at the shipped default `dx = 12 um` and at the hosted-instance
`dx = 30 um`, on the flagship `xenium_319` section.

**Nothing in this file was optimised.** It is a measurement and a ranked list of
candidates. Every proposed change in section 6 is written as a diff and NOT
applied, because the files they touch belong to other work in flight.

---

## 0. Read this before quoting any number here

### 0.1 A performance claim is a claim

This repository has published a performance claim that was wrong in its
**frame**, not merely in its digits: "60x cheaper, 4.4 s vs ~250 s" was really
48.6x (0.90 s vs 43.70 s), and the full run it was meant to speed up got **4.2 %
slower**. So, explicitly:

* Every timing below comes from a command that was executed, and the command is
  pasted in section 8.
* No speedup is quoted without its effect on the end-to-end time a user waits.
* Where a change would trade time for memory or memory for time, both
  directions are stated. One of the five candidates in section 6 makes the run
  **slower** (#3), and it is still on the list, because the constraint it
  relieves is an OOM kill rather than a wait.

### 0.2 The machine was loaded, and that dominates the absolutes

Measured on an Apple-silicon macOS host, 10 logical CPUs, 16 GiB. Other agents
were running the whole time. `os.getloadavg()` during this session ranged from
**~4.6 to ~48 on 10 CPUs**.

The same RunSpec, on the same code, measured five different wall-clocks:

| arm | dx | wall (s) | cpu u+s (s) | loadavg (1 min) |
|---|---|---|---|---|
| `stages`, cold, first measurement of the session | 12 | **47.083** | not recorded | not recorded |
| `stages`, warm cache, immediately after | 12 | 58.131 | not recorded | not recorded |
| `instrument` arm A, cold | 12 | 69.576 | not recorded | not recorded |
| `instrument` arm B, warm, shimmed | 12 | 103.479 | not recorded | not recorded |
| `sweep` (fresh subprocess) | 12 | **107.909** | 99.846 | 6.81 |
| `redundant-solve` instrumented run | 12 | 174.270 | 142.44 | 26.67 |

"not recorded" is literal: the CPU-time and loadavg stamps were added to the
harness after those first four runs, and they are not being back-filled from
memory. Every arm run after that change carries both.

A **2.3x spread on identical work**, and 3.7x once the box hit load 27. Note
that **CPU seconds moved too** — 105.8 s of CPU for two `dx = 12` runs early on
(shell `time`: 93.98 s user + 11.79 s system across a 105.2 s two-run command),
versus 99.8 s of CPU for a single one later. On an Apple-silicon host under
load the scheduler migrates work onto efficiency cores, so the same
instructions cost more CPU-seconds. Neither wall-clock nor CPU-time is a stable
absolute here.

**What is stable, and what this file therefore leads with, is the SHAPE** — the
fraction of the run each stage takes, the ratio between the factorisation and
the back-substitutions, the count of redundant solves, and the bytes-per-voxel
slope of peak RSS. Those are reproducible across every arm above. The absolutes
are quoted with the loadavg attached, and should be re-measured on a quiet box
before anyone plans around them.

This is also why `tests/test_performance_invariants.py` contains no wall-clock
assertion at all. See its module docstring.

### 0.3 Unit calibration for the memory numbers

`resource.getrusage(RUSAGE_SELF).ru_maxrss` is **bytes on Darwin** and
kilobytes on Linux. Calibrated here rather than assumed — allocating a
400,000,000-byte array moved it by 400,031,744:

```
$ .venv/bin/python -c "
import resource, numpy as np
b=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
a=np.ones(50_000_000); a[::1000]=2
c=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print('delta', c-b)"
delta 400031744
```

All RSS figures below are **peak** RSS (high-water mark), which is the quantity
an OOM killer acts on, not current RSS.

### 0.4 Code state

Branch `main` at `f9d443d`, on a working tree with several other agents editing
concurrently. SHA-256 of the four files that carry every number below, **taken
at the time of the measurements**:

```
233b9664272953c4f5b94594c4842ad1c9f296f5e42e83cfb1ff4a4f3e3f50ea  tmesim/physics/diffusion.py
2d5964258a14f195fb09782fff41d8c0f5f5345cd6cd92017a4bbbd909a589b4  tmesim/sim/engine.py
828089444e0828b0350451e01e38356362cba061c3c7c24e5a8bd7fe455ee0d9  app/server/jobs.py
0e4d9e17324febe220cb11535d6f8f8ece501a8ef072e6a3691c33b25fb523b2  app/server/cache.py
```

`tmesim/physics/diffusion.py` and `tmesim/sim/engine.py` still hash to those
values at the time of writing — they were never touched.

**`app/server/jobs.py` and `app/server/cache.py` DID change after the
measurements**, by other work in flight, to
`46860b3ad74f8857740b5e25db5adaa68eda16c432a5d36835fbbf5ba9870586` and
`b250318e15230f5d596639b82af11cdc098c858f03d7f41ca06862d3445c351c`. The change
adds `"cisplatin"` to `VERDICT_STABLE_COMPOUNDS` and its accompanying
disclosure block: **a frozenset literal and comments, on no hot path, for a
compound this profile did not run.** It does not move any number here — but it
is recorded rather than glossed, because "the hashes matched" is the kind of
sentence that is true when it is written and false when it is read.

`app/server/api.py` (which owns `default_spec()`) was likewise being edited.
The spec it returned is printed verbatim by every arm of the harness and is
reproduced in section 1, so a later change to it is visible as a changed spec
rather than as a silent divergence.

`calibrate_alpha()` asserted equal to `22.223477711172563` at the start and the
end of this work (section 9).

---

## 1. The spec, and the reconciliation with `docs/DEPLOYMENT.md`

The harness calls `app.server.api.default_spec()` and overrides only `dx_um`
(and, for the section sweep, `section`). At `dx = 12` that is verbatim:

```json
{"attribute_mechanism": true, "drug": "temozolomide", "dx_um": 12.0,
 "mgmt": null, "oxygen_mode": "pde",
 "pk": {"dose_interval_s": 86400.0, "half_life_s": 14400.0,
        "n_doses": 1, "plasma_peak": 50.0},
 "section": "xenium_319", "sensitivity_mode": "uniform",
 "state_mode": "oxygen", "total_time_s": 21600.0,
 "vascular_mode": "discrete"}
```

`docs/DEPLOYMENT.md` line 277 records **`runtime_s` 47.13 cold** for this spec,
on a 4-CPU VM-backed daemon in a container, with
`kill_fraction = 0.11221232991768859` on 40,887 cells.

**The first cold measurement of this session, on this host:**

```
grid: 424x675 = 286,200 voxels
kill_fraction: 0.11221232991768859
TOTAL 47.083 s
```

* The grid matches `docs/DEPLOYMENT.md`'s own 286,200-voxel figure exactly.
* The `kill_fraction` is **bit-identical** to the containerised one, to all 17
  digits. That is the reconciliation that matters: the two runs computed the
  same thing.
* The wall-clock agrees to 0.1 %, **and that agreement is partly luck.** Per
  §0.2 the same spec measured 107.9 s later in the same session. The honest
  statement is: *47.13 s in that container and 47.08 s here are consistent, and
  this host's spread for the identical spec is 47–108 s.* Do not read either
  number as a benchmark.

---

## 2. Stage breakdown

Stage boundaries are the server's own `ctx.stage()` transitions, timestamped by
a `_Ctx` stand-in; they are not a re-implementation of the stage list.
`pre-stage` is the block of function-local imports at the top of `execute_run`.

### 2.1 `dx = 12 um` — the shipped default (286,200 voxels)

Clean (no profiler, no shims), cold reconstruction cache:

| stage | seconds | fraction |
|---|---|---|
| pre-stage (imports) | 0.713 | 1.5 % |
| **reconstruct** | 3.243 | 6.9 % |
| sensitivity | 0.002 | 0.0 % |
| compile_drug | 0.894 | 1.9 % |
| **simulate** | **40.733** | **86.5 %** |
| metrics | 1.468 | 3.1 % |
| payloads | 0.029 | 0.1 % |
| **TOTAL** | **47.083** | |

Warm (reconstruction served from cache), immediately after, same process:

| stage | seconds | fraction |
|---|---|---|
| reconstruct | 0.000 | 0.0 % |
| compile_drug | 0.823 | 1.4 % |
| **simulate** | 55.606 | 95.7 % |
| metrics | 1.670 | 2.9 % |
| payloads | 0.030 | 0.1 % |
| **TOTAL** | 58.131 | |

The warm run being *slower* than the cold one is machine noise, not an effect:
see §0.2. What the warm column does establish is that the cache works — the
reconstruction stage goes to zero — and that **caching the reconstruction
removes at most 6.9 % of a `dx = 12` run.** The reconstruction cache is a
correctness and memory device far more than a speed one.

### 2.2 `dx = 30 um` — the hosted-instance resolution (46,070 voxels)

Clean, cold:

| stage | seconds | fraction |
|---|---|---|
| pre-stage (imports) | 0.560 | 8.3 % |
| reconstruct | 0.695 | 10.3 % |
| sensitivity | 0.001 | 0.0 % |
| compile_drug | 0.718 | 10.6 % |
| **simulate** | **4.650** | **68.5 %** |
| metrics | 0.156 | 2.3 % |
| payloads | 0.004 | 0.1 % |
| **TOTAL** | **6.783** | |

`kill_fraction = 0.11212833865278016`, grid 170x271.

At the coarse resolution the fixed costs (imports, drug compilation) are a
third of the run; at `dx = 12` they are under 4 %. **`dx` moves the run by
6.9x (6.78 s -> 47.08 s) for a 6.2x voxel increase** — very close to linear,
which §3 explains.

### 2.3 The six stages the brief asked for

The server's stage list does not separate transport from kill from mechanism —
they are all inside `simulate`. `cProfile` does. Cumulative time by function,
from the profiled `dx = 12` run (profiled total 92.181 s; the tracing overhead
is real and is why the *fractions*, not the seconds, are the result):

These are **cumulative** times and some of them nest, so they do not sum to
100 %; the indentation says which contains which. The six top-level rows
(`build_reconstruction`, `compile_drug`, `simulate_response`, the metrics
block, payloads) are disjoint and do.

| requested stage | function | cum s | % of profiled run |
|---|---|---|---|
| **reconstruction** | `cache.build_reconstruction` | 4.806 | 5.2 % |
| ├─ *field building* | `reconstruct_real` -> `solve_oxygen` #1 | 3.815 | 4.1 % |
| └─ *field building* | `_max_darcy_speed` -> `solve_ifp` #1 | 0.985 | 1.1 % |
| drug compilation | `compile_drug` stage (`calibrate_alpha` 0.655) | 1.459 | 1.6 % |
| **simulate** | `sim.engine.simulate_response` | 81.941 | 88.9 % |
| ├─ *field building* | `solve_oxygen` #2 (re-solve) | ~1.4 | ~1.5 % |
| ├─ *field building* | `solve_ifp` #2 | ~1.1 | ~1.2 % |
| ├─ **transport solve** | `simulate_drug_transport` | **78.109** | **84.7 %** |
| ├─ **kill evaluation** | `apply_drug_kill` (x7) + `_viable_history` | 0.033 | **0.04 %** |
| └─ **mechanism attribution** | `sim.mechanism.attribute` | 0.003 | **0.003 %** |
| metrics | `direct_vessel_distance_um` + `compute_all` + `compute_reach` | 2.892 | 3.1 % |
| └─ *field building* | `solve_ifp` #3, inside `compute_reach` | ~1.2 | ~1.3 % |
| **payload encoding** | `_build_payloads` + `encode_*` | 0.059 | **0.06 %** |

(`solve_oxygen` totals 5.230 s cum over its 2 calls and `solve_ifp` 3.282 s over
its 3; the per-call splits marked `~` are those totals apportioned by the
directly-measured per-call times in §6 #4, since `cProfile` reports a function's
calls in aggregate. **Field building, added up across all five calls, is 8.512 s
cum — 9.2 % of the profiled run — and three of those five calls are
provably redundant.** See §6 #4.)

Two of these are worth stating out loud:

* **The transport solve is the run.** Everything else together is under 15 %.
* **Kill evaluation, mechanism attribution and payload encoding together are
  0.1 % of the run.** `attribute_mechanism` defaults to `True` and its
  docstring calls it "a handful of numpy passes" — measured, it is 3
  milliseconds against 47 seconds. There is no performance argument for
  turning any of the three off, and none should be made.

Note `solve_oxygen` is called **twice** and `solve_ifp` **three times** per cold
run. That is not a typo; see §5.

---

## 3. Inside the transport solve

`cProfile` cannot see inside `factorized(A)` — SuperLU's `solve` is a C method —
so the split below comes from the `instrument` arm, which wraps
`tmesim.physics.diffusion.factorized` and `_advective_divergence` with
`time.perf_counter()` calls and changes no arithmetic. The harness asserts
that: arm A (clean) and arm B (shimmed) both returned
`kill_fraction = 0.11221232991768859`, **`BIT-EXACT vs arm A: True`**.

`dx = 12`, warm cache, `simulate` stage = 98.138 s (loadavg not recorded on this
arm; its total, 103.479 s, is the fourth row of §0.2's table, so read the
percentages and not the seconds):

| | seconds | % of simulate |
|---|---|---|
| sparse LU factorisation (once) | 1.335 | 1.4 % |
| **2,160 triangular solves** | **68.443** | **69.7 %** |
| 2,160 upwind advection passes | 15.548 | 15.8 % |
| everything else in `simulate` | 12.813 | 13.1 % |

Operator `A = I - dt*L + dt*K`: **286,200 x 286,200, nnz = 1,041,638.**
Per-solve mean **31.69 ms**.

The same shape from the unshimmed `cProfile` run, by self time:

```
   ncalls  tottime  percall  filename:lineno(function)
     2161   59.005    0.027  {method 'solve' of 'SuperLU' objects}
     2160   10.951    0.005  physics/diffusion.py:152(_advective_divergence)
        5    8.261    1.652  {scipy.sparse.linalg._dsolve._superlu.gssv}     <- spsolve: oxygen + IFP
        1    6.348    6.348  physics/diffusion.py:186(simulate_drug_transport)
        2    2.238    1.119  {scipy.sparse.linalg._dsolve._superlu.gstrf}    <- the transport factorisation
```

So the run is, to first order, **2,160 back-substitutions against one sparse
LU**. `PDE_DT_S = 10.0` and `total_time_s = 21600` fix the step count at 2,160
regardless of `dx`; `dx` sets the cost of each step. That is why §2.2's
"6.9x time for 6.2x voxels" is nearly linear: the step *count* does not move.

---

## 4. Memory, part 1 — the floor before any grid exists

`experiments/profile_pipeline.py baseline` imports the production modules in
the order `execute_run` reaches them and weighs each one:

```
  step                                      peak RSS    rose by    secs
  ---------------------------------------------------------------------
  (interpreter + stdlib)                       16.7 M
  numpy                                        26.7 M      10.0    0.07
  scipy.sparse.linalg                          49.3 M      22.6    0.22
  pandas                                       89.8 M      40.5    0.74
  tmesim.core                                  90.1 M       0.2    0.01
  tmesim.physics.diffusion                     90.1 M       0.0    0.00
  tmesim.sim.engine                            94.2 M       4.1    0.05
  tmesim.reconstruct.real_pipeline            107.8 M      13.6    0.17
  tmesim.molecular.library                    133.9 M      26.0    0.28
  tmesim.validation.calibration               134.8 M       1.0    0.01
  tmesim.validation.delivery_metrics          191.1 M      56.2    1.17
  app.server.cache                            191.1 M       0.0    0.00
  app.server.jobs                             191.1 M       0.0    0.00
  app.server.api                              202.0 M      10.9    0.28
  api.list_sections()                         293.2 M      91.2    0.17
  api.default_spec()                          293.2 M       0.0    0.06
  calibrate_alpha()                           300.7 M       7.4    0.61
  ---------------------------------------------------------------------
  FLOOR before any RunSpec executes: 300.7 MiB
```

**300.7 MiB is gone before `dx` has been read.** On a 512 MB instance that is
59 % of the budget. Two entries deserve naming:

* **`list_sections()` costs 91.2 MiB of high-water in 0.17 s.** Enumerating
  sections one at a time, in isolated processes, shows it is not spread across
  the registry — it is **one section**:

  Each row below is `describe_section(id)` as the **first and only** call in a
  fresh interpreter, so every one is charged its true cost rather than
  inheriting the previous section's high-water mark:

  | section | +peak RSS | secs | n_cells |
  |---|---|---|---|
  | `gbmspace_at10_bra_5_fo_1_1` (ingested) | **+173.5 MiB** | 0.37 | **2,956** |
  | `visium_gbm1` (built-in) | **+71.2 MiB** | 0.32 | 57,050 |
  | `visium_10x_parent_gbm` (built-in) | **+71.3 MiB** | 0.27 | 41,560 |
  | `xenium_319`, `xenium_io`, `gse324268` | +1.5, +1.1, +1.1 MiB | 0.00 | 40,887 / 816,769 / 851,509 |

  The sections memoised in `app/runs/_section_meta.json` cost nothing — note
  `gse324268` reports **851,509 cells for 0.00 s and 1.1 MiB**. The ones that
  miss the memo read a matrix to report a cell count, and **a section with
  2,956 cells costs 173.5 MiB to describe**. Because RSS is a high-water mark,
  a *sequential* enumeration pays the largest single miss once, not the sum:
  all 38 sections together cost 174.7 MiB, essentially the gbmspace row alone.

  `GET /api/sections` and `default_spec()` both call `list_sections()`.

  **Scope.** This host has 38 registered sections (7 built-in + 31 ingested
  `gbmspace_*`). The hosted container has only the 7 built-ins, so it does not
  pay the 173.5 MiB. With the ingested registry suppressed to mirror the
  container, enumerating all 7 built-ins still costs a **72.5 MiB** high-water
  spike, essentially all of it `visium_10x_parent_gbm`. So this is a real cost
  on the hosted instance too — and it is a cost that **grows with every section
  an operator ingests**, on a code path that runs at startup and on every
  `/api/sections` poll.

* **`tmesim.validation.delivery_metrics` costs 56.2 MiB and 1.17 s to import**
  and pulls in the whole validation stack behind it.

---

## 5. Memory, part 2 — the curve, and where the cliff is

### 5.1 The peak is reached in `reconstruct`, not in the solve

This is the finding that reframes the whole `TMESIM_MAX_GRID_CELLS` question.
Peak RSS snapshotted at every stage exit, one cold subprocess per point:

| dx | voxels | pre-stage | **reconstruct** | compile_drug | simulate | metrics | payloads |
|---|---|---|---|---|---|---|---|
| 96 | 4,644 | 311.4 | 471.4 | 477.1 | 477.2 | 477.8 | 477.8 |
| 60 | 11,696 | 311.7 | 447.0 | 448.8 | 448.8 | 448.8 | 448.8 |
| 36 | 32,092 | 312.0 | 512.6 | 518.6 | 551.7 | 551.7 | 551.7 |
| 30 | 46,070 | 312.8 | 432.0 | 435.1 | 461.1 | 461.1 | 461.1 |
| 24 | 71,994 | 312.3 | 588.0 | 588.7 | 610.4 | 610.4 | 610.4 |
| 20 | 103,530 | 314.0 | 615.4 | 616.5 | 616.5 | 616.5 | 616.5 |
| 16 | 161,226 | 313.2 | 742.5 | 742.5 | 742.5 | 742.5 | 742.5 |
| 14 | 210,177 | 311.9 | 729.9 | 729.9 | **781.3** | 781.3 | 781.3 |
| **12** | **286,200** | 314.2 | **1011.6** | 1011.6 | 1011.6 | 1011.6 | 1011.6 |

At `dx = 12` the entire 697 MiB rise happens inside `reconstruct`, and nothing
afterwards — not the 2,160-step transport solve, not the 24-frame time-lapse,
not the 31 MB of payloads — moves the high-water mark again.

### 5.2 The single allocation that sets the cliff

Tracing the reconstruction at `dx = 12`, checkpoint by checkpoint:

```
checkpoint                                     peak RSS     rose
imports                                          113.8M
default_spec (list_sections)                     229.4M    115.6
  -> reconstruct_real ENTER                      229.4M      0.0
  -> spsolve #1 A=286200 nnz=1428802             934.0M    704.6     <-- HERE
  -> reconstruct_real EXIT                       934.0M      0.0
  -> _max_darcy_speed ENTER                      934.0M      0.0
  -> spsolve #2 A=286200 nnz=1038592             934.0M      0.0
  -> solve_ifp #1                                934.0M      0.0
build_reconstruction DONE (6.5s)                 934.0M      0.0
```

**One `scipy.sparse.linalg.spsolve` on the steady-state oxygen system —
286,200 x 286,200, 1,428,802 nonzeros — costs 704.6 MiB.** That is the memory
cliff. It is a sparse LU: the factors fill in far beyond the 5-point stencil,
and SuperLU builds them, uses them once, and frees them.

The second `spsolve` (the IFP system, 1,038,592 nnz) reuses the freed pages, so
it does not move the high-water mark — but it would allocate comparably if it
ran first.

Isolated, on the captured matrices:

| solver | time | +peak RSS | bit-exact vs `spsolve` default |
|---|---|---|---|
| `spsolve` (SuperLU, COLAMD) — **current** | 5.18 s | **+501.4 MiB** | — |
| `spsolve(permc_spec="MMD_AT_PLUS_A")` @ dx=20 | 0.97 s | +107.9 MiB (vs 130.6) | **No** — rel 1.28e-12 |
| `spsolve(permc_spec="NATURAL")` @ dx=20 | 29.77 s | +927.2 MiB | No — rel 3.43e-11 |
| Jacobi-preconditioned CG, rtol 1e-10 (2,493 its) | 18.35 s | **+65.9 MiB** | No — rel 3.87e-11 |
| Jacobi-preconditioned CG, rtol 1e-12 (2,797 its) | 20.04 s | +74.6 MiB | No — rel 8.06e-13 |

So the memory lever exists and is large (**501 -> 75 MiB, 6.7x**), and it costs
**3.9x more wall-clock on that solve** and is **not bit-exact**. Reordering is
not the answer: at `dx = 20` a better symmetric ordering buys only 17 % of the
memory, and it is not bit-exact either.

### 5.3 The curve, fitted

Peak RSS against voxels, 12 cold subprocesses, `xenium_319`,
loadavg 4.6–6.8 throughout:

| dx (um) | grid | voxels | peak RSS (MiB) | run delta (MiB) | wall (s) | cpu (s) |
|---|---|---|---|---|---|---|
| 96 | 54x86 | 4,644 | 477.8 | 248.9 | 2.51 | 2.64 |
| 72 | 72x114 | 8,208 | 529.9 | 301.2 | 2.68 | 2.91 |
| 60 | 86x136 | 11,696 | 448.8 | 219.9 | 3.74 | 3.81 |
| 48 | 107x170 | 18,190 | 520.3 | 289.3 | 3.85 | 4.11 |
| 40 | 128x204 | 26,112 | 505.9 | 276.9 | 4.80 | 4.95 |
| 36 | 142x226 | 32,092 | 551.7 | 322.6 | 5.61 | 5.58 |
| 30 | 170x271 | 46,070 | 461.1 | 231.8 | 9.85 | 9.62 |
| 24 | 213x338 | 71,994 | 610.4 | 381.4 | 14.98 | 14.52 |
| 20 | 255x406 | 103,530 | 616.5 | 385.8 | 24.80 | 23.96 |
| 16 | 318x507 | 161,226 | 742.5 | 512.7 | 32.26 | 32.00 |
| 14 | 363x579 | 210,177 | 781.3 | 552.8 | 55.02 | 52.38 |
| 12 | 424x675 | 286,200 | **1011.6** | 781.0 | 107.91 | 99.85 |

Least squares on the grid-dominated regime (voxels >= 46,070, n = 6):

```
peak RSS  =  407.1 MiB  +  2.074 KiB/voxel        R^2 = 0.954
run delta =  177.8 MiB  +  2.071 KiB/voxel        R^2 = 0.954
```

Below ~46,000 voxels the curve is **flat and noisy** (448–552 MiB across a 7x
voxel range, non-monotone). That is the floor of §4 plus a reconstruction
transient that does not depend on the mesh, plus allocator variance. **There is
no useful `dx` signal at the coarse end**: making the mesh coarser than
`dx = 30` buys essentially no memory on this host.

`rss_before_run` was 228.5–231.0 MiB on every one of the 12 points — the floor
is extremely reproducible; the run delta is what moves.

### 5.4 Where the cliff is

Two facts are already in `docs/DEPLOYMENT.md`, measured on the real 512 MB
Render instance and **not by me**:

* `dx = 30` (46,070 voxels) peaks at **220 MiB** — survives.
* `dx = 12` (286,200 voxels) is **OOM-killed**.

`TMESIM_MAX_GRID_CELLS=auto` currently adopts **46,070** — the larger of the two
known-good points — and the page says, correctly, that it "is deliberately
conservative about everything between them". This profile supplies the missing
middle: a **slope**.

Anchoring this host's measured slope (2.074 KiB/voxel) on that container's
measured point (46,070 voxels -> 220 MiB):

| budget | predicted ceiling | dx on `xenium_319` |
|---|---|---|
| 512 MiB (no headroom) | ~190,000 voxels | ~14.7 um |
| 480 MiB (~6 % headroom) | ~174,000 voxels | ~15.4 um |
| **450 MiB (~12 % headroom)** | **~160,000 voxels** | **~16 um** |

The model is internally consistent with both known container points: it puts
46,070 at 220 MiB by construction, and 286,200 at **~752 MiB**, which is above
512 and therefore OOM — which is what was observed.

**This is a prediction, not a measurement, and it is falsifiable with one run.**
The slope comes from macOS/arm64 and the anchor from Linux/x86; the slope is
dominated by SuperLU fill-in, which is a property of the matrix and the
ordering rather than of the OS, which is the reason to expect it to transfer —
but that is an argument, and the point of this file is not to make arguments.

> **The test.** Set `TMESIM_MAX_GRID_CELLS=200000` on the 512 MB instance and
> submit `dx = 16` (161,226 voxels). This model says it peaks near **460 MiB**
> and completes. If it OOMs, the slope does not transfer and `auto` should stay
> at 46,070. If it completes, `auto` is **3.5x more conservative than the
> hardware requires** and is refusing resolutions the instance can serve.

Either outcome replaces a guess with a measurement, which is the point.

**Scope, and it is narrow.** Everything in §5.3 and §5.4 is `xenium_319` only.
The slope is a slope *for that section*, and the next subsection shows that on
four of the seven built-in sections the grid is not the dominant term at all —
so **this model must not be used to set a global ceiling.** A voxel ceiling
tuned on `xenium_319` will pass a section that needs six times the memory at a
twenty-eighth of the voxel count.

### 5.5 Against section size — and this inverts everything above

Peak RSS at a **fixed `dx = 30 um`**, one cold subprocess per section,
loadavg 24–38 throughout (so read the memory, not the time):

| section | source | voxels | n_cells | **peak RSS (MiB)** | run delta (MiB) | wall (s) |
|---|---|---|---|---|---|---|
| `visium_gbm3` | Visium | **1,665** | 21,279 | **3004.1** | 2772.0 | 21.2 |
| `visium_gbm1` | Visium | 3,190 | 57,050 | **2475.2** | 2243.6 | 35.2 |
| `visium_gbm2` | Visium | 2,970 | 34,707 | **2242.3** | 2010.8 | 33.9 |
| `visium_10x_parent_gbm` | Visium | 2,850 | 41,560 | **2216.3** | 1985.1 | 30.8 |
| `xenium_319` | Xenium | 46,070 | 40,887 | **553.2** | 321.3 | 41.9 |
| `xenium_io` | Xenium | 139,029 | 816,769 | **2024.9** | 1790.8 | 145.8 |
| `gse324268` | Xenium | 169,087 | 851,509 | **2087.5** | 1855.2 | 177.1 |

Read the first row again. **`visium_gbm3` is 1,665 voxels — 28x SMALLER than
`xenium_319` at the same `dx` — and it costs 5.4x the peak memory.** The §5.3
fit (2.07 KiB/voxel) predicts 3.4 MiB for that grid. It cost 2,772 MiB, a
factor of **800x** out.

**Peak RSS is not primarily a function of the grid.** It is a function of the
SECTION. `TMESIM_MAX_GRID_CELLS` counts voxels, and voxels are not the
dominant term for four of the seven built-in sections.

Tracing `visium_gbm3` end to end:

```
  imports                                    111.6 M
  deconvolve.load_reference (7.1s)          2034.2 M   +1922.6     <-- HERE
  reconstruct_from_visium DONE              2830.1 M    +795.8
  ... spsolve N=1665 nnz=8161                2838.0 M      +0.0    (the grid is free)
```

`tmesim/data/real/scref_gbm.py::load_reference()`, measured on its own:

```
load_reference(): 6.04 s   peak RSS 96.0 -> 1708.5 MiB (+1612.5)
  expression DataFrame: 23,460 genes x 3,589 cells, dtype int64
  resident size (memory_usage deep): 643.8 MiB
  nonzero fraction: 0.1100        -> CSR float64 would be ~106.0 MiB
  max value 77,122 -> fits int32  -> dense int32 would be 321.2 MiB
```

It is `pd.read_csv(..., sep=r"\s+")` on a 20.5 MB gzip, producing a **dense
int64 genes-x-cells DataFrame** whose parsing transient is 2.5x its own
resident size. It is loaded **on every Visium reconstruction**, is never
cached, and everything downstream needs from it is a small per-cell-type
signature matrix over the ~few-thousand genes shared with the Visium panel.

The two large Xenium sections (`xenium_io`, `gse324268`, ~850k cells) cost
~1.8–1.9 GB for a different reason again — cell count, not voxel count.

**Consequence for the hosted instance.** `docs/DEPLOYMENT.md` records that
mounting `data_cache` read-only makes **all 7 built-in sections available**.
Four of those seven are Visium. On this host each of them needs 2.2–3.0 GB at
`dx = 30`, which is 4–6x a 512 MB instance, and the shortfall is a data-loading
cost (numpy/pandas array sizes) rather than anything OS-specific. **The
prediction is that all four Visium sections OOM-kill a 512 MB container at any
`dx`, and that no setting of `TMESIM_MAX_GRID_CELLS` prevents it**, because
1,665 voxels passes every ceiling anyone would configure. That prediction has
not been run in the container and should be, before it is believed — but it is
one `POST /api/runs` away from being settled either way, and if it holds, the
Studio is currently offering four sections that take the service down.

---

## 6. The ranked top five

Ordering is by **the size of the thing measured**, not by how easy it looks.
Three of the five are time and two are memory; the column says which, because
ranking a second against a megabyte is not a ranking.

| # | cost | kind | size | bit-exact fix possible? |
|---|---|---|---|---|
| 1 | 2,160 implicit triangular solves | time | 69.7 % of `simulate` | **No** |
| 2 | scRNA reference densified per Visium run | memory | **+1,922.6 MiB**, 6–7 s | **Yes** |
| 3 | oxygen sparse LU (the Xenium cliff) | memory | +704.6 MiB | **No** |
| 4 | elliptic solves recomputed bit-identically | time | ~4.6 s of 47 s (~10 %) | **Yes, verified** |
| 5 | 2,160 upwind advection passes | time | 15.8 % of `simulate` | Yes, conditionally |

Two smaller items that are also bit-exact and cheap are in §6.6.

---

### #1 — The 2,160 implicit triangular solves. 69.7 % of `simulate`, ~60 % of the run.

**What it is.** `physics/diffusion.py::simulate_drug_transport` marches the
transport PDE with an implicit diffusion step. The operator `A` is
time-invariant, so it is factorised once (1.335 s) and back-substituted once per
step. `PDE_DT_S = 10.0` over `total_time_s = 21600` fixes the count at 2,160;
`dx = 12` makes each back-substitution 31.69 ms on 286,200 unknowns.

**Why it costs that.** 2,160 x 31.69 ms = 68.4 s. There is no waste in it — the
factorisation is already hoisted, the operator is already the minimal 5-point
stencil (`nnz = 1,041,638` on 286,200 rows = 3.6/row), and `SuperLU.solve` is C.
It is simply the arithmetic the chosen discretisation asks for.

**What could be done.**

* **(a) A larger `dt`.** `results/numerical_rigor.md` (quoted verbatim in every
  run's caveats) says the prediction is "dt-converged to <1e-2 for
  dt <~ 20 s". `dt = 20` halves the step count: ~34 s off a 47 s run. The
  advection CFL is the binding constraint and must be re-checked
  (`max_speed*dt/dx <= 1`, and `dx` shrinks as the mesh refines).
* **(b) The steady-state field instead of the march.** This repo has already
  measured that a steady field is 48.6x cheaper to solve than the march for the
  reach map, and `app/server/reach.py` already uses it. If the *reported*
  quantities (AUC, coverage, penetration length) can be expressed on a steady
  field, the march disappears entirely. That is a large scientific question,
  not a performance one, and it is out of scope here.
* **(c) Nothing inside the loop.** Do not micro-optimise around
  `SuperLU.solve`; it is 69.7 % of `simulate` and it is already C.

**Risk.** (a) is a change to the discretisation: every stored calibration,
every regression baseline and every `app/runs/*/result.json` was produced at
`dt = 10`. (b) changes what is being computed.

**Could a fix be bit-exact?** **No.** Both (a) and (b) change the arithmetic.
Per the standing rule, either is a **model change**, is pre-registered and
measured as one, or is not done.

---

### #2 — The scRNA reference is densified on every Visium reconstruction. +1,922.6 MiB.

**The largest single allocation measured anywhere in this product** — 2.7x the
oxygen LU of #3, on sections whose grids are 1,665 to 3,190 voxels.

**What it is.** `tmesim/reconstruct/deconvolve.py::deconvolve_sample` calls
`tmesim/data/real/scref_gbm.py::load_reference()`, which does

```python
expr = pd.read_csv(counts_path, sep=r"\s+", index_col=0)
```

on `data_cache/scref_gbm/GSE84465_GBM_All_data.csv.gz` (20.5 MB on disk) and
returns a **dense int64 DataFrame, 23,460 genes x 3,589 cells**.

**Why it costs that.**

| | |
|---|---|
| file on disk | 20.5 MB gzip |
| resident DataFrame (`memory_usage(deep=True)`) | **643.8 MiB** |
| peak RSS during the load | **+1,612.5 MiB** (measured alone) / **+1,922.6 MiB** (in the pipeline) |
| time | **6.04 s** alone, 7.1 s in the pipeline |
| nonzero fraction of the matrix | **0.1100** |
| same data as CSR float64 | ~106.0 MiB (**6.1x smaller**) |
| max count 77,122, so int32 suffices | 321.2 MiB dense (**2x smaller**) |

Three separate multipliers stack: the matrix is stored **dense** when it is
89 % zeros; it is stored **int64** when the largest count is 77,122; and
`read_csv` on a whitespace-delimited gzip peaks at **2.5x** the size of the
object it finally produces.

And it is **not cached**. Every Visium reconstruction — every cache miss, for
each of the four built-in Visium sections — pays the full 6 s and 1.9 GB
again. What the pipeline actually consumes from it is
`build_signatures(shared_genes, reference=...)`: a per-cell-type mean over the
few thousand genes shared with the Visium panel, i.e. a `(G_shared, K)` matrix
with `K` around seven. **A ~1 MB result is being produced from a 1.9 GB
transient, from first principles, on every run.**

**What could be done.** In increasing order of surgery:

**(D2-a) Memoise the signatures on disk, and stop loading the reference in the
server at all.** The signature matrix is a deterministic function of the
reference file and the requested gene list (`max_cells=None`, so the
`random_state` subsample never runs).

```diff
--- a/tmesim/reconstruct/deconvolve.py
@@ deconvolve_sample
-    expr_df, labels_df, _ = load_reference()
-    ref_syms = set(expr_df.index[~expr_df.index.duplicated()])
+    # The reference is a 23,460 x 3,589 dense int64 frame: 643.8 MiB resident
+    # and a +1,612.5 MiB parse transient, for a (G_shared, K) signature matrix
+    # of a few hundred kB (results/performance_profile_2026-09-20.md §6 #2).
+    # `ref_symbols()` reads only the gene index; `cached_signatures()` returns
+    # the memoised means and falls back to load_reference() on a miss, so the
+    # first run after a reference change still pays and every later one does
+    # not.
+    ref_syms = ref_symbols()
@@
-    signatures, type_names = build_signatures(
-        shared_genes, reference=(expr_df, labels_df, None)
-    )
+    signatures, type_names = cached_signatures(shared_genes)
```

*Bit-exact:* **yes, by construction, and it must be proven that way** — the
cache-populating run and a `load_reference()` run must be asserted
`np.array_equal` on `signatures`, on a real sample, with the cache key
including a fingerprint of the two reference files so a changed reference
invalidates rather than serves a stale signature. Without that fingerprint this
becomes a correctness bug rather than an optimisation.

**(D2-b) If the reference must be loaded, load it sparse and narrow.** Read the
counts with `dtype=np.int32` and convert to CSR immediately: 643.8 -> ~106 MiB
resident, and a smaller parse transient. *Bit-exact:* **int32 yes** (every value
is representable exactly — max 77,122 — and `astype` of an exact integer is
lossless), **CSR yes** for the sums and means `build_signatures` takes.
The parse itself is the harder part: `read_csv(dtype=...)` still materialises
chunks, so the transient does not fall proportionally and needs its own
measurement rather than an assumption.

**(D2-c) `functools.lru_cache` on `load_reference()`.** One line, and **the
wrong fix on the machine that matters**: it removes the repeat *time* but pins
643.8 MiB resident for the life of the process. On a 512 MB instance that is
worse than the status quo. Listed so it is not reached for by reflex.

**Risk.** Low for (a) with the fingerprint, moderate without it. Medium for (b)
— the dtype claim is exact but the parse-transient claim is not yet measured.

**Could a fix be bit-exact?** **Yes** — this is the largest item on this list
that can be, and the only one of the two memory items that can be. It is worth
doing before anything in #1 or #3 is contemplated.

**And it may be urgent rather than merely large.** See §5.5: four of the seven
built-in sections are Visium, and on this host each needs 2.2–3.0 GB at
`dx = 30`. If that transfers to the 512 MB instance — untested, and one run
settles it — the Studio is currently offering four sections that OOM-kill the
container, at a voxel count no ceiling would ever refuse.

---

### #3 — The oxygen sparse LU: +704.6 MiB, and it is the OOM cliff.

**What it is.** §5.2. One `spsolve` on a 286,200-row system inside
`reconstruct_real`. It is 100 % of the memory rise at `dx = 12` and it is the
reason `TMESIM_MAX_GRID_CELLS` exists.

**Why it costs that.** A sparse direct LU of a 2D 5-point system fills in far
beyond the original 1.43M nonzeros; SuperLU builds the factors, uses them for
one right-hand side, and frees them. The right-hand side count is **one**,
which is the whole problem: a direct factorisation amortises over many solves
and here there are none to amortise over.

**What could be done.** Replace the direct solve with a
Jacobi-preconditioned CG (the matrix is SPD: `-D*L` positive definite plus a
non-negative diagonal). Measured on the real `dx = 12` oxygen system:

```
spsolve (current)            5.18 s   +501.4 MiB
cg(jacobi, rtol=1e-12)      20.04 s    +74.6 MiB    2,797 iterations
                                                    max|d| vs direct = 3.40e-13
```

**6.7x less memory, 3.9x more time on that solve.** End-to-end that is roughly
**+15 s on a 47 s run (+32 %)** in exchange for removing ~430 MiB from the peak.

**Risk.** High, and in an unusual direction: it makes the product *slower* in
order to make it *possible* on a small host. It also introduces a solver
tolerance where there was none — a knob that changes the answer and will need
its own convergence study.

**Could a fix be bit-exact?** **No.** Measured `rel = 8.06e-13` at
`rtol = 1e-12`, and 1.28e-12 even for a mere reordering
(`permc_spec="MMD_AT_PLUS_A"`). **There is no memory fix here that is bit-exact.**
Any of them is a model change and must be pre-registered as one — including the
one-word `permc_spec` change, which looks like a pure optimisation and is not.

A genuinely bit-exact partial: the IFP system is solved 2–3 times (see #4) and
the oxygen system twice. Removing the duplicates does **not** lower the peak,
because the first solve already sets the high-water mark — this is worth
stating because it is the natural wrong conclusion to draw from #2.

---

### #4 — Elliptic solves recomputed bit-identically. 1 redundant oxygen + 2 redundant IFP per cold run.

**This is the only item on this list that can be bit-exact, and it was verified
to be, not assumed.**

**What it is.** One cold `dx = 12` run performs `solve_oxygen` **twice** and
`solve_ifp` **three times**, on inputs that do not change between the calls.
Instrumented (every namespace that binds the two names was patched, so no call
site is invisible):

```
solve_oxygen calls: 2   solve_ifp calls: 3
  ox:  #1  6.00s   #2  5.07s   bit-identical to #1: True   max|d|=0.000e+00
  ifp: #1  2.53s   #2  2.07s   bit-identical to #1: True   max|d|=0.000e+00
                   #3  1.83s   bit-identical to #1: True   max|d|=0.000e+00
```

The call sites, with stacks:

| # | site | what it is for |
|---|---|---|
| O1 | `real_fields.add_oxygen` <- `real_pipeline._install_oxygen` | the reconstruction's `fields.oxygen` |
| O2 | `sim/engine.py:306` in `simulate_response` | re-solved and overwrites `work.fields.oxygen` |
| I1 | `app/server/cache.py:1214` `_max_darcy_speed` | solves IFP purely to report `max_darcy_speed_um_s` for the CFL banner; discards both outputs |
| I2 | `sim/engine.py:314` in `simulate_response` | the velocity the transport step actually uses |
| I3 | `app/server/reach.py:945` in `compute_reach` | re-solves because `compute_reach` was not handed a velocity |

`add_oxygen` defaults to `OxygenConfig()` and `SimConfig.oxygen` defaults to
`OxygenConfig()`; the server never overrides either, which is why O1 and O2 are
byte-identical. The three IFP solves are byte-identical for the same reason.

**What it costs.** Two independent measurements, because neither alone is
clean:

* *Directly, under load 26.67* (the instrumented run printed above): the five
  elliptic solves cost **17.50 s** in total, of which the three duplicates
  (O2 5.07 s + I2 2.07 s + I3 1.83 s) are **8.97 s**, on a 174.27 s run.
* *From the clean `cProfile` run*: all five direct solves appear as
  `_superlu.gssv` self time, **8.261 s over 5 calls** — C time, which tracing
  barely inflates. Apportioning by the per-call ratio measured above
  (oxygen ~2.6x an IFP solve) gives ~2.63 s per oxygen solve and ~1.00 s per
  IFP solve, so the three duplicates are **~4.6 s — about 10 % of a 47 s cold
  run.**

On a *warm*-cache run, O1 and I1 do not execute, so the remaining duplicate is
I3 alone — and it is why `metrics` is 3.1 % of the run (`compute_reach` cum
2.676 s in the profile, of which the IFP solve is the bulk).

**What could be done** — three independent diffs, smallest first.

**(D4a) Hand `compute_reach` the velocity the engine already solved.** Removes
I3. Additive, and the arrays were measured byte-identical.

```diff
--- a/tmesim/sim/engine.py
@@ simulate_response, the SimulationResult construction
         fields={
             "oxygen": oxygen,
             "ifp": ifp,
+            # The Darcy velocity solve_ifp returned alongside `ifp`. Carried so
+            # a downstream consumer (app/server/reach.compute_reach) does not
+            # have to re-solve the same system to get it: measured
+            # bit-identical, max|d| = 0.0 (results/performance_profile_2026-09-20.md §6 #4).
+            "velocity": velocity,
             "auc": transport.auc,
             "vascular_density": work.fields.vascular_density,
         },

--- a/app/server/reach.py
@@ compute_reach, the convection branch
         if not convection:
             velocity = None
         elif velocity is None:
-            from tmesim.physics.pressure import solve_ifp
-            _ifp, velocity = solve_ifp(micro)
+            velocity = (sim_result.fields or {}).get("velocity")
+        if convection and velocity is None:
+            from tmesim.physics.pressure import solve_ifp
+            _ifp, velocity = solve_ifp(micro)
```

*Bit-exact:* **yes**, and measured (`max|d| = 0.000e+00` across all three
solves). The engine solves on `_working_copy(micro)` and `compute_reach` solves
on `micro`; the copy preserves every field the IFP solve reads, which is why
the outputs are identical rather than merely close. A guard belongs with this
diff, asserting `np.array_equal` between the carried velocity and a fresh solve
on one real section, so the claim is tested rather than reasoned.

**(D4b) Stop `_max_darcy_speed` re-solving what the engine is about to solve.**
Removes I1. `_max_darcy_speed` exists to report the advection CFL. The cheapest
honest version keeps the number and stops paying for it on every build:

```diff
--- a/app/server/cache.py
@@ build_reconstruction
-    built["max_darcy_speed_um_s"] = _max_darcy_speed(micro)
+    # Lazily: the CFL banner is read by /api/validate and by the run report,
+    # neither of which is on the hot path. A full solve_ifp on the flagship
+    # grid is 2.5 s and ~500 MiB of transient LU
+    # (results/performance_profile_2026-09-20.md §5.2), paid on every cache
+    # miss whether or not anybody asks for the number.
+    built["max_darcy_speed_um_s"] = None
+    built["_max_darcy_speed_fn"] = lambda: _max_darcy_speed(micro)
```

*Bit-exact:* **yes for the run** — the run never reads this number. But it is a
**contract change** for `report["max_darcy_speed_um_s"]` (always-present
becomes computed-on-demand), and `/api/validate`'s CFL check would pay the 2.5 s
instead. Prefer D4c below if both land.

**(D4c) Let the engine reuse the reconstruction's oxygen, opt-in.** Removes O2,
the largest single duplicate (5.07 s measured, ~2.6 s at low load).

```diff
--- a/tmesim/core/config.py
@@ class SimConfig
     oxygen: OxygenConfig = field(default_factory=OxygenConfig)
+    #: Reuse ``micro.fields.oxygen`` instead of re-solving it, when the caller
+    #: KNOWS the field was produced by solve_oxygen(micro, this same cfg.oxygen).
+    #: Default False keeps every stored result bit-identical; the server sets
+    #: it only for oxygen_mode == 'pde', which is exactly the case in which the
+    #: two solves were measured byte-identical (max|d| = 0.0).
+    reuse_field_oxygen: bool = False

--- a/tmesim/sim/engine.py
@@ simulate_response step 2
-    oxygen = solve_oxygen(work, cfg.oxygen)
+    if getattr(cfg, "reuse_field_oxygen", False) and work.fields.oxygen is not None:
+        oxygen = np.asarray(work.fields.oxygen, dtype=np.float64)
+    else:
+        oxygen = solve_oxygen(work, cfg.oxygen)

--- a/app/server/jobs.py
@@ execute_run, the SimConfig construction
-    cfg = SimConfig(drug=dp, pk=pk_cfg, pde=pde, seed=SIM_SEED, **cfg_kwargs)
+    # oxygen_mode == 'pde' means the reconstruction's fields.oxygen came from
+    # solve_oxygen(micro, OxygenConfig()) -- the identical call simulate_response
+    # would make. Measured byte-identical on the flagship section
+    # (results/performance_profile_2026-09-20.md §6 #4). Any other mode leaves
+    # the engine re-solving, which is what the existing caveat describes.
+    cfg = SimConfig(drug=dp, pk=pk_cfg, pde=pde, seed=SIM_SEED,
+                    reuse_field_oxygen=(str(oxy_mode) == "pde"),
+                    **cfg_kwargs)
```

*Bit-exact:* **yes, conditionally, and the condition is the whole thing.** It
holds only while `add_oxygen`'s default config and `SimConfig.oxygen`'s default
stay equal. That equality must be asserted in a test, not left as a comment —
otherwise this is exactly the class of change that is bit-exact on the day it
lands and silently a model change six weeks later. It also interacts with the
standing caveat in `jobs.py` that says the engine "RE-SOLVES the steady-state
oxygen PDE internally and overwrites fields.oxygen"; that caveat becomes wrong
for `oxygen_mode == 'pde'` and must be rewritten in the same commit.

**Risk.** Low for D4a, medium for D4c (the invariant above), and D4b changes a
reported field's availability. **Combined expected saving: ~4.6 s off a 47 s
cold run (~10 %) and ~1–2 s off a warm run, with no number changing.** No
memory saving: the peak is a high-water mark and the *first* solve already sets
it (§5.1).

---

### #5 — 2,160 upwind advection passes. 15.8 % of `simulate`.

**What it is.** `_advective_divergence` is called once per step and is pure
Python + numpy: 10.951 s of self time over 2,160 calls at `dx = 12` (5.07 ms
each) in the clean profile; 15.548 s in the instrumented arm.

**Why it costs that.** Each call allocates a fresh set of full-grid
temporaries — measured **49.3 bytes per voxel of peak traced allocation** at
64x64 and 48.7 at 128x128 (i.e. ~6 float64 grids per call). At 286,200 voxels
that is ~13 MiB allocated and freed **2,160 times**, plus the arithmetic.

**What could be done.** Hoist the per-call temporaries out of the loop and
write into them with `out=`:

```diff
--- a/tmesim/physics/diffusion.py
@@ simulate_drug_transport, before the step loop
+    # Scratch for the advection pass. _advective_divergence allocates ~6
+    # grid-sized temporaries per call and is called n_steps times (2,160 on the
+    # shipped window): 49.3 bytes/voxel x 286,200 voxels x 2,160 calls of
+    # allocate-and-free (results/performance_profile_2026-09-20.md §6 #4).
+    adv_scratch = _AdvectionScratch(ny, nx) if use_convection else None
@@ inside the loop
-            rhs -= dt * _advective_divergence(c, vx, vy, tissue, dx)
+            rhs -= dt * _advective_divergence(c, vx, vy, tissue, dx,
+                                              scratch=adv_scratch)
```

with `_advective_divergence` gaining an optional `scratch` argument that
supplies pre-allocated `div`, `vf`, `c_up` and `flux` buffers and uses
`np.multiply(..., out=...)` / `np.where(..., out=...)` in place of the current
expression temporaries.

**Risk.** Moderate, and specific: `np.where(vf > 0.0, a, b)` into a preallocated
buffer is the same operation in the same order; but the moment anyone "tidies"
`flux = np.where(active, flux, 0.0) / dx` into `flux /= dx; flux[~active] = 0`
the **order of the divide and the mask changes**, and on a field with any
non-finite or denormal values that is a different answer.

**Could a fix be bit-exact?** **Yes, if and only if the operation order is
preserved exactly** — same ufuncs, same sequence, only the destination changes.
That is provable and must be proven: run the unoptimised and optimised paths on
a real section and assert `np.array_equal` on `auc` and
`final_concentration`, not `allclose`. If the proof fails, it is a model change
and stops. Realistic saving: a fraction of 15.8 % of `simulate` — call it
3–6 s of a 47 s run — and it should not be attempted before #4, which is
larger, safer and already verified.

---

## 6.6 Also measured, smaller, and both bit-exact

### 6.6a — The 300.7 MiB floor, of which `list_sections()` is a 72–175 MiB transient.

**What it is.** §4. Before `dx` is read: 300.7 MiB on this host (38 sections),
~202 MiB before `list_sections()` runs. On a container-shaped registry
(7 built-ins) the enumeration still costs a **72.5 MiB** high-water spike, all
of it one section whose metadata memo misses.

**Why it costs that.** `describe_section` reads a real matrix to report a cell
count when `app/runs/_section_meta.json` has no valid fingerprint for it. A
section with **2,956 cells** costs **173.5 MiB and 0.37 s** to describe;
memoised sections with **851,509 cells** cost **0.0 MiB and 0.00 s**. The memo
works; the misses are what cost.

**What could be done.**

* Ship the memo, or populate it at image-build time, so no first request pays a
  file read. One prebuilt `_section_meta.json` covering the 7 built-ins removes
  the spike from the container entirely.
* Make `describe_section` read only the metadata it reports (n_cells, extent,
  panel) rather than materialising a matrix — a `describe` path that streams a
  header instead of a body.
* Import `tmesim.validation.delivery_metrics` (56.2 MiB, 1.17 s) lazily inside
  `execute_run` rather than pulling the validation stack into the server's
  resident set. It already *is* a function-local import in `jobs.py`; the cost
  above is the import itself, paid on the first run and kept forever.

**Risk.** Low. None of this touches arithmetic.

**Could a fix be bit-exact?** **Yes, trivially** — it is enumeration and import
order, and no number in any result depends on it. This is the only item on the
list with a good risk/benefit ratio *and* a bit-exact guarantee, which is why
it is worth doing even though it is fifth by size.

---

### 6.6b — `calibrate_alpha()` is recomputed on every run.

0.655 s per run at low load (cum, from the profile), 3.165 s under load, and
7.4 MiB. It runs `molecular.brownian.simulate_effective_diffusivity` ten times.
It is a pure function of frozen constants — its value is
`22.223477711172563`, which this repo already treats as an invariant and
asserts on. `functools.lru_cache` on it is bit-exact **by construction**, saves
~1.4 % of a `dx = 12` run and ~10 % of a `dx = 30` run (where `compile_drug` is
10.6 % of the total), and carries essentially no risk. It is listed here rather
than in the top five only because 0.65 s is small next to 68 s.

---

## 7. `tests/test_performance_invariants.py`

16 tests, **no wall-clock assertion anywhere**, because §0.2 shows a wall-clock
threshold on this repo's machines would be a measurement of the neighbours.
Three kinds of assertion, all exact and load-independent:

1. **Counts** — one sparse factorisation per transport solve regardless of step
   count; exactly `n_steps` triangular solves; exactly `n_steps` advection
   passes with convection on and **zero** with it off; exactly one
   vessel-distance transform per `execute_run`.
2. **Dtypes** — `auc`, `final_concentration`, `times`,
   `mean_concentration_history`, `d_eff`, the operator's data and the
   vessel-distance map must all be `float64`. A float32 downcast is the single
   most tempting "memory optimisation" available on the transport path and it
   is a model change, not an optimisation.
3. **Allocation scaling** — peak traced bytes at `N` versus `4N`, via numpy's
   `tracemalloc` domain (numpy array data is `malloc`'d and is invisible to
   `tracemalloc` without it). Linear gives ~4x; quadratic gives ~16x; the
   bounds sit at 8.0, which no constant-factor change reaches and no class
   change survives. Plus the exact five-point identity
   `nnz == 5*ny*nx - 2*(ny+nx)`, and "peak memory does not grow with the
   window".

**Every one of them was planted and watched go red.** Thirteen defects were
introduced into a throwaway copy of the tree (never the repo), the relevant
test was run, and the copy was restored and re-verified green:

| # | planted defect | test that caught it | how it failed |
|---|---|---|---|
| P1 | `factorized(a)` moved inside the step loop | `test_exactly_one_factorisation_per_transport_solve` | 50/200 factorisations, expected 1 |
| P2 | a second `solve(rhs)` per step | same | solve count doubled |
| P3 | advection computed then scaled by zero when convection is off | `test_advection_cost_is_one_pass_per_step_and_zero_when_off` | `50 == 0` failed |
| P4a | `auc`/`c` allocated as float32 | `test_transport_outputs_are_float64` | dtype |
| P4b | `_effective_diffusivity` returns float32 | `test_operator_and_diffusivity_are_float64` | dtype |
| P5 | every step's field appended to a list | `test_transport_peak_memory_does_not_grow_with_the_window` | ratio > 1.5 |
| P6 | 24 time-lapse frames allocated when none requested | `test_no_timelapse_allocates_no_frames` | `frames is not None` |
| P7 | one explicit structural zero stored per row | `test_diffusion_operator_nnz_is_the_five_point_stencil_exactly` | nnz != identity |
| P8 | `cKDTree` replaced by the one-line pairwise broadcast | `test_vessel_distance_is_not_quadratic_in_cell_count` | **ratio 3.988 < 2.0 failed** |
| P9 | dense NxN adjacency scratch in the operator assembly | `test_diffusion_operator_allocation_is_linear_in_grid_size` | **ratio 15.255 < 8.0 failed** |
| P11 | dense NxN scratch in the advection pass | `test_advective_divergence_allocation_is_linear_in_grid_size` | **ratio 15.858 < 8.0 failed** |
| P12 | dense NxN scratch in the distance map | `test_vessel_distance_is_linear_in_grid_size` | **ratio 15.851 < 8.0 failed** |
| P13 | the call counter itself perturbs the solve by 1e-15 | `test_counting_does_not_change_the_answer` | `array_equal` failed |
| P10 | `distance_um=dist` dropped from the `compute_all` call in `jobs.py` | `test_one_run_factorises_once_and_maps_vessel_distance_once` | `2 == 1` failed |

P8's ratio of **3.988** is worth pausing on: that is the predicted 4x for an
O(N*M) broadcast when M is quadrupled, arriving within 0.3 % of theory. The
allocation-scaling tests are not heuristics.

One test went red on its first run against *unmodified* code and was corrected
rather than loosened: `test_saved_history_is_o_steps_over_save_every_not_o_steps`
was written with the obvious formula `n_steps // save_every + 1`, which is
wrong whenever the window does not divide evenly — and the shipped 2,160/500
does not. The real count is in the test.

```
$ .venv/bin/python -m pytest tests/test_performance_invariants.py -q
................                                                         [100%]
16 passed in 11.31s
```

---

## 8. Commands

```bash
# the harness
.venv/bin/python experiments/profile_pipeline.py stages     --dx 12 --repeats 2
.venv/bin/python experiments/profile_pipeline.py stages     --dx 30 --repeats 1
.venv/bin/python experiments/profile_pipeline.py profile    --dx 12 --top 22
.venv/bin/python experiments/profile_pipeline.py profile    --dx 30 --top 12
.venv/bin/python experiments/profile_pipeline.py instrument --dx 12
.venv/bin/python experiments/profile_pipeline.py baseline
.venv/bin/python experiments/profile_pipeline.py sweep \
    --dx 96,72,60,48,40,36,30,24,20,16,14,12 --section xenium_319
.venv/bin/python experiments/profile_pipeline.py sweep --dx 30 \
    --section visium_gbm3,visium_gbm2,visium_gbm1,visium_10x_parent_gbm,\
xenium_319,xenium_io,gse324268

# the tests
.venv/bin/python -m pytest tests/test_performance_invariants.py -q
```

`python` does not exist in this environment; `.venv/bin/python` is the
interpreter. The `sweep` arm launches one fresh subprocess per point on purpose
— a peak-RSS figure from a process that has already held another grid is not a
peak-RSS figure for the grid you asked about.

---

## 9. `calibrate_alpha()`

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; \
  v=calibrate_alpha(); print(repr(v)); assert v == 22.223477711172563"
22.223477711172563          # start of this work
22.223477711172563          # end of this work
```

---

## 10. What this file does not claim

* It does not claim any of the five candidates is worth doing. It measures
  them.
* It does not claim a runtime for this product. §0.2 measured a 2.3x spread for
  one spec on one host in one hour.
* The `dx = 16` memory prediction in §5.4 is a **prediction**, stated with the
  single run that falsifies it.
* The claim in §5.5 and §6 #2 that the four Visium sections would OOM-kill a
  512 MB container is a **prediction from a macOS measurement**, not a
  container measurement. The measured part — 2,216 to 3,004 MiB peak on this
  host, and +1,922.6 MiB inside one named function — is solid; the transfer to
  the instance is not, and is settled by one `POST /api/runs`.
* The wall-clock column of §5.5 was taken at loadavg 24–38 and is quoted only
  so the table is complete. Do not compare those seven numbers to each other.
* Nothing here touches the scientific claim surface. No number in
  `results/` other than this file changes, no verdict moves, and the
  `kill_fraction` for the shipped default spec was bit-identical
  (`0.11221232991768859`) across every arm of every measurement in this
  document — including the instrumented one, which is how the instrumentation
  is known to be instrumentation.
