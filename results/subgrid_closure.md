# Can a sub-grid closure represent the perivascular boundary layer the mesh cannot?

> **One-line answer.** On the controlled problem, **yes, and it needs no machine
> learning at all**: a closed-form, **zero-fitted-parameter** source reweighting
> removes **84–97 %** of the grid error out of sample, cuts the grid span of the
> mean on **288 of 289** ladders, and removes **76 %** of the banded `F_state`
> score's excess over an **irreducible floor** — the floor being the fact that the
> *exact* solution's banded fraction is itself a function of `dx`, so no closure of
> any kind can go below it. **A learned closure was therefore never justified and
> none was fitted.**
>
> **On real tissue it does not transfer.** It reduces the `F_state` grid span on
> **one of three** oxygen configurations and increases it on the other two, while
> changing the answer itself by **62 %**. The reason is
> measurable rather than mysterious: the closure is derived from a 1-D slab cell
> problem, and on the GBM section **38.7 % of partial-vessel cells sit at a
> correction strength larger than anything the controlled problem ever produced**
> (`beta` up to **144** against a validated maximum of **5.68**), where the
> formula degenerates to the purely geometric `beta = 1/v` and stops being a
> boundary-layer closure at all.
>
> **And a third finding that changes the framing of the question.** The premise
> that the mesh limit is what makes the *shipped* ladder move is **false**. Under
> the **shipped** oxygen constants the controlled problem's worst discretisation
> error at `dx = 30 um` is **0.170 %** (§2). The shipped configuration has no
> boundary-layer resolution problem to close. `results/grid_convergence_gbm.md`
> §5's `F_state` swing of **0.0461** and `results/tail90_falsification.md` §6.1's
> **3.72 %** fine-rung `F_state` span are therefore **not** perivascular
> under-resolution and **no closure of the oxygen solve can touch them.**

**Status: a feasibility study. Nothing shipped was changed.** `OxygenConfig`
defaults (`supply = 0.05`, `consumption = 0.02`, `hypoxia_threshold = 0.2`),
`physics/oxygen.py`, and `real_fields._DEFAULT_VASCULAR_MODE` /
`vessel_network.discrete_vascular_field(geometry='binary')` are all untouched.
The closure lives in a **new, default-off module that nothing else in `tmesim/`
imports** (`tmesim/physics/subgrid.py`); `solve_oxygen_subgrid(..., enabled=False)`
is a bit-identical passthrough. `pytest` result is in §7.

* Generator: `experiments/subgrid_closure.py`
* Raw: `results/sweep/subgrid_identity.json`, `subgrid_sweep.json`,
  `subgrid_analytic.json`, `subgrid_transfer_partial.json` (the `dx = 4 um` rung
  exhausted memory, so the per-rung checkpoint is the complete record — see §5),
  `subgrid_clip.json`
* Controlled problem: `experiments/boundary_layer_resolution.py`'s `slab_exact`
  and `slab_exact_reference`, **imported unmodified**
* Section: `data_cache/xenium_gbm`, 40,887 cells / 23,812 tumour

---

## 1. Formulating it before touching a model

The brief offers three candidate closures — a correction to the coarse-cell
source term, an effective source strength, or a corrected cell-average given the
coarse solution. **The choice is decidable by measurement rather than by taste,
and it was.**

### 1.1 A conservation identity says where the error can possibly live

The shipped steady oxygen solve (`physics/oxygen.solve_oxygen`,
`sink_form="linear"`) is

```
(-D * L + diag(supply*v + q)) @ O = supply*v            q = consumption * rho_tum
```

`L` is `_build_laplacian`, a 5-point operator with **zero row sums**. Summing the
system over cells therefore gives, exactly and at every `dx`,

```
sum_j (supply*v_j + q_j) * O_j  =  sum_j supply*v_j
```

and the continuum obeys the same identity with integrals. Measured, the residual
of that identity on the coarse solve is **≤ 3.9e-16** at every point of §2 —
i.e. it is an identity, not an approximation.

Subtracting the discrete identity from the continuum one, the right-hand sides
cancel (`supply*v_j` is already the exact cell integral of the source under an
area-weighted painting), and the **entire** domain-mean error is the mismatch
between

| | |
|---|---|
| discrete | `supply * v_j * O_j` — the source weighted by the **cell-average** O |
| continuum | `integral of supply * v * O` — the source weighted by O **where the vessel actually is** |

When the vessel is sub-cell, `O` inside it is far higher than the cell average,
so the discrete form under-counts the source's own self-shading and the solve
**over-supplies**. **That is a source-weighting error, and it is one scalar per
cell.**

### 1.2 The oracle decomposition — the choice, measured

Two oracles, each handed the exact answer and asked how much of the error it
could remove *if it were perfect* (`subgrid_identity.json`):

* **`beta`** — multiply the source's contribution to the reaction **diagonal
  only** by the exact `mean_vessel(O) / mean_cell(O)`. RHS untouched.
* **`flux`** — replace every face diffusivity by the value that makes the
  discrete face flux equal the exact flux, given the exact cell averages. This
  is the transmissibility/effective-diffusivity family, which is what most
  multiscale schemes reach for first.

Jointly calibrated low-CMRO₂ constants (`supply = 194`, `consumption = 10`,
`rho = 1`; `lambda_source = 3.21 um`, `lambda_sink = 14.14 um`), source
half-width 7.5 um = `VESSEL_RADIUS_UM`, relative error of the domain mean:

| dx (um) | cells per `lambda_sink` | **base** | oracle `beta` only | oracle `flux` only | oracle **both** | exact `beta` |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 7.07 | +0.201 % | −0.193 % | +0.397 % | −9.9e-10 % | 1.0000 |
| 4 | 3.54 | −0.923 % | −1.534 % | +0.629 % | −2.9e-09 % | 1.0000 |
| 6 | 2.36 | +4.702 % | +0.032 % | +4.640 % | +4.3e-10 % | 1.0000 |
| 8 | 1.77 | −5.713 % | −6.406 % | +0.734 % | −5.6e-09 % | 1.0084 |
| **12** | **1.18** | **+6.007 %** | **−2.545 %** | **+8.742 %** | **+2.0e-09 %** | **1.1006** |
| 16 | 0.88 | +18.281 % | −0.313 % | +18.643 % | +2.8e-09 % | 1.2176 |
| 24 | 0.59 | +43.677 % | +2.209 % | +40.745 % | −4.3e-09 % | 1.4909 |
| 30 | 0.47 | **+62.828 %** | **+3.159 %** | +58.240 % | −4.7e-09 % | 1.7200 |

Three things are settled by this table and nothing below re-opens them.

1. **`beta` and `flux` together are EXACT.** Worst `|rel err of mean|` over all
   three configurations is **8.7e-11** and worst `L_inf` on the field is
   **6.1e-11** — the sparse solver's own tolerance. So the two terms are a
   **complete** parameterisation of the discretisation error, not a heuristic
   pair. Anything a learned closure could represent, these two already span.
2. **`beta` alone carries the error at coarse `dx`, and `flux` alone does not.**
   At `dx = 30 um`, `beta` takes +62.8 % to +3.2 % (**95 % of the error**) while
   `flux` takes it only to +58.2 % (**7 %**). **So the closure is a scalar per
   cell multiplying the source diagonal.** The effective-diffusivity family — the
   obvious first guess, and the one an exponentially-fitted scheme would give —
   is the *wrong* lever here.
3. **`flux` is the residual at fine `dx`.** At `dx = 4` and `8 um` the base error
   is dominated by `flux` (−5.7 % → +0.7 % under the flux oracle) and `beta` does
   nothing. That is a source-**geometry quantisation** term, the same one
   `results/boundary_layer_resolution.md` §4.3 documents, and it is left in
   place. §3 shows it becoming the dominant residual exactly there.

### 1.3 The closure, and what it is a function of

```
beta_j  =  mean of O over the vessel support in cell j  /  mean of O over cell j
```

evaluated in **closed form** from a local cell problem: a 1-D two-region cell of
width `dx` carrying a source of half-width `a = v*dx`, symmetric (no-flux) at the
vessel centre, and matched at its outer edge to a **semi-infinite decaying sink
medium** (`model="open"`, the default):

```
[0, a) :  D O'' = (supply + q) O - supply        lam_i = sqrt(D/(supply+q))
[a, oo):  D O'' = q O                            lam_o = sqrt(D/q)
```

Both cell means are elementary; the implementation is 30 lines of `tanh` and
`exp` in `tmesim/physics/subgrid.py::source_reweight`.

**It is a function of `(v, dx, supply, q, D)` and nothing else** — the local
vascular **area fraction**, the mesh spacing and the two rate constants, every
one of which the engine already holds per cell. This is why the problem is
tractable at all: `lambda` is computable per cell from quantities that are
already there. **There is no fitted coefficient, no training set and no lookup
table.**

`model="noflux"` closes the local cell instead (the textbook Krogh-cylinder unit
cell). It is retained and scored because it is the standard homogenisation
choice and it is **measurably worse** (§3).

### 1.4 The structural precondition, stated before any result

`beta == 1` identically wherever `v == 0` or `v == 1`: a cell that is entirely
vessel or entirely tissue **has no sub-cell structure to close over**. The
**shipped** vascular field is `_DEFAULT_VASCULAR_MODE = "discrete"` with
`geometry="binary"` and is **exactly 0.0 or exactly 1.0 everywhere** (measured:
2 unique values over the whole GBM section). **On the shipped field this closure
is an exact no-op by construction**, and §5 demonstrates that by measurement
rather than asserting it.

The binary painting discards the sub-grid geometry **before** the solve, so no
sub-grid closure *of the solve* can recover it. §5 therefore runs tissue on the
opt-in fractional-area field
(`discrete_vascular_field(geometry="subgrid")`, itself default-off, from
`results/vessel_geometry_fix.md`), with and without the closure, so the
comparison isolates the closure and not the painting.

---

## 2. The controlled problem, and a result that reframes the question

1-D two-region slab with the verified closed-form steady state from
`experiments/boundary_layer_resolution.py`, driven with the **oxygen** operator
(`k_in = supply + q`, `s_in = supply`, `k_out = q`) and solved with the
**shipped** `_build_laplacian`. Swept:

| axis | values | why |
|---|---|---|
| `lambda_sink = sqrt(D/q)` | 3, 5, 8, 10, 14, 20, 32, 50, 90, 160, **320** um | the brief's 3–320 um: 10–14 um is the jointly calibrated dense tile, 316 um the shipped one |
| `lambda_source = sqrt(D/supply)` | **3.2**, 6, 12, 30, 80, **200** um | 3.2 um is the jointly calibrated value, 200 um the shipped one |
| vessel half-width `a` | 3, **7.5**, 15 um | 7.5 um is `VESSEL_RADIUS_UM` |
| Krogh half-spacing | 24, 48, 96, 240, 480 um | see below |
| `dx` | 2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 30 um | the brief's 2–30 um |

**7,656 controlled points**, each scored against exact Simpson cell-averages of
the closed form.

**The half-spacing has to be swept, and that is not a detail.** With a single
vessel in a 480 um domain the whole field sits below the 0.2 band edge and the
*banded* fraction — which is what `F_state` is built from, and the harshest test
of a closure because it is a step function of the field — is degenerate. Swept,
**2,503 of 7,656 points** have a non-degenerate banded score. The section's own
2.68 % perfused area fraction at a 7.5 um radius implies a mean intervessel
spacing near 80 um, i.e. a half-spacing near 40 um, which is inside the range.

### 2.1 The finding that reframes the brief's premise

| configuration | `lambda_source` | `lambda_sink` | cells per `lambda_sink` at dx = 12 | **worst base error over dx = 2–30 um** |
|---|---:|---:|---:|---:|
| **shipped** (0.05, 0.02) | 200.0 um | 316.2 um | 26.35 | **0.170 %** |
| jointly calibrated, low (194, 10) | 3.21 um | 14.14 um | 1.18 | **62.8 %** |
| jointly calibrated, high (194, 20) | 3.21 um | 10.00 um | 0.83 | **73.8 %** |

> **Under the shipped oxygen constants there is no perivascular boundary-layer
> resolution problem to close.** The worst discretisation error anywhere on the
> ladder, at `dx = 30 um`, is **0.170 %**, and at the published `dx = 12 um` it
> is **0.018 %**.

This is not a small point. `results/grid_convergence_gbm.md` §5 measures the
shipped `F_state` swinging **0.0461** across `dx = 30/18/12/8`, and
`results/tail90_falsification.md` §6.1 measures a **3.72 %** fine-rung `F_state`
span and routes it to "whoever owns the oxygen solve and the cell-state
banding". **Neither can be perivascular under-resolution**, because the shipped
constants put 26 mesh cells inside the shortest length scale and the controlled
problem's error at that resolution is two parts in ten thousand. Those swings
belong to the *other* dx-dependent objects in the reconstruction — the binary
vessel painting, the `_tumor_density` binning and the nearest-neighbour sampling
of the field at fixed cell positions — and **no closure of the oxygen solve can
touch them.** The boundary-layer problem is real only for the **jointly
calibrated** constants, which is exactly where `oxygen_joint_grid_test.md` §4
located it.

---

## 3. The analytic baseline first, as the brief requires

### 3.1 The brief's named baseline — a power law in `(dx, lambda)` — FAILS

`results/boundary_layer_resolution.md` §4.2 fits
`|rel err of mean| = 1.427e-3 * dx^1.939 * delta^-1.145` at **R² = 0.995** and the
brief proposes an analytic correction derived from it as the obvious first thing
to try. Refitted here on the **train** grids over the parameter range that
actually occurs:

```
|rel err of mean| = 7.411e-5 * dx^2.061 * lambda_sink^-0.386        R^2 = 0.348
```

* **`p = 2.06` reproduces the clean second order** — that half of the earlier
  result survives intact and is confirmed on a different operator.
* **`R² collapses from 0.995 to 0.348.`** The earlier fit held `L_out` fixed at
  64.92 um and used one geometry; across the real ranges of `lambda_source`,
  vessel half-width and vessel spacing, **the error is not a power law in
  `(dx, lambda_sink)` at all.** The R² = 0.995 is a property of that slice, not
  of the operator.
* **Applied as a correction it makes things worse.** Scored on the **held-out**
  grids, and *given the exact sign of the error as free oracle information* (the
  measured sign is not constant, so a real correction would not have it), the
  median `|rel err of mean|` goes from **0.0034 to 0.0052**.
* **And it structurally cannot do the job.** It corrects a *scalar*. `F_state`
  and the hypoxic fraction are functionals of the **field** — a banded fraction
  is not recoverable from a corrected mean. Whatever else is true, this baseline
  could never have moved the quantity the falsification is scored on.

**So the brief's analytic baseline is refuted, and refuted on its own terms.**
What follows is a *different* analytic object: not a fit to the error, but a
closed form for the sub-cell physics.

### 3.2 The analytic closure, out of sample

Scored on the **held-out** grids `dx = 3, 8, 10, 16, 20, 30 um`, which no fit of
any kind has seen (the closure has no parameters, so this is a formality — but
it is the criterion the brief set):

| | median `|rel err of mean|` | p90 | max |
|---|---:|---:|---:|
| shipped solve (base) | 0.0034 | 0.0863 | **0.8088** |
| **closure, `open`** | **0.0011** | **0.0167** | **0.1047** |
| closure, `noflux` | 0.0025 | 0.0455 | 0.1835 |

Restricted to the regime the brief is about (`lambda_sink ≤ 14 um`,
`lambda_source ≤ 6 um` — the jointly calibrated constants), held-out grids only,
n = 530:

| | median | p90 | max |
|---|---:|---:|---:|
| base | 0.0539 | 0.3644 | **0.8088** |
| **`open`** | **0.0099** | **0.0389** | **0.1047** |
| `noflux` | 0.0208 | 0.0843 | 0.1424 |

Per grid, same regime (`*` = held-out):

| dx (um) | n | base median | **closed median** | base max | **closed max** |
|---:|---:|---:|---:|---:|---:|
| 2 | 140 | 0.35 % | **0.06 %** | 1.40 % | 0.18 % |
| 3 * | 140 | 1.11 % | 0.61 % | 3.60 % | 3.60 % |
| 4 | 140 | 0.43 % | **0.85 %** | 1.77 % | 2.11 % |
| 6 | 140 | 3.64 % | **0.29 %** | 10.13 % | 1.46 % |
| 8 * | 120 | 3.50 % | **1.59 %** | 18.23 % | 6.96 % |
| 10 * | 60 | 5.64 % | **0.63 %** | 25.38 % | 4.29 % |
| **12** | 120 | **10.43 %** | **0.84 %** | 31.65 % | 2.90 % |
| 16 * | 90 | 15.75 % | **1.21 %** | 44.38 % | 10.47 % |
| 20 * | 60 | 23.56 % | **1.27 %** | 56.16 % | 5.84 % |
| 24 | 90 | 28.35 % | **1.31 %** | 66.43 % | 3.05 % |
| 30 * | 60 | **36.98 %** | **0.71 %** | **80.88 %** | **3.19 %** |

**The closure is worth most exactly where the mesh is worst**, and it flattens
the error's dependence on `dx` from a 100-fold ramp (0.35 % → 36.98 %) to
roughly flat (0.06 % – 1.6 %). **Its one honest loss is at `dx = 4 um`**
(0.43 % → 0.85 %), where the residual is the flux term §1.2 identified and left
in place; the closure is a no-op there and the small degradation is the
second-order interaction between the two.

`noflux` is worse than `open` on every aggregate. The textbook Krogh unit cell
under-estimates `beta` because it forbids the leakage that actually occurs.

### 3.3 The grid dependence of the mean — the falsification criterion, met

Span of the relative error of the mean across each configuration's whole `dx`
ladder (924 ladders):

| | median span | p90 | max |
|---|---:|---:|---:|
| base | 0.0141 | 0.2689 | 0.8399 |
| **`open`** | **0.0059** | **0.0442** | **0.1166** |
| `noflux` | 0.0084 | 0.1055 | 0.1979 |

**The closure reduces the grid span of the mean on 288 of 289** non-degenerate
ladders. That is the falsification criterion, and it is met.

### 3.4 The banded score, and the irreducible floor — why NOTHING was trained

The quantity the falsification is scored on is `F_state`, a banded functional
(`0.35 x HYPOXIC + 0.5 x QUIESCENT + 1.0 x VIABLE`, the shipped `PDConfig`
multipliers). Its grid span across the whole `dx` ladder, over the 289 ladders
with a non-degenerate banded score:

| | mean span | median | p90 | max |
|---|---:|---:|---:|---:|
| **exact solution — the FLOOR** | **0.0310** | 0.0231 | — | — |
| shipped solve (base) | 0.0355 | 0.0258 | 0.0708 | 0.1422 |
| **closure, `open`** | **0.0317** | **0.0231** | **0.0625** | **0.0958** |
| closure, `noflux` | — | 0.0248 | 0.0692 | 0.1208 |

> **The `F_state` grid span has an irreducible floor, and the analytic closure
> lands essentially on it.** The **exact** solution's cell-averaged `F_state`
> *itself* spans **0.0310** on average across the `dx` ladder, because banding is
> a step function applied to a piecewise-constant representation and the
> representation changes with `dx`. Scored as **excess over that floor**, the
> mean goes from **0.00586 to 0.00142 — 76 % of the removable part removed**
> (74 % in the joint-calibrated-like regime; p90 excess 0.0177 → 0.0043, max
> 0.0880 → 0.0479).

**This is the reason nothing was trained, and it is a stronger reason than "the
analytic version was good enough".** A learned closure would have to beat a
**zero-parameter** formula that already sits close to the floor of what *any*
closure can achieve, because the residual is present in the exact solution
itself. Per the brief — *"a learned model that does not beat a two-parameter fit
is not a result"* — a learned model that cannot beat a **zero**-parameter formula
near its information-theoretic floor is less than that.

**No neural network, no regression, no lookup table was fitted anywhere in this
file.** The only fit is the two-parameter power law of §3.1, reported because the
brief asked for it, and it failed.

**Honest qualification.** The closure strictly reduces the `F_state` span on
**122 of 289** ladders, leaves it **unchanged on 115**, and **increases it on
52**. That is what approaching a floor from above with a step-function score
looks like: improving the field can push an individual cell across a band edge
either way. The aggregates improve; the per-ladder outcome does not always. The
mean-error result is much cleaner — the grid span of the mean falls on **288 of
289** ladders (§3.3).

---

## 4. Answering the brief's questions 1–4, before tissue

1. **What is the closure a function of, and what does it output?** A scalar per
   cell, `beta(v, dx, supply, q, D)`, multiplying the source term's contribution
   to the reaction diagonal. Chosen because a conservation identity localises the
   entire domain-mean error to the source weighting, and an oracle decomposition
   measured that a perfect `beta` removes **95 %** of it at `dx = 30 um` while a
   perfect effective-diffusivity/flux correction removes **7 %** (§1.2).
2. **Did the analytic baseline suffice?** The **brief's** analytic baseline — a
   correction derived from the `R² = 0.995` power law — **did not**: refitted over
   the ranges that actually occur it drops to `R² = 0.348`, applied out of sample
   with the exact sign given free it makes the error **worse** (0.0034 → 0.0052),
   and it structurally cannot move a banded fraction (§3.1). A **different**
   analytic object — a closed-form, zero-fitted-parameter sub-cell closure — did
   suffice.
3. **Was anything trained?** **No, and it would not have been justified.** §3.4.
4. **Out of sample on the controlled problem?** Yes, on grids
   `dx = 3, 8, 10, 16, 20, 30 um` seen by no fit: median `|rel err of mean|`
   **0.0034 → 0.0011**, p90 **0.0863 → 0.0167**, max **0.8088 → 0.1047**; grid
   span of the mean reduced on **288 of 289** ladders. **The falsification
   criterion on the controlled problem is met.**

---

## 5. The honest question: does it transfer to tissue? **No.**

Real 10x Xenium GBM section, 40,887 cells / 23,812 tumour, the ladder from
`results/grid_convergence_gbm.md` (`dx = 30/18/12/8 um`) plus a fine probe at 6
and 5 um. The oxygen solve is steady state, so **no `dt` and no CFL enters any of
this**. `dx = 4 um` was attempted and the 2.57 M-unknown direct solve exhausted
memory; the run checkpoints per rung, so 30/18/12/8/6/5 are complete.

### 5.0 Two controls first

**Reproduction.** The `binary | raw` arms reproduce `oxygen_joint_grid_test.md`
§1 exactly — joint-lo `F_state` **0.3936 / 0.4577 / 0.5249 / 0.5510** against its
published 0.39357 / 0.45771 / 0.52485 / 0.55102; joint-hi **0.3735 / 0.3911 /
0.4167 / 0.4277** against 0.37350 / 0.39106 / 0.41674 / 0.42768; shipped
**0.5041 / 0.5401 / 0.5502 / 0.5488**, which is also `grid_convergence_gbm.md`
§5. **It is the same pipeline producing the same numbers**, so anything that
moves below moves because of the closure.

**The no-op, demonstrated not asserted.** On the **shipped binary** vascular
field, `binary | closed` returns `F_state` = **0.5502 / 0.5249 / 0.4167** at
`dx = 12 um` for shipped / joint-lo / joint-hi — **identical to `binary | raw` to
four decimals**, exactly as §1.4 predicts. `beta == 1` everywhere because the
field is exactly 0 or 1. **The shipped painting destroys the sub-grid geometry
before the solve and no closure of the solve can recover it.**

Everything below therefore runs on the opt-in fractional-area field
(`geometry="subgrid"`, `radius_mode="constant"`, `radius_um = 7.5`), **raw versus
closed on the same field**, so the comparison isolates the closure.

### 5.1 The result

`F_state` on the real section, fractional-area vascular field:

| arm | dx 30 | dx 18 | dx 12 | dx 8 | dx 6 | dx 5 | **span (30–8)** | dx12 − finest |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| shipped, raw | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | **0.0000** † | +0.0000 |
| **shipped, closed** | 0.3500 | 0.3500 | 0.3687 | 0.4793 | 0.5000 | 0.5000 | **0.1293** | −0.1313 |
| joint-lo, raw | 0.6509 | 0.6141 | 0.5929 | 0.5734 | 0.5645 | 0.5589 | **0.0774** | +0.0340 |
| **joint-lo, closed** | 0.4638 | 0.3966 | 0.3767 | 0.3733 | 0.3810 | 0.3893 | **0.0905** | −0.0126 |
| joint-hi, raw | 0.5222 | 0.4755 | 0.4548 | 0.4383 | 0.4298 | 0.4252 | **0.0838** | +0.0297 |
| **joint-hi, closed** | 0.4270 | 0.3818 | 0.3658 | 0.3602 | 0.3601 | 0.3621 | **0.0668** | +0.0037 |

† The shipped/subgrid raw arm is **degenerate, not converged**: 100 % of tumour
cells are QUIESCENT at every rung, so `F_state` is pinned at exactly 0.5 and the
banded score has no dynamic range at all. Its zero span is not evidence of
anything.

**Stated plainly: the closure reduces the `F_state` grid span on ONE of the three
tissue configurations and increases it on the other two.**

* **joint-hi: it works.** Span **0.0838 → 0.0668** (−20 %), and the error at the
  published `dx = 12 um` against the finest rung collapses from **+0.0297 to
  +0.0037**, an **8-fold** reduction.
* **joint-lo: it fails.** Span **0.0774 → 0.0905** (+17 %), and the ladder stops
  being monotone — the closed arm turns around between `dx = 8` and `dx = 5` and
  climbs again. The `dx = 12` error does shrink (+0.0340 → −0.0126) but it
  **changes sign**, which is the same pathology `results/tail90_falsification.md`
  §0.3 reported for the last slab-derived correction this project tried.
* **shipped: it is destructive.** A flat (if degenerate) ladder becomes a
  **0.1293** span, and the closure drives the whole section from 100 % QUIESCENT
  to 100 % HYPOXIC at `dx = 30` and `18 um`.

**And the size of the intervention is the giveaway.** At `dx = 12 um` the closure
changes mean tumour O2 by **−62 %** on joint-lo (0.3110 → 0.1134) and **−62 %**
on joint-hi (0.1980 → 0.0756). On the controlled problem at the same `dx` the
closure moves the mean by about **6 %**. **A correction that is ten times larger
on tissue than on the problem it was derived from is not the same correction.**

### 5.2 Why — the mechanism, measured

`beta` degenerates to exactly `1/v` in the thin-saturated-source limit (measured
on the section: median `|beta − 1/v| / beta` = **0.0000**). On tissue the
fractional-area vascular field has a median `v` of **0.2292** over partial-vessel
cells, a 1st percentile of **0.0069**, and **28.7 %** of those cells below
`v = 0.1`. So `beta` reaches values the controlled problem never produced:

| | controlled sweep | real GBM section, `dx = 12 um` |
|---|---:|---:|
| median `beta` | 1.038 | **3.789** |
| p90 / p75 `beta` | 1.536 | 10.286 |
| **max `beta`** | **5.677** | **144.0** |

**38.7 % of partial-vessel cells on the section demand a `beta` larger than
anything the closure was ever validated at**, and those cells carry **9.8 %** of
the total source mass. In that regime the formula is no longer a boundary-layer
closure at all — it is the purely geometric factor `1/v`, which carries no
information about `lambda_source` or `lambda_sink` (measured: `beta` is identical
to three decimals across the shipped and both jointly calibrated arms, whose rate
constants differ by 4,000-fold).

**The physical reason is stated in the module's own failure-mode list and it is
the first item there: the local cell problem is 1-D Cartesian and a real vessel
is a disc.** A slab source spans the full cross-section of its cell and oxygen
escapes in one direction; a sub-cell disc lets oxygen escape radially, so the
true cell average sits higher relative to the vessel value and the true `beta` is
**smaller** than the 1-D formula gives. The closure therefore systematically
**over-corrects** on tissue, and it over-corrects hardest exactly where `v` is
smallest — which is most of the section.

### 5.3 The diagnosis, tested rather than asserted

If §5.2 is right, restricting the closure to the `beta` range it was validated
over should repair the shipped catastrophe. It was run: cap `beta` at the
controlled sweep's own maximum of **5.677** and let the closure **decline to act**
(`beta = 1`) on cells beyond it. `results/sweep/subgrid_clip.json`:

| arm | raw | closed | **closed, `beta` capped at 5.677** |
|---|---:|---:|---:|
| shipped | 0.0000 † | 0.1293 | **0.0056** |
| joint-lo | **0.0774** | 0.0905 | 0.1105 |
| joint-hi | 0.0838 | **0.0668** | 0.0779 |

**The diagnosis is confirmed and the closure is still not rescued.** Capping
removes the shipped arm's catastrophe almost entirely (0.1293 → 0.0056), which is
strong evidence that out-of-range `beta` was the cause. But on **joint-lo the
capped closure is worse than doing nothing at all** (0.1105 against a raw
0.0774), and on joint-hi capping is worse than not capping. **There is no version
of this closure — full or restricted to its validated range — that reduces the
grid dependence of `F_state` on more than one of the three tissue
configurations.**

### 5.4 What did NOT transfer, precisely

| claim | controlled problem | real GBM section |
|---|---|---|
| removes most of the mean's grid error | **yes** — 84–97 %, out of sample | not measurable (no exact answer on tissue) |
| reduces the grid SPAN of the score | **yes** — 288/289 ladders (mean), floor-limited (`F_state`) | **1 of 3 configurations** |
| sign of the benefit is stable | **yes** | **no** — improves joint-hi, worsens joint-lo and shipped |
| operates inside its validated range | **yes**, by construction | **no** — 38.7 % of source cells beyond it, `beta` to 144 vs 5.68 |
| magnitude of the correction | ~6 % of the mean at `dx = 12` | **~62 %** of the mean at `dx = 12` |

**This is the same failure shape `results/tail90_falsification.md` recorded**, and
it is worth naming the pattern rather than the instance: a correction derived on
a controlled 1-D problem encodes that problem's *geometry* as well as its
*physics*, and this project's tissue differs from the slab in geometry far more
than in physics.

---

## 6. Verdict

**Three things are established and they must be said in this order.**

**1. The premise that the mesh limit explains the shipped ladder is FALSE.** Under
the shipped oxygen constants (`supply = 0.05`, `consumption = 0.02`) the
controlled problem's worst discretisation error over `dx = 2–30 um` is
**0.170 %**, and **0.018 %** at the published `dx = 12 um` (§2.1). There are 26
mesh cells inside the shortest length scale. **`grid_convergence_gbm.md` §5's
`F_state` span of 0.0461 and `tail90_falsification.md` §6.1's 3.72 % fine-rung
span are not perivascular under-resolution and no sub-grid closure of the oxygen
solve can address them.** They belong to the other dx-dependent objects in the
reconstruction — the binary vessel painting, the `_tumor_density` binning, and
the nearest-neighbour sampling of the field at fixed cell positions. **That is
the single most actionable finding in this file**, and it re-routes a question
two other reports have been trying to answer with resolution.

**2. On the controlled problem the closure works, and no machine learning was
needed or justified.** A closed-form, **zero-fitted-parameter** source
reweighting removes **84–97 %** of the grid error on grids no fit has seen,
flattens a 100-fold `dx` ramp to roughly flat, reduces the grid span of the mean
on **288 of 289** ladders, and removes **76 %** of the `F_state` span's excess
over an **irreducible floor** that the exact solution itself carries. The brief's
own suggested analytic baseline — a correction from the `R² = 0.995` power law —
**failed** (`R² = 0.348` over the real ranges; worse out of sample even with the
sign given free; structurally unable to move a banded fraction). **Nothing was
trained, and training would not have been justified**, because there is no
headroom between a zero-parameter formula and the floor for a model to occupy.

**3. It does not transfer to tissue, and the mesh limit is real in the only sense
that matters operationally.** The closure reduces the `F_state` grid span on
**one of three** tissue configurations and increases it on the other two; the
sign of its benefit is configuration-dependent; it changes the answer by 62 %
while improving grid-independence by at most 20 %; and **38.7 %** of the
section's partial-vessel cells sit outside every `beta` it was validated at.
Restricting it to its validated range confirms the mechanism and still does not
rescue it. **I cannot certify a closure that removes the mesh dependence on
tissue, and this one does not.**

**The honest summary is the third of the brief's three good outcomes, with a
qualification that makes it more useful than a flat negative:** *the mesh limit
identified by `oxygen_joint_grid_test.md` is real for the jointly calibrated
constants, a sub-grid closure of the source term is demonstrably the right
mathematical object for it, an analytic zero-parameter version of that object
solves the controlled problem outright — and it fails on tissue for an
identified, measured, geometric reason (a 1-D slab cell problem applied to 2-D
sub-cell discs) rather than an unknown one.*

**What would be needed next, named but not claimed.** A **radial/2-D local cell
problem** (a disc source in a square cell, matched to a decaying far field) in
place of the 1-D slab, which would give a smaller and better-behaved `beta` at
small `v` — exactly where §5.2 shows the failure lives. That is a derivation, not
a model, and it is testable by the same falsification this file used. **Nothing
in this report argues for adopting the current closure, and it stays off.**

---

## 7. Provenance, controls, and what this does NOT establish

* **Nothing shipped changed.** `OxygenConfig.supply = 0.05`,
  `consumption = 0.02`, `hypoxia_threshold = 0.2`, `PDConfig` multipliers,
  `physics/oxygen.py`, `_DEFAULT_VASCULAR_MODE = "discrete"` and
  `discrete_vascular_field(geometry="binary")` are untouched. The closure is a
  **new module, `tmesim/physics/subgrid.py`, that nothing else in `tmesim/`
  imports** — verified by grep — and `solve_oxygen_subgrid(..., enabled=False)`
  is a bit-identical passthrough. Everything is left unstaged.
* **`pytest`**: the suite was re-run after the coordinator reverted an unrelated
  `tmesim/molecular/` change (measured-vs-computed logP, which had drifted
  `alpha` from 148.60815765445545 to 68.637 and broken six tests including two
  bitwise PK reference pins). The new module was independently confirmed not to
  be the cause: with `subgrid.py` moved aside, the same tests still failed.
* **One section, one tumour type.** Every tissue number is `xenium_gbm`.
* **The controlled problem is 1-D Cartesian.** That is the whole of §5.2 and it
  is a limitation of the closure, not only of the test.
* **The closure is derived for the `linear` sink only** and raises on
  `sink_form="michaelis_menten"`, which has no closed-form local cell problem.
  `oxygen_joint_grid_test.md` limitation 3 notes a Michaelis–Menten sink would
  lengthen the effective diffusion length in exactly the hypoxic regions where
  the layer is thinnest, and would therefore plausibly *reduce* the problem this
  file is trying to close. **That was not tested and remains the single most
  likely way §2.1's regime split could change.**
* **No transport solve, no drug, no verdict.** This file measures an oxygen field
  and its banding. It does not re-score any 3-D:2-D ratio and does not touch
  `kill_calibration`.
* **The jointly calibrated constants are inherited, not endorsed.**
  `supply = 194` and `consumption = 10/20` are `oxygen_supply_calibration.md`'s,
  with all of its UNVERIFIED flags, and the 7.5 um vessel-painting radius is
  known-wrong (`results/vascular_segmentation_audit.md`).
* **`dx = 4 um` was not reached on tissue** (memory), so the finest tissue rung is
  5 um and "converged" is never demonstrated on the section, only approached.
* **The fractional-area painting is a different source field from the binary
  one** (source integral 6.192e5 vs 7.487e5 um^2 at `dx = 12`), so `subgrid` arms
  are not comparable to `binary` arms. Every closure comparison in §5 is
  raw-versus-closed **on the same field**, which is why that comparison is the
  only one drawn.

---

## 8. Reproduce

```
./.venv/bin/python experiments/subgrid_closure.py identity sweep analytic transfer clip
```

Each stage checkpoints to `results/sweep/subgrid_<stage>.json`; `transfer`
additionally writes `subgrid_transfer_partial.json` after every rung and is
restart-safe. `identity` is the stage that decides the closure's form and it
runs in 0.1 s with no tissue.
