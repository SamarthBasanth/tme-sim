# What actually makes `F_state` move with `dx`? Not the boundary layer, and not the sampler either.

> **One-line answer.** The `F_state` grid swing is **the binary vessel-source
> rasterisation on the WRITE side, amplified by a band edge that cuts through
> the densest part of the population.** Holding the oxygen field completely
> fixed and varying only how cells read it reproduces **1.3 %** of the movement;
> varying only the dx-rasterised **source** while solving on a fixed fine grid
> reproduces **104 %** of it. Formal variance attribution:
> **source rasterisation +1.057, PDE discretisation −0.056, cell-to-grid
> sampling −0.001** of `var(F_state)`.
>
> **Both previously named mechanisms are wrong.** It is not perivascular
> under-resolution (`results/subgrid_closure.md` §2 already killed that: 0.170 %
> discretisation error at `dx = 30 um`), and it is **not** the cell-to-grid
> sampling rule either — which was this brief's leading suspicion and is
> **refuted here by direct experiment**.
>
> **Is it fixable? Partly, and the two halves have opposite answers.**
> The *driver* is a defect and is cheaply fixable: the shipped binary
> centre-in-disc painting does not conserve the painted vessel area under `dx`
> (**4.73 %** span), while the already-present, opt-in `geometry='subgrid'`
> partial-area painting does (**0.077 %** span, **61x** more stable).
> The *amplifier* is **not** a defect and is **not** fixable by refining
> anything: **10.10 %** of tumour cells sit within **0.01** of the `0.4`
> QUIESCENT ceiling and the whole tumour population spans only
> **0.2924–0.4237** in O2, so `dq/d(mean O2)` is **6.72 per unit O2**. Any
> residual field movement is multiplied by that. **That half is the honest
> resolution limit of banding a continuous field, and the response to it is to
> stop using hard bands — not to refine the grid.**

**Status: a diagnosis. Nothing shipped was changed.** No default, band edge or
constant was touched. `OxygenConfig` (`supply = 0.05`, `consumption = 0.02`,
`hypoxia_threshold = 0.2`), `PDConfig`, `real_fields._DEFAULT_VASCULAR_MODE` and
`discrete_vascular_field`'s `geometry='binary'` default are all as shipped.
`calibrate_alpha()` returns **148.60815765445545** (§7). Nothing in `tmesim/`
was edited; the `geometry='subgrid'` arm of §6 is an existing opt-in branch
reached by passing a keyword.

* Generator: `experiments/fstate_mechanism.py`
* Raw: `results/sweep/fstate_identity.json`, `fstate_ladder.json`,
  `fstate_fixed_field.json`, `fstate_raster.json`, `fstate_band_edges.json`,
  `fstate_source_fix.json`, `fstate_o2_cells.npz` (per-cell O2 at every rung)
* Section: `data_cache/xenium_gbm`, 40,887 cells / **23,812 tumour**, the same
  section `results/kill_calibration.md` and `results/grid_convergence_gbm.md` use

---

## 0. The question, and why the two standing answers are both wrong

`results/grid_convergence_gbm.md` §5 measures `F_state` — the oxygen-banded
cell-state multiplier that multiplies straight into every scored 3-D:2-D ratio —
swinging **0.5041–0.5502** across the dx ladder. `results/tail90_falsification.md`
§6.1 measures a **3.72 %** fine-rung span carrying 25–51 % of the whole ladder's
variance. Both attribute it to perivascular under-resolution.

`results/subgrid_closure.md` §2 falsified that: under the shipped oxygen
constants the worst discretisation error on a controlled problem with an exact
solution is **0.170 %** at `dx = 30 um`. A field that accurate cannot move a
banded fraction by 4 points on its own.

So the field is essentially right and the label still moves. This file finds out
why by decomposition rather than by argument.

### 0.1 `F_state` is a one-variable problem, which makes this tractable

The measured HYPOXIC fraction is **exactly 0.0000 at every rung** (§4 shows why:
no tumour cell anywhere on this section is within 0.09 of the `0.2` hypoxia
edge). With no hypoxic cells,

```
F_state  =  1 - 0.5 * q          q = QUIESCENT fraction of tumour cells
```

exactly — `0.8996` quiescent gives `0.5502`, the published number. So the entire
question reduces to: **what fraction of tumour cells have `local_o2 < 0.4`, and
why does it move with `dx`?**

---

## 1. The decomposition, stated before it was run

`F_state` is a functional of one thing only: the per-tumour-cell oxygen value
`o_i` that `cells.agents.update_cell_states` bands. Exactly four separable
things stand between the continuum problem and `o_i`, and each is switched on
one at a time:

| arm | field solved on | cells read it by |
|---|---|---|
| **D** | reference grid (`dx = 5 um`), reference source | bilinear at the cell's true `(x_i, y_i)` |
| **E(dx)** | reference grid, **dx-rasterised source** | bilinear at the cell's true `(x_i, y_i)` |
| **B(dx)** | **dx grid**, dx-rasterised source | bilinear at the cell's true `(x_i, y_i)` |
| **A(dx)** | dx grid, dx-rasterised source | **production floor index** (`Domain.to_index`) |

Reading down the table each row switches exactly one mechanism on, so the four
arms telescope into an exact additive identity:

```
A(dx) - D  =  [E(dx) - D]  +  [B(dx) - E(dx)]  +  [A(dx) - B(dx)]
              source raster    discretisation      cell-to-grid sampling
```

**`A` is the shipped number.** `D` is what every arm would return if nothing
depended on the mesh.

### 1.1 Preconditions (`--stage identity`)

Everything below depends on the harness reproducing the shipped pipeline, so
that is asserted first, not assumed:

| check | result |
|---|---|
| local `solve_from_fields` vs shipped `physics.oxygen.solve_oxygen` | **bit-identical**, `max abs diff = 0.000e+00` |
| local `band_states` + floor sampler vs shipped `update_cell_states` | **identical labels, all 23,812 cells** |
| `F_state` at `dx = 12 um` | **0.5502**, the published value |
| shipped `fields.vascular_density` vs `geometry='binary'` | **equal at all 10 rungs** (§6) |

The local solver exists only so an arbitrary `(v, rho)` pair can be fed to it —
it is `solve_oxygen`'s `sink_form='linear'` branch line for line.

---

## 2. THE DECISIVE CONTROL: hold the field fixed, vary only the read

The brief's leading suspicion was the cell-to-grid sampling. `Domain.to_index`
is a **floor**, so a cell is handed the O2 stored at the centre of whichever
grid cell contains it — a point up to `dx/sqrt(2)` away. Cells sit at fixed
measured positions; the grid does not. That is dx-dependent and has nothing to
do with the PDE, so it is a clean and cheap thing to test decisively.

**The test.** One oxygen field, solved once at `dx = 5 um`, **never re-solved**,
re-read at every rung through the production floor rule. `C_avg` block-averages
it onto the dx grid first (the finite-volume restriction, i.e. what a coarse
cell-centred value *means*); `C_inj` point-injects. Nothing about the PDE, the
source, or mesh convergence can enter these numbers.

| `dx` (um) | mean cell→node displacement (um) | mean per-cell \|ΔO2\| | `C_avg` `F_state` | `C_inj` `F_state` |
|---:|---:|---:|---:|---:|
| 30 | 11.5 | 0.00036 | 0.5443 | 0.5442 |
| 24 | 9.1 | 0.00029 | 0.5442 | 0.5441 |
| 18 | 6.9 | 0.00023 | 0.5442 | 0.5442 |
| 16 | 6.1 | 0.00020 | 0.5440 | 0.5440 |
| 12 | 4.6 | 0.00016 | 0.5440 | 0.5439 |
| 10 | 3.8 | 0.00012 | 0.5438 | 0.5438 |
| 8 | 3.1 | 0.00010 | 0.5439 | 0.5439 |
| 7 | 2.7 | 0.00009 | 0.5438 | 0.5438 |
| 6 | 2.3 | 0.00009 | 0.5438 | 0.5438 |

> **`C_avg` span = 0.0006. `C_inj` span = 0.0005. The production span is 0.0461.
> Sampling reproduces 1.3 % of the movement. THE LEADING SUSPICION IS REFUTED.**

The independent cross-check agrees: on the **real** ladder, re-reading each
rung's own field with a dx-independent bilinear sampler instead of the
production floor sampler changes `F_state` by at most **0.00048**, and the
bilinear ladder's span (**0.0461**) is identical to the production ladder's
(**0.0461**) to four decimals.

**Why sampling is negligible, mechanistically.** The mean `|grad O2|` at tumour
cells is **4.75e-05 per um**. Over even the `dx = 30 um` mean displacement of
11.5 um that is an O2 error of ~5e-4 — three orders of magnitude below the
0.0123 level shift §3 finds. The shipped oxygen field is *smooth on the scale of
the grid*, which is the same fact `subgrid_closure.md` §2 reports as a 0.170 %
discretisation error. A smooth field is forgiving of where you sample it.

---

## 3. What it IS: the source rasterisation, on the write side

Arm **E** solves on the fixed `dx = 5 um` reference grid with **only the
dx-rasterised source** injected, and reads with the dx-independent sampler. It
therefore contains the write-side rasterisation and nothing else.

| `dx` (um) | **A** production | **B** bilinear read | **E** dx-source, fine solve | **C** fixed field | source `E−D` | discretisation `B−E` | sampling `A−B` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 0.5041 | 0.5042 | 0.5046 | 0.5443 | **−0.0391** | −0.0004 | −0.0001 |
| 24 | 0.5392 | 0.5397 | 0.5395 | 0.5442 | **−0.0041** | +0.0002 | −0.0005 |
| 18 | 0.5401 | 0.5403 | 0.5397 | 0.5442 | **−0.0039** | +0.0006 | −0.0002 |
| 16 | 0.5333 | 0.5334 | 0.5305 | 0.5440 | **−0.0132** | +0.0030 | −0.0001 |
| 12 | 0.5502 | 0.5504 | 0.5527 | 0.5440 | **+0.0091** | −0.0024 | −0.0001 |
| 10 | 0.5305 | 0.5308 | 0.5310 | 0.5438 | **−0.0126** | −0.0002 | −0.0003 |
| 8 | 0.5488 | 0.5490 | 0.5523 | 0.5439 | **+0.0086** | −0.0033 | −0.0001 |
| 7 | 0.5407 | 0.5408 | 0.5479 | 0.5438 | **+0.0042** | −0.0071 | −0.0001 |
| 6 | 0.5390 | 0.5392 | 0.5334 | 0.5438 | **−0.0102** | +0.0058 | −0.0002 |

`D = 0.5437` (dx-independent baseline, `q = 0.9127`).

**Attribution.** Because the identity of §1 is additive and exact, each term's
share of the production variance is `cov(term, A)/var(A)`, which sums to exactly
1 and — unlike a span ratio — prices cancellation:

| term | span | % of `A` span | **`cov/var(A)`** | mean \|term\| |
|---|---:|---:|---:|---:|
| **source rasterisation** `E−D` | 0.0481 | 104.4 % | **+1.057** | 0.01167 |
| PDE discretisation `B−E` | 0.0128 | 27.8 % | **−0.056** | 0.00254 |
| cell-to-grid sampling `A−B` | 0.0004 | 0.8 % | **−0.001** | 0.00021 |

> **The source rasterisation alone over-explains the movement (+1.057); the PDE
> discretisation slightly cancels it (−0.056); the sampler does nothing
> (−0.001).**

Restricted to the fine rungs (`dx <= 18 um`), which is what
`tail90_falsification.md` §6.1's **3.72 %** refers to: production span **3.65 %**,
source-only span **4.11 %**, fixed-field sampling span **0.08 %**.

### 3.1 The write-side defect, named exactly

`vessel_network.discrete_vascular_field(geometry='binary')` paints a grid cell as
vessel **iff its centre lies within `VESSEL_RADIUS_UM = 7.5 um` of an endothelial
cell**, and the painted value is exactly 1 or 0. A 7.5 um disc has area
**176.71 um²**. A grid cell has area `dx²`:

| `dx` (um) | cell area (um²) | disc area / cell area | painted vessel area (um²) |
|---:|---:|---:|---:|
| 30 | 900 | **0.196** | **713,700** |
| 24 | 576 | 0.307 | 748,800 |
| 18 | 324 | 0.545 | 746,496 |
| 16 | 256 | 0.690 | 745,728 |
| 12 | 144 | **1.227** | 748,656 |
| 10 | 100 | 1.767 | 739,200 |
| 8 | 64 | 2.761 | 746,816 |
| 7 | 49 | 3.606 | 741,517 |
| 6 | 36 | 4.909 | 743,904 |
| 5 | 25 | 7.069 | 744,000 |

At `dx = 30 um` an entire vessel is **one fifth of one grid cell**, so it is
quantised to 0 or 1 cells: either missed entirely, or painted **5.1x too large**.
At the published `dx = 12 um` it is 1.2 cells. The result is that the **total
painted vessel area is not conserved under `dx`** — it spans **4.73 %** — and it
moves **non-monotonically**, which is precisely the non-monotone, oscillating
character `grid_convergence_gbm.md` §0 flags and cannot explain with a
convergence argument.

**This is the whole chain, and it correlates almost perfectly** (over the full
10-rung ladder `dx = 30…5 um` of the table above; the 9-rung decomposition set of
§3, which drops `dx = 5 um` as the reference, gives `corr(area, meanO2) = 0.9872`):

| quantity | span across the ladder | correlation |
|---|---|---|
| painted vessel area | 713,700–748,800 um² (**4.73 %**) | — |
| mean tumour O2 | 0.3364–0.3486 (**0.0123**) | `corr(area, meanO2)` = **0.9868** |
| `F_state` | 0.5041–0.5502 (**0.0461**) | `corr(meanO2, F_state)` = **0.9357** |

**The sink is not implicated.** `_tumor_density`'s areal-density normalisation
works exactly as intended: the integrated tumour sink is
**23812.00000000** at every single rung, conserved to **2.9e-11**. Only its
spatial spread changes with `dx`, and that is the small `B−E` term. The
grid-independence fix documented in `physics/oxygen.py` is doing its job; the
*source* never received the equivalent fix.

---

## 4. The amplifier: a band edge through the densest part of the population

A 0.0123 shift in mean O2 is 3.5 % of the mean. It produces a **0.0922** swing
in the quiescent fraction — a 7.5x amplification. That factor is not a bug in
anything; it is where the band edge sits relative to the population.

Measured on the reference field, over all 23,812 tumour cells:

| statistic | value |
|---|---:|
| mean tumour O2 | **0.3476** |
| median | 0.3453 |
| **min / max** | **0.2924 / 0.4237** |
| p1 / p99 | 0.2957 / 0.4164 |
| p90 / p95 | **0.3976 / 0.4046** |
| fraction below the `0.4` QUIESCENT ceiling | **0.9127** |
| fraction below the `0.2` HYPOXIC edge | **0.0000** |

> **The entire tumour population lives in `O2 ∈ [0.2924, 0.4237]` — a full range
> of just 0.1313 — and the `0.4` band edge cuts it between p90 and p95.**

How hard the population is piled against that edge:

| within ±ε of the `0.4` ceiling | fraction of tumour cells |
|---:|---:|
| 0.001 | 1.18 % |
| 0.002 | 2.24 % |
| 0.005 | 6.27 % |
| **0.01** | **10.10 %** |
| 0.02 | 15.81 % |
| 0.05 | 41.22 % |

The corresponding local population density at the ceiling is **5.05 per unit
O2**, and the empirical `dq/d(mean O2)` regressed across the ladder is **−6.72
per unit O2**. The two agree to the accuracy one should expect given the shift
is not perfectly uniform, and they close the loop numerically:

```
predicted q swing = 0.0123  x  6.72  =  0.0824
observed  q swing =                     0.0922      (89 % accounted for)
```

**Nothing is near the `0.2` edge** (the closest tumour cell sits at O2 = 0.2924, **0.0924** above it), which is why
the HYPOXIC band is empty at every rung and `hypoxia_resistance` acts on no cell
at all — consistent with `results/cell_state_impact.md` and
`results/state_validation.md`.

---

## 5. So: fixable, or irreducible? Both, in different halves

**The driver is a defect and is fixable.** The binary centre-in-disc painting
fails to conserve the one quantity the oxygen solve is most sensitive to. That
is not a resolution limit, it is a rasterisation that was never made
area-conserving — the exact counterpart of the `dx²` sink fix that
`physics/oxygen.py` already documents, applied to the source instead of the sink.

**The amplifier is not a defect and is not fixable by refining anything.** With
10 % of the population within 0.01 of a hard band edge and the whole population
spanning 0.13 in O2, *any* residual perturbation of the field — from a better
source, a finer mesh, a different sink form, or a genuinely better physical model
— is multiplied by ~6.7 into the quiescent fraction and by ~3.4 into `F_state`.
Refining the grid does not reduce that gain; it is a property of putting a step
function through a smooth unimodal distribution at its densest point.

> **The honest framing: fixing the source rasterisation would remove the
> *forcing*, not the *sensitivity*. A banded cell-state model on this section is
> intrinsically a high-gain amplifier of any change to the oxygen field, and the
> response to that is to stop using hard bands — a graded, continuous
> `resistance(O2)` in place of the `VIABLE / QUIESCENT / HYPOXIC` step — not to
> refine the mesh.** That is a modelling decision with its own validation
> burden, it is out of scope for a diagnosis, and **nothing here implements it.**

---

## 6. Evidence that the driver really is fixable (measured, not asserted)

`discrete_vascular_field` already carries an opt-in `geometry='subgrid'`
partial-area painting, whose docstring *claims* the painted total is "computed
without reference to `domain.dx`, so it is grid-independent by construction".
That is a claim, so it was measured. **No default was changed** — `geometry` is
passed as a keyword to an existing branch, and the shipped binary arm is
re-measured on the identical reconstruction at every rung (and confirmed equal
to `micro.fields.vascular_density` at all 10 rungs).

| `dx` (um) | **binary** painted area (um²) | **subgrid** painted area (um²) |
|---:|---:|---:|
| 30 | 713,700 | 347,327 |
| 24 | 748,800 | 347,349 |
| 18 | 746,496 | 347,499 |
| 16 | 745,728 | 347,468 |
| 12 | 748,656 | 347,340 |
| 10 | 739,200 | 347,230 |
| 8 | 746,816 | 347,381 |
| 7 | 741,517 | 347,481 |
| 6 | 743,904 | 347,396 |
| 5 | 744,000 | 347,307 |

> **Span across the ladder: binary 4.731 %, subgrid 0.077 % — 61x more stable.**

**What this does and does not establish.** It establishes that the dx-dependence
of the *painted source area* — the quantity §3 shows carries `cov/var = +1.057`
of the `F_state` variance — is removable, cheaply, with code already in the
repository. It does **not** establish the resulting `F_state` ladder, and this
file does not run one: the two paintings deposit **different total areas**
(347k vs 744k um², because `subgrid` paints one footprint per segmented vessel
*object* at a data-derived radius rather than one 7.5 um disc per endothelial
cell), so switching would move the absolute oxygen level and re-open
`OxygenConfig.supply`, every downstream `F_state`, and every published verdict.
**That is a deliberate re-run, not a bug fix, and it is not made here.** The
experiment that would settle it is the §3 ladder re-run under
`geometry='subgrid'` with `supply` re-anchored — named here so it is not
mistaken for done.

---

## 7. Preconditions and hygiene

* `calibrate_alpha()` returns **148.60815765445545** — checked before and after
  this work.
* **`tmesim/` was not modified.** Every shipped default, band edge and constant
  is untouched; this file is a measurement of the pipeline as it ships.
* The oxygen solve is unconditionally stable (steady, sparse-direct, no `dt`), so
  the ladder is not bound by the advection CFL floor that stops
  `grid_convergence_gbm.md` at `dx = 7 um`; rungs at 6 and 5 um are included
  here and no drug transport solve was run.
* `results/sweep/fstate_o2_cells.npz` carries the per-cell O2 vector at every
  rung under both samplers, so §2–§4 can be recomputed without re-solving.

---

## 8. What this changes in the standing record

Two documents name the wrong mechanism. **Their measurements stand and are
reproduced here exactly** — `grid_convergence_gbm.md` §5's `0.0461` swing is
reproduced to four decimals, and its `dx = 12 um` `F_state = 0.5502` is
reproduced exactly. **Only the mechanism attribution is wrong**, and a
correction note has been added to each:

* `results/grid_convergence_gbm.md` §5 — attributes the swing to the oxygen PDE
  being solved on the same grid. **This file is machine-generated by
  `experiments/grid_convergence_gbm.py` (`write_report`, §5 emitted at line
  1588), so the hand-added note there WILL BE DESTROYED on the next run of that
  script.** The note says so in its own text and names its own re-insertion
  point. `experiments/fstate_mechanism.py` does not own that generator and did
  not edit it.
* `results/tail90_falsification.md` §6.1 — attributes the `F_state` span to
  "something dx-dependent in the reconstruction", following
  `results/boundary_layer_resolution.md` §6. That is *half right and was the
  closest anyone got*: it is the reconstruction, specifically the vessel source
  rasterisation — but it is not the boundary layer. Despite its header claiming
  "Generated by `experiments/tail90_falsification.py`", that script writes only
  JSON and **no markdown**, so that file is hand-maintained and the note there
  is durable.
