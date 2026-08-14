# The vessel SOURCE never got the fix the SINK already has. Here it is, and it removes 98.8 % of the `F_state` grid drift.

> ### One-line answer
>
> **`physics/oxygen.py` fixed the oxygen SINK by refusing to feed a raw
> per-grid-cell COUNT into the PDE — it divides by the grid-cell area and
> re-expresses the result against a fixed physical tile. The SOURCE has the
> identical defect in the identical place and never received the identical
> treatment.** The shipped binary painting asks each grid cell a YES/NO question
> ("is your centre within 7.5 µm of an endothelial cell?") and then deposits a
> WHOLE grid cell of vessel, area `dx²`, whenever the answer is yes. That is a
> count dressed as a field, and the total painted vessel area consequently
> **spans 4.7312 %** across `dx = 30…5 µm` — reproduced here independently,
> row-for-row against `results/fstate_mechanism.md` §3.1.
>
> **The fix is the same move, applied to the source: paint the covered-area
> FRACTION instead of a binary token.** Same discs, same `VESSEL_RADIUS_UM =
> 7.5 µm`, same endothelial cells, same tissue mask — only the *measurement* of
> the disc changes. Painted area then becomes **742,713.0 µm² at every single
> rung, a span of 6.3e-14 % (exact to floating point, not merely convergent)**,
> because the sub-pixel lattice it is measured on is anchored in physical space
> and the grid only decides how its points are bucketed.
>
> **It removes 98.82 % of the `F_state` drift.** `F_state` span
> **0.04609 → 0.000546**; mean tumour O2 span **0.012269 → 0.000117**. The
> `binary` control arm reproduces the published ladder exactly (0.5041–0.5502,
> `dq/d(meanO2) = −6.72`, `corr(area, meanO2) = 0.9868`, `corr(meanO2, F_state)
> = 0.9357` — all four to the published digits).
>
> **IT DOES CHANGE THE TOTAL, AND THAT IS THE FINDING THAT DECIDES THE MERGE.**
> At the production `dx = 12 µm` the shipped total is 748,656 µm² and the
> conserved total is 742,713 µm² — **−0.794 %**. It cannot be otherwise: on this
> ladder the shipped total is a *different number at every rung*, so "conserve
> the shipped total" is not a well-posed constraint, and any single conserved
> value must disagree with nine of the ten rungs. **Conserving the total and
> removing the drift are in direct conflict, and this is not a matter of
> implementation.** The `area_matched` arm proves it by experiment: holding the
> tissue-integrated total EXACTLY at the shipped per-rung value and changing only
> the distribution removes just **24.10 %** of the drift, because the shipped
> total *is* the defect.
>
> **The second half of the drift is untouched and is not fixable by this or
> anything else.** Under the fixed source **9.78 %** of tumour cells still sit
> within 0.01 of the 0.4 band ceiling (9.87 % before), the local population
> density at the ceiling is **4.89 per unit O2** (4.94 before), and the whole
> tumour population still spans only **0.1309** in O2. The fix removes the
> *forcing*; the *gain* is undiminished (marginally higher, if anything — §5).
> The residual 0.000546 is that undiminished gain acting on the 0.000117 of field
> movement the PDE discretisation and the sampler still contribute.
>
> **The nine published verdicts do not move.** Re-scored through
> `delivery_source_geometry.py`'s own blast harness at `dx = 12 µm`: **zero of
> nine verdicts change, zero of the seven scored rows flip**, largest ratio move
> 2.62 %. The project's only falsification (breast docetaxel FAIL) survives, and
> so does the flagship (GBM temozolomide PASS). For contrast, the same harness
> measured `geometry='subgrid'` flipping three of seven.
>
> **RECOMMENDATION: land it DEFAULT-OFF, as shipped here.** `geometry='binary'`
> is untouched and bit-identical; `geometry='area'` is opt-in and reachable only
> by keyword. Promoting it to the default is a **deliberate −0.794 % source
> recalibration** and must be run as one. §10 names the three things that re-run
> needs; §8 records why the ~0.8 % dose perturbation is not the reason to
> hesitate.

**Status: a fix, landed opt-in, plus the measurement of what promoting it would
cost.** No band edge, no `OxygenConfig.supply`, no `consumption`, no
`hypoxia_threshold`, no `VESSEL_RADIUS_UM`, no `LINK_DISTANCE_UM`, no
`_DEFAULT_VASCULAR_MODE` and no `DX_UM` was changed. `calibrate_alpha()` returns
**148.60815765445545** — checked before and after (§9).

* Generator: `experiments/vessel_source_conservation.py`
* Code: `tmesim/reconstruct/vessel_network.py`
  (`_area_coverage_field`, `SOURCE_SUBPIXEL_PITCH_UM`, `geometry='area'`)
* Tests: `tests/test_vessel_source.py`
* Raw: `results/sweep/vsrc_area_ladder.json`, `vsrc_fstate_ladder.json`,
  `vsrc_band_edge.json`, `vsrc_verdicts.json`, `vsrc_hygiene.json`,
  `vsrc_verdicts_dx12.jsonl`
* Section: `data_cache/xenium_gbm`, 40,887 cells / 23,812 tumour /
  **4,295 endothelial → 3,091 vessel objects** — the same section
  `results/fstate_mechanism.md`, `kill_calibration.md` and
  `grid_convergence_gbm.md` use. Breast rows from
  `data_cache/xenium_breast_janesick/rep1`.

---

## 0. The model: what the SINK's fix actually does

Read this first, because the source needs the same treatment and not a new
invention.

`physics/oxygen.py` solves

```
D * lap(O) + supply * vasc * (1 - O)  -  consumption * rho_tum * O  =  0
                └── SOURCE ──┘             └──── SINK ────┘
```

`_tumor_density` builds `rho_tum`. An earlier revision binned tumour cells with
`np.add.at` and fed the raw per-grid-cell COUNT straight into the sink. A count
scales as `dx²`, so halving the grid quartered the consumption term: the module's
own docstring records the hypoxic tumour fraction moving 0.937 / 0.364 / 0.134 at
`dx = 30 / 20 / 12 µm`. The fix is two lines and one idea:

```python
np.add.at(density, (row, col), 1.0)
density *= TUMOR_REFERENCE_AREA_UM2 / (float(domain.dx) ** 2)
```

**Divide the count by the grid-cell area to make it an areal DENSITY, then
re-express it against a fixed physical reference tile** (`TUMOR_REFERENCE_AREA_UM2
= 144 µm²`, chosen as the production `dx²` so the change is a pure
re-parameterisation *at production* and removes the spurious `dx²` factor
everywhere else). The conserved quantity is the integral: `sum(rho_tum) * dx² /
144` is the tumour cell count, which is a property of the tissue and cannot
depend on the mesh.

It works. `results/fstate_mechanism.md` §3 measures the integrated tumour sink at
**23812.00000000 at every rung, conserved to 2.9e-11**, and attributes the
sink term only the small `B−E` residual.

**The source is the same shape of quantity — an integral over tissue that
describes the tissue — and it received none of this.**

---

## 1. The diagnosis, reproduced independently before anything was changed

`vessel_network.discrete_vascular_field(geometry='binary')`:

```python
hit = tree.query_ball_point(grid_pts, r=float(radius_um), return_length=True) > 0
field = hit.reshape(ny, nx).astype(np.float64)
```

A grid cell is 1.0 iff its **centre** lies within 7.5 µm of an endothelial cell,
and 0.0 otherwise. The painted vessel area is therefore `(number of cells whose
centre won that lottery) × dx²`. A 7.5 µm disc is **176.71 µm²**; a `dx = 30 µm`
cell is 900 µm². Re-measured here on the real GBM section, `--stage area`:

| `dx` (µm) | cell area (µm²) | disc / cell | **binary** painted area (µm²) | published (`fstate_mechanism.md` §3.1) |
|---:|---:|---:|---:|---:|
| 30 | 900 | **0.196** | **713,700** | 713,700 |
| 24 | 576 | 0.307 | 748,800 | 748,800 |
| 18 | 324 | 0.545 | 746,496 | 746,496 |
| 16 | 256 | 0.690 | 745,728 | 745,728 |
| 12 | 144 | **1.227** | **748,656** | 748,656 |
| 10 | 100 | 1.767 | 739,200 | 739,200 |
| 8 | 64 | 2.761 | 746,816 | 746,816 |
| 7 | 49 | 3.606 | 741,517 | 741,517 |
| 6 | 36 | 4.909 | 743,904 | 743,904 |
| 5 | 25 | 7.069 | 744,000 | 744,000 |

> **Span 4.7312 % over the full ladder; 4.7307 % over the brief's own
> `dx = 30/18/12/8 µm` subset. Every row reproduces the published table
> exactly. THE DIAGNOSIS IS CONFIRMED — proceed.**

The shipped reconstruction was asserted to BE this arm at every rung
(`shipped_is_binary: true` for all 10, i.e. `micro.fields.vascular_density`
equals `geometry='binary'`), so this is a measurement of production and not of a
harness.

The movement is non-monotone, which is the tell: a discretisation error
converges, a quantisation lottery does not. At `dx = 30 µm` an entire vessel is
one fifth of one grid cell, so it is rounded to 0 or 1 cells — missed, or painted
5.1x too large.

---

## 2. The fix: measure the disc instead of voting on it

`geometry='area'` (`vessel_network._area_coverage_field`) paints **the same
discs, at the same radius, around the same endothelial cells**, and gives each
grid cell the **fraction of its area** those discs cover. That is the source-side
form of the sink's move: a per-cell token becomes a per-unit-area density, and
the conserved quantity becomes the integral.

**How the `dx`-independence is obtained — by construction, not by convergence.**
Coverage is measured on a sub-pixel lattice of pitch `h = dx / n_sub`,
`n_sub = ceil(dx / SOURCE_SUBPIXEL_PITCH_UM)`. Sub-point `k` of grid cell `j` sits
at `j*dx + (k+0.5)*h`, i.e. at `(m + 0.5)*h` with `m = j*n_sub + k` — **a single
global lattice anchored at the domain origin.** The grid decides only how those
points are bucketed into cells, never where they are. `SOURCE_SUBPIXEL_PITCH_UM =
1.0 µm` is chosen so that **every rung of the production ladder is an exact
integer multiple of it**, in which case `h == 1.0 µm` at every rung and the
lattice is literally the same point set. Discs around neighbouring endothelial
cells overlap, so the painted set is taken as a UNION (`np.unique` on the packed
sub-pixel index), never a sum.

| `dx` (µm) | binary (µm²) | **area (µm²)** | sub-grid (µm², context) |
|---:|---:|---:|---:|
| 30 | 713,700 | **742,713.0** | 347,327 |
| 24 | 748,800 | **742,713.0** | 347,349 |
| 18 | 746,496 | **742,713.0** | 347,499 |
| 16 | 745,728 | **742,713.0** | 347,468 |
| 12 | 748,656 | **742,713.0** | 347,340 |
| 10 | 739,200 | **742,713.0** | 347,230 |
| 8 | 746,816 | **742,713.0** | 347,381 |
| 7 | 741,517 | **742,713.0** | 347,481 |
| 6 | 743,904 | **742,713.0** | 347,396 |
| 5 | 744,000 | **742,713.0** | 347,307 |

> **binary 4.7312 % — area 6.3e-14 % — sub-grid 0.0774 %.**
> The fix is not "61x more stable" like the sub-grid painting; on the production
> ladder it is **exact**.

### 2.1 It is exact on the ladder, and still 65x better off it

The exactness above is bought by the pitch dividing `dx`. That is a real
condition, so it was tested rather than assumed, on rungs deliberately chosen NOT
to divide 1.0 µm:

| `dx` (µm) | binary (µm²) | area (µm²) |
|---:|---:|---:|
| 45 | 710,775 | 741,945 |
| 22.5 | 773,044 | 742,684 |
| 11.25 | 739,758 | 742,907 |
| 7.5 | 742,500 | 742,907 |
| 3.75 | 744,370 | 742,907 |

> **Off-lattice span: binary 8.391 %, area 0.1295 % — 65x better.** Off the
> lattice the guarantee degrades from "identical" to a Riemann error of order
> `perimeter x h`, which is what a fixed-pitch estimator should give and is still
> two orders below the defect.

### 2.2 The absolute value is the true disc area, and the estimator's own error is measured

The conserved number has to be the physically right one, or the fix is a
re-scaling. Sweeping the sub-pixel pitch on the real section (identical results
at `dx = 12` and `dx = 30`, which is the global-lattice property again):

| pitch `h` (µm) | painted area (µm²) at `dx = 12` | at `dx = 30` |
|---:|---:|---:|
| 4.0 | 742,512.0 | 744,370.3 |
| 2.0 | 742,480.0 | 742,480.0 |
| **1.0 (shipped)** | **742,713.0** | **742,713.0** |
| 0.5 | 742,605.0 | 742,605.0 |
| 0.25 | 742,643.8 | 742,643.8 |

The shipped `h = 1.0 µm` sits **+0.0093 %** from the `h = 0.25 µm` value, i.e.
the estimator's own error is ~120x smaller than the 1.12 % worth of movement it
removes. (On an isolated disc the same estimator returns 177.0000 µm² against an
exact 176.7146 µm², +0.16 %, at every `dx` — `tests/test_vessel_source.py`.)
`h = 4.0 µm` no longer divides `dx = 30`, which is why that one cell of the table
breaks ranks: it is the §2.1 effect, visible in miniature.

---

## 3. What the fix does to `F_state`

`--stage fstate` re-solves the shipped linear oxygen PDE at every rung under each
source arm and bands the result exactly as `cells.agents.update_cell_states`
does, reusing `experiments/fstate_mechanism.py`'s own solver, sampler and banding
helpers so these numbers are comparable to that report rather than merely similar
to it.

| `dx` (µm) | **binary** `F_state` | published §3 | **area** `F_state` | **area_matched** |
|---:|---:|---:|---:|---:|
| 30 | 0.5041 | 0.5041 | 0.5424 | 0.5124 |
| 24 | 0.5392 | 0.5392 | 0.5425 | 0.5474 |
| 18 | 0.5401 | 0.5401 | 0.5422 | 0.5455 |
| 16 | 0.5333 | 0.5333 | 0.5421 | 0.5449 |
| **12** | **0.5502** | 0.5502 | **0.5422** | 0.5471 |
| 10 | 0.5305 | 0.5305 | 0.5420 | 0.5387 |
| 8 | 0.5488 | 0.5488 | 0.5420 | 0.5454 |
| 7 | 0.5407 | 0.5407 | 0.5419 | 0.5407 |
| 6 | 0.5390 | 0.5390 | 0.5420 | 0.5428 |
| 5 | 0.5435 | — | 0.5419 | 0.5427 |

| arm | `F_state` span | mean-O2 span | quiescent span | source-area span | **drift removed** |
|---|---:|---:|---:|---:|---:|
| `binary` (shipped) | **0.04609** | 0.012269 | 0.09218 | 4.7312 % | — |
| **`area` (the fix)** | **0.000546** | **0.000117** | 0.001092 | **0.0000 %** | **98.82 %** |
| `area_matched` | 0.034982 | 0.010649 | 0.069965 | 4.7312 % | 24.10 % |

The `binary` arm reproduces `results/fstate_mechanism.md` to its published
digits, including the whole causal chain, so nothing below is the harness:

| quantity | this run | published |
|---|---:|---:|
| `F_state` span | 0.04609 | 0.0461 |
| mean tumour O2 span | 0.012269 | 0.0123 |
| `corr(painted area, mean O2)` | **0.9868** | 0.9868 |
| `corr(mean O2, F_state)` | **0.9357** | 0.9357 |
| `dq/d(mean O2)` | **−6.72** | −6.72 |
| `F_state` at `dx = 12 µm` | 0.5502 | 0.5502 |

**A corroboration worth naming.** `fstate_mechanism.md` §3 defines `D = 0.5437`
as "what every arm would return if nothing depended on the mesh" — the reference
grid with the reference source. The fixed source lands the entire ladder at
**0.5419–0.5425**, i.e. within **0.0015** of that mesh-independent baseline,
having been derived by a completely different route.

---

## 4. THE HARD PART: does the fix change the total, or only its distribution?

It changes the total. **−0.794 % at the production spacing.** And it could not
have done otherwise.

| | painted area (µm²) |
|---|---:|
| shipped, at production `dx = 12 µm` | 748,656 |
| shipped, ladder minimum (`dx = 30`) | 713,700 |
| shipped, ladder maximum (`dx = 24`) | 748,800 |
| shipped, ladder mean | 741,882 |
| **conserved (`geometry='area'`), every rung** | **742,713** |

* against the production rung: **−0.794 %**
* against the ladder mean: **+0.112 %**
* against the finest shipped rung (`dx = 5`): **−0.173 %**

> **"Conserve the shipped total" is not a well-posed constraint here.** The
> shipped total is a *different number at every rung*, spanning 4.73 %. Any
> single conserved value necessarily disagrees with nine of the ten rungs. The
> −0.794 % is not a design choice inside the fix; it is the distance from the
> production rung to the value that has no `dx` in it, and the production rung
> happens to sit near the TOP of the quantisation lottery (748,656 against a
> ladder mean of 741,882).

### 4.1 The `area_matched` arm settles it by experiment

`match_total_to` rescales the fixed field so its tissue-integrated total is
EXACTLY the shipped per-rung total. That is the strictest possible reading of
"conserve the shipped total while removing the grid dependence" — and it is
measured, not argued:

> **It removes 24.10 % of the drift, against the raw fix's 98.82 %.**

The reason is structural. Rescaling to a shipped total that itself spans 4.73 %
re-imports the entire defect through the gauge: `area_matched`'s source-area span
is 4.7312 %, identical to `binary`'s by construction, and its
`corr(mean O2, F_state)` is 0.9993. At `dx = 30 µm` it drops to `F_state =
0.5124` purely because the binary total it is being matched to is 3.9 % low
there.

> **Conserving the total and removing the drift are in direct conflict, because
> THE TOTAL IS THE DEFECT.** A fix that preserved the shipped number at every
> rung would preserve the thing that carries `cov/var = +1.057` of the variance.
> That is the honest negative half of this result, and it is why the fix ships
> default-off.

### 4.2 Where the production-rung change actually goes

At `dx = 12 µm`, decomposing the total `F_state` move:

| | `F_state` | Δ vs shipped |
|---|---:|---:|
| `binary` (shipped) | 0.5502 | — |
| `area_matched` (distribution only, total held exactly) | 0.5471 | **−0.0031** |
| `area` (the fix) | 0.5422 | **−0.0080** |

**39 % of the production-rung change is placement; 61 % is the −0.794 % of
dose.** Both are small next to the 0.0461 of grid artefact removed.

---

## 5. The other half of the drift is untouched, and is not fixable by anything

`results/fstate_mechanism.md` §4 measured the amplifier: a hard band edge cutting
through the densest part of a population that spans only 0.13 in O2. That is not
a defect and no rasterisation can address it. Measured on the tumour population
at `dx = 12 µm` under both sources (`results/sweep/vsrc_band_edge.json`):

| statistic | **binary** (shipped) | **area** (fixed) |
|---|---:|---:|
| mean tumour O2 | 0.34863 | 0.34728 |
| min / max | 0.2853 / 0.4289 | 0.2905 / 0.4214 |
| **full population range** | **0.1436** | **0.1309** |
| p90 / p95 | 0.4001 / 0.4075 | 0.3968 / 0.4042 |
| fraction below the 0.4 ceiling | 0.8996 | 0.9157 |
| fraction below the 0.2 hypoxia edge | 0.0000 | 0.0000 |
| within 0.005 of the ceiling | 4.86 % | 6.13 % |
| **within 0.01 of the ceiling** | **9.87 %** | **9.78 %** |
| within 0.02 of the ceiling | 15.63 % | 15.52 % |
| **local density at the ceiling (per unit O2)** | **4.94** | **4.89** |

> **The amplifier is intact.** At the ±0.01 and ±0.02 scales the pile-up is
> unchanged to within 1 % relative (9.87 → 9.78 %, 15.63 → 15.52 %) and the
> ceiling density moves 4.94 → 4.89 per unit O2. At tighter windows it is
> slightly WORSE under the fix (±0.005: 4.86 → 6.13 %; ±0.001: 0.84 → 1.15 %),
> because the fixed field is marginally narrower (range 0.1436 → 0.1309) and
> packs the same population into less O2. The fix removes the *forcing* and
> removes **none** of the *gain* — exactly as `fstate_mechanism.md` §5 predicted
> it would, and if anything the gain is very slightly higher.

**So how is the 0.0461 split?** Not into two additive halves — that framing does
not survive contact with the measurement, and saying otherwise would be
inventing a decomposition. The drift is a product:

```
F_state drift  =  (field movement)  x  (band-edge gain)
     0.04609   =      0.012269      x        3.76
```

(3.76 is the measured span ratio; the regressed `dq/d(mean O2) = −6.72` implies
3.36 into `F_state`, and the two agree to the accuracy the non-uniformity of the
shift allows — the same 89 %-closure `fstate_mechanism.md` §4 reports.)

The fix attacks the left factor and removes **99.05 %** of it (mean-O2 span
0.012269 → 0.000117). The right factor is unchanged: 9.78 % of tumour cells
remain within 0.01 of the ceiling, the ceiling density is still ~4.9 per unit O2,
and the population still spans 0.13. **The residual 0.000546 of `F_state` span is
not a leftover piece of the source defect — it is the same undiminished gain
acting on the 0.000117 of movement that the PDE discretisation and the
cell-to-grid sampler still contribute** (`fstate_mechanism.md` §3 priced those at
`cov/var = −0.056` and `−0.001`).

Concretely: **98.82 % of the 0.0461 is removed; the 1.18 % that remains is the
band-edge amplifier operating on an irreducible residual field movement, and
refining anything further will not shrink it below what the amplifier does to
whatever movement is left.** The response to the amplifier is still the one
`fstate_mechanism.md` §5 named — stop using hard bands — and **nothing here
implements that.**

---

## 6. What the fix does NOT fix (read before promoting it)

**The painted area is conserved. The painted SUPPORT is not, and cannot be.**

| `dx` (µm) | binary support (frac. of tissue) | **area support** |
|---:|---:|---:|
| 30 | 0.0256 | **0.2105** |
| 18 | 0.0268 | 0.1216 |
| 12 | 0.0268 | 0.0838 |
| 8 | 0.0268 | 0.0613 |
| 5 | 0.0267 | 0.0454 |

A grid cell that contains any part of a vessel is nonzero, so at `dx = 30 µm` the
source is spread — at low amplitude — over 21 % of the tissue against 4.5 % at
`dx = 5 µm`. This is not a defect of the fix; it is what a cell-averaged
(finite-volume) representation of a sub-grid object *means*, and it is the
honest form of the information the binary field destroyed by rounding. But it
does mean the fix conserves the source's **integral**, not its **placement**,
and any distance-resolved claim about perivascular gradients still needs
`dx <= 8 µm` exactly as `results/vessel_promotion.md` requires today. At `dx =
12 µm` only **0.0067 %** of tissue grid cells reach a value of 1.0 — the field is
essentially entirely partial, which is the point.

**Two knock-on consequences, named so they are not discovered later:**

1. **`results/subgrid_closure.md`'s no-op result is a property of the BINARY
   field, not of the closure.** That report found `beta` identically 1 because
   the field is exactly 0 or 1 — "the binary painting destroys the sub-grid
   geometry before the solve runs". Under `geometry='area'` the field is
   continuous in [0, 1], so the sub-grid closure is no longer identically a
   no-op. **What it then does is unmeasured. This report does not run it**, and
   nothing in `tmesim/physics/subgrid.py` was touched.
2. **`results/oxygen_supply_calibration.md` §1 derives `supply`'s units from the
   premise that `vasc` is "exactly 1.0 inside a vessel pixel and exactly 0
   outside".** Under the fix `vasc` becomes the vessel AREA FRACTION of the grid
   cell, which is the quantity a per-unit-tissue-volume transvascular transfer
   coefficient should actually multiply — so the fix makes that derivation *more*
   dimensionally honest, not less. **`OxygenConfig.supply` was not touched and
   that derivation was not re-run.**

---

## 7. The nine published verdicts

`--stage verdicts` re-scores the nine published 3-D:2-D rows through
`experiments/delivery_source_geometry.py`'s existing blast harness — the same
`_run_panel` / `kill_calibration.score_ratio_row` path
`results/delivery_source_geometry.md` used — so these arms are directly
comparable to that report's `binary` / `subgrid` arms. The `binary` control is
**replayed from that report's own checkpoint** rather than re-solved, so the
control cannot drift: it is literally the published one.
`experiments/delivery_source_geometry.py` was NOT edited; `_source_field` is
patched in-process to teach it the two new arms and every other arm delegates
straight back to the original.

> **ZERO of the nine rows change verdict. Zero of the seven SCORED rows flip.**
> The largest ratio move anywhere is **2.62 %**. For contrast, the same harness
> measured `geometry='subgrid'` flipping **three of seven** — GBM temozolomide
> PASS→FAIL, GBM cisplatin PASS→FAIL, breast docetaxel FAIL→PASS.

| row | published ratio | published verdict | **`area` ratio** | Δ | **`area` verdict** | `area_matched` | `subgrid` (for contrast) |
|---|---:|---|---:|---:|---|---|---|
| GBM temozolomide | 0.2626 | PASS | **0.2581** | −1.73 % | **PASS** | 0.2607 PASS | 0.1550 **FAIL** |
| GBM cisplatin | 0.2791 | PASS | **0.2718** | −2.62 % | **PASS** | 0.2745 PASS | 0.1042 **FAIL** |
| GBM doxorubicin | 0.2077 | PASS | **0.2076** | −0.06 % | **PASS** | 0.2096 PASS | 0.0660 PASS |
| breast doxorubicin | 0.2301 | PASS | **0.2303** | +0.07 % | **PASS** | 0.2304 PASS | 0.1307 PASS |
| breast paclitaxel | 0.1503 | NOT TESTABLE | **0.1534** | +2.07 % | NOT TESTABLE | 0.1535 | 0.0763 |
| breast docetaxel | 0.1484 | **FAIL** | **0.1516** | +2.18 % | **FAIL** | 0.1518 FAIL | 0.0752 **PASS** |
| breast epirubicin | 0.2301 | PASS | **0.2303** | +0.07 % | **PASS** | 0.2304 PASS | 0.1307 PASS |
| breast 5-fluorouracil | 0.3379 | NOT TESTABLE | **0.3358** | −0.61 % | NOT TESTABLE | 0.3360 | 0.2530 |
| breast cisplatin | 0.3085 | PASS | **0.3051** | −1.10 % | **PASS** | 0.3053 PASS | 0.2160 PASS |

The `binary` control reproduces all nine published ratios to `< 5e-5`
(`binary_control_reproduces_published: true`), and `calibrate_alpha()` returned
148.60815765445545 inside the run.

**Three things this table establishes.**

1. **The project's one falsification survives.** `breast docetaxel` is the only
   FAIL in the published record, and `results/delivery_source_geometry.md`
   rejected the sub-grid geometry partly because it *removed* it. Under this fix
   it stays FAIL, and its margin moves the *safe* way (0.1484 → 0.1516 against a
   band it is still outside).
2. **The flagship survives.** GBM temozolomide — the only drug with a matched
   patient-derived ex-vivo anchor — stays PASS at 0.2581 against a published
   0.2626.
3. **The fix is a different kind of perturbation from the sub-grid geometry, and
   the ~2 % scale is why.** Both change the source. `subgrid` changes the
   geometric object, the weighting law and the total (744k → 347k µm²) and moves
   the nine ratios to **0.32–0.75x** of published. `area` changes only the
   rasterisation of the shipped object and moves them by at most 2.6 %.

**Source and state under each arm at `dx = 12 µm`:**

| panel | source Σ over tissue | support | mean O2 | `F_state` |
|---|---:|---:|---:|---:|
| GBM `binary` | 5199.0 | 0.0268 | 0.35205 | 0.5502 |
| GBM `area` | 5157.7 (**−0.79 %**) | 0.0838 | 0.35066 | 0.5422 |
| GBM `area_matched` | 5199.0 (exact) | 0.0838 | 0.35245 | 0.5471 |
| breast `binary` | 15191.0 | 0.0567 | 0.34572 | 0.5063 |
| breast `area` | 15164.5 (**−0.17 %**) | 0.1533 | 0.34532 | 0.5059 |
| breast `area_matched` | 15191.0 (exact) | 0.1533 | 0.34569 | 0.5063 |

The dose change is **section-dependent** — −0.79 % on GBM, −0.17 % on breast —
which is exactly what it should be, since it is the distance from each section's
own `dx = 12 µm` quantisation draw to its own conserved value, and those draws
are independent. Breast `F_state` barely moves (0.5063 → 0.5059), consistent with
its much smaller dose change; the breast section's own band-edge distribution was
not measured here, so no mechanism is claimed for it.

---

## 8. Why the ~0.8 % dose change is not the reason to hesitate

`results/oxygen_supply_calibration.md` establishes that `OxygenConfig.supply` is
**4,570x too small**, and that it is a *compensating* error against a sink
~140x too small. The source term the PDE actually sees is `supply * vasc`, and
`vasc`'s tissue integral is exactly the quantity this fix conserves. So:

> **This fix moves `supply x ∫vasc` by a factor of 0.99206. That same product
> carries a standing, documented, sourced miscalibration of a factor of 4,570.**
> On a log scale the perturbation introduced here is ~1,060x smaller than the
> error already known to be sitting in the same term.

That is not an argument that 0.794 % is negligible — it changes `F_state` at
production by 0.0049 and it is reported as a change, not smuggled. It is an
argument that **the `supply` question this fix "re-opens" was never closed**, and
that declining a fix which removes 98.8 % of a grid artefact in order to protect
the eighth-of-a-percent stability of a constant known to be wrong by three and a
half orders of magnitude would be the wrong trade. What must NOT happen is
`supply` being quietly re-tuned to absorb the 0.794 %; that would recreate
exactly the compensating-error structure that report documents. **`supply` is
untouched here and the recommendation is that it stay untouched until it is
recalibrated on its own terms.**

---

## 9. Preconditions and hygiene

Every constant this work was forbidden to move, checked rather than claimed
(`results/sweep/vsrc_hygiene.json`):

| | value | required |
|---|---|---|
| `calibrate_alpha()` | **148.60815765445545** | unchanged ✓ |
| `OxygenConfig.supply` | 0.05 | unchanged ✓ |
| `OxygenConfig.consumption` | 0.02 | unchanged ✓ |
| `OxygenConfig.hypoxia_threshold` | 0.2 | unchanged ✓ |
| QUIESCENT ceiling | 0.4 | unchanged ✓ |
| `VESSEL_RADIUS_UM` | 7.5 | unchanged ✓ |
| `LINK_DISTANCE_UM` | 20.0 | unchanged ✓ |
| `real_fields._DEFAULT_VASCULAR_MODE` | `discrete` | unchanged ✓ |
| `discrete_vascular_field(geometry=)` default | **`binary`** | unchanged ✓ |
| shipped field == `geometry='binary'` | true at all 10 rungs | ✓ |

`tmesim/physics/oxygen.py` was read and **not modified**. Nothing under
`tmesim/molecular/`, `tmesim/validation/`, `app/` or `landing-3d/` was touched,
and no existing `results/` file was edited.

**Test suite.** `./.venv/bin/python -m pytest -q` → **950 passed, 0 failed**
(baseline before this work, same command, same machine: **928 passed, 0
failed**). This work adds the **13** tests of `tests/test_vessel_source.py`; the
other 9 arrive from concurrent work under `tmesim/molecular/` by another agent
(`tests/test_two_route.py`, `tmesim/molecular/lipophilicity.py`) and are not
this change's. No pre-existing test was modified, and the known ordering flake
`test_app_server.py::test_failed_runs_surface_their_error` did not fire in either
run.

**What changed under `tmesim/`, in full:** `reconstruct/vessel_network.py` gains
`SOURCE_SUBPIXEL_PITCH_UM`, `_area_coverage_field`, a `subpixel_pitch_um`
keyword, and an `'area'` case in `discrete_vascular_field`'s geometry dispatch.
The pre-existing binary rasterisation is byte-for-byte the same code, now guarded
by `if geom_mode == "binary":`, and `tests/test_vessel_source.py` pins that the
default is still `'binary'`, still binary-valued, and still exactly what
`vascular_density_from_cells` returns.

---

## 10. Merge recommendation

**LAND AS SHIPPED — opt-in, default-off.** The branch is additive, the default is
bit-identical, and the pinned nine verdicts are unmoved because nothing routes
through it.

**Then run the promotion as a deliberate, separately-scored recalibration**, not
as a bug fix, because it moves the production source total by −0.794 % and
production `F_state` by −0.0080. §4.1 is the reason this cannot be finessed: the
only version of the fix that preserves the shipped total is `area_matched`, and
that one recovers just 24 % of the drift because it re-imports the defect through
the gauge. The promotion re-run needs exactly three things, none of which is done
here:

1. the published `results/kill_calibration.md` / `p99_convention.md` numbers
   REGENERATED under `geometry='area'`, not merely re-scored. §7 already
   establishes the decisive part — no verdict moves and the largest ratio move is
   2.62 % — so this is a bookkeeping re-run rather than a risk, but the published
   ratios in those files would still be stale by up to 2.62 % until it is done;
2. `results/subgrid_closure.md` re-run, because §6 shows its no-op result was a
   property of the binary field and is void under a continuous one;
3. `OxygenConfig.supply` left ALONE, so the 0.794 % is visible as a change rather
   than absorbed into a constant already known to be 4,570x wrong (§8).

**Do not conflate this with `geometry='subgrid'`.** That painting deposits
347k µm² against 744k, weights vessels by a `k²` cross-section law instead of the
wall-area law, and `results/delivery_source_geometry.md` rejected it for the
delivery source on physical grounds after flipping three of seven scored
verdicts. `geometry='area'` keeps the shipped object, the shipped radius and the
shipped per-endothelial-cell weighting, and fixes only the rasterisation. They
are not alternatives to each other and the argument against one is not an
argument against the other.
