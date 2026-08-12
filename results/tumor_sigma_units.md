# The tumour-occupancy smoothing length was in pixels, so refining the mesh changed the model

**Status: landed.** `tmesim/physics/pressure.py` now expresses the tumour-occupancy
smoothing length in micrometres. Four bitwise reference values in
`tests/data/static_source_reference.json` were re-pinned as a consequence; the
superseded values are kept in the same file. The production grid `dx = 12 um` is
bit-identical and gained a bitwise pin it did not previously have.

Suite: **772 -> 775 passed** (three new tests, no test removed or weakened).

---

## 1. The defect

`_tumor_field` rasterises tumour cells and smooths them into a `[0, 1]`
occupancy field `phi`. That field is not cosmetic — it sets both terms that make
a tumour a tumour in this model:

```python
lp_eff = _LP_NORMAL * (1.0 + _LEAK_GAIN * phi)   # leaky tumour vessels
d      = _LAMBDA_LYMPH * (1.0 - phi)             # lymphatic drainage -> 0 in tumour
```

So the length over which occupancy is smeared is a **model** parameter: it is
the width of the leaky-and-undrained region, i.e. the tumour margin. It was
written as

```python
_TUMOR_SIGMA = 5.0   # gaussian smoothing (grid cells) for the tumor-density field
...
phi = ndimage.gaussian_filter(phi, sigma=_TUMOR_SIGMA, mode="nearest")
```

`scipy.ndimage.gaussian_filter` takes `sigma` in **pixels**. The comment is
therefore accurate about the units and wrong about the consequence: a model
parameter had been expressed in mesh units, so its physical value was whatever
the current grid happened to be.

| `dx` (um) | 8 | **12 (production)** | 15 | 18 | 24 | 80 |
|---|---:|---:|---:|---:|---:|---:|
| physical smoothing, before (um) | 40 | **60** | 75 | 90 | 120 | 400 |

**Why this matters more than a factor-of-two in a constant.** Refining the mesh
silently redefined the model. Every grid-convergence study on an IFP quantity in
this repo was therefore measuring two things superimposed — discretisation error
*and* a change of physics — and could not converge even in principle. The defect
undermines the method used to check other defects, which is why it is worth a
document of its own.

Originally reported in `results/ifp_validation.md` §5.3; diagnosed, implemented
and then deliberately reverted in `results/small_fixes.md` §1, because landing it
requires re-pinning a bitwise guard. This document is that landing.

## 2. The fix

Two lines. Make the constant a physical length; convert at the point of use.

```diff
-_TUMOR_SIGMA = 5.0   # gaussian smoothing (grid cells) for the tumor-density field
+_TUMOR_SIGMA_UM = 60.0  # gaussian smoothing (MICROMETRES) for the tumor-density field
...
-    phi = ndimage.gaussian_filter(phi, sigma=_TUMOR_SIGMA, mode="nearest")
+    sigma_px = _TUMOR_SIGMA_UM / float(micro.domain.dx)
+    phi = ndimage.gaussian_filter(phi, sigma=sigma_px, mode="nearest")
```

**`60.0` is not a new number.** It is what the old constant already produced at
the production grid: `5.0 cells x 12 um/cell = 60.0 um`. Nothing was calibrated,
tuned or fitted; no other constant was touched. The change is a units
correction, and it is deliberately anchored so that the grid everything ships on
does not move.

`dx = 12` is consequently the exact fixed point:

| `dx` (um) | 8 | **12** | 15 | 18 | 24 | 80 |
|---|---:|---:|---:|---:|---:|---:|
| sigma before (px) | 5.0 | **5.0** | 5.0 | 5.0 | 5.0 | 5.0 |
| sigma after (px) = 60/dx | 7.5 | **5.0** | 4.0 | 3.333 | 2.5 | 0.75 |
| physical, before (um) | 40 | **60** | 75 | 90 | 120 | 400 |
| physical, after (um) | 60 | **60** | 60 | 60 | 60 | 60 |

## 3. What it does and does not fix

Solving IFP on the real Xenium GBM section at four grids, with the old
pixel-space constant reproduced for comparison:

| `dx` (um) | physical sigma, before | mean IFP in tissue, before | physical sigma, after | mean IFP in tissue, after |
|---:|---:|---:|---:|---:|
| 8 | 40 um | 0.096335 | 60 um | 0.121386 |
| 12 | 60 um | **0.128387** | 60 um | **0.128387** |
| 16 | 80 um | 0.150662 | 60 um | 0.127746 |
| 24 | 120 um | 0.181899 | 60 um | 0.132363 |
| | | **1.89x spread** | | **1.09x spread** |

Before, the mean rose monotonically with `dx` — which is the signature of a
model changing with the mesh, not of a solver converging. After, it sits in a
narrow band and the remaining variation is discretisation error, which is the
thing a convergence study is supposed to measure.

**This does not make the IFP field grid-converged, and nothing here should be
read as claiming it does.** The peak still moves substantially with `dx` (11.74
-> 18.34 mmHg across the same ladder). That is consistent with what
`landing-3d/public/models/ifp_contours.json` already states in its
`resolution_caveat`: the screening length is ~7.1 um, below the 12 um grid at
99.8 % of tissue nodes, so a peak statistic on this field is not resolved. What
the fix changes is that convergence is now a *measurable* property rather than a
confounded one.

## 4. Production is bit-identical — verified, not asserted

`60.0 / 12.0 == 5.0` is exact in IEEE-754 (`0x1.4000000000000p+2` both ways), so
at `dx = 12` `gaussian_filter` receives a bit-identical argument.

Checked at four levels, before and after the change:

* **Unit.** `np.array_equal(gaussian_filter(phi, 5.0), gaussian_filter(phi, 60.0/12.0))`
  -> `True`.
* **Field.** On the real Xenium GBM section at `dx = 12` via the production
  reconstruction call, SHA-256 of the solved pressure field, of the Darcy
  velocity, and of the 250 um display coarse-grain are unchanged
  (`7969ca5d...`, `5969c952...`, `d5be310c...`).
* **Published figures.** `raw_field_max_mmhg = 12.93283`,
  `raw_field_mean_in_tissue_mmhg = 0.128387`, `displayed_field_max_mmhg = 0.72608`
  — all reproduced exactly, and all **5** exported iso-contours hash identically
  vertex for vertex.
* **End to end.** A full engine run on that section with convection ON gives an
  identical `kill_fraction`, `mean_auc`, `auc_sum`, `final_c_sum`,
  `penetration_index` and `viable_history`, to the last bit.

No shipped result, figure or landing asset changes.

## 5. The re-pin, and why it is legitimate

With the fix in place, four bitwise pins in `tests/test_time_resolved_pk.py`
failed against `tests/data/static_source_reference.json`. They have been
regenerated. Because a bitwise guard is only worth anything if re-pinning is
visible, the change is recorded three times: here, in the test module's
docstring, and in the reference file's own `_README._repin_history`, with the
superseded values kept verbatim under
`_superseded_2026-08-12_grid_cell_tumor_sigma`.

**The failure pattern confirms the mechanism rather than merely correlating with
it.** `phi` reaches a transport result by exactly one route — the Darcy velocity
`v = -K grad p` — so it can bite only when convection is enabled:

| case | `dx` (um) | convection | outcome |
|---|---:|---|---|
| `default_1dose` | 15 | **on** | re-pinned |
| `weak_1dose` | 15 | **on** | re-pinned |
| `weak_3dose` | 15 | **on** | re-pinned |
| `real_xenium_dx80_tmz_3dose` | 80 | **on** | re-pinned |
| `weak_10dose_noconv` | 15 | off | **byte-for-byte unchanged** |
| `solver_flat` (raw solver, no IFP) | 15 | off | **byte-for-byte unchanged** |

Every convection-on case moved; **every convection-off case did not**. This was
predicted before the values were regenerated and was the stated stop condition:
had a convection-off case or a `dx = 12` value moved, the mechanism would not
have been what we thought and the fix would not have landed. Neither did.

The magnitudes are also consistent with a small change to the velocity field
rather than with something structural: the transport sums moved by 1-3 ULP
(e.g. `auc_sum` for `default_1dose`, `0x1.9127422a88cfbp+38` ->
`0x1.9127422a88cf9p+38`), and `plasma_mean` — a pure-PK quantity with no spatial
dependence — is unchanged in all four. `penetration_index` moved most, as a
shape statistic should. `kill_fraction` is unchanged in the three hand-built
cases (it is quantised at 1/300 cells there) and moved at `dx = 80`.

**Why the new values are right and not merely different.** The old values encode
a model whose smoothing length was an artefact of the grid it happened to be run
on: 75 um at `dx = 15`, 400 um at `dx = 80`. Both now use 60 um — the same
physical length the production grid always used and the one every published
result was computed under. The new numbers are the model these cases were always
meant to be running, evaluated on their own meshes.

## 6. The guard is stronger than it was

The old reference pinned `dx = 15` and `dx = 80`. Nothing ships on either. The
production grid `dx = 12` — the grid behind every number in `results/` and behind
the landing asset — had **no bitwise pin at all**. So the grids that were pinned
were exactly the ones whose values were wrong, and the grid that mattered was
unguarded.

Three tests were added in the same commit that re-pins:

1. `test_production_dx12_ifp_solve_matches_stored_reference_bitwise` — pins the
   IFP solve itself on the real section at `dx = 12` (max and mean pressure in
   tissue, max Darcy speed). This is the level the defect acts at, so a
   regression shows up without needing a transport solve to reveal it.
2. `test_production_dx12_static_path_matches_stored_reference_bitwise` — the full
   engine on that section, convection ON, pinning `kill_fraction`, `mean_auc`,
   `plasma_mean`, `penetration_index`, `auc_sum`, `final_c_sum` and the whole
   `viable_history`.
3. `test_tumor_occupancy_smoothing_is_a_physical_length` — asserts that
   `_TUMOR_SIGMA` no longer exists, that the pixel sigma scales as `1/dx`, and
   that `_TUMOR_SIGMA_UM / 12.0 == 5.0` exactly. It asserts the *conversion*,
   not the value, so a future deliberate recalibration of the 60 um remains
   possible; a silent return to grid-cell units does not.

The values behind (1) and (2) were generated **before** the fix was applied and
verified byte-identical **after** it. The production pin therefore does not
merely record current behaviour — it pins the invariance claim this whole change
rests on.

## 7. Known loose end (not owned here)

`tmesim/validation/scorecard.py` still describes the defect as live, in prose,
at two places (lines ~239 and ~2127: "`_TUMOR_SIGMA = 5.0` is specified in GRID
CELLS"). Those strings are now stale. They are documentation, not behaviour —
nothing imports the constant outside `pressure.py` — but the scorecard's owner
should update them. Similarly `results/ifp_validation.md` §5.3 and
`results/small_fixes.md` §1 describe the defect as unfixed; both are historical
records of the diagnosis and are left as written, with this document as the
landing record.
