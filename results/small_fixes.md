# Three small independent fixes: two landed, one blocked

Three defects found by other agents this session, each cheap, each verified
independently before it was acted on. All three are real. **Two landed. The
third is correct, is implemented, and is reverted** — because landing it
requires re-pinning a bitwise regression reference that is not mine to re-pin,
and forcing it past that guard would have been the wrong call.

Test suite: **754 passed**, matching the baseline.

| # | Defect | Confirmed independently? | Status | Shipped `dx = 12` behaviour |
|---|---|---|---|---|
| 1 | `_TUMOR_SIGMA` is in grid cells, so the model changes with the mesh | **yes**, as described | **reverted — blocked, see §1** | would have been bitwise identical |
| 2 | Exported IFP metadata states a screening length ~8.5x too large | **yes**, two independent ways | **landed** | bitwise identical (every contour vertex) |
| 3 | Grid ladders discard the per-cell exposure vectors | yes for GBM; the breast half of the claim was wrong in its details (§3) | **landed** | unchanged — nothing re-run |

---

## 1. `_TUMOR_SIGMA` is a mesh-dependent *model* — confirmed, fixed, and then reverted

**File:** `tmesim/physics/pressure.py` — **left unmodified, at `HEAD`.**

### The defect is real
`results/ifp_validation.md` §5.3 reported it; I confirmed it directly. The
constant is fed straight to `scipy.ndimage.gaussian_filter`, whose `sigma` is in
**pixels**:

```python
_TUMOR_SIGMA = 5.0   # gaussian smoothing (grid cells) for the tumor-density field
...
phi = ndimage.gaussian_filter(phi, sigma=_TUMOR_SIGMA, mode="nearest")
```

The tumour-occupancy field `phi` sets both the vessel-leak gain
(`lp_eff = _LP_NORMAL * (1 + _LEAK_GAIN * phi)`) and the lymphatic suppression
(`d = _LAMBDA_LYMPH * (1 - phi)`), so this is a **model** parameter, not a
numerical one. Its physical value is currently:

| `dx` (um) | 8 | 12 (production) | 15 | 18 | 24 | 80 |
|---|---:|---:|---:|---:|---:|---:|
| physical smoothing (um) | 40 | 60 | 75 | 90 | 120 | 400 |

A grid-convergence study on any IFP quantity is therefore measuring two things
at once — refinement *and* a redefinition of the model — and can never converge.
This is the most serious of the three defects, because it silently undermines
the method used to check the others.

### The fix, which works
Two lines. Make the constant a physical length and convert at the point of use:

```diff
-_TUMOR_SIGMA = 5.0   # gaussian smoothing (grid cells) for the tumor-density field
+_TUMOR_SIGMA_UM = 60.0   # gaussian smoothing (MICROMETRES)
...
-    phi = ndimage.gaussian_filter(phi, sigma=_TUMOR_SIGMA, mode="nearest")
+    sigma_px = _TUMOR_SIGMA_UM / float(micro.domain.dx)
+    phi = ndimage.gaussian_filter(phi, sigma=sigma_px, mode="nearest")
```

`60.0` is derived from what the old constant produces at production resolution:
`5.0 cells x 12 um/cell = 60.0 um`. So `dx = 12` is the fixed point of the
change:

| `dx` (um) | 8 | **12 (production)** | 15 | 18 | 24 | 80 |
|---|---:|---:|---:|---:|---:|---:|
| old sigma (px) | 5.0 | **5.0** | 5.0 | 5.0 | 5.0 | 5.0 |
| new sigma (px) = 60/dx | 7.5 | **5.0** | 4.0 | 3.333 | 2.5 | 0.75 |
| physical sigma, old (um) | 40 | **60** | 75 | 90 | 120 | 400 |
| physical sigma, new (um) | 60 | **60** | 60 | 60 | 60 | 60 |

`60.0 / 12.0 == 5.0` is exact in IEEE double (checked, not assumed), so
`gaussian_filter` receives an identical argument and returns an identical array.
Verified three ways:

* **Unit level.** `np.array_equal(gaussian_filter(phi, 5.0),
  gaussian_filter(phi, 60.0/12.0))` → `True`.
* **End to end on real tissue.** With the fix applied, the full GBM landing
  export at `dx = 12` (§2) reproduced `raw_field_max_mmhg = 12.93283`,
  `raw_field_mean_in_tissue_mmhg = 0.128387` and
  `displayed_field_max_mmhg = 0.72608`, and **every contour vertex** matched the
  pre-existing file. Re-running the same export after reverting produced the
  same values again — so `dx = 12` really is invariant to this change, checked
  from both sides.
* **Self-test.** `python -m tmesim.physics.pressure` (which runs at `dx = 15`)
  passed all three Jain-effect assertions with the fix in place.

### Why it is reverted: the repo pins engine behaviour bitwise at `dx != 12`
With the fix applied the suite goes **754 → 750 passed, 4 failed**, all in
`tests/test_time_resolved_pk.py`:

```
FAILED test_static_path_matches_stored_reference_bitwise[default_1dose]
FAILED test_static_path_matches_stored_reference_bitwise[weak_1dose]
FAILED test_static_path_matches_stored_reference_bitwise[weak_3dose]
FAILED test_real_tissue_static_path_matches_stored_reference_bitwise
```

These compare against hex-exact values stored in
`tests/data/static_source_reference.json`.

**The attribution is certain, and the failure pattern confirms the mechanism
rather than merely correlating with it.** I restored the old pixel-space sigma
in place and re-ran the file: **16 passed**. The only route from `_tumor_field`
into a transport result is the Darcy velocity `v = -K grad p`, so the fix can
only bite when convection is enabled — and that is exactly the split observed:

| case | `dx` (um) | `convection` | result under the fix |
|---|---:|---|---|
| `default_1dose` | 15 | **True** | **FAIL** |
| `weak_1dose` | 15 | **True** | **FAIL** |
| `weak_3dose` | 15 | **True** | **FAIL** |
| `weak_10dose_noconv` | 15 | False | pass |
| `real_xenium_dx80_tmz_3dose` | 80 | **True** | **FAIL** |
| `solver_flat` (raw solver, no IFP) | 15 | False | pass |

Every convection-on case fails; every convection-off case passes. Nothing else
in the run moved.

**Why I did not just regenerate the reference.** The values are wrong under the
old model, so re-pinning them is the *correct* end state — but:

* the reference lives under `tests/`, which is not among the files I was given;
* it is a bitwise guard whose entire value is that nobody quietly regenerates it
  to make their own change pass. Regenerating it as part of the change it is
  catching destroys the guard and leaves no independent check that the new
  numbers are right rather than merely different;
* the brief is explicit — if a fix is riskier than described, do the other two
  and report rather than force it. The report I was given said the fix would
  preserve shipped `dx = 12` behaviour, which is true, but did not anticipate
  that the repo pins engine behaviour bitwise at `dx = 15` and `dx = 80`.

So `tmesim/physics/pressure.py` is left byte-identical to `HEAD`.

### To land it (about ten minutes, for whoever owns the reference)
1. Apply the two-line diff above.
2. Regenerate `tests/data/static_source_reference.json` for the five
   convection-enabled entries, and record in that file (or the test module) that
   the re-pin was caused by the units correction — otherwise the next reader
   cannot tell a deliberate re-pin from a silent regression.
3. Sanity-check the direction: at `dx = 15` the tumour-occupancy smoothing
   narrows 75 → 60 um and at `dx = 80` it narrows 400 → 60 um, so the margin
   pressure gradient should get *sharper*, not smoother, at both.
4. Nothing at `dx = 12` needs regenerating — that is the point of choosing 60.

### Also stale once it lands (files I do not own)
* `tmesim/validation/scorecard.py` (~L239, ~L1751) describes
  `_TUMOR_SIGMA = 5.0` as "specified in GRID CELLS" and lists it as a live
  defect. Correct today; must be updated with the fix.
* `results/ifp_validation.md` §5.3 and its parameter table (L146), likewise.

### What is *not* affected either way
No result on disk changes, because nothing at `dx = 12` changes. But every IFP
number computed at another resolution — the `dx = 24 / 18 / 8` rows of
`results/sweep/ifp_validation.json`, and the non-`dx = 12` rungs of both grid
ladders, which take their convective velocity from `solve_ifp` — was computed
under the mesh-dependent model and will move when the fix lands. I did not
re-run them.

---

## 2. The exported IFP metadata stated a screening length ~8.5x too large — LANDED

**Files:** `experiments/export_landing_asset.py`,
`landing-3d/public/models/ifp_contours.json`

### Confirmed twice, independently
The JSON's `display_coarse_grain_note` read "screening length ~60 um". The
screened-Poisson operator's decay length is `L = sqrt(K/(a+d))`.

1. **From the measurement.** `results/sweep/ifp_validation.json`, GBM at
   `dx = 12` — the exact configuration this asset is exported at — gives median
   `L = 7.081 um`, p5–p99 `6.61–9.79 um`. All eight section x resolution rows
   lie in `6.1–7.6 um`.
2. **From the shipped coefficients, by hand.** Away from the tumour `phi -> 0`,
   so `d -> _LAMBDA_LYMPH = 3.0e-2` and `a = vasc * 1.0e-3` is negligible beside
   it, while `K -> _K0 = 1.2`. That gives `L = sqrt(1.2 / 0.03) = 6.3 um`
   directly. 60 um would need `K/(a+d)` to be 100x larger than the constants at
   the top of `pressure.py` permit. The published figure is not a rounding
   disagreement; it is unreachable.

The measured value is right and the published one is wrong. (In passing: 60 um
is exactly `_TUMOR_SIGMA` x `dx` at production — the two were most likely
conflated, which is how §1 and §2 came to share a number without sharing a
cause.)

### Changed
Regenerated through the established export path — **not hand-edited**:

```
./.venv/bin/python experiments/export_landing_asset.py --no-glb
```

so the file and its generator cannot drift. The generator now carries the
measured diagnostics as named constants restated from
`results/sweep/ifp_validation.json` — restated, not recomputed, so this exporter
does not add a second implementation of the same quantity that could disagree
with the first.

**Before:**

> "The raw field's structure is per-vessel (screening length ~60 um) and lies
> below one screen pixel ..."

**After:** the corrected number, plus a new `screening_length_um` block carrying
the definition, median, p5–p99, source, and an explicit `correction_note`
recording what it supersedes.

**And the part that matters more than the number.** Because `L ~ 7 um` is below
*every* grid spacing the engine runs at — including this export's `dx = 12` — a
new top-level `resolution_caveat` states what the mesh actually resolves:

* `L < dx` at **99.8 %** of tissue nodes;
* the reaction term dominates the diagonal, median `(a+d)*dx^2/K = 2.87`;
* each node therefore relaxes toward its own local algebraic equilibrium
  `a*p_v/(a+d)`, and the solved field correlates with that pointwise limit at
  **r = 0.90** on this section;
* so this **is** a real solve of the stated equation — nothing authored, nothing
  fitted to this patient — but at this resolution it is closer to a pointwise
  nonlinear map of the vascular-density field than to a resolved transport
  solve, and the margin pressure gradient carrying the Jain effect is not
  resolved.

A provenance block that names the solver but not what the mesh resolves
overstates the artefact. Correcting 60 → 7.1 alone would have fixed the typo and
left the misleading impression intact.

### Shipped behaviour: metadata only
Verified by diffing the regenerated file against the previous one:

* `levels` — every iso-level, every path, every vertex — **identical**;
* `raw_field_max_mmhg`, `raw_field_mean_in_tissue_mmhg`,
  `displayed_field_max_mmhg` — identical;
* keys added: `screening_length_um`, `resolution_caveat`; keys removed: none;
* only pre-existing key whose value changed: `display_coarse_grain_note`.

**The landing page's rendering and its annotation were not touched.** The
numeric range was already removed from the label in `a5b96bf` and that label is
deliberate. Adding two keys is safe for the build because
`landing-3d/lib/ifpContours.ts` reads the JSON through
`raw as unknown as { ... }`, which does not excess-property-check.

The export in the tree was produced against `pressure.py` at `HEAD`, and
separately reproduced identically with §1's fix applied — so it is valid either
way and §2 does not depend on §1.

### Flagged, not touched (file I do not own)
`landing-3d/lib/ifpContours.ts` (module docstring, ~L27) repeats the same wrong
figure — "the field's screening length is of order 60 um" — in a source comment.
It does not reach the rendered page, but it is the same error in a file outside
my scope.

---

## 3. The grid ladders discarded the data that makes re-scoring free — LANDED

**Files:** `experiments/grid_convergence_gbm.py`,
`experiments/grid_convergence_breast.py`

### Confirmed for GBM; the breast half of the claim was wrong in its details
* **GBM — exactly as reported.** `run_ladder` built `e_unit_by_key`, populated
  it, returned it, and no caller or code path ever wrote it.
* **Breast — same defect, different mechanism.** `grid_convergence_breast.py`
  has no `e_unit_by_key` at all and never did. It computes `e_unit` as a local,
  reduces it, and lets it fall out of scope. The consequence is identical and
  the fix is the same, so I fixed both; but the specific claim that *both*
  scripts "build `e_unit_by_key` in memory, return it, and never write it" is
  true of one script only. Recorded because the description is what prompted the
  look.

Either way the rung records kept only `mean / median / p01 / p99 / max` (plus
`_e_unit_sum` / `_e_unit_first` as no-op guards), so a tail mean was not
recoverable and re-scoring under a different normalising convention cost a
1.5–2.5 h re-run.

### Changed
Both scripts now persist the per-cell vector the moment it exists, before it is
scored. New in each: `E_UNIT_STORE`, `store_e_unit`, `load_e_unit`,
`e_unit_manifest`, and an `e_unit_store` manifest embedded in the assembled
sweep JSON. Both the GBM ladder (whose `dt`-control points share `run_ladder`)
and the breast `dt`-control loop are covered.

**Format and cost.** One compressed `.npz` per script, one array per rung,
member name == the rung key verbatim so a record and its vector cannot be
mismatched. `float32`.

| store | rungs x cells | on disk |
|---|---|---:|
| `results/sweep/gbmgrid_e_unit.npz` | 30 x 23,812 | **2.4 MB** |
| `results/sweep/gridconv_breast_e_unit.npz` | 42 x 89,785 | **12.8 MB** |
| | **total** | **~15 MB** |

Measured, not estimated: compressing a real cached exposure vector
(`data_cache/kill_calibration_runs.json`) costs **3.40 bytes/element** in
`float32` and 7.40 in `float64`. The alternatives were `float64` (~33 MB total)
and inlining into the existing sweep JSONs (~90 MB of text, in files other tools
already parse). `float32` round-trips to `5.9e-8` relative on nonzero entries
and `8.9e-11` on the sum — orders of magnitude below anything a percentile or
tail mean of this quantity can resolve — and the `float64` `_e_unit_sum` in each
record is untouched, so the existing byte-level control against the published
cache is not weakened.

**Write discipline.** Rewritten after every solve, mirroring the JSONL
checkpoint's restart-safety, so an interrupted ladder keeps every vector it has
already paid for. Measured merge-and-rewrite cost: 0.08 s (GBM, 2.4 MB) and
0.4 s (breast, 13 MB), against multi-minute solves. Writes go to a temp file
followed by `os.replace`, so a kill mid-write cannot leave a corrupt archive.

**Purely additive.** No existing output changes meaning; the sweep JSONs and
JSONL checkpoints stay byte-compatible with everything already written against
them. The brief preferred extending the existing sweep JSONs — the manifest key
does exactly that, but the *vectors* cannot be back-filled into them, because
the data was never written and recovering it is the 1.5–2.5 h re-run the fix
exists to avoid.

**The ladders were not re-run.** The two `.npz` files therefore do not exist
yet; the next run of either script creates its own. The repo still holds only
the 25 per-cell exposure vectors it held before, all at `dx = 12`. What changed
is that the next run is the last one that has to pay for them.

### Verified
Both modules import cleanly, and `store_e_unit` / `load_e_unit` /
`e_unit_manifest` were exercised end to end on a temporary path at realistic
array sizes: incremental merge across two keys, exact round-trip within
`float32`, key names matching `_key`, and a JSON-serialisable manifest. No
existing store was touched (both target paths were asserted absent first).

---

## Test suite

`./.venv/bin/python -m pytest -q` → **754 passed**, matching the baseline taken
before any of these changes. No test was added, removed, skipped or edited.

For the record, the intermediate state with §1 applied was **750 passed, 4
failed**; that is the measurement §1 rests on, and it is why §1 is reverted
rather than shipped.
