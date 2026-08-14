# Audit: `vessel_mask_from_density`'s percentile threshold — why it disagrees with the direct cell-based distance, and whether a percentile fix exists

Research-only pass. No files under `tmesim/` or `experiments/` were modified. All numbers below were produced by importing the **unmodified** `tmesim.validation.vessel_profile.vessel_mask_from_density` and `distance_from_vessel_um`, and the **unmodified** `tmesim.validation.vascular_validation.distance_to_nearest`, against one in-memory reconstruction of the flagship 10x Xenium GBM section (`reconstruct_real('xenium', dx_um=12)`, 40,887 cells, 23,812 tumour, 4,295 endothelial). Scratch driver scripts used to produce these numbers are not part of the repo.

## Verdict up front

**The percentile is not the real bug — the field it's applied to is.** `vessel_mask_from_density`'s 95th-percentile threshold produces a ~3× inflated distance **only when it is run on the smooth KDE field** (`vascular_mode='kde'`, the historical default). Under the **current production default** (`vascular_mode='discrete'`, since `results/vessel_promotion.md`), the same function's own degenerate-threshold fallback (`vasc > 0`) makes the percentile choice almost irrelevant and already reproduces the ground-truth cell-based distance to within 1%. So there are really two separate findings to act on, not one:

1. On the **KDE field**, no single percentile is a robust fix — the "correct" percentile itself shifts by ~15 points as the KDE bandwidth moves over its own documented sweep range (25–60 µm), because percentile and bandwidth both control the same thing (how much of the smoothed field counts as "vessel") and trade off against each other.
2. The **discrete field already gets the right answer**, but only through a fallback branch whose original purpose (guard against a degenerate all-tissue mask) is incidental to why it happens to also be geometrically correct. This is fragile: it works because the discrete field's nonzero footprint (~2.7% of tissue) already approximates true vessel geometry, not because `vessel_mask_from_density`'s logic was designed to measure "nearest vessel."

**Recommendation:** replace the percentile-threshold construction with a direct distance calculation from clustered endothelial-cell/vessel positions (option (b) below), matching `vascular_validation.py`'s already-correct method, rather than trying to pick a better percentile. Details and reasoning in §7–8.

---

## 1. What `vessel_mask_from_density` actually does (code trace)

`tmesim/validation/vessel_profile.py:173-218`:

```python
def vessel_mask_from_density(vascular_density, tissue_mask, percentile=VESSEL_PERCENTILE):
    vasc = np.asarray(vascular_density, dtype=np.float64)
    tissue = np.asarray(tissue_mask, dtype=bool)
    thr = float(np.percentile(vasc[tissue], float(percentile)))
    if thr > 0.0:
        mask = (vasc >= thr) & tissue
        rule = f"p{float(percentile):g} threshold on the field"
    else:
        mask = (vasc > 0.0) & tissue
        rule = "field > 0 (sparse/segmented source; p-threshold degenerate at 0)"
    ...
```

`VESSEL_PERCENTILE = 95.0` is the module-wide default, used by `run_arm`, `threshold_sensitivity`, `bandwidth_sensitivity`, and by every external caller that omits the `percentile` kwarg (see §6).

Two completely different regimes fall out of this one function depending on what `vascular_density` field it's handed:

- **KDE field** (`vascular_mode='kde'`, 40 µm Gaussian bandwidth on measured endothelial-cell positions): only 0.22% of tissue pixels are exactly zero (measured directly), so `np.percentile(vasc[tissue], 95)` returns a strictly positive threshold and the function takes the `thr > 0.0` branch — a genuine percentile threshold on a smooth, mostly-nonzero field.
- **Discrete field** (`vascular_mode='discrete'`, the production default since `vessel_promotion.md`): 97.32% of tissue pixels are exactly zero (a discrete, painted vessel-radius source), so `np.percentile(vasc[tissue], p)` is 0.0 for every `p` up to ~97.3, and the function falls into the `thr <= 0.0` branch — `vasc > 0`, i.e. "wherever the segmented vessel source is nonzero at all." **The percentile argument is then not doing anything** for any `p ≤ ~97.3` (see the flat rows in §4 Table A) — the mask is identical regardless of whether you asked for p50 or p95.

## 2. Why the KDE-mode percentile threshold overshoots — the mechanism

A percentile threshold on a max-normalised, smoothly-varying field answers "where is the field in its own top X%?", not "where is a vessel?". Two effects compound:

1. **A smooth KDE has no zero baseline to anchor a threshold to.** With 40 µm bandwidth, 99.78% of tissue has *some* positive endothelial-KDE value — the field is a continuous hill-and-valley surface, not a sparse indicator of vessel footprint. A percentile threshold therefore always returns *some* nonzero area fraction, no matter how the threshold is chosen, and that area fraction is placed wherever multiple nearby endothelial cells' kernels constructively overlap — i.e. **the biggest, most crowded vessel clusters**, not "any pixel within X µm of any vessel." True vessel geometry (many small, fairly evenly spread objects — median vessel-vessel spacing ≈ 41–60 µm per `results/vascular_validation.md`) gets replaced by a small number of large "hotspot" blobs.
2. **The blob-vs-scatter mismatch grows with the percentile.** Lower percentiles keep more area and therefore approach (but never equal) the true vessel footprint; higher percentiles keep only the hottest hotspots, which are spatially sparser than real vessels, so the distance from an average tumour cell to the nearest surviving blob grows faster than the area fraction shrinks. Measured on the flagship section (bandwidth = 40 µm, the default):

| percentile | vessel area fraction | median tumour→"vessel" distance (µm) | ratio to ground truth (42.85 µm) |
|---:|---:|---:|---:|
| 50 | 50.0% | 0.0 | 0.00 |
| 70 | 30.0% | 26.8 | 0.63 |
| 75 | 25.0% | 36.0 | 0.84 |
| **78** | **22.0%** | **43.3** | **1.01** |
| 80 | 20.0% | 48.0 | 1.12 |
| 85 | 15.0% | 61.2 | 1.43 |
| 90 | 10.0% | 86.5 | 2.02 |
| **95 (current default)** | **5.0%** | **136.8** | **3.19** |
| 97 | 3.0% | 192.0 | 4.48 |
| 99 | 1.0% | 429.3 | 10.02 |

(Ground truth = `vascular_validation.py`'s `distance_to_nearest(tumour_xy, endothelial_xy)` median, 42.85 µm — the same anchor `results/vascular_validation.md` reports as "44 µm" to one significant figure with a slightly different vessel definition.)

The p95 number here, **136.8 µm**, matches the ~134 µm figure that motivated this audit (the small difference is almost certainly just which reconstruction pass / rounding the original number came from — this rerun is internally consistent to itself).

## 3. Does grid resolution (dx) or the field's own smoothing (bandwidth) compound this?

**Bandwidth compounds it heavily; grid resolution does not.**

Fixed-percentile, bandwidth sweep (KDE field rebuilt at each bandwidth via `build_real_fields(..., bandwidth_um=bw, vascular_mode='kde')`, same dx = 12 µm, same reconstruction):

| bandwidth (µm) | p80 median d (µm) | ratio | p90 median d (µm) | ratio | p95 median d (µm) | ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | 26.8 | 0.63 | 60.0 | 1.40 | 91.4 | 2.13 |
| 40 (default) | 48.0 | 1.12 | 86.5 | 2.02 | 136.8 | 3.19 |
| 60 | 72.0 | 1.68 | 129.2 | 3.02 | 204.4 | 4.77 |

At a **fixed percentile**, moving the bandwidth alone from 25→60 µm moves the answer by up to **2.7×** (p80: 26.8 → 72.0 µm). This is corroborated by the independent bandwidth-sweep pass in `results/sweep/vasc_bandwidth_12um_2026-07-29T05-54-00Z.md`, which measured the same quantity ("field-based hotspot distance, median (p95 threshold)") at bandwidth = 12 µm and got **37.9 µm** — consistent with the monotonic bandwidth→distance relationship found here (smaller bandwidth ⇒ smaller hotspot distance).

Because bandwidth and percentile both control the same underlying thing (how selective the threshold is on a field whose *shape* itself changes with bandwidth), **the percentile that best matches ground truth also moves with bandwidth** — it is not a fixed correction:

| bandwidth (µm) | percentile giving median d ≈ 42.85 µm (ratio ≈ 1.0) |
|---:|---:|
| 25 | **~p85** (43.3 µm, ratio 1.01) |
| 40 (default) | **~p78** (43.3 µm, ratio 1.01) |
| 60 | below p78 (p78 already overshoots at ratio 1.43); extrapolating the trend, roughly **~p68–72** |

A ~15–17 percentile-point swing in "the correct answer" across the module's own already-documented bandwidth sweep range means **there is no percentile that is simultaneously correct and bandwidth-robust.** Any fixed percentile chosen to match today's 40 µm bandwidth will silently drift wrong if the bandwidth default ever changes (or if a caller passes a different bandwidth, which several already do — see `bandwidth_sensitivity` in §6).

Grid resolution, by contrast, has a much smaller effect. Fixed p95, KDE field, dx sweep:

| dx (µm) | median tumour→vessel distance (µm) |
|---:|---:|
| 30 | 153.0 |
| 20 | 141.4 |
| 12 | 136.8 |

A ~12% drift, monotonically converging as dx shrinks — second-order next to the ~170%+ swing bandwidth produces. **Grid resolution is not the compounding factor here; the field's own smoothing kernel is.**

## 4. The discrete field's accidental correctness (Table A, full data)

Under the current default field (`vascular_mode='discrete'`), the percentile is inert for any value ≤ ~97.3 (the field's own zero fraction), because every such percentile hits the `thr <= 0.0` fallback and collapses to the same `vasc > 0` mask:

| percentile | threshold | rule taken | vessel area % | median d (µm) | ratio to ground truth |
|---:|---:|---|---:|---:|---:|
| 50–97 | 0.0000 | `field > 0` fallback | 2.683 | 43.3 | 1.01 |
| 97.4–99.5 | 1.0000 | genuine p-threshold | 2.683 | 43.3 | 1.01 |

(Above ~97.4 the genuine percentile threshold also happens to land back on `vasc==1` exactly, i.e. the same discrete-source footprint, because the discrete field only takes two values inside tissue: 0 and 1. So even the "genuine" branch reproduces the fallback's mask here — not a coincidence to rely on, since it depends on the discrete field being binary, which is an implementation detail of `discrete_vascular_field`, not a documented contract.)

**This means the production default is already within 1% of the cell-based ground truth today** — but for a reason that has nothing to do with percentile selection being correct. It's correct because (a) the discrete field's nonzero footprint already *is* a segmented vessel mask (built by `reconstruct.vessel_network.segment_vessels` + `discrete_vascular_field`, independently of any percentile logic), and (b) the fallback branch happens to select exactly that footprint. If a future discrete field construction ever made the "zero fraction" smaller than 100 − 95 = 5% (i.e. more than 5% of pixels nonzero), the `thr > 0.0` branch would activate and the *real* percentile-driven hotspot bias documented in §2–3 would silently reappear, indistinguishable at a glance from today's correct behavior.

## 5. The sibling bug this audit surfaced (not `vessel_mask_from_density`, but its blind reimplementation)

`tmesim/validation/vascular_validation.py:499-530`, `field_based_distance()`, is a **separate, hand-rolled reimplementation** of the same percentile-threshold idea, used only as an internal cross-check inside that module (`grid_invariance`, and the "opposing bias" section of `results/vascular_validation.md`):

```python
thr = float(np.percentile(vasc[tm], percentile))
mask = (vasc >= thr) & tm
```

This is **missing `vessel_mask_from_density`'s degenerate-threshold fallback**. Under the current default discrete field, `thr` comes out `0.0` (same as above), but `mask = (vasc >= 0.0) & tm` selects the **entire tissue mask** (since the field is never negative) instead of falling back to `vasc > 0`. The result — visible today in the live `results/vascular_validation.md` (§5, "The opposing bias") — is a nonsensical **median field-based distance of 0 µm**, printed as "0.0× larger" and "nan×" across the dx sweep. That report's "two incompatible answers" framing (44 µm vs. "0 µm") is consequently **not** a faithful representation of what `vessel_mask_from_density` itself currently returns (43.3 µm, ratio 1.01, per §4) — it's an artifact of `vascular_validation.py` not calling the real function and instead re-deriving a buggy variant of it. This should be fixed by either (i) calling `vessel_mask_from_density` directly from `vascular_validation.py` instead of hand-rolling `field_based_distance`, or (ii) adding the same `thr > 0.0` guard there. This is a one-line, low-risk fix and is flagged here only because it was found in the course of this audit; it is not itself the KDE-percentile hotspot bug this audit was asked to investigate, and whoever centrally implements the recommendation in §8 should sweep it up in the same pass since it touches the identical logic.

## 6. Callers — what depends on `vessel_mask_from_density`'s current behaviour

```
grep -rn "vessel_mask_from_density" --include="*.py" .
```

| caller | field passed | percentile | exposure to the KDE-mode bias |
|---|---|---|---|
| `experiments/regimen.py` | `micro.fields.vascular_density` (default = discrete) | default (95) | **none today** — inert percentile, discrete fallback active |
| `experiments/combination.py:402` | `micro.fields.vascular_density` (default = discrete) | explicit `VESSEL_PERCENTILE` (95) | **none today** |
| `experiments/dose_schedule.py:638` | `micro.fields.vascular_density` (default = discrete) | default (95) | **none today** |
| `tmesim/validation/delivery_metrics.py:451,738` | `micro.fields.vascular_density` (default = discrete) | default / configurable `vessel_percentile` kwarg | **none today**, but the kwarg exists — a caller passing a percentile > ~97.3 on the discrete field, or switching to KDE, re-exposes the bias |
| `tmesim/validation/heldout_patient.py:341` | `micro.fields.vascular_density` (default = discrete) | `DM.VESSEL_PERCENTILE` (95) | **none today** |
| `tmesim/validation/vessel_profile.py` (`run_arm`, `threshold_sensitivity`) | `micro.fields.vascular_density` (default = discrete, since `build_arm_microenvironment` never overrides `vascular_mode`) | 95, and the 90/95/98/99 sweep in `threshold_sensitivity` | **none today** — but see next row |
| `tmesim/validation/vessel_profile.py` (`bandwidth_sensitivity`, §6c "THE DECISIVE CONTROL") | rebuilds fields via `build_real_fields(dom, base.cells, bandwidth_um=bw)` **without overriding `vascular_mode`**, so it stays on `discrete` | 95 (via `VESSEL_PERCENTILE`) | **This sweep is currently measuring nothing.** `vascular_density_from_cells`'s `discrete` branch does not accept or use `bandwidth_um` at all (see `real_fields.py:186-195`, "`bandwidth_um` is ignored in this mode"). So the "KDE bandwidth sweep" in `vessel_profile.md` §"THE DECISIVE CONTROL" rebuilds an **identical vessel mask and distance map at all three bandwidths** — confirmed directly: `vessel_mask_from_density` is called on `vasc` from a field that literally cannot change with `bw` in discrete mode. The small residual drift already visible in that table (doxorubicin L = 34.7 / 34.0 / 33.0 µm at bw = 25/40/60) is coming entirely from `ecm_density_from_cells`'s stromal-KDE bandwidth (which *does* use `bandwidth_um`) changing tissue diffusivity/tortuosity slightly — **not** from any change in vessel geometry. The module's own docstring ("this is the decisive test of whether the fitted decay length is drug transport or our own smoothing kernel") is currently not testing what it claims to test on the vessel-mask/distance axis. This is a second, adjacent bug worth flagging to whoever synthesizes the parallel bandwidth-sweep investigations. |
| `tmesim/validation/vessel_network_impact.py:457` | **explicit** `fields["A_kde"]` | `VESSEL_PERCENTILE` (95), labelled `"LEGACY axis"` | **exposed, but intentionally** — this call exists specifically to reproduce the historical KDE-era 26–28 µm numbers for comparison and is documented as such |
| `tmesim/validation/vessel_network_impact.py:930` (`stale_baseline_check`) | `vascular_density_from_cells(m.domain, micro.cells)` — **note: this variable is named `kde` but, since the function's default mode flipped to `discrete`, it is actually building the discrete field**, a naming staleness bug adjacent to but distinct from the main finding | `VESSEL_PERCENTILE` | Currently inert (discrete fallback), but the surrounding narrative text still describes it as "the KDE" |
| `app/server/api.py`, `app/server/cache.py`, `app/server/jobs.py` | user-selectable `vascular_mode` request parameter, default `"discrete"` but **`"kde"` is an exposed, documented option** (`api.py:373`, `cache.py:806-878`) | whatever the downstream validation/experiment module uses (95, typically) | **This is the one live path by which an ordinary user (not just a legacy-comparison test) can hit the full 3.19× hotspot bias in production** — anyone requesting `vascular_mode='kde'` through the app gets a vessel mask built by exactly the mechanism audited in §2–3. |
| `tests/test_vessel_network.py:163-190` | both `discrete` and `kde`, explicitly | default | **Tests, and protects, exactly the fallback behaviour described in §1 and §4** (`test_sparse_source_does_not_degenerate_the_vessel_mask`). Any fix must not break this test's assertions: `vd["threshold"] == 0.0`, `"field > 0" in vd["rule"]`, `vd["mask"] == (disc > 0) & tissue` for the discrete case; `vk["threshold"] > 0.0`, percentile rule, for the KDE case. |

**Summary of caller exposure:** every production caller that goes through `micro.fields.vascular_density` under the current default is *currently* unaffected by the KDE-mode hotspot bias, purely because the default field switched to `discrete` (a change made for unrelated reasons in `vessel_promotion.md`) and the fallback branch happens to catch it. The bias is live and reachable today only via (a) the app server's exposed `vascular_mode='kde'` option, and (b) the two intentional/adjacent legacy-comparison call sites in `vessel_network_impact.py`. It would become live again everywhere, silently, if the discrete field's zero-fraction ever dropped below ~5% (activating the genuine-percentile branch on a field it wasn't designed for) or if `vascular_mode='kde'` were ever promoted back to default.

## 7. Two fix options

### (a) A corrected percentile value

Rejected as a real fix, for reasons already quantified above:

- There is no percentile that is correct across the module's own documented bandwidth range: ~p85 at bw=25, ~p78 at bw=40, ~p68–72 (extrapolated) at bw=60. Hard-coding one of these only fixes the specific (bandwidth, percentile) pair it was tuned on.
- It does nothing to address the discrete-mode fallback being "correct by accident" (§4) — a future change to the discrete field's zero-fraction could silently reintroduce the bias even with a "corrected" percentile, since the corrected value was tuned against the KDE field, not the discrete one.
- It leaves `vessel_mask_from_density` conceptually measuring "top X% of a smoothed density field" when the question actually being asked, everywhere it's used, is "distance to the nearest vessel" — a category mismatch that a better-tuned threshold papers over rather than closes.

If this were pursued anyway despite the above, **p78 (bandwidth=40, the current default)** is the closest available match (ratio 1.01 at exactly that bandwidth), against the current `VESSEL_PERCENTILE = 95.0` — but this number should be treated as disposable, not load-bearing, given how much it moves with bandwidth.

### (b) Replace the percentile-threshold approach with a direct distance calculation

Compute distance-to-nearest-vessel directly from clustered endothelial-cell (or segmented-vessel-object) positions, the way `vascular_validation.py`'s `distance_to_nearest` / `cluster_vessels` already does, instead of thresholding a density field at all.

**Concrete shape of the change**, scoped to `tmesim/validation/vessel_profile.py`:

1. Add a new function, e.g. `vessel_mask_from_cells(cells, domain, link_um=VESSEL_LINK_UM)`, that either (i) reuses `reconstruct.vessel_network.segment_vessels` (already used by the discrete-field pipeline itself — no new clustering logic needed) to get vessel objects, then paints them onto the grid at their true footprint, or (ii) skips the mask/EDT step entirely for the tumour-cell-specific use case and calls `scipy.spatial.cKDTree` / `vascular_validation.distance_to_nearest` directly on cell coordinates, matching `vascular_validation.py` exactly.
2. `distance_from_vessel_um` (the EDT wrapper) stays as-is; it just needs a mask that is already geometrically a vessel footprint rather than a percentile-selected subset of a density field. Since `reconstruct.vessel_network.segment_vessels` + `discrete_vascular_field` already build exactly this for the `discrete` vascular mode, **the discrete field's own construction already contains the fix** — the cleanest version of option (b) may be as small as: always build the mask from `discrete_vascular_field`'s vessel geometry (or `vasc > 0` on it) regardless of which `vascular_density` field is being used for the *drug source* term. That decouples "what field the drug is delivered through" (KDE or discrete, a legitimate, still-open modelling choice) from "what counts as a vessel for the distance axis" (should always be the actual segmented vessel objects).
3. `VESSEL_PERCENTILE` / `VESSEL_PERCENTILE_SWEEP` and `percentile=` kwargs across the 6 caller sites in §6 would be removed or repurposed as a no-op / deprecated argument, with a clear migration note, since there would no longer be a percentile to choose.
4. `tests/test_vessel_network.py`'s `test_sparse_source_does_not_degenerate_the_vessel_mask` would need to be rewritten (not just left passing) — it currently asserts on percentile/threshold/rule fields that wouldn't exist in the new function. This is expected churn, not a regression.
5. `vessel_network_impact.py`'s intentional `"LEGACY axis"` KDE comparison (§6) should probably be **kept** as an explicit, clearly-labelled special case for reproducing historical numbers — it already understands it's measuring something different ("legacy p95 axis") and doesn't need this fix applied to it.
6. `bandwidth_sensitivity` (§6, "THE DECISIVE CONTROL") would need a decision either way: if the vessel/distance axis is decoupled from the drug-source field's bandwidth (per point 2 above), the module's bandwidth sweep should be re-scoped to explicitly test *only* what it can actually affect under the new design (ECM tortuosity, and the KDE-mode drug-source field if `vascular_mode='kde'` is exercised) rather than implicitly promising to test vessel-mask sensitivity it can no longer produce.

## 8. Recommendation

**Go with (b).** The percentile-threshold approach is not just currently mistuned, it is the wrong tool for the question being asked. `vessel_mask_from_density` is used everywhere as an answer to "how far is this tissue from a vessel", but a percentile on a density field can only ever answer "how far is this tissue from the densest few percent of the field" — those coincide only by luck (as they currently do, via the discrete-mode fallback) and diverge badly and unpredictably otherwise (§2–3). A corrected percentile is a parameter tuned to today's specific (field, bandwidth) pair; it will look right in this report and drift wrong the next time the bandwidth default, the discrete field's sparsity, or the vascular mode changes — and nothing in the code would flag that drift, because the function would still run and return a number.

The engineering cost of (b) is genuinely small relative to that risk: the codebase already has both ingredients built and tested — `reconstruct.vessel_network.segment_vessels` (vessel clustering) and `vascular_validation.distance_to_nearest` / `cluster_vessels` (the exact distance computation this audit used as ground truth throughout). The fix is substantially "point `vessel_mask_from_density`'s callers at code that already exists and is already validated against the literature in `vascular_validation.md`," not "write new geometry logic." The caller list in §6 shows the blast radius is six files plus one test, all of which currently pass a plain `micro.fields.vascular_density` with no exotic overrides — a mechanical, low-risk migration.

The one thing worth deciding explicitly in the centralized fix (not something this audit should pre-empt): whether the *drug source term* (`delivery = vascular_density × bbb_penetrance × plasma`) stays free to be KDE or discrete (a legitimate, separately-debated modelling question, per `vessel_promotion.md`), while the *distance-to-vessel axis* always uses the direct cell/segmented-vessel geometry regardless of which source field is active. That decoupling is what fixes both the percentile bug and the currently-inert bandwidth sweep (§6) in one change, and it is consistent with how `vascular_validation.py` already treats these as two separate questions.

---

*Scope note:* this file only reports findings; `tmesim/validation/vessel_profile.py`, `vascular_validation.py`, and all other project files were read but not modified. No test suite was run. Driver scripts used to produce the numbers above were scratch files outside the repository and are not retained here.
