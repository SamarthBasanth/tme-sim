# Is the percentile-threshold vessel-distance discrepancy a property of the method, or of one section?

Read-only cross-check. No file under `tmesim/` was modified; nothing was run through the test suite; this is a new file only.

## What was asked

`results/vascular_validation.md` is reported to have found the percentile-threshold vessel mask in `tmesim/validation/vessel_profile.py` (`vessel_mask_from_density`, p95 threshold + `distance_transform_edt`) giving a distance-to-vessel ~3x larger (134 µm) than the direct cell-position measurement (44 µm) on the GBM flagship section. This check runs **both** distance methods, exactly as implemented in `vessel_profile.py` and `vascular_validation.py`, on all three real sections this project has cached, to see whether that gap is a fixed property of the percentile-threshold method (present everywhere) or specific to one section/field.

**Important caveat up front: the specific 134 µm / 44 µm / ~3x figures could not be located or reproduced.** They do not appear anywhere in the current `results/vascular_validation.md`, `results/vessel_profile.md`, or `results/sweep/*.md` (checked by `grep -rn "134"` across `results/`). The currently-saved `vascular_validation.md` (timestamped after this check started, presumably from one of the two agents already working the problem) instead reports the *opposite* pathology on the incumbent section — a **0 µm** field-based distance, not 134 µm (its §5, "The opposing bias"). What follows is what actually happens when both methods are executed today, on all three sections, so the two agents auditing root cause have an independent, reproducible data point instead of a possibly-stale number.

## Method definitions (read from source, not assumed)

- **DIRECT (cell-position)**: `vascular_validation.distance_to_nearest(tumour_xy, endothelial_xy)` — Euclidean distance from every tumour cell to the nearest *measured* endothelial cell coordinate (`cKDTree` query, k=1). Computed from cell coordinates; independent of grid/dx by construction.
- **PERCENTILE-threshold**: `vessel_profile.vessel_mask_from_density(vascular_density, tissue_mask, percentile=95)` → `vessel_profile.distance_from_vessel_um` (`scipy.ndimage.distance_transform_edt`, `sampling=dx`). Thresholds `vascular_density` at its 95th percentile *over tissue*; **if that threshold comes out ≤ 0** (true when the field is sparse/discrete, i.e. the segmented-vessel default promoted in `results/vessel_promotion.md`), the function already falls back to `vasc > 0` rather than degenerately selecting the whole tissue mask. This fallback is existing code in `vessel_profile.py`, not something added by this check.
- Both methods were evaluated at the **same point set** (grid cell containing each tumour cell, via `Domain.to_index`) so the comparison is apples-to-apples, and cross-checked against the median over *all* tissue pixels (matching what `vessel_profile.run_arm` itself prints).

## Reconstruction calls used (verbatim from the named reports / their driving experiment scripts)

- GBM flagship: `reconstruct_real("xenium", sample="data_cache/gse324268_xenium", dx_um=12.0, state_mode="markers")` — the exact call in `experiments/flagship_markers.py:497`.
- GBM incumbent (the section `vascular_validation.md`'s own headline numbers, e.g. 44.2 µm, are computed from): `reconstruct_real("xenium", dx_um=12.0)` — included because it is the section the discrepancy claim traces back to.
- Breast: `reconstruct_real(source="xenium", sample="data_cache/xenium_breast_janesick/rep1", dx_um=12.0)` — verbatim from `results/breast_cancer_flagship.md`.

## Results

| section | n tumour cells | direct-method median (µm) | percentile-method median (µm) | ratio (pct / direct) |
|---|---:|---:|---:|---:|
| GBM flagship (GSE324268, 366-gene) | 515,887 | 20.5 | 17.0 | **0.83×** |
| GBM incumbent (brain_319, default xenium sample) | 23,812 | 42.8 | 43.3 | **1.01×** |
| Breast (Xenium Janesick, Rep1) | 89,785 | 31.2 | 26.8 | **0.86×** |

(All medians sampled at tumour-cell grid positions; the all-tissue-pixel basis gives materially the same ratios — 0.88×, 0.98×, 1.01× respectively — so the point-set choice doesn't drive the answer.)

**None of the three sections show anything close to a 3x gap.** All three cluster tightly around parity (0.83–1.01×), and where they diverge from 1.0 at all, the percentile method reads *shorter* than the direct method (on the two larger, denser-panel sections), not 3x longer.

## Where a real, large, and genuinely section-specific discrepancy DOES exist

It is not in `vessel_profile.vessel_mask_from_density` (above) — it is in `vascular_validation.py`'s own, separate `field_based_distance()` helper, which implements the identical p95-threshold logic **without** the `vasc > 0` fallback. Run directly on all three sections:

| section | direct-method median (µm) | `field_based_distance` (no fallback) median (µm) | vessel-mask threshold |
|---|---:|---:|---:|
| GBM flagship (GSE324268) | 20.5 | 17.0 | 1.0 (p95 lands at the field's max) |
| GBM incumbent (brain_319) | 42.8 | **0.0** (degenerate) | **0.0** (p95 lands at exactly 0) |
| Breast (Janesick Rep1) | 31.2 | 26.8 | 1.0 (p95 lands at the field's max) |

On the incumbent section only, the 95th-percentile threshold lands at exactly 0.0 — this section's `vascular_density` is the sparse/discrete segmented-vessel field (~97.5% of tissue pixels are exact zero; see `results/sweep/vasc_discrete_baseline_2026-07-29T05-54-00Z.md`), so `vasc >= 0` is trivially true almost everywhere and the "vessel mask" swallows the whole tissue, collapsing every distance to 0. On the two richer/larger sections (GSE324268, breast) the field is dense enough that the p95 threshold saturates at 1.0 (the field's own max, since it's max-normalised), so the unpatched function never hits the degenerate branch and agrees with the already-patched `vessel_profile.py` version exactly (0.83×/0.86× both places). This matches `results/sweep/vasc_discrete_baseline_2026-07-29T05-54-00Z.md`'s line calling the p95-threshold-degeneracy bug "now-fixed" — that fix exists in `vessel_profile.vessel_mask_from_density`'s `field>0` fallback, but has not been applied to `vascular_validation.field_based_distance`, which still exhibits it today.

## What this means for the two parallel fix efforts

The ~3x-overestimate framing does not generalize and is not reproducible against either distance function as currently written, on any of the three cached real sections. The reproducible, section-dependent failure runs in the *opposite* direction (collapse toward 0, not inflation toward 3x) and triggers on exactly one condition — a vascular-density field sparse enough that its own 95th percentile is ≤ 0 — which is a property of how discrete/segmented that section's vessel field is, not a property of the percentile-threshold *method* in the abstract. Among the three sections tested this condition is true only for the incumbent brain_319 section and false for both of the richer/denser flagship candidates (GSE324268 and the breast Xenium section), where the p95 threshold saturates at the field's own max instead. Practically: **the fix that matters is porting `vessel_profile.vessel_mask_from_density`'s existing `field > 0` fallback into `vascular_validation.field_based_distance`**, which is section-specific in exactly the same way (it will only change the incumbent section's number, not the two flagship candidates') — not a systematic ~3x correction that needs to be applied everywhere.
