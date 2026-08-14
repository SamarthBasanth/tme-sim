# Discrete-vessel baseline (no KDE at all) — does the already-default field solve the clustering problem?

One data point in the vascular-clustering sweep: does `mode='discrete'` (`tmesim/reconstruct/vessel_network.py`, promoted to default in `results/vessel_promotion.md`) already fix what `results/vascular_validation.md` found — that the reconstructed vasculature is statistically close to a random (Poisson/CSR) point field? Measured in-memory on `reconstruct_real('xenium', dx_um=12)`; no existing file edited, nothing run through the test suite.

**Verdict up front: NO — the discrete field is exactly as close to CSR as the number `vascular_validation.md` already published, because that published number was never computed from the KDE field in the first place.** This is not a new failure; it is confirmation that the whole KDE-bandwidth framing of the problem was a red herring.

## Provenance check (do this before reading anything else)

| check | result |
|---|---|
| `real_fields._DEFAULT_VASCULAR_MODE` | `'discrete'` |
| `reconstruct_real('xenium', dx_um=12).fields.vascular_density` bit-identical to explicit `mode='discrete'` rebuild | **True** |
| same field bit-identical to explicit `mode='kde'` rebuild | False |
| `vascular_density` nonzero fraction of tissue / mean over tissue | 2.68% / 0.0268 (discrete signature; KDE would be ~99.7% / ~0.084) |
| `vascular_validation.py` ever passes `vascular_mode` to `reconstruct_real`, or reads it anywhere | **Never** — 0 occurrences in the source |
| `vascular_validation.measure()`'s headline clustering test reads `micro.fields.vascular_density` | Only for the two things explicitly flagged as secondary in the report (`kde_area_fractions`, `field_based_distance`) — **not** for the headline CSR ratio |

**The load-bearing fact:** `vascular_validation.py`'s headline "0.99× / 1.31×" finding is computed by `cluster_vessels()` — DBSCAN/single-linkage directly on measured `ENDOTHELIAL`-cell **coordinates** — never touching `vascular_density` at all. `reconstruct.vessel_network.segment_vessels` (the literal source of the `discrete` field) runs the *same* single-linkage clustering at the *same* 20 µm link distance on the *same* coordinates. So the discrete field's own vessel segmentation and `vascular_validation.md`'s headline diagnostic are, to within a small denominator/implementation difference, **the same computation**. Neither was ever a property of the KDE's Gaussian bandwidth; both are a property of the raw endothelial point pattern.

## The clustering diagnostic, computed directly on the discrete field's own vessel objects

| quantity | value | published baseline (`vascular_validation.md`) |
|---|---:|---:|
| MVD (20 µm linking) | 110.8 /mm² | 111 /mm² |
| median tumour→vessel distance (observed) | 42.85 µm | 44.2 µm (42.8 µm in its own §3.4) |
| CSR median at this MVD (Poisson null) | 44.63 µm | 44.6 µm |
| **CSR clustering index, median (obs/CSR)** | **0.960×** | **0.99×** |
| tumour cells >100 µm from a vessel (observed) | 4.02% | 4.0% |
| CSR fraction >100 µm at this MVD | 3.08% | 3.1% |
| **CSR clustering index, far tail (obs/CSR)** | **1.307×** | **1.31×** |
| pixel-level field distance cross-check (median / p90), using the fixed `field>0` fallback mask on the actual `vascular_density` field | 43.3 / 84.9 µm | — (not previously computed on this field due to the now-fixed p95-threshold-degeneracy bug) |

Both the direct vessel-object computation and the independent pixel-level distance-transform on the actual field agree to within ~1 µm. **The discrete field reproduces the near-random arrangement to three significant figures.**

## Heterogeneity, 500 µm tiles, two ways

| quantity | value | note |
|---|---:|---|
| endothelial-**cell-count** based p90/p10 (mode-independent) | **2.42×** | reuses `vascular_validation.local_heterogeneity`; matches the published 2.4× exactly, as expected since it never reads the field |
| `vascular_density` **field-value** based p90/p10 (discrete field, this mode) | **2.25×** | computed here directly on the discrete field's own tile-mean values |
| real published CD34 vascular area fraction, across GBM cases (PMID 38928515) | 0.35–16.89% (**~48×** spread) | the target this heterogeneity is meant to approach |

The field-value heterogeneity is *slightly lower* than the cell-count heterogeneity, not higher — going discrete does not manufacture heterogeneity the raw point pattern doesn't already have. Both are roughly two orders of magnitude short of the ~48× spread real GBM shows across cases. Vessel area fraction where the engine's source is actually nonzero (`field>0`) is 2.68% of tissue, bracketing the published whole-section mean (2.1 ± 0.3%, PMID 36611806) reasonably — the *mean* amplitude is fine; it is exclusively the *spatial* unevenness that is missing.

## Downstream sanity check — cisplatin, dx = 20 µm, discrete default

| quantity | this run (dx=20) | `results/real_drug_panel.md` (dx=12) |
|---|---:|---:|
| alpha | 148.6082 | 148.6082 (same closed-form solve) |
| mean tumour AUC | 1.801e6 µM·s | 1.75e6 µM·s |
| coverage @1×IC50 | 0.9938 | 0.9936 (`vessel_promotion.md` §5, discrete arm) |
| kill fraction | 0.4239 | 0.4362 |

`reconstruct_real()` has **no `vascular_mode` parameter** — it always inherits whatever `real_fields._DEFAULT_VASCULAR_MODE` currently is, which is `'discrete'`. Per `vessel_promotion.md` §7.4, `real_drug_panel.md` was explicitly regenerated under the discrete default. **This is therefore not a new datapoint — it is a confirmation.** The dx=20 vs dx=12 numbers agree to within the grid-drift already characterized in `vessel_promotion.md` §2 (worst-drug absolute drift ≤0.007 for coverage, ≤0.03 relative for AUC), and nothing here moves the ranking or the amplitude story.

## Bottom line

The discrete field was already the default, and computing the clustering diagnostic on it (both via its own vessel-segmentation function and via the fixed pixel-level distance transform) reproduces `vascular_validation.md`'s "indistinguishable from a random vessel field" finding almost exactly (0.96–0.96× vs 0.99×, 1.31× vs 1.31×), because that finding was never a KDE artifact to begin with — it is a property of the measured PECAM1-gated endothelial-cell point pattern itself, which both fields are built from. Bypassing the Gaussian smoothing bypasses smoothing, not the underlying near-randomness; there was no clustering signal for the KDE to have erased. **The bandwidth question this sweep was framed around is moot: `results/vascular_validation.md` §7 item 4's diagnosis ("the smoothing kernel is the same size as the structure it is meant to resolve") is the wrong mechanism.** The real open question is why this section's PECAM1+ endothelial cells form a near-CSR point pattern at all — segmentation/gating resolution, single-section sampling, or a genuine biological feature of this tumour — and no change to `vascular_density_from_cells` (KDE bandwidth or discrete-vs-smooth) can fix that, because neither construction adds information the cell-typing step didn't already produce.

---
*Reproduce (read-only): `reconstruct_real('xenium', dx_um=12)`, `reconstruct.vessel_network.segment_vessels`/`segmentation_report`, `validation.vascular_validation.local_heterogeneity`, `validation.vessel_profile.vessel_mask_from_density`/`distance_from_vessel_um`, `validation.delivery_metrics.compile_real_drug`/`compute_all` — all called as-is, no source file modified.*
