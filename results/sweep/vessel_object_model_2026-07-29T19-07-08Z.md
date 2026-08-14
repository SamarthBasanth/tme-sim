# Does the DBSCAN vessel-object model explain the near-random CSR finding? (lever b)

**A genuine attempt to move `results/vascular_validation.md` §6a's "too close to random" verdict by changing the vessel-clustering OBJECT MODEL itself, not the density field.** §6a already ruled out field-construction (9 KDE bandwidths, discrete, adaptive-kNN — none move the CSR ratio because the diagnostic never reads `vascular_density`) and named two untried levers: (a) richer marker-based endothelial identification (a parallel task) and (b) changing the DBSCAN vessel-clustering step in `measure()`. This is lever (b), tested directly rather than argued.

**Bottom line up front: lever (b) does not move the needle for GBM, and the honest reason why is now understood rather than assumed.** Every object-model alternative tried — a deliberate DBSCAN eps×min_samples sweep, HDBSCAN (a genuinely different, variable-density algorithm), and a spatially-varying empirical CSR null in place of the homogeneous-Poisson one — reproduces GBM's near-CSR result at every parameterisation that keeps the vessel count in a literature-plausible range. Several of them push the result to *look even closer to random* than the module's current default. For breast, the previously-flagged instability of the tail-clustering ratio across DBSCAN linking distances (1.2×–21.5× in `results/vascular_validation_breast.md` §3.3) is substantially explained here: it is largely an artefact of counting isolated single endothelial cells as their own one-cell "vessel," and stabilises once that convention is changed — but the underlying direction (breast shows more real clustering signal than GBM) survives and is not reversed by any candidate.

`tmesim/validation/vascular_validation.py` was read in full and **not modified**. `vv.measure()`, `vv.verdict()`, `vv.cluster_vessels()`, `vv.poisson_distance_stats()` and `vv.distance_to_nearest()` are called verbatim from a standalone prototype script; nothing under `tmesim/` was edited.

## 1. Reproduce: current parameterisation, read from the source

`tmesim/validation/vascular_validation.py`:

```python
VESSEL_LINK_UM = 20.0
VESSEL_LINK_SWEEP: Sequence[float] = (15.0, 20.0, 25.0, 30.0)
...
def cluster_vessels(endo_xy, link_um):
    labels = DBSCAN(eps=float(link_um), min_samples=1).fit_predict(endo_xy)
    ...
```

`min_samples=1` means DBSCAN degenerates to **single-linkage agglomeration**: every endothelial cell is automatically a core point, so no point is ever discarded as noise — an isolated cell more than `eps` from every other endothelial cell becomes its own one-cell "vessel" rather than being dropped. The module already documents (§3.2 of `results/vascular_validation.md`) that at the default `eps=20`, 81% of GBM "vessels" are exactly this: a single isolated cell.

`verdict()`'s CSR test compares `primary["tumor_to_vessel_median"]` (i.e. distance to the nearest DBSCAN-object centroid, at `eps=20`) against `poisson_distance_stats(q_rec)`, a **closed-form homogeneous-Poisson null at the object density** `q_rec = n_objects / area_mm2`. That is the whole CSR machinery this task targets.

## 2. Hypothesis

Two concrete, literature-grounded ways the current parameterisation could be masking real structure:

1. **`eps` too small relative to true inter-vessel spacing, `min_samples=1` too permissive.** A vessel lumen's endothelial lining cells sit ~10–30 µm apart (this section's own measured endothelial-cell nearest-neighbour distance is 22.9 µm median for GBM, 13.0 µm for breast — `results/vascular_validation.md` §3.2, `results/vascular_validation_breast.md` §2.2). Genuinely separate vessels sit far apart — this project's own directly-measured anchor (Yoshii & Sugiyama 1988, PMID 3359450) puts glioma intercapillary spacing at 77.8 ± 15.9 µm. At `eps=20` and `min_samples=1`, an endothelial cell more than 20 µm from any neighbour becomes its OWN vessel object even though it sits well inside the ~78–124 µm envelope of a real capillary's exclusion zone — i.e. the object count is inflated toward the raw cell count, and the CSR null is computed at a density much closer to "endothelial cell density" than "vessel density." Because `min_samples=1` can never denoise this (every point is a core point by construction), no `eps` choice alone can fix it — a `min_samples ≥ 2` requirement is a structurally different lever.
2. **The null itself may be wrong, not just its density.** `poisson_distance_stats` assumes one *global*, spatially homogeneous vessel density. Real tissue cellularity is heterogeneous for *every* cell type, not just endothelium (the module's own `pair_enrichment` control already shows STROMAL and IMMUNE cells enriched at multiple scales relative to an all-cells null). If some of the apparent "1.3–1.4× tail excess" over the homogeneous-Poisson null is actually generic tissue-density heterogeneity rather than vessel-specific clustering, a null that itself preserves that heterogeneity (drawn from real cell positions, the way `pair_enrichment` already works) should shrink the excess, not grow it — and that is a testable, falsifiable prediction, not a hand-wave.

## 3. Candidates tried

Standalone prototype (not committed under `tmesim/`), calling `vv.measure()`/`vv.verdict()`/`vv.cluster_vessels()`/`vv.poisson_distance_stats()`/`vv.distance_to_nearest()` unmodified:

- **(A) Deliberate DBSCAN eps × min_samples sweep.** `eps ∈ {8, 10, 12, 15, 18, 20, 25, 30, 40}` µm (finer near the ~10–30 µm endothelial lining spacing and extending to the ~40–125 µm real inter-vessel spacing), crossed with `min_samples ∈ {1, 2, 3}` — `min_samples ≥ 2` is the structurally different lever hypothesis 1 above calls for: it lets isolated cells become **noise** (dropped) rather than their own vessel.
- **(B) HDBSCAN** (`sklearn.cluster.HDBSCAN`, variable-density clustering — no single fixed radius) swept over `min_cluster_size ∈ {2,3,4,5,8}` and `min_samples ∈ {1, None}`, tested both as "clusters only" (noise points dropped entirely) and "clusters + singleton noise" (noise points kept as one-cell objects, matching the module's own `min_samples=1` DBSCAN convention) — this directly tests whether DBSCAN's fixed-eps assumption, rather than eps itself, is the problem.
- **(C) Spatially-varying (inhomogeneous) empirical CSR null.** Instead of `poisson_distance_stats`'s closed-form homogeneous-Poisson formula, draw the same number of "null vessel" points from the empirical positions of ALL measured cells (20 resamples, same philosophy as `vv.pair_enrichment`'s own null, which preserves the tissue outline and local cellularity gradients and destroys only which points are "vessels"), at each `eps` in the module's own `VESSEL_LINK_SWEEP`.

Run against the real sections exactly as the module's own `main()` loads them: `reconstruct_real("xenium", dx_um=12)` for GBM (40,887 cells, 4,295 endothelial, 27.90 mm² footprint) and `reconstruct_real("xenium", sample="data_cache/xenium_breast_janesick/rep1", dx_um=12)` for breast (167,780 cells, 13,677 endothelial, 38.57 mm² footprint) — matching `results/vascular_validation.md` and `results/vascular_validation_breast.md` exactly, confirmed by reproducing the module's own baseline verdict first: GBM `csr_median_ratio=0.991, csr_tail_ratio=1.307`, matching the on-disk report's "0.99×"/"~1.3×" to rounding.

**Note on a mid-session, external change to `verdict()` (not made by this task):** the breast baseline reproduced here as `csr_median_ratio=0.907, csr_tail_ratio=4.801` does NOT match the `0.9068`/`4.5361` recorded in the on-disk `results/vascular_validation_breast.md` §3.1, even though `csr_median_ratio` matches to rounding. Investigated directly (not assumed): `tmesim/validation/vascular_validation.py` changed by 66 lines partway through this task — **not edited by this task**, `cluster_vessels()`'s DBSCAN step (`eps=link_um, min_samples=1`, the actual subject of this report) is byte-for-byte unchanged — but `verdict()` picked up an unrelated fix, tagged in its own docstring as a "DENOMINATOR-BASIS FIX," that makes the tail ratio consistently use the vessel-OBJECT tail fraction (`prim["tumor_to_vessel_frac_gt_100"]`, 4.42%) instead of the previous code's mix of the raw endothelial-CELL tail fraction (4.17%, from `tumor_to_endo_stats`) divided by a vessel-object-density Poisson null — exactly the apples-to-oranges mismatch `results/vascular_validation_breast.md` §3.1 already flagged as a finding, not a bug it fixed. `4.17% / 0.92% = 4.536` (the old, on-disk number); `4.42% / 0.92% = 4.80` (the new, live number, reproduced here). This is a legitimate fix of a real inconsistency — likely made by the parallel lever-(a) task or a reviewer, not a problem with this report's methodology — and it does not touch the DBSCAN parameterization this report investigates. Every candidate result in §4–§5 below was computed with this task's own standalone `csr_against_homogeneous_poisson`/`csr_against_empirical_null` helpers, which always use the vessel-OBJECT distance distribution (never the raw-cell one) for both numerator and denominator, so they were unaffected by either version of `verdict()` and require no correction.

## 4. Results — GBM

### 4a. Candidate A: DBSCAN eps × min_samples

| eps (µm) | min_samples | n objects | MVD (/mm²) | spacing 1/√MVD (µm) | median ratio | tail ratio |
|---|---|---|---|---|---|---|
| 8 | 1 | 4,238 | 151.9 | 81 | 1.12 | **4.75** |
| 15 | 1 | 3,698 | 132.5 | 87 | 1.06 | 2.64 |
| 20 ← (module default) | 1 | 3,091 | 110.8 | 95 | 0.99 | 1.44 |
| 25 | 1 | 2,665 | 95.5 | 102 | 0.95 | 0.98 |
| 30 | 1 | 2,421 | 86.8 | 107 | 0.93 | 0.80 |
| 40 | 1 | 2,060 | 73.8 | 116 | 0.91 | 0.63 |
| 20 | **2** | 579 | 20.8 | 220 | 1.02 | 1.02 |
| 20 | **3** | 226 | 8.1 | 351 | 1.09 | 1.01 |

**The apparent "tail ratio blows up at small eps" is a density-mismatch artefact, not real clustering signal.** The observed tumour→vessel distance distribution barely changes with `eps` at `min_samples=1` (the module's own report already notes this: "insensitive to the linking distance, median 43.2–46.7 µm across the whole sweep") — because merging a handful of already-adjacent cells barely moves a centroid. What changes is the **object count fed into the Poisson-null denominator**, which at `eps=8` is 151.9/mm², i.e. essentially the raw endothelial-cell density (154/mm²), not a vessel density. A Poisson null at cell-like density predicts an almost-zero tail (`frac_gt_100 ≈ 0.85%` at that density), so dividing the (essentially fixed) observed tail by that artificially tiny null number inflates the ratio to 4.75× without reflecting anything about real spatial arrangement — this is exactly the "cells are not lumens" confound the module already documents, now shown quantitatively to be what drives the small-eps tail blow-up.

**At the eps range that keeps MVD in the literature-plausible band (65 ± 30 /mm², i.e. eps 25–40), both ratios sit at or below 1.0** (median 0.91–0.95×, tail 0.63–0.98×) — if anything a *slightly* less clustered result than the module's own default `eps=20`, not more.

**`min_samples ≥ 2` (denoising isolated cells rather than keeping them as one-cell vessels) collapses the object count far below any literature anchor** (20.8/mm² and 8.1/mm² vs. 65 ± 30/mm² published) — an over-correction, since a genuine capillary cut transversely near its edge can legitimately present as a single endothelial cell in a 2D section, so `min_samples ≥ 2` throws away real single-cell vessel profiles along with spurious ones. At those (too-sparse) densities the ratios sit almost exactly at 1.0 (median 1.02–1.09×, tail 1.01–1.02×) — a third, independent confirmation of near-CSR, now from the opposite (maximally conservative) end of the object-count spectrum.

### 4b. Candidate B: HDBSCAN (variable-density clustering)

Across `min_cluster_size ∈ {2,3,4,5,8}` × `min_samples ∈ {1, None}` × both noise-handling conventions (20 configurations): **median ratio 0.86–1.12×, tail ratio 0.63–1.23×** (one outlier configuration, `min_cluster_size=8, min_samples=None, clusters_plus_singleton_noise`, reaches 1.23× — still unremarkable). HDBSCAN's variable-density assumption — the direct test of whether DBSCAN's *fixed*-eps assumption was itself the problem — produces the same qualitative answer as DBSCAN at every setting tried. **This rules out "DBSCAN's fixed-radius assumption" as the explanation.**

### 4c. Candidate C: empirical (spatially-varying) CSR null vs. the homogeneous-Poisson null

| eps (µm) | n objects | homogeneous-null median ratio | homogeneous-null tail ratio | empirical-null median ratio | empirical-null tail ratio |
|---|---|---|---|---|---|
| 15 | 3,698 | 1.06 | 2.64 | 1.17 | 1.17 |
| 20 | 3,091 | 0.99 | 1.44 | 1.07 | **0.81** |
| 25 | 2,665 | 0.95 | 0.98 | 1.01 | 0.63 |
| 30 | 2,421 | 0.93 | 0.80 | 0.98 | 0.56 |

**The more careful, spatially-varying null makes GBM look *less* clustered than the homogeneous null does, not more**, at every eps tested — the tail ratio drops by 30–43% when the null is allowed to reflect real tissue-density heterogeneity instead of assuming one flat rate. This directly confirms the hypothesis in §2.2: part of the ~1.3–1.4× tail excess the module's default reports is generic tissue-cellularity heterogeneity (shared by all cell types, per `pair_enrichment`'s STROMAL/IMMUNE controls), not vessel-specific clustering. Correcting for it does not reveal *hidden* real clustering — it removes some of the *appearance* of clustering that the naive null was creating. This is the opposite of what would be needed to overturn §6a's verdict.

## 5. Results — breast

### 5a. Candidate A

| eps (µm) | min_samples | n objects | MVD (/mm²) | median ratio | tail ratio |
|---|---|---|---|---|---|
| 15 | 1 | 7,640 | 198.1 | 0.99 | 21.53 |
| 20 ← (module default) | 1 | 5,756 | 149.2 | 0.91 | 4.80 |
| 30 | 1 | 3,884 | 100.7 | 0.86 | 1.20 |
| 40 | 1 | 2,609 | 67.6 | 0.87 | 0.64 |
| 20 | **2** | 2,171 | 56.3 | 0.95 | 1.17 |
| 20 | **3** | 1,288 | 33.4 | 0.95 | 0.96 |
| 15 | **2** | 2,147 | 55.7 | 1.01 | 1.37 |
| 30 | **2** | 1,918 | 49.7 | 0.88 | 0.85 |

**This directly explains the instability `results/vascular_validation_breast.md` §3.3 already flagged** (tail ratio swinging 1.2×–21.5× across the 15–30 µm linking sweep at `min_samples=1`): breast's endothelial cells sit much closer together on average (13.0 µm median NN vs. GBM's 22.9 µm), so more of them remain "singleton" one-cell vessels at small `eps`, which pushes the object density toward the raw cell density (355/mm²) exactly as in GBM — the same density-mismatch artefact, just larger here because breast's clustering correction is doing more real work (§2.2 of the breast report already noted this). **Once `min_samples ≥ 2` removes singleton "vessels" from the count, the tail ratio stabilises dramatically** — 0.85×–1.37× across the same 15–30 µm eps range at `min_samples=2` (vs. 1.2×–21.5× at `min_samples=1`), and 0.96× (a single value, since `min_samples=3` collapses to a near-fixed object count) at `eps=20`. **The instability was substantially an artefact of the "isolated cell = one-cell vessel" convention, not evidence about real vascular clustering** — but note the *filtered* object density (33–56/mm²) is now below even the low end of breast's own (admittedly hotspot-biased) literature anchor (61–136/mm²), so this trades one bias (density inflation from singletons) for another (density deflation from over-filtering); it is not a free improvement.

### 5b. Candidate B: HDBSCAN

Across the same 20 configurations: **median ratio 0.86–1.04×, tail ratio 0.78–4.30×** (one outlier, `min_cluster_size=8, min_samples=None, clusters_plus_singleton_noise` = 4.30×, driven by the same singleton-inflation mechanism at a low `min_cluster_size`/permissive setting). Excluding that outlier, the bulk of configurations sit at 0.78–2.02× — modest, not dramatic, excess — consistent with, not contradicting, breast's already-reported "genuinely different, more textured than GBM but not real-GBM-pathology-level clustering" reading.

### 5c. Candidate C: empirical null

| eps (µm) | homogeneous-null median ratio | homogeneous-null tail ratio | empirical-null median ratio | empirical-null tail ratio |
|---|---|---|---|---|
| 15 | 0.99 | 21.53 | 1.31 | 12.60 |
| 20 | 0.91 | 4.80 | 1.19 | **6.77** |
| 25 | 0.87 | 2.12 | 1.13 | 3.94 |
| 30 | 0.86 | 1.20 | 1.11 | 2.45 |

**Breast diverges from GBM here, and in a way that reinforces rather than undermines the cross-tissue contrast already reported.** At 3 of 4 eps values (20/25/30) the empirical null makes breast's excess tail LARGER, not smaller (opposite of GBM's direction at every eps). Once genuine tissue-density heterogeneity is accounted for, breast's vasculature still shows more excess clustering than that heterogeneity alone would predict — while GBM's does not. This is exactly the asymmetry `results/vascular_validation_breast.md`'s verdict already describes ("both signs point toward more genuine vascular heterogeneity in this breast section than GBM's reconstruction showed") — now supported by a fairer null instead of only the homogeneous one, which is a real strengthening of that claim, not a new one.

## 6. Honest conclusion

**For GBM: no object-model choice tried — 8 fine-grained DBSCAN eps values, DBSCAN with noise-filtering (`min_samples` 2 and 3), a genuinely different variable-density algorithm (HDBSCAN, 20 configurations), or a spatially-varying empirical null in place of the homogeneous-Poisson one — overturns the near-CSR finding.** Every configuration that keeps the vessel count in a literature-plausible range (density comparable to 35–95/mm², spacing comparable to 78–124 µm) lands at median ratio 0.86–1.12× and tail ratio 0.56–1.44×, i.e. statistically unremarkable departure from a random point field. The apparent large tail-ratio numbers at small `eps` (up to 4.75×) are shown here to be a measurement artefact — comparing an almost-fixed observed distance distribution against an artificially tiny Poisson null caused by counting individual cells as vessels — not a signal that a different clustering choice reveals hidden structure. If anything, the more principled corrections tried (noise-filtering via `min_samples`, and the spatially-varying empirical null) push the result *closer to or below* the CSR baseline, not further from it.

**This means §6a's conclusion can now be stated more strongly and more specifically than before: the near-randomness of GBM's reconstructed vasculature is a property of the actual endothelial cell POSITIONS in this section, not an artefact of any post-hoc clustering or null-model choice tested — across two independent axes (field construction, per §6a, and object-model construction, per this report).** Combined with lever (a) (richer marker-based endothelial identification, tested in parallel), the two concrete remaining levers named in §6a have both now been tried on the object/identification side of the pipeline; the field-smoothing and vessel-clustering axes are exhausted.

**For breast: the finding is more nuanced and genuinely different from GBM, exactly as the existing breast report already found — but the specific number `vv.verdict()` currently reports as "the" tail ratio (4.54×, or 4.80× reproduced here) is now understood to be substantially inflated by counting isolated endothelial cells as one-cell vessel objects.** A `min_samples ≥ 2` requirement stabilises the previously-flagged eps-instability (1.2×–21.5× down to roughly 0.85×–1.37×) but trades density-inflation bias for density-deflation bias (33–56/mm² vs. the 61–136/mm² hotspot-biased literature floor), so it is not offered here as a strict improvement — only as a diagnosis of *why* the instability exists. The direction of breast's excess clustering (real, if modest) survives every candidate tried, including the more careful empirical null, which if anything strengthens it relative to GBM's null result.

## 7. Proposed diff — NOT applied, for review only

No change to `tmesim/validation/vascular_validation.py` is proposed as an unambiguous fix, because every candidate that changes the headline ratio does so by trading one form of bias for another (density inflation from singleton-counting vs. density deflation from over-filtering), and because the core finding — GBM is near-CSR — is unchanged regardless. If a reviewer wants to act on the singleton-instability diagnosis in §5a specifically, the smallest well-scoped change would be to report the `min_samples=1` and `min_samples=2` CSR ratios side by side (not replace one with the other), e.g. adding a `min_samples` axis to `VESSEL_LINK_SWEEP`'s reporting in `measure()`/`comparison_rows()` so the instability is visible in the report rather than only in the primary `eps=20, min_samples=1` number — a documentation/transparency change, not a "fix," and left for explicit review rather than applied here per the task's instructions.

## 8. Reproduce

Standalone research script, **not committed to the repo** (same convention as `results/sweep/vasc_bandwidth_breast_default_vs_15um_2026-07-29T05-54-00Z.md`): `vessel_object_model_prototype.py`, imports `tmesim.validation.vascular_validation` unmodified and calls `reconstruct_real` the same way `vv.main()` and `experiments/vascular_validation_breast.py` do. Run from the repo root with the project venv:

```
./.venv/bin/python /path/to/vessel_object_model_prototype.py           # GBM
./.venv/bin/python /path/to/vessel_object_model_prototype.py --breast  # breast
```

Full per-configuration JSON output for both sections (every DBSCAN eps×min_samples combination, every HDBSCAN configuration, and every empirical-null resample summarised in §4–§5) is on disk alongside the script; re-running reproduces the module's own baseline verdict first as a sanity check (GBM `csr_median_ratio=0.991, csr_tail_ratio=1.307`; breast `csr_median_ratio=0.907, csr_tail_ratio=4.801`) before running the candidates.
