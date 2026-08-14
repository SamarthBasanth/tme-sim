# Vascular KDE bandwidth sweep — 16 µm

One point in a ~20-point bandwidth sweep (context: `results/vascular_validation.md`, `tmesim/reconstruct/real_fields.py::_kde_on_grid`, default `_DEFAULT_BANDWIDTH_UM = 40.0`). Field built in-memory via `vascular_density_from_cells(domain, cells, bandwidth_um=16.0, mode='kde')` on `reconstruct_real('xenium', dx_um=12)`; no existing file edited. A `bandwidth_um=40.0` KDE field was built the same way, in the same run, as a directly comparable within-script baseline (see note below on why this is more trustworthy here than the published KDE figures).

## Field sanity (dx=12µm build)

| check | value |
|---|---|
| NaN / Inf | 0 / 0 |
| range | [0.0, 1.0] |
| mean over tissue | 0.0233 |
| nonzero fraction over tissue | 56.7% (vs 67.5% measured here for bw=40, vs ~99.7% claimed for the 40µm default in `real_fields.py`'s own docstring — see caveat below) |
| all-zero? | no |

## Clustering diagnostic

Two families. The module's **cell-coordinate** diagnostic (vessel clustering + CSR null, endothelial-count tile heterogeneity) is computed from raw endothelial-cell positions, never from the KDE field, so it is **bandwidth-invariant by construction** and reproduces `vascular_validation.md` exactly. The bandwidth-*sensitive* signal is the KDE field's own spatial structure.

| diagnostic | bw=16µm | bw=40µm (this run) | note |
|---|---:|---:|---|
| median tumour→vessel dist (cell-based) | 44.2 µm | 44.2 µm | identical — bandwidth-invariant |
| ratio vs CSR null (median) | 0.99× | 0.99× | identical, matches `vascular_validation.md` |
| ratio vs CSR null (>100µm tail) | 1.44× | 1.44× | identical (published 1.31×; small sweep/rounding difference, not a bandwidth effect) |
| endothelial-count heterogeneity, p90/p10 (cell-based) | 2.42× | 2.42× | identical |
| **field-value heterogeneity, p90/p10 (500µm tiles)** | **2.58×** | **2.55×** | essentially unchanged — 500µm tiles are ≫16–40µm, so the kernel-width difference washes out at that scale; both far below CD34 area-fraction's ~48× across cases |
| **field hotspot distance @p95 threshold, median (p90)** | **67.9 µm (146.5)** | **136.8 µm (303.6)** | bw=16's percentile-selected "hotspots" are *more numerous and closer together*, not fewer/further apart — narrowing the kernel breaks large blobs into many small local maxima, so a relative (percentile) threshold picks up more scattered small peaks. This is a threshold-selection effect, not evidence of "recovered clustering." |

**Correction on the field-distance baseline.** `vascular_validation.md` §5 reports a degenerate "0 µm" field-based distance — but that number was computed on the pipeline's *default* `vascular_density` field, which is `mode='discrete'` in the current codebase (zero over ~97.5% of tissue, so its p95 threshold collapses to 0 and the fallback selects the whole tissue). It is **not** a true 40µm-KDE number despite the module's docstring describing the field as a 40µm KDE. The bw=40µm figure above was built explicitly via `mode='kde'` in this script and gives a non-degenerate 136.8 µm, consistent with `vessel_promotion.md`'s own claim that "the KDE's p95 is strictly positive, so no existing number moves." Use this run's bw=40 row, not the published "0 µm," as the reference point.

## Downstream sanity — cisplatin, dx=20µm (rough, not grid-converged)

| metric | bw=16µm @ dx=20 | bw=40µm @ dx=20 (this run) | 40µm-KDE reference, dx=12 (`vessel_promotion.md` §2.1) | discrete-default, dx=12 (`real_drug_panel.md`) |
|---|---:|---:|---:|---:|
| mean tumour AUC (µM·s) | 2.356e6 | 5.074e6 | 5.758e6 | 1.75e6 |
| therapeutic_coverage@1×IC50 | 0.9983 | 0.9988 | 0.9986 | not reported there (`delivery_metrics.md` discrete: 0.9936) |
| kill_fraction | 0.663 | 0.929 | 0.912 | 0.4362 |
| transport field NaN/Inf | 0 / 0 | 0 / 0 | — | — |

My own bw=40 dx=20 run (5.07e6 AUC, kill 0.929, coverage 0.9988) sits close to the published dx=12 40µm-KDE reference (5.76e6, 0.912, 0.9986), which cross-checks the script's methodology; the residual gap is attributable to dx=20 vs dx=12 (not grid-converged, as instructed) rather than to the bandwidth swap.

## Bottom line

The bw=16µm field is numerically sound (no NaN/Inf, not degenerate, 56.7% nonzero). As with every other point in this sweep, the module's headline "close to random" clustering finding (median tumour→vessel distance 0.99× the CSR null) is computed from raw endothelial-cell coordinates and literally cannot move with bandwidth — that verdict is unaffected by any bandwidth choice in this sweep, 16µm included. The field-derived diagnostics that *are* bandwidth-sensitive give a mixed and partly counter-intuitive picture at 16µm: 500µm-tile heterogeneity is essentially unchanged from bw=40 (2.58× vs 2.55×, both far short of the ~48× real CD34-area-fraction spread across cases), and the percentile-threshold hotspot distance actually *shrinks* relative to bw=40 (68µm vs 137µm) because a narrower kernel fragments large smoothed blobs into many small, closely-spaced local maxima rather than recovering fewer, more separated true clusters. Downstream, narrowing 40→16µm lowers predicted cisplatin delivery substantially (mean AUC 5.07e6→2.36e6 µM·s, kill 0.929→0.663 in this run's matched pair), while `therapeutic_coverage@1×IC50` stays pinned near ceiling (~0.998–0.999) at both bandwidths and is therefore not a useful discriminator for this drug on this metric. Net: 16µm does not resolve the heterogeneity gap `vascular_validation.md` flags, and its main visible effect is to redistribute the field's mass into more, smaller local peaks rather than into larger clustered/avascular regions — the synthesis pass should weigh the 500µm-tile heterogeneity ratio (not the hotspot-distance number) as the more meaningful signal for ranking bandwidths, and should discard the published "0 µm" field-distance figure as a discrete-mode artifact rather than a true 40µm-KDE baseline.
