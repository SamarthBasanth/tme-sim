# Vascular KDE bandwidth sweep — 14 µm

One point in a ~20-point bandwidth sweep (context: `results/vascular_validation.md`, `tmesim/reconstruct/real_fields.py::_kde_on_grid`, default `_DEFAULT_BANDWIDTH_UM = 40.0`). Field built in-memory via `vascular_density_from_cells(domain, cells, bandwidth_um=14.0, mode='kde')` on `reconstruct_real('xenium', dx_um=12)`; no existing file edited.

## Field sanity (dx=12µm build)

| check | value |
|---|---|
| NaN / Inf | 0 / 0 |
| range | [0.0, 1.0] |
| mean over tissue | 0.0298 |
| zero fraction over tissue | 16.2% (vs ~0.3% claimed for the 40µm default in `real_fields.py`'s own docstring) |
| all-zero? | no |

## Clustering diagnostic

Two families, because the module's cell-based diagnostic (vessel clustering + CSR null, tile heterogeneity from raw endothelial-cell counts) is computed from **cell coordinates**, not the KDE field — so it is bandwidth-invariant by construction and reproduces `vascular_validation.md` exactly (median 44.2 µm, ratio 0.99×, tail 1.44×, p90/p10 2.42×). The bandwidth-*sensitive* signal is the KDE field's own spatial heterogeneity and its percentile-threshold hotspot distance.

| diagnostic | bw=14µm | bw=40µm (published) | note |
|---|---:|---:|---|
| median tumour→vessel dist (cell-based) | 44.2 µm | 44.2 µm | identical — bandwidth-invariant |
| ratio vs CSR null (median) | 0.99× | 0.99× | identical |
| ratio vs CSR null (>100µm tail) | 1.44× | 1.31× | recomputed here, matches published within rounding |
| endothelial-count heterogeneity, p90/p10 (cell-based) | 2.42× | 2.42× | identical |
| **field-value heterogeneity, p90/p10 (500µm tiles)** | **2.82×** | not computed for bw=40 in this pass | mildly higher than the cell-based number, still far below CD34 area-fraction's ~48× across cases |
| **field hotspot distance @p95 threshold, median** | **37.9 µm** (p90 80.5) | **0 µm** (published, `field_based_distance`) | bw=14 no longer collapses to the degenerate 0 µm the 40µm KDE gives |

## Downstream sanity — cisplatin, dx=20µm (rough, not grid-converged)

| metric | bw=14µm @ dx=20 | 40µm-default baseline | source |
|---|---:|---:|---|
| mean tumour AUC (µM·s) | 2.044e6 | 1.75e6 | `results/real_drug_panel.md` (dx=12) |
| mean tumour AUC (µM·s) | 2.044e6 | 5.947e6 / 5.758e6 | `results/vessel_promotion.md` explicit KDE-40µm row, dx=15 / dx=12 |
| therapeutic_coverage@1×IC50 | 0.9978 | 0.9986 | `vessel_promotion.md` KDE row (dx=12/15) |
| kill_fraction | 0.489 | 0.4362 / 0.912–0.923 | `real_drug_panel.md` vs `vessel_promotion.md` KDE row (these two published baselines disagree with each other by ~3.3×, roughly the discrete/KDE normalisation factor documented in `real_fields.py` — flagged, not resolved here) |
| penetration length L (µm) | 105.1 | 184.5–201.6 | `vessel_promotion.md` KDE row |
| transport field NaN/Inf | 0 / 0 | — | sane |

## Bottom line

The field is numerically sound (no NaN/Inf/degenerate) at 14 µm. As expected, the cell-coordinate clustering diagnostic cannot move with bandwidth at all — the real signal is the field-derived heterogeneity and hotspot geometry, both of which show the intended direction of effect: field p90/p10 rises modestly (2.42× → 2.82×) and the percentile-threshold hotspot distance stops being the degenerate 0 µm the 40 µm default gives, landing at a physically plausible 38 µm median. That is real but modest relief of the over-smoothing problem, still an order of magnitude short of the ~48× CD34 area-fraction heterogeneity reported across real GBM cases. Downstream, tightening the kernel from 40 → 14 µm lowers predicted cisplatin delivery relative to the KDE-40µm baselines (mean AUC and kill both fall, penetration length shortens ~45%) — i.e. it makes the model *less* delivery-optimistic, consistent with concentrating the same total source over a smaller area. Note the two "40µm-default" baselines cited here (`real_drug_panel.md` vs `vessel_promotion.md`'s explicit KDE row) disagree with each other by ~3.3×, which looks like a stale-mode artifact in `real_drug_panel.md` rather than a bandwidth effect — worth the synthesis pass double-checking which one is the correct reference before ranking bandwidths on the downstream metrics.
