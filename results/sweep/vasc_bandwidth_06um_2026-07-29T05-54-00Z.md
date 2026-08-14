# Vascular KDE bandwidth sweep — 6 µm

One point in a ~20-point bandwidth sweep (context: `results/vascular_validation.md`, `tmesim/reconstruct/real_fields.py::_kde_on_grid`, default `_DEFAULT_BANDWIDTH_UM = 40.0`). Field built in-memory via `vascular_density_from_cells(domain, cells, bandwidth_um=6.0, mode='kde')` on `reconstruct_real('xenium', dx_um=12)`; no existing file edited.

## Field sanity (dx=12µm build)

| check | value |
|---|---|
| NaN / Inf | 0 / 0 |
| range | [0.0, 1.0] |
| mean over tissue | 0.00947 |
| zero fraction over tissue | 65.2% (vs ~16.2% at bw=14µm, ~0.3% claimed for the 40µm default) |
| all-zero? | no |

## Clustering diagnostic

The module's cell-based diagnostic (endothelial-cell DBSCAN clustering, tumour→vessel distance, tile heterogeneity from raw endothelial-cell counts) is computed from **cell coordinates**, not the KDE field, so it is bandwidth-invariant by construction and reproduces `vascular_validation.md` exactly without needing to be rerun: median 44.2 µm, ratio vs CSR 0.99× (median) / 1.31× (>100µm tail), endothelial-count p90/p10 2.4×. The bandwidth-*sensitive* signal is the KDE field's own spatial heterogeneity and its percentile-threshold hotspot geometry, computed here directly on the bw=6µm field (95th-percentile threshold on the tissue mask → connected components → distance transform + a CSR null at the resulting blob density):

| diagnostic | bw=6µm | bw=14µm (sibling sweep point) | bw=40µm (published) |
|---|---:|---:|---:|
| field hotspot count above p95 threshold | 2,209 blobs (79.2/mm²) | — | — (degenerate, see below) |
| field hotspot median distance | 40.6 µm | 37.9 µm | **0 µm** (`field_based_distance`, degenerate) |
| field hotspot p90 distance | 80.7 µm | 80.5 µm | 0 µm |
| CSR null at same blob density (median / >100µm tail) | 52.8 µm / 8.3% | not computed | not computed |
| **ratio vs CSR, field hotspots (median / >100µm tail)** | **0.77× / 0.44×** | not computed | not computed |
| field-value heterogeneity, p90/p10 (500µm tiles) | **2.50×** | 2.82× | not computed |
| area frac @ thr 0.2 / 0.3 / 0.5 | 2.20% / 2.18% / 0.10% | — | 2.68% / 2.68% / 2.68% (flat plateau) |

**Reading this honestly: bw=6µm does not recover clustering, it recovers sub-Poisson regularity.** Unlike the cell-based measure (which sits close to CSR, 0.99–1.31×), the field-hotspot geometry at bw=6µm has a *shorter* median distance and a *thinner* far tail than a random point field at the same blob density (0.77× / 0.44×) — the many small, sharp single-cell spikes this bandwidth produces are more evenly spread than CSR would predict, not more clustered. It does fix the degenerate 0 µm the 40 µm default gives (thresholding a field that's smooth over ~97% of tissue trivially selects contiguous cores everywhere), landing at a plausible-looking 41 µm median — but that is a side-effect of un-smoothing the field back toward the near-CSR raw cell pattern, not evidence of recovered vascular heterogeneity. The 500 µm-tile field heterogeneity (2.50×) is essentially unchanged from bw=14's 2.82× and from the cell-count-based 2.4× at any bandwidth — all three orders of magnitude short of the ~48× CD34 area-fraction spread across real GBM cases (PMID 38928515) that motivates this sweep.

## Downstream sanity — cisplatin, dx=20µm (rough, not grid-converged)

| metric | bw=6µm @ dx=20 | discrete-mode default (dx=12) | bw=40µm KDE (dx=12) | source |
|---|---:|---:|---:|---|
| mean tumour AUC (µM·s) | 1.058e6 | 1.69–1.75e6 | 5.758e6 | `real_drug_panel.md` / `delivery_metrics.md` / `vessel_promotion.md` KDE row |
| therapeutic_coverage@1×IC50 | 0.968 | 0.9936 | 0.9986 | `delivery_metrics.md` / `vessel_promotion.md` KDE row |
| kill_fraction | 0.328 | 0.375–0.436 | 0.912 | `vessel_promotion.md` (discrete/KDE rows) / `real_drug_panel.md` |
| transport field NaN/Inf | 0 / 0 | — | — | sane |

IC50 used: 24.55 µM (matches the panel's 24.6 µM). No NaN, no degenerate field; the solve is numerically sound.

## Bottom line

The field is numerically sound at 6 µm (no NaN/Inf/degenerate), but tightening the kernel this far does **not** buy the clustering this sweep is chasing — it undershoots CSR on the field-hotspot distance test (0.77×/0.44×, i.e. *more* regular than random) rather than overshooting it, and its 500 µm-tile heterogeneity (2.50×) is statistically indistinguishable from bw=14µm's (2.82×) and the bandwidth-invariant cell-count number (2.4×), all far short of real GBM's ~48× case-to-case spread. Downstream, cisplatin delivery collapses relative to both baselines — mean AUC falls to 1.06e6 µM·s (0.6× the current discrete-mode default, 0.18× the legacy 40 µm KDE) and kill_fraction to 0.328 (below the discrete default's 0.37–0.44 and far below the 40 µm KDE's 0.91) — though this comparison is confounded by dx (20 µm here vs 12 µm for both baselines) as well as by field construction, so read it as "same order of magnitude, lower delivery," not a precise ratio. Net: 6 µm looks too aggressive — it fixes the degenerate 0 µm field-hotspot distance the 40 µm default gives, but by regularising the field toward the near-CSR raw cell pattern rather than by revealing real mesoscale clustering, and it pays for that with a large drop in predicted delivery.
