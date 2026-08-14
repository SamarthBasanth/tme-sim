# Vascular KDE bandwidth sweep — 26 µm

One point in a ~20-point bandwidth sweep (context: `results/vascular_validation.md`, `tmesim/reconstruct/real_fields.py::_kde_on_grid`, default `_DEFAULT_BANDWIDTH_UM = 40.0`; current pipeline default vascular mode is actually `'discrete'`, not KDE — see caveat below). Field built in-memory via `vascular_density_from_cells(domain, cells, bandwidth_um=26.0, mode='kde')` on `reconstruct_real('xenium', dx_um=12)`; no existing file edited. A `bandwidth_um=40.0` KDE field was built the same way, in the same run, as a directly comparable within-script baseline.

## Field sanity (dx=12µm build)

| check | bw=26µm | bw=40µm (this run) |
|---|---:|---:|
| NaN / Inf | 0 / 0 | 0 / 0 |
| range | [0.0, 1.0] | [0.0, 1.0] |
| mean over tissue mask | 0.0508 | 0.0747 |
| nonzero fraction over tissue mask (193,752 px) | 98.9% | 99.8% |
| all-zero? | no | no |

Both fields are numerically sane — no degeneracy, no NaN/Inf, full [0,1] range after max-normalisation. Narrowing the kernel from 40→26µm leaves the field nonzero over almost all of tissue either way (98.9% vs 99.8%): at dx=12µm with ~40,887 measured cells over 27.9 mm², the median inter-cell gap (~16.7µm) is well under either bandwidth, so the KDE stays essentially fully connected at both settings — it does not fragment into isolated blobs at 26µm.

## Clustering diagnostic

Two families. The module's **cell-coordinate** diagnostic (DBSCAN vessel clustering + CSR null, endothelial-count tile heterogeneity) is computed straight from raw endothelial-cell positions — it never touches the KDE field — so it is **bandwidth-invariant by construction** and reproduces `vascular_validation.md` exactly regardless of which bandwidth is tested. The bandwidth-*sensitive* signal is the KDE field's own spatial structure, reported separately below.

| diagnostic | value | note |
|---|---:|---|
| median tumour→vessel distance (cell-based, link=20µm) | 44.21 µm | identical at every bandwidth |
| CSR (Poisson) null at same MVD (110.8/mm²) | 44.63 µm | — |
| ratio vs CSR null | **0.991×** | matches `vascular_validation.md`'s reported 0.99× |
| endothelial-count heterogeneity, p90/p10 (500µm tiles, cell-based) | **2.42×** | matches published 2.4× |

| field diagnostic (bandwidth-sensitive) | bw=26µm | bw=40µm (this run) |
|---|---:|---:|
| field hotspot distance @p95 threshold, median (p90), µm | 93.7 (205.1) | 136.8 (303.6) |
| hotspot area (5% of tissue by construction — percentile threshold) | 1.395 mm² | 1.395 mm² |
| ratio vs cell-based CSR null (44.63 µm) | 2.10× | 3.07× |
| tile-mean-field heterogeneity, p90/p10 (500µm tiles) | **2.62×** | **2.63×** |

**Correction on the field-distance baseline.** `vascular_validation.md` §5 reports a degenerate "0 µm" field-based distance, but that number was computed against the pipeline's *default* `vascular_density` field, which is `mode='discrete'` in the current codebase (exactly zero over most of tissue by construction, so a p95 threshold on it behaves very differently from a true KDE). It is not a genuine 40µm-KDE figure despite the module's docstring describing the field that way. Both rows above were built explicitly with `mode='kde'` in this script and give non-degenerate numbers — use the bw=40 column here as the reference, not the published "0 µm".

**Reading the bandwidth effect honestly.** Narrowing 40→26µm shrinks the field-hotspot distance substantially (136.8→93.7µm median, i.e. the top-5%-density patches sit closer to tumour cells), because a tighter kernel produces smaller, more numerous local peaks rather than a few broad blobs — this is a real, bandwidth-driven change in the field's fine structure. But at the 500µm tile scale used for the headline heterogeneity number, the effect nearly vanishes (2.62× vs 2.63×): 500µm is >12× either bandwidth, so total density mass per tile is roughly conserved under blurring at these kernel widths, and neither setting comes close to the ~48× CD34-area-fraction spread measured across real GBM cases (PMID 38928515). **Bandwidth=26µm changes the field's local geometry but does not fix the tile-scale heterogeneity gap this validation flags** — and it cannot move the cell-based 0.99×-CSR verdict at all, since that number is computed upstream of any KDE.

## Downstream sanity — cisplatin, dx=20µm (rough, not grid-converged)

| metric | bw=26µm @ dx=20 (this run) | bw=40µm @ dx=20 (matched, this run) | discrete-default, dx=12 (`real_drug_panel.md`) | KDE(40µm), dx=12 (`vessel_promotion.md` §2.1) |
|---|---:|---:|---:|---:|
| mean tumour AUC (µM·s) | 3.498e6 | 5.074e6 | 1.75e6 | 5.758e6 |
| therapeutic_coverage@1×IC50 | 0.9987 | 0.9988 | not reported | 0.9986 |
| kill_fraction | 0.8934 | 0.9290 | 0.4362 | 0.9120 |
| transport field NaN/Inf | 0 / 0 | 0 / 0 | — | — |

The bw=26µm field is numerically well-behaved end to end: no NaN/Inf in the AUC field, finite bounded output, no crash. My own bw=40µm dx=20 run (5.074e6 AUC, coverage 0.9988, kill 0.929) sits close to the published dx=12 40µm-KDE reference (5.758e6, 0.9986, 0.9120), cross-checking this script's methodology; the residual gap is attributable to dx=20 vs dx=12 (deliberately unconverged here, as instructed), not to any bandwidth error. Against that matched pair, narrowing to 26µm lowers predicted cisplatin delivery (AUC 5.074e6→3.498e6 µM·s, kill 0.929→0.893) while `therapeutic_coverage@1×IC50` stays pinned near ceiling at both (0.9988 vs 0.9987) and is not a useful discriminator here. The comparison against `real_drug_panel.md`'s discrete-default numbers (1.75e6 AUC, kill 0.4362) is **confounded three ways** — KDE vs discrete field construction (real_fields.py documents a ~3.36× normalisation difference between the two source conventions), dx=20 (rough) vs dx=12 (converged), and possibly small alpha/CFL differences at the coarser mesh — so it should not be read as evidence about bandwidth specifically; the bw=40-matched column is the correct like-for-like anchor.

## Bottom line

The bw=26µm field is numerically sound (no NaN/Inf, [0,1] range, 98.9% nonzero over tissue, not degenerate) and runs cleanly through a full transport solve. As at every other point in this sweep, the module's headline "close to random" finding — median tumour→vessel distance 0.99× the CSR null — is computed from raw endothelial-cell coordinates and **cannot move with bandwidth at all**; that verdict is unaffected by 26µm or any other kernel width tested here. Among the diagnostics that *are* bandwidth-sensitive, 26µm meaningfully tightens the field's local hotspot structure relative to 40µm (median hotspot distance 93.7 vs 136.8µm) but leaves the 500µm-tile heterogeneity essentially unchanged (2.62× vs 2.63×, both far short of the ~48× real cross-case CD34 spread) and lowers downstream cisplatin delivery in the matched dx=20 pair (AUC 3.498e6 vs 5.074e6 µM·s, kill 0.893 vs 0.929). Net: 26µm does not resolve the heterogeneity gap `vascular_validation.md` flags, buys a real but modest improvement in local hotspot compactness, and costs some predicted delivery relative to the 40µm default — the synthesis pass should weigh the tile-heterogeneity ratio (the diagnostic closest to the report's own §3.6 metric) as the primary ranking signal across bandwidths, not the hotspot-distance number, which is a threshold-selection artifact as much as a clustering measurement.
