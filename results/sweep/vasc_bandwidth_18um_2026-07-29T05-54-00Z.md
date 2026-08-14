# Vascular KDE bandwidth sweep — 18 µm

One point in a ~20-point bandwidth sweep (context: `results/vascular_validation.md`, `tmesim/reconstruct/real_fields.py::_kde_on_grid`, default `_DEFAULT_BANDWIDTH_UM = 40.0`). Field built in-memory via `vascular_density_from_cells(domain, cells, bandwidth_um=18.0, mode='kde')` on `reconstruct_real('xenium', dx_um=12)`; no existing file edited. A `bandwidth_um=40.0` KDE field was built the same way, in the same run, as a directly comparable within-script control — the pipeline's *actual current default* is `mode='discrete'` (since `vessel_promotion.md`), so `results/real_drug_panel.md` / `results/delivery_metrics.md` are **not** true 40 µm-KDE baselines despite older docstrings describing the field that way; `results/vessel_promotion.md` §2.1's cisplatin dx=12 KDE ladder (coverage 0.9986 · AUC 5.758e6 · kill 0.9120) is the one legitimate external KDE anchor and is used below to cross-check this run's own bw=40 control.

## Field sanity (dx=12µm build)

| check | bw=18µm | bw=40µm (matched control) |
|---|---:|---:|
| NaN / Inf | 0 / 0 | 0 / 0 |
| range | [0.0, 1.0] | [0.0, 1.0] |
| mean over tissue | 0.0392 | (see below) |
| nonzero fraction over tissue | 91.3% | 99.8% (docstring claims ~99.7% for the default — matches) |
| all-zero? | no | no |

## Clustering diagnostic

Two families. The module's **cell-coordinate** diagnostic (vessel clustering + CSR null, endothelial-count tile heterogeneity) is computed from raw endothelial-cell positions and never touches the KDE field, so it is **bandwidth-invariant by construction**. Reproduced here using `vascular_validation.py`'s own `verdict()` definition — median ratio from vessel-CENTROID distances, tail ratio from raw-endothelial-CELL distances against a CSR null at the vessel-object density (a real asymmetry in the published methodology, matched exactly rather than approximated):

| diagnostic | bw=18µm | bw=40µm (this run) | published (`vascular_validation.md`) |
|---|---:|---:|---:|
| median tumour→vessel dist (cell-based) | 44.2 µm | 44.2 µm | 44 µm |
| ratio vs CSR null (median) | 0.991× | 0.991× | 0.99× |
| ratio vs CSR null (>100µm tail) | 1.307× | 1.307× | 1.31× |
| endothelial-count heterogeneity, p90/p10 (cell-based) | 2.42× | 2.42× | ~2.4× |
| MVD (link=20µm) /mm² | 111 | 111 | 111 |

The bandwidth-*sensitive* signal is the KDE field's own spatial structure:

| diagnostic | bw=18µm | bw=40µm (this run) | note |
|---|---:|---:|---|
| field-value heterogeneity, p90/p10 (500µm tiles) | 5.00× | 5.19× | essentially unchanged — 500µm tiles are ≫18–40µm, so the kernel-width difference mostly washes out at that scale; both far below CD34 area-fraction's ~48× spread across published cases |
| hotspot CSR median-ratio @ p90/p95/p98 threshold | 0.908 / 0.951 / 1.013 | 0.949 / 1.001 / 1.105 | all within ±11% of 1.0 (CSR) at every threshold, both bandwidths — no clean "heavier tail" clustering signature at either |
| hotspot CSR tail-ratio (>100µm) @ p90/p95/p98 | 0.552 / 0.929 / 0.985 | 0.967 / 0.985 / 0.993 | bw=18 is *more sub-Poisson* (more regular, not more clustered) than bw=40 at the loosest threshold |
| n hotspot objects @ p95 / median inter-hotspot distance | 697 objects, 89.4 µm | 191 objects, 179.7 µm | narrowing the kernel fragments large blobs into more, closer, smaller local maxima — a threshold-selection artefact, not recovered biological clustering |

## Downstream sanity — cisplatin, dx=20µm (rough, not grid-converged)

| metric | bw=18µm @ dx=20 | bw=40µm @ dx=20 (this run) | 40µm-KDE reference, dx=12 (`vessel_promotion.md` §2.1) | discrete-default, dx=12 (`real_drug_panel.md`) |
|---|---:|---:|---:|---:|
| mean tumour AUC (µM·s) | 2.658e6 | 5.074e6 | 5.758e6 | 1.75e6 |
| therapeutic_coverage@1×IC50 | 0.9983 | 0.9988 | 0.9986 | not reported there (`delivery_metrics.md` discrete: 0.9936) |
| kill_fraction | 0.804 | 0.929 | 0.912 | 0.4362 |
| transport field NaN/Inf | 0 / 0 | 0 / 0 | — | — |

This run's own bw=40 control (5.07e6 AUC, kill 0.929, coverage 0.9988) sits close to the external `vessel_promotion.md` dx=12 KDE reference (5.76e6, 0.912, 0.9986), cross-checking the script's methodology; the residual gap is attributable to dx=20 vs dx=12 (deliberately not grid-converged) rather than to anything bandwidth-related.

## Bottom line

The bw=18µm field is numerically sound (no NaN/Inf, not degenerate, 91.3% nonzero over tissue). As at every other sweep point, the module's headline "close to random" finding (median tumour→vessel distance 0.99× the CSR null, tail 1.31×) is computed from raw endothelial-cell coordinates and cannot move with bandwidth — reproduced here to 3 significant figures and identical at bw=18 and the bw=40 control. The genuinely bandwidth-sensitive field-hotspot diagnostic gives a mixed, non-monotone picture: at 18µm the field is somewhat less smeared than at 40µm (91.3% vs 99.8% nonzero, tile heterogeneity 5.00× vs 5.19×), but its percentile-threshold "hotspots" stay within ±11% of the CSR null at every threshold tested for both bandwidths — narrowing to 18µm does **not** recover the heavier-than-CSR far tail that would signal real vascular clustering, and if anything the 18µm hotspots read slightly *more* regular (sub-Poisson) than the 40µm ones at loose thresholds, because a narrower kernel fragments large smoothed blobs into more numerous, more closely-spaced small peaks rather than into fewer, farther-apart true clusters. Downstream, narrowing 40→18µm nearly halves predicted cisplatin delivery in this run's matched pair (mean AUC 5.07e6→2.66e6 µM·s, kill 0.929→0.804), while `therapeutic_coverage@1×IC50` stays pinned near ceiling (~0.998–0.999) at both bandwidths and both are still roughly 2× the converged discrete-mode kill (0.436) — a mode effect, not evidence for or against this bandwidth. Net: 18µm does not resolve the clustering/heterogeneity gap `vascular_validation.md` flags; the synthesis pass should track the 500µm-tile heterogeneity ratio and the hotspot-CSR-ratio spread (not the published discrete-mode drug numbers, which are the wrong baseline) when ranking candidate bandwidths.
