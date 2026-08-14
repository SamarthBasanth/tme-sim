# Vascular KDE bandwidth sweep — 22 µm

One point in a ~20-point bandwidth sweep (context: `results/vascular_validation.md`, `tmesim/reconstruct/real_fields.py::_kde_on_grid`, default `_DEFAULT_BANDWIDTH_UM = 40.0`). Field built in-memory via `vascular_density_from_cells(domain, cells, bandwidth_um=22.0, mode='kde')` on `reconstruct_real('xenium', dx_um=12)`; no existing file edited. A `bandwidth_um=40.0` KDE field was built the same way, in the same run family, as a directly comparable within-script baseline (`real_fields`'s current pipeline default is `mode='discrete'`, not KDE, so the published KDE numbers elsewhere need the same caveat other sweep points already flagged — see below).

## Field sanity (dx=12µm build)

| check | value |
|---|---|
| NaN / Inf | 0 / 0 |
| range | [0.0, 1.0] |
| mean over tissue | 0.0451 |
| nonzero fraction over tissue | 95.64% (vs 99.78% measured here for bw=40, cross-checks the ~99.7% the `real_fields.py` docstring claims for the 40µm default) |
| all-zero? | no |

## Clustering diagnostic

Two families, as in every other point of this sweep. The **cell-coordinate** diagnostic (vessel clustering + CSR null, endothelial-count tile heterogeneity) is computed from raw endothelial-cell positions, never from the KDE field, so it is **bandwidth-invariant by construction** and reproduces `vascular_validation.md` exactly — not recomputed here, just carried through. The bandwidth-*sensitive* signal is the KDE field's own spatial structure.

| diagnostic | bw=22µm | bw=40µm (this run) | note |
|---|---:|---:|---|
| median tumour→vessel dist (cell-based) | 44.2 µm | 44.2 µm | identical — bandwidth-invariant (`vascular_validation.md`) |
| ratio vs CSR null (median) | 0.99× | 0.99× | identical, matches `vascular_validation.md` |
| endothelial-count heterogeneity, p90/p10 (cell-based) | 2.4× | 2.4× | identical, published value, not recomputed |
| **field-value heterogeneity, p90/p10 (500µm tiles)** | **2.55×** | **2.55×** | essentially unchanged — 500µm tiles are ≫22–40µm, so this bandwidth difference washes out at that scale; matches the bw=16µm sweep point's bw=40 control exactly (2.55×), cross-validating the method; both far below CD34 area-fraction's ~48× spread across cases |
| field hotspot distance @p95 threshold, median (p90) | **80.5 µm (181.6)** | **136.8 µm (303.6)** | n_components 517 vs 191 at matched q=18.5/mm² vs 6.9/mm² |
| ratio vs CSR null at matched hotspot density | **0.738×** | **0.762×** | both **<1×** — the percentile-selected "hotspots" sit *closer together* than a random field of the same density would, not further apart. Narrowing the kernel fragments large smoothed blobs into more, smaller, closely-spaced local maxima; it is not evidence of recovered clustering (real clustering would push this ratio *above* 1×). |

## Downstream sanity — cisplatin, dx=20µm (rough, not grid-converged)

| metric | bw=22µm @ dx=20 | bw=40µm @ dx=20 (control, this run) | 40µm-KDE reference, dx=12 (`vessel_promotion.md`, cisplatin) | discrete-default, dx=12 (`real_drug_panel.md`, cisplatin) |
|---|---:|---:|---:|---:|
| mean tumour AUC (µM·s) | 3.084e6 | 5.074e6 | 5.758e6 | 1.75e6* |
| therapeutic_coverage@1×IC50 | 0.9986 | 0.9988 | 0.9986 | not reported there |
| kill_fraction | 0.878 | 0.929 | 0.912 | 0.4362 |
| transport field NaN/Inf | 0 / 0 | 0 / 0 | — | — |

\*`real_drug_panel.md`'s "mean AUC" column is the **tissue-wide** mean (`res.extras['mean_auc']`), not the tumour-cell mean quoted above, and it runs under the CURRENT pipeline default (`mode='discrete'`, segmented vessel objects), not a 40µm KDE, since `results/vessel_promotion.md` flipped the default. It is not a valid same-mode comparator for a bandwidth sweep; `vessel_promotion.md`'s own KDE-arm cisplatin row (dx=12, 5.758e6/0.9986/0.912) is the better reference and sits close to this run's bw=40µm/dx=20 control, so the residual gap there is attributable to dx=20 vs dx=12 rather than to methodology.

## Bottom line

The bw=22µm field is numerically sound (no NaN/Inf, not degenerate, 95.6% nonzero over tissue). As with every other sweep point, the module's headline "close to random" verdict (median tumour→vessel distance 0.99× the CSR null) is computed from raw endothelial-cell coordinates and cannot move with bandwidth. Of the two bandwidth-*sensitive* field diagnostics tested, neither shows 22µm recovering real clustering relative to 40µm: 500µm-tile heterogeneity is unchanged (2.55× at both, still ~19× short of the ~48× CD34-area-fraction spread across real GBM cases), and the percentile-threshold hotspot-distance ratio to its own matched-density CSR null actually stays below 1× at both bandwidths (0.738× at 22µm vs 0.762× at 40µm) — narrowing the kernel makes the selected hotspots more numerous and closer together, not more separated/clustered. Downstream on a rough dx=20µm cisplatin solve, narrowing 40→22µm lowers predicted delivery and kill (mean tumour AUC 5.07e6→3.08e6 µM·s, kill 0.929→0.878), while `therapeutic_coverage@1×IC50` stays pinned near ceiling (0.9988→0.9986) and is not a useful discriminator for this drug. Net for the synthesis pass: 22µm does not fix the heterogeneity gap by either field-based measure used here, and its main visible effect — like 16µm's — is redistributing the field's mass into more, smaller local peaks rather than into larger clustered/avascular regions; weigh the 500µm-tile heterogeneity ratio over the hotspot-distance ratio, and discount `real_drug_panel.md`'s discrete-mode numbers as a bandwidth comparator in favor of `vessel_promotion.md`'s KDE-arm figures.
