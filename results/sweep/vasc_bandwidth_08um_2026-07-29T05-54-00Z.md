# Bandwidth sweep point: 8 µm

One point in a 20-point KDE-bandwidth sweep (baseline: 40 µm, `real_fields._kde_on_grid`). Field built in-memory via `vascular_density_from_cells(domain, cells, bandwidth_um=8.0, mode='kde')` on `reconstruct_real('xenium', dx_um=12)`; no source file modified. Reproduce: `./.venv/bin/python <script above>`.

| check | bw=8 µm | bw=40 µm (default, same run) |
|---|---:|---:|
| sanity (NaN / all-zero / range) | 0 NaN, not all-zero, [0,1] — **PASS** | — |
| zero-fraction of tissue mask | **45.2%** | 0.2% |
| mean field value (tissue mask) | 0.0181 | 0.0747 |
| 500 µm-tile heterogeneity, p90/p10 of *field intensity* | **2.58×** | 2.63× |
| tumour→nearest-nonzero-pixel distance, median / p90 | 0 / **38 µm** | 0 / 0 µm |
| nonzero-tissue fraction | 54.8% | 99.8% |

Cell-based CSR diagnostic (bandwidth-**invariant** by construction — DBSCAN-clustered endothelial cell coordinates, not the KDE field; reproduces `vascular_validation.md`'s own numbers as a sanity check on the reused code): median tumour→vessel distance is **0.99×** the Poisson-CSR-null at the same MVD (110.8/mm²), tail beyond 100 µm is **1.44×** the null (report's own mixed convention gives 1.31× — both self-consistent, same conclusion). Endothelial-cell-density heterogeneity across 500 µm tiles: **2.42×** p90/p10.

Downstream sanity check — one rough cisplatin solve, dx=20 µm (not grid-converged), bw=8 field substituted for `vascular_density` before `simulate_response`:

| metric | this run (bw=8, dx=20, unconverged) | published (dx=12) |
|---|---:|---:|
| mean tumour AUC (µM·s) | 1.16e+06 | 1.75e+06 (`real_drug_panel.md`) / 1.69e+06 (`delivery_metrics.md`) |
| therapeutic_coverage @1×IC50 | 0.984 | 0.9936 (`delivery_metrics.md`; not reported in `real_drug_panel.md`) |
| kill_fraction / penetration_index | 0.341 / 0.421 | 0.4362 / 0.428 (`real_drug_panel.md`, `mode='uniform'`) |

Field is numerically sane (no NaN, no degeneracy, bounded [0,1], not all-zero). Narrowing the bandwidth to 8 µm makes a real, large difference to the field's own *footprint*: the zero-fraction of tissue jumps from 0.2% (40 µm) to 45.2% (8 µm), and the tumour→nonzero-field p90 distance opens from 0 to 38 µm — so the "smooth-everywhere hotspot" artifact documented in `vascular_validation.md` §5 is substantially reduced at this bandwidth. But the **500 µm-tile heterogeneity is essentially unchanged** (2.58× vs 2.63×, both ≈ the cell-based 2.42×): tile-scale heterogeneity is set by where the endothelial cells actually are, not by the smoothing kernel, so bandwidth alone will not fix the "no vascular hotspots / no avascular zones" verdict at that scale — only the sub-tile structure changes. The cell-based CSR result (~1.0× median, ~1.3–1.4× tail) is reproduced unchanged, confirming it is genuinely bandwidth-independent, as the diagnostic never touches the KDE field. **Caveat on the transport comparison**: by file mtime, both `real_drug_panel.md` (2026-07-28 19:30) and `delivery_metrics.md` (2026-07-28 12:48) post-date `vessel_promotion.md` (2026-07-28 04:18), which flipped the pipeline default from `mode='kde'` (40 µm) to `mode='discrete'` — so those published numbers most likely reflect the **discrete** vessel field, not a 40 µm KDE, meaning this transport check compares bandwidth *and* field construction at once (plus an unconverged dx=20 vs the published dx=12), not bandwidth alone; read the ~0.67× AUC gap and the near-identical, near-saturated coverage (0.984 vs 0.9936) accordingly.
