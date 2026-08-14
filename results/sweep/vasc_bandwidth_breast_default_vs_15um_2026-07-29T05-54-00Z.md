# Does the GBM "vasculature looks random" finding generalize to breast? (40µm default vs 15µm candidate)

`reconstruct_real(source="xenium", sample="data_cache/xenium_breast_janesick/rep1", dx_um=12.0)` — the exact call `results/breast_cancer_flagship.md` used (167,780 cells, 89,785 tumor, 13,677 endothelial via `BREAST_PANEL_MARKER_MAP` = PECAM1/VWF/CLDN5). Current `real_fields._DEFAULT_BANDWIDTH_UM = 40.0`; `vascular_density_from_cells(..., mode='kde')` compared at 40µm (current default) vs 15µm (candidate). Not run through the repo — a standalone research script (not committed).

| diagnostic | breast @ default (40µm) | breast @ 15µm | GBM published (`vascular_validation.md`) |
|---|---:|---:|---:|
| median tumor→vessel dist, **raw cell positions** (bandwidth-invariant) | 31.2 µm | 31.2 µm (unchanged) | 42.8–44.2 µm |
| ratio vs CSR null (median, raw) | **1.25×** | 1.25× (unchanged) | **0.99×** |
| fraction tumor cells >100µm from a vessel (raw) | **4.17%** | 4.17% (unchanged) | **4.02%** |
| endothelial cell-count p90/p10, 500µm tiles (raw, bandwidth-invariant) | 6.87× | 6.87× (unchanged) | 2.42× |
| field-threshold (≥0.3×max) occupied area | 3.81% | 1.69% | — (n/a, different method) |
| median dist to field-threshold "vessel" pixel vs CSR | 5.55× | 2.77× | — |
| field-value p90/p10, 500µm tiles | 6.32× | 5.95× | 2.42× (bw40) → 2.82× (bw14) |
| p95-hotspot median distance | 135.1 µm | 73.2 µm | ~0 µm (bw40, degenerate) → 37.9 µm (bw14) |

**Yes, the same problem generalizes.** The diagnostic that actually measures real-vasculature clustering — nearest-vessel distance from raw, measured endothelial-cell positions vs. a CSR (complete-spatial-randomness) null at the same density — is by construction bandwidth-invariant, and on breast it gives ratio 1.25× and a >100µm tumor tail of 4.17%, essentially the same near-random signature GBM reported (0.99×, 4.02%). This holds despite breast's richer 13-endothelial-gene Xenium panel; the production marker map (`marker_panels.BREAST_PANEL_MARKER_MAP`) only actually types ENDOTHELIAL from 3 of those 13 genes (PECAM1/VWF/CLDN5), so the extra panel richness isn't yet exploited by vessel identification and doesn't explain or fix the near-random result. Whether 15µm "helps" is inconclusive and possibly moot here: the bandwidth-sensitive field diagnostics (threshold-occupied area, p95-hotspot distance) do shrink going 40µm→15µm, in the same direction the parallel GBM sweep reports, but the field's own tile-level p90/p10 heterogeneity barely moves (6.32×→5.95×, essentially flat, unlike GBM's reported rise) — plausibly because the underlying raw cell positions are already this close to CSR-random, so there is little real clustering signal left for any bandwidth to reveal. A couple of quick searches found no breast-cancer equivalent of the Marien et al. 2016 whole-slide stereological MVD anchor used for GBM (that paper covers GBM/colorectal/renal/ovarian, not breast); one secondary hit mentioned breast "immature vasculature clusters" near ~42±8 CD31+ vessels/mm² co-locating with perivascular hypoxia hotspots, but it wasn't traced to primary literature in the time budgeted, so it is noted only, not relied on.
