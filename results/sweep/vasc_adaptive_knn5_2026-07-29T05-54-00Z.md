# Adaptive (k-NN) bandwidth KDE — does it recover the clustering a fixed bandwidth erases?

Follow-up to `results/vascular_validation.md` and the fixed-bandwidth sweep in this directory
(`vasc_bandwidth_{08,12,14,16}um_*.md`, which already found that narrowing the *fixed* Gaussian
KDE bandwidth does not fix the 500 µm-tile heterogeneity gap). This point tests the other
well-established KDE variant instead: **adaptive (variable) bandwidth**, scaled per-cell by local
k-nearest-neighbour spacing. Built entirely in a standalone script — `real_fields.py` was not
modified — on `reconstruct_real('xenium', dx_um=12)` (4,295 real endothelial cells, grid
(424, 675)). No repo file was edited; no test suite run.

## Formula

For each endothelial cell *i*, let `d_k(i)` be the Euclidean distance to its k-th nearest
**other** endothelial cell (k = 5). The local bandwidth is

```
h_i = c · d_k(i),   c = 0.5, k = 5
```

and the field is the sum of per-cell Gaussians, each using its own bandwidth, then masked to
tissue and max-normalised — identical post-processing convention to `real_fields._kde_on_grid`,
so the two fields are a like-for-like comparison:

```
vascular_density(x, y) = clip( Σ_i exp( -|(x,y) - (x_i,y_i)|² / (2·h_i²) )  ·  tissue_mask(x,y)  /  max(·),  0, 1 )
```

This is the simplest member of the **k-NN adaptive / "balloon estimator" bandwidth family**
(Breiman, Meisel & Purcell 1977, *Technometrics* 19:135–144; see also Silverman, *Density
Estimation for Statistics and Data Analysis*, 1986, §5.3, for the general variable-kernel
construction) — bandwidth directly proportional to k-NN spacing, no pilot-density iteration.
This is **not** the Abramson (1982) two-pass estimator, which instead sets `h_i` from a pilot
density evaluated at the point with a fixed −1/2 sensitivity exponent; that variant was not
implemented here (flagged per the deviation instruction).

`h` on this section: median **45.5 µm**, mean 46.1, p10 17.5, p90 **72.6**, max **205.2 µm** —
i.e. mostly *larger* than the fixed 40 µm bandwidth it's compared against, and heavily
right-skewed.

## Results

| diagnostic | fixed 40 µm KDE | adaptive k=5, c=0.5 KDE | bandwidth-sensitive? |
|---|---:|---:|---|
| cell-based median tumour→vessel dist. vs CSR null | 0.99× | 0.99× (identical) | **no** — computed from raw endothelial-cell / DBSCAN-vessel coordinates, never reads `vascular_density` |
| cell-based tail >100 µm from a vessel | 4.4% | 4.4% (identical) | **no** |
| endothelial-cell-count tile heterogeneity, p90/p10 (500 µm, ≥100 cells/tile) | 2.42× | 2.42× (identical) | **no** |
| **field-value** tile heterogeneity, p90/p10 (500 µm, sampled at cell positions, ≥100 cells/tile) | **2.78×** | **1.27×** | yes — **adaptive is WORSE** |
| field coefficient of variation over tissue | **1.005** | **0.403** | yes — **adaptive is WORSE** |
| zero-fraction of tissue | 0.2% | 0.04% | adaptive is *more* uniformly nonzero |
| sanity (NaN/Inf, range) | 0/0, [0,1] | 0/0, [0,1] | **PASS**, both fields numerically sane |

`c`/`k` sensitivity (same section, field-value tile p90/p10, no transport re-run): every
combination of k ∈ {3, 5, 10} × c ∈ {0.15, 0.25, 0.5, 1.0} gives **1.5×–1.8×**, i.e. below the
fixed-bandwidth family's 2.4×–3.4× (this run and `vasc_bandwidth_{08,12,14,16}um`) at *every*
setting tested — this is not a c = 0.5 / k = 5 artifact.

### Downstream sanity check — cisplatin, dx = 20 µm (rough, not grid-converged)

Same-run comparator (fixed-40 µm KDE built and solved in the same script, same dx = 20 grid) is
used as the trustworthy baseline, since `real_drug_panel.md`'s published numbers now come from
the unrelated `mode='discrete'` field (per `vessel_promotion.md`), not a 40 µm KDE.

| metric | adaptive k=5,c=0.5 | fixed 40 µm (same run) | published discrete-mode, dx=12 (`real_drug_panel.md`) |
|---|---:|---:|---:|
| mean tumour AUC (µM·s) | **2.54e+07** | 5.07e+06 | 1.75e+06 |
| therapeutic_coverage@1×IC50 | 0.9993 | 0.9988 | not reported there (0.9936 under discrete, `delivery_metrics.md`) |
| kill_fraction | 0.9896 | 0.9290 | 0.4362 |
| transport field NaN/Inf | 0/0 | 0/0 | — |

## Bottom line

Adaptive k-NN bandwidth does **not** preserve clustering better than a fixed bandwidth on this
section — it is worse by every field-level measure tested (CV 0.40 vs 1.01, tile p90/p10 1.27×
vs 2.78×), robustly across a full c/k sweep. The mechanism is specific to this data: `vascular_validation.md` §3.2 already found that 81% of "putative vessels" here are single isolated
endothelial cells, so most cells' 5th-nearest-*other*-endothelial-cell distance is large and
right-skewed (median 45.5 µm, p90 72.6, max 205 µm) — instead of shrinking in dense clusters and
expanding only in genuinely sparse gaps, the adaptive bandwidth mostly just inflates across the
(majority) sparse regions, flooding them with broad Gaussian blobs that erase the very avascular
pockets it was meant to preserve. The cell-position CSR/heterogeneity diagnostics
`vascular_validation.py` actually implements are unaffected by any of this — they never read
`vascular_density` at all, fixed or adaptive, so the section's headline "0.99× CSR" verdict is a
property of the raw endothelial point pattern, not fixable by *any* smoothing-kernel choice.
Downstream, the adaptive field is numerically sane (no NaN/Inf, bounded [0,1], finite
kill/coverage) but delivers ~5× more drug than the matched fixed-40 µm run and ~14× more than the
published discrete-mode panel, because its flooded, near-uniform geometry removes the
perivascular boundary layer the transport model relies on — on this real GBM section, k-NN
adaptive bandwidth is a worse fix than doing nothing, not a better one.

---
*Reproduce:* standalone script (not checked into this repo), calling only existing, unmodified
`tmesim` functions: `reconstruct.real_pipeline.reconstruct_real`,
`reconstruct.real_fields.vascular_density_from_cells(..., mode='kde')` and `add_oxygen`,
`validation.vascular_validation.{cluster_vessels, distance_to_nearest, local_heterogeneity,
poisson_distance_stats}`, `validation.calibration.{calibrate_alpha, calibrated_drugparams}`,
`validation.delivery_metrics.{tumor_cell_exposure, therapeutic_coverage}`,
`sim.engine.simulate_response`.
