# Vessel promotion (`ceb96c0`) — inputs and provenance

Filed after the fact. The promotion of the discrete segmented-vessel source to
the default vascular field was decided and committed; this records **what it was
run on**, which had never been written down. The raw drivers lived in a
session-scoped temp directory and **have already been deleted** — they were not
recoverable by the time this was filed, so what follows is the record, not a
re-derivation. Anything below that matters should be re-measured rather than
cited.

## Data

Real 10x Xenium GBM section via `reconstruct_real('xenium')`: 40,887 measured
cells, 30,174 tumour, 3,984 `CellType.ENDOTHELIAL`. Tissue mask `'geometric'`
(union of 30 µm discs). Grids 339×541 (dx=15), 424×675 (dx=12), 635×1012
(dx=8), 846×1349 (dx=6); tissue pixels 124,074 / 193,752 / 435,924 / 775,035.

## Solver config

`dt = 5.0 s` **fixed at every dx**, not the historical 10 s. `total_time`
21,600 s (6 h), `convection = True`, `seed = 0`. Engine, core and physics
unmodified.

Stated reason for the reduced `dt`: diffusion and reaction are implicit and
unconditionally stable, so only the explicit upwind advection is CFL-bound.
Holding `dt` fixed keeps it unconfounded with `dx`.

## Parameters

`alpha = 148.60815765445545` from `calibrate_alpha()`, unmodified.
`LINK_DISTANCE_UM = 20.0`, `VESSEL_RADIUS_UM = 7.5`, `amplitude = 'unit'`.
BBB `L_PARACELLULAR_LEAK = 0.238`; `REGION_DISRUPTION` core 0.50 / margin 0.00.
Convergence criteria taken from `delivery_metrics.py` rather than invented:
`CONVERGENCE_TOL_ABS = 0.03`, `CONVERGENCE_TOL = 0.05`, `DISCRIMINATION_MIN_RATIO = 3.0`.

## Drugs

Convergence set: temozolomide, cisplatin, doxorubicin (doxorubicin is the
resolution-limited case, shortest L). Panel adds lomustine, methotrexate.
Out-of-sample: mitoxantrone, topotecan — compiled from verified SMILES through
the same pipeline, neither in the GBM library.

## Literature anchors

| use | value | source |
|---|---|---|
| doxorubicin fit anchor | 40–50 µm half-distance → 64.92 µm e-fold | PMID 16361566 |
| out-of-sample band | 25–75 µm half-distance → 36.1–108.2 µm e-fold | PMID 23680920 |
| vascular area fraction | 2.1 ± 0.3 % CD34 | PMID 36611806 |
| vessel caliber for r = 7.5 µm | — | PMID 19856441, PMID 24390982 |

Note the area-fraction anchor has since been shown **not commensurable** with
what the model computes — CD34 stains endothelium, so the lumen contributes no
stained pixels, while the model paints the whole profile. See
`results/vascular_segmentation_audit.md`. The commensurable replacement
(Marien 2016, PMID 27583442, areal fraction 4.70 % ± 2.38 %, counting lumen)
was already in the repo.

## Choices worth auditing

* `L` fitted on a **common** Euclidean-distance-from-vessel axis rebuilt at each
  dx, so it means the same thing in both arms. The legacy p95-on-the-KDE axis is
  also reported for the KDE arm to reproduce published numbers.
* Kill reported as **seed-free expected kill** (mean per-cell probability) to
  remove Bernoulli noise from the convergence question; realised 2-seed kill
  agrees to < 0.003.
* Empirical alpha solved by 2-probe power law plus a verification solve rather
  than bisection, saving ~20 solves; verification residual < 1 %.
* **The KDE arm at dx = 6 is doxorubicin only**, on cost. Its non-convergence
  therefore rests on one drug at the finest step, though the causal source
  statistic (KDE tissue-mean, −6.4 % from dx=8 to 6) is drug-independent.

## Cross-check

At dx=12 / dt=5 both arms reproduce `vessel_network_impact.md`'s dt=10 arms A
and B to four significant figures (discrete coverage 0.3836 / 0.9936 / 0.0000,
AUC 1.083e4 / 1.691e6 / 7.063e5, L 39.6 / 152.5 / 273.3; KDE 0.9692, 3.712e4,
95.9) — evidence the `dt` change is not doing the work.

Grid ladder compute ~5,700 s, plus three other drivers.

---

## ⚠ A CFL discrepancy this record exposes, and it bears on a published ladder

This record states measured maximum interstitial speeds for the **discrete**
field of **0.730 / 0.779 / 0.843 / 0.868 µm/s** at dx = 15/12/8/6, giving
`dt_max` of 20.5 / 15.4 / **9.5** / 6.9 s — i.e. **`dt = 10 s` violates the
advection CFL condition for the discrete field at dx ≤ 8**, which is why it
used dt = 5 s.

`results/grid_convergence_gbm.md` ran its whole ladder at **dt = 10 s** and
measured 0.7513 µm/s at dx=8 (CFL 0.939) and 0.6684 at dx=7 (CFL 0.955) — under
1, so the solver accepted them, but with almost no margin. The two speed
measurements **straddle the stability limit at dx = 8**: 0.7513 gives CFL 0.939,
0.843 gives 1.054.

Both cannot be right for the same configuration, and the difference is not
academic — the GBM ladder's conclusions rest on its fine rungs, and those are
the rungs at issue. Possible causes, none yet tested: a different tissue mask,
a different vascular amplitude, or drift in the vessel field between the two
runs.

**What should be done:** re-measure max interstitial speed on the shipped
configuration at dx = 8 and dx = 7, and if the higher figure is right, re-run
those rungs at a CFL-satisfying `dt` and state whether the GBM ladder's
non-monotonicity survives. Until then the two finest rungs of that ladder should
be read as provisional.

This also revises the correction in `dea609e`, which established the true grid
floor at dt = 10 s as dx = 7 (GBM) / dx = 8 (breast). That correction replaced a
**false** claim (dx = 10 violates CFL — it runs at 0.663) and its direction was
right, but the floor it quotes was measured on one field and may be optimistic
for the discrete default.

## Related open item

The promotion changed the **unit** of every absolute exposure by 3.36×, which is
a normalisation change and not a physics result (geometry-only effect ≤ 0.0166).
That was flagged at the time. What was not known then is that `supply = 0.05`
was never re-derived for the new referent —
see `results/oxygen_supply_calibration.md` and
`results/oxygen_level_shift.md`.
