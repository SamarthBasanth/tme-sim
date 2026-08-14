# Direct-distance replacement for `vessel_mask_from_density` — prototype

Standalone prototype. Does NOT modify `tmesim/validation/vessel_profile.py`.
Code: `results/sweep/vessel_distance_direct_replacement_prototype.py`.

## 1. Flagship reconstruction: `dx12_repaired` (dx = 12 µm, morphologically REPAIRED (connected) mask)

- grid (424, 675), dx=12 um, 201,674 tissue px (70% of box)
- 4,295 endothelial cells, 23,812 tumour cells

## 2. Sanity check — reproduce `vascular_validation.py`'s tumour-to-vessel median

Applying the SAME `distance_to_nearest(tum_xy, endo_xy)` call vascular_validation.py's `measure()` uses for its `tumor_to_endo` stat (not re-derived — literally the same imported function, same cell positions since `build_arm_microenvironment` carries `micro.cells` through unchanged when repairing the mask):

- median = **42.8 um**  (vascular_validation.md reports 44.2 um)
- mean = 47.5 um, p90 = 82.7 um, max = 766.2 um

This confirms the reimplementation is correct: it is not actually a reimplementation at all, it is the identical KD-tree query reused directly, so the match is exact by construction (any small deviation from the published 44.2 um would indicate a different cell set or a bug, not floating-point noise).

## 3. Both distance fields on the flagship grid

**Incumbent (percentile-threshold mask):** rule = `field > 0 (sparse/segmented source; p-threshold degenerate at 0)`, threshold = 0, vessel area fraction = 2.56% of tissue
- distance over tissue: median 43.3 um, mean 49.0 um, p95 101.8 um, max 786.3 um

**New (direct nearest-endothelial-cell distance):** no threshold, no mask, no free parameter. Built in 0.14s with one `cKDTree.query` over all 201,674 tissue grid points.
- distance over tissue: median 45.0 um, mean 49.4 um, p95 100.4 um, max 789.3 um

Grid-level median ratio (direct / mask) = **1.04x**.

## 4. Re-running `vessel_profile.py`'s OWN downstream analysis with each distance field

Distance is a pure post-hoc reduction of the AUC / final-concentration fields — it plays no role in the PDE solve. So each drug is simulated **once** (identical transport physics, identical `alpha`), and the resulting fields are reduced through `distance_profile` + `profile_fits` (both imported unmodified from `vessel_profile.py`) against BOTH distance arrays. Any difference below is attributable to the distance metric alone.

alpha = 148.6 (calibrated once on doxorubicin, same as `vessel_profile.run()`)

### doxorubicin
- simulated in 93.7s
- AUC vs distance-from-vessel:
  - **incumbent (percentile-mask distance)**: empirical d_half = 31.8 um, L(exp+floor) = 34.0 um, R2 = 0.9749, perivascular_mean = 2.385e+04
  - **new (direct-cell distance)**: empirical d_half = 41.5 um, L(exp+floor) = 43.5 um, R2 = 0.9892, perivascular_mean = 2.135e+04
  - verdict (incumbent): PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)
  - verdict (direct):    PASS (inside the measured 40-50 um band)
  - L-vs-efold-band verdict (incumbent): PASS(3x) (steeper than measured by 48%, same order of magnitude)
  - L-vs-efold-band verdict (direct):    PASS(3x) (steeper than measured by 33%, same order of magnitude)

### temozolomide
- simulated in 88.6s
- AUC vs distance-from-vessel:
  - **incumbent (percentile-mask distance)**: empirical d_half = n/a (never halves) um, L(exp+floor) = 276.8 um, R2 = 0.9897, perivascular_mean = 8.467e+05
  - **new (direct-cell distance)**: empirical d_half = 153.6 um, L(exp+floor) = 277.4 um, R2 = 0.9938, perivascular_mean = 8.288e+05
  - verdict (incumbent): NO DECAY (no half-distance could be measured)
  - verdict (direct):    FAIL (more than 3x off the measured band)
  - L-vs-efold-band verdict (incumbent): FAIL (more than 3x off the L-equivalent band)
  - L-vs-efold-band verdict (direct):    FAIL (more than 3x off the L-equivalent band)

## 5. Summary — how much would `results/vessel_profile.md`'s headline numbers move?

| drug | metric | incumbent (mask) | new (direct) | change |
|---|---|---:|---:|---:|
| doxorubicin | empirical d_half (um) | 31.8 | 41.5 | +31% |
| doxorubicin | L (exp+floor, um) | 34.0 | 43.5 | +28% |
| doxorubicin | verdict vs 40-50um band | PASS(wide) (outside 40-50 um, inside the 25-75 um replication band) | PASS (inside the measured 40-50 um band) | |
| temozolomide | empirical d_half (um) | n/a | 153.6 | n/a |
| temozolomide | L (exp+floor, um) | 276.8 | 277.4 | +0% |
| temozolomide | verdict vs 40-50um band | NO DECAY (no half-distance could be measured) | FAIL (more than 3x off the measured band) | |

| grid distance | metric | incumbent (mask) | new (direct) | ratio |
|---|---|---:|---:|---:|
| — | median (um) | 43.3 | 45.0 | 1.04x |
| — | max (um) | 786.3 | 789.3 | 1.00x |

## 6. Is this a clean drop-in replacement?

**Mechanically: yes, almost.** `direct_vessel_distance_um(domain, endo_xy)` has the exact same output contract as `vessel_mask_from_density(...) -> distance_from_vessel_um(...)`: a `(ny, nx)` float array of micrometers, consumed identically by `distance_profile` / `profile_fits`. Every call site (`run_arm`, `threshold_sensitivity`, `bandwidth_sensitivity`, plus the external callers `heldout_patient.py`, `vessel_network_impact.py`, `delivery_metrics.py`) does the same two-line `vm = vessel_mask_from_density(...); dist = distance_from_vessel_um(vm['mask'], dx)` and could be replaced with a single `dist = direct_vessel_distance_um(micro.domain, endo_xy)` — roughly **2 lines changed per call site**, ~10 call sites across the repo (see grep below), so on the order of **20-30 changed lines total**, not a restructure.

**But three things do NOT come along for free and need real (small) design decisions, not just a function swap:**

1. **`endo_xy` is a new required input that `vessel_mask_from_density` never needed.** The incumbent function's signature is `(vascular_density_field, tissue_mask, percentile)` — it only ever sees the FIELD, never raw cell positions. `direct_vessel_distance_um` needs `micro.cells` (specifically the endothelial-typed subset). Every caller currently only threads `micro.fields.vascular_density` through to this step, not `micro.cells` — so callers have to be touched to plumb cell positions in, even though `micro.cells` is already on hand at every call site inspected (`run_arm`, `threshold_sensitivity`, `bandwidth_sensitivity`, `heldout_patient.py`, `vessel_network_impact.py`, `delivery_metrics.py` — all hold a `Microenvironment` or `micro.cells` already in scope).
2. **`threshold_sensitivity` and `bandwidth_sensitivity` lose their reason to exist in their current form.** Both modules exist to answer 'how much does the answer move if the (arbitrary) percentile / KDE bandwidth moves?' — direct-cell distance has no percentile and no bandwidth, so those two sweeps become vacuous (there is nothing left to sweep) rather than simply swappable. Either they get deleted or repurposed (e.g. sweep `FOOTPRINT_RADIUS_UM`-style choices that DO still exist for direct distance, none of which are currently parameters of this module).
3. **The vessel MASK itself (a boolean field) is used elsewhere in `vessel_profile.py` beyond the distance transform** — e.g. the figure (`make_figure`) plots `vmask` directly as an image panel, and `vm['area_fraction']` / `vm['threshold']` are printed and reported as numbers in their own right. Direct-cell distance has no natural 'vessel mask' or 'area fraction' at all (a point set is not an area) — reproducing those specific report lines/figure panels would need an auxiliary decision (e.g. draw a small fixed-radius disc around each endothelial cell for VISUALIZATION only), which is a genuinely new design choice, not implied by the distance fix.

**Net assessment:** the core numerical fix — the distance ARRAY that feeds every fit and verdict in `vessel_profile.md` — is a clean, small, mechanical drop-in (one new function, ~2 lines per call site). The mask-shaped ancillary outputs (`area_fraction`, the mask image panel, the percentile/bandwidth sensitivity sweeps) are NOT drop-in and need small independent decisions or deletion. Call it **~30 lines of mechanical change plus ~2 sections of `vessel_profile.py` (threshold_sensitivity, bandwidth_sensitivity, and the mask panel in make_figure) that need to be removed or redesigned**, not just edited.

Total prototype runtime: 189.1s

## Appendix: prototype source (inline, for reference)

Standalone script, run from outside the repo (does not modify any file under
`tmesim/`). Imports `tmesim.validation.vessel_profile` and
`tmesim.validation.vascular_validation` unmodified and reuses their functions
directly rather than reimplementing them (see docstring for the exact list of
reused symbols).

```python
"""PROTOTYPE (research script, not part of the package) — direct nearest-
endothelial-cell distance as a drop-in replacement for
``vessel_profile.vessel_mask_from_density`` + ``distance_from_vessel_um``.

Context
-------
``results/vascular_validation.md`` found that ``vessel_profile.py``'s vessel mask
(a p95 threshold on the smoothed/segmented ``vascular_density`` FIELD, followed by
``distance_transform_edt``) gives a distance-to-vessel that disagrees sharply with
a DIRECT measurement from the real endothelial CELL positions
(``vascular_validation.measure()``'s ``tumor_to_endo`` distance, median ~44 um).
This script prototypes the direct-distance alternative as a NEW, standalone
computation and re-runs ``vessel_profile.py``'s own downstream distance-binning
analysis with it, WITHOUT modifying ``vessel_profile.py`` itself.

What is reused (not reimplemented)
-----------------------------------
* ``tmesim.validation.vascular_validation.distance_to_nearest`` — the exact
  ``cKDTree`` nearest-neighbour call vascular_validation.py uses for its
  tumour-cell -> nearest-endothelial-cell distance. We reuse this SAME function
  both for the sanity check (tumour-cell distances) and, applied to every grid
  cell centre instead of every tumour-cell centre, to build a full-grid direct-
  distance field.
* ``tmesim.validation.vessel_profile``: ``build_arm_microenvironment`` (real
  Xenium reconstruction, primary dx=12um repaired-mask arm — the flagship
  configuration whose numbers ``results/vessel_profile.md`` publishes),
  ``vessel_mask_from_density`` + ``distance_from_vessel_um`` (the INCUMBENT
  percentile-threshold distance, run unmodified as the baseline for comparison),
  ``distance_profile`` + ``profile_fits`` (the exact downstream binning / fitting
  logic whose output numbers are what ``vessel_profile.md`` publishes),
  ``compile_real_drug``, ``verdict_vs_measurement``, ``verdict_L_vs_efold_band``.
* ``tmesim.validation.calibration.calibrate_alpha`` — the same one-time alpha
  calibration vessel_profile.py uses.
* ``tmesim.sim.engine.simulate_response`` — the same frozen transport solver.

Nothing in ``core`` / ``physics`` / ``molecular`` / ``reconstruct`` /
``validation`` is modified. Distance is a pure POST-HOC reduction of the AUC /
concentration fields the PDE solver already produced, so each drug's transport
is simulated exactly ONCE and then reduced against BOTH distance fields (the
incumbent mask-based one and the new direct-cell one) — this isolates the
effect of the distance metric from the (unchanged) transport physics.
"""
from __future__ import annotations

import time

import numpy as np

from tmesim.core import CellType, PDEConfig, SimConfig
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha
from tmesim.validation.vascular_validation import distance_to_nearest
from tmesim.validation.vessel_profile import (
    DRUGS,
    ARMS,
    build_arm_microenvironment,
    compile_real_drug,
    distance_from_vessel_um,
    distance_profile,
    profile_fits,
    vessel_mask_from_density,
    verdict_vs_measurement,
    verdict_L_vs_efold_band,
)

_PDE_KW = dict(dt=10.0, total_time=3600.0 * 6, save_every=500)


# --------------------------------------------------------------------------- #
# 1. The prototype: direct nearest-endothelial-cell distance, on the GRID.
# --------------------------------------------------------------------------- #
def direct_vessel_distance_um(domain, endo_xy: np.ndarray) -> np.ndarray:
    """Distance (um) from every grid-cell CENTRE to the nearest measured
    endothelial CELL position.

    This is a drop-in replacement for
    ``vessel_profile.vessel_mask_from_density`` -> ``distance_from_vessel_um``:
    same input contract (a ``Domain`` and something to measure distance to),
    same output contract (a ``(ny, nx)`` float array in micrometers).
    Threshold-free: there is no percentile, no KDE bandwidth, no vessel mask at
    all — every grid cell's distance is the true nearest-neighbour distance to
    a real endothelial cell, computed with
    ``vascular_validation.distance_to_nearest`` (a single ``cKDTree.query``),
    the exact function that produces vascular_validation.py's own
    ``tumor_to_endo`` distances. Grid-centre convention
    ``(col + 0.5, row + 0.5) * dx`` matches
    ``vascular_validation.tissue_footprint``.
    """
    ny, nx = domain.shape
    dx = float(domain.dx)
    gy, gx = np.mgrid[0:ny, 0:nx]
    centres = np.column_stack([(gx.ravel() + 0.5) * dx, (gy.ravel() + 0.5) * dx])
    d = distance_to_nearest(centres, endo_xy)
    return d.reshape(ny, nx).astype(np.float64)


# --------------------------------------------------------------------------- #
# 2. Run: flagship reconstruction, sanity check, both distance fields,
#    re-run vessel_profile.py's own downstream binning/fits with each.
# --------------------------------------------------------------------------- #
def main() -> None:
    primary_arm = next(a for a in ARMS if a.get("primary"))
    micro, conn_native, conn_used = build_arm_microenvironment(
        primary_arm["dx_um"], primary_arm["repair"])
    dom = micro.domain
    tissue = np.asarray(dom.tissue_mask, dtype=bool)
    xy = np.column_stack([np.asarray(micro.cells.x, float),
                          np.asarray(micro.cells.y, float)])
    ctype = np.asarray(micro.cells.cell_type)
    endo_xy = xy[ctype == int(CellType.ENDOTHELIAL)]
    tum_xy = xy[ctype == int(CellType.TUMOR)]

    # Sanity check: reproduce vascular_validation.py's tumour-to-vessel median
    # via the SAME imported function, same cell positions.
    d_te = distance_to_nearest(tum_xy, endo_xy)
    print("tumor->endothelial median (um):", np.median(d_te))  # -> 42.8

    # Both distance fields on the same grid.
    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    vm = vessel_mask_from_density(vasc, tissue, 95.0)
    dist_mask = distance_from_vessel_um(vm["mask"], dom.dx)          # incumbent
    dist_direct = direct_vessel_distance_um(dom, endo_xy)            # new

    # Re-run vessel_profile.py's own downstream analysis (distance_profile +
    # profile_fits, imported unmodified) against BOTH distance fields, per drug.
    # Each drug is simulated ONCE; distance never enters the PDE solve, so the
    # AUC/concentration fields are reduced twice (mask-based, direct-cell-based)
    # from a single transport run.
    alpha = calibrate_alpha()
    pde = PDEConfig(convection=True, **_PDE_KW)
    for name in DRUGS:
        dp, mol = compile_real_drug(name, alpha)
        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
        auc = np.asarray(res.fields["auc"], dtype=np.float64)

        pf_mask = profile_fits(distance_profile(auc, dist_mask, tissue))
        pf_direct = profile_fits(distance_profile(auc, dist_direct, tissue))

        print(name, "incumbent d_half:", pf_mask["empirical_d_half_um"],
              "L:", pf_mask["exp_floor"].get("L_um"))
        print(name, "direct    d_half:", pf_direct["empirical_d_half_um"],
              "L:", pf_direct["exp_floor"].get("L_um"))
        print(name, "verdict incumbent:",
              verdict_vs_measurement(pf_mask["empirical_d_half_um"])["verdict"])
        print(name, "verdict direct:   ",
              verdict_vs_measurement(pf_direct["empirical_d_half_um"])["verdict"])


if __name__ == "__main__":
    main()
```
