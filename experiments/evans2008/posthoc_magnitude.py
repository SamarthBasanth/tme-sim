"""POST-HOC diagnostic. Not pre-registered, and it CANNOT move the verdict.

The pre-registered primary (``results/evans2008_perivascular_sign_2026-10-10.md`` 4.3)
returned branch R3. This script asks only *why*: how large is the per-vessel perivascular
pO2 change whose sign the primary statistic counts? A sign is only meaningful if the
underlying change is resolvable.

Everything printed here is explicitly labelled post-hoc in the report and is used to
explain the R3 result, never to reclassify it.

Run::

    ./.venv/bin/python experiments/evans2008/posthoc_magnitude.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perivascular_sign import (ALPHA_EXPECTED, BAND_CENTRES,  # noqa: E402
                               XENIUM_SECTIONS, assert_alpha, band_table,
                               clopper_pearson, ols_slope)

SECTIONS = ("xenium_319", "xenium_io", "gse324268")
DX = 8.0
OUT = os.path.join("experiments", "evans2008", "out", "posthoc_magnitude.json")


def main() -> None:
    from tmesim.core import OxygenConfig
    from tmesim.physics.oxygen import REFERENCE_PO2_MMHG, solve_oxygen
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.reconstruct.vessel_network import VESSEL_RADIUS_UM, segment_vessels

    assert_alpha("posthoc start")
    out = {"alpha": ALPHA_EXPECTED, "dx_um": DX, "sections": {}}
    pooled = {"slope": [], "drop": [], "cells": []}

    for sec in SECTIONS:
        micro = reconstruct_real(source="xenium", sample=XENIUM_SECTIONS[sec], dx_um=DX)
        O = solve_oxygen(micro, OxygenConfig())
        assert np.array_equal(O, np.asarray(micro.fields.oxygen))
        vessels = segment_vessels(micro.cells, micro.domain)
        bt = band_table(O, micro, vessels, VESSEL_RADIUS_UM, 3)

        mm = bt["band_means_mmhg"]
        slope = ols_slope(BAND_CENTRES, bt["band_means_o2"]) * float(REFERENCE_PO2_MMHG)
        drop = mm[:, 0] - mm[:, 4]            # >0 means pO2 falls away from the vessel
        cells = bt["cells_per_vessel"]
        pooled["slope"].append(slope)
        pooled["drop"].append(drop)
        pooled["cells"].append(cells)

        exp_m, anti_m = slope < 0, slope > 0
        rec = {
            "n_analysable": int(bt["n_analysable"]),
            "median_abs_slope_mmhg_per_um": float(np.median(np.abs(slope))),
            "median_abs_slope_expected_group": float(np.median(np.abs(slope[exp_m]))),
            "median_abs_slope_anti_group": float(np.median(np.abs(slope[anti_m]))),
            "abs_drop_mmhg_percentiles": {
                p: float(np.percentile(np.abs(drop), p)) for p in (50, 90, 99, 100)
            },
            "frac_abs_drop_below_0p1_mmhg": float((np.abs(drop) < 0.1).mean()),
            "frac_abs_drop_below_0p01_mmhg": float((np.abs(drop) < 0.01).mean()),
        }
        out["sections"][sec] = rec
        print(f"\n{sec}: n={rec['n_analysable']}")
        print(f"  median |slope| = {rec['median_abs_slope_mmhg_per_um']:.6f} mmHg/um "
              f"(expected group {rec['median_abs_slope_expected_group']:.6f}, "
              f"anti group {rec['median_abs_slope_anti_group']:.6f})")
        print(f"  |pO2(band1) - pO2(band5)| mmHg: p50 "
              f"{rec['abs_drop_mmhg_percentiles'][50]:.4f}  p90 "
              f"{rec['abs_drop_mmhg_percentiles'][90]:.4f}  p99 "
              f"{rec['abs_drop_mmhg_percentiles'][99]:.4f}  max "
              f"{rec['abs_drop_mmhg_percentiles'][100]:.4f}")
        print(f"  fraction with |drop| < 0.1 mmHg: "
              f"{rec['frac_abs_drop_below_0p1_mmhg']:.4f};  < 0.01 mmHg: "
              f"{rec['frac_abs_drop_below_0p01_mmhg']:.4f}")

    slope = np.concatenate(pooled["slope"])
    drop = np.concatenate(pooled["drop"])
    cells = np.concatenate(pooled["cells"])

    print("\n=== POOLED, post-hoc ===")
    rows = []
    for thr in (0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0):
        sel = np.abs(drop) >= thr
        if sel.sum() == 0:
            rows.append({"threshold_mmhg": thr, "n": 0})
            continue
        e = int((slope[sel] < 0).sum())
        a = int((slope[sel] > 0).sum())
        lo, hi = clopper_pearson(e, e + a)
        rows.append({"threshold_mmhg": thr, "n": int(sel.sum()),
                     "frac_of_all": float(sel.mean()),
                     "n_expected": e, "n_anti": a,
                     "p_expected": e / (e + a), "cp95": [lo, hi]})
        print(f"  |drop| >= {thr:5.2f} mmHg: n = {int(sel.sum()):6d} "
              f"({sel.mean():6.2%} of vessels)  P_exp = {e/(e+a):.4f} "
              f"[{lo:.4f}, {hi:.4f}]")
    out["pooled_sign_split_vs_resolvability"] = rows
    out["pooled"] = {
        "n": int(slope.size),
        "median_abs_slope_mmhg_per_um": float(np.median(np.abs(slope))),
        "median_abs_drop_mmhg": float(np.median(np.abs(drop))),
        "frac_abs_drop_below_0p1_mmhg": float((np.abs(drop) < 0.1).mean()),
        "singleton_fraction_of_analysable": float((cells == 1).mean()),
        "median_cells_per_analysable_vessel": float(np.median(cells)),
    }
    print(f"\n  pooled median |slope| = "
          f"{out['pooled']['median_abs_slope_mmhg_per_um']:.6f} mmHg/um")
    print(f"  pooled median |drop over 40 um| = "
          f"{out['pooled']['median_abs_drop_mmhg']:.4f} mmHg")
    print(f"  pooled fraction of analysable vessels that are single endothelial "
          f"cells: {out['pooled']['singleton_fraction_of_analysable']:.4f}")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    assert_alpha("posthoc end")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
