"""DOES THE P90/P10 DELIVERY SPREAD SEPARATE COMPOUNDS ABOVE ITS OWN NOISE FLOOR?

WHY THIS EXISTS
---------------
``results/kyle2004_penetration.md`` §8 measured a noise floor for the map's
*penetration-length* readout and found it FAILS on reconstructed tissue: the
empirical half-distances of doxorubicin (39.4-63.1 um) and daunorubicin
(32.7-53.9 um) overlap once the grid, the tissue-mask repair, the distance-bin
width and the distance metric are allowed to vary. The compound difference sat
inside the method's own spread.

The product has since been repositioned onto a DIFFERENT delivery readout: the
**P90/P10 delivery spread** of per-tumour-cell exposure
(``delivery_metrics.exposure_percentiles``, key ``E_p90_over_p10``), headlined on
the claim that it discriminates compounds — temozolomide ~1.99x against
doxorubicin ~3.80x on a matched section, grid, window and dose.

**That separation has never been run through a nuisance sweep.** This script does
exactly that, and nothing else. It is the same question asked of a different
metric.

WHICH NUISANCE AXES ACTUALLY ENTER THIS METRIC — read before comparing to §8
----------------------------------------------------------------------------
The §8 sweep varied four things. Only some of them can touch P90/P10, and saying
which is part of the result rather than a caveat on it:

* **distance-bin width — DOES NOT ENTER.** P90/P10 is a percentile ratio over
  tumour cells. There is no distance binning anywhere in it.
* **distance metric (direct-endothelial-cell vs percentile-mask) — DOES NOT
  ENTER.** No distance map is used.
* **exposure vs AUC — ALGEBRAICALLY DEGENERATE.** ``exposure_uM = auc_uM_s / T``
  with ``T`` a single scalar per run, and percentiles are homogeneous of degree
  one, so ``E_p90/E_p10`` and ``auc_p90/auc_p10`` are the SAME NUMBER. Asserted
  in ``--stage degenerate`` rather than argued.
* **grid resolution and tissue-mask repair — DO ENTER**, through the field solve
  and through which grid cell each measured tumour cell samples.
* **exposure window — DOES ENTER**, through the shape of the finite-window AUC
  field.

So P90/P10 is structurally exposed to FEWER reduction choices than the
half-distance readout was. That is a real point in its favour and it is why this
sweep is narrower than §8's — it is not a weaker sweep, it is the complete set of
axes the metric can see.

``delivery_metrics.tumor_cell_exposure`` states the design intent outright: *"The
tumor cell population is a FIXED set of measured Xenium cells, identical at every
dx -- which is what makes these statistics grid-comparable."* This script tests
that claim rather than accepting it.

WHAT IS CHANGED
---------------
Nothing. ``calibrate_alpha()`` is asserted live at ``22.223477711172563``. Drugs
are compiled through the flagship path the real panel uses
(``calibrated_drugparams`` -> ``apply_barrier_to_drugparams``). Nothing under
``tmesim/`` or ``app/`` is written to or modified; this script only drives them.
**Nothing is tuned to make the compounds separate.**

Run:
    ./.venv/bin/python experiments/delivery_spread_noise.py --stage degenerate
    ./.venv/bin/python experiments/delivery_spread_noise.py --stage all
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Sequence

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.core import Domain, Microenvironment, PDEConfig, SimConfig
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.reconstruct.real_fields import build_real_fields
from tmesim.sim.engine import simulate_response
from tmesim.validation import vessel_profile as VP
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.validation.delivery_metrics import (
    exposure_percentiles,
    tumor_cell_exposure,
)

SHIPPED_ALPHA = 22.223477711172563

#: The headline pair the repositioned panel rests on.
PAIR: Sequence[str] = ("temozolomide", "doxorubicin")

#: The published headline values this sweep is testing the robustness of.
CLAIMED = {"temozolomide": 1.99, "doxorubicin": 3.80}

#: Grid / mask configurations. The first is the shipped flagship configuration.
CONFIGS: List[Dict] = [
    {"key": "dx12_repair2", "dx_um": 12.0, "repair": 2, "shipped": True},
    {"key": "dx12_native", "dx_um": 12.0, "repair": 0, "shipped": False},
    {"key": "dx16_repair2", "dx_um": 16.0, "repair": 2, "shipped": False},
    {"key": "dx20_native", "dx_um": 20.0, "repair": 0, "shipped": False},
    {"key": "dx30_native", "dx_um": 30.0, "repair": 0, "shipped": False},
]

#: Exposure windows (s). 6 h is the shipped window.
WINDOWS_S: Sequence[float] = (3600.0 * 6, 3600.0 * 3, 3600.0 * 12)

_RESULTS = os.path.join(_ROOT, "results")
_DATA = os.path.join(_RESULTS, "delivery_spread_noise.json")


def _assert_alpha() -> float:
    a = calibrate_alpha()
    if a != SHIPPED_ALPHA:
        raise AssertionError(
            f"calibrate_alpha() = {a!r}, expected {SHIPPED_ALPHA!r}")
    return a


def compile_flagship(name: str, alpha: float):
    """The pipeline the shipped real-tissue panel uses."""
    dp = calibrated_drugparams(name, alpha)
    return apply_barrier_to_drugparams(dp, get_molecule(name))


def build(dx_um: float, repair: int):
    micro, _cn, _cr = VP.build_arm_microenvironment(dx_um=dx_um, repair=repair)
    return micro


def spread_for(micro, dp, window_s: float) -> Dict:
    pde = PDEConfig(dt=10.0, total_time=float(window_s), convection=True,
                    save_every=500)
    res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
    exp = tumor_cell_exposure(micro, res)
    pct = exposure_percentiles(exp["exposure_uM"], exp["auc_uM_s"])
    lo, hi = pct["auc_p10_uM_s"], pct["auc_p90_uM_s"]
    return {
        "E_p90_over_p10": float(pct["E_p90_over_p10"]),
        "auc_p90_over_p10": float(hi / lo) if lo > 0 else float("nan"),
        "E_p10_uM": float(pct["E_p10_uM"]),
        "E_p50_uM": float(pct["E_p50_uM"]),
        "E_p90_uM": float(pct["E_p90_uM"]),
        "n_tumor": int(exp["n_tumor"]),
        "T_s": float(exp["T_s"]),
    }


def stage_degenerate() -> None:
    """Prove the E-vs-AUC axis is degenerate rather than asserting it in prose."""
    alpha = _assert_alpha()
    micro = build(30.0, 0)                     # cheapest grid; the claim is exact
    dp = compile_flagship("temozolomide", alpha)
    r = spread_for(micro, dp, 3600.0 * 6)
    rel = abs(r["E_p90_over_p10"] - r["auc_p90_over_p10"]) / r["E_p90_over_p10"]
    print(f"E_p90/E_p10   = {r['E_p90_over_p10']:.10f}")
    print(f"auc_p90/p10   = {r['auc_p90_over_p10']:.10f}")
    print(f"relative diff = {rel:.3e}   -> "
          f"{'DEGENERATE, as claimed' if rel < 1e-12 else 'NOT degenerate'}")


def main() -> Dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["all", "degenerate"])
    args = ap.parse_args()
    if args.stage == "degenerate":
        stage_degenerate()
        return {}

    alpha = _assert_alpha()
    t0 = time.time()
    drugs = {n: compile_flagship(n, alpha) for n in PAIR}
    rows: List[Dict] = []

    micros = {}
    for cfg in CONFIGS:
        micros[cfg["key"]] = build(cfg["dx_um"], cfg["repair"])
        print(f"  built {cfg['key']}")
    for w in WINDOWS_S:
        for cfg in CONFIGS:
            micro = micros[cfg["key"]]
            for name in PAIR:
                t1 = time.time()
                r = spread_for(micro, drugs[name], w)
                rows.append({"drug": name, "config": cfg["key"],
                             "dx_um": cfg["dx_um"], "repair": cfg["repair"],
                             "window_s": w, "shipped_config": bool(
                                 cfg["shipped"] and w == 3600.0 * 6),
                             **r})
                print(f"  {cfg['key']:<14} {int(w/3600):>2}h  {name:<14} "
                      f"P90/P10 = {r['E_p90_over_p10']:7.3f}   "
                      f"n_tumor={r['n_tumor']}  ({time.time()-t1:.0f}s)")

    out: Dict = {"alpha": alpha, "rows": rows, "claimed": CLAIMED,
                 "elapsed_s": time.time() - t0}

    print("\n" + "=" * 78)
    print("NOISE FLOOR — P90/P10 delivery spread")
    print("=" * 78)
    spans = {}
    for name in PAIR:
        v = [r["E_p90_over_p10"] for r in rows if r["drug"] == name
             and np.isfinite(r["E_p90_over_p10"])]
        ship = [r["E_p90_over_p10"] for r in rows if r["drug"] == name
                and r["shipped_config"]]
        spans[name] = {"min": min(v), "max": max(v), "n": len(v),
                       "shipped": ship[0] if ship else None,
                       "claimed": CLAIMED[name]}
        print(f"  {name:<14} shipped {spans[name]['shipped']:.3f} "
              f"(claimed {CLAIMED[name]})   span over {len(v)} cells: "
              f"{min(v):.3f} - {max(v):.3f}")
    a, b = PAIR
    disjoint = (spans[a]["max"] < spans[b]["min"]
                or spans[b]["max"] < spans[a]["min"])
    out["spans"] = spans
    out["disjoint"] = bool(disjoint)
    gap = (spans[b]["min"] - spans[a]["max"] if spans[a]["max"] < spans[b]["min"]
           else spans[a]["min"] - spans[b]["max"])
    out["gap"] = float(gap)
    print(f"\n  DISJOINT: {disjoint}    "
          f"{'gap' if disjoint else 'overlap'} = {abs(gap):.3f}")
    print(f"  -> the P90/P10 separation "
          f"{'SURVIVES' if disjoint else 'DOES NOT SURVIVE'} the nuisance sweep")

    with open(_DATA, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {_DATA}\ntotal {time.time()-t0:.1f}s")
    return out


if __name__ == "__main__":
    main()
