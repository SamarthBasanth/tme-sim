"""Facts about the sweep data that decide what a surrogate can even be asked.

Three things, all measured from ``sweep_data.json``:

1.  WHICH AXES ARE NULL. The sweep varies ``half_life_s`` at a fixed window and
    ``plasma_peak`` on a few sections. If a target does not move, the axis
    carries no information and rows along it are duplicates -- which matters
    because duplicated rows shrink a cross-validated error without the model
    having learned anything.
2.  WHETHER THE REACH LADDER IS A USABLE TARGET AT ALL. A target with no
    variance is predicted perfectly by a constant, so an impressive-looking
    error on it means nothing.
3.  THE SPREAD OF EACH TARGET ACROSS DONORS, which is the denominator any
    reported error has to be read against.

Prints and writes ``experiments/surrogate_feasibility/analyse_results.json``.
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from typing import Dict, List

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_IN = os.path.join(_HERE, "sweep_data.json")
_OUT = os.path.join(_HERE, "analyse_results.json")


def _f(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def targets_of(r: Dict) -> Dict:
    e10, e90 = _f(r.get("E_p10_uM")), _f(r.get("E_p90_uM"))
    med, pk = _f(r.get("median_tumor_auc")), _f(r.get("plasma_peak_uM"))
    return {
        "E_p90_over_p10": (e90 / e10) if (e10 and e90 and e10 > 0) else None,
        "penetration_length_um": _f(r.get("penetration_length_um")),
        "median_auc_over_peak": (med / pk) if (med is not None and pk) else None,
        "penetration_index": _f(r.get("penetration_index")),
        "reach_headline_F": _f(r.get("reach_headline_F")),
        "E_p50_uM": _f(r.get("E_p50_uM")),
    }


def main() -> int:
    with open(_IN) as fh:
        blob = json.load(fh)
    rows = [r for r in blob["rows"] if r.get("status") == "ok"]
    out: Dict = {"n_rows_ok": len(rows), "n_rows_total": len(blob["rows"]),
                 "alpha": blob.get("alpha"), "dx_um": blob.get("dx_um")}

    # ---------------------------------------------------------------- 1. axes
    # Group rows that are identical EXCEPT for one axis, and report the largest
    # relative change the axis produces in each target.
    axis_report: Dict[str, Dict] = {}
    for axis, held in (("half_life_s", ("section", "drug", "n_doses",
                                        "total_time_s", "plasma_peak")),
                       ("plasma_peak", ("section", "drug", "n_doses",
                                        "total_time_s", "half_life_s")),
                       ("total_time_s", ("section", "drug", "n_doses",
                                         "half_life_s", "plasma_peak"))):
        groups: Dict[tuple, List[Dict]] = defaultdict(list)
        for r in rows:
            groups[tuple(r.get(k) for k in held)].append(r)
        per_target: Dict[str, List[float]] = defaultdict(list)
        n_cmp = 0
        for g in groups.values():
            if len({r.get(axis) for r in g}) < 2:
                continue
            n_cmp += 1
            tg = [targets_of(r) for r in g]
            for k in tg[0]:
                vals = [t[k] for t in tg if t[k] is not None and t[k] > 0]
                if len(vals) >= 2:
                    per_target[k].append(max(vals) / min(vals))
        axis_report[axis] = {
            "n_comparison_groups": n_cmp,
            "max_fold_change_per_target": {
                k: {"max": float(np.max(v)), "median": float(np.median(v)),
                    "n": len(v)}
                for k, v in per_target.items()},
        }
    out["axis_sensitivity"] = axis_report

    # -------------------------------------------------------- 2. reach ladder
    Fs = [t["reach_headline_F"] for t in map(targets_of, rows)]
    Fok = [v for v in Fs if v is not None]
    nulls = sum(1 for r in rows if r.get("reach_null_reason"))
    ladders = [tuple(r.get("reach_ladder_F") or []) for r in rows]
    out["reach_ladder"] = {
        "n_rows": len(rows),
        "n_null_reason": nulls,
        "n_headline_F_finite": len(Fok),
        "n_headline_F_exactly_zero": sum(1 for v in Fok if v == 0.0),
        "frac_headline_F_exactly_zero": (
            sum(1 for v in Fok if v == 0.0) / len(Fok)) if Fok else None,
        "headline_F_min": min(Fok) if Fok else None,
        "headline_F_max": max(Fok) if Fok else None,
        "headline_F_median": float(np.median(Fok)) if Fok else None,
        "n_distinct_ladders": len(set(ladders)),
        # Does the headline move with the COMPOUND at fixed tissue+config?
        "compound_fold_change": None,
    }
    per_cfg: Dict[tuple, set] = defaultdict(set)
    for r in rows:
        t = targets_of(r)
        if t["reach_headline_F"] is not None:
            per_cfg[(r["section"], r["total_time_s"], r["half_life_s"],
                     r["plasma_peak"])].add(round(t["reach_headline_F"], 12))
        # noqa
    moved = sum(1 for v in per_cfg.values() if len(v) > 1)
    out["reach_ladder"]["n_configs_where_F_moves_with_compound"] = moved
    out["reach_ladder"]["n_configs_checked"] = len(per_cfg)

    # ------------------------------------------------- 3. across-donor spread
    spreads = {}
    for k in ("E_p90_over_p10", "penetration_length_um",
              "median_auc_over_peak", "penetration_index"):
        shipped = [r for r in rows
                   if r["total_time_s"] == 21600.0 and r["half_life_s"] == 14400.0
                   and r["plasma_peak"] == 50.0 and r["drug"] == "temozolomide"]
        by_donor: Dict[str, List[float]] = defaultdict(list)
        for r in shipped:
            v = targets_of(r)[k]
            if v is not None:
                by_donor[r["donor"]].append(v)
        med = {d: float(np.median(v)) for d, v in by_donor.items() if v}
        valid_med = {d: float(np.median([
            targets_of(r)[k] for r in shipped
            if r["donor"] == d and targets_of(r)[k] is not None]))
            for d in med
            if any(r["spatial_spread_valid"] for r in shipped
                   if r["donor"] == d)}
        vals = list(med.values())
        spreads[k] = {
            "config": ("temozolomide, 50 uM peak, 14400 s half-life, "
                       "1 dose, 21600 s window, dx=30 um"),
            "n_donors": len(med),
            "min": min(vals) if vals else None,
            "max": max(vals) if vals else None,
            "fold_spread": (max(vals) / min(vals))
            if vals and min(vals) > 0 else None,
            "abs_range": (max(vals) - min(vals)) if vals else None,
            "per_donor_median": med,
            "valid_only_per_donor_median": valid_med,
            "valid_only_fold_spread": (
                max(valid_med.values()) / min(valid_med.values()))
            if valid_med and min(valid_med.values()) > 0 else None,
        }
    out["across_donor_spread"] = spreads

    # ---------------------------------------------------------------- 4. cost
    ok = [r for r in rows]
    walls = [r["wall_clock_s"] for r in ok if _f(r.get("wall_clock_s"))]
    warm = [r["wall_clock_s"] for r in ok if r.get("recon_cache_hit")]
    out["solver_cost_dx30"] = {
        "n": len(walls),
        "wall_median_s": float(np.median(walls)) if walls else None,
        "wall_p90_s": float(np.percentile(walls, 90)) if walls else None,
        "wall_max_s": float(np.max(walls)) if walls else None,
        "warm_median_s": float(np.median(warm)) if warm else None,
        "warm_max_s": float(np.max(warm)) if warm else None,
        "recon_build_s_median": float(np.median(
            [r["recon_build_s"] for r in ok if _f(r.get("recon_build_s"))])),
    }

    with open(_OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str)[:7000])
    print(f"\n[analyse] -> {_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
