"""The derived quantities the written record quotes, persisted so they trace.

Two things in the report are not fields of any stage's output:

1.  The coarse-vs-fine comparison restricted to pairs where ``dx`` ACTUALLY
    CHANGED THE GRID. ``coarse_baseline.py``'s own aggregate includes the three
    Visium sections, where the reconstruction pins the grid to the spot lattice
    so ``dx = 30`` and ``dx = 12`` are the identical solve -- 18 of 30 pairs
    comparing a run with itself, which is why that aggregate reads
    ``median_rel_err = 0.000 %``. The honest comparison is the 12 pairs on the
    two Xenium sections, and it is computed here rather than in the report's
    prose.
2.  The across-donor spread of each target on the ``spatial_spread``-VALID
    subset in ABSOLUTE units. ``analyse.py`` stores the fold spread and the
    per-donor medians for that subset but not the absolute range, and the
    absolute range is the number an error has to be read against.

Writes ``experiments/surrogate_feasibility/derived.json``.
"""
from __future__ import annotations

import json
import math
import os
import statistics as st
from typing import Dict, List

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_HERE, "derived.json")

TARGETS = ("E_p90_over_p10", "penetration_length_um", "penetration_index",
           "median_auc_over_peak")


def _f(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _pct(sorted_vals: List[float], q: float) -> float:
    if not sorted_vals:
        return float("nan")
    i = min(int(round(q * (len(sorted_vals) - 1))), len(sorted_vals) - 1)
    return sorted_vals[i]


def coarse_vs_fine() -> Dict:
    with open(os.path.join(_HERE, "coarse_baseline.json")) as fh:
        cb = json.load(fh)
    real, vac = [], []
    for p in cb["pairs"]:
        c, f = p.get("coarse") or {}, p.get("fine") or {}
        if c.get("status") != "ok" or f.get("status") != "ok":
            continue
        (real if c.get("grid") != f.get("grid") else vac).append(p)

    out: Dict = {
        "note": ("only pairs where dx changed the grid are compared; the rest "
                 "are the identical solve twice"),
        "n_pairs_grid_differs": len(real),
        "n_pairs_grid_identical": len(vac),
        "sections_grid_differs": sorted({p["section"] for p in real}),
        "sections_grid_identical": sorted({p["section"] for p in vac}),
        "grids": {p["section"]: {"coarse": p["coarse"]["grid"],
                                 "fine": p["fine"]["grid"]}
                  for p in real + vac},
        "targets": {},
    }
    for t in TARGETS:
        rel, abs_ = [], []
        for p in real:
            cv, fv = _f(p["coarse"].get(t)), _f(p["fine"].get(t))
            if cv is None or fv is None or fv == 0:
                continue
            rel.append(abs(cv - fv) / abs(fv))
            abs_.append(abs(cv - fv))
        rel.sort(); abs_.sort()
        if not rel:
            continue
        out["targets"][t] = {
            "n": len(rel),
            "median_rel_err": st.median(rel),
            "p90_rel_err": _pct(rel, 0.90),
            "max_rel_err": rel[-1],
            "median_rel_err_pct": st.median(rel) * 100.0,
            "p90_rel_err_pct": _pct(rel, 0.90) * 100.0,
            "max_rel_err_pct": rel[-1] * 100.0,
            "median_abs_err": st.median(abs_),
            "max_abs_err": abs_[-1],
        }
    for key in ("coarse", "fine"):
        w = sorted(p[key]["wall_clock_s"] for p in real)
        la = [p[key]["loadavg_1m"] for p in real
              if p[key].get("loadavg_1m") is not None]
        out[f"wall_{key}"] = {
            "n": len(w), "median_s": st.median(w), "min_s": w[0],
            "max_s": w[-1],
            "loadavg_1m_min": min(la) if la else None,
            "loadavg_1m_max": max(la) if la else None,
        }
    out["cost_ratio_fine_over_coarse_medians"] = (
        out["wall_fine"]["median_s"] / out["wall_coarse"]["median_s"])
    return out


def valid_only_spread() -> Dict:
    with open(os.path.join(_HERE, "sweep_data.json")) as fh:
        sw = json.load(fh)
    rows = [r for r in sw["rows"] if r.get("status") == "ok"
            and r.get("spatial_spread_valid")
            and r["total_time_s"] == 21600.0 and r["half_life_s"] == 14400.0
            and r["plasma_peak"] == 50.0 and r["drug"] == "temozolomide"]

    def target(r, t):
        if t == "E_p90_over_p10":
            a, b = _f(r.get("E_p90_uM")), _f(r.get("E_p10_uM"))
            return (a / b) if (a and b) else None
        if t == "median_auc_over_peak":
            a, b = _f(r.get("median_tumor_auc")), _f(r.get("plasma_peak_uM"))
            return (a / b) if (a is not None and b) else None
        return _f(r.get(t))

    out: Dict = {"config": ("temozolomide, 50 uM peak, 14400 s half-life, "
                            "1 dose, 21600 s window, dx = 30 um"),
                 "n_rows": len(rows),
                 "sections": sorted({r["section"] for r in rows}),
                 "targets": {}}
    for t in TARGETS:
        per = {r["section"]: target(r, t) for r in rows}
        per = {k: v for k, v in per.items() if v is not None}
        if len(per) < 2:
            continue
        vals = sorted(per.values())
        out["targets"][t] = {
            "n_donors": len(per),
            "per_section": per,
            "min": vals[0], "max": vals[-1],
            "abs_range": vals[-1] - vals[0],
            "fold_spread": (vals[-1] / vals[0]) if vals[0] > 0 else None,
        }
    return out


def main() -> int:
    out = {"coarse_vs_fine_grid_differs": coarse_vs_fine(),
           "valid_only_across_donor_spread": valid_only_spread()}
    with open(_OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
    print(f"\n[derived] -> {_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
