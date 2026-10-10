#!/usr/bin/env python3
"""Re-derive, from the surrogate study's OWN rows, which readouts are invariant
under a pure dose / half-life change and which scale -- and by what exponent.

Nothing here trusts ``results/surrogate_feasibility_2026-10-09.md`` section 3.
It groups rows that are identical on every axis EXCEPT the one under test and
reports, per readout, whether the two values are BIT-IDENTICAL, exactly a
stated factor, or neither.
"""
from __future__ import annotations

import json
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP = os.path.normpath(
    os.path.join(HERE, "..", "surrogate_feasibility", "sweep_data.json"))

#: Everything that identifies a run other than the axis under test.
IDENTITY = ("section", "drug", "dx_um", "total_time_s", "n_doses",
            "plasma_peak", "half_life_s")


def numeric_keys(rows):
    keys = set()
    for r in rows:
        for k, v in r.items():
            if isinstance(v, bool) or v is None:
                continue
            if isinstance(v, (int, float)):
                keys.add(k)
    return sorted(keys)


def group(rows, axis):
    """{identity-without-axis: {axis value: row}}"""
    out = defaultdict(dict)
    rest = tuple(k for k in IDENTITY if k != axis)
    for r in rows:
        if r.get("status") != "ok":
            continue
        ident = tuple(r.get(k) for k in rest)
        out[ident][r.get(axis)] = r
    return {k: v for k, v in out.items() if len(v) >= 2}


def classify(axis, rows, keys):
    groups = group(rows, axis)
    verdict = {}
    for k in keys:
        bit = 0
        folds = []
        ratio_norm = []          # fold / axis-ratio, i.e. exponent-1 residual
        n = 0
        for by_val in groups.values():
            vals = sorted(by_val)
            for i in range(len(vals) - 1):
                a, b = by_val[vals[i]], by_val[vals[i + 1]]
                va, vb = a.get(k), b.get(k)
                if va is None or vb is None:
                    continue
                if isinstance(va, bool) or isinstance(vb, bool):
                    continue
                if not (math.isfinite(va) and math.isfinite(vb)):
                    continue
                n += 1
                if va == vb:
                    bit += 1
                    folds.append(1.0)
                    ratio_norm.append(1.0 / (vals[i + 1] / vals[i]))
                elif va != 0:
                    folds.append(vb / va)
                    ratio_norm.append((vb / va) / (vals[i + 1] / vals[i]))
                else:
                    folds.append(float("inf"))
                    ratio_norm.append(float("inf"))
        if not n:
            continue
        verdict[k] = {
            "n_pairs": n,
            "bit_identical": bit,
            "bit_identical_all": bit == n,
            "min_fold": min(folds),
            "max_fold": max(folds),
            "max_abs_dev_from_1": max(abs(f - 1.0) for f in folds),
            # how far the fold is from being EXACTLY the axis ratio
            "max_abs_dev_from_axis_ratio": max(abs(r - 1.0) for r in ratio_norm),
        }
    return groups, verdict


def axis_ratios(axis, groups):
    rs = set()
    for by_val in groups.values():
        vals = sorted(by_val)
        for i in range(len(vals) - 1):
            rs.add(vals[i + 1] / vals[i])
    return sorted(rs)


def main():
    with open(SWEEP) as fh:
        d = json.load(fh)
    rows = d["rows"]
    keys = numeric_keys(rows)
    out = {"n_rows": len(rows), "alpha": d.get("alpha")}
    for axis in ("plasma_peak", "half_life_s", "total_time_s"):
        groups, verdict = classify(axis, rows, keys)
        out[axis] = {
            "n_groups": len(groups),
            "axis_ratios": axis_ratios(axis, groups),
            "readouts": verdict,
        }
        print(f"\n=== axis {axis}: {len(groups)} groups, "
              f"ratios {axis_ratios(axis, groups)} ===")
        for k, v in sorted(verdict.items(),
                           key=lambda kv: kv[1]["max_abs_dev_from_1"]):
            tag = ("BIT-IDENTICAL" if v["bit_identical_all"] else
                   f"dev1 {v['max_abs_dev_from_1']:.4e} "
                   f"devR {v['max_abs_dev_from_axis_ratio']:.4e}")
            print(f"  {k:34s} n={v['n_pairs']:4d} fold "
                  f"[{v['min_fold']!r}, {v['max_fold']!r}]  {tag}")
    dest = os.path.join(HERE, "axis_exponents.json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote", dest)


if __name__ == "__main__":
    sys.exit(main())
