#!/usr/bin/env python3
"""Classify every RunResult leaf from ``probe_scaling.py`` per axis, ACROSS all
cases, into

  BIT-IDENTICAL         identical doubles in every arm of that axis, every case
  SCALES EXACTLY        fold == the arm's realised plasma-scalar ratio exactly
  INVARIANT to <= N ulp fold == 1.0 to within N units in the last place
  SCALES to <= N ulp    fold == the ratio to within N ulp
  NEITHER               anything else -- a candidate for REFUSAL

The scale a magnitude readout is PREDICTED to take is the ratio of the
engine's own PK scalar ``plasma_mean``, because in the shipped
(``time_resolved=False``) path that single float is the only way any PK number
reaches the transport solve (``tmesim/sim/engine.py``: ``source =
delivery_source(work, drug, plasma_mean)``).
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

DOSE_ARMS = ("dose_x4", "dose_x1.5", "dose_x3")
HL_ARMS = ("hl_x4", "hl_x1.37")
BOTH_ARMS = ("both",)


def ulps(a: float, b: float) -> float:
    """|a-b| in units in the last place of max(|a|,|b|)."""
    if a == b:
        return 0.0
    m = max(abs(a), abs(b))
    if m == 0.0 or not math.isfinite(m):
        return float("inf")
    return abs(a - b) / math.ulp(m)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, "probe_scaling.json")
    with open(path) as fh:
        rep = json.load(fh)

    cases = rep["cases"]
    leaves = sorted({k for c in cases for f in c["leaves"].values() for k in f})

    axes = {"dose": DOSE_ARMS, "half_life": HL_ARMS, "both": BOTH_ARMS}
    verdict = {}

    for leaf in leaves:
        verdict[leaf] = {}
        for axis, arms in axes.items():
            n = 0
            bit = True
            inv_ulps = 0.0
            scale_exact = True
            scale_ulps = 0.0
            folds = []
            for c in cases:
                base = c["leaves"]["base"].get(leaf)
                pm0 = c["plasma_mean_uM"]["base"]
                for arm in arms:
                    v = c["leaves"][arm].get(leaf)
                    if base is None or v is None:
                        continue
                    if not (math.isfinite(base) and math.isfinite(v)):
                        continue
                    n += 1
                    if v != base:
                        bit = False
                    if base != 0.0:
                        fold = v / base
                        folds.append(fold)
                        inv_ulps = max(inv_ulps, ulps(fold, 1.0))
                        ratio = c["plasma_mean_uM"][arm] / pm0
                        if fold != ratio:
                            scale_exact = False
                        scale_ulps = max(scale_ulps, ulps(fold, ratio))
                    else:
                        folds.append(None if v != 0 else 1.0)
                        if v != 0:
                            bit = False
                            scale_exact = False
                            inv_ulps = float("inf")
                            scale_ulps = float("inf")
            if not n:
                continue
            klass = "NEITHER"
            if bit:
                klass = "BIT_IDENTICAL"
            elif scale_exact:
                klass = "SCALES_EXACTLY"
            elif inv_ulps <= 4:
                klass = f"INVARIANT_{int(math.ceil(inv_ulps))}ulp"
            elif scale_ulps <= 4:
                klass = f"SCALES_{int(math.ceil(scale_ulps))}ulp"
            verdict[leaf][axis] = {
                "class": klass, "n": n,
                "max_ulp_from_1": inv_ulps, "max_ulp_from_ratio": scale_ulps,
                "min_fold": min((f for f in folds if f is not None), default=None),
                "max_fold": max((f for f in folds if f is not None), default=None),
            }

    for axis in axes:
        print(f"\n################ axis {axis} ################")
        by_class = {}
        for leaf, v in verdict.items():
            if axis not in v:
                continue
            by_class.setdefault(v[axis]["class"], []).append((leaf, v[axis]))
        for klass in sorted(by_class):
            rows = by_class[klass]
            print(f"\n-- {klass}  ({len(rows)} leaves)")
            for leaf, rec in sorted(rows):
                if klass == "BIT_IDENTICAL":
                    print(f"     {leaf}")
                else:
                    print(f"     {leaf:48s} n={rec['n']:2d} "
                          f"ulp(1)={rec['max_ulp_from_1']:.3g} "
                          f"ulp(ratio)={rec['max_ulp_from_ratio']:.3g} "
                          f"folds [{rec['min_fold']!r}, {rec['max_fold']!r}]")

    dest = os.path.join(HERE, "probe_classification.json")
    with open(dest, "w") as fh:
        json.dump({"source": os.path.basename(path),
                   "n_cases": len(cases),
                   "cases": [(c["section"], c["drug"]) for c in cases],
                   "verdict": verdict}, fh, indent=1)
    print("\nwrote", dest)


if __name__ == "__main__":
    sys.exit(main())
