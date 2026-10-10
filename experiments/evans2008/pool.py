"""Pool the per-section records and apply the pre-registered decision rule.

Pre-registration: ``results/evans2008_perivascular_sign_2026-10-10.md`` section 4.3.
The rule is applied exactly as written; no branch is reinterpreted here.
"""
from __future__ import annotations

import glob
import json
import os

from scipy.stats import beta

EVANS_K, EVANS_N = 40, 69
EVANS_CI = (0.4548, 0.6976)
PRIMARY = ("xenium_319", "xenium_io", "gse324268")


def cp(k: int, n: int):
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return lo, hi


def load():
    runs = {}
    for p in sorted(glob.glob(os.path.join("experiments", "evans2008", "out", "*.json"))):
        with open(p) as fh:
            runs.update(json.load(fh).get("runs", {}))
    return runs


def branch(p: float) -> str:
    if p >= 0.90 and p > EVANS_CI[1]:
        return "R1 REFUTED-TOO-CLEAN"
    if EVANS_CI[0] <= p <= EVANS_CI[1]:
        return "R2 CONSISTENT"
    if EVANS_CI[1] < p < 0.90:
        return "R3 DISCREPANT-NOT-DETERMINISTIC"
    if p < EVANS_CI[0]:
        return "R4 REFUTED-TOO-NOISY"
    return "unclassified"


def main() -> None:
    runs = load()
    print(f"Evans patient F: {EVANS_K}/{EVANS_N} = {EVANS_K/EVANS_N:.4f}  "
          f"CP95 = [{EVANS_CI[0]:.4f}, {EVANS_CI[1]:.4f}]\n")

    print(f"{'run':44s} {'n_an':>6s} {'n_exp':>6s} {'n_anti':>6s} {'P_exp':>7s} "
          f"{'CP95 lo':>8s} {'CP95 hi':>8s} {'classes':>8s}")
    for k in sorted(runs):
        r = runs[k]
        if r.get("status") != "ok":
            print(f"{k:44s} {r.get('status')}  {r.get('reason') or r.get('void_reason','')}")
            continue
        b = r["bands_1_5"]
        lo, hi = b["p_expected_cp95"]
        print(f"{k:44s} {r['n_analysable']:6d} {b['n_expected_sign']:6d} "
              f"{b['n_anti_gradient']:6d} {b['p_expected']:7.4f} {lo:8.4f} {hi:8.4f} "
              f"{r['n_distinct_evans_classes']:8d}")

    # --- the pre-registered primary: three Xenium sections, dx=8, linear, mc3 --- #
    print("\n=== PRIMARY (prereg 4.3): 3 Xenium sections, dx=8 um, linear sink, "
          "all vessels, min_count=3 ===")
    tot_e = tot_a = 0
    classes = set()
    for sec in PRIMARY:
        k = f"{sec}@dx8@linear@mc3"
        r = runs.get(k)
        if r is None or r.get("status") != "ok":
            print(f"  {sec}: MISSING/not ok -> primary cannot be pooled as registered")
            return
        b = r["bands_1_5"]
        tot_e += b["n_expected_sign"]
        tot_a += b["n_anti_gradient"]
        classes |= set(r["class_counts"].keys())
        print(f"  {sec:12s} {b['n_expected_sign']:5d} expected / "
              f"{b['n_anti_gradient']:5d} anti  P_exp = {b['p_expected']:.4f}  "
              f"[{b['p_expected_cp95'][0]:.4f}, {b['p_expected_cp95'][1]:.4f}]  "
              f"-> {branch(b['p_expected'])}")
    n = tot_e + tot_a
    p = tot_e / n
    lo, hi = cp(tot_e, n)
    print(f"\n  POOLED  {tot_e}/{n} = P_exp = {p:.4f}  CP95 = [{lo:.4f}, {hi:.4f}]")
    print(f"  PRE-REGISTERED BRANCH: {branch(p)}")
    print(f"  distinct Evans Table 1 classes spanned, pooled: {len(classes)} "
          f"-> {sorted(classes)}")
    print(f"  R6 (class stratification void): "
          f"{'FIRES' if len(classes) == 1 else 'does not fire'}")

    # bands 2-5 arm, pooled
    e2 = sum(runs[f"{s}@dx8@linear@mc3"]["bands_2_5"]["n_expected_sign"] for s in PRIMARY)
    a2 = sum(runs[f"{s}@dx8@linear@mc3"]["bands_2_5"]["n_anti_gradient"] for s in PRIMARY)
    print(f"\n  bands 2-5 pooled: {e2}/{e2+a2} = {e2/(e2+a2):.4f} -> {branch(e2/(e2+a2))}")

    # >=4 endothelial cells arm, pooled
    e4 = a4 = 0
    for s in PRIMARY:
        g = runs[f"{s}@dx8@linear@mc3"].get("bands_1_5_vessels_ge4_cells", {})
        e4 += g.get("n_expected_sign", 0)
        a4 += g.get("n_anti_gradient", 0)
    if e4 + a4:
        print(f"  vessels with >=4 endothelial cells, pooled: {e4}/{e4+a4} = "
              f"{e4/(e4+a4):.4f} -> {branch(e4/(e4+a4))}")

    # perivascular dynamic range, pooled view
    print("\n=== perivascular pO2 range on the primary arm (mmHg) ===")
    for sec in PRIMARY:
        r = runs[f"{sec}@dx8@linear@mc3"]
        pv = r["perivascular_mmhg"]
        sl = r["slope_mmhg_per_um"]
        print(f"  {sec:12s} field {r['field_mmhg_min']:.3f}-{r['field_mmhg_max']:.3f}  "
              f"band1 med {pv['band1_median']:.3f}  band5 med {pv['band5_median']:.3f}  "
              f"median(b1-b5) {pv['band1_minus_band5_median']:+.4f}  "
              f"max|b1-b5| {pv['band1_minus_band5_max_abs']:.4f}  "
              f"median|slope| {sl['median_abs']:.6f} mmHg/um")


if __name__ == "__main__":
    main()
