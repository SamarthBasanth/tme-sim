#!/usr/bin/env python3
"""Score every algebraically-scaled readout against the run it replaces.

Reads the store ``measure_speedup.py`` left behind in ``_runs`` (which holds a
1x and a 3x dose solve of the same spec), asks ``scale_lookup`` for the 3x
answer WITHOUT letting it see the 3x run, and reports the realised error per
readout in units in the last place. So run ``measure_speedup.py`` first;
``_runs`` is a real ``JobManager`` root, is gitignored, and can be deleted.

3.0 is NOT a power of two, so this is the case where the scaled answer is exact
in exact arithmetic and only round-off-close in floating point. The point of
the table is to show how close, measured, rather than asserted.
"""
from __future__ import annotations

import json
import math
import os
import sys

_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

HERE = os.path.dirname(os.path.abspath(__file__))


def ulps(a: float, b: float) -> float:
    if a == b:
        return 0.0
    m = max(abs(a), abs(b))
    if m == 0.0 or not math.isfinite(m):
        return float("inf")
    return abs(a - b) / math.ulp(m)


def main():
    from tmesim.validation.calibration import calibrate_alpha
    alpha0 = calibrate_alpha()
    assert alpha0 == 22.223477711172563, alpha0

    from app.server.jobs import JobManager
    from app.server.memo import scale_lookup

    root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "_runs")
    mgr = JobManager(root=root, start_worker=False)
    try:
        peaks = {}
        for job in mgr.list():
            if job.status != "done":
                continue
            peaks.setdefault(
                float((job.spec.get("pk") or {}).get("plasma_peak")), job)
        if len(peaks) < 2:
            print("need a store holding two doses of one spec", file=sys.stderr)
            return 1
        base_peak = min(peaks)
        target_peak = max(peaks)
        target = peaks[target_peak]

        # Hide the target run so the scaler cannot answer it from the memo.
        hidden = target.status
        target.status = "hidden-for-this-measurement"
        try:
            out = scale_lookup(mgr, target.spec)
        finally:
            target.status = hidden
        assert out["mode"] == "algebraic_scale", out.get("reason")

        rows = []
        for block, leaves in out["readouts"].items():
            for leaf, rec in leaves.items():
                truth = (target.body.get(block) or {}).get(leaf)
                if rec.get("value") is None or truth is None:
                    continue
                rows.append({
                    "readout": f"{block}.{leaf}",
                    "basis": rec["basis"],
                    "declared_max_ulp": rec["max_measured_ulp"],
                    "scaled": rec["value"],
                    "solved": float(truth),
                    "ulp": ulps(rec["value"], float(truth)),
                    "rel": (abs(rec["value"] - float(truth))
                            / abs(float(truth)) if truth else 0.0),
                })
        rows.sort(key=lambda r: -r["ulp"])
        report = {
            "section": target.spec["section"],
            "dx_um": target.spec["dx_um"],
            "total_time_s": target.spec["total_time_s"],
            "base_plasma_peak_uM": base_peak,
            "scaled_to_plasma_peak_uM": target_peak,
            "dose_ratio": target_peak / base_peak,
            "dose_ratio_is_power_of_two": bool(
                math.frexp(target_peak / base_peak)[0] == 0.5),
            "bit_exact_claimed": out["derivation"]["bit_exact"],
            "source_run_id": out["source_run_id"],
            "n_readouts_served": len(rows),
            "worst_ulp": max((r["ulp"] for r in rows), default=0.0),
            "worst_relative": max((r["rel"] for r in rows), default=0.0),
            "n_bit_identical": sum(1 for r in rows if r["ulp"] == 0.0),
            "rows": rows,
            "refused": sorted(out["refused"]),
            "alpha": alpha0,
        }
        print(f"{'readout':46s} {'basis':14s} {'ulp':>6s} {'declared':>9s}  relative")
        for r in rows:
            print(f"{r['readout']:46s} {r['basis']:14s} {r['ulp']:6.1f} "
                  f"{r['declared_max_ulp']:9.1f}  {r['rel']:.3e}")
        print(f"\n{report['n_bit_identical']} of {len(rows)} served readouts "
              f"are BIT-IDENTICAL to the solve; worst "
              f"{report['worst_ulp']:.1f} ulp "
              f"({report['worst_relative']:.3e} relative) at a dose ratio of "
              f"{report['dose_ratio']} "
              f"(power of two: {report['dose_ratio_is_power_of_two']})")
        print(f"refused: {len(report['refused'])} -> {report['refused']}")
        with open(os.path.join(HERE, "verify_against_solve.json"), "w") as fh:
            json.dump(report, fh, indent=1)
    finally:
        mgr.shutdown()
    assert calibrate_alpha() == alpha0
    return 0


if __name__ == "__main__":
    sys.exit(main())
