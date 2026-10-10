#!/usr/bin/env python3
"""What the exact memo and the exact algebra actually save, on real runs.

A performance claim is a claim, and ``results/performance_profile_2026-09-20.md``
records one RunSpec on unchanged code measuring five different wall-clocks on
this machine. So this reports the RATIO and the load, and solves the same spec
TWICE to show how much of the spread is the machine rather than the change.

The memo's saving is the whole run: a hit does no reconstruction and no solve.
The algebra's saving is also the whole run, and it applies to a spec that has
NEVER been solved -- which is the part a cache cannot do.

Usage: measure_speedup.py [section] [dx_um] [total_time_s]
"""
from __future__ import annotations

import json
import os
import sys
import time

_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    section = sys.argv[1] if len(sys.argv) > 1 else "xenium_breast"
    dx = float(sys.argv[2]) if len(sys.argv) > 2 else 30.0
    window = float(sys.argv[3]) if len(sys.argv) > 3 else 21600.0

    from tmesim.validation.calibration import calibrate_alpha
    alpha0 = calibrate_alpha()
    assert alpha0 == 22.223477711172563, alpha0

    from app.server.jobs import JobManager
    from app.server.memo import answer_spec, engine_fingerprint

    engine_fingerprint()          # one-off, so it is not billed to a lookup

    root = os.path.join(HERE, "_runs")
    os.makedirs(root, exist_ok=True)
    mgr = JobManager(root=root)

    def spec(peak):
        return {
            "section": section, "dx_um": dx, "drug": "temozolomide",
            "sensitivity_mode": "uniform", "vascular_mode": "discrete",
            "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": None,
            "pk": {"plasma_peak": peak, "half_life_s": 14400.0,
                   "n_doses": 1, "dose_interval_s": 86400.0},
            "total_time_s": window, "attribute_mechanism": True,
            "timelapse_frames": 0,
        }

    def solve(s, label):
        t0 = time.time()
        job = mgr.submit(s)
        while job.status not in ("done", "failed", "cancelled"):
            time.sleep(0.02)
        assert job.status == "done", job.error
        dt = time.time() - t0
        print(f"  solve {label:22s} {dt:8.3f}s  (run {job.run_id})")
        return dt, job.run_id

    def lookup(s, label, allow_scaling=False, n=20):
        best = None
        for _ in range(n):
            t0 = time.time()
            out = answer_spec(mgr, s, allow_scaling=allow_scaling)
            dt = time.time() - t0
            best = dt if best is None else min(best, dt)
        print(f"  {label:28s} {best * 1000:8.3f}ms  mode={out['mode']}")
        return best, out

    print(f"section {section}  dx {dx} um  window {window} s  "
          f"loadavg {os.getloadavg()}")
    try:
        first, _ = solve(spec(50.0), "cold (reconstruct)")
        second, _ = solve(spec(50.0), "warm (recon cached)")
        third, _ = solve(spec(50.0), "warm again")
        memo_s, memo_out = lookup(spec(50.0), "exact memo hit")
        scale_s, scale_out = lookup(spec(150.0), "algebraic 3x dose",
                                    allow_scaling=True)
        dose_solve, _ = solve(spec(150.0), "3x dose, solved")

        out = {
            "section": section, "dx_um": dx, "total_time_s": window,
            "alpha": alpha0, "loadavg": list(os.getloadavg()),
            "solve_cold_s": first,
            "solve_warm_s": [second, third],
            "memo_hit_s": memo_s,
            "memo_speedup_vs_warm_solve": min(second, third) / memo_s,
            "memo_speedup_vs_cold_solve": first / memo_s,
            "scaled_answer_s": scale_s,
            "scaled_speedup_vs_its_own_solve": dose_solve / scale_s,
            "unsolved_dose_solve_s": dose_solve,
            "memo_mode": memo_out["mode"],
            "scale_mode": scale_out["mode"],
            "scale_bit_exact": scale_out.get("derivation", {}).get("bit_exact"),
        }
        print("\n" + json.dumps(out, indent=1))
        dest = os.path.join(HERE, f"speedup_dx{int(dx)}.json")
        with open(dest, "w") as fh:
            json.dump(out, fh, indent=1)
        print("wrote", dest)
    finally:
        mgr.shutdown()
    assert calibrate_alpha() == alpha0


if __name__ == "__main__":
    sys.exit(main())
