#!/usr/bin/env python3
"""Measure, on REAL solves, exactly which RunResult leaves are invariant under a
pure dose / half-life change and which scale -- and whether "exact" survives a
ratio that is NOT a power of two.

Why this exists and the surrogate study's sweep does not answer it:

* that sweep varied ``plasma_peak`` 50 -> 200 only, i.e. a ratio of 4.0, which
  is a POWER OF TWO. Scaling an IEEE-754 double by a power of two is exact, so
  a bit-identical shape readout there does not establish that the shape is
  bit-identical under, say, a 1.5x dose change. This probe runs 1.5x and 3x.
  (``app/server/reach.py::SCALE_INVARIANCE_NOTE`` already records this exact
  lesson for ``bbb_penetrance``: "EXACT IN EXACT ARITHMETIC ... an earlier
  wording of this clause said 'bit-identical' about all of them, which is
  false for two".)
* that sweep recorded 44 columns. A RunResult has ~145 numeric leaves, and the
  ones it left out include every THRESHOLD readout
  (``therapeutic_coverage_*``, ``median_E_over_ic50``, ``kill_fraction``,
  ``viable_tumor_fraction``), which are exactly the ones a scaler must refuse.
  Unmeasured is not the same as invariant.

Every run is a real ``app.server.jobs.execute_run`` call on one caller-owned
reconstruction cache, so the reconstruction is built once per section and the
only thing differing between arms is the PK block.

Usage:  probe_scaling.py [dx_um] [total_time_s]
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

HERE = os.path.dirname(os.path.abspath(__file__))

#: (section, compound). Two tissues and three compounds, so no claim below
#: rests on one reconstruction or one descriptor set. ``xenium_breast`` is the
#: one section here on which this repository holds per-cell spatial statistics
#: VALID, so the headline spread readout is measured somewhere it is defended.
CASES = (
    ("gbmspace_at12_bra_4_fo_4_1", "temozolomide"),
    ("gbmspace_at12_bra_4_fo_4_1", "doxorubicin"),
    ("gbmspace_at12_bra_4_fo_4_1", "vincristine"),
    ("xenium_breast", "temozolomide"),
    ("xenium_breast", "doxorubicin"),
)

BASE_PEAK, BASE_HL = 50.0, 14400.0

#: Each arm is a PURE change on one axis (plus one deliberate both-axes arm).
#: 4.0 is a power of two and 1.5 / 3.0 / 1.37 are not -- that contrast is the
#: whole point of the probe.
ARMS = {
    "base":      dict(peak=BASE_PEAK,       hl=BASE_HL),
    "dose_x4":   dict(peak=BASE_PEAK * 4.0, hl=BASE_HL),          # power of 2
    "dose_x1.5": dict(peak=BASE_PEAK * 1.5, hl=BASE_HL),          # NOT
    "dose_x3":   dict(peak=BASE_PEAK * 3.0, hl=BASE_HL),          # NOT
    "hl_x4":     dict(peak=BASE_PEAK,       hl=BASE_HL * 4.0),
    "hl_x1.37":  dict(peak=BASE_PEAK,       hl=BASE_HL * 1.37),
    "both":      dict(peak=BASE_PEAK * 3.0, hl=BASE_HL * 1.37),
}


def flatten(obj, prefix=""):
    """Every numeric leaf of a RunResult body, by dotted path."""
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            if str(k).startswith("_"):
                continue
            out.update(flatten(v, f"{prefix}.{k}" if prefix else str(k)))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            out.update(flatten(v, f"{prefix}[{i}]"))
    elif isinstance(obj, bool) or obj is None or isinstance(obj, str):
        pass
    elif isinstance(obj, (int, float)):
        out[prefix] = float(obj)
    return out


def main():
    dx = float(sys.argv[1]) if len(sys.argv) > 1 else 30.0
    total_time = float(sys.argv[2]) if len(sys.argv) > 2 else 3600.0

    from tmesim.validation.calibration import calibrate_alpha
    alpha0 = calibrate_alpha()
    assert alpha0 == 22.223477711172563, alpha0

    from app.server.cache import ReconstructionCache
    from app.server.jobs import _Ctx, execute_run
    sys.path.insert(0, os.path.join(HERE, "..", "surrogate_feasibility"))
    from solver_runner import _StdoutJob, make_spec  # noqa: E402

    def solve(spec, run_id, cache):
        """The RAW RunResult body -- solver_runner.run_one() returns a
        flattened summary that drops most of the leaves this probe is about."""
        res = execute_run(spec, _Ctx(_StdoutJob(run_id, spec)), cache=cache)
        res.pop("_payloads", None)
        return res

    cases = []
    for section, drug in CASES:
        cache = ReconstructionCache(maxsize=1)
        flat, pmean = {}, {}
        for name, a in ARMS.items():
            spec = make_spec(section, drug, dx, total_time,
                             a["peak"], a["hl"], 1, 86400.0)
            t0 = time.time()
            body = solve(spec, f"{section}:{drug}:{name}", cache)
            flat[name] = flatten(body)
            pmean[name] = flat[name].get("extra_metrics.plasma_mean_uM")
            print(f"  {section:30s} {drug:13s} {name:10s} "
                  f"peak={a['peak']:<7g} hl={a['hl']:<9g} "
                  f"{time.time() - t0:6.2f}s  plasma_mean={pmean[name]!r}")
        cases.append({"section": section, "drug": drug,
                      "plasma_mean_uM": pmean, "leaves": flat})

    alpha1 = calibrate_alpha()
    assert alpha1 == alpha0, (alpha0, alpha1)

    out = {"dx_um": dx, "total_time_s": total_time, "alpha": alpha0,
           "base_peak_uM": BASE_PEAK, "base_half_life_s": BASE_HL,
           "arms": ARMS, "cases": cases}
    dest = os.path.join(HERE, "probe_scaling.json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote", dest, f"({len(cases)} cases x {len(ARMS)} arms)")
    assert math.isfinite(alpha1)


if __name__ == "__main__":
    sys.exit(main())
