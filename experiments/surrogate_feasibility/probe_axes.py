"""Scratch probe 2: what does the solve actually cost, and what moves the targets?

Answers, before any model is fitted:

* WARM solve cost. A surrogate is only worth anything if the thing it replaces
  is slow. The reconstruction is cached per section in the Studio, so the cost
  a preview would actually hide is the solve on a WARM cache -- measured here at
  the shipped flagship resolution (dx = 12 um), not just at the cheap dx = 30.
* Which axes move the candidate targets at all. Probe 1 showed ``plasma_peak``
  is exactly scalar on the exposure field, so the remaining candidate axes are
  the COMPOUND's transport descriptors and the PK time shape (half-life, number
  of doses) and the exposure window.

Writes nothing. Prints.
"""
from __future__ import annotations

import sys
import time

from solver_runner import assert_alpha, make_spec, run_one

SECTION = sys.argv[1] if len(sys.argv) > 1 else "xenium_319"
DX = float(sys.argv[2]) if len(sys.argv) > 2 else 12.0


def show(tag, out):
    if out["status"] != "ok":
        print(f"  {tag:34s} FAILED {out['error_type']}: {out['error'][:160]}")
        return
    r = out["E_p90_uM"] / out["E_p10_uM"] if out["E_p10_uM"] else float("nan")
    print(f"  {tag:34s} wall={out['wall_clock_s']:7.2f}s "
          f"hit={str(out['recon_cache_hit'])[:5]:5s} "
          f"E_p90/p10={r:9.4f} L={out['penetration_length_um']!s:>20.20s} "
          f"pen_idx={out['penetration_index']!s:>8.6s} "
          f"med/peak={out['median_tumor_auc'] / out['plasma_peak_uM']:12.4g} "
          f"F={out['reach_headline_F']}")


def main() -> int:
    from app.server.cache import ReconstructionCache

    assert_alpha("probe2 start")
    cache = ReconstructionCache(maxsize=1)

    print(f"\n### {SECTION} @ dx={DX:g}  — cold then warm, default compound")
    base = dict(section=SECTION, drug="temozolomide", dx_um=DX,
                total_time_s=21600.0, plasma_peak=50.0, half_life_s=14400.0,
                n_doses=1, dose_interval_s=86400.0)
    t0 = time.time()
    out = run_one(SECTION, make_spec(**base), cache, "cold")
    show("cold (recon + solve)", out)
    if out["status"] != "ok":
        print(out.get("traceback", "")[-1500:])
        return 1
    print(f"    recon_build_s={out['recon_build_s']:.2f} "
          f"grid={out['grid']} n_cells={out['n_cells']}")
    cold_total = time.time() - t0

    # Warm solves, repeated, so the solve cost is a measurement not one sample.
    warm = []
    for i in range(3):
        o = run_one(SECTION, make_spec(**base), cache, f"warm{i}")
        warm.append(o["wall_clock_s"])
    print(f"    WARM solve wall: {['%.2f' % w for w in warm]} "
          f"min={min(warm):.2f}s  (cold total {cold_total:.1f}s)")

    print(f"\n### compound axis (same section, warm cache)")
    for drug in ("temozolomide", "cisplatin", "doxorubicin", "lomustine",
                 "carmustine", "paclitaxel", "bevacizumab"):
        spec = make_spec(**{**base, "drug": drug})
        o = run_one(SECTION, spec, cache, drug)
        if o["status"] == "ok":
            print(f"    D={o['drug_diffusivity_um2_s']:9.4g} "
                  f"bbb={o['drug_bbb_penetrance']:7.4f} "
                  f"clr={o['drug_clearance_1_per_s']:9.4g} "
                  f"bind={o['drug_binding_1_per_s']:9.4g}", end="")
        show(drug, o)

    print(f"\n### PK time shape + window (temozolomide, warm cache)")
    for hl, nd, tt in ((3600.0, 1, 21600.0), (14400.0, 1, 21600.0),
                       (57600.0, 1, 21600.0), (14400.0, 3, 21600.0),
                       (14400.0, 1, 3600.0), (14400.0, 1, 86400.0)):
        spec = make_spec(**{**base, "half_life_s": hl, "n_doses": nd,
                            "total_time_s": tt})
        o = run_one(SECTION, spec, cache, f"hl{hl:g}_n{nd}_t{tt:g}")
        show(f"hl={hl:g} n={nd} T={tt:g}", o)

    assert_alpha("probe2 end")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
