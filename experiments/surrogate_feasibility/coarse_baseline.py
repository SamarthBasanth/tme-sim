"""The baseline a learned surrogate actually has to beat: a COARSER SOLVE.

A surrogate is one way to answer "show the user something now and solve
properly behind it". The other way is to run the SAME solver on a coarser grid.
That option has properties no learned model can have:

* its output is the engine's, so it carries the engine's own caveats and
  refusals, and cannot disagree with the real run about what is defensible;
* it needs no training set, so it cannot be stale, and no new section, compound
  or engine change invalidates it;
* it makes no AI claim.

Its cost is a discretisation error, and that error is MEASURABLE. This script
measures it directly, in the readouts a preview would show: the same section
and compound solved at the preview resolution (dx = 30 um) and at the shipped
flagship resolution (dx = 12 um), compared pair by pair.

If the coarse solve already matches the fine solve to within a few per cent,
then the ceiling on a useful surrogate is not "how well can a model learn
this?" but "can a model beat a few per cent, out of sample, on a new donor?" --
and that is the question the study has to answer.

WHAT "AGREES" MEANS HERE, EXACTLY
---------------------------------
NOT "the coarse solve approximates the true answer". ``results/
CLAIM_INVENTORY.md`` sec. 1f RETRACTS claim N2: dx = 12 um is NOT grid
converged -- the GBM ladder is non-monotone with observed orders 0.96 and
-3.33, and "neither panel is in an asymptotic regime". There is no converged
answer on offer for either resolution to approximate.

The question a preview has to answer is narrower and does not need one: does
the number shown now match THE NUMBER THE PRODUCT WILL PRINT when the real run
finishes? That is a comparison between two things that both exist. A preview
that disagrees with the run it precedes is broken whether or not either is
right about the tissue.

TIMING IS A CLAIM
-----------------
``results/performance_profile_2026-09-20.md`` sec. 0.1-0.2: one RunSpec on
unchanged code measured FIVE different wall-clocks on this machine (47.1 s to
174.3 s at dx = 12, under loadavg 4.6-48). So every timing here is recorded
with the load average it was taken under, and no speedup is quoted without it.

Writes ``experiments/surrogate_feasibility/coarse_baseline.json``.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from typing import Dict, List, Optional

from solver_runner import assert_alpha, make_spec, run_one, section_record

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_HERE, "coarse_baseline.json")

PREVIEW_DX = 30.0
SHIPPED_DX = 12.0
COMPOUNDS = ["temozolomide", "cisplatin", "doxorubicin", "lomustine",
             "methotrexate", "vincristine"]
PK = dict(plasma_peak=50.0, half_life_s=14400.0, n_doses=1,
          dose_interval_s=86400.0, total_time_s=21600.0)

#: xenium_319 and xenium_breast are the two sections here where the project
#: holds per-cell spatial statistics VALID and whose dx = 12 um solve is
#: affordable. The other two valid sections (xenium_io 816,769 cells,
#: gse324268 851,509 cells) are ~20x the flagship's grid cells at dx = 12 and
#: are not run: that is a budget limit on this comparison, stated rather than
#: worked around.
DEFAULT_SECTIONS = [
    "xenium_319",                    # flagship, spatial_spread VALID
    "gbmspace_at10_bra_5_fo_1_1",
    "gbmspace_at13_bra_3_fo_1_1",
    "visium_gbm2",
    "xenium_breast",                 # spatial_spread VALID, 2.9x flagship
]

TARGETS = ("E_p90_over_p10", "penetration_length_um", "penetration_index",
           "median_auc_over_peak")


def _f(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def readouts(o: Dict) -> Dict:
    e10, e90 = _f(o.get("E_p10_uM")), _f(o.get("E_p90_uM"))
    med, pk = _f(o.get("median_tumor_auc")), _f(o.get("plasma_peak_uM"))
    return {
        "E_p90_over_p10": (e90 / e10) if (e10 and e90 and e10 > 0) else None,
        "penetration_length_um": _f(o.get("penetration_length_um")),
        "penetration_index": _f(o.get("penetration_index")),
        "median_auc_over_peak": (med / pk) if (med is not None and pk) else None,
    }


def main(argv: Optional[List[str]] = None) -> int:
    from app.server.cache import ReconstructionCache

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sections", default=",".join(DEFAULT_SECTIONS))
    args = p.parse_args(argv)
    sections = [s for s in args.sections.split(",") if s]

    assert_alpha("coarse start")
    blob: Dict = {"alpha": None, "preview_dx_um": PREVIEW_DX,
                  "shipped_dx_um": SHIPPED_DX, "pk": PK,
                  "compounds": COMPOUNDS, "pairs": [], "timing": {}}
    if os.path.exists(_OUT):
        with open(_OUT) as fh:
            blob = json.load(fh)
    blob["alpha"] = 22.223477711172563

    done = {(r["section"], r["drug"]) for r in blob["pairs"]}
    t_all = time.time()
    for section in sections:
        rec = section_record(section)
        todo = [d for d in COMPOUNDS if (section, d) not in done]
        if not todo:
            print(f"=== {section}: complete", flush=True)
            continue
        print(f"\n=== {section}  src={rec.get('source')}  "
              f"({len(todo)} compounds x 2 resolutions)", flush=True)
        for dx in (PREVIEW_DX, SHIPPED_DX):
            # One fresh cache per (section, dx): dx is part of the cache key, so
            # a shared cache of size 1 would thrash between the two.
            cache = ReconstructionCache(maxsize=1)
            for drug in todo:
                spec = make_spec(section, drug, dx, PK["total_time_s"],
                                 PK["plasma_peak"], PK["half_life_s"],
                                 PK["n_doses"], PK["dose_interval_s"])
                out = run_one(section, spec, cache, f"{section[:10]}/{drug[:5]}")
                key = "coarse" if dx == PREVIEW_DX else "fine"
                hit = next((r for r in blob["pairs"]
                            if r["section"] == section and r["drug"] == drug),
                           None)
                if hit is None:
                    hit = {"section": section, "drug": drug,
                           "source": rec.get("source"),
                           "spatial_spread_valid": None}
                    blob["pairs"].append(hit)
                if out["status"] != "ok":
                    hit[key] = {"status": "FAILED",
                                "error": f"{out['error_type']}: "
                                         f"{out['error'][:200]}"}
                    print(f"    !! {drug} dx={dx:g}: {out['error'][:120]}",
                          flush=True)
                    continue
                hit[key] = {
                    "status": "ok",
                    "wall_clock_s": out["wall_clock_s"],
                    "recon_build_s": out["recon_build_s"],
                    "recon_cache_hit": out["recon_cache_hit"],
                    "grid": out["grid"],
                    # performance_profile_2026-09-20.md sec. 0.2: a wall-clock
                    # on this machine is meaningless without the load it was
                    # taken under.
                    "loadavg_1m": os.getloadavg()[0],
                    "n_cpu": os.cpu_count(),
                    **readouts(out),
                }
                print(f"    {drug:14s} dx={dx:>4.0f} "
                      f"{out['wall_clock_s']:7.2f}s grid={out['grid']} "
                      f"E_p90/p10={readouts(out)['E_p90_over_p10']} "
                      f"L={out['penetration_length_um']}", flush=True)
        with open(_OUT, "w") as fh:
            json.dump(blob, fh, indent=1, default=str)

    # ---- the comparison --------------------------------------------------- #
    summary: Dict = {}
    for t in TARGETS:
        errs, rel, rows = [], [], 0
        for pr in blob["pairs"]:
            c, f = pr.get("coarse") or {}, pr.get("fine") or {}
            if c.get("status") != "ok" or f.get("status") != "ok":
                continue
            cv, fv = _f(c.get(t)), _f(f.get(t))
            if cv is None or fv is None or fv == 0:
                continue
            rows += 1
            errs.append(abs(cv - fv))
            rel.append(abs(cv - fv) / abs(fv))
        if rows:
            errs_s, rel_s = sorted(errs), sorted(rel)
            summary[t] = {
                "n_pairs": rows,
                "mae": sum(errs) / rows,
                "median_ae": errs_s[rows // 2],
                "max_ae": errs_s[-1],
                "median_rel_err": rel_s[rows // 2],
                "max_rel_err": rel_s[-1],
                "mean_rel_err": sum(rel) / rows,
            }
    blob["coarse_vs_fine"] = summary
    # Cost of each resolution, measured on the runs just done.
    for key in ("coarse", "fine"):
        w = [pr[key]["wall_clock_s"] for pr in blob["pairs"]
             if (pr.get(key) or {}).get("status") == "ok"]
        warm = [pr[key]["wall_clock_s"] for pr in blob["pairs"]
                if (pr.get(key) or {}).get("status") == "ok"
                and pr[key].get("recon_cache_hit")]
        if w:
            ws = sorted(w)
            la = [pr[key].get("loadavg_1m") for pr in blob["pairs"]
                  if (pr.get(key) or {}).get("status") == "ok"
                  and pr[key].get("loadavg_1m") is not None]
            blob["timing"][key] = {
                "n": len(w), "median_s": ws[len(ws) // 2],
                "min_s": ws[0], "max_s": ws[-1],
                "warm_median_s": (sorted(warm)[len(warm) // 2]
                                  if warm else None),
                "warm_n": len(warm),
                "loadavg_1m_min": min(la) if la else None,
                "loadavg_1m_max": max(la) if la else None,
                "n_cpu": os.cpu_count(),
            }
    blob["elapsed_s"] = time.time() - t_all
    blob["alpha_end"] = assert_alpha("coarse end")
    with open(_OUT, "w") as fh:
        json.dump(blob, fh, indent=1, default=str)

    print("\n=== coarse (dx=30) vs fine (dx=12), same solver ===")
    for t, s in summary.items():
        print(f"  {t:26s} n={s['n_pairs']:3d} "
              f"median_rel_err={s['median_rel_err'] * 100:7.3f}%  "
              f"max_rel_err={s['max_rel_err'] * 100:7.3f}%  "
              f"MAE={s['mae']:.4g}")
    print(f"  timing: {json.dumps(blob['timing'])}")
    print(f"\n[coarse] -> {_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
