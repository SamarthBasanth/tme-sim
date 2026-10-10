"""The reach LADDER across every section, at one PK shape — and why one is enough.

The main sweep (``sweep.py``) recorded the reach HEADLINE ``F(0.1)`` correctly
but read the full ladder off the wrong level of the payload, so 822 rows carry
an empty ladder. ``solver_runner`` is fixed; this re-collects the ladder rather
than re-running the whole 822-row design.

ONE PK SHAPE IS SUFFICIENT, AND THAT IS A MEASUREMENT
-----------------------------------------------------
Not a shortcut. From ``analyse_results.json`` on the main sweep, the maximum
fold change in ``reach_headline_F`` is EXACTLY 1.0 across all three PK axes
varied -- ``half_life_s`` (n=204 comparison groups, 16x range), ``plasma_peak``
(n=6, 4x range) and ``total_time_s`` (n=143, 6x range). The reason is in
``app/server/reach.py``: the shipped map is the STEADY-STATE one
(``window_mode='steady_state'``), and a steady-state, self-normalised ratio
cannot depend on the dosing schedule. So the ladder is a function of
(section, compound) alone, and sweeping PK again would add 3 duplicate rows per
cell.

That invariance is RE-ASSERTED here rather than assumed: ``--check-invariance``
runs a second PK shape on a few sections and fails if any ladder rung moves.

Writes ``experiments/surrogate_feasibility/ladder_data.json``.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from typing import Dict, List, Optional

from solver_runner import (
    assert_alpha,
    make_spec,
    run_one,
    section_record,
    spatial_spread_valid,
)
from sweep import COMPOUNDS, DX_UM, donor_of

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_HERE, "ladder_data.json")

#: The shipped PKConfig default.
SHIPPED = dict(plasma_peak=50.0, half_life_s=14400.0, n_doses=1,
               dose_interval_s=86400.0, total_time_s=21600.0)
#: A deliberately different shape, for the invariance re-assertion only.
ALT = dict(plasma_peak=200.0, half_life_s=3600.0, n_doses=1,
           dose_interval_s=86400.0, total_time_s=3600.0)


def load() -> Dict:
    if os.path.exists(_OUT):
        with open(_OUT) as fh:
            return json.load(fh)
    return {"alpha": None, "dx_um": DX_UM, "shipped_pk": SHIPPED,
            "alt_pk": ALT, "compounds": COMPOUNDS, "rows": [],
            "invariance": []}


def save(b: Dict) -> None:
    tmp = _OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(b, fh, default=str)
    os.replace(tmp, _OUT)


def _row(section: str, rec: Dict, drug: str, pk: Dict, out: Dict) -> Dict:
    r = {
        "section": section, "donor": donor_of(section),
        "source": rec.get("source"),
        "spatial_spread_valid": spatial_spread_valid(rec),
        "drug": drug, "dx_um": DX_UM, **{f"pk_{k}": v for k, v in pk.items()},
    }
    for k in ("status", "wall_clock_s", "reach_null_reason",
              "reach_headline_x", "reach_headline_F", "reach_ladder_x",
              "reach_ladder_F", "reach_window_mode", "reach_E_ref",
              "n_viable_reach", "error_type", "error",
              "n_cells", "n_tumor", "grid", "extent_um",
              "frac_TUMOR", "frac_STROMAL", "frac_IMMUNE",
              "frac_ENDOTHELIAL", "frac_NECROTIC", "mean_vascular_density",
              "mean_ecm_density", "mean_oxygen", "tissue_coverage_fraction",
              "max_darcy_speed_um_s", "drug_diffusivity_um2_s",
              "drug_binding_1_per_s", "drug_clearance_1_per_s",
              "drug_bbb_penetrance", "drug_mw", "drug_logp",
              "penetration_length_analytic_um"):
        r[k] = out.get(k)
    return r


def main(argv: Optional[List[str]] = None) -> int:
    from app.server.cache import ReconstructionCache, all_sections

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sections", default="")
    p.add_argument("--check-invariance", default="xenium_319,visium_gbm2",
                   help="comma list; a second PK shape is run on these and "
                        "every ladder rung must be identical")
    args = p.parse_args(argv)

    blob = load()
    blob["alpha"] = assert_alpha("ladder start")
    avail = list(all_sections().keys())
    sections = ([s for s in args.sections.split(",") if s] if args.sections
                else avail)
    inv = {s for s in args.check_invariance.split(",") if s}

    t0 = time.time()
    done = {(r["section"], r["drug"]) for r in blob["rows"]}
    for section in sections:
        if section not in avail:
            print(f"=== {section}: not in registry", flush=True)
            continue
        todo = [d for d in COMPOUNDS if (section, d) not in done]
        if not todo and section not in inv:
            print(f"=== {section}: complete", flush=True)
            continue
        rec = section_record(section)
        print(f"\n=== {section}  donor={donor_of(section)}  "
              f"valid={spatial_spread_valid(rec)}", flush=True)
        cache = ReconstructionCache(maxsize=1)
        for drug in todo:
            spec = make_spec(section, drug, DX_UM, SHIPPED["total_time_s"],
                             SHIPPED["plasma_peak"], SHIPPED["half_life_s"],
                             SHIPPED["n_doses"], SHIPPED["dose_interval_s"])
            out = run_one(section, spec, cache, f"{section[:12]}/{drug[:5]}")
            row = _row(section, rec, drug, SHIPPED, out)
            blob["rows"].append(row)
            print(f"    {drug:14s} {out.get('wall_clock_s', 0):6.2f}s "
                  f"F={row['reach_headline_F']} "
                  f"ladder={row['reach_ladder_F']}", flush=True)
        if section in inv:
            for drug in COMPOUNDS[:2]:
                spec = make_spec(section, drug, DX_UM, ALT["total_time_s"],
                                 ALT["plasma_peak"], ALT["half_life_s"],
                                 ALT["n_doses"], ALT["dose_interval_s"])
                out = run_one(section, spec, cache, f"inv/{drug[:5]}")
                blob["invariance"].append(
                    _row(section, rec, drug, ALT, out))
                print(f"    [inv] {drug:12s} F={out.get('reach_headline_F')} "
                      f"ladder={(out.get('reach_ladder_F'))}", flush=True)
        save(blob)

    # ---- the invariance re-assertion ------------------------------------- #
    checks = []
    for a in blob["invariance"]:
        base = next((r for r in blob["rows"]
                     if r["section"] == a["section"] and r["drug"] == a["drug"]),
                    None)
        if base is None:
            continue
        same_F = repr(base["reach_headline_F"]) == repr(a["reach_headline_F"])
        same_ladder = (repr(base["reach_ladder_F"])
                       == repr(a["reach_ladder_F"]))
        checks.append({
            "section": a["section"], "drug": a["drug"],
            "shipped_F": base["reach_headline_F"], "alt_F": a["reach_headline_F"],
            "shipped_ladder": base["reach_ladder_F"],
            "alt_ladder": a["reach_ladder_F"],
            "headline_identical": same_F, "ladder_identical": same_ladder,
        })
    blob["invariance_checks"] = checks
    n_bad = sum(1 for c in checks if not (c["headline_identical"]
                                          and c["ladder_identical"]))
    blob["invariance_violations"] = n_bad
    blob["elapsed_s"] = time.time() - t0
    blob["alpha_end"] = assert_alpha("ladder end")
    save(blob)

    ok = sum(1 for r in blob["rows"] if r.get("status") == "ok")
    filled = sum(1 for r in blob["rows"] if r.get("reach_ladder_F"))
    print(f"\n[ladder] {ok}/{len(blob['rows'])} rows ok, "
          f"{filled} with a populated ladder, "
          f"{len(checks)} invariance pairs, {n_bad} violations, "
          f"{blob['elapsed_s'] / 60:.1f} min -> {_OUT}")
    for c in checks:
        print(f"    inv {c['section'][:22]:22s} {c['drug'][:12]:12s} "
              f"headline_same={c['headline_identical']} "
              f"ladder_same={c['ladder_identical']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
