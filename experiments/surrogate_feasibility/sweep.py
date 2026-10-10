"""Generate the surrogate feasibility study's training set, by running the SOLVER.

A surrogate's ground truth is the solver, so every row here is a real
``app.server.jobs.execute_run`` call -- the same function the Studio queue
calls. Nothing is predicted, interpolated or assumed.

DESIGN
------
Section axis
    Every section the backend can see: the 7 built-ins plus the 31 ingested
    ones (30 GBM-Space Visium + 1 Xenium breast). On 34 of the 38 the project
    REFUSES per-cell spatial statistics as readings of that tissue (the Visium
    spot refusal in ``app.server.api``), and those rows are flagged
    ``spatial_spread_valid=False``. They are kept anyway and that is deliberate:
    this study measures whether a model reproduces THE SOLVER's output, and the
    solver produces a number on those sections. No tissue claim is made from
    them. The 4 sections where the project does hold the statistic valid are
    reported separately, because that is the only subset that can answer "would
    a preview be usable on a donor".

Dose axis: DELIBERATELY ABSENT
    Measured first (``probe.py``): doubling ``pk.plasma_peak`` doubles E_p10 and
    E_p90 and leaves ``E_p90/E_p10`` and ``median_tumor_auc / plasma_peak``
    BIT-IDENTICAL. The transport chain is exactly linear in the plasma
    amplitude, so the dose axis carries no information for a relative readout
    and adding it would only inflate n with duplicate rows -- a model trained on
    it would report a smaller error for having been handed the same row twice.
    A handful of confirmation rows are run at a second amplitude
    (``--invariance-sections``) and recorded, so the claim is in the data.

Compound + PK axes
    The same configuration grid on every section, so the SECTION is the only
    thing that differs between CV folds and the design is balanced.

Writes ``experiments/surrogate_feasibility/sweep_data.json``, appending one
section at a time so a killed run keeps its work. Never touches ``results/``
or ``app/runs/``.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from typing import Dict, List, Optional

from solver_runner import (
    assert_alpha,
    make_spec,
    run_one,
    section_record,
    spatial_spread_valid,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_HERE, "sweep_data.json")

DX_UM = 30.0
BASE_PEAK = 50.0          # the shipped PKConfig default
ALT_PEAK = 200.0          # the dose-invariance confirmation amplitude

#: Six compounds spanning the transport descriptor space the solver actually
#: uses (diffusivity, BBB penetrance, clearance, binding). This is a FEATURE
#: axis, not a comparison: no ordering between compounds is produced anywhere
#: in this study (that claim was withdrawn 2026-08-21).
COMPOUNDS: List[str] = [
    "temozolomide",
    "cisplatin",
    "doxorubicin",
    "lomustine",
    "methotrexate",
    "vincristine",
]

#: (half_life_s, n_doses, total_time_s). The shipped default is the second.
#:
#: Chosen AFTER measuring which axes move anything (``probe_axes.py`` on
#: xenium_319 at dx=12). At a fixed window, ``half_life_s`` over 3600/14400/
#: 57600 s and ``n_doses`` over 1/3 leave ``E_p90/E_p10`` and the fitted
#: penetration length identical to 8 significant figures -- they are EXACTLY
#: null axes for the shape readouts, and only move the exposure MAGNITUDE.
#: ``total_time_s`` moves everything. So the window is varied, the half-life is
#: varied for the magnitude target, and ``n_doses`` is held at 1: with the
#: shipped 86400 s dose interval, doses 2 and 3 fall outside every window here
#: and changed nothing at all. Spending runs on an axis measured to be null
#: would only pad n with duplicate rows and shrink a reported error.
PK_SHAPES = [
    (14400.0, 1, 3600.0),
    (14400.0, 1, 21600.0),
    (3600.0, 1, 21600.0),
    (57600.0, 1, 21600.0),
]

#: Sections that get the second dose amplitude as an invariance check.
INVARIANCE_SECTIONS = ["xenium_319", "gbmspace_at10_bra_5_fo_1_1",
                       "visium_gbm2"]


def donor_of(section_id: str) -> str:
    """The DONOR a section came from.

    GBM-Space ships up to five sections per patient (``AT10-BRA-5-FO-1_1``,
    ``...-1_2``, ``...-2_1``, ...). Holding out one of them while training on
    another section of the SAME patient is a leak: the model has already seen
    that tissue. The CV in ``fit.py`` groups on this, not on the section id.
    """
    m = re.match(r"gbmspace_(at\d+)_", section_id)
    return m.group(1) if m else section_id


def load() -> Dict:
    if os.path.exists(_OUT):
        with open(_OUT) as fh:
            return json.load(fh)
    return {"alpha": None, "dx_um": DX_UM, "rows": [], "sections": {},
            "compounds": COMPOUNDS, "pk_shapes": PK_SHAPES,
            "base_peak_uM": BASE_PEAK, "alt_peak_uM": ALT_PEAK}


def save(blob: Dict) -> None:
    tmp = _OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(blob, fh, default=str)
    os.replace(tmp, _OUT)


def configs_for(section: str):
    """(drug, peak, half_life, n_doses, total_time) for one section."""
    for drug in COMPOUNDS:
        for hl, nd, tt in PK_SHAPES:
            yield (drug, BASE_PEAK, hl, nd, tt)
    if section in INVARIANCE_SECTIONS:
        for hl, nd, tt in PK_SHAPES[1:2]:
            for drug in COMPOUNDS[:3]:
                yield (drug, ALT_PEAK, hl, nd, tt)


def run_section(section: str, blob: Dict, quiet: bool = False) -> None:
    from app.server.cache import ReconstructionCache

    rec = section_record(section)
    valid = spatial_spread_valid(rec)
    done = {(r["drug"], r["plasma_peak"], r["half_life_s"], r["n_doses"],
             r["total_time_s"]) for r in blob["rows"] if r["section"] == section}
    cfgs = [c for c in configs_for(section) if c not in done]
    if not cfgs:
        print(f"=== {section}: already complete", flush=True)
        return

    print(f"\n=== {section}  donor={donor_of(section)}  "
          f"src={rec.get('source')}  spatial_spread_valid={valid}  "
          f"({len(cfgs)} runs)", flush=True)
    # A FRESH cache per section: the process-global one would hold an 850k-cell
    # reconstruction resident for the rest of the sweep.
    cache = ReconstructionCache(maxsize=1)
    t_sec = time.time()
    for i, (drug, peak, hl, nd, tt) in enumerate(cfgs):
        spec = make_spec(section, drug, DX_UM, tt, peak, hl, nd, 86400.0)
        out = run_one(section, spec, cache, run_id=f"{section[:14]}/{drug[:6]}")
        row = {
            "section": section,
            "donor": donor_of(section),
            "source": rec.get("source"),
            "label": rec.get("label"),
            "spatial_spread_valid": valid,
            "dx_um": DX_UM,
            "drug": drug,
            "plasma_peak": peak,
            "half_life_s": hl,
            "n_doses": nd,
            "total_time_s": tt,
        }
        row.update(out)
        row.pop("log_tail", None)
        if out["status"] != "ok":
            print(f"    !! {drug} hl={hl:g} n={nd} T={tt:g}: "
                  f"{out['error_type']}: {out['error'][:140]}", flush=True)
        elif not quiet:
            print(f"    [{i + 1:3d}/{len(cfgs)}] {drug:14s} hl={hl:>7.0f} "
                  f"n={nd} T={tt:>6.0f} peak={peak:>5.0f} "
                  f"{out['wall_clock_s']:6.2f}s "
                  f"E_p90/p10="
                  f"{(out['E_p90_uM'] / out['E_p10_uM']) if out['E_p10_uM'] else float('nan'):8.4f} "
                  f"L={out['penetration_length_um']}", flush=True)
        blob["rows"].append(row)
    blob["sections"][section] = {
        "donor": donor_of(section), "source": rec.get("source"),
        "label": rec.get("label"), "spatial_spread_valid": valid,
        "caveats": [{"capability": c.get("capability"),
                     "status": c.get("status")}
                    for c in (rec.get("caveats") or [])],
        "wall_clock_s": time.time() - t_sec,
    }
    save(blob)
    print(f"    section total {time.time() - t_sec:.1f}s", flush=True)


def main(argv: Optional[List[str]] = None) -> int:
    from app.server.cache import all_sections

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sections", default="", help="comma list; default all")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)

    blob = load()
    blob["alpha"] = assert_alpha("sweep start")
    avail = list(all_sections().keys())
    sections = ([s for s in args.sections.split(",") if s] if args.sections
                else avail)
    print(f"[sweep] {len(sections)} sections, "
          f"{len(COMPOUNDS)} compounds x {len(PK_SHAPES)} PK shapes, "
          f"dx={DX_UM:g} um", flush=True)

    t0 = time.time()
    for s in sections:
        if s not in avail:
            print(f"=== {s}: NOT in the registry, skipped", flush=True)
            continue
        try:
            run_section(s, blob, quiet=args.quiet)
        except KeyboardInterrupt:
            save(blob)
            raise
        except Exception as exc:  # noqa: BLE001
            print(f"=== {s}: SECTION FAILED {type(exc).__name__}: {exc}",
                  flush=True)
            blob["sections"][s] = {"donor": donor_of(s), "status": "FAILED",
                                   "error": f"{type(exc).__name__}: {exc}"}
            save(blob)
    blob["elapsed_s"] = time.time() - t0
    blob["alpha_end"] = assert_alpha("sweep end")
    save(blob)
    ok = sum(1 for r in blob["rows"] if r.get("status") == "ok")
    print(f"\n[sweep] {ok}/{len(blob['rows'])} rows ok in "
          f"{blob['elapsed_s'] / 60:.1f} min -> {_OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
