"""Does the vessel SOURCE conserve its painted area under ``dx``? (No. Fix it.)

THE DEFECT, AND WHY IT IS THE SINK'S DEFECT IN A DIFFERENT PLACE
----------------------------------------------------------------
``results/fstate_mechanism.md`` decomposed the ``F_state`` grid drift with an
exact additive identity and attributed it as ``cov(term, F_state)/var(F_state)``:

    vessel-source rasterisation (write side)   +1.057
    PDE discretisation                         -0.056
    cell-to-grid sampling (read side)          -0.001

i.e. essentially 100 % of it is the binary vessel painting.
``vessel_network.discrete_vascular_field(geometry='binary')`` marks a grid cell
as vessel iff its CENTRE falls within ``VESSEL_RADIUS_UM = 7.5 um`` of an
endothelial cell, and then deposits a whole grid cell of vessel -- area ``dx^2``.
A 7.5 um disc is 176.71 um^2 against a ``dx = 30 um`` cell's 900 um^2, so a
vessel quantises to 0 or 1 cells and the TOTAL PAINTED VESSEL AREA IS NOT
CONSERVED under ``dx``.

``physics/oxygen.py`` already documents, and already fixes, exactly this defect
on the other side of the same PDE. ``_tumor_density`` refused to feed a raw
per-grid-cell COUNT into the sink -- a count carries a ``dx^2`` factor -- and
divides by the grid-cell area, re-expressing the result against a fixed physical
reference tile. That fix works: ``fstate_mechanism.md`` sec. 3 measures the
integrated tumour sink conserved to 2.9e-11 over the same ladder. THE SOURCE
NEVER RECEIVED THE EQUIVALENT TREATMENT.

WHAT THIS SCRIPT DOES
---------------------
``--stage area``      Reproduce the diagnosis independently (painted vessel area
                      across the dx ladder, binary), then measure the
                      area-conserving painting against it, plus the sub-pixel
                      pitch convergence and the sub-grid painting for context.
``--stage fstate``    The ``F_state`` ladder under three source arms, so "how
                      much of the 0.0461 does the fix remove?" is a measurement.
``--stage verdicts``  The nine published 3-D:2-D rows under the fixed source,
                      through ``experiments.delivery_source_geometry``'s existing
                      blast harness (the binary control arm is replayed from that
                      report's own checkpoint, so the control is not re-derived
                      here and cannot drift).
``--stage alpha``     ``calibrate_alpha()`` must return 148.60815765445545.

WHAT IS CHANGED UNDER ``tmesim/``
---------------------------------
One thing: ``reconstruct/vessel_network.py`` gains an OPT-IN
``geometry='area'`` branch (:func:`~reconstruct.vessel_network._area_coverage_field`)
and the constant ``SOURCE_SUBPIXEL_PITCH_UM``. The ``'binary'`` default is
untouched and bit-identical. No band edge, no ``OxygenConfig.supply``, no
``DX_UM``, no alpha, no ``_DEFAULT_VASCULAR_MODE`` moves.

Usage
-----
    ./.venv/bin/python -u experiments/vessel_source_conservation.py --stage area
    ./.venv/bin/python -u experiments/vessel_source_conservation.py --stage fstate
    ./.venv/bin/python -u experiments/vessel_source_conservation.py --stage alpha
    ./.venv/bin/python -u -m experiments.vessel_source_conservation --stage verdicts
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.core.config import OxygenConfig, PDConfig                 # noqa: E402
from tmesim.physics.oxygen import _tumor_density                      # noqa: E402
from tmesim.reconstruct.real_pipeline import reconstruct_real         # noqa: E402
from tmesim.reconstruct.vessel_network import (                       # noqa: E402
    SOURCE_SUBPIXEL_PITCH_UM,
    VESSEL_RADIUS_UM,
    _area_coverage_field,
    discrete_vascular_field,
    segment_vessels,
)

OUT_DIR = os.path.join(_ROOT, "results", "sweep")

#: The same ladder ``experiments/fstate_mechanism.py`` uses, so the numbers here
#: are directly comparable to ``results/fstate_mechanism.md``.
DX_LADDER = (30.0, 24.0, 18.0, 16.0, 12.0, 10.0, 8.0, 7.0, 6.0, 5.0)

#: The four rungs the brief names explicitly.
DX_BRIEF = (30.0, 18.0, 12.0, 8.0)

#: Production spacing (``results/numerical_rigor.md``). NOT changed anywhere.
PRODUCTION_DX_UM = 12.0

#: Sub-pixel pitches swept to show the estimator's own convergence.
PITCHES = (4.0, 2.0, 1.0, 0.5, 0.25)

#: Rungs that are NOT integer multiples of :data:`SOURCE_SUBPIXEL_PITCH_UM`, so
#: the sub-lattice genuinely differs between them and the conservation guarantee
#: weakens from "identical" to "a Riemann error of order perimeter x pitch". The
#: claim that conservation does not DEPEND on the convenient ladder is worth an
#: experiment, not a footnote.
DX_OFF_LATTICE = (45.0, 22.5, 11.25, 7.5, 3.75)

SAMPLE_GBM = os.path.join(_ROOT, "data_cache", "xenium_gbm")


def log(msg: str) -> None:
    print(f"[vsrc {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _dump(obj, name: str) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=float)
    log(f"written -> results/sweep/{name}")


def _span_pct(a: Sequence[float]) -> float:
    x = np.asarray(list(a), dtype=np.float64)
    return float(100.0 * (x.max() - x.min()) / x.mean()) if x.size else float("nan")


# --------------------------------------------------------------------------- #
# The three source arms
# --------------------------------------------------------------------------- #
#: ``binary``       -- production, unchanged. Whole grid cell per centre-in-disc.
#: ``area``         -- the fix. SAME discs, SAME radius, SAME endothelial cells;
#:                    the grid cell carries the covered-area FRACTION.
#: ``area_matched`` -- the fix with ``match_total_to`` set to the production
#:                    field, i.e. the tissue-integrated total held EXACTLY at the
#:                    shipped value and only the distribution allowed to move.
#:                    This is the arm that separates "the total changed" from
#:                    "the distribution changed", which is the whole question of
#:                    whether this is a fix or a recalibration.
ARMS = ("binary", "area", "area_matched")


def source_field(domain, vessels: Dict, arm: str) -> np.ndarray:
    base = discrete_vascular_field(domain, vessels)
    if arm == "binary":
        return base
    if arm == "area":
        return discrete_vascular_field(domain, vessels, geometry="area")
    if arm == "area_matched":
        return discrete_vascular_field(domain, vessels, geometry="area",
                                       match_total_to=base)
    if arm == "subgrid":
        return discrete_vascular_field(domain, vessels, geometry="subgrid")
    raise ValueError(f"unknown arm {arm!r}")


# --------------------------------------------------------------------------- #
# Stage 1 -- painted area conservation
# --------------------------------------------------------------------------- #
def stage_area(dx_ladder: Sequence[float] = DX_LADDER,
               pitches: Sequence[float] = PITCHES) -> Dict:
    """Reproduce the 4.73 % binary span, then measure the fix against it.

    Nothing here solves a PDE: this is purely "how much vessel did the painting
    deposit, and does that depend on the mesh?".
    """
    rows: List[Dict] = []
    pitch_rows: List[Dict] = []
    for dx in dx_ladder:
        dx = float(dx)
        t0 = time.time()
        micro = reconstruct_real(source="xenium", dx_um=dx)
        dom, cells = micro.domain, micro.cells
        v = segment_vessels(cells, dom)
        tissue = np.asarray(dom.tissue_mask, dtype=bool)
        f_bin = source_field(dom, v, "binary")
        f_area = source_field(dom, v, "area")
        f_sub = source_field(dom, v, "subgrid")
        # The painting BEFORE the tissue mask, so the mask's own (separately
        # owned, separately documented) dx-dependence cannot be mistaken for the
        # painting's.
        f_area_raw = _area_coverage_field(dom, v["xy"], VESSEL_RADIUS_UM)
        cell_area = dx * dx
        rows.append({
            "dx_um": dx,
            "grid": list(dom.shape),
            "n_endothelial_cells": int(v["n_endothelial_cells"]),
            "n_vessels": int(v["n_vessels"]),
            "cell_area_um2": cell_area,
            "disc_over_cell_area": float(np.pi * VESSEL_RADIUS_UM ** 2 / cell_area),
            "binary_area_um2": float(f_bin.sum() * cell_area),
            "area_area_um2": float(f_area.sum() * cell_area),
            "area_unmasked_um2": float(f_area_raw.sum() * cell_area),
            "subgrid_area_um2": float(f_sub.sum() * cell_area),
            "binary_support_frac": float((f_bin > 0)[tissue].mean()),
            "area_support_frac": float((f_area > 0)[tissue].mean()),
            "area_saturated_frac": float((f_area >= 1.0)[tissue].mean()),
            "tissue_frac": float(tissue.mean()),
            # the shipped reconstruction must BE the binary arm, or nothing below
            # is about production
            "shipped_is_binary": bool(np.array_equal(
                np.clip(np.asarray(micro.fields.vascular_density, float), 0, 1),
                np.clip(f_bin, 0, 1))),
            "secs": time.time() - t0,
        })
        log(f"dx={dx:g}: binary {rows[-1]['binary_area_um2']:.0f} um2 | "
            f"area {rows[-1]['area_area_um2']:.0f} um2 | "
            f"subgrid {rows[-1]['subgrid_area_um2']:.0f} um2 "
            f"({rows[-1]['secs']:.0f}s)")
        _dump({"rows": rows}, "vsrc_area_ladder.json")

    # Sub-pixel pitch convergence, at the production spacing and at the coarsest
    # rung: the estimator's OWN error, which is what sets the absolute value.
    for dx in (PRODUCTION_DX_UM, float(max(dx_ladder))):
        micro = reconstruct_real(source="xenium", dx_um=dx)
        v = segment_vessels(micro.cells, micro.domain)
        for p in pitches:
            f = _area_coverage_field(micro.domain, v["xy"], VESSEL_RADIUS_UM, float(p))
            pitch_rows.append({"dx_um": dx, "pitch_um": float(p),
                               "area_um2": float(f.sum() * dx * dx)})
            log(f"pitch dx={dx:g} h={p:g}: {pitch_rows[-1]['area_um2']:.1f} um2")

    # Off-lattice rungs: does conservation depend on dx dividing the pitch?
    off_rows: List[Dict] = []
    for dx in DX_OFF_LATTICE:
        dx = float(dx)
        micro = reconstruct_real(source="xenium", dx_um=dx)
        v = segment_vessels(micro.cells, micro.domain)
        f_b = source_field(micro.domain, v, "binary")
        f_a = source_field(micro.domain, v, "area")
        off_rows.append({"dx_um": dx,
                         "binary_area_um2": float(f_b.sum() * dx * dx),
                         "area_area_um2": float(f_a.sum() * dx * dx)})
        log(f"off-lattice dx={dx:g}: binary "
            f"{off_rows[-1]['binary_area_um2']:.0f} | area "
            f"{off_rows[-1]['area_area_um2']:.0f} um2")

    binary = [r["binary_area_um2"] for r in rows]
    area = [r["area_area_um2"] for r in rows]
    sub = [r["subgrid_area_um2"] for r in rows]
    brief = {d: r for d, r in zip([r["dx_um"] for r in rows], rows)}
    prod = brief[PRODUCTION_DX_UM]
    out = {
        "rows": rows, "pitch_rows": pitch_rows, "off_lattice_rows": off_rows,
        "off_lattice_binary_span_pct": _span_pct(
            [r["binary_area_um2"] for r in off_rows]),
        "off_lattice_area_span_pct": _span_pct(
            [r["area_area_um2"] for r in off_rows]),
        "vessel_radius_um": float(VESSEL_RADIUS_UM),
        "subpixel_pitch_um": float(SOURCE_SUBPIXEL_PITCH_UM),
        "binary_span_pct": _span_pct(binary),
        "area_span_pct": _span_pct(area),
        "subgrid_span_pct": _span_pct(sub),
        "binary_span_pct_brief_rungs": _span_pct(
            [brief[d]["binary_area_um2"] for d in DX_BRIEF]),
        "area_span_pct_brief_rungs": _span_pct(
            [brief[d]["area_area_um2"] for d in DX_BRIEF]),
        # THE question the brief asks: does the fix move the SHIPPED total?
        "production_binary_area_um2": prod["binary_area_um2"],
        "production_area_area_um2": prod["area_area_um2"],
        "production_total_change_pct": float(
            100.0 * (prod["area_area_um2"] / prod["binary_area_um2"] - 1.0)),
        "binary_ladder_mean_um2": float(np.mean(binary)),
        "area_vs_binary_ladder_mean_pct": float(
            100.0 * (np.mean(area) / np.mean(binary) - 1.0)),
        "all_shipped_is_binary": all(r["shipped_is_binary"] for r in rows),
    }
    log(f"SPAN binary {out['binary_span_pct']:.4f} % | area "
        f"{out['area_span_pct']:.4f} % | subgrid {out['subgrid_span_pct']:.4f} %")
    log(f"total change at production dx={PRODUCTION_DX_UM:g}: "
        f"{out['production_total_change_pct']:+.3f} %")
    _dump(out, "vsrc_area_ladder.json")
    return out


# --------------------------------------------------------------------------- #
# Stage 2 -- what it does to F_state
# --------------------------------------------------------------------------- #
def stage_fstate(dx_ladder: Sequence[float] = DX_LADDER,
                 arms: Sequence[str] = ARMS) -> Dict:
    """The ``F_state`` ladder under each source arm.

    Reuses ``experiments.fstate_mechanism``'s solver / sampler / banding helpers
    verbatim, so this ladder is directly comparable to ``results/
    fstate_mechanism.md`` sec. 3 rather than merely similar to it. The ``binary``
    arm MUST reproduce that report's published ladder (0.5041-0.5502); that is
    asserted, not hoped for.
    """
    from experiments.fstate_mechanism import f_state, sample_floor, solve_from_fields

    cfg, pd_cfg = OxygenConfig(), PDConfig()
    rows: List[Dict] = []
    for dx in dx_ladder:
        dx = float(dx)
        t0 = time.time()
        micro = reconstruct_real(source="xenium", dx_um=dx)
        dom = micro.domain
        v = segment_vessels(micro.cells, dom)
        rho = _tumor_density(micro)
        tum = micro.cells.tumor_mask()
        xt, yt = micro.cells.x[tum], micro.cells.y[tum]
        tissue = np.asarray(dom.tissue_mask, dtype=bool)
        row: Dict = {"dx_um": dx, "grid": list(dom.shape),
                     "n_tumour": int(tum.sum())}
        for arm in arms:
            field = source_field(dom, v, arm)
            o = solve_from_fields(field, rho, dx, cfg)
            fs = f_state(sample_floor(o, dx, xt, yt), cfg, pd_cfg)
            row[arm] = {
                **fs,
                "source_area_um2": float(field.sum() * dx * dx),
                "source_sum_over_tissue": float(field[tissue].sum()),
                "oxygen_mean_tissue": float(o[tissue].mean()),
            }
        row["secs"] = time.time() - t0
        rows.append(row)
        log("dx={:g}: ".format(dx) + " | ".join(
            f"{a} F={row[a]['f_state']:.4f} q={row[a]['quiescent_frac']:.4f} "
            f"O2={row[a]['mean_o2']:.4f}" for a in arms)
            + f" ({row['secs']:.0f}s)")
        _dump({"rows": rows}, "vsrc_fstate_ladder.json")

    out: Dict = {"rows": rows, "arms": list(arms), "summary": {}}
    for arm in arms:
        fs = [r[arm]["f_state"] for r in rows]
        o2 = [r[arm]["mean_o2"] for r in rows]
        q = [r[arm]["quiescent_frac"] for r in rows]
        out["summary"][arm] = {
            "f_state_min": float(np.min(fs)), "f_state_max": float(np.max(fs)),
            "f_state_span": float(np.max(fs) - np.min(fs)),
            "mean_o2_span": float(np.max(o2) - np.min(o2)),
            "quiescent_span": float(np.max(q) - np.min(q)),
            "f_state_at_production": float(
                [r[arm]["f_state"] for r in rows
                 if r["dx_um"] == PRODUCTION_DX_UM][0]),
            "source_area_span_pct": _span_pct(
                [r[arm]["source_area_um2"] for r in rows]),
        }
        log(f"{arm}: F_state span {out['summary'][arm]['f_state_span']:.4f} "
            f"({out['summary'][arm]['f_state_min']:.4f}-"
            f"{out['summary'][arm]['f_state_max']:.4f}), "
            f"production {out['summary'][arm]['f_state_at_production']:.4f}")
    base = out["summary"]["binary"]["f_state_span"]
    for arm in arms:
        s = out["summary"][arm]["f_state_span"]
        out["summary"][arm]["drift_removed_frac"] = float(1.0 - s / base) if base else 0.0
    _dump(out, "vsrc_fstate_ladder.json")
    return out


# --------------------------------------------------------------------------- #
# Stage 3 -- the nine published verdicts
# --------------------------------------------------------------------------- #
def stage_verdicts(dx_um: float = PRODUCTION_DX_UM,
                   arms: Sequence[str] = ("area", "area_matched"),
                   tissues: Sequence[str] = ("gbm", "breast")) -> Dict:
    """Re-score the nine published rows under the fixed source.

    Runs through ``experiments.delivery_source_geometry``'s existing blast
    harness -- the same ``_run_panel`` / ``score_ratio_row`` path
    ``results/delivery_source_geometry.md`` used -- so the arms here are directly
    comparable to that report's ``binary`` / ``subgrid`` arms. The ``binary``
    control is REPLAYED from that report's own checkpoint rather than re-solved,
    which is both cheaper and stricter: the control cannot drift because it is
    literally the published one.

    ``dsg._source_field`` is patched IN THIS PROCESS ONLY to teach it the two new
    arms. ``experiments/delivery_source_geometry.py`` is NOT edited (it is owned
    elsewhere), and the patch delegates every arm it does not own straight back
    to the original function.
    """
    import experiments.delivery_source_geometry as dsg
    from tmesim.reconstruct.vessel_network import LINK_DISTANCE_UM

    original = dsg._source_field

    def patched(micro, arm: str) -> np.ndarray:
        if arm in ("area", "area_matched"):
            v = dsg._segment_safe(micro.cells, micro.domain, float(LINK_DISTANCE_UM))
            return source_field(micro.domain, v, arm)
        return original(micro, arm)

    dsg._source_field = patched
    ck = os.path.join(OUT_DIR, f"vsrc_verdicts_dx{dx_um:g}.jsonl")
    published_ck = os.path.join(
        OUT_DIR, f"delivery_source_geometry_blast_dx{dx_um:g}.jsonl")

    panels: List[Dict] = []
    # the published binary control, replayed
    control: Dict[str, Dict] = {}
    if os.path.exists(published_ck):
        with open(published_ck) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                p = json.loads(line)
                if p.get("arm") == "binary" and float(p.get("dx_um", -1)) == float(dx_um):
                    control[p["tissue"]] = p
        log(f"replayed binary control panels: {sorted(control)}")

    done: Dict[str, Dict] = {}
    if os.path.exists(ck):
        with open(ck) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                p = json.loads(line)
                if float(p.get("dx_um", -1)) == float(dx_um):
                    done[f"{p['tissue']}|{p['arm']}"] = p
        log(f"checkpoint: {len(done)} panel(s) already complete")

    from tmesim.validation.calibration import calibrate_alpha
    alpha = float(calibrate_alpha())
    log(f"alpha = {alpha!r}")

    handle = open(ck, "a")
    try:
        for tissue in tissues:
            if tissue in control:
                panels.append(control[tissue])
            drugs = (dsg.BREAST_DRUGS if tissue == "breast" else dsg.GBM_DRUGS)
            for arm in arms:
                key = f"{tissue}|{arm}"
                if key in done:
                    log(f"{key}: cached")
                    panels.append(done[key])
                    continue
                log(f"=== panel {key} (alpha={alpha:.6f}) ===")
                p = dsg._run_panel(tissue, arm, alpha, drugs, float(dx_um))
                panels.append(p)
                handle.write(json.dumps(p, default=float) + "\n")
                handle.flush()
    finally:
        handle.close()
        dsg._source_field = original

    for p in panels:
        n = dsg.rescore_panel(p)
        log(f"{p['tissue']}|{p['arm']}: {n} verdict(s) re-derived")

    # assemble the nine-row comparison
    by = {(p["tissue"], p["arm"]): p for p in panels}
    table: List[Dict] = []
    for (tissue, drug), pub_ratio in dsg.PUBLISHED_RATIO.items():
        rec: Dict = {"tissue": tissue, "drug": drug,
                     "published_ratio": float(pub_ratio),
                     "published_verdict": dsg.PUBLISHED_VERDICT[(tissue, drug)]}
        for arm in ("binary",) + tuple(arms):
            p = by.get((tissue, arm))
            if p is None:
                continue
            r = next((x for x in p["rows"] if x["drug"] == drug), None)
            if r is None:
                continue
            rec[arm] = {"ratio": float(r["ratio_exact"]),
                        "verdict": r["verdict"],
                        "verdict_scored": bool(r["verdict_scored"]),
                        "f_penetration": float(r["f_penetration"]),
                        "f_state": float(r["f_state"]),
                        "e_unit_mean": float(r["e_unit_mean"])}
        table.append(rec)

    flips = {}
    for arm in arms:
        n_flip = sum(1 for rec in table
                     if arm in rec and "binary" in rec
                     and rec[arm]["verdict_scored"]
                     and rec[arm]["verdict"] != rec["binary"]["verdict"])
        flips[arm] = int(n_flip)
    control_reproduces = all(
        abs(rec["binary"]["ratio"] - rec["published_ratio"]) < 5e-5
        for rec in table if "binary" in rec)

    out = {"dx_um": float(dx_um), "alpha": alpha, "arms": list(arms),
           "panels": panels, "table": table, "scored_verdict_flips": flips,
           "binary_control_reproduces_published": bool(control_reproduces)}
    for rec in table:
        line = f"{rec['tissue']}/{rec['drug']}: pub {rec['published_ratio']:.4f} " \
               f"[{rec['published_verdict']}]"
        for arm in ("binary",) + tuple(arms):
            if arm in rec:
                line += f" | {arm} {rec[arm]['ratio']:.4f} {rec[arm]['verdict']}"
        log(line)
    log(f"scored verdict flips: {flips}; binary control reproduces published: "
        f"{control_reproduces}")
    _dump(out, "vsrc_verdicts.json")
    return out


# --------------------------------------------------------------------------- #
# Stage 4 -- hygiene
# --------------------------------------------------------------------------- #
#: The value ``calibrate_alpha()`` must return, before and after this work.
ALPHA_EXPECTED = 148.60815765445545


def stage_alpha() -> Dict:
    """Every constant this work is forbidden to move, checked rather than claimed."""
    from tmesim.core.config import OxygenConfig as _Ox
    from tmesim.reconstruct import real_fields
    from tmesim.validation.calibration import calibrate_alpha

    cfg = _Ox()
    out = {
        "alpha": float(calibrate_alpha()),
        "alpha_expected": ALPHA_EXPECTED,
        "alpha_unchanged": bool(float(calibrate_alpha()) == ALPHA_EXPECTED),
        "oxygen_supply": float(cfg.supply),
        "oxygen_consumption": float(cfg.consumption),
        "hypoxia_threshold": float(cfg.hypoxia_threshold),
        "quiescent_ceiling": 2.0 * float(cfg.hypoxia_threshold),
        "vessel_radius_um": float(VESSEL_RADIUS_UM),
        "default_vascular_mode": str(real_fields._DEFAULT_VASCULAR_MODE),
        "discrete_default_geometry": _default_geometry(),
    }
    log(f"alpha={out['alpha']!r} unchanged={out['alpha_unchanged']}; "
        f"supply={out['oxygen_supply']}, mode={out['default_vascular_mode']}, "
        f"geometry default={out['discrete_default_geometry']}")
    _dump(out, "vsrc_hygiene.json")
    return out


def _default_geometry() -> str:
    import inspect
    sig = inspect.signature(discrete_vascular_field)
    return str(sig.parameters["geometry"].default)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="area",
                    choices=("area", "fstate", "verdicts", "alpha"))
    ap.add_argument("--dx", type=float, default=PRODUCTION_DX_UM)
    ap.add_argument("--arms", default="")
    ap.add_argument("--tissues", default="gbm,breast")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    if args.stage == "area":
        stage_area()
    elif args.stage == "fstate":
        arms = (tuple(a.strip() for a in args.arms.split(",") if a.strip())
                if args.arms else ARMS)
        stage_fstate(arms=arms)
    elif args.stage == "verdicts":
        arms = (tuple(a.strip() for a in args.arms.split(",") if a.strip())
                if args.arms else ("area", "area_matched"))
        stage_verdicts(dx_um=args.dx, arms=arms,
                       tissues=tuple(t.strip() for t in args.tissues.split(",")
                                     if t.strip()))
    else:
        stage_alpha()


if __name__ == "__main__":
    main()
