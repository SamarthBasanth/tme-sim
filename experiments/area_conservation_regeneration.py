"""Regenerate the published results under ``geometry='area'``, and re-run the
sub-grid closure on the field the promotion would actually ship.

WHY THIS EXISTS
---------------
``results/vessel_source_conservation.md`` §10 lands ``geometry='area'`` opt-in
and default-off, and names three things a promotion to the default would need.
Two of them are measurements and are made here:

1. **The published ratios REGENERATED under ``geometry='area'``**, not merely
   re-scored, so the -0.794 % source change is separately visible rather than
   absorbed. Both arms are re-solved in the SAME process off the SAME
   reconstruction, so ``binary`` is a live control that must reproduce
   ``delivery_source_geometry.PUBLISHED_RATIO`` rather than a replayed cache.
   The per-tumour-cell exposure field is captured on every solve, so
   ``results/p99_convention.md``'s normaliser family is scored on those same
   solves without a second set.

2. **``results/subgrid_closure.md`` re-run on the continuous field.** That
   report's §5.0 control -- "on the shipped binary vascular field the closure is
   an exact no-op, because ``beta == 1`` identically when the field is exactly 0
   or 1" -- is a property of the BINARY painting and is VOID the moment the
   default becomes continuous. §5 itself ran on ``geometry='subgrid'``, which is
   a DIFFERENT continuous field (one footprint per vessel OBJECT at a
   data-derived radius, ~347k um^2 against ~744k). Nobody has measured the
   closure on ``geometry='area'``, which is the field the promotion would ship.

The third thing §10 names -- leaving ``OxygenConfig.supply`` alone -- is a
prohibition, not a measurement, and is checked in :func:`stage_alpha` alongside
every other constant this work is forbidden to move.

NOTHING IN ``tmesim/`` IS TOUCHED. ``geometry='binary'`` stays the default;
``experiments/delivery_source_geometry.py`` is not edited (it is owned
elsewhere) and is instead patched IN THIS PROCESS ONLY, exactly as
``experiments/vessel_source_conservation.py`` already does.

Run::

    ./.venv/bin/python -u -m experiments.area_conservation_regeneration --stage alpha
    ./.venv/bin/python -u -m experiments.area_conservation_regeneration --stage closure
    ./.venv/bin/python -u -m experiments.area_conservation_regeneration --stage ratios
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

SWEEP = os.path.join(_ROOT, "results", "sweep")
_T0 = time.time()

#: The value ``calibrate_alpha()`` must return, before and after this work.
ALPHA_EXPECTED = 148.60815765445545

#: The published GBM ladder (``results/grid_convergence_gbm.md`` §5), which is
#: also the ladder ``results/subgrid_closure.md`` §5.1 spans, plus its fine probe.
DX_LADDER = (30.0, 18.0, 12.0, 8.0)
DX_FINE = (6.0, 5.0)
PRODUCTION_DX_UM = 12.0


def log(msg: str) -> None:
    print(f"[areagen +{time.time() - _T0:7.1f}s] {msg}", flush=True)


def _dump(obj, name: str) -> str:
    os.makedirs(SWEEP, exist_ok=True)
    path = os.path.join(SWEEP, f"areagen_{name}")
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=float)
    log(f"wrote {os.path.relpath(path, _ROOT)}")
    return path


# --------------------------------------------------------------------------- #
# Stage 0 -- the constants this work is forbidden to move
# --------------------------------------------------------------------------- #
def stage_alpha() -> Dict:
    """Every constant this work is forbidden to move, CHECKED rather than claimed."""
    from tmesim.core.config import OxygenConfig
    from tmesim.reconstruct import real_fields
    from tmesim.reconstruct.vessel_network import (
        LINK_DISTANCE_UM, VESSEL_RADIUS_UM, discrete_vascular_field)
    from tmesim.validation.calibration import calibrate_alpha
    from tmesim.validation.heldout_patient import DX_UM as _HELDOUT_DX_UM
    from tmesim.validation.mgmt_scalar import DX_UM as _MGMT_DX_UM

    ox = OxygenConfig()
    import inspect
    sig = inspect.signature(discrete_vascular_field)
    out = {
        "alpha": float(calibrate_alpha()),
        "alpha_expected": ALPHA_EXPECTED,
        "oxygen_supply": float(ox.supply),
        "oxygen_consumption": float(ox.consumption),
        "hypoxia_threshold": float(ox.hypoxia_threshold),
        "default_vascular_mode": real_fields._DEFAULT_VASCULAR_MODE,
        "default_geometry": str(sig.parameters["geometry"].default),
        "vessel_radius_um": float(VESSEL_RADIUS_UM),
        "link_distance_um": float(LINK_DISTANCE_UM),
        "dx_um_heldout_patient": float(_HELDOUT_DX_UM),
        "dx_um_mgmt_scalar": float(_MGMT_DX_UM),
    }
    out["dx_um_ok"] = bool(_HELDOUT_DX_UM == 12.0 and _MGMT_DX_UM == 12.0)
    out["alpha_ok"] = bool(out["alpha"] == ALPHA_EXPECTED)
    out["supply_ok"] = bool(out["oxygen_supply"] == 0.05)
    out["default_geometry_ok"] = bool(out["default_geometry"] == "binary")
    out["default_mode_ok"] = bool(out["default_vascular_mode"] == "discrete")
    out["all_ok"] = bool(out["alpha_ok"] and out["supply_ok"]
                         and out["default_geometry_ok"] and out["default_mode_ok"]
                         and out["dx_um_ok"])
    for k, v in out.items():
        log(f"  {k} = {v!r}")
    _dump(out, "alpha.json")
    return out


# --------------------------------------------------------------------------- #
# Stage 1 -- the sub-grid closure on the field the promotion would ship
# --------------------------------------------------------------------------- #
def _beta_diagnostics(vasc: np.ndarray, dx: float, supply: float,
                      consumption_rho: np.ndarray) -> Dict:
    """The §5.2 table of ``results/subgrid_closure.md``, on THIS field."""
    from experiments.subgrid_closure import BETA_VALIDATED_MAX
    from tmesim.physics.subgrid import source_reweight

    beta = source_reweight(vasc, float(dx), float(supply), consumption_rho,
                           2000.0, model="open")
    partial = (vasc > 0) & (vasc < 1)
    v = vasc[partial]
    b = beta[partial]
    if v.size == 0:
        return {"n_partial_cells": 0}
    src = vasc * float(supply)
    tot = float(src.sum())
    return {
        "n_partial_cells": int(v.size),
        "frac_source_cells_partial": float(partial.sum() / max((vasc > 0).sum(), 1)),
        "v_median": float(np.median(v)),
        "v_p01": float(np.percentile(v, 1.0)),
        "v_frac_below_0p1": float((v < 0.1).mean()),
        "beta_median": float(np.median(b)),
        "beta_p75": float(np.percentile(b, 75.0)),
        "beta_p90": float(np.percentile(b, 90.0)),
        "beta_max": float(b.max()),
        "beta_validated_max": float(BETA_VALIDATED_MAX),
        "frac_partial_beyond_validated": float((b > BETA_VALIDATED_MAX).mean()),
        "source_mass_frac_beyond_validated": float(
            src[partial][b > BETA_VALIDATED_MAX].sum() / tot) if tot > 0 else 0.0,
        "median_rel_dev_from_inv_v": float(
            np.median(np.abs(b - 1.0 / np.maximum(v, 1e-300)) / b)),
    }


def stage_closure(rungs: Sequence[float] = DX_LADDER + DX_FINE,
                  paintings: Sequence[str] = ("binary", "area"),
                  beta_cap: bool = True) -> Dict:
    """``results/subgrid_closure.md`` §5, re-run on ``geometry='area'``.

    Every arm is a fresh steady oxygen solve; the closure and the band decision
    both happen after the solve, but ``beta`` enters the reaction diagonal, so
    the closed arms CANNOT be re-scored off a cached field and are re-solved.
    The ``binary`` arms ARE the published control and are re-solved too, so the
    reproduction of ``grid_convergence_gbm.md`` §5 is live rather than replayed.
    """
    from experiments.subgrid_closure import (
        BETA_VALIDATED_MAX, JOINT_CONSUMPTION_HI, JOINT_CONSUMPTION_LO,
        JOINT_SUPPLY, SHIPPED_CONSUMPTION, SHIPPED_SUPPLY, VESSEL_R, _tissue_arm)
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    arms = (("shipped", SHIPPED_SUPPLY, SHIPPED_CONSUMPTION),
            ("joint_lo", JOINT_SUPPLY, JOINT_CONSUMPTION_LO),
            ("joint_hi", JOINT_SUPPLY, JOINT_CONSUMPTION_HI))

    ck = os.path.join(SWEEP, "areagen_closure_partial.json")
    out: Dict = {"section": "data_cache/xenium_gbm", "dx_ladder": list(rungs),
                 "published_ladder": list(DX_LADDER),
                 "arms": {}, "beta": {}, "per_cell": {}}
    done = set()
    if os.path.exists(ck):
        with open(ck) as fh:
            prev = json.load(fh)
        out["arms"] = prev.get("arms", {})
        out["beta"] = prev.get("beta", {})
        for key, rows in out["arms"].items():
            for r in rows:
                done.add((key, float(r["dx_um"])))
        log(f"checkpoint: {len(done)} arm/rung measurement(s) already complete")

    for dx in rungs:
        t0 = time.time()
        micro = None
        for tag, sup, con in arms:
            for painting in paintings:
                for closure in (False, True):
                    if painting == "binary" and closure and dx != PRODUCTION_DX_UM:
                        # beta == 1 identically on a binary field. Run it ONCE at
                        # the published grid to reproduce §5.0's no-op control.
                        continue
                    key = f"{tag}|{painting}|{'closed' if closure else 'raw'}"
                    if (key, float(dx)) in done:
                        continue
                    if micro is None:
                        micro = reconstruct_real(source="xenium", dx_um=float(dx))
                        log(f"dx={dx:g}: reconstructed {micro.domain.shape} "
                            f"in {time.time() - t0:.1f}s")
                    r = _tissue_arm(micro, sup, con, painting=painting,
                                    closure=closure)
                    r.pop("o2_per_tumour_cell", None)
                    r["dx_um"] = float(dx)
                    out["arms"].setdefault(key, []).append(r)
                    log(f"  {key:26s} dx={dx:5g}  F_state={r['f_state']:.4f}  "
                        f"meanO2={r['mean_o2_tumour']:.4f}  "
                        f"src={r['source_integral']:.1f} um^2  "
                        f"binary_field={r['vasc_is_binary']}")

        # the beta table, on each painting, at this rung
        for painting in paintings:
            bkey = f"{painting}|{dx:g}"
            if bkey in out["beta"]:
                continue
            if micro is None:
                micro = reconstruct_real(source="xenium", dx_um=float(dx))
            import dataclasses
            from tmesim.core import Microenvironment
            from tmesim.physics.oxygen import _tumor_density
            from tmesim.reconstruct.real_fields import vascular_density_from_cells
            if painting == "binary":
                v = np.asarray(micro.fields.vascular_density, dtype=np.float64)
            else:
                v = vascular_density_from_cells(
                    micro.domain, micro.cells, mode="discrete",
                    geometry=painting, radius_um=VESSEL_R,
                    **({"radius_mode": "constant"} if painting == "subgrid" else {}))
            work = Microenvironment(
                domain=micro.domain,
                fields=dataclasses.replace(micro.fields, vascular_density=v),
                cells=micro.cells.copy(), patient_id=micro.patient_id,
                metadata=micro.metadata)
            rho = _tumor_density(work)
            for tag, sup, con in arms:
                out["beta"][f"{painting}|{dx:g}|{tag}"] = _beta_diagnostics(
                    v, float(micro.domain.dx), float(sup), float(con) * rho)

        with open(ck, "w") as fh:
            json.dump({k: v for k, v in out.items() if k != "per_cell"}, fh,
                      indent=2, default=float)
        log(f"dx={dx:g}: rung complete ({time.time() - t0:.1f}s)")
        micro = None

    # --- the capped-beta diagnosis arm, at the published ladder only ---------
    if beta_cap:
        for dx in DX_LADDER:
            need = [(t, s, c) for t, s, c in arms
                    if (f"{t}|area|capped", float(dx)) not in done]
            if not need:
                continue
            micro = reconstruct_real(source="xenium", dx_um=float(dx))
            for tag, sup, con in need:
                r = _tissue_arm(micro, sup, con, painting="area", closure=True,
                                beta_cap=BETA_VALIDATED_MAX)
                r.pop("o2_per_tumour_cell", None)
                r["dx_um"] = float(dx)
                out["arms"].setdefault(f"{tag}|area|capped", []).append(r)
                log(f"  {tag}|area|capped dx={dx:5g}  F_state={r['f_state']:.4f}")
            with open(ck, "w") as fh:
                json.dump({k: v for k, v in out.items() if k != "per_cell"},
                          fh, indent=2, default=float)

    out.pop("per_cell", None)
    out["summary"] = _closure_summary(out["arms"])
    for key in sorted(out["summary"]):
        s = out["summary"][key]
        span = s["f_state_span_ladder"]
        log(f"SUMMARY {key:28s} span(ladder)="
            f"{'n/a' if span is None else format(span, '.4f')} "
            f"dx12={s['f_state_at_dx12']} finest={s['f_state_finest']}")
    _dump(out, "closure.json")
    return out


def _closure_summary(arms: Dict) -> Dict:
    summary = {}
    for key, rows in arms.items():
        rows = sorted(rows, key=lambda r: -r["dx_um"])
        lad = [r for r in rows if r["dx_um"] in DX_LADDER]
        if len(lad) < 2:
            summary[key] = {
                "f_state_by_dx": {str(r["dx_um"]): r["f_state"] for r in rows},
                "f_state_span_ladder": None, "f_state_at_dx12": next(
                    (r["f_state"] for r in rows
                     if r["dx_um"] == PRODUCTION_DX_UM), None),
                "f_state_finest": rows[-1]["f_state"],
                "n_rungs": len(rows)}
            continue
        fs = [r["f_state"] for r in lad]
        allfs = [r["f_state"] for r in rows]
        d = np.diff(fs)
        summary[key] = {
            "f_state_by_dx": {str(r["dx_um"]): r["f_state"] for r in rows},
            "mean_o2_by_dx": {str(r["dx_um"]): r["mean_o2_tumour"] for r in rows},
            "source_integral_by_dx": {str(r["dx_um"]): r["source_integral"]
                                      for r in rows},
            "f_state_span_ladder": float(max(fs) - min(fs)),
            "f_state_span_all_rungs": float(max(allfs) - min(allfs)),
            "mean_o2_span_ladder": float(
                max(r["mean_o2_tumour"] for r in lad)
                - min(r["mean_o2_tumour"] for r in lad)),
            "monotone_in_dx": bool(np.all(d > 0) or np.all(d < 0)),
            "n_direction_reversals": int(np.sum(np.sign(d[1:]) != np.sign(d[:-1]))),
            "f_state_at_dx12": next((r["f_state"] for r in rows
                                     if r["dx_um"] == PRODUCTION_DX_UM), None),
            "f_state_finest": rows[-1]["f_state"],
            "err_at_dx12_vs_finest": (
                next((r["f_state"] for r in rows
                      if r["dx_um"] == PRODUCTION_DX_UM), float("nan"))
                - rows[-1]["f_state"]),
            "vasc_is_binary": bool(lad[0]["vasc_is_binary"]),
            "n_rungs": len(rows),
        }
    return summary


# --------------------------------------------------------------------------- #
# Stage 2 -- the published ratios, regenerated
# --------------------------------------------------------------------------- #
def stage_ratios(dx_um: float = PRODUCTION_DX_UM,
                 tissues: Sequence[str] = ("gbm", "breast"),
                 arm_list: Sequence[str] = ("binary", "area")) -> Dict:
    """The nine published rows, RE-SOLVED under both source paintings.

    ``experiments/vessel_source_conservation.py`` §7 replayed the ``binary``
    control from ``delivery_source_geometry``'s checkpoint. Here it is re-solved
    in the same process, off the same reconstruction, with the same alpha, so
    (a) the reproduction of ``PUBLISHED_RATIO`` is a live check rather than a
    tautology and (b) the per-tumour-cell exposure field exists for BOTH arms,
    which is what ``results/p99_convention.md``'s normaliser family needs.

    Nothing here is re-scored off a cache: the source field enters the transport
    PDE, so the exposure must be re-solved. What IS re-scored off the solve is
    the normaliser family -- ``decompose`` reads only ``e_unit``, so all
    fourteen conventions come from one set of solves.
    """
    import experiments.delivery_source_geometry as dsg
    import tmesim.validation.kill_calibration as kc
    from experiments.p99_convention import (
        HOMOGENEOUS_TAGS, SHIPPED_TAG, decompose_general)
    from experiments.vessel_source_conservation import source_field
    from tmesim.reconstruct.vessel_network import LINK_DISTANCE_UM
    from tmesim.validation.calibration import calibrate_alpha

    alpha = float(calibrate_alpha())
    if alpha != ALPHA_EXPECTED:
        raise AssertionError(f"alpha drifted: {alpha!r} != {ALPHA_EXPECTED!r}")
    log(f"alpha = {alpha!r} (guard OK)")

    original_source = dsg._source_field
    original_decompose = kc.decompose
    captured: List[Dict] = []

    def patched_source(micro, arm: str) -> np.ndarray:
        if arm in ("area", "area_matched"):
            v = dsg._segment_safe(micro.cells, micro.domain,
                                  float(LINK_DISTANCE_UM))
            return source_field(micro.domain, v, arm)
        return original_source(micro, arm)

    def capturing_decompose(e_unit, sens, resist, c_bath, ic50, hill,
                            max_kill=1.0, **kw):
        """Record the scoring inputs, then delegate to the SHIPPED function.

        This changes nothing: the return value is the shipped one, byte for
        byte. It exists so the normaliser family can be scored on exactly the
        arrays the published verdict was taken on.
        """
        captured.append({
            "e_unit": np.asarray(e_unit, dtype=np.float64).copy(),
            "sens": np.asarray(sens, dtype=np.float64).copy(),
            "resist": np.asarray(resist, dtype=np.float64).copy(),
            "c_bath": float(c_bath), "ic50": float(ic50), "hill": float(hill),
            "max_kill": float(max_kill),
        })
        return original_decompose(e_unit, sens, resist, c_bath, ic50, hill,
                                  max_kill, **kw)

    ck = os.path.join(SWEEP, f"areagen_ratios_dx{dx_um:g}.jsonl")
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

    dsg._source_field = patched_source
    kc.decompose = capturing_decompose
    panels: List[Dict] = []
    try:
        for tissue in tissues:
            drugs = (dsg.BREAST_DRUGS if tissue == "breast" else dsg.GBM_DRUGS)
            for arm in arm_list:
                key = f"{tissue}|{arm}"
                if key in done:
                    log(f"{key}: cached (normaliser family NOT available)")
                    panels.append(done[key])
                    continue
                log(f"=== panel {key} at dx={dx_um:g} (alpha={alpha:.6f}) ===")
                t0 = time.time()
                captured.clear()
                p = dsg._run_panel(tissue, arm, alpha, drugs, float(dx_um))
                # one capture per drug, in order
                if len(captured) == len(p["rows"]):
                    for row, cap in zip(p["rows"], captured):
                        row["normalisers"] = _norm_family(
                            cap, decompose_general, HOMOGENEOUS_TAGS, SHIPPED_TAG,
                            shipped_ratio=row["ratio_exact"])
                else:  # pragma: no cover - defensive
                    log(f"WARNING: {len(captured)} captures for "
                        f"{len(p['rows'])} rows; normaliser family skipped")
                panels.append(p)
                with open(ck, "a") as fh:
                    fh.write(json.dumps(p, default=float) + "\n")
                log(f"=== panel {key} done in {time.time() - t0:.0f}s ===")
    finally:
        dsg._source_field = original_source
        kc.decompose = original_decompose

    for p in panels:
        n = dsg.rescore_panel(p)
        log(f"{p['tissue']}|{p['arm']}: {n} verdict(s) re-derived through the "
            f"shared band decision")

    out = _ratio_table(panels, dsg, alpha, dx_um, arm_list)
    _dump(out, f"ratios_dx{dx_um:g}.json")
    return out


def _norm_family(cap: Dict, decompose_general, tags, shipped_tag: str,
                 shipped_ratio: float) -> Dict:
    """Score one drug row under every homogeneous normaliser convention."""
    fam = {}
    for tag in tags:
        d = decompose_general(cap["e_unit"], cap["sens"], cap["resist"],
                              cap["c_bath"], cap["ic50"], cap["hill"],
                              cap["max_kill"], tag=tag)
        fam[tag] = {"ratio_exact": float(d["ratio_exact"]),
                    "f_penetration": float(d["f_penetration"]),
                    "ref_unit": float(d["ref_unit_uM_per_uM_plasma"])}
    # the shipped convention must reproduce the shipped decompose exactly
    fam["_shipped_reproduced"] = bool(
        abs(fam[shipped_tag]["ratio_exact"] - float(shipped_ratio)) < 1e-12)
    return fam


def _ratio_table(panels, dsg, alpha: float, dx_um: float,
                 arm_list: Sequence[str]) -> Dict:
    by = {(p["tissue"], p["arm"]): p for p in panels}
    table: List[Dict] = []
    for (tissue, drug), pub_ratio in dsg.PUBLISHED_RATIO.items():
        rec: Dict = {"tissue": tissue, "drug": drug,
                     "published_ratio": float(pub_ratio),
                     "published_verdict": dsg.PUBLISHED_VERDICT[(tissue, drug)]}
        for arm in arm_list:
            p = by.get((tissue, arm))
            if p is None:
                continue
            r = next((x for x in p["rows"] if x["drug"] == drug), None)
            if r is None:
                continue
            rec[arm] = {
                "ratio": float(r["ratio_exact"]),
                "verdict": r["verdict"],
                "verdict_scored": bool(r["verdict_scored"]),
                "f_penetration": float(r["f_penetration"]),
                "f_state": float(r["f_state"]),
                "e_unit_mean": float(r["e_unit_mean"]),
                "e_unit_p99": float(r["e_unit_p99"]),
                "implied_ic50_fold": float(r["implied_ic50_fold"]),
                "normalisers": r.get("normalisers"),
            }
        if "binary" in rec and "area" in rec:
            rec["ratio_pct_move"] = 100.0 * (rec["area"]["ratio"]
                                             - rec["binary"]["ratio"]) / rec["binary"]["ratio"]
            rec["verdict_flip"] = bool(
                rec["area"]["verdict_scored"]
                and rec["area"]["verdict"] != rec["binary"]["verdict"])
        table.append(rec)

    control_reproduces = all(
        abs(rec["binary"]["ratio"] - rec["published_ratio"]) < 5e-5
        for rec in table if "binary" in rec)
    worst_control = max((abs(rec["binary"]["ratio"] - rec["published_ratio"])
                         for rec in table if "binary" in rec), default=float("nan"))
    flips = sum(1 for rec in table if rec.get("verdict_flip"))
    worst_move = max((abs(rec["ratio_pct_move"]) for rec in table
                      if "ratio_pct_move" in rec), default=float("nan"))

    # normaliser family: does the area change move any convention's verdict?
    norm_moves: Dict[str, Dict] = {}
    for rec in table:
        if not (rec.get("binary", {}).get("normalisers")
                and rec.get("area", {}).get("normalisers")):
            continue
        for tag in rec["binary"]["normalisers"]:
            if tag.startswith("_"):
                continue
            b = rec["binary"]["normalisers"][tag]["ratio_exact"]
            a = rec["area"]["normalisers"][tag]["ratio_exact"]
            d = norm_moves.setdefault(tag, {"worst_pct_move": 0.0, "rows": {}})
            pct = 100.0 * (a - b) / b if b else float("nan")
            d["rows"][f"{rec['tissue']}|{rec['drug']}"] = {
                "binary": b, "area": a, "pct_move": pct}
            if abs(pct) > abs(d["worst_pct_move"]):
                d["worst_pct_move"] = pct

    out = {"dx_um": float(dx_um), "alpha": alpha, "arms": list(arm_list),
           "table": table, "panels": panels,
           "scored_verdict_flips": int(flips),
           "worst_abs_ratio_pct_move": float(worst_move),
           "binary_control_reproduces_published": bool(control_reproduces),
           "worst_binary_control_deviation": float(worst_control),
           "normaliser_family": norm_moves}
    for rec in table:
        line = (f"{rec['tissue']}/{rec['drug']}: pub {rec['published_ratio']:.4f} "
                f"[{rec['published_verdict']}]")
        for arm in arm_list:
            if arm in rec:
                line += f" | {arm} {rec[arm]['ratio']:.6f} {rec[arm]['verdict']}"
        if "ratio_pct_move" in rec:
            line += f" | move {rec['ratio_pct_move']:+.3f} %"
        log(line)
    log(f"scored verdict flips: {flips}; worst |ratio move| {worst_move:.3f} %; "
        f"binary control reproduces published: {control_reproduces} "
        f"(worst dev {worst_control:.2e})")
    for tag in sorted(norm_moves):
        log(f"  normaliser {tag:8s} worst ratio move "
            f"{norm_moves[tag]['worst_pct_move']:+.3f} %")
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", required=True,
                    choices=("alpha", "closure", "ratios"))
    ap.add_argument("--dx", type=float, default=PRODUCTION_DX_UM)
    ap.add_argument("--rungs", type=str, default="")
    ap.add_argument("--no-beta-cap", action="store_true")
    a = ap.parse_args(argv)
    if a.stage == "alpha":
        stage_alpha()
    elif a.stage == "closure":
        rungs = (tuple(float(x) for x in a.rungs.split(","))
                 if a.rungs else DX_LADDER + DX_FINE)
        stage_closure(rungs=rungs, beta_cap=not a.no_beta_cap)
    elif a.stage == "ratios":
        stage_ratios(dx_um=float(a.dx))
    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
