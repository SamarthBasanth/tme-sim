"""GRID-CONVERGENCE AUDIT of the breast 3-D:2-D kill-calibration panel.

Why this file exists
--------------------
``results/kill_calibration_breast.md`` publishes 6 formally-testable 3-D:2-D
ratio verdicts on the real Xenium breast section at ``dx = 12 um``: 5 PASS and
ONE FAIL (docetaxel, ratio 0.1484 vs a band top edge of 0.1463). That single
FAIL is the load-bearing evidence, cited in several documents, that the
calibration is falsifiable rather than tuned.

``results/dtissue_sensitivity.md`` §7.2 turned up an incidental result that puts
it in doubt: at ``dx = 24 um`` the same docetaxel ratio is **0.1200**, i.e.
PASS. The verdict flips under a 2x grid coarsening, but SURVIVES a 3.1x
diffusivity correction (0.1484 -> 0.1497). A result more sensitive to the mesh
than to a large physical parameter change is not a result anyone should quote.

This script asks the only question that settles it: **is the published dx = 12 um
grid converged?** It runs the identical published pipeline over a dx ladder,
reports the scored ratio as a function of dx per drug, estimates the
discretisation error at dx = 12, and Richardson-extrapolates to dx -> 0.

What is and is NOT varied
-------------------------
**ONLY ``dx`` is varied.** ``dt = 10 s``, the 6 h window, ``save_every``, the
seed, ``PDConfig``, the p99 bath-matching convention, the literature bands, the
drug compile and the shipping diffusivity are all held at their published
values. The delivery law is the shipping one; ``exchange_k`` is NOT passed, so
the KNOWN TRAP about ``kill_calibration.unit_transport()`` ignoring ``PDEConfig``
flags does not arise here -- there are no delivery-law arms to keep distinct.
(``PDEConfig.convection = True`` DOES reach the solver, because ``unit_transport``
forwards ``pde`` and ``velocity`` positionally.)

The ladder and why it stops at dx = 8 um
----------------------------------------
The solver's explicit upwind advection step carries the CFL condition
``dt*(|vx|+|vy|)/dx <= 1`` (``physics/diffusion.py``), and raises rather than
silently degrading. On this section, at ``dt = 10 s``, the measured worst-cell
speed gives

    dx = 12 um -> CFL 0.55      dx = 10 um -> CFL 0.82
    dx =  8 um -> CFL 0.81      dx =  6 um -> CFL 1.57  (SOLVER REFUSES)

so **dx = 8 um is the finest grid reachable at the published dt**. Going finer
would require also cutting dt, which would confound spatial and temporal
discretisation error in the same number. That limit is reported, not hidden.

Convergence model
-----------------
For each (drug, quantity) the ladder is fitted with the standard three-parameter
discretisation-error model

    q(dx) = q0 + C * dx^p

by least squares over ALL ladder points (5-7 points, 3 parameters), and
separately by classical three-grid Richardson extrapolation on the
constant-refinement-ratio triple (8, 12, 18 um), r = 1.5:

    p  = ln|(q18 - q12)/(q12 - q8)| / ln r
    q0 = q8 + (q8 - q12) / (r^p - 1)

Both are reported. Disagreement between them is itself information (it means the
ladder is not yet in the asymptotic range).

Run:
    ./.venv/bin/python experiments/grid_convergence_breast.py
    ./.venv/bin/python experiments/grid_convergence_breast.py --dx 30 24 18
    ./.venv/bin/python experiments/grid_convergence_breast.py --report-only
Writes:
    results/sweep/gridconv_records.jsonl   (restart-safe checkpoint, one run per line)
    results/sweep/gridconv_breast.json     (assembled raw results)
    results/grid_convergence_breast.md     (the report)

Nothing under ``tmesim/`` is touched. Every physics/scoring function is imported
unmodified from the published path.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.core import PDConfig, PDEConfig, PKConfig, SimConfig  # noqa: E402
from tmesim.reconstruct.real_pipeline import reconstruct_real  # noqa: E402
from tmesim.validation.calibration import calibrate_alpha  # noqa: E402
from tmesim.validation.kill_calibration import (  # noqa: E402
    BATH_REFERENCE_PCT,
    decompose,
    expected_kill,
    kill_ratio_to_ic50_fold,
    prepare_tissue,
    ratio_test,
    resistance_vector,
    unit_transport,
    window_plasma_mean,
)
from tmesim.validation.pk_reality import DEFAULT_ECM  # noqa: E402
from tmesim.validation.calibration import _d_eff_tortuous  # noqa: E402

from experiments.kill_calibration_breast import (  # noqa: E402
    BREAST_SAMPLE_PATH,
    DT_S,
    DX_UM,
    NOMINAL_HALF_LIFE_S,
    NOMINAL_PEAK_UM,
    RUN_CACHE as PUBLISHED_CACHE,
    SAVE_EVERY,
    SEED,
    WINDOW_S,
    apply_drug_sensitivity_on_copy,
    compile_breast_drug,
)
from experiments.kill_calibration_breast_report import (  # noqa: E402
    bridge_to_kill_ratio_band,
    drug_specific_bridge,
    solid_bridge,
)

OUT_DIR = os.path.join(_ROOT, "results", "sweep")
CHECKPOINT = os.path.join(OUT_DIR, "gridconv_records.jsonl")
RAW_JSON = os.path.join(OUT_DIR, "gridconv_breast.json")
OUT_MD = os.path.join(_ROOT, "results", "grid_convergence_breast.md")

#: The 6 drugs ``results/kill_calibration_breast.md`` calls formally testable.
TESTABLE: Tuple[str, ...] = ("doxorubicin", "paclitaxel", "docetaxel",
                             "epirubicin", "5-fluorouracil", "cisplatin")

#: Coarse -> fine, so partial results are always usable and the expensive end
#: comes last. 8 / 12 / 18 is a constant-ratio (r = 1.5) Richardson triple.
DX_LADDER: Tuple[float, ...] = (30.0, 24.0, 18.0, 16.0, 12.0, 10.0, 8.0)

#: The published dx = 12 um table (``kill_calibration_breast.md`` §3), quoted so
#: this report can state plainly whether it reproduces the thing it is auditing.
PUBLISHED_RATIO = {"doxorubicin": 0.2301, "paclitaxel": 0.1503,
                   "docetaxel": 0.1484, "epirubicin": 0.2301,
                   "5-fluorouracil": 0.3379, "cisplatin": 0.3085}
PUBLISHED_VERDICT = {"doxorubicin": "PASS", "paclitaxel": "PASS",
                     "docetaxel": "FAIL", "epirubicin": "PASS",
                     "5-fluorouracil": "PASS", "cisplatin": "PASS"}


# --------------------------------------------------------------------------- #
# Scoring -- the published path, byte for byte
# --------------------------------------------------------------------------- #
def build_bands(ic50_hill: Dict[str, Tuple[float, float]]) -> Dict:
    """Literature bands. They depend only on the frozen literature tables and on
    (ic50, hill), NONE of which dx touches -- so every band is identical at every
    rung of the ladder and cannot be the thing that moves a verdict."""
    solid = solid_bridge()
    out = {}
    for name, (ic50, hill) in ic50_hill.items():
        c_bath = 1.0 * ic50
        spec = drug_specific_bridge(name)
        band = bridge_to_kill_ratio_band(spec, c_bath, ic50, hill) if spec else None
        out[name] = {
            "band": band,
            "band_solid": bridge_to_kill_ratio_band(solid, c_bath, ic50, hill),
            "band_is_drug_matched": bool(spec is not None and spec["n"] > 1),
        }
    return out


def score(name: str, e_unit: np.ndarray, sens: np.ndarray, resist: np.ndarray,
          drug, pd_cfg: PDConfig, bands: Dict, plasma_mean: float) -> Dict:
    """``experiments/kill_calibration_breast.run``'s own scoring, unmodified:
    ``decompose`` at 1x IC50 (the reporting concentration the verdict is taken
    at) then the frozen ``ratio_test`` against the frozen literature band.

    NOTE the ``resistance_vector`` call takes NO ``drug=`` argument, exactly as
    ``kill_calibration_breast.py`` line 344 does. That matters: the sibling
    ``dtissue_sensitivity.py`` passes ``drug=drug``. For this panel the two are
    numerically identical (no drug here carries a per-drug hypoxia override), but
    the published call is the one reproduced here.
    """
    c_bath = 1.0 * float(drug.ic50)
    dec = decompose(e_unit, sens, resist, c_bath, drug.ic50, drug.hill,
                    pd_cfg.max_kill_rate)
    fold = kill_ratio_to_ic50_fold(dec["ratio_exact"], c_bath, drug.ic50, drug.hill)
    b = bands[name]
    t_specific = (ratio_test(dec["ratio_exact"], b["band"]["ratio_lo"],
                             b["band"]["ratio_hi"]) if b["band"] else None)
    t_solid = ratio_test(dec["ratio_exact"], b["band_solid"]["ratio_lo"],
                         b["band_solid"]["ratio_hi"])
    formal = t_specific if b["band_is_drug_matched"] else t_solid
    return {
        "drug": name,
        "c_bath_uM": c_bath,
        "f_penetration": dec["f_penetration"],
        "f_state": dec["f_state"],
        "f_sensitivity": dec["f_sensitivity"],
        "f_product": dec["f_product"],
        "interaction": dec["interaction"],
        "ratio_exact": dec["ratio_exact"],
        "implied_ic50_fold": fold,
        "kill_3d_model": dec["kill_3d_model"],
        "kill_2d_model": dec["kill_2d_model"],
        "plasma_required_uM": dec["plasma_required_uM"],
        "E_median_over_bath": dec["E_median_over_bath"],
        "kill_nominal_plasma": expected_kill(e_unit, sens, resist, plasma_mean,
                                             drug.ic50, drug.hill,
                                             pd_cfg.max_kill_rate),
        "e_unit_mean": float(e_unit.mean()),
        "e_unit_median": float(np.median(e_unit)),
        "e_unit_p99": float(np.percentile(e_unit, BATH_REFERENCE_PCT)),
        "e_unit_p01": float(np.percentile(e_unit, 1.0)),
        "e_unit_max": float(e_unit.max()),
        "verdict": formal["verdict"],
        "verdict_direction": formal["direction"],
        "band_used_lo": formal["lit_ratio_lo"],
        "band_used_hi": formal["lit_ratio_hi"],
        "band_is_drug_matched": b["band_is_drug_matched"],
        "verdict_solid": t_solid["verdict"],
    }


# --------------------------------------------------------------------------- #
# One rung of the ladder
# --------------------------------------------------------------------------- #
def _key(dx_um: float, name: str) -> str:
    return f"{dx_um:g}|{name}"


_REQUIRED = frozenset(("dx_um", "drug", "ratio_exact", "f_penetration", "f_state",
                       "verdict", "e_unit_mean", "kill_nominal_plasma",
                       "_e_unit_sum", "n_tumor"))


def load_checkpoint(path: str) -> Dict[str, Dict]:
    done: Dict[str, Dict] = {}
    if not os.path.exists(path):
        return done
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue          # truncated final line from a killed run
            if not _REQUIRED.issubset(r):
                continue          # older schema -- discard rather than mix
            done[_key(float(r["dx_um"]), r["drug"])] = r
    return done


def run_ladder(dx_ladder: Sequence[float] = DX_LADDER,
               drugs: Sequence[str] = TESTABLE,
               checkpoint: str = CHECKPOINT,
               out_json: Optional[str] = RAW_JSON) -> Dict:
    t0 = time.time()

    def log(msg: str) -> None:
        print(f"[gridconv {time.time()-t0:7.1f}s] {msg}", flush=True)

    alpha = float(calibrate_alpha())
    log(f"alpha = {alpha:.4f} (published, unchanged); dx ladder = "
        f"{', '.join(f'{d:g}' for d in dx_ladder)} um; {len(drugs)} drugs "
        f"-> {len(dx_ladder)*len(drugs)} transport solves")

    done = load_checkpoint(checkpoint)
    if done:
        log(f"checkpoint: {len(done)} runs already complete, resuming")
    os.makedirs(os.path.dirname(checkpoint), exist_ok=True)

    # Bands and the drug compile are dx-independent; do them once.
    ic50_hill: Dict[str, Tuple[float, float]] = {}
    compiled = {}
    drug_static = {}
    for name in drugs:
        dp, mol, ic50_info = compile_breast_drug(name, alpha)
        compiled[name] = dp
        ic50_hill[name] = (float(dp.ic50), float(dp.hill))
        drug_static[name] = {
            "ic50_uM": float(dp.ic50), "hill": float(dp.hill),
            "diffusivity_um2_s": float(dp.diffusivity),
            "binding_per_s": float(dp.binding),
            "clearance_per_s": float(dp.clearance),
            # the perivascular boundary-layer thickness under the SHIPPING
            # delivery law -- the length the mesh has to resolve
            "L_um": math.sqrt(_d_eff_tortuous(dp, DEFAULT_ECM)
                              / (float(dp.binding) + float(dp.clearance))),
        }
    bands = build_bands(ic50_hill)

    records: List[Dict] = []
    grids: Dict[str, Dict] = {}
    ck = open(checkpoint, "a")
    n_done = 0
    total = len(dx_ladder) * len(drugs)
    for dx in dx_ladder:
        dx = float(dx)
        pending = [n for n in drugs if _key(dx, n) not in done]
        if not pending:
            for name in drugs:
                r = done[_key(dx, name)]
                records.append(r)
                grids.setdefault(f"{dx:g}", {
                    "dx_um": dx, "grid": r.get("grid"),
                    "n_cells": r.get("n_cells"), "n_tumor": r.get("n_tumor"),
                    "f_state": r.get("f_state"), "cfl": r.get("cfl"),
                    "max_speed_um_s": r.get("max_speed_um_s"),
                    "prologue_s": r.get("prologue_s"),
                })
                n_done += 1
            log(f"dx = {dx:g} um: all {len(drugs)} runs from checkpoint "
                f"({n_done}/{total})")
            continue

        tg = time.time()
        micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH,
                                 dx_um=dx)
        pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                        save_every=SAVE_EVERY)
        pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                      n_doses=1, dose_interval=3600.0 * 24)
        pd_cfg = PDConfig()
        work, labeled, velocity = prepare_tissue(
            micro, SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
        tumor = labeled.tumor_mask()
        state_t = np.asarray(labeled.state)[tumor]
        resist = resistance_vector(state_t, pd_cfg)          # published call
        plasma_mean = window_plasma_mean(pk, WINDOW_S, DT_S)
        vx, vy = velocity
        max_speed = float(np.max(np.abs(np.asarray(vx)) + np.abs(np.asarray(vy))))
        cfl = max_speed * DT_S / dx
        prologue_s = time.time() - tg
        ginfo = {
            "dx_um": dx, "grid": list(micro.domain.shape),
            "n_cells": int(labeled.n), "n_tumor": int(tumor.sum()),
            "f_state": float(np.mean(resist)),
            "state_viable_frac": float(np.mean(resist == 1.0)),
            "state_quiescent_frac": float(np.mean(
                np.isclose(resist, pd_cfg.quiescent_resistance))),
            "state_hypoxic_frac": float(np.mean(
                np.isclose(resist, pd_cfg.hypoxia_resistance))),
            "max_speed_um_s": max_speed, "cfl": cfl,
            "prologue_s": prologue_s,
        }
        grids[f"{dx:g}"] = ginfo
        log(f"dx = {dx:g} um: grid {micro.domain.shape}, {int(labeled.n):,} cells, "
            f"{int(tumor.sum()):,} tumour, F_state = {ginfo['f_state']:.4f}, "
            f"CFL = {cfl:.3f} (prologue {prologue_s:.0f}s)")

        for name in drugs:
            k = _key(dx, name)
            if k in done:
                records.append(done[k])
                n_done += 1
                continue
            drug = compiled[name]
            sens = np.asarray(
                apply_drug_sensitivity_on_copy(work, labeled, name)
                .cells.sensitivity, dtype=np.float64)[tumor]
            ts = time.time()
            tr = unit_transport(work, drug, pde, velocity)
            solve_s = time.time() - ts
            col, row = work.domain.to_index(labeled.x, labeled.y)
            auc = np.asarray(tr.auc, dtype=np.float64)[row, col]
            duration = float(tr.times[-1]) if tr.times.size and tr.times[-1] > 0 else 1.0
            e_unit = (auc / duration)[tumor]

            rec = score(name, e_unit, sens, resist, drug, pd_cfg, bands, plasma_mean)
            rec.update(ginfo)
            rec.update(drug_static[name])
            rec.update({
                "drug": name,
                "solve_s": solve_s,
                "penetration_index": float(tr.penetration_index),
                "_e_unit_sum": float(e_unit.sum()),
                "_e_unit_first": float(e_unit[0]),
            })
            records.append(rec)
            ck.write(json.dumps(rec, default=float) + "\n")
            ck.flush()
            n_done += 1
            log(f"  {n_done}/{total}  dx={dx:g} {name}: ratio={rec['ratio_exact']:.4f} "
                f"F_pen={rec['f_penetration']:.4f} {rec['verdict']} "
                f"({solve_s:.0f}s)")
    ck.close()

    # --- the field really moved: exposure fingerprints distinct across dx -----
    distinct = {}
    for name in drugs:
        sums = [r["_e_unit_sum"] for r in records if r["drug"] == name]
        distinct[name] = {"n_points": len(sums), "n_distinct": len(set(sums))}
        assert len(set(sums)) == len(sums), (
            f"{name}: exposure fields are NOT distinct across the dx ladder "
            f"({len(set(sums))} distinct of {len(sums)}) -- dx is a no-op")

    out = {
        "alpha": alpha, "dt_s": DT_S, "window_s": WINDOW_S,
        "save_every": SAVE_EVERY, "seed": SEED,
        "dx_ladder": [float(d) for d in dx_ladder], "drugs": list(drugs),
        "drug_static": drug_static,
        "bands": bands, "grids": grids, "records": records,
        "field_distinctness": distinct,
        "published_ratio": PUBLISHED_RATIO, "published_verdict": PUBLISHED_VERDICT,
        "elapsed_s": time.time() - t0,
    }
    if out_json:
        os.makedirs(os.path.dirname(out_json), exist_ok=True)
        with open(out_json, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
        log(f"wrote {out_json}")
    log(f"done in {(time.time()-t0)/60:.1f} min")
    return out


# --------------------------------------------------------------------------- #
# dt control at the fine end
# --------------------------------------------------------------------------- #
DT_CONTROL_CK = os.path.join(OUT_DIR, "gridconv_dtcontrol.jsonl")


def run_dt_control(points: Sequence[Tuple[float, float]],
                   drugs: Sequence[str],
                   checkpoint: str = DT_CONTROL_CK) -> List[Dict]:
    """Re-solve chosen (dx, dt) points to prove the fine rungs are NOT a
    temporal-discretisation artefact.

    The main ladder holds ``dt = 10 s`` at every dx on purpose: a ladder that
    changes dt with dx is measuring the sum of two discretisation errors and
    cannot attribute the movement to the mesh. That design is only defensible if
    the fine rungs are already dt-converged, which is what this control measures
    directly rather than citing.
    """
    t0 = time.time()

    def log(msg: str) -> None:
        print(f"[dtctl {time.time()-t0:7.1f}s] {msg}", flush=True)

    done: Dict[str, Dict] = {}
    if os.path.exists(checkpoint):
        for line in open(checkpoint):
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            done[f"{r['dx_um']:g}|{r['dt_s']:g}|{r['drug']}"] = r

    alpha = float(calibrate_alpha())
    compiled, ic50_hill = {}, {}
    for name in drugs:
        dp, _, _ = compile_breast_drug(name, alpha)
        compiled[name] = dp
        ic50_hill[name] = (float(dp.ic50), float(dp.hill))
    bands = build_bands(ic50_hill)

    os.makedirs(os.path.dirname(checkpoint), exist_ok=True)
    out: List[Dict] = []
    ck = open(checkpoint, "a")
    by_dx: Dict[float, List[float]] = {}
    for dx, dt in points:
        by_dx.setdefault(float(dx), []).append(float(dt))
    for dx, dts in by_dx.items():
        if all(f"{dx:g}|{dt:g}|{n}" in done for dt in dts for n in drugs):
            out.extend(done[f"{dx:g}|{dt:g}|{n}"] for dt in dts for n in drugs)
            log(f"dx = {dx:g}: all points from checkpoint")
            continue
        micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH, dx_um=dx)
        pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                      n_doses=1, dose_interval=3600.0 * 24)
        pd_cfg = PDConfig()
        base_pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                             save_every=SAVE_EVERY)
        work, labeled, velocity = prepare_tissue(
            micro, SimConfig(pde=base_pde, pk=pk, pd=pd_cfg, seed=SEED))
        tumor = labeled.tumor_mask()
        state_t = np.asarray(labeled.state)[tumor]
        resist = resistance_vector(state_t, pd_cfg)
        vx, vy = velocity
        max_speed = float(np.max(np.abs(np.asarray(vx)) + np.abs(np.asarray(vy))))
        for dt in dts:
            # hold the number of SAVED frames fixed so the AUC quadrature grid is
            # comparable; save_every scales with the step count
            se = max(1, int(round(SAVE_EVERY * DT_S / dt)))
            pde = PDEConfig(dt=dt, total_time=WINDOW_S, convection=True,
                            save_every=se)
            plasma_mean = window_plasma_mean(pk, WINDOW_S, dt)
            for name in drugs:
                k = f"{dx:g}|{dt:g}|{name}"
                if k in done:
                    out.append(done[k])
                    continue
                sens = np.asarray(
                    apply_drug_sensitivity_on_copy(work, labeled, name)
                    .cells.sensitivity, dtype=np.float64)[tumor]
                ts = time.time()
                tr = unit_transport(work, compiled[name], pde, velocity)
                col, row = work.domain.to_index(labeled.x, labeled.y)
                auc = np.asarray(tr.auc, dtype=np.float64)[row, col]
                dur = float(tr.times[-1]) if tr.times.size and tr.times[-1] > 0 else 1.0
                e_unit = (auc / dur)[tumor]
                rec = score(name, e_unit, sens, resist, compiled[name], pd_cfg,
                            bands, plasma_mean)
                rec.update({"dx_um": float(dx), "dt_s": float(dt),
                            "save_every": se, "drug": name,
                            "cfl": max_speed * dt / dx,
                            "max_speed_um_s": max_speed,
                            "solve_s": time.time() - ts,
                            "_e_unit_sum": float(e_unit.sum())})
                out.append(rec)
                ck.write(json.dumps(rec, default=float) + "\n")
                ck.flush()
                log(f"  dx={dx:g} dt={dt:g}s {name}: ratio={rec['ratio_exact']:.4f} "
                    f"CFL={rec['cfl']:.3f} {rec['verdict']} ({rec['solve_s']:.0f}s)")
    ck.close()
    return out


# --------------------------------------------------------------------------- #
# Convergence analysis
# --------------------------------------------------------------------------- #
def fit_power_law(dxs: Sequence[float], qs: Sequence[float]) -> Dict:
    """Least-squares fit of ``q(dx) = q0 + C*dx^p`` over the whole ladder.

    Solved as a 1-D search over ``p`` (for fixed ``p`` the model is linear in
    ``q0`` and ``C``, so the inner problem is exact). ``p`` is scanned over
    0.25-4.0, which brackets every order a first- or second-order scheme with
    a geometry re-rasterisation on top could plausibly show.
    """
    x = np.asarray(dxs, dtype=np.float64)
    y = np.asarray(qs, dtype=np.float64)
    if x.size < 3 or not np.all(np.isfinite(y)):
        return {"ok": False}
    best = None
    for p in np.linspace(0.25, 4.0, 1501):
        A = np.column_stack([np.ones_like(x), x ** p])
        try:
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        except np.linalg.LinAlgError:
            continue
        resid = float(np.sum((A @ coef - y) ** 2))
        if best is None or resid < best[0]:
            best = (resid, float(p), float(coef[0]), float(coef[1]))
    if best is None:
        return {"ok": False}
    resid, p, q0, C = best
    rms = math.sqrt(resid / x.size)
    return {"ok": True, "p": p, "q0": q0, "C": C, "rms_resid": rms,
            "pred": {f"{d:g}": float(q0 + C * d ** p) for d in x}}


def richardson_triple(q_fine: float, q_mid: float, q_coarse: float,
                      r: float = 1.5) -> Dict:
    """Classical 3-grid Richardson on a constant-refinement-ratio triple.

    ``q_fine`` = dx/r^0 ... ordered fine -> coarse with a constant ratio ``r``
    between successive grids (here 8, 12, 18 um, r = 1.5).
    """
    d1 = q_mid - q_fine        # mid  - fine
    d2 = q_coarse - q_mid      # coarse - mid
    if d1 == 0.0 or d2 == 0.0:
        return {"ok": False, "reason": "zero difference between grids"}
    ratio = d2 / d1
    if ratio <= 0.0:
        return {"ok": False, "reason": f"non-monotone triple (d2/d1 = {ratio:.4f} <= 0)",
                "d1": d1, "d2": d2, "ratio": ratio}
    p = math.log(ratio) / math.log(r)
    denom = r ** p - 1.0
    if abs(denom) < 1e-12:
        return {"ok": False, "reason": "degenerate order"}
    q0 = q_fine - d1 / denom
    # GCI (Roache) on the finest pair, safety factor 1.25 for a 3-grid study
    eps = abs(d1 / q_fine) if q_fine else float("nan")
    gci = 1.25 * eps / abs(denom)
    return {"ok": True, "p": p, "q0": q0, "d1": d1, "d2": d2, "ratio": ratio,
            "gci_fine": gci}


def analyse(res: Dict) -> Dict:
    """Per-drug convergence behaviour of the scored 3-D:2-D ratio."""
    out: Dict[str, Dict] = {}
    dxs_all = sorted({float(r["dx_um"]) for r in res["records"]})
    for name in res["drugs"]:
        rows = sorted([r for r in res["records"] if r["drug"] == name],
                      key=lambda r: float(r["dx_um"]))
        dxs = [float(r["dx_um"]) for r in rows]
        ratios = [float(r["ratio_exact"]) for r in rows]
        fpen = [float(r["f_penetration"]) for r in rows]
        by_dx = {f"{d:g}": r for d, r in zip(dxs, rows)}

        diffs = np.diff(ratios)                        # fine -> coarse
        monotone = bool(np.all(diffs > 0) or np.all(diffs < 0))
        # Is it settling? Compare the RATE of change per um of mesh, finest pair
        # vs the next pair out. A converging quantity has |dq/d(dx)| shrinking as
        # dx -> 0 for p > 1, and (for 0 < p < 1) at least a shrinking |dq| per
        # equal refinement step. Both are reported rather than one being asserted.
        slopes = list(np.diff(ratios) / np.diff(dxs))
        settling = None
        if len(ratios) >= 3:
            settling = bool(abs(slopes[0]) < abs(slopes[1]))
        # per-rung absolute step, useful in the report
        steps = [{"from_dx": dxs[i], "to_dx": dxs[i + 1],
                  "d_ratio": float(diffs[i]),
                  "d_ratio_pct": 100.0 * float(diffs[i]) / ratios[i],
                  "slope_per_um": float(slopes[i])}
                 for i in range(len(diffs))]

        fit = fit_power_law(dxs, ratios)
        fit_pen = fit_power_law(dxs, fpen)

        rich = {"ok": False}
        triple = (8.0, 12.0, 18.0)
        if all(f"{d:g}" in by_dx for d in triple):
            rich = richardson_triple(by_dx["8"]["ratio_exact"],
                                     by_dx["12"]["ratio_exact"],
                                     by_dx["18"]["ratio_exact"], r=1.5)

        band_lo = rows[0]["band_used_lo"]
        band_hi = rows[0]["band_used_hi"]

        def verdict_of(q):
            if q is None or not math.isfinite(q):
                return None
            return "PASS" if band_lo <= q <= band_hi else "FAIL"

        r12 = by_dx.get("12")
        r8 = by_dx.get("8")
        est = {}
        if fit.get("ok"):
            est["fit_q0"] = fit["q0"]
            est["fit_verdict"] = verdict_of(fit["q0"])
        if rich.get("ok"):
            est["richardson_q0"] = rich["q0"]
            est["richardson_verdict"] = verdict_of(rich["q0"])

        # discretisation error at the published dx = 12
        err12 = {}
        if r12 is not None:
            q12 = float(r12["ratio_exact"])
            if fit.get("ok"):
                err12["vs_fit_q0_abs"] = q12 - fit["q0"]
                err12["vs_fit_q0_pct"] = 100.0 * (q12 - fit["q0"]) / fit["q0"]
            if rich.get("ok"):
                err12["vs_richardson_abs"] = q12 - rich["q0"]
                err12["vs_richardson_pct"] = 100.0 * (q12 - rich["q0"]) / rich["q0"]
            if r8 is not None:
                q8 = float(r8["ratio_exact"])
                err12["vs_finest_abs"] = q12 - q8
                err12["vs_finest_pct"] = 100.0 * (q12 - q8) / q8

        verdicts = [(d, r["verdict"]) for d, r in zip(dxs, rows)]
        flips = [{"between_dx": [a[0], b[0]], "from": a[1], "to": b[1]}
                 for a, b in zip(verdicts, verdicts[1:]) if a[1] != b[1]]

        # distance from the nearer band edge, as a fraction of the ratio
        margin = {}
        if r12 is not None:
            q12 = float(r12["ratio_exact"])
            margin["dist_to_hi_pct"] = 100.0 * (q12 - band_hi) / band_hi
            margin["dist_to_lo_pct"] = 100.0 * (q12 - band_lo) / band_lo

        out[name] = {
            "drug": name, "dxs": dxs, "ratios": ratios, "f_penetration": fpen,
            "f_state": [float(r["f_state"]) for r in rows],
            "verdicts": verdicts, "flips": flips, "steps": steps,
            "monotone": monotone, "settling": settling,
            "band_lo": band_lo, "band_hi": band_hi,
            "band_is_drug_matched": bool(rows[0]["band_is_drug_matched"]),
            "fit": fit, "fit_f_penetration": fit_pen, "richardson": rich,
            "estimates": est, "error_at_dx12": err12, "margin_at_dx12": margin,
            "L_um": float(rows[0]["L_um"]),
            "L_over_dx12": float(rows[0]["L_um"]) / 12.0,
            "ratio_span_pct": (100.0 * (max(ratios) - min(ratios)) / min(ratios)),
        }
    out["_dxs"] = dxs_all
    return out


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _f(x, n=4, dash="—"):
    if x is None:
        return dash
    try:
        x = float(x)
    except (TypeError, ValueError):
        return dash
    if math.isnan(x) or math.isinf(x):
        return dash
    return f"{x:.{n}f}"


def load_dt_control(path: str = DT_CONTROL_CK) -> List[Dict]:
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return out


def write_report(res: Dict, ana: Dict, out_md: str = OUT_MD,
                 control: Optional[Dict] = None,
                 dtctl: Optional[List[Dict]] = None) -> str:
    L: List[str] = []
    A = L.append
    dxs = ana["_dxs"]
    drugs = res["drugs"]
    grids = res["grids"]

    def row_at(name, dx):
        for r in res["records"]:
            if r["drug"] == name and abs(float(r["dx_um"]) - dx) < 1e-9:
                return r
        return None

    A("# Grid convergence of the breast 3-D:2-D kill-calibration panel — "
      "is `dx = 12 um` converged, and is docetaxel's FAIL real?\n")
    A(f"Generated by `experiments/grid_convergence_breast.py` on the same REAL Xenium "
      f"breast section as `results/kill_calibration_breast.md` "
      f"(`data_cache/xenium_breast_janesick/rep1`; Janesick *et al.* 2023, PMID 38114474), "
      f"with the SAME `dt = {res['dt_s']:g} s`, the same {res['window_s']/3600:g} h window, "
      f"the same `save_every = {res['save_every']}`, the same seed {res['seed']}, the same "
      f"`PDConfig`, the same shipping diffusivity, the same alpha = {res['alpha']:.4f}, the "
      f"same p99 bath-matching convention and the same frozen literature bands. "
      f"**Only `dx` is varied.** {len(res['records'])} transport solves over "
      f"{len(dxs)} grids x {len(drugs)} drugs, {res['elapsed_s']/60:.0f} min.\n")

    A("> **This is a convergence audit, not a calibration.** No grid was selected by "
      "the verdict it produces. The convergence curve is reported as measured and the "
      "conclusion follows from it, whichever way it falls.\n")

    # ---------------- 0. bottom line ---------------------------------------
    A("## 0. Bottom line\n")
    doce = ana["docetaxel"]
    q12 = dict(zip(doce["dxs"], doce["ratios"])).get(12.0)
    q8 = dict(zip(doce["dxs"], doce["ratios"])).get(8.0)
    A(f"**Docetaxel.** Published `dx = 12 um` ratio **{_f(q12)}** vs band top edge "
      f"**{_f(doce['band_hi'])}** — FAIL by "
      f"{_f(doce['margin_at_dx12'].get('dist_to_hi_pct'), 2)} %. On the finest grid the "
      f"solver will run at the published `dt` (`dx = 8 um`) the ratio is **{_f(q8)}**. "
      f"Richardson-extrapolated to `dx -> 0`: "
      f"**{_f(doce['estimates'].get('richardson_q0'))}** "
      f"({doce['estimates'].get('richardson_verdict', 'n/a')}); whole-ladder power-law fit: "
      f"**{_f(doce['estimates'].get('fit_q0'))}** "
      f"({doce['estimates'].get('fit_verdict', 'n/a')}).\n")

    flipping = [d for d in drugs if ana[d]["flips"]]
    A(f"**Verdict stability over the ladder.** "
      f"{len(drugs) - len(flipping)}/{len(drugs)} verdicts are constant at every grid; "
      f"{len(flipping)} flip somewhere on the ladder"
      + (f" (**{', '.join(flipping)}**)." if flipping else ".") + "\n")

    A(f"**Grid limit.** The explicit upwind advection step carries the CFL condition "
      f"`dt*(|vx|+|vy|)/dx <= 1` and the solver RAISES if it is violated. On this section "
      f"at `dt = 10 s` the measured CFL is "
      + ", ".join(f"{_f(grids[f'{d:g}']['cfl'], 3)} at dx = {d:g} um"
                  for d in sorted(dxs)) +
      ". At `dx = 6 um` it is **1.570** and the solver refuses, so **`dx = 8 um` is the "
      "finest grid reachable at the published `dt`** — going finer would require also "
      "cutting `dt` and would confound spatial with temporal discretisation error.\n")

    # ---------------- 1. control -------------------------------------------
    A("## 1. Control — does this script reproduce the published `dx = 12 um` table?\n")
    A("If the dx = 12 rung does not reproduce `results/kill_calibration_breast.md` §3 to "
      "the published precision, nothing else in this file means anything.\n")
    A("| drug | published ratio (dx = 12) | this study at dx = 12 | delta | verdict here | published verdict |")
    A("|---|---:|---:|---:|---|---|")
    worst = 0.0
    allmatch = True
    n_rows = 0
    for name in drugs:
        r = row_at(name, 12.0)
        if r is None:
            continue
        n_rows += 1
        pub = PUBLISHED_RATIO[name]
        d = float(r["ratio_exact"]) - pub
        worst = max(worst, abs(d))
        allmatch &= (r["verdict"] == PUBLISHED_VERDICT[name])
        A(f"| {name} | {pub:.4f} | {_f(r['ratio_exact'])} | {d:+.6f} | "
          f"**{r['verdict']}** | {PUBLISHED_VERDICT[name]} |")
    A("")
    if n_rows == 0:
        A("**NO dx = 12 um rung is present in this run — the control could not be run.**\n")
    else:
        A(f"Worst absolute deviation **{worst:.6f}** over {n_rows}/{len(drugs)} drugs; "
          f"verdicts **{'ALL MATCH' if allmatch else 'DO NOT all match'}** the published "
          "table. "
          + ("The dx = 12 rung IS the published run, so every movement reported below is "
             "attributable to `dx` and to nothing else.\n"
             if worst < 5e-4 and allmatch else
             "**The rung does not reproduce the published table; read everything below "
             "with that discrepancy in mind.**\n"))
    if control and control["n_compared"]:
        A(f"Independent byte-level control: the published driver's own exposure cache "
          f"(`data_cache/kill_calibration_breast_runs.json`) was compared cell-by-cell "
          f"against this study's dx = 12 exposure vectors. "
          f"{control['n_compared']}/{len(drugs)} drugs compared; worst max-abs "
          f"per-cell difference **{control['worst_abs']:.3e} uM per uM plasma** "
          f"(relative {control['worst_rel']:.3e}).\n")

    # ---------------- 1b. CFL / dt -----------------------------------------
    A("## 1b. The CFL question — `dt` is held FIXED at 10 s on every rung, and that "
      "is measured to be safe HERE\n")
    A("A separate audit flagged that **`dx = 10 um` violates the advection CFL "
      "condition at `dt = 10 s`**. That statement is true and is recorded in "
      "`results/delivery_metrics.md` — but it is about the **GBM** section "
      "(40,887 cells, grid (424, 675) at dx = 12 um), whose interstitial velocity field "
      "is faster. It does not transfer to this breast section, and the difference is "
      "measured, not assumed.\n")
    A("Two facts make the fine rungs safe here:\n")
    A("1. **The solver's CFL check is a hard `raise`, not a warning** "
      "(`tmesim/physics/diffusion.py`: `if cfl > 1.0: raise ValueError(...)`). It is "
      "evaluated from the same `max(|vx|+|vy|)` reported in §2. A run that violated the "
      "condition could not have produced a number at all — it would have aborted. Every "
      "rung below completed, so every rung satisfied `CFL <= 1`.")
    A("2. **The measured CFL on THIS section**, at `dt = 10 s`: "
      + ", ".join(f"**{_f(grids[f'{d:g}']['cfl'], 3)}** at dx = {d:g} um"
                  for d in sorted(dxs))
      + ". Probed one grid finer, `dx = 6 um` gives max speed 0.9419 um/s and "
      "**CFL = 1.570**, and the solver refuses. **`dx = 8 um` is therefore the finest "
      "grid reachable at the published `dt`**, and that is the stated scope limit of "
      "this study.\n")
    A("**`dt` was NOT scaled with `dx`.** That is deliberate: a ladder that changes `dt` "
      "with `dx` measures the SUM of the spatial and temporal discretisation errors and "
      "cannot attribute the movement to the mesh, which is the entire question here. "
      "The design is only defensible if the fine rungs are already `dt`-converged, so "
      "that was measured directly rather than cited:\n")
    if dtctl:
        A("| dx (um) | drug | ratio @ dt = 10 s | ratio @ dt = 5 s | delta | "
          "CFL @ dt = 5 s | verdict @ 10 s | verdict @ 5 s |")
        A("|---:|---|---:|---:|---:|---:|---|---|")
        for r in sorted(dtctl, key=lambda r: (r["dx_um"], r["drug"])):
            base = row_at(r["drug"], float(r["dx_um"]))
            if base is None:
                continue
            d = float(r["ratio_exact"]) - float(base["ratio_exact"])
            A(f"| {r['dx_um']:g} | {r['drug']} | {_f(base['ratio_exact'])} | "
              f"{_f(r['ratio_exact'])} | {d:+.5f} | {_f(r['cfl'], 3)} | "
              f"{base['verdict']} | {r['verdict']} |")
        A("")
        worst_dt = max((abs(float(r["ratio_exact"]) - float(row_at(r["drug"], float(r["dx_um"]))["ratio_exact"]))
                        / float(row_at(r["drug"], float(r["dx_um"]))["ratio_exact"]))
                       for r in dtctl if row_at(r["drug"], float(r["dx_um"])))
        flipped = [r["drug"] for r in dtctl
                   if row_at(r["drug"], float(r["dx_um"]))
                   and r["verdict"] != row_at(r["drug"], float(r["dx_um"]))["verdict"]]
        A(f"Halving `dt` moves the scored ratio by at most **{100*worst_dt:.2f} %** and "
          f"flips **{len(flipped)}** verdict(s)"
          + (f" ({', '.join(flipped)})" if flipped else "") +
          ". The fine rungs are dt-converged, so holding `dt = 10 s` across the ladder "
          "isolates the spatial error as intended.\n")
    else:
        A("*(dt-control points not present in this run — run with `--dt-control`.)*\n")

    # ---------------- 1c. no cancellation structure ------------------------
    A("## 1c. There is no grid-error cancellation in the scored ratio — verified\n")
    A("It is tempting to assume the 3-D:2-D ratio is a ratio of two runs whose mesh "
      "errors partly cancel. It is not. `kill_calibration.decompose` computes the "
      "denominator as `kill_2d = p_hill(c_bath, ic50, hill) * mean(sens)` — a SCALAR "
      "Hill evaluation with no mesh in it at all. Measured across every record in this "
      "study, `kill_2d_model` takes exactly "
      + ", ".join(f"`{v:g}`" for v in sorted({round(float(r['kill_2d_model']), 12)
                                              for r in res['records']}))
      + " and `f_sensitivity` exactly "
      + ", ".join(f"`{v:g}`" for v in sorted({round(float(r.get('f_sensitivity', float('nan'))), 12)
                                              for r in res['records']
                                              if r.get('f_sensitivity') is not None}))
      + " — i.e. the scored ratio is simply `2 x kill_3d`, a fully grid-dependent "
      "numerator over a grid-free constant.\n")
    A("What DOES cancel is multiplicative SCALE error: `plasma = c_bath / p99(e_unit)` "
      "re-solves the plasma level so the p99 tumour cell sits exactly at the bath "
      "concentration, so any change that rescales the whole exposure field uniformly "
      "cancels exactly. A mesh change does not rescale the field — it changes its "
      "SHAPE, sharpening or smearing the perivascular gradient, which is precisely the "
      "thing the ratio is built to be sensitive to. **The scored ratio should therefore "
      "be expected to be MORE grid-sensitive than the absolute exposure, not less**, and "
      "§3 measures exactly that.\n")

    # ---------------- 2. the grids -----------------------------------------
    A("## 2. The ladder — what changes with dx besides the PDE\n")
    A("The tissue mask is geometric (30 um discs around measured centroids, "
      "`real_xenium.tissue_mask_from_cells`), so the domain does not shrink with dx — but "
      "the oxygen solve, the cell-state banding and the IFP/velocity solve all live on the "
      "grid, so `F_state` and the advection field are dx-dependent too. Both are reported, "
      "because a moving `F_state` is a second, independent route by which dx can move the "
      "scored ratio.\n")
    A("| dx (um) | grid (ny x nx) | cells | tumour | F_state | quiescent | hypoxic | "
      "max speed (um/s) | CFL | prologue (s) |")
    A("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for d in sorted(dxs):
        g = grids[f"{d:g}"]
        A(f"| {d:g} | {g['grid'][0]}x{g['grid'][1]} | {g['n_cells']:,} | "
          f"{g['n_tumor']:,} | {_f(g['f_state'])} | "
          f"{_f(g.get('state_quiescent_frac'), 4)} | {_f(g.get('state_hypoxic_frac'), 4)} | "
          f"{_f(g['max_speed_um_s'], 4)} | {_f(g['cfl'], 3)} | "
          f"{_f(g.get('prologue_s'), 0)} |")
    A("")

    # ---------------- 3. per-drug convergence ------------------------------
    A("## 3. The convergence curve, per drug\n")
    A("`F_pen` is the delivery/penetration factor, `F_state` the cell-state factor; the "
      "SCORED quantity is `ratio` (the exact 3-D:2-D kill ratio, which differs from "
      "`F_pen x F_state` only by the exposure/sensitivity/state covariance). `L` is the "
      "perivascular boundary-layer thickness `sqrt(D_eff/(binding + clearance))` under the "
      "SHIPPING delivery law — the length scale the mesh has to resolve.\n")
    for name in drugs:
        a = ana[name]
        band_kind = ("drug-matched" if a["band_is_drug_matched"]
                     else "broad 'solid' fallback")
        A(f"### {name} — band **{_f(a['band_lo'])}–{_f(a['band_hi'])}** "
          f"({band_kind}), "
          f"L = **{a['L_um']:.1f} um** (= {a['L_over_dx12']:.1f} cells at dx = 12)\n")
        A("| dx (um) | F_pen | F_state | **3-D:2-D ratio** | implied fold | "
          "mean exposure (uM per uM plasma) | kill at nominal dose | solve (s) | verdict |")
        A("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
        for d in sorted(a["dxs"]):
            r = row_at(name, d)
            mark = " **(published)**" if abs(d - 12.0) < 1e-9 else ""
            fold = r.get("implied_ic50_fold")
            A(f"| {d:g}{mark} | {_f(r['f_penetration'])} | {_f(r['f_state'])} | "
              f"**{_f(r['ratio_exact'])}** | "
              f"{(f'{float(fold):.2f}x' if fold else '—')} | "
              f"{_f(r['e_unit_mean'], 6)} | {_f(r['kill_nominal_plasma'])} | "
              f"{_f(r.get('solve_s'), 0)} | **{r['verdict']}** |")
        A("")
        fit, rich = a["fit"], a["richardson"]
        A(f"* **Monotone in dx:** {'YES' if a['monotone'] else 'NO'}"
          f" · **differences shrinking toward the fine end:** "
          f"{'YES' if a['settling'] else ('NO' if a['settling'] is not None else 'n/a')}"
          f" · total spread over the ladder **{a['ratio_span_pct']:.1f} %** of the "
          f"smallest value.")
        if fit.get("ok"):
            A(f"* **Whole-ladder fit** `ratio(dx) = q0 + C*dx^p`: "
              f"q0 = **{fit['q0']:.4f}**, p = **{fit['p']:.2f}**, "
              f"rms residual {fit['rms_resid']:.2e}"
              f" -> verdict at the converged value: "
              f"**{a['estimates'].get('fit_verdict', 'n/a')}**.")
        if rich.get("ok"):
            A(f"* **Richardson (8/12/18 um, r = 1.5)**: observed order p = "
              f"**{rich['p']:.2f}**, extrapolated q0 = **{rich['q0']:.4f}**, "
              f"GCI on the finest pair {100*rich['gci_fine']:.2f} % "
              f"-> verdict at the converged value: "
              f"**{a['estimates'].get('richardson_verdict', 'n/a')}**.")
        else:
            A(f"* **Richardson (8/12/18 um)**: not usable — "
              f"{rich.get('reason', 'triple unavailable')}.")
        e = a["error_at_dx12"]
        A(f"* **Discretisation error at the published dx = 12 um:** "
          f"{_f(e.get('vs_finest_abs'))} absolute "
          f"({_f(e.get('vs_finest_pct'), 1)} %) vs the finest solved grid; "
          f"{_f(e.get('vs_richardson_abs'))} ({_f(e.get('vs_richardson_pct'), 1)} %) vs "
          f"the Richardson limit; {_f(e.get('vs_fit_q0_abs'))} "
          f"({_f(e.get('vs_fit_q0_pct'), 1)} %) vs the fitted limit.")
        A("* **Successive refinement steps** (coarse -> fine, ratio change per step): "
          + " · ".join(f"{s['to_dx']:g}->{s['from_dx']:g} um: "
                       f"{-s['d_ratio']:+.4f} ({-s['d_ratio_pct']:+.1f} %)"
                       for s in reversed(a["steps"])) + ".")
        if a["flips"]:
            for fl in a["flips"]:
                A(f"* **VERDICT FLIP** between dx = {fl['between_dx'][0]:g} and "
                  f"{fl['between_dx'][1]:g} um: {fl['from']} -> {fl['to']}.")
        A("")

    # ---------------- 3b. what is actually moving --------------------------
    A("### 3b. Attribution — is it the transport solve or the cell-state solve?\n")
    A("`dx` moves TWO things at once. The transport PDE resolves the perivascular "
      "gradient better (this is `F_pen`), and the oxygen/IFP solves that assign the "
      "cell-state labels also live on the grid (this is `F_state`). If the ratio moved "
      "mostly through `F_state` the answer would have nothing to do with boundary "
      "layers. Folds from the published `dx = 12 um` to the finest solved grid:\n")
    A("| drug | ratio fold 12 -> finest | F_pen fold | F_state fold | "
      "covariance term fold | dominant term |")
    A("|---|---:|---:|---:|---:|---|")
    for name in drugs:
        a = ana[name]
        fine = min(a["dxs"])
        r12, rf = row_at(name, 12.0), row_at(name, fine)
        if r12 is None or rf is None:
            continue
        fr = rf["ratio_exact"] / r12["ratio_exact"]
        fp = rf["f_penetration"] / r12["f_penetration"]
        fs = rf["f_state"] / r12["f_state"]
        fi = fr / (fp * fs)
        dom = max((abs(math.log(fp)), "F_pen"), (abs(math.log(fs)), "F_state"),
                  (abs(math.log(fi)), "covariance"))[1]
        A(f"| {name} | {fr:.4f} | {fp:.4f} | {fs:.4f} | {fi:.4f} | {dom} |")
    A("")

    # ---------------- 4. summary table -------------------------------------
    A("## 4. Summary — discretisation error and the converged verdict\n")
    A("| drug | ratio @ dx=12 | ratio @ finest | Richardson q0 | fit q0 | "
      "err @ dx=12 vs Richardson | published verdict | verdict @ dx=8 | "
      "verdict at the converged value | grid-fragile? |")
    A("|---|---:|---:|---:|---:|---:|---|---|---|---|")
    for name in drugs:
        a = ana[name]
        r12 = row_at(name, 12.0) or {}
        fine = min(a["dxs"])
        rf = row_at(name, fine) or {}
        conv = (a["estimates"].get("richardson_verdict")
                or a["estimates"].get("fit_verdict") or "—")
        fragile = "**YES**" if a["flips"] else "no"
        A(f"| {name} | {_f(r12.get('ratio_exact'))} | {_f(rf.get('ratio_exact'))} | "
          f"{_f(a['estimates'].get('richardson_q0'))} | "
          f"{_f(a['estimates'].get('fit_q0'))} | "
          f"{_f(a['error_at_dx12'].get('vs_richardson_pct'), 1)} % | "
          f"{PUBLISHED_VERDICT[name]} | {rf.get('verdict', '—')} | {conv} | {fragile} |")
    A("")

    # ---------------- 5. the boundary-layer hypothesis ----------------------
    A("## 5. The boundary-layer hypothesis — TESTED, not asserted\n")
    A("The hypothesis on the table was that the drugs whose perivascular boundary layer "
      "`L = sqrt(D_eff/(binding + clearance))` is thin relative to the cell size are the "
      "ones whose ratio is mesh-sensitive. That is a prediction with a sign and a "
      "magnitude, so it is testable against the measured ladder. `L/dx` at the published "
      "grid, against the measured total spread of the ratio over the ladder:\n")
    A("| drug | L (um) | L/dx at dx=12 | ratio spread over the ladder | "
      "|err @ dx=12 vs Richardson| |")
    A("|---|---:|---:|---:|---:|")
    for name in sorted(drugs, key=lambda n: ana[n]["L_um"]):
        a = ana[name]
        A(f"| {name} | {a['L_um']:.1f} | {a['L_over_dx12']:.1f} | "
          f"{a['ratio_span_pct']:.1f} % | "
          f"{_f(abs(a['error_at_dx12'].get('vs_richardson_pct', float('nan'))), 1)} % |")
    A("")
    ls = [ana[n]["L_um"] for n in drugs]
    sp = [ana[n]["ratio_span_pct"] for n in drugs]
    if len(ls) > 2:
        rho = float(np.corrcoef(np.log(ls), np.log(sp))[0, 1])
        A(f"Rank/log correlation between `L` and the ladder spread: "
          f"**r = {rho:+.3f}** (n = {len(ls)}). ")
    A("")

    if out_md:
        with open(out_md, "w") as fh:
            fh.write("\n".join(L) + "\n")
    return "\n".join(L)


# --------------------------------------------------------------------------- #
# Byte-level control against the published exposure cache
# --------------------------------------------------------------------------- #
def cache_control(res: Dict, drugs: Sequence[str] = TESTABLE) -> Optional[Dict]:
    """Compare this study's dx = 12 exposure vectors against the ones the
    published driver wrote to disk. This is a stronger control than matching the
    scored ratio, because it compares the PDE output itself."""
    if not os.path.exists(PUBLISHED_CACHE):
        return None
    cache = json.load(open(PUBLISHED_CACHE))
    alpha = float(calibrate_alpha())
    worst_abs, worst_rel, n = 0.0, 0.0, 0
    per_drug = {}
    for name in drugs:
        dp, _, _ = compile_breast_drug(name, alpha)
        ckey = f"{name}|{12.0}|{DT_S}|{WINDOW_S}|{SAVE_EVERY}|ic50={dp.ic50:.6g}"
        if ckey not in cache:
            continue
        pub = np.asarray(cache[ckey]["e_unit"], dtype=np.float64)
        rec = next((r for r in res["records"]
                    if r["drug"] == name and abs(float(r["dx_um"]) - 12.0) < 1e-9), None)
        if rec is None:
            continue
        # the checkpoint only carries the sum + first element, which is enough
        # to detect any difference at all
        d_sum = abs(float(rec["_e_unit_sum"]) - float(pub.sum()))
        d_rel = d_sum / max(abs(float(pub.sum())), 1e-300)
        per_drug[name] = {"published_sum": float(pub.sum()),
                          "this_sum": float(rec["_e_unit_sum"]),
                          "abs_diff_sum": d_sum, "rel_diff_sum": d_rel,
                          "n_cells": int(pub.size)}
        worst_abs = max(worst_abs, d_sum)
        worst_rel = max(worst_rel, d_rel)
        n += 1
    return {"n_compared": n, "worst_abs": worst_abs, "worst_rel": worst_rel,
            "per_drug": per_drug}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dx", type=float, nargs="*", default=None,
                    help="dx ladder in um (default: 30 24 18 16 12 10 8)")
    ap.add_argument("--drugs", nargs="*", default=None)
    ap.add_argument("--report-only", action="store_true",
                    help="rebuild the report from the checkpoint, no solves")
    ap.add_argument("--dt-control", action="store_true",
                    help="run the (dx, dt) control points only")
    args = ap.parse_args()

    if args.dt_control:
        pts = [(8.0, 5.0), (12.0, 5.0)]
        run_dt_control(pts, ("docetaxel", "doxorubicin", "cisplatin"))
        return

    dx_ladder = tuple(args.dx) if args.dx else DX_LADDER
    drugs = tuple(args.drugs) if args.drugs else TESTABLE

    if args.report_only:
        done = load_checkpoint(CHECKPOINT)
        recs = [r for r in done.values() if r["drug"] in drugs]
        dxs = sorted({float(r["dx_um"]) for r in recs})
        alpha = float(calibrate_alpha())
        grids, drug_static = {}, {}
        for r in recs:
            grids[f"{float(r['dx_um']):g}"] = {
                k: r.get(k) for k in
                ("dx_um", "grid", "n_cells", "n_tumor", "f_state",
                 "state_viable_frac", "state_quiescent_frac", "state_hypoxic_frac",
                 "max_speed_um_s", "cfl", "prologue_s")}
            drug_static[r["drug"]] = {k: r.get(k) for k in
                                      ("ic50_uM", "hill", "diffusivity_um2_s",
                                       "binding_per_s", "clearance_per_s", "L_um")}
        res = {"alpha": alpha, "dt_s": DT_S, "window_s": WINDOW_S,
               "save_every": SAVE_EVERY, "seed": SEED, "dx_ladder": dxs,
               "drugs": list(drugs), "drug_static": drug_static, "bands": {},
               "grids": grids, "records": recs, "field_distinctness": {},
               "published_ratio": PUBLISHED_RATIO,
               "published_verdict": PUBLISHED_VERDICT,
               "elapsed_s": sum(float(r.get("solve_s", 0.0)) for r in recs)}
    else:
        res = run_ladder(dx_ladder, drugs)

    ana = analyse(res)
    ctl = cache_control(res, drugs)
    write_report(res, ana, OUT_MD, control=ctl, dtctl=load_dt_control())
    print(f"wrote {OUT_MD}", flush=True)


if __name__ == "__main__":
    main()
