"""Mechanism attribution under BOTH transvascular delivery laws, side by side.

Why this script exists
======================
``results/mechanism_attribution.md`` §3 was measured under the historical
ONE-WAY vascular source, ``S = PS*Cp``, which has no path back from the
interstitium to the blood and therefore lets tissue concentration exceed the
plasma level. ``PDEConfig(transvascular_exchange=True)`` replaces it with
gradient-driven exchange, ``S = PS*(Cp - C)``, whose steady state
``PS*Cp/(PS + clearance + binding)`` is bounded by ``Cp`` by construction.

Which law SHIPS is an open decision; the flag currently defaults to False. So
the report has to carry both arms as first-class, and that needs both arms
MEASURED rather than one measured and one inferred. This script produces the
full §3 table for each drug under each law from one code path, in one process,
with everything but the flag held identical.

The propagation trap
====================
``PDEConfig(transvascular_exchange=True)`` is a silent no-op on any path that
calls ``simulate_drug_transport`` directly (see
``tmesim/validation/propagation.py``). This script goes through
``sim.engine.simulate_response``, which is the one declared PROPAGATED path,
and it does not take that on trust: every run is wrapped in
``propagation.instrument()`` so the solver boundary is observed, and the ON arm
is required to have been handed ``exchange_k`` on every solver call while the
OFF arm was handed it on none. The two arms' AUC fields are additionally
required to differ bit-wise before any comparison is reported.

Run::

    ./.venv/bin/python experiments/mechanism_attribution_arms.py

Writes a JSON summary to stdout and to the path given by ``--out`` (default:
none, stdout only). It writes nothing into ``results/`` itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from typing import Dict, List

import numpy as np

# --- production constants, copied from app/server/jobs.py so this script is
#     the same run the API would have executed -------------------------------
PDE_DT_S = 10.0
PDE_SAVE_EVERY = 500
PDE_CONVECTION = True
SIM_SEED = 0
SECTION = "xenium"
DX_UM = 12.0
TOTAL_TIME_S = 21600.0
PLASMA_PEAK_UM = 50.0
HALF_LIFE_S = 14400.0
N_DOSES = 1
DOSE_INTERVAL_S = 86400.0
SENS_MODE = "uniform"
TIME_RESOLVED = True
TIMELAPSE_FRAMES = 24
DRUGS = ("temozolomide", "cisplatin")


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def _digest(arr: np.ndarray) -> str:
    return hashlib.blake2b(
        np.ascontiguousarray(arr, dtype=np.float64).tobytes(), digest_size=16
    ).hexdigest()


def run_arm(micro, entry_sensitivity, drug_name: str, exchange: bool) -> Dict:
    """One production-equivalent run. ``exchange`` picks the delivery law."""
    from tmesim.core import PDEConfig, PKConfig, SimConfig
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.library import get_molecule
    from tmesim.reconstruct.real_pipeline import apply_drug_sensitivity
    from tmesim.sim.engine import plasma_profile_safe, simulate_response
    from tmesim.validation import delivery_metrics, propagation
    from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

    # Reset per-cell sensitivity exactly as the job runner does, so one arm can
    # never inherit the other's.
    micro.cells.sensitivity = np.array(entry_sensitivity, copy=True)
    apply_drug_sensitivity(micro, drug_name, mode=SENS_MODE)

    alpha = float(calibrate_alpha())
    dp = apply_barrier_to_drugparams(
        calibrated_drugparams(drug_name, alpha), get_molecule(drug_name)
    )

    pde = PDEConfig(
        dt=PDE_DT_S,
        total_time=TOTAL_TIME_S,
        convection=PDE_CONVECTION,
        save_every=PDE_SAVE_EVERY,
        time_resolved=TIME_RESOLVED,
        transvascular_exchange=exchange,
    )
    pk_cfg = PKConfig(
        plasma_peak=PLASMA_PEAK_UM,
        half_life=HALF_LIFE_S,
        n_doses=N_DOSES,
        dose_interval=DOSE_INTERVAL_S,
    )
    cfg = SimConfig(drug=dp, pk=pk_cfg, pde=pde, seed=SIM_SEED)

    _log(f"  run {drug_name} exchange={exchange} ...")
    t0 = time.time()
    with propagation.instrument() as probe:
        res = simulate_response(
            micro, cfg, attribute=True, capture_frames=TIMELAPSE_FRAMES
        )
    dt_run = time.time() - t0
    _log(f"  done in {dt_run:.1f}s  kill_fraction={res.kill_fraction!r}")

    dm = delivery_metrics.compute_all(micro, res, dp)
    head = dm["headline"]
    mech = res.mechanism
    counts = mech["counts"]
    sb = mech["survival_breakdown"]
    rec = mech["reconciliation"]

    # --- tissue / plasma exposure ratio ---------------------------------- #
    times = np.arange(0.0, TOTAL_TIME_S, PDE_DT_S)
    plasma = plasma_profile_safe(cfg, times)
    plasma_auc = float(plasma.sum() * PDE_DT_S)
    auc = np.asarray(res.fields["auc"], dtype=np.float64)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    auc_tissue = auc[tissue]

    # cell-wise, which is what results/transvascular_exchange_revalidation.md §2
    # quotes ("max cellwise AUC/plasma-AUC"). domain.to_index is the same
    # nearest-grid-cell map the pharmacodynamic kill uses.
    col, row = micro.domain.to_index(micro.cells.x, micro.cells.y)
    auc_cells = auc[row, col]
    tumor_mask = micro.cells.tumor_mask()
    auc_tumor_cells = auc[row[tumor_mask], col[tumor_mask]]

    frames = res.transport.frames
    peak_by_frame = (
        [float(f.max()) for f in np.asarray(frames)] if frames is not None else None
    )

    return {
        "drug": drug_name,
        "transvascular_exchange": bool(exchange),
        "runtime_s": dt_run,
        # --- propagation evidence ---------------------------------------- #
        "solver_calls": probe.solver_calls,
        "solver_calls_with_exchange_k": probe.exchange_k_calls,
        "exchange_coefficient_built": probe.exchange_coefficient_calls,
        "auc_field_digest": _digest(auc),
        # --- headline ----------------------------------------------------- #
        "ic50_uM": float(dp.ic50),
        "bbb_penetrance": float(dp.bbb_penetrance),
        "kill_fraction": float(res.kill_fraction),
        "kill_fraction_repr": repr(float(res.kill_fraction)),
        "mean_tumor_auc_uM_s": head["mean_tumor_auc_uM_s"],
        "therapeutic_coverage_1x": head["therapeutic_coverage_1x_ic50"],
        "penetration_index": float(res.extras["penetration_index"]),
        # --- exposure vs plasma ------------------------------------------- #
        "plasma_auc_uM_s": plasma_auc,
        "mean_tissue_auc_over_plasma_auc": float(auc_tissue.mean() / plasma_auc),
        "max_tissue_auc_over_plasma_auc": float(auc_tissue.max() / plasma_auc),
        "mean_cell_auc_over_plasma_auc": float(auc_cells.mean() / plasma_auc),
        "max_cell_auc_over_plasma_auc": float(auc_cells.max() / plasma_auc),
        "mean_tumor_cell_auc_over_plasma_auc": float(auc_tumor_cells.mean() / plasma_auc),
        "max_tumor_cell_auc_over_plasma_auc": float(auc_tumor_cells.max() / plasma_auc),
        # --- mechanism ----------------------------------------------------- #
        "n_tumor": mech["n_tumor"],
        "n_survived": mech["n_survived"],
        "counts": {k: int(v) for k, v in counts.items()},
        "n_killed": int(counts["KILLED"]),
        "survival_breakdown": {k: float(v) for k, v in sb.items()},
        "unattributed": int(counts["UNATTRIBUTED"]),
        "joint_limited": {
            "count": int(mech["joint_limited"]["count"]),
            "fraction_of_exposure_limited": float(
                mech["joint_limited"]["fraction_of_exposure_limited"]),
        },
        "exposure_ratio_percentiles_tumor": {
            k: (float(v) if v is not None else None)
            for k, v in mech["exposure_ratio_percentiles_tumor"].items()
        },
        "exposure_shortfall_fold": (
            float(mech["exposure_shortfall_fold"])
            if mech["exposure_shortfall_fold"] is not None else None),
        "thresholds": {k: (float(v) if isinstance(v, (int, float)) else v)
                       for k, v in mech["thresholds"].items()},
        # --- reconciliation (the correctness property, law-independent) ---- #
        "reconciliation": {
            "bit_exact": bool(rec["bit_exact"]),
            "kill_fraction_from_labels": repr(float(rec["kill_fraction_from_labels"])),
            "kill_fraction_from_engine": repr(float(rec["kill_fraction_from_engine"])),
            "abs_difference": float(rec["abs_difference"]),
        },
        "timelapse_peak_concentration_uM": peak_by_frame,
        "timelapse_frame_times_s": (
            [float(t) for t in np.asarray(res.transport.frame_times)]
            if res.transport.frame_times is not None else None),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="write JSON summary here")
    args = ap.parse_args()

    from tmesim.reconstruct.real_pipeline import reconstruct_real

    _log(f"reconstruct_real({SECTION!r}, dx_um={DX_UM}) ...")
    t0 = time.time()
    micro = reconstruct_real(SECTION, dx_um=DX_UM)
    _log(f"reconstruction done in {time.time() - t0:.1f}s  "
         f"grid={micro.domain.shape}  cells={len(micro.cells.x)}")
    pristine_sensitivity = np.array(micro.cells.sensitivity, copy=True)

    rows: List[Dict] = []
    for drug in DRUGS:
        for exchange in (False, True):
            rows.append(run_arm(micro, pristine_sensitivity, drug, exchange))

    # ---- the guard: the two arms must not be the same run ----------------- #
    checks: Dict[str, Dict] = {}
    for drug in DRUGS:
        off = next(r for r in rows if r["drug"] == drug and not r["transvascular_exchange"])
        on = next(r for r in rows if r["drug"] == drug and r["transvascular_exchange"])
        checks[drug] = {
            "off_solver_calls": off["solver_calls"],
            "off_solver_calls_with_exchange_k": off["solver_calls_with_exchange_k"],
            "on_solver_calls": on["solver_calls"],
            "on_solver_calls_with_exchange_k": on["solver_calls_with_exchange_k"],
            "auc_fields_differ": off["auc_field_digest"] != on["auc_field_digest"],
            "off_digest": off["auc_field_digest"],
            "on_digest": on["auc_field_digest"],
        }
        assert off["solver_calls"] >= 1 and on["solver_calls"] >= 1, (
            f"{drug}: the transport solver never ran; there is nothing to compare")
        assert off["solver_calls_with_exchange_k"] == 0, (
            f"{drug}: the OFF arm was handed exchange_k, which it must not be")
        assert on["solver_calls_with_exchange_k"] == on["solver_calls"], (
            f"{drug}: the ON arm did NOT receive exchange_k on every solver "
            f"call ({on['solver_calls_with_exchange_k']}/{on['solver_calls']}). "
            f"The flag is INERT on this path and the comparison would be a run "
            f"against itself.")
        assert off["auc_field_digest"] != on["auc_field_digest"], (
            f"{drug}: the two arms produced BIT-IDENTICAL AUC fields")
        assert off["reconciliation"]["bit_exact"] and on["reconciliation"]["bit_exact"], (
            f"{drug}: attribution failed to reconcile with kill_fraction")

    out = {
        "spec": {
            "section": "xenium_319",
            "dx_um": DX_UM,
            "sensitivity_mode": SENS_MODE,
            "vascular_mode": "discrete",
            "oxygen_mode": "pde",
            "state_mode": "oxygen",
            "pk": {"plasma_peak": PLASMA_PEAK_UM, "half_life_s": HALF_LIFE_S,
                   "n_doses": N_DOSES, "dose_interval_s": DOSE_INTERVAL_S},
            "total_time_s": TOTAL_TIME_S,
            "dt_s": PDE_DT_S,
            "convection": PDE_CONVECTION,
            "time_resolved": TIME_RESOLVED,
            "timelapse_frames": TIMELAPSE_FRAMES,
            "seed": SIM_SEED,
        },
        "propagation_checks": checks,
        "runs": rows,
    }
    text = json.dumps(out, indent=2)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")
        _log(f"wrote {args.out}")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
