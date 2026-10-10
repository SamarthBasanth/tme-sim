"""One solver run, the SHIPPED way, with nothing persisted to the run store.

This module exists so the surrogate feasibility study's training data is
produced by exactly the function the Studio queue calls --
``app.server.jobs.execute_run`` -- and not by a reimplementation of the
transport chain. ``_StdoutJob`` is copied in shape from
``experiments/multisection_generalisation.py``: ``execute_run`` only ever calls
``ctx.log`` / ``ctx.stage`` / ``ctx.check_cancelled``.

Nothing here writes to ``results/`` or ``app/runs/``.
"""
from __future__ import annotations

import os
import sys
import threading
import time
import traceback
from typing import Dict, List

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

ALPHA_EXPECTED = 22.223477711172563


def assert_alpha(where: str) -> float:
    from tmesim.validation.calibration import calibrate_alpha

    a = float(calibrate_alpha())
    if repr(a) != repr(ALPHA_EXPECTED):
        raise SystemExit(
            f"[{where}] calibrate_alpha() = {a!r}, expected {ALPHA_EXPECTED!r}. "
            f"The engine is not the one this study claims to measure."
        )
    print(f"[alpha] {where}: {a!r} OK", flush=True)
    return a


class _StdoutJob:
    def __init__(self, run_id: str, spec: Dict, quiet: bool = True):
        self.run_id = run_id
        self.spec = spec
        self.stage = "queued"
        self.progress = 0.0
        self.cancel_event = threading.Event()
        self.logs: List[str] = []
        self._quiet = quiet

    def log(self, msg: str) -> None:
        self.logs.append(str(msg))
        if not self._quiet:
            print(f"    [{self.run_id}] {msg}", flush=True)


def make_spec(section: str, drug: str, dx_um: float, total_time_s: float,
              plasma_peak: float, half_life_s: float, n_doses: int,
              dose_interval_s: float) -> Dict:
    """The shipped Studio default RunSpec with the study's axes varying."""
    from tmesim.reconstruct.real_pipeline import (
        OXYGEN_MODE_DEFAULT,
        STATE_MODE_DEFAULT,
    )

    return {
        "section": section,
        "dx_um": float(dx_um),
        "drug": str(drug),
        "sensitivity_mode": "uniform",
        "vascular_mode": "discrete",
        "oxygen_mode": OXYGEN_MODE_DEFAULT,
        "state_mode": STATE_MODE_DEFAULT,
        "mgmt": None,
        "pk": {
            "plasma_peak": float(plasma_peak),
            "half_life_s": float(half_life_s),
            "n_doses": int(n_doses),
            "dose_interval_s": float(dose_interval_s),
        },
        "total_time_s": float(total_time_s),
        # Mechanism attribution is OFF: the run is bit-identical without it and
        # none of this study's targets are mechanism counts.
        "attribute_mechanism": False,
        "timelapse_frames": 0,
    }


def section_record(section: str) -> Dict:
    from app.server.api import _with_caveats
    from app.server.cache import lookup_section

    rec = dict(lookup_section(section) or {})
    rec.setdefault("id", section)
    try:
        return _with_caveats(rec)
    except Exception:  # noqa: BLE001
        return rec


def spatial_spread_valid(rec: Dict) -> bool:
    """Does the project declare per-cell spatial statistics valid here?

    The Visium (spot) refusal in ``app.server.api`` says the delivery map, its
    falloff and the delivery spread are NOT readings of that tissue, because
    the endothelial positions driving the source are synthesised. That does not
    stop the solver producing a number, and a surrogate's ground truth is the
    SOLVER -- so those sections are kept for the surrogate-error measurement,
    and flagged so no tissue claim is built on them.
    """
    for c in (rec.get("caveats") or []):
        if str(c.get("capability")) == "spatial_spread":
            return str(c.get("status")) == "ok"
    return rec.get("source") == "xenium"


def run_one(section: str, spec: Dict, cache, run_id: str = "x") -> Dict:
    """Execute one spec against a caller-owned reconstruction cache."""
    from app.server.jobs import _Ctx, execute_run

    job = _StdoutJob(run_id, spec)
    ctx = _Ctx(job)
    t0 = time.time()
    try:
        result = execute_run(spec, ctx, cache=cache)
    except Exception as exc:  # noqa: BLE001
        return {
            "status": "FAILED",
            "wall_clock_s": time.time() - t0,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc()[-2000:],
            "log_tail": job.logs[-20:],
        }
    wall = time.time() - t0
    result.pop("_payloads", None)
    rep = result.get("report") or {}
    met = result.get("metrics") or {}
    xm = result.get("extra_metrics") or {}
    reach = result.get("reach") or {}
    dp = result.get("drug_params") or {}
    prov = result.get("provenance") or {}
    cf = rep.get("celltype_fractions") or {}
    return {
        "status": "ok",
        "wall_clock_s": wall,
        "recon_build_s": rep.get("reconstruction_build_time_s"),
        "recon_cache_hit": rep.get("reconstruction_cache_hit"),
        # ---- cheap tissue features (available once a section is loaded) ----
        "n_cells": rep.get("n_cells"),
        "n_tumor": rep.get("n_tumor"),
        "grid": rep.get("grid"),
        "extent_um": rep.get("extent_um"),
        "frac_TUMOR": cf.get("TUMOR"),
        "frac_STROMAL": cf.get("STROMAL"),
        "frac_IMMUNE": cf.get("IMMUNE"),
        "frac_ENDOTHELIAL": cf.get("ENDOTHELIAL"),
        "frac_NECROTIC": cf.get("NECROTIC"),
        "mean_vascular_density": rep.get("mean_vascular_density"),
        "mean_ecm_density": rep.get("mean_ecm_density"),
        "mean_oxygen": rep.get("mean_oxygen"),
        "tissue_coverage_fraction": rep.get("tissue_coverage_fraction"),
        "max_darcy_speed_um_s": rep.get("max_darcy_speed_um_s"),
        # ---- compound transport descriptors (free, no solve) ----
        "drug_diffusivity_um2_s": dp.get("diffusivity_um2_s"),
        "drug_binding_1_per_s": dp.get("binding_1_per_s"),
        "drug_clearance_1_per_s": dp.get("clearance_1_per_s"),
        "drug_bbb_penetrance": dp.get("bbb_penetrance"),
        "drug_mw": dp.get("molecular_weight"),
        "drug_logp": dp.get("logp"),
        "drug_ic50_uM": dp.get("ic50_uM"),
        # ---- SOLVER OUTPUT: the candidate surrogate targets ----
        "E_p10_uM": xm.get("E_p10_uM"),
        "E_p50_uM": xm.get("E_p50_uM"),
        "E_p90_uM": xm.get("E_p90_uM"),
        "penetration_length_um": met.get("penetration_length_um"),
        "penetration_length_analytic_um": xm.get("penetration_length_analytic_um"),
        "penetration_index": xm.get("penetration_index"),
        "mean_tumor_auc": met.get("mean_tumor_auc"),
        "median_tumor_auc": xm.get("median_tumor_auc"),
        "mean_tissue_auc": xm.get("mean_tissue_auc"),
        "plasma_mean_uM": xm.get("plasma_mean_uM"),
        "plasma_peak_uM": xm.get("plasma_peak_uM"),
        "low_o2_over_high_o2_exposure": xm.get("low_o2_over_high_o2_exposure"),
        "reach_null_reason": reach.get("null_reason"),
        "reach_headline_x": reach.get("headline_x"),
        "reach_headline_F": reach.get("headline_under_reached_fraction"),
        # The LADDER is nested under ``curve``, not flat on the reach block:
        # app/server/reach.py returns {"curve": {"x": [...],
        # "under_reached_fraction": [...]}}. Reading it off the top level
        # silently yields None for every row, which is how the first pass of
        # this sweep recorded 822 empty ladders and no error.
        "reach_ladder_x": (reach.get("curve") or {}).get("x"),
        "reach_ladder_F": (reach.get("curve") or {}).get(
            "under_reached_fraction"),
        "reach_window_mode": reach.get("window_mode"),
        "reach_E_ref": reach.get("E_ref_uM_s") or reach.get("E_ref"),
        "n_viable_reach": reach.get("n_viable"),
        "metric_reasons": result.get("metric_reasons") or {},
        "alpha_at_run": prov.get("alpha"),
        "log_tail": job.logs[-12:],
    }
