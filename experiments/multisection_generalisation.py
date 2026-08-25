"""Does ANY of it generalise off xenium_319? — the shipped pipeline, every runnable section.

Every quantitative claim in this repository is "on this section": ONE 10x Xenium
GBM section (``xenium_319``, 40,887 cells). The reach ladder, the delivery
spread, the falloff, the fitted penetration length, the per-cell mechanism split
— all argued on that one specimen. Seven sections are runnable. This driver runs
the SHIPPED pipeline (``app.server.jobs.execute_run``, the exact function the
Studio queue calls) on each of them, with ONE compound and ONE configuration, so
that the SECTION is the only thing that varies.

WHAT IT DOES NOT DO
-------------------
* It does not tune anything. Nothing under ``tmesim/`` or ``app/`` is touched.
  ``calibrate_alpha()`` is asserted live at ``22.223477711172563`` before the
  first section and after the last.
* It does not work around a refusal. A section that cannot produce a quantity
  gets ``null_reason`` / ``metric_reasons`` recorded verbatim, never a
  substitute.
* It does not build a cross-section comparison on absolute kill fraction.
  ``kill_fraction`` is recorded because the shipped payload carries it and
  suppressing a number the product prints would be its own dishonesty, but it is
  flagged ``not_defensible`` in the payload and the report must not rank on it.
* It compares ONE compound across SECTIONS. It does not rank compounds — that
  claim was withdrawn 2026-08-21 and the prohibition is absolute.

THE SPOT REFUSAL IS LOAD-BEARING
--------------------------------
``app.server.api.SPOT_SPATIAL_REFUSAL`` declares that on a Visium (spot) section
the per-cell spatial statistics AND the delivery map, its falloff, its
perivascular reference and the delivery spread are all INVALID, because the
delivery source is a density over ENDOTHELIAL CELL POSITIONS and those positions
are synthesised (12 pseudo-cells per spot at the spot centroid with uniform
jitter), not measured. This driver runs the Visium sections anyway — the engine
does not raise — and records both the numbers the pipeline emits AND the refusal
that says they are not readings of that tissue. The report quarantines them.

Run::

    ./.venv/bin/python experiments/multisection_generalisation.py \
        --sections xenium_319,visium_gbm1 --dx 30

Writes (append-safe, one section at a time, so a killed run keeps its work):
    results/multisection_generalisation.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
import traceback
from typing import Dict, List, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_RESULTS = os.path.join(_ROOT, "results")
_DATA = os.path.join(_RESULTS, "multisection_generalisation.json")

#: The single value ``calibrate_alpha()`` must return, asserted at the top and
#: the bottom of every stage. If this moves, the run is describing a different
#: engine and every number in it is void.
ALPHA_EXPECTED = 22.223477711172563

#: ONE compound, everywhere. Named in the same sentence as every number.
DRUG = "temozolomide"

#: The seven built-in sections, in registry order.
BUILTIN_SECTIONS: List[str] = [
    "xenium_319",
    "xenium_io",
    "gse324268",
    "visium_10x_parent_gbm",
    "visium_gbm1",
    "visium_gbm2",
    "visium_gbm3",
]


def _assert_alpha(where: str) -> float:
    from tmesim.validation.calibration import calibrate_alpha

    a = float(calibrate_alpha())
    if repr(a) != repr(ALPHA_EXPECTED):
        raise SystemExit(
            f"[{where}] calibrate_alpha() = {a!r}, expected {ALPHA_EXPECTED!r}. "
            f"The engine is not the one this experiment claims to measure; "
            f"refusing to produce numbers."
        )
    print(f"[alpha] {where}: calibrate_alpha() = {a!r}  OK", flush=True)
    return a


# --------------------------------------------------------------------------- #
# The smallest object app.server.jobs.execute_run will accept as a job context.
# It is NOT a reimplementation of anything: execute_run only ever calls
# ctx.log / ctx.stage / ctx.check_cancelled, and _Ctx forwards those to a Job.
# --------------------------------------------------------------------------- #
class _StdoutJob:
    def __init__(self, run_id: str, spec: Dict, quiet: bool = False):
        self.run_id = run_id
        self.spec = spec
        self.stage = "queued"
        self.progress = 0.0
        self.cancel_event = threading.Event()
        self.logs: List[str] = []
        self._quiet = quiet

    def log(self, msg: str) -> None:
        line = f"    [{self.run_id}] {msg}"
        self.logs.append(str(msg))
        if not self._quiet:
            print(line, flush=True)


def _spec_for(section: str, dx_um: float, total_time_s: float) -> Dict:
    """The shipped Studio default RunSpec, with ONLY the section varying."""
    from tmesim.core import PKConfig
    from tmesim.reconstruct.real_pipeline import (
        OXYGEN_MODE_DEFAULT,
        STATE_MODE_DEFAULT,
    )

    pk = PKConfig()
    return {
        "section": section,
        "dx_um": float(dx_um),
        "drug": DRUG,
        "sensitivity_mode": "uniform",
        "vascular_mode": "discrete",
        "oxygen_mode": OXYGEN_MODE_DEFAULT,
        "state_mode": STATE_MODE_DEFAULT,
        "mgmt": None,
        "pk": {
            "plasma_peak": float(pk.plasma_peak),
            "half_life_s": float(pk.half_life),
            "n_doses": int(pk.n_doses),
            "dose_interval_s": float(pk.dose_interval),
        },
        "total_time_s": float(total_time_s),
        # The mechanism split IS one of the readouts under test, and it is a
        # read-out of arrays the kill model already computed (the run is
        # bit-identical with it on or off). No time-lapse: it changes nothing
        # and costs memory.
        "attribute_mechanism": True,
        "timelapse_frames": 0,
    }


def _section_record(section: str) -> Dict:
    """Registry record + the platform caveats the shipped API attaches."""
    from app.server.api import _with_caveats
    from app.server.cache import lookup_section

    rec = lookup_section(section) or {}
    out = dict(rec)
    out.setdefault("id", section)
    try:
        return _with_caveats(out)
    except Exception:  # noqa: BLE001
        return out


def _extract(result: Dict, section: str, rec: Dict, spec: Dict,
             wall_s: float) -> Dict:
    """Pull the product's actual shipped readouts out of a RunResult body."""
    report = result.get("report") or {}
    prov = result.get("provenance") or {}
    metrics = result.get("metrics") or {}
    extra = result.get("extra_metrics") or {}
    reach = result.get("reach") or {}
    mech = result.get("mechanism") or {}

    e_p90 = extra.get("E_p90_uM")
    e_p10 = extra.get("E_p10_uM")
    spread = None
    if isinstance(e_p90, (int, float)) and isinstance(e_p10, (int, float)) \
            and e_p10 not in (0, None):
        spread = float(e_p90) / float(e_p10)

    breakdown = (mech.get("survival_breakdown") or {}) if mech else {}
    dominant = None
    if breakdown:
        dominant = max(breakdown.items(), key=lambda kv: kv[1])

    curve = (reach.get("curve") or {}) if reach else {}
    falloff = (reach.get("falloff") or {}) if reach else {}

    caveat_rows = rec.get("caveats") or []

    return {
        "section": section,
        "label": rec.get("label"),
        "platform": rec.get("platform") or rec.get("source"),
        "spatial_spread_valid": rec.get("spatial_spread_valid"),
        "section_caveats": caveat_rows,
        "spec": spec,
        "wall_clock_s": wall_s,
        "status": "ok",
        # --- geometry / population ------------------------------------------ #
        "n_cells": report.get("n_cells"),
        "grid": report.get("grid") or report.get("grid_shape"),
        "extent_um": report.get("extent_um"),
        "dx_um": spec["dx_um"],
        "celltype_counts": report.get("celltype_counts"),
        "celltype_fractions": report.get("celltype_fractions"),
        "tissue_coverage_fraction": report.get("tissue_coverage_fraction"),
        "tissue_mask_mode": report.get("tissue_mask_mode"),
        "oxygen_block": report.get("oxygen"),
        "cell_state_block": report.get("cell_state"),
        "vascular_source": report.get("vascular_source"),
        "mean_vascular_density": report.get("mean_vascular_density"),
        "mean_ecm_density": report.get("mean_ecm_density"),
        "mean_oxygen": report.get("mean_oxygen"),
        "n_tumor": extra.get("n_tumor"),
        "n_viable_reach": reach.get("n_viable"),
        "n_tumor_reach": reach.get("n_tumor"),
        "viable_tumor_fraction": extra.get("viable_tumor_fraction"),
        "hypoxic_tumor_fraction": extra.get("hypoxic_tumor_fraction"),
        # --- delivery spread (the headline that needs no threshold) --------- #
        "E_p10_uM": e_p10,
        "E_p50_uM": extra.get("E_p50_uM"),
        "E_p90_uM": e_p90,
        "E_p90_over_p10": spread,
        # --- reach ladder --------------------------------------------------- #
        "reach_null_reason": reach.get("null_reason"),
        "reach_headline_x": reach.get("headline_x"),
        "reach_headline_F": reach.get("headline_under_reached_fraction"),
        "reach_ladder_x": curve.get("x"),
        "reach_ladder_F": curve.get("under_reached_fraction"),
        "reach_window_mode": reach.get("window_mode"),
        "E_ref": reach.get("E_ref_uM_s"),
        "E_ref_units": reach.get("E_ref_units"),
        "falloff_distance_um": falloff.get("distance_um"),
        "falloff_normalised_mean_auc": falloff.get("normalised_mean_auc"),
        "falloff_count": falloff.get("count"),
        "penetration_fit_error": reach.get("penetration_fit_error"),
        # --- penetration length --------------------------------------------- #
        "penetration_length_um": metrics.get("penetration_length_um"),
        "penetration_length_analytic_um": extra.get(
            "penetration_length_analytic_um"),
        # --- exposure / coverage -------------------------------------------- #
        "mean_tumor_auc": metrics.get("mean_tumor_auc"),
        "median_tumor_auc": extra.get("median_tumor_auc"),
        "mean_tissue_auc": extra.get("mean_tissue_auc"),
        "therapeutic_coverage_1x": metrics.get("therapeutic_coverage_1x"),
        "median_E_over_ic50": extra.get("median_E_over_ic50"),
        "p10_E_over_ic50": extra.get("p10_E_over_ic50"),
        "plasma_peak_uM": extra.get("plasma_peak_uM"),
        "plasma_peak_over_ic50": extra.get("plasma_peak_over_ic50"),
        "tissue_median_over_plasma_peak": extra.get(
            "tissue_median_over_plasma_peak"),
        "low_o2_over_high_o2_exposure": extra.get("low_o2_over_high_o2_exposure"),
        # --- per-cell mechanism split ---------------------------------------- #
        "mechanism_reason": result.get("mechanism_reason"),
        "mech_counts": mech.get("counts") if mech else None,
        "mech_survival_breakdown": breakdown or None,
        "mech_dominant_limiter": dominant[0] if dominant else None,
        "mech_dominant_share": dominant[1] if dominant else None,
        "mech_joint_limited": mech.get("joint_limited") if mech else None,
        "mech_reconciliation": mech.get("reconciliation") if mech else None,
        "mech_thresholds": mech.get("thresholds") if mech else None,
        # --- NOT a defensible cross-section output --------------------------- #
        "kill_fraction": metrics.get("kill_fraction"),
        "kill_fraction_not_defensible": (
            "Absolute kill fraction is never a defensible output and no "
            "cross-section comparison in this experiment is built on it. "
            "Recorded because the shipped payload carries it."),
        # --- every refusal the pipeline emitted ------------------------------ #
        "metric_reasons": result.get("metric_reasons") or {},
        "run_caveats": prov.get("caveats") or [],
        "drug_params": result.get("drug_params"),
    }


def _load() -> Dict:
    if os.path.exists(_DATA):
        try:
            with open(_DATA) as fh:
                return json.load(fh)
        except Exception:  # noqa: BLE001
            pass
    return {"alpha": ALPHA_EXPECTED, "drug": DRUG, "runs": {}}


def _save(blob: Dict) -> None:
    os.makedirs(_RESULTS, exist_ok=True)
    tmp = _DATA + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(blob, fh, indent=1, default=str)
    os.replace(tmp, _DATA)


def run_section(section: str, dx_um: float, total_time_s: float,
                quiet: bool = False) -> Dict:
    from app.server.cache import ReconstructionCache
    from app.server.jobs import _Ctx, execute_run

    rec = _section_record(section)
    spec = _spec_for(section, dx_um, total_time_s)
    key = f"{section}@dx{dx_um:g}"
    print(f"\n=== {key} — {rec.get('label') or section} "
          f"[{rec.get('platform') or rec.get('source') or 'platform unknown'}]",
          flush=True)
    for c in (rec.get("caveats") or []):
        print(f"    [section caveat] {c.get('capability')}: {c.get('status')}",
              flush=True)

    job = _StdoutJob(key, spec, quiet=quiet)
    ctx = _Ctx(job)
    # A FRESH cache per section: the process-global RECON_CACHE would hold a
    # 20x-cells reconstruction resident for the rest of the sweep.
    cache = ReconstructionCache(maxsize=1)
    t0 = time.time()
    try:
        result = execute_run(spec, ctx, cache=cache)
    except Exception as exc:  # noqa: BLE001
        wall = time.time() - t0
        print(f"    !! {type(exc).__name__}: {exc}", flush=True)
        return {
            "section": section,
            "label": rec.get("label"),
            "platform": rec.get("platform") or rec.get("source"),
            "section_caveats": rec.get("caveats") or [],
            "spec": spec,
            "wall_clock_s": wall,
            "status": "FAILED",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
            "log_tail": job.logs[-40:],
        }
    wall = time.time() - t0
    result.pop("_payloads", None)
    out = _extract(result, section, rec, spec, wall)
    out["log_tail"] = job.logs[-40:]
    print(f"    -> done in {wall:.1f}s  "
          f"n_tumor={out['n_tumor']} "
          f"E_p90/p10={out['E_p90_over_p10']} "
          f"F({out['reach_headline_x']})={out['reach_headline_F']} "
          f"L={out['penetration_length_um']} "
          f"limiter={out['mech_dominant_limiter']}", flush=True)
    return out


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sections", default=",".join(BUILTIN_SECTIONS),
                   help="comma-separated section ids (default: the 7 built-ins)")
    p.add_argument("--dx", type=float, default=30.0,
                   help="grid spacing in um, IDENTICAL for every section")
    p.add_argument("--total-time", type=float, default=3600.0 * 6,
                   help="exposure window in s (default: the shipped 6 h)")
    p.add_argument("--tag", default="", help="suffix for the run key")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)

    _assert_alpha("start")
    blob = _load()
    sections = [s.strip() for s in args.sections.split(",") if s.strip()]
    print(f"[plan] {len(sections)} section(s) at dx = {args.dx:g} um, "
          f"drug = {DRUG}, window = {args.total_time:g} s", flush=True)

    t_all = time.time()
    for sid in sections:
        key = f"{sid}@dx{args.dx:g}" + (f"|{args.tag}" if args.tag else "")
        out = run_section(sid, args.dx, args.total_time, quiet=args.quiet)
        out["alpha_at_run"] = _assert_alpha(f"after {sid}")
        blob["runs"][key] = out
        blob["updated"] = time.strftime("%Y-%m-%d %H:%M:%S")
        _save(blob)
    print(f"\n[done] {len(sections)} section(s) in {time.time() - t_all:.1f}s "
          f"-> {_DATA}", flush=True)
    _assert_alpha("end")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
