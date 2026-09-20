"""Is the spot-data refusal over-applied? — a jitter-reseed test on the four Visium GBM sections.

``app/server/api.py::SPOT_SPATIAL_REFUSAL`` is attached, by platform alone, to
every Visium section. Its stated reason is specific and it is a good reason: a
~55 um Visium spot is not a cell, the pipeline SYNTHESISES
``round(frac * 12)`` pseudo-cells per spot at the spot centroid with a seeded
uniform-disc jitter (``reconstruct.real_visium._place_cells``), and a per-cell
spatial statistic computed on those positions "would measure that jitter, not
the tissue".

THIS EXPERIMENT DOES NOT ARGUE WITH THAT REASON. It takes the reason literally
and asks, quantity by quantity, WHICH of the shipped run's outputs are in fact
functions of the jitter. A quantity that is provably not a function of the
jitter is not defended by that reason, and a refusal that lands on it is landing
outside its own stated subject. A quantity that IS a function of the jitter
stays refused, and nothing below weakens that.

Run::

    ./.venv/bin/python experiments/visium_section_level.py

Writes its raw record to ``--out`` (default: a scratchpad path). It writes
NOTHING under ``results/`` and touches nothing under ``tmesim/`` or ``app/``.


PRE-REGISTERED CRITERION — WRITTEN BEFORE THE FIRST RUN
=======================================================
A scalar emitted by the shipped pipeline on a Visium section is called
**RECOVERABLE** if and only if it passes ALL THREE of the following. Failing any
one leaves the refusal in force for that quantity. The thresholds and the
per-quantity class assignments below were fixed before any number was computed.

**C1 — JITTER-INVARIANCE (measured).** Re-run the shipped pipeline with the
pseudo-cell placement seed changed from the shipped ``0`` to ``1`` and ``2``,
changing nothing else. For every one of the four sections, the quantity must
satisfy

  * ``RATIO`` rule (quantities that are dimensionless ratios or have physical
    units): ``max_s |x_s - x_0| / |x_0| <= 0.01`` (1 %);
  * ``FRACTION`` rule (quantities that are fractions in [0, 1], where a relative
    test near zero is unfair): ``max_s |x_s - x_0| <= 0.01`` absolute.

  Three seeds is a small sample and is stated as such. C1 is a test of "is this
  a function of the jitter DRAW", which is exactly what the refusal's stated
  reason asserts. It is NOT a test of "is this a function of the synthesis", and
  C2 exists because it is not.

**C2 — MULTIPLICITY-INDEPENDENCE (judged from code, declared in advance).** The
quantity must be definable without the constant ``_CELLS_PER_SPOT = 12``. A cell
COUNT is ``12 x (spot occupancy)`` and is an artefact of that constant no matter
how jitter-invariant it is; it fails C2. A quantity computed over the grid or
over the tissue mask does not involve the constant at all and passes. A quantity
computed over the pseudo-cell population but as a FRACTION of it is recorded as
``WEIGHTED`` — it is a spot-level statistic weighted by the deconvolved
fractions, which is a real thing, but it is NOT a per-cell statistic and must
never be reported as one; ``WEIGHTED`` is treated as a C2 FAIL for the purposes
of the verdict, and the reasoning is in the report rather than hidden in a
threshold.

**C3 — RESOLUTION HONESTY (judged from the definition, declared in advance).**
The quantity's definition must contain no length scale below the Visium spot
pitch (~100 um). Nothing below that pitch was measured on these sections, so a
statistic whose definition reaches under it is reporting a resolution the data
does not have. Passing C3 does NOT make the resolution limit go away: every
recovered quantity is a ~100 um-resolution quantity and the report says so on
every line.

**FALSIFIER, stated in advance.** If ``mean_vascular_density``,
``mean_ecm_density`` or ``mean_tissue_auc`` moves by more than the C1 threshold
under a reseed, then the delivery source on a Visium section IS downstream of
the synthesised positions, the refusal's stated mechanism is right as written,
and the answer to this whole experiment is "no, the refusal is correctly
scoped". That is a complete result and it is the one this file is prepared to
report.

**TWO STRUCTURAL CHECKS, with a positive control.**

1. ``reconstruct.real_fields.vascular_density_from_cells`` — the function that
   builds a vessel source FROM ENDOTHELIAL CELL POSITIONS — is instrumented with
   a call counter. On the four Visium sections the expected count is **0**; on
   ``xenium_319`` (the positive control) it is expected to be **> 0**. If the
   Visium count is non-zero the refusal's mechanism sentence is correct and this
   experiment stops there.
2. The three arrays the transport solve actually consumes
   (``domain.tissue_mask``, ``fields.vascular_density``, ``fields.ecm_density``)
   are compared BITWISE across seeds. The AUC field is a deterministic function
   of those three plus the drug and the PDE config, so bitwise-identical inputs
   are a stronger statement than any tolerance on a summary of the output.

**WHAT IS NOT UNDER TEST AND IS NOT RECOVERED BY ANYTHING HERE.** Cross-compound
ranking (withdrawn 2026-08-21, absolutely). Absolute kill fraction. Any clinical
or diagnostic reading. This file runs ONE compound on several sections; it does
not compare compounds. ``kill_fraction`` is recorded because the shipped payload
carries it and deleting a number the product prints would be its own dishonesty;
it is flagged and no verdict rests on it.
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

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

ALPHA_EXPECTED = 22.223477711172563

#: ONE compound, every section. Not a comparison between compounds.
DRUG = "temozolomide"

#: The four Visium GBM sections R12 names, and the single-cell positive control.
VISIUM: List[str] = [
    "visium_10x_parent_gbm",
    "visium_gbm1",
    "visium_gbm2",
    "visium_gbm3",
]
CONTROL = "xenium_319"

#: Shipped placement seed is 0; 1 and 2 are the probes.
SEEDS = (0, 1, 2)

#: Same configuration the 2026-08-23 sweep used, so its rows are comparable.
DX_UM = 30.0
TOTAL_TIME_S = 3600.0 * 6

# --------------------------------------------------------------------------- #
# The pre-registered classification. FIXED BEFORE THE RUN.
#
#   key -> (C1 rule, C2 class, C3 verdict, one-line reason)
#
# C1 rule   : "ratio" (relative <= 1%) or "fraction" (absolute <= 0.01)
# C2 class  : "grid"      - computed over the grid/mask; _CELLS_PER_SPOT absent
#             "weighted"  - computed over the pseudo-cell population (FAIL)
#             "count"     - a cell count; 12 x occupancy by construction (FAIL)
#             "scalar"    - a compiled constant or a PK input; no tissue at all
# C3        : True if no length below the ~100 um spot pitch enters its
#             definition.
# --------------------------------------------------------------------------- #
CLASSIFY: Dict[str, tuple] = {
    # --- grid / mask only: no cell enters the definition --------------------
    "tissue_coverage_fraction": ("fraction", "grid", True,
        "fraction of the domain inside the spot-support mask"),
    "mean_vascular_density": ("ratio", "grid", True,
        "mean of the 100 um KDE of the per-spot ENDOTHELIAL fraction"),
    "mean_ecm_density": ("ratio", "grid", True,
        "mean of the 100 um KDE of the per-spot STROMAL fraction"),
    "mean_tissue_auc": ("ratio", "grid", True,
        "area-weighted mean of the solved AUC field over the tissue mask"),
    "mean_oxygen": ("ratio", "grid", True,
        "mean of the solved oxygen field; its SINK reads tumour-cell density"),
    "penetration_length_analytic_um": ("ratio", "grid", True,
        "sqrt(D_eff/k_eff) from mean_ecm_density; closed form, no field solve"),
    "plasma_peak_uM": ("ratio", "scalar", True, "a PK input, not a tissue read"),
    "plasma_peak_over_ic50": ("ratio", "scalar", True,
        "PK input over a compiled IC50; no tissue"),
    # --- read the solved field AT PSEUDO-CELL POSITIONS ---------------------
    "mean_tumor_auc": ("ratio", "weighted", True,
        "AUC sampled at pseudo-cell positions, averaged over them"),
    "median_tumor_auc": ("ratio", "weighted", True, "as above, median"),
    "E_p10_uM": ("ratio", "weighted", True,
        "10th percentile of per-pseudo-cell exposure"),
    "E_p50_uM": ("ratio", "weighted", True, "50th percentile, same population"),
    "E_p90_uM": ("ratio", "weighted", True, "90th percentile, same population"),
    "E_p90_over_p10": ("ratio", "weighted", True,
        "THE HEADLINE. A percentile ratio over the pseudo-cell population"),
    "therapeutic_coverage_1x": ("fraction", "weighted", True,
        "fraction of pseudo-tumour-cells above 1x IC50"),
    "median_E_over_ic50": ("ratio", "weighted", True, "same population"),
    "p10_E_over_ic50": ("ratio", "weighted", True, "same population"),
    "tissue_median_over_plasma_peak": ("ratio", "grid", True,
        "median of the AUC field over the tissue mask, over a PK input"),
    "viable_tumor_fraction": ("fraction", "weighted", True,
        "state labels banded at pseudo-cell positions"),
    "hypoxic_tumor_fraction": ("fraction", "weighted", True, "as above"),
    "low_o2_over_high_o2_exposure": ("ratio", "weighted", True,
        "an exposure ratio between two pseudo-cell bands"),
    # --- explicitly sub-spot-pitch: C3 FAIL by construction -----------------
    "penetration_length_um": ("ratio", "grid", False,
        "fitted against distance to the nearest endothelial PSEUDO-CELL; the "
        "fit's own length scale is tens of um, under the spot pitch"),
    "reach_headline_F": ("fraction", "weighted", False,
        "binned on the same sub-pitch vessel-distance transform"),
    # --- counts: 12 x occupancy by construction -----------------------------
    "n_cells": ("ratio", "count", True, "12 x spot occupancy"),
    "n_tumor": ("ratio", "count", True, "12 x spot occupancy x tumour fraction"),
    "n_viable_reach": ("ratio", "count", True, "a pseudo-cell count"),
    # --- never a defensible output, recorded only ---------------------------
    "kill_fraction": ("fraction", "weighted", True,
        "NOT DEFENSIBLE as an output at all; recorded, never read"),
}

#: The three arrays the transport solve consumes. Compared BITWISE.
TRANSPORT_INPUTS = ("tissue_mask", "vascular_density", "ecm_density")


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
# Instrumentation. Both hooks are installed in THIS process only and neither
# changes any file under tmesim/ or app/.
# --------------------------------------------------------------------------- #
_VDC_CALLS = {"n": 0}


def install_hooks(seed_box: Dict[str, int]) -> None:
    """Reseed ``_place_cells`` and count ``vascular_density_from_cells`` calls.

    ``reconstruct_from_visium`` calls the module-global ``_place_cells`` with no
    ``seed``, so it takes the default 0. Rebinding the global is the smallest
    possible probe: the cell COUNTS are ``round(frac * 12)`` and do not depend on
    the RNG at all, so the only thing a different seed changes is WHERE inside
    its spot disc each pseudo-cell lands — which is precisely the jitter the
    refusal names.
    """
    import tmesim.reconstruct.real_visium as rv
    import tmesim.reconstruct.real_pipeline as rp

    if getattr(rv, "_vsl_hooked", False):
        return

    orig_place = rv._place_cells

    def place(*a, **kw):
        kw["seed"] = int(seed_box["seed"])
        return orig_place(*a, **kw)

    rv._place_cells = place
    rv._vsl_hooked = True

    orig_vdc = rp.vascular_density_from_cells

    def vdc(*a, **kw):
        _VDC_CALLS["n"] += 1
        return orig_vdc(*a, **kw)

    rp.vascular_density_from_cells = vdc


# --------------------------------------------------------------------------- #
# Stage 1 — the transport inputs, bitwise, across seeds.
# --------------------------------------------------------------------------- #
def transport_input_check(section: str, seed_box: Dict[str, int]) -> Dict:
    from app.server.cache import lookup_section
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    rec = lookup_section(section) or {}
    source = rec.get("source")
    sample = rec.get("sample")

    arrays: Dict[int, Dict[str, np.ndarray]] = {}
    coords: Dict[int, np.ndarray] = {}
    vdc_by_seed: Dict[int, int] = {}
    meta: Dict = {}
    for s in SEEDS:
        seed_box["seed"] = s
        before = _VDC_CALLS["n"]
        micro = reconstruct_real(source=source, sample=sample, dx_um=DX_UM)
        vdc_by_seed[s] = _VDC_CALLS["n"] - before
        arrays[s] = {
            "tissue_mask": np.asarray(micro.domain.tissue_mask).copy(),
            "vascular_density": np.asarray(
                micro.fields.vascular_density, dtype=np.float64).copy(),
            "ecm_density": np.asarray(
                micro.fields.ecm_density, dtype=np.float64).copy(),
        }
        if micro.fields.oxygen is not None:
            arrays[s]["oxygen"] = np.asarray(
                micro.fields.oxygen, dtype=np.float64).copy()
        coords[s] = np.column_stack(
            [np.asarray(micro.cells.x), np.asarray(micro.cells.y)]).copy()
        if s == SEEDS[0]:
            rp_meta = (micro.metadata or {}).get("real_pipeline") or {}
            meta = {
                "vascular_source": rp_meta.get("vascular_source"),
                "ecm_source": rp_meta.get("ecm_source"),
                "source": rp_meta.get("source"),
                "n_spots": (micro.metadata or {}).get("n_spots"),
                "n_cells": int(micro.cells.n),
                "grid": list(micro.domain.shape),
                "dx_um": float(micro.domain.dx),
            }
        del micro

    out: Dict = {"section": section, "meta": meta,
                 "vascular_density_from_cells_calls": vdc_by_seed}

    # Did the probe actually move the cells? A reseed that changed nothing would
    # make every invariance below vacuous, so this is checked first and the
    # whole section is marked VACUOUS if it fails.
    moved = {}
    for s in SEEDS[1:]:
        if coords[s].shape != coords[SEEDS[0]].shape:
            moved[s] = {"same_shape": False}
            continue
        d = np.abs(coords[s] - coords[SEEDS[0]])
        moved[s] = {
            "same_shape": True,
            "n_cells": int(coords[s].shape[0]),
            "frac_moved": float((d.max(axis=1) > 0).mean()),
            "median_shift_um": float(np.median(np.hypot(
                coords[s][:, 0] - coords[SEEDS[0]][:, 0],
                coords[s][:, 1] - coords[SEEDS[0]][:, 1]))),
            "max_shift_um": float(np.max(np.hypot(
                coords[s][:, 0] - coords[SEEDS[0]][:, 0],
                coords[s][:, 1] - coords[SEEDS[0]][:, 1]))),
        }
    out["probe_moved_cells"] = moved

    fields: Dict[str, Dict] = {}
    for name in list(TRANSPORT_INPUTS) + ["oxygen"]:
        if name not in arrays[SEEDS[0]]:
            continue
        a0 = arrays[SEEDS[0]][name]
        row = {"bitwise_identical": True, "max_abs_diff": 0.0,
               "max_rel_diff": 0.0}
        for s in SEEDS[1:]:
            a = arrays[s][name]
            same = bool(a.shape == a0.shape and np.array_equal(a, a0))
            row["bitwise_identical"] &= same
            if a.shape == a0.shape:
                d = np.abs(a.astype(np.float64) - a0.astype(np.float64))
                row["max_abs_diff"] = max(row["max_abs_diff"], float(d.max()))
                den = np.abs(a0.astype(np.float64))
                nz = den > 0
                if nz.any():
                    row["max_rel_diff"] = max(
                        row["max_rel_diff"], float((d[nz] / den[nz]).max()))
        fields[name] = row
    out["transport_inputs"] = fields
    return out


# --------------------------------------------------------------------------- #
# Stage 2 — the shipped pipeline, per seed.
# --------------------------------------------------------------------------- #
class _StdoutJob:
    def __init__(self, run_id: str, spec: Dict):
        self.run_id = run_id
        self.spec = spec
        self.stage = "queued"
        self.progress = 0.0
        self.cancel_event = threading.Event()
        self.logs: List[str] = []

    def log(self, msg: str) -> None:
        self.logs.append(str(msg))


def _spec_for(section: str) -> Dict:
    from tmesim.core import PKConfig
    from tmesim.reconstruct.real_pipeline import (
        OXYGEN_MODE_DEFAULT, STATE_MODE_DEFAULT)

    pk = PKConfig()
    return {
        "section": section,
        "dx_um": DX_UM,
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
        "total_time_s": TOTAL_TIME_S,
        "attribute_mechanism": True,
        "timelapse_frames": 0,
    }


def _scalars(result: Dict) -> Dict:
    report = result.get("report") or {}
    metrics = result.get("metrics") or {}
    extra = result.get("extra_metrics") or {}
    reach = result.get("reach") or {}
    mech = result.get("mechanism") or {}
    e90, e10 = extra.get("E_p90_uM"), extra.get("E_p10_uM")
    spread = (float(e90) / float(e10)
              if isinstance(e90, (int, float)) and isinstance(e10, (int, float))
              and e10 else None)
    breakdown = (mech.get("survival_breakdown") or {}) if mech else {}
    dominant = max(breakdown.items(), key=lambda kv: kv[1]) if breakdown else None
    out = {
        "n_cells": report.get("n_cells"),
        "tissue_coverage_fraction": report.get("tissue_coverage_fraction"),
        "mean_vascular_density": report.get("mean_vascular_density"),
        "mean_ecm_density": report.get("mean_ecm_density"),
        "mean_oxygen": report.get("mean_oxygen"),
        "vascular_source": report.get("vascular_source"),
        "grid": report.get("grid") or report.get("grid_shape"),
        "n_tumor": extra.get("n_tumor"),
        "n_viable_reach": reach.get("n_viable"),
        "E_p10_uM": e10,
        "E_p50_uM": extra.get("E_p50_uM"),
        "E_p90_uM": e90,
        "E_p90_over_p10": spread,
        "mean_tumor_auc": metrics.get("mean_tumor_auc"),
        "median_tumor_auc": extra.get("median_tumor_auc"),
        "mean_tissue_auc": extra.get("mean_tissue_auc"),
        "therapeutic_coverage_1x": metrics.get("therapeutic_coverage_1x"),
        "median_E_over_ic50": extra.get("median_E_over_ic50"),
        "p10_E_over_ic50": extra.get("p10_E_over_ic50"),
        "tissue_median_over_plasma_peak": extra.get(
            "tissue_median_over_plasma_peak"),
        "plasma_peak_uM": extra.get("plasma_peak_uM"),
        "plasma_peak_over_ic50": extra.get("plasma_peak_over_ic50"),
        "viable_tumor_fraction": extra.get("viable_tumor_fraction"),
        "hypoxic_tumor_fraction": extra.get("hypoxic_tumor_fraction"),
        "low_o2_over_high_o2_exposure": extra.get("low_o2_over_high_o2_exposure"),
        "penetration_length_um": metrics.get("penetration_length_um"),
        "penetration_length_analytic_um": extra.get(
            "penetration_length_analytic_um"),
        "reach_headline_x": reach.get("headline_x"),
        "reach_headline_F": reach.get("headline_under_reached_fraction"),
        "reach_null_reason": reach.get("null_reason"),
        "kill_fraction": metrics.get("kill_fraction"),
        "mech_dominant_limiter": dominant[0] if dominant else None,
        "mech_dominant_share": dominant[1] if dominant else None,
        "metric_reasons": result.get("metric_reasons") or {},
    }
    return out


def run_one(section: str, seed: int, seed_box: Dict[str, int]) -> Dict:
    from app.server.cache import ReconstructionCache
    from app.server.jobs import _Ctx, execute_run

    seed_box["seed"] = seed
    spec = _spec_for(section)
    job = _StdoutJob(f"{section}|seed{seed}", spec)
    ctx = _Ctx(job)
    cache = ReconstructionCache(maxsize=1)  # never reuse another seed's build
    t0 = time.time()
    try:
        result = execute_run(spec, ctx, cache=cache)
    except Exception as exc:  # noqa: BLE001
        return {"status": "FAILED", "error_type": type(exc).__name__,
                "error": str(exc), "traceback": traceback.format_exc(),
                "wall_clock_s": time.time() - t0}
    result.pop("_payloads", None)
    out = _scalars(result)
    out["status"] = "ok"
    out["wall_clock_s"] = time.time() - t0
    return out


# --------------------------------------------------------------------------- #
# Stage 3 — apply the pre-registered criterion.
# --------------------------------------------------------------------------- #
def deviation(values: List[Optional[float]], rule: str) -> Optional[float]:
    """Max deviation of the probe seeds from the shipped seed, per ``rule``."""
    if not values or values[0] is None:
        return None
    try:
        base = float(values[0])
    except (TypeError, ValueError):
        return None
    worst = 0.0
    for v in values[1:]:
        if v is None:
            return None
        d = abs(float(v) - base)
        if rule == "ratio":
            if base == 0.0:
                return float("inf") if d > 0 else 0.0
            d = d / abs(base)
        worst = max(worst, d)
    return worst


def classify(per_section: Dict[str, Dict[int, Dict]]) -> Dict:
    rows: Dict[str, Dict] = {}
    for key, (rule, c2, c3, why) in CLASSIFY.items():
        thr = 0.01
        devs: Dict[str, Optional[float]] = {}
        vals: Dict[str, List] = {}
        for sec, by_seed in per_section.items():
            v = [by_seed.get(s, {}).get(key) for s in SEEDS]
            vals[sec] = v
            devs[sec] = deviation(v, rule)
        finite = [d for d in devs.values() if d is not None]
        c1_pass = bool(finite) and len(finite) == len(devs) and max(finite) <= thr
        c2_pass = c2 == "grid" or c2 == "scalar"
        verdict = "RECOVERABLE" if (c1_pass and c2_pass and c3) else "REFUSED"
        fails = []
        if not c1_pass:
            fails.append("C1")
        if not c2_pass:
            fails.append("C2")
        if not c3:
            fails.append("C3")
        rows[key] = {
            "c1_rule": rule, "c1_threshold": thr,
            "c1_max_deviation_over_sections": (max(finite) if finite else None),
            "c1_deviation_by_section": devs,
            "c1_pass": c1_pass,
            "c2_class": c2, "c2_pass": c2_pass,
            "c3_pass": bool(c3),
            "verdict": verdict, "fails": fails,
            "why": why,
            "values_by_section": vals,
        }
    return rows


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="jitter-reseed test on Visium GBM")
    default_out = os.environ.get("VSL_OUT") or os.path.join(
        _HERE, "..", "..", "visium_section_level.json")
    p.add_argument("--out", default=os.path.abspath(default_out))
    p.add_argument("--sections", default=",".join(VISIUM))
    args = p.parse_args(argv)

    alpha_start = _assert_alpha("start")
    sections = [s.strip() for s in args.sections.split(",") if s.strip()]
    seed_box = {"seed": 0}
    install_hooks(seed_box)

    blob: Dict = {
        "experiment": "visium_section_level",
        "alpha_start": alpha_start,
        "drug": DRUG,
        "dx_um": DX_UM,
        "total_time_s": TOTAL_TIME_S,
        "seeds": list(SEEDS),
        "criterion": {
            "C1": "reseed jitter; ratio rule rel<=0.01, fraction rule abs<=0.01,"
                  " on EVERY section",
            "C2": "definable without _CELLS_PER_SPOT=12; 'weighted'/'count' FAIL",
            "C3": "no length scale below the ~100 um spot pitch",
        },
        "transport_input_checks": {},
        "runs": {},
        "control": {},
    }

    print("\n=== STAGE 1: transport inputs, bitwise, across jitter seeds ===",
          flush=True)
    for sec in sections:
        t0 = time.time()
        chk = transport_input_check(sec, seed_box)
        blob["transport_input_checks"][sec] = chk
        vs = chk["meta"].get("vascular_source")
        print(f"  {sec}: vascular_source={vs!r} "
              f"vdc_calls={chk['vascular_density_from_cells_calls']} "
              f"({time.time() - t0:.1f}s)", flush=True)
        for name, row in chk["transport_inputs"].items():
            print(f"      {name:18s} bitwise_identical="
                  f"{row['bitwise_identical']} max_abs_diff="
                  f"{row['max_abs_diff']:.3e}", flush=True)
        for s, m in chk["probe_moved_cells"].items():
            print(f"      probe seed {s}: moved {m.get('frac_moved')} of cells, "
                  f"median {m.get('median_shift_um')} um", flush=True)

    print("\n=== STAGE 1b: positive control on the single-cell section ===",
          flush=True)
    t0 = time.time()
    ctrl = transport_input_check(CONTROL, seed_box)
    blob["control"][CONTROL] = ctrl
    print(f"  {CONTROL}: vascular_source="
          f"{ctrl['meta'].get('vascular_source')!r} "
          f"vdc_calls={ctrl['vascular_density_from_cells_calls']} "
          f"({time.time() - t0:.1f}s)", flush=True)

    print("\n=== STAGE 2: the shipped pipeline, per seed ===", flush=True)
    per_section: Dict[str, Dict[int, Dict]] = {}
    for sec in sections:
        per_section[sec] = {}
        for s in SEEDS:
            out = run_one(sec, s, seed_box)
            per_section[sec][s] = out
            blob["runs"].setdefault(sec, {})[str(s)] = out
            print(f"  {sec} seed={s} status={out.get('status')} "
                  f"E_p90/p10={out.get('E_p90_over_p10')} "
                  f"mean_tissue_auc={out.get('mean_tissue_auc')} "
                  f"({out.get('wall_clock_s', 0):.1f}s)", flush=True)
        _assert_alpha(f"after {sec}")

    print("\n=== STAGE 3: the pre-registered criterion ===", flush=True)
    blob["classification"] = classify(per_section)
    for key, row in blob["classification"].items():
        d = row["c1_max_deviation_over_sections"]
        print(f"  {key:34s} {row['verdict']:12s} "
              f"C1max={'—' if d is None else f'{d:.3e}'} "
              f"C2={row['c2_class']:8s} C3={row['c3_pass']} "
              f"{'fails ' + ','.join(row['fails']) if row['fails'] else ''}",
              flush=True)

    blob["alpha_end"] = _assert_alpha("end")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(blob, fh, indent=1, default=str)
    print(f"\n[done] -> {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
