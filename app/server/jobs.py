"""Run execution + the job queue.

The API must never block on a simulation, so every ``POST /api/runs`` enqueues a
:class:`Job` and returns immediately. A single background **worker thread**
drains the queue.

Why exactly one worker, and why a thread rather than a process:

* **Thread.** The expensive object is the reconstruction, which lives in
  :data:`app.server.cache.RECON_CACHE` as a graph of large numpy arrays. A
  process pool would have to rebuild or pickle it per run, throwing away the
  entire point of the cache. The heavy inner loops (numpy, ``scipy.sparse``
  factorisation and solves) release the GIL, so the API's event loop stays
  responsive while a run is in flight.
* **One.** Runs on the same key share the cached ``Microenvironment``, and
  ``apply_drug_sensitivity`` writes ``cells.sensitivity`` on it. Serialising the
  worker makes that safe by construction. It also stops N simultaneous runs from
  fighting over the same cores and making every one of them slower.

Everything scientific here is a call into the production entry points —
``reconstruct_real`` (via the cache), ``apply_drug_sensitivity``,
``calibrate_alpha`` -> ``calibrated_drugparams`` -> ``apply_barrier_to_drugparams``,
``simulate_response``, and ``validation.delivery_metrics.compute_all``. That is
the same chain ``experiments/patient_report.py`` runs, in the same order, with
the same production constants (see :data:`PDE_DT_S`), which is what makes an API
run numerically identical to a CLI run.
"""
from __future__ import annotations

import json
import math
import os
import queue
import threading
import time
import traceback
from collections import deque
from typing import Callable, Deque, Dict, List, Optional

import numpy as np

from .cache import (
    RECON_CACHE,
    ReconstructionCache,
    make_key,
    runs_dir,
    state_route_support,
)
from .payloads import (
    FIELD_NAMES,
    MAX_TIMELAPSE_FRAMES,
    encode_field,
    encode_geometry,
    encode_mechanism,
    encode_timelapse,
)

__all__ = [
    "PDE_DT_S",
    "PDE_SAVE_EVERY",
    "SIM_SEED",
    "DEFAULT_TIMELAPSE_FRAMES",
    "TIMELAPSE_BUDGET_BYTES",
    "RunCancelled",
    "Job",
    "JobManager",
    "execute_run",
    "jsonable",
]

# --------------------------------------------------------------------------- #
# Production simulation constants.
#
# The frozen RunSpec carries no solver settings, so these are fixed at the
# values every tmesim experiment uses (experiments/patient_report.py,
# validation/delivery_metrics.py: PDEConfig(dt=10.0, convection=True,
# save_every=500) and seed=0). They are reported in every RunResult's
# provenance so the choice is visible, never implicit.
# --------------------------------------------------------------------------- #
PDE_DT_S = 10.0
PDE_SAVE_EVERY = 500
PDE_CONVECTION = True
SIM_SEED = 0

#: Time-lapse frames captured per run unless the RunSpec says otherwise. 24 is
#: the largest count that keeps the flagship 677x425 grid under 30 MB on the
#: wire (1.15 MB/frame), and is dense enough that a 6 h window is sampled every
#: 15 min. ``MAX_TIMELAPSE_FRAMES`` (40) is the hard ceiling; 0 disables capture.
DEFAULT_TIMELAPSE_FRAMES = 24

#: Wire/disk budget for one run's time-lapse. If the requested frame count would
#: exceed it, the count is REDUCED to the largest that fits and the reason is
#: returned with the run — the payload is never silently truncated and the
#: capture is never silently skipped.
TIMELAPSE_BUDGET_BYTES = 48 * 1024 * 1024

#: Terminal statuses — a job in one of these is never re-run.
TERMINAL = ("done", "failed", "cancelled")

#: Longest log kept in memory per job; ``/status`` returns a tail of it.
_LOG_MAXLEN = 500


class RunCancelled(Exception):
    """Raised inside the worker when a job's cancel flag is seen."""


# --------------------------------------------------------------------------- #
# JSON sanitising
# --------------------------------------------------------------------------- #
#: Arrays longer than this are summarised rather than serialised — a RunResult
#: is a JSON document, and the bulk data has a binary endpoint.
_MAX_INLINE_ARRAY = 64


def jsonable(obj, _depth: int = 0):
    """Convert engine output into something ``json.dumps`` accepts.

    Non-finite floats become ``None`` rather than the ``NaN`` / ``Infinity``
    tokens Python emits, which are not valid JSON and which no browser will
    parse. A ``None`` produced this way always travels with a reason elsewhere
    in the payload — see :func:`_metric` — so a null is never mistaken for a
    measured zero.
    """
    if obj is None or isinstance(obj, (str, bool)):
        return obj
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return v if math.isfinite(v) else None
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        if obj.size > _MAX_INLINE_ARRAY:
            return {
                "_array_omitted": True,
                "shape": list(obj.shape),
                "dtype": str(obj.dtype),
                "reason": (
                    f"array of {obj.size} elements omitted from JSON; bulk "
                    "grids are served by /api/runs/{id}/field/{name}"
                ),
            }
        return [jsonable(v, _depth + 1) for v in obj.tolist()]
    if isinstance(obj, dict):
        return {str(k): jsonable(v, _depth + 1) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [jsonable(v, _depth + 1) for v in obj]
    return str(obj)


def _metric(value, reasons: Dict[str, str], name: str, reason: str):
    """Record a metric, or ``None`` plus ``reason`` when it is not a number.

    ``penetration_length_fit`` returns ``nan`` when the exponential-plus-floor
    fit does not converge, and the contract forbids inventing a substitute.
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        reasons[name] = reason
        return None
    if not math.isfinite(v):
        reasons[name] = reason
        return None
    return v


# --------------------------------------------------------------------------- #
# The run itself
# --------------------------------------------------------------------------- #
#: (stage, progress-at-stage-START). Progress is stage-granular: the engine
#: entry points are single blocking calls with no callback, and this backend
#: does not modify them, so a within-stage percentage would be fabricated.
STAGES = (
    ("queued", 0.00),
    ("reconstruct", 0.05),
    ("sensitivity", 0.45),
    ("compile_drug", 0.50),
    ("simulate", 0.55),
    ("metrics", 0.85),
    ("payloads", 0.93),
    ("done", 1.00),
)
_STAGE_PROGRESS = dict(STAGES)


class _Ctx:
    """Progress / logging / cancellation channel handed to :func:`execute_run`."""

    def __init__(self, job: "Job"):
        self.job = job

    def log(self, msg: str) -> None:
        self.job.log(msg)

    def stage(self, name: str) -> None:
        self.check_cancelled()
        self.job.stage = name
        self.job.progress = _STAGE_PROGRESS.get(name, self.job.progress)
        self.job.log(f"[stage] {name}")

    def check_cancelled(self) -> None:
        if self.job.cancel_event.is_set():
            raise RunCancelled(
                f"cancelled by request during stage {self.job.stage!r}"
            )


def execute_run(spec: Dict, ctx: _Ctx,
                cache: ReconstructionCache = RECON_CACHE) -> Dict:
    """Run one ``RunSpec`` end to end and return the ``RunResult`` body.

    Returns the dict that becomes ``RunResult`` minus the envelope fields
    (``run_id``, ``status``, ``runtime_s``) that :class:`JobManager` owns, plus
    a private ``_payloads`` entry holding the binary blobs to persist.

    A spec carrying ``kind='protocol'`` is a TPL program rather than a RunSpec
    and is handed to :func:`app.server.protocol_api.execute_protocol_job`, which
    runs a whole sweep as ONE job on this same worker thread and this same
    reconstruction cache. The dispatch is additive: a spec without ``kind``
    takes the original path, unchanged.
    """
    if spec.get("kind") == "protocol":
        from .protocol_api import execute_protocol_job

        return execute_protocol_job(spec, ctx, cache=cache)

    from tmesim.core import PDConfig, PDEConfig, PKConfig, SimConfig
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.library import get_molecule
    from tmesim.molecular.potency import resolve_ic50
    from tmesim.reconstruct.real_pipeline import (
        apply_drug_sensitivity,
        reconstruction_report,
    )
    from tmesim.sim.engine import simulate_response
    from tmesim.validation import delivery_metrics
    from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
    from tmesim.validation.pk_reality import penetration_length_um

    section = spec["section"]
    dx_um = float(spec["dx_um"])
    drug_name = spec.get("drug")
    sens_mode = spec.get("sensitivity_mode") or "uniform"
    vasc_mode = spec.get("vascular_mode") or "discrete"
    oxy_mode = spec.get("oxygen_mode") or "pde"
    state_mode = spec.get("state_mode") or "oxygen"
    mgmt = spec.get("mgmt")
    pk = spec.get("pk") or {}
    total_time_s = float(spec["total_time_s"])
    # The COUNTERFACTUAL arm of the rank attribution: an otherwise bit-identical
    # run with PDConfig(hypoxia_resistance=1, quiescent_resistance=1), i.e. this
    # tumour's oxygen-driven resistance switched off. Dividing the two kills
    # gives rho, the O2-RESISTANCE term in tmesim.report.attribution. It is a
    # RunSpec field rather than a private hook so that both arms of a report go
    # through this same queue and this same reconstruction cache — the arm is
    # visible in /api/runs, cancellable, and reproducible from its own spec.
    resistance_off = bool(spec.get("resistance_off", False))
    # --- mechanism attribution + time-lapse (both additive, both optional) --- #
    # attribute_mechanism is on by default: it is a read-out of arrays the kill
    # model already computed, costs a handful of numpy passes, and the run is
    # bit-identical with it on or off (tmesim.sim.mechanism). time_resolved is
    # OFF by default because it changes the physics — it is the only way tissue
    # concentration can fall, and therefore the only way a time-lapse can show
    # the drug CLEAR rather than plateau, but it is not the historical path.
    attribute_mechanism = bool(spec.get("attribute_mechanism", True))
    timelapse_frames = spec.get("timelapse_frames", DEFAULT_TIMELAPSE_FRAMES)
    timelapse_frames = int(timelapse_frames) if timelapse_frames is not None else 0
    time_resolved = bool(spec.get("time_resolved", False))

    caveats: List[str] = []
    metric_reasons: Dict[str, str] = {}

    # --- 1. reconstruction (cached) --------------------------------------- #
    ctx.stage("reconstruct")
    key = make_key(section, dx_um, vasc_mode, oxy_mode, state_mode)
    was_cached = cache.peek(key) is not None
    entry = cache.get_or_build(key, log=ctx.log)
    micro = entry["micro"]
    caveats.extend(entry.get("notes") or [])

    # Reset the shared cached microenvironment to the state reconstruct_real
    # left it in, so a previous run's per-cell sensitivity — or cell-state
    # labelling — can never leak into this one. (The worker is single-threaded,
    # so this is a full reset.) ``state_mode`` is part of the cache key, so a
    # differently-stated reconstruction is a different entry and is never reused
    # here; this reset guards the remaining case of two runs sharing one entry.
    micro.cells.sensitivity = np.array(entry["pristine_sensitivity"], copy=True)
    pristine_state = entry.get("pristine_state")
    if pristine_state is not None:
        micro.cells.state = np.array(pristine_state, copy=True)
    micro.metadata = dict(entry["pristine_metadata"])

    ctx.check_cancelled()

    # --- 2. per-cell sensitivity ------------------------------------------ #
    ctx.stage("sensitivity")
    sens_info: Optional[Dict] = None
    if drug_name:
        ctx.log(f"apply_drug_sensitivity({drug_name!r}, mode={sens_mode!r}"
                + (f", mgmt={mgmt!r}" if mgmt is not None else "") + ")")
        sens_info = apply_drug_sensitivity(
            micro, drug_name, mode=sens_mode,
            **({"mgmt": mgmt} if mgmt is not None else {}),
        )
        ctx.log(f"sensitivity mean={sens_info['mean']:.4f} "
                f"upper_bound={sens_info.get('upper_bound')}")
    else:
        ctx.log("no drug in spec: reconstruction-only run, sensitivity left at "
                "the drug-agnostic reconstruction default (uniform 1.0)")

    report = reconstruction_report(micro)
    report["vascular_mode"] = vasc_mode
    report["vascular_source"] = entry["report"].get("vascular_source")
    report["reconstruction_build_time_s"] = entry["build_time_s"]
    report["reconstruction_cache_hit"] = bool(was_cached)
    report["max_darcy_speed_um_s"] = entry.get("max_darcy_speed_um_s")

    oxy_block = report.get("oxygen") or {}
    if oxy_block.get("caveat"):
        caveats.append(f"oxygen ({oxy_block.get('mode')}): {oxy_block['caveat']}")

    # --- cell state: which route labelled cells.state, and what it cost ------ #
    state_block = report.get("cell_state") or {}
    if state_block.get("caveat"):
        caveats.append(
            f"cell_state ({state_block.get('mode')}): {state_block['caveat']}")
    unmeasured = state_block.get("unmeasured_states") or []
    if unmeasured:
        caveats.append(
            f"state_mode={state_mode!r}: {sorted(unmeasured)} is UNMEASURED on "
            "this panel and is assigned to NO cell. That is NOT the same "
            "statement as 'computed and found to be zero', and it must not be "
            "reported as if it were (results/cell_state_impact.md §3.1)."
        )
    # results/cell_state_impact.md §5: on a panel that carries NO hypoxia-response
    # transcript, no state route clears the bar — the default's HYPOXIC band is
    # empty for a bookkeeping reason, the marker route declares hypoxia UNMEASURED,
    # and the vessel-distance route is a geometric proxy. The honest output of such
    # a section is a RANKING, a coverage map and a delivery metric; the absolute
    # level is not earned under ANY route. Derived from the panel's gene list, not
    # hardcoded to a section id.
    support = state_route_support(section)
    if support["genes_known"] and not support["hypoxia_present"]:
        caveats.append(
            "ABSOLUTE KILL MUST NOT BE QUOTED FROM THIS SECTION UNDER ANY STATE "
            f"ROUTE. Its panel carries none of {support['hypoxia_missing']}, so "
            "hypoxia is UNMEASURED here and NO state route clears the bar "
            "(results/cell_state_impact.md §5): kill moves by up to 0.3252 "
            "absolute (1.77x relative) purely with the choice of state route, "
            "with tissue, drug, transport and RNG held identical. The CROSS-DRUG "
            "RANKING is invariant across all three routes and is the defensible "
            "output; kill_fraction and any kill ceiling on this section carry an "
            "unresolved uncertainty of that size."
        )
    if str(state_mode) != "oxygen":
        caveats.append(
            "sim.engine.simulate_response RE-BANDS every living cell from its "
            "OWN internal oxygen solve (cells.agents.update_cell_states), so the "
            f"state_mode={state_mode!r} labelling reported here is the "
            "RECONSTRUCTION's and did NOT drive this run's kill_fraction. It "
            "changes the reported band occupancies and their provenance, not the "
            "simulated kill. Propagating a state route into a kill means "
            "re-running cells.pkpd.apply_drug_kill on the engine's own transport "
            "result, which is what validation.cell_state_impact does."
        )
    if str(oxy_mode) != "pde":
        caveats.append(
            "sim.engine.simulate_response RE-SOLVES the steady-state oxygen PDE "
            "internally and overwrites fields.oxygen, so the "
            f"oxygen_mode={oxy_mode!r} field served by /field/oxygen is the "
            "RECONSTRUCTION's oxygen and did NOT drive this run's kill. See the "
            "note in reconstruct.real_pipeline.reconstruct_real."
        )

    # --- 3. drug compilation ---------------------------------------------- #
    ctx.stage("compile_drug")
    alpha = float(calibrate_alpha())
    ctx.log(f"calibrate_alpha() -> {alpha:.4f}")
    provenance: Dict = {
        "alpha": alpha,
        "alpha_source": (
            "validation.calibration.calibrate_alpha(): the single perivascular-"
            "sink constant, fit ONCE on the doxorubicin 40-50 um penetration "
            "half-distance (Primeau 2005, PMID 16361566) converted to an "
            "e-folding length L = d_half/ln2 ~= 64.9 um, and revalidated "
            "out-of-sample on mitoxantrone."
        ),
        "ic50_uM": None,
        "ic50_source": None,
        "sensitivity_note": None,
        "caveats": caveats,
        # --- additive, and mandatory per the backend brief ---------------- #
        "sensitivity_mode": sens_mode if drug_name else None,
        "sensitivity_upper_bound": None,
        "sensitivity_determinant_gene": None,
        # --- cell-state route, surfaced verbatim from the engine ------------ #
        # reconstruction_report()['cell_state'] is the pipeline's own provenance
        # block: mode, source, is_measured, caveat, the marker-coverage detail
        # and the resulting tumour band fractions. It is copied here whole rather
        # than summarised, because the coverage block is exactly what has to
        # travel with any number derived from a marker route
        # (results/cell_state_impact.md §1).
        "state_mode": state_mode,
        "cell_state": jsonable(state_block),
        "state_is_measured": bool(state_block.get("is_measured", False)),
        "state_source": state_block.get("source"),
        "state_unmeasured": sorted(unmeasured) if unmeasured else [],
        "state_fractions_tumor": jsonable(state_block.get("fractions_tumor")),
        "state_drives_kill": False,
        "state_drives_kill_reason": (
            "sim.engine.simulate_response re-bands cells.state from its own "
            "oxygen solve before applying the kill, so state_mode changes the "
            "REPORTED band occupancies and their provenance, not this run's "
            "kill_fraction."
        ),
        "pde_dt_s": PDE_DT_S,
        "pde_save_every": PDE_SAVE_EVERY,
        "pde_convection": PDE_CONVECTION,
        "resistance_off": resistance_off,
        "pd_hypoxia_resistance": 1.0 if resistance_off else float(
            PDConfig().hypoxia_resistance),
        "pd_quiescent_resistance": 1.0 if resistance_off else float(
            PDConfig().quiescent_resistance),
        # PDEConfig.time_resolved: whether the transport source follows the real
        # plasma(t) curve or its window average. It defaults to the engine's own
        # default (False, the historical Stage-0 path) and is only raised when
        # the RunSpec asks for it — which is what a time-lapse needs if the drug
        # is to be seen CLEARING rather than plateauing.
        "pde_time_resolved": time_resolved,
        "pde_time_resolved_default": getattr(PDEConfig(), "time_resolved", None),
        "seed": SIM_SEED,
        "engine_entrypoints": [
            "reconstruct.real_pipeline.reconstruct_real",
            "reconstruct.real_pipeline.apply_drug_sensitivity",
            "validation.calibration.calibrate_alpha",
            "validation.calibration.calibrated_drugparams",
            "molecular.binding_site_barrier.apply_barrier_to_drugparams",
            "sim.engine.simulate_response",
            "validation.delivery_metrics.compute_all",
        ],
    }

    if not drug_name:
        provenance["ic50_reason"] = (
            "no drug in this RunSpec: this is a reconstruction-only run, so "
            "there is no potency to resolve."
        )
        provenance["sensitivity_note"] = (
            "No drug was simulated. The reconstruction's per-cell sensitivity "
            "is the drug-agnostic default (uniform 1.0); no intrinsic "
            "resistance is asserted for any drug."
        )
        caveats.append(
            "This tissue is an untreated resection with no outcome label; "
            "nothing in this result has been scored against a measured "
            "treatment response."
        )
        ctx.stage("payloads")
        payloads = _build_payloads(micro, report, auc_field=None, log=ctx.log)
        for name in ("kill_fraction", "mean_tumor_auc",
                     "therapeutic_coverage_1x", "penetration_length_um"):
            metric_reasons[name] = (
                "requires a drug; this RunSpec has drug=null "
                "(reconstruction-only run)"
            )
        no_drug = ("requires a drug; this RunSpec has drug=null "
                   "(reconstruction-only run), so there is no kill to attribute "
                   "and no drug field to film")
        return {
            "report": jsonable(report),
            "metrics": {k: None for k in metric_reasons},
            "metric_reasons": metric_reasons,
            "provenance": jsonable(provenance),
            "drug_params": None,
            "mechanism": None,
            "mechanism_reason": no_drug,
            "timelapse": {"available": False, "reason": no_drug},
            "_payloads": payloads,
        }

    pot = resolve_ic50(drug_name)
    try:
        mol = get_molecule(drug_name)
    except KeyError:
        # "palbociclib" is curated in molecular.breast_library, not here (the
        # two drug libraries are kept disjoint), but is a valid GBM-context
        # drug too -- calibration._compile already special-cases it. See
        # app/server/api.py's drug_rows() for the fuller comment.
        if drug_name == "palbociclib":
            from tmesim.molecular.breast_library import get_molecule_breast

            mol = get_molecule_breast(drug_name)
        else:
            raise
    dp = apply_barrier_to_drugparams(
        calibrated_drugparams(drug_name, alpha), mol
    )
    ctx.log(f"compiled {drug_name}: ic50={dp.ic50:.4g} uM  D={dp.diffusivity:.1f} "
            f"um^2/s  bbb={dp.bbb_penetrance:.4f}  binding={dp.binding:.4g}")

    n_lines = pot.get("n_cell_lines")
    if pot["confidence"] == "measured":
        ic50_source = f"GDSC measured (n={n_lines})"
    elif pot["confidence"] == "literature":
        ic50_source = "curated literature"
    else:
        ic50_source = "flagged default (LOW confidence)"
    provenance.update({
        "ic50_uM": float(dp.ic50),
        "ic50_source": ic50_source,
        "ic50_confidence": pot["confidence"],
        "ic50_note": pot["note"],
        "ic50_n_cell_lines": n_lines,
        "ic50_resolved_uM": float(pot["ic50_uM"]),
    })
    if pot["confidence"] != "measured":
        caveats.append(
            f"IC50 for {drug_name} is NOT a GDSC measurement: {ic50_source}. "
            f"{pot['note']}"
        )

    if sens_info is not None:
        provenance["sensitivity_note"] = sens_info.get("comparability")
        provenance["sensitivity_upper_bound"] = bool(
            sens_info.get("upper_bound", True))
        provenance["sensitivity_determinant_gene"] = sens_info.get(
            "determinant_gene")
        provenance["sensitivity_kind"] = sens_info.get("kind")
        provenance["sensitivity_rationale"] = sens_info.get("rationale")
        provenance["sensitivity_mean"] = sens_info.get("mean")
        for k in ("mgmt_applied", "mgmt_scalar", "mgmt_input", "mgmt_model",
                  "mgmt_scalar_ci95", "mgmt_calibration_n", "mgmt_note",
                  "mgmt_resolution"):
            if k in sens_info:
                provenance[k] = sens_info[k]
        if provenance["sensitivity_upper_bound"]:
            caveats.append(
                "ABSOLUTE kill_fraction is a no-intrinsic-resistance UPPER "
                f"BOUND for {drug_name} in sensitivity_mode={sens_mode!r}. "
                + str(sens_info.get("comparability") or "")
            )
        if sens_info.get("mgmt_supplied") and not sens_info.get("mgmt_applied"):
            caveats.append(str(sens_info.get("mgmt_note")))

    # --- 4. simulate ------------------------------------------------------- #
    ctx.stage("simulate")
    pde = PDEConfig(dt=PDE_DT_S, total_time=total_time_s,
                    convection=PDE_CONVECTION, save_every=PDE_SAVE_EVERY,
                    time_resolved=time_resolved)
    n_frames, timelapse_reason = _plan_timelapse(
        timelapse_frames, micro.domain.shape, total_time_s)
    if timelapse_reason:
        ctx.log("timelapse: " + timelapse_reason)
    pk_cfg = PKConfig(
        plasma_peak=float(pk.get("plasma_peak", PKConfig().plasma_peak)),
        half_life=float(pk.get("half_life_s", PKConfig().half_life)),
        n_doses=int(pk.get("n_doses", PKConfig().n_doses)),
        dose_interval=float(pk.get("dose_interval_s", PKConfig().dose_interval)),
    )
    cfg_kwargs = {}
    if resistance_off:
        cfg_kwargs["pd"] = PDConfig(hypoxia_resistance=1.0,
                                    quiescent_resistance=1.0)
        caveats.append(
            "COUNTERFACTUAL ARM: this run was simulated with "
            "PDConfig(hypoxia_resistance=1, quiescent_resistance=1), i.e. this "
            "tumour's oxygen-driven resistance switched OFF. Its kill_fraction "
            "is NOT a prediction for this tissue — it is the denominator of "
            "rho = kill / kill(resistance off), the O2-RESISTANCE term of the "
            "rank attribution (tmesim.report.attribution). Do not read it, "
            "quote it, or compare it as a therapy result."
        )
    cfg = SimConfig(drug=dp, pk=pk_cfg, pde=pde, seed=SIM_SEED, **cfg_kwargs)
    ctx.log(f"simulate_response(dt={PDE_DT_S}s, total_time={total_time_s}s, "
            f"convection={PDE_CONVECTION}, seed={SIM_SEED}, "
            f"plasma_peak={pk_cfg.plasma_peak} uM, "
            f"half_life={pk_cfg.half_life} s, n_doses={pk_cfg.n_doses}, "
            f"resistance_off={resistance_off})")
    t_sim = time.time()
    res = simulate_response(
        micro, cfg,
        attribute=attribute_mechanism,
        capture_frames=(n_frames if n_frames >= 2 else None),
    )
    ctx.log(f"simulate_response done in {time.time() - t_sim:.1f}s: "
            f"kill_fraction={res.kill_fraction:.6f}")

    ctx.check_cancelled()

    # --- 5. metrics -------------------------------------------------------- #
    ctx.stage("metrics")
    dm = delivery_metrics.compute_all(micro, res, dp)
    head = dm["headline"]
    ctx.log(f"delivery_metrics: coverage_1x={head['therapeutic_coverage_1x_ic50']} "
            f"mean_tumor_auc={head['mean_tumor_auc_uM_s']} "
            f"L={head['penetration_length_L_um']}")

    metrics = {
        "kill_fraction": _metric(
            res.kill_fraction, metric_reasons, "kill_fraction",
            "simulate_response returned a non-finite kill fraction"),
        "mean_tumor_auc": _metric(
            head["mean_tumor_auc_uM_s"], metric_reasons, "mean_tumor_auc",
            "no TUMOR cells in this reconstruction, so there is no "
            "tumor-cell exposure to average"),
        "therapeutic_coverage_1x": _metric(
            head["therapeutic_coverage_1x_ic50"], metric_reasons,
            "therapeutic_coverage_1x",
            "no TUMOR cells in this reconstruction, so coverage is undefined"),
        "penetration_length_um": _metric(
            head["penetration_length_L_um"], metric_reasons,
            "penetration_length_um",
            "the exp+floor fit of AUC vs distance-from-vessel did not "
            "converge on this field (delivery_metrics.penetration_length_fit); "
            "no substitute value is reported"),
    }
    try:
        l_analytic = float(penetration_length_um(
            dp, ecm_density=float(report["mean_ecm_density"])))
    except Exception:  # noqa: BLE001
        l_analytic = float("nan")

    extra_metrics = {
        "median_tumor_auc": head["median_tumor_auc_uM_s"],
        "mean_tissue_auc": dm["absolute_exposure"]["mean_tissue_auc"],
        "therapeutic_coverage_0_5x": head["therapeutic_coverage_0.5x_ic50"],
        "therapeutic_coverage_2x": head["therapeutic_coverage_2x_ic50"],
        "median_E_over_ic50": head["median_E_over_ic50"],
        "p10_E_over_ic50": head["p10_E_over_ic50"],
        "E_p10_uM": head["E_p10_uM"],
        "E_p50_uM": head["E_p50_uM"],
        "E_p90_uM": head["E_p90_uM"],
        "penetration_length_analytic_um": l_analytic,
        "penetration_index": float(res.extras["penetration_index"]),
        "viable_tumor_fraction": float(res.viable_tumor_fraction),
        "hypoxic_tumor_fraction": float(res.extras["hypoxic_tumor_fraction"]),
        "plasma_mean_uM": float(res.extras["plasma_mean"]),
        "plasma_peak_uM": float(res.extras["plasma_peak"]),
        "n_tumor": int(res.extras["n_tumor"]),
        "exposure_window_s": float(dm["T_s"]),
        "low_o2_over_high_o2_exposure": head["low_o2_over_high_o2_exposure"],
    }
    caveats.append(
        "penetration_index is resolution-sensitive and is OVER-ESTIMATED on "
        "grids coarser than the ~45 um perivascular boundary layer "
        f"(results/numerical_rigor.md section 2); this run used dx = {dx_um:g} um."
        + ("" if dx_um <= 20.0 else
           " kill_fraction is grid-converged only for dx <= 20 um, so this "
           "run's kill fraction is NOT in the converged regime either.")
    )
    if time_resolved:
        caveats.append(
            "Systemic PK was driven TIME-RESOLVED for this run: the real "
            "plasma(t) curve rescales the delivery source at every solver step "
            "(PDEConfig.time_resolved=True), so tissue concentration genuinely "
            "rises and falls. Both paths deliver the same TOTAL drug over the "
            "window, so a difference against a static-source run is a SCHEDULE "
            "effect, not a dose effect. Stored calibrations and regressions in "
            "this repo were produced on the static path, so this run's absolute "
            "numbers are not directly comparable with them."
        )
    else:
        caveats.append(
            "Systemic PK is the documented Stage-0 simplification: a single "
            "time-averaged sustained plasma level drives the transport solve "
            "(tmesim/sim/engine.py), not a time-resolved infusion."
        )
    caveats.append(
        "This tissue is an untreated resection with no outcome label. Every "
        "number here is a PREDICTION of the physical model and has never been "
        "scored against a measured treatment response. NOT FOR CLINICAL USE."
    )

    # --- 6. payloads ------------------------------------------------------- #
    ctx.stage("payloads")
    payloads = _build_payloads(
        micro, report, auc_field=res.fields["auc"], log=ctx.log)

    # --- mechanism attribution -------------------------------------------- #
    mechanism_block, mechanism_reason = _mechanism_block(
        res, attribute_mechanism)
    if mechanism_block is not None:
        payloads["mechanism"] = encode_mechanism(
            mechanism_block.pop("_labels"), mechanism_block)
        rec = mechanism_block["reconciliation"]
        ctx.log(
            f"mechanism: {mechanism_block['counts']['KILLED']} killed / "
            f"{mechanism_block['n_tumor']} tumour cells; reconciles bit-exact="
            f"{rec['bit_exact']}; survival "
            + ", ".join(f"{k} {v:.3f}"
                        for k, v in mechanism_block["survival_breakdown"].items())
        )
        if mechanism_block["counts"]["UNATTRIBUTED"]:
            caveats.append(
                f"MECHANISM: {mechanism_block['counts']['UNATTRIBUTED']} "
                "surviving tumour cells are UNATTRIBUTED — their survival is not "
                "explained by the exposure/resistance/chance bands, and they "
                "have NOT been assigned to the nearest plausible category "
                "(tmesim.sim.mechanism)."
            )

    # --- time-lapse --------------------------------------------------------- #
    timelapse_block = _timelapse_block(
        res, float(micro.domain.dx), n_frames, timelapse_reason, time_resolved)
    if timelapse_block.get("available"):
        payloads["timelapse"] = encode_timelapse(
            res.transport.frames,
            res.transport.frame_times,
            float(micro.domain.dx),
            interval_s=res.transport.frame_interval_s,
        )
        timelapse_block["payload_bytes"] = len(payloads["timelapse"])
        ctx.log(
            f"timelapse: {timelapse_block['n_frames']} frames every "
            f"{timelapse_block['interval_s']:.1f} s over "
            f"{timelapse_block['window_s']:.1f} s "
            f"({len(payloads['timelapse']) / 1e6:.2f} MB)"
        )

    drug_params = {
        "name": dp.name,
        "ic50_uM": float(dp.ic50),
        "hill": float(dp.hill),
        "diffusivity_um2_s": float(dp.diffusivity),
        "binding_1_per_s": float(dp.binding),
        "clearance_1_per_s": float(dp.clearance),
        "bbb_penetrance": float(dp.bbb_penetrance),
        "molecular_weight": float(dp.molecular_weight),
        "mw_from_structure": float(mol.get("mw")) if mol.get("mw") else None,
        "logp": mol.get("logp"),
    }
    return {
        "report": jsonable(report),
        "metrics": metrics,
        "metric_reasons": metric_reasons,
        "extra_metrics": jsonable(extra_metrics),
        "provenance": jsonable(provenance),
        "drug_params": jsonable(drug_params),
        "mechanism": mechanism_block,
        "mechanism_reason": mechanism_reason,
        "timelapse": timelapse_block,
        "_payloads": payloads,
    }


# --------------------------------------------------------------------------- #
# Mechanism + time-lapse helpers
# --------------------------------------------------------------------------- #
def _plan_timelapse(requested: int, grid_shape, total_time_s: float):
    """Decide how many concentration frames to capture, and say why.

    Returns ``(n_frames, reason_or_None)``. The count is clamped to
    ``MAX_TIMELAPSE_FRAMES`` and then reduced, if necessary, to the largest
    number that fits :data:`TIMELAPSE_BUDGET_BYTES`. A capture is never silently
    skipped or silently shortened: whenever the answer differs from what was
    asked for, the reason travels with the run.
    """
    if requested <= 0:
        return 0, "disabled by RunSpec (timelapse_frames <= 0)"
    if requested < 2:
        return 0, (f"timelapse_frames={requested} is below the minimum of 2 "
                   f"(a time-lapse needs a first and a last frame)")

    reason = None
    n = int(requested)
    if n > MAX_TIMELAPSE_FRAMES:
        reason = (f"requested {n} frames, capped at MAX_TIMELAPSE_FRAMES="
                  f"{MAX_TIMELAPSE_FRAMES}")
        n = MAX_TIMELAPSE_FRAMES

    rows, cols = int(grid_shape[0]), int(grid_shape[1])
    per_frame = rows * cols * 4  # f32 on the wire
    if per_frame <= 0:
        return 0, "grid has no cells"
    fits = int(TIMELAPSE_BUDGET_BYTES // per_frame)
    if fits < n:
        if fits < 2:
            return 0, (
                f"a single {rows}x{cols} frame is {per_frame / 1e6:.1f} MB, so "
                f"even two frames exceed the {TIMELAPSE_BUDGET_BYTES / 1e6:.0f} "
                f"MB time-lapse budget; no frames were captured"
            )
        reason = (
            f"requested {requested} frames but {rows}x{cols} is "
            f"{per_frame / 1e6:.2f} MB/frame, so the count was reduced to "
            f"{fits} to stay inside the {TIMELAPSE_BUDGET_BYTES / 1e6:.0f} MB "
            f"budget"
        )
        n = fits

    steps = int(round(float(total_time_s) / PDE_DT_S))
    if n > steps + 1:
        reason = (f"requested {n} frames but the solve takes only {steps} steps; "
                  f"capturing {steps + 1}")
        n = steps + 1
    return n, reason


def _mechanism_block(res, requested: bool):
    """Split the engine's attribution into a JSON block + the label array.

    Returns ``(block_or_None, reason_or_None)``. The block carries the labels
    under the private key ``_labels``; the caller pops them into the binary
    payload so the JSON never duplicates 850k bytes.
    """
    if not requested:
        return None, "disabled by RunSpec (attribute_mechanism=false)"
    mech = getattr(res, "mechanism", None)
    if mech is None:
        return None, (
            "simulate_response returned no attribution; this engine build "
            "predates tmesim.sim.mechanism"
        )
    block = jsonable({k: v for k, v in mech.items() if k != "labels"})
    block["_labels"] = np.asarray(mech["labels"], dtype=np.uint8)
    block["payload"] = "GET /api/runs/{id}/mechanism (TMSM)"
    block["order"] = (
        "one byte per cell, in the same order as the geometry payload "
        "(GET /api/runs/{id}/geometry)"
    )
    return block, None


def _timelapse_block(res, dx_um: float, n_frames: int,
                     plan_reason: Optional[str], time_resolved: bool) -> Dict:
    """Frame metadata for the RunResult (the frames themselves are binary)."""
    transport = getattr(res, "transport", None)
    frames = getattr(transport, "frames", None) if transport is not None else None
    block: Dict = {
        "available": frames is not None and len(frames) > 0,
        "requested_frames": int(n_frames),
        "reason": plan_reason,
        "time_resolved": bool(time_resolved),
        "quantity": "drug concentration C (uM)",
        "quantity_note": (
            "snapshots are the CONCENTRATION field, not the AUC field the "
            "/field/drug_auc overlay serves. AUC only accumulates, so it can "
            "never show drug leaving."
        ),
        "payload": "GET /api/runs/{id}/timelapse (TMSL)",
    }
    if not block["available"]:
        block.setdefault("reason", None)
        if block["reason"] is None:
            block["reason"] = "no frames were captured for this run"
        return block

    times = np.asarray(transport.frame_times, dtype=np.float64)
    block.update({
        "n_frames": int(len(frames)),
        "frame_times_s": [float(t) for t in times],
        "interval_s": float(transport.frame_interval_s or 0.0),
        "realised_interval_s": (
            float(np.diff(times).mean()) if times.size > 1 else 0.0),
        "window_s": float(times[-1]) if times.size else 0.0,
        "grid": [int(frames.shape[1]), int(frames.shape[2])],
        "dx_um": float(dx_um),
        "max_concentration_uM": float(np.nanmax(frames)) if frames.size else None,
    })
    if not time_resolved:
        block["plateau_note"] = (
            "this run used the static (window-averaged) plasma source, so the "
            "delivery rate is held constant for the whole window by "
            "construction: the drug is seen to ENTER and SPREAD, and it "
            "approaches a quasi-steady state rather than clearing. Set "
            "time_resolved=true in the RunSpec to drive the solver with the "
            "real plasma(t) curve and see it fall."
        )
    return block


def _build_payloads(micro, report: Dict, auc_field, log: Callable[[str], None]
                    ) -> Dict[str, bytes]:
    """Encode the geometry + every available field into their wire payloads."""
    cells = micro.cells
    extent = report.get("extent_um") or micro.domain.extent_um
    out: Dict[str, bytes] = {
        "geometry": encode_geometry(cells.x, cells.y, cells.cell_type, extent),
    }
    dx = float(micro.domain.dx)
    sources = {
        "vascular_density": micro.fields.vascular_density,
        "ecm_density": micro.fields.ecm_density,
        "oxygen": micro.fields.oxygen,
        "drug_auc": auc_field,
    }
    for name in FIELD_NAMES:
        field = sources.get(name)
        if field is None:
            continue
        out[f"field_{name}"] = encode_field(np.asarray(field), dx)
    log(f"payloads: geometry {len(out['geometry'])} B for {cells.n} cells; "
        f"fields {[k for k in out if k != 'geometry']}")
    return out


# --------------------------------------------------------------------------- #
# Job + manager
# --------------------------------------------------------------------------- #
class Job:
    """One queued / running / finished run."""

    def __init__(self, run_id: str, spec: Dict):
        self.run_id = run_id
        self.spec = spec
        self.status = "queued"
        self.stage = "queued"
        self.progress = 0.0
        self.created_at = time.time()
        self.started_at: Optional[float] = None
        self.finished_at: Optional[float] = None
        self.error: Optional[str] = None
        self.traceback: Optional[str] = None
        self.body: Optional[Dict] = None
        self.cancel_event = threading.Event()
        self._log: Deque[str] = deque(maxlen=_LOG_MAXLEN)
        self._lock = threading.Lock()
        self.log(f"queued: {json.dumps(spec, sort_keys=True, default=str)}")

    # -- logging ----------------------------------------------------------- #
    def log(self, msg: str) -> None:
        stamp = time.strftime("%H:%M:%S")
        with self._lock:
            self._log.append(f"{stamp} {msg}")

    def log_tail(self, n: int = 40) -> List[str]:
        with self._lock:
            return list(self._log)[-int(n):]

    # -- views ------------------------------------------------------------- #
    @property
    def elapsed_s(self) -> float:
        if self.started_at is None:
            return 0.0
        end = self.finished_at if self.finished_at is not None else time.time()
        return float(end - self.started_at)

    def status_view(self, log_lines: int = 40) -> Dict:
        return {
            "run_id": self.run_id,
            "status": self.status,
            "progress": round(float(self.progress), 4),
            "stage": self.stage,
            "elapsed_s": round(self.elapsed_s, 3),
            "log_tail": self.log_tail(log_lines),
        }

    def result_view(self) -> Dict:
        out = {
            "run_id": self.run_id,
            "spec": self.spec,
            "status": self.status,
            "report": None,
            "metrics": None,
            "provenance": None,
            "runtime_s": round(self.elapsed_s, 3) if self.started_at else None,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "stage": self.stage,
            "progress": round(float(self.progress), 4),
            "error": self.error,
        }
        if self.body:
            for k, v in self.body.items():
                if not k.startswith("_"):
                    out[k] = v
        return out


class JobManager:
    """Queue, worker thread, and on-disk history under ``app/runs/``.

    Completed runs (result JSON + every binary payload) are written to
    ``app/runs/<run_id>/`` so history — including the geometry and field
    payloads the UI renders — survives a server restart.
    """

    def __init__(self, root: Optional[str] = None,
                 cache: ReconstructionCache = RECON_CACHE,
                 start_worker: bool = True):
        self.root = root or runs_dir()
        os.makedirs(self.root, exist_ok=True)
        self.cache = cache
        self._jobs: "Dict[str, Job]" = {}
        self._order: List[str] = []          # submission order, oldest first
        self._lock = threading.RLock()
        self._queue: "queue.Queue[str]" = queue.Queue()
        self._counter = 0
        self._worker: Optional[threading.Thread] = None
        self._stop = threading.Event()
        self._restore()
        if start_worker:
            self.start()

    # -- lifecycle ---------------------------------------------------------- #
    def start(self) -> None:
        if self._worker is not None and self._worker.is_alive():
            return
        self._stop.clear()
        self._worker = threading.Thread(
            target=self._run_worker, name="tmesim-run-worker", daemon=True)
        self._worker.start()

    def shutdown(self, timeout: float = 5.0) -> None:
        self._stop.set()
        self._queue.put("")  # unblock the worker
        if self._worker is not None:
            self._worker.join(timeout=timeout)

    # -- persistence -------------------------------------------------------- #
    def _dir(self, run_id: str) -> str:
        return os.path.join(self.root, run_id)

    def _restore(self) -> None:
        """Load previously-completed runs from disk into the in-memory index."""
        try:
            names = sorted(os.listdir(self.root))
        except OSError:
            return
        for name in names:
            path = os.path.join(self.root, name, "result.json")
            if not (name.startswith("r_") and os.path.exists(path)):
                continue
            try:
                with open(path) as fh:
                    stored = json.load(fh)
            except Exception:  # noqa: BLE001 - skip an unreadable run
                continue
            job = Job(stored.get("run_id", name), stored.get("spec") or {})
            job.status = stored.get("status", "done")
            job.stage = stored.get("stage", job.status)
            job.progress = float(stored.get("progress") or
                                 (1.0 if job.status == "done" else 0.0))
            job.created_at = stored.get("created_at") or job.created_at
            job.started_at = stored.get("started_at")
            job.finished_at = stored.get("finished_at")
            job.error = stored.get("error")
            job.body = {k: v for k, v in stored.items()
                        if k in ("report", "metrics", "metric_reasons",
                                 "extra_metrics", "provenance", "drug_params",
                                 "mechanism", "mechanism_reason", "timelapse",
                                 # a TPL sweep's grid, rows and provenance —
                                 # without this a restored protocol run would
                                 # come back with an empty result grid
                                 "protocol")}
            job._log.clear()
            for line in (stored.get("log_tail") or []):
                job._log.append(line)
            job.log("restored from disk (this server was restarted)")
            with self._lock:
                self._jobs[job.run_id] = job
                self._order.append(job.run_id)
                try:
                    self._counter = max(
                        self._counter, int(job.run_id.split("_")[-1]))
                except ValueError:
                    pass

    def _persist(self, job: Job, payloads: Optional[Dict[str, bytes]],
                 status: Optional[str] = None) -> None:
        """Write every binary payload, then the result document, atomically.

        ``status`` stamps the document with a terminal state the in-memory job
        has not been given yet. That is what lets :meth:`_execute` finish ALL
        the writes before it announces ``done`` — see the note there.
        """
        d = self._dir(job.run_id)
        os.makedirs(d, exist_ok=True)
        for name, blob in (payloads or {}).items():
            tmp = os.path.join(d, name + ".bin.tmp")
            with open(tmp, "wb") as fh:
                fh.write(blob)
            os.replace(tmp, os.path.join(d, name + ".bin"))
        doc = job.result_view()
        if status is not None:
            doc["status"] = status
            doc["stage"] = status
            doc["progress"] = 1.0
        doc["log_tail"] = job.log_tail(_LOG_MAXLEN)
        tmp = os.path.join(d, "result.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(jsonable(doc), fh, indent=1)
        os.replace(tmp, os.path.join(d, "result.json"))

    def payload_path(self, run_id: str, name: str) -> str:
        """Absolute path of a persisted binary payload (may not exist)."""
        return os.path.join(self._dir(run_id), f"{name}.bin")

    # -- queue -------------------------------------------------------------- #
    def submit(self, spec: Dict) -> Job:
        with self._lock:
            self._counter += 1
            run_id = f"r_{self._counter:05d}"
            job = Job(run_id, spec)
            self._jobs[run_id] = job
            self._order.append(run_id)
        self._queue.put(run_id)
        return job

    def get(self, run_id: str) -> Optional[Job]:
        with self._lock:
            return self._jobs.get(run_id)

    def list(self, limit: Optional[int] = None) -> List[Job]:
        """Newest first."""
        with self._lock:
            ids = list(reversed(self._order))
            jobs = [self._jobs[i] for i in ids if i in self._jobs]
        return jobs[:limit] if limit else jobs

    def cancel(self, run_id: str) -> Optional[Job]:
        """Cancel a queued or running job.

        A queued job is cancelled immediately. A running job is cancelled
        COOPERATIVELY: the flag is checked at every stage boundary, so a job
        already inside ``simulate_response`` finishes that solve first — the
        engine is not modified to accept an interrupt.
        """
        job = self.get(run_id)
        if job is None or job.status in TERMINAL:
            return job
        job.cancel_event.set()
        job.log("cancel requested")
        if job.status == "queued":
            job.status = "cancelled"
            job.stage = "cancelled"
            job.finished_at = time.time()
            job.started_at = job.started_at or job.finished_at
            job.error = "cancelled before it started running"
            self._persist(job, None)
        return job

    # -- worker ------------------------------------------------------------- #
    def _run_worker(self) -> None:
        while not self._stop.is_set():
            try:
                run_id = self._queue.get(timeout=0.25)
            except queue.Empty:
                continue
            if not run_id:
                continue
            job = self.get(run_id)
            if job is None or job.status in TERMINAL:
                continue
            if job.cancel_event.is_set():
                job.status = "cancelled"
                job.stage = "cancelled"
                job.started_at = job.started_at or time.time()
                job.finished_at = time.time()
                job.error = "cancelled before it started running"
                self._persist(job, None)
                continue
            self._execute(job)

    def _execute(self, job: Job) -> None:
        job.status = "running"
        job.started_at = time.time()
        job.log("running")
        ctx = _Ctx(job)
        try:
            body = execute_run(job.spec, ctx, cache=self.cache)
            payloads = body.pop("_payloads", None)
            job.body = body
            job.finished_at = time.time()
            job.log(f"done in {job.elapsed_s:.1f}s")
            # Persist BEFORE announcing done. ``done`` is the signal the UI (and
            # a restarting server) treat as "the payloads are there": the poller
            # sees it and immediately GETs /geometry, /mechanism, /timelapse. If
            # the flag were flipped first, every one of those would race the
            # writes and 409 on a run that had in fact succeeded — a window that
            # was narrow while a run wrote ~1 MB and is not narrow now that a
            # time-lapse is tens of megabytes. The document is stamped with the
            # terminal status it is about to be given, so what lands on disk is
            # a finished run, not a snapshot mid-write.
            self._persist(job, payloads, status="done")
            job.status = "done"
            job.stage = "done"
            job.progress = 1.0
        except RunCancelled as exc:
            job.status = "cancelled"
            job.stage = "cancelled"
            job.finished_at = time.time()
            job.error = str(exc)
            job.log(f"cancelled: {exc}")
            self._persist(job, None)
        except Exception as exc:  # noqa: BLE001 - surfaced, never swallowed
            job.status = "failed"
            job.stage = "failed"
            job.finished_at = time.time()
            job.error = f"{type(exc).__name__}: {exc}"
            job.traceback = traceback.format_exc()
            job.log(f"FAILED: {job.error}")
            for line in job.traceback.rstrip().splitlines()[-12:]:
                job.log(line)
            self._persist(job, None)
