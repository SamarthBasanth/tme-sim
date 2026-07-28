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

from .cache import RECON_CACHE, ReconstructionCache, make_key, runs_dir
from .payloads import FIELD_NAMES, encode_field, encode_geometry

__all__ = [
    "PDE_DT_S",
    "PDE_SAVE_EVERY",
    "SIM_SEED",
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
    """
    from tmesim.core import PDEConfig, PKConfig, SimConfig
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
    mgmt = spec.get("mgmt")
    pk = spec.get("pk") or {}
    total_time_s = float(spec["total_time_s"])

    caveats: List[str] = []
    metric_reasons: Dict[str, str] = {}

    # --- 1. reconstruction (cached) --------------------------------------- #
    ctx.stage("reconstruct")
    key = make_key(section, dx_um, vasc_mode, oxy_mode)
    was_cached = cache.peek(key) is not None
    entry = cache.get_or_build(key, log=ctx.log)
    micro = entry["micro"]
    caveats.extend(entry.get("notes") or [])

    # Reset the shared cached microenvironment to the state reconstruct_real
    # left it in, so a previous run's per-cell sensitivity can never leak into
    # this one. (The worker is single-threaded, so this is a full reset.)
    micro.cells.sensitivity = np.array(entry["pristine_sensitivity"], copy=True)
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
        "pde_dt_s": PDE_DT_S,
        "pde_save_every": PDE_SAVE_EVERY,
        "pde_convection": PDE_CONVECTION,
        # PDEConfig.time_resolved (whether the source follows plasma(t) or its
        # window average) is not in the frozen RunSpec, so this backend always
        # runs the engine's own default. Read off the dataclass rather than
        # retyped, and reported as None if this engine build has no such field.
        "pde_time_resolved": getattr(PDEConfig(), "time_resolved", None),
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
        return {
            "report": jsonable(report),
            "metrics": {k: None for k in metric_reasons},
            "metric_reasons": metric_reasons,
            "provenance": jsonable(provenance),
            "drug_params": None,
            "_payloads": payloads,
        }

    pot = resolve_ic50(drug_name)
    mol = get_molecule(drug_name)
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
                    convection=PDE_CONVECTION, save_every=PDE_SAVE_EVERY)
    pk_cfg = PKConfig(
        plasma_peak=float(pk.get("plasma_peak", PKConfig().plasma_peak)),
        half_life=float(pk.get("half_life_s", PKConfig().half_life)),
        n_doses=int(pk.get("n_doses", PKConfig().n_doses)),
        dose_interval=float(pk.get("dose_interval_s", PKConfig().dose_interval)),
    )
    cfg = SimConfig(drug=dp, pk=pk_cfg, pde=pde, seed=SIM_SEED)
    ctx.log(f"simulate_response(dt={PDE_DT_S}s, total_time={total_time_s}s, "
            f"convection={PDE_CONVECTION}, seed={SIM_SEED}, "
            f"plasma_peak={pk_cfg.plasma_peak} uM, "
            f"half_life={pk_cfg.half_life} s, n_doses={pk_cfg.n_doses})")
    t_sim = time.time()
    res = simulate_response(micro, cfg)
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
        "_payloads": payloads,
    }


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
                                 "extra_metrics", "provenance", "drug_params")}
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

    def _persist(self, job: Job, payloads: Optional[Dict[str, bytes]]) -> None:
        d = self._dir(job.run_id)
        os.makedirs(d, exist_ok=True)
        for name, blob in (payloads or {}).items():
            tmp = os.path.join(d, name + ".bin.tmp")
            with open(tmp, "wb") as fh:
                fh.write(blob)
            os.replace(tmp, os.path.join(d, name + ".bin"))
        doc = job.result_view()
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
            job.status = "done"
            job.stage = "done"
            job.progress = 1.0
            job.finished_at = time.time()
            job.log(f"done in {job.elapsed_s:.1f}s")
            self._persist(job, payloads)
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
