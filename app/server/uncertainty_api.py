"""``/api/uncertainty/*`` — the interval, as a first-class output.

===========  ==========================================  =====================
Method       Path                                        Returns
===========  ==========================================  =====================
GET          /api/uncertainty/spans                      the recorded-span
                                                         registry + the two
                                                         KINDS, with provenance
GET          /api/runs/{run_id}/uncertainty              that finished run's
                                                         intervals and refusals
===========  ==========================================  =====================

Why this router exists
----------------------
``results/kill_calibration.md`` section 5.1 says it plainly — *"Absolute kill is
defensible as an INTERVAL, not as a number, and the interval is wide"* —
temozolomide 0.0048-0.1442, a 30x span. Until now that sentence lived only in a
prose caveat in a results document, while every surface a reader actually looks
at reported a single-valued ``kill_fraction``. This router puts the interval
where the number is.

**It does not make absolute kill fraction a defensible output.** Every interval
comes back with ``quotable: false`` and the reason attached. The interval is
published as the SIZE OF THE REFUSAL: a reader who sees ``kill_fraction =
0.1122`` and nothing else has no way to know that the characterised uncertainty
around it spans a factor of thirty, and that is the failure this fixes.

**The two bars are different KINDS and the payload keeps them apart.** The
potency bar is BIOLOGICAL — patients really do differ that much, and no better
model narrows it; it narrows only by measuring this patient's own sensitivity.
The cell-state bar is EPISTEMIC — it is our ignorance, and one Ki-67 stain on a
serial section deletes it. A payload that merged them would tell a reader the
opposite of the truth about both, so
:func:`tmesim.validation.uncertainty.assert_kinds_separated` runs over every
response before it is returned.

**No cross-compound anything.** Both routes take at most ONE run for ONE
compound. The ranking was withdrawn on 2026-08-21 and is not reinstated here in
any form; there is no entry point that accepts two compounds.

How a finished run is re-read
-----------------------------
The band re-evaluates the run's own kill at the recorded IC50 band edges, which
needs the run's per-cell exposure and per-cell state. Both are recovered from
the payloads the run already persisted — ``geometry.bin`` (positions + types),
``field_drug_auc.bin`` (the exposure field) and ``field_oxygen.bin`` (what the
state bands are cut from). **Nothing is re-simulated and no PDE is solved.**

That recovery is then GATED, not trusted: the recomputed expectation must
reproduce the run's own recorded ``kill_fraction`` to within
:data:`~tmesim.validation.uncertainty.REPRODUCTION_SIGMA_GATE` sigma of the
engine's Bernoulli sampling scale, or the band is refused with the residual
quoted. A recovery that cannot reproduce the run is describing a different run,
and a band around the wrong number is worse than no band.

Refusals are ``200``, not errors. "No measured span exists for this metric" is a
*result* — it is the honest answer and the one a UI needs in order to show a
bare number without implying it is tight. Only an unknown run (404) or an
unfinished one (409) is an error.
"""
from __future__ import annotations

import os
import re
import time
from functools import lru_cache
from typing import Dict, Tuple

import numpy as np
from fastapi import APIRouter, HTTPException

from .cache import redact_exception
from .jobs import JobManager
from .payloads import decode_field, decode_geometry

from tmesim.validation.uncertainty import (
    BANDABLE_GRADES,
    COHORT_MEDIAN_SPAN,
    F_STATE_BY_ROUTE,
    F_STATE_DEFAULT_ROUTE,
    KIND_MEANING,
    NOT_COVERED,
    POTENCY_SPAN,
    REPRODUCTION_SIGMA_GATE,
    STATE_ROUTE_SPAN,
    STATE_SPAN_FACTOR,
    STATE_SPAN_SECTIONS,
    UNCHARACTERISED_METRICS,
    FinishedRun,
    assert_kinds_separated,
    run_interval_report,
)

__all__ = ["build_uncertainty_router", "load_finished_run", "RunNotBandable"]

#: Run ids this server mints (``jobs.JobManager.submit``). Anything else is
#: rejected before it can reach a filesystem join — ``jobs.get()`` already only
#: knows ids it minted, but a path is built from this string and a second,
#: cheaper guard costs nothing.
_RUN_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")

#: Only ``sensitivity_mode='uniform'`` is recoverable from the persisted
#: payloads: under any other mode the per-cell sensitivity is a marker-derived
#: array that geometry.bin does not carry. Refused rather than approximated.
_RECOVERABLE_SENSITIVITY_MODES = ("uniform",)


@lru_cache(maxsize=64)
def _drug_constants(compound: str, alpha: float) -> Tuple[float, float]:
    """``(hill, effective_hypoxia_resistance)`` for one compound at one alpha.

    MEASURED: ``calibration.calibrated_drugparams`` costs ~0.55 s per call — it
    recompiles the molecule — and it was 100% of this route's cost, against
    ~0.3 ms to decode all three payloads. Only two floats are ever read out of
    the result, so only two floats are cached: handing a shared mutable
    ``DrugParams`` back to callers would be a footgun for a saving this does not
    need.
    """
    from tmesim.cells.pkpd import effective_hypoxia_resistance
    from tmesim.core.config import PDConfig
    from tmesim.validation.calibration import calibrated_drugparams

    drug = calibrated_drugparams(compound, float(alpha))
    return (float(drug.hill),
            float(effective_hypoxia_resistance(drug, PDConfig())))


class RunNotBandable(Exception):
    """A finished run whose per-cell state cannot be recovered from disk.

    Carries a ``detail`` dict that is returned to the client verbatim as a
    refusal — it names what is missing, never a guess at a value.
    """

    def __init__(self, detail: Dict) -> None:
        super().__init__(str(detail.get("reason") or "not bandable"))
        self.detail = detail


# --------------------------------------------------------------------------- #
# Recovering one finished run's per-cell state from its persisted payloads
# --------------------------------------------------------------------------- #
def _read(path: str) -> bytes:
    with open(path, "rb") as fh:
        return fh.read()


def _cell_indices(x, y, dx: float, ny: int, nx: int) -> Tuple[np.ndarray, np.ndarray]:
    """``tmesim.core.Domain.to_index`` on the persisted geometry.

    Reproduced rather than imported because the domain object is gone once the
    run finished; it is three lines and the reproduction gate checks the result.
    """
    col = np.clip((np.asarray(x, dtype=np.float64) / dx).astype(int), 0, nx - 1)
    row = np.clip((np.asarray(y, dtype=np.float64) / dx).astype(int), 0, ny - 1)
    return col, row


def load_finished_run(jobs: JobManager, run_id: str, doc: Dict) -> FinishedRun:
    """Rebuild the per-cell (exposure, modulation) arrays of a finished run.

    Raises :class:`RunNotBandable` — never a fabricated fallback — when a
    payload is missing, a field's grid disagrees with another's, or the run used
    a configuration whose per-cell inputs the payloads do not carry.
    """
    from tmesim.core import CellState, CellType
    from tmesim.core.config import OxygenConfig, PDConfig

    spec = doc.get("spec") or {}
    prov = doc.get("provenance") or {}
    metrics = doc.get("metrics") or {}

    compound = spec.get("drug")
    if not compound:
        raise RunNotBandable({
            "status": "REFUSED",
            "reason": (
                "this is a reconstruction-only run (spec.drug is null): it has "
                "no drug, no exposure field and no kill fraction, so there is "
                "nothing to put an interval around."),
        })

    # NOT ``metrics.get(...) or {}``: a kill_fraction of exactly 0.0 is a real
    # measured value (a compound that killed nothing is a finding, not a
    # missing metric) and ``or`` would swallow it into the refusal branch.
    kill = metrics.get("kill_fraction")
    kill_value = kill.get("value") if isinstance(kill, dict) else kill
    if kill_value is None or not np.isfinite(float(kill_value)):
        raise RunNotBandable({
            "status": "REFUSED",
            "reason": ("this run reports no finite kill_fraction, so there is "
                       "no point for an interval to bracket."),
        })

    sens_mode = str(prov.get("sensitivity_mode") or "uniform")
    if sens_mode not in _RECOVERABLE_SENSITIVITY_MODES:
        raise RunNotBandable({
            "status": "REFUSED",
            "reason": (
                f"sensitivity_mode={sens_mode!r}: per-cell sensitivity is a "
                f"marker-derived array that the persisted geometry payload does "
                f"not carry, so the run's own modulation cannot be recovered "
                f"from disk. Re-running with the interval computed in-process "
                f"would band it; approximating the array here would not."),
        })

    needed = {"geometry": "geometry", "auc": "field_drug_auc",
              "oxygen": "field_oxygen"}
    paths = {k: jobs.payload_path(run_id, v) for k, v in needed.items()}
    missing = sorted(k for k, p in paths.items() if not os.path.exists(p))
    if missing:
        raise RunNotBandable({
            "status": "REFUSED",
            "reason": (
                f"payload(s) {missing} are not on disk for this run (pruned, or "
                f"written by a build that predates them). The band re-evaluates "
                f"this run's OWN per-cell exposure and state; without them there "
                f"is nothing to re-evaluate, and substituting another run's "
                f"fields would produce a band around a different tumour."),
            "missing_payloads": missing,
        })

    geo = decode_geometry(_read(paths["geometry"]))
    auc = decode_field(_read(paths["auc"]))
    oxy = decode_field(_read(paths["oxygen"]))
    if (auc["rows"], auc["cols"]) != (oxy["rows"], oxy["cols"]) or \
            abs(auc["dx_um"] - oxy["dx_um"]) > 1e-9:
        raise RunNotBandable({
            "status": "REFUSED",
            "reason": (
                f"the exposure field ({auc['rows']}x{auc['cols']} @ "
                f"{auc['dx_um']} um) and the oxygen field ({oxy['rows']}x"
                f"{oxy['cols']} @ {oxy['dx_um']} um) are not on the same grid, "
                f"so a cell cannot be given both an exposure and a state."),
        })

    dx = float(auc["dx_um"])
    ny, nx = int(auc["rows"]), int(auc["cols"])
    col, row = _cell_indices(geo["x"], geo["y"], dx, ny, nx)
    auc_cell = np.asarray(auc["values"], dtype=np.float64)[row, col]
    o2_cell = np.asarray(oxy["values"], dtype=np.float64)[row, col]

    duration = float(spec.get("total_time_s") or 0.0)
    if duration <= 0.0:
        raise RunNotBandable({
            "status": "REFUSED",
            "reason": ("spec.total_time_s is missing or non-positive, so the "
                       "AUC field cannot be reduced to the time-averaged "
                       "concentration the Hill term takes."),
        })
    E_uM = auc_cell / duration

    # --- state banding, exactly as cells.agents.update_cell_states cuts it ---
    thr = float(OxygenConfig().hypoxia_threshold)
    state = np.full(o2_cell.size, int(CellState.VIABLE), dtype=np.int64)
    state[o2_cell < thr] = int(CellState.HYPOXIC)
    state[(o2_cell >= thr) & (o2_cell < 2.0 * thr)] = int(CellState.QUIESCENT)

    pd_cfg = PDConfig()
    hill, r_hypoxia = _drug_constants(str(compound), float(prov.get("alpha")))
    resistance = np.ones(state.size, dtype=np.float64)
    if not bool(prov.get("resistance_off")):
        resistance[state == int(CellState.HYPOXIC)] = r_hypoxia
        resistance[state == int(CellState.QUIESCENT)] = float(
            pd_cfg.quiescent_resistance)
    sensitivity = float(prov.get("sensitivity_mean") or 1.0)

    tumor = np.asarray(geo["cell_type"]) == int(CellType.TUMOR)
    if not tumor.any():
        raise RunNotBandable({
            "status": "REFUSED",
            "reason": "this reconstruction has no TUMOR cells to band.",
        })

    ic50 = prov.get("ic50_resolved_uM", prov.get("ic50_uM"))
    return FinishedRun(
        compound=str(compound),
        section_id=str(spec.get("section") or ""),
        kill_fraction=float(kill_value),
        E_uM=E_uM[tumor],
        modulation=(sensitivity * resistance)[tumor],
        ic50_uM=float(ic50),
        hill=float(hill),
        max_kill_rate=float(pd_cfg.max_kill_rate),
        state_mode=str(spec.get("state_mode") or "oxygen"),
        dx_um=dx,
        run_id=run_id,
    )


# --------------------------------------------------------------------------- #
# The span registry, as a response
# --------------------------------------------------------------------------- #
def span_registry() -> Dict:
    """Every recorded span this build can band with, and what kind each one is.

    Cheap — it touches no run and no data cache — so a client can show the two
    kinds and their provenance before any simulation is paid for.
    """
    return {
        "kinds": KIND_MEANING,
        "the_distinction": (
            "The two bars are different KINDS and conflating them would make "
            "this feature worse than no feature. The potency bar is largely "
            "BIOLOGICAL variation between patients: it does not shrink with a "
            "better model, it shrinks by measuring this patient's own "
            "sensitivity. The cell-state bar is EPISTEMIC: it is our ignorance, "
            "and a better model — or one Ki-67 stain on a serial section — "
            "removes it."
        ),
        "potency_spans": {k: v.as_dict() for k, v in POTENCY_SPAN.items()},
        "cohort_median_spans": {k: v.as_dict()
                                for k, v in COHORT_MEDIAN_SPAN.items()},
        "cell_state_route_span": STATE_ROUTE_SPAN.as_dict(),
        "cell_state_route_applies_to_sections": list(STATE_SPAN_SECTIONS),
        "cell_state_route_default": F_STATE_DEFAULT_ROUTE,
        "f_state_by_route": dict(F_STATE_BY_ROUTE),
        "cell_state_factor": {"lo": STATE_SPAN_FACTOR[0],
                              "hi": STATE_SPAN_FACTOR[1]},
        "bandable_grades": list(BANDABLE_GRADES),
        "uncharacterised_metrics": dict(UNCHARACTERISED_METRICS),
        "not_covered": list(NOT_COVERED),
        "reproduction_gate_sigma": REPRODUCTION_SIGMA_GATE,
        "cross_compound": (
            "REFUSED. No route here accepts more than one compound, and no "
            "response ranks, orders or compares compounds in any form "
            "(withdrawn 2026-08-21)."
        ),
        "record": "results/uncertainty_band_2026-09-20.md",
    }


# --------------------------------------------------------------------------- #
# Router
# --------------------------------------------------------------------------- #
def build_uncertainty_router(jobs: JobManager) -> APIRouter:
    """``/api/uncertainty/*`` bound to the SAME job manager as ``/api/runs``.

    It owns no queue, no cache and no second copy of a span: the registry is
    read straight out of :mod:`tmesim.validation.uncertainty`, so a span cannot
    be corrected in the library and stay stale on the wire.
    """
    router = APIRouter()

    @router.get("/api/uncertainty/spans")
    def get_spans() -> Dict:
        return span_registry()

    @router.get("/api/runs/{run_id}/uncertainty")
    def get_run_uncertainty(run_id: str) -> Dict:
        t0 = time.perf_counter()
        if not _RUN_ID_RE.match(run_id or ""):
            raise HTTPException(status_code=404, detail={
                "message": "unknown run", "run_id": run_id,
                "reason": "run id is not in this server's id format"})
        job = jobs.get(run_id)
        if job is None:
            raise HTTPException(status_code=404, detail={
                "message": "unknown run", "run_id": run_id})
        if job.status != "done":
            raise HTTPException(status_code=409, detail={
                "message": "this run has not finished",
                "run_id": run_id, "status": job.status, "stage": job.stage,
                "reason": ("the interval re-reads the run's persisted per-cell "
                           "exposure and state; those are written when a run "
                           "completes."),
            })

        doc = job.result_view()
        base = {
            "run_id": run_id,
            "compound": (doc.get("spec") or {}).get("drug"),
            "section_id": (doc.get("spec") or {}).get("section"),
            "kinds": KIND_MEANING,
            "not_covered": list(NOT_COVERED),
            "record": "results/uncertainty_band_2026-09-20.md",
        }
        try:
            run = load_finished_run(jobs, run_id, doc)
        except RunNotBandable as exc:
            # A refusal is a RESULT, not a server error: 200 with the reason.
            out = dict(base)
            out.update(exc.detail)
            out["intervals"] = {}
            out["cost_s"] = time.perf_counter() - t0
            return out
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(status_code=500, detail={
                "message": "could not re-read this run's payloads",
                "run_id": run_id,
                "reason": redact_exception(
                    exc, f"uncertainty band for {run_id}"),
            }) from None

        point_metrics = {
            name: (m.get("value") if isinstance(m, dict) else m)
            for name, m in (doc.get("metrics") or {}).items()
        }
        point_metrics.update({
            k: v for k, v in (doc.get("extra_metrics") or {}).items()
            if isinstance(v, (int, float))
        })
        report = run_interval_report(run, point_metrics={
            k: v for k, v in point_metrics.items()
            if isinstance(v, (int, float)) and np.isfinite(float(v))
        })
        for interval in report["intervals"].values():
            if interval.get("band") is not None:
                # The guard runs on the way OUT, over the object the client will
                # actually receive — not over an intermediate that a later edit
                # could stop being the same object.
                assert_kinds_separated(interval)
        out = dict(base)
        out.update(report)
        out["recomputed_from"] = [
            "geometry.bin", "field_drug_auc.bin", "field_oxygen.bin"]
        out["resimulated"] = False
        out["cost_s"] = time.perf_counter() - t0
        return out

    return router
