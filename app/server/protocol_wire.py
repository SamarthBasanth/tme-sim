"""Editor-facing projection of a TPL plan/result.

Two clients read these routes and they were specified independently:

* **API_CONTRACT.md** (this backend's own contract) — ``diagnostics[]`` with
  **1-based** ``col`` / ``end_col``, ``n_runs``, ``estimated_cost``, ``grid``,
  ``rows``.
* **``app/web/js/protocol/api.js``** (the editor's frozen wire doc) —
  ``errors[]`` / ``warnings[]`` with **0-based, half-open** ``col`` /
  ``end_col``, ``run_count``, ``estimated_cost_s``, ``axes``, ``refusal``, and a
  sweep snapshot under ``/api/protocol/sweeps/{id}``.

Rather than pick one and break the other, every response carries **both**, under
disjoint key names, and this module owns the second projection. The two column
conventions never collide because they never share a key: ``diagnostics[].col``
is 1-based, ``errors[].col`` is 0-based, and both say so here and in
API_CONTRACT.md. If the two are ever reconciled, deleting this module and its
three call sites is the whole change.

Nothing here computes anything new — it is a pure reshaping of
:class:`tmesim.protocol.executor.Plan` and the executor's result.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

__all__ = [
    "editor_diagnostics", "editor_validate", "editor_plan", "sweep_snapshot",
    "estimate_cost_s", "RECON_S_PER_MCELL", "SOLVE_S_PER_MSTEPCELL",
]

# --------------------------------------------------------------------------- #
# Cost estimate.
#
# The contract's `estimated_cost.wall_time_s` is deliberately null — this
# backend does not benchmark in advance. The editor's panel wants a number, so
# these two constants exist to produce one, and they are LABELLED an estimate
# everywhere they surface (`cost_note`). They are back-computed from one
# measured run on this repository:
#
#   xenium_319, dx = 30 um -> 170 x 271 = 45,630 cells; reconstruction 3.6 s
#   24 h window at dt = 10 s -> 8,640 steps; simulate + metrics ~110 s
#
# They are a property of THIS MACHINE and this section, not of the model, and a
# different host will be off by whatever its single-core numpy/scipy throughput
# differs by. That is exactly why the contract's own field stays null.
# --------------------------------------------------------------------------- #
#: seconds of reconstruction per million grid cells (measured, see above)
RECON_S_PER_MCELL = 3.6 / 0.045630
#: seconds of transport solve per million (timestep x grid cell)
SOLVE_S_PER_MSTEPCELL = 110.0 / (8640 * 45630 / 1e6)


def estimate_cost_s(plan) -> "tuple[float, str]":
    """Rough wall time for a plan, in seconds, plus the sentence that qualifies it."""
    cost = plan.estimated_cost or {}
    cells = cost.get("approx_grid_cells") or {}
    mean_cells = (sum(cells.values()) / len(cells)) if cells else 0.0
    steps = float(cost.get("total_pde_timesteps") or 0.0)
    n_recon = int(cost.get("reconstructions_to_build") or 0)
    total = (n_recon * RECON_S_PER_MCELL * mean_cells / 1e6
             + SOLVE_S_PER_MSTEPCELL * steps * mean_cells / 1e6)
    note = (
        f"ESTIMATE, not a measurement. {int(steps)} timesteps at dt = "
        f"{cost.get('pde_dt_s')} s over ~{int(mean_cells)} grid cells, plus "
        f"{n_recon} reconstruction(s) still to build, scaled by throughput "
        "constants back-computed from ONE run on this repository "
        "(app/server/protocol_wire.py). It is a property of this machine, not "
        "of the model; the contract's own estimated_cost.wall_time_s stays "
        "null for that reason."
    )
    return (round(total, 1) if math.isfinite(total) else 0.0), note


# --------------------------------------------------------------------------- #
# Diagnostics
# --------------------------------------------------------------------------- #
def editor_diagnostics(diags: List[Dict]) -> "tuple[List[Dict], List[Dict]]":
    """Split 1-based diagnostics into the editor's 0-based errors / warnings."""
    errors: List[Dict] = []
    warnings: List[Dict] = []
    for d in diags:
        out: Dict = {
            "line": d["line"] or 1,
            "message": d["message"],
            "severity": d["severity"],
            "code": d["code"],
        }
        if d.get("col"):
            out["col"] = max(int(d["col"]) - 1, 0)
            out["end_col"] = max(int(d.get("end_col") or d["col"]) - 1,
                                 out["col"] + 1)
        if d.get("hint"):
            out["remedy"] = d["hint"]
        if d.get("remedy_text"):
            out["remedy_text"] = d["remedy_text"]
            out["remedy_line"] = d.get("remedy_line") or out["line"]
        (errors if d["severity"] == "error" else warnings).append(out)
    return errors, warnings


def editor_validate(body: Dict, plan) -> Dict:
    """Add the editor's `errors` / `warnings` / `symbols` to a /validate body."""
    errors, warnings = editor_diagnostics(body.get("diagnostics") or [])
    prog = getattr(plan, "program", None)
    body["errors"] = errors
    body["warnings"] = warnings
    body["symbols"] = {
        "drugs": [d.name for d in (prog.drugs if prog else [])],
        "protocols": [p.name for p in (prog.protocols if prog else [])],
        "sweeps": [s.param for s in (prog.sweeps if prog else [])],
    }
    body["run_count"] = body.get("n_runs", 0)
    return body


# --------------------------------------------------------------------------- #
# Plan
# --------------------------------------------------------------------------- #
def _metric_rows(plan) -> List[Dict]:
    from tmesim.protocol.ast import METRICS

    meta = {
        "coverage": {"label": "Coverage @ 1x IC50", "unit": "%", "fmt": "pct"},
        "kill": {"label": "Kill fraction", "unit": "%", "fmt": "pct"},
        "mean_auc": {"label": "Mean tumour AUC", "unit": "uM.s", "fmt": "sci"},
        "penetration_length": {"label": "Penetration length", "unit": "um",
                               "fmt": "num1"},
    }
    keys = plan.runs[0].report if plan.runs else list(METRICS)
    return [{"key": k, "spec": METRICS[k], "better": "higher", **meta.get(k, {})}
            for k in keys if k in METRICS]


def _axes(plan) -> List[Dict]:
    """Sweep axes as the editor's grid dimensions: name, unit, numeric values."""
    out = []
    for axis in plan.sweep_axes:
        if axis["kind"] != "range":
            continue          # a fixed override is one point, not a dimension
        values, unit = [], None
        for text in axis["values"]:
            parts = str(text).split()
            try:
                values.append(float(parts[0]))
            except (ValueError, IndexError):
                values.append(0.0)
            if len(parts) > 1:
                unit = parts[1]
        out.append({"name": axis["param"], "unit": unit, "values": values})
    return out


def _coords(plan, run, axes: List[Dict]) -> List[int]:
    """Index of this run along each axis (what the editor's grid keys on)."""
    coords = []
    for axis in axes:
        text = run.axes.get(axis["name"])
        try:
            value = float(str(text).split()[0])
        except (ValueError, IndexError, AttributeError):
            coords.append(0)
            continue
        best, best_d = 0, float("inf")
        for i, v in enumerate(axis["values"]):
            d = abs(v - value)
            if d < best_d:
                best, best_d = i, d
        coords.append(best)
    return coords


def _assumption_lines(plan) -> List[str]:
    """The dose->PK assumptions, as the sentences the editor renders verbatim."""
    seen, out = set(), []
    for run in plan.runs:
        conv = run.dose_conversion
        if conv is None:
            continue
        for line in conv.conversions:
            if line not in seen:
                seen.add(line)
                out.append(line)
        for a in conv.assumptions:
            line = (f"{a.quantity} = {a.value:g} {a.unit} [{a.provenance}] — "
                    f"{a.source}")
            if line not in seen:
                seen.add(line)
                out.append(line)
        for note in conv.notes + [conv.time_resolved_reason]:
            if note and note not in seen:
                seen.add(note)
                out.append(note)
    return out


def _caveat_lines(plan) -> List[str]:
    out = [d["message"] for d in plan.diagnostics.to_dicts()
           if d["severity"] == "warning"]
    return out


def _refusal(plan) -> Optional[Dict]:
    """The single refusal the editor pins to a line, if there is one."""
    budget = [d for d in plan.diagnostics.errors
              if d.code in ("over-budget", "over-allowance")]
    chosen = budget[0] if budget else (plan.diagnostics.errors[0]
                                       if plan.diagnostics.errors else None)
    if chosen is None:
        return None
    return {
        "code": ("max_runs_exceeded" if chosen.code in
                 ("over-budget", "over-allowance") else chosen.code),
        "line": chosen.line or 1,
        "message": chosen.message,
        "remedy": chosen.hint,
        "remedy_text": chosen.remedy_text,
        "remedy_line": chosen.remedy_line or chosen.line or 1,
        "run_count": plan.n_runs,
        "max_runs": plan.allow_runs if plan.allow_runs is not None
        else plan.max_runs,
    }


def editor_plan(body: Dict, plan) -> Dict:
    """Add the editor's plan projection to a /plan body."""
    body = editor_validate(body, plan)
    axes = _axes(plan)
    cost_s, cost_note = estimate_cost_s(plan)
    per_run = cost_s / plan.n_runs if plan.n_runs else 0.0
    body["run_count"] = plan.n_runs
    body["estimated_cost_s"] = cost_s
    body["cost_note"] = cost_note
    body["axes"] = axes
    body["metrics"] = _metric_rows(plan)
    body["assumptions"] = _assumption_lines(plan)
    body["caveats"] = _caveat_lines(plan)
    body["refusal"] = _refusal(plan)
    # `runs` keeps every contract field and gains the editor's.
    for entry, run in zip(body.get("runs") or [], plan.runs):
        entry["coords"] = _coords(plan, run, axes)
        entry["drug"] = run.drug_alias
        entry["params"] = dict(run.axes)
        entry["est_s"] = round(per_run, 1)
        entry["report"] = list(run.report)
        entry["rank"] = run.rank_by
    return body


# --------------------------------------------------------------------------- #
# Sweep snapshot
# --------------------------------------------------------------------------- #
#: Why a sweep cell has no `run_id`. Stated in the payload rather than left as a
#: silent null.
NO_RUN_ID_REASON = (
    "a TPL sweep is ONE job on the shared queue, not N: its runs share a cached "
    "reconstruction and must be serialised on the single worker anyway. There is "
    "therefore no per-point run_id to open in the viewer — the point's metrics, "
    "provenance and caveats are carried in this payload instead."
)


def sweep_snapshot(job, body: Dict) -> Dict:
    """Project a protocol job into the editor's sweep payload."""
    grid = body.get("grid") or getattr(job, "protocol_grid", None) or []
    rows = {r["index"]: r for r in (body.get("rows") or [])}
    plan_doc = body.get("plan") or {}
    counts = {"queued": 0, "running": 0, "done": 0, "failed": 0, "cancelled": 0}
    for cell in grid:
        counts[cell["status"]] = counts.get(cell["status"], 0) + 1

    status = job.status
    if status == "done" and counts.get("failed"):
        status = "partial" if counts.get("done") else "failed"

    runs = []
    for cell in grid:
        row = rows.get(cell["index"], {})
        runs.append({
            "index": cell["index"],
            "coords": (plan_doc.get("runs") or [{}] * (cell["index"] + 1))[
                cell["index"]].get("coords", [])
            if cell["index"] < len(plan_doc.get("runs") or []) else [],
            "label": cell["label"],
            "drug": cell["drug"],
            "params": cell.get("axes") or {},
            "spec": row.get("spec"),
            "run_id": None,
            "run_id_reason": NO_RUN_ID_REASON,
            "status": cell["status"],
            "progress": 1.0 if cell["status"] == "done" else 0.0,
            "stage": cell["status"],
            "metrics": _editor_metrics(cell.get("metrics")),
            "result": row or None,
            "error": cell.get("error"),
            "runtime_s": cell.get("runtime_s"),
            "est_s": None,
            "report": row.get("reported"),
        })

    return {
        "sweep_id": job.run_id,
        "status": status,
        "run_count": len(grid) or body.get("n_runs") or 0,
        "counts": counts,
        "axes": plan_doc.get("axes") or [],
        "metrics": plan_doc.get("metrics") or [],
        "caveats": plan_doc.get("caveats") or [],
        "assumptions": plan_doc.get("assumptions") or [],
        "expects": body.get("expects") or [],
        "estimated_cost_s": plan_doc.get("estimated_cost_s"),
        "elapsed_s": round(job.elapsed_s, 3),
        "runs": runs,
        "error": job.error,
    }


def _editor_metrics(metrics: Optional[Dict]) -> Optional[Dict]:
    """Engine metric keys -> the TPL metric names the editor's grid selects on."""
    if not metrics:
        return None
    from tmesim.protocol.ast import METRICS

    out = {tpl: metrics.get(engine) for tpl, engine in METRICS.items()}
    out.update(metrics)          # the engine keys stay available too
    return out
