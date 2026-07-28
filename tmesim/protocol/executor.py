"""Plan and execute a TPL program.

Two halves, deliberately separable:

:func:`plan`
    Everything that can be decided WITHOUT compute — resolve the drugs, expand
    the sweep cross-product, convert dose -> PK, check the specimen's
    capability, count the runs and refuse an over-budget sweep. This is what
    ``POST /api/protocol/plan`` returns.
:func:`execute_plan`
    Runs the plan through the EXISTING engine entry points, in the same order
    and with the same production constants ``app/server/jobs.py`` uses, so a
    protocol run is numerically the same kind of object as an ``/api/runs`` run.
    Reconstructions come from ``app.server.cache.RECON_CACHE``, so a sweep over
    dose on one specimen reconstructs ONCE.

Nothing here reimplements science. The chain is::

    RECON_CACHE.get_or_build            (reconstruct.real_pipeline)
    apply_drug_sensitivity              (reconstruct.real_pipeline)
    molecule -> molecule_to_drugparams -> calibrated_drugparams(alpha)
             -> apply_barrier_to_drugparams        (molecular / validation)
    simulate_response                   (sim.engine)
    delivery_metrics.compute_all        (validation)

The one thing this module DOES set that the ``/api/runs`` path does not is
``PDEConfig.time_resolved = True``, and it says so in every result — without it
the solver collapses the plasma curve to its window mean and `every` / `repeat`
/ `cycles` / `over` become structurally inert (``results/dose_schedule.md``).
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from . import ast as A
from .dosing import DEFAULT_MAX_WINDOW_S, DoseConversion, convert_dose
from .errors import BudgetError, Diagnostic, ResolutionError, TplError, collect
from .parser import parse
from .resolve import ResolvedDrug, check_capability, resolve_drug
from .units import Quantity, _fmt

__all__ = [
    "DEFAULT_MAX_RUNS", "RunPlan", "Plan", "plan", "execute_plan",
    "validate_source",
]

#: The sweep budget. A two-axis sweep is trivially a hundred PDE solves, so the
#: planner refuses above this unless the script says `allow runs N`.
DEFAULT_MAX_RUNS = 64

#: Production solver constants, taken from app/server/jobs.py so a protocol run
#: and an /api/runs run are the same computation.
PDE_DT_S = 10.0
PDE_SAVE_EVERY = 500
PDE_CONVECTION = True
SIM_SEED = 0

#: Fixed reconstruction modes. TPL exposes `state_mode` (the spec's run option)
#: but not `oxygen_mode` / `vascular_mode`, so both stay at the production
#: default rather than being silently variable.
VASCULAR_MODE = "discrete"
OXYGEN_MODE = "pde"


# --------------------------------------------------------------------------- #
# Plan containers
# --------------------------------------------------------------------------- #
@dataclass
class RunPlan:
    """One point of the cross-product: one drug, one parameter combination."""

    index: int = 0
    drug_alias: str = ""
    protocol_name: str = ""
    section: str = ""
    axes: Dict[str, str] = field(default_factory=dict)
    dx_um: float = 0.0
    state_mode: str = "oxygen"
    sensitivity_mode: str = "uniform"
    report: List[str] = field(default_factory=list)
    rank_by: Optional[str] = None
    dose_conversion: Optional[DoseConversion] = None
    recon_key: Tuple = ()
    line: int = 0
    col: int = 0

    @property
    def label(self) -> str:
        axes = " ".join(f"{k}={v}" for k, v in self.axes.items())
        return f"{self.drug_alias}/{self.protocol_name}@{self.section}" + (
            f" [{axes}]" if axes else "")

    def spec(self) -> Dict:
        """A ``RunSpec``-shaped dict (API_CONTRACT.md), for the UI and the log."""
        dc = self.dose_conversion
        return {
            "section": self.section,
            "dx_um": self.dx_um,
            "drug": self.drug_alias,
            "sensitivity_mode": self.sensitivity_mode,
            "vascular_mode": VASCULAR_MODE,
            "oxygen_mode": OXYGEN_MODE,
            "state_mode": self.state_mode,
            "mgmt": None,
            "pk": (dc.to_dict()["pk"] if dc else None),
            "total_time_s": (dc.window_s if dc else None),
        }

    def to_dict(self) -> Dict:
        return {
            "index": self.index,
            "label": self.label,
            "drug": self.drug_alias,
            "protocol": self.protocol_name,
            "section": self.section,
            "axes": dict(self.axes),
            "dx_um": self.dx_um,
            "state_mode": self.state_mode,
            "sensitivity_mode": self.sensitivity_mode,
            "report": list(self.report),
            "rank_by": self.rank_by,
            "spec": self.spec(),
            "dose_conversion": (self.dose_conversion.to_dict()
                                if self.dose_conversion else None),
            "reconstruction_key": list(self.recon_key),
            "line": self.line,
            "col": self.col,
        }


@dataclass
class Plan:
    """The whole expanded, checked, still-unexecuted request."""

    program: Optional[A.Program] = None
    diagnostics: collect = field(default_factory=collect)
    runs: List[RunPlan] = field(default_factory=list)
    drugs: Dict[str, ResolvedDrug] = field(default_factory=dict)
    capability: Dict[str, Dict] = field(default_factory=dict)
    sweep_axes: List[Dict] = field(default_factory=list)
    expects: List[Dict] = field(default_factory=list)
    max_runs: int = DEFAULT_MAX_RUNS
    allow_runs: Optional[int] = None
    estimated_cost: Dict = field(default_factory=dict)
    source: str = ""

    @property
    def ok(self) -> bool:
        return self.diagnostics.ok

    @property
    def n_runs(self) -> int:
        return len(self.runs)

    def to_dict(self) -> Dict:
        return {
            "ok": self.ok,
            "diagnostics": self.diagnostics.to_dicts(),
            "n_runs": self.n_runs,
            "max_runs": self.max_runs,
            "allow_runs": self.allow_runs,
            "runs": [r.to_dict() for r in self.runs],
            "drugs": {k: v.to_dict() for k, v in self.drugs.items()},
            "capability": dict(self.capability),
            "sweep_axes": list(self.sweep_axes),
            "expects": list(self.expects),
            "estimated_cost": dict(self.estimated_cost),
        }


# --------------------------------------------------------------------------- #
# validate (parse + cheap semantics only)
# --------------------------------------------------------------------------- #
def validate_source(source: str) -> Dict:
    """``POST /api/protocol/validate`` — parse and check without compiling.

    Returns ``{ok, diagnostics[]}``. Every diagnostic carries a line and a
    column. Nothing here runs the molecular pipeline or touches the solver, so
    it is fast enough for keystroke-level use.
    """
    p = plan(source, compile_drugs=False, check_budget=True)
    return {"ok": p.ok, "diagnostics": p.diagnostics.to_dicts(),
            "n_runs": p.n_runs, "max_runs": p.max_runs}


# --------------------------------------------------------------------------- #
# plan
# --------------------------------------------------------------------------- #
def plan(source: str, compile_drugs: bool = True, check_budget: bool = True,
         max_runs: int = DEFAULT_MAX_RUNS,
         max_window_s: Optional[float] = None) -> Plan:
    """Parse, resolve, expand and check a TPL program without executing it."""
    if max_window_s is None:
        from . import dosing as _dosing

        max_window_s = _dosing.DEFAULT_MAX_WINDOW_S
    program, diags = parse(source)
    out = Plan(program=program, diagnostics=diags, max_runs=int(max_runs),
               source=source or "")
    if program is None:
        return out
    # A statement that failed to PARSE never made it into the program, so a
    # later `run X with Y` would then also report "unknown protocol Y". That
    # cascade buries the real error, so name-resolution failures are suppressed
    # (and said to be suppressed) when the file already has a syntax error.
    syntax_broken = bool(diags.errors)
    if program.allow_runs is not None:
        out.allow_runs = program.allow_runs.count

    # ---- duplicate definitions ------------------------------------------- #
    drug_defs: Dict[str, A.DrugDef] = {}
    for d in program.drugs:
        if d.name in drug_defs:
            diags.error(f"drug {d.name!r} is defined twice (first on line "
                        f"{drug_defs[d.name].line}).", node=d,
                        code="duplicate-drug")
            continue
        drug_defs[d.name] = d
    proto_defs: Dict[str, A.ProtocolDef] = {}
    for pr in program.protocols:
        if pr.name in proto_defs:
            diags.error(f"protocol {pr.name!r} is defined twice (first on line "
                        f"{proto_defs[pr.name].line}).", node=pr,
                        code="duplicate-protocol")
            continue
        proto_defs[pr.name] = pr

    if not program.runs:
        diags.error(
            "this program defines no run: add a `run <drug> with <protocol> on "
            "<specimen>` (or a `compare`) block. Definitions alone execute "
            "nothing.", line=1, col=1, code="no-run")

    # ---- protocols are complete ------------------------------------------ #
    for pr in proto_defs.values():
        if pr.dose is None:
            diags.error(f"protocol {pr.name}: no `dose` statement. A protocol "
                        "without a dose delivers nothing.", node=pr,
                        code="incomplete-protocol")
        if pr.every is None and (pr.repeat and pr.repeat.count > 1):
            diags.error(
                f"protocol {pr.name}: `repeat {pr.repeat.count}` needs an "
                "`every <T>` interval, otherwise every dose lands at t = 0.",
                node=pr.repeat, code="incomplete-protocol")
        if pr.rest is not None and (pr.cycles is None or pr.cycles.count <= 1):
            diags.warn(
                f"protocol {pr.name}: `then rest` only has an effect with "
                "`cycles > 1` — with a single cycle nothing follows the rest.",
                node=pr.rest, code="inert-statement")

    # ---- resolve the drugs ------------------------------------------------ #
    alpha: Optional[float] = None
    for name, node in drug_defs.items():
        try:
            if compile_drugs and alpha is None:
                from ..validation.calibration import calibrate_alpha

                alpha = float(calibrate_alpha())
            out.drugs[name] = resolve_drug(
                node, compile_params=compile_drugs, alpha=alpha)
        except TplError as exc:
            for d in exc.diagnostics:
                if not d.line:
                    d.line, d.col = node.line, node.col
                diags.add(d)
        except Exception as exc:  # noqa: BLE001 - never a bare traceback
            diags.error(
                f"drug {name}: compilation failed — {type(exc).__name__}: {exc}",
                node=node, code="compile-failed")

    # ---- sweeps ----------------------------------------------------------- #
    axes: List[Tuple[str, List[Quantity], A.Sweep]] = []
    seen_params: Dict[str, int] = {}
    for sw in program.sweeps:
        if sw.param in seen_params:
            diags.error(
                f"sweep {sw.param!r} is declared twice (first on line "
                f"{seen_params[sw.param]}); a parameter has one axis.",
                node=sw, code="duplicate-sweep")
            continue
        seen_params[sw.param] = sw.line
        try:
            values = _expand_sweep(sw)
        except TplError as exc:
            for d in exc.diagnostics:
                diags.add(d)
            continue
        axes.append((sw.param, values, sw))
        out.sweep_axes.append({
            "param": sw.param, "kind": sw.kind, "n_points": len(values),
            "values": [str(v) for v in values],
            "line": sw.line, "col": sw.col,
        })

    # ---- expand the cross-product ----------------------------------------- #
    combos = _cross_product(axes)
    index = 0
    for stmt in program.runs:
        drugs = list(getattr(stmt, "drugs", []) or [])
        proto_name = stmt.protocol
        section = stmt.section
        opts = stmt.options

        missing_drug = [d for d in drugs if d not in drug_defs]
        for d in missing_drug:
            if syntax_broken:
                continue
            diags.error(
                f"unknown drug {d!r}: no `drug {d} from ...` block in this "
                f"program. Defined: {', '.join(drug_defs) or '(none)'}.",
                node=stmt, code="unknown-drug",
                col=getattr(stmt, "drug_col", 0) or stmt.col)
        if proto_name not in proto_defs:
            if not syntax_broken:
                diags.error(
                    f"unknown protocol {proto_name!r}: no `protocol "
                    f"{proto_name}` block. Defined: "
                    f"{', '.join(proto_defs) or '(none)'}.",
                    node=stmt, code="unknown-protocol",
                    col=getattr(stmt, "protocol_col", 0) or stmt.col)
            continue
        if missing_drug:
            continue
        proto = proto_defs[proto_name]
        if proto.dose is None:
            continue

        state_mode = opts.state_mode or "oxygen"
        sens_mode = opts.sensitivity or "uniform"
        cap, cap_diags = check_capability(
            section, state_mode, sens_mode, OXYGEN_MODE,
            line=(opts.state_mode_line or stmt.line),
            col=(opts.state_mode_col or getattr(stmt, "section_col", 0)
                 or stmt.col))
        out.capability[f"{section}:{state_mode}:{sens_mode}"] = cap
        diags.extend(cap_diags)
        if not cap_diags.ok:
            continue

        for combo in combos:
            dx_q = combo.get("dx") or opts.dx
            if dx_q is None:
                diags.error(
                    f"run {proto_name} on {section}: no grid spacing. Add "
                    "`dx <N> um` to the run block or `sweep dx <N> um`.",
                    node=stmt, code="missing-dx")
                break
            dx_um, _ = dx_q.canonical()
            if dx_um <= 0:
                diags.error("dx must be > 0.", node=stmt, code="out-of-range")
                break

            for alias in drugs:
                resolved = out.drugs.get(alias)
                try:
                    conv = _dose_to_pk(proto, combo, resolved, compile_drugs,
                                       max_window_s)
                except TplError as exc:
                    for d in exc.diagnostics:
                        diags.add(d)
                    conv = None
                except Exception as exc:  # noqa: BLE001
                    diags.error(
                        f"protocol {proto_name}: dose -> PK conversion failed — "
                        f"{exc}", node=proto.dose or proto,
                        code="dose-conversion-failed")
                    conv = None

                rp = RunPlan(
                    index=index, drug_alias=alias, protocol_name=proto_name,
                    section=section, dx_um=float(dx_um),
                    state_mode=state_mode, sensitivity_mode=sens_mode,
                    report=list(opts.report) or list(A.METRICS),
                    rank_by=opts.rank_by,
                    dose_conversion=conv,
                    axes={k: str(v) for k, v in combo.items()},
                    recon_key=(section, round(float(dx_um), 6), VASCULAR_MODE,
                               OXYGEN_MODE, state_mode),
                    line=stmt.line, col=stmt.col,
                )
                out.runs.append(rp)
                index += 1

    # ---- expectations ------------------------------------------------------ #
    for e in program.expects:
        out.expects.append({
            "metric": e.metric, "comparator": e.comparator,
            "value": e.value.value if e.value else None,
            "unit": e.value.unit if e.value else None,
            "line": e.line, "col": e.col,
            "note": ("reported as pass/fail; a failed expectation MARKS the "
                     "run, it never alters or blocks the result"),
        })

    # ---- cost + budget ----------------------------------------------------- #
    out.estimated_cost = _estimate_cost(out)
    if check_budget:
        _check_budget(out)
    return out


def _expand_sweep(sw: A.Sweep) -> List[Quantity]:
    """Materialise one sweep axis into its list of values."""
    if sw.kind == "fixed":
        return [sw.fixed]
    start, stop = float(sw.start.value), float(sw.stop.value)
    unit = sw.stop.unit
    if sw.n_steps is not None:
        n = int(sw.n_steps)
        if n == 1:
            vals = [start]
        else:
            step = (stop - start) / (n - 1)
            vals = [start + i * step for i in range(n)]
    else:
        step = float(sw.step)
        if (stop - start) * step < 0:
            raise BudgetError(
                f"sweep {sw.param}: `step {_fmt(step)}` walks away from "
                f"{_fmt(stop)} (from {_fmt(start)}). Reverse the range or the "
                "sign of the step.",
                line=sw.line, col=sw.col, end_col=sw.end_col,
                code="empty-sweep")
        n = int(math.floor(abs(stop - start) / abs(step) + 1e-9)) + 1
        if n > 10_000:
            raise BudgetError(
                f"sweep {sw.param}: `from {_fmt(start)} to {_fmt(stop)} step "
                f"{_fmt(step)}` is {n} points. That is not a sweep, it is a "
                "typo — check the step size.",
                line=sw.line, col=sw.col, end_col=sw.end_col,
                code="absurd-sweep")
        vals = [start + i * step for i in range(n)]
    return [Quantity(v, unit, sw.line, sw.col, sw.end_col) for v in vals]


def _cross_product(axes) -> List[Dict[str, Quantity]]:
    combos: List[Dict[str, Quantity]] = [{}]
    for param, values, _sw in axes:
        combos = [dict(c, **{param: v}) for c in combos for v in values]
    return combos


def _dose_to_pk(proto: A.ProtocolDef, combo: Dict[str, Quantity],
                resolved: Optional[ResolvedDrug], compiled: bool,
                max_window_s: float) -> Optional[DoseConversion]:
    """Apply the sweep coordinates to a protocol and convert it to PKConfig."""
    dose_q = combo.get("dose") or proto.dose.amount
    over_q = combo.get("over") or proto.dose.over
    every_q = combo.get("every") or (proto.every.interval if proto.every else None)
    repeat = int(combo["repeat"].value) if "repeat" in combo else (
        proto.repeat.count if proto.repeat else 1)
    cycles = int(combo["cycles"].value) if "cycles" in combo else (
        proto.cycles.count if proto.cycles else 1)
    rest_s = proto.rest.duration.canonical()[0] if proto.rest else 0.0
    every_s = every_q.canonical()[0] if every_q is not None else 86400.0

    route = proto.route.value if proto.route else "iv"

    mw = None
    key = None
    if resolved is not None:
        key = resolved.identity
        if resolved.drug_params is not None:
            mw = float(resolved.drug_params.molecular_weight)
        elif resolved.molecule is not None:
            mw = float(resolved.molecule.get("mw") or 0.0)
    if not mw:
        if not compiled:
            # validate-only: the molecular weight is not needed to check syntax,
            # units, budget or capability, and compiling is what /plan is for.
            return None
        raise ResolutionError(
            f"protocol {proto.name}: the drug has no molecular weight, so a "
            "milligram dose cannot be converted to a molar plasma "
            "concentration. Supply `molecular_weight <N>` in the drug block.",
            line=proto.dose.line, col=proto.dose.col, code="missing-mw")

    return convert_dose(
        dose=dose_q, route=route, every_s=every_s, repeat=repeat,
        rest_s=rest_s, cycles=cycles, molecular_weight=mw,
        drug_key=key, over=over_q, max_window_s=max_window_s,
    )


def _estimate_cost(p: Plan) -> Dict:
    """What this plan will cost, in units that are actually knowable."""
    from app.server.cache import RECON_CACHE, SECTIONS, describe_section

    keys = []
    for r in p.runs:
        if r.recon_key not in keys:
            keys.append(r.recon_key)
    cached = [k for k in keys if RECON_CACHE.peek(k) is not None]

    total_steps = 0
    grids: Dict[str, List[int]] = {}
    for r in p.runs:
        if r.dose_conversion:
            total_steps += int(math.ceil(r.dose_conversion.window_s / PDE_DT_S))
        row = describe_section(r.section) if r.section in SECTIONS else {}
        extent = row.get("extent_um")
        if extent and r.dx_um > 0:
            grids.setdefault(
                f"{r.section}@dx={_fmt(r.dx_um)}um",
                [int(math.ceil(extent[0] / r.dx_um)),
                 int(math.ceil(extent[1] / r.dx_um))])

    return {
        "n_runs": p.n_runs,
        "n_reconstructions": len(keys),
        "reconstructions_cached": len(cached),
        "reconstructions_to_build": len(keys) - len(cached),
        "reconstruction_keys": [list(k) for k in keys],
        "pde_dt_s": PDE_DT_S,
        "total_pde_timesteps": total_steps,
        "approx_grid_cells": {k: v[0] * v[1] for k, v in grids.items()},
        "wall_time_s": None,
        "wall_time_reason": (
            "not predicted. Wall time is dominated by the reconstruction "
            "(1 s to ~90 s depending on section and dx) and by the sparse "
            "solve, neither of which this backend benchmarks in advance. The "
            "timestep count and grid size above are the honest cost figures."),
        "sweep_saves_reconstructions": (
            f"{p.n_runs} runs share {len(keys)} reconstruction(s): a sweep over "
            "dose on one specimen at one dx and one state_mode reconstructs "
            "ONCE and reuses it."),
    }


def _check_budget(p: Plan) -> None:
    """Refuse a sweep above ``max_runs`` without an explicit ``allow runs N``."""
    n = p.n_runs
    if n == 0:
        return
    limit = p.allow_runs if p.allow_runs is not None else p.max_runs
    if n <= limit:
        return
    axes = ", ".join(f"{a['param']}({a['n_points']})" for a in p.sweep_axes)
    node_line = p.program.allow_runs.line if p.program.allow_runs else (
        p.program.sweeps[0].line if p.program.sweeps else 1)
    node_col = p.program.allow_runs.col if p.program.allow_runs else (
        p.program.sweeps[0].col if p.program.sweeps else 1)
    if p.allow_runs is not None:
        msg = (f"this program expands to {n} runs, above the {limit} you "
               f"allowed with `allow runs {limit}`. Axes: {axes or '(none)'}. "
               f"Raise the allowance to at least {n} or shrink the sweep.")
        code = "over-allowance"
    else:
        msg = (f"this program expands to {n} runs, above the default budget of "
               f"{p.max_runs}. Axes: {axes or '(none)'}. Each run is a full PDE "
               "solve, so a two-axis sweep is trivially a hundred of them. Add "
               f"`allow runs {n}` to say you meant it, or shrink the sweep.")
        code = "over-budget"
    p.diagnostics.error(msg, line=node_line, col=node_col, code=code,
                        hint=f"allow runs {n}")


# --------------------------------------------------------------------------- #
# execute
# --------------------------------------------------------------------------- #
def execute_plan(p: Plan, log: Optional[Callable[[str], None]] = None,
                 on_progress: Optional[Callable[[List[Dict]], None]] = None,
                 check_cancelled: Optional[Callable[[], None]] = None,
                 cache=None) -> Dict:
    """Execute an already-checked plan, sequentially, on the shared cache.

    Sequential on purpose: the cached ``Microenvironment`` is shared, and
    ``apply_drug_sensitivity`` writes ``cells.sensitivity`` on it. One run at a
    time makes that safe by construction — the same reason
    ``app/server/jobs.py`` uses exactly one worker.
    """
    def _log(msg: str) -> None:
        if log is not None:
            log(msg)

    if not p.ok:
        raise TplError("refusing to execute a plan with errors",
                       diagnostics=list(p.diagnostics.items))

    import numpy as np

    from app.server.cache import RECON_CACHE
    from tmesim.core import PDEConfig, PKConfig, SimConfig
    from tmesim.reconstruct.real_pipeline import (
        apply_drug_sensitivity,
        reconstruction_report,
    )
    from tmesim.sim.engine import simulate_response
    from tmesim.validation import delivery_metrics

    cache = cache if cache is not None else RECON_CACHE

    grid: List[Dict] = [
        {"index": r.index, "label": r.label, "axes": dict(r.axes),
         "drug": r.drug_alias, "section": r.section, "status": "queued",
         "metrics": None, "error": None, "runtime_s": None,
         "reconstruction_cache_hit": None}
        for r in p.runs
    ]
    if on_progress:
        on_progress(grid)

    rows: List[Dict] = []
    recon_built: List[Tuple] = []
    recon_hits = 0
    t_all = time.time()

    for cell, rp in zip(grid, p.runs):
        if check_cancelled:
            check_cancelled()
        cell["status"] = "running"
        if on_progress:
            on_progress(grid)
        t0 = time.time()
        try:
            resolved = p.drugs[rp.drug_alias]
            dp = resolved.drug_params
            conv = rp.dose_conversion
            caveats: List[str] = []

            # --- reconstruction (shared, cached) ------------------------- #
            was_cached = cache.peek(rp.recon_key) is not None
            entry = cache.get_or_build(rp.recon_key, log=log)
            if was_cached:
                recon_hits += 1
            elif rp.recon_key not in recon_built:
                recon_built.append(rp.recon_key)
            cell["reconstruction_cache_hit"] = bool(was_cached)
            micro = entry["micro"]
            caveats.extend(entry.get("notes") or [])

            micro.cells.sensitivity = np.array(entry["pristine_sensitivity"],
                                               copy=True)
            if entry.get("pristine_state") is not None:
                micro.cells.state = np.array(entry["pristine_state"], copy=True)
            micro.metadata = dict(entry["pristine_metadata"])

            # --- per-cell sensitivity ------------------------------------ #
            sens_name = resolved.identity or rp.drug_alias
            sens_info = apply_drug_sensitivity(
                micro, sens_name, mode=rp.sensitivity_mode)
            report = reconstruction_report(micro)
            report["reconstruction_cache_hit"] = bool(was_cached)
            report["reconstruction_build_time_s"] = entry["build_time_s"]
            report["vascular_mode"] = VASCULAR_MODE
            report["vascular_source"] = entry["report"].get("vascular_source")

            state_block = report.get("cell_state") or {}
            if state_block.get("caveat"):
                caveats.append(f"cell_state ({state_block.get('mode')}): "
                               f"{state_block['caveat']}")
            if sens_info.get("upper_bound"):
                caveats.append(
                    "ABSOLUTE kill_fraction is a NO-INTRINSIC-RESISTANCE UPPER "
                    f"BOUND in sensitivity_mode={rp.sensitivity_mode!r}. "
                    + str(sens_info.get("comparability") or ""))

            # --- simulate ------------------------------------------------- #
            pk = PKConfig(plasma_peak=conv.plasma_peak_uM,
                          half_life=conv.half_life_s,
                          n_doses=conv.n_doses,
                          dose_interval=conv.dose_interval_s)
            pde = PDEConfig(dt=PDE_DT_S, total_time=conv.window_s,
                            convection=PDE_CONVECTION,
                            save_every=PDE_SAVE_EVERY, time_resolved=True)
            if not getattr(pde, "time_resolved", False):  # pragma: no cover
                raise RuntimeError(
                    "this engine build has no PDEConfig.time_resolved; a TPL "
                    "schedule cannot be resolved without it.")
            cfg = SimConfig(drug=dp, pk=pk, pde=pde, seed=SIM_SEED)
            _log(f"[{rp.index}] {rp.label}: simulate_response("
                 f"dt={PDE_DT_S}s, total_time={_fmt(conv.window_s)}s, "
                 f"time_resolved=True, plasma_peak={_fmt(conv.plasma_peak_uM)} uM, "
                 f"half_life={_fmt(conv.half_life_s)} s, n_doses={conv.n_doses})")
            res = simulate_response(micro, cfg)

            dm = delivery_metrics.compute_all(micro, res, dp)
            head = dm["headline"]
            metrics = {
                "kill_fraction": _finite(res.kill_fraction),
                "therapeutic_coverage_1x": _finite(
                    head["therapeutic_coverage_1x_ic50"]),
                "mean_tumor_auc": _finite(head["mean_tumor_auc_uM_s"]),
                "penetration_length_um": _finite(
                    head["penetration_length_L_um"]),
            }
            caveats.append(conv.time_resolved_reason)
            caveats.extend(conv.notes)
            caveats.append(
                "penetration_index and kill_fraction are resolution-sensitive; "
                f"this run used dx = {_fmt(rp.dx_um)} um"
                + ("" if rp.dx_um <= 20.0 else
                   ", which is COARSER than the 20 um at which kill_fraction is "
                   "grid-converged (results/numerical_rigor.md §2), so both are "
                   "upper bounds."))
            caveats.append(
                "This tissue is an untreated resection with no outcome label. "
                "Every number here is a PREDICTION of the physical model and "
                "has never been scored against a measured treatment response. "
                "NOT FOR CLINICAL USE.")

            row = {
                "index": rp.index,
                "label": rp.label,
                "axes": dict(rp.axes),
                "drug": rp.drug_alias,
                "protocol": rp.protocol_name,
                "section": rp.section,
                "spec": rp.spec(),
                "metrics": metrics,
                "reported": {m: metrics.get(A.METRICS[m]) for m in rp.report},
                "extra_metrics": {
                    "median_tumor_auc": head["median_tumor_auc_uM_s"],
                    "therapeutic_coverage_0_5x":
                        head["therapeutic_coverage_0.5x_ic50"],
                    "therapeutic_coverage_2x":
                        head["therapeutic_coverage_2x_ic50"],
                    "viable_tumor_fraction": float(res.viable_tumor_fraction),
                    "plasma_peak_uM": float(res.extras["plasma_peak"]),
                    "plasma_mean_uM": float(res.extras["plasma_mean"]),
                    "time_resolved": float(res.extras.get("time_resolved", 0.0)),
                    "n_tumor": int(res.extras["n_tumor"]),
                },
                "provenance": {
                    "drug": resolved.to_dict(),
                    "dose_to_pk": conv.to_dict(),
                    "sensitivity": {
                        "mode": rp.sensitivity_mode,
                        "drug_used_for_lookup": sens_name,
                        "upper_bound": bool(sens_info.get("upper_bound", True)),
                        "comparability": sens_info.get("comparability"),
                        "determinant_gene": sens_info.get("determinant_gene"),
                    },
                    "capability": p.capability.get(
                        f"{rp.section}:{rp.state_mode}:{rp.sensitivity_mode}"),
                    "cell_state": state_block,
                    "solver": {
                        "pde_dt_s": PDE_DT_S,
                        "pde_save_every": PDE_SAVE_EVERY,
                        "pde_convection": PDE_CONVECTION,
                        "pde_time_resolved": True,
                        "seed": SIM_SEED,
                        "vascular_mode": VASCULAR_MODE,
                        "oxygen_mode": OXYGEN_MODE,
                        "state_mode": rp.state_mode,
                    },
                    "engine_entrypoints": [
                        "app.server.cache.ReconstructionCache.get_or_build",
                        "reconstruct.real_pipeline.apply_drug_sensitivity",
                        "molecular.compile_drug.molecule_to_drugparams",
                        "validation.calibration.calibrated_drugparams",
                        "molecular.binding_site_barrier.apply_barrier_to_drugparams",
                        "sim.engine.simulate_response",
                        "validation.delivery_metrics.compute_all",
                    ],
                    "caveats": caveats,
                },
                "report": {k: report.get(k) for k in
                           ("n_cells", "grid", "extent_um", "mean_ecm_density",
                            "mean_vascular_density", "reconstruction_cache_hit",
                            "vascular_source")},
                "runtime_s": round(time.time() - t0, 3),
            }
            row["expectations"] = _score_expectations(p.expects, metrics)
            rows.append(row)
            cell.update(status="done", metrics=metrics,
                        runtime_s=row["runtime_s"])
            _log(f"[{rp.index}] done in {row['runtime_s']:.1f}s: "
                 f"kill={metrics['kill_fraction']}")
        except Exception as exc:  # noqa: BLE001 - surfaced per run, never swallowed
            if type(exc).__name__ == "RunCancelled":
                cell.update(status="cancelled", error=str(exc))
                if on_progress:
                    on_progress(grid)
                raise
            cell.update(status="failed", error=f"{type(exc).__name__}: {exc}",
                        runtime_s=round(time.time() - t0, 3))
            rows.append({"index": rp.index, "label": rp.label,
                         "axes": dict(rp.axes), "drug": rp.drug_alias,
                         "section": rp.section, "metrics": None,
                         "error": f"{type(exc).__name__}: {exc}"})
            _log(f"[{rp.index}] FAILED: {type(exc).__name__}: {exc}")
        if on_progress:
            on_progress(grid)

    ranking = _rank(p, rows)
    return {
        "n_runs": len(p.runs),
        "grid": grid,
        "rows": rows,
        "ranking": ranking,
        "sweep_axes": list(p.sweep_axes),
        "expects": list(p.expects),
        "diagnostics": p.diagnostics.to_dicts(),
        "reconstructions": {
            "built": [list(k) for k in recon_built],
            "n_built": len(recon_built),
            "n_cache_hits": recon_hits,
            "note": (f"{len(p.runs)} run(s) needed {len(recon_built)} "
                     "reconstruction(s); the rest were served from "
                     "app.server.cache.RECON_CACHE."),
        },
        "runtime_s": round(time.time() - t_all, 3),
    }


def _finite(value) -> Optional[float]:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _score_expectations(expects: List[Dict], metrics: Dict) -> List[Dict]:
    """Score every ``expect`` against this run. Reported, never enforced."""
    out = []
    for e in expects:
        key = A.METRICS.get(e["metric"])
        actual = metrics.get(key)
        target = e["value"]
        if e.get("unit") == "%":
            target = float(target) / 100.0
        if actual is None or target is None:
            out.append({**e, "actual": actual, "passed": None,
                        "reason": "metric unavailable for this run"})
            continue
        cmp = e["comparator"]
        passed = {
            ">": actual > target, "<": actual < target,
            ">=": actual >= target, "<=": actual <= target,
            "==": actual == target, "!=": actual != target,
        }[cmp]
        out.append({**e, "actual": actual, "passed": bool(passed)})
    return out


def _rank(p: Plan, rows: List[Dict]) -> Optional[Dict]:
    """``rank by <metric>`` for a ``compare`` statement."""
    metric = next((r.rank_by for r in p.runs if r.rank_by), None)
    if not metric:
        return None
    key = A.METRICS[metric]
    scored = [(r["drug"], r.get("metrics", {}).get(key) if r.get("metrics")
               else None, r["index"]) for r in rows]
    usable = [s for s in scored if s[1] is not None]
    usable.sort(key=lambda s: s[1], reverse=True)
    return {
        "metric": metric,
        "engine_key": key,
        "order": [{"drug": d, "value": v, "run_index": i}
                  for d, v, i in usable],
        "unranked": [{"drug": d, "run_index": i}
                     for d, v, i in scored if v is None],
        "note": ("ranking is the defensible output of a section whose panel "
                 "cannot support an absolute kill; it is invariant across the "
                 "three state routes (results/cell_state_impact.md §3.3)."),
    }
