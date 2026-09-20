#!/usr/bin/env python
"""Where the wall-clock and the peak RSS of one Studio run actually go.

This is a MEASUREMENT harness, not a physics experiment. It drives exactly the
code path ``app/server/jobs.py`` drives — ``execute_run(spec, ctx, cache)`` —
with a ``_Ctx`` stand-in that timestamps every ``ctx.stage()`` transition, so
the stage split is the SERVER'S OWN stage split and not a re-implementation of
it. Nothing here re-derives a number the pipeline already computes.

Three arms, because they answer three different questions and cannot be run in
one process without contaminating each other:

``stages``
    One run, clean (no profiler, no instrumentation). Authoritative wall-clock
    per server stage, from ``ctx.stage()`` boundaries. This is the number to
    quote.

``profile``
    The same run under ``cProfile``. Attributes time to named functions —
    which is how ``simulate`` gets split into oxygen / IFP / transport / kill /
    mechanism, none of which the server's own stage list distinguishes. The
    profiled total is ALWAYS reported beside the clean total so the tracing
    overhead is visible rather than folded into a claim.

``instrument``
    The same run with two timing shims installed: one around
    ``tmesim.physics.diffusion.factorized`` (splitting the sparse LU
    factorisation from the per-step triangular solves, which cProfile cannot
    see because ``splu(A).solve`` is a C method) and one around
    ``_advective_divergence``. The shims add timing calls and change NO
    arithmetic; the run's ``kill_fraction`` is printed so it can be compared
    against the clean arm's.

``memory``
    One run in THIS process, reporting ``resource.getrusage`` peak RSS. Meant
    to be launched one-per-process by ``sweep`` so each point is a cold
    interpreter with an empty reconstruction cache — a peak-RSS number from a
    process that has already held another grid is not a peak-RSS number for the
    grid you asked about.

``sweep``
    Spawns ``memory`` in a fresh subprocess per (section, dx) point and prints
    the curve.

Usage
-----
    .venv/bin/python experiments/profile_pipeline.py stages --dx 12
    .venv/bin/python experiments/profile_pipeline.py profile --dx 30 --top 30
    .venv/bin/python experiments/profile_pipeline.py instrument --dx 12
    .venv/bin/python experiments/profile_pipeline.py memory --dx 30 --json
    .venv/bin/python experiments/profile_pipeline.py sweep \
        --dx 60,48,36,30,24,20,16,12 --section xenium_319

Every arm prints the full RunSpec it executed. A runtime without a spec is not
a measurement on this product (docs/DEPLOYMENT.md), and neither is a peak RSS.
"""
from __future__ import annotations

import argparse
import cProfile
import json
import os
import pstats
import resource
import subprocess
import sys
import time
from typing import Dict, List, Optional

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)


# --------------------------------------------------------------------------- #
# RSS helpers
# --------------------------------------------------------------------------- #
def _maxrss_bytes() -> int:
    """Peak resident set size of this process, in BYTES.

    ``ru_maxrss`` is bytes on Darwin and kilobytes on Linux; the unit is not
    portable and getting it wrong is a 1024x error in the one number an
    operator sizes an instance from. Calibrated on this machine at the top of
    results/performance_profile_2026-09-20.md.
    """
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(raw) if sys.platform == "darwin" else int(raw) * 1024


def _mib(n: float) -> float:
    return float(n) / (1024.0 * 1024.0)


# --------------------------------------------------------------------------- #
# The spec, and the stage-timing ctx
# --------------------------------------------------------------------------- #
def build_spec(section: str, dx_um: float, drug: str, total_time_s: float,
               timelapse_frames: Optional[int],
               attribute_mechanism: bool) -> Dict:
    """The stock ``GET /api/defaults`` spec (api.default_spec) with dx / drug / window overridden.

    Built by CALLING ``app.server.api.defaults()`` rather than by retyping its
    fields, so a change to the shipped default shows up here as a different
    measurement instead of as a silent divergence between what the Studio runs
    and what this file profiles.
    """
    from app.server.api import default_spec

    spec = dict(default_spec())
    spec["section"] = section
    spec["dx_um"] = float(dx_um)
    spec["drug"] = drug
    spec["total_time_s"] = float(total_time_s)
    if timelapse_frames is not None:
        spec["timelapse_frames"] = int(timelapse_frames)
    spec["attribute_mechanism"] = bool(attribute_mechanism)
    return spec


class StageTimer:
    """A ``_Ctx`` stand-in that records the wall-clock of each server stage.

    Duck-types the three methods ``execute_run`` calls on its ctx (``log``,
    ``stage``, ``check_cancelled``). It cannot cancel, which is the point: a
    profiling harness that could abort a run halfway would report a partial
    time as a whole one.
    """

    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.marks: List = []          # (stage_name, seconds_in_stage)
        # Peak RSS at the END of each stage. Peak, not current: the number an
        # operator needs is the high-water mark the OOM killer sees, and a
        # stage that allocates 300 MiB and frees it before its boundary is
        # still a stage that needed 300 MiB.
        self.rss: List = []            # (stage_name, peak_rss_bytes_at_exit)
        self.logs: List[str] = []
        self.t0 = time.perf_counter()
        self._t_last = self.t0
        self._cur = "pre-stage"

    def log(self, msg: str) -> None:
        self.logs.append(msg)
        if self.verbose:
            print(f"    [log] {msg}", file=sys.stderr)

    def stage(self, name: str) -> None:
        now = time.perf_counter()
        self.marks.append((self._cur, now - self._t_last))
        self.rss.append((self._cur, _maxrss_bytes()))
        self._cur = name
        self._t_last = now
        if self.verbose:
            print(f"  [stage] {name} (t+{now - self.t0:.2f}s)", file=sys.stderr)

    def check_cancelled(self) -> None:
        return None

    def finish(self) -> None:
        now = time.perf_counter()
        self.marks.append((self._cur, now - self._t_last))
        self.rss.append((self._cur, _maxrss_bytes()))
        self._t_last = now

    @property
    def total_s(self) -> float:
        return sum(d for _n, d in self.marks)

    def table(self) -> List[Dict]:
        total = self.total_s or 1.0
        merged: Dict[str, float] = {}
        order: List[str] = []
        for name, d in self.marks:
            if name not in merged:
                merged[name] = 0.0
                order.append(name)
            merged[name] += d
        return [{"stage": n, "seconds": merged[n],
                 "fraction": merged[n] / total} for n in order]


def _fresh_cache():
    """A private, empty reconstruction cache.

    ``RECON_CACHE`` is a process-global LRU; profiling against it would report
    a cache hit as a reconstruction cost of zero on the second call and as 40 s
    on the first, in the same process, for the same spec.
    """
    from app.server.cache import ReconstructionCache

    return ReconstructionCache(maxsize=4)


def _cpu_s() -> float:
    """User+system CPU seconds consumed by this process so far.

    Reported beside every wall-clock because this repo's machines are shared.
    Wall-clock on a box at load 14 is a measurement of the box; CPU time is a
    measurement of the code. Neither one alone is honest: CPU time hides
    contention an operator genuinely pays for, and wall-clock attributes a
    neighbour's build to the solver.
    """
    r = resource.getrusage(resource.RUSAGE_SELF)
    return float(r.ru_utime) + float(r.ru_stime)


def _loadavg() -> str:
    try:
        one, five, fifteen = os.getloadavg()
        return f"{one:.2f}/{five:.2f}/{fifteen:.2f}"
    except (OSError, AttributeError):  # pragma: no cover - platform dependent
        return "unavailable"


def run_once(spec: Dict, cache=None, verbose: bool = False):
    """Execute one run through the server's own entry point. Returns
    ``(body, StageTimer, wall_s, cpu_s, loadavg_str)``."""
    from app.server.jobs import execute_run

    ctx = StageTimer(verbose=verbose)
    cache = cache if cache is not None else _fresh_cache()
    load = _loadavg()
    t0 = time.perf_counter()
    c0 = _cpu_s()
    body = execute_run(spec, ctx, cache=cache)
    wall = time.perf_counter() - t0
    cpu = _cpu_s() - c0
    ctx.finish()
    return body, ctx, wall, cpu, load


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def _grid_of(body: Dict) -> str:
    rep = body.get("report") or {}
    shape = rep.get("grid") or rep.get("grid_shape") or rep.get("shape")
    if shape:
        try:
            return f"{int(shape[0])}x{int(shape[1])} = {int(shape[0]) * int(shape[1]):,} voxels"
        except Exception:  # noqa: BLE001
            return str(shape)
    return "grid_shape not in report"


def _print_header(spec: Dict, body: Optional[Dict] = None) -> None:
    print("RunSpec:", json.dumps(spec, sort_keys=True))
    if body is not None:
        met = body.get("metrics") or {}
        print(f"grid: {_grid_of(body)}")
        print(f"kill_fraction: {met.get('kill_fraction')!r}")


def _print_stage_table(ctx: StageTimer, wall: float, cpu: Optional[float] = None,
                       load: Optional[str] = None) -> None:
    print(f"\n  {'stage':<16}{'seconds':>10}{'fraction':>11}")
    print("  " + "-" * 37)
    for row in ctx.table():
        if row["seconds"] < 5e-4 and row["stage"] == "pre-stage":
            continue
        print(f"  {row['stage']:<16}{row['seconds']:>10.3f}"
              f"{row['fraction'] * 100:>10.1f}%")
    print("  " + "-" * 37)
    print(f"  {'TOTAL':<16}{wall:>10.3f}{100.0:>10.1f}%")
    if ctx.rss:
        print(f"\n  {'stage':<16}{'peak RSS at exit':>20}{'rose by':>12}")
        print("  " + "-" * 48)
        prev = None
        for name, rss in ctx.rss:
            rose = "" if prev is None else f"{_mib(rss - prev):>11.1f}"
            print(f"  {name:<16}{_mib(rss):>19.1f} MiB{rose}")
            prev = rss
    if cpu is not None:
        print(f"  {'cpu (u+s)':<16}{cpu:>10.3f}"
              f"{cpu / wall if wall else float('nan'):>9.2f}x")
    if load is not None:
        print(f"  loadavg at start: {load}   (1/5/15 min; "
              f"{os.cpu_count()} logical CPUs)")


# --------------------------------------------------------------------------- #
# Arms
# --------------------------------------------------------------------------- #
def arm_stages(args) -> int:
    spec = build_spec(args.section, args.dx, args.drug, args.total_time,
                      args.timelapse_frames, not args.no_mechanism)
    cache = _fresh_cache()
    print("=" * 72)
    print(f"STAGES (clean, no profiler)  dx={args.dx} um  section={args.section}")
    print("=" * 72)
    _print_header(spec)

    rows = []
    for i in range(args.repeats):
        label = "COLD (empty recon cache)" if i == 0 else f"WARM #{i} (cache hit)"
        body, ctx, wall, cpu, load = run_once(spec, cache=cache,
                                              verbose=args.verbose)
        print(f"\n--- {label} ---")
        print(f"grid: {_grid_of(body)}")
        print(f"kill_fraction: {(body.get('metrics') or {}).get('kill_fraction')!r}")
        _print_stage_table(ctx, wall, cpu, load)
        print(f"  peak RSS so far: {_mib(_maxrss_bytes()):.1f} MiB")
        rows.append({"label": label, "wall_s": wall, "cpu_s": cpu,
                     "load": load,
                     **{r["stage"]: r["seconds"] for r in ctx.table()}})
    if len(rows) > 1:
        print("\n  repeat spread (the machine, not the code)")
        print(f"  {'':<26}{'min':>9}{'max':>9}{'max/min':>9}")
        for field in ("wall_s", "cpu_s", "simulate", "reconstruct", "metrics"):
            vals = [r[field] for r in rows if r.get(field) is not None]
            vals = [v for v in vals if v > 0]
            if len(vals) < 2:
                continue
            print(f"  {field:<26}{min(vals):>9.3f}{max(vals):>9.3f}"
                  f"{max(vals) / min(vals):>9.2f}x")
    return 0


_PROFILE_TARGETS = (
    # (label, module-substring, function name)
    ("reconstruct_real", "reconstruct/real_pipeline.py", "reconstruct_real"),
    ("build_reconstruction", "app/server/cache.py", "build_reconstruction"),
    ("_max_darcy_speed", "app/server/cache.py", "_max_darcy_speed"),
    ("apply_drug_sensitivity", "reconstruct/real_pipeline.py",
     "apply_drug_sensitivity"),
    ("reconstruction_report", "reconstruct/real_pipeline.py",
     "reconstruction_report"),
    ("calibrate_alpha", "validation/calibration.py", "calibrate_alpha"),
    ("simulate_response", "sim/engine.py", "simulate_response"),
    ("solve_oxygen", "physics/oxygen.py", "solve_oxygen"),
    ("solve_ifp", "physics/pressure.py", "solve_ifp"),
    ("simulate_drug_transport", "physics/diffusion.py",
     "simulate_drug_transport"),
    ("_build_diffusion_operator", "physics/diffusion.py",
     "_build_diffusion_operator"),
    ("_advective_divergence", "physics/diffusion.py", "_advective_divergence"),
    ("apply_drug_kill", "cells/pkpd.py", "apply_drug_kill"),
    ("mechanism.attribute", "sim/mechanism.py", "attribute"),
    ("_viable_history", "sim/engine.py", "_viable_history"),
    ("direct_vessel_distance_um", "validation/vessel_profile.py",
     "direct_vessel_distance_um"),
    ("delivery_metrics.compute_all", "validation/delivery_metrics.py",
     "compute_all"),
    ("compute_reach", "app/server/reach.py", "compute_reach"),
    ("_build_payloads", "app/server/jobs.py", "_build_payloads"),
    ("encode_field", "app/server/payloads.py", "encode_field"),
    ("encode_geometry", "app/server/payloads.py", "encode_geometry"),
    ("encode_timelapse", "app/server/payloads.py", "encode_timelapse"),
    ("encode_mechanism", "app/server/payloads.py", "encode_mechanism"),
    ("jsonable", "app/server/jobs.py", "jsonable"),
)


def _cum_for(stats: pstats.Stats, mod_sub: str, fname: str):
    """(cumtime, tottime, ncalls) summed over every profiled entry whose file
    ends in ``mod_sub`` and whose function is ``fname``. Returns None if the
    function was never called — absence is reported as absence, never as 0.0."""
    cum = tot = 0.0
    ncalls = 0
    found = False
    for (fn, _lineno, name), (cc, nc, tt, ct, _callers) in stats.stats.items():
        if name != fname:
            continue
        if not fn.replace("\\", "/").endswith(mod_sub):
            continue
        found = True
        cum += ct
        tot += tt
        ncalls += nc
    return (cum, tot, ncalls) if found else None


def arm_profile(args) -> int:
    spec = build_spec(args.section, args.dx, args.drug, args.total_time,
                      args.timelapse_frames, not args.no_mechanism)
    cache = _fresh_cache()
    print("=" * 72)
    print(f"cPROFILE  dx={args.dx} um  section={args.section}")
    print("=" * 72)
    _print_header(spec)

    if args.warm:
        print("\n(priming the reconstruction cache first; the profiled run is "
              "therefore a CACHE HIT)")
        run_once(spec, cache=cache)

    prof = cProfile.Profile()
    ctx = StageTimer(verbose=args.verbose)
    from app.server.jobs import execute_run

    load = _loadavg()
    t0 = time.perf_counter()
    c0 = _cpu_s()
    prof.enable()
    body = execute_run(spec, ctx, cache=cache)
    prof.disable()
    wall = time.perf_counter() - t0
    cpu = _cpu_s() - c0
    ctx.finish()

    print(f"grid: {_grid_of(body)}")
    print(f"kill_fraction: {(body.get('metrics') or {}).get('kill_fraction')!r}")
    print(f"\nPROFILED wall-clock: {wall:.3f} s  "
          "(includes cProfile tracing overhead; compare against the `stages` "
          "arm before quoting any absolute time)")
    _print_stage_table(ctx, wall, cpu, load)

    stats = pstats.Stats(prof)
    print(f"\n  {'target':<30}{'cum s':>9}{'self s':>9}{'calls':>10}"
          f"{'cum %':>8}")
    print("  " + "-" * 66)
    for label, mod_sub, fname in _PROFILE_TARGETS:
        got = _cum_for(stats, mod_sub, fname)
        if got is None:
            print(f"  {label:<30}{'not called':>9}")
            continue
        cum, tot, nc = got
        print(f"  {label:<30}{cum:>9.3f}{tot:>9.3f}{nc:>10}"
              f"{cum / wall * 100:>7.1f}%")

    print(f"\n  top {args.top} by SELF time")
    print("  " + "-" * 66)
    stats.sort_stats("tottime")
    stats.stream = sys.stdout
    stats.print_stats(args.top)

    if args.dump:
        prof.dump_stats(args.dump)
        print(f"\nprofile written to {args.dump}")
    return 0


def arm_instrument(args) -> int:
    """Split the transport solve into factorisation / triangular solves /
    advection, which cProfile cannot do (``splu(A).solve`` is a C method).

    The shims wrap calls and change no arithmetic. The clean arm is run FIRST
    in the same process on a primed cache so the two kill fractions can be
    compared bit-for-bit — if the instrumentation changed the answer it would
    not be instrumentation.
    """
    import numpy as np  # noqa: F401  (kept: asserts the import cost is paid early)
    from tmesim.physics import diffusion as diff

    spec = build_spec(args.section, args.dx, args.drug, args.total_time,
                      args.timelapse_frames, not args.no_mechanism)
    cache = _fresh_cache()
    print("=" * 72)
    print(f"INSTRUMENTED TRANSPORT  dx={args.dx} um  section={args.section}")
    print("=" * 72)
    _print_header(spec)

    print("\n--- arm A: clean (cold) ---")
    body_a, ctx_a, wall_a, cpu_a, load_a = run_once(spec, cache=cache,
                                                    verbose=args.verbose)
    print(f"grid: {_grid_of(body_a)}")
    kill_a = (body_a.get("metrics") or {}).get("kill_fraction")
    print(f"kill_fraction: {kill_a!r}")
    _print_stage_table(ctx_a, wall_a, cpu_a, load_a)

    acc = {"factorize_s": 0.0, "solve_s": 0.0, "n_solves": 0,
           "advect_s": 0.0, "n_advect": 0}
    real_factorized = diff.factorized
    real_advect = diff._advective_divergence

    def timed_factorized(a):
        t = time.perf_counter()
        inner = real_factorized(a)
        acc["factorize_s"] += time.perf_counter() - t
        acc["nnz_A"] = int(a.nnz)
        acc["n_A"] = int(a.shape[0])

        def solve(b):
            t0 = time.perf_counter()
            out = inner(b)
            acc["solve_s"] += time.perf_counter() - t0
            acc["n_solves"] += 1
            return out

        return solve

    def timed_advect(*a, **kw):
        t0 = time.perf_counter()
        out = real_advect(*a, **kw)
        acc["advect_s"] += time.perf_counter() - t0
        acc["n_advect"] += 1
        return out

    diff.factorized = timed_factorized
    diff._advective_divergence = timed_advect
    try:
        print("\n--- arm B: instrumented (warm recon cache) ---")
        body_b, ctx_b, wall_b, cpu_b, load_b = run_once(
            spec, cache=cache, verbose=args.verbose)
    finally:
        diff.factorized = real_factorized
        diff._advective_divergence = real_advect

    kill_b = (body_b.get("metrics") or {}).get("kill_fraction")
    print(f"kill_fraction: {kill_b!r}")
    print(f"BIT-EXACT vs arm A: {kill_a == kill_b}"
          + ("" if kill_a == kill_b else "   <-- INSTRUMENTATION CHANGED THE "
                                         "ANSWER; every number below is void"))
    _print_stage_table(ctx_b, wall_b, cpu_b, load_b)

    sim_s = next((r["seconds"] for r in ctx_b.table() if r["stage"] == "simulate"),
                 float("nan"))
    print(f"\n  transport internals (instrumented, simulate stage = {sim_s:.3f} s)")
    print("  " + "-" * 60)
    print(f"  {'sparse LU factorisation':<34}{acc['factorize_s']:>9.3f} s"
          f"{acc['factorize_s'] / sim_s * 100:>8.1f}% of simulate")
    print(f"  {'triangular solves (n=' + str(acc['n_solves']) + ')':<34}"
          f"{acc['solve_s']:>9.3f} s{acc['solve_s'] / sim_s * 100:>8.1f}%")
    print(f"  {'upwind advection (n=' + str(acc['n_advect']) + ')':<34}"
          f"{acc['advect_s']:>9.3f} s{acc['advect_s'] / sim_s * 100:>8.1f}%")
    rest = sim_s - acc["factorize_s"] - acc["solve_s"] - acc["advect_s"]
    print(f"  {'everything else in simulate':<34}{rest:>9.3f} s"
          f"{rest / sim_s * 100:>8.1f}%")
    if "n_A" in acc:
        print(f"\n  operator A: {acc['n_A']:,} x {acc['n_A']:,}, "
              f"nnz = {acc['nnz_A']:,}")
    if acc["n_solves"]:
        print(f"  per-solve mean: {acc['solve_s'] / acc['n_solves'] * 1e3:.2f} ms")
    return 0


def arm_baseline(args) -> int:
    """The RSS FLOOR: what the process costs before any grid exists.

    The 512 MB ceiling is usually discussed as if it were spent on the grid.
    It is not spent only on the grid, and this arm establishes how much of it
    is gone before ``dx`` has been read — import by import, each one timed and
    weighed on its own. Every step is a real import of a real production
    module, in the order ``execute_run`` reaches them.
    """
    steps = [
        ("numpy", "import numpy"),
        ("scipy.sparse.linalg", "import scipy.sparse.linalg"),
        ("pandas", "import pandas"),
        ("tmesim.core", "import tmesim.core"),
        ("tmesim.physics.diffusion", "import tmesim.physics.diffusion"),
        ("tmesim.sim.engine", "import tmesim.sim.engine"),
        ("tmesim.reconstruct.real_pipeline",
         "import tmesim.reconstruct.real_pipeline"),
        ("tmesim.molecular.library", "import tmesim.molecular.library"),
        ("tmesim.validation.calibration", "import tmesim.validation.calibration"),
        ("tmesim.validation.delivery_metrics",
         "import tmesim.validation.delivery_metrics"),
        ("app.server.cache", "import app.server.cache"),
        ("app.server.jobs", "import app.server.jobs"),
        ("app.server.api", "import app.server.api"),
        ("api.list_sections()",
         "import app.server.api as _a; _a.list_sections()"),
        ("api.default_spec()", "import app.server.api as _a; _a.default_spec()"),
        ("calibrate_alpha()",
         "from tmesim.validation.calibration import calibrate_alpha; "
         "calibrate_alpha()"),
    ]
    print("=" * 72)
    print("RSS FLOOR — cost before any grid exists")
    print("=" * 72)
    print(f"  {'step':<38}{'peak RSS':>12}{'rose by':>11}{'secs':>8}")
    print("  " + "-" * 69)
    prev = _maxrss_bytes()
    print(f"  {'(interpreter + stdlib)':<38}{_mib(prev):>11.1f} M")
    ns: Dict = {}
    for label, code in steps:
        t0 = time.perf_counter()
        exec(compile(code, f"<{label}>", "exec"), ns)  # noqa: S102
        dt = time.perf_counter() - t0
        now = _maxrss_bytes()
        print(f"  {label:<38}{_mib(now):>11.1f} M{_mib(now - prev):>10.1f}"
              f"{dt:>8.2f}")
        prev = now
    print("  " + "-" * 69)
    print(f"  FLOOR before any RunSpec executes: {_mib(prev):.1f} MiB")
    return 0


def arm_memory(args) -> int:
    """One run, this process, peak RSS. Launch one per process."""
    rss_import = _maxrss_bytes()
    spec = build_spec(args.section, args.dx, args.drug, args.total_time,
                      args.timelapse_frames, not args.no_mechanism)
    rss_ready = _maxrss_bytes()
    t0 = time.perf_counter()
    body, ctx, wall, cpu, load = run_once(spec, verbose=args.verbose)
    rss_peak = _maxrss_bytes()

    rep = body.get("report") or {}
    shape = rep.get("grid") or rep.get("grid_shape") or rep.get("shape") or [0, 0]
    try:
        voxels = int(shape[0]) * int(shape[1])
    except Exception:  # noqa: BLE001
        voxels = 0
    payload_bytes = sum(len(v) for v in (body.get("_payloads") or {}).values())
    out = {
        "section": args.section,
        "dx_um": float(args.dx),
        "drug": args.drug,
        "total_time_s": float(args.total_time),
        "grid": [int(shape[0]), int(shape[1])],
        "voxels": voxels,
        "n_cells": int(rep.get("n_cells") or 0),
        "wall_s": round(wall, 3),
        "cpu_s": round(cpu, 3),
        "loadavg_at_start": load,
        "rss_after_import_MiB": round(_mib(rss_import), 1),
        "rss_before_run_MiB": round(_mib(rss_ready), 1),
        "rss_peak_MiB": round(_mib(rss_peak), 1),
        "rss_run_delta_MiB": round(_mib(rss_peak - rss_ready), 1),
        "payload_bytes": payload_bytes,
        "rss_by_stage_MiB": {n: round(_mib(v), 1) for n, v in ctx.rss},
        "kill_fraction": (body.get("metrics") or {}).get("kill_fraction"),
        "stages": {r["stage"]: round(r["seconds"], 3) for r in ctx.table()},
    }
    del body
    if args.json:
        print("@@JSON@@" + json.dumps(out))
    else:
        print(json.dumps(out, indent=2))
    return 0


def arm_sweep(args) -> int:
    """One fresh subprocess per point, so every peak RSS is a cold one."""
    dxs = [float(x) for x in str(args.dx).split(",") if x.strip()]
    sections = [s.strip() for s in str(args.section).split(",") if s.strip()]
    rows: List[Dict] = []
    print("=" * 96)
    print("PEAK-RSS SWEEP — one cold subprocess per point")
    print("=" * 96)
    for sec in sections:
        for dx in dxs:
            cmd = [sys.executable, os.path.abspath(__file__), "memory",
                   "--section", sec, "--dx", str(dx), "--drug", args.drug,
                   "--total-time", str(args.total_time), "--json"]
            if args.timelapse_frames is not None:
                cmd += ["--timelapse-frames", str(args.timelapse_frames)]
            if args.no_mechanism:
                cmd += ["--no-mechanism"]
            t0 = time.perf_counter()
            proc = subprocess.run(cmd, cwd=_REPO, capture_output=True,
                                  text=True, timeout=args.timeout)
            elapsed = time.perf_counter() - t0
            line = next((ln for ln in proc.stdout.splitlines()
                         if ln.startswith("@@JSON@@")), None)
            if line is None:
                print(f"  {sec} dx={dx}: FAILED rc={proc.returncode} "
                      f"after {elapsed:.1f}s")
                tail = (proc.stderr or "").strip().splitlines()[-4:]
                for t in tail:
                    print(f"      {t}")
                rows.append({"section": sec, "dx_um": dx, "failed": True,
                             "rc": proc.returncode})
                continue
            row = json.loads(line[len("@@JSON@@"):])
            row["subprocess_wall_s"] = round(elapsed, 2)
            rows.append(row)
            print(f"  {sec:<16} dx={dx:>5} um  "
                  f"{row['voxels']:>9,} voxels  "
                  f"peak {row['rss_peak_MiB']:>7.1f} MiB  "
                  f"(run delta {row['rss_run_delta_MiB']:>7.1f})  "
                  f"wall {row['wall_s']:>7.2f} s")
    print("\nCSV")
    print("section,dx_um,voxels,n_cells,rss_peak_MiB,rss_run_delta_MiB,"
          "rss_before_run_MiB,wall_s,cpu_s,simulate_s,payload_bytes,"
          "kill_fraction,loadavg_at_start")
    for r in rows:
        if r.get("failed"):
            print(f"{r['section']},{r['dx_um']},FAILED,,,,,,,,,,")
            continue
        print(f"{r['section']},{r['dx_um']},{r['voxels']},{r['n_cells']},"
              f"{r['rss_peak_MiB']},{r['rss_run_delta_MiB']},"
              f"{r['rss_before_run_MiB']},{r['wall_s']},{r['cpu_s']},"
              f"{r['stages'].get('simulate')},{r['payload_bytes']},"
              f"{r['kill_fraction']},{r['loadavg_at_start']}")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(rows, fh, indent=2)
        print(f"\nrows written to {args.out}")
    return 0


# --------------------------------------------------------------------------- #
def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="arm", required=True)

    def common(sp):
        sp.add_argument("--section", default="xenium_319")
        sp.add_argument("--dx", default="12")
        sp.add_argument("--drug", default="temozolomide")
        sp.add_argument("--total-time", dest="total_time", type=float,
                        default=3600.0 * 6)
        sp.add_argument("--timelapse-frames", dest="timelapse_frames",
                        type=int, default=None)
        sp.add_argument("--no-mechanism", action="store_true")
        sp.add_argument("--verbose", "-v", action="store_true")
        return sp

    sp = common(sub.add_parser("stages"))
    sp.add_argument("--repeats", type=int, default=2)
    sp.set_defaults(func=arm_stages, dx_is_float=True)

    sp = common(sub.add_parser("profile"))
    sp.add_argument("--top", type=int, default=25)
    sp.add_argument("--warm", action="store_true")
    sp.add_argument("--dump", default=None)
    sp.set_defaults(func=arm_profile)

    sp = common(sub.add_parser("instrument"))
    sp.set_defaults(func=arm_instrument)

    sp = common(sub.add_parser("baseline"))
    sp.set_defaults(func=arm_baseline)

    sp = common(sub.add_parser("memory"))
    sp.add_argument("--json", action="store_true")
    sp.set_defaults(func=arm_memory)

    sp = common(sub.add_parser("sweep"))
    sp.add_argument("--timeout", type=float, default=3600.0)
    sp.add_argument("--out", default=None)
    sp.set_defaults(func=arm_sweep)

    args = p.parse_args(argv)
    if args.arm != "sweep":
        args.dx = float(args.dx)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
