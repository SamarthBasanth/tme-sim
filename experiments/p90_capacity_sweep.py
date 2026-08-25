"""DOES ``E_p90/E_p10`` — THE HEADLINE — SURVIVE A TISSUE-CAPACITY TERM?

WHY THIS EXISTS
---------------
``E_p90/E_p10`` (``tmesim.validation.delivery_metrics.exposure_percentiles``) is
the product's headline delivery number and **its capacity stability has never
been measured**. Two of its supporting arguments failed in the two days before
this file was written (``results/reach_invariances_measure.md`` sec. 7.1 retracted
the "~10x grid sensitivity" of the viable denominator, measured at 0.6 %;
``results/spatial_noise_decomposition.md`` showed distance binning was never the
mechanism, so ``results/delivery_spread_noise.md``'s "blind to binning" argument
was qualified). Meanwhile its nearest MEASURED relatives on the same section are
severely capacity-sensitive: the normalised falloff moves 0.3772 and the fitted
penetration length collapses 279.715 -> 76.740 um by R = 39
(``results/reach_invariances_measure.md`` sec. 4), and that file states outright
that nothing in it licenses any claim about ``E_p90/E_p10``.

The pre-declaration — expectations, failure thresholds, and the verdict rule —
was written to ``results/p90_capacity_sweep.md`` sec. 0 BEFORE this script was
run. The thresholds are inherited verbatim from ``delivery_spread_noise.md``'s
own measured nuisance envelope; none was chosen after seeing capacity data.

THE TWO ARMS, AND THE TRAP BETWEEN THEM
---------------------------------------
``R`` multiplies the whole AUC field and therefore cancels EXACTLY from a
percentile ratio. That is not the question. Capacity rescales TIME, so at a fixed
wall-clock window a capacity-R tissue has only reached the state the R = 1 tissue
reached at ``T/R``; the finite-window SHAPE changes and P90/P10 reads the shape.

* **ARM A — truncated window (CORRECT).** Run the shipping R = 1 model at
  ``(dt/R, T/R)``. The solver's step is ``c_new = A^-1[c + h S - h adv(c)]`` with
  ``A = I - h L + h k`` (``tmesim/physics/diffusion.py`` :377, :447, :452), so the
  capacity equation ``R dC/dt = div(D grad C) - div(v C) + S - k C`` discretised
  at timestep ``dt`` IS the R = 1 solver at ``h = dt/R``, iterate for iterate —
  **and the advection term is retarded by construction**, because the same ``h``
  multiplies it. Step count unchanged (2,160), so no time resolution is traded
  away.
* **ARM B — the published parameter-division emulation, CONTAMINATED WITH
  CONVECTION ON.** ``capacity_ranking_robustness.retard`` divides ``diffusivity``,
  ``binding``, ``clearance`` and ``bbb_penetrance`` by R and **cannot divide the
  velocity field**, which is built inside ``sim.engine.simulate_response``.
  Multiply the emulated equation through by R and the advection term keeps its R:
  with convection ON this is a capacity change PLUS an R-fold Peclet boost, not a
  capacity change. ``capacity_ranking_robustness.md`` sec. 2d runs convection OFF
  for exactly that reason. The shipping config is convection ON.

Both are measured. Every printed number names its arm.

WHAT IS CHANGED
---------------
Nothing. ``calibrate_alpha()`` is asserted live at ``22.223477711172563`` at the
TOP and BOTTOM of every stage. Nothing under ``tmesim/`` or ``app/`` is written
to or modified. No parameter is tuned toward any anchor.

Run:
    ./.venv/bin/python experiments/p90_capacity_sweep.py --stage guards
    ./.venv/bin/python experiments/p90_capacity_sweep.py --stage plant
    ./.venv/bin/python experiments/p90_capacity_sweep.py --stage identity
    ./.venv/bin/python experiments/p90_capacity_sweep.py --stage truncated
    ./.venv/bin/python experiments/p90_capacity_sweep.py --stage emulation
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.core import OxygenConfig, PDEConfig, PKConfig, SimConfig
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.physics import diffusion as _diffusion
from tmesim.sim.engine import simulate_response
from tmesim.validation import vessel_profile as VP
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.validation.delivery_metrics import (
    exposure_percentiles,
    tumor_cell_exposure,
)


def _load_published_retard():
    """``capacity_ranking_robustness.retard``, loaded from its file by path.

    Imported rather than re-implemented so that Arm B is literally the published
    emulation — the arm this task's brief says has already caught two people —
    and not a paraphrase of it that could quietly differ. ``experiments/`` is not
    a package, so it is loaded by explicit path rather than by name.
    """
    import importlib.util
    p = os.path.join(_ROOT, "experiments", "capacity_ranking_robustness.py")
    spec = importlib.util.spec_from_file_location(
        "_published_capacity_ranking_robustness", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.retard


retard = _load_published_retard()

# --------------------------------------------------------------------------- #
# Frozen state. Asserted against LIVE code, never adjusted.
# --------------------------------------------------------------------------- #
SHIPPED_ALPHA = 22.223477711172563
FROZEN_PDE_DEFAULTS = dict(dt=0.5, total_time=43200.0, convection=True,
                           save_every=200)

#: The flagship configuration ``results/delivery_spread_noise.md`` sec. 2 reports
#: 2.111 / 3.942 at: dx = 12 um, mask repair = 2, dt = 10 s, 6 h, convection ON.
DX_UM = 12.0
REPAIR = 2
DT_S = 10.0
WINDOW_S = 3600.0 * 6
CONVECTION = True
SAVE_EVERY = 10_000              # history is not read here; keep it out of memory
SEED = 0

PAIR: Sequence[str] = ("temozolomide", "doxorubicin")

#: The ladder the task fixes and the sibling capacity files use. 3.2 / 39 / 4184
#: are the minimum, best-available constant and maximum of the 33 measured rat
#: brain-slice Vu,brain values (``results/brain_slice_anchor.md``). 1 is shipping.
R_LADDER: Sequence[float] = (1.0, 3.2, 8.0, 20.0, 39.0, 60.0, 80.0, 120.0,
                             620.0, 4184.0)
R_ANCHOR = 39.0
R_MEASURED_MIN = 3.2
R_MEASURED_MAX = 4184.0

#: Capacities at which the two routes are cross-checked against each other.
R_IDENTITY: Sequence[float] = (39.0,)

# --------------------------------------------------------------------------- #
# INHERITED thresholds. Transcribed from results/delivery_spread_noise.md
# secs. 0 and 2 (30 solves: 5 grid/mask reconstructions x 3 windows x 2
# compounds). NOT chosen for this study, and pre-declared in
# results/p90_capacity_sweep.md sec. 0.3 before this script was run.
# --------------------------------------------------------------------------- #
NUISANCE = {
    "temozolomide": {"shipped": 2.111, "lo": 1.969, "hi": 2.171},
    "doxorubicin": {"shipped": 3.942, "lo": 3.770, "hi": 3.996},
}
for _n, _d in NUISANCE.items():
    _d["width"] = _d["hi"] - _d["lo"]
#: The worst-case doxorubicin/temozolomide separation anywhere on that nuisance
#: grid (``delivery_spread_noise.md`` sec. 0). C2's bar.
SEPARATION_FLOOR = 1.74

#: Tolerance for "the two routes agree to solver precision" in the identity
#: cross-check. 1e-9 relative is ~7 orders looser than the 3.331e-16 the sibling
#: study measured on the reach falloff and ~2 orders tighter than any effect this
#: study reports, so it can neither pass a real disagreement nor fail on rounding.
ROUTE_REL_TOL = 1e-9

_RESULTS = os.path.join(_ROOT, "results")
_DATA = os.path.join(_RESULTS, "p90_capacity_sweep.json")


# --------------------------------------------------------------------------- #
# Guards. Each is planted against in --stage plant.
# --------------------------------------------------------------------------- #
def _assert_alpha(where: str) -> float:
    """The calibration, checked against LIVE code. Top AND bottom of every stage."""
    a = calibrate_alpha()
    if a != SHIPPED_ALPHA:
        raise AssertionError(
            f"[{where}] calibrate_alpha() = {a!r}, expected {SHIPPED_ALPHA!r}")
    return a


def _assert_frozen_pde() -> None:
    d = PDEConfig()
    for k, v in FROZEN_PDE_DEFAULTS.items():
        if getattr(d, k) != v:
            raise AssertionError(
                f"PDEConfig().{k} = {getattr(d, k)!r}, expected {v!r}")
    if d.time_resolved is not False:
        raise AssertionError(
            "PDEConfig().time_resolved is not False: the source would no longer "
            "be a single scalar and the truncated-window identity would fail")
    if d.transvascular_exchange is not False:
        raise AssertionError("PDEConfig().transvascular_exchange is not False")


def _assert_population_invariant(rows: Sequence[Dict]) -> None:
    """The tumour population must be identical at every R.

    ``E_p90/E_p10`` is a percentile over a cell population. If the population
    moved, the values at two R are percentiles of different things and nothing
    may be compared elementwise.
    """
    if not rows:
        return
    n0 = rows[0]["n_tumor"]
    for r in rows:
        if r["n_tumor"] != n0:
            raise AssertionError(
                f"tumour population moved at R = {r['R']} ({r['arm']}): "
                f"n_tumor = {r['n_tumor']} vs {n0} — the percentiles are of "
                "different populations and are not comparable")


def _assert_steps_invariant(rows: Sequence[Dict]) -> None:
    """The truncation must not trade time resolution for window length."""
    if not rows:
        return
    s0 = rows[0]["n_steps"]
    for r in rows:
        if r["n_steps"] != s0:
            raise AssertionError(
                f"step count moved at R = {r['R']} ({r['arm']}): "
                f"{r['n_steps']} vs {s0} — the arms differ in time resolution "
                "as well as in capacity, so the comparison is confounded")


def _assert_ratio_defined(rows: Sequence[Dict]) -> None:
    """``E_p10`` must stay strictly positive and the ratio finite.

    At R = 4184 the window is ~5 s. If the worst-served decile underflows to
    zero the ratio silently becomes ``inf``/``nan`` and a broken run would be
    reported as an enormous spread.
    """
    for r in rows:
        if not (r["E_p10_uM"] > 0.0):
            raise AssertionError(
                f"E_p10 = {r['E_p10_uM']!r} at R = {r['R']} ({r['arm']}): the "
                "ratio is undefined and must not be reported as a spread")
        if not np.isfinite(r["E_p90_over_p10"]):
            raise AssertionError(
                f"E_p90/E_p10 = {r['E_p90_over_p10']!r} at R = {r['R']} "
                f"({r['arm']})")


def _routes_agree(a: Dict, b: Dict, tol: float = ROUTE_REL_TOL) -> Dict:
    """Do two routes to the same capacity give the same P90/P10?

    Returns the comparison; ``agree`` is the boolean the identity stage asserts
    on. This is the check the whole study's instrument choice rests on, so it is
    planted against directly: fed the convection-ON emulation it MUST report
    ``agree = False``, and fed the convection-OFF one it MUST report True.
    """
    x, y = a["E_p90_over_p10"], b["E_p90_over_p10"]
    rel = abs(x - y) / abs(y) if y else float("inf")
    return {"a": float(x), "b": float(y), "rel_diff": float(rel),
            "tol": float(tol), "agree": bool(rel <= tol)}


# --------------------------------------------------------------------------- #
# The non-negativity clamp watcher — the one non-linearity a rescale cannot see
# --------------------------------------------------------------------------- #
class ClampWatch:
    """Count firings of ``np.maximum(c_new, 0.0, out=c_new)``.

    ``max(a*x, 0) = a*max(x, 0)`` for ``a > 0``, so the clamp is positively
    homogeneous and a rescale test survives it intact whether it fires or not.
    It has to be instrumented directly. The shim forwards everything else to the
    real numpy and is removed unconditionally.
    """

    def __init__(self) -> None:
        self.clamp_calls = 0
        self.steps_with_negatives = 0
        self.negative_entries = 0
        self.worst_undershoot = 0.0
        self._real = None

    def __enter__(self) -> "ClampWatch":
        real = _diffusion.np
        self._real = real
        watch = self

        class _Shim:
            def __getattr__(self, k):
                return getattr(real, k)

            def maximum(self, a, b, out=None, **kw):
                # Only the clamp is an in-place maximum against scalar 0.0;
                # _effective_diffusivity's maximum has out=None and a non-zero
                # floor, so it is not counted.
                if out is not None and real.ndim(b) == 0 and float(b) == 0.0:
                    watch.clamp_calls += 1
                    mn = float(real.min(a))
                    if mn < 0.0:
                        watch.steps_with_negatives += 1
                        watch.negative_entries += int(
                            real.count_nonzero(real.asarray(a) < 0.0))
                        watch.worst_undershoot = min(watch.worst_undershoot, mn)
                return real.maximum(a, b, out=out, **kw)

        _diffusion.np = _Shim()
        return self

    def __exit__(self, *exc) -> None:
        _diffusion.np = self._real

    def as_dict(self) -> Dict:
        return {"clamp_calls": self.clamp_calls,
                "steps_with_negatives": self.steps_with_negatives,
                "negative_entries": self.negative_entries,
                "worst_undershoot": self.worst_undershoot,
                "clamp_ever_fired": bool(self.negative_entries)}


# --------------------------------------------------------------------------- #
# One solve
# --------------------------------------------------------------------------- #
_MICRO: Dict[str, object] = {}


def micro() -> object:
    if "m" not in _MICRO:
        t0 = time.time()
        m, _cn, _cr = VP.build_arm_microenvironment(dx_um=DX_UM, repair=REPAIR)
        _MICRO["m"] = m
        print(f"  built dx = {DX_UM:g} um repair = {REPAIR}, grid "
              f"{m.domain.shape}, {m.cells.n} cells ({time.time() - t0:.0f}s)",
              flush=True)
    return _MICRO["m"]


def compile_flagship(name: str, alpha: float):
    """The pipeline the shipped real-tissue panel uses."""
    return apply_barrier_to_drugparams(calibrated_drugparams(name, alpha),
                                       get_molecule(name))


def spread(dp, dt: float, total_time: float, convection: bool,
           watch: Optional[ClampWatch] = None) -> Dict:
    """One P90/P10 delivery spread, through the same calls the panel makes."""
    m = micro()
    pde = PDEConfig(dt=dt, total_time=total_time, convection=convection,
                    save_every=SAVE_EVERY)
    t0 = time.time()
    if watch is not None:
        with watch:
            res = simulate_response(m, SimConfig(drug=dp, pde=pde, seed=SEED))
    else:
        res = simulate_response(m, SimConfig(drug=dp, pde=pde, seed=SEED))
    exp = tumor_cell_exposure(m, res)
    pct = exposure_percentiles(exp["exposure_uM"], exp["auc_uM_s"])
    lo, hi = pct["auc_p10_uM_s"], pct["auc_p90_uM_s"]
    return {
        "dt_s": float(dt), "total_time_s": float(total_time),
        "convection": bool(convection),
        "n_steps": int(round(total_time / dt)),
        "E_p90_over_p10": float(pct["E_p90_over_p10"]),
        "auc_p90_over_p10": float(hi / lo) if lo > 0 else float("nan"),
        "E_p10_uM": float(pct["E_p10_uM"]),
        "E_p50_uM": float(pct["E_p50_uM"]),
        "E_p90_uM": float(pct["E_p90_uM"]),
        "auc_p10_uM_s": float(pct["auc_p10_uM_s"]),
        "auc_p50_uM_s": float(pct["auc_p50_uM_s"]),
        "auc_p90_uM_s": float(pct["auc_p90_uM_s"]),
        "n_tumor": int(exp["n_tumor"]), "T_s": float(exp["T_s"]),
        "elapsed_s": time.time() - t0,
    }


def arm_a(dp, R: float) -> Dict:
    """ARM A — truncated window ``(dt/R, T/R)``, convection ON. The CORRECT
    capacity change: the advection term is retarded by the same ``h``."""
    r = spread(dp, DT_S / R, WINDOW_S / R, CONVECTION)
    r.update({"R": float(R), "arm": "A_truncated_window_convection_ON",
              "arm_is_capacity": True})
    return r


def arm_b(dp, R: float, convection: bool = CONVECTION) -> Dict:
    """ARM B — the published ``retard`` parameter-division emulation at the
    shipping ``(dt, T)``. With convection ON this is a capacity change PLUS an
    R-fold Peclet boost, because the velocity field cannot be divided by R."""
    r = spread(retard(dp, R), DT_S, WINDOW_S, convection)
    r.update({"R": float(R),
              "arm": ("B_parameter_division_convection_ON" if convection
                      else "B_parameter_division_convection_OFF"),
              "arm_is_capacity": not convection})
    return r


# --------------------------------------------------------------------------- #
# Stages
# --------------------------------------------------------------------------- #
def stage_plant() -> Dict:
    """Plant the defect each guard claims to catch and watch it bite.

    Nothing here is a simulation. These are the cheap guards; the expensive one
    (``_routes_agree``) is planted against with real solves in ``stage_identity``.
    """
    import tmesim.validation.calibration as _cal

    out: List[Dict] = []

    def _bit(name: str, fn) -> None:
        try:
            fn()
        except AssertionError as exc:
            out.append({"guard": name, "bit": True, "message": str(exc)})
            print(f"  PLANT [{name}]: BIT -> {exc}", flush=True)
            return
        out.append({"guard": name, "bit": False, "message": ""})
        print(f"  PLANT [{name}]: DID NOT BITE — the guard is vacuous", flush=True)

    ok = {"R": 1.0, "arm": "plant", "n_tumor": 23812, "n_steps": 2160,
          "E_p10_uM": 1.0, "E_p90_over_p10": 2.0}

    # 1. alpha
    real_alpha = _cal.calibrate_alpha
    try:
        _cal.calibrate_alpha = lambda *a, **k: SHIPPED_ALPHA + 1e-9
        globals()["calibrate_alpha"] = _cal.calibrate_alpha
        _bit("_assert_alpha", lambda: _assert_alpha("plant"))
    finally:
        _cal.calibrate_alpha = real_alpha
        globals()["calibrate_alpha"] = real_alpha
    # and it must NOT bite once restored
    _assert_alpha("plant/restored")
    print("  PLANT [_assert_alpha]: passes again once restored", flush=True)

    # 2. population
    bad = dict(ok, n_tumor=23811)
    _bit("_assert_population_invariant",
         lambda: _assert_population_invariant([ok, bad]))
    _assert_population_invariant([ok, dict(ok)])

    # 3. step count
    bad = dict(ok, n_steps=2159)
    _bit("_assert_steps_invariant",
         lambda: _assert_steps_invariant([ok, bad]))
    _assert_steps_invariant([ok, dict(ok)])

    # 4. E_p10 underflow
    bad = dict(ok, E_p10_uM=0.0)
    _bit("_assert_ratio_defined", lambda: _assert_ratio_defined([bad]))
    _assert_ratio_defined([ok])

    # 5. _routes_agree, on synthetic inputs. The real, solved plant is in
    #    stage_identity; this only proves the comparator is not a rubber stamp.
    same = _routes_agree({"E_p90_over_p10": 2.0}, {"E_p90_over_p10": 2.0})
    diff = _routes_agree({"E_p90_over_p10": 2.0}, {"E_p90_over_p10": 2.001})
    rubber = (not same["agree"]) or diff["agree"]
    out.append({"guard": "_routes_agree(synthetic)", "bit": not rubber,
                "message": f"identical -> agree={same['agree']}, "
                           f"0.05% apart -> agree={diff['agree']}"})
    print(f"  PLANT [_routes_agree]: identical agree={same['agree']}, "
          f"0.05%-apart agree={diff['agree']} -> "
          f"{'BIT' if not rubber else 'VACUOUS'}", flush=True)

    vacuous = [r["guard"] for r in out if not r["bit"]]
    if vacuous:
        raise AssertionError(f"guards that did not bite: {vacuous}")
    return {"plants": out, "all_bit": True}


def stage_identity(alpha: float) -> Dict:
    """Validate Arm A against Arm B where Arm B is exact, and MEASURE the trap.

    (a) The clamp, at both ends of Arm A's sweep.
    (b) convection OFF: the two routes MUST agree to solver precision. This is
        what validates this implementation of the truncated route.
    (c) convection ON: the two routes MUST NOT agree — this is the planted
        defect (using Arm B as a capacity change at the shipping config), and
        the size of the disagreement is the size of the Peclet contamination.
    """
    out: Dict = {"clamp": [], "routes_off": [], "routes_on": []}
    dp = compile_flagship("temozolomide", alpha)

    for label, dt, T in (("ArmA R=1 (dt 10 s, 6 h) = shipping", DT_S, WINDOW_S),
                         (f"ArmA R={R_MEASURED_MAX:g} truncated",
                          DT_S / R_MEASURED_MAX, WINDOW_S / R_MEASURED_MAX)):
        w = ClampWatch()
        r = spread(dp, dt, T, CONVECTION, watch=w)
        rec = {"label": label, "dt_s": dt, "total_time_s": T, **w.as_dict(),
               "E_p90_over_p10": r["E_p90_over_p10"]}
        out["clamp"].append(rec)
        print(f"  clamp [{label}]: {rec['clamp_calls']} calls, "
              f"{rec['negative_entries']} negative entries, worst "
              f"{rec['worst_undershoot']:.3e} -> "
              f"{'FIRED' if rec['clamp_ever_fired'] else 'never fired'}",
              flush=True)

    for R in R_IDENTITY:
        # (b) convection OFF — Arm B is exact here; the routes must agree.
        a_off = spread(dp, DT_S / R, WINDOW_S / R, False)
        b_off = arm_b(dp, R, convection=False)
        cmp_off = _routes_agree(a_off, b_off)
        cmp_off["R"] = float(R)
        out["routes_off"].append(cmp_off)
        print(f"  identity R={R:g} convection OFF: truncated "
              f"{cmp_off['a']:.9f} vs emulation {cmp_off['b']:.9f}  rel "
              f"{cmp_off['rel_diff']:.3e} -> agree={cmp_off['agree']}",
              flush=True)
        if not cmp_off["agree"]:
            raise AssertionError(
                f"the two routes DISAGREE at R = {R} with convection OFF "
                f"(rel {cmp_off['rel_diff']:.3e}), where the parameter-division "
                "emulation is exact. Arm A's implementation is not validated "
                "and no number in this study may be reported.")

        # (c) convection ON — THE PLANTED DEFECT. Using Arm B as if it were a
        #     capacity change at the shipping config. The check must bite.
        a_on = arm_a(dp, R)
        b_on = arm_b(dp, R, convection=True)
        cmp_on = _routes_agree(a_on, b_on)
        cmp_on["R"] = float(R)
        cmp_on["plant_bit"] = not cmp_on["agree"]
        out["routes_on"].append(cmp_on)
        print(f"  PLANT R={R:g} convection ON (Arm B used as capacity): "
              f"truncated {cmp_on['a']:.6f} vs emulation {cmp_on['b']:.6f}  rel "
              f"{cmp_on['rel_diff']:.3e} -> "
              f"{'BIT' if cmp_on['plant_bit'] else 'DID NOT BITE'}", flush=True)
        if not cmp_on["plant_bit"]:
            raise AssertionError(
                f"the planted defect did not bite at R = {R}: the convection-ON "
                "emulation agreed with the correct route to within "
                f"{ROUTE_REL_TOL:g}, so this study's central instrument claim "
                "is untested.")
    return out


def stage_sweep(alpha: float, arm: str) -> Dict:
    """The ladder, one arm, both compounds."""
    fn = arm_a if arm == "A" else arm_b
    drugs = {n: compile_flagship(n, alpha) for n in PAIR}
    rows: List[Dict] = []
    for R in R_LADDER:
        for name in PAIR:
            r = fn(drugs[name], R)
            r["drug"] = name
            rows.append(r)
            print(f"  [arm {arm}] R = {R:<8g} {name:<14} "
                  f"E_p90/E_p10 = {r['E_p90_over_p10']:10.4f}  "
                  f"E_p10 = {r['E_p10_uM']:.4g} uM  n_tumor={r['n_tumor']}  "
                  f"({r['elapsed_s']:.0f}s)", flush=True)
            _dump({f"partial_arm_{arm}": rows})
    for name in PAIR:
        sub = [r for r in rows if r["drug"] == name]
        _assert_population_invariant(sub)
        _assert_steps_invariant(sub)
        _assert_ratio_defined(sub)
    return {"arm": arm, "rows": rows,
            "r_ladder": [float(v) for v in R_LADDER],
            "route": ("truncated window (dt/R, T/R), convection ON — a genuine "
                      "capacity change, advection retarded by construction"
                      if arm == "A" else
                      "parameter division (D, binding, clearance, "
                      "bbb_penetrance)/R at the shipping (dt, T), convection ON "
                      "— a capacity change PLUS an R-fold Peclet boost, NOT a "
                      "capacity change")}


def analyse(rows: List[Dict]) -> Dict:
    """C1, C2 and C3 exactly as pre-declared in results/p90_capacity_sweep.md."""
    by = {n: {r["R"]: r for r in rows if r["drug"] == n} for n in PAIR}
    per_R: List[Dict] = []
    for R in R_LADDER:
        rec: Dict = {"R": float(R)}
        for n in PAIR:
            v = by[n][R]["E_p90_over_p10"]
            rec[n] = float(v)
            rec[f"{n}_move_vs_R1"] = float(v - by[n][1.0]["E_p90_over_p10"])
            rec[f"{n}_fold_vs_R1"] = float(v / by[n][1.0]["E_p90_over_p10"])
        rec["separation"] = float(rec[PAIR[1]] / rec[PAIR[0]])
        per_R.append(rec)

    # C1 — value stability at R = 39 against the inherited nuisance width.
    c1: Dict = {"R": R_ANCHOR, "per_drug": {}, "fail": False}
    for n in PAIR:
        move = abs(by[n][R_ANCHOR]["E_p90_over_p10"]
                   - by[n][1.0]["E_p90_over_p10"])
        d = {"R1": by[n][1.0]["E_p90_over_p10"],
             "R39": by[n][R_ANCHOR]["E_p90_over_p10"],
             "abs_move": float(move),
             "nuisance_width": NUISANCE[n]["width"],
             "times_nuisance_width": float(move / NUISANCE[n]["width"]),
             "fail": bool(move > NUISANCE[n]["width"])}
        c1["per_drug"][n] = d
        c1["fail"] = c1["fail"] or d["fail"]

    # C2 — separation over the MEASURED capacity range only.
    meas = [p for p in per_R if R_MEASURED_MIN <= p["R"] <= R_MEASURED_MAX]
    worst = min(meas, key=lambda p: p["separation"])
    c2 = {"floor": SEPARATION_FLOOR,
          "min_separation": worst["separation"], "at_R": worst["R"],
          "range_lo": R_MEASURED_MIN, "range_hi": R_MEASURED_MAX,
          "fail": bool(worst["separation"] < SEPARATION_FLOOR)}

    # C3 — descriptive only, declared as such before running.
    c3: Dict = {}
    for n in PAIR:
        vals = [by[n][R]["E_p90_over_p10"] for R in R_LADDER
                if R_MEASURED_MIN <= R <= R_MEASURED_MAX]
        c3[n] = {"min": float(min(vals)), "max": float(max(vals)),
                 "fold_span": float(max(vals) / min(vals)),
                 "shipped": NUISANCE[n]["shipped"],
                 "nuisance_lo": NUISANCE[n]["lo"],
                 "nuisance_hi": NUISANCE[n]["hi"]}

    if c1["fail"] and c2["fail"]:
        verdict = ("THE HEADLINE IS CAPACITY-SENSITIVE AND SO IS THE CLAIM "
                   "BUILT ON IT")
    elif c1["fail"]:
        verdict = ("THE HEADLINE'S VALUE IS CAPACITY-SENSITIVE; THE SEPARATION "
                   "IS NOT")
    elif c2["fail"]:
        verdict = "THE SEPARATION IS CAPACITY-SENSITIVE"
    else:
        verdict = "THE HEADLINE SURVIVES CAPACITY"
    return {"per_R": per_R, "C1": c1, "C2": c2, "C3": c3, "verdict": verdict}


def check_reproduction(rows: List[Dict]) -> Dict:
    """R = 1 on Arm A IS the shipping config and must reproduce
    delivery_spread_noise.md sec. 2's flagship row. A check, not a threshold —
    and nothing is adjusted if it fails."""
    out: Dict = {}
    for n in PAIR:
        got = next(r["E_p90_over_p10"] for r in rows
                   if r["drug"] == n and r["R"] == 1.0)
        pub = NUISANCE[n]["shipped"]
        out[n] = {"measured": float(got), "published": pub,
                  "abs_diff": float(abs(got - pub)),
                  "rounds_to_published": bool(round(got, 3) == pub)}
        print(f"  reproduction {n:<14} measured {got:.4f} vs published "
              f"{pub} -> {'MATCHES' if out[n]['rounds_to_published'] else 'DIFFERS'}",
              flush=True)
    return out


# --------------------------------------------------------------------------- #
_OUT: Dict = {}


def _dump(extra: Optional[Dict] = None) -> None:
    if extra:
        _OUT.update(extra)
    with open(_DATA, "w") as fh:
        json.dump(_OUT, fh, indent=1)


def _load_existing() -> None:
    """Stages ACCUMULATE rather than clobber, so a partial file never looks
    like a complete record of a run that was never made."""
    if os.path.exists(_DATA):
        try:
            with open(_DATA) as fh:
                _OUT.update(json.load(fh))
            print(f"  loaded existing {os.path.basename(_DATA)}: {sorted(_OUT)}",
                  flush=True)
        except (OSError, ValueError) as exc:               # pragma: no cover
            print(f"  existing {_DATA} unreadable ({exc}); starting fresh")


def main() -> Dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "guards", "plant", "identity",
                             "truncated", "emulation", "analyse"])
    args = ap.parse_args()

    alpha = _assert_alpha(f"{args.stage}/top")
    _assert_frozen_pde()
    print(f"calibrate_alpha() = {alpha!r}  (asserted, top of stage)")
    print(f"PDEConfig() defaults asserted; PKConfig() = {PKConfig()}")
    print(f"OxygenConfig().hypoxia_threshold = "
          f"{OxygenConfig().hypoxia_threshold!r}")
    print(f"config: dx = {DX_UM:g} um, repair = {REPAIR}, dt = {DT_S:g} s, "
          f"window = {WINDOW_S:g} s, convection = {CONVECTION}, seed = {SEED}")
    if args.stage == "guards":
        _assert_alpha("guards/bottom")
        return {}

    t0 = time.time()
    _load_existing()
    _OUT.update({"alpha": alpha,
                 "config": {"dx_um": DX_UM, "repair": REPAIR, "dt_s": DT_S,
                            "window_s": WINDOW_S, "convection": CONVECTION,
                            "seed": SEED, "drugs": list(PAIR),
                            "section": "xenium GBM, build_arm_microenvironment",
                            "pk": dataclasses.asdict(PKConfig()),
                            "hypoxia_threshold":
                                OxygenConfig().hypoxia_threshold},
                 "r_ladder": [float(v) for v in R_LADDER],
                 "inherited_thresholds": {
                     "nuisance": NUISANCE,
                     "separation_floor": SEPARATION_FLOOR,
                     "source": "results/delivery_spread_noise.md secs. 0 and 2"},
                 })
    _dump()

    if args.stage in ("all", "plant"):
        print("\n--- PLANT: every guard, against the defect it claims to catch ---")
        _dump({"plant": stage_plant()})
        _assert_alpha("plant/bottom")

    if args.stage in ("all", "identity"):
        print("\n--- IDENTITY: Arm A validated, and the trap measured ---")
        _dump({"identity": stage_identity(alpha)})
        _assert_alpha("identity/bottom")

    if args.stage in ("all", "truncated"):
        print("\n--- ARM A: truncated window (dt/R, T/R), convection ON "
              "— the CORRECT capacity change ---")
        s = stage_sweep(alpha, "A")
        s["reproduction"] = check_reproduction(s["rows"])
        s["analysis"] = analyse(s["rows"])
        _dump({"arm_A": s})
        _assert_alpha("truncated/bottom")
        _report(s)

    if args.stage in ("all", "emulation"):
        print("\n--- ARM B: parameter division at (dt, T), convection ON "
              "— a capacity change PLUS an R-fold Peclet boost, NOT capacity ---")
        s = stage_sweep(alpha, "B")
        s["analysis"] = analyse(s["rows"])
        _dump({"arm_B": s})
        _assert_alpha("emulation/bottom")
        _report(s)

    if args.stage == "analyse":
        for key in ("arm_A", "arm_B"):
            if key in _OUT and _OUT[key].get("rows"):
                _OUT[key]["analysis"] = analyse(_OUT[key]["rows"])
                _report(_OUT[key])
        _dump()

    _OUT["elapsed_s"] = _OUT.get("elapsed_s", 0.0) + (time.time() - t0)
    _OUT.pop("partial_arm_A", None)
    _OUT.pop("partial_arm_B", None)
    _dump()
    _assert_alpha(f"{args.stage}/bottom")
    print(f"\ncalibrate_alpha() = {alpha!r}  (asserted, bottom of stage)")
    print(f"wrote {_DATA}  ({time.time() - t0:.0f}s this invocation)")
    return _OUT


def _report(s: Dict) -> None:
    a = s["analysis"]
    arm = ("ARM A (truncated window, convection ON — genuine capacity)"
           if s["arm"] == "A" else
           "ARM B (parameter division, convection ON — capacity PLUS an "
           "R-fold Peclet boost, NOT capacity)")
    print("\n" + "=" * 78)
    print(f"  {arm}")
    print("=" * 78)
    print(f"  {'R':>8} {PAIR[0]:>16} {PAIR[1]:>16} {'sep':>8}")
    for p in a["per_R"]:
        print(f"  {p['R']:>8g} {p[PAIR[0]]:>16.4f} {p[PAIR[1]]:>16.4f} "
              f"{p['separation']:>8.3f}")
    for n in PAIR:
        d = a["C1"]["per_drug"][n]
        print(f"  C1 {n:<14} R1 {d['R1']:.4f} -> R39 {d['R39']:.4f}, move "
              f"{d['abs_move']:.4f} = {d['times_nuisance_width']:.1f}x its "
              f"nuisance width {d['nuisance_width']:.3f} -> "
              f"{'FAIL' if d['fail'] else 'pass'}")
    c2 = a["C2"]
    print(f"  C2 min separation {c2['min_separation']:.3f} at R = {c2['at_R']:g} "
          f"(range {c2['range_lo']}-{c2['range_hi']}), floor "
          f"{c2['floor']} -> {'FAIL' if c2['fail'] else 'pass'}")
    print(f"\n  VERDICT [{s['arm']}]: {a['verdict']}")


if __name__ == "__main__":
    main()
