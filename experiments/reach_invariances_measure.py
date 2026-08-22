"""THE TWO REACH-MAP INVARIANCE MEASUREMENTS THAT ARE SHIPPED BUT UNRECORDED.

WHY THIS EXISTS
---------------
Two measurements are quoted on the Studio reach panel
(``app/web/js/ui/reach.js``) and in the payload note
``app/server/reach.SCALE_INVARIANCE_NOTE``:

1. **The reach map is exactly invariant to the transvascular barrier term.**
   ``bbb_penetrance`` is a pure global multiplier on the delivery source
   (``treatment.drugs.delivery_source``: ``S = vascular_density *
   bbb_penetrance * plasma``), so it divides out of every within-map ratio. The
   shipped words are *"verified bit-identical at bbb_penetrance 0.1 and 1.0"*.
2. **The falloff SHAPE is not capacity-invariant while the under-reached ladder
   nearly is.** The shipped words are *"At R = 39 this curve moves 0.3772
   absolute and the fitted penetration length collapses from 279.7 to 76.7
   um"*, and *"over R = 1 to 39 the whole ladder moves by at most 0.0011, and
   it is exactly flat out to R = 60 ... across the full measured range of 3.2
   to 4184 mL/g it moves by up to 0.80"*.

Neither had a ``results/`` record, so by this repo's convention nothing else may
cite them. This script re-measures both from scratch. It transcribes no number:
every figure in ``results/reach_invariances_measure.md`` comes out of
``results/reach_invariances_measure.json``, which comes out of this file.

WHAT IS CHANGED
---------------
Nothing. ``calibrate_alpha()`` is asserted live at ``22.223477711172563`` and
the four ``PDEConfig`` defaults are asserted at every stage. No default,
constant or shipped parameter is modified; ``tmesim/`` and ``app/`` are only
driven, never written to. The two ``bbb_penetrance`` values in stage ``bbb``
are per-run ``dataclasses.replace`` copies of a compiled ``DrugParams`` — the
library value is untouched.

THE CAPACITY IDENTITY, AND WHY IT NEEDS NO CAPACITY TERM
--------------------------------------------------------
A tissue-capacity (retardation) term ``R`` makes the transport equation

    R dC/dt = div(D grad C) - div(v C) + S - k C

i.e. ``dC/dt = (1/R)[...]``, a pure rescale of TIME: ``C_R(x, t) =
C_1(x, t/R)`` and ``AUC_R(x, T) = R * AUC_1(x, T/R)``.

The DISCRETE statement is stronger than the continuous one, and it is the one
used here. The solver's step is
``c_new = A^-1 [c + h S - h adv(c)]`` with ``A = I + h (L + k)``. Running it
with the capacity-R coefficients ``(L/R, k/R, S/R)`` at timestep ``dt`` gives
``A = I + (dt/R)(L + k)`` and ``rhs = c + (dt/R) S - (dt/R) adv(c)`` — which is
**the R = 1 solver at timestep ``h = dt/R``**, iterate for iterate. So:

    the capacity-R reach map at window T, timestep dt
      ==  the R = 1 reach map at window T/R, timestep dt/R

with the AUC field differing by the global factor R, which every reach quantity
(``falloff = bin mean / E_ref``, ``r_i = AUC_i / E_ref``, the fitted decay
length ``L``) divides out exactly. The step COUNT is unchanged, so no
resolution is traded away for the truncation.

**This route is not an approximation of the parameter-division route — it is
the other way round.** The velocity field is built inside
``sim.engine.simulate_response`` and cannot be divided by R from outside
(``results/capacity_ranking_robustness.md`` sec. 2d says so and runs
convection OFF for that reason). The ``dt/R`` route retards the advection term
correctly because ``h`` multiplies it too. Stage ``identity`` measures both:
the two routes agree to solver precision with convection OFF, and the size of
the term the parameter-division route gets wrong with convection ON is
reported rather than assumed small.

PRECONDITIONS, VERIFIED RATHER THAN ASSUMED (stage ``identity``)
----------------------------------------------------------------
* **The non-negativity clamp never fires.** ``physics/diffusion.py`` clamps
  ``np.maximum(c_new, 0.0, out=c_new)`` every step. A plasma-rescale test
  CANNOT establish this — ``max(a*x, 0) = a*max(x, 0)`` for ``a > 0``, so the
  clamp is positively homogeneous and an exact rescale result survives it
  intact. It is instrumented directly instead: a counting shim on the module's
  ``np`` records, for every step of every solve in this stage, how many entries
  were negative before the clamp and the worst undershoot.
* **The source is static.** ``PDEConfig().time_resolved`` is False, so the
  plasma level enters as one scalar ``plasma_mean``, asserted here.
* **``plasma_mean`` is window-dependent and that is harmless — measured, not
  argued.** ``simulate_response`` averages the decaying plasma profile over
  ``arange(0, total_time, dt)``, so a truncated window changes the source by a
  global factor. That factor is exactly what measurement 1 shows the reach map
  divides out, and stage ``identity`` re-measures it on the truncated windows
  rather than inheriting it.
* **The viable band does not move.** Oxygen is solved independently of
  ``PDEConfig``, so ``n_viable`` must be identical at every R; asserted.
* **The usable distance bins do not move.** They are set by per-bin grid-point
  counts, i.e. pure geometry; asserted identical at every R.

Run:
    ./.venv/bin/python experiments/reach_invariances_measure.py --stage guards
    ./.venv/bin/python experiments/reach_invariances_measure.py --stage all
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
from tmesim.reconstruct.real_pipeline import reconstruct_real
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.validation.delivery_metrics import penetration_length_fit
from tmesim.validation.vessel_profile import (
    direct_vessel_distance_um,
    endothelial_xy,
)

# --------------------------------------------------------------------------- #
# THE REACH MODULE IS PINNED TO A COMMIT, ON PURPOSE
# --------------------------------------------------------------------------- #
# A sibling agent is converting the reach map from a FINITE-WINDOW object to a
# STEADY-STATE one and owns ``app/server/reach.py``. Mid-run its working-tree
# copy raised ``NameError: name 'FINITE_WINDOW_NOTE' is not defined``, which is
# a partially-applied edit, not a bug in this study — but it is also the reason
# a measurement of the CURRENT disclosures must not read a file that is being
# rewritten underneath it.
#
# So the module is loaded from the git blob at the commit whose disclosures this
# study is documenting, verified by SHA-256, and the version identifier is
# written into the output. Every number here describes THAT reach map. If the
# steady-state rewrite lands, these finite-window numbers describe the previous
# configuration and this pin is how a reader can tell.
PINNED_REACH_COMMIT = "23e51f6"
PINNED_REACH_PATH = "app/server/reach.py"
PINNED_REACH_SHA256 = (
    "aaf5d5893b697780fd88af5de555dd11a78b4a241f8aeca07cf03fca7bbe7a66")


def _load_pinned_reach():
    """``app/server/reach.py`` as of :data:`PINNED_REACH_COMMIT`, hash-checked.

    Executed as an anonymous module rather than imported, so nothing in the
    working tree is read, written or shadowed.
    """
    import hashlib
    import subprocess
    import types

    src = subprocess.check_output(
        ["git", "-C", _ROOT, "show",
         f"{PINNED_REACH_COMMIT}:{PINNED_REACH_PATH}"])
    got = hashlib.sha256(src).hexdigest()
    if got != PINNED_REACH_SHA256:
        raise AssertionError(
            f"{PINNED_REACH_COMMIT}:{PINNED_REACH_PATH} hashes to {got}, "
            f"expected {PINNED_REACH_SHA256}")
    name = f"reach_pinned_{PINNED_REACH_COMMIT}"
    mod = types.ModuleType(name)
    mod.__file__ = f"<git:{PINNED_REACH_COMMIT}:{PINNED_REACH_PATH}>"
    exec(compile(src, mod.__file__, "exec"), mod.__dict__)
    return mod


def _worktree_reach_sha256() -> Optional[str]:
    """What the working tree's copy hashes to right now — recorded, not used."""
    import hashlib
    p = os.path.join(_ROOT, PINNED_REACH_PATH)
    try:
        with open(p, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError:                                        # pragma: no cover
        return None


_REACH = _load_pinned_reach()
compute_reach = _REACH.compute_reach
REACH_LADDER = _REACH.REACH_LADDER
HEADLINE_X = _REACH.HEADLINE_X

# --------------------------------------------------------------------------- #
# Frozen state this study is only valid against. Asserted, never adjusted.
# --------------------------------------------------------------------------- #
SHIPPED_ALPHA = 22.223477711172563
FROZEN_PDE_DEFAULTS = dict(dt=0.5, total_time=43200.0, convection=True,
                           save_every=200)

#: The configuration measurement 2 is quoted at. dx = 15 um is the
#: ``results/ranking_robustness.md`` / ``capacity_ranking_robustness.md``
#: Monte-Carlo grid; dt = 10 s, 6 h and convection ON are the SHIPPING studio
#: values (``app/server/jobs.py`` PDE_DT_S / PDE_CONVECTION,
#: ``app/server/api.py`` STUDIO_DEFAULT_TOTAL_TIME_S).
DX_UM = 15.0
DT_S = 10.0
WINDOW_S = 3600.0 * 6
CONVECTION = True
SAVE_EVERY = 10_000            # history is not read here; keep it out of memory
DRUG = "temozolomide"

#: Grids the barrier measurement is repeated on. 15 um is measurement 2's grid,
#: 12 um is the studio default (``app/server/api.default_spec``), 30 um is the
#: coarse arm ``results/delivery_spread_noise.md`` also runs. Three resolutions
#: because an exactness claim that holds at one grid and is only asserted at the
#: others is not the claim the panel makes. It also puts ``n_viable`` on the
#: record at three grids, which the panel's own grid-sensitivity remark needs.
STUDIO_DX_UM = 12.0
COARSE_DX_UM = 30.0
BBB_DX: Sequence[float] = (15.0, 12.0, 30.0)

#: bbb_penetrance values for measurement 1. 0.1 and 1.0 are the two the panel
#: names; the compiled value is included so the shipping point is on the record.
BBB_VALUES: Sequence[float] = (0.1, 1.0)

#: Capacity grid. 1 is the shipping value; 3.2 / 39 / 4184 are the minimum,
#: best-available constant and maximum of the 33 measured rat brain-slice
#: Vu,brain values (results/brain_slice_anchor.md). 50-120 is deliberately dense
#: because the shipped claim locates a "flat to R = 60" boundary there.
R_GRID: Sequence[float] = (1.0, 2.0, 3.2, 5.0, 8.0, 13.0, 20.0, 30.0, 39.0,
                           50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 270.0, 620.0,
                           1420.0, 4184.0)

#: The convection-OFF arm's grid. Deliberately a SUBSET of :data:`R_GRID`, not a
#: different one: it keeps R = 1 and R = 39 (the claim's two endpoints), the
#: 60/80 region where the shipping arm's ladder stops being exactly flat, and
#: the ends of the measured capacity range. A cross-arm comparison has to be
#: made at shared R values or it is not a comparison.
R_GRID_OFF: Sequence[float] = (1.0, 3.2, 20.0, 39.0, 60.0, 80.0, 120.0, 620.0,
                               4184.0)

#: Capacities at which the truncated-window route is checked against a genuine
#: parameter-division retarded run.
R_IDENTITY: Sequence[float] = (4.0, 39.0)

_RESULTS = os.path.join(_ROOT, "results")
_DATA = os.path.join(_RESULTS, "reach_invariances_measure.json")


# --------------------------------------------------------------------------- #
# Guards
# --------------------------------------------------------------------------- #
def assert_frozen_state() -> float:
    """Every constant this study leans on, checked against LIVE code."""
    a = calibrate_alpha()
    if a != SHIPPED_ALPHA:
        raise AssertionError(
            f"calibrate_alpha() = {a!r}, expected {SHIPPED_ALPHA!r}")
    d = PDEConfig()
    for k, v in FROZEN_PDE_DEFAULTS.items():
        if getattr(d, k) != v:
            raise AssertionError(f"PDEConfig().{k} = {getattr(d, k)!r}, "
                                 f"expected {v!r}")
    if getattr(d, "time_resolved") is not False:
        raise AssertionError("PDEConfig().time_resolved is not False: the "
                             "source would no longer be a single scalar and "
                             "the truncated-window identity would not hold")
    if getattr(d, "transvascular_exchange") is not False:
        raise AssertionError("PDEConfig().transvascular_exchange is not False")
    if tuple(REACH_LADDER) != (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
        raise AssertionError(f"REACH_LADDER moved: {tuple(REACH_LADDER)!r}")
    if float(HEADLINE_X) != 0.1:
        raise AssertionError(f"HEADLINE_X moved: {HEADLINE_X!r}")
    return a


def compile_flagship(name: str, alpha: float):
    """The pipeline the shipped real-tissue panel uses."""
    return apply_barrier_to_drugparams(calibrated_drugparams(name, alpha),
                                       get_molecule(name))


# --------------------------------------------------------------------------- #
# The non-negativity clamp watcher
# --------------------------------------------------------------------------- #
class ClampWatch:
    """Count every firing of ``np.maximum(c_new, 0.0, out=c_new)``.

    The clamp is the ONE non-linearity in the transport solve, and it is exactly
    the one a rescale test cannot see (it is positively homogeneous). This
    swaps the diffusion module's ``np`` for a forwarding shim that inspects the
    array before the clamp is applied, and restores it unconditionally.
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
            def __getattr__(self, k):          # everything else is numpy's
                return getattr(real, k)

            def maximum(self, a, b, out=None, **kw):
                # The clamp is the only in-place maximum against scalar 0.0.
                # ``_effective_diffusivity``'s maximum has out=None and a
                # non-zero floor, so it is not counted.
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
# One reach map
# --------------------------------------------------------------------------- #
_MICROS: Dict[float, object] = {}
_DIST: Dict[float, np.ndarray] = {}


def micro_for(dx_um: float):
    if dx_um not in _MICROS:
        t0 = time.time()
        m = reconstruct_real(source="xenium", dx_um=dx_um)
        _MICROS[dx_um] = m
        _DIST[dx_um] = direct_vessel_distance_um(m.domain,
                                                 endothelial_xy(m.cells))
        print(f"  reconstructed dx = {dx_um:g} um, grid {m.domain.shape}, "
              f"{m.cells.n} cells ({time.time() - t0:.0f}s)")
    return _MICROS[dx_um]


def reach_map(dx_um: float, dp, dt: float, total_time: float,
              convection: bool = CONVECTION,
              watch: Optional[ClampWatch] = None) -> Dict:
    """One reach map, through the SAME calls ``app/server/jobs.py`` makes."""
    micro = micro_for(dx_um)
    dist = _DIST[dx_um]
    pde = PDEConfig(dt=dt, total_time=total_time, convection=convection,
                    save_every=SAVE_EVERY)
    t0 = time.time()
    if watch is not None:
        with watch:
            res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
    else:
        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
    pen = penetration_length_fit(res.fields["auc"], micro, distance_um=dist,
                                 include_source_null=False)
    reach = compute_reach(micro, res, distance_um=dist, penetration=pen)
    if reach["null_reason"]:
        raise AssertionError(f"reach refused: {reach['null_reason']}")
    return {
        "dx_um": float(dx_um), "dt_s": float(dt),
        "total_time_s": float(total_time), "convection": bool(convection),
        "bbb_penetrance": float(dp.bbb_penetrance),
        "diffusivity": float(dp.diffusivity), "binding": float(dp.binding),
        "clearance": float(dp.clearance),
        "n_steps": int(round(total_time / dt)),
        "E_ref_uM_s": float(reach["E_ref_uM_s"]),
        "headline": float(reach["headline_under_reached_fraction"]),
        "ladder_x": [float(v) for v in reach["curve"]["x"]],
        "F": [float(v) for v in reach["curve"]["under_reached_fraction"]],
        "falloff_distance_um": [float(v) for v in
                                reach["falloff"]["distance_um"]],
        "falloff": [float(v) for v in reach["falloff"]["normalised_mean_auc"]],
        "falloff_count": [int(v) for v in reach["falloff"]["count"]],
        "L_um": float(pen["drug"]["L_um"]),
        "near_exp_L_um": float(pen["drug"]["near_exp_L_um"]),
        "n_viable": int(reach["n_viable"]), "n_tumor": int(reach["n_tumor"]),
        "elapsed_s": time.time() - t0,
    }


def _bitwise(a: Sequence[float], b: Sequence[float]) -> Dict:
    """Bit-identical or not, and by how much if not."""
    x = np.asarray(a, dtype=np.float64)
    y = np.asarray(b, dtype=np.float64)
    if x.shape != y.shape:
        return {"identical": False, "max_abs": float("nan"),
                "max_rel": float("nan"), "n_differing": -1,
                "note": f"shape {x.shape} vs {y.shape}"}
    d = np.abs(x - y)
    denom = np.where(np.abs(y) > 0, np.abs(y), 1.0)
    return {"identical": bool(np.array_equal(x, y)),
            "max_abs": float(d.max()) if d.size else 0.0,
            "max_rel": float((d / denom).max()) if d.size else 0.0,
            "n_differing": int(np.count_nonzero(x != y))}


# --------------------------------------------------------------------------- #
# Stage 1 — the barrier invariance
# --------------------------------------------------------------------------- #
def stage_bbb(alpha: float, dx_list: Sequence[float]) -> Dict:
    base = compile_flagship(DRUG, alpha)
    out: Dict = {"drug": DRUG, "compiled_bbb_penetrance":
                 float(base.bbb_penetrance), "arms": [], "comparisons": []}
    for dx in dx_list:
        arms = []
        for bbb in BBB_VALUES:
            dp = dataclasses.replace(base, bbb_penetrance=float(bbb))
            r = reach_map(dx, dp, DT_S, WINDOW_S)
            r["arm"] = f"dx{dx:g}_bbb{bbb:g}"
            arms.append(r)
            print(f"  bbb {bbb:<4g} dx {dx:g}: F(0.1) = {r['headline']:.6f}  "
                  f"E_ref = {r['E_ref_uM_s']:.6g}  L = {r['L_um']:.4f} um  "
                  f"({r['elapsed_s']:.0f}s)")
        lo, hi = arms[0], arms[1]
        cmp_ = {
            "dx_um": float(dx),
            "bbb_lo": float(BBB_VALUES[0]), "bbb_hi": float(BBB_VALUES[1]),
            "F": _bitwise(lo["F"], hi["F"]),
            "falloff": _bitwise(lo["falloff"], hi["falloff"]),
            "headline_identical": bool(lo["headline"] == hi["headline"]),
            "headline_lo": lo["headline"], "headline_hi": hi["headline"],
            "L_um_lo": lo["L_um"], "L_um_hi": hi["L_um"],
            "L_um_identical": bool(lo["L_um"] == hi["L_um"]),
            "L_um_rel_diff": abs(lo["L_um"] - hi["L_um"]) / hi["L_um"],
            "E_ref_lo": lo["E_ref_uM_s"], "E_ref_hi": hi["E_ref_uM_s"],
            "E_ref_ratio": lo["E_ref_uM_s"] / hi["E_ref_uM_s"],
            "E_ref_ratio_expected": float(BBB_VALUES[0]) / float(BBB_VALUES[1]),
        }
        cmp_["E_ref_ratio_rel_err"] = abs(
            cmp_["E_ref_ratio"] - cmp_["E_ref_ratio_expected"]
        ) / cmp_["E_ref_ratio_expected"]
        out["arms"].extend(arms)
        out["comparisons"].append(cmp_)
        print(f"  -> dx {dx:g}: ladder identical = {cmp_['F']['identical']}, "
              f"falloff identical = {cmp_['falloff']['identical']} "
              f"(max abs {cmp_['falloff']['max_abs']:.3e}), "
              f"E_ref ratio {cmp_['E_ref_ratio']:.17g}")
    return out


# --------------------------------------------------------------------------- #
# Stage 2 — preconditions for the capacity identity
# --------------------------------------------------------------------------- #
def stage_identity(alpha: float) -> Dict:
    base = compile_flagship(DRUG, alpha)
    out: Dict = {"clamp": [], "routes": [], "convection_term": []}

    # (a) the clamp, at the shipping configuration and at the most extreme
    #     truncated window the capacity sweep uses.
    for label, dt, T in (("shipping R=1 (dt 10 s, 6 h)", DT_S, WINDOW_S),
                         (f"truncated R={max(R_GRID):g}",
                          DT_S / max(R_GRID), WINDOW_S / max(R_GRID))):
        w = ClampWatch()
        r = reach_map(DX_UM, base, dt, T, watch=w)
        rec = {"label": label, "dt_s": dt, "total_time_s": T, **w.as_dict(),
               "headline": r["headline"], "E_ref_uM_s": r["E_ref_uM_s"]}
        out["clamp"].append(rec)
        print(f"  clamp [{label}]: {rec['clamp_calls']} calls, "
              f"{rec['negative_entries']} negative entries, worst "
              f"{rec['worst_undershoot']:.3e} -> "
              f"{'FIRED' if rec['clamp_ever_fired'] else 'never fired'}")

    # (b) truncated-window route vs genuine parameter-division retarded run.
    #     Convection OFF is where the parameter-division route is exact
    #     (the velocity field cannot be divided by R from outside).
    for R in R_IDENTITY:
        trunc = reach_map(DX_UM, base, DT_S / R, WINDOW_S / R,
                          convection=False)
        ret = dataclasses.replace(
            base,
            diffusivity=base.diffusivity / R,
            binding=base.binding / R,
            clearance=base.clearance / R,
            bbb_penetrance=base.bbb_penetrance / R,
        )
        param = reach_map(DX_UM, ret, DT_S, WINDOW_S, convection=False)
        rec = {"R": float(R), "convection": False,
               "F": _bitwise(trunc["F"], param["F"]),
               "falloff": _bitwise(trunc["falloff"], param["falloff"]),
               "L_um_truncated": trunc["L_um"], "L_um_param": param["L_um"],
               "headline_truncated": trunc["headline"],
               "headline_param": param["headline"]}
        rec["L_um_rel_diff"] = abs(
            rec["L_um_truncated"] - rec["L_um_param"]) / rec["L_um_param"]
        out["routes"].append(rec)
        print(f"  identity R={R:g} (convection OFF): falloff max abs "
              f"{rec['falloff']['max_abs']:.3e}, ladder max abs "
              f"{rec['F']['max_abs']:.3e}, L {rec['L_um_truncated']:.4f} vs "
              f"{rec['L_um_param']:.4f} um")

        # (c) how big is the term the parameter-division route gets wrong when
        #     convection is ON — i.e. why the truncated route is the correct one.
        trunc_on = reach_map(DX_UM, base, DT_S / R, WINDOW_S / R,
                             convection=True)
        param_on = reach_map(DX_UM, ret, DT_S, WINDOW_S, convection=True)
        rec2 = {"R": float(R), "convection": True,
                "F": _bitwise(trunc_on["F"], param_on["F"]),
                "falloff": _bitwise(trunc_on["falloff"], param_on["falloff"]),
                "L_um_truncated": trunc_on["L_um"],
                "L_um_param": param_on["L_um"]}
        out["convection_term"].append(rec2)
        print(f"  convection term R={R:g}: falloff max abs "
              f"{rec2['falloff']['max_abs']:.3e}, L "
              f"{rec2['L_um_truncated']:.4f} vs {rec2['L_um_param']:.4f} um")
    return out


# --------------------------------------------------------------------------- #
# Stage 3 — the capacity sweep
# --------------------------------------------------------------------------- #
def stage_capacity(alpha: float, convection: bool = CONVECTION,
                   r_grid: Sequence[float] = R_GRID) -> Dict:
    """The capacity sweep, by the truncated-window route.

    ``convection`` selects the ARM, and the arm has to be stated with every
    number that comes out of it:

    * ``True`` — the shipping configuration. The truncated-window route retards
      the advection term correctly (``h = dt/R`` multiplies it), so this arm is
      a genuine capacity change with ``v -> v/R``.
    * ``False`` — the protocol ``results/capacity_ranking_robustness.md`` sec.
      2d runs, and the protocol under which a parameter-division emulation
      (divide ``D``, ``binding``, ``clearance``, ``bbb_penetrance`` by ``R``)
      is exact. Run so that a claim made under THAT protocol can be checked
      against a like-for-like arm rather than against the shipping one.
    """
    base = compile_flagship(DRUG, alpha)
    rows: List[Dict] = []
    key = "partial_capacity" + ("" if convection else "_off")
    for R in r_grid:
        r = reach_map(DX_UM, base, DT_S / R, WINDOW_S / R,
                      convection=convection)
        r["R"] = float(R)
        rows.append(r)
        print(f"  R = {R:<8g} F(0.1) = {r['headline']:.6f}  "
              f"L = {r['L_um']:8.3f} um  E_ref = {r['E_ref_uM_s']:.6g}  "
              f"n_viable = {r['n_viable']}  ({r['elapsed_s']:.0f}s)")
        _dump({key: rows})
    return {"rows": rows, "convection": bool(convection),
            "r_grid": [float(v) for v in r_grid],
            "route": "truncated window (dt/R, T/R); no capacity term is "
                     "implemented and no parameter is divided"}


def analyse_capacity(rows: List[Dict]) -> Dict:
    base = next(r for r in rows if r["R"] == 1.0)
    f0 = np.asarray(base["falloff"])
    F0 = np.asarray(base["F"])
    d0 = base["falloff_distance_um"]

    # Preconditions that make the comparison elementwise-meaningful.
    for r in rows:
        if r["n_viable"] != base["n_viable"] or r["n_tumor"] != base["n_tumor"]:
            raise AssertionError(
                f"viable band moved at R = {r['R']}: {r['n_viable']} vs "
                f"{base['n_viable']} — oxygen is supposed to be independent of "
                "PDEConfig and the ladder would not be comparable")
        if r["falloff_distance_um"] != d0:
            raise AssertionError(f"usable distance bins moved at R = {r['R']}")
        if r["n_steps"] != base["n_steps"]:
            raise AssertionError(f"step count moved at R = {r['R']}: "
                                 f"{r['n_steps']} vs {base['n_steps']}")

    per_r = []
    for r in rows:
        f = np.asarray(r["falloff"])
        F = np.asarray(r["F"])
        per_r.append({
            "R": r["R"],
            "falloff_max_abs_move": float(np.abs(f - f0).max()),
            "falloff_argmax_distance_um": float(
                d0[int(np.argmax(np.abs(f - f0)))]),
            "falloff_min": float(f.min()), "falloff_max": float(f.max()),
            "ladder_max_abs_move": float(np.abs(F - F0).max()),
            "ladder_exactly_flat": bool(np.array_equal(F, F0)),
            "headline": r["headline"],
            "headline_move": float(r["headline"] - base["headline"]),
            "L_um": r["L_um"], "near_exp_L_um": r["near_exp_L_um"],
            "E_ref_uM_s": r["E_ref_uM_s"],
        })

    def _upto(limit: float, key: str) -> float:
        return max(p[key] for p in per_r if p["R"] <= limit)

    flat = [p["R"] for p in per_r if p["ladder_exactly_flat"]]
    first_move = [p["R"] for p in per_r if not p["ladder_exactly_flat"]]
    r39 = next(p for p in per_r if p["R"] == 39.0)
    return {
        "per_R": per_r,
        "falloff_R1_range": [float(f0.min()), float(f0.max())],
        "ladder_max_move_R1_to_39": _upto(39.0, "ladder_max_abs_move"),
        "ladder_max_move_full_range": max(p["ladder_max_abs_move"]
                                          for p in per_r),
        "ladder_exactly_flat_upto_R": max(flat) if flat else None,
        "ladder_first_R_that_moves": min(first_move) if first_move else None,
        "falloff_move_at_R39": r39["falloff_max_abs_move"],
        "falloff_move_at_R39_distance_um": r39["falloff_argmax_distance_um"],
        "L_um_R1": base["L_um"], "L_um_R39": r39["L_um"],
        "headline_R1": base["headline"],
        "n_viable": base["n_viable"], "n_tumor": base["n_tumor"],
        "n_steps": base["n_steps"],
    }


# --------------------------------------------------------------------------- #
_OUT: Dict = {}


def _dump(extra: Optional[Dict] = None) -> None:
    if extra:
        _OUT.update(extra)
    with open(_DATA, "w") as fh:
        json.dump(_OUT, fh, indent=1)


def _load_existing() -> None:
    """Stages ACCUMULATE into one record rather than clobbering it.

    Without this, running a stage on its own would silently drop every stage
    already measured, and the file would look like a complete record of a run
    that was never made.
    """
    if os.path.exists(_DATA):
        try:
            with open(_DATA) as fh:
                _OUT.update(json.load(fh))
            print(f"  loaded existing {os.path.basename(_DATA)}: "
                  f"{sorted(_OUT)}")
        except (OSError, ValueError) as exc:            # pragma: no cover
            print(f"  existing {_DATA} unreadable ({exc}); starting fresh")


def main() -> Dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "guards", "bbb", "identity", "capacity",
                             "capacity_off"])
    args = ap.parse_args()

    alpha = assert_frozen_state()
    print(f"calibrate_alpha() = {alpha!r}  (asserted)")
    print(f"PDEConfig() defaults asserted; PKConfig() = {PKConfig()}")
    print(f"OxygenConfig().hypoxia_threshold = "
          f"{OxygenConfig().hypoxia_threshold!r}")
    if args.stage == "guards":
        return {}

    t0 = time.time()
    _load_existing()
    _OUT.update({
        "alpha": alpha,
        "reach_module": {
            "pinned_commit": PINNED_REACH_COMMIT,
            "path": PINNED_REACH_PATH,
            "sha256": PINNED_REACH_SHA256,
            "worktree_sha256": _worktree_reach_sha256(),
            "worktree_matches_pin": (_worktree_reach_sha256()
                                     == PINNED_REACH_SHA256),
        },
        "config": {"dx_um": DX_UM, "dt_s": DT_S, "window_s": WINDOW_S,
                   "convection": CONVECTION, "drug": DRUG,
                   "section": "xenium GBM (reconstruct_real source='xenium')",
                   "studio_dx_um": STUDIO_DX_UM,
                   "pk": dataclasses.asdict(PKConfig()),
                   "hypoxia_threshold": OxygenConfig().hypoxia_threshold},
        "ladder": [float(v) for v in REACH_LADDER],
        "headline_x": float(HEADLINE_X),
        "r_grid": [float(v) for v in R_GRID],
    })
    _dump()

    if args.stage in ("all", "bbb"):
        print("\n--- MEASUREMENT 1: barrier (bbb_penetrance) invariance ---")
        _dump({"bbb": stage_bbb(alpha, list(BBB_DX))})

    if args.stage in ("all", "identity"):
        print("\n--- PRECONDITIONS: clamp, and the capacity identity ---")
        _dump({"identity": stage_identity(alpha)})

    if args.stage in ("all", "capacity"):
        print("\n--- MEASUREMENT 2: capacity sweep via truncated windows ---")
        cap = stage_capacity(alpha)
        cap["analysis"] = analyse_capacity(cap["rows"])
        _dump({"capacity": cap})
        a = cap["analysis"]
        print("\n" + "=" * 78)
        print(f"  falloff at R = 1 spans {a['falloff_R1_range'][1]:.4f} -> "
              f"{a['falloff_R1_range'][0]:.4f}")
        print(f"  falloff max absolute move at R = 39: "
              f"{a['falloff_move_at_R39']:.4f} "
              f"(at {a['falloff_move_at_R39_distance_um']:.0f} um)")
        print(f"  fitted L: {a['L_um_R1']:.1f} um (R = 1) -> "
              f"{a['L_um_R39']:.1f} um (R = 39)")
        print(f"  ladder max move R = 1..39: "
              f"{a['ladder_max_move_R1_to_39']:.6f}")
        print(f"  ladder exactly flat up to R = "
              f"{a['ladder_exactly_flat_upto_R']}; first R that moves it: "
              f"{a['ladder_first_R_that_moves']}")
        print(f"  ladder max move over the full grid: "
              f"{a['ladder_max_move_full_range']:.4f}")

    if args.stage in ("all", "capacity_off"):
        print("\n--- MEASUREMENT 2, CONVECTION-OFF ARM ---")
        print("    (the protocol capacity_ranking_robustness.md sec. 2d runs, "
              "and the only protocol under which a parameter-division "
              "emulation is exact)")
        cap = stage_capacity(alpha, convection=False, r_grid=R_GRID_OFF)
        cap["analysis"] = analyse_capacity(cap["rows"])
        _dump({"capacity_convection_off": cap})
        a = cap["analysis"]
        print(f"  [conv OFF] falloff move at R = 39: "
              f"{a['falloff_move_at_R39']:.4f}; L {a['L_um_R1']:.1f} -> "
              f"{a['L_um_R39']:.1f} um; ladder max move R = 1..39 "
              f"{a['ladder_max_move_R1_to_39']:.6f}; exactly flat to R = "
              f"{a['ladder_exactly_flat_upto_R']}; full-grid ladder move "
              f"{a['ladder_max_move_full_range']:.4f}")

    _OUT["elapsed_s"] = time.time() - t0
    _OUT.pop("partial_capacity", None)
    _OUT.pop("partial_capacity_off", None)
    _dump()
    print(f"\nwrote {_DATA}  ({_OUT['elapsed_s']:.0f}s)")
    return _OUT


if __name__ == "__main__":
    main()
