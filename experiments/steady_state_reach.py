"""IS THE REACH MAP BETTER AS A STEADY-STATE OBJECT THAN A FINITE-WINDOW ONE?

WHY THIS EXISTS
---------------
``app/server/reach.py`` publishes a scale-free reach map: a ratio between two
points of one exposure field. Its capacity-invariance disclosure has been
rewritten twice in two days, and the second version rests on a claim this script
exists to test:

    *"THIS WINDOW IS NOT CONVERGED: the gap to the infinite-window limit is
    0.01302 for temozolomide at the shipping 6 h and decays as 281/T with no
    plateau, so reaching 0.001 would need about 78 h."*

**That is a statement about an ESTIMATOR, not about the field**, and the
distinction decides the whole design. Write ``AUC(x, T) = C_ss(x)*T + J(x, T)``
where ``J`` is the integrated transient. If ``C`` converges then ``J`` tends to a
finite constant, so the AUC-based ratio's error tends to ``const/T`` — an
``O(1/T)`` tail that **can never plateau, no matter how converged the field is**.
The field ``C(x, t)`` itself converges at the reaction rate, i.e.
*exponentially*. Those are completely different situations and only one of them
is a reason to distrust the map.

Sec. 1 measures which it is, rather than deriving it.

WHAT THE STEADY-STATE ROUTE IS
------------------------------
Not a longer window and not an integration to a tolerance. The shipped stepper
(``physics.diffusion.simulate_drug_transport``) advances

    (I - dt*L + dt*K) c_new  =  c + dt*S - dt*Adv(c)

so its own fixed point ``c_new = c = c*`` satisfies, with every ``dt`` cancelling
identically,

    (K - L + Adv) c*  =  S

— one sparse linear system, built from the engine's OWN ``L``
(``_build_diffusion_operator``, imported not copied) and its own ``K``. That is
``app.server.reach.steady_state_concentration``. It costs ONE factorisation and
ONE triangular solve where the 6 h march costs one factorisation and 2,160
triangular solves, so it is strictly cheaper than the object it replaces.

Sec. 2 measures the cost. Sec. 3 checks the solve against a long march rather
than trusting the algebra.

WHY IT SHOULD BE EXACTLY CAPACITY-INVARIANT
-------------------------------------------
A tissue-capacity term ``R`` enters as ``R*dC/dt = ...``: it multiplies the TIME
DERIVATIVE AND NOTHING ELSE, so it is absent from the steady problem by
inspection. ``results/capacity_ranking_robustness.md`` sec. 2b already measured
that at 1e-9 on a homogeneous probe. Sec. 4 here measures it on the flagship
section, on the two quantities the product actually publishes, against the
finite-window residual it is meant to remove (falloff movement 0.3772 at R = 39).

WHAT IS CHANGED
---------------
Nothing under ``tmesim/``. ``calibrate_alpha()`` is asserted live at
``22.223477711172563``. Drugs are compiled through the flagship path the shipped
real-tissue panel uses (``calibrated_drugparams`` -> ``apply_barrier_to_drugparams``).
Nothing is tuned.

Run:
    ./.venv/bin/python experiments/steady_state_reach.py --stage operator
    ./.venv/bin/python experiments/steady_state_reach.py --stage estimator
    ./.venv/bin/python experiments/steady_state_reach.py --stage cost
    ./.venv/bin/python experiments/steady_state_reach.py --stage capacity
    ./.venv/bin/python experiments/steady_state_reach.py --stage spread
    ./.venv/bin/python experiments/steady_state_reach.py --stage all
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

from app.server.reach import (  # noqa: E402
    REACH_LADDER,
    _advection_matrix,
    steady_state_concentration,
)
from tmesim.core import (  # noqa: E402
    DrugParams,
    OxygenConfig,
    PDEConfig,
    SimConfig,
)
from tmesim.molecular.binding_site_barrier import (  # noqa: E402
    apply_barrier_to_drugparams,
)
from tmesim.molecular.library import get_molecule  # noqa: E402
from tmesim.physics.diffusion import (  # noqa: E402
    _advective_divergence,
    simulate_drug_transport,
)
from tmesim.physics.oxygen import solve_oxygen  # noqa: E402
from tmesim.physics.pressure import solve_ifp  # noqa: E402
from tmesim.sim.engine import simulate_response  # noqa: E402
from tmesim.treatment.drugs import delivery_source  # noqa: E402
from tmesim.validation import vessel_profile as VP  # noqa: E402
from tmesim.validation.calibration import (  # noqa: E402
    calibrate_alpha,
    calibrated_drugparams,
)
from tmesim.validation.delivery_metrics import (  # noqa: E402
    _sample_at_cells,
    exposure_percentiles,
    tumor_cell_exposure,
)
from tmesim.validation.vessel_profile import (  # noqa: E402
    direct_vessel_distance_um,
    distance_profile,
    endothelial_xy,
    usable_bins,
)

SHIPPED_ALPHA = 22.223477711172563

#: The flagship shipping configuration: 10x Xenium GBM section, dx = 12 um,
#: mask repair 2 — the same one ``results/delivery_spread_noise.md`` calls
#: ``dx12_repair2``.
FLAGSHIP = {"dx_um": 12.0, "repair": 2}

#: The shipped exposure window, for every finite-window arm below.
WINDOW_S = 3600.0 * 6

#: The headline pair the delivery half of the product rests on.
PAIR: Sequence[str] = ("temozolomide", "doxorubicin")

#: The nuisance grid ``results/delivery_spread_noise.md`` sec. 2 swept. The
#: WINDOW axis of that grid is deliberately absent here: under steady state
#: there is no window to vary, which is the whole point.
CONFIGS: List[Dict] = [
    {"key": "dx12_repair2", "dx_um": 12.0, "repair": 2, "shipped": True},
    {"key": "dx12_native", "dx_um": 12.0, "repair": 0, "shipped": False},
    {"key": "dx16_repair2", "dx_um": 16.0, "repair": 2, "shipped": False},
    {"key": "dx20_native", "dx_um": 20.0, "repair": 0, "shipped": False},
    {"key": "dx30_native", "dx_um": 30.0, "repair": 0, "shipped": False},
]

#: The finite-window P90/P10 pair the steady-state values have to be compared
#: against. ``results/delivery_spread_noise.md`` sec. 0.
FINITE_SPREAD_SHIPPED = {"temozolomide": 2.111, "doxorubicin": 3.942}
FINITE_SPREAD_CLAIMED = {"temozolomide": 1.99, "doxorubicin": 3.80}

#: The capacity value ``results/capacity_ranking_robustness.md`` treats as the
#: central estimate, and the one the 0.3772 falloff residual is quoted at.
CAPACITY_R = 39.0
CAPACITY_GRID = (1.0, 3.2, 8.0, 20.0, 39.0, 120.0, 4184.0)

_RESULTS = os.path.join(_ROOT, "results")
_DATA = os.path.join(_RESULTS, "steady_state_reach.json")


def _p(*a):
    print(*a, flush=True)


def _assert_alpha() -> float:
    a = calibrate_alpha()
    if a != SHIPPED_ALPHA:
        raise AssertionError(
            f"calibrate_alpha() = {a!r}, expected {SHIPPED_ALPHA!r}")
    return a


def compile_flagship(name: str, alpha: float) -> DrugParams:
    """The pipeline the shipped real-tissue panel uses."""
    return apply_barrier_to_drugparams(calibrated_drugparams(name, alpha),
                                       get_molecule(name))


def retard(dp: DrugParams, R: float) -> DrugParams:
    """A capacity term ``R``, applied EXACTLY as
    ``experiments/capacity_ranking_robustness.py`` applies it (its ``retard``),
    so the residuals here are comparable to the 0.3772 that file reports."""
    return dataclasses.replace(
        dp,
        diffusivity=dp.diffusivity / R,
        binding=dp.binding / R,
        clearance=dp.clearance / R,
        bbb_penetrance=dp.bbb_penetrance / R,
    )


class Tissue:
    """One reconstruction, with the oxygen / IFP stages the engine runs before
    transport already done — so a steady-state solve can reuse them instead of
    re-deriving them, and so the velocity handed to the solve is the SAME array
    ``simulate_response`` hands its own transport call."""

    def __init__(self, dx_um: float, repair: int):
        t0 = time.time()
        micro, _cn, _cr = VP.build_arm_microenvironment(dx_um=dx_um,
                                                        repair=repair)
        self.micro = micro
        self.build_s = time.time() - t0
        t0 = time.time()
        self.oxygen = solve_oxygen(micro, OxygenConfig())
        micro.fields.oxygen = self.oxygen
        self.ifp, self.velocity = solve_ifp(micro)
        micro.fields.ifp = self.ifp
        micro.fields.velocity = self.velocity
        self.pre_s = time.time() - t0
        self.tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
        t0 = time.time()
        self.distance_um = direct_vessel_distance_um(micro.domain,
                                                     endothelial_xy(micro.cells))
        self.dist_s = time.time() - t0

    def source(self, dp: DrugParams, plasma: float = 1.0) -> np.ndarray:
        return delivery_source(self.micro, dp, plasma)

    def steady(self, dp: DrugParams, convection: bool = True,
               plasma: float = 1.0) -> np.ndarray:
        return steady_state_concentration(
            self.micro, dp, self.source(dp, plasma),
            velocity=(self.velocity if convection else None))

    def march(self, dp: DrugParams, total_time: float, convection: bool = True,
              dt: float = 10.0):
        pde = PDEConfig(dt=dt, total_time=float(total_time),
                        convection=convection, save_every=10 ** 9)
        return simulate_drug_transport(
            self.micro, dp, self.source(dp),
            pde, velocity=(self.velocity if convection else None))


# --------------------------------------------------------------------------- #
# The two reductions the product publishes, taken off an arbitrary field
# --------------------------------------------------------------------------- #
def falloff_of(field: np.ndarray, ts: Tissue):
    """The normalised distance falloff: bin means divided by the first bin.

    Identical binning to ``reach.compute_reach`` — the same distance map, the
    same ``vessel_profile.distance_profile`` / ``usable_bins`` — so the numbers
    here are the numbers the payload would carry.
    """
    prof = distance_profile(np.asarray(field, dtype=np.float64),
                            ts.distance_um, ts.tissue)
    d, m, n = usable_bins(prof)
    if m.size == 0 or float(m[0]) <= 0.0:
        return None
    return {"distance_um": [float(v) for v in d],
            "normalised": [float(v) for v in (m / float(m[0]))],
            "count": [int(v) for v in n], "E_ref": float(m[0])}


def ladder_of(field: np.ndarray, ts: Tissue, e_ref: float):
    """``F(x)`` over :data:`REACH_LADDER`, on the VIABLE tumour cells, exactly
    as ``reach.compute_reach`` defines them."""
    res = _CellProbe(ts, field)
    r = np.sort(res.viable_values / e_ref)
    x = np.asarray(list(REACH_LADDER), dtype=np.float64)
    F = np.searchsorted(r, x, side="left").astype(np.float64) / float(r.size)
    return {"x": [float(v) for v in x], "F": [float(v) for v in F],
            "n_viable": int(r.size)}


class _CellProbe:
    """Per-tumour-cell values of an arbitrary grid field, with the engine's own
    'viable' band (o >= 2 * hypoxia_threshold), i.e. the population
    ``delivery_metrics.hypoxic_vs_normoxic_exposure_ratio`` uses."""

    def __init__(self, ts: Tissue, field: np.ndarray):
        micro = ts.micro
        tumor = micro.cells.tumor_mask()
        self.values = np.asarray(
            _sample_at_cells(np.asarray(field, dtype=np.float64), micro, tumor),
            dtype=np.float64)
        o2 = np.asarray(_sample_at_cells(ts.oxygen, micro, tumor),
                        dtype=np.float64)
        self.viable = o2 >= 2.0 * OxygenConfig().hypoxia_threshold
        self.viable_values = self.values[self.viable]
        self.n_tumor = int(tumor.sum())


def p90_p10(values: np.ndarray) -> float:
    lo = float(np.percentile(values, 10.0))
    hi = float(np.percentile(values, 90.0))
    return float(hi / lo) if lo > 0.0 else float("nan")


# --------------------------------------------------------------------------- #
# 0. The operator check: the ONE piece of the engine re-expressed here
# --------------------------------------------------------------------------- #
def stage_operator() -> Dict:
    _p("=" * 78)
    _p("0. THE ADVECTION MATRIX IS THE ENGINE'S OWN ADVECTION OPERATOR")
    _p("=" * 78)
    rng = np.random.default_rng(0)
    worst = 0.0
    rows = []
    for trial in range(6):
        ny, nx = int(rng.integers(4, 20)), int(rng.integers(4, 20))
        dx = float(rng.uniform(5.0, 40.0))
        tissue = rng.random((ny, nx)) > 0.25
        vx = rng.normal(size=(ny, nx))
        vy = rng.normal(size=(ny, nx))
        c = rng.random((ny, nx))
        c[~tissue] = 0.0
        direct = _advective_divergence(c, vx, vy, tissue, dx)
        viamat = (_advection_matrix(vx, vy, tissue, dx)
                  @ c.ravel()).reshape(ny, nx)
        scale = max(float(np.abs(direct).max()), 1e-300)
        err = float(np.abs(direct - viamat).max()) / scale
        worst = max(worst, err)
        rows.append({"ny": ny, "nx": nx, "dx": dx, "rel_err": err})
        _p(f"  {ny:>3}x{nx:<3} dx={dx:6.2f}   max rel err = {err:.3e}")
    _p(f"\n  WORST over 6 random cases: {worst:.3e}")
    _p("  -> the matrix form APPLIES the same operator "
       f"({'PASS' if worst < 1e-13 else 'FAIL'})")
    return {"worst_rel_err": worst, "cases": rows}


# --------------------------------------------------------------------------- #
# 1. Is the field converged, or is 281/T a property of the estimator?
# --------------------------------------------------------------------------- #
def stage_estimator(ts: Tissue, dp: DrugParams) -> Dict:
    _p("=" * 78)
    _p("1. WHAT 281/T IS ACTUALLY MEASURING")
    _p("=" * 78)
    k = float(dp.clearance) + float(dp.binding)
    tau = 1.0 / k
    _p(f"  k = clearance + binding = {k:.6g} 1/s   ->   tau = 1/k = {tau:.1f} s "
       f"= {tau/3600:.2f} h")
    _p(f"  shipping window {WINDOW_S:.0f} s = {WINDOW_S/tau:.2f} tau\n")

    css = ts.steady(dp, convection=True)
    ref_fall = falloff_of(css, ts)
    ref_norm = np.asarray(ref_fall["normalised"])
    tis = ts.tissue
    css_scale = float(css[tis].max())

    rows = []
    windows = (1800.0, 3600.0, WINDOW_S, 3600.0 * 12, 3600.0 * 24, 3600.0 * 48)
    for T in windows:
        t0 = time.time()
        tr = ts.march(dp, T, convection=True)
        el = time.time() - t0
        # (a) the FIELD's own distance from the steady field
        field_err = float(np.abs(tr.final_concentration[tis] - css[tis]).max()
                          ) / css_scale
        # (b) AUC(T) = C_ss*T + J(T): is J bounded?
        j = tr.auc[tis] - css[tis] * T
        j_mean = float(j.mean())
        j_max = float(np.abs(j).max())
        # (c) the published REDUCTION's gap: the AUC falloff against the steady
        #     falloff, which is the quantity 0.01302 was reported for
        auc_fall = falloff_of(tr.auc, ts)
        gap = float(np.abs(np.asarray(auc_fall["normalised"]) - ref_norm).max())
        rows.append({"T_s": T, "field_rel_err": field_err,
                     "J_mean_uM_s": j_mean, "J_absmax_uM_s": j_max,
                     "auc_falloff_gap": gap, "gap_times_T": gap * T,
                     "march_s": el})
        _p(f"  T = {T/3600:5.1f} h  ({el:5.1f}s):  "
           f"|C(T)-C_ss|/max = {field_err:9.3e}   "
           f"AUC falloff gap = {gap:.6f}   gap*T = {gap*T:8.1f}   "
           f"mean J = {j_mean:.5g} uM*s")

    _p("\n  READ THIS COLUMNWISE:")
    ratios = []
    for i in range(1, len(rows)):
        a, b = rows[i - 1], rows[i]
        if a["field_rel_err"] > 0 and b["field_rel_err"] > 0:
            # exponential convergence predicts err(T2)/err(T1) = exp(-k*(T2-T1))
            pred = float(np.exp(-k * (b["T_s"] - a["T_s"])))
            obs = b["field_rel_err"] / a["field_rel_err"]
            ratios.append({"from_s": a["T_s"], "to_s": b["T_s"],
                           "observed": obs, "exp_minus_k_dT": pred})
            _p(f"    field error {a['T_s']/3600:4.1f}h -> {b['T_s']/3600:4.1f}h: "
               f"observed x{obs:.3e}   exp(-k*dT) = {pred:.3e}")
    gt = [r["gap_times_T"] for r in rows]
    _p(f"    gap*T over the whole ladder: {min(gt):.1f} .. {max(gt):.1f}  "
       f"(a CONSTANT would mean the gap is exactly const/T)")
    jm = [r["J_mean_uM_s"] for r in rows]
    _p(f"    mean J over the ladder:      {min(jm):.6g} .. {max(jm):.6g} uM*s  "
       f"(BOUNDED => AUC = C_ss*T + const)")
    return {"k_per_s": k, "tau_s": tau, "rows": rows, "ratios": ratios}


# --------------------------------------------------------------------------- #
# 2 + 3. Cost, and agreement with a long march
# --------------------------------------------------------------------------- #
def stage_cost(ts: Tissue, dp: DrugParams) -> Dict:
    _p("=" * 78)
    _p("2. COST, AND AGREEMENT WITH A LONG MARCH")
    _p("=" * 78)
    out: Dict = {"grid": [int(ts.micro.domain.ny), int(ts.micro.domain.nx)],
                 "dx_um": float(ts.micro.domain.dx),
                 "build_s": ts.build_s, "pre_s": ts.pre_s,
                 "distance_s": ts.dist_s}
    _p(f"  grid {ts.micro.domain.ny}x{ts.micro.domain.nx} at dx = "
       f"{ts.micro.domain.dx:g} um "
       f"({ts.micro.domain.ny*ts.micro.domain.nx:,} unknowns)")

    for conv in (True, False):
        ts.steady(dp, convection=conv)          # warm any lazy import
        t0 = time.time()
        css = ts.steady(dp, convection=conv)
        el = time.time() - t0
        out[f"steady_s_convection_{conv}"] = el
        out[f"steady_min_uM_convection_{conv}"] = float(css.min())
        out[f"steady_max_uM_convection_{conv}"] = float(css[ts.tissue].max())
        _p(f"  steady solve, convection={str(conv):<5}: {el:6.2f} s   "
           f"C_ss in [{css.min():.4g}, {css[ts.tissue].max():.4g}] uM")

    t0 = time.time()
    tr6 = ts.march(dp, WINDOW_S, convection=True)
    out["march_6h_s"] = time.time() - t0
    _p(f"  6 h march (the object being replaced): {out['march_6h_s']:6.2f} s   "
       f"= {out['march_6h_s'] / out['steady_s_convection_True']:.1f}x the steady "
       f"solve")

    _p("\n  AGREEMENT WITH A LONG MARCH (the algebra, checked rather than trusted)")
    css_c = ts.steady(dp, convection=True)
    css_n = ts.steady(dp, convection=False)
    tis = ts.tissue
    agree = []
    for conv, css in ((True, css_c), (False, css_n)):
        scale = float(css[tis].max())
        for T in (3600.0 * 24, 3600.0 * 72):
            t0 = time.time()
            tr = ts.march(dp, T, convection=conv)
            el = time.time() - t0
            err = float(np.abs(tr.final_concentration[tis] - css[tis]).max()
                        ) / scale
            agree.append({"convection": conv, "T_s": T, "rel_err": err,
                          "march_s": el})
            _p(f"    convection={str(conv):<5} T = {T/3600:4.0f} h "
               f"({el:6.1f}s): max |C(T) - C_ss| / max C_ss = {err:.3e}")
    out["agreement"] = agree
    out["worst_agreement_rel_err"] = max(a["rel_err"] for a in agree)
    return out


# --------------------------------------------------------------------------- #
# 4. Capacity invariance, re-measured on the flagship section
# --------------------------------------------------------------------------- #
def stage_capacity(ts: Tissue, dp: DrugParams) -> Dict:
    _p("=" * 78)
    _p("4. CAPACITY INVARIANCE, MEASURED NOT ASSERTED")
    _p("=" * 78)
    _p("  Convection is OFF in every arm here, for the reason "
       "capacity_ranking_robustness.md sec. 2(d) gives: R is emulated by "
       "dividing D, binding, clearance and bbb_penetrance, and the engine's "
       "velocity field cannot be divided from outside. Same protocol, so the "
       "numbers are comparable to that file's 0.3772.\n")

    base_ss = ts.steady(dp, convection=False)
    base_ss_fall = np.asarray(falloff_of(base_ss, ts)["normalised"])
    base_ss_ref = falloff_of(base_ss, ts)["E_ref"]
    base_ss_lad = np.asarray(ladder_of(base_ss, ts, base_ss_ref)["F"])

    base_auc = ts.march(dp, WINDOW_S, convection=False).auc
    base_auc_fall = np.asarray(falloff_of(base_auc, ts)["normalised"])
    base_auc_ref = falloff_of(base_auc, ts)["E_ref"]
    base_auc_lad = np.asarray(ladder_of(base_auc, ts, base_auc_ref)["F"])

    rows = []
    for R in CAPACITY_GRID:
        dpr = retard(dp, R)
        ss = ts.steady(dpr, convection=False)
        f_ss = falloff_of(ss, ts)
        ss_fall = np.asarray(f_ss["normalised"])
        ss_lad = np.asarray(ladder_of(ss, ts, f_ss["E_ref"])["F"])
        auc = ts.march(dpr, WINDOW_S, convection=False).auc
        f_auc = falloff_of(auc, ts)
        auc_fall = np.asarray(f_auc["normalised"])
        auc_lad = np.asarray(ladder_of(auc, ts, f_auc["E_ref"])["F"])
        row = {
            "R": R,
            "steady_falloff_absmove": float(np.abs(ss_fall - base_ss_fall).max()),
            "steady_ladder_absmove": float(np.abs(ss_lad - base_ss_lad).max()),
            "finite_falloff_absmove": float(np.abs(auc_fall - base_auc_fall).max()),
            "finite_ladder_absmove": float(np.abs(auc_lad - base_auc_lad).max()),
        }
        rows.append(row)
        _p(f"  R = {R:>7.1f}   STEADY falloff move {row['steady_falloff_absmove']:.3e}  "
           f"ladder {row['steady_ladder_absmove']:.3e}   |   "
           f"FINITE(6h) falloff move {row['finite_falloff_absmove']:.4f}  "
           f"ladder {row['finite_ladder_absmove']:.4f}")

    at39 = [r for r in rows if r["R"] == CAPACITY_R][0]
    _p(f"\n  AT R = {CAPACITY_R:g}, the number the disclosure quotes:")
    _p(f"    finite-window falloff movement = {at39['finite_falloff_absmove']:.4f}  "
       f"(capacity_ranking_robustness.md reports 0.3772)")
    _p(f"    STEADY-STATE falloff movement  = {at39['steady_falloff_absmove']:.3e}")
    _p(f"    worst steady falloff movement over R = 1..4184: "
       f"{max(r['steady_falloff_absmove'] for r in rows):.3e}")
    _p(f"    worst steady ladder  movement over R = 1..4184: "
       f"{max(r['steady_ladder_absmove'] for r in rows):.3e}")
    return {"rows": rows,
            "steady_worst_falloff_absmove": float(
                max(r["steady_falloff_absmove"] for r in rows)),
            "steady_worst_ladder_absmove": float(
                max(r["steady_ladder_absmove"] for r in rows)),
            "finite_falloff_absmove_at_R39": float(at39["finite_falloff_absmove"]),
            "steady_falloff_absmove_at_R39": float(at39["steady_falloff_absmove"])}


def stage_capacity_convection(ts: Tissue, dp: DrugParams) -> Dict:
    """The counter-case the capacity protocol hides: R with convection ON.

    The published emulation divides ``D``, ``binding``, ``clearance`` and
    ``bbb_penetrance`` by ``R`` and leaves the velocity field alone, because it
    is built inside ``simulate_response`` and cannot be reached from outside
    (``capacity_ranking_robustness.md`` sec. 2(d)). Multiply the emulated steady
    equation through by ``R``:

    ``div(D grad C) - R*div(v C) + S - k*C = 0``

    The advection term keeps its ``R``. So with convection ON the emulation is
    NOT a capacity change — it is a capacity change plus a factor-``R`` boost to
    the Peclet number, and the steady map moves. That is a property of the
    EMULATION, not of the map: in the real model ``R`` multiplies ``dC/dt`` and
    the steady problem has none.

    Three arms, so the claim is identified rather than asserted:
      (a) convection OFF, published emulation      -> should be ~0
      (b) convection ON,  published emulation      -> moves; measure how much
      (c) convection ON,  emulation + ``v/R`` too  -> should return to ~0
    """
    _p("=" * 78)
    _p("4b. THE SAME SWEEP WITH CONVECTION ON — WHERE THE CLAIM COULD BE MISREAD")
    _p("=" * 78)
    rows = []
    for R in (39.0, 4184.0):
        dpr = retard(dp, R)
        arms = {}
        for label, conv, vel in (
                ("a_conv_off", False, None),
                ("b_conv_on_published", True, ts.velocity),
                ("c_conv_on_velocity_scaled", True, ts.velocity / R)):
            base_vel = None if not conv else ts.velocity
            base = steady_state_concentration(
                ts.micro, dp, ts.source(dp, 1.0), velocity=base_vel)
            got = steady_state_concentration(
                ts.micro, dpr, ts.source(dpr, 1.0), velocity=vel)
            fb, fg = falloff_of(base, ts), falloff_of(got, ts)
            fall = float(np.abs(np.asarray(fg["normalised"])
                                - np.asarray(fb["normalised"])).max())
            lad = float(np.abs(
                np.asarray(ladder_of(got, ts, fg["E_ref"])["F"])
                - np.asarray(ladder_of(base, ts, fb["E_ref"])["F"])).max())
            arms[label] = {"falloff_absmove": fall, "ladder_absmove": lad}
            _p(f"  R = {R:>7.1f}  {label:<26} falloff {fall:.3e}   "
               f"ladder {lad:.3e}")
        rows.append({"R": R, **arms})
    return {"rows": rows}


# --------------------------------------------------------------------------- #
# 5. Does P90/P10 still discriminate?  THE STOP-AND-REPORT TEST
# --------------------------------------------------------------------------- #
def stage_spread(alpha: float, tissues: Dict[str, Tissue]) -> Dict:
    _p("=" * 78)
    _p("5. THE DELIVERY SPREAD P90/P10 UNDER STEADY STATE")
    _p("=" * 78)
    drugs = {n: compile_flagship(n, alpha) for n in PAIR}
    rows = []
    for cfg in CONFIGS:
        ts = tissues[cfg["key"]]
        for name in PAIR:
            dp = drugs[name]
            t0 = time.time()
            css = ts.steady(dp, convection=True)
            el = time.time() - t0
            probe = _CellProbe(ts, css)
            ss_ratio = p90_p10(probe.values)
            t1 = time.time()
            res = simulate_response(
                ts.micro,
                SimConfig(drug=dp,
                          pde=PDEConfig(dt=10.0, total_time=WINDOW_S,
                                        convection=True, save_every=500),
                          seed=0))
            fin_el = time.time() - t1
            exp = tumor_cell_exposure(ts.micro, res)
            pct = exposure_percentiles(exp["exposure_uM"], exp["auc_uM_s"])
            rows.append({
                "config": cfg["key"], "dx_um": cfg["dx_um"],
                "repair": cfg["repair"], "shipped": bool(cfg["shipped"]),
                "drug": name,
                "steady_p90_over_p10": ss_ratio,
                "finite_p90_over_p10": float(pct["E_p90_over_p10"]),
                "n_tumor": int(probe.n_tumor),
                "steady_s": el, "finite_s": fin_el,
            })
            _p(f"  {cfg['key']:<14} {name:<14}  STEADY {ss_ratio:7.3f} "
               f"({el:5.2f}s)   FINITE-6h {pct['E_p90_over_p10']:7.3f} "
               f"({fin_el:5.1f}s)   n_tumor={probe.n_tumor}")

    out: Dict = {"rows": rows}
    _p("\n  SEPARATION")
    for key, label in (("steady_p90_over_p10", "STEADY-STATE"),
                       ("finite_p90_over_p10", "FINITE 6 h  ")):
        spans = {}
        for name in PAIR:
            v = [r[key] for r in rows if r["drug"] == name and np.isfinite(r[key])]
            sh = [r[key] for r in rows if r["drug"] == name and r["shipped"]]
            spans[name] = {"min": min(v), "max": max(v), "n": len(v),
                           "shipped": sh[0] if sh else None}
        a, b = PAIR
        disjoint = (spans[a]["max"] < spans[b]["min"]
                    or spans[b]["max"] < spans[a]["min"])
        gap = (spans[b]["min"] - spans[a]["max"]
               if spans[a]["max"] < spans[b]["min"]
               else spans[a]["min"] - spans[b]["max"])
        worst = (spans[b]["min"] / spans[a]["max"]
                 if spans[a]["max"] < spans[b]["min"] else float("nan"))
        out[key] = {"spans": spans, "disjoint": bool(disjoint),
                    "gap": float(gap), "worst_case_ratio": float(worst),
                    "shipped_ratio": float(spans[b]["shipped"]
                                           / spans[a]["shipped"])}
        _p(f"    {label}: tmz {spans[a]['shipped']:.3f} "
           f"[{spans[a]['min']:.3f}-{spans[a]['max']:.3f}]   "
           f"dox {spans[b]['shipped']:.3f} "
           f"[{spans[b]['min']:.3f}-{spans[b]['max']:.3f}]   "
           f"DISJOINT={disjoint}  gap={gap:.3f}  "
           f"worst-case dox/tmz = {worst:.3f}x  "
           f"shipped dox/tmz = {out[key]['shipped_ratio']:.3f}x")
    return out


# --------------------------------------------------------------------------- #
# 6. LOOKING FOR THE COUNTER-CASE ON PURPOSE
# --------------------------------------------------------------------------- #
def stage_pathology(alpha: float) -> Dict:
    """Where a ONE-SHOT DIRECT SOLVE should break if it is going to.

    Three measurements pointed at steady state and all three were about
    invariance and signal. None of them was about whether a sparse direct solve
    is well behaved on the tissues this product actually gets. So this section
    goes looking for the failure rather than waiting for it:

    * **fine grids** — the sparse factorisation is the whole cost, and fill-in
      on a 2-D stencil grows faster than the unknown count. If anything makes
      the route unshippable it is dx = 8 um, not dx = 12.
    * **near-zero clearance** — ``k -> 0`` sends the penetration length
      ``L = sqrt(D/k)`` to infinity and the system towards a pure Neumann
      Laplacian, which is SINGULAR. That is the conditioning cliff, and the
      question is how close the real panel gets to it.
    * **sparse vasculature** — a section with few endothelial cells has a
      near-empty source. The map should REFUSE (no usable first bin / a
      non-positive reference), not return a ratio of two round-off values.
    * **non-negativity** — the time-stepper clamps ``C`` at zero every step. A
      direct solve does not clamp anything, so if the discrete steady state
      were not non-negative the two routes would disagree by construction.
    """
    _p("=" * 78)
    _p("6. THE COUNTER-CASE HUNT")
    _p("=" * 78)
    out: Dict = {}

    _p("  (a) grid refinement — cost and solve accuracy")
    rows = []
    for dx in (30.0, 20.0, 16.0, 12.0, 8.0):
        try:
            ts = Tissue(dx, 2 if dx in (12.0, 16.0) else 0)
        except Exception as exc:                      # pragma: no cover
            rows.append({"dx_um": dx, "error": repr(exc)})
            _p(f"    dx = {dx:4.1f} um   BUILD FAILED: {exc!r}")
            continue
        dp = compile_flagship("temozolomide", alpha)
        t0 = time.time()
        try:
            css = ts.steady(dp, convection=True)
            el = time.time() - t0
        except Exception as exc:                      # pragma: no cover
            rows.append({"dx_um": dx, "error": repr(exc)})
            _p(f"    dx = {dx:4.1f} um   SOLVE FAILED: {exc!r}")
            del ts
            continue
        # Residual of the system actually solved, relative to the source norm:
        # the direct-solve equivalent of a convergence criterion.
        res = _residual(ts, dp, css, convection=True)
        n = ts.micro.domain.ny * ts.micro.domain.nx
        rows.append({"dx_um": dx, "n_unknowns": int(n), "solve_s": el,
                     "rel_residual": res, "min_uM": float(css.min()),
                     "max_uM": float(css[ts.tissue].max())})
        _p(f"    dx = {dx:4.1f} um  n = {n:>8,}   solve {el:6.2f} s   "
           f"rel residual {res:.3e}   min C_ss = {css.min():.3e}")
        del ts, css
    out["grid"] = rows

    _p("\n  (b) clearance -> 0: how close is the conditioning cliff?")
    ts = Tissue(FLAGSHIP["dx_um"], FLAGSHIP["repair"])
    dp = compile_flagship("temozolomide", alpha)
    k0 = float(dp.clearance) + float(dp.binding)
    krows = []
    for fac in (1.0, 1e-2, 1e-4, 1e-6, 1e-8):
        dpk = dataclasses.replace(dp, clearance=dp.clearance * fac,
                                  binding=dp.binding * fac)
        try:
            css = ts.steady(dpk, convection=True)
            res = _residual(ts, dpk, css, convection=True)
            L = float(np.sqrt(dpk.diffusivity / (k0 * fac))) if fac > 0 else np.inf
            krows.append({"k_factor": fac, "k_per_s": k0 * fac, "L_um": L,
                          "rel_residual": res,
                          "max_uM": float(css[ts.tissue].max()),
                          "min_uM": float(css.min())})
            _p(f"    k x {fac:<7.0e}  k = {k0*fac:.3e} 1/s  "
               f"L = sqrt(D/k) = {L:11.1f} um   rel residual {res:.3e}   "
               f"max C_ss = {css[ts.tissue].max():.4g} uM  "
               f"min = {css.min():.3e}")
        except Exception as exc:
            krows.append({"k_factor": fac, "error": repr(exc)})
            _p(f"    k x {fac:<7.0e}  SOLVE FAILED: {exc!r}")
    out["clearance"] = krows
    out["k_shipping_per_s"] = k0

    _p("\n  (c) sparse / absent vasculature — does it REFUSE?")
    from app.server.reach import compute_reach
    srows = []
    base_vasc = np.array(ts.micro.fields.vascular_density, copy=True)
    res_run = simulate_response(
        ts.micro,
        SimConfig(drug=dp, pde=PDEConfig(dt=10.0, total_time=WINDOW_S,
                                         convection=True, save_every=500),
                  seed=0))
    for frac, label in ((1.0, "full"), (1e-3, "0.1% of vessel density"),
                        (0.0, "no vasculature at all")):
        ts.micro.fields.vascular_density = base_vasc * frac
        try:
            r = compute_reach(ts.micro, res_run, distance_um=ts.distance_um,
                              drug=dp, velocity=ts.velocity, convection=True)
            srows.append({"vascular_factor": frac, "label": label,
                          "null_reason": r["null_reason"],
                          "E_ref": r["E_ref_uM_s"],
                          "headline": r["headline_under_reached_fraction"]})
            _p(f"    {label:<24} null_reason = {r['null_reason']!r}  "
               f"E_ref = {r['E_ref_uM_s']!r}")
        except Exception as exc:
            srows.append({"vascular_factor": frac, "label": label,
                          "raised": repr(exc)})
            _p(f"    {label:<24} RAISED {exc!r}")
    ts.micro.fields.vascular_density = base_vasc
    out["sparse_vasculature"] = srows

    _p("\n  (d) non-negativity of the unclamped direct solve")
    neg = []
    for name in PAIR:
        dpn = compile_flagship(name, alpha)
        css = ts.steady(dpn, convection=True)
        neg.append({"drug": name, "min_uM": float(css.min())})
        _p(f"    {name:<14} min C_ss over the whole grid = {css.min():.3e} uM")
    out["non_negativity"] = neg
    return out


def _residual(ts: Tissue, dp: DrugParams, css: np.ndarray,
              convection: bool) -> float:
    """``||A c - S|| / ||S||`` (max norm) for the system that was solved.

    Rebuilt from the same pieces the solve used, so this measures the linear
    solve's own accuracy — the direct-method analogue of the convergence
    criterion a time march would need.
    """
    import scipy.sparse as sp

    from app.server.reach import _advection_matrix
    from tmesim.physics.diffusion import (
        _build_diffusion_operator,
        _effective_diffusivity,
    )

    micro = ts.micro
    tissue = ts.tissue
    ny, nx = micro.domain.ny, micro.domain.nx
    n = ny * nx
    dx = float(micro.domain.dx)
    d_eff = _effective_diffusivity(dp, micro.fields.ecm_density)
    k = np.full((ny, nx), float(dp.clearance) + float(dp.binding))
    k[~tissue] = 0.0
    diag = k.ravel() + (~tissue).ravel().astype(np.float64)
    a = (-_build_diffusion_operator(d_eff, tissue, dx)
         + sp.csc_matrix((diag, (np.arange(n), np.arange(n))), shape=(n, n)))
    if convection:
        a = a + _advection_matrix(ts.velocity[0], ts.velocity[1], tissue, dx)
    src = np.where(tissue, ts.source(dp, 1.0), 0.0).ravel()
    r = a @ css.ravel() - src
    denom = max(float(np.abs(src).max()), 1e-300)
    return float(np.abs(r).max()) / denom


# --------------------------------------------------------------------------- #
def main() -> Dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "operator", "estimator", "cost",
                             "capacity", "capacity_convection", "spread",
                             "pathology"])
    args = ap.parse_args()
    stage = args.stage

    alpha = _assert_alpha()
    _p(f"calibrate_alpha() = {alpha!r}  (asserted)")
    out: Dict = {"alpha": alpha, "stage": stage, "window_s": WINDOW_S,
                 "flagship": FLAGSHIP}
    t_all = time.time()

    if stage in ("all", "operator"):
        out["operator"] = stage_operator()
    if stage == "operator":
        return out

    tissues: Dict[str, Tissue] = {}
    need = ([c["key"] for c in CONFIGS] if stage in ("all", "spread")
            else ["dx12_repair2"])
    for cfg in CONFIGS:
        if cfg["key"] in need:
            _p(f"building {cfg['key']} ...")
            tissues[cfg["key"]] = Tissue(cfg["dx_um"], cfg["repair"])
            _p(f"  {cfg['key']}: grid "
               f"{tissues[cfg['key']].micro.domain.shape} "
               f"build {tissues[cfg['key']].build_s:.1f}s "
               f"oxygen+ifp {tissues[cfg['key']].pre_s:.1f}s "
               f"distance {tissues[cfg['key']].dist_s:.1f}s")
    flag = tissues["dx12_repair2"]
    dp_tmz = compile_flagship("temozolomide", alpha)

    if stage in ("all", "estimator"):
        out["estimator"] = stage_estimator(flag, dp_tmz)
    if stage in ("all", "cost"):
        out["cost"] = stage_cost(flag, dp_tmz)
    if stage in ("all", "capacity"):
        out["capacity"] = stage_capacity(flag, dp_tmz)
    if stage in ("all", "capacity_convection"):
        out["capacity_convection"] = stage_capacity_convection(flag, dp_tmz)
    if stage in ("all", "spread"):
        out["spread"] = stage_spread(alpha, tissues)
    if stage in ("all", "pathology"):
        out["pathology"] = stage_pathology(alpha)

    out["elapsed_s"] = time.time() - t_all
    # MERGED, not overwritten. Each section here is expensive enough to be run
    # on its own (the estimator ladder alone marches 91.5 h of simulated time),
    # so a single-stage run must not silently discard the sections it did not
    # recompute. Whatever this run measured replaces the same key; everything
    # else is left exactly as it was, and ``stage`` records which is which.
    merged: Dict = {}
    if os.path.exists(_DATA):
        try:
            with open(_DATA) as fh:
                merged = json.load(fh)
        except (OSError, ValueError):                 # pragma: no cover
            merged = {}
    stages_run = merged.get("stages_run", [])
    merged.update(out)
    merged["stages_run"] = sorted(set(stages_run) | {stage})
    with open(_DATA, "w") as fh:
        json.dump(merged, fh, indent=1)
    _p(f"\nwrote {_DATA}  (sections present: "
       f"{', '.join(sorted(k for k in merged if isinstance(merged[k], dict)))})")
    _p(f"total {out['elapsed_s']:.1f}s")
    return out


if __name__ == "__main__":
    main()
