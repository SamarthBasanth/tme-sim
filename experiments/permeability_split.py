"""Separating ENTRY ("how much gets in") from SPREAD ("where it goes once in").

Driver for ``results/permeability_split.md``. Four stages, each of which can be
run alone::

    python experiments/permeability_split.py --stage entanglement
    python experiments/permeability_split.py --stage lambda
    python experiments/permeability_split.py --stage tissue
    python experiments/permeability_split.py --stage scope
    python experiments/permeability_split.py            # all four + JSON

WHAT IT MEASURES, AND WHY IN THIS ORDER
=======================================
**entanglement** -- the exact property. Under the shipped one-way law the
solved field is positively homogeneous of degree 1 in the entry term ``PS``, so
``PS`` scales the LEVEL and contributes exactly nothing to the SHAPE. Measured
by sweeping ``PS`` over twelve decades on the arm-A slab and watching the
near:far ratio not move. Then the same sweep with the split ON, which is the
only thing that lets ``PS`` touch shape at all -- and which, in this geometry,
moves it by 0.15 % in the wrong direction.

**lambda** -- the by-product, and the most useful number here. The arm-A
near:far ratio is (very nearly) a function of the single length
``lambda = sqrt(D_eff / (clearance + binding))``. Mapping ratio against lambda
by two independent routes -- vary ``k`` at fixed ``D``, vary ``D`` at fixed
``k`` -- converts the "~14x too shallow" absolute gradient failure into an
error in ONE length, which is a much smaller and much more findable thing.
**This stage adopts nothing.** It maps the model's own readout onto the model's
own coordinate and reports where the measurement falls on that map. No constant
is changed, no value is fitted, and the compiled compounds are re-derived from
``calibrate_alpha()`` at its asserted live value in every stage.

**tissue** -- the verdict test, on the flagship section. Does the split change
a per-cell DOSE / DELIVERY / RESISTANCE call, and does it change what a sweep
of the unmeasured entry scale can do to one?

**scope** -- the honest accounting of what the split costs. ``PS [1/s] =
P [um/s] * (S/V) [1/um]``; the repo has neither factor, and this stage prints
what the shipped law would have to be claiming about ``P`` if it were a
permeability law.

CONFIGURATION IS PRINTED WITH EVERY NUMBER. Nothing here is quoted without the
dx, the window, the convection setting and the exchange setting it came from.

WIRING. Every solve goes through ``tmesim.sim.engine.simulate_response``, which
is the ONE production place that turns ``PDEConfig.transvascular_exchange``
into an ``exchange_k`` argument. This file therefore calls
``simulate_drug_transport`` nowhere and is not a wiring point; see
``tmesim/validation/propagation.py`` for why doing it the other way is how a
"flag on vs flag off" comparison silently becomes a comparison of a run to
itself.

NOTHING SHIPPED IS CHANGED. ``PDEConfig.transvascular_exchange`` and
``transvascular.SPLIT_DEFAULT`` are both asserted ``False`` at the top of every
run. The ``PS`` sweeps construct HYPOTHETICAL compounds by scaling
``bbb_penetrance`` on a copy of the compiled ``DrugParams`` -- which is exactly
scaling ``PS``, since ``bbb_penetrance`` enters the engine through
``delivery_source`` and nowhere else -- and no swept value is adopted anywhere.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import sys
import time
from collections import Counter
from typing import Dict, List, Optional, Sequence

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
if os.path.dirname(os.path.abspath(__file__)) not in sys.path:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tmesim.core import PDEConfig, PKConfig, SimConfig                # noqa: E402
from tmesim.physics.diffusion import ALPHA_TORTUOSITY                 # noqa: E402
from tmesim.sim.engine import simulate_response                       # noqa: E402
from tmesim.sim.mechanism import LABELS                               # noqa: E402
from tmesim.treatment.drugs import PS_SCALE_S_INV                     # noqa: E402
from tmesim.treatment.transvascular import (                          # noqa: E402
    SPLIT_DEFAULT,
    boundary_layer_um,
    entry_elasticity,
    entry_transport_regime,
    implied_permeability_um_s,
    one_way_plasma_multiple,
    saturation_fraction,
)
from tmesim.validation import vessel_profile as VP                    # noqa: E402
from tmesim.validation.calibration import calibrate_alpha             # noqa: E402

import isomer_penetration_map as IPM                                  # noqa: E402

ALPHA_LIVE = 22.223477711172563

#: Arm-A geometry, taken UNCHANGED from ``experiments/isomer_penetration_map.py``
#: so this file's 2.9977 is comparable with that file's published 3.00 without
#: any conversion. Kyle's own windows and thickness.
ECM = 0.0875
SLAB_DX_UM = 2.0
SLAB_THICKNESS_UM = 140.0
SLAB_WINDOW_UM = 30.0
SLAB_EXPOSURE_S = 3600.0

#: Kyle 2004 Table 1, fluorescence column, as ``isomer_penetration_map`` records
#: it. Quoted, never fitted to.
KYLE_NEAR_FAR = {"doxorubicin": (11.0, 5.0), "epirubicin": (12.0, 3.0),
                 "daunorubicin": (50.0, 15.0)}

#: Flagship section arm, unchanged from ``results/kyle2004_penetration.md`` §7.
TISSUE_DX_UM = 12.0
TISSUE_REPAIR = 1
TISSUE_DT_S = 10.0
TISSUE_TOTAL_S = 3600.0 * 6
TISSUE_CONVECTION = True

#: Entry-scale sweep. A multiplier on ``bbb_penetrance``, i.e. on ``PS``. None
#: of these is proposed as a value; the sweep exists to show what the LAW does
#: with an unmeasured scale, which is the whole question.
PS_SWEEP = (1e-6, 1e-3, 1.0, 1e3, 1e6)
TISSUE_PS_SWEEP = (1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3)

#: Lambda coordinates for the ratio map, um. Chosen to bracket both the shipped
#: value (72.5) and both measured targets; NOT chosen to hit anything.
LAMBDA_GRID = (20.0, 25.0, 27.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 72.487, 90.0)

_DATA_PATH = os.path.join(_ROOT, "results", "permeability_split_data.json")


# --------------------------------------------------------------------------- #
def _assert_frozen() -> float:
    """Every stage re-asserts that nothing shipped moved under it."""
    alpha = calibrate_alpha()
    if alpha != ALPHA_LIVE:
        raise AssertionError(
            f"calibrate_alpha() is {alpha!r}, not the asserted live value "
            f"{ALPHA_LIVE!r}. A change that moves it is not additive and this "
            "experiment's numbers are not comparable with results/.")
    if PDEConfig().transvascular_exchange is not False:
        raise AssertionError("PDEConfig.transvascular_exchange is not False; "
                             "the shipped default has been changed")
    if SPLIT_DEFAULT is not False:
        raise AssertionError("transvascular.SPLIT_DEFAULT is not False")
    if PS_SCALE_S_INV != 1.0:
        raise AssertionError(
            f"drugs.PS_SCALE_S_INV is {PS_SCALE_S_INV!r}, not the shipped 1.0; "
            "every absolute number below would be on a different scale")
    return alpha


def _scaled(dp, c: float):
    """A HYPOTHETICAL compound whose entry term is ``c`` times the real one.

    ``bbb_penetrance`` reaches the engine only through ``delivery_source``
    (verified: no other consumer in ``tmesim/`` outside the CLI's printing and
    the library's own definitions), so scaling it is exactly scaling ``PS`` and
    touches neither ``D_eff`` nor ``clearance`` nor ``binding``. Nothing here is
    adopted.
    """
    return dataclasses.replace(dp, bbb_penetrance=float(dp.bbb_penetrance) * c)


def _slab_ratio(dp, exchange: bool, dx_um: float = SLAB_DX_UM,
                exposure_s: float = SLAB_EXPOSURE_S) -> Dict:
    """Kyle's near:far exposure ratio on the arm-A slab, one delivery law."""
    micro, nx = IPM._slab(ECM, SLAB_THICKNESS_UM, dx_um)
    d_eff = float(dp.diffusivity) * max(1.0 - ALPHA_TORTUOSITY * ECM, 1e-6)
    dt = min(2.0, IPM.MCC_CFL_SAFETY * dx_um * dx_um / max(d_eff, 1e-9))
    pde = PDEConfig(convection=False, dt=dt, total_time=float(exposure_s),
                    save_every=100000, transvascular_exchange=bool(exchange))
    res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
    auc = np.asarray(res.fields["auc"], dtype=np.float64).mean(axis=0)
    x = np.arange(nx) * dx_um
    near = auc[x <= SLAB_WINDOW_UM]
    far = auc[x >= SLAB_THICKNESS_UM - SLAB_WINDOW_UM]
    return {"ratio": float(near.mean() / far.mean()),
            "auc_wall": float(auc[0]), "auc_far": float(auc[-1]),
            "profile": auc, "x_um": x, "dt_s": float(dt)}


def _cfg(exchange: bool, dx_um: float = SLAB_DX_UM) -> str:
    return (f"slab {SLAB_THICKNESS_UM:.0f} um, dx = {dx_um:g} um, ecm "
            f"{ECM:g}, {SLAB_EXPOSURE_S / 3600:.0f} h window, convection OFF, "
            f"exchange {'ON' if exchange else 'OFF'}")


# --------------------------------------------------------------------------- #
# STAGE 1 -- the entanglement, measured
# --------------------------------------------------------------------------- #
def stage_entanglement(verbose: bool = True) -> Dict:
    alpha = _assert_frozen()
    comp = IPM.compile_all(alpha, two_route=True)
    dp = comp[IPM.DOX]["drug_params"]
    out: Dict = {"alpha": alpha, "config_off": _cfg(False),
                 "config_on": _cfg(True), "one_way": [], "split": []}

    if verbose:
        print("=" * 82)
        print("STAGE 1 — the entanglement: does the ENTRY term touch the SHAPE?")
        print("=" * 82)
        print(f"\nreproduction check ({_cfg(False)}):")
    base = _slab_ratio(dp, exchange=False)
    out["reproduction_near_far"] = base["ratio"]
    out["published_near_far"] = 3.00
    if verbose:
        print(f"  doxorubicin near:far = {base['ratio']:.4f}  against "
              f"results/kyle2004_penetration.md §6.1's published 3.00")
    if not (2.99 <= base["ratio"] <= 3.01):
        raise AssertionError(
            f"arm A did not reproduce: {base['ratio']!r} is not 3.00. The rest "
            "of this file is not comparable with the published result.")

    if verbose:
        print(f"\n(a) ONE-WAY law ({_cfg(False)}). PS scaled by c:")
        print(f"    {'c':>10}{'near:far':>26}{'AUC at the wall':>20}")
    ref_shape = base["profile"] / base["profile"][0]
    for c in PS_SWEEP:
        r = _slab_ratio(_scaled(dp, c), exchange=False)
        shape_err = float(np.max(np.abs(r["profile"] / r["profile"][0] - ref_shape)))
        out["one_way"].append({"c": c, "ratio": r["ratio"],
                               "auc_wall": r["auc_wall"],
                               "max_shape_deviation": shape_err})
        if verbose:
            print(f"    {c:>10.0e}{r['ratio']!r:>26}{r['auc_wall']:>20.6e}")
    ratios = [row["ratio"] for row in out["one_way"]]
    out["one_way_ratio_spread"] = max(ratios) / min(ratios)
    out["one_way_max_shape_deviation"] = max(
        row["max_shape_deviation"] for row in out["one_way"])
    if verbose:
        print(f"    -> over {len(PS_SWEEP)} values spanning "
              f"{PS_SWEEP[-1] / PS_SWEEP[0]:.0e}, the ratio spread is "
              f"{out['one_way_ratio_spread'] - 1.0:.3e} and the worst profile-"
              f"shape deviation is {out['one_way_max_shape_deviation']:.3e}.")
        print("       The entry term sets the LEVEL and contributes NOTHING to "
              "the shape.")
    if out["one_way_ratio_spread"] - 1.0 > 1e-10:
        raise AssertionError(
            "the one-way law was expected to be exactly homogeneous in PS; the "
            f"near:far ratio moved by {out['one_way_ratio_spread'] - 1.0:.3e}")

    if verbose:
        print(f"\n(b) SPLIT law ({_cfg(True)}). Same sweep:")
        print(f"    {'c':>10}{'PS_wall (1/s)':>16}{'near:far':>12}"
              f"{'AUC at the wall':>20}")
    for c in TISSUE_PS_SWEEP:
        r = _slab_ratio(_scaled(dp, c), exchange=True)
        ps_wall = float(dp.bbb_penetrance) * c * PS_SCALE_S_INV  # phi = 1 at the wall
        out["split"].append({"c": c, "ps_wall_s_inv": ps_wall,
                             "ratio": r["ratio"], "auc_wall": r["auc_wall"]})
        if verbose:
            print(f"    {c:>10.0e}{ps_wall:>16.4e}{r['ratio']:>12.4f}"
                  f"{r['auc_wall']:>20.6e}")
    sr = [row["ratio"] for row in out["split"]]
    out["split_ratio_range"] = [min(sr), max(sr)]
    # The LARGEST move the split can buy, signed: positive would be toward the
    # measurement (steeper), negative away from it.
    moves = [(r - base["ratio"]) / base["ratio"] for r in sr]
    out["split_largest_move"] = max(moves, key=abs)
    out["split_moves_toward_measurement"] = out["split_largest_move"] > 0.0
    if verbose:
        dox_gap = KYLE_NEAR_FAR["doxorubicin"][0] / base["ratio"]
        print(f"    -> the split moves the ratio from {base['ratio']:.4f} to "
              f"{min(sr):.4f}-{max(sr):.4f}: the largest move is "
              f"{100 * out['split_largest_move']:+.2f} %, i.e. "
              f"{'TOWARD' if out['split_moves_toward_measurement'] else 'AWAY FROM'}"
              f" the measurement, against a gap of {dox_gap:.2f}x for "
              f"doxorubicin (11/{base['ratio']:.2f}) and 13.7x for "
              f"daunorubicin (50/3.65, results/kyle2004_penetration.md §6.1).")
        print("       PS is supported only where vascular_density > 0. Kyle's")
        print("       gradient is measured where it is not. The split cannot")
        print("       reach the tissue interior — that is structural, not small.")
    return out


# --------------------------------------------------------------------------- #
# STAGE 2 -- the lambda map
# --------------------------------------------------------------------------- #
def stage_lambda(verbose: bool = True) -> Dict:
    alpha = _assert_frozen()
    comp = IPM.compile_all(alpha, two_route=True)
    tort = 1.0 - ALPHA_TORTUOSITY * ECM
    out: Dict = {"alpha": alpha, "config": _cfg(False), "shipped": {}, "map": []}

    for n in IPM.ORDER:
        dp = comp[n]["drug_params"]
        d_eff = float(dp.diffusivity) * tort
        k = float(dp.clearance) + float(dp.binding)
        out["shipped"][n] = {"d_eff_um2_s": d_eff, "k_loss_s_inv": k,
                             "lambda_um": math.sqrt(d_eff / k)}
    if verbose:
        print("\n" + "=" * 82)
        print("STAGE 2 — the near:far ratio is a function of ONE length")
        print("=" * 82)
        print(f"\n{_cfg(False)}\n")
        print(f"    {'compound':<15}{'D_eff':>10}{'k_loss':>12}{'lambda_um':>11}")
        for n in IPM.ORDER:
            s = out["shipped"][n]
            print(f"    {n:<15}{s['d_eff_um2_s']:>10.3f}"
                  f"{s['k_loss_s_inv']:>12.4e}{s['lambda_um']:>11.3f}")
        print("\n  ratio(lambda), by two INDEPENDENT routes to the same lambda:")
        print(f"    {'lambda_um':>10}{'via k':>12}{'via D':>12}{'disagreement':>14}")

    dp0 = comp[IPM.DOX]["drug_params"]
    d0 = float(dp0.diffusivity) * tort
    k0 = float(dp0.clearance) + float(dp0.binding)
    for lam in LAMBDA_GRID:
        via_k = _slab_ratio(dataclasses.replace(
            dp0, clearance=d0 / (lam * lam), binding=0.0), exchange=False)["ratio"]
        via_d = _slab_ratio(dataclasses.replace(
            dp0, diffusivity=k0 * lam * lam / tort), exchange=False)["ratio"]
        row = {"lambda_um": lam, "ratio_via_k": via_k, "ratio_via_d": via_d,
               "disagreement": max(via_k, via_d) / min(via_k, via_d)}
        out["map"].append(row)
        if verbose:
            print(f"    {lam:>10.3f}{via_k:>12.4f}{via_d:>12.4f}"
                  f"{row['disagreement']:>13.3f}x")

    lam_grid = np.array([r["lambda_um"] for r in out["map"]])
    r_grid = np.array([r["ratio_via_k"] for r in out["map"]])
    order = np.argsort(r_grid)

    def lam_for(ratio: float) -> float:
        return float(np.interp(ratio, r_grid[order], lam_grid[order]))

    out["implied_lambda"] = {}
    if verbose:
        print("\n  where the MEASUREMENT falls on the model's own coordinate:")
        print(f"    {'compound':<15}{'measured':>18}{'implied lambda':>16}"
              f"{'shipped lambda':>16}{'fold error':>12}")
    for n, (mean, sd) in KYLE_NEAR_FAR.items():
        shipped = out["shipped"][n]["lambda_um"]
        lo, hi = lam_for(mean + sd), lam_for(max(mean - sd, r_grid.min()))
        rec = {"measured_near_far": mean, "measured_sd": sd,
               "implied_lambda_um": lam_for(mean),
               "implied_lambda_um_band": [lo, hi],
               "shipped_lambda_um": shipped,
               "fold_error": shipped / lam_for(mean)}
        out["implied_lambda"][n] = rec
        if verbose:
            print(f"    {n:<15}{f'{mean:.0f} +- {sd:.0f}':>18}"
                  f"{rec['implied_lambda_um']:>13.2f} um"
                  f"{shipped:>13.2f} um{rec['fold_error']:>11.2f}x")
    if verbose:
        worst = max(v["fold_error"] for v in out["implied_lambda"].values())
        print(f"\n  -> the ~14x absolute-gradient failure is an error of at most "
              f"{worst:.2f}x in ONE length.")
        print("     Nothing is adopted here. lambda alone CANNOT say whether the")
        print("     error is in D_eff or in k; the two routes above are the only")
        print("     handle and they agree to within "
              f"{max(r['disagreement'] for r in out['map']):.2f}x.")
    return out


# --------------------------------------------------------------------------- #
# STAGE 3 -- the verdict test on real tissue
# --------------------------------------------------------------------------- #
def _tissue_cfg() -> str:
    return (f"10x Xenium GBM flagship section, dx = {TISSUE_DX_UM:g} um, "
            f"repair = {TISSUE_REPAIR}, dt = {TISSUE_DT_S:g} s, "
            f"{TISSUE_TOTAL_S / 3600:.0f} h window, convection "
            f"{'ON' if TISSUE_CONVECTION else 'OFF'}")


def stage_tissue(verbose: bool = True, sweep: bool = True,
                 drugs: Optional[Sequence[str]] = None) -> Dict:
    alpha = _assert_frozen()
    micro, _, _ = VP.build_arm_microenvironment(dx_um=TISSUE_DX_UM,
                                                repair=TISSUE_REPAIR)
    comp = IPM.compile_all(alpha, two_route=True)
    names = list(drugs) if drugs else list(IPM.ORDER)
    pde_kw = dict(dt=TISSUE_DT_S, total_time=TISSUE_TOTAL_S,
                  convection=TISSUE_CONVECTION, save_every=500)
    out: Dict = {"alpha": alpha, "config": _tissue_cfg(),
                 "compounds": {}, "ps_sweep": []}

    if verbose:
        print("\n" + "=" * 82)
        print("STAGE 3 — does the split change a PER-CELL verdict on real tissue?")
        print("=" * 82)
        print(f"\n{_tissue_cfg()}\n")

    def _run(dp, exchange: bool):
        pde = PDEConfig(transvascular_exchange=exchange, **pde_kw)
        return simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0),
                                 attribute=True)

    for n in names:
        dp = comp[n]["drug_params"]
        labels, rec = {}, {}
        for exchange in (False, True):
            t0 = time.time()
            res = _run(dp, exchange)
            counts = {k: int(v) for k, v in res.mechanism["counts"].items()}
            bd = res.mechanism["survival_breakdown"]
            n_tumor = int(res.mechanism["n_tumor"])
            labels[exchange] = np.asarray(res.mechanism["labels"])
            rec["ON" if exchange else "OFF"] = {
                "kill_fraction": float(res.kill_fraction),
                "n_tumor": n_tumor,
                "counts": counts,
                "unreached_fraction_of_tumor":
                    counts["SURVIVED_UNREACHED"] / n_tumor if n_tumor else None,
                "survival_breakdown": {k: float(v) for k, v in bd.items()},
                "dominant_limiter": max(bd, key=bd.get),
                "mean_auc": float(np.asarray(res.fields["auc"]).mean()),
                "solve_s": time.time() - t0,
            }
            if verbose:
                a = rec["ON" if exchange else "OFF"]
                print(f"  {n:<15} exchange {'ON ' if exchange else 'OFF'}  "
                      f"kill = {a['kill_fraction']:.6f}  unreached/tumour = "
                      f"{a['unreached_fraction_of_tumor']:.4f}  dominant = "
                      f"{a['dominant_limiter']:<9} ({a['solve_s']:.0f}s)")
        changed = labels[False] != labels[True]
        trans = Counter(zip(labels[False][changed].tolist(),
                            labels[True][changed].tolist()))
        rec["per_cell"] = {
            "n_cells": int(changed.size),
            "n_changed": int(changed.sum()),
            "fraction_changed": float(changed.mean()),
            "transitions": {f"{LABELS[u]}->{LABELS[v]}": int(c)
                            for (u, v), c in trans.items()},
        }
        rec["dominant_limiter_inverts"] = (
            rec["OFF"]["dominant_limiter"] != rec["ON"]["dominant_limiter"])
        out["compounds"][n] = rec
        if verbose:
            top = ", ".join(
                f"{LABELS[u]}->{LABELS[v]}: {c}"
                for (u, v), c in sorted(trans.items(), key=lambda kv: -kv[1])[:4])
            print(f"    -> {rec['per_cell']['n_changed']:,} of "
                  f"{rec['per_cell']['n_cells']:,} cells change label "
                  f"({100 * rec['per_cell']['fraction_changed']:.3f} %)"
                  + (f": {top}" if top else ""))
            print(f"    -> dominant limiter "
                  f"{rec['OFF']['dominant_limiter']} -> "
                  f"{rec['ON']['dominant_limiter']}"
                  + ("  *** INVERTS ***" if rec["dominant_limiter_inverts"] else ""))

    if sweep:
        if verbose:
            print(f"\n  ENTRY-SCALE SWEEP on {names[0]}, both laws. PS x c:")
            print(f"    {'c':>8}{'law':>6}{'kill':>10}{'mean AUC':>13}"
                  f"{'KILLED':>9}{'UNREACHED':>11}{'UNDERDOSED':>12}")
        dp = comp[names[0]]["drug_params"]
        for c in TISSUE_PS_SWEEP:
            for exchange in (False, True):
                res = _run(_scaled(dp, c), exchange)
                counts = {k: int(v) for k, v in res.mechanism["counts"].items()}
                row = {"c": c, "exchange": exchange, "drug": names[0],
                       "kill_fraction": float(res.kill_fraction),
                       "mean_auc": float(np.asarray(res.fields["auc"]).mean()),
                       "counts": counts}
                out["ps_sweep"].append(row)
                if verbose:
                    print(f"    {c:>8.0e}{'ON' if exchange else 'OFF':>6}"
                          f"{row['kill_fraction']:>10.6f}{row['mean_auc']:>13.4e}"
                          f"{counts['KILLED']:>9,}"
                          f"{counts['SURVIVED_UNREACHED']:>11,}"
                          f"{counts['SURVIVED_UNDERDOSED']:>12,}")
        off = [r for r in out["ps_sweep"] if not r["exchange"]]
        on = [r for r in out["ps_sweep"] if r["exchange"]]
        out["sweep_auc_span_off"] = (max(r["mean_auc"] for r in off)
                                     / min(r["mean_auc"] for r in off))
        out["sweep_auc_span_on"] = (max(r["mean_auc"] for r in on)
                                    / min(r["mean_auc"] for r in on))
        out["sweep_c_span"] = TISSUE_PS_SWEEP[-1] / TISSUE_PS_SWEEP[0]
        if verbose:
            print(f"\n    over a {out['sweep_c_span']:.0e} sweep of the entry "
                  f"scale, mean exposure spans {out['sweep_auc_span_off']:.4g}x "
                  f"under the one-way law and {out['sweep_auc_span_on']:.4g}x "
                  f"under the split.")
            print("    The one-way span EQUALS the sweep, exactly: elasticity 1.")
    return out


# --------------------------------------------------------------------------- #
# STAGE 4 -- honest scope
# --------------------------------------------------------------------------- #
def stage_scope(verbose: bool = True,
                radii_um: Sequence[float] = (2.0, 4.0, 8.0)) -> Dict:
    alpha = _assert_frozen()
    comp = IPM.compile_all(alpha, two_route=True)
    tort = 1.0 - ALPHA_TORTUOSITY * ECM
    out: Dict = {"alpha": alpha, "ps_scale_s_inv": PS_SCALE_S_INV,
                 "compounds": {}, "elasticity": []}

    if verbose:
        print("\n" + "=" * 82)
        print("STAGE 4 — honest scope: what parameter does the split trade for?")
        print("=" * 82)
        print(f"\n  PS [1/s] = P [um/s] x (S/V) [1/um], and S/V = 2*phi/R.")
        print(f"  The shipped law is PS = phi * bbb_penetrance * "
              f"{PS_SCALE_S_INV:g} 1/s, so phi cancels and the IMPLIED wall")
        print("  permeability is a per-compound constant that depends on a "
              "vessel radius\n  this repository does not measure.\n")
        print(f"    {'compound':<15}{'bbb_pen':>10}"
              + "".join(f"{'P @ R=' + f'{r:g}um':>14}" for r in radii_um))
    for n in IPM.ORDER:
        dp = comp[n]["drug_params"]
        rec = {"bbb_penetrance": float(dp.bbb_penetrance),
               "implied_permeability_um_s": {
                   str(r): implied_permeability_um_s(dp.bbb_penetrance, r,
                                                     PS_SCALE_S_INV)
                   for r in radii_um}}
        out["compounds"][n] = rec
        if verbose:
            print(f"    {n:<15}{rec['bbb_penetrance']:>10.5f}"
                  + "".join(f"{rec['implied_permeability_um_s'][str(r)]:>14.5g}"
                            for r in radii_um))
    if verbose:
        print("\n  These are NOT measurements and NOT predictions. They are what "
              "the shipped\n  law would have to be claiming if it were a "
              "permeability law. Every factor\n  is contested: bbb_penetrance is "
              "a dimensionless PARTITION coefficient,\n  the 1 s^-1 is undocumented"
              " (results/transvascular_exchange_revalidation.md\n  §6.2), and R is "
              "the caller's.\n")

    dp = comp[IPM.DOX]["drug_params"]
    k = float(dp.clearance) + float(dp.binding)
    d_eff = float(dp.diffusivity) * tort
    if verbose:
        print(f"  Entry elasticity d ln(exposure)/d ln(PS), doxorubicin, "
              f"k_loss = {k:.5g} 1/s:")
        print(f"    {'PS (1/s)':>12}{'one-way':>10}{'split':>10}"
              f"{'C_ss/Cp one-way':>18}{'C_ss/Cp split':>16}{'delta_um':>11}")
    for ps in (1e-6, 1e-4, k, 1e-2, 1e-1, 1.0, 1e2):
        row = {"ps_s_inv": ps,
               "elasticity_one_way": entry_elasticity(ps, k, split=False),
               "elasticity_split": entry_elasticity(ps, k, split=True),
               "plasma_multiple_one_way": one_way_plasma_multiple(ps, k),
               "plasma_fraction_split": saturation_fraction(ps, k),
               "boundary_layer_um_split": boundary_layer_um(d_eff, ps, k,
                                                            split=True)}
        out["elasticity"].append(row)
        if verbose:
            print(f"    {ps:>12.4g}{row['elasticity_one_way']:>10.4f}"
                  f"{row['elasticity_split']:>10.4f}"
                  f"{row['plasma_multiple_one_way']:>18.4g}"
                  f"{row['plasma_fraction_split']:>16.4f}"
                  f"{row['boundary_layer_um_split']:>11.2f}")
    ps_shipped = float(dp.bbb_penetrance) * PS_SCALE_S_INV
    out["regime_at_shipped_ps"] = entry_transport_regime(
        d_eff, ps_shipped, k, spacing_um=SLAB_THICKNESS_UM)
    out["boundary_layer_um_at_shipped_ps_split"] = boundary_layer_um(
        d_eff, ps_shipped, k, split=True)
    out["boundary_layer_um_one_way"] = boundary_layer_um(d_eff, ps_shipped, k,
                                                         split=False)
    if verbose:
        print(f"\n  At the shipped wall PS = {ps_shipped:.4e} 1/s (phi = 1), "
              f"the perivascular decay\n  length is "
              f"{out['boundary_layer_um_one_way']:.2f} um under the one-way law "
              f"and {out['boundary_layer_um_at_shipped_ps_split']:.2f} um under "
              f"the split;\n  the Damkohler number at a {SLAB_THICKNESS_UM:.0f} um"
              f" spacing is "
              f"{out['regime_at_shipped_ps']['damkohler']:.3g} -> "
              f"{out['regime_at_shipped_ps']['regime']}.")
        delta = out["boundary_layer_um_at_shipped_ps_split"]
        print(f"  The split's {delta:.2f} um is "
              f"{delta / TISSUE_DX_UM:.1f}x the flagship section's "
              f"dx = {TISSUE_DX_UM:g} um, so it is resolved on that grid at the "
              f"shipped PS.\n  It stops being resolved above PS ~ "
              f"{d_eff / (TISSUE_DX_UM ** 2) - k:.3g} 1/s. The discretisation "
              f"question is owned by\n  experiments/boundary_layer_resolution.py, "
              "not by this file.")
    return out


# --------------------------------------------------------------------------- #
def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def main(argv: Optional[List[str]] = None) -> Dict:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", default="all",
                    choices=("all", "entanglement", "lambda", "tissue", "scope"))
    ap.add_argument("--no-sweep", action="store_true",
                    help="stage 'tissue': skip the entry-scale sweep (faster)")
    ap.add_argument("--drug", action="append", default=None,
                    help="stage 'tissue': restrict to these compounds")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)

    t0 = time.time()
    res: Dict = {"alpha_asserted": ALPHA_LIVE,
                 "shipped_defaults": {
                     "PDEConfig.transvascular_exchange":
                         PDEConfig().transvascular_exchange,
                     "transvascular.SPLIT_DEFAULT": SPLIT_DEFAULT,
                     "drugs.PS_SCALE_S_INV": PS_SCALE_S_INV,
                     "PKConfig": dataclasses.asdict(PKConfig())}}
    if args.stage in ("all", "entanglement"):
        res["entanglement"] = stage_entanglement()
    if args.stage in ("all", "lambda"):
        res["lambda"] = stage_lambda()
    if args.stage in ("all", "tissue"):
        res["tissue"] = stage_tissue(sweep=not args.no_sweep, drugs=args.drug)
    if args.stage in ("all", "scope"):
        res["scope"] = stage_scope()
    res["runtime_s"] = time.time() - t0
    _assert_frozen()

    if not args.no_write and args.stage == "all":
        with open(_DATA_PATH, "w") as fh:
            json.dump(_jsonable(res), fh, indent=1)
        print(f"\nwrote {os.path.relpath(_DATA_PATH, _ROOT)} "
              f"({res['runtime_s']:.0f}s)")
    return res


if __name__ == "__main__":
    main()
