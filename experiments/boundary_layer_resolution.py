"""Is there a length scale the shipped dx = 12 um grid cannot resolve, and where?

The hypothesis under test
-------------------------
``results/transvascular_exchange_revalidation.md`` section 5 records a
perivascular boundary layer introduced by the gradient-driven delivery law::

    delta = sqrt(D_eff / (PS + clearance + binding))

and states that it "collapses to roughly 15-17 um for bbb_penetrance >~ 0.6".
That same section is ANNOTATED AS REFUTED on its directional consequence: it
predicted finer grids would push exchange-ON values DOWN, that was measured, and
the values went UP by 1.3-2.6%. So the mechanism may be real while the direction
inferred from it is not. This script measures rather than argues.

Three separate questions, three separate stages
-----------------------------------------------
STAGE 1 -- the length-scale census (no PDE, seconds).
    For the REAL GBM section and the REAL breast section, at the shipped
    dx = 12 um, compute every length scale in the problem from the actual
    compiled DrugParams and the actual reconstructed fields:

      * ``delta``  = sqrt(D_eff/(PS + clearance + binding))  -- exchange ON,
        evaluated ONLY where PS > 0 (PS is exactly zero away from a vessel under
        the shipped discrete vascular field, so delta is not defined there).
      * ``L``      = sqrt(D_eff/(clearance + binding))       -- exchange OFF.
        This is the length scale the SHIPPED verdicts actually run on; see the
        propagation note below.
      * ``R``      = VESSEL_RADIUS_UM = 7.5 um, the support radius of the
        source itself, which is a GEOMETRIC scale and is drug-independent.

    Reported as cells-per-layer at dx = 8/12/16/24/30.

STAGE 2 -- controlled numerics (a 1-D two-region slab with a CLOSED-FORM steady
    state). The shipped solver is run on that slab at a range of dx and a range
    of delta, and its discretisation error is measured against the exact
    solution. Two variants:

      * interface ALIGNED with a cell face at every dx -> isolates the cost of
        under-resolving the exponential layer;
      * interface OFFSET by half a cell -> adds the geometric quantisation
        error that the real 7.5 um vessel discs suffer.

    Output: error vs delta/dx, which is a property of the SCHEME and transfers
    to every future run.

STAGE 3 -- does delta predict which drugs move with dx? The per-drug movement
    numbers already measured in this repo are read from the results files (they
    are not recomputed here; two sibling agents own those sweeps) and correlated
    against the stage-1 length scales.

THE KNOWN TRAP (tmesim/validation/propagation.py)
--------------------------------------------------
``experiments/kill_calibration_breast.py`` line 327 and
``tmesim/validation/kill_calibration.py`` both go through
``kill_calibration.unit_transport``, which calls ``simulate_drug_transport``
DIRECTLY with no ``exchange_k``. ``kill_calibration.py`` is DECLARED ``INERT``
in ``propagation.DECLARATIONS`` for exactly this reason. Consequence for this
brief: the shipped breast and GBM 3-D:2-D verdicts -- including the docetaxel
verdict that moves with dx -- are EXCHANGE-OFF numbers, so ``delta`` is not even
an active length scale on that path. That is checked mechanically in stage 0
rather than asserted.

Stage 2 injects ``exchange_k`` by hand, which is the pattern
``propagation.WIRING_POINTS`` documents as correct for a validation experiment
whose subject path is INERT. It asserts the two arms are NOT bit-identical
before any comparison is drawn.

Nothing under ``tmesim/`` is modified or re-implemented. The transport operator
exercised in stage 2 is the shipped ``simulate_drug_transport``.

Run:
    ./.venv/bin/python experiments/boundary_layer_resolution.py --stage all
    ./.venv/bin/python experiments/boundary_layer_resolution.py --stage 1
Writes:
    results/sweep/blayer_lengthscales.json
    results/sweep/blayer_controlled.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.core import (  # noqa: E402
    CellPopulation,
    CellState,
    CellType,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    TissueFields,
)
from tmesim.physics.diffusion import (  # noqa: E402
    ALPHA_TORTUOSITY,
    simulate_drug_transport,
)
from tmesim.reconstruct.vessel_network import VESSEL_RADIUS_UM  # noqa: E402
from tmesim.treatment.drugs import PS_SCALE_S_INV, exchange_coefficient  # noqa: E402

OUT_DIR = os.path.join(_ROOT, "results", "sweep")

#: The grids the brief asks for, plus the shipped one.
DX_GRID_UM: Tuple[float, ...] = (8.0, 12.0, 16.0, 24.0, 30.0)
SHIPPED_DX_UM = 12.0

BREAST_SAMPLE_PATH = os.path.join(_ROOT, "data_cache", "xenium_breast_janesick",
                                  "rep1")


def log(msg: str) -> None:
    print(f"[blayer {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# STAGE 0 -- the propagation trap, checked mechanically
# --------------------------------------------------------------------------- #
def stage0_trap_check() -> Dict:
    """Confirm from the repo's own declaration table which paths are INERT.

    This is a READ of ``propagation.DECLARATIONS`` plus a static AST scan, not a
    claim of this script's own.
    """
    from tmesim.validation.propagation import (
        Category,
        declaration_for,
        scan_call_sites,
        wiring_points,
    )

    out: Dict = {"wiring_points_all": sorted(wiring_points(False)),
                 "wiring_points_production": sorted(wiring_points(True))}

    kc = declaration_for("tmesim/validation/kill_calibration.py")
    out["kill_calibration_declared"] = {"category": kc.category.value,
                                        "reason": kc.reason}

    # The breast panel script is not itself declared; scan it directly.
    for rel in ("experiments/kill_calibration_breast.py",
                "tmesim/validation/kill_calibration.py"):
        cs = scan_call_sites(os.path.join(_ROOT, rel))
        out[rel] = {
            "direct_solver_calls": list(cs.direct_calls),
            "direct_calls_with_exchange_k": list(cs.direct_calls_with_exchange_k),
            "bare_direct_calls": list(cs.bare_direct_calls),
            "engine_calls": list(cs.engine_calls),
            "static_prediction": cs.predicted.value,
        }
    out["conclusion"] = (
        "kill_calibration.unit_transport is the transport call behind BOTH the "
        "GBM and breast 3-D:2-D panels and it passes no exchange_k, so the "
        "shipped verdicts are exchange-OFF and `delta` is not an active length "
        "scale on that path."
    )
    _ = Category  # imported for the reader; category strings above are its values
    return out


# --------------------------------------------------------------------------- #
# STAGE 1 -- length-scale census on the two real sections
# --------------------------------------------------------------------------- #
def _compile_gbm_panel(alpha: float) -> Dict[str, DrugParams]:
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.library import GBM_DRUGS, get_molecule
    from tmesim.validation.calibration import calibrated_drugparams

    out: Dict[str, DrugParams] = {}
    for name in GBM_DRUGS:
        try:
            dp = calibrated_drugparams(name, alpha)
            dp = apply_barrier_to_drugparams(dp, get_molecule(name))
        except Exception as exc:  # a molecule with no potency entry
            log(f"  GBM: skipping {name} ({type(exc).__name__}: {exc})")
            continue
        out[name] = dp
    return out


def _compile_breast_panel(alpha: float) -> Dict[str, DrugParams]:
    from experiments.kill_calibration_breast import DRUGS as BREAST_DRUGS
    from experiments.kill_calibration_breast import compile_breast_drug

    out: Dict[str, DrugParams] = {}
    for name in BREAST_DRUGS:
        dp, _mol, _ic50 = compile_breast_drug(name, alpha)
        out[name] = dp
    return out


def _field_stats(vals: np.ndarray) -> Dict:
    if vals.size == 0:
        return {"n": 0}
    return {
        "n": int(vals.size),
        "min": float(np.min(vals)),
        "p05": float(np.percentile(vals, 5)),
        "median": float(np.median(vals)),
        "mean": float(np.mean(vals)),
        "p95": float(np.percentile(vals, 95)),
        "max": float(np.max(vals)),
    }


def length_scales_for_section(tissue_name: str, micro, drugs: Dict[str, DrugParams]
                              ) -> Dict:
    """Every length scale in the problem, from the ACTUAL reconstructed fields."""
    ecm = np.asarray(micro.fields.ecm_density, dtype=np.float64)
    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)

    vasc_t = vasc[tissue]
    perfused = tissue & (vasc > 0.0)

    section = {
        "tissue": tissue_name,
        "dx_um": float(micro.domain.dx),
        "grid": [int(micro.domain.ny), int(micro.domain.nx)],
        "n_tissue_cells": int(tissue.sum()),
        "n_perfused_cells": int(perfused.sum()),
        "perfused_area_fraction": float(perfused.sum()) / max(1, int(tissue.sum())),
        "vascular_density_over_tissue": _field_stats(vasc_t),
        "vascular_density_where_positive": _field_stats(vasc[perfused]),
        "ecm_density_over_tissue": _field_stats(ecm[tissue]),
        "ps_scale_s_inv": float(PS_SCALE_S_INV),
        "vessel_radius_um": float(VESSEL_RADIUS_UM),
        "drugs": [],
    }

    for name in sorted(drugs):
        dp = drugs[name]
        # exchange_coefficient IS delivery_source at unit plasma: the shipped PS.
        ps = np.asarray(exchange_coefficient(micro, dp), dtype=np.float64)
        d_raw = float(dp.diffusivity) * (1.0 - ALPHA_TORTUOSITY * np.clip(ecm, 0, 1))
        floor = max(1e-6, 1e-3 * float(dp.diffusivity))
        d_eff = np.maximum(d_raw, floor)

        k_off = float(dp.clearance) + float(dp.binding)          # scalar, 1/s
        k_on = k_off + ps                                        # field, 1/s

        # delta is only defined where PS > 0 -- the shipped discrete vascular
        # field is EXACTLY zero away from a vessel, so "the boundary layer" has
        # no meaning outside the perfused voxels.
        delta = np.sqrt(d_eff[perfused] / k_on[perfused])
        L_field = np.sqrt(d_eff[tissue] / k_off)

        row = {
            "drug": name,
            "bbb_penetrance": float(dp.bbb_penetrance),
            "diffusivity_um2_s": float(dp.diffusivity),
            "clearance_per_s": float(dp.clearance),
            "binding_per_s": float(dp.binding),
            "k_off_per_s": k_off,
            "PS_max_per_s": float(ps.max()),
            "PS_median_where_positive": (float(np.median(ps[perfused]))
                                         if perfused.any() else 0.0),
            "PS_over_k_off_max": float(ps.max()) / k_off,
            "D_eff_median_um2_s": float(np.median(d_eff[tissue])),
            "delta_um": _field_stats(delta),
            "L_um": _field_stats(L_field),
        }
        # Cells-per-layer at each dx, using the MEDIAN delta / L over the
        # relevant support (a distribution reduced to one number only for the
        # table; the full distribution is above).
        dmed = float(np.median(delta)) if delta.size else float("nan")
        lmed = float(np.median(L_field)) if L_field.size else float("nan")
        row["cells_per_layer_delta"] = {f"dx{int(d)}": dmed / d for d in DX_GRID_UM}
        row["cells_per_layer_L"] = {f"dx{int(d)}": lmed / d for d in DX_GRID_UM}
        row["cells_per_vessel_radius"] = {
            f"dx{int(d)}": float(VESSEL_RADIUS_UM) / d for d in DX_GRID_UM}
        section["drugs"].append(row)

    return section


def stage1(dx_um: float = SHIPPED_DX_UM, sections: Sequence[str] = ("gbm", "breast")
           ) -> Dict:
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.validation.calibration import calibrate_alpha

    alpha = float(calibrate_alpha())
    log(f"stage 1: alpha = {alpha:.6f}")
    out: Dict = {"alpha": alpha, "dx_um": float(dx_um),
                 "dx_grid_um": list(DX_GRID_UM), "sections": []}

    if "gbm" in sections:
        log("stage 1: reconstructing REAL Xenium GBM section ...")
        t = time.time()
        micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
        log(f"  GBM reconstructed in {time.time()-t:.1f} s")
        drugs = _compile_gbm_panel(alpha)
        log(f"  GBM panel: {len(drugs)} drugs compiled")
        out["sections"].append(length_scales_for_section("gbm", micro, drugs))

    if "breast" in sections:
        log("stage 1: reconstructing REAL Xenium breast section ...")
        t = time.time()
        micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH,
                                 dx_um=float(dx_um))
        log(f"  breast reconstructed in {time.time()-t:.1f} s")
        drugs = _compile_breast_panel(alpha)
        log(f"  breast panel: {len(drugs)} drugs compiled")
        out["sections"].append(length_scales_for_section("breast", micro, drugs))

    return out


def stage1b_source_vs_dx(dx_list: Sequence[float] = DX_GRID_UM,
                         sections: Sequence[str] = ("gbm", "breast")) -> Dict:
    """How much of a dx sweep on REAL tissue is a change to the PHYSICS?

    The shipped vascular field is DISCRETE: a grid cell is a source iff its
    CENTRE lies within ``VESSEL_RADIUS_UM = 7.5 um`` of a segmented endothelial
    cell. A 7.5 um radius disc has area 177 um^2; one dx = 12 um cell has area
    144 um^2 and one dx = 30 um cell has 900 um^2. So changing dx changes the
    painted vessel AREA and hence the total delivered mass -- a dx sweep on real
    tissue is therefore NOT a pure convergence study, and any interpretation of
    one has to net this out. Measured here rather than assumed.
    """
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    out: Dict = {"vessel_radius_um": float(VESSEL_RADIUS_UM),
                 "disc_area_um2": math.pi * VESSEL_RADIUS_UM ** 2,
                 "sections": {}}
    for tissue_name in sections:
        rows = []
        for dx in dx_list:
            t0 = time.time()
            kw = {} if tissue_name == "gbm" else {"sample": BREAST_SAMPLE_PATH}
            micro = reconstruct_real(source="xenium", dx_um=float(dx), **kw)
            vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
            tis = np.asarray(micro.domain.tissue_mask, dtype=bool)
            cell_area = dx * dx
            rows.append({
                "dx_um": float(dx),
                "grid": [int(micro.domain.ny), int(micro.domain.nx)],
                "n_tissue_cells": int(tis.sum()),
                "tissue_area_um2": float(tis.sum()) * cell_area,
                "n_source_cells": int((tis & (vasc > 0)).sum()),
                "source_area_um2": float((tis & (vasc > 0)).sum()) * cell_area,
                "source_area_fraction": float((tis & (vasc > 0)).sum())
                                        / max(1, int(tis.sum())),
                "source_integral_um2": float(vasc[tis].sum()) * cell_area,
                "seconds": time.time() - t0,
            })
            log(f"  {tissue_name} dx={dx:g}: source area frac "
                f"{rows[-1]['source_area_fraction']:.5f}, "
                f"integral {rows[-1]['source_integral_um2']:.4g} um^2 "
                f"({rows[-1]['seconds']:.1f}s)")
        ref = [r for r in rows if r["dx_um"] == SHIPPED_DX_UM][0]
        for r in rows:
            r["source_integral_rel_to_dx12"] = (
                r["source_integral_um2"] / ref["source_integral_um2"])
        out["sections"][tissue_name] = rows
    return out


# --------------------------------------------------------------------------- #
# STAGE 2 -- the controlled problem: a 1-D two-region slab, exact steady state
# --------------------------------------------------------------------------- #
#
#   region A   0 <= x < a      k_in  = PS + kappa      source = PS * Cp
#   region B   a <= x <= Lx    k_out = kappa           source = 0
#   no-flux at x = 0 and x = Lx, constant D.
#
# Steady state of  D C'' = k C - S :
#   A:  C = Ceq + A cosh(x/delta),         Ceq = PS*Cp/(PS+kappa), delta = sqrt(D/(PS+kappa))
#   B:  C = B cosh((Lx-x)/Lout),           Lout = sqrt(D/kappa)
# matched by continuity of C and of D C' at x = a. Both branches already satisfy
# their own no-flux end condition by construction.
#
def slab_exact(x: np.ndarray, *, d: float, k_in: float, k_out: float,
               s_in: float, a: float, lx: float) -> np.ndarray:
    """Exact steady state of ``D C'' = k C - S`` on the two-region slab.

    ``k_in``/``s_in`` apply on ``[0, a)``; ``k_out`` and zero source on
    ``[a, Lx]``. Written generally in ``k_in`` so the SAME closed form serves
    the exchange-ON case (``k_in = PS + kappa``) and the exchange-OFF one-way
    leak (``k_in = kappa``, same source): the two differ only in ``k_in``.
    """
    lam_a = math.sqrt(d / k_in)
    lam_b = math.sqrt(d / k_out)
    ceq = s_in / k_in

    # A cosh(x/lam_a) + Ceq  on [0,a);  B cosh((Lx-x)/lam_b) on [a,Lx]
    # continuity of C and of C' at a (D is constant, so flux continuity is C').
    ca, sa = math.cosh(a / lam_a), math.sinh(a / lam_a)
    cb, sb = math.cosh((lx - a) / lam_b), math.sinh((lx - a) / lam_b)
    # [ca, -cb; sa/lam_a, sb/lam_b] [A; B] = [-Ceq; 0]   (Cramer)
    det = ca * (sb / lam_b) + cb * (sa / lam_a)
    aa = (-ceq * (sb / lam_b)) / det
    bb = (ceq * (sa / lam_a)) / det

    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    inA = x < a
    out[inA] = ceq + aa * np.cosh(x[inA] / lam_a)
    out[~inA] = bb * np.cosh((lx - x[~inA]) / lam_b)
    return out


def _simpson(y: np.ndarray, h: float) -> float:
    n = y.size
    w = np.ones(n)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return float((h / 3.0) * np.dot(w, y))


def slab_exact_reference(*, d: float, k_in: float, k_out: float, s_in: float,
                         a: float, lx: float, dx: float,
                         sub: int = 512) -> Dict:
    """Exact domain mean AND exact per-cell averages on the ``dx`` grid.

    The numerical scheme is finite-volume, so its unknowns approximate CELL
    AVERAGES. Comparing them against point samples of the exact solution would
    inject an O(dx^2) error that belongs to the comparison, not to the scheme.
    Every cell average below is a Simpson quadrature of the CLOSED FORM over
    that cell with ``sub`` sub-intervals, i.e. exact to round-off.
    """
    nx = int(round(lx / dx))
    # (nx, sub+1) sub-sample grid: one Simpson panel set per cell.
    t = np.linspace(0.0, 1.0, sub + 1)
    xs = (np.arange(nx)[:, None] + t[None, :]) * dx
    vals = slab_exact(xs.ravel(), d=d, k_in=k_in, k_out=k_out, s_in=s_in,
                      a=a, lx=lx).reshape(nx, sub + 1)
    w = np.ones(sub + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    h = dx / sub
    # (vals @ w) would be identical but routes through BLAS, which raises
    # spurious FP-flag RuntimeWarnings on this build; the explicit reduction is
    # bit-stable and quiet. Verified insensitive to ``sub`` from 64 to 16384
    # (mean agrees to 1e-16), so the quadrature is exact for this purpose.
    cell_avg = (h / 3.0) * (vals * w[None, :]).sum(axis=1) / dx
    return {"cell_avg": cell_avg, "mean": float(cell_avg.mean()),
            "peak_cell": float(cell_avg[0])}


def slab_exact_penetration(*, d: float, k_in: float, k_out: float, s_in: float,
                           a: float, lx: float, n: int = 400_001) -> float:
    """The CONTINUUM limit of ``TransportResult.penetration_index``.

    The shipped index is ``mean(AUC over tissue) / percentile(AUC, 99)``. Its
    continuum analogue is the length-weighted mean of C over the domain divided
    by the length-weighted 99th percentile of C -- both computed here from the
    closed form on a 4e5-point grid, which is the dx -> 0 limit the discrete
    index is converging to. At steady state AUC is proportional to C, and the
    index is a RATIO, so the proportionality cancels exactly.
    """
    x = np.linspace(0.0, lx, n)
    c = slab_exact(x, d=d, k_in=k_in, k_out=k_out, s_in=s_in, a=a, lx=lx)
    return float(np.clip(c.mean() / np.percentile(c, 99.0), 0.0, 1.0))


def _slab_micro(nx: int, dx: float) -> Microenvironment:
    """A 1-cell-tall strip domain. ecm = 0 so D_eff == drug.diffusivity exactly."""
    domain = Domain(nx=int(nx), ny=1, dx=float(dx))
    fields = TissueFields(ecm_density=domain.zeros(),
                          vascular_density=domain.zeros())
    cells = CellPopulation(
        x=np.array([dx * 0.5]), y=np.array([dx * 0.5]),
        cell_type=np.array([int(CellType.TUMOR)]),
        state=np.array([int(CellState.VIABLE)]),
        sensitivity=np.array([1.0]))
    return Microenvironment(domain=domain, fields=fields, cells=cells)


def slab_numeric(*, d: float, ps: float, kappa: float, cp: float, a: float,
                 lx: float, dx: float, exchange: bool,
                 t_mult: float = 40.0) -> Dict:
    """Run the SHIPPED solver on the slab and return its steady state.

    ``exchange=True`` passes ``exchange_k = PS`` (the honest gradient-driven law,
    the thing whose boundary layer is under test). ``exchange=False`` is the
    one-way leak with the SAME source, kept so the two can be compared.

    Integrated to ``t_mult`` times the slowest relaxation time so the result is
    the discrete STEADY state (checked, and the residual drift is reported).
    """
    nx = int(round(lx / dx))
    if abs(nx * dx - lx) > 1e-9:
        raise ValueError(f"Lx={lx} is not an integer number of dx={dx} cells")
    micro = _slab_micro(nx, dx)
    xc = (np.arange(nx) + 0.5) * dx           # cell centres

    # Cell-averaged region indicator: the fraction of each cell lying in region A.
    left = np.arange(nx) * dx
    right = left + dx
    frac_a = np.clip((np.minimum(right, a) - left) / dx, 0.0, 1.0)

    drug = DrugParams(name="slab", diffusivity=float(d), clearance=float(kappa),
                      binding=0.0, bbb_penetrance=1.0)
    source = (ps * cp * frac_a).reshape(1, nx)
    ex = (ps * frac_a).reshape(1, nx) if exchange else None

    # relaxation time of the slowest mode ~ 1/kappa (plus diffusive time)
    tau = max(1.0 / kappa, lx * lx / d)
    total = t_mult * tau
    dt = min(tau / 50.0, total / 200.0)
    n_steps = int(round(total / dt))
    if n_steps > 200_000:
        dt = total / 200_000.0
    pde = PDEConfig(dt=float(dt), total_time=float(total), convection=False,
                    save_every=max(1, int(round(total / dt)) // 4))

    tr = simulate_drug_transport(micro, drug, source, pde, exchange_k=ex)
    c = np.asarray(tr.final_concentration, dtype=np.float64).ravel()
    hist = np.asarray(tr.mean_concentration_history, dtype=np.float64)
    drift = (abs(hist[-1] - hist[-2]) / max(abs(hist[-1]), 1e-300)
             if hist.size >= 2 else float("nan"))
    peak99 = float(np.percentile(c, 99.0))
    pen = float(np.clip(c.mean() / peak99, 0.0, 1.0)) if peak99 > 0 else 0.0
    return {"x": xc, "c": c, "mean": float(c.mean()), "peak": float(c[0]),
            "penetration_index": pen,
            "steady_drift": float(drift), "n_steps": int(round(total / dt)),
            "dt_s": float(dt), "total_time_s": float(total)}


#: Slab geometries. ``a_um`` is the half-width of the source region.
#:
#: ``wide_aligned``  a = 48 um, an integer number of cells at EVERY dx tested,
#:     so the region interface sits exactly on a cell face and the ONLY error
#:     measured is the failure to resolve the exponential layer. This is the
#:     clean error-vs-(delta/dx) curve the brief asks for.
#: ``wide_offset``   a = 54 um, deliberately NOT a cell multiple at dx = 12/16/
#:     24, so the geometric quantisation error of the interface is added on top.
#: ``vessel``        a = 7.5 um = ``VESSEL_RADIUS_UM``, the ACTUAL half-width of
#:     the shipped discrete vascular source. This is the geometry the real runs
#:     have, and it is deliberately included because the stage-1 census shows
#:     delta > 2*a for every drug in both panels -- i.e. on the real problem the
#:     source's own support, not the boundary layer, may be the binding scale.
SLAB_GEOMETRIES: Tuple[Tuple[str, float], ...] = (
    ("wide_aligned", 48.0),
    ("wide_offset", 54.0),
    ("vessel", float(VESSEL_RADIUS_UM)),
)


def stage2(delta_targets_um: Sequence[float] = (5.0, 10.0, 15.0, 20.0, 30.0,
                                                45.0),
           dx_list: Sequence[float] = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0,
                                       16.0, 24.0, 30.0),
           lx: float = 480.0, cp: float = 1.0) -> Dict:
    """Error vs (delta/dx) for the shipped scheme on the closed-form slab.

    ``kappa`` is fixed so the OUTER length ``L = sqrt(D/kappa) = 64.92 um`` (the
    measured doxorubicin e-fold anchor) is the same in every series, and ``PS``
    is set to hit each target ``delta``, so within a series the only thing that
    changes is ``delta/dx``.

    Both delivery laws are run at every point, each scored against ITS OWN exact
    solution, so "does exchange-ON need a finer grid than exchange-OFF" is
    answered by a difference of two measured error curves rather than by
    argument.
    """
    d = 140.41            # um^2/s -- the shipping D_tissue (dtissue_sensitivity.md)
    l_out = 64.92         # um -- the measured doxorubicin e-fold anchor
    kappa = d / (l_out ** 2)

    out: Dict = {"D_um2_s": d, "L_out_um": l_out, "kappa_per_s": kappa,
                 "Lx_um": lx, "Cp_uM": cp,
                 "geometries": {t: a for t, a in SLAB_GEOMETRIES},
                 "series": []}

    for tag, a_um in SLAB_GEOMETRIES:
        for dtarget in delta_targets_um:
            ps = d / (dtarget ** 2) - kappa
            if ps <= 0:
                continue
            series = {"tag": tag, "a_um": a_um, "delta_um": dtarget,
                      "L_out_um": l_out, "PS_per_s": ps,
                      "PS_over_kappa": ps / kappa, "points": []}
            for dx in dx_list:
                if abs(round(lx / dx) * dx - lx) > 1e-9:
                    continue
                t0 = time.time()
                rec: Dict = {
                    "dx_um": dx,
                    "delta_over_dx": dtarget / dx,
                    "L_over_dx": l_out / dx,
                    "a_over_dx": a_um / dx,
                    "interface_aligned": abs(round(a_um / dx) * dx - a_um) < 1e-9,
                }
                for law, exchange in (("on", True), ("off", False)):
                    k_in = (ps + kappa) if exchange else kappa
                    ref = slab_exact_reference(d=d, k_in=k_in, k_out=kappa,
                                               s_in=ps * cp, a=a_um, lx=lx,
                                               dx=dx)
                    num = slab_numeric(d=d, ps=ps, kappa=kappa, cp=cp, a=a_um,
                                       lx=lx, dx=dx, exchange=exchange)
                    prof_err = float(
                        np.sqrt(np.mean((num["c"] - ref["cell_avg"]) ** 2))
                        / np.sqrt(np.mean(ref["cell_avg"] ** 2)))
                    rec[f"{law}_mean"] = num["mean"]
                    rec[f"{law}_exact_mean"] = ref["mean"]
                    rec[f"{law}_rel_err_mean"] = (
                        (num["mean"] - ref["mean"]) / ref["mean"])
                    rec[f"{law}_rel_err_peakcell"] = (
                        (num["peak"] - ref["peak_cell"]) / ref["peak_cell"])
                    rec[f"{law}_rel_l2_profile"] = prof_err
                    rec[f"{law}_steady_drift"] = num["steady_drift"]
                    # The shipped penetration_index and its continuum limit.
                    pen_exact = slab_exact_penetration(
                        d=d, k_in=k_in, k_out=kappa, s_in=ps * cp, a=a_um,
                        lx=lx)
                    rec[f"{law}_pen_index"] = num["penetration_index"]
                    rec[f"{law}_pen_index_exact"] = pen_exact
                    rec[f"{law}_pen_index_rel_err"] = (
                        (num["penetration_index"] - pen_exact) / pen_exact)
                    # dt is reported per point: the coordinator's CFL note.
                    # convection is OFF on this slab (velocity is None), so the
                    # upwind CFL does not apply and dt is chosen purely to reach
                    # steady state; ``*_steady_drift`` proves it got there, and
                    # backward Euler's fixed point is dt-EXACT, so every error
                    # below is purely SPATIAL.
                    rec[f"{law}_dt_s"] = num["dt_s"]
                    rec[f"{law}_n_steps"] = num["n_steps"]
                    rec[f"{law}_convection"] = False
                rec["seconds"] = time.time() - t0
                series["points"].append(rec)
            out["series"].append(series)
            log(f"  stage2 {tag:13s} delta={dtarget:5g} um  ON err: "
                + " ".join(f"dx{p['dx_um']:g}:{p['on_rel_err_mean']*100:+.2f}%"
                           for p in series["points"]
                           if p["dx_um"] in (4.0, 8.0, 12.0, 24.0, 30.0)))

    return out


def stage2_trap_assert() -> Dict:
    """Prove the two arms are NOT bit-identical before any comparison is drawn.

    ``simulate_drug_transport`` is called DIRECTLY here, so a ``PDEConfig`` flag
    would be a silent no-op; this script injects ``exchange_k`` itself, which is
    the pattern ``propagation.WIRING_POINTS`` documents for an INERT subject.
    The assertion below is the guard that pattern requires.
    """
    d, l_out = 140.41, 64.92
    kappa = d / l_out ** 2
    ps = d / (15.0 ** 2) - kappa
    off = slab_numeric(d=d, ps=ps, kappa=kappa, cp=1.0, a=48.0, lx=480.0,
                       dx=12.0, exchange=False)
    on = slab_numeric(d=d, ps=ps, kappa=kappa, cp=1.0, a=48.0, lx=480.0,
                      dx=12.0, exchange=True)
    identical = bool(np.array_equal(off["c"], on["c"]))
    assert not identical, (
        "exchange ON and OFF produced BIT-IDENTICAL fields -- exchange_k did not "
        "reach the solver and every stage-2 number would be a comparison of a "
        "run to itself. Refusing to continue.")
    return {"exchange_off_mean": off["mean"], "exchange_on_mean": on["mean"],
            "bit_identical": identical,
            "ratio_on_over_off": on["mean"] / off["mean"]}


# --------------------------------------------------------------------------- #
# STAGE 3 -- what does a finer grid actually COST, as a (dx, dt) pair?
# --------------------------------------------------------------------------- #
#: (dx, dt) pairs to price. dt is NOT free to choose: ``results/vessel_promotion.md``
#: records that the discrete vascular source steepens the IFP gradient enough
#: that the explicit upwind advection CFL gives dt_max = 20.5 / 15.4 / 9.5 / 6.9 s
#: at dx = 15 / 12 / 8 / 6, i.e. the shipped dt = 10 s VIOLATES CFL at dx <= 8.
#: A recommendation of "use dx = 8" therefore silently also costs a halved dt, so
#: the cost is measured as a PAIR and the CFL number is recomputed here, not
#: quoted.
COST_POINTS: Tuple[Tuple[float, float], ...] = (
    (24.0, 10.0),
    (16.0, 10.0),
    (12.0, 10.0),   # shipped
    (8.0, 5.0),
    (6.0, 3.0),
)


def stage3_cost(tissue_name: str = "breast",
                points: Sequence[Tuple[float, float]] = COST_POINTS) -> Dict:
    """Measured wall-clock and memory for ONE transport solve at each (dx, dt).

    Uses the SHIPPED path: ``reconstruct_real`` -> ``prepare_tissue`` (the
    engine's own oxygen/state/IFP prologue) -> ``unit_transport``. The measured
    max interstitial speed gives the true CFL limit at each dx, so "is dt = 10 s
    legal here" is answered from this run rather than quoted.
    """
    import tracemalloc

    from tmesim.core import PDConfig, PKConfig, SimConfig
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.validation.calibration import calibrate_alpha
    from tmesim.validation.kill_calibration import prepare_tissue, unit_transport

    alpha = float(calibrate_alpha())
    if tissue_name == "breast":
        from experiments.kill_calibration_breast import compile_breast_drug
        drug = compile_breast_drug("docetaxel", alpha)[0]
        rec_kw = {"sample": BREAST_SAMPLE_PATH}
    else:
        from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
        from tmesim.molecular.library import get_molecule
        from tmesim.validation.calibration import calibrated_drugparams
        drug = apply_barrier_to_drugparams(
            calibrated_drugparams("doxorubicin", alpha),
            get_molecule("doxorubicin"))
        rec_kw = {}

    out: Dict = {"tissue": tissue_name, "drug": drug.name, "alpha": alpha,
                 "window_s": 3600.0 * 6, "points": []}

    for dx, dt in points:
        t0 = time.time()
        micro = reconstruct_real(source="xenium", dx_um=float(dx), **rec_kw)
        t_recon = time.time() - t0

        pde = PDEConfig(dt=float(dt), total_time=3600.0 * 6, convection=True,
                        save_every=10)
        cfg = SimConfig(pde=pde, pk=PKConfig(), pd=PDConfig(), seed=0)
        t0 = time.time()
        work, _labeled, velocity = prepare_tissue(micro, cfg)
        t_prep = time.time() - t0

        v = np.asarray(velocity, dtype=np.float64)
        max_speed = float(np.max(np.abs(v[0]) + np.abs(v[1])))
        dt_max_cfl = (float(dx) / max_speed) if max_speed > 0 else float("inf")

        tracemalloc.start()
        t0 = time.time()
        tr = unit_transport(work, drug, pde, velocity)
        t_solve = time.time() - t0
        _cur, peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        tis = np.asarray(micro.domain.tissue_mask, dtype=bool)
        rec = {
            "dx_um": float(dx), "dt_s": float(dt),
            "grid": [int(micro.domain.ny), int(micro.domain.nx)],
            "n_grid_cells": int(micro.domain.ny * micro.domain.nx),
            "n_tissue_cells": int(tis.sum()),
            "n_steps": int(round(pde.total_time / dt)),
            "max_interstitial_speed_um_s": max_speed,
            "cfl_dt_max_s": dt_max_cfl,
            "cfl_ok_at_this_dt": bool(dt <= dt_max_cfl),
            "cfl_ok_at_dt10": bool(10.0 <= dt_max_cfl),
            "t_reconstruct_s": t_recon,
            "t_prologue_s": t_prep,
            "t_one_solve_s": t_solve,
            "peak_solve_mem_MiB": peak_bytes / 2 ** 20,
            "mean_auc": float(np.asarray(tr.auc)[tis].mean()),
            "penetration_index": float(tr.penetration_index),
        }
        out["points"].append(rec)
        log(f"  cost dx={dx:g} dt={dt:g}: solve {t_solve:.1f}s "
            f"({rec['n_steps']} steps, {rec['n_tissue_cells']} tissue cells), "
            f"peak {rec['peak_solve_mem_MiB']:.0f} MiB, "
            f"CFL dt_max {dt_max_cfl:.1f}s (dt=10 legal: {rec['cfl_ok_at_dt10']}), "
            f"pen {rec['penetration_index']:.4f}")

    ref = [p for p in out["points"] if p["dx_um"] == SHIPPED_DX_UM]
    if ref:
        r = ref[0]
        for p in out["points"]:
            p["solve_cost_vs_dx12"] = p["t_one_solve_s"] / r["t_one_solve_s"]
            p["mem_vs_dx12"] = p["peak_solve_mem_MiB"] / r["peak_solve_mem_MiB"]
    return out


# --------------------------------------------------------------------------- #
# STAGE 4 -- does a length scale PREDICT which drugs move with dx?
# --------------------------------------------------------------------------- #
#: Per-drug dx movement ALREADY MEASURED in this repo. Nothing here is recomputed
#: (two sibling agents own those sweeps); every value is transcribed from the
#: named file, with the delivery law that run used.
#:
#: GBM: ``results/transvascular_exchange_revalidation.md`` section 5, the swept
#: dx = 12/10/8 um table, run with exchange ON (that section's whole subject).
MEASURED_GBM_ON: Dict[str, Dict] = {
    # drug -> ratio at dx = 12 and dx = 8
    "lomustine":    {"r12": 0.4584, "r8": 0.4671},
    "carmustine":   {"r12": 0.6041, "r8": 0.6123},
    "procarbazine": {"r12": 0.4374, "r8": 0.4453},
    "ribociclib":   {"r12": 0.2609, "r8": 0.2677},
    "cisplatin":    {"r12": 0.6184, "r8": 0.6266},
    "doxorubicin":  {"r12": 0.0165, "r8": 0.0165},
}
#: Breast: ``results/dtissue_sensitivity.md`` section 7.2, the dx = 12 (published)
#: vs dx = 24 grid control, run through ``kill_calibration.unit_transport`` and
#: therefore with exchange OFF (stage 0 establishes this).
MEASURED_BREAST_OFF: Dict[str, Dict] = {
    "doxorubicin":     {"r12": 0.2301, "r24": 0.2191},
    "paclitaxel":      {"r12": 0.1503, "r24": 0.1223},
    "docetaxel":       {"r12": 0.1484, "r24": 0.1200},
    "epirubicin":      {"r12": 0.2301, "r24": 0.2191},
    "5-fluorouracil":  {"r12": 0.3379, "r24": 0.3317},
    "cisplatin":       {"r12": 0.3085, "r24": 0.3048},
}


def _spearman(a: Sequence[float], b: Sequence[float]) -> float:
    from scipy.stats import rankdata
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = math.sqrt(float((ra ** 2).sum() * (rb ** 2).sum()))
    return float((ra * rb).sum() / denom) if denom else float("nan")


def stage4_predicts() -> Dict:
    """Correlate the stage-1 length scales against the MEASURED dx movement."""
    with open(os.path.join(OUT_DIR, "blayer_lengthscales.json")) as fh:
        ls = json.load(fh)
    by_tissue = {s["tissue"]: {r["drug"]: r for r in s["drugs"]}
                 for s in ls["sections"]}

    out: Dict = {"panels": []}
    for tissue, measured, law, key_hi in (
            ("gbm", MEASURED_GBM_ON, "exchange ON", "r8"),
            ("breast", MEASURED_BREAST_OFF, "exchange OFF", "r24")):
        rows = []
        for drug, m in measured.items():
            src = by_tissue.get(tissue, {}).get(drug)
            if src is None:
                log(f"  stage4: {tissue}/{drug} absent from stage 1 -- skipped")
                continue
            move = abs(m[key_hi] - m["r12"]) / m["r12"]
            rows.append({
                "drug": drug,
                "delta_median_um": src["delta_um"]["median"],
                "L_median_um": src["L_um"]["median"],
                "bbb_penetrance": src["bbb_penetrance"],
                "ratio_dx12": m["r12"], f"ratio_{key_hi}": m[key_hi],
                "abs_rel_movement": move,
                "signed_rel_movement": (m[key_hi] - m["r12"]) / m["r12"],
            })
        if len(rows) >= 3:
            mv = [r["abs_rel_movement"] for r in rows]
            panel = {
                "tissue": tissue, "law": law, "n": len(rows),
                "comparison": f"dx12 vs {key_hi[1:]} um",
                "spearman_movement_vs_delta": _spearman(
                    [r["delta_median_um"] for r in rows], mv),
                "spearman_movement_vs_L": _spearman(
                    [r["L_median_um"] for r in rows], mv),
                "spearman_movement_vs_bbb": _spearman(
                    [r["bbb_penetrance"] for r in rows], mv),
                "rows": sorted(rows, key=lambda r: -r["abs_rel_movement"]),
            }
            out["panels"].append(panel)
            log(f"  stage4 {tissue} ({law}): rho(move, delta) = "
                f"{panel['spearman_movement_vs_delta']:+.3f}, "
                f"rho(move, L) = {panel['spearman_movement_vs_L']:+.3f}")
    return out


def stage4_fit() -> Dict:
    """Fit ``|rel err of mean| = C * dx^p * delta^q`` on the ALIGNED slab series.

    dx = 30 um is EXCLUDED: a = 48 um is not an integer multiple of it, so those
    points carry interface-quantisation error and belong to the other curve.
    """
    with open(os.path.join(OUT_DIR, "blayer_controlled.json")) as fh:
        ct = json.load(fh)
    xs, ys, zs = [], [], []
    for s in ct["series"]:
        if s["tag"] != "wide_aligned":
            continue
        for p in s["points"]:
            if not p["interface_aligned"] or p["on_rel_err_mean"] == 0:
                continue
            xs.append(math.log(p["dx_um"]))
            ys.append(math.log(s["delta_um"]))
            zs.append(math.log(abs(p["on_rel_err_mean"])))
    a = np.column_stack([np.ones(len(xs)), xs, ys])
    coef, *_ = np.linalg.lstsq(a, np.asarray(zs), rcond=None)
    pred = a @ coef
    ss_res = float(((np.asarray(zs) - pred) ** 2).sum())
    ss_tot = float(((np.asarray(zs) - np.mean(zs)) ** 2).sum())
    return {"n_points": len(xs), "log_C": float(coef[0]),
            "C": float(math.exp(coef[0])), "p_dx": float(coef[1]),
            "q_delta": float(coef[2]), "r2": 1.0 - ss_res / ss_tot,
            "form": "|rel_err_mean| = C * dx_um^p * delta_um^q  (exchange ON, "
                    "interface aligned, L_out = 64.92 um fixed)"}


# --------------------------------------------------------------------------- #
def _json_safe(obj):
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items() if k != "x" and k != "c"}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "0", "1", "1b", "2", "3", "4"])
    ap.add_argument("--dx", type=float, default=SHIPPED_DX_UM)
    ap.add_argument("--sections", default="gbm,breast")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    if args.stage in ("all", "0"):
        r0 = stage0_trap_check()
        with open(os.path.join(OUT_DIR, "blayer_trap_check.json"), "w") as fh:
            json.dump(_json_safe(r0), fh, indent=2)
        log("stage 0 written")
        print(json.dumps(_json_safe(r0), indent=2))

    if args.stage in ("all", "1"):
        r1 = stage1(args.dx, tuple(s.strip() for s in args.sections.split(",")))
        with open(os.path.join(OUT_DIR, "blayer_lengthscales.json"), "w") as fh:
            json.dump(_json_safe(r1), fh, indent=2)
        log("stage 1 written -> results/sweep/blayer_lengthscales.json")

    if args.stage in ("all", "1b"):
        r1b = stage1b_source_vs_dx(
            DX_GRID_UM, tuple(s.strip() for s in args.sections.split(",")))
        with open(os.path.join(OUT_DIR, "blayer_source_vs_dx.json"), "w") as fh:
            json.dump(_json_safe(r1b), fh, indent=2)
        log("stage 1b written -> results/sweep/blayer_source_vs_dx.json")

    if args.stage in ("all", "2"):
        trap = stage2_trap_assert()
        log(f"stage 2 trap check: on/off mean ratio = "
            f"{trap['ratio_on_over_off']:.6f}, bit_identical="
            f"{trap['bit_identical']}")
        r2 = stage2()
        r2["trap_check"] = trap
        with open(os.path.join(OUT_DIR, "blayer_controlled.json"), "w") as fh:
            json.dump(_json_safe(r2), fh, indent=2)
        log("stage 2 written -> results/sweep/blayer_controlled.json")

    if args.stage in ("all", "3"):
        r3 = stage3_cost()
        with open(os.path.join(OUT_DIR, "blayer_cost.json"), "w") as fh:
            json.dump(_json_safe(r3), fh, indent=2)
        log("stage 3 written -> results/sweep/blayer_cost.json")

    if args.stage in ("all", "4"):
        r4 = stage4_predicts()
        r4["error_law_fit"] = stage4_fit()
        with open(os.path.join(OUT_DIR, "blayer_analysis.json"), "w") as fh:
            json.dump(_json_safe(r4), fh, indent=2)
        log("stage 4 written -> results/sweep/blayer_analysis.json")


if __name__ == "__main__":
    main()
