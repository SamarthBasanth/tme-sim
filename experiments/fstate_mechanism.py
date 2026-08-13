"""What actually makes ``F_state`` move with ``dx``?

THE QUESTION, AND WHY THE OBVIOUS ANSWER IS ALREADY DEAD
-------------------------------------------------------
``results/grid_convergence_gbm.md`` sec. 5 measures the oxygen-banded cell-state
multiplier ``F_state`` swinging **0.5041-0.5502** (0.0461) across the GBM dx
ladder, and ``results/tail90_falsification.md`` sec. 6.1 measures a **3.72 %**
fine-rung span that carries 25-51 % of the whole ladder's variance. Both
documents attribute that movement to perivascular under-resolution -- the mesh
being unable to resolve the oxygen boundary layer around vessels.

``results/subgrid_closure.md`` sec. 2 killed that attribution: under the SHIPPED
oxygen constants the worst discretisation error of the oxygen solve on a
controlled problem with an exact solution is **0.170 % at dx = 30 um**. A field
that accurate cannot, by itself, move a banded fraction by 4 points.

So the field is right and the label still moves. This module finds out why, by
decomposition rather than by argument.

THE DECOMPOSITION
-----------------
``F_state`` is a functional of one thing only: the per-tumour-cell oxygen value
``o_i`` that ``cells.agents.update_cell_states`` bands. Four separable things
stand between the continuum problem and ``o_i``, and each is switched on one at
a time:

===========  =====================================  ===========================
arm          field solved on                        cells read it by
===========  =====================================  ===========================
``D``        reference grid, reference source       bilinear at (x_i, y_i)
``E(dx)``    reference grid, **dx-rasterised src**  bilinear at (x_i, y_i)
``B(dx)``    **dx grid**, dx-rasterised source      bilinear at (x_i, y_i)
``A(dx)``    dx grid, dx-rasterised source          **production floor index**
===========  =====================================  ===========================

``A`` is the shipped number and is asserted equal to the production pipeline.
Reading down the table, each row switches exactly one mechanism on:

* ``E - D``  = **source rasterisation** (the vessel/sink painting changes with dx)
* ``B - E``  = **discretisation error** of the oxygen PDE itself
* ``A - B``  = **cell-to-grid sampling** (``Domain.to_index`` is a floor, so a
  cell is given the O2 of the grid cell CONTAINING it, i.e. of a point up to
  ``dx/sqrt(2)`` away)
* ``A - D``  = the total, which is what the two documents measured

THE DECISIVE CONTROL
--------------------
If sampling is the cause then it must reproduce itself with the oxygen field
**held completely fixed**. Arms ``C_avg`` and ``C_inj`` do that: one single
reference field, never re-solved, re-read at every dx through the production
floor rule (``C_avg`` block-averages it onto the dx grid first, as a finite
volume scheme would; ``C_inj`` point-injects). Nothing about the PDE, the
source or the mesh convergence can enter those two numbers -- only the sampling.

Everything here is READ-ONLY with respect to ``tmesim/``. No default, band edge
or constant is changed; ``solve_oxygen``'s linear branch is re-implemented
locally ONLY so that an arbitrary (v, rho) pair can be fed to it, and that
re-implementation is asserted bit-identical to the shipped function (sec.
``--stage identity``).

Usage
-----
    ./.venv/bin/python experiments/fstate_mechanism.py --stage ladder
    ./.venv/bin/python experiments/fstate_mechanism.py --stage fixed
    ./.venv/bin/python experiments/fstate_mechanism.py --stage raster
    ./.venv/bin/python experiments/fstate_mechanism.py --stage edges
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Sequence

import numpy as np
import scipy.sparse as sp
from scipy.ndimage import map_coordinates
from scipy.sparse.linalg import spsolve

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.cells.agents import update_cell_states                  # noqa: E402
from tmesim.core import CellState, PDConfig                          # noqa: E402
from tmesim.core.config import OxygenConfig                          # noqa: E402
from tmesim.physics.oxygen import (                                  # noqa: E402
    TUMOR_REFERENCE_AREA_UM2,
    _build_laplacian,
    _tumor_density,
    solve_oxygen,
)
from tmesim.reconstruct.real_pipeline import reconstruct_real        # noqa: E402
from tmesim.validation.kill_calibration import resistance_vector     # noqa: E402

OUT_DIR = os.path.join(_ROOT, "results", "sweep")

#: The published GBM ladder (``results/grid_convergence_gbm.md`` sec. 2), plus
#: three finer rungs that the oxygen solve (unconditionally stable, no dt) can
#: reach even though the ADVECTION CFL stops the drug ladder at 7 um. Nothing
#: below depends on a drug transport solve, so the CFL floor does not apply.
DX_LADDER = (30.0, 24.0, 18.0, 16.0, 12.0, 10.0, 8.0, 7.0, 6.0, 5.0)

#: ``results/grid_convergence_gbm.md`` sec. 0 splits coarse from fine here.
FINE_MAX_DX = 18.0

#: The reference grid for the fixed-field control and for arms D / E. Chosen as
#: the finest rung the direct sparse solve completes (see --stage identity).
REF_DX = 4.0


def log(msg: str) -> None:
    print(f"[fstate {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# F_state, exactly as the production pipeline computes it
# --------------------------------------------------------------------------- #
def band_states(o2: np.ndarray, cfg: OxygenConfig) -> np.ndarray:
    """``cells.agents.update_cell_states``'s banding, on a per-cell O2 vector.

    Duplicated in closed form (not called) because the shipped function insists
    on sampling the field itself through ``Domain.to_index`` -- and the whole
    point here is to vary how that sample is taken. Asserted equal to the
    shipped function on the production arm (:func:`assert_matches_production`).
    """
    t = float(cfg.hypoxia_threshold)
    o2 = np.asarray(o2, dtype=np.float64)
    st = np.full(o2.shape, int(CellState.VIABLE), dtype=np.int8)
    st[o2 < t] = int(CellState.HYPOXIC)
    st[(o2 >= t) & (o2 < 2.0 * t)] = int(CellState.QUIESCENT)
    return st


def f_state(o2_tumour: np.ndarray, cfg: OxygenConfig,
            pd_cfg: PDConfig) -> Dict:
    """``F_state = mean_i[resistance_i]`` over TUMOUR cells, plus the bands."""
    st = band_states(o2_tumour, cfg)
    r = resistance_vector(st, pd_cfg)
    return {
        "f_state": float(np.mean(r)),
        "quiescent_frac": float(np.mean(st == int(CellState.QUIESCENT))),
        "hypoxic_frac": float(np.mean(st == int(CellState.HYPOXIC))),
        "viable_frac": float(np.mean(st == int(CellState.VIABLE))),
        "mean_o2": float(np.mean(o2_tumour)),
    }


# --------------------------------------------------------------------------- #
# The two samplers
# --------------------------------------------------------------------------- #
def sample_floor(field: np.ndarray, dx: float, x: np.ndarray,
                 y: np.ndarray) -> np.ndarray:
    """PRODUCTION rule: ``Domain.to_index`` -- floor to the containing cell.

    Reproduces ``core.domain.Domain.to_index`` exactly (clip included). This is
    a piecewise-constant read: cell ``i`` is handed the value stored at the
    CENTRE of whichever grid cell contains it, a point up to ``dx/sqrt(2)``
    away from where the cell actually is.
    """
    ny, nx = field.shape
    col = np.clip((np.asarray(x) / dx).astype(int), 0, nx - 1)
    row = np.clip((np.asarray(y) / dx).astype(int), 0, ny - 1)
    return field[row, col]


def sample_bilinear(field: np.ndarray, dx: float, x: np.ndarray,
                    y: np.ndarray) -> np.ndarray:
    """dx-INDEPENDENT rule: bilinear interpolation at the cell's true position.

    Grid values live at cell centres ``((j+0.5)*dx, (i+0.5)*dx)``
    (``Domain.coords_um``), so the fractional index of a point is
    ``x/dx - 0.5``. ``mode='nearest'`` extends the edge value, matching the
    no-flux boundary the solve was run with.
    """
    coords = np.vstack([np.asarray(y) / dx - 0.5, np.asarray(x) / dx - 0.5])
    return map_coordinates(field, coords, order=1, mode="nearest")


# --------------------------------------------------------------------------- #
# A local copy of solve_oxygen's linear branch, so (v, rho) can be swapped
# --------------------------------------------------------------------------- #
def solve_from_fields(vasc: np.ndarray, rho_tum: np.ndarray, dx: float,
                      cfg: OxygenConfig) -> np.ndarray:
    """``(-D*L + diag(supply*v + consumption*rho)) @ O = supply*v``.

    Line-for-line the ``sink_form == 'linear'`` branch of
    ``physics.oxygen.solve_oxygen``, differing ONLY in taking ``vasc`` and
    ``rho_tum`` as arguments instead of deriving them from a Microenvironment.
    Asserted bit-identical on the real section by :func:`stage_identity`.
    """
    ny, nx = vasc.shape
    vasc = np.clip(np.asarray(vasc, dtype=np.float64), 0.0, 1.0)
    D = float(cfg.diffusivity)
    L = _build_laplacian(ny, nx, float(dx))
    source = (float(cfg.supply) * vasc).ravel()
    reaction_bulk = (float(cfg.consumption) * np.asarray(rho_tum, float)).ravel()
    if not np.any(source) and not np.any(reaction_bulk):
        return np.zeros((ny, nx), dtype=np.float64)
    A = (-D * L) + sp.diags(source + reaction_bulk, format="csr")
    o = spsolve(A.tocsc(), source)
    o = np.asarray(o, dtype=np.float64).reshape(ny, nx)
    np.clip(o, 0.0, 1.0, out=o)
    return o


# --------------------------------------------------------------------------- #
# Grid-to-grid transfer, in PHYSICAL coordinates (the grids do not nest)
# --------------------------------------------------------------------------- #
def inject_to(field: np.ndarray, dx_src: float, shape_dst, dx_dst: float
              ) -> np.ndarray:
    """Piecewise-constant reconstruction of a coarse field on a finer grid.

    For every destination cell centre, take the value of whichever SOURCE cell
    contains it -- i.e. the exact piecewise-constant function the coarse field
    represents, re-tabulated. Conservative for an intensive quantity (a density
    or an occupancy), which is what both ``vasc`` and ``rho_tum`` are.
    """
    ny_d, nx_d = shape_dst
    ny_s, nx_s = field.shape
    xs = (np.arange(nx_d) + 0.5) * dx_dst
    ys = (np.arange(ny_d) + 0.5) * dx_dst
    col = np.clip((xs / dx_src).astype(int), 0, nx_s - 1)
    row = np.clip((ys / dx_src).astype(int), 0, ny_s - 1)
    return field[np.ix_(row, col)]


def block_average(field: np.ndarray, dx_src: float, shape_dst, dx_dst: float
                  ) -> np.ndarray:
    """Average a fine field over each coarse cell -- the finite-volume restriction.

    This is what a coarse-grid cell-centred value MEANS: the mean of the
    continuum over that cell. Destination cells with no source cell inside them
    (impossible here, ``dx_dst >> dx_src``) would fall back to injection.
    """
    ny_d, nx_d = shape_dst
    ny_s, nx_s = field.shape
    xs = (np.arange(nx_s) + 0.5) * dx_src
    ys = (np.arange(ny_s) + 0.5) * dx_src
    col = np.clip((xs / dx_dst).astype(int), 0, nx_d - 1)
    row = np.clip((ys / dx_dst).astype(int), 0, ny_d - 1)
    flat = (row[:, None] * nx_d + col[None, :]).ravel()
    tot = np.bincount(flat, weights=field.ravel(), minlength=ny_d * nx_d)
    cnt = np.bincount(flat, minlength=ny_d * nx_d)
    out = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    return out.reshape(ny_d, nx_d)


def block_centre_of(x: np.ndarray, y: np.ndarray, dx: float, shape):
    """Where the production sampler EFFECTIVELY reads each cell's oxygen.

    The centre of the ``dx`` grid cell containing ``(x_i, y_i)``. Feeding these
    displaced coordinates to a fixed fine field isolates the relocation alone.
    """
    ny, nx = shape
    col = np.clip((np.asarray(x) / dx).astype(int), 0, nx - 1)
    row = np.clip((np.asarray(y) / dx).astype(int), 0, ny - 1)
    return (col + 0.5) * dx, (row + 0.5) * dx


# --------------------------------------------------------------------------- #
# Stages
# --------------------------------------------------------------------------- #
def _dump(obj, name: str) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=float)
    log(f"written -> results/sweep/{name}")


def stage_identity() -> Dict:
    """Preconditions: the local solver is the shipped one, the local banding is
    the shipped one, and the reference grid actually solves."""
    cfg, pd_cfg = OxygenConfig(), PDConfig()
    out: Dict = {}
    dx = 12.0
    micro = reconstruct_real(source="xenium", dx_um=dx)
    o_ship = solve_oxygen(micro, cfg)
    rho = _tumor_density(micro)
    o_loc = solve_from_fields(micro.fields.vascular_density, rho, dx, cfg)
    out["solver_max_abs_diff"] = float(np.max(np.abs(o_ship - o_loc)))
    out["solver_bit_identical"] = bool(np.array_equal(o_ship, o_loc))

    # banding: shipped update_cell_states vs the closed form + floor sampler
    micro.fields.oxygen = o_ship
    update_cell_states(micro, o_ship, cfg)
    tum = micro.cells.tumor_mask()
    ship_states = np.asarray(micro.cells.state)[tum]
    mine = band_states(
        sample_floor(o_ship, dx, micro.cells.x[tum], micro.cells.y[tum]), cfg)
    out["band_states_identical"] = bool(np.array_equal(ship_states, mine))
    out["f_state_shipped"] = float(np.mean(resistance_vector(ship_states, pd_cfg)))
    out["f_state_local"] = float(np.mean(resistance_vector(mine, pd_cfg)))
    out["dx_um"] = dx
    log(f"identity: solver bit-identical={out['solver_bit_identical']} "
        f"(max|diff|={out['solver_max_abs_diff']:.3e}), bands identical="
        f"{out['band_states_identical']}, F_state={out['f_state_local']:.4f}")
    return out


def _rung(dx: float, cfg: OxygenConfig) -> Dict:
    """Reconstruct + solve one rung, returning everything later stages need."""
    t0 = time.time()
    micro = reconstruct_real(source="xenium", dx_um=dx)
    vasc = np.clip(np.asarray(micro.fields.vascular_density, float), 0.0, 1.0)
    rho = _tumor_density(micro)
    o = solve_from_fields(vasc, rho, dx, cfg)
    tum = micro.cells.tumor_mask()
    return {
        "dx": dx, "micro": micro, "vasc": vasc, "rho": rho, "oxygen": o,
        "tum": tum, "x": micro.cells.x, "y": micro.cells.y,
        "shape": tuple(micro.domain.shape), "secs": time.time() - t0,
    }


def stage_ladder(dx_ladder: Sequence[float] = DX_LADDER) -> Dict:
    """Arms A and B at every rung, plus the raw per-cell O2 vectors."""
    cfg, pd_cfg = OxygenConfig(), PDConfig()
    rows: List[Dict] = []
    os.makedirs(OUT_DIR, exist_ok=True)
    npz: Dict[str, np.ndarray] = {}
    for dx in dx_ladder:
        r = _rung(float(dx), cfg)
        xt, yt = r["x"][r["tum"]], r["y"][r["tum"]]
        o_floor = sample_floor(r["oxygen"], dx, xt, yt)
        o_bilin = sample_bilinear(r["oxygen"], dx, xt, yt)
        A = f_state(o_floor, cfg, pd_cfg)
        B = f_state(o_bilin, cfg, pd_cfg)
        row = {
            "dx_um": float(dx), "grid": list(r["shape"]),
            "n_tumour": int(r["tum"].sum()),
            "A_production": A, "B_bilinear": B,
            # how much source the dx rasterisation actually paints
            "vasc_integral_um2": float(r["vasc"].sum() * dx * dx),
            "vasc_frac_nonzero": float(np.mean(r["vasc"] > 0)),
            "rho_integral": float(r["rho"].sum() * dx * dx
                                  / TUMOR_REFERENCE_AREA_UM2),
            "field_mean": float(r["oxygen"].mean()),
            "secs": r["secs"],
        }
        rows.append(row)
        npz[f"floor_{dx:g}"] = o_floor.astype(np.float32)
        npz[f"bilin_{dx:g}"] = o_bilin.astype(np.float32)
        log(f"dx={dx:g}: grid={r['shape']} A(prod) F={A['f_state']:.4f} "
            f"q={A['quiescent_frac']:.4f} | B(bilin) F={B['f_state']:.4f} "
            f"q={B['quiescent_frac']:.4f} | vasc_int="
            f"{row['vasc_integral_um2']:.0f} ({r['secs']:.0f}s)")
        _dump({"rows": rows}, "fstate_ladder.json")
        np.savez_compressed(os.path.join(OUT_DIR, "fstate_o2_cells.npz"), **npz)
    return {"rows": rows}


def stage_fixed(ref_dx: float = REF_DX,
                dx_ladder: Sequence[float] = DX_LADDER) -> Dict:
    """THE DECISIVE CONTROL: one field, never re-solved, read at every dx."""
    cfg, pd_cfg = OxygenConfig(), PDConfig()
    log(f"solving the reference field once at dx = {ref_dx:g} um")
    ref = _rung(float(ref_dx), cfg)
    o_ref = ref["oxygen"]
    xt, yt = ref["x"][ref["tum"]], ref["y"][ref["tum"]]
    log(f"reference grid {ref['shape']} solved in {ref['secs']:.0f}s")

    # D: the dx-independent baseline -- the reference field read at the cells'
    # true positions. This is the number every arm would return if nothing
    # depended on the mesh.
    D = f_state(sample_bilinear(o_ref, ref_dx, xt, yt), cfg, pd_cfg)
    log(f"D (baseline, dx-independent): F_state={D['f_state']:.4f} "
        f"q={D['quiescent_frac']:.4f} mean_o2={D['mean_o2']:.4f}")

    rows: List[Dict] = []
    for dx in dx_ladder:
        if dx <= ref_dx:
            continue
        ny = int(np.ceil(ref["shape"][0] * ref_dx / dx))
        nx = int(np.ceil(ref["shape"][1] * ref_dx / dx))
        # C_avg: finite-volume restriction of the FIXED field, then floor-read
        avg = block_average(o_ref, ref_dx, (ny, nx), dx)
        avg = np.where(np.isnan(avg), 0.0, avg)
        C_avg = f_state(sample_floor(avg, dx, xt, yt), cfg, pd_cfg)
        # C_inj: the FIXED field read at the centre of the containing dx block
        cx, cy = block_centre_of(xt, yt, dx, (ny, nx))
        C_inj = f_state(sample_bilinear(o_ref, ref_dx, cx, cy), cfg, pd_cfg)
        disp = np.hypot(cx - xt, cy - yt)
        rows.append({
            "dx_um": float(dx), "grid": [ny, nx],
            "C_avg": C_avg, "C_inj": C_inj,
            "mean_displacement_um": float(disp.mean()),
            "max_displacement_um": float(disp.max()),
            "mean_abs_o2_shift_avg": float(np.mean(np.abs(
                sample_floor(avg, dx, xt, yt)
                - sample_bilinear(o_ref, ref_dx, xt, yt)))),
            "mean_abs_o2_shift_inj": float(np.mean(np.abs(
                sample_bilinear(o_ref, ref_dx, cx, cy)
                - sample_bilinear(o_ref, ref_dx, xt, yt)))),
        })
        log(f"dx={dx:g}: C_avg F={C_avg['f_state']:.4f} "
            f"q={C_avg['quiescent_frac']:.4f} | C_inj F={C_inj['f_state']:.4f} "
            f"q={C_inj['quiescent_frac']:.4f} | mean|dO2|="
            f"{rows[-1]['mean_abs_o2_shift_avg']:.4f} "
            f"mean disp={disp.mean():.1f} um")
        _dump({"ref_dx_um": float(ref_dx), "ref_grid": list(ref["shape"]),
               "D_baseline": D, "rows": rows}, "fstate_fixed_field.json")
    return {"ref_dx_um": float(ref_dx), "ref_grid": list(ref["shape"]),
            "D_baseline": D, "rows": rows}


def stage_raster(ref_dx: float = REF_DX,
                 dx_ladder: Sequence[float] = DX_LADDER) -> Dict:
    """Arm E: the dx-rasterised SOURCE, solved on the reference grid.

    Isolates 'how much source landed where' from 'how well the PDE was solved'.
    """
    cfg, pd_cfg = OxygenConfig(), PDConfig()
    ref = _rung(float(ref_dx), cfg)
    xt, yt = ref["x"][ref["tum"]], ref["y"][ref["tum"]]
    D = f_state(sample_bilinear(ref["oxygen"], ref_dx, xt, yt), cfg, pd_cfg)
    rows: List[Dict] = []
    for dx in dx_ladder:
        if dx <= ref_dx:
            continue
        c = _rung(float(dx), cfg)
        v_f = inject_to(c["vasc"], dx, ref["shape"], ref_dx)
        r_f = inject_to(c["rho"], dx, ref["shape"], ref_dx)
        o_e = solve_from_fields(v_f, r_f, ref_dx, cfg)
        E = f_state(sample_bilinear(o_e, ref_dx, xt, yt), cfg, pd_cfg)
        rows.append({
            "dx_um": float(dx), "E": E,
            "vasc_integral_coarse_um2": float(c["vasc"].sum() * dx * dx),
            "vasc_integral_injected_um2": float(v_f.sum() * ref_dx * ref_dx),
            "rho_integral_injected": float(r_f.sum() * ref_dx * ref_dx
                                           / TUMOR_REFERENCE_AREA_UM2),
        })
        log(f"dx={dx:g}: E(coarse source, fine solve) F={E['f_state']:.4f} "
            f"q={E['quiescent_frac']:.4f} mean_o2={E['mean_o2']:.4f}")
        _dump({"ref_dx_um": float(ref_dx), "D_baseline": D, "rows": rows},
              "fstate_raster.json")
    return {"ref_dx_um": float(ref_dx), "D_baseline": D, "rows": rows}


def stage_edges(ref_dx: float = REF_DX) -> Dict:
    """Is the population piled against the band edge, and how hard?

    Measures the O2 distribution near ``t`` and ``2t`` on the reference field,
    the local gradient the sampler has to fight, and -- the load-bearing number
    -- what fraction of cells sit within the per-cell O2 shift a dx step
    actually produces.
    """
    cfg, pd_cfg = OxygenConfig(), PDConfig()
    t = float(cfg.hypoxia_threshold)
    ref = _rung(float(ref_dx), cfg)
    o_ref = ref["oxygen"]
    xt, yt = ref["x"][ref["tum"]], ref["y"][ref["tum"]]
    o = sample_bilinear(o_ref, ref_dx, xt, yt)
    n = o.size

    gy, gx = np.gradient(o_ref, ref_dx)
    gmag = np.hypot(gx, gy)
    g_cells = sample_bilinear(gmag, ref_dx, xt, yt)

    out: Dict = {
        "ref_dx_um": float(ref_dx), "n_tumour": int(n),
        "band_edges": {"hypoxia_t": t, "quiescent_ceiling": 2.0 * t},
        "mean_o2": float(o.mean()), "median_o2": float(np.median(o)),
        "percentiles": {str(p): float(np.percentile(o, p))
                        for p in (1, 5, 10, 25, 50, 75, 90, 95, 99)},
        "frac_below_ceiling": float(np.mean(o < 2.0 * t)),
        "frac_below_t": float(np.mean(o < t)),
        "mean_grad_um": float(g_cells.mean()),
        "median_grad_um": float(np.median(g_cells)),
        # how many cells sit within eps of the 0.4 ceiling
        "within_of_ceiling": {
            f"{eps:g}": float(np.mean(np.abs(o - 2.0 * t) < eps))
            for eps in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)},
        "within_of_t": {
            f"{eps:g}": float(np.mean(np.abs(o - t) < eps))
            for eps in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)},
        # local density of the O2 distribution AT the ceiling (per unit O2),
        # which is the sensitivity dq/d(shift) in the small-shift limit
        "density_at_ceiling_per_unit_o2": float(
            np.mean(np.abs(o - 2.0 * t) < 0.01) / 0.02),
        "density_at_t_per_unit_o2": float(
            np.mean(np.abs(o - t) < 0.01) / 0.02),
    }
    # A histogram tight around the ceiling, so the shape is on the record.
    edges = np.linspace(2.0 * t - 0.05, 2.0 * t + 0.05, 21)
    out["ceiling_histogram"] = {
        "edges": [float(e) for e in edges],
        "counts": [int(c) for c in np.histogram(o, bins=edges)[0]],
    }
    log(f"edges: mean_o2={out['mean_o2']:.4f}, {out['frac_below_ceiling']:.4f} "
        f"below the 0.4 ceiling, {out['within_of_ceiling']['0.01']*100:.2f} % "
        f"within 0.01 of it, mean|grad O2| = {out['mean_grad_um']:.5f} /um")
    _dump(out, "fstate_band_edges.json")
    return out


def stage_fix(dx_ladder: Sequence[float] = DX_LADDER) -> Dict:
    """Is the dominant term fixable? Measure, do not assert.

    The dominant term is the BINARY centre-in-disc vessel painting, whose
    integral is not conserved under dx. ``vessel_network.discrete_vascular_field``
    already carries an opt-in ``geometry='subgrid'`` partial-area painting whose
    docstring claims the painted total is 'computed without reference to
    domain.dx, so it is grid-independent by construction'. That is a claim, so
    it is measured here rather than believed.

    **NOTHING SHIPPED IS CHANGED.** ``geometry`` is passed as a keyword to an
    existing opt-in branch; ``_DEFAULT_VASCULAR_MODE`` and the ``'binary'``
    default are untouched, and the shipped arm is re-measured alongside so the
    two are compared on identical reconstructions.
    """
    from tmesim.reconstruct.real_fields import vascular_density_from_cells
    cfg = OxygenConfig()
    rows: List[Dict] = []
    for dx in dx_ladder:
        dx = float(dx)
        micro = reconstruct_real(source="xenium", dx_um=dx)
        dom, cells = micro.domain, micro.cells
        v_bin = vascular_density_from_cells(dom, cells, geometry="binary")
        v_sub = vascular_density_from_cells(dom, cells, geometry="subgrid")
        rows.append({
            "dx_um": dx, "grid": list(dom.shape),
            "binary_integral_um2": float(v_bin.sum() * dx * dx),
            "subgrid_integral_um2": float(v_sub.sum() * dx * dx),
            "shipped_matches_binary": bool(np.array_equal(
                np.clip(np.asarray(micro.fields.vascular_density, float), 0, 1),
                np.clip(v_bin, 0, 1))),
        })
        log(f"dx={dx:g}: binary integral={rows[-1]['binary_integral_um2']:.0f} "
            f"um2 | subgrid integral={rows[-1]['subgrid_integral_um2']:.0f} um2")
        _dump({"rows": rows}, "fstate_source_fix.json")
    b = np.array([r["binary_integral_um2"] for r in rows])
    s = np.array([r["subgrid_integral_um2"] for r in rows])
    out = {"rows": rows,
           "binary_span_pct": float(100 * (b.max() - b.min()) / b.mean()),
           "subgrid_span_pct": float(100 * (s.max() - s.min()) / s.mean())}
    log(f"painted-area span across the ladder: binary "
        f"{out['binary_span_pct']:.2f} % vs subgrid {out['subgrid_span_pct']:.2f} %")
    _dump(out, "fstate_source_fix.json")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="ladder",
                    choices=["identity", "ladder", "fixed", "raster", "edges",
                             "fix"])
    ap.add_argument("--ref-dx", type=float, default=REF_DX)
    args = ap.parse_args()
    if args.stage == "identity":
        _dump(stage_identity(), "fstate_identity.json")
    elif args.stage == "ladder":
        stage_ladder()
    elif args.stage == "fixed":
        stage_fixed(ref_dx=args.ref_dx)
    elif args.stage == "raster":
        stage_raster(ref_dx=args.ref_dx)
    elif args.stage == "edges":
        stage_edges(ref_dx=args.ref_dx)
    elif args.stage == "fix":
        stage_fix()


if __name__ == "__main__":
    main()
