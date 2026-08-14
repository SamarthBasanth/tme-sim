"""Can a sub-grid closure represent the perivascular boundary layer the mesh cannot?

THE QUESTION
------------
``results/oxygen_joint_grid_test.md`` §4 measured, from the constants alone, that
with the jointly calibrated oxygen configuration the source equilibration length
``sqrt(D/supply)`` is **3.2 um** and the sink diffusion length
``sqrt(D/(consumption*rho))`` is **10-14 um**, against a production mesh of
``dx = 12 um`` -- **0.47 to 1.77 mesh cells per length scale** -- and concluded
that *"a physically correct vessel is a boundary-layer object this project's
meshes cannot represent"*. This file asks whether a sub-grid closure can
represent it instead, and holds itself to a falsification criterion.

WHY THIS IS TRACTABLE
---------------------
The training data can be GENERATED, not found. The controlled problem already
exists: ``experiments/boundary_layer_resolution.py`` builds a 1-D two-region slab
with a verified closed-form steady state (``slab_exact``) and exact cell averages
(``slab_exact_reference``), and this file imports both UNMODIFIED. The oxygen
operator is the same two-region reaction-diffusion problem with
``k_in = supply + q``, ``s_in = supply``, ``k_out = q``, so one closed form scores
every point below.

WHAT IS BEING CLOSED, AND WHY THAT AND NOT SOMETHING ELSE (stage `identity`)
----------------------------------------------------------------------------
The shipped steady solve is ``(-D L + diag(supply*v + q)) O = supply*v``. ``L``
has zero row sums, so summing over cells gives, exactly and at every ``dx``,
``sum_j (supply*v_j + q_j) O_j = sum_j supply*v_j``; the continuum obeys the same
identity with integrals. Therefore the ENTIRE domain-mean error is the mismatch
between ``supply*v_j*O_j`` and ``integral(supply*v*O)`` over the cell -- a
SOURCE-WEIGHTING error. Stage `identity` verifies that claim by construction: it
computes the exact per-cell factor
``beta = mean_vessel(O_exact) / mean_cell(O_exact)`` and shows the coarse solve
with that factor reproduces the exact cell averages to machine precision, while
an exact face-flux ("transmissibility") oracle alone does not. That is what
picks the closure's OUTPUT: one scalar per cell multiplying the source's
contribution to the reaction diagonal.

THE CLOSURE ITSELF (``tmesim/physics/subgrid.py``, opt-in, default OFF)
----------------------------------------------------------------------
``beta`` is then evaluated in CLOSED FORM from a local cell problem with ZERO
fitted parameters -- a function of ``(v, dx, supply, q, D)``, all of which the
engine already holds per cell. Nothing is trained. Per the brief's ordering, the
analytic baseline is fitted and reported FIRST; a learned closure is only
justified if it beats a parameter-free formula, and stage `analytic` reports
whether it needs to exist.

THE STRUCTURAL PRECONDITION
---------------------------
``beta == 1`` identically where ``v == 0`` or ``v == 1``. The SHIPPED vascular
field is binary (``geometry="binary"``), so on it the closure is an exact no-op.
Stage `transfer` therefore runs the tissue ladder on the opt-in fractional-area
field (``geometry="subgrid"``, ``results/vessel_geometry_fix.md``) and reports
the binary arm alongside it -- the closure is measured against the SAME source
field with and without it, so the comparison isolates the closure.

STAGES (each checkpoints to results/sweep/subgrid_<stage>.json)
---------------------------------------------------------------
  identity  -- the conservation identity, the oracle decomposition, and the
               machine-precision proof that the source reweighting is a COMPLETE
               parameterisation of the error. No fitting.
  sweep     -- generate the controlled data: lambda_sink 3-320 um, lambda_source
               3-200 um, dx 2-30 um, source half-width swept. Split into TRAIN
               and HELD-OUT dx.
  analytic  -- score the analytic closure out of sample, against the shipped
               solve and against the power-law error-law baseline the brief
               names. Decides whether anything needs training.
  transfer  -- the honest question: the REAL GBM section, dx = 30/18/12/8 (+ a
               fine probe), F_state span with and without the closure.

Run:
    ./.venv/bin/python experiments/subgrid_closure.py identity sweep analytic transfer

NOTHING SHIPPED IS CHANGED. ``OxygenConfig`` defaults, the vascular painting
default and ``physics/oxygen.py`` are untouched; the closure lives in a new
default-off module that nothing else imports.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# The controlled problem and its verified closed form, imported UNMODIFIED.
from boundary_layer_resolution import slab_exact, slab_exact_reference  # noqa: E402

from tmesim.physics.oxygen import _build_laplacian  # noqa: E402
from tmesim.physics.subgrid import source_reweight  # noqa: E402

SWEEP = os.path.join(_ROOT, "results", "sweep")
_T0 = time.time()

#: Shipped OxygenConfig.diffusivity.
D_O2 = 2000.0
#: Shipped VESSEL_RADIUS_UM.
VESSEL_R = 7.5
#: Controlled-problem domain half-length (um) -- the Krogh HALF-SPACING between
#: vessels, with no-flux at both ends (symmetry at the vessel centre and at the
#: midpoint between two vessels). It is SWEPT, not fixed: with a single vessel in
#: a 480 um domain the whole field sits below the 0.2 band edge and the banded
#: fraction -- the quantity F_state is built from, and the harshest test of a
#: closure because it is a step function of the field -- is degenerate. The real
#: section's perfused area fraction of 2.68 % at a 7.5 um radius implies a mean
#: intervessel spacing near 80 um, i.e. a half-spacing near 40 um, which is
#: inside the swept range.
HALF_SPACING_UM: Tuple[float, ...] = (24.0, 48.0, 96.0, 240.0, 480.0)
#: Default used by the stage-`identity` oracle, which is about the OPERATOR and
#: not about a banded fraction.
LX = 480.0

#: Shipped and jointly calibrated oxygen constants (nothing re-derived here;
#: taken from results/oxygen_supply_calibration.md and oxygen_joint_grid_test.md).
SHIPPED_SUPPLY, SHIPPED_CONSUMPTION = 0.05, 0.02
JOINT_SUPPLY = 194.0
JOINT_CONSUMPTION_LO, JOINT_CONSUMPTION_HI = 10.0, 20.0

#: dx the closure is DERIVED for -- none, it has no fitted parameter. These are
#: the grids the fit of the power-law BASELINE is allowed to see.
DX_TRAIN: Tuple[float, ...] = (2.0, 4.0, 6.0, 12.0, 24.0)
#: Grids no fit of any kind may see. The falsification is scored here.
DX_HELDOUT: Tuple[float, ...] = (3.0, 8.0, 10.0, 16.0, 20.0, 30.0)
DX_ALL: Tuple[float, ...] = tuple(sorted(DX_TRAIN + DX_HELDOUT))

#: Band edges of the shipped cell-state labelling (OxygenConfig.hypoxia_threshold
#: = 0.2 and the QUIESCENT/VIABLE edge at 0.4, cells/agents.py) and the
#: PDConfig resistance multipliers they map to. The banded fraction is the
#: SCORED quantity: it is what F_state is built from and it is a step function
#: of the field, so it is the harshest test of a closure.
BAND_EDGES = (0.2, 0.4)
BAND_MULT = (0.35, 0.5, 1.0)    # HYPOXIC, QUIESCENT, VIABLE (PDConfig defaults)


def log(msg: str) -> None:
    print(f"[subgrid +{time.time() - _T0:7.1f}s] {msg}", flush=True)


def save(stage: str, payload) -> str:
    os.makedirs(SWEEP, exist_ok=True)
    path = os.path.join(SWEEP, f"subgrid_{stage}.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=float)
    log(f"checkpoint -> {os.path.relpath(path, _ROOT)}")
    return path


# --------------------------------------------------------------------------- #
# The controlled problem: the shipped OXYGEN discretisation on a 1-D strip.
# --------------------------------------------------------------------------- #
def _slab_operator(nx: int, dx: float, reaction: np.ndarray,
                   d: float = D_O2) -> sp.csr_matrix:
    """``-D*L + diag(reaction)`` using the SHIPPED ``_build_laplacian``."""
    return (-d * _build_laplacian(1, nx, dx)) + sp.diags(reaction, format="csr")


def slab_vasc(a_um: float, dx: float, painting: str = "area",
              lx: float = LX) -> np.ndarray:
    """Vessel fraction per coarse cell for a source of half-width ``a_um``.

    ``'area'`` is the fractional-area painting (``geometry='subgrid'``);
    ``'binary'`` is the shipped centre-in-radius painting, which is 0 or 1.
    """
    nx = int(round(lx / dx))
    left = np.arange(nx) * dx
    right = left + dx
    if painting == "area":
        return np.clip((np.minimum(right, a_um) - left) / dx, 0.0, 1.0)
    if painting == "binary":
        return ((left + right) * 0.5 < a_um).astype(np.float64)
    raise ValueError(painting)


def slab_solve(*, supply: float, q: float, a_um: float, dx: float,
               painting: str = "area", beta: Optional[np.ndarray] = None,
               d: float = D_O2, lx: float = LX) -> np.ndarray:
    """Shipped oxygen solve on the slab, optionally with the closure applied."""
    nx = int(round(lx / dx))
    v = slab_vasc(a_um, dx, painting, lx=lx)
    src = supply * v
    b = np.ones(nx) if beta is None else np.asarray(beta, dtype=np.float64)
    o = spsolve(_slab_operator(nx, dx, src * b + q, d).tocsc(), src)
    return np.clip(np.asarray(o, dtype=np.float64), 0.0, 1.0)


def _simpson_mean(f, x0: float, x1: float, sub: int = 512) -> float:
    t = np.linspace(x0, x1, sub + 1)
    w = np.ones(sub + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return float(((x1 - x0) / sub / 3.0) * np.dot(w, f(t)) / (x1 - x0))


def _scores(o: np.ndarray) -> Dict[str, float]:
    """The quantities a verdict is actually read off."""
    lo, hi = BAND_EDGES
    f_hyp = float((o < lo).mean())
    f_qui = float(((o >= lo) & (o < hi)).mean())
    f_via = float((o >= hi).mean())
    return {
        "mean": float(o.mean()),
        "median": float(np.median(o)),
        "hypoxic_frac": f_hyp,
        "quiescent_frac": f_qui,
        "f_state": (f_hyp * BAND_MULT[0] + f_qui * BAND_MULT[1]
                    + f_via * BAND_MULT[2]),
    }


def _exact_scores(*, supply: float, q: float, a_um: float, dx: float,
                  lx: float = LX) -> Tuple[Dict[str, float], np.ndarray]:
    ref = slab_exact_reference(d=D_O2, k_in=supply + q, k_out=q, s_in=supply,
                               a=a_um, lx=lx, dx=dx)
    return _scores(ref["cell_avg"]), ref["cell_avg"]


# --------------------------------------------------------------------------- #
# STAGE identity -- what is the closure a function OF, and what does it output?
# --------------------------------------------------------------------------- #
def stage_identity() -> Dict:
    """Decide the closure's form by measurement, not by preference.

    Three oracles, each given the EXACT answer and asked how much of the error it
    can remove if it were perfect:

      * ``beta``  -- reweight the source term's diagonal by the exact
        ``mean_vessel(O)/mean_cell(O)``.
      * ``flux``  -- replace every face diffusivity by the value that makes the
        discrete face flux equal the exact flux, given the exact cell averages.
      * both.

    If ``beta`` alone is near-complete, the closure is a scalar per cell and
    nothing more is needed. If it is not, the closure has to carry a flux term
    too and is a harder object.
    """
    out: Dict = {"note": "oracle decomposition; every beta here is EXACT, taken "
                         "from the closed form, not modelled",
                 "arms": []}
    for tag, supply, cons in (("joint_lo", JOINT_SUPPLY, JOINT_CONSUMPTION_LO),
                              ("joint_hi", JOINT_SUPPLY, JOINT_CONSUMPTION_HI),
                              ("shipped", SHIPPED_SUPPLY, SHIPPED_CONSUMPTION)):
        q = cons * 1.0                       # rho = 1 (the section's p90 tile)
        arm = {"tag": tag, "supply": supply, "consumption": cons, "rho": 1.0,
               "lambda_source_um": math.sqrt(D_O2 / supply),
               "lambda_sink_um": math.sqrt(D_O2 / q),
               "a_um": VESSEL_R, "points": []}
        for dx in DX_ALL:
            if abs(round(LX / dx) * dx - LX) > 1e-9:
                continue
            nx = int(round(LX / dx))
            ref = slab_exact_reference(d=D_O2, k_in=supply + q, k_out=q,
                                       s_in=supply, a=VESSEL_R, lx=LX, dx=dx)
            exact = ref["cell_avg"]

            def fex(t, supply=supply, q=q):
                return slab_exact(t, d=D_O2, k_in=supply + q, k_out=q,
                                  s_in=supply, a=VESSEL_R, lx=LX)

            v = slab_vasc(VESSEL_R, dx, "area")
            # ORACLE beta: exact mean of O over the vessel support in each cell.
            beta = np.ones(nx)
            left = np.arange(nx) * dx
            for j in np.nonzero(v > 0)[0]:
                hi = min(left[j] + dx, VESSEL_R)
                if hi > left[j]:
                    beta[j] = _simpson_mean(fex, left[j], hi) / exact[j]
            # ORACLE flux: face diffusivity reproducing the exact face flux.
            faces = (np.arange(nx - 1) + 1.0) * dx
            eps = 1e-4
            grad = (fex(faces + eps) - fex(faces - eps)) / (2 * eps)
            dif = (exact[1:] - exact[:-1]) / dx
            with np.errstate(divide="ignore", invalid="ignore"):
                dface = np.where(np.abs(dif) > 1e-300, D_O2 * grad / dif, D_O2)
            dface = np.where(np.isfinite(dface), dface, D_O2)

            def solve_face(bet):
                rows, cols, vals = [], [], []
                diag = np.zeros(nx)
                i = np.arange(nx - 1)
                f = dface / dx ** 2
                rows += [i, i + 1]
                cols += [i + 1, i]
                vals += [f, f]
                np.add.at(diag, i, -f)
                np.add.at(diag, i + 1, -f)
                rows.append(np.arange(nx))
                cols.append(np.arange(nx))
                vals.append(diag)
                lap = sp.coo_matrix(
                    (np.concatenate(vals),
                     (np.concatenate(rows), np.concatenate(cols))),
                    shape=(nx, nx)).tocsr()
                src = supply * v
                a = (-lap) + sp.diags(src * bet + q, format="csr")
                return np.clip(spsolve(a.tocsc(), src), 0.0, 1.0)

            base = slab_solve(supply=supply, q=q, a_um=VESSEL_R, dx=dx)
            o_beta = slab_solve(supply=supply, q=q, a_um=VESSEL_R, dx=dx,
                                beta=beta)
            o_flux = solve_face(np.ones(nx))
            o_both = solve_face(beta)

            def err(o):
                return float((o.mean() - exact.mean()) / exact.mean())

            def linf(o):
                return float(np.max(np.abs(o - exact)))

            arm["points"].append({
                "dx_um": dx,
                "cells_per_lambda_sink": math.sqrt(D_O2 / q) / dx,
                "beta_at_vessel_cell": float(beta[0]),
                "rel_err_mean_base": err(base),
                "rel_err_mean_oracle_beta": err(o_beta),
                "rel_err_mean_oracle_flux": err(o_flux),
                "rel_err_mean_oracle_both": err(o_both),
                "linf_base": linf(base),
                "linf_oracle_both": linf(o_both),
                # the discrete conservation identity, checked not asserted
                "conservation_residual_base": float(abs(
                    ((supply * v + q) * base).sum() - (supply * v).sum())
                    / max((supply * v).sum(), 1e-300)),
            })
        worst_both = max(abs(p["rel_err_mean_oracle_both"]) for p in arm["points"])
        arm["worst_abs_rel_err_oracle_both"] = worst_both
        arm["worst_linf_oracle_both"] = max(p["linf_oracle_both"]
                                            for p in arm["points"])
        out["arms"].append(arm)
        log(f"identity {tag}: oracle(beta+flux) worst |rel err of mean| = "
            f"{worst_both:.3e}, worst Linf = {arm['worst_linf_oracle_both']:.3e}")
    return out


# --------------------------------------------------------------------------- #
# STAGE sweep -- generate the controlled training/held-out data
# --------------------------------------------------------------------------- #
#: Sink lengths spanning what actually occurs: the jointly calibrated 10-14 um at
#: a dense tile up to the shipped 316 um, and down to 3 um.
LAMBDA_SINK_UM: Tuple[float, ...] = (3.0, 5.0, 8.0, 10.0, 14.0, 20.0, 32.0,
                                     50.0, 90.0, 160.0, 320.0)
#: Source equilibration lengths: the jointly calibrated 3.2 um to the shipped 200.
LAMBDA_SOURCE_UM: Tuple[float, ...] = (3.2, 6.0, 12.0, 30.0, 80.0, 200.0)
#: Vessel half-widths. 7.5 um is the shipped VESSEL_RADIUS_UM.
A_UM: Tuple[float, ...] = (3.0, 7.5, 15.0)


def stage_sweep() -> Dict:
    """Every controlled point: exact cell averages vs the coarse solve."""
    out: Dict = {
        "D_um2_s": D_O2,
        "half_spacing_um": list(HALF_SPACING_UM),
        "lambda_sink_um": list(LAMBDA_SINK_UM),
        "lambda_source_um": list(LAMBDA_SOURCE_UM),
        "a_um": list(A_UM),
        "dx_train_um": list(DX_TRAIN), "dx_heldout_um": list(DX_HELDOUT),
        "band_edges": list(BAND_EDGES),
        "min_cells_per_domain": 4,
        "records": [],
    }
    n = 0
    for lam_q in LAMBDA_SINK_UM:
        q = D_O2 / lam_q ** 2
        for lam_s in LAMBDA_SOURCE_UM:
            supply = D_O2 / lam_s ** 2
            for a_um in A_UM:
                for lx in HALF_SPACING_UM:
                    if a_um >= lx * 0.5:
                        continue     # the vessel must not fill the half-cell
                    for dx in DX_ALL:
                        if abs(round(lx / dx) * dx - lx) > 1e-9:
                            continue
                        nx = int(round(lx / dx))
                        if nx < 4:
                            continue
                        ex, exact = _exact_scores(supply=supply, q=q, a_um=a_um,
                                                  dx=dx, lx=lx)
                        v = slab_vasc(a_um, dx, "area", lx=lx)
                        base = slab_solve(supply=supply, q=q, a_um=a_um, dx=dx,
                                          lx=lx)
                        b_open = source_reweight(v, dx, supply,
                                                 np.full_like(v, q), D_O2,
                                                 model="open")
                        b_nofl = source_reweight(v, dx, supply,
                                                 np.full_like(v, q), D_O2,
                                                 model="noflux")
                        o_open = slab_solve(supply=supply, q=q, a_um=a_um,
                                            dx=dx, beta=b_open, lx=lx)
                        o_nofl = slab_solve(supply=supply, q=q, a_um=a_um,
                                            dx=dx, beta=b_nofl, lx=lx)
                        rec = {
                            "lambda_sink_um": lam_q, "lambda_source_um": lam_s,
                            "a_um": a_um, "half_spacing_um": lx, "dx_um": dx,
                            "n_cells": nx,
                            "split": "train" if dx in DX_TRAIN else "heldout",
                            "cells_per_lambda_sink": lam_q / dx,
                            "cells_per_lambda_source": lam_s / dx,
                            "beta_max": float(b_open.max()),
                            "exact": ex,
                        }
                        for tag, o in (("base", base), ("open", o_open),
                                       ("noflux", o_nofl)):
                            sc = _scores(o)
                            rec[tag] = sc
                            rec[f"{tag}_rel_err_mean"] = (
                                (sc["mean"] - ex["mean"])
                                / max(ex["mean"], 1e-300))
                            rec[f"{tag}_abs_err_f_state"] = (sc["f_state"]
                                                             - ex["f_state"])
                            rec[f"{tag}_abs_err_hypoxic"] = (
                                sc["hypoxic_frac"] - ex["hypoxic_frac"])
                            rec[f"{tag}_rel_l2"] = float(
                                np.sqrt(np.mean((o - exact) ** 2))
                                / max(np.sqrt(np.mean(exact ** 2)), 1e-300))
                        out["records"].append(rec)
                        n += 1
        log(f"sweep lambda_sink = {lam_q:g} um done ({n} points)")
    # How many ladders have a NON-degenerate banded score, i.e. the exact field
    # actually straddles a band edge somewhere on the ladder. A closure cannot be
    # tested on a saturated statistic and this reports how many points can test it.
    nd = sum(1 for r in out["records"]
             if 0.0 < r["exact"]["hypoxic_frac"] < 1.0)
    out["n_points"] = n
    out["n_points_nondegenerate_band"] = nd
    log(f"sweep: {n} controlled points, {nd} with a non-degenerate banded score")
    return out


# --------------------------------------------------------------------------- #
# STAGE analytic -- the baseline the brief names, then the closure, out of sample
# --------------------------------------------------------------------------- #
def _fit_power_law(recs: Sequence[Dict]) -> Dict:
    """``|rel err of mean| = C * dx^p * lambda_sink^q``, the brief's baseline.

    Fitted on the TRAIN dx only, exactly as ``boundary_layer_resolution.md`` §4.2
    fits its own law, and then used as a multiplicative correction to the mean.
    """
    xs, ys, zs = [], [], []
    for r in recs:
        e = r["base_rel_err_mean"]
        if not np.isfinite(e) or abs(e) < 1e-14:
            continue
        xs.append(math.log(r["dx_um"]))
        ys.append(math.log(r["lambda_sink_um"]))
        zs.append(math.log(abs(e)))
    a = np.column_stack([np.ones(len(xs)), xs, ys])
    coef, *_ = np.linalg.lstsq(a, np.asarray(zs), rcond=None)
    pred = a @ coef
    ss_res = float(((np.asarray(zs) - pred) ** 2).sum())
    ss_tot = float(((np.asarray(zs) - np.mean(zs)) ** 2).sum())
    return {"n_points": len(xs), "C": float(math.exp(coef[0])),
            "p_dx": float(coef[1]), "q_lambda": float(coef[2]),
            "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
            "form": "|rel_err_mean| = C * dx^p * lambda_sink^q"}


def _grid_span(recs: Sequence[Dict], key: str) -> float:
    """Max-min of a quantity across the dx ladder -- the grid dependence."""
    vals = [r[key] for r in recs]
    return float(max(vals) - min(vals)) if vals else float("nan")


def stage_analytic(sweep: Dict) -> Dict:
    recs = sweep["records"]
    train = [r for r in recs if r["split"] == "train"]
    held = [r for r in recs if r["split"] == "heldout"]

    out: Dict = {"n_train": len(train), "n_heldout": len(held)}

    # --- the brief's named baseline: a two-parameter power law in (dx, lambda) --
    fit = _fit_power_law(train)
    out["power_law_baseline"] = fit
    # Applied as a correction it needs a SIGN, and the measured sign is not
    # constant, so the honest scoring is to give it the exact sign (a free extra
    # bit of oracle information) and see whether it still helps.
    corr_err = []
    for r in held:
        e = r["base_rel_err_mean"]
        mag = fit["C"] * r["dx_um"] ** fit["p_dx"] \
            * r["lambda_sink_um"] ** fit["q_lambda"]
        corr_err.append(abs(e - math.copysign(mag, e)))
    out["power_law_baseline_heldout"] = {
        "median_abs_rel_err_after": float(np.median(corr_err)),
        "p90_abs_rel_err_after": float(np.percentile(corr_err, 90)),
        "median_abs_rel_err_before": float(np.median(
            [abs(r["base_rel_err_mean"]) for r in held])),
        "note": "given the EXACT sign of the error as free oracle information; "
                "it still only corrects a scalar and cannot touch a banded "
                "fraction, which is a functional of the FIELD",
        "can_correct_banded_fraction": False,
    }

    # --- the analytic closure, scored on train and held-out separately ---------
    for split, rs in (("train", train), ("heldout", held)):
        blk = {}
        for tag in ("base", "open", "noflux"):
            blk[tag] = {
                "median_abs_rel_err_mean": float(np.median(
                    [abs(r[f"{tag}_rel_err_mean"]) for r in rs])),
                "p90_abs_rel_err_mean": float(np.percentile(
                    [abs(r[f"{tag}_rel_err_mean"]) for r in rs], 90)),
                "max_abs_rel_err_mean": float(np.max(
                    [abs(r[f"{tag}_rel_err_mean"]) for r in rs])),
                "median_abs_err_f_state": float(np.median(
                    [abs(r[f"{tag}_abs_err_f_state"]) for r in rs])),
                "p90_abs_err_f_state": float(np.percentile(
                    [abs(r[f"{tag}_abs_err_f_state"]) for r in rs], 90)),
                "median_rel_l2": float(np.median(
                    [r[f"{tag}_rel_l2"] for r in rs])),
            }
        out[f"scores_{split}"] = blk

    # --- the falsification criterion: GRID DEPENDENCE of a scored quantity ----
    # For each (lambda_sink, lambda_source, a) configuration, the span of the
    # scored quantity across the whole dx ladder. That is the thing the closure
    # has to shrink, and it is scored on the HELD-OUT grids as well as all.
    ladders: List[Dict] = []
    keys = sorted({(r["lambda_sink_um"], r["lambda_source_um"], r["a_um"],
                    r["half_spacing_um"]) for r in recs})
    for lam_q, lam_s, a_um, lx in keys:
        grp = [r for r in recs
               if (r["lambda_sink_um"], r["lambda_source_um"], r["a_um"],
                   r["half_spacing_um"]) == (lam_q, lam_s, a_um, lx)]
        if len(grp) < 3:
            continue
        grp_h = [r for r in grp if r["split"] == "heldout"]
        row = {"lambda_sink_um": lam_q, "lambda_source_um": lam_s, "a_um": a_um,
               "half_spacing_um": lx, "n_rungs": len(grp),
               "band_nondegenerate": bool(
                   any(0.0 < r["exact"]["hypoxic_frac"] < 1.0 for r in grp))}
        for tag in ("base", "open", "noflux"):
            row[f"{tag}_span_mean"] = _grid_span(grp, f"{tag}_rel_err_mean")
            row[f"{tag}_span_f_state"] = _grid_span(
                [{"k": r[tag]["f_state"]} for r in grp], "k")
            row[f"{tag}_span_f_state_heldout"] = _grid_span(
                [{"k": r[tag]["f_state"]} for r in grp_h], "k")
            row[f"{tag}_span_hypoxic"] = _grid_span(
                [{"k": r[tag]["hypoxic_frac"]} for r in grp], "k")
        row["exact_f_state"] = grp[0]["exact"]["f_state"]
        ladders.append(row)
    out["ladders"] = ladders

    nd = [l for l in ladders if l["band_nondegenerate"]]
    out["n_ladders"] = len(ladders)
    out["n_ladders_band_nondegenerate"] = len(nd)

    def _agg(field: str, subset=None) -> Dict:
        rows = ladders if subset is None else subset
        return {tag: {
            "median": float(np.median([l[f"{tag}_{field}"] for l in rows])),
            "p90": float(np.percentile([l[f"{tag}_{field}"] for l in rows], 90)),
            "max": float(np.max([l[f"{tag}_{field}"] for l in rows])),
        } for tag in ("base", "open", "noflux")}

    out["grid_span_summary"] = {
        "rel_err_mean": _agg("span_mean"),
        "f_state_all_ladders": _agg("span_f_state"),
        "f_state_nondegenerate_only": _agg("span_f_state", nd),
        "f_state_heldout_grids_only": _agg("span_f_state_heldout", nd),
        "hypoxic_frac_nondegenerate_only": _agg("span_hypoxic", nd),
    }

    # Fraction of ladders on which the closure strictly reduces the F_state span.
    for tag in ("open", "noflux"):
        wins = sum(1 for l in nd
                   if l[f"{tag}_span_f_state"] < l["base_span_f_state"] - 1e-12)
        ties = sum(1 for l in nd
                   if abs(l[f"{tag}_span_f_state"] - l["base_span_f_state"]) <= 1e-12)
        out[f"{tag}_f_state_span_wins"] = {
            "reduced": wins, "unchanged": ties,
            "increased": len(nd) - wins - ties, "n": len(nd),
            "scope": "ladders with a non-degenerate banded score"}
        rw = sum(1 for l in nd if l[f"{tag}_span_mean"] < l["base_span_mean"] - 1e-12)
        out[f"{tag}_mean_span_wins"] = {
            "reduced": rw, "n": len(nd)}

    log("analytic: median |rel err of mean| on HELD-OUT grids  base "
        f"{out['scores_heldout']['base']['median_abs_rel_err_mean']:.4f} -> open "
        f"{out['scores_heldout']['open']['median_abs_rel_err_mean']:.4f}")
    g = out["grid_span_summary"]["f_state_nondegenerate_only"]
    log(f"analytic: {out['n_ladders_band_nondegenerate']}/{out['n_ladders']} "
        f"ladders have a non-degenerate banded score; median F_state grid SPAN "
        f"base {g['base']['median']:.4f} -> open {g['open']['median']:.4f}")
    return out


# --------------------------------------------------------------------------- #
# STAGE transfer -- the honest question: does it move on real tissue?
# --------------------------------------------------------------------------- #
#: The ladder from results/grid_convergence_gbm.md, plus a fine probe. The
#: oxygen solve is steady state, so no dt and no CFL enters any of this.
DX_LADDER: Tuple[float, ...] = (30.0, 18.0, 12.0, 8.0)
DX_FINE: Tuple[float, ...] = (6.0, 5.0, 4.0)


def _tissue_arm(micro, supply: float, consumption: float, *,
                painting: str, closure: bool, beta_cap: Optional[float] = None
                ) -> Dict:
    """One (section, oxygen configuration, painting, closure) measurement.

    ``painting`` selects the vascular field the closure is measured on:

    ``'binary'``   the SHIPPED field (``geometry='binary'``), exactly 0 or 1.
    ``'subgrid'``  the per-vessel-OBJECT fractional-area field
                   (``geometry='subgrid'``, ``results/vessel_geometry_fix.md``).
                   This is what §5 of ``results/subgrid_closure.md`` used.
    ``'area'``     the AREA-CONSERVING rasterisation of the SHIPPED discs
                   (``geometry='area'``, ``results/vessel_source_conservation.md``)
                   — same objects, same radius, same per-endothelial-cell
                   weighting as ``'binary'``; only the measurement of the disc
                   changes. This is the promotion candidate, and it is a
                   DIFFERENT continuous field from ``'subgrid'``: it deposits
                   ~744k um^2 against ~347k. Added by
                   ``results/area_conservation_regeneration.md``.
    """
    import dataclasses

    from tmesim.cells.agents import update_cell_states
    from tmesim.core import CellState, CellType, Microenvironment, OxygenConfig, PDConfig
    from tmesim.physics.oxygen import REFERENCE_PO2_MMHG, solve_oxygen
    from tmesim.physics.subgrid import solve_oxygen_subgrid
    from tmesim.reconstruct.real_fields import vascular_density_from_cells

    cfg = dataclasses.replace(OxygenConfig(), supply=supply,
                              consumption=consumption)
    if painting == "binary":
        vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    elif painting == "area":
        # radius_mode is IGNORED by geometry='area' (it paints the shipped discs
        # at radius_um), so this is the production call with one keyword changed.
        vasc = vascular_density_from_cells(
            micro.domain, micro.cells, mode="discrete", geometry="area",
            radius_um=VESSEL_R)
    elif painting == "subgrid":
        vasc = vascular_density_from_cells(
            micro.domain, micro.cells, mode="discrete", geometry="subgrid",
            radius_mode="constant", radius_um=VESSEL_R)
    else:
        raise ValueError(f"unknown painting {painting!r}")
    work = Microenvironment(
        domain=micro.domain,
        fields=dataclasses.replace(micro.fields, vascular_density=vasc),
        cells=micro.cells.copy(), patient_id=micro.patient_id,
        metadata=micro.metadata)

    beta_stats = None
    if closure:
        v_for_beta = vasc
        if beta_cap is not None:
            # Restrict the closure to the range it was VALIDATED over on the
            # controlled problem. beta -> 1/v in the thin-source limit, so
            # capping beta is exactly capping 1/v; cells below the floor keep
            # beta = 1, i.e. the closure declines to act rather than
            # extrapolating. This tests the diagnosis, it does not tune it.
            from tmesim.physics.oxygen import _tumor_density as _td
            b_full = source_reweight(vasc, micro.domain.dx, supply,
                                     float(consumption) * _td(work), 2000.0,
                                     model="open")
            v_for_beta = np.where(b_full > beta_cap, 1.0, vasc)
            beta_stats = {"beta_cap": float(beta_cap),
                          "frac_cells_declined": float(
                              (b_full[(vasc > 0) & (vasc < 1)] > beta_cap).mean())
                          if np.any((vasc > 0) & (vasc < 1)) else 0.0}
        oxy = solve_oxygen_subgrid(work, cfg, enabled=True, model="open",
                                   vascular_fraction=v_for_beta)
    else:
        oxy = solve_oxygen(work, cfg)

    pdc = PDConfig()
    mult = {"HYPOXIC": pdc.hypoxia_resistance,
            "QUIESCENT": pdc.quiescent_resistance, "VIABLE": 1.0}
    scored = Microenvironment(
        domain=work.domain, fields=work.fields, cells=micro.cells.copy(),
        patient_id=micro.patient_id, metadata=micro.metadata)
    update_cell_states(scored, oxy, OxygenConfig())
    st = np.asarray(scored.cells.state)
    tum = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
    n_tum = int(tum.sum())
    fr = {s.name: float((st[tum] == int(s)).sum()) / n_tum
          for s in (CellState.HYPOXIC, CellState.QUIESCENT, CellState.VIABLE)}
    col, row = micro.domain.to_index(np.asarray(micro.cells.x)[tum],
                                     np.asarray(micro.cells.y)[tum])
    o_tum = oxy[row, col]
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    return {
        "f_state": float(sum(fr[k] * mult[k] for k in fr)),
        "band_fractions": fr,
        "mean_o2_tumour": float(o_tum.mean()),
        "median_o2_tumour_mmhg": float(np.median(o_tum) * REFERENCE_PO2_MMHG),
        "mean_o2_tissue": float(oxy[tissue].mean()),
        "source_integral": float((vasc[tissue]).sum()) * float(micro.domain.dx) ** 2,
        "n_source_cells": int((vasc[tissue] > 0).sum()),
        "vasc_is_binary": bool(np.all((vasc == 0) | (vasc == 1))),
        "beta_stats": beta_stats,
        "o2_per_tumour_cell": o_tum,           # stripped before JSON
    }


def stage_transfer(fine: bool = True) -> Dict:
    """Real GBM section, dx ladder, with and without the closure."""
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.physics.subgrid import source_reweight as _sr

    rungs = list(DX_LADDER) + (list(DX_FINE) if fine else [])
    arms = (
        ("shipped", SHIPPED_SUPPLY, SHIPPED_CONSUMPTION),
        ("joint_lo", JOINT_SUPPLY, JOINT_CONSUMPTION_LO),
        ("joint_hi", JOINT_SUPPLY, JOINT_CONSUMPTION_HI),
    )
    out: Dict = {"section": "data_cache/xenium_gbm", "dx_ladder": rungs,
                 "arms": {}, "per_cell": {}}
    ck = os.path.join(SWEEP, "subgrid_transfer_partial.json")

    for dx in rungs:
        t0 = time.time()
        micro = reconstruct_real(source="xenium", dx_um=float(dx))
        log(f"transfer dx={dx:g}: reconstructed {micro.domain.shape} in "
            f"{time.time() - t0:.1f}s")
        for tag, sup, con in arms:
            for painting in ("binary", "subgrid"):
                for closure in (False, True):
                    if painting == "binary" and closure:
                        # beta == 1 identically on a binary field; run it anyway
                        # ONCE at the published grid to demonstrate the no-op
                        # rather than assert it.
                        if dx != 12.0:
                            continue
                    key = f"{tag}|{painting}|{'closed' if closure else 'raw'}"
                    r = _tissue_arm(micro, sup, con, painting=painting,
                                    closure=closure)
                    percell = r.pop("o2_per_tumour_cell")
                    out["per_cell"].setdefault(key, {})[str(dx)] = percell
                    r["dx_um"] = float(dx)
                    out["arms"].setdefault(key, []).append(r)
        # beta diagnostics on the subgrid field at this dx
        with open(ck, "w") as fh:
            json.dump({k: v for k, v in out.items() if k != "per_cell"}, fh,
                      indent=2, default=float)
        log(f"transfer dx={dx:g}: done ({time.time() - t0:.1f}s total)")

    # --- spans, and the band-change count that oxygen_joint_grid_test.md uses --
    summary = {}
    for key, rows in out["arms"].items():
        rows.sort(key=lambda r: -r["dx_um"])
        lad = [r for r in rows if r["dx_um"] in DX_LADDER]
        if len(lad) < 2:
            continue
        fs = [r["f_state"] for r in lad]
        allfs = [r["f_state"] for r in rows]
        d = np.diff(fs)
        summary[key] = {
            "f_state_by_dx": {str(r["dx_um"]): r["f_state"] for r in rows},
            "f_state_span_ladder": float(max(fs) - min(fs)),
            "f_state_span_all_rungs": float(max(allfs) - min(allfs)),
            "monotone_in_dx": bool(np.all(d > 0) or np.all(d < 0)),
            "n_direction_reversals": int(np.sum(np.sign(d[1:])
                                                != np.sign(d[:-1]))),
            "f_state_at_dx12": next((r["f_state"] for r in rows
                                     if r["dx_um"] == 12.0), None),
            "f_state_finest": rows[-1]["f_state"],
            "err_at_dx12_vs_finest": (
                next((r["f_state"] for r in rows if r["dx_um"] == 12.0), np.nan)
                - rows[-1]["f_state"]),
            "vasc_is_binary": lad[0]["vasc_is_binary"],
        }
        # cells changing band across the ladder
        pc = out["per_cell"].get(key, {})
        have = [str(x) for x in DX_LADDER if str(x) in pc]
        if len(have) == len(DX_LADDER):
            bands = np.stack([np.digitize(pc[h], BAND_EDGES) for h in have])
            changed = float((bands != bands[0]).any(axis=0).mean())
            summary[key]["frac_tumour_cells_changing_band"] = changed
            summary[key]["median_per_cell_grid_sigma"] = float(
                np.median(np.std(np.stack([pc[h] for h in have]), axis=0)))
    out["summary"] = summary
    out.pop("per_cell")
    _ = _sr
    for key in sorted(summary):
        s = summary[key]
        log(f"transfer {key:28s} F_state span(ladder) {s['f_state_span_ladder']:.4f}"
            f"  dx12 {s['f_state_at_dx12']}  finest {s['f_state_finest']:.4f}")
    return out


#: The largest ``beta`` reached ANYWHERE in the controlled sweep. Beyond this the
#: closure is being extrapolated outside everything it was validated on.
BETA_VALIDATED_MAX = 5.677


def stage_clip() -> Dict:
    """Does restricting the closure to its VALIDATED range make it transfer?

    ``beta`` degenerates to ``1/v`` for a thin saturated source, so on tissue --
    where the fractional-area vascular field has a median ``v`` of 0.229 and a
    1st percentile of 0.0069 -- it reaches values two orders of magnitude beyond
    anything the controlled problem produced. This stage caps it at the
    controlled sweep's own maximum and lets the closure DECLINE to act
    (``beta = 1``) on cells outside that range. It tests the diagnosis in §5;
    it does not tune the physics, and the cap is the measured
    ``BETA_VALIDATED_MAX``, not a fitted number.
    """
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    arms = (("shipped", SHIPPED_SUPPLY, SHIPPED_CONSUMPTION),
            ("joint_lo", JOINT_SUPPLY, JOINT_CONSUMPTION_LO),
            ("joint_hi", JOINT_SUPPLY, JOINT_CONSUMPTION_HI))
    out: Dict = {"beta_cap": BETA_VALIDATED_MAX, "dx_ladder": list(DX_LADDER),
                 "arms": {}}
    for dx in DX_LADDER:
        t0 = time.time()
        micro = reconstruct_real(source="xenium", dx_um=float(dx))
        for tag, sup, con in arms:
            r = _tissue_arm(micro, sup, con, painting="subgrid", closure=True,
                            beta_cap=BETA_VALIDATED_MAX)
            r.pop("o2_per_tumour_cell", None)
            r["dx_um"] = float(dx)
            out["arms"].setdefault(f"{tag}|subgrid|clipped", []).append(r)
        log(f"clip dx={dx:g}: done ({time.time() - t0:.1f}s)")
    summary = {}
    for k, rows in out["arms"].items():
        rows.sort(key=lambda r: -r["dx_um"])
        fs = [r["f_state"] for r in rows]
        summary[k] = {"f_state_by_dx": {str(r["dx_um"]): r["f_state"]
                                        for r in rows},
                      "f_state_span_ladder": float(max(fs) - min(fs)),
                      "frac_cells_declined": rows[0]["beta_stats"][
                          "frac_cells_declined"]}
        log(f"clip {k:28s} span {summary[k]['f_state_span_ladder']:.4f}  "
            + " ".join(f"{v:.4f}" for v in fs))
    out["summary"] = summary
    return out


# --------------------------------------------------------------------------- #
def main(argv: Optional[Sequence[str]] = None) -> int:
    stages = list(argv if argv is not None else sys.argv[1:])
    if not stages:
        stages = ["identity", "sweep", "analytic", "transfer"]
    os.makedirs(SWEEP, exist_ok=True)

    sweep = None
    if "identity" in stages:
        save("identity", stage_identity())
    if "sweep" in stages:
        sweep = stage_sweep()
        save("sweep", sweep)
    if "analytic" in stages:
        if sweep is None:
            with open(os.path.join(SWEEP, "subgrid_sweep.json")) as fh:
                sweep = json.load(fh)
        save("analytic", stage_analytic(sweep))
    if "transfer" in stages:
        save("transfer", stage_transfer())
    if "clip" in stages:
        save("clip", stage_clip())
    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
