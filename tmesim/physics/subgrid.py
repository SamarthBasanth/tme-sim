"""Opt-in sub-grid closure for the perivascular boundary layer of the O2 solve.

**NOTHING IN THIS MODULE IS ON BY DEFAULT AND NOTHING ELSE IN ``tmesim/``
IMPORTS IT.** ``solve_oxygen_subgrid(..., enabled=False)`` is a bit-identical
passthrough to :func:`tmesim.physics.oxygen.solve_oxygen`. It exists so
``experiments/subgrid_closure.py`` can measure a closure rather than assert one,
and so anyone who wants to re-measure it can import one auditable function.

THE PROBLEM
-----------
``results/oxygen_joint_grid_test.md`` §4 established that with the jointly
calibrated oxygen constants (``supply = 194``, ``consumption = 10-20``) the
source equilibration length ``sqrt(D/supply)`` falls to **3.2 um** and the sink
diffusion length ``sqrt(D/(consumption*rho))`` to **10-14 um**, against a
production mesh of ``dx = 12 um`` -- **0.47 to 1.77 mesh cells per length
scale** -- and concluded that "a physically correct vessel is a boundary-layer
object this project's meshes cannot represent".

WHAT IS ACTUALLY WRONG, DERIVED NOT ASSUMED
--------------------------------------------
The shipped steady solve (``physics/oxygen.solve_oxygen``, ``sink_form="linear"``)
is::

    (-D * L + diag(supply*v + q)) @ O = supply*v          q = consumption*rho

``L`` is the 5-point Laplacian with **zero row sums**, so summing the discrete
system over cells gives, exactly and at every ``dx``::

    sum_j (supply*v_j + q_j) * O_j  =  sum_j supply*v_j

and the continuum satisfies the same identity with integrals. Subtracting the
two, the **entire** domain-mean discretisation error is the mismatch between

    supply * v_j * O_j            (discrete: cell-average O in the source term)
    integral of supply * v * O    (continuum: O weighted by where the vessel IS)

When the vessel is sub-cell, ``O`` inside it is far higher than the cell average,
so the discrete form under-weights the source's own self-shading and the solve
**over-supplies**. That is a *source-weighting* error, not a flux error, and it
is removed by ONE scalar per cell.

THE CLOSURE
-----------
Multiply the **source term's contribution to the reaction diagonal only** by

    beta_j  =  mean of O over the vessel support in cell j
               ----------------------------------------------
               mean of O over the whole of cell j

leaving the right-hand side ``supply*v`` untouched (it is already the exact cell
integral of the source). With the exact ``beta`` the coarse solve reproduces the
exact cell averages to **machine precision** -- measured, see
``results/subgrid_closure.md`` §2 -- so this is a complete parameterisation of
the error, not a heuristic.

``beta`` is then evaluated from a **closed-form local cell problem** with **zero
fitted parameters**: a 1-D two-region cell of width ``dx`` carrying a source of
half-width ``a = v*dx``, symmetric (no-flux) at the vessel centre and, in the
default ``model="open"``, matched to a semi-infinite decaying sink medium at the
outer edge::

    [0, a) :  D O'' = (supply + q) O - supply       lam_i = sqrt(D/(supply+q))
    [a, oo):  D O'' = q O                           lam_o = sqrt(D/q)

Its two cell means are elementary. ``model="noflux"`` closes the cell instead
(a Krogh-cylinder unit cell); it is retained because it is the textbook
homogenisation choice and it is **measurably worse** here, which is worth being
able to reproduce.

WHAT IT IS A FUNCTION OF
------------------------
``(v, dx, supply, q, D)`` -- i.e. the local vascular AREA FRACTION, the mesh
spacing, and the two rate constants. Every one is already a per-cell quantity the
engine holds. There is no fitted coefficient, no training set and no lookup
table.

THE ONE STRUCTURAL PRECONDITION -- READ THIS BEFORE USING IT
-------------------------------------------------------------
``beta == 1`` identically whenever ``v == 0`` or ``v == 1``: a cell that is
entirely vessel or entirely tissue has no sub-cell structure to close over. The
**shipped** vascular field is
``reconstruct.real_fields._DEFAULT_VASCULAR_MODE = "discrete"`` with
``geometry="binary"``, which is **exactly 1.0 or exactly 0.0 everywhere**. On
that field this closure is an exact no-op by construction. It requires the
fractional-area field
(``vessel_network.discrete_vascular_field(..., geometry="subgrid")``, itself
opt-in and default-off, from ``results/vessel_geometry_fix.md``). The binary
painting discards the sub-grid geometry BEFORE the solve, so no sub-grid closure
of the solve can recover it.

KNOWN FAILURE MODES (measured in ``results/subgrid_closure.md`` §3-§4)
----------------------------------------------------------------------
1. The local cell problem is **1-D Cartesian**; a real vessel is a disc. The
   radial ``1/r`` term steepens the near field, so ``beta`` is under-estimated
   for a curved boundary.
2. It closes the **source weighting only**. The face-flux (transmissibility)
   error is left in place; it is sub-dominant at ``dx >= 10 um`` but is the
   dominant residual at ``dx = 4-8 um``, where it is non-monotone in ``dx``.
3. It is derived for the ``"linear"`` sink. It is **not** valid for
   ``sink_form="michaelis_menten"``, where the local problem has no closed form.
4. It assumes a single vessel per cell. Two vessels in one coarse cell share one
   ``v`` and get one ``beta``.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from ..core import Microenvironment, OxygenConfig
from .oxygen import _build_laplacian, _tumor_density

__all__ = ["source_reweight", "solve_oxygen_subgrid", "SUBGRID_MODELS"]

#: Local cell problems available for the closed form. ``"open"`` is the default
#: and the one that measures better; see the module docstring.
SUBGRID_MODELS = ("open", "noflux")


def source_reweight(
    vascular_fraction: np.ndarray,
    dx: float,
    supply: float,
    sink: np.ndarray,
    diffusivity: float,
    model: str = "open",
) -> np.ndarray:
    """Per-cell source-reweighting factor ``beta`` (>= 1), closed form.

    Parameters
    ----------
    vascular_fraction : ndarray
        Vessel AREA FRACTION per cell, in [0, 1]. **A binary field gives
        ``beta == 1`` everywhere** -- see the module docstring.
    dx : float
        Grid spacing (um).
    supply : float
        ``OxygenConfig.supply`` (1/s).
    sink : ndarray
        ``consumption * rho_tum`` per cell (1/s).
    diffusivity : float
        ``OxygenConfig.diffusivity`` (um^2/s).
    model : {'open', 'noflux'}
        Outer boundary of the local cell problem.

    Returns
    -------
    ndarray
        ``beta``, same shape as ``vascular_fraction``, exactly 1.0 wherever
        ``v <= 0`` or ``v >= 1``, and 1.0 wherever the formula is not finite.
    """
    m = str(model).strip().lower()
    if m not in SUBGRID_MODELS:
        raise ValueError(f"unknown model {model!r}; expected one of {SUBGRID_MODELS}")

    v = np.clip(np.asarray(vascular_fraction, dtype=np.float64), 0.0, 1.0)
    d = float(diffusivity)
    dx = float(dx)
    # A cell with no sink at all has an infinite outer length; the floor keeps
    # lam_o finite and makes beta -> 1 there, which is the correct limit (a
    # sink-free cell equilibrates to O = 1 uniformly and has no sub-cell shape).
    q = np.maximum(np.asarray(sink, dtype=np.float64) + np.zeros_like(v), 1e-300)
    k_src = float(supply) * v

    a = v * dx                       # source half-width inside the cell
    b = dx - a                       # sink-only remainder
    lam_i = np.sqrt(d / np.maximum(k_src + q, 1e-300))
    lam_o = np.sqrt(d / q)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        u = np.clip(a / lam_i, 0.0, 300.0)
        w = np.clip(b / lam_o, 1e-300, 300.0)
        t = np.tanh(u)
        if m == "open":
            # outer region is semi-infinite: O ~ O_a * exp(-(x-a)/lam_o)
            denom = 1.0 + (lam_o / lam_i) * t
            o_a = (lam_o / lam_i) * t / denom
            mean_b = o_a * lam_o * (1.0 - np.exp(-w)) / np.maximum(b, 1e-300)
        else:
            # outer region closed by a no-flux face at x = dx (Krogh unit cell)
            denom = 1.0 + (lam_o / lam_i) * t / np.tanh(w)
            mean_b = (lam_o ** 2 / lam_i) * t / (np.maximum(b, 1e-300) * denom)
        mean_a = 1.0 - (lam_i / np.maximum(a, 1e-300)) * t / denom
        mean_c = (a * mean_a + b * mean_b) / dx
        beta = np.where(mean_c > 0.0, mean_a / mean_c, 1.0)

    beta = np.where(np.isfinite(beta), beta, 1.0)
    return np.where((v <= 0.0) | (v >= 1.0), 1.0, beta)


def solve_oxygen_subgrid(
    micro: Microenvironment,
    cfg: OxygenConfig,
    *,
    enabled: bool = False,
    model: str = "open",
    vascular_fraction: Optional[np.ndarray] = None,
) -> np.ndarray:
    """``solve_oxygen`` with the sub-grid source closure, OFF by default.

    ``enabled=False`` (the default) forwards to the unmodified shipped
    :func:`tmesim.physics.oxygen.solve_oxygen` and is bit-identical to it.

    ``vascular_fraction`` overrides ``micro.fields.vascular_density`` for the
    closure's ``v`` ONLY (the source term itself still uses the field on the
    microenvironment), so a binary source field can be closed against a
    fractional-area geometry without changing the source. When it is ``None``
    the field's own values are used, and on the shipped binary field that makes
    the closure an exact no-op.
    """
    from .oxygen import solve_oxygen

    if not enabled:
        return solve_oxygen(micro, cfg)

    sink_form = str(getattr(cfg, "sink_form", "linear")).strip().lower()
    if sink_form != "linear":
        raise ValueError(
            "the sub-grid closure is derived for the linear sink only; "
            f"cfg.sink_form={sink_form!r} has no closed-form local cell problem"
        )

    domain = micro.domain
    ny, nx = domain.shape
    vasc = np.clip(np.asarray(micro.fields.vascular_density, dtype=np.float64),
                   0.0, 1.0)
    if vasc.shape != (ny, nx):
        raise ValueError(f"vascular_density shape {vasc.shape} != {(ny, nx)}")
    rho_tum = _tumor_density(micro)

    d = float(cfg.diffusivity)
    supply = float(cfg.supply)
    sink = float(cfg.consumption) * rho_tum

    v_for_beta = vasc if vascular_fraction is None else np.clip(
        np.asarray(vascular_fraction, dtype=np.float64), 0.0, 1.0)
    if v_for_beta.shape != (ny, nx):
        raise ValueError("vascular_fraction shape does not match the domain")

    beta = source_reweight(v_for_beta, domain.dx, supply, sink, d, model=model)

    source = (supply * vasc).ravel()
    if not np.any(source) and not np.any(sink):
        return domain.zeros()

    # THE ONE CHANGE: the source's contribution to the reaction DIAGONAL is
    # scaled by beta. The right-hand side is untouched -- it is already the
    # exact cell integral of the source.
    reaction = source * beta.ravel() + sink.ravel()
    a = (-d * _build_laplacian(ny, nx, domain.dx)) + sp.diags(reaction,
                                                             format="csr")
    oxygen = np.asarray(spsolve(a.tocsc(), source),
                        dtype=np.float64).reshape(ny, nx)
    np.clip(oxygen, 0.0, 1.0, out=oxygen)
    return oxygen
