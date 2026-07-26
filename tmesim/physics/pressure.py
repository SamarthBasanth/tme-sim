"""Interstitial fluid pressure (IFP) and interstitial flow — the Jain effect.

Solid tumors famously develop *elevated interstitial fluid pressure*. Two
coupled pathologies cause it:

1.  **Hyperpermeable ("leaky") vessels.** Tumor microvessels have grossly
    elevated hydraulic conductivity, so plasma filtrate pours into the
    interstitium (a distributed fluid *source*).
2.  **Absent/compressed lymphatics.** The tumor interior has no functional
    lymphatic drainage, so that filtrate has nowhere to go (the fluid *sink*
    that normal tissue relies on is missing).

Fluid therefore accumulates until the interstitial pressure rises to nearly the
microvascular pressure. Across the tumor *interior* the pressure is high and
*flat* (a plateau), so its gradient — and hence the interstitial fluid velocity
`v = -K grad(p)` (Darcy's law) — is near zero. All the pressure drop is
squeezed into a thin rim at the **tumor margin**, where a steep outward gradient
drives a radially **outward** convective flow.

That outward margin flow is the crux of the **Jain effect** (Jain 1987; Baxter &
Jain 1989): it physically *opposes* the inward transport of a blood-borne drug,
washing molecules back out of the rim and starving the core. A drug-transport
model that ignores it will over-predict delivery.

Model (poroelastic-lite, steady state)
--------------------------------------
We solve a screened-Poisson (Helmholtz-type) equation for `p` (mmHg) on the 2D
grid. Combining Darcy interstitial flow with distributed Starling source/sink
terms and linearising about the interstitial pressure gives

    -div( K grad p ) + (a + d) * p = a * p_v

where, per grid node,

* `K`   — interstitial hydraulic conductivity  [um^2 / (mmHg * s)]
          (mildly reduced by dense ECM: stiff matrix resists flow).
* `a`   — vessel-wall filtration conductance = vascular_density * Lp_eff, with
          `Lp_eff` strongly ELEVATED inside the tumor (leaky vessels)
          [1 / (mmHg * s)].  This is both the source strength and, since the
          source is driven by `(p_v - p)`, part of the reaction coefficient.
* `d`   — lymphatic drainage conductance = lambda_lymph * (1 - phi_tumor);
          FUNCTIONAL in normal tissue, ~0 inside the tumor core
          [1 / (mmHg * s)].
* `p_v` — effective microvascular driving pressure (mmHg).

Behaviour that falls out of this (and is asserted in the self-test):

* **Tumor core:** `d -> 0`, `a` large  =>  `p -> a*p_v/(a+d) -> p_v`, a high
  flat plateau  =>  `grad p ~ 0`  =>  `v ~ 0`.
* **Tumor margin:** `d` and `a` swap dominance over a thin rim  =>  steep
  pressure drop  =>  large `-K grad p` pointing radially OUTWARD.
* **Normal tissue:** strong drainage `d` pulls `p -> 0`; a Dirichlet boundary
  pins the domain edge at 0 mmHg (normal tissue reference).

Numerics: 5-point finite-volume Laplacian with harmonic-mean face
conductances (handles spatially varying `K` conservatively), assembled as a
`scipy.sparse` matrix and solved with a sparse direct solver. Dirichlet p=0 on
the domain border and outside the tissue mask; the screening term `(a+d)`
provides the interior closure (no-flux is implicit where `a=d=0`).

Parameters are order-of-magnitude biological values, not patient-calibrated —
consistent with the rest of `tmesim` (see `core/config.py`). What matters for
Stage 0 is that the *qualitative* Jain effect emerges robustly.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.linalg import spsolve

from ..core.microenvironment import Microenvironment

# --- Model parameters (order-of-magnitude; see module docstring) -------------
_K0 = 1.2            # baseline interstitial hydraulic conductivity, um^2/(mmHg*s)
_K_ECM_ATTEN = 0.6   # dense ECM reduces K by up to this fraction
_PV = 25.0           # effective microvascular driving pressure, mmHg
_LP_NORMAL = 1.0e-3  # baseline vessel-wall filtration conductance, 1/(mmHg*s)
_LEAK_GAIN = 12.0    # tumor vessels are ~(1+gain)x leakier at phi_tumor=1
_LAMBDA_LYMPH = 3.0e-2  # lymphatic drainage conductance in normal tissue, 1/(mmHg*s)
_TUMOR_SIGMA = 5.0   # gaussian smoothing (grid cells) for the tumor-density field


def _tumor_field(micro: Microenvironment) -> np.ndarray:
    """Rasterise tumor cells into a smooth [0,1] tumor-occupancy field.

    Used to spatially modulate vessel leakiness (up) and lymphatic drainage
    (down) inside the tumor. Degrades gracefully to all-zeros when the
    microenvironment contains no tumor cells.
    """
    ny, nx = micro.domain.shape
    phi = np.zeros((ny, nx), dtype=np.float64)
    cells = micro.cells
    tum = cells.tumor_mask()
    if not np.any(tum):
        return phi
    col, row = micro.domain.to_index(cells.x[tum], cells.y[tum])
    np.add.at(phi, (row, col), 1.0)
    phi = ndimage.gaussian_filter(phi, sigma=_TUMOR_SIGMA, mode="nearest")
    pmax = float(phi.max())
    if pmax > 0.0:
        phi /= pmax
    return phi


def solve_ifp(micro: Microenvironment) -> Tuple[np.ndarray, np.ndarray]:
    """Solve for interstitial fluid pressure and the interstitial fluid velocity.

    Poroelastic-lite Darcy model (see module docstring): a screened-Poisson
    solve with leaky-vessel fluid sources and tumor-suppressed lymphatic
    drainage, so an elevated flat IFP plateau builds in the tumor core and the
    pressure drops steeply across the margin. This reproduces the Jain effect:
    near-zero interstitial velocity in the core and an outward convective flow
    at the tumor margin that opposes inward drug delivery.

    Parameters
    ----------
    micro : Microenvironment
        Provides the grid (`domain`), the `vascular_density` and `ecm_density`
        fields (source strength / interstitial conductivity), and the tumor
        `cells` (used to locate the leaky, poorly-drained tumor region).

    Returns
    -------
    ifp : np.ndarray, shape (ny, nx)
        Interstitial fluid pressure in mmHg (>= 0), elevated and flat in the
        tumor core, ~0 in normal tissue and on the boundary.
    velocity : np.ndarray, shape (2, ny, nx)
        Interstitial (Darcy) fluid velocity in um/s, `v = -K grad(p)`.
        `velocity[0]` is the x-component (columns), `velocity[1]` the
        y-component (rows). Near-zero in the core, points radially OUTWARD at
        the margin.
    """
    domain = micro.domain
    ny, nx = domain.shape
    dx = float(domain.dx)
    dx2 = dx * dx
    N = ny * nx

    fields = micro.fields
    vasc = np.clip(np.asarray(fields.vascular_density, dtype=np.float64), 0.0, None)
    ecm = np.asarray(fields.ecm_density, dtype=np.float64)
    ecm = np.clip(ecm, 0.0, 1.0) if ecm.shape == (ny, nx) else np.zeros((ny, nx))

    phi = _tumor_field(micro)  # [0,1] tumor occupancy

    # Spatially varying coefficients.
    K = _K0 * (1.0 - _K_ECM_ATTEN * ecm)                 # interstitial conductivity
    lp_eff = _LP_NORMAL * (1.0 + _LEAK_GAIN * phi)       # leaky tumor vessels
    a = vasc * lp_eff                                    # vessel-wall source/coupling
    d = _LAMBDA_LYMPH * (1.0 - phi)                      # lymphatic drainage (tumor->0)
    react = a + d                                        # screening (reaction) term
    rhs_src = a * _PV                                    # distributed fluid source

    # Dirichlet nodes: domain border + anything outside the tissue mask (normal
    # tissue reference, p = 0 mmHg).
    dirichlet = np.zeros((ny, nx), dtype=bool)
    dirichlet[0, :] = dirichlet[-1, :] = True
    dirichlet[:, 0] = dirichlet[:, -1] = True
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    dirichlet |= ~tissue
    interior = ~dirichlet

    kidx = np.arange(N).reshape(ny, nx)

    # Harmonic-mean face conductances (divided by dx^2 for the FV stencil).
    def _hmean(k1: np.ndarray, k2: np.ndarray) -> np.ndarray:
        return 2.0 * k1 * k2 / (k1 + k2 + 1e-30)

    face_e = np.zeros((ny, nx))            # face between (i,j) and (i,j+1)
    face_e[:, :-1] = _hmean(K[:, :-1], K[:, 1:]) / dx2
    face_w = np.zeros((ny, nx))            # face between (i,j) and (i,j-1)
    face_w[:, 1:] = face_e[:, :-1]
    face_s = np.zeros((ny, nx))            # face between (i,j) and (i+1,j)
    face_s[:-1, :] = _hmean(K[:-1, :], K[1:, :]) / dx2
    face_n = np.zeros((ny, nx))            # face between (i,j) and (i-1,j)
    face_n[1:, :] = face_s[:-1, :]

    diag = np.zeros((ny, nx))
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []

    # For each direction: add the face conductance to the source node's diagonal
    # (always — a Dirichlet neighbour acts as a p=0 ghost), and add an
    # off-diagonal coupling only when the neighbour is itself an unknown.
    for face, nbr_off in (
        (face_e, +1),      # east  (j+1)
        (face_w, -1),      # west  (j-1)
        (face_s, +nx),     # south (i+1)
        (face_n, -nx),     # north (i-1)
    ):
        # Build a mask of source nodes that actually have a neighbour in-domain.
        has = np.zeros((ny, nx), dtype=bool)
        if nbr_off == +1:
            has[:, :-1] = True
        elif nbr_off == -1:
            has[:, 1:] = True
        elif nbr_off == +nx:
            has[:-1, :] = True
        else:
            has[1:, :] = True

        src = interior & has
        diag[src] += face[src]

        # Neighbour indices for source nodes; couple only to non-Dirichlet ones.
        src_k = kidx[src]
        nbr_k = src_k + nbr_off
        nbr_interior = interior.ravel()[nbr_k]
        keep = nbr_interior
        if np.any(keep):
            rows.append(src_k[keep])
            cols.append(nbr_k[keep])
            data.append(-face[src][keep])

    # Reaction term on interior diagonals.
    diag[interior] += react[interior]

    # Interior diagonal entries.
    rows.append(kidx[interior])
    cols.append(kidx[interior])
    data.append(diag[interior])

    # Dirichlet rows: identity.
    dk = kidx[dirichlet]
    rows.append(dk)
    cols.append(dk)
    data.append(np.ones(dk.size))

    A = sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(N, N),
    ).tocsr()

    b = np.zeros(N)
    b[kidx[interior]] = rhs_src[interior]

    p = spsolve(A, b).reshape(ny, nx)
    p = np.clip(p, 0.0, None)  # pressure is non-negative (max principle)
    ifp = np.ascontiguousarray(p)

    # Darcy velocity v = -K grad(p), central differences (um/s).
    grad_y = np.gradient(ifp, dx, axis=0)
    grad_x = np.gradient(ifp, dx, axis=1)
    vx = -K * grad_x
    vy = -K * grad_y
    velocity = np.stack([vx, vy], axis=0)

    return ifp, velocity


if __name__ == "__main__":
    # Self-test: a central, high-vascular-leak tumor in otherwise normal tissue.
    # Assert the two signatures of the Jain effect.
    import numpy as np

    from ..core.domain import Domain
    from ..core.microenvironment import (
        CellPopulation,
        CellState,
        CellType,
        Microenvironment,
        TissueFields,
    )

    rng = np.random.default_rng(0)
    nx = ny = 120
    dx = 15.0
    domain = Domain(nx=nx, ny=ny, dx=dx)
    W, H = domain.extent_um
    cx, cy = W / 2.0, H / 2.0
    R_tumor = 500.0  # um

    X, Y = domain.coords_um()
    r = np.hypot(X - cx, Y - cy)

    # Vascular density: perfused everywhere, elevated (leaky) toward the core.
    vasc = 0.30 + 0.20 * np.exp(-(r ** 2) / (2 * (0.6 * R_tumor) ** 2))
    vasc = np.clip(vasc, 0.0, 1.0)
    ecm = 0.2 + 0.3 * np.exp(-(r ** 2) / (2 * (0.7 * R_tumor) ** 2))
    fields = TissueFields(ecm_density=ecm, vascular_density=vasc)

    # Tumor cells filling a central disk; a few stromal cells scattered outside.
    n_tumor = 4000
    ang = rng.uniform(0, 2 * np.pi, n_tumor)
    rad = R_tumor * np.sqrt(rng.uniform(0, 1, n_tumor))
    tx = cx + rad * np.cos(ang)
    ty = cy + rad * np.sin(ang)
    n_str = 800
    sx = rng.uniform(0, W, n_str)
    sy = rng.uniform(0, H, n_str)
    x = np.concatenate([tx, sx])
    y = np.concatenate([ty, sy])
    ctype = np.concatenate([
        np.full(n_tumor, int(CellType.TUMOR)),
        np.full(n_str, int(CellType.STROMAL)),
    ])
    state = np.full(x.size, int(CellState.VIABLE))
    sens = rng.uniform(0.2, 0.9, x.size)
    cells = CellPopulation(x=x, y=y, cell_type=ctype, state=state, sensitivity=sens)

    micro = Microenvironment(domain=domain, fields=fields, cells=cells,
                             patient_id="pressure-selftest")

    ifp, velocity = solve_ifp(micro)

    assert ifp.shape == (ny, nx), ifp.shape
    assert velocity.shape == (2, ny, nx), velocity.shape

    # (a) Core IFP >> margin IFP.
    core_mask = r < 0.35 * R_tumor
    margin_mask = (r > 0.8 * R_tumor) & (r < 1.1 * R_tumor)
    core_ifp = float(ifp[core_mask].mean())
    margin_ifp = float(ifp[margin_mask].mean())

    # (b) Radial interstitial velocity at the margin points OUTWARD.
    eps = 1e-9
    rhat_x = (X - cx) / (r + eps)
    rhat_y = (Y - cy) / (r + eps)
    v_radial = velocity[0] * rhat_x + velocity[1] * rhat_y  # >0 == outward
    margin_vrad = float(v_radial[margin_mask].mean())

    # Core interstitial velocity is near zero (flat plateau).
    core_speed = float(np.hypot(velocity[0], velocity[1])[core_mask].mean())
    margin_speed = float(np.hypot(velocity[0], velocity[1])[margin_mask].mean())

    print("=== IFP self-test (Jain effect) ===")
    print(f"core   IFP  mean = {core_ifp:8.3f} mmHg")
    print(f"margin IFP  mean = {margin_ifp:8.3f} mmHg")
    print(f"peak   IFP       = {float(ifp.max()):8.3f} mmHg")
    print(f"core   speed     = {core_speed:8.4f} um/s (should be ~0)")
    print(f"margin speed     = {margin_speed:8.4f} um/s")
    print(f"margin v_radial  = {margin_vrad:8.4f} um/s (>0 == OUTWARD, opposes delivery)")

    assert core_ifp > margin_ifp, (
        f"core IFP {core_ifp:.3f} should exceed margin IFP {margin_ifp:.3f}")
    assert margin_vrad > 0.0, (
        f"margin radial velocity {margin_vrad:.4f} should be outward (>0)")
    assert core_speed < margin_speed, (
        f"core flow {core_speed:.4f} should be slower than margin {margin_speed:.4f}")

    print("\nOK: elevated flat core IFP -> ~0 core flow; outward margin flow (Jain effect).")
