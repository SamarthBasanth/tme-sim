"""Steady-state oxygen reaction-diffusion field and hypoxia labeling.

Oxygen matters for two reasons in the Stage 0 kill-test:

1. Tumor cores are poorly vascularized (see ``reconstruct.vasculature``), so O2
   diffuses in from vessels but is consumed by dense tumor tissue. The balance
   produces the classic hypoxic core.
2. Hypoxic cells are drug/therapy resistant, which the pharmacodynamic module
   (``cells.pkpd``) uses to modulate kill. This module supplies the O2 field and
   the per-cell HYPOXIC label that resistance keys off of.

We solve the *steady-state* reaction-diffusion balance

    D * lap(O) + supply * vasc * (1 - O) - consumption * tumor_density * O = 0

on the (ny, nx) grid with no-flux (Neumann) boundaries. Rearranged into a linear
system ``A @ O = b`` it is symmetric, sparse and diagonally dominant, so a single
sparse direct solve gives the exact steady state (no iteration/CFL concerns).

    diffusion:  D * lap(O)
    source:     supply * vasc          (turned off as O -> 1 via the (1-O) term)
    sink:       consumption * rho_tum  (proportional to LOCAL tumor-cell density)

Units follow the project convention: length um, diffusivity um^2/s. The oxygen
field is dimensionless relative saturation in [0, 1] (1 = fully oxygenated).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from ..core import CellState, Microenvironment, OxygenConfig


def _tumor_density(micro: Microenvironment) -> np.ndarray:
    """Bin tumor cells onto the grid, returning per-cell tumor counts (ny, nx).

    Uses ``Domain.to_index`` to map each tumor cell's (x, y) in micrometers to its
    (col, row) grid index and accumulates counts with ``np.add.at`` (vectorized).
    The raw count per grid cell is the local tumor-cell density that drives the
    oxygen sink.
    """
    domain = micro.domain
    cells = micro.cells
    density = domain.zeros()
    tum = cells.tumor_mask()
    if not np.any(tum):
        return density
    col, row = domain.to_index(cells.x[tum], cells.y[tum])
    np.add.at(density, (row, col), 1.0)
    return density


def _build_laplacian(ny: int, nx: int, dx: float) -> sp.csr_matrix:
    """5-point finite-difference Laplacian on (ny, nx) with no-flux boundaries.

    No-flux (Neumann) boundaries are imposed by ghost cells equal to the interior
    cell, i.e. a missing neighbor simply drops out of the stencil. The returned
    operator ``L`` approximates ``lap(O)`` such that ``(L @ O.ravel())`` has the
    same sign convention as ``(neighbors - k*center)/dx^2``.
    """
    n = ny * nx
    idx = np.arange(n).reshape(ny, nx)
    inv = 1.0 / (dx * dx)

    rows = []
    cols = []
    vals = []
    diag = np.zeros((ny, nx), dtype=np.float64)

    # Four cardinal neighbors; each existing neighbor contributes +inv off-diagonal
    # and -inv to the center. Missing (boundary) neighbors are dropped => no flux.
    for src, dst, dslice in (
        (idx[1:, :], idx[:-1, :], (slice(1, None), slice(None))),   # north
        (idx[:-1, :], idx[1:, :], (slice(0, -1), slice(None))),     # south
        (idx[:, 1:], idx[:, :-1], (slice(None), slice(1, None))),   # west
        (idx[:, :-1], idx[:, 1:], (slice(None), slice(0, -1))),     # east
    ):
        rows.append(src.ravel())
        cols.append(dst.ravel())
        vals.append(np.full(src.size, inv))
        diag[dslice] -= inv

    rows.append(idx.ravel())
    cols.append(idx.ravel())
    vals.append(diag.ravel())

    L = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
    ).tocsr()
    return L


def solve_oxygen(micro: Microenvironment, cfg: OxygenConfig) -> np.ndarray:
    """Solve the steady-state reaction-diffusion oxygen field.

    Solves ``D*lap(O) + supply*vasc*(1-O) - consumption*rho_tum*O = 0`` on the
    domain grid with no-flux boundaries, where ``vasc`` is
    ``micro.fields.vascular_density`` (the O2 source) and ``rho_tum`` is the local
    tumor-cell density binned onto the grid (the O2 sink). The linear system

        (-D*L + diag(supply*vasc + consumption*rho_tum)) @ O = supply*vasc

    is symmetric positive-definite and solved directly with a sparse solver.

    Parameters
    ----------
    micro : Microenvironment
        Provides the domain, ``fields.vascular_density`` and the cell population.
    cfg : OxygenConfig
        ``diffusivity`` (um^2/s), ``supply`` and ``consumption`` rate constants.

    Returns
    -------
    np.ndarray
        Oxygen saturation field, shape (ny, nx), clipped to [0, 1]. The caller is
        responsible for storing it into ``micro.fields.oxygen``.
    """
    domain = micro.domain
    ny, nx = domain.shape
    n = ny * nx

    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    if vasc.shape != (ny, nx):
        raise ValueError(
            f"vascular_density shape {vasc.shape} != domain shape {(ny, nx)}"
        )
    vasc = np.clip(vasc, 0.0, 1.0)
    rho_tum = _tumor_density(micro)

    D = float(cfg.diffusivity)
    supply = float(cfg.supply)
    consumption = float(cfg.consumption)

    L = _build_laplacian(ny, nx, domain.dx)

    source = (supply * vasc).ravel()                          # RHS b
    reaction = (supply * vasc + consumption * rho_tum).ravel()  # O-proportional loss

    # A = -D * L + diag(reaction). (-L is positive-definite; reaction >= 0.)
    A = (-D * L) + sp.diags(reaction, format="csr")

    # If there is no source and no reaction anywhere the system is singular; guard.
    if not np.any(source) and not np.any(reaction):
        return domain.zeros()

    oxygen = spsolve(A.tocsc(), source)
    oxygen = np.asarray(oxygen, dtype=np.float64).reshape(ny, nx)
    np.clip(oxygen, 0.0, 1.0, out=oxygen)
    return oxygen


def label_hypoxia(
    micro: Microenvironment, oxygen: np.ndarray, cfg: OxygenConfig
) -> None:
    """Mark non-dead cells sitting in low-oxygen grid cells as HYPOXIC (in place).

    Each cell's local oxygen is sampled at its (col, row) grid index. Cells that
    are not DEAD and sit where ``oxygen < cfg.hypoxia_threshold`` have their state
    set to ``CellState.HYPOXIC``. Modifies ``micro.cells.state`` in place; returns
    nothing.

    Parameters
    ----------
    micro : Microenvironment
        Population whose ``cells.state`` is updated.
    oxygen : np.ndarray
        Oxygen field (ny, nx), e.g. from :func:`solve_oxygen`.
    cfg : OxygenConfig
        Uses ``hypoxia_threshold``.
    """
    domain = micro.domain
    cells = micro.cells
    oxygen = np.asarray(oxygen, dtype=np.float64)
    if oxygen.shape != domain.shape:
        raise ValueError(
            f"oxygen shape {oxygen.shape} != domain shape {domain.shape}"
        )

    col, row = domain.to_index(cells.x, cells.y)
    local_o2 = oxygen[row, col]

    alive = cells.state != int(CellState.DEAD)
    hypoxic = alive & (local_o2 < float(cfg.hypoxia_threshold))
    cells.state[hypoxic] = np.int8(CellState.HYPOXIC)


if __name__ == "__main__":
    from ..core import (
        CellPopulation,
        CellType,
        Domain,
        TissueFields,
    )

    # ------------------------------------------------------------------
    # Build a synthetic Domain: vasculature is sparse in a central tumor
    # region (well-perfused rim, avascular core) with tumor cells packed
    # into that central region. Expect a hypoxic core.
    # ------------------------------------------------------------------
    rng = np.random.default_rng(0)
    ny = nx = 80
    dx = 15.0
    domain = Domain(nx=nx, ny=ny, dx=dx)

    X, Y = domain.coords_um()
    cx, cy = X.mean(), Y.mean()
    r = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    core_radius = 0.30 * min(domain.extent_um)

    # Background vasculature ~0.6, suppressed toward the core (avascular center).
    vasc = 0.6 - 0.59 * np.exp(-(r ** 2) / (2 * (0.5 * core_radius) ** 2))
    vasc = np.clip(vasc, 0.005, 1.0)

    ecm = domain.ones() * 0.4
    fields = TissueFields(ecm_density=ecm, vascular_density=vasc)

    # Tumor cells packed densely into the central core (strong O2 sink); a few
    # scattered stromal cells elsewhere.
    n_tumor = 12000
    ang = rng.uniform(0, 2 * np.pi, n_tumor)
    rad = core_radius * np.sqrt(rng.uniform(0, 1, n_tumor))
    tx = np.clip(cx + rad * np.cos(ang), 0, domain.extent_um[0] - 1e-3)
    ty = np.clip(cy + rad * np.sin(ang), 0, domain.extent_um[1] - 1e-3)

    n_stroma = 300
    sx = rng.uniform(0, domain.extent_um[0], n_stroma)
    sy = rng.uniform(0, domain.extent_um[1], n_stroma)

    x = np.concatenate([tx, sx])
    y = np.concatenate([ty, sy])
    ctype = np.concatenate([
        np.full(n_tumor, int(CellType.TUMOR)),
        np.full(n_stroma, int(CellType.STROMAL)),
    ])
    state = np.zeros(x.size, dtype=np.int8)
    sensitivity = rng.uniform(0, 1, x.size)
    cells = CellPopulation(x=x, y=y, cell_type=ctype, state=state,
                           sensitivity=sensitivity)

    micro = Microenvironment(domain=domain, fields=fields, cells=cells,
                             patient_id="oxygen-selftest")

    cfg = OxygenConfig()
    oxygen = solve_oxygen(micro, cfg)

    assert oxygen.shape == domain.shape
    assert oxygen.min() >= 0.0 and oxygen.max() <= 1.0, "oxygen out of [0,1]"

    # Determinism: solving again gives an identical field.
    oxygen2 = solve_oxygen(micro, cfg)
    assert np.array_equal(oxygen, oxygen2), "solve_oxygen is not deterministic"

    # Core vs periphery: mean O2 in the central core disc vs an outer annulus.
    core_mask = r < 0.5 * core_radius
    periph_mask = r > 1.3 * core_radius
    core_o2 = float(oxygen[core_mask].mean())
    periph_o2 = float(oxygen[periph_mask].mean())

    print(f"grid = {domain.shape}, dx = {dx} um")
    print(f"oxygen range = [{oxygen.min():.4f}, {oxygen.max():.4f}]")
    print(f"tumor-core mean O2  = {core_o2:.4f}")
    print(f"periphery mean O2   = {periph_o2:.4f}")
    assert core_o2 < periph_o2, (
        f"tumor core ({core_o2:.4f}) is not more hypoxic than periphery "
        f"({periph_o2:.4f})"
    )

    # Hypoxia labeling marks some (core) cells.
    label_hypoxia(micro, oxygen, cfg)
    n_hypoxic = int((micro.cells.state == int(CellState.HYPOXIC)).sum())
    print(f"hypoxia_threshold   = {cfg.hypoxia_threshold}")
    print(f"n hypoxic cells     = {n_hypoxic} / {micro.cells.n}")
    assert n_hypoxic > 0, "label_hypoxia marked no cells"

    print("\nself-test PASSED")
