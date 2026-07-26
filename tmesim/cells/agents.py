"""Agent-based cell-state relabeling from the local oxygen microenvironment.

This module implements the *cheap, per-step* agent update: given a steady-state
oxygen field on the domain grid, relabel each living cell's :class:`CellState`
according to how much oxygen it sits in. Hypoxia is the biologically load-bearing
signal here — hypoxic glioma cells are quiescent and notoriously drug-resistant,
which is exactly the effect the downstream pharmacodynamic kill model
(:mod:`tmesim.cells.pkpd`) leans on.

Oxygen -> state banding
-----------------------
Let ``t = cfg.hypoxia_threshold`` (default 0.2) and let ``o`` be the local,
relative oxygen (in ``[0, 1]``; 1 == well-oxygenated) sampled at the cell's grid
index. For every cell whose current state is NOT ``DEAD``::

    o <  t          -> HYPOXIC     (severe O2 debt; resistant/quiescent biology)
    t <= o < 2*t    -> QUIESCENT   (marginal O2; alive but not proliferating)
    o >= 2*t        -> VIABLE      (well perfused)

The intermediate ``[t, 2t)`` QUIESCENT band is a deliberately simple, documented
choice (an "at-risk" perfusion shoulder just above the hypoxia line). ``DEAD``
cells are never touched — the update cannot resurrect a killed cell.

The relabel is fully vectorized (one boolean-masked fancy-index into the oxygen
field) and deterministic — no randomness is involved.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from ..core.config import OxygenConfig
from ..core.microenvironment import CellState, CellType, Microenvironment


def update_cell_states(
    micro: Microenvironment, oxygen: np.ndarray, cfg: OxygenConfig
) -> None:
    """Relabel VIABLE / HYPOXIC / QUIESCENT from local oxygen, in place.

    For each non-``DEAD`` cell, the local oxygen value is sampled at the cell's
    grid index (via :meth:`Domain.to_index`) and the cell's ``state`` is
    overwritten according to the banding documented in the module docstring::

        o <  t        -> HYPOXIC
        t <= o < 2*t  -> QUIESCENT
        o >= 2*t      -> VIABLE

    where ``t = cfg.hypoxia_threshold``. ``DEAD`` cells are left untouched (no
    resurrection). Operates IN PLACE on ``micro.cells.state``; returns ``None``.

    Parameters
    ----------
    micro : Microenvironment
        The microenvironment whose ``cells`` are relabeled in place.
    oxygen : np.ndarray, shape (ny, nx)
        Relative oxygen field in ``[0, 1]`` on ``micro.domain``'s grid.
    cfg : OxygenConfig
        Supplies ``hypoxia_threshold``; the QUIESCENT band top is ``2 *``
        that threshold.
    """
    cells = micro.cells
    domain = micro.domain

    oxygen = np.asarray(oxygen, dtype=np.float64)
    if oxygen.shape != domain.shape:
        raise ValueError(
            f"oxygen shape {oxygen.shape} != domain shape {domain.shape}"
        )

    # Sample local oxygen at each cell's grid index. to_index returns (col, row);
    # fields are indexed [row, col] == [y, x].
    col, row = domain.to_index(cells.x, cells.y)
    local_o2 = oxygen[row, col]

    t = cfg.hypoxia_threshold

    # Only relabel living cells; DEAD cells are frozen.
    alive = cells.state != int(CellState.DEAD)

    hypoxic = alive & (local_o2 < t)
    quiescent = alive & (local_o2 >= t) & (local_o2 < 2.0 * t)
    viable = alive & (local_o2 >= 2.0 * t)

    cells.state[hypoxic] = int(CellState.HYPOXIC)
    cells.state[quiescent] = int(CellState.QUIESCENT)
    cells.state[viable] = int(CellState.VIABLE)


def state_fractions(micro: Microenvironment) -> Dict[str, float]:
    """Fractions of tumor cells in each :class:`CellState` (diagnostics).

    Considers only ``CellType.TUMOR`` cells (the simulation's readout population)
    and returns the fraction of them in each state. Keys are the lowercase state
    names ``{"viable", "hypoxic", "quiescent", "dead"}`` and the values sum to
    ``1.0`` (up to floating point) when at least one tumor cell exists; if there
    are no tumor cells, all fractions are ``0.0``.

    Parameters
    ----------
    micro : Microenvironment
        Source of the cell population.

    Returns
    -------
    dict[str, float]
        Per-state fraction of the tumor cells.
    """
    cells = micro.cells
    tumor = cells.tumor_mask()
    n_tumor = int(tumor.sum())

    names = {
        CellState.VIABLE: "viable",
        CellState.HYPOXIC: "hypoxic",
        CellState.QUIESCENT: "quiescent",
        CellState.DEAD: "dead",
    }
    if n_tumor == 0:
        return {name: 0.0 for name in names.values()}

    tumor_states = cells.state[tumor]
    return {
        name: float((tumor_states == int(state)).sum()) / n_tumor
        for state, name in names.items()
    }


if __name__ == "__main__":
    # Self-test: build a Domain + cells + a horizontal gradient oxygen field and
    # assert that cells in the low-O2 region become HYPOXIC while high-O2 cells
    # stay VIABLE, and that the QUIESCENT band lands in the middle.
    from ..core.domain import Domain
    from ..core.microenvironment import CellPopulation, TissueFields

    nx, ny, dx = 40, 10, 10.0
    domain = Domain(nx=nx, ny=ny, dx=dx)

    # Oxygen increases linearly with x from ~0 (left) to ~1 (right).
    X, _ = domain.coords_um()
    oxygen = X / (nx * dx)  # in [0, 1)
    fields = TissueFields(
        ecm_density=domain.zeros(),
        vascular_density=domain.zeros(),
        oxygen=oxygen,
    )

    cfg = OxygenConfig()  # hypoxia_threshold = 0.2 -> QUIESCENT band [0.2, 0.4)
    t = cfg.hypoxia_threshold

    # Place tumor cells along a single row, spanning the whole x-gradient.
    xs = (np.arange(nx) + 0.5) * dx
    ys = np.full(nx, (ny // 2 + 0.5) * dx)
    n = nx
    cells = CellPopulation(
        x=xs,
        y=ys,
        cell_type=np.full(n, int(CellType.TUMOR), dtype=np.int8),
        state=np.full(n, int(CellState.VIABLE), dtype=np.int8),
        sensitivity=np.full(n, 0.5),
    )
    # Mark one cell in the low-O2 region DEAD to prove it is not resurrected.
    dead_idx = 0
    cells.state[dead_idx] = int(CellState.DEAD)

    micro = Microenvironment(domain=domain, fields=fields, cells=cells)

    update_cell_states(micro, oxygen, cfg)

    # Recompute the oxygen each cell was sampled at, to check bands exactly.
    col, row = domain.to_index(cells.x, cells.y)
    local_o2 = oxygen[row, col]

    st = cells.state
    for i in range(n):
        if i == dead_idx:
            assert st[i] == int(CellState.DEAD), "DEAD cell was resurrected!"
            continue
        o = local_o2[i]
        if o < t:
            assert st[i] == int(CellState.HYPOXIC), (i, o, st[i])
        elif o < 2.0 * t:
            assert st[i] == int(CellState.QUIESCENT), (i, o, st[i])
        else:
            assert st[i] == int(CellState.VIABLE), (i, o, st[i])

    # High-level assertions the task asks for explicitly.
    low_o2_cells = local_o2 < t
    # exclude the DEAD cell from the "became hypoxic" claim
    became_hypoxic = (st == int(CellState.HYPOXIC))
    assert np.all(became_hypoxic[low_o2_cells & (np.arange(n) != dead_idx)]), \
        "some low-O2 living cell did not become HYPOXIC"

    high_o2_cells = local_o2 >= 2.0 * t
    assert np.all(st[high_o2_cells] == int(CellState.VIABLE)), \
        "some high-O2 cell is not VIABLE"

    fracs = state_fractions(micro)
    assert abs(sum(fracs.values()) - 1.0) < 1e-9, fracs

    print("self-test PASSED")
    print(f"  hypoxia_threshold t = {t}; QUIESCENT band = [{t}, {2*t})")
    print(f"  n tumor cells      = {n}  (1 pre-marked DEAD)")
    print(f"  local_o2 range     = [{local_o2.min():.3f}, {local_o2.max():.3f}]")
    print(f"  state_fractions    = "
          + ", ".join(f"{k}={v:.3f}" for k, v in fracs.items()))
