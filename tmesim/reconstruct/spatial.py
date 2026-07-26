"""Continuum tissue fields from a sparse cell population (`reconstruct` stage).

Given a discrete `CellPopulation` (positions + phenotypes) and the simulation
`Domain`, we build the continuum `TissueFields` the physics solvers need:

- **ECM density** — a stroma-weighted kernel-density estimate of cell positions.
  Stromal (reactive astrocyte / CAF-like) cells lay down the most matrix, tumor
  cells some, immune/endothelial little. Dense stromal/tumor regions therefore
  get high ECM (raising tortuosity -> lowering effective drug diffusivity).
  Normalized to [0, 1].
- **Vascular density** — normally supplied by `reconstruct.vasculature`. If the
  caller does not pass one, we synthesize a simple plausible field (perfused at
  the tissue rim, sparser in dense-tumor interior) so the object is usable, but
  callers SHOULD prefer the real vasculature module.

Oxygen / IFP / velocity are left `None` here — they are solved later by the
`physics` package. The estimator is fully vectorized and deterministic.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.ndimage import gaussian_filter

from ..core.domain import Domain
from ..core.microenvironment import CellPopulation, CellType, TissueFields

# Per-CellType weight for how much extracellular matrix a cell contributes.
# Stroma dominates ECM deposition; tumor contributes moderately; immune,
# endothelial and necrotic debris contribute little structured matrix.
_ECM_WEIGHTS: dict[int, float] = {
    int(CellType.TUMOR): 0.6,
    int(CellType.STROMAL): 1.0,
    int(CellType.IMMUNE): 0.1,
    int(CellType.ENDOTHELIAL): 0.1,
    int(CellType.NECROTIC): 0.2,
}

# KDE bandwidth in micrometers. Chosen so single cells smear over a biologically
# reasonable neighbourhood (~a few tens of um) rather than a single grid cell.
_KDE_BANDWIDTH_UM: float = 40.0


def _bin_cells(
    domain: Domain,
    cells: CellPopulation,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Accumulate (optionally weighted) cell counts into a (ny, nx) grid map.

    Uses `Domain.to_index` so binning matches the grid convention the physics
    solvers assume, then scatters with `np.add.at` (vectorized, handles repeated
    indices correctly).
    """
    grid = domain.zeros()
    if cells.n == 0:
        return grid
    col, row = domain.to_index(cells.x, cells.y)
    if weights is None:
        np.add.at(grid, (row, col), 1.0)
    else:
        np.add.at(grid, (row, col), np.asarray(weights, dtype=np.float64))
    return grid


def _normalize01(arr: np.ndarray) -> np.ndarray:
    """Scale an array to [0, 1]. A flat (or all-zero) array maps to all zeros."""
    lo = float(arr.min())
    hi = float(arr.max())
    if hi - lo <= 0.0:
        return np.zeros_like(arr)
    return (arr - lo) / (hi - lo)


def _synthesize_vascular_density(
    domain: Domain, cells: CellPopulation
) -> np.ndarray:
    """Fallback vascular field when none is supplied.

    Deliberately simple and deterministic: uniformly perfused baseline, made
    SPARSER where tumor cells are dense (mimicking the poorly-vascularized,
    hypoxia-prone tumor interior that motivates the whole delivery problem).
    Callers should prefer `reconstruct.vasculature.build_vascular_network`.
    """
    base = np.full(domain.shape, 0.35, dtype=np.float64)
    if cells.n > 0:
        tum = cells.tumor_mask()
        if tum.any():
            sigma = _KDE_BANDWIDTH_UM / domain.dx
            tumor_map = _bin_cells(
                domain,
                CellPopulation(
                    x=cells.x[tum], y=cells.y[tum],
                    cell_type=cells.cell_type[tum], state=cells.state[tum],
                    sensitivity=cells.sensitivity[tum],
                ),
            )
            tumor_density = _normalize01(gaussian_filter(tumor_map, sigma=sigma))
            # Suppress vasculature by up to 80% in the densest tumor core.
            base = base * (1.0 - 0.8 * tumor_density)
    # Confine to tissue and clip into [0, 1].
    base = np.where(domain.tissue_mask, base, 0.0)
    return np.clip(base, 0.0, 1.0)


def build_tissue_fields(
    domain: Domain,
    cells: CellPopulation,
    vascular_density: Optional[np.ndarray] = None,
) -> TissueFields:
    """Build continuum `TissueFields` from a discrete cell population.

    Kernel-density-estimates the (stroma-weighted) cell positions onto the grid
    via a Gaussian filter on a binned cell-count map to produce ECM density, and
    either uses the supplied vascular field or synthesizes a simple one.

    Parameters
    ----------
    domain : Domain
        The simulation grid. Fields are returned with shape ``domain.shape``.
    cells : CellPopulation
        Discrete cells with positions in micrometers. Zero cells is handled
        gracefully (ECM is all-zero).
    vascular_density : np.ndarray, optional
        A precomputed (ny, nx) vascular-density field in [0, 1] — normally from
        `reconstruct.vasculature`. If ``None``, a simple field is synthesized.

    Returns
    -------
    TissueFields
        With ``ecm_density`` and ``vascular_density`` populated (both in [0, 1]),
        and ``oxygen``/``ifp``/``velocity`` left ``None`` for the physics stage.

    Raises
    ------
    ValueError
        If a supplied ``vascular_density`` does not match ``domain.shape``.
    """
    sigma = _KDE_BANDWIDTH_UM / domain.dx

    # --- ECM density: stroma-weighted KDE of cell positions ------------------
    if cells.n == 0:
        ecm = domain.zeros()
    else:
        weights = np.array(
            [_ECM_WEIGHTS.get(int(t), 0.1) for t in cells.cell_type],
            dtype=np.float64,
        )
        weighted_map = _bin_cells(domain, cells, weights=weights)
        ecm = gaussian_filter(weighted_map, sigma=sigma)
        ecm = _normalize01(ecm)
    # Only tissue carries matrix.
    ecm = np.where(domain.tissue_mask, ecm, 0.0)

    # --- Vascular density: supplied or synthesized ---------------------------
    if vascular_density is None:
        vasc = _synthesize_vascular_density(domain, cells)
    else:
        vasc = np.asarray(vascular_density, dtype=np.float64)
        if vasc.shape != domain.shape:
            raise ValueError(
                f"vascular_density shape {vasc.shape} != domain.shape {domain.shape}"
            )
        vasc = np.clip(vasc, 0.0, 1.0)

    return TissueFields(ecm_density=ecm, vascular_density=vasc)


if __name__ == "__main__":
    # Self-test: build fields for a random population and for the empty edge case.
    rng = np.random.default_rng(0)
    dom = Domain(nx=120, ny=120, dx=15.0)
    n = 4000
    xmax, ymax = dom.extent_um
    ctype = rng.integers(0, 5, size=n).astype(np.int8)
    pop = CellPopulation(
        x=rng.uniform(0, xmax, size=n),
        y=rng.uniform(0, ymax, size=n),
        cell_type=ctype,
        state=np.zeros(n, dtype=np.int8),
        sensitivity=rng.uniform(0, 1, size=n),
    )
    fields = build_tissue_fields(dom, pop)
    print("ecm  range [%.3f, %.3f] mean %.3f shape %s"
          % (fields.ecm_density.min(), fields.ecm_density.max(),
             fields.ecm_density.mean(), fields.ecm_density.shape))
    print("vasc range [%.3f, %.3f] mean %.3f shape %s"
          % (fields.vascular_density.min(), fields.vascular_density.max(),
             fields.vascular_density.mean(), fields.vascular_density.shape))
    print("oxygen is None:", fields.oxygen is None)

    empty = CellPopulation(
        x=np.array([]), y=np.array([]), cell_type=np.array([]),
        state=np.array([]), sensitivity=np.array([]),
    )
    ef = build_tissue_fields(dom, empty)
    print("zero-cell ecm all-zero:", bool((ef.ecm_density == 0).all()))
    print("spatial.py self-test OK")
