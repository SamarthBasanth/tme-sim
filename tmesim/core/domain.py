"""Spatial domain for the tumor microenvironment.

Everything in the pipeline lives on a 2D finite-difference grid. We work in 2D
on purpose for the Stage 0 kill-test: reconstructing a faithful 3D
microenvironment from sparse patient samples is an unsolved inference problem,
so we deliberately dodge it here and prove the machinery in 2D (a whole-slide
section or an organoid cross-section is naturally 2D anyway).

Coordinates are in micrometers (um). Grid arrays are indexed [row=y, col=x] with
shape (ny, nx), matching image / numpy convention.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass
class Domain:
    """A 2D finite-difference grid describing the physical extent of a sample.

    Attributes
    ----------
    nx, ny : int
        Number of grid columns (x) and rows (y).
    dx : float
        Grid spacing in micrometers (isotropic; dy == dx).
    tissue_mask : np.ndarray of bool, shape (ny, nx)
        True where tissue is present (inside the simulated region). Defaults to
        all-True. Used to apply no-flux boundaries at tissue edges.
    """

    nx: int
    ny: int
    dx: float = 10.0
    tissue_mask: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        if self.tissue_mask is None:
            self.tissue_mask = np.ones((self.ny, self.nx), dtype=bool)
        else:
            self.tissue_mask = np.asarray(self.tissue_mask, dtype=bool)
            if self.tissue_mask.shape != (self.ny, self.nx):
                raise ValueError(
                    f"tissue_mask shape {self.tissue_mask.shape} != ({self.ny},{self.nx})"
                )

    @property
    def shape(self) -> Tuple[int, int]:
        """(ny, nx) — numpy array shape for fields on this domain."""
        return (self.ny, self.nx)

    @property
    def extent_um(self) -> Tuple[float, float]:
        """(width_x_um, height_y_um) physical extent."""
        return (self.nx * self.dx, self.ny * self.dx)

    def zeros(self) -> np.ndarray:
        """A fresh (ny, nx) float64 array of zeros on this grid."""
        return np.zeros((self.ny, self.nx), dtype=np.float64)

    def ones(self) -> np.ndarray:
        return np.ones((self.ny, self.nx), dtype=np.float64)

    def coords_um(self) -> Tuple[np.ndarray, np.ndarray]:
        """Meshgrid (X, Y) of cell-center coordinates in micrometers.

        Returns arrays of shape (ny, nx).
        """
        xs = (np.arange(self.nx) + 0.5) * self.dx
        ys = (np.arange(self.ny) + 0.5) * self.dx
        return np.meshgrid(xs, ys)

    def to_index(self, x_um: np.ndarray, y_um: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Map physical (x, y) in um to integer grid indices (col, row), clamped."""
        col = np.clip((np.asarray(x_um) / self.dx).astype(int), 0, self.nx - 1)
        row = np.clip((np.asarray(y_um) / self.dx).astype(int), 0, self.ny - 1)
        return col, row
