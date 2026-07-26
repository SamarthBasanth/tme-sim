"""The central data objects: cells, tissue fields, and the reconstructed
patient-specific tumor microenvironment.

`Microenvironment` is the hinge of the whole project. It is the OUTPUT of the
`reconstruct` package (sample -> digital microenvironment) and the INPUT to the
`sim` package (microenvironment -> predicted drug response). Every module either
produces one, consumes one, or transforms fields inside one.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Dict, Optional

import numpy as np

from .domain import Domain


class CellType(IntEnum):
    """Coarse phenotype classes in the glioblastoma microenvironment."""

    TUMOR = 0        # malignant glioma cells
    STROMAL = 1      # reactive astrocytes / CAF-like stroma
    IMMUNE = 2       # tumor-associated macrophages / T cells
    ENDOTHELIAL = 3  # vessel-lining (drug source locations)
    NECROTIC = 4     # necrotic core debris


class CellState(IntEnum):
    """Dynamic per-cell state, updated by the agent model each step."""

    VIABLE = 0
    HYPOXIC = 1     # low O2, often drug-resistant / quiescent
    QUIESCENT = 2   # non-proliferating but alive
    DEAD = 3        # killed (by drug or necrosis)


@dataclass
class CellPopulation:
    """Vectorized agent-based cells. Parallel numpy arrays, all length n_cells.

    Using parallel arrays (not per-cell objects) keeps the agent update fast and
    vectorizable for tens of thousands of cells.

    Attributes
    ----------
    x, y : (n,) float — cell-center positions in micrometers.
    cell_type : (n,) int — CellType values.
    state : (n,) int — CellState values.
    sensitivity : (n,) float in [0,1] — intrinsic drug sensitivity
        (1 = fully sensitive, 0 = fully resistant). Encodes the biology a pure
        transport model cannot: MGMT status, resistance clones, etc.
    exposure : (n,) float — accumulated local drug exposure (AUC) per cell,
        filled in during simulation. Initialized to zeros.
    """

    x: np.ndarray
    y: np.ndarray
    cell_type: np.ndarray
    state: np.ndarray
    sensitivity: np.ndarray
    exposure: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=np.float64)
        self.y = np.asarray(self.y, dtype=np.float64)
        self.cell_type = np.asarray(self.cell_type, dtype=np.int8)
        self.state = np.asarray(self.state, dtype=np.int8)
        self.sensitivity = np.asarray(self.sensitivity, dtype=np.float64)
        n = len(self.x)
        for name in ("y", "cell_type", "state", "sensitivity"):
            if len(getattr(self, name)) != n:
                raise ValueError(f"CellPopulation array '{name}' length != {n}")
        if self.exposure is None:
            self.exposure = np.zeros(n, dtype=np.float64)
        else:
            self.exposure = np.asarray(self.exposure, dtype=np.float64)

    @property
    def n(self) -> int:
        return len(self.x)

    def tumor_mask(self) -> np.ndarray:
        return self.cell_type == CellType.TUMOR

    def alive_mask(self) -> np.ndarray:
        return self.state != CellState.DEAD

    def viable_tumor_fraction(self) -> float:
        """Fraction of the ORIGINAL tumor cells still alive. This is the primary
        simulation readout (1.0 = no kill, 0.0 = complete kill)."""
        tum = self.tumor_mask()
        n_tumor = int(tum.sum())
        if n_tumor == 0:
            return 0.0
        alive_tumor = int((tum & self.alive_mask()).sum())
        return alive_tumor / n_tumor

    def copy(self) -> "CellPopulation":
        return CellPopulation(
            x=self.x.copy(), y=self.y.copy(), cell_type=self.cell_type.copy(),
            state=self.state.copy(), sensitivity=self.sensitivity.copy(),
            exposure=self.exposure.copy(),
        )


@dataclass
class TissueFields:
    """Continuum physical fields on the Domain grid, each shape (ny, nx).

    ecm_density : [0,1] relative extracellular-matrix density (raises tortuosity,
        lowers effective drug diffusivity).
    vascular_density : [0,1] fraction of grid cell that is perfused vessel. Acts
        as the drug SOURCE and the oxygen source.
    oxygen : relative [0,1] (1 = well-oxygenated). Filled by physics.oxygen.
    ifp : interstitial fluid pressure (mmHg). Filled by physics.pressure.
    velocity : optional (2, ny, nx) convective velocity (um/s) derived from -grad(ifp).
    """

    ecm_density: np.ndarray
    vascular_density: np.ndarray
    oxygen: Optional[np.ndarray] = None
    ifp: Optional[np.ndarray] = None
    velocity: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        self.ecm_density = np.asarray(self.ecm_density, dtype=np.float64)
        self.vascular_density = np.asarray(self.vascular_density, dtype=np.float64)


@dataclass
class Microenvironment:
    """Patient-specific reconstructed TME — the central object of the pipeline.

    metadata may carry non-spatial patient features used by baselines and by the
    per-cell sensitivity model, e.g. {'subtype': 'mesenchymal', 'mgmt_methylated':
    True, 'expression': {...}}.
    """

    domain: Domain
    fields: TissueFields
    cells: CellPopulation
    patient_id: str = "synthetic"
    metadata: Dict = field(default_factory=dict)

    def summary(self) -> Dict:
        return {
            "patient_id": self.patient_id,
            "grid": self.domain.shape,
            "dx_um": self.domain.dx,
            "n_cells": self.cells.n,
            "n_tumor": int(self.cells.tumor_mask().sum()),
            "mean_ecm": float(self.fields.ecm_density.mean()),
            "mean_vascular": float(self.fields.vascular_density.mean()),
            "metadata_keys": list(self.metadata.keys()),
        }
