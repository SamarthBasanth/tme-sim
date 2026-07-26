"""Core data model for tmesim — shared, dependency-free objects.

Import from here everywhere else, e.g.::

    from tmesim.core import Domain, Microenvironment, CellType, DrugParams
"""
from .domain import Domain
from .microenvironment import (
    CellPopulation,
    CellState,
    CellType,
    Microenvironment,
    TissueFields,
)
from .config import (
    DRUG_LIBRARY,
    DrugParams,
    OxygenConfig,
    PDConfig,
    PDEConfig,
    PKConfig,
    SimConfig,
)
from .results import (
    ExperimentResult,
    PatientRecord,
    SimulationResult,
    TransportResult,
)

__all__ = [
    "Domain",
    "CellPopulation", "CellState", "CellType", "Microenvironment", "TissueFields",
    "DRUG_LIBRARY", "DrugParams", "OxygenConfig", "PDConfig", "PDEConfig",
    "PKConfig", "SimConfig",
    "ExperimentResult", "PatientRecord", "SimulationResult", "TransportResult",
]
