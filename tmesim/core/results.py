"""Shared result containers passed between physics, simulation, and validation.

Defining these here (rather than letting each module invent its own return type)
is what lets the simulation engine and the validation harness be written against
a fixed contract while the physics/cell modules are built in parallel.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np


@dataclass
class TransportResult:
    """Output of a drug-transport PDE solve on the microenvironment grid.

    auc : (ny, nx) accumulated exposure = integral of concentration over time
        (uM * s). This is the field the pharmacodynamic kill model consumes.
    final_concentration : (ny, nx) drug concentration at the last step (uM).
    times : (n_saved,) time stamps of saved snapshots (s).
    mean_concentration_history : (n_saved,) tissue-mean concentration over time.
    penetration_index : scalar summary in [0,1] — how far the drug reached
        relative to a well-mixed ideal (1 = perfect penetration).
    """

    auc: np.ndarray
    final_concentration: np.ndarray
    times: np.ndarray
    mean_concentration_history: np.ndarray
    penetration_index: float = 0.0


@dataclass
class SimulationResult:
    """Output of a full microenvironment simulation for one patient/drug.

    predicted_response : float in [0,1] — the model's prediction that this tumor
        RESPONDS to this drug (1 = strong response/kill, 0 = no response). This is
        the number the validation harness scores against observed outcomes.
    viable_tumor_fraction : final fraction of tumor cells still alive.
    kill_fraction : 1 - viable_tumor_fraction.
    viable_history : (n_saved,) viable tumor fraction over the treatment course.
    times : (n_saved,) time stamps (s).
    transport : the TransportResult that drove the kill (for inspection/plots).
    fields : dict of named (ny,nx) fields to visualize (oxygen, ifp, auc, ...).
    extras : free-form scalars/diagnostics.
    """

    predicted_response: float
    viable_tumor_fraction: float
    kill_fraction: float
    viable_history: np.ndarray
    times: np.ndarray
    transport: Optional[TransportResult] = None
    fields: Dict[str, np.ndarray] = field(default_factory=dict)
    extras: Dict[str, float] = field(default_factory=dict)


@dataclass
class PatientRecord:
    """One patient/sample for the Stage 0 cohort experiment.

    micro : the reconstructed Microenvironment.
    observed_response : ground-truth label (1 = responder, 0 = non-responder) OR
        a continuous sensitivity (e.g. -log(IC50)); the harness handles both.
    features : flat dict of scalar features for the non-spatial baselines.
    """

    patient_id: str
    micro: object  # tmesim.core.Microenvironment (avoid circular import)
    observed_response: float
    features: Dict[str, float] = field(default_factory=dict)


@dataclass
class ExperimentResult:
    """Output of the Stage 0 kill-test across a cohort.

    The headline number is `verdict`: did the mechanistic microenvironment model
    beat the cheap baseline at predicting observed drug response?
    """

    model_scores: Dict[str, float]        # metric_name -> value for the sim model
    baseline_scores: Dict[str, Dict[str, float]]  # baseline_name -> {metric: value}
    per_patient: List[Dict] = field(default_factory=list)
    primary_metric: str = "roc_auc"
    verdict: str = ""                     # human-readable pass/fail summary
    beat_baseline: bool = False
