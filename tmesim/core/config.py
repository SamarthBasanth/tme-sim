"""Configuration dataclasses and the drug library schema.

All physical parameters carry explicit units in comments. Defaults are chosen to
be biologically plausible for a small-molecule chemotherapeutic in brain tumor
tissue (order-of-magnitude, not patient-calibrated — calibration is exactly the
Stage-1+ work this scaffold is built to enable).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict


@dataclass
class DrugParams:
    """Transport + pharmacodynamic parameters for one agent.

    diffusivity : effective free diffusion coefficient in tissue, um^2/s.
        (Small molecules ~ 100-500 um^2/s; antibodies ~ 1-10 um^2/s.)
    molecular_weight : Da (informational / for scaling heuristics).
    ic50 : concentration giving half-maximal kill effect, uM.
    hill : Hill coefficient for the dose-response curve.
    binding : first-order tissue binding/uptake rate, 1/s (retention).
    clearance : first-order clearance/metabolism rate in tissue, 1/s.
    bbb_penetrance : [0,1] fraction of vascular concentration that crosses into
        tissue (blood-brain-barrier factor — central to glioblastoma).
    """

    name: str
    diffusivity: float = 250.0
    molecular_weight: float = 300.0
    ic50: float = 5.0
    hill: float = 1.5
    binding: float = 5e-4
    clearance: float = 2e-4
    bbb_penetrance: float = 0.3


@dataclass
class PKConfig:
    """Systemic pharmacokinetics setting the vascular (blood) concentration that
    drives tissue delivery. A simple 1-compartment exponential by default.

    plasma_peak : peak plasma concentration after a dose, uM.
    half_life : plasma elimination half-life, seconds.
    n_doses : number of doses in the treatment course.
    dose_interval : seconds between doses.
    """

    plasma_peak: float = 50.0
    half_life: float = 3600.0 * 4
    n_doses: int = 1
    dose_interval: float = 3600.0 * 24


@dataclass
class PDEConfig:
    """Numerical settings for the drug-transport PDE solve."""

    dt: float = 0.5              # time step, seconds
    total_time: float = 3600.0 * 12  # simulated exposure window, seconds
    convection: bool = True      # include interstitial-flow advection
    save_every: int = 200        # store history every N steps


@dataclass
class OxygenConfig:
    """Reaction-diffusion parameters for the steady-state oxygen field."""

    diffusivity: float = 2000.0  # um^2/s
    consumption: float = 0.02    # per-tumor-cell consumption rate (relative)
    supply: float = 0.05         # vascular supply rate (relative)
    hypoxia_threshold: float = 0.2  # oxygen below this => HYPOXIC state


@dataclass
class PDConfig:
    """Pharmacodynamics: how accumulated exposure maps to cell kill."""

    hypoxia_resistance: float = 0.35  # multiplier on kill for hypoxic cells (<1)
    quiescent_resistance: float = 0.5
    max_kill_rate: float = 1.0        # scaling on the Hill kill term


@dataclass
class SimConfig:
    """Top-level knobs bundling the per-stage configs."""

    drug: DrugParams = field(default_factory=lambda: DrugParams(name="temozolomide-like"))
    pk: PKConfig = field(default_factory=PKConfig)
    pde: PDEConfig = field(default_factory=PDEConfig)
    oxygen: OxygenConfig = field(default_factory=OxygenConfig)
    pd: PDConfig = field(default_factory=PDConfig)
    seed: int = 0


# A small default glioblastoma-relevant drug library. Agents may extend this.
DRUG_LIBRARY: Dict[str, DrugParams] = {
    "temozolomide-like": DrugParams(
        name="temozolomide-like", diffusivity=300.0, molecular_weight=194.0,
        ic50=25.0, hill=1.2, bbb_penetrance=0.4, clearance=3e-4),
    "doxorubicin-like": DrugParams(
        name="doxorubicin-like", diffusivity=150.0, molecular_weight=543.0,
        ic50=0.5, hill=1.5, bbb_penetrance=0.05, binding=2e-3, clearance=1e-4),
    "cisplatin-like": DrugParams(
        name="cisplatin-like", diffusivity=350.0, molecular_weight=300.0,
        ic50=8.0, hill=1.3, bbb_penetrance=0.15, clearance=2e-4),
}
