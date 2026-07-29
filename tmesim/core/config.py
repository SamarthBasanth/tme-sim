"""Configuration dataclasses and the drug library schema.

All physical parameters carry explicit units in comments. Defaults are chosen to
be biologically plausible for a small-molecule chemotherapeutic in brain tumor
tissue (order-of-magnitude, not patient-calibrated — calibration is exactly the
Stage-1+ work this scaffold is built to enable).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


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
    hypoxia_resistance : optional PER-DRUG override of
        ``PDConfig.hypoxia_resistance`` (see that field's docstring for the
        general unvalidated-prior caveat). ``None`` (the default for every drug
        except temozolomide, see :data:`DRUG_LIBRARY`) means "use the global
        ``PDConfig`` default" — i.e. this field changes nothing for a drug that
        does not set it. It exists because the literature basis found in
        ``results/sweep/state_multiplier_literature_check_*.md`` is
        drug-specific (real GBM-cell-line hypoxia-vs-normoxia TMZ data), not a
        generic chemotherapy constant, and ``PDConfig`` has no drug axis of its
        own. ``tmesim.cells.pkpd.apply_drug_kill`` (and every module that
        recomputes the same resistance vector directly, e.g.
        ``validation.kill_calibration.resistance_vector``) resolves the
        effective multiplier as ``drug.hypoxia_resistance if drug.hypoxia_resistance
        is not None else pd.hypoxia_resistance``.
    """

    name: str
    diffusivity: float = 250.0
    molecular_weight: float = 300.0
    ic50: float = 5.0
    hill: float = 1.5
    binding: float = 5e-4
    clearance: float = 2e-4
    bbb_penetrance: float = 0.3
    hypoxia_resistance: Optional[float] = None


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
    #: Drive the transport source with the full time-resolved plasma curve
    #: ``plasma(t)`` instead of its window time-average. ``False`` (default)
    #: preserves the historical Stage-0 plasma-averaging behaviour BIT-FOR-BIT;
    #: ``True`` lets tissue concentration genuinely rise and fall, which is what
    #: makes dosing SCHEDULE (fractionation, infusion duration) resolvable at
    #: all. See :func:`tmesim.sim.engine.simulate_response` and the
    #: ``source_profile`` argument of
    #: :func:`tmesim.physics.diffusion.simulate_drug_transport`.
    time_resolved: bool = False


@dataclass
class OxygenConfig:
    """Reaction-diffusion parameters for the steady-state oxygen field.

    ``consumption`` is UNCALIBRATED against any measured pO2 (see
    ``tmesim/physics/oxygen.py`` module docstring and
    ``results/oxygen_sensitivity_band.md`` for the quantified, literature-
    bracketed uncertainty range this default sits inside, rather than at the
    center of). ``sink_form`` selects between two functional forms for how the
    sink depends on local oxygen -- see :func:`tmesim.physics.oxygen.solve_oxygen`.
    Both ``consumption`` and ``sink_form`` default to the ORIGINAL behavior so
    nothing existing changes silently; the alternative sink form and the
    consumption sensitivity band are opt-in / reported-not-adopted.
    """

    diffusivity: float = 2000.0  # um^2/s
    consumption: float = 0.02    # per-tumor-cell consumption rate (relative)
    supply: float = 0.05         # vascular supply rate (relative)
    hypoxia_threshold: float = 0.2  # oxygen below this => HYPOXIC state

    #: {"linear", "michaelis_menten"}. "linear" (default) is the ORIGINAL sink,
    #: ``consumption * rho_tum * O`` -- proportional to O everywhere, so it
    #: self-throttles as O falls and never dominates the (bounded) supply term.
    #: "michaelis_menten" is ``consumption * rho_tum * O / (O + mm_half_saturation)``
    #: (Grimes, Fletcher & Partridge 2014, PMID 26064525) -- saturates toward a
    #: CONSTANT rate as O falls, so it does not self-throttle. See
    #: ``tmesim.physics.oxygen.solve_oxygen`` and ``results/oxygen_sensitivity_band.md``
    #: for why this axis is NOT a minor correction on this project's real flagship
    #: reconstruction (hypoxic tumor fraction swings 0.0 -> 0.69-1.0 across the
    #: literature-plausible half-saturation range at the current ``consumption``).
    sink_form: str = "linear"

    #: Dimensionless Michaelis-Menten half-saturation constant K (only used when
    #: ``sink_form == "michaelis_menten"``). Default = Grimes' own literature
    #: base case, ``Km = 1 mmHg`` (their sensitivity ceiling is 10 mmHg), expressed
    #: against this project's own O=1 anchor (``oxygen.REFERENCE_PO2_MMHG = 35``
    #: mmHg, Erecinska & Silver 2001, PMID 11718758): ``K = 1 / 35 = 0.0286``.
    #: This is A defensible choice, not THE resolved value -- see
    #: ``results/oxygen_sensitivity_band.md``.
    mm_half_saturation: float = 1.0 / 35.0

    #: Newton-Raphson iteration cap for the Michaelis-Menten sink (unused by
    #: "linear", which is a single direct sparse solve). Typical convergence is
    #: a handful of iterations; see ``tmesim/physics/oxygen.py``.
    mm_max_iter: int = 50

    #: Newton-Raphson convergence tolerance (max absolute change in O between
    #: iterations) for the Michaelis-Menten sink.
    mm_tol: float = 1e-9


@dataclass
class PDConfig:
    """Pharmacodynamics: how accumulated exposure maps to cell kill.

    ``hypoxia_resistance`` (this GLOBAL default) and ``quiescent_resistance`` are
    still UNVALIDATED PRIORS with no citation (flagged in-repo at
    ``validation/kill_calibration.py:1954,1997``). A literature search
    (``results/sweep/state_multiplier_literature_check_*.md``) found real,
    GBM-cell-line, hypoxia-vs-normoxia temozolomide data (Musah-Eroje & Watson
    2019, PMID 30694423; Wang et al. 2017, PMID 28801626; Ahmed et al. 2018, PMID
    29492900) that converges on a much SMALLER hypoxia effect for TMZ
    specifically than 0.35 implies -- IC50 shifts of only ~1.1-1.4x, and one clean
    null result in the matched-2D-oxygen comparison closest to what this
    multiplier represents.

    **TMZ-specific recalibration (done).** That IC50-fold evidence has been
    converted, Hill-curve-consistently, into an equivalent kill-probability
    multiplier at this project's own documented operating point (temozolomide,
    ``results/dose_schedule.md``: median tumour-cell exposure/IC50 = 0.056 at 1x
    nominal dose, drug Hill coefficient h=1.2 -- the "not saturated" sub-clinical
    regime this project's real runs actually occupy). For a hypoxia-only IC50
    shift by factor f at fixed absolute exposure E (ratio r = E/IC50):

        R_equivalent(r, f, h) = (1 + r^h) / (f^h + r^h)

    which as r -> 0 (this project's regime) tends to the closed form
    ``R_equivalent ~ 1 / f^h``. Using a literature-representative f = 1.2 (within
    the reported 1.1-1.4x TMZ range; Wang et al. 2017) gives R_equivalent ~
    0.81 at r=0.056, essentially unchanged (0.806-0.811) across this project's
    whole documented sub-saturating dose band (r = 0.034-0.078, the 0.6x-1.4x
    "approximate clinical band" in ``dose_schedule.md``). The chosen value,
    **0.80**, is stored as a PER-DRUG override on temozolomide's ``DrugParams``
    (see ``DrugParams.hypoxia_resistance`` and ``DRUG_LIBRARY``), NOT here --
    this global default stays at the original, still-uncited 0.35 for cisplatin
    and doxorubicin, for which no comparably quantitative literature basis was
    found in that search (cisplatin: real, GBM-specific, statistically
    significant hypoxia effect confirmed in direction, magnitude unextracted;
    doxorubicin: a real 4-10x fold-change exists but is prostate tissue, not
    glioma, and a mismatched exposure protocol).

    For ``quiescent_resistance``, no quantitative fold-change was found anywhere
    in the literature for any of this project's three drugs -- only mechanistic
    support for the direction (<1). It has NOT been changed: inventing a
    precise number with no measurement behind it would be worse than leaving
    the honestly-labeled unvalidated prior in place.
    """

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


#: TMZ's literature-informed hypoxia-resistance override -- see
#: ``PDConfig.hypoxia_resistance`` docstring for the full operating-point
#: derivation (R_equivalent = (1+r^h)/(f^h+r^h), f~1.2 from Wang et al. 2017
#: PMID 28801626, r=0.056 from ``results/dose_schedule.md``). Applied both to
#: the synthetic "temozolomide-like" library entry below AND to the real,
#: GDSC-calibrated molecular-compile pathway
#: (``validation.calibration._compile``), which is what every real-data report
#: (``kill_calibration.md``, ``real_drug_panel.md``, ``dose_schedule.md``,
#: ``cell_state_impact.md``) actually uses for temozolomide.
TMZ_HYPOXIA_RESISTANCE: float = 0.80

# A small default glioblastoma-relevant drug library. Agents may extend this.
DRUG_LIBRARY: Dict[str, DrugParams] = {
    "temozolomide-like": DrugParams(
        name="temozolomide-like", diffusivity=300.0, molecular_weight=194.0,
        ic50=25.0, hill=1.2, bbb_penetrance=0.4, clearance=3e-4,
        hypoxia_resistance=TMZ_HYPOXIA_RESISTANCE),
    "doxorubicin-like": DrugParams(
        name="doxorubicin-like", diffusivity=150.0, molecular_weight=543.0,
        ic50=0.5, hill=1.5, bbb_penetrance=0.05, binding=2e-3, clearance=1e-4),
    "cisplatin-like": DrugParams(
        name="cisplatin-like", diffusivity=350.0, molecular_weight=300.0,
        ic50=8.0, hill=1.3, bbb_penetrance=0.15, clearance=2e-4),
}
