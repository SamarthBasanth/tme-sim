"""dose -> PK: the conversion most likely to be quietly wrong, made explicit.

A TPL protocol says ``dose 150 mg/m2 / every 24 h / repeat 5 / then rest 23 d /
cycles 6``. The engine consumes
``core.PKConfig(plasma_peak_uM, half_life_s, n_doses, dose_interval_s)``. Getting
from one to the other requires four quantities that are NOT in the script —
body-surface area, body weight, bioavailability and volume of distribution —
plus the drug's own plasma half-life. Every one of them is written into the
returned provenance with a source and a confidence, and anything without
evidence is labelled ``assumed``. Nothing here is silent.

The arithmetic
--------------
::

    amount_mg = dose x BSA            (mg/m2)      BSA = 1.73 m2
              = dose x weight         (mg/kg)      weight = 70 kg
              = dose                  (mg)

    plasma_peak_uM = 1000 x amount_mg x F / (Vd_L x MW_g_per_mol)

The 1000 converts mmol/L to umol/L: ``mg / (L x (mg/mmol)) = mmol/L``, and
``MW`` in g/mol is numerically ``mg/mmol``.

This is the standard one-compartment bolus relation ``C0 = F x D / Vd``. It is an
*initial* concentration, which is exactly what
:func:`tmesim.treatment.drugs.plasma_profile` wants: that function administers
each dose instantaneously to ``plasma_peak`` and decays it at ``half_life``.

`over T` is an infusion, and is NOT a bolus
-------------------------------------------
Two mechanisms, chosen by whether the resulting administration train stays
uniform (which is all ``PKConfig`` can express):

* **Sub-dose expansion (exact, preferred).** The infusion is discretised into
  ``k`` equal micro-boluses of ``amount/k`` spaced ``T/k`` apart. Drug genuinely
  ARRIVES over ``T`` instead of at ``t = 0``: the plasma curve rises to a plateau
  and then decays. Used whenever the full administration train — every dose of
  every cycle inside the observation window, each expanded into its ``k``
  sub-doses — is uniformly spaced, which covers the single infusion and the
  continuous infusion.
* **Peak attenuation (approximate, and labelled).** When the train is not
  uniform (an infusion repeated on a longer cycle), ``PKConfig`` cannot represent
  it. The peak is then scaled by the exact one-compartment infusion factor
  ``(1 - exp(-lambda T)) / (lambda T)``, which reproduces the true Cmax of a
  constant-rate infusion and preserves the AUC exactly, but starts at that Cmax
  rather than rising to it. The approximation, its direction and its size are
  reported.

Either way the executor sets ``PDEConfig.time_resolved = True``, because with
the default (``False``) the solver replaces the whole plasma curve by its window
mean and every schedule becomes the same input — the artefact
``results/dose_schedule.md`` was regenerated to remove.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .units import Quantity, _fmt

__all__ = [
    "BSA_M2", "BODY_WEIGHT_KG", "PK_LITERATURE", "DEFAULT_MAX_WINDOW_S",
    "Assumption", "DoseConversion", "schedule_times", "convert_dose",
    "WindowPlan", "window_plan",
]

# --------------------------------------------------------------------------- #
# The two patient constants. Named in PROTOCOL_LANGUAGE.md; fixed here so a
# script never silently depends on a body size it did not state.
# --------------------------------------------------------------------------- #
#: Body-surface area used to turn mg/m2 into milligrams. 1.73 m2 is the standard
#: adult normalisation (the same one renal-function reporting uses) and is the
#: value the spec names.
BSA_M2 = 1.73

#: Body weight used to turn mg/kg into milligrams.
BODY_WEIGHT_KG = 70.0

#: Generic volume of distribution when the drug has no curated entry: 0.6 L/kg,
#: i.e. total body water. An ASSUMPTION, flagged as one everywhere it is used.
GENERIC_VD_L_PER_KG = 0.6

#: How long a transport solve may cover. A Stupp course is ~168 days; a
#: PDE solve at dt = 10 s over that is 1.45e6 steps and is not a thing this
#: backend will start. The window is truncated to this and the truncation is
#: reported (with how many doses fell outside it) rather than hidden.
#: 24 h covers a full daily dosing interval, so a `every 24 h` schedule is
#: resolved for one complete cycle of its own period. Overridable with
#: ``TMESIM_PROTOCOL_MAX_WINDOW_S`` (seconds) — the tests use a short window so
#: the suite stays fast, and the value in force is reported in every result.
DEFAULT_MAX_WINDOW_S = float(
    os.environ.get("TMESIM_PROTOCOL_MAX_WINDOW_S", 24.0 * 3600.0))

#: Target spacing of the micro-boluses an infusion is discretised into, and the
#: clamp on their count. 15 min is well below the shortest half-life in the
#: curated table (carmustine, ~22 min), so the discretisation never aliases the
#: elimination it is riding on.
_SUBDOSE_TARGET_S = 900.0
_SUBDOSE_MIN = 2
_SUBDOSE_MAX = 64


@dataclass
class Assumption:
    """One quantity the conversion needed and the script did not supply."""

    quantity: str
    value: float
    unit: str
    provenance: str          # measured | literature | computed | assumed | flagged-default
    source: str

    def to_dict(self) -> Dict:
        return {"quantity": self.quantity, "value": self.value,
                "unit": self.unit, "provenance": self.provenance,
                "source": self.source}


# --------------------------------------------------------------------------- #
# Curated per-drug plasma PK. Only drugs whose numbers can be cited are here;
# everything else falls to the generic defaults and is labelled `assumed`.
# `confidence` is the provenance tag that travels with each number.
# --------------------------------------------------------------------------- #
PK_LITERATURE: Dict[str, Dict] = {
    "temozolomide": {
        "half_life_s": 1.8 * 3600.0,
        "half_life_source": (
            "plasma elimination half-life ~1.8 h in malignant-glioma patients "
            "(Ostermann et al., Clin Cancer Res 2004;10:3728-3736, "
            "PMID 15173079)"),
        "vd_L": 0.4 * BODY_WEIGHT_KG,
        "vd_source": (
            "apparent volume of distribution ~0.4 L/kg (~17 L/m2) x 70 kg "
            "(Ostermann et al., Clin Cancer Res 2004, PMID 15173079)"),
        "f_oral": 1.0,
        "f_source": (
            "oral bioavailability ~100%; temozolomide is completely absorbed "
            "(Baker et al., Clin Cancer Res 1999;5:309-317, PMID 10037179)"),
        "confidence": "literature",
    },
    "cisplatin": {
        "half_life_s": 0.5 * 3600.0,
        "half_life_source": (
            "FREE (ultrafilterable) platinum disappears from plasma with a "
            "half-life of ~20-30 min; the ~5-day terminal half-life belongs to "
            "protein-BOUND platinum, which is not the pharmacologically active "
            "species this model transports (Himmelstein et al., Cancer Treat "
            "Rep 1981;65:1123-1126)"),
        "vd_L": 25.0,
        "vd_source": (
            "volume of distribution of free platinum ~15-25 L in adults "
            "(clinical cisplatin PK literature); 25 L used"),
        "f_oral": None,   # not orally administered
        "f_source": "cisplatin is not given orally; `route oral` is refused",
        "confidence": "literature",
    },
    "doxorubicin": {
        "half_life_s": 26.0 * 3600.0,
        "half_life_source": (
            "terminal plasma half-life ~20-30 h; 26 h used (Robert, Clin "
            "Pharmacokinet 1982;7:428-441)"),
        "vd_L": 1500.0,
        "vd_source": (
            "very large volume of distribution, ~20-25 L/kg (extensive tissue "
            "and DNA binding); 1500 L used"),
        "f_oral": None,
        "f_source": "doxorubicin is not orally bioavailable; `route oral` is refused",
        "confidence": "literature",
    },
}


@dataclass
class DoseConversion:
    """The full record of one dose -> PKConfig conversion."""

    # engine-facing outputs (the four PKConfig fields)
    plasma_peak_uM: float = 0.0
    half_life_s: float = 0.0
    n_doses: int = 0
    dose_interval_s: float = 0.0

    # window bookkeeping
    window_s: float = 0.0
    protocol_duration_s: float = 0.0
    n_administrations_total: int = 0
    n_administrations_outside_window: int = 0
    truncated: bool = False

    # arithmetic
    dose_input: str = ""
    amount_mg: float = 0.0
    route: str = "iv"
    bioavailability: float = 1.0
    volume_of_distribution_L: float = 0.0
    molecular_weight_g_per_mol: float = 0.0
    formula: str = ""
    conversions: List[str] = field(default_factory=list)
    assumptions: List[Assumption] = field(default_factory=list)

    # infusion
    infusion: Optional[Dict] = None

    # solver
    time_resolved: bool = True
    time_resolved_reason: str = ""

    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            "pk": {
                "plasma_peak": self.plasma_peak_uM,
                "half_life_s": self.half_life_s,
                "n_doses": self.n_doses,
                "dose_interval_s": self.dose_interval_s,
            },
            "window_s": self.window_s,
            "protocol_duration_s": self.protocol_duration_s,
            "n_administrations_total": self.n_administrations_total,
            "n_administrations_outside_window":
                self.n_administrations_outside_window,
            "window_truncated": self.truncated,
            "dose_input": self.dose_input,
            "amount_mg": self.amount_mg,
            "route": self.route,
            "bioavailability": self.bioavailability,
            "volume_of_distribution_L": self.volume_of_distribution_L,
            "molecular_weight_g_per_mol": self.molecular_weight_g_per_mol,
            "formula": self.formula,
            "conversions": list(self.conversions),
            "assumptions": [a.to_dict() for a in self.assumptions],
            "infusion": self.infusion,
            "time_resolved": self.time_resolved,
            "time_resolved_reason": self.time_resolved_reason,
            "notes": list(self.notes),
        }


# --------------------------------------------------------------------------- #
# Schedule expansion
# --------------------------------------------------------------------------- #
def schedule_times(every_s: float, repeat: int, rest_s: float,
                   cycles: int) -> Tuple[List[float], float]:
    """Administration times (s from t=0) and the total protocol duration.

    One cycle is ``repeat`` doses spaced ``every_s`` apart, followed by
    ``rest_s`` of no drug, so its length is ``repeat*every_s + rest_s`` — for
    Stupp adjuvant (5 x 24 h + 23 d) that is exactly the 28-day cycle.
    """
    every_s = float(every_s)
    repeat = max(int(repeat), 1)
    rest_s = max(float(rest_s), 0.0)
    cycles = max(int(cycles), 1)
    cycle_len = repeat * every_s + rest_s
    times = [c * cycle_len + r * every_s
             for c in range(cycles) for r in range(repeat)]
    return times, cycles * cycle_len


@dataclass
class WindowPlan:
    """What the observation window does to a schedule.

    Computed from the SCHEDULE ALONE — no compound, no molecular weight, no
    compilation. That is the point of splitting it out of
    :func:`convert_dose`: ``convert_dose`` is reached only once a drug has been
    compiled, so ``POST /api/protocol/validate`` and ``POST /api/protocol/run``
    — both of which plan with ``compile_drugs=False`` — never reached it and
    said NOTHING about a truncated schedule. The operator committed the compute
    at ``/run`` without being told that 29 of a Stupp course's 30
    administrations would be dropped. Truncation does not depend on the
    compound, so it does not have to wait for one.

    One arithmetic, in one place. A second copy in the planner would be free to
    drift from this one, and then the number the operator is shown stops
    matching the number the solver used.
    """

    window_s: float = 0.0
    protocol_duration_s: float = 0.0
    administrations: List[float] = field(default_factory=list)
    in_window: List[float] = field(default_factory=list)
    truncated: bool = False

    @property
    def n_administrations_total(self) -> int:
        return len(self.administrations)

    @property
    def n_administrations_outside_window(self) -> int:
        return len(self.administrations) - len(self.in_window)

    @property
    def note(self) -> str:
        """The sentence stating the truncation; ``""`` when nothing was cut."""
        if not self.truncated:
            return ""
        return (
            f"This protocol spans {_fmt(self.protocol_duration_s / 86400.0)} d "
            f"but the transport solve covers the first "
            f"{_fmt(self.window_s / 3600.0)} h. "
            f"{self.n_administrations_outside_window} of "
            f"{self.n_administrations_total} administrations fall OUTSIDE the "
            "observed window and deliver no drug in this result. The window is "
            "a compute limit (a full course at dt = 10 s is ~1.5e6 timesteps), "
            "not a claim that the later doses do not matter.")

    def to_dict(self) -> Dict:
        return {
            "window_s": self.window_s,
            "protocol_duration_s": self.protocol_duration_s,
            "n_administrations_total": self.n_administrations_total,
            "n_administrations_outside_window":
                self.n_administrations_outside_window,
            "truncated": self.truncated,
            "note": self.note,
        }


def window_plan(every_s: float, repeat: int, rest_s: float, cycles: int,
                max_window_s: float) -> WindowPlan:
    """Apply the observation window to a schedule and record what it cut."""
    admin, duration = schedule_times(every_s, repeat, rest_s, cycles)
    window = min(duration, float(max_window_s))
    # Always cover at least the first administration and one interval, so a
    # single-dose protocol still has something to observe.
    window = max(window, min(float(max_window_s), max(every_s, 3600.0)))
    return WindowPlan(
        window_s=window,
        protocol_duration_s=duration,
        administrations=list(admin),
        in_window=[t for t in admin if t < window - 1e-9],
        truncated=duration > window + 1e-9,
    )


def _uniform_spacing(times: List[float], rel_tol: float = 1e-9
                     ) -> Optional[float]:
    """The common spacing of ``times`` if it is uniform, else ``None``."""
    if len(times) <= 1:
        return None
    gaps = [b - a for a, b in zip(times, times[1:])]
    first = gaps[0]
    if first <= 0:
        return None
    for g in gaps[1:]:
        if abs(g - first) > max(rel_tol * first, 1e-9):
            return None
    return first


# --------------------------------------------------------------------------- #
# The conversion
# --------------------------------------------------------------------------- #
def convert_dose(
    dose: Quantity,
    route: str,
    every_s: float,
    repeat: int,
    rest_s: float,
    cycles: int,
    molecular_weight: float,
    drug_key: Optional[str] = None,
    over: Optional[Quantity] = None,
    max_window_s: float = DEFAULT_MAX_WINDOW_S,
) -> DoseConversion:
    """Convert a TPL protocol into ``PKConfig`` fields plus its full provenance.

    Parameters
    ----------
    dose : Quantity
        The ``dose`` statement's amount, in ``mg/m2``, ``mg/kg`` or ``mg``.
    route : {'oral', 'iv', 'ip'}
        Sets the bioavailability term.
    every_s, repeat, rest_s, cycles
        The expanded schedule (seconds / counts).
    molecular_weight : float
        g/mol, from the compiled molecule. Required: without it a milligram
        amount cannot become a molar concentration.
    drug_key : str, optional
        Generic drug name used to look up :data:`PK_LITERATURE`. ``None`` (or an
        unlisted name) falls to the generic defaults, all flagged ``assumed``.
    over : Quantity, optional
        Infusion duration. ``None`` = bolus.
    max_window_s : float
        Longest transport window this run may cover.

    Returns
    -------
    DoseConversion
    """
    out = DoseConversion()
    key = (drug_key or "").strip().lower()
    lit = PK_LITERATURE.get(key, {})
    out.route = route

    if not molecular_weight or molecular_weight <= 0:
        raise ValueError(
            "dose -> PK needs a molecular weight to turn milligrams into a "
            "molar concentration, and the compiled molecule has none "
            f"(mw={molecular_weight!r}). No default is invented."
        )
    out.molecular_weight_g_per_mol = float(molecular_weight)

    # --- 1. dose -> milligrams ------------------------------------------- #
    out.dose_input = str(dose)
    if dose.unit == "mg/m2":
        out.amount_mg = float(dose.value) * BSA_M2
        out.conversions.append(
            f"{_fmt(dose.value)} mg/m2 x {BSA_M2} m2 = "
            f"{_fmt(out.amount_mg)} mg")
        out.assumptions.append(Assumption(
            "body_surface_area", BSA_M2, "m2", "assumed",
            f"standard adult body-surface-area normalisation ({BSA_M2:g} m2), "
            "the value named in PROTOCOL_LANGUAGE.md. This script states no "
            "patient, so no patient-specific BSA exists to use."))
    elif dose.unit == "mg/kg":
        out.amount_mg = float(dose.value) * BODY_WEIGHT_KG
        out.conversions.append(
            f"{_fmt(dose.value)} mg/kg x {BODY_WEIGHT_KG} kg = "
            f"{_fmt(out.amount_mg)} mg")
        out.assumptions.append(Assumption(
            "body_weight", BODY_WEIGHT_KG, "kg", "assumed",
            f"generic {BODY_WEIGHT_KG:g} kg adult. This script states no "
            "patient."))
    else:  # mg
        out.amount_mg = float(dose.value)
        out.conversions.append(f"{_fmt(out.amount_mg)} mg (flat amount)")

    # --- 2. bioavailability ---------------------------------------------- #
    if route in ("iv", "ip"):
        out.bioavailability = 1.0
        out.assumptions.append(Assumption(
            "bioavailability", 1.0, "fraction",
            "assumed" if route == "ip" else "computed",
            "F = 1 by definition for an intravenous dose."
            if route == "iv" else
            "F = 1 ASSUMED for intraperitoneal administration. Peritoneal "
            "absorption is neither instantaneous nor complete for every "
            "compound; this is an upper bound."))
    else:  # oral
        f_oral = lit.get("f_oral")
        if f_oral is None and key in PK_LITERATURE:
            raise ValueError(
                f"`route oral` is not available for {key}: {lit.get('f_source')}"
            )
        if f_oral is None:
            out.bioavailability = 1.0
            out.assumptions.append(Assumption(
                "bioavailability", 1.0, "fraction", "assumed",
                f"no curated oral bioavailability for {drug_key or 'this drug'}; "
                "F = 1 assumed, which is an UPPER BOUND — a real oral F below 1 "
                "would lower every concentration below proportionally."))
        else:
            out.bioavailability = float(f_oral)
            out.assumptions.append(Assumption(
                "bioavailability", float(f_oral), "fraction",
                lit.get("confidence", "literature"), lit["f_source"]))

    # --- 3. volume of distribution ---------------------------------------- #
    vd = lit.get("vd_L")
    if vd:
        out.volume_of_distribution_L = float(vd)
        out.assumptions.append(Assumption(
            "volume_of_distribution", float(vd), "L",
            lit.get("confidence", "literature"), lit["vd_source"]))
    else:
        out.volume_of_distribution_L = GENERIC_VD_L_PER_KG * BODY_WEIGHT_KG
        out.assumptions.append(Assumption(
            "volume_of_distribution", out.volume_of_distribution_L, "L",
            "assumed",
            f"no curated Vd for {drug_key or 'this drug'}; "
            f"{GENERIC_VD_L_PER_KG} L/kg x {BODY_WEIGHT_KG} kg = total body "
            "water. This is a generic prior, not a measurement, and it scales "
            "the peak concentration inversely — a drug with a large Vd "
            "(doxorubicin, ~1500 L) would be over-predicted 36x by it."))

    # --- 4. plasma half-life ---------------------------------------------- #
    if lit.get("half_life_s"):
        out.half_life_s = float(lit["half_life_s"])
        out.assumptions.append(Assumption(
            "plasma_half_life", out.half_life_s, "s",
            lit.get("confidence", "literature"), lit["half_life_source"]))
    else:
        from ..core import PKConfig

        out.half_life_s = float(PKConfig().half_life)
        out.assumptions.append(Assumption(
            "plasma_half_life", out.half_life_s, "s", "flagged-default",
            f"no curated plasma half-life for {drug_key or 'this drug'}; using "
            "the engine's own tmesim.core.PKConfig default "
            f"({out.half_life_s / 3600.0:g} h). This is a placeholder, not a "
            "property of the molecule."))

    # --- 5. bolus peak ----------------------------------------------------- #
    bolus_peak = (1000.0 * out.amount_mg * out.bioavailability
                  / (out.volume_of_distribution_L
                     * out.molecular_weight_g_per_mol))
    out.formula = ("plasma_peak_uM = 1000 x amount_mg x F / (Vd_L x MW_g_per_mol)"
                   "   [one-compartment bolus, C0 = F*D/Vd]")
    out.conversions.append(
        f"1000 x {_fmt(out.amount_mg)} mg x {_fmt(out.bioavailability)} / "
        f"({_fmt(out.volume_of_distribution_L)} L x "
        f"{_fmt(out.molecular_weight_g_per_mol)} g/mol) = "
        f"{_fmt(bolus_peak)} uM per dose")

    # --- 6. schedule expansion + window ------------------------------------ #
    # The arithmetic and the sentence both live in `window_plan` so the PLANNER
    # can state the same truncation before any drug is compiled — see
    # :class:`WindowPlan`.
    wp = window_plan(every_s, repeat, rest_s, cycles, max_window_s)
    admin = list(wp.administrations)
    duration = wp.protocol_duration_s
    window = wp.window_s
    out.protocol_duration_s = duration
    out.n_administrations_total = wp.n_administrations_total
    out.window_s = window
    out.truncated = wp.truncated
    in_window = list(wp.in_window)
    out.n_administrations_outside_window = wp.n_administrations_outside_window
    if out.truncated:
        out.notes.append(wp.note)
    if not in_window:
        in_window = [0.0]
        out.notes.append(
            "no administration falls inside the window; the first dose is "
            "placed at t = 0 so the run is not silently empty.")

    # --- 7. infusion vs bolus ---------------------------------------------- #
    lam = math.log(2.0) / out.half_life_s
    if over is not None:
        t_inf, rec = over.canonical()
        out.conversions.append(f"infusion duration {rec}")
        k = int(min(max(round(t_inf / _SUBDOSE_TARGET_S), _SUBDOSE_MIN),
                    _SUBDOSE_MAX))
        sub_dt = t_inf / k
        expanded: List[float] = []
        for a in in_window:
            expanded.extend(a + j * sub_dt for j in range(k))
        expanded = sorted(t for t in expanded if t < window - 1e-9)
        spacing = _uniform_spacing(expanded) if len(expanded) > 1 else sub_dt
        if spacing is not None and abs(spacing - sub_dt) <= 1e-6 * max(sub_dt, 1.0):
            out.plasma_peak_uM = bolus_peak / k
            out.n_doses = len(expanded)
            out.dose_interval_s = sub_dt
            out.infusion = {
                "duration_s": t_inf,
                "mode": "sub-dose expansion (exact)",
                "n_subdoses_per_administration": k,
                "subdose_interval_s": sub_dt,
                "per_subdose_peak_uM": out.plasma_peak_uM,
                "explanation": (
                    f"the infusion is discretised into {k} equal micro-boluses "
                    f"of {_fmt(bolus_peak / k)} uM spaced {_fmt(sub_dt)} s "
                    "apart, so drug ARRIVES over the infusion duration and the "
                    "plasma curve rises to a plateau before decaying. This is a "
                    "SUSTAINED profile, not a bolus."),
            }
            dropped = k * len(in_window) - len(expanded)
            if dropped > 0:
                out.infusion["subdoses_outside_window"] = dropped
                out.notes.append(
                    f"{dropped} of {k * len(in_window)} infusion micro-boluses "
                    f"fall outside the {_fmt(window / 3600.0)} h window, so "
                    f"only {len(expanded)}/{k} of each administration's drug is "
                    "delivered in this result. The infusion is longer than the "
                    "observed window.")
        else:
            atten = ((1.0 - math.exp(-lam * t_inf)) / (lam * t_inf)
                     if lam * t_inf > 0 else 1.0)
            out.plasma_peak_uM = bolus_peak * atten
            out.n_doses = len(in_window)
            out.dose_interval_s = float(every_s)
            out.infusion = {
                "duration_s": t_inf,
                "mode": "peak attenuation (APPROXIMATE)",
                "attenuation_factor": atten,
                "explanation": (
                    "this infusion repeats on a longer cycle, so the full "
                    "administration train is NOT uniformly spaced and "
                    "core.PKConfig — which takes exactly one interval — cannot "
                    "represent it. The peak is instead scaled by the exact "
                    "one-compartment infusion factor "
                    f"(1 - exp(-lambda*T))/(lambda*T) = {_fmt(atten)}, which "
                    "reproduces the true Cmax of a constant-rate infusion and "
                    "preserves the AUC exactly. What it does NOT reproduce is "
                    "the RISE: the modelled curve starts at Cmax instead of "
                    "climbing to it, so early-window exposure is OVER-stated "
                    f"and late-window exposure under-stated by at most {_fmt(t_inf)} s "
                    "of lead time."),
            }
    else:
        out.plasma_peak_uM = bolus_peak
        out.n_doses = len(in_window)
        out.dose_interval_s = float(every_s)
        out.infusion = None

    out.time_resolved = True
    out.time_resolved_reason = (
        "PDEConfig.time_resolved = True. With the engine default (False) the "
        "solver replaces the whole plasma curve by its window MEAN before the "
        "tissue sees it, which makes fractionation and infusion duration "
        "bit-identical inputs — the artefact results/dose_schedule.md was "
        "regenerated to remove. Every TPL run sets it True so `every`, "
        "`repeat`, `cycles` and `over` are structurally able to change the "
        "answer."
    )
    return out
