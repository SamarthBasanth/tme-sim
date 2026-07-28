"""Carrier-mediated UPTAKE transport (OATP / OCT / OAT) as a tissue-differential factor.

Why this module exists
----------------------
Three measured-data validations in this repo converge, independently, on the
same named missing mechanism — and it is **not** efflux:

* ``results/gist_imatinib.md`` — the transport-only model is *structurally*
  falsified for tissue partition. The measured GIST tumour:normal imatinib
  ratios (0.023-0.142) fall **below the model's mathematically reachable floor**
  (``R_normal/R_tumour`` = 0.208 nominal, 0.083 at the most favourable end of the
  literature intervessel band). Falling below that floor requires MORE drug in
  normal liver than passive transport can put there. The report names the
  mechanism: *"active hepatic uptake (OCT1/OATP) plus lysosomal ion-trapping of a
  weak base in hepatocytes concentrates imatinib in normal liver."*
* ``results/moncayo_explant.md`` §4.3 — the *other* half of that diagnosis is
  ruled out. A spatially uniform ion-trapping retardation ``R`` gives
  ``C = R * C_bath`` **everywhere** at equilibrium, so a scalar trapping term can
  buy a time constant but never a spatial pattern. The trapping half of the GIST
  diagnosis cannot be the answer at the scalar level.
* ``results/efflux_impact.md`` — the tractable alternative was implemented,
  cited and scored, and it FAILED on measured data. It makes the glioma
  scorecard worse (GM 1.995x -> 2.424x, 0/8 configurations improve), it is
  *provably inert* on GIST (a scalar per-drug factor cancels exactly out of a
  within-drug tissue ratio, verified to machine precision), and its
  physiologically-correct placement — canalicular ABCB1/ABCG2 exporting imatinib
  out of normal hepatocytes into bile — has the **wrong sign**, taking the GIST
  fold error from 7.92x to 298x.

So this module implements the one named candidate left standing. It is not a
guess: it is the mechanism that both measured datasets point at, and the efflux
run's own conclusion was *"the missing mechanism is uptake transport, not
efflux"*.

Why uptake can do what efflux structurally could not
-----------------------------------------------------
Efflux failed on GIST for an algebraic reason, not an empirical one: it is a
**per-drug scalar**, and a scalar cancels out of a within-drug tissue ratio.
Uptake is different in kind, and the difference is the entire point:

**Uptake transporters are tissue-specific, so the factor is per-(drug, tissue),
not per-drug.** OATP1B1, OATP1B3 and OCT1 are hepatocyte transporters — OATP1B1
and OATP1B3 are essentially liver-exclusive — and a GIST liver metastasis is not
made of hepatocytes. A term that is large in normal liver and 1.0 in the tumour
does not cancel; it raises the *denominator* of the tumour:normal ratio, which is
exactly the direction the measurement demands and exactly the direction efflux
could not go.

It also lifts a ceiling the repo has documented from the start.
:func:`tmesim.reconstruct.barrier.effective_penetrance` clips penetrance at 1.0
and says why: *"a partition above 1 requires tissue binding or active uptake,
neither of which this term models."* Adavosertib's measured human-glioma value is
**2.46**, above that ceiling, and no passive model can reach it at any parameter
value. :func:`uptake_adjusted_penetrance` raises the ceiling to the uptake factor
itself, which is the physically correct bound: a pump can concentrate unbound
drug up to its own uptake ratio and no further.

The functional form — and it is UNIQUE, unlike efflux
------------------------------------------------------
Passive permeability ``P`` across the membrane in both directions, plus an
active influx clearance ``CL_up`` in parallel on the entry face. At steady state
the unbound tissue-to-blood ratio is

.. math::

    f_\\mathrm{uptake} = K_{p,uu}
      = \\frac{P + CL_\\mathrm{up}}{P}
      = 1 + b_\\mathrm{up}
      = \\mathrm{UR}

where ``UR`` is the **uptake ratio**, the standard published observable. Note
what this is NOT: it is not a second knob. The efflux module had to choose
between two forms (``1/ER`` and ``2/(1+ER)``, differing by 2x at high ER) because
a pump on the apical face of a polarised monolayer is read differently depending
on whether the cell interior is resolved. **For uptake the ambiguity does not
arise.** All three ways the quantity is measured collapse to the same number:

===============================================  ================================
observable                                       value
===============================================  ================================
cell accumulation fold, transfected / control    ``1 + CL_up/P``
directional ``Papp(A->B) / Papp(B->A)``          ``1 + CL_up/P``
``Papp(A->B)`` vs the passive-only monolayer     ``1 + CL_up/P``
===============================================  ================================

This is verified from the three-compartment mass balance in this module's
self-test, not asserted. **One form, no form-shopping, one fewer discretionary
choice than the efflux term had.** An in vivo tissue ``Kp,uu`` is the same
quantity read in vivo, so it needs no conversion either.

Combining transporters
----------------------
Several carriers can act on the same molecule at the same membrane, and their
clearances ADD: ``CL_tot = sum_i CL_i``. Since ``CL_i/P = UR_i - 1``,

.. math::

    \\mathrm{UR}_\\mathrm{combined} = \\sum_i \\mathrm{UR}_i - (n - 1)

which is :func:`combine_uptake_ratios` — the exact mirror of the efflux module's
parallel-pump rule, and legitimate for the same reason and under the same
condition (the ratios must share a passive permeability ``P``, i.e. same molecule,
comparable system).

.. _uptake-circularity:

The circularity guard — the most important thing in this module
----------------------------------------------------------------
An uptake ratio can be published in two very different forms, and treating them
as interchangeable would invalidate the entire test:

* an **in vitro** fold uptake in transporter-transfected vs vector-control cells.
  This is a *mechanistic input*: it measures the transporter, in a cell line,
  against a defined control. Independent of any tissue outcome.
* an **in vivo tissue** ``Kp,uu``. This is an *outcome*, and for adavosertib in
  human glioma it is the same class of measurement, in the same tissue, as the
  scorecard target itself. Using a measured human-glioma unbound partition of 3.2
  to "predict" a measured human-glioma tissue:plasma ratio of 2.46 is not a
  prediction. It is copying the answer across with a units change.

So :data:`PREFER_IN_VITRO` is the default and a record may declare an observable
``circular_against`` a named scorecard. When the only usable number for a drug is
one flagged circular, :func:`uptake_factor` returns **exactly 1.0** with reason
``"only_circular_observable"`` — the drug is treated as having no number at all.
The circular variant is still reachable via ``allow_circular=True``, but only so
that ``tmesim/validation/uptake_impact.py`` can report it as an explicitly
labelled upper bound. It is never the headline. This is the same discipline the
efflux module applied when it refused to fit the in vivo scaling factor ``SF``.

Honesty / scope — read this before using it
--------------------------------------------
* **Tissue expression is BINARY here, and that is a modelling choice.**
  :data:`TISSUE_TRANSPORTERS` says a transporter is either functionally present
  in a tissue or it is not, from published expression/localisation data. The
  alternative — scaling transport by measured protein abundance in fmol/ug —
  needs a transport-per-fmol conversion factor, and the only data we could
  calibrate it against are the two scorecards under test. That is circular for
  exactly the reason the efflux module refused ``SF``, so it is not done. The
  binary map is a literature fact; a scaled map would be a fitted parameter.
* **Substrate status is looked up, never inferred.** No structure-based OATP or
  OCT substrate classifier exists here. A drug absent from :data:`UPTAKE_TABLE`,
  or present with ``"unknown"`` status, returns 1.0.
* **A confirmed substrate with no published ratio also returns 1.0.** We do not
  invent numbers; :func:`uptake_detail` reports
  ``reason="substrate_but_no_published_ratio"`` so the gap is visible.
* **Substrate is not inhibitor.** Several drugs here (rucaparib, imatinib in some
  reports) INHIBIT carriers they may or may not be transported by. Inhibition is
  binding, not transport, and is never coded as substrate evidence.
* **The factor is >= 1.0 always.** This module models influx only. It must never
  reduce penetrance — that is the efflux module's job, and the two are
  deliberately not combined here.
* **Nothing is fitted.** Every ratio carries a PMID. Zero free parameters.

Architecture: this module reads nothing from the engine and writes nothing back.
``core``, ``physics``, ``sim.engine`` and ``treatment.drugs`` are NOT modified —
this module is purely additive and nothing in the shipped pipeline calls it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

__all__ = [
    "UptakeEntry",
    "UptakeRecord",
    "UPTAKE_TABLE",
    "TRANSPORTERS",
    "TISSUE_TRANSPORTERS",
    "TISSUES",
    "PREFER_IN_VITRO",
    "PREFER_IN_VIVO",
    "SOURCE_PREFERENCES",
    "combine_uptake_ratios",
    "transporters_in_tissue",
    "uptake_ratio_for",
    "uptake_factor",
    "uptake_detail",
    "uptake_adjusted_penetrance",
    "substrate_summary",
]


# ===========================================================================
# 0. Combination rule
# ===========================================================================
def combine_uptake_ratios(ratios: Sequence[Optional[float]]) -> Optional[float]:
    """Combine uptake ratios for carriers acting in PARALLEL on one membrane.

    Their influx clearances add, ``CL_tot = sum_i CL_i``, and ``CL_i/P = UR_i - 1``,
    hence ``UR_combined = sum_i UR_i - (n - 1)``.

    ``None`` entries mean *no usable number* (not "ratio 1") and are dropped.
    Returns ``None`` if nothing usable remains. Individual ratios below 1.0 are
    clamped to 1.0: this module models influx and never reduces exposure.
    """
    parts = [max(1.0, float(r)) for r in ratios if r is not None]
    if not parts:
        return None
    return float(sum(parts) - (len(parts) - 1))


# ===========================================================================
# 1. Transporters and where they are expressed
# ===========================================================================
#: Canonical tissue keys. ``TISSUE_GENERIC`` is the safe default: it carries NO
#: transporters, so :func:`uptake_factor` returns exactly 1.0 for it, whatever
#: the drug. A caller who does not know the tissue therefore gets no effect
#: rather than a guessed one.
TISSUE_GENERIC = "generic"
TISSUE_NORMAL_LIVER = "normal_liver"
TISSUE_LIVER_METASTASIS = "liver_metastasis"
TISSUE_NORMAL_INTESTINE = "normal_intestine"
TISSUE_INTESTINAL_TUMOUR = "intestinal_tumour"
TISSUE_NORMAL_BRAIN = "normal_brain"
TISSUE_GLIOMA = "glioma"
TISSUE_KIDNEY = "kidney"

TISSUES = (
    TISSUE_GENERIC,
    TISSUE_NORMAL_LIVER,
    TISSUE_LIVER_METASTASIS,
    TISSUE_NORMAL_INTESTINE,
    TISSUE_INTESTINAL_TUMOUR,
    TISSUE_NORMAL_BRAIN,
    TISSUE_GLIOMA,
    TISSUE_KIDNEY,
)


@dataclass(frozen=True)
class Transporter:
    """One uptake carrier, its gene, and where it is functionally expressed.

    ``tissues`` is the binary expression map discussed in the module docstring.
    A transporter contributes to :func:`uptake_factor` in a tissue if and only if
    that tissue is listed here — which is the mechanism by which the term becomes
    tissue-differential and therefore able to move a tumour:normal ratio.
    """

    name: str
    gene: str
    tissues: Tuple[str, ...]
    tissue_citation: str
    notes: str = ""

    def __post_init__(self) -> None:
        for t in self.tissues:
            if t not in TISSUES:
                raise ValueError(f"{self.name}: unknown tissue {t!r}")
        if self.tissues and not self.tissue_citation.strip():
            raise ValueError(f"{self.name}: tissue expression is a factual claim, cite it")


TRANSPORTERS: Dict[str, Transporter] = {}


def _add_transporter(t: Transporter) -> None:
    TRANSPORTERS[t.name] = t


def transporters_in_tissue(tissue: str) -> Tuple[str, ...]:
    """Names of the carriers functionally expressed in ``tissue``."""
    if tissue not in TISSUES:
        raise ValueError(f"unknown tissue {tissue!r}; expected one of {TISSUES}")
    return tuple(n for n, t in TRANSPORTERS.items() if tissue in t.tissues)


# ===========================================================================
# 2. The record types
# ===========================================================================
PREFER_IN_VITRO = "in_vitro"
PREFER_IN_VIVO = "in_vivo"
SOURCE_PREFERENCES = (PREFER_IN_VITRO, PREFER_IN_VIVO)


@dataclass(frozen=True)
class UptakeEntry:
    """One drug's relationship to ONE uptake transporter.

    Attributes
    ----------
    transporter : str
        Key into :data:`TRANSPORTERS`.
    status : {"yes", "no", "unknown"}
        Substrate status. ``"yes"`` requires a citation. ``"unknown"`` means the
        literature was searched and nothing usable was found; it is treated as a
        non-substrate downstream but reported separately.
    citation : str
        PMID + one-line justification for ``status``.
    uptake_ratio : float or None
        **In vitro** fold uptake, transporter-expressing cells over vector
        control (or the equivalent inhibitor-sensitive component). ``None`` means
        no usable published number, which is NOT the same as 1.0.
    ratio_assay, ratio_citation : str
        The system the ratio was measured in, and its PMID.
    kpuu_in_vivo : float or None
        An **in vivo tissue** unbound partition attributable to this carrier.
        Same quantity, read in vivo, so it needs no conversion.
    kpuu_tissue, kpuu_assay, kpuu_citation : str
        Which tissue and species, the design, and the PMID.
    circular_against : tuple of str
        Names of scorecards against which ``kpuu_in_vivo`` would be circular —
        i.e. the observable IS (near enough) the thing being predicted. See
        :ref:`the circularity guard <uptake-circularity>`. Non-empty means the
        in vivo number is withheld from the default configuration entirely.
    notes : str
        Everything a reader needs in order to distrust the row appropriately.
    """

    transporter: str
    status: str
    citation: str = ""
    uptake_ratio: Optional[float] = None
    ratio_assay: str = ""
    ratio_citation: str = ""
    kpuu_in_vivo: Optional[float] = None
    kpuu_tissue: str = ""
    kpuu_assay: str = ""
    kpuu_citation: str = ""
    circular_against: Tuple[str, ...] = ()
    notes: str = ""

    def __post_init__(self) -> None:
        if self.transporter not in TRANSPORTERS:
            raise ValueError(f"unknown transporter {self.transporter!r}")
        if self.status not in ("yes", "no", "unknown"):
            raise ValueError(f"{self.transporter}: status={self.status!r} must be "
                             "'yes', 'no' or 'unknown'")
        if self.status != "unknown" and not self.citation.strip():
            raise ValueError(f"{self.transporter}: status={self.status!r} asserts a "
                             "fact and must carry a citation")
        if self.uptake_ratio is not None:
            if not (self.uptake_ratio > 0.0):
                raise ValueError(f"{self.transporter}: uptake_ratio must be positive")
            if self.status != "yes":
                raise ValueError(f"{self.transporter}: has a ratio but status="
                                 f"{self.status!r}")
            if not self.ratio_citation.strip():
                raise ValueError(f"{self.transporter}: uptake_ratio is a numeric "
                                 "claim and needs ratio_citation")
        if self.kpuu_in_vivo is not None:
            if not (self.kpuu_in_vivo > 0.0):
                raise ValueError(f"{self.transporter}: kpuu_in_vivo must be positive")
            if self.status != "yes":
                raise ValueError(f"{self.transporter}: has a Kp,uu but status="
                                 f"{self.status!r}")
            if not self.kpuu_citation.strip():
                raise ValueError(f"{self.transporter}: kpuu_in_vivo is a numeric "
                                 "claim and needs kpuu_citation")

    def ratio(self, prefer: str = PREFER_IN_VITRO,
              allow_circular: bool = False) -> Optional[float]:
        """Best usable uptake ratio for this entry, honouring the circularity guard.

        Returns ``None`` when the only available observable is one flagged
        :attr:`circular_against` and ``allow_circular`` is False.
        """
        if self.status != "yes":
            return None
        vivo = self.kpuu_in_vivo
        if vivo is not None and self.circular_against and not allow_circular:
            vivo = None
        a, b = (self.uptake_ratio, vivo)
        if prefer == PREFER_IN_VIVO:
            a, b = b, a
        return a if a is not None else b

    def source(self, prefer: str = PREFER_IN_VITRO,
               allow_circular: bool = False) -> Optional[str]:
        """Which observable :meth:`ratio` actually used."""
        if self.ratio(prefer, allow_circular) is None:
            return None
        vivo_ok = (self.kpuu_in_vivo is not None
                   and (allow_circular or not self.circular_against))
        primary = self.uptake_ratio if prefer == PREFER_IN_VITRO else (
            self.kpuu_in_vivo if vivo_ok else None)
        if primary is not None:
            return prefer
        return PREFER_IN_VIVO if prefer == PREFER_IN_VITRO else PREFER_IN_VITRO

    @property
    def suppressed_circular(self) -> bool:
        """True if this entry HAS a number that the guard is withholding."""
        return (self.status == "yes" and self.uptake_ratio is None
                and self.kpuu_in_vivo is not None and bool(self.circular_against))


@dataclass(frozen=True)
class UptakeRecord:
    """One drug's literature-sourced uptake-transporter assignments."""

    drug: str
    entries: Tuple[UptakeEntry, ...] = ()
    notes: str = ""

    def __post_init__(self) -> None:
        seen = set()
        for e in self.entries:
            if e.transporter in seen:
                raise ValueError(f"{self.drug}: duplicate entry for {e.transporter}")
            seen.add(e.transporter)
        if not self.notes.strip():
            raise ValueError(f"{self.drug} must carry notes")

    def entries_in_tissue(self, tissue: str) -> Tuple[UptakeEntry, ...]:
        """Entries whose transporter is functionally expressed in ``tissue``."""
        live = set(transporters_in_tissue(tissue))
        return tuple(e for e in self.entries if e.transporter in live)

    @property
    def is_substrate(self) -> bool:
        """True if ANY carrier is a confirmed ``"yes"``."""
        return any(e.status == "yes" for e in self.entries)


UPTAKE_TABLE: Dict[str, UptakeRecord] = {}


def _add(rec: UptakeRecord) -> None:
    UPTAKE_TABLE[rec.drug] = rec


# ===========================================================================
# 3. TABLE_PLACEHOLDER
# ===========================================================================


# ===========================================================================
# 4. Lookup and the public factor
# ===========================================================================
def _resolve_name(drug: Union[str, Mapping]) -> Optional[str]:
    """Extract a lower-case lookup key from a name or a descriptor mapping."""
    if isinstance(drug, str):
        return drug.strip().lower()
    if isinstance(drug, Mapping):
        name = drug.get("name")
        return str(name).strip().lower() if name else None
    name = getattr(drug, "name", None)
    return str(name).strip().lower() if name else None


def _is_antibody(drug: Union[str, Mapping]) -> bool:
    if isinstance(drug, Mapping):
        return bool(drug.get("is_antibody", False))
    return bool(getattr(drug, "is_antibody", False))


def uptake_ratio_for(drug: Union[str, Mapping],
                     tissue: str = TISSUE_GENERIC,
                     record: Optional[UptakeRecord] = None,
                     prefer: str = PREFER_IN_VITRO,
                     allow_circular: bool = False) -> Optional[float]:
    """Combined uptake ratio for ``drug`` in ``tissue``, or ``None`` if there is none.

    Only carriers expressed in ``tissue`` contribute; their ratios are combined
    in parallel by :func:`combine_uptake_ratios`.
    """
    if prefer not in SOURCE_PREFERENCES:
        raise ValueError(f"prefer must be one of {SOURCE_PREFERENCES}, got {prefer!r}")
    rec = record if record is not None else UPTAKE_TABLE.get(_resolve_name(drug) or "")
    if rec is None:
        return None
    ratios = [e.ratio(prefer, allow_circular) for e in rec.entries_in_tissue(tissue)]
    return combine_uptake_ratios(ratios)


def uptake_detail(drug: Union[str, Mapping],
                  tissue: str = TISSUE_GENERIC,
                  prefer: str = PREFER_IN_VITRO,
                  allow_circular: bool = False,
                  uptake_ratio: Optional[float] = None) -> Dict[str, object]:
    """Full, auditable explanation of the factor from :func:`uptake_factor`.

    The ``reason`` field is the one to read when the factor is 1.0, because 1.0
    has SEVEN distinct meanings here:

    =======================================  ====================================
    reason                                   meaning
    =======================================  ====================================
    ``"not_in_table"``                       never looked up
    ``"antibody"``                           macromolecule; carriers N/A
    ``"non_substrate"``                      MEASURED not to be a substrate
    ``"unknown_substrate_status"``           searched, no usable source found
    ``"no_transporter_in_this_tissue"``      IS a substrate, but of a carrier
                                             this tissue does not express — THE
                                             tissue-differential case
    ``"substrate_but_no_published_ratio"``   confirmed substrate, no number
    ``"only_circular_observable"``           a number exists but using it would
                                             be circular against a scorecard
    ``"applied"``                            a real published ratio was applied
    =======================================  ====================================
    """
    if tissue not in TISSUES:
        raise ValueError(f"unknown tissue {tissue!r}; expected one of {TISSUES}")
    key = _resolve_name(drug)
    base: Dict[str, object] = {
        "drug": key, "tissue": tissue, "factor": 1.0, "uptake_ratio": None,
        "applied": False, "transporters_applied": (), "transporters_in_tissue":
            transporters_in_tissue(tissue), "citations": "", "notes": "",
        "source": None, "suppressed_circular": False,
    }
    if uptake_ratio is not None:
        r = max(1.0, float(uptake_ratio))
        return {**base, "factor": r, "uptake_ratio": float(uptake_ratio),
                "applied": True, "reason": "explicit_override",
                "notes": "uptake_ratio supplied by caller, table bypassed"}
    if _is_antibody(drug):
        return {**base, "reason": "antibody",
                "notes": "a ~149 kDa IgG is not a solute-carrier substrate; "
                         "antibody tissue entry is FcRn/convection-mediated"}
    rec = UPTAKE_TABLE.get(key or "")
    if rec is None:
        return {**base, "reason": "not_in_table",
                "notes": "no literature assignment exists in UPTAKE_TABLE; treated "
                         "as a non-substrate, which may be wrong"}
    base["notes"] = rec.notes
    base["citations"] = "; ".join(
        c for e in rec.entries
        for c in (e.citation, e.ratio_citation, e.kpuu_citation) if c)
    here = rec.entries_in_tissue(tissue)
    if not rec.is_substrate:
        any_unknown = any(e.status == "unknown" for e in rec.entries) or not rec.entries
        return {**base,
                "reason": "unknown_substrate_status" if any_unknown else "non_substrate"}
    if not any(e.status == "yes" for e in here):
        # The drug IS a substrate -- but of a carrier this tissue does not
        # express. This branch is the whole mechanism: it is what makes the term
        # tissue-differential and therefore able to move a tumour:normal ratio.
        return {**base, "reason": "no_transporter_in_this_tissue"}
    r = uptake_ratio_for(drug, tissue, record=rec, prefer=prefer,
                         allow_circular=allow_circular)
    if r is None:
        suppressed = any(e.suppressed_circular for e in here)
        return {**base, "suppressed_circular": suppressed,
                "reason": ("only_circular_observable" if suppressed
                           else "substrate_but_no_published_ratio")}
    applied = tuple(e.transporter for e in here
                    if e.ratio(prefer, allow_circular) is not None)
    srcs = {e.source(prefer, allow_circular) for e in here
            if e.ratio(prefer, allow_circular) is not None}
    return {**base, "factor": max(1.0, float(r)), "uptake_ratio": float(r),
            "applied": True, "reason": "applied", "transporters_applied": applied,
            "source": "+".join(sorted(s for s in srcs if s)),
            "suppressed_circular": any(e.suppressed_circular for e in here)}


def uptake_factor(drug: Union[str, Mapping],
                  tissue: str = TISSUE_GENERIC,
                  prefer: str = PREFER_IN_VITRO,
                  allow_circular: bool = False,
                  uptake_ratio: Optional[float] = None) -> float:
    """Exposure multiplier ``>= 1.0`` for carrier-mediated uptake. **1.0 = no effect.**

    Multiply a drug's passive tissue exposure by this to obtain its exposure in
    the presence of the uptake carriers that ``tissue`` expresses::

        net_exposure = passive_exposure * uptake_factor(drug, tissue)

    The factor is the **uptake ratio** ``UR = 1 + CL_up/P`` (see the module
    docstring — this form is unique, there is no second candidate), combined
    across carriers in parallel by :func:`combine_uptake_ratios`.

    **Everything unproven returns exactly 1.0**, and there are seven distinct
    ways that can happen — call :func:`uptake_detail` to find out which. In
    particular a drug that IS a confirmed substrate returns 1.0 in a tissue that
    does not express its carrier, which is the point of the whole module.

    Parameters
    ----------
    drug : str, Mapping, or object with ``.name``
        Drug name, a descriptor mapping (honouring ``"is_antibody"``), or any
        object with a ``.name`` attribute such as a
        :class:`tmesim.core.DrugParams`.
    tissue : str
        One of :data:`TISSUES`. Defaults to :data:`TISSUE_GENERIC`, which
        expresses no carriers and therefore always returns 1.0 — a caller who
        does not specify a tissue gets no effect rather than a guessed one.
    prefer : {"in_vitro", "in_vivo"}
        Which observable wins when a carrier has both. Default is the
        pre-committed ``"in_vitro"``.
    allow_circular : bool
        Permit observables flagged :attr:`UptakeEntry.circular_against`. Default
        False. See :ref:`the circularity guard <uptake-circularity>`; this exists
        for labelled sensitivity analysis only and must never be the headline.
    uptake_ratio : float, optional
        Bypass the table with an explicit ratio, for sweeps. Values below 1 are
        clamped to 1 (this module never reduces exposure).

    Returns
    -------
    float
        Multiplier ``>= 1.0``.
    """
    return float(uptake_detail(drug, tissue=tissue, prefer=prefer,
                               allow_circular=allow_circular,
                               uptake_ratio=uptake_ratio)["factor"])


def uptake_adjusted_penetrance(passive_penetrance: float,
                               drug: Union[str, Mapping],
                               tissue: str = TISSUE_GENERIC,
                               bbb_disruption: float = 0.0,
                               leak: Optional[float] = None,
                               prefer: str = PREFER_IN_VITRO,
                               allow_circular: bool = False,
                               uptake_ratio: Optional[float] = None) -> float:
    """Compose uptake with the regional paracellular leak, in the right order.

    ::

        U     = uptake_factor(drug, tissue)
        P_eff = clip( P_passive * U  +  L * bbb_disruption , 0, U )

    Two physical claims are encoded here, both deliberate.

    **1. Uptake multiplies the transcellular term only.** OATP/OCT/OAT carriers
    sit in the plasma membrane of the endothelial or parenchymal cell, so they
    act only on drug taking the transcellular route. A disrupted tight junction
    is a paracellular shunt *around* the carrier, exactly as it is around an
    efflux pump — the same composition
    :func:`tmesim.molecular.efflux.efflux_adjusted_penetrance` uses, for the same
    reason.

    **2. The ceiling rises from 1.0 to U.**
    :func:`tmesim.reconstruct.barrier.effective_penetrance` clips at 1.0 and its
    docstring states the reason: *"a partition above 1 requires tissue binding or
    active uptake, neither of which this term models."* This function models
    active uptake, so the ceiling becomes the largest unbound partition a carrier
    can sustain, which is its own uptake ratio ``U``. That is not a free choice:
    at steady state ``Kp,uu -> (P + CL_up)/P = U`` as passive equilibration
    completes, and no pump can exceed its own ratio.

    With ``U = 1`` (any non-substrate, or any tissue lacking the carrier) this
    reduces **exactly** to :func:`tmesim.reconstruct.barrier.effective_penetrance`
    — verified bit-for-bit in this module's self-test and in
    ``tmesim/validation/uptake_impact.py``.
    """
    u = uptake_factor(drug, tissue=tissue, prefer=prefer,
                      allow_circular=allow_circular, uptake_ratio=uptake_ratio)
    transcellular = float(passive_penetrance) * u
    d = float(bbb_disruption)
    if d == 0.0:
        return float(min(u, max(0.0, transcellular)))
    if leak is None:
        from ..reconstruct.barrier import L_PARACELLULAR_LEAK
        leak = L_PARACELLULAR_LEAK
    if not (-1e-12 <= d <= 1.0 + 1e-12):
        raise ValueError(f"bbb_disruption must lie in [0, 1]; got {d}")
    return float(min(u, max(0.0, transcellular + float(leak) * d)))


def substrate_summary() -> Dict[str, int]:
    """Counts of assignments in :data:`UPTAKE_TABLE`, for the report."""
    out = {"n_drugs": len(UPTAKE_TABLE), "substrate_any": 0, "non_substrate": 0,
           "unknown_only": 0, "substrate_with_number": 0,
           "substrate_without_number": 0, "circular_suppressed": 0}
    for rec in UPTAKE_TABLE.values():
        if rec.is_substrate:
            out["substrate_any"] += 1
            has = any(e.uptake_ratio is not None or e.kpuu_in_vivo is not None
                      for e in rec.entries)
            out["substrate_with_number" if has else "substrate_without_number"] += 1
            if any(e.suppressed_circular for e in rec.entries):
                out["circular_suppressed"] += 1
        elif any(e.status == "unknown" for e in rec.entries) or not rec.entries:
            out["unknown_only"] += 1
        else:
            out["non_substrate"] += 1
    for name, t in TRANSPORTERS.items():
        out[f"n_yes_{name}"] = sum(
            1 for rec in UPTAKE_TABLE.values()
            for e in rec.entries if e.transporter == name and e.status == "yes")
    return out
