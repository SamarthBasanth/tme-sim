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
    "ANTIBODY_DRUGS",
    "HEPATIC_TUMOUR_DOWNREGULATION",
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
# 3a. The transporters, and the binary tissue map
# ===========================================================================
_add_transporter(Transporter(
    name="OATP1B1", gene="SLCO1B1",
    tissues=(TISSUE_NORMAL_LIVER,),
    tissue_citation="PMID 10644574 (Abe 2000): LIVER-EXCLUSIVE, hepatocyte "
                    "basolateral (sinusoidal) membrane; cloning PMID 10358072.",
    notes="The cleanest tissue-restriction in the table, and the reason the GIST "
          "test is possible at all: a transporter found in hepatocytes and "
          "essentially nowhere else gives a normal-liver-vs-tumour contrast for "
          "free, with no fitted asymmetry.",
))

_add_transporter(Transporter(
    name="OATP1B3", gene="SLCO1B3",
    tissues=(TISSUE_NORMAL_LIVER,),
    tissue_citation="PMID 10779507 (Konig 2000): liver-restricted in NORMAL "
                    "tissue, hepatocyte basolateral, perivenous-zonated.",
    notes="CARRIES A REAL THREAT TO THE MODEL'S PREMISE and it is recorded here "
          "rather than buried. Unlike OATP1B1, OATP1B3 is reported to be "
          "abundantly expressed in gastric, colon and pancreatic CANCER (PMID "
          "11375950, Abe 2001), i.e. it is liver-restricted in normal tissue but "
          "NOT tumour-absent in general. It is coded as absent from the GIST "
          "liver metastasis because GIST is a mesenchymal (KIT-driven interstitial "
          "cell of Cajal) tumour, not an adenocarcinoma — but no measurement of "
          "OATP1B3 in GIST tissue exists in the literature, so this is an "
          "inference from lineage and it is flagged as such in the report.",
))

_add_transporter(Transporter(
    name="OCT1", gene="SLC22A1",
    tissues=(TISSUE_NORMAL_LIVER,),
    tissue_citation="PMID 9187257 (Zhang 1997): HUMAN OCT1 is liver-predominant. "
                    "Species caveat stated in the source itself: rodent Oct1 is "
                    "also kidney, intestine and colon, so mouse data do not "
                    "transfer tissue-for-tissue.",
    notes="Quantitatively liver-skewed: SLC22A1 mRNA in human liver is ~500-fold "
          "higher than in CML mononuclear/CD34+ cells (PMID 24352644, primary "
          "text). NOTE: a dedicated hepatocyte-sinusoidal-membrane "
          "immunolocalisation paper was NOT retrieved, so the sub-cellular "
          "localisation claim is deliberately not made here — only the "
          "tissue-level predominance, which is what this module uses.",
))

_add_transporter(Transporter(
    name="OATP1A2", gene="SLCO1A2",
    tissues=(TISSUE_NORMAL_BRAIN, TISSUE_GLIOMA, TISSUE_KIDNEY),
    tissue_citation="PMID 15632119 (Lee 2005): human BRAIN CAPILLARY ENDOTHELIUM, "
                    "luminal membrane, plus renal distal nephron. In glioma "
                    "specifically PMID 16357150 (Bronger 2005 Cancer Res), 61 human "
                    "glioma specimens: protein in the luminal membrane of the "
                    "endothelial cells forming the BBB and the blood-TUMOUR "
                    "barrier, but NOT in the glioma cells themselves.",
    notes="Two things to hold. (1) Bronger 2005 places OATP1A2 on the BARRIER, not "
          "on tumour cells — which is exactly where a penetrance multiplier "
          "belongs, so the placement in `uptake_adjusted_penetrance` is supported "
          "rather than assumed. (2) CONTRADICTED in part: PMID 25132355 (Gao 2015) "
          "localises OATP1A2 to NEURONS and puts OATP2B1, not OATP1A2, in brain "
          "capillary endothelium. (3) INTESTINE IS DELIBERATELY EXCLUDED. OATP1A2 "
          "is widely asserted to sit in the enterocyte apical membrane and to "
          "mediate duodenal absorption of imatinib (PMID 21508937), but no "
          "verified human intestinal localisation citation could be retrieved for "
          "this table, and the module forbids an uncited tissue claim. This "
          "exclusion is not neutral: including intestine would move the GIST "
          "intestinal control pair. It is reported as an open item, not hidden.",
))

_add_transporter(Transporter(
    name="OAT1", gene="SLC22A6",
    tissues=(TISSUE_KIDNEY,),
    tissue_citation="PMID 11912245 (Motohashi 2002): kidney proximal tubule, "
                    "basolateral membrane, cortex.",
    notes="Renal. Present in this table only so that a drug whose real carrier "
          "biology is organic-anion secretion (cefazolin) is not silently recorded "
          "as having no carrier at all. It contributes nothing in brain or liver, "
          "which is the correct behaviour and not a limitation.",
))

_add_transporter(Transporter(
    name="OAT3", gene="SLC22A8",
    tissues=(TISSUE_KIDNEY,),
    tissue_citation="PMID 11912245 (Motohashi 2002): kidney proximal tubule, "
                    "basolateral membrane, cortex.",
    notes="Renal, as OAT1. Same rationale for inclusion.",
))

_add_transporter(Transporter(
    name="OCT2", gene="SLC22A2",
    tissues=(TISSUE_KIDNEY,),
    tissue_citation="PMID 11912245 (Motohashi 2002): kidney proximal tubule "
                    "basolateral membrane, cortex and medullary ray; the same "
                    "paper reports hOCT1 protein UNDETECTABLE in kidney, which is "
                    "the human OCT1/OCT2 division of labour.",
    notes="Also demonstrated in COCHLEAR HAIR CELLS (PMID 20110413) — the "
          "cisplatin ototoxicity route, and the reason OCT1/2 double-knockout mice "
          "show no ototoxicity. Cochlea is not a tissue this repo models, so that "
          "localisation is recorded and not encoded.",
))

_add_transporter(Transporter(
    name="RFC", gene="SLC19A1",
    tissues=(TISSUE_NORMAL_LIVER, TISSUE_LIVER_METASTASIS,
             TISSUE_NORMAL_INTESTINE, TISSUE_INTESTINAL_TUMOUR,
             TISSUE_NORMAL_BRAIN, TISSUE_GLIOMA, TISSUE_KIDNEY),
    tissue_citation="PMID 12144527 (Whetstine 2002 Biochem J): hRFC transcripts "
                    "'ubiquitously expressed' across a 76-tissue array, highest in "
                    "placenta and lowest in skeletal muscle.",
    notes="DELIBERATELY LISTED IN EVERY TISSUE, and that is the point. RFC is the "
          "canonical folate/antifolate importer and it is UBIQUITOUS, so it is "
          "present on both sides of any tumour-versus-normal comparison and "
          "CANCELS out of the ratio — the same algebra that made a scalar efflux "
          "factor inert on GIST. A ubiquitous carrier is not a tissue-differential "
          "mechanism no matter how strong it is. Including it makes that visible "
          "in the table rather than leaving it as an unstated assumption.",
))

_add_transporter(Transporter(
    name="CTR1", gene="SLC31A1",
    tissues=(),
    tissue_citation="",
    notes="NO TISSUES ENCODED, ON PURPOSE. CTR1 is universally described as the "
          "principal cisplatin uptake route, but no verifiable primary human "
          "tissue-distribution citation could be retrieved, and this module "
          "forbids an uncited tissue claim. It therefore contributes in no tissue "
          "and cisplatin's factor is 1.0 everywhere. That is the correct behaviour "
          "for an unverified premise, and it is reported rather than papered over. "
          "The substrate evidence itself (PMID 12370430) is also accumulation- and "
          "resistance-based rather than a direct flux measurement, with no "
          "Km/Vmax.",
))

_add_transporter(Transporter(
    name="SLC22A16", gene="SLC22A16",
    tissues=(),
    tissue_citation="",
    notes="NO TISSUES ENCODED. Cloned as a carnitine transporter 'specific to "
          "human testis' — Sertoli cells and epididymal luminal membrane (PMID "
          "12089149) — and testis is not a tissue this repo models. Recorded "
          "because it carries the only clean Km for doxorubicin uptake (5.2 +/- "
          "0.4 uM, PMID 15963465), which would otherwise look like a usable number "
          "for a tumour compartment. CONTRADICTION FLAGGED: PMID 15963465 detects "
          "it in primary acute-leukaemia patient samples, so 'testis-specific' may "
          "be too strong — but not strong enough in the other direction to justify "
          "encoding a tissue.",
))

_add_transporter(Transporter(
    name="GLUT3", gene="SLC2A3",
    tissues=(TISSUE_NORMAL_BRAIN, TISSUE_GLIOMA),
    tissue_citation="PMID 1505464 (Nagamatsu 1992): neuron-predominant in human "
                    "cerebral cortex by double labelling. Glioblastoma / "
                    "brain-tumour-initiating cells PMID 23995067 (Flavahan 2013 "
                    "Nat Neurosci).",
    notes="Do NOT cite the original cloning paper (PMID 3170580) for "
          "brain-specificity — it reports GLUT3 in most tissues. The "
          "brain-predominance claim rests on the two PMIDs above. GLUT3 is in "
          "this table for one drug only (temozolomide) and on one unreplicated "
          "2025 report; see that row.",
))


#: Measured downregulation of the hepatic uptake carriers in a liver tumour vs
#: its own adjacent non-tumour liver. Used by
#: ``validation/uptake_impact.py::hcc_downregulation_check`` to test the model's
#: central premise against data from outside this repo's two scorecards.
#:
#: HCC is chosen because it is the HARDEST case for the premise: it is a tumour
#: that actually derives from hepatocytes. A GIST metastasis has no hepatocyte
#: lineage at all, so a premise that survives HCC survives GIST a fortiori.
#: **No number from these papers enters the model** — this is a directional
#: check, not a calibration.
HEPATIC_TUMOUR_DOWNREGULATION = (
    {"gene": "SLC22A1", "transporter": "OCT1", "direction": "down",
     "finding": "protein significantly DOWN in HCC vs paired adjacent normal "
                "liver, n = 71 paired, tissue-microarray IHC, P < 0.0001",
     "citation": "PMID 22196450 (Schaeffeler 2011 Genome Med)"},
    {"gene": "SLC22A1", "transporter": "OCT1", "direction": "down",
     "finding": "39 advanced HCC biopsies: absent 10%, weak 38%; and of the "
                "tumours that DO express it only ~1/3 have it at the plasma "
                "membrane (15% PM vs 36% cytosolic) — mislocalisation, not just "
                "loss of message",
     "citation": "PMID 28178663 (Geier/Marin 2017 Oncotarget, TRANSFER)"},
    {"gene": "SLC22A1", "transporter": "OCT1", "direction": "down",
     "finding": "aberrant splice variants that fail to reach the plasma membrane; "
                ">= 1 inactivating variant in 48% of HCC and 40% of "
                "cholangiocarcinoma",
     "citation": "PMID 23532667 (Herraez/Marin 2013 Hepatology)"},
    {"gene": "SLCO1B3", "transporter": "OATP1B3", "direction": "down",
     "finding": "decreased in 60% of 13 paired HCC / adjacent non-tumour samples, "
                "mRNA and protein, via HNF3-beta repression",
     "citation": "PMID 14739090 (Vavricka 2004 J Hepatol)"},
    {"gene": "SLCO1B1", "transporter": "OATP1B1", "direction": "not_significant",
     "finding": "CONTRADICTS THE CLEAN STORY: OATP1B1 (OATP-C) was NOT "
                "significantly decreased in the same 13 paired samples",
     "citation": "PMID 14739090 (Vavricka 2004 J Hepatol)"},
)


# ===========================================================================
# 3b. The drug table
# ===========================================================================
# --- the GIST anchor, and the row that decides the whole exercise -----------
_add(UptakeRecord(
    drug="imatinib",
    entries=(
        UptakeEntry(
            transporter="OCT1", status="no",
            citation="PMID 24352644 (Nies 2014 Clin Cancer Res): NO increase over "
                     "vector control in FIVE OCT1-expressing lines (HEK, MDCK, "
                     "V79, CHO) and in oocytes, with the MPP+ positive control up "
                     "8-fold in the SAME oocytes; decynium-22 and prazosin failed "
                     "to inhibit; shRNA knockdown no effect; and Oct1(-/-) mice "
                     "show NO change in plasma or HEPATIC imatinib. Corroborated "
                     "by PMID 18483382 (Hu/Sparreboom 2008): oocytes not a "
                     "substrate (P = 0.21), HEK293 only 1.20-fold, and "
                     "reduced-function SLC22A1 variants had no effect on oral "
                     "clearance in 74 patients (P = 0.99). Also PMID 30734690 "
                     "(Blanc Mettral 2019): not a substrate of OCT1, OCT2 or OCT3.",
            notes="THIS ROW FALSIFIES THE MECHANISM results/gist_imatinib.md "
                  "NAMED, and it is the single most consequential entry in the "
                  "module. The large clinical 'OCT-1 activity' literature (PMID "
                  "16597591 White 2006, PMID 17761829 White 2007, PMID 17568400 "
                  "Wang 2008, PMID 18559609 Hiwase 2008 with a prazosin-inhibitable "
                  "IUR of 6.38 ng/200,000 cells) is real and clinically predictive "
                  "— but every one of those studies measures an "
                  "INHIBITOR-DIFFERENCE in non-transfected cells, never transport "
                  "by OCT1. PMID 26108972 (Burger/Sparreboom 2015 Mol Pharmacol) "
                  "shows in HEK293, K562, SD-1 and GIST-T1 that the "
                  "prazosin-sensitive component is LYSOSOMAL SEQUESTRATION, not "
                  "SLC22A1. And imatinib is a potent OCT1 INHIBITOR (IC50 0.095 "
                  "uM, PMID 24352644), which is exactly the property that makes "
                  "inhibitor-difference assays read as if it were a substrate. "
                  "Inhibitor is not substrate; that distinction is the whole "
                  "difference between this row being 'no' and being 'yes'.",
        ),
        UptakeEntry(
            transporter="OATP1B1", status="unknown",
            notes="GENUINELY DISPUTED between two direct transfected-cell studies, "
                  "so coded unknown rather than resolved in either direction. FOR: "
                  "PMID 23340295 (Zimmerman 2013 Clin Cancer Res), HEK293 stably "
                  "transfected, 0.1 uM, 30 min — imatinib among the TKIs with "
                  "statistically significant uptake over vector in BOTH OATP1B1 "
                  "and OATP1B3. AGAINST: PMID 18483382 (Hu 2008), oocytes at 0.2 "
                  "uM — OATP1B1 negative while OATP1B3 was positive.",
        ),
        UptakeEntry(
            transporter="OATP1B3", status="yes",
            citation="PMID 23340295 (Zimmerman 2013 Clin Cancer Res): HEK293 "
                     "stably transfected with OATP1B3, [3H]imatinib 0.1 uM, 30 min, "
                     "statistically significant uptake over vector control; the "
                     "same paper states imatinib was 'previously identified as a "
                     "substrate of OATP1B3'. Concordant with PMID 18483382 (Hu "
                     "2008), oocytes 0.2 uM.",
            notes="CONFIRMED SUBSTRATE, NO USABLE NUMBER — the exact situation "
                  "doxorubicin was in for the efflux table. Zimmerman publishes "
                  "imatinib's uptake only as bars in Figure 1A with asterisks; no "
                  "fold appears in the text, the tables or the figure caption, and "
                  "no Km/Vmax exists for imatinib on any OATP. Reading a value off "
                  "the figure would be a guess, so the factor is 1.0. What the "
                  "figure DOES license is an upper BOUND, and the report uses it "
                  "as such rather than as a value: sorafenib was the highest TKI "
                  "in that panel at ~1.5-fold, so imatinib is necessarily below "
                  "1.5-fold. A bound that falsifies is legitimate where a bound "
                  "used as a value would not be.",
        ),
        UptakeEntry(
            transporter="OATP1A2", status="yes",
            citation="PMID 21508937 (Eechoute/Sparreboom 2011 Clin Pharmacol Ther): "
                     "Xenopus oocytes, [3H]imatinib, uptake over water-injected "
                     "control; rosuvastatin-inhibitable at both pH values tested. "
                     "Independently PMID 21633340 (Yamakawa 2011 Clin Pharmacol "
                     "Ther): significantly enhanced uptake in OATP1A2-transfected "
                     "HEK293 (P = 0.002), naringin-inhibitable in HEK293, Caco-2 "
                     "and K562.",
            uptake_ratio=1.5,
            ratio_assay="Xenopus oocyte uptake of [3H]imatinib vs water-injected "
                        "control at pH 7.4 (1.5-fold, P = 0.0008); rises to "
                        "2.4-fold at pH 5.0 (P = 0.0001), i.e. the carrier is "
                        "pH-sensitive and the physiological-pH value is recorded "
                        "here as the conservative one",
            ratio_citation="PMID 21508937 (Eechoute 2011 Clin Pharmacol Ther)",
            notes="THE ONLY NUMERIC UPTAKE RATIO ANYWHERE FOR IMATINIB — and it is "
                  "for a transporter that is NOT hepatic, so it cannot rescue the "
                  "GIST prediction. OATP1A2 sits in brain capillary endothelium "
                  "and renal distal nephron (PMID 15632119), not in hepatocytes. "
                  "Two further warnings travel with it. (1) THE CLINICAL "
                  "TRANSLATION FAILED, in the source paper itself: no association "
                  "with OATP1A2 variants in 94 patients, and rosuvastatin had NO "
                  "effect on imatinib PK in 12 GIST patients. The authors' own "
                  "conclusion is that 'this transporter by itself is unlikely to "
                  "contribute substantially to the absorption profiles of imatinib "
                  "in humans'. (2) The 1.5-fold is an oocyte number at 0.2-1 uM.",
        ),
    ),
    notes="THE ANCHOR CASE, AND THE MECHANISM IT WAS SUPPOSED TO SUPPLY IS "
          "FALSIFIED. results/gist_imatinib.md named 'active hepatic uptake "
          "(OCT1/OATP)' as the missing physics. Verified rather than assumed, the "
          "OCT1 half is contradicted by five transfected cell lines, oocytes with "
          "a working positive control, shRNA knockdown and Oct1(-/-) mice whose "
          "HEPATIC imatinib is indistinguishable from wild-type (PMID 24352644). "
          "The OATP half survives as substrate status but carries no publishable "
          "number and is bounded below 1.5-fold, against a required ~11.7-fold. "
          "Imatinib's real intracellular-accumulation mechanism, per the group "
          "that dismantled the OCT1 story, is LYSOSOMAL SEQUESTRATION (PMID "
          "26108972) — tested in GIST-T1 among other lines. That is the "
          "redirection this run produces and it is developed in "
          "results/uptake_impact.md. Note also that no measurement of ANY uptake "
          "carrier in GIST tumour tissue exists: PubMed 'SLCO1A2 AND GIST' returns "
          "zero, and the one GIST transporter review (PMID 21163869) covers only "
          "ABCB1/ABCC1. Human liver:plasma and human tumour Kp,uu for imatinib are "
          "both UNKNOWN; the available figures are mouse (liver:plasma AUC 2.29, "
          "PMID 21635257), and Oct1(-/-) mice show that number is not "
          "OCT1-dependent.",
))

# --- the glioma anchor ------------------------------------------------------
_add(UptakeRecord(
    drug="adavosertib",
    entries=(
        UptakeEntry(
            transporter="OATP1A2", status="yes",
            citation="PMID 28928160 (Li 2017 Clin Cancer Res): cellular uptake of "
                     "AZD1775 in HEK293 overexpressing OATP1A2 significantly higher "
                     "than vector control (P = 0.001), reduced to control level by "
                     "100 uM sulfobromophthalein (P = 0.005), with an "
                     "estrone-3-sulfate positive control (P < 0.001). OATP2B1 and "
                     "OAT3 negative in the same experiment.",
            kpuu_in_vivo=3.2,
            kpuu_tissue="human glioblastoma tumour",
            kpuu_assay="unbound tumour-to-plasma concentration ratio, median 3.2 "
                       "(range 1.3-24.4) in n = 20 patients, phase 0 single oral "
                       "dose with resection at 2-26 h; fu MEASURED in the patients' "
                       "own specimens (fu,brain 0.067 +/- 0.013, fu,plasma 0.182 "
                       "+/- 0.024), Kp,total median 9.1",
            kpuu_citation="PMID 28928160 (Li 2017 Clin Cancer Res)",
            circular_against=("glioma_regional_pk",),
            notes="A PROPERLY CONTROLLED POSITIVE WITH NO PUBLISHABLE NUMBER. The "
                  "uptake fold appears ONLY as bars in Figure 4A; no value is given "
                  "in the text, the tables or the abstract, and no Km/Vmax is "
                  "reported. PubMed 'AZD1775 AND OATP1A2' returns exactly one "
                  "record, this one, and it publishes only P-values. Under this "
                  "module's rule a figure is not a number, so uptake_ratio stays "
                  "None — the same treatment doxorubicin got in the efflux table. "
                  "THE Kp,uu IS WITHHELD AS CIRCULAR: it is a human glioblastoma "
                  "tumour-to-plasma ratio, i.e. the same quantity in the same "
                  "tissue as the scorecard target. Worse, it is not even an "
                  "independent measurement of the mechanism: Li 2017's PBPK uptake "
                  "clearance of 44 L/h was 'assigned based on sensitivity "
                  "analysis', i.e. TUNED so that the model reproduces 3.2 (passive "
                  "plus efflux alone give Kp,uu 0.78-0.92). And no pH experiment "
                  "was ever run on OATP1A2 — all the pH work in that paper is on "
                  "MDCKII efflux — so the widely repeated 'acidic tumour "
                  "microenvironment activates OATP1A2 uptake' story is a PBPK "
                  "inference, not a measurement.",
        ),
    ),
    notes="THE PRE-REGISTERED TARGET OF THIS RUN, AND IT CANNOT BE SCORED. The "
          "efflux run identified adavosertib as the drug where efflux was most "
          "clearly the wrong mechanism and named uptake as the right one. Uptake "
          "IS the right mechanism — the substrate status is demonstrated with "
          "inhibitor and positive controls — but the model needs a number and no "
          "number exists that is both published and non-circular. Two independent "
          "human cohorts agree the drug accumulates (Li 2017 Kp,uu median 3.2, "
          "n = 20; Lee 2025 PMID 39820282 intratumoral:plasma 4.18 +/- 3.36 "
          "contrast-enhancing and 0.74 +/- 0.63 non-enhancing), and note that the "
          "non-enhancing value BELOW 1.0 means the accumulation is a "
          "blood-TUMOUR-barrier phenomenon in leaky core, not a BBB phenomenon. "
          "Contradicting mouse data are on record and unresolved: brain:plasma 1.0 "
          "in wild-type mice (PMID 29147815) and normal-brain:blood 0.04-0.05 "
          "(PMID 25609063), a 9-25x interspecies gap neither human paper explains. "
          "Placement note: PMID 16357150 puts OATP1A2 in the glioma BBB/blood-tumour "
          "BARRIER endothelium and explicitly NOT in glioma cells, which is exactly "
          "where a penetrance multiplier acts — so the composition used here is "
          "supported rather than assumed.",
))

# --- first-line GBM alkylators ----------------------------------------------
_add(UptakeRecord(
    drug="temozolomide",
    entries=(
        UptakeEntry(
            transporter="GLUT3", status="yes",
            citation="PMID 40813374 (Diao 2025 Cell Death Discov): GLUT3 knockdown "
                     "in U343-MG reduced TMZ uptake by 54.3%; stable GLUT3 "
                     "overexpression in U251-MG increased uptake ~70% vs "
                     "empty-vector control, LC-MS/MS with a 2H3-TMZ internal "
                     "standard; site-directed mutagenesis of 8 residues; mouse "
                     "xenograft confirmation.",
            uptake_ratio=1.70,
            ratio_assay="cellular TMZ uptake, stable GLUT3-overexpressing U251-MG "
                        "vs stable empty-vector control (+70%), LC-MS/MS with "
                        "2H3-TMZ internal standard — the canonical "
                        "transfected-versus-vector form of the observable",
            ratio_citation="PMID 40813374 (Diao 2025 Cell Death Discov)",
            notes="THE WEAKEST NUMERIC ROW IN THE TABLE, INCLUDED WITH ITS "
                  "PROBLEMS ON DISPLAY, and leave-one-out tested in the report. "
                  "(1) Single 2025 paper, no independent replication found. (2) The "
                  "companion affinity data are internally suspect: the paper "
                  "reports GLUT3 glucose Kd 1-10 uM against a well-established "
                  "glucose Km of ~1.4 mM, a ~1000x discrepancy, so the MST numbers "
                  "should not be trusted — the uptake numbers are a different "
                  "measurement and are the only thing used here. (3) Kd is binding, "
                  "not transport; no Km/Vmax for TMZ transport exists anywhere. "
                  "(4) The knockdown arm implies a larger endogenous ratio of "
                  "1/(1-0.543) = 2.19; the smaller overexpression figure is used "
                  "because it matches the canonical observable and is the "
                  "conservative one. (5) The paper's own framing claim that TMZ "
                  "uptake 'is regulated by Wnt signalling' could not be traced to a "
                  "primary source and is not relied on.",
        ),
    ),
    notes="Corrects a common assumption in the other direction from the efflux "
          "table. TMZ is routinely described as entering cells purely by passive "
          "diffusion; a 2025 report contradicts that with a carrier. Everything "
          "else is UNKNOWN rather than negative: no OATP, OCT, OAT, MCT, ENT, PEPT "
          "or OCTN assay for TMZ exists. LAT1 is coded nowhere because the only "
          "LAT1 result (PMID 38052257) concerns an ENGINEERED amino-acid prodrug, "
          "not temozolomide. The active species MTIC is generated intracellularly "
          "by spontaneous, pH-driven ring opening (PMID 8049205) rather than by a "
          "localised enzyme, and no MTIC transport study exists. The human CSF "
          "anchor, AUC_CSF/AUC_plasma = 20% in n = 35 (PMID 15173079), is total "
          "drug and is not attributable to any transporter, so it is NOT used here.",
))

_add(UptakeRecord(
    drug="lomustine",
    entries=(
        UptakeEntry(
            transporter="OATP1B1", status="no",
            citation="PMID 557367 (Begleiter 1977 Cancer Res): 'uptake of CCNU also "
                     "occurs by passive diffusion' — temperature-independent, "
                     "Na-independent, unaffected by NaF/NaCN/2,4-DNP, and NOT "
                     "self-inhibited by excess unlabelled CCNU or BCNU, i.e. "
                     "non-saturable. A non-saturable, energy-independent uptake "
                     "excludes carrier mediation generically.",
            notes="Recorded against OATP1B1 as the representative carrier because "
                  "the source excludes carrier mediation as a class rather than "
                  "testing named transporters. TRAP IN THIS PAPER: total "
                  "RADIOACTIVITY uptake for ring-labelled CCNU WAS "
                  "temperature-dependent, which the authors attribute to drug "
                  "decomposition preceding uptake. Anyone citing "
                  "'temperature-dependent CCNU uptake' has misread it.",
        ),
    ),
    notes="An improvement on the efflux table, where lomustine was 'unknown' in "
          "both cells: here there is a genuine, well-designed MEASURED NEGATIVE "
          "from 1977 that excludes carrier-mediated uptake as a class. Corroborated "
          "by PMID 6227421 (intracellular CCNU constant at 85 pmol/uL in L1210 over "
          "10-60 min). Discard the DrugBank-derived claim that lomustine is a major "
          "CYP2D6 substrate: it is unverifiable and contradicts the GLEOSTINE label, "
          "which states the elimination pathways have not been characterised.",
))

_add(UptakeRecord(
    drug="carmustine",
    entries=(
        UptakeEntry(
            transporter="OATP1B1", status="no",
            citation="PMID 557367 (Begleiter 1977 Cancer Res): 'uptake of intact "
                     "BCNU was by passive diffusion'; cell/medium distribution ratio "
                     "0.2-0.6, i.e. BELOW 1 and therefore not concentrative, and no "
                     "saturation. Replicated by PMID 1692254: steady-state "
                     "cell:medium 0.52 +/- 0.12 in SK-MG-1 human glioma.",
            notes="As for lomustine, recorded against a representative carrier "
                  "because the source excludes carrier mediation as a class. A "
                  "cell:medium ratio below 1 is a particularly clean negative for an "
                  "UPTAKE term: a concentrative carrier cannot produce it.",
        ),
    ),
    notes="CRITICAL MIS-ATTRIBUTION RISK, AVOIDED AND RECORDED. PMID 8341272 (Noe "
          "1993 Mol Pharmacol) reports the only Km/Vmax in the entire nitrosourea "
          "literature — Km 2.39 +/- 0.37 mM, Vmax 236 +/- 53 pmol/uL/s, epinephrine "
          "Ki 163 +/- 15 uM, via the extraneuronal monoamine transporter "
          "(EMT/OCT3, SLC22A3) — but it is for SarCNU, an EXPERIMENTAL ANALOGUE, "
          "not carmustine. The same paper explicitly contrasts SarCNU with BCNU, "
          "'which enters cells primarily by passive diffusion'. Importing that "
          "Km/Vmax would attribute a different molecule's carrier to this one. "
          "(Even within the SarCNU series it is contested: PMID 11181897 found no "
          "correlation between EMT expression and SarCNU efficacy across nine "
          "xenografts.) Also do not miscite PMID 6631721: MeCCNU/semustine "
          "impairing renal PAH transport in rats is transporter INJURY, not "
          "substrate status.",
))

_add(UptakeRecord(
    drug="procarbazine",
    entries=(),
    notes="GENUINELY UNKNOWN for every carrier, and honestly so. The single "
          "on-topic primary source is PMID 708471 (Lam 1978 Biochem Pharmacol), "
          "'On the mechanism of uptake of procarbazine by L5178Y lymphoblasts in "
          "vitro' — PubMed carries no abstract and the full text could not be "
          "retrieved, so its conclusion cannot be stated. The paper exists and is "
          "exactly on topic; that is all that can be said. Do NOT infer 'passive by "
          "analogy' with the nitrosoureas: the same group published "
          "carrier-mediated transport for melphalan in the same cell line (PMID "
          "837375). Procarbazine is small, lipophilic and equilibrates rapidly with "
          "CSF, which is consistent with passive diffusion — but that is "
          "physicochemical REASONING, not a transport measurement, and it does not "
          "exclude a carrier. Modelling caveat: procarbazine is a prodrug requiring "
          "CYP/MAO activation (PMID 2908840, PMID 1316811) and transport of the "
          "active azoxy metabolites has never been studied at all. Its MAO "
          "INHIBITION is enzyme pharmacology and must not be read as transporter "
          "substrate status.",
))

# --- the rest of GBM_DRUGS ---------------------------------------------------
_add(UptakeRecord(
    drug="methotrexate",
    entries=(
        UptakeEntry(
            transporter="OATP1A2", status="yes",
            citation="PMID 16702441 (Badagnani 2006 JPET): saturable "
                     "OATP1A2-mediated MTX transport in Xenopus oocytes, Km 457 +/- "
                     "118 uM, Vmax 17.5 +/- 4.9 pmol/oocyte/60 min; acidic pH "
                     "stimulated MTX uptake up to 7-fold.",
            notes="A MEASURED Km IS NOT AN UPTAKE RATIO, and the distinction is "
                  "load-bearing. This module's factor is `UR = 1 + CL_up/P`, which "
                  "needs the transporter clearance AND the passive permeability in "
                  "the same system. Km/Vmax alone gives neither, and no "
                  "transfected-versus-vector fold is published for MTX on OATP1A2. "
                  "The 'up to 7-fold' figure is a pH RATIO within the transfected "
                  "cells, not a transfected-over-control ratio, and using it would "
                  "be a category error. (The 30-fold value in that abstract is for "
                  "the estrone-3-sulfate probe, not MTX — a trap.) Factor stays 1.0.",
        ),
        UptakeEntry(
            transporter="OATP1B1", status="yes",
            citation="PMID 26474710 (Kim 2015 Mol Pharm): MTX uptake in "
                     "HEK293-OATP1B1, characterised via inhibitor IC50s "
                     "(rifampicin 0.9 uM at 1B1, 0.3 uM at 1B3; telmisartan 6.7 / "
                     "7.9 uM).",
            notes="Substrate status is implied by the assay design rather than by a "
                  "reported Km or fold; no usable ratio. Population-PK support "
                  "exists but is NOT transport evidence: SLCO1B1 variants explain "
                  "10.7% of MTX clearance variance in childhood ALL (PMID "
                  "22147369).",
        ),
        UptakeEntry(
            transporter="OATP1B3", status="yes",
            citation="PMID 11375950 (Abe 2001 Gastroenterology): saturable, "
                     "dose-dependent MTX transport in Xenopus oocytes; "
                     "OATP1B3/LST-2 transfection increased MTX sensitivity. No Km "
                     "reported.",
        ),
        UptakeEntry(
            transporter="OAT1", status="yes",
            citation="PMID 12130730 (Takeda 2002 JPET): hOAT1-mediated MTX "
                     "transport, Km 553.8 uM.",
        ),
        UptakeEntry(
            transporter="OAT3", status="yes",
            citation="PMID 12130730 (Takeda 2002 JPET): hOAT3-mediated MTX "
                     "transport, Km 21.1 uM (hOAT4 17.8 uM).",
        ),
        UptakeEntry(
            transporter="RFC", status="yes",
            citation="PMID 36071163 (Dang 2022 Nature): cryo-EM structure of human "
                     "RFC WITH METHOTREXATE BOUND; 'major importer of folates ... "
                     "as well as chemotherapeutic agents such as methotrexate'.",
            notes="No Km could be verified from any primary source read for this "
                  "table — the commonly quoted low-micromolar value did not appear "
                  "in any retrievable paper, so it is marked unverified rather than "
                  "hard-coded. In any case RFC is UBIQUITOUS and therefore cancels "
                  "out of any tumour:normal ratio (see the RFC transporter note).",
        ),
    ),
    notes="THE BEST-COVERED DRUG IN THE TABLE AND STILL NO USABLE FACTOR. "
          "Methotrexate is a confirmed substrate of six distinct carriers spanning "
          "liver, kidney, BBB and every tissue — and not one of them yields the "
          "quantity this module needs. Km/Vmax and inhibitor IC50s are abundant; "
          "transfected-over-control uptake FOLDS are absent. That is a precise "
          "statement of the coverage gap: the field measures affinity and "
          "inhibition, while a penetrance multiplier needs a ratio against passive "
          "permeability in the same system. PCFT/SLC46A1 is left out entirely: 27 "
          "PubMed records link MTX and PCFT but no primary transport paper could be "
          "verified, and an unverified carrier is not encoded. Also do not cite "
          "PMID 20448812 for OATP1B1/1B3 MTX transport — its probe is "
          "fluorescein-methotrexate, a conjugate, not MTX.",
))

_add(UptakeRecord(
    drug="doxorubicin",
    entries=(
        UptakeEntry(
            transporter="OATP1A2", status="yes",
            citation="PMID 24554572 (Durmus 2014 Int J Cancer): doxorubicin uptake "
                     "in OATP1A2-overexpressing HEK293 vs control. Independently "
                     "confirmed in a HeLa vaccinia single-transporter system, PMID "
                     "27777271 (Lee 2016 Mol Pharmacol), where four SLCO1A2 variants "
                     "showed impaired transport.",
            uptake_ratio=2.0,
            ratio_assay="~2-fold cellular uptake in OATP1A2-transfected HEK293 over "
                        "vector control — the canonical transfected-versus-control "
                        "observable, in a single-transporter background",
            ratio_citation="PMID 24554572 (Durmus 2014 Int J Cancer)",
            notes="A QUIETLY IMPORTANT ROW. Doxorubicin was the anchor case the "
                  "efflux table could NOT put a number on — the textbook P-gp "
                  "substrate and the molecule this repo's perivascular constant "
                  "alpha was fitted to, left at exactly 1.0 because no defensible "
                  "efflux ratio existed. Here it is one of the few drugs that DOES "
                  "carry a usable published uptake ratio. Doxorubicin is not on "
                  "either scorecard, so this changes no score; it is recorded "
                  "because it shows the coverage gap is drug-specific and "
                  "mechanism-specific rather than uniform.",
        ),
        UptakeEntry(
            transporter="OATP1B1", status="unknown",
            notes="GENUINELY ASSAY-FORMAT DEPENDENT, so coded unknown. NEGATIVE in "
                  "vitro: both PMID 24554572 and PMID 27777271 found no doxorubicin "
                  "uptake by OATP1B1 in HEK293 / HeLa single-transporter systems. "
                  "POSITIVE in the SAME paper (PMID 27777271) in MDCKII TRANSWELL "
                  "assays: 'OATP1B1, OATP1B3, and OATP1A2, either alone or in "
                  "combination with MDR1, significantly transported doxorubicin'. "
                  "POSITIVE in vivo: liver-specific transgenic OATP1B1 partially "
                  "rescued hepatic uptake in Oatp1a/1b-/- mice (PMID 24554572), "
                  "where Oatp1a/1b deficiency lowered the liver-to-plasma ratio up "
                  "to 4.1-fold.",
        ),
        UptakeEntry(transporter="OATP1B3", status="unknown"),
        UptakeEntry(
            transporter="SLC22A16", status="yes",
            citation="PMID 15963465 (Okabe 2005 BBRC): saturable, dose-dependent "
                     "doxorubicin transport by SLC22A16 (OCT6/CT2) in Xenopus "
                     "oocytes, apparent Km 5.2 +/- 0.4 uM; Jurkat transfectants "
                     "significantly more doxorubicin-sensitive at 2 uM.",
            notes="Contributes nothing because SLC22A16 is encoded with no tissues "
                  "(testis-specific, PMID 12089149) — see that transporter's note. "
                  "The Km is real and is the cleanest doxorubicin uptake kinetic "
                  "anywhere, which is exactly why it is recorded and neutralised "
                  "rather than silently omitted.",
        ),
        UptakeEntry(
            transporter="OCT1", status="unknown",
            notes="PubMed returns ZERO records for doxorubicin uptake by SLC22A1 or "
                  "SLC22A3. This is a true null, not a negative.",
        ),
    ),
    notes="The mirror image of its efflux row. Note the in vivo mouse numbers are "
          "attributable to Oatp1a/1b as a family, not to any single human "
          "orthologue: Oatp1a/1b-/- gave plasma AUC 1.3-fold higher and liver AUC "
          "1.4-fold lower with a liver-to-plasma ratio up to 4.1-fold lower (PMID "
          "24554572), and Slco1a/1b-/- gave plasma AUC 1.7-fold higher with a "
          "2.3-fold lower liver-to-plasma ratio and 40% reduced clearance (PMID "
          "27777271). Those are mouse-family numbers and are NOT used as a human "
          "single-transporter uptake ratio.",
))

_add(UptakeRecord(
    drug="cisplatin",
    entries=(
        UptakeEntry(
            transporter="CTR1", status="yes",
            citation="PMID 12370430 (Ishida 2002 PNAS): yeast CTR1 deletion confers "
                     "resistance with reduced intracellular cisplatin; mouse cells "
                     "lacking one or two mCtr1 alleles show resistance and decreased "
                     "accumulation in parallel with gene dosage.",
            notes="ACCUMULATION AND RESISTANCE, NOT A DIRECT FLUX MEASUREMENT, and "
                  "no Km/Vmax exists. Combined with the unverified tissue "
                  "distribution (see the CTR1 transporter note) this row "
                  "contributes exactly nothing, which is the honest outcome for a "
                  "mechanism everybody cites and nobody has localised in a "
                  "retrievable human paper.",
        ),
        UptakeEntry(
            transporter="OCT2", status="yes",
            citation="PMID 16914559 (Yonezawa 2006 JPET): HEK293 transient "
                     "transfection — both cytotoxicity and CELLULAR ACCUMULATION of "
                     "cisplatin enhanced by hOCT2; OCTN1 and OCTN2 explicitly "
                     "negative, and carboplatin and nedaplatin were not substrates "
                     "of any of the seven carriers tested (a useful internal "
                     "control).",
            notes="Renal, so it contributes nothing in brain or liver. In vivo "
                  "corroboration: OCT1/2 double-knockout mice show no ototoxicity "
                  "and only mild nephrotoxicity, and cimetidine protects wild-type "
                  "mice (PMID 20110413). Pharmacogenetics (association only, not "
                  "transport): SLC22A2 rs316019 protective against ototoxicity.",
        ),
        UptakeEntry(
            transporter="OCT1", status="unknown",
            notes="CONTRADICTORY and left unresolved. PMID 16914559 reports "
                  "accumulation 'weakly' enhanced by hOCT1, i.e. a weak substrate; "
                  "PMID 16314463 (Ciarimboli 2005) found cisplatin at 100 uM "
                  "inhibited transport via hOCT2-HEK293 but NOT hOCT1-HEK293, and "
                  "no effect in human hepatocyte couplets. Note also that "
                  "Ciarimboli measured inhibition of a fluorescent cation probe "
                  "plus apoptosis rather than cisplatin flux, so its claims are "
                  "inferential in both directions.",
        ),
    ),
    notes="Expectation verified rather than assumed, as in the efflux table, and "
          "with the same conclusion from the other side: platinum handling runs "
          "through copper transporters and glutathione, not through the OATP/OCT "
          "carriers this module models. The net effect on both scorecards is "
          "exactly 1.0, and cisplatin is on neither.",
))

_add(UptakeRecord(
    drug="etoposide",
    entries=(
        UptakeEntry(
            transporter="OATP1B1", status="unknown",
            notes="DISPUTED WITHIN A SINGLE PAPER, AND THE ABSTRACT IS THE PART "
                  "THAT IS WRONG. PMID 26663398 (Drenberg 2016 Clin Pharmacol Ther) "
                  "abstract asserts 'cytarabine, daunorubicin, etoposide, and "
                  "mitoxantrone are substrates for OATP1B1 and the mouse ortholog "
                  "Oatp1b2', but its own Results state verbatim that intracellular "
                  "uptake was significantly enhanced over vector control 'with the "
                  "exception of VP16' — VP16 being etoposide. Read the number, not "
                  "the claim. The paper's Discussion then cites OTHERS for "
                  "etoposide's OATP1B1 substrate status, not its own null arm.",
        ),
        UptakeEntry(
            transporter="OATP1B3", status="unknown",
            notes="Tested only as an INHIBITOR: Ki 13.5 uM against "
                  "OATP1B3-mediated [3H]CCK-8 (PMID 25618019). No transport of "
                  "etoposide was demonstrated, so this must not be upgraded to "
                  "substrate.",
        ),
    ),
    notes="A SPECIES SPLIT, recorded because it is exactly the kind of thing that "
          "gets flattened into a wrong 'yes'. In PMID 26663398 the MOUSE ortholog "
          "Oatp1b2 did significantly transport etoposide while human OATP1B1 did "
          "not reach significance in the same experiment; in vivo, Oatp1b2-/- mice "
          "showed 1.2-2.4-fold reduced liver accumulation and 1.3-1.5-fold higher "
          "plasma AUC. Monkey OATP1B1/1B3/2B1 also transport it (PMID 30652318). "
          "The human single-transporter evidence is the null arm. Etoposide is on "
          "neither scorecard.",
))

_add(UptakeRecord(
    drug="vincristine",
    entries=(
        UptakeEntry(
            transporter="OATP1B1", status="yes",
            citation="PMID 28828541 (Lancaster 2017 Pharm Res): "
                     "OATP1B1-transfected CHO cells showed significantly increased "
                     "vincristine uptake versus wild-type, reduced by OATP "
                     "inhibitors.",
            notes="REAL BUT QUANTITATIVELY MINOR, and the source itself says so: in "
                  "suspended human hepatocytes 'only 10% of the total vincristine "
                  "uptake process could be attributed to OATP-mediated transport'. "
                  "So even if a fold existed, the correct factor would be near 1.1, "
                  "not near the value a transfected CHO line would suggest. No "
                  "fold, no Km, no Vmax is published; single 2017 paper; "
                  "heterologous non-hepatic background.",
        ),
        UptakeEntry(
            transporter="OATP1B3", status="yes",
            citation="PMID 28828541 (Lancaster 2017 Pharm Res): same experiment, "
                     "OATP1B3-transfected CHO cells.",
            notes="Also an INHIBITOR of OATP1B3, Ki 30.6 uM (PMID 25618019) — not "
                  "substrate evidence. And beware the vinca substitution trap that "
                  "the efflux table also flags: much of the canonical literature is "
                  "VINBLASTINE (Ki 10.2 uM at OATP1B1, 15.9 uM at OATP1B3), which "
                  "is not interchangeable with vincristine for a quantitative term.",
        ),
    ),
    notes="Hepatic disposition of vincristine is dominated by a "
          "benzbromarone-sensitive EFFLUX mechanism (possibly MRP3) and by "
          "MDR1/P-gp, not by uptake. 'Vincristine is an OATP1B1 substrate' is true "
          "and should not carry a large uptake term in any PBPK model. Factor 1.0 "
          "for want of a number, which here happens to also be roughly the right "
          "answer.",
))

# --- macromolecule -----------------------------------------------------------
_add(UptakeRecord(
    drug="bevacizumab",
    entries=(),
    notes="NOT APPLICABLE rather than measured-negative, and the module "
          "short-circuits antibodies to 1.0 before this row is read. Solute "
          "carriers translocate small molecules through a substrate-binding pore; a "
          "~149 kDa IgG1 cannot be an SLC substrate. That is mechanistic REASONING, "
          "not an experimental negative — nobody has run a bevacizumab OATP1B1 "
          "assay and nobody should. The relevant transport route is FcRn (FCGRT), a "
          "RECEPTOR and not an SLC: IgG catabolism is 10x faster in FcRn-deficient "
          "mice (PMID 8643606), it mediates bidirectional transcytosis across human "
          "epithelium (PMID 10510331), and it is present in endothelium, "
          "macrophages, enterocytes and choroid plexus epithelium (PMID 28402755, "
          "PMID 12067234). IMPORTANT NUANCE THAT CUTS AGAINST THE USUAL NARRATIVE: "
          "PMID 19636712 measured a blood-corrected brain:plasma AUC ratio of "
          "0.0022 for IgG1 in wild-type mice and 0.0021 in FcRn-KO (P = 0.33), "
          "concluding FcRn is NOT responsible for poor brain exposure — so FcRn "
          "should not be cited as the cause of an antibody's BBB exclusion. This "
          "repo already handles bevacizumab in a size/permeability term "
          "(bbb_penetrance = 0 via the is_antibody and MW gates), which is the "
          "right place for it.",
))

# --- the glioma-scorecard drug the efflux run was built for ------------------
_add(UptakeRecord(
    drug="eribulin",
    entries=(
        UptakeEntry(transporter="OATP1B1", status="unknown"),
        UptakeEntry(transporter="OATP1B3", status="unknown"),
        UptakeEntry(transporter="OCT1", status="unknown"),
        UptakeEntry(transporter="OAT1", status="unknown"),
        UptakeEntry(transporter="OAT3", status="unknown"),
    ),
    notes="CODED UNKNOWN, NOT 'NO', AND THE PRECEDENT IS THIS REPO'S OWN. The "
          "HALAVEN US Prescribing Information section 12.3 states that eribulin at "
          "clinically relevant concentrations 'is not a substrate of ... OATP1B1, "
          "OATP1B3 ... OCT1, OCT2 ... OAT1, OAT3 ... MATE1 ... or BCRP', and that "
          "it IS a substrate of P-gp. That is a label claim with NO PMID and an "
          "unpublished underlying study — exactly the situation in which "
          "results/efflux_impact.md coded eribulin's BCRP status 'unknown' rather "
          "than 'no'. The same standard is applied here rather than a more "
          "convenient one. Note also the label's qualifier 'at clinically relevant "
          "concentrations', which leaves higher concentrations untested. OATP1A2, "
          "LAT1, MCT, ENT, PEPT and OCTN have never been tested. Either way the "
          "factor is 1.0, so eribulin's 6.94x over-prediction — the largest single "
          "error on the glioma scorecard — is untouched by this mechanism.",
))

_add(UptakeRecord(
    drug="cefazolin",
    entries=(
        UptakeEntry(
            transporter="OAT3", status="yes",
            citation="PMID 14984259 (Sakurai 2004 Pharm Res): 'the uptake study "
                     "using HEK293 transfectants revealed that cefazolin and "
                     "phenolsulfonphthalein were transported by hOAT3'. Corroborated "
                     "by PMID 16098483 (Ueo 2005 Biochem Pharmacol): accumulation "
                     "significantly greater in HEK-hOAT3 than control.",
            notes="No cefazolin-specific fold could be verified (paywalled); "
                  "companion cephalosporins in the same panel ran 35-50x over "
                  "control, which is suggestive but is a DIFFERENT MOLECULE'S "
                  "number and is not imported. No Km for cefazolin at hOAT3 exists. "
                  "Renal anyway, so it contributes nothing in brain.",
        ),
        UptakeEntry(
            transporter="OAT1", status="no",
            citation="PMID 16098483 (Ueo 2005 Biochem Pharmacol): uptake of these "
                     "antibiotics by HEK-hOAT1 was 'within two-fold of that by "
                     "control cells'; the open-access review PMC3153278 (PMID "
                     "22069672) reads the primary data and states hOAT1-dependent "
                     "uptake was non-significant for cefazolin and that 'cefazolin "
                     "is also not an hOAT1 substrate'.",
            notes="THIS CORRECTS results/efflux_impact.md, WHICH RECORDED CEFAZOLIN "
                  "AS 'a confirmed OAT1 substrate (PMID 15618660, Ki 72 uM vs PAH)'. "
                  "Two errors are compounded in that note and both are fixed here. "
                  "(1) PMID 15618660 is RAT Oat1, and the human orthologue does not "
                  "transport cefazolin — a species split that a 2025 review (PMID "
                  "40563504) already mis-propagates into a human OAT1 substrate "
                  "table. (2) A Ki is INHIBITION of PAH transport, not substrate "
                  "status. The rat result stands on its own terms (rOAT1 "
                  "significantly stimulated cefazolin uptake, IC50 72 +/- 8 uM vs "
                  "PAH) and is recorded as rat.",
        ),
        UptakeEntry(
            transporter="OATP1B1", status="yes",
            citation="PMID 18974612 (Nakakariya 2008 Drug Metab Pharmacokinet): "
                     "OATP1B1-mediated cefazolin transport, Km 20.8 +/- 3.6 mM, "
                     "Vmax 31.7 pmol/hr/oocyte.",
            notes="MILLIMOLAR, i.e. two orders of magnitude above clinical cefazolin "
                  "plasma concentrations (~100-200 uM total), and the value is "
                  "EXTRAPOLATED BEYOND THE TESTED RANGE (0.1-10 mM). Recorded for "
                  "completeness; almost certainly not clinically meaningful. Factor "
                  "1.0 regardless, since a Km is not an uptake ratio.",
        ),
        UptakeEntry(
            transporter="OATP1B3", status="yes",
            citation="PMID 18974612 (Nakakariya 2008 Drug Metab Pharmacokinet): "
                     "OATP1B3-mediated cefazolin transport, Km 3.89 +/- 2.21 mM, "
                     "Vmax 10.5 pmol/hr/oocyte.",
            notes="Millimolar, as OATP1B1. Same caveat.",
        ),
    ),
    notes="PEPT1 is deliberately coded nowhere and the reason is a good one: "
          "cefazolin lacks the alpha-amino group that PEPT1 substrates carry, and "
          "this is verified rather than assumed — GlySar uptake was 'unaffected by "
          "cefazolin' (PMID 22490229), with no intracellular pH drop in Caco-2 "
          "(PMID 15618692) and a negative rabbit in vivo result (PMID 11765153), "
          "while cefadroxil and cephradine (which DO have the alpha-amino group) "
          "are substrates. PEPT2 has never been directly tested. HUMAN IN VIVO "
          "CAVEAT: probenecid prolongs cefazolin serum half-life by 20-40% (PMID "
          "4758839, n = 12, replicated), but probenecid is non-selective across "
          "OAT1, OAT3, MRP2 and MRP4, so the interaction does not localise to OAT3 "
          "— and PMID 4758839 itself concludes cefazolin is 'cleared primarily by "
          "the glomerulus', a flagged contradiction with the modern "
          "OAT3-secretion picture. An expression-phenotype correlation exists "
          "(renal hOAT3 mRNA vs cefazolin elimination rate constant, r = 0.536, "
          "PMID 14984259) but no genotype-PK association does. Everything relevant "
          "here is renal; in glioma cefazolin's factor is 1.0.",
))

_add(UptakeRecord(
    drug="levetiracetam",
    entries=(
        UptakeEntry(
            transporter="OCT1", status="no",
            citation="PMID 35955563 (Sachkova 2022 Int J Mol Sci): OCT1-HEK293 vs "
                     "mock uptake ratio = 1.00 against the authors' own substrate "
                     "threshold of > 1.5.",
            notes="A clean measured negative in a properly controlled screen.",
        ),
    ),
    notes="THE ORPHAN-CARRIER HYPOTHESIS WAS TESTED HEAD-ON AND FAILED, which is a "
          "stronger result than 'unknown'. The proton-coupled organic-cation "
          "antiporter is the long-postulated unidentified influx carrier at the "
          "BBB; PMID 35955563 screened 72 drugs in hCMEC/D3 and levetiracetam gave "
          "antiport ratios of 1.17 +/- 0.28 and 1.19 +/- 0.24 against the authors' "
          "criterion that ratios near 1.0 indicate NON-substrates. No saturable "
          "brain uptake exists in any species: serum and CSF rise linearly with "
          "dose in rat (PMID 10210031) and brain-ECF entry is linear with no "
          "regional heterogeneity (PMID 11454660). CORRECTION TO A COMMON PREMISE, "
          "including one this repo could easily have inherited: levetiracetam's "
          "passive permeability is MODEST, not high — PS = 0.015 mL/min/g, 21x "
          "lower than brivaracetam (PMID 26663401) — so its slow brain entry is low "
          "passive permeability, not a missing carrier. Two traps avoided: the only "
          "Michaelis-Menten constant in the levetiracetam literature (Km 435 uM, "
          "PMID 14530892) belongs to a blood ESTERASE, not a transporter; and the "
          "metabolite ucb L057 is actively secreted and probenecid-sensitive while "
          "the parent is not, so 'levetiracetam is actively secreted' conflates the "
          "two. Every other carrier is UNKNOWN (never tested), which is a genuine "
          "null rather than a negative.",
))

# --- the PARP inhibitors (Moncayo) ------------------------------------------
_add(UptakeRecord(
    drug="niraparib",
    entries=(
        UptakeEntry(
            transporter="OATP1B1", status="no",
            citation="PMID 34730366 (Martins/Loos/Schinkel 2021 Mol Pharm): "
                     "Oatp1a/1b deletion did not significantly affect niraparib "
                     "oral bioavailability or LIVER DISTRIBUTION in vivo (it did "
                     "reduce liver uptake of the M1 metabolite). Concordant with "
                     "the ZEJULA US PI section 12.3 and the EMA EPAR "
                     "(EMA/648982/2017), both of which state niraparib is not a "
                     "substrate of OATP1B1, OATP1B3 or OCT1 — label claims with no "
                     "PMID, cited here only as corroboration of the peer-reviewed "
                     "in vivo negative.",
            notes="An in vivo hepatic negative in the transporter-knockout design "
                  "is a stronger negative than a monolayer screen, and it is "
                  "peer-reviewed, which is why this row is 'no' where eribulin's is "
                  "'unknown'.",
        ),
        UptakeEntry(transporter="OATP1B3", status="unknown"),
        UptakeEntry(transporter="OCT1", status="unknown"),
    ),
    notes="Distinguish substrate from inhibitor, which matters for this drug: "
          "niraparib is a potent MATE1/MATE2-K INHIBITOR (IC50 0.18 uM / <= 0.14 "
          "uM, ZEJULA label) and that is not transport. UNRESOLVED CONTRADICTION "
          "recorded rather than smoothed: the FDA label says niraparib 'did not "
          "inhibit ... OCT1' while the EMA EPAR says it 'weakly inhibits OCT1 with "
          "an IC50 = 34.4 uM' and advises caution with metformin. The measured "
          "6-7x brain:plasma rise in Abcb1a/1b-KO mice (PMID 34730366) is LOSS OF "
          "EFFLUX, not gain of uptake, and belongs to the efflux module. No "
          "SLCO1B1/SLC22A1 genotype-exposure study exists for any PARP inhibitor.",
))

_add(UptakeRecord(
    drug="rucaparib",
    entries=(
        UptakeEntry(
            transporter="OATP1B1", status="no",
            citation="PMID 32129697 (Liao 2020 Xenobiotica): 'rucaparib was a "
                     "substrate for P-gp and BCRP, but not for OATP1B1, OATP1B3, "
                     "OAT1, OAT3, or OCT2.' Peer-reviewed, unlike the label-only "
                     "negatives elsewhere in this table.",
        ),
        UptakeEntry(
            transporter="OATP1B3", status="no",
            citation="PMID 32129697 (Liao 2020 Xenobiotica): same experiment.",
        ),
        UptakeEntry(
            transporter="OAT1", status="no",
            citation="PMID 32129697 (Liao 2020 Xenobiotica): same experiment.",
        ),
        UptakeEntry(
            transporter="OAT3", status="no",
            citation="PMID 32129697 (Liao 2020 Xenobiotica): same experiment.",
        ),
        UptakeEntry(
            transporter="OCT1", status="unknown",
            notes="TESTED ONLY AS AN INHIBITOR, NEVER AS A SUBSTRATE. Rucaparib "
                  "inhibits OCT1 with IC50 4.3 uM (PMID 32129697), which is not "
                  "evidence either way about transport.",
        ),
    ),
    notes="THE DRUG MOST LIKELY TO BE MISCODED IN THIS WHOLE TABLE, because it is a "
          "potent inhibitor of the very carriers it is not a substrate of: MATE2-K "
          "IC50 0.19 uM, MATE1 0.63 uM, OCT1 4.3 uM, OCT2 31 uM (PMID 32129697). "
          "UNRESOLVED REGULATORY CONTRADICTION recorded: the RUBRACA US PI states "
          "flatly that 'rucaparib inhibited OATP1B1, OATP1B3, OAT1, OAT3 ...' while "
          "the EMA EPAR states that at clinical exposures it did NOT inhibit them; "
          "PMID 32129697 explains the divergence (only 42-66% inhibition at 300 uM, "
          "far above clinical exposure) but the US PI omits the concentration "
          "qualifier. The clinical creatinine rise is attributed to MATE1/MATE2-K "
          "INHIBITION blocking renal creatinine secretion — again inhibition, not "
          "substrate status.",
))

_add(UptakeRecord(
    drug="olaparib",
    entries=(
        UptakeEntry(
            transporter="OATP1B1", status="unknown",
            notes="THE SINGLE MOST CONSEQUENTIAL UNRESOLVED ITEM IN THE PARP SET, "
                  "and the contradiction is INTERNAL TO ONE DOCUMENT — the EMA "
                  "Lynparza EPAR (EMA/CHMP/789139/2014) section 2.3.3. Study "
                  "13ASTRUKP7S2, in HEK293 overexpressing OATP1B1/1B3/OCT1 at 1-100 "
                  "uM with working positive controls, found 'no evidence that "
                  "olaparib was a substrate for any of OATP1B1, OATP1B3 or OCT1'. "
                  "Study KMX042 in the same document suggested possible OATP1B1 "
                  "involvement in hepatic uptake, but the assessor records that it "
                  "'was assessed qualitatively and not quantitatively' and that "
                  "'no conclusion can be drawn due to the variability of the test "
                  "system'. Coding this 'no' would over-state the evidence.",
        ),
        UptakeEntry(transporter="OATP1B3", status="unknown"),
        UptakeEntry(transporter="OCT1", status="unknown"),
    ),
    notes="You cannot read uptake-substrate status off the FDA label for this drug: "
          "the LYNPARZA US PI section 12.3 lists olaparib only as an INHIBITOR of "
          "BCRP, OATP1B1, OCT1, OCT2, OAT3, MATE1 and MATE2-K, and as a "
          "substrate+inhibitor of P-gp, without stating uptake-substrate status at "
          "all. PMID 27684210 (McCormick & Swaisland 2016 Xenobiotica) says only "
          "that olaparib 'was found to be a substrate of MDR1 but not of several "
          "other transporters' without enumerating them. Inhibitor IC50s, for "
          "completeness and NOT as substrate evidence: OATP1B1 20.3-27.1 uM, OCT1 "
          "37.9 uM, OCT2 19.9 uM, OAT3 18.4 uM, MATE1 5.50 uM.",
))



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


#: Drugs in the table that are macromolecules, for which a solute carrier is
#: NOT APPLICABLE rather than unmeasured. Needed because a caller passing a bare
#: name string carries no ``is_antibody`` flag, and reporting bevacizumab as
#: ``"unknown_substrate_status"`` would wrongly imply nobody had looked.
ANTIBODY_DRUGS = frozenset({"bevacizumab"})


def _is_antibody(drug: Union[str, Mapping]) -> bool:
    if isinstance(drug, Mapping):
        if drug.get("is_antibody", False):
            return True
        name = drug.get("name")
        return str(name).strip().lower() in ANTIBODY_DRUGS if name else False
    if getattr(drug, "is_antibody", False):
        return True
    return (_resolve_name(drug) or "") in ANTIBODY_DRUGS


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


# ===========================================================================
# 5. Self-test
# ===========================================================================
if __name__ == "__main__":
    # --- the functional form is what the docstring derives, and is UNIQUE ---
    # Three-compartment polarised monolayer with an uptake carrier on the apical
    # face, simulated directly from the mass balance. All three published
    # observables must collapse to the SAME number -- this is the property that
    # spares this module the two-forms problem the efflux module had.
    for P in (0.4, 1.0, 3.0):
        for U in (0.0, 1.0, 5.0, 50.0):
            c_ab = (P + U) / (2.0 * P)          # unit apical donor, sink basal
            papp_ab = P * c_ab
            c_ba = P / (2.0 * P)                # unit basal donor, sink apical
            papp_ba = P * c_ba
            ur = 1.0 + U / P
            assert abs(papp_ab / papp_ba - ur) < 1e-12, "directional ratio"
            assert abs(papp_ab / (P / 2.0) - ur) < 1e-12, "vs passive-only monolayer"
            assert abs((P + U) / P - ur) < 1e-12, "cell accumulation fold"

    # --- parallel-carrier combination --------------------------------------
    assert combine_uptake_ratios([]) is None
    assert combine_uptake_ratios([None, None]) is None
    assert combine_uptake_ratios([5.0, None]) == 5.0
    assert combine_uptake_ratios([3.0, 5.0]) == 7.0          # 3 + 5 - 1
    assert combine_uptake_ratios([1.0, 1.0]) == 1.0          # two dead carriers
    assert combine_uptake_ratios([0.4]) == 1.0               # never reduces

    # --- the factor never drops below 1.0, in any tissue, for any drug -----
    for name in list(UPTAKE_TABLE) + ["no_such_drug"]:
        for t in TISSUES:
            assert uptake_factor(name, tissue=t) >= 1.0, (name, t)

    # --- everything unproven is exactly 1.0, with the RIGHT reason ---------
    for probe, tissue, expect in (
            ("no_such_drug", TISSUE_GENERIC, "not_in_table"),
            ("procarbazine", TISSUE_NORMAL_LIVER, "unknown_substrate_status"),
            ("lomustine", TISSUE_NORMAL_LIVER, "non_substrate"),
            ("levetiracetam", TISSUE_NORMAL_LIVER, "non_substrate"),
            ("methotrexate", TISSUE_GLIOMA, "substrate_but_no_published_ratio"),
            ("adavosertib", TISSUE_GLIOMA, "only_circular_observable"),
            # THE tissue-differential branch: a real substrate of a real carrier,
            # in a tissue that does not express it.
            ("imatinib", TISSUE_LIVER_METASTASIS, "no_transporter_in_this_tissue"),
            ("methotrexate", TISSUE_NORMAL_LIVER, "substrate_but_no_published_ratio"),
            # doxorubicin IS a substrate (OATP1A2, SLC22A16) but of carriers no
            # hepatocyte expresses, so it too lands in the differential branch
            ("doxorubicin", TISSUE_NORMAL_LIVER, "no_transporter_in_this_tissue"),
    ):
        d = uptake_detail(probe, tissue=tissue)
        assert d["factor"] == 1.0, (probe, tissue, d["factor"])
        assert d["reason"] == expect, (probe, tissue, d["reason"], expect)
    assert uptake_factor({"name": "bevacizumab", "is_antibody": True}) == 1.0
    assert uptake_detail({"name": "x", "is_antibody": True})["reason"] == "antibody"
    # a BARE NAME must also resolve to "antibody", not to "unknown": for a
    # macromolecule an SLC is not applicable, which is not the same as unmeasured
    assert uptake_detail("bevacizumab")["reason"] == "antibody"

    # --- the drugs that DO carry a number ----------------------------------
    assert abs(uptake_factor("temozolomide", tissue=TISSUE_GLIOMA) - 1.70) < 1e-12
    assert abs(uptake_factor("doxorubicin", tissue=TISSUE_GLIOMA) - 2.00) < 1e-12
    assert abs(uptake_factor("imatinib", tissue=TISSUE_GLIOMA) - 1.50) < 1e-12
    # ... and are inert in tissues lacking the carrier
    assert uptake_factor("temozolomide", tissue=TISSUE_NORMAL_LIVER) == 1.0
    assert uptake_factor("doxorubicin", tissue=TISSUE_LIVER_METASTASIS) == 1.0

    # --- the circularity guard actually guards -----------------------------
    assert uptake_factor("adavosertib", tissue=TISSUE_GLIOMA) == 1.0
    assert abs(uptake_factor("adavosertib", tissue=TISSUE_GLIOMA,
                             allow_circular=True) - 3.2) < 1e-12
    assert UPTAKE_TABLE["adavosertib"].entries[0].suppressed_circular

    # --- every asserted fact carries a citation ----------------------------
    for rec in UPTAKE_TABLE.values():
        assert rec.notes.strip(), rec.drug
        for e in rec.entries:
            if e.status != "unknown":
                assert e.citation.strip(), (rec.drug, e.transporter)
            if e.uptake_ratio is not None:
                assert "PMID" in e.ratio_citation, (rec.drug, e.transporter)
            if e.kpuu_in_vivo is not None:
                assert "PMID" in e.kpuu_citation, (rec.drug, e.transporter)
    for t in TRANSPORTERS.values():
        assert t.notes.strip(), t.name
        if t.tissues:
            assert "PMID" in t.tissue_citation, t.name

    # --- composition reduces EXACTLY when the factor is 1.0 ----------------
    from ..reconstruct.barrier import L_PARACELLULAR_LEAK, effective_penetrance
    for name, tissue in (("cisplatin", TISSUE_GLIOMA),
                         ("no_such_drug", TISSUE_GENERIC),
                         ("imatinib", TISSUE_LIVER_METASTASIS)):
        assert uptake_factor(name, tissue=tissue) == 1.0
        for p in (0.05, 0.5, 0.9, 1.0):
            for d in (0.0, 0.5, 1.0):
                assert (uptake_adjusted_penetrance(p, name, tissue=tissue,
                                                   bbb_disruption=d)
                        == float(effective_penetrance(p, d))), (name, p, d)
    # uptake hits the transcellular route only: the paracellular leak survives
    hit = uptake_adjusted_penetrance(0.2, "temozolomide", tissue=TISSUE_GLIOMA,
                                     bbb_disruption=1.0)
    assert abs(hit - (0.2 * 1.70 + L_PARACELLULAR_LEAK)) < 1e-12
    # and the free-drug ceiling rises to U, not to 1.0
    assert abs(uptake_adjusted_penetrance(1.0, "doxorubicin", tissue=TISSUE_GLIOMA)
               - 2.0) < 1e-12

    # --- THE STRUCTURAL POINT: this term does NOT cancel out of a ratio ----
    # The efflux module's self-test proves the opposite property for a scalar.
    e_t, e_n = 0.5572, 0.9607                    # GIST Krogh delivery efficiencies
    u_t = uptake_factor("imatinib", tissue=TISSUE_LIVER_METASTASIS)
    u_n = uptake_factor("imatinib", tissue=TISSUE_NORMAL_LIVER)
    scalar = uptake_factor("imatinib", tissue=TISSUE_GLIOMA)     # a per-drug value
    assert abs((scalar * e_t) / (scalar * e_n) - e_t / e_n) < 1e-15, \
        "a SCALAR must cancel -- that was efflux's problem"
    assert (u_t * e_t) / (u_n * e_n) == (u_t / u_n) * (e_t / e_n), \
        "a tissue-differential factor does not cancel; it scales the ratio by u_t/u_n"

    # --- report -------------------------------------------------------------
    s = substrate_summary()
    print(f"UPTAKE_TABLE: {s['n_drugs']} drugs, {len(TRANSPORTERS)} transporters")
    print(f"  substrate of >= 1 carrier          : {s['substrate_any']}")
    print(f"  .. WITH a usable published number  : {s['substrate_with_number']}")
    print(f"  .. WITHOUT one                     : {s['substrate_without_number']}")
    print(f"  .. only number withheld as circular: {s['circular_suppressed']}")
    print(f"  measured non-substrate             : {s['non_substrate']}")
    print(f"  unknown throughout                 : {s['unknown_only']}")
    print()
    tis = [t for t in TISSUES if t != TISSUE_GENERIC]
    hdr = f"{'drug':14s} " + " ".join(f"{t[:9]:>9s}" for t in tis) + "   reason(glioma)"
    print(hdr)
    print("-" * len(hdr))
    for name in UPTAKE_TABLE:
        cells = " ".join(f"{uptake_factor(name, tissue=t):9.2f}" for t in tis)
        print(f"{name:14s} {cells}   {uptake_detail(name, tissue=TISSUE_GLIOMA)['reason']}")
    print()
    print("tissue -> carriers:")
    for t in TISSUES:
        print(f"  {t:20s} {', '.join(transporters_in_tissue(t)) or '(none)'}")
