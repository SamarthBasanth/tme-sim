"""Active efflux (P-glycoprotein / ABCB1, BCRP / ABCG2) as a penetrance factor.

Why this module exists
----------------------
Two independent measured-data validations in this repo name **active efflux** as
a missing mechanism, and one of them names it as *the* largest residual on the
delivery axis:

* ``results/glioma_regional_pk.md`` — eribulin, a 730 Da macrocycle and a
  documented P-gp substrate, is over-predicted 6.94x even after the CNS-MPO veto
  fix. The report states plainly: *"Even a correct passive score would be too
  high, because the real limit on eribulin is active efflux, which the
  function's own docstring concedes it 'cannot capture'."*
* ``results/gist_imatinib.md`` — the transport model is *structurally* falsified
  for tissue partition, with **hepatic uptake transporters** and **lysosomal
  ion-trapping** named as the missing physics.

Of the mechanisms those two runs name, efflux is the one that is **tractable
today**. It is a well-characterised, per-drug, published property with a
standard assay and a standard reported number (the efflux ratio), and it acts at
the blood-brain barrier and the cell membrane in a way a first-order term can
legitimately represent. The other named mechanism — a *spatial* acidic-
compartment field for ion-trapping — is deliberately NOT implemented here: the
Moncayo explant run (``results/moncayo_explant.md`` §4.3) proved that a SCALAR
ion-trapping term mathematically cannot produce spatial dispersion (uniform
retardation gives ``C = R * C_bath`` everywhere at equilibrium), and two
literature passes established that no measured acidic-phospholipid content in
usable units exists for tumour tissue. See :ref:`the honesty section
<efflux-honesty>` below for what THIS module can and cannot do, which turns out
to be subject to a closely related structural limit.

The functional form
-------------------
The published number for a transporter substrate is the **efflux ratio**

.. math::

    \\mathrm{ER} = \\frac{P_\\mathrm{app}(B \\rightarrow A)}{P_\\mathrm{app}(A \\rightarrow B)}

measured in a polarised monolayer (MDCKII-ABCB1, MDCKII-ABCG2, Caco-2) with the
transporter on the apical (``A``) face. We need the factor by which that pump
reduces net penetration relative to the same molecule with the pump switched
off. Both candidate forms below parameterise the same underlying quantity, the
**efflux burden** ``b = CL_efflux / P_passive >= 0``, and both give ``ER = 1 + b``.
They differ only in how the barrier is idealised.

**1. Single-barrier / unbound-partition form —** :data:`FORM_KPUU` **, the
default.** Treat the barrier as one membrane with efflux clearance in parallel
with passive permeability, which is the standard derivation behind the published
definition ``ER = (P_pass + V_max/K_m) / P_pass`` (Desai et al. 2013, PMID
23421687). Then

.. math::

    f_\\mathrm{efflux} = K_{p,uu}
      = \\frac{P_\\mathrm{pass}}{P_\\mathrm{pass} + CL_\\mathrm{eff}}
      = \\frac{1}{1 + b} = \\frac{1}{\\mathrm{ER}}

**2. Two-membrane monolayer form —** :data:`FORM_MONOLAYER` **, the conservative
alternative.** Resolve the cell interior explicitly: apical ``A`` — cell ``C`` —
basolateral ``B``, passive permeability ``P`` across BOTH membranes, efflux
clearance ``E`` at the apical membrane only, sink receiver. Quasi-steady state on
the cell gives ``P_app(A->B) = P^2/(2P+E)`` and ``P_app(B->A) = P(P+E)/(2P+E)``,
so ``ER = 1 + E/P`` as before, while the passive-only two-membranes-in-series
permeability is ``P/2``, hence

.. math::

    f_\\mathrm{efflux} = \\frac{P_\\mathrm{app}(A\\!\\to\\!B)}{P/2}
                       = \\frac{2}{1 + \\mathrm{ER}}

which is exactly **2x more permissive** than form 1 at large ER.

**Why form 1 is the default, and how that was decided.** Both return exactly 1.0
at ``ER = 1``, so neither penalises a non-substrate — which rules out the plain
``1/(1+ER)`` starting form, since that would silently halve every molecule that
has no efflux at all. Between the two, three things favour ``1/ER``:

* It is the form the field actually uses and the one with published empirical
  support against measured ``Kp,uu,brain``: in vitro-to-in vivo r^2 = 0.813 for a
  high-expressing MDCK-MDR1 line (Kikuchi et al. 2013, PMID 24009309), good
  transwell-to-knockout agreement across 32 CNS compounds (Feng et al. 2008,
  PMID 17962372), and ``Kp,uu`` inversely tracking ER above ER ~ 3 (Jiang et al.
  2022, PMID 35156331). No comparable empirical support exists for ``2/(1+ER)``.
* It is what the published *definition* of ER already encodes; the ``+1`` is
  inside ER, so dividing by ER does not double-count.
* On the one drug in our own table carrying BOTH observables it fits better —
  see the consistency test below.

Form 2 is retained and reported as a sensitivity precisely because it is the
weaker penalty: if a conclusion flips between the two, that fact belongs in the
report rather than in a footnote.

**Rejected: an in vivo scaling factor.** The best-supported practical
formulation in the literature is not a closed form at all but ``Kp,uu =
1/(1 + SF*(ER-1))`` with an assay- and species-specific scaling factor ``SF``
calibrated to neuro-PK data (Trapa et al. 2019, PMID 30683809; transporter
abundance differences quantified by Uchida et al. 2011, PMID 21291474 — human
BBB P-gp is 2.3x LOWER than mouse mdr1a, human BCRP 1.9x higher). ``SF`` is
deliberately **not** implemented. There is exactly one place we could calibrate
it — the six-drug measured scorecard this term is about to be tested against —
and fitting it there would make the test circular. The shipped term therefore has
**zero free parameters** and is correspondingly allowed to fail.

The other published observable is the same model, not a second model
--------------------------------------------------------------------
Several drugs have no in vitro ER but DO have an in vivo transporter-knockout
brain-to-plasma fold change, ``fold_KO = Kp(knockout)/Kp(wild-type)``. That is
the same model read out a different way: deleting the pump sets ``b = 0``, so
``fold_KO = 1/f_efflux`` in either form, and therefore

.. math::

    f_\\mathrm{efflux} = \\frac{1}{\\mathrm{fold}_\\mathrm{KO}}
    \\qquad\\text{(identical in both forms)}

with the equivalent ER being ``ER = fold_KO`` under :data:`FORM_KPUU` and
``ER = 2*fold_KO - 1`` under :data:`FORM_MONOLAYER`. (Both relations are verified
numerically in this module's self-test.) So the module is not form-shopping when
it uses whichever observable a given drug happens to have — there is one model
and two windows onto it, and the in vivo window gives the same answer either way.

It also gives a real **consistency test** wherever both windows exist. Under
:data:`FORM_KPUU` the in vitro ER predicts ``fold_KO = ER``; under
:data:`FORM_MONOLAYER` it predicts ``fold_KO = (1+ER)/2``. Exactly two drugs in
:data:`EFFLUX_TABLE` carry both:

=============  ==============  ===========  ==========================  ==========
drug           in vitro ER     measured     FORM_KPUU discrepancy       FORM_MONO
                               fold_KO      (pred = ER)                 (pred = (1+ER)/2)
=============  ==============  ===========  ==========================  ==========
adavosertib    14.1            25           1.77x (in vivo stronger)    3.31x
imatinib       39              28           0.72x (in vivo weaker)      1.40x
=============  ==============  ===========  ==========================  ==========

Geometric mean absolute discrepancy: **1.57x for FORM_KPUU, 2.15x for
FORM_MONOLAYER**. So the default form is the better-calibrated of the two on our
own rows as well as in the published correlations. **n = 2** — this is a
tie-breaker between two pre-specified forms, not a calibration, and it is
reported rather than hidden.

**Which observable is preferred, decided in advance.** In vitro ER is the primary
source when it exists: it is per-pump, it is the standard reported quantity, and
it is measured on human transporters. The in vivo knockout fold is the fallback,
used only where no in vitro ER is published. That rule is fixed here so it cannot
be chosen per drug to flatter a result; the opposite preference is run as a
sensitivity in ``results/efflux_impact.md``.

Combining two pumps
-------------------
P-gp and BCRP sit on the same apical membrane and act in parallel, so their
efflux clearances ADD: ``E_total = E_Pgp + E_BCRP``. Because ``E/P = ER - 1``,

.. math::

    \\mathrm{ER}_\\mathrm{combined} = \\mathrm{ER}_\\mathrm{Pgp}
                                    + \\mathrm{ER}_\\mathrm{BCRP} - 1

which is :func:`combine_efflux_ratios`. Note this is only valid when both ratios
were measured against the *same* passive permeability ``P`` — i.e. same
molecule, comparable monolayer. Ratios from an in vivo knockout study are NOT
combined this way and are flagged as such in the table.

.. _efflux-honesty:

Honesty / scope — read this before using it
--------------------------------------------
* **It is a scalar per-drug property.** :func:`efflux_factor` returns one number
  per molecule. That number therefore **cancels exactly** out of any
  within-drug, between-tissue RATIO — the GIST tumour:normal imatinib target is
  exactly such a ratio, so this term is *structurally incapable* of moving it,
  for the same algebraic reason the Moncayo run gave for scalar ion-trapping.
  This is demonstrated numerically in ``tmesim/validation/efflux_impact.py``,
  not merely asserted. Anything that WOULD move a tissue ratio requires a
  per-tissue transporter-expression field, which is a different (and much less
  well measured) object.
* **Substrate status is looked up, never inferred.** There is no structure-based
  P-gp-substrate classifier here, because a wrong classification is worse than
  no classification. A drug not in :data:`EFFLUX_TABLE`, or in it as
  ``"unknown"``, is treated as a **non-substrate (factor 1.0)** and the reason is
  reported by :func:`efflux_detail`.
* **A confirmed substrate with no published ER also gets factor 1.0.** We will
  not invent a ratio to fill a gap; :func:`efflux_detail` reports
  ``reason="substrate_but_no_published_er"`` so the omission is visible rather
  than silent.
* **ER is an in vitro number with real scatter against in vivo.** It is a
  monolayer measurement in an overexpressing cell line; the same drug can differ
  several-fold between labs and assay formats. Every value in the table carries
  its assay type so a user can see what they are trusting.
* **Efflux is transcellular.** Where the regional paracellular leak of
  :mod:`tmesim.reconstruct.barrier` is active, efflux must multiply the
  *transcellular* penetrance only — a broken tight junction bypasses the pumps.
  :func:`efflux_adjusted_penetrance` enforces this composition.
* **Nothing here is fitted.** Every ER in the table is a published measurement
  with a PMID/DOI. The module introduces **no free parameter**.

Architecture: this module reads nothing from the engine and writes nothing back.
``core``, ``physics`` and ``sim.engine`` are NOT modified.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Union

__all__ = [
    "EffluxRecord",
    "EFFLUX_TABLE",
    "FORM_MONOLAYER",
    "FORM_KPUU",
    "EFFLUX_FORMS",
    "PREFER_IN_VITRO",
    "PREFER_IN_VIVO",
    "SOURCE_PREFERENCES",
    "combine_efflux_ratios",
    "knockout_fold_predicted_from_er",
    "efflux_ratio_for",
    "efflux_factor",
    "efflux_detail",
    "efflux_adjusted_penetrance",
    "substrate_summary",
]


# ===========================================================================
# 2. The functional forms
# ===========================================================================
#: Single-barrier unbound-partition form, ``f = Kp,uu = 1 / ER``. **The default.**
#:
#: Efflux clearance in parallel with passive permeability across one barrier —
#: the derivation behind the published definition of ER itself (Desai et al.
#: 2013, PMID 23421687) and the form with actual empirical support against
#: measured brain ``Kp,uu`` (PMIDs 24009309, 17962372, 35156331).
#: ``f(1) = 1`` exactly.
FORM_KPUU = "kpuu"

#: Two-membrane monolayer form, ``f = 2 / (1 + ER)``. Conservative alternative.
#:
#: Resolves the cell interior explicitly (apical / cell / basolateral, efflux at
#: the apical membrane). ``f(1) = 1`` exactly and ``f -> 2/ER`` for large ER, so
#: it is exactly 2x MORE PERMISSIVE than :data:`FORM_KPUU`. Retained for
#: sensitivity analysis: if a conclusion depends on which of these two is used,
#: that is a fact about the conclusion.
FORM_MONOLAYER = "monolayer"

EFFLUX_FORMS = (FORM_KPUU, FORM_MONOLAYER)


def _apply_form(er: float, form: str) -> float:
    """Map an efflux ratio to a penetrance multiplier in ``(0, 1]``."""
    if form not in EFFLUX_FORMS:
        raise ValueError(f"unknown efflux form {form!r}; expected one of {EFFLUX_FORMS}")
    if er < 1.0:
        # An ER below 1 means net ABSORPTIVE transport (an uptake transporter, or
        # assay noise). This module models efflux only; it must never boost
        # penetrance, so the factor is clamped at 1.
        er = 1.0
    f = (1.0 / er) if form == FORM_KPUU else (2.0 / (1.0 + er))
    return float(min(1.0, max(0.0, f)))


def _er_from_knockout_fold(fold: float, form: str) -> float:
    """Express an in vivo knockout fold change on the ER scale OF THE GIVEN FORM.

    Both forms satisfy ``fold_KO = 1/f_efflux`` — deleting the pump sets the
    efflux burden to zero — so inverting each form's ``f(ER)`` gives

    * :data:`FORM_KPUU`      ``ER = fold``
    * :data:`FORM_MONOLAYER` ``ER = 2*fold - 1``

    and consequently ``efflux_factor`` computed from a knockout fold is
    ``1/fold`` in EITHER form. The form only matters for in vitro ratios.
    """
    if form not in EFFLUX_FORMS:
        raise ValueError(f"unknown efflux form {form!r}; expected one of {EFFLUX_FORMS}")
    fold = max(1.0, float(fold))
    return fold if form == FORM_KPUU else (2.0 * fold - 1.0)


def knockout_fold_predicted_from_er(er: float, form: str = FORM_KPUU) -> float:
    """The knockout brain-to-plasma fold change an in vitro ER implies.

    The inverse of :func:`_er_from_knockout_fold`: ``fold = ER`` under
    :data:`FORM_KPUU` and ``fold = (1 + ER)/2`` under :data:`FORM_MONOLAYER`.
    Comparing this to a MEASURED knockout fold is the module's one available
    internal check on the functional form.
    """
    return 1.0 / _apply_form(float(er), form)


def combine_efflux_ratios(er_pgp: Optional[float],
                          er_bcrp: Optional[float]) -> Optional[float]:
    """Combine P-gp and BCRP efflux ratios measured on the same molecule.

    The two pumps sit on the same apical membrane and act in PARALLEL, so their
    efflux clearances add: ``E_tot = E_Pgp + E_BCRP``. With ``E/P = ER - 1``
    this gives

    ::

        ER_combined = ER_Pgp + ER_BCRP - 1

    Passing ``None`` for either pump means "no usable number", not "ratio 1", so
    it is simply omitted from the sum. Returns ``None`` if both are ``None``.

    This addition is only legitimate when both ratios were measured against the
    same passive permeability, i.e. same molecule in comparable monolayers. It is
    NOT applied to in vivo knockout fold-changes, which already report a
    whole-barrier effect.
    """
    parts = [e for e in (er_pgp, er_bcrp) if e is not None]
    if not parts:
        return None
    return float(sum(max(1.0, p) for p in parts) - (len(parts) - 1))


# ===========================================================================
# 1. The record type
# ===========================================================================
@dataclass(frozen=True)
class EffluxRecord:
    """One drug's literature-sourced ABC-transporter substrate assignment.

    Every field that asserts a fact about the world carries its own citation.
    ``"unknown"`` is a first-class value and is treated downstream exactly like
    ``"no"`` (factor 1.0), but is reported differently so the distinction between
    *measured to be a non-substrate* and *never measured* is never lost.

    Attributes
    ----------
    drug : str
        Lower-case drug name, matching ``molecular.library.GBM_DRUGS`` where the
        drug is in the library.
    pgp : {"yes", "no", "unknown"}
        P-glycoprotein (ABCB1 / MDR1) substrate status.
    pgp_citation : str
        PMID/DOI + one-line justification for ``pgp``. Required whenever ``pgp``
        is not ``"unknown"``.
    bcrp : {"yes", "no", "unknown"}
        BCRP (ABCG2) substrate status.
    bcrp_citation : str
        PMID/DOI + one-line justification for ``bcrp``.
    er_pgp, er_bcrp : float or None
        Published **in vitro** bidirectional monolayer efflux ratios,
        ``Papp(B->A)/Papp(A->B)``. ``None`` means *no usable published number*,
        which is NOT the same as 1.0 and is handled separately.
    er_assay : str
        The monolayer system ``er_pgp`` / ``er_bcrp`` were measured in
        (MDCKII-ABCB1, MDCKII-ABCG2, Caco-2, ...).
    er_citation : str
        PMID/DOI for the in vitro ratio(s).
    fold_ko : float or None
        Published **in vivo** transporter-knockout brain-to-plasma fold change,
        ``Kp(knockout)/Kp(wild-type)``. Under the module's model this equals
        ``(1 + ER)/2`` and is therefore a direct measurement of ``1/f_efflux``.
        Used only when no in vitro ER exists; otherwise it serves as the
        consistency check.
    fold_ko_assay, fold_ko_citation : str
        Which knockout strain / route, and its PMID/DOI.
    notes : str
        Anything a reader needs in order to distrust the row appropriately:
        disputes in the literature, assay caveats, species splits,
        non-transporter resistance mechanisms, and so on.
    """

    drug: str
    pgp: str
    pgp_citation: str
    bcrp: str
    bcrp_citation: str
    er_pgp: Optional[float] = None
    er_bcrp: Optional[float] = None
    er_assay: str = ""
    er_citation: str = ""
    fold_ko: Optional[float] = None
    fold_ko_assay: str = ""
    fold_ko_citation: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        for field, val in (("pgp", self.pgp), ("bcrp", self.bcrp)):
            if val not in ("yes", "no", "unknown"):
                raise ValueError(
                    f"{self.drug}: {field}={val!r} must be 'yes', 'no' or 'unknown'"
                )
            if val != "unknown":
                cit = self.pgp_citation if field == "pgp" else self.bcrp_citation
                if not cit.strip():
                    raise ValueError(
                        f"{self.drug}: {field}={val!r} asserts a fact and must carry "
                        f"a citation. Use 'unknown' if there is no source."
                    )
        for field, val in (("er_pgp", self.er_pgp), ("er_bcrp", self.er_bcrp)):
            if val is None:
                continue
            if not (val > 0.0):
                raise ValueError(f"{self.drug}: {field}={val!r} must be positive")
            if not self.er_citation.strip():
                raise ValueError(
                    f"{self.drug}: {field} is a numeric claim and needs er_citation"
                )
        if self.er_pgp is not None and self.pgp != "yes":
            raise ValueError(f"{self.drug}: has an er_pgp but pgp={self.pgp!r}")
        if self.er_bcrp is not None and self.bcrp != "yes":
            raise ValueError(f"{self.drug}: has an er_bcrp but bcrp={self.bcrp!r}")
        if self.fold_ko is not None:
            if not (self.fold_ko > 0.0):
                raise ValueError(f"{self.drug}: fold_ko={self.fold_ko!r} must be positive")
            if not self.fold_ko_citation.strip():
                raise ValueError(f"{self.drug}: fold_ko is a numeric claim, cite it")
            if not self.is_substrate:
                raise ValueError(
                    f"{self.drug}: has a knockout fold change but neither pump is 'yes'"
                )

    @property
    def is_substrate(self) -> bool:
        """True if EITHER pump is a confirmed ``"yes"``."""
        return self.pgp == "yes" or self.bcrp == "yes"

    @property
    def er_in_vitro(self) -> Optional[float]:
        """Combined in vitro efflux ratio across both pumps, or ``None``."""
        return combine_efflux_ratios(self.er_pgp, self.er_bcrp)

    def er_from_knockout(self, form: str = FORM_KPUU) -> Optional[float]:
        """The in vivo knockout fold change expressed on ``form``'s ER scale."""
        if self.fold_ko is None:
            return None
        return _er_from_knockout_fold(self.fold_ko, form)

    @property
    def has_number(self) -> bool:
        """True if at least one usable published quantity exists (either window)."""
        return self.er_in_vitro is not None or self.fold_ko is not None

    def consistency(self, form: str = FORM_KPUU) -> Optional[Dict[str, float]]:
        """In vitro vs in vivo agreement, where BOTH observables exist.

        The in vitro ER implies a knockout fold change of
        :func:`knockout_fold_predicted_from_er`. Returns that prediction, the
        measured fold, and their ratio — greater than 1 means the monolayer assay
        UNDER-states the intact barrier, which is the direction the literature
        generally reports. ``None`` when one of the two observables is missing.
        """
        er = self.er_in_vitro
        if er is None or self.fold_ko is None:
            return None
        pred = knockout_fold_predicted_from_er(er, form)
        return {"fold_ko_predicted_from_er": pred,
                "fold_ko_measured": float(self.fold_ko),
                "in_vivo_over_in_vitro": float(self.fold_ko) / pred}


# ===========================================================================
# 3. The table
# ===========================================================================
#: LITERATURE-SOURCED substrate assignments for this repo's drug library.
#:
#: Covers all ten agents in :data:`tmesim.molecular.library.GBM_DRUGS`, the four
#: extra drugs the human glioma scorecard scores
#: (``tmesim/validation/glioma_regional_pk.py``), imatinib (the GIST validation),
#: and the three PARP inhibitors from the Moncayo explant validation.
#:
#: **Every non-``"unknown"`` cell carries a PMID.** ``"unknown"`` means the
#: literature was searched and nothing usable was found; it is treated as a
#: non-substrate downstream but is reported separately, because "measured not to
#: be a substrate" and "nobody has looked" are different states of knowledge.
#:
#: **Numeric ratios are only recorded when a specific value could be attributed
#: to a specific paper.** Several confirmed substrates therefore carry no number
#: at all (doxorubicin, methotrexate, vincristine, etoposide) — see each row's
#: ``notes`` for why. That is a real coverage gap and it is what the
#: ``substrate_but_no_published_er`` reason code exists to surface.
EFFLUX_TABLE: Dict[str, EffluxRecord] = {}


def _add(rec: EffluxRecord) -> None:
    EFFLUX_TABLE[rec.drug] = rec


# --- first-line GBM alkylators ---------------------------------------------
_add(EffluxRecord(
    drug="temozolomide",
    pgp="yes",
    pgp_citation="PMID 29852323 (de Gooijer 2018 Neoplasia): WEAK. No clear "
                 "transport in LLC-Mdr1a monolayers, but Abcb1a/b-/- mouse brain "
                 "accumulated ~20% more TMZ than WT.",
    bcrp="yes",
    bcrp_citation="PMID 29852323: TMZ is a GOOD Bcrp1 substrate — clear vectorial "
                  "transport in MDCKII-Bcrp1, abolished by elacridar 5 uM.",
    fold_ko=1.5,
    fold_ko_assay="in vivo brain-to-plasma fold change, Abcb1a/b;Abcg2-/- vs WT "
                  "FVB mice, 50 mg/kg i.v.; reproduced pharmacologically by elacridar",
    fold_ko_citation="PMID 29852323 (de Gooijer 2018 Neoplasia), p < 0.001, "
                     "plasma clearance unchanged",
    notes="NO in vitro ER exists and none can: the authors state TMZ's ~30 min "
          "half-life at pH 7.4 makes Papp in cm/s uncalculable. The combined "
          "knockout effect is only 1.5x (1.2x for each pump alone) — i.e. TMZ is "
          "a real but very weak efflux substrate, which corrects the common "
          "assumption that it is not one at all. Goldwirt 2014 (PMID 24867782) "
          "asserts ABCB1 transport in its abstract but reports a mdr1a-/- vs +/+ "
          "brain:plasma ratio of 1.1, i.e. no effect; read its number, not its claim.",
))

_add(EffluxRecord(
    drug="lomustine",
    pgp="unknown",
    pgp_citation="",
    bcrp="unknown",
    bcrp_citation="",
    notes="No transport assay of any kind exists — no bidirectional monolayer, no "
          "transfected-cell transport study, no Mdr1a/b-/- or Bcrp1-/- mouse PK. "
          "The one relevant primary result is a RESISTANCE readout, not a transport "
          "measurement: PMID 36252775 (Radtke 2022 Pharmacol Res) found that "
          "CRISPR ABCB1 knockout in glioblastoma cells improved the temozolomide "
          "response and altered the carmustine response but did NOT alter the "
          "lomustine response. That is suggestive of non-substrate status and is "
          "recorded here, but a viability assay is not a substrate assay, so the "
          "status stays UNKNOWN rather than being upgraded to 'no'.",
))

_add(EffluxRecord(
    drug="carmustine",
    pgp="unknown",
    pgp_citation="",
    bcrp="unknown",
    bcrp_citation="",
    notes="No primary substrate study exists. Searches for carmustine/BCNU against "
          "P-glycoprotein, ABCB1, ABCG2/BCRP, efflux transporters and knockout-mouse "
          "PK returned zero substrate assays; every carmustine+MDR1 hit is a "
          "gene-therapy chemoprotection study in which MDR1 is transduced to protect "
          "marrow against a DIFFERENT drug while MGMT protects against BCNU (e.g. "
          "PMID 10933194). Mechanistically BCNU is small, lipophilic and hydrolyses "
          "within minutes, which makes it both a poor efflux candidate and hard to "
          "assay — but that is reasoning, not a citation.",
))

_add(EffluxRecord(
    drug="procarbazine",
    pgp="unknown",
    pgp_citation="",
    bcrp="unknown",
    bcrp_citation="",
    notes="No primary literature found. The broadest productive query "
          "('procarbazine transport') returns ~13 records, none addressing "
          "ABC-transporter substrate status. Circumstantially procarbazine is small "
          "(MW 221), lipophilic, and used in the CNS-directed PCV regimen, which is "
          "consistent with low efflux liability — an inference, not data, and not a "
          "basis for a numeric term.",
))

# --- antimetabolite ---------------------------------------------------------
_add(EffluxRecord(
    drug="methotrexate",
    pgp="no",
    pgp_citation="PMID 38981917 (Volpe 2024, FDA): MTX efflux ratio measured in "
                 "C2BBe1 wild-type vs single/double CRISPR knockouts falls in "
                 "BCRP-KO and MRP2-KO cells but NOT in P-gp-KO cells.",
    bcrp="yes",
    bcrp_citation="PMID 14500392 (Volk & Schneider 2003 Cancer Res): ATP-dependent "
                  "vesicular transport by WILD-TYPE (R482) BCRP, Km 680 uM, Vmax "
                  "2400 pmol/mg/min; PMID 12874005 (Chen 2003): Km 1.34 mM, Vmax "
                  "687 pmol/mg/min. MTX-Glu2/Glu3 also transported. In vivo "
                  "confirmation PMID 15313923, PMID 19996279.",
    notes="A confirmed BCRP substrate with NO usable number for this purpose. The "
          "available in vivo figure (1.8-1.9x, PMID 15313923) is a PLASMA CLEARANCE "
          "change in Bcrp1-/- mice, not a brain- or tissue-to-plasma penetration "
          "fold change, so it does not measure the quantity this module needs and is "
          "deliberately NOT used as fold_ko. Note also the mirror-image of "
          "doxorubicin: the R482T/R482G drug-selected BCRP mutants LOSE MTX "
          "transport while gaining doxorubicin transport. MTX is additionally an "
          "OAT/OATP/RFC substrate, which is out of this module's scope and is part "
          "of why an efflux-only term is incomplete for it. One contradicting P-gp "
          "claim exists (PMID 30977678) from Caco-2 accumulation rather than "
          "bidirectional flux; not weighted.",
))

# --- anthracycline / platinum / topo-II / vinca ------------------------------
_add(EffluxRecord(
    drug="doxorubicin",
    pgp="yes",
    pgp_citation="PMID 3472246 (Ueda 1987 PNAS): MDR1 cDNA transfection confers "
                 "doxorubicin resistance — the canonical demonstration. In vivo: "
                 "PMID 10681369 (van Asperen 2000), biliary secretion of unchanged "
                 "doxorubicin falls 13.3% -> 2.4% of dose in mdr1a-/- mice (~5.5x).",
    bcrp="no",
    bcrp_citation="PMID 11559526 (Honjo 2001) and PMID 14612912 (Robey 2003): "
                  "WILD-TYPE R482 human BCRP does NOT transport doxorubicin; only "
                  "the R482T/R482G drug-selected mutants do (mutant transfectants "
                  "13-71x resistant vs 3-4x for wild-type).",
    notes="THE ANCHOR CASE, AND IT CARRIES NO NUMBER. Doxorubicin is the textbook "
          "P-gp substrate and the drug this repo's perivascular-uptake constant "
          "alpha was originally fitted to, yet no defensible single in vitro ER "
          "could be attributed: PMID 39669004 (Kotze 2024) ran it bidirectionally "
          "in MDCK-II-MDR1 and states ER >= 2, but publishes ER only as pH curves "
          "in figures with values in unretrieved supplementary material. Recording "
          "'ER ~ 2' from a figure would be a guess, so the factor stays 1.0. The "
          "widespread belief that doxorubicin is a BCRP substrate traces entirely "
          "to mutant cell lines and is coded here as 'no' for wild-type BCRP.",
))

_add(EffluxRecord(
    drug="cisplatin",
    pgp="no",
    pgp_citation="PMID 22792399 (Stordal 2012 PLoS One): elacridar does not reverse "
                 "cisplatin resistance in a P-gp overexpressor; PMID 24403508 "
                 "(Doherty 2014): P-gp-overexpressing KB-8-5-11 cells are "
                 "COLLATERALLY SENSITIVE to cisplatin; PMID 16534613 (Takara 2006): "
                 "cisplatin-selected HeLa show DECREASED MDR1 expression/activity.",
    bcrp="no",
    bcrp_citation="PMID 18820913, PMID 19638344 — negative in BCRP substrate screens.",
    notes="Expectation verified rather than assumed: platinum resistance runs "
          "through copper transporters (CTR1/ATP7A/ATP7B) and glutathione/GST, not "
          "ABC efflux. One trap avoided: PMID 15139851 (Alqawi 2004) shows cisplatin "
          "competes for a photoaffinity label on ABCG2 — that is BINDING, not "
          "transport, and must not be coded as substrate evidence.",
))

_add(EffluxRecord(
    drug="etoposide",
    pgp="yes",
    pgp_citation="PMID 20028753 (Lagas 2010 Clin Cancer Res): P-gp restricts oral "
                 "(re)uptake and mediates gut-wall excretion; PMID 12649196 (Allen "
                 "2003): GF120918 in WT mice raised oral etoposide plasma exposure "
                 "4-5x; PMID 11245453: mdr1a/1b;mrp1 triple-KO mice are 3-5x more "
                 "etoposide-sensitive.",
    bcrp="unknown",
    bcrp_citation="",
    notes="Genuinely split by species/system, hence BCRP unknown: mouse Bcrp1 "
          "transports etoposide across MDCKII monolayers and confers resistance "
          "reversible by Ko143 and GF120918 (PMID 12649196), while human ABCG2 in "
          "HEK293 conferred no appreciable etoposide resistance (PMID 14612912). "
          "NO NUMBER IS RECORDED. A candidate ER of 5.38 (C2BBe1 wild-type, falling "
          "to 3.65 in MDR1-KO and 1.24 in MDR1/MRP2 double-KO, PMID 38981917) could "
          "not be verified against the source table, and in any case that wild-type "
          "ratio is dominated by MRP2/ABCC2 rather than P-gp — etoposide is "
          "primarily an MRP2 substrate, so attributing its whole efflux to P-gp "
          "would over-state the term this module models.",
))

_add(EffluxRecord(
    drug="vincristine",
    pgp="yes",
    pgp_citation="PMID 20501632 (Wang 2010 Neuro-Oncology); PMID 11245453 (Johnson "
                 "2001): mdr1a/1b;mrp1 triple-KO mice show 128x increased vincristine "
                 "toxicity; PMID 11489261.",
    bcrp="no",
    bcrp_citation="PMID 18820913 — negative in a BCRP substrate screen.",
    notes="Substrate status is not in doubt but the QUANTITATIVE brain evidence is "
          "weaker than its reputation. The only brain-specific figure found (PMID "
          "20501632) is a 1.3-3.6x normal-brain partition increase in a TRIPLE "
          "Pgp+Mrp1 knockout, so it is not attributable to P-gp and BCRP alone and "
          "is not recorded as fold_ko. No numeric monolayer ER was located. Beware "
          "also that much of the canonical 'vinca + P-gp' literature is VINBLASTINE "
          "(e.g. PMID 7910522), which is not interchangeable with vincristine for a "
          "quantitative term.",
))

# --- macromolecule marker ---------------------------------------------------
_add(EffluxRecord(
    drug="bevacizumab",
    pgp="no",
    pgp_citation="PMID 23421687 (Desai 2013 ACS Chem Neurosci): ER = (P_pass + "
                 "Vmax/Km)/P_pass is undefined as P_pass -> 0; PMID 27572805 "
                 "(Pardridge 2016): therapeutic antibodies do not cross the BBB and "
                 "must be re-engineered for receptor-mediated transcytosis.",
    bcrp="no",
    bcrp_citation="PMID 23421687; PMID 27572805 — same argument.",
    notes="NOT APPLICABLE rather than measured-negative. ABC efflux transporters act "
          "on solutes that partition into the lipid bilayer; a ~149 kDa IgG has "
          "essentially no passive transcellular permeability, so the efflux ratio is "
          "ill-posed, not merely unmeasured. Bevacizumab's barrier restriction "
          "belongs in a size/permeability term, and this repo already handles it "
          "there (bbb_penetrance = 0 via the is_antibody and MW >= 1000 gates). "
          "efflux_factor short-circuits antibodies to 1.0 before this row is read.",
))

# --- glioma-scorecard drugs outside the library ------------------------------
_add(EffluxRecord(
    drug="levetiracetam",
    pgp="unknown",
    pgp_citation="",
    bcrp="no",
    bcrp_citation="PMID 25645391 (Roemermann 2015 Neuropharmacology): no transport "
                  "in MDCKII-Bcrp1 or MDCKII-BCRP concentration-equilibrium assays "
                  "at therapeutic concentrations, in a paper whose lamotrigine "
                  "positive control DID transport — a well-controlled negative.",
    notes="P-gp status is genuinely DISPUTED, and the dispute looks like a species "
          "and assay-format artefact, so 'unknown' is the honest code. FOR: PMID "
          "18824002 (Luna-Tortos 2008) found directional transport by human P-gp in "
          "a concentration-equilibrium transport assay — used precisely because "
          "levetiracetam's high passive permeability masks low-affinity efflux in a "
          "conventional bidirectional assay — though it was also MK571-inhibitable, "
          "so MRPs contribute. AGAINST: PMID 17045309 (Baltes 2006) found transport "
          "by mouse mdr1a but NOT human MDR1; PMID 15120740 (Potschka 2004) found "
          "verapamil and probenecid did not raise cortical levetiracetam in rat "
          "microdialysis; PMID 17916018 rated it only a weak canine P-gp substrate. "
          "No numeric ER exists either way — a CETA reports AUC differences, not a "
          "ratio.",
))

_add(EffluxRecord(
    drug="cefazolin",
    pgp="unknown",
    pgp_citation="",
    bcrp="yes",
    bcrp_citation="PMID 18339814 (Kato 2008 Drug Metab Dispos): ATP-dependent uptake "
                  "of cefazolin into RAT Bcrp-expressing Sf9 membrane vesicles.",
    notes="Two caveats make this row weak in both cells. P-gp: there is NO substrate "
          "study at all — four PubMed hits for cefazolin + P-glycoprotein, none "
          "testing substrate status — hence unknown, not 'no'. BCRP: the positive is "
          "RAT Bcrp, measured as vesicular uptake in a biliary-excretion context, "
          "not a human-BCRP monolayer ratio, and human BCRP was not tested. No ER "
          "exists, so the factor is 1.0 regardless. Cefazolin's real transporter "
          "biology is organic-anion, not ABC-efflux: a confirmed OAT1 substrate "
          "(PMID 15618660, Ki 72 uM vs PAH) and an MRP4 substrate (PMID 17344354, "
          "Km 80 uM, kidney:plasma up 2.7x in Mrp4-/- mice).",
))

_add(EffluxRecord(
    drug="adavosertib",
    pgp="yes",
    pgp_citation="PMID 28928160 (Li 2017 Clin Cancer Res); independently confirmed "
                 "PMID 29147815 (de Gooijer 2018 Invest New Drugs).",
    bcrp="yes",
    bcrp_citation="PMID 28928160; PMID 29147815.",
    er_pgp=10.6,
    er_bcrp=4.5,
    er_assay="bidirectional MDCKII-ABCB1 (ER 10.6) and MDCKII-ABCG2 (ER 4.5) "
             "monolayers at basolateral pH 7.4; abolished by cyclosporine A 10 uM "
             "and Ko143 1 uM respectively; loperamide/gefitinib positive controls, "
             "85-115% mass recovery",
    er_citation="PMID 28928160 (Li 2017 Clin Cancer Res)",
    fold_ko=25.0,
    fold_ko_assay="in vivo brain-to-plasma ratio, wild-type vs Abcb1a/b;Abcg2-/- mice",
    fold_ko_citation="PMID 29147815 (de Gooijer 2018 Invest New Drugs)",
    notes="The best-evidenced drug in the table — two assay types, two independent "
          "groups — and, for our purposes, the most dangerous. Two warnings. (1) The "
          "ERs are pH-DEPENDENT: ABCB1 falls 9.2 -> 3.7 and ABCG2 4.7 -> 2.6 as "
          "basolateral pH drops 7.4 -> 6.5, so the acidic tumour microenvironment "
          "roughly halves both. The pH 7.4 values are recorded, which is the "
          "efflux-unfavourable end for a tumour. (2) DIRECTLY CONTRADICTED IN OUR "
          "OWN TISSUE: despite dual efflux, Li 2017 measured human GBM tumour Kp,uu "
          "at a MEDIAN OF 3.2 (range 1.3-24.4) — i.e. adavosertib ACCUMULATES in "
          "human glioma, because OATP1A2-mediated uptake dominates in the acidic "
          "tumour microenvironment. The mouse and human results genuinely disagree "
          "and the authors say so explicitly. Applying this row to human glioma is "
          "contradicted by a direct human measurement of the same drug in the same "
          "tissue.",
))

_add(EffluxRecord(
    drug="eribulin",
    pgp="yes",
    pgp_citation="PMID 21162698 (Taur 2010 Xenobiotica).",
    bcrp="unknown",
    bcrp_citation="",
    fold_ko=30.0,
    fold_ko_assay="in vivo brain penetration, abcb1a-deficient CF-1 mice vs CF-1 "
                  "wild-type (oral bioavailability also 62.3% vs 7.6%)",
    fold_ko_citation="PMID 21162698 (Taur 2010 Xenobiotica)",
    notes="The strongest single efflux number in the table and the drug this whole "
          "exercise was motivated by. Caveats: it rests on ONE primary paper (an "
          "Eisai in-house study) and the CF-1 model knocks out abcb1a ONLY, not "
          "abcb1a/1b. No in vitro ER exists — the Caco-2 work in the same paper used "
          "eribulin as an INHIBITOR of digoxin efflux (IC50 > 10 uM), which is "
          "inhibitor data and must not be conflated with substrate data. BCRP is "
          "coded unknown, not 'no', because the only source for the negative is the "
          "FDA Halaven label section 12.3, whose underlying in vitro study is "
          "unpublished and carries no PMID. Counterpoint worth carrying: the current "
          "label reports no clinically meaningful DDI with ketoconazole or rifampin, "
          "so the 30x mouse brain effect does not translate to a clinical IV "
          "interaction.",
))

# --- GIST validation ---------------------------------------------------------
_add(EffluxRecord(
    drug="imatinib",
    pgp="yes",
    pgp_citation="PMID 12604685 (Dai 2003 JPET): B->A flux 39x greater than A->B in "
                 "MDCKII-MDR1, reversed by zosuquidar.",
    bcrp="yes",
    bcrp_citation="PMID 15251980 (Burger 2004 Blood): established in "
                  "BCRP-overexpressing lines with mitoxantrone competition and Ko143 "
                  "reversal. IN VITRO ONLY — see notes.",
    er_pgp=39.0,
    er_bcrp=None,
    er_assay="bidirectional MDCKII-MDR1, Papp(B->A)/Papp(A->B) = 39 (parental "
             "MDCKII background 8, so the transfectant-specific component is real "
             "but the background is not negligible)",
    er_citation="PMID 12604685 (Dai 2003 JPET)",
    fold_ko=28.0,
    fold_ko_assay="in vivo brain-to-plasma, Mdr1a/1b;Bcrp TRIPLE knockout vs WT",
    fold_ko_citation="PMID 19225039 (Bihorel/Oberoi 2009 Drug Metab Dispos)",
    notes="P-gp is unambiguous and replicated: 11x brain:plasma at 1 h after oral "
          "dosing and ~7x after IV in mdr1a/b-/- mice (PMID 12604685), 4x in "
          "Mdr1a/1b-/- (PMID 19225039), 5.5x by in situ brain perfusion (PMID "
          "17696988), 6.4x IV (PMID 17380257). BCRP is WEAKER than its reputation "
          "and the two windows disagree: in vitro yes, but in vivo Bcrp(-/-) mice "
          "show brain:plasma comparable to wild-type (PMID 19225039) and Bcrp1(-/-) "
          "uptake matches WT at non-saturating concentrations (PMID 17696988). The "
          "28x figure is the P-gp/Bcrp SYNERGY in the triple knockout and must not "
          "be attributed to BCRP alone. NOTE FOR THE GIST VALIDATION: none of this "
          "helps there. A scalar per-drug factor cancels exactly out of a "
          "tumour:normal ratio, and hepatic P-gp/BCRP sit on the CANALICULAR "
          "membrane exporting imatinib into bile, which lowers NORMAL liver and "
          "pushes the predicted ratio the wrong way.",
))

# --- Moncayo explant validation (PARP inhibitors) ----------------------------
_add(EffluxRecord(
    drug="niraparib",
    pgp="yes",
    pgp_citation="PMID 34730366 (Martins/Loos/Schinkel 2021 Mol Pharm): human ABCB1 "
                 "transported niraparib moderately in MDCK-II cells; oral elacridar "
                 "in WT mice reproduced the knockout phenotype.",
    bcrp="yes",
    bcrp_citation="PMID 34730366 (in vitro, mouse Abcg2); PMID 41998207 (Macdonald "
                  "2026 Br J Cancer, polarised-cell efflux assays). IN VITRO ONLY: "
                  "single Abcg2-/- mice were indistinguishable from wild-type, so "
                  "BCRP is not BBB-limiting for niraparib in vivo.",
    fold_ko=6.5,
    fold_ko_assay="in vivo brain-to-plasma fold change, Abcb1a/1b-/- and "
                  "Abcb1a/1b;Abcg2-/- vs WT mice (reported as a 6- to 7-fold range; "
                  "6.5 is the midpoint of the authors' own range, not a fitted value)",
    fold_ko_citation="PMID 34730366 (Martins 2021 Mol Pharm)",
    notes="No numeric in vitro transwell ER is published for ANY of the three PARP "
          "inhibitors in accessible sources — the source papers report transport "
          "qualitatively or place ratios in paywalled supplements. Niraparib is the "
          "brain-penetrant member of the class: measured brain:plasma AUC 0.29-0.33 "
          "vs ~0.03 for olaparib in two xenograft models (PMID 30647846), "
          "corroborated by NHP and mouse imaging (PMID 41023116). That ~10x measured "
          "separation is a useful external check on this module and it does not pass "
          "it — see results/efflux_impact.md.",
))

_add(EffluxRecord(
    drug="rucaparib",
    pgp="yes",
    pgp_citation="PMID 24962512 (Durmus/Schinkel 2015 Pharm Res); PMID 26438157 "
                 "(Parrish/Elmquist/Sarkaria 2015 Mol Cancer Ther): intracellular "
                 "rucaparib reduced to 32.0+/-7.4% of control in MDR1 cells, "
                 "reversed by zosuquidar.",
    bcrp="yes",
    bcrp_citation="PMID 24962512 (very efficiently transported by mouse Abcg2); "
                  "PMID 26438157 (intracellular rucaparib 14.5+/-5.2% of control in "
                  "BCRP1 cells, reversed by Ko143).",
    fold_ko=14.6,
    fold_ko_assay="steady-state osmotic-pump brain-to-plasma, Mdr1a/1b-/-Bcrp1-/- "
                  "triple knockout (1.61+/-0.25) vs WT (0.11+/-0.08)",
    fold_ko_citation="PMID 26438157 (Parrish 2015 Mol Cancer Ther), abstract values",
    notes="INTERNAL INCONSISTENCY IN THE SOURCE, recorded rather than smoothed: the "
          "abstract gives 1.61 vs 0.11 (14.6x), the Results section of the same "
          "paper gives 2.0+/-0.37 vs 0.12+/-0.11 (~17x), and the Discussion says "
          "'almost 20-fold'. The abstract value is used because it is the most "
          "conservative and the most citable; the true uncertainty on this row is at "
          "least 14.6-20x. Two independent groups agree on dual P-gp/BCRP substrate "
          "status.",
))

_add(EffluxRecord(
    drug="olaparib",
    pgp="yes",
    pgp_citation="PMID 35341895 (Song 2022 Eur J Pharm Sci): brain-to-plasma UNBOUND "
                 "ratio rose 5.6x in Abcb1a/1b-/- mice vs WT; elacridar reproduced "
                 "the effect.",
    bcrp="yes",
    bcrp_citation="PMID 30396915 (Song 2018 JPET): efficiently transported by human "
                  "ABCG2 in vitro — the component is masked by endogenous canine "
                  "ABCB1 in MDCKII and only appears once that is suppressed "
                  "(corrected ER rose 3x under verapamil, 5x under shRNA knockdown). "
                  "IN VIVO MINOR: single Abcg2-/- was indistinguishable from WT.",
    fold_ko=8.1,
    fold_ko_assay="in vivo brain-to-plasma UNBOUND concentration ratio (Kp,uu), "
                  "Abcb1a/1b;Abcg2-/- vs WT mice",
    fold_ko_citation="PMID 35341895 (Song 2022 Eur J Pharm Sci)",
    notes="Abcb1 dominates (5.6x alone) and Abcg2 adds little in vivo (8.1x "
          "combined) despite efficient human ABCG2 transport in vitro — a clean "
          "example of why in vitro substrate status does not imply an in vivo "
          "barrier effect. Note this quantity is a Kp,uu ratio, i.e. already "
          "unbound-corrected, which is exactly the quantity FORM_KPUU models.",
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


def _is_antibody(drug: Union[str, Mapping]) -> bool:
    if isinstance(drug, Mapping):
        return bool(drug.get("is_antibody", False))
    return bool(getattr(drug, "is_antibody", False))


#: Which published observable to prefer when a drug has both. Fixed in advance
#: (see the module docstring) so the choice cannot be made per drug to flatter a
#: result. ``"in_vitro"`` = combined monolayer ER; ``"in_vivo"`` = knockout fold
#: change re-expressed as ``ER = 2*fold - 1``.
PREFER_IN_VITRO = "in_vitro"
PREFER_IN_VIVO = "in_vivo"
SOURCE_PREFERENCES = (PREFER_IN_VITRO, PREFER_IN_VIVO)


def efflux_ratio_for(drug: Union[str, Mapping],
                     record: Optional[EffluxRecord] = None,
                     prefer: str = PREFER_IN_VITRO,
                     form: str = FORM_KPUU) -> Optional[float]:
    """Best available efflux ratio for a drug, or ``None`` if there is none.

    Returns an ER on ``form``'s scale in both cases: an in vivo knockout fold
    change is converted by :func:`_er_from_knockout_fold`, the exact inverse of
    that form's ``f(ER)``. ``prefer`` selects which observable wins when both
    exist; the fallback is always the other one.
    """
    if prefer not in SOURCE_PREFERENCES:
        raise ValueError(f"prefer must be one of {SOURCE_PREFERENCES}, got {prefer!r}")
    rec = record if record is not None else EFFLUX_TABLE.get(_resolve_name(drug) or "")
    if rec is None or not rec.is_substrate:
        return None
    a, b = rec.er_in_vitro, rec.er_from_knockout(form)
    if prefer == PREFER_IN_VIVO:
        a, b = b, a
    return a if a is not None else b


def efflux_detail(drug: Union[str, Mapping],
                  form: str = FORM_KPUU,
                  efflux_ratio: Optional[float] = None,
                  prefer: str = PREFER_IN_VITRO) -> Dict[str, object]:
    """Full, auditable explanation of the factor returned by :func:`efflux_factor`.

    Returns a dict with ``factor``, ``efflux_ratio``, ``reason``, ``applied``,
    ``pgp``, ``bcrp``, ``citations`` and ``notes``. The ``reason`` field is the
    one to read when a factor is 1.0, because 1.0 has four different meanings:

    ==================================  ======================================
    reason                              meaning
    ==================================  ======================================
    ``"not_in_table"``                  never looked up; treated as non-substrate
    ``"antibody"``                      macromolecule; ABC transporters N/A
    ``"non_substrate"``                 MEASURED not to be a substrate
    ``"unknown_substrate_status"``      searched, no usable source found
    ``"substrate_but_no_published_er"`` confirmed substrate, no number to use
    ``"applied"``                       a real published ratio was applied
    ==================================  ======================================
    """
    key = _resolve_name(drug)
    base = {
        "drug": key, "factor": 1.0, "efflux_ratio": None, "form": form,
        "applied": False, "pgp": None, "bcrp": None, "citations": "", "notes": "",
        "source": None, "consistency": None,
    }
    if efflux_ratio is not None:
        f = _apply_form(float(efflux_ratio), form)
        return {**base, "factor": f, "efflux_ratio": float(efflux_ratio),
                "applied": True, "reason": "explicit_override",
                "notes": "efflux_ratio supplied by caller, table bypassed"}
    if _is_antibody(drug):
        return {**base, "reason": "antibody",
                "notes": "a ~149 kDa IgG does not cross membranes passively, so an "
                         "apical membrane transporter has nothing to act on"}
    rec = EFFLUX_TABLE.get(key or "")
    if rec is None:
        return {**base, "reason": "not_in_table",
                "notes": "no literature assignment exists in EFFLUX_TABLE; treated "
                         "as a non-substrate, which may be wrong"}
    base.update({"pgp": rec.pgp, "bcrp": rec.bcrp, "notes": rec.notes,
                 "consistency": rec.consistency(form),
                 "citations": "; ".join(
                     c for c in (rec.pgp_citation, rec.bcrp_citation,
                                 rec.er_citation, rec.fold_ko_citation) if c)})
    if not rec.is_substrate:
        # "unknown" wins over "no": if EITHER pump is unresolved then the drug's
        # substrate status as a whole is unresolved, even where the other pump
        # has a clean measured negative (levetiracetam is exactly this case).
        any_unknown = rec.pgp == "unknown" or rec.bcrp == "unknown"
        return {**base,
                "reason": "unknown_substrate_status" if any_unknown else "non_substrate"}
    er = efflux_ratio_for(drug, record=rec, prefer=prefer, form=form)
    if er is None:
        return {**base, "reason": "substrate_but_no_published_er"}
    other = PREFER_IN_VIVO if prefer == PREFER_IN_VITRO else PREFER_IN_VITRO
    primary = rec.er_in_vitro if prefer == PREFER_IN_VITRO else rec.er_from_knockout(form)
    source = prefer if primary is not None else other
    return {**base, "factor": _apply_form(er, form), "efflux_ratio": er,
            "applied": True, "reason": "applied", "source": source}


def efflux_factor(drug: Union[str, Mapping],
                  form: str = FORM_KPUU,
                  efflux_ratio: Optional[float] = None,
                  prefer: str = PREFER_IN_VITRO) -> float:
    """Penetrance multiplier in ``(0, 1]`` for active efflux. **1.0 = no effect.**

    Multiply a drug's passive/transcellular penetrance by this to obtain the
    net penetrance in the presence of P-gp / BCRP:

    ::

        net_penetrance = passive_penetrance * efflux_factor(drug)

    with ``efflux_factor = 2 / (1 + ER_combined)`` by default (see the module
    docstring for the derivation, and :func:`combine_efflux_ratios` for how two
    pumps are combined).

    **Everything unproven returns exactly 1.0.** A drug that is absent from
    :data:`EFFLUX_TABLE`, is recorded as a non-substrate, has ``"unknown"``
    status, or is a confirmed substrate with no published ratio, all return 1.0.
    Use :func:`efflux_detail` to find out which of those happened — the four
    cases are very different and the report must not conflate them.

    Parameters
    ----------
    drug : str, Mapping, or object with ``.name``
        A drug name (e.g. ``"eribulin"``), a descriptor mapping as produced by
        :func:`tmesim.molecular.descriptors.compute_descriptors` (with a
        ``"name"`` key, and honouring ``"is_antibody"``), or any object carrying
        a ``.name`` attribute such as a :class:`tmesim.core.DrugParams`.
    form : {"kpuu", "monolayer"}
        Functional form. ``"kpuu"`` (default) is ``1/ER``, the single-barrier
        unbound-partition relation the field uses. ``"monolayer"`` is
        ``2/(1+ER)``, which resolves the cell interior and is exactly 2x more
        permissive at high ER; provided for sensitivity analysis.
    efflux_ratio : float, optional
        Bypass the table with an explicit ratio. Intended for sweeps and for
        drugs whose ratio a caller has from a source not carried here. Values
        below 1 are clamped to 1 (this module never boosts penetrance).
    prefer : {"in_vitro", "in_vivo"}
        Which published observable wins when a drug has both a monolayer ER and
        a knockout fold change. Defaults to the pre-committed ``"in_vitro"``.

    Returns
    -------
    float
        Multiplier in ``(0, 1]``.
    """
    return float(efflux_detail(drug, form=form, efflux_ratio=efflux_ratio,
                               prefer=prefer)["factor"])


def efflux_adjusted_penetrance(passive_penetrance: float,
                               drug: Union[str, Mapping],
                               bbb_disruption: float = 0.0,
                               leak: Optional[float] = None,
                               form: str = FORM_KPUU,
                               efflux_ratio: Optional[float] = None,
                               prefer: str = PREFER_IN_VITRO) -> float:
    """Compose efflux with the regional paracellular leak, in the right order.

    ::

        P_eff = clip( P_passive * f_efflux  +  L * bbb_disruption , 0, 1 )

    **Efflux multiplies the transcellular term only.** P-gp and BCRP sit in the
    luminal membrane of the endothelial cell, so they can only act on drug that
    takes the transcellular route. A disrupted tight junction is a paracellular
    shunt *around* the cell and therefore around the pump — which is precisely
    why efflux inhibition helps least where the barrier is already broken.
    Multiplying the whole bracket instead would assert the opposite, so the
    ordering here is a physical claim, not a coding convenience.

    With ``bbb_disruption = 0`` (the default) this reduces exactly to
    ``passive_penetrance * efflux_factor(...)``, and with a non-substrate it
    reduces exactly to :func:`tmesim.reconstruct.barrier.effective_penetrance`.

    Parameters
    ----------
    passive_penetrance : float
        Structure-only penetrance in [0, 1], e.g. from
        :func:`tmesim.molecular.permeability.bbb_penetrance`.
    drug : str or Mapping
        Passed to :func:`efflux_factor`.
    bbb_disruption : float
        Regional barrier disruption in [0, 1]; 0 = intact.
    leak : float, optional
        Additive paracellular leak constant. Defaults to
        :data:`tmesim.reconstruct.barrier.L_PARACELLULAR_LEAK` (0.238), imported
        lazily so this module does not depend on ``reconstruct`` unless the
        regional term is actually in play.
    form, efflux_ratio
        Passed to :func:`efflux_factor`.

    Returns
    -------
    float
        Effective penetrance in [0, 1].
    """
    f = efflux_factor(drug, form=form, efflux_ratio=efflux_ratio, prefer=prefer)
    transcellular = float(passive_penetrance) * f
    d = float(bbb_disruption)
    if d == 0.0:
        return float(min(1.0, max(0.0, transcellular)))
    if leak is None:
        from ..reconstruct.barrier import L_PARACELLULAR_LEAK
        leak = L_PARACELLULAR_LEAK
    if not (-1e-12 <= d <= 1.0 + 1e-12):
        raise ValueError(f"bbb_disruption must lie in [0, 1]; got {d}")
    return float(min(1.0, max(0.0, transcellular + float(leak) * d)))


def substrate_summary() -> Dict[str, int]:
    """Counts of substrate assignments in :data:`EFFLUX_TABLE`, for the report."""
    out = {"n_drugs": len(EFFLUX_TABLE), "pgp_yes": 0, "pgp_no": 0, "pgp_unknown": 0,
           "bcrp_yes": 0, "bcrp_no": 0, "bcrp_unknown": 0,
           "substrate_with_number": 0, "substrate_without_number": 0,
           "non_substrate": 0}
    for rec in EFFLUX_TABLE.values():
        out[f"pgp_{rec.pgp}"] += 1
        out[f"bcrp_{rec.bcrp}"] += 1
        if rec.is_substrate:
            out["substrate_with_number" if rec.has_number
                else "substrate_without_number"] += 1
        else:
            out["non_substrate"] += 1
    return out


# ===========================================================================
# 6. Self-test
# ===========================================================================
if __name__ == "__main__":
    # --- the two functional forms are what the docstring derives -----------
    # Three-compartment monolayer, simulated directly from the mass balance,
    # to check that FORM_MONOLAYER really is the exact solution.
    for P in (0.4, 1.0, 3.0):
        for E in (0.0, 1.0, 5.0, 50.0):
            c_cell = P / (2.0 * P + E)          # unit donor, sink receiver
            papp_ab = P * c_cell
            papp_ba = (P + E) * c_cell
            er = papp_ba / papp_ab
            assert abs(er - (1.0 + E / P)) < 1e-12, "ER = 1 + E/P"
            # net vs passive-only (two membranes in series => P/2)
            assert abs(papp_ab / (P / 2.0) - _apply_form(er, FORM_MONOLAYER)) < 1e-12
            # knockout fold is 1/f in BOTH forms, by construction
            for form in EFFLUX_FORMS:
                f = _apply_form(er, form)
                fold = 1.0 / f
                assert abs(_er_from_knockout_fold(fold, form) - er) < 1e-9, form
                assert abs(knockout_fold_predicted_from_er(er, form) - fold) < 1e-9

    # both forms are exactly 1.0 for a non-substrate, and never boost
    for form in EFFLUX_FORMS:
        assert _apply_form(1.0, form) == 1.0, f"{form} must be 1.0 at ER=1"
        assert _apply_form(0.2, form) == 1.0, f"{form} must clamp ER<1 to no effect"
        assert _apply_form(100.0, form) < 0.05
    assert _apply_form(50.0, FORM_MONOLAYER) > _apply_form(50.0, FORM_KPUU), \
        "monolayer form must be the more permissive one"

    # --- parallel-pump combination ----------------------------------------
    assert combine_efflux_ratios(None, None) is None
    assert combine_efflux_ratios(5.0, None) == 5.0
    assert combine_efflux_ratios(10.6, 4.5) == 14.1                # adavosertib
    assert combine_efflux_ratios(1.0, 1.0) == 1.0                  # two non-pumps

    # --- everything unproven is exactly 1.0 --------------------------------
    for probe, expect in (("no_such_drug", "not_in_table"),
                          ("cisplatin", "non_substrate"),
                          ("carmustine", "unknown_substrate_status"),
                          ("doxorubicin", "substrate_but_no_published_er")):
        d = efflux_detail(probe)
        assert d["factor"] == 1.0, (probe, d)
        assert d["reason"] == expect, (probe, d["reason"])
    assert efflux_factor({"name": "bevacizumab", "is_antibody": True}) == 1.0
    assert efflux_detail({"name": "x", "is_antibody": True})["reason"] == "antibody"

    # --- the drugs that DO carry a number ----------------------------------
    assert abs(efflux_factor("eribulin") - 1.0 / 30.0) < 1e-12
    assert abs(efflux_factor("adavosertib") - 1.0 / 14.1) < 1e-12
    assert abs(efflux_factor("temozolomide") - 1.0 / 1.5) < 1e-12
    # an in vivo fold gives the SAME factor in both forms
    for form in EFFLUX_FORMS:
        assert abs(efflux_factor("eribulin", form=form) - 1.0 / 30.0) < 1e-12, form
    # an in vitro ER does not
    assert efflux_factor("adavosertib", form=FORM_MONOLAYER) > \
        efflux_factor("adavosertib", form=FORM_KPUU)

    # --- every asserted fact carries a citation ----------------------------
    for rec in EFFLUX_TABLE.values():
        if rec.pgp != "unknown":
            assert rec.pgp_citation.strip(), rec.drug
        if rec.bcrp != "unknown":
            assert rec.bcrp_citation.strip(), rec.drug
        if rec.er_in_vitro is not None:
            assert "PMID" in rec.er_citation or "doi" in rec.er_citation.lower()
        if rec.fold_ko is not None:
            assert "PMID" in rec.fold_ko_citation or "doi" in rec.fold_ko_citation.lower()
        assert rec.notes.strip(), f"{rec.drug} must carry notes"

    # --- composition with the paracellular leak ----------------------------
    from ..reconstruct.barrier import L_PARACELLULAR_LEAK, effective_penetrance
    # a non-substrate reduces EXACTLY to the barrier module's expression
    for p in (0.05, 0.5, 0.9):
        for d in (0.0, 0.5, 1.0):
            assert (efflux_adjusted_penetrance(p, "cisplatin", d)
                    == float(effective_penetrance(p, d))), (p, d)
    # efflux hits the transcellular route only: at full disruption the
    # paracellular leak survives untouched
    hit = efflux_adjusted_penetrance(0.5, "eribulin", 1.0)
    assert abs(hit - (0.5 / 30.0 + L_PARACELLULAR_LEAK)) < 1e-12

    # --- the structural limit this module is subject to --------------------
    # A SCALAR per-drug factor cancels EXACTLY out of a tumour:normal ratio.
    f = efflux_factor("imatinib")
    e_tumour, e_normal = 0.5572, 0.9607           # GIST Krogh delivery efficiencies
    assert abs((f * e_tumour) / (f * e_normal) - e_tumour / e_normal) < 1e-15, \
        "efflux must cancel out of a within-drug tissue ratio -- it is scalar"

    # --- report -------------------------------------------------------------
    s = substrate_summary()
    print(f"EFFLUX_TABLE: {s['n_drugs']} drugs")
    print(f"  P-gp   yes/no/unknown : {s['pgp_yes']}/{s['pgp_no']}/{s['pgp_unknown']}")
    print(f"  BCRP   yes/no/unknown : {s['bcrp_yes']}/{s['bcrp_no']}/{s['bcrp_unknown']}")
    print(f"  substrate WITH a usable number    : {s['substrate_with_number']}")
    print(f"  substrate WITHOUT a usable number : {s['substrate_without_number']}")
    print(f"  non-substrate                     : {s['non_substrate']}")
    print()
    hdr = f"{'drug':14s} {'pgp':8s} {'bcrp':8s} {'ER':>7s} {'factor':>8s}  source     reason"
    print(hdr)
    print("-" * len(hdr))
    for name, rec in EFFLUX_TABLE.items():
        d = efflux_detail(name)
        er = d["efflux_ratio"]
        print(f"{name:14s} {rec.pgp:8s} {rec.bcrp:8s} "
              f"{(f'{er:.2f}' if er else '-'):>7s} {d['factor']:8.5f}  "
              f"{str(d['source'] or '-'):10s} {d['reason']}")
    print()
    print("in vitro / in vivo consistency (the only internal check on the form):")
    for form in EFFLUX_FORMS:
        rows = [(n, r.consistency(form)) for n, r in EFFLUX_TABLE.items()
                if r.consistency(form)]
        gm = 1.0
        for _, c in rows:
            v = c["in_vivo_over_in_vitro"]
            gm *= (v if v >= 1.0 else 1.0 / v)
        gm **= 1.0 / len(rows)
        for n, c in rows:
            print(f"  {form:10s} {n:12s} pred fold {c['fold_ko_predicted_from_er']:6.2f}"
                  f"  measured {c['fold_ko_measured']:6.2f}"
                  f"  ratio {c['in_vivo_over_in_vitro']:.2f}x")
        print(f"  {form:10s} -> geometric-mean |discrepancy| over n={len(rows)}: {gm:.2f}x")
    print("\nAll self-tests passed.")
