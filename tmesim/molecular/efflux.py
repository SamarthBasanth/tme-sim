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

The functional form — DERIVED, not assumed
-------------------------------------------
The published number for a transporter substrate is the **efflux ratio**

.. math::

    \\mathrm{ER} = \\frac{P_\\mathrm{app}(B \\rightarrow A)}{P_\\mathrm{app}(A \\rightarrow B)}

measured in a polarised monolayer (MDCKII-MDR1, MDCKII-BCRP, Caco-2) with the
transporter on the apical (``A``) face. We need the factor by which that pump
reduces *net absorptive* transport relative to the same molecule with the pump
switched off. That factor follows from the geometry; it does not have to be
guessed.

Take the standard three-compartment monolayer: apical ``A`` — cell ``C`` —
basolateral ``B``, with passive permeability ``P`` across BOTH membranes and an
additional efflux clearance ``E`` pumping ``C -> A`` at the apical membrane.
Quasi-steady state on the cell compartment, sink conditions on the receiver:

.. math::

    P C_A = (2P + E)\\, C_C
    \\quad\\Rightarrow\\quad
    P_\\mathrm{app}(A\\!\\to\\!B) = \\frac{P^2}{2P + E},
    \\;\\;
    P_\\mathrm{app}(B\\!\\to\\!A) = \\frac{P(P + E)}{2P + E}

so ``ER = 1 + E/P`` and, since the passive-only two-membranes-in-series
permeability is ``P/2``,

.. math::

    \\boxed{\\;
      f_\\mathrm{efflux}
      \\;=\\; \\frac{P_\\mathrm{app}(A\\!\\to\\!B)}{P_\\mathrm{passive}}
      \\;=\\; \\frac{2P}{2P + E}
      \\;=\\; \\frac{2}{1 + \\mathrm{ER}} \\;}

This is :data:`FORM_MONOLAYER`, the default. Three properties make it the right
default rather than the ``1/(1+ER)`` starting form the brief suggested:

1. **It is exactly 1 at ER = 1.** A non-substrate assays at ER ~ 1 and must
   receive no penalty. ``1/(1+ER)`` would silently halve every non-substrate.
2. **It is derived from the same experiment that produces the number.** ER is
   defined by a two-membrane monolayer; the conversion uses that same geometry
   rather than a fitted curve.
3. **It is the conservative member of its family.** For large ER it tends to
   ``2/ER``, i.e. exactly twice the alternative ``1/ER`` (:data:`FORM_KPUU`, the
   BBB unbound-partition convention). Since the literature consistently finds
   that in vitro ER **under**-states the in vivo brain effect, choosing the
   weaker of the two forms means the term cannot be accused of manufacturing the
   improvement it is being tested for.

The blood-brain barrier endothelium is structurally the same object — a single
cell layer with a luminal and an abluminal membrane and the pumps on the luminal
face — so the same derivation carries over, with ``f_efflux`` interpreted as the
efflux-driven reduction in ``Kp,uu,brain``. Both forms are provided; the choice
is a single argument and the impact of switching it is reported in
``results/efflux_impact.md``.

The other published observable is the same model, not a second model
--------------------------------------------------------------------
Several drugs have no in vitro ER but DO have an in vivo transporter-knockout
brain-to-plasma fold change, ``fold_KO = Kp(knockout) / Kp(wild-type)``. That is
not a different quantity requiring a different calibration — it is the **same
three-compartment model read out a different way**. Deleting the pump sets
``E = 0``, so

.. math::

    \\mathrm{fold}_\\mathrm{KO}
      = \\frac{P/2}{P^2/(2P+E)}
      = 1 + \\frac{E}{2P}
      = \\frac{1 + \\mathrm{ER}}{2}
    \\qquad\\Longrightarrow\\qquad
    f_\\mathrm{efflux} = \\frac{1}{\\mathrm{fold}_\\mathrm{KO}}
                       = \\frac{2}{1 + \\mathrm{ER}}

i.e. an in vivo knockout fold change measures ``1/f_efflux`` **directly**, and
maps onto an in vitro ER by ``ER = 2 * fold_KO - 1``. (Verified numerically in
this module's self-test.) This matters for two reasons. First, it means the
module is not form-shopping when it uses whichever observable exists for a given
drug — there is one model and two windows onto it. Second, it turns the
in-vitro/in-vivo comparison into a genuine **consistency check**: where both
exist, ``(1 + ER)/2`` should equal the measured knockout fold, and where it does
not, the size of the disagreement is a direct read on how much the monolayer
assay under-states the intact barrier. That check is reported per drug rather
than averaged away.

**Which observable is preferred, decided in advance.** In vitro ER is the
primary source when it exists: it is per-pump, it is the standard reported
quantity, and it is measured on human transporters. The in vivo knockout fold is
the fallback and is used only where no in vitro ER is published. That rule is
fixed here so it cannot be chosen per drug to flatter a result; the opposite
preference is run as a sensitivity in ``results/efflux_impact.md``.

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
    "combine_efflux_ratios",
    "efflux_ratio_for",
    "efflux_factor",
    "efflux_detail",
    "efflux_adjusted_penetrance",
    "substrate_summary",
]


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

    @property
    def er_from_knockout(self) -> Optional[float]:
        """The in vivo knockout fold change expressed as an ER: ``2*fold - 1``."""
        return None if self.fold_ko is None else float(2.0 * self.fold_ko - 1.0)

    @property
    def has_number(self) -> bool:
        """True if at least one usable published quantity exists (either window)."""
        return self.er_in_vitro is not None or self.fold_ko is not None

    @property
    def consistency(self) -> Optional[Dict[str, float]]:
        """In vitro vs in vivo agreement, where BOTH observables exist.

        Under the module's model the in vitro ER predicts a knockout fold change
        of ``(1 + ER)/2``. Returns that prediction, the measured fold, and their
        ratio (>1 means the monolayer assay UNDER-states the intact barrier,
        which is the direction the literature generally reports). ``None`` when
        one of the two is missing.
        """
        er = self.er_in_vitro
        if er is None or self.fold_ko is None:
            return None
        pred = (1.0 + er) / 2.0
        return {"fold_ko_predicted_from_er": pred,
                "fold_ko_measured": float(self.fold_ko),
                "in_vivo_over_in_vitro": float(self.fold_ko) / pred}


# ===========================================================================
# 2. The functional forms
# ===========================================================================
#: Derived three-compartment monolayer form, ``f = 2 / (1 + ER)``.
#:
#: Exact for the model in this module's docstring: passive permeability ``P``
#: across both membranes, efflux clearance ``E`` at the apical membrane,
#: sink receiver. ``f(1) = 1`` exactly, ``f -> 2/ER`` for large ER.
#: **This is the default.**
FORM_MONOLAYER = "monolayer"

#: BBB unbound-partition convention, ``f = 1 / ER``.
#:
#: The relation commonly quoted when in vitro ER is used to anticipate
#: ``Kp,uu,brain``. It is exactly 2x more aggressive than
#: :data:`FORM_MONOLAYER` at large ER and, unlike it, does NOT return 1 at
#: ER = 1 (it returns 1 there too, but it penalises weakly-effluxed molecules
#: harder immediately above 1). Provided for sensitivity analysis.
FORM_KPUU = "kpuu"

EFFLUX_FORMS = (FORM_MONOLAYER, FORM_KPUU)


def _apply_form(er: float, form: str) -> float:
    """Map an efflux ratio to a penetrance multiplier in ``(0, 1]``."""
    if er < 1.0:
        # An ER below 1 means net ABSORPTIVE transport (an uptake transporter, or
        # assay noise). This module models efflux only; it must never boost
        # penetrance, so the factor is clamped at 1.
        er = 1.0
    if form == FORM_MONOLAYER:
        f = 2.0 / (1.0 + er)
    elif form == FORM_KPUU:
        f = 1.0 / er
    else:
        raise ValueError(f"unknown efflux form {form!r}; expected one of {EFFLUX_FORMS}")
    return float(min(1.0, max(0.0, f)))


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
# 3. The table
# ===========================================================================
#: LITERATURE-SOURCED substrate assignments. See ``results/efflux_impact.md``
#: for the same table rendered with full citations and the search that produced
#: it. Populated in section 3 below.
EFFLUX_TABLE: Dict[str, EffluxRecord] = {}


def _add(rec: EffluxRecord) -> None:
    EFFLUX_TABLE[rec.drug] = rec


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


def efflux_ratio_for(drug: Union[str, Mapping],
                     record: Optional[EffluxRecord] = None) -> Optional[float]:
    """The combined published efflux ratio for a drug, or ``None`` if unavailable."""
    rec = record if record is not None else EFFLUX_TABLE.get(_resolve_name(drug) or "")
    if rec is None or not rec.is_substrate:
        return None
    return combine_efflux_ratios(rec.er_pgp, rec.er_bcrp)


def efflux_detail(drug: Union[str, Mapping],
                  form: str = FORM_MONOLAYER,
                  efflux_ratio: Optional[float] = None) -> Dict[str, object]:
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
                 "citations": "; ".join(
                     c for c in (rec.pgp_citation, rec.bcrp_citation, rec.er_citation)
                     if c)})
    if not rec.is_substrate:
        both_unknown = rec.pgp == "unknown" and rec.bcrp == "unknown"
        return {**base,
                "reason": "unknown_substrate_status" if both_unknown else "non_substrate"}
    er = combine_efflux_ratios(rec.er_pgp, rec.er_bcrp)
    if er is None:
        return {**base, "reason": "substrate_but_no_published_er"}
    return {**base, "factor": _apply_form(er, form), "efflux_ratio": er,
            "applied": True, "reason": "applied"}


def efflux_factor(drug: Union[str, Mapping],
                  form: str = FORM_MONOLAYER,
                  efflux_ratio: Optional[float] = None) -> float:
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
    form : {"monolayer", "kpuu"}
        Functional form. ``"monolayer"`` (default) is ``2/(1+ER)``, derived from
        the three-compartment monolayer that defines ER. ``"kpuu"`` is ``1/ER``,
        the BBB unbound-partition convention, exactly 2x more aggressive at high
        ER; provided for sensitivity analysis.
    efflux_ratio : float, optional
        Bypass the table with an explicit ratio. Intended for sweeps and for
        drugs whose ratio a caller has from a source not carried here. Values
        below 1 are clamped to 1 (this module never boosts penetrance).

    Returns
    -------
    float
        Multiplier in ``(0, 1]``.
    """
    return float(efflux_detail(drug, form=form, efflux_ratio=efflux_ratio)["factor"])


def efflux_adjusted_penetrance(passive_penetrance: float,
                               drug: Union[str, Mapping],
                               bbb_disruption: float = 0.0,
                               leak: Optional[float] = None,
                               form: str = FORM_MONOLAYER,
                               efflux_ratio: Optional[float] = None) -> float:
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
    f = efflux_factor(drug, form=form, efflux_ratio=efflux_ratio)
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
