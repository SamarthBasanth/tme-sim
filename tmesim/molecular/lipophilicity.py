"""MEASURED octanol/water lipophilicity (logP) resolution for the molecular layer.

Closing the "logP gap", exactly the way :mod:`tmesim.molecular.potency` closes the
potency gap. Everything downstream of ``molecular.descriptors`` that touches
lipophilicity -- plasma free fraction, tissue binding, clearance, passive membrane
permeability and BBB penetrance -- was being driven by RDKit's **computed Crippen
cLogP**, an atom-contribution ESTIMATE, not a measurement. For several drugs that
estimate is wrong by more than two log units, which is more than the entire
dynamic range some of those estimators have.

The documented case (``results/binding_sink_refuted.md`` section 3):

    methotrexate   Crippen  +0.268   vs   experimental  -1.85   (2.1 log units)
    5-fluorouracil Crippen  -0.798   vs   experimental  -0.89   (0.09, fine)

and the model's own methotrexate:5-FU tissue-binding ratio collapses from 4.83x
to 1.51x on that substitution alone, against a MEASURED breast ratio of 1.02x.
Two other modules had already flagged the same pathology independently:
``permeability.membrane_permeability`` documents Crippen 4.59 vs experimental
~3.0 for imatinib, and ``permeability.bbb_penetrance`` documents Crippen +0.20
vs experimental ~-2.2 for cisplatin.

THIS IS A DATA-QUALITY FIX, NOT A TUNING KNOB
---------------------------------------------
Every value in :data:`EXPERIMENTAL_LOGP` is a published MEASUREMENT with a
citation. None of them was chosen, adjusted or selected to move any anchor, and
no value here was fitted to anything. Where a drug has no verifiable measured
logP, it is LEFT ON RDKIT and said so out loud (see :data:`NO_EXPERIMENTAL_LOGP`)
-- a partial table with honest gaps is the correct outcome; a complete table with
guessed numbers is not. Some of these substitutions make downstream predictions
WORSE; that is reported rather than hidden, because the input is still the right
input.

Three tiers, in strict order of trust (mirrors ``potency.resolve_ic50``)
-----------------------------------------------------------------------
  1. EXPERIMENTAL (:data:`EXPERIMENTAL_LOGP`) -- a published, measured
     octanol/water partition coefficient. tier = 'experimental'.
  2. COMPUTED -- RDKit's Crippen ``MolLogP`` (whatever ``descriptors`` produced).
     tier = 'rdkit_crippen'.
  3. FAIL LOUDLY -- if there is neither a measured value nor a finite computed
     one, :func:`resolve_logp` raises ``ValueError``. There is deliberately no
     silent numeric default: a wrong logP is exactly the defect this module
     exists to remove, so an absent one must not be papered over.

Evidence bar for tier 1 (applied uniformly; drugs that failed it are listed in
:data:`NO_EXPERIMENTAL_LOGP` with the reason)
---------------------------------------------------------------------------
ACCEPTED as measured:
  * Hansch C, Leo A, Hoekman D. *Exploring QSAR - Hydrophobic, Electronic, and
    Steric Constants.* Washington, DC: American Chemical Society, 1995 -- the
    standard experimental logP compilation. Surfaced through PubChem's
    Hazardous Substances Data Bank (HSDB) records as ``log Kow = X`` with a page
    number, and independently re-imported by DrugBank citing the same source.
    Each entry below carries its HSDB record id AND its page in Hansch.
  * A regulatory assessment report stating a measured partition coefficient for
    the active substance (used for exactly one drug, temozolomide, below).

REJECTED, and NOT used anywhere in this table:
  * Human Metabolome Database (HMDB) "LogP" -- these are ALOGPS *predictions*.
    PubChem lists them beside experimental values in the same LogP section,
    which is precisely how a computed number gets mistaken for a measured one.
  * XLogP3 / ChemAxon / Crippen / any other calculated descriptor. Substituting
    one computed logP for another computed logP is not a data-quality fix.
  * Database entries whose underlying source could not be established.

Units: logP is the base-10 octanol/water partition coefficient of the NEUTRAL
species, dimensionless. Note the corollary: it is a logP, not a logD -- none of
these numbers carries a pH, so the ionisation limitation already documented in
``permeability`` (no pKa anywhere in this layer) is unchanged by this module.
"""
from __future__ import annotations

import math
from typing import Dict, Optional

# ---------------------------------------------------------------------------
# Tier 1: MEASURED experimental octanol/water logP.
#
# Keys are the lowercased library drug names used by
# ``molecular.library.GBM_DRUGS`` and ``molecular.breast_library.BREAST_DRUGS``.
# Every entry carries: the value, the PubChem CID the value was read from, the
# primary source, and a note. NOTHING here is interpolated, averaged, rounded
# from a different compound, or inferred from a structural analogue.
# ---------------------------------------------------------------------------
_HANSCH = (
    "Hansch C, Leo A, Hoekman D. Exploring QSAR - Hydrophobic, Electronic, and "
    "Steric Constants. Washington, DC: American Chemical Society, 1995"
)

EXPERIMENTAL_LOGP: Dict[str, Dict] = {
    # --- GBM library (molecular.library.GBM_DRUGS) -------------------------
    "temozolomide": {
        "logp": -1.153,
        "pubchem_cid": 5394,
        "source": (
            "European Medicines Agency, Temozolomide Teva / Temomedac EPAR "
            "public assessment report, active-substance section: 'Its "
            "dissociation constant pKa has been found 0 and partition "
            "coefficient has been determined logP = -1.153.' "
            "https://www.ema.europa.eu/en/documents/assessment-report/"
            "temozolomide-teva-epar-public-assessment-report_en.pdf "
            "(also carried by DrugBank DB00853, which cites this same EPAR)"
        ),
        "note": (
            "Regulatory dossier measurement, not a compilation entry -- the "
            "only non-Hansch value in this table. Quoted verbatim from the "
            "EPAR PDF; the leading minus sign sits at a line break in that "
            "document and is easy to lose in text extraction, so it was "
            "re-checked against the layout-preserved text. RDKit Crippen gives "
            "-2.08, i.e. the estimate was already roughly right in sign and "
            "about 0.9 log units too polar."
        ),
    },
    "lomustine": {
        "logp": 2.83,
        "pubchem_cid": 3950,
        "source": f"{_HANSCH}, p. 64 (via PubChem/HSDB record 6519, 'log Kow = 2.83')",
        "note": "Independently re-imported by DrugBank DB01206 citing Hansch & Leo.",
    },
    "carmustine": {
        "logp": 1.53,
        "pubchem_cid": 2578,
        "source": f"{_HANSCH}, p. 14 (via PubChem/HSDB record 7761, 'log kow = 1.53')",
        "note": "Independently re-imported by DrugBank DB00262 citing Hansch et al. 1995.",
    },
    "procarbazine": {
        "logp": 0.06,
        "pubchem_cid": 4915,
        "source": f"{_HANSCH}, p. 106 (via PubChem/HSDB record 3250, 'log Kow = 0.06')",
        "note": (
            "Independently re-imported by DrugBank DB01168 citing Hansch et al. "
            "1995. RDKit Crippen gives 1.05, ~1 log unit too lipophilic."
        ),
    },
    "methotrexate": {
        "logp": -1.85,
        "pubchem_cid": 126941,
        "source": f"{_HANSCH}, p. 1071 (via PubChem/HSDB record 3123, 'log Kow = -1.85')",
        "note": (
            "THE case that motivated this module. RDKit Crippen gives +0.268 -- "
            "2.1 log units too lipophilic. Crippen's atom-contribution scheme "
            "cannot see that methotrexate's two carboxylates and pteridine "
            "make it extremely hydrophilic. Independently re-imported by "
            "DrugBank DB00563 citing Hansch et al. 1995. See "
            "results/binding_sink_refuted.md section 3."
        ),
    },
    "doxorubicin": {
        "logp": 1.27,
        "pubchem_cid": 31703,
        "source": f"{_HANSCH}, p. 187 (via PubChem/HSDB record 3070, 'log Kow = 1.27')",
        "note": (
            "Independently re-imported by DrugBank DB00997 citing Hansch et al. "
            "1995. NOT to be confused with the -0.5 PubChem also lists for this "
            "compound, which comes from HMDB and is an ALOGPS prediction. "
            "CONSEQUENCE TO BE AWARE OF: doxorubicin is the drug "
            "validation.calibration fits `alpha` on, so correcting its logP "
            "moves alpha. That is a downstream fact, not a reason to withhold "
            "the measurement -- see this module's docstring."
        ),
    },
    "etoposide": {
        "logp": 0.60,
        "pubchem_cid": 36462,
        "source": f"{_HANSCH}, p. 188 (via PubChem/HSDB record 6517, 'log Kow = 0.60')",
        "note": "Independently re-imported by DrugBank DB00773 citing Hansch et al. 1995.",
    },
    "vincristine": {
        "logp": 2.82,
        "pubchem_cid": 5978,
        "source": f"{_HANSCH}, p. 192 (via PubChem/HSDB record 3199, 'log Kow = 2.82')",
        "note": "Independently re-imported by DrugBank DB00541 citing Hansch et al. 1995.",
    },
    # --- Breast library (molecular.breast_library.BREAST_DRUGS) ------------
    "5-fluorouracil": {
        "logp": -0.89,
        "pubchem_cid": 3385,
        "source": f"{_HANSCH}, p. 7 (via PubChem/HSDB record 3228, 'log Kow = -0.89')",
        "note": (
            "The methotrexate comparator in results/binding_sink_refuted.md. "
            "RDKit Crippen gives -0.798, i.e. this one the estimator already "
            "had essentially right (0.09 log units) -- which is why the "
            "methotrexate:5-FU ratio was almost entirely a methotrexate error. "
            "Independently re-imported by DrugBank DB00544 citing Hansch 1995."
        ),
    },
    "cyclophosphamide": {
        "logp": 0.63,
        "pubchem_cid": 2907,
        "source": f"{_HANSCH}, p. 35 (via PubChem/HSDB record 3047, 'log Kow = 0.63')",
        "note": (
            "Corroborated independently by the ILO-WHO International Chemical "
            "Safety Card 0689, which also gives 0.63."
        ),
    },
}


# ---------------------------------------------------------------------------
# The honest gaps. Every library drug NOT in EXPERIMENTAL_LOGP is listed here
# with the reason it failed the evidence bar. These drugs keep RDKit's computed
# Crippen logP and are reported as tier 'rdkit_crippen'.
#
# This list is not decoration: it is the difference between "we checked and
# there is nothing" and "we didn't check". A drug that is in neither dict is a
# drug nobody looked at, and the self-test below enforces that both libraries
# are fully accounted for.
# ---------------------------------------------------------------------------
NO_EXPERIMENTAL_LOGP: Dict[str, str] = {
    "cisplatin": (
        "PubChem CID 5702198 carries NO experimental LogP section at all. The "
        "nearest measured source is Screnci D, McKeage MJ, Galettis P, et al., "
        "Br J Cancer 2000;82(4):966-72 (PMID 10732773, PMC2374394), which "
        "measured log P by octanol-aqueous shake flask for eight platinum "
        "drugs but publishes only the RANGE (-2.53 to -0.16) in the accessible "
        "text; the per-drug value for cisplatin is in a table that could not be "
        "retrieved. permeability.py already notes an 'experimental ~-2.2' for "
        "cisplatin, but that figure has no citation attached in this repo and "
        "is NOT promoted to a measured value here on that basis. Left on RDKit "
        "Crippen (+0.20), which is very likely ~2 log units too lipophilic -- "
        "a KNOWN, unfixed defect, not a resolved one."
    ),
    "irinotecan": (
        "PubChem CID 60838's only LogP entry (3.2) is from HMDB, i.e. an "
        "ALOGPS prediction, which this module rejects. No Hansch entry, no "
        "regulatory partition-coefficient statement found. Left on RDKit."
    ),
    "ribociclib": (
        "PubChem CID 44631912 carries no experimental LogP section. A modern "
        "(2017) approval with no entry in the Hansch compilation. Left on RDKit."
    ),
    "bevacizumab": (
        "A ~149 kDa IgG. An octanol/water partition coefficient is not a "
        "meaningful quantity for a folded protein and none is published. This "
        "drug resolves through descriptors._FALLBACK_TABLE (is_antibody=True), "
        "and every logP-consuming estimator downstream already hard-gates on "
        "is_antibody, so the value is inert. Left as-is."
    ),
    "paclitaxel": (
        "PubChem CID 36314 lists 3 from DrugBank DB01229 and 3 from HMDB. The "
        "HMDB value is an ALOGPS prediction; the DrugBank entry carries no "
        "citation, and DrugBank's page could not be retrieved to establish "
        "whether its logP is the experimental or the predicted property. A "
        "bare '3' agreeing exactly with a known prediction is not evidence of "
        "measurement. Left on RDKit rather than promoted on a guess."
    ),
    "docetaxel": (
        "PubChem CID 148124 lists 2.4 from DrugBank DB01248 and 2.4 from HMDB "
        "(ALOGPS prediction). Same problem as paclitaxel: the DrugBank entry "
        "is uncited and agrees exactly with the prediction, so its "
        "experimental status could not be established. Left on RDKit."
    ),
    "epirubicin": (
        "PubChem CID 41867 lists -0.5 from DrugBank DB00445 and -0.5 from "
        "HMDB (ALOGPS). Note that -0.5 is ALSO the HMDB prediction for "
        "doxorubicin, whose MEASURED Hansch value is +1.27; epirubicin is "
        "doxorubicin's 4'-epimer, so a genuine measurement should sit near "
        "+1.27, not -0.5. Both listed values therefore look predicted. "
        "Deliberately NOT assigned doxorubicin's 1.27 by analogy -- a "
        "stereoisomer's measured value is not this compound's measured value. "
        "Left on RDKit."
    ),
    "tamoxifen": (
        "PubChem CID 2733526's only LogP entry (7.1) is from HMDB, i.e. an "
        "ALOGPS prediction. No Hansch entry and no measured shake-flask value "
        "with a retrievable citation was found. Left on RDKit (6.00), which is "
        "the drug in these libraries where the estimator is most exposed -- it "
        "sits far above the logP ~3 clamp in permeability.membrane_permeability "
        "and has no measurement to check it against."
    ),
    "vinorelbine": (
        "PubChem CID 5311497 lists a bare 4 from DrugBank DB00361 with no "
        "citation and no corroborating source. A single-significant-figure, "
        "uncited database number does not meet the bar. Note the contrast with "
        "vincristine, the same drug class, which DOES have a Hansch value "
        "(2.82) and is in the measured tier. Left on RDKit."
    ),
    "palbociclib": (
        "PubChem CID 5330286 lists 0.99 from DrugBank DB09073 attributed to "
        "the 'IBRANCE (palbociclib) monograph'. That attribution could not be "
        "confirmed against the monograph text, and 0.99 for a compound whose "
        "basic piperazine has pKa 7.4 is much more consistent with a "
        "distribution coefficient (logD at ~pH 7.4) than with a neutral-species "
        "logP. A logD is NOT a logP and must not be fed to estimators that "
        "assume the neutral species. Left on RDKit."
    ),
}


def experimental_logp(drug_name: str) -> Optional[Dict]:
    """Return the curated experimental logP record for a drug, or ``None``.

    Case-insensitive on ``drug_name``. The returned dict is a copy, so callers
    cannot mutate the table.
    """
    key = (drug_name or "").strip().lower()
    rec = EXPERIMENTAL_LOGP.get(key)
    return dict(rec) if rec is not None else None


def resolve_logp(drug_name: str, computed_logp: Optional[float] = None) -> Dict:
    """Resolve a logP for a drug using the best available evidence tier.

    Order of trust: measured experimental -> RDKit computed -> raise.

    Parameters
    ----------
    drug_name : str
        Library drug name; case-insensitive. May be empty/unknown, in which case
        only the computed tier is available.
    computed_logp : float or None
        RDKit's Crippen ``MolLogP`` for this molecule, if it could be computed.

    Returns
    -------
    dict
        Always ``{logp, tier, source, confidence, note, computed_logp,
        experimental_logp, delta}``.

        * ``tier`` is ``'experimental'`` or ``'rdkit_crippen'``.
        * ``computed_logp`` is the RDKit value that was superseded (or used),
          kept so the substitution is always inspectable after the fact.
        * ``experimental_logp`` is the measured value or ``None``.
        * ``delta`` is ``experimental - computed`` when both exist, else ``None``.

    Raises
    ------
    ValueError
        If the drug has no experimental entry AND no finite computed logP. This
        is the deliberate "fail loudly" tier: there is no silent default,
        because an unnoticed wrong logP is the exact defect this module removes.
    """
    comp: Optional[float]
    try:
        comp = float(computed_logp) if computed_logp is not None else None
    except (TypeError, ValueError):
        comp = None
    if comp is not None and not math.isfinite(comp):
        comp = None

    rec = experimental_logp(drug_name)
    if rec is not None:
        exp = float(rec["logp"])
        return {
            "logp": exp,
            "tier": "experimental",
            "confidence": "measured",
            "source": rec["source"],
            "note": rec.get("note", ""),
            "pubchem_cid": rec.get("pubchem_cid"),
            "computed_logp": comp,
            "experimental_logp": exp,
            "delta": (exp - comp) if comp is not None else None,
        }

    if comp is not None:
        key = (drug_name or "").strip().lower()
        gap = NO_EXPERIMENTAL_LOGP.get(key)
        return {
            "logp": comp,
            "tier": "rdkit_crippen",
            "confidence": "computed_ESTIMATE",
            "source": "RDKit Crippen MolLogP (atom-contribution estimate)",
            "note": (
                gap if gap else
                "no curated experimental logP was sought or found for this "
                "molecule; it is not one of the GBM/breast library drugs the "
                "EXPERIMENTAL_LOGP table covers."
            ),
            "pubchem_cid": None,
            "computed_logp": comp,
            "experimental_logp": None,
            "delta": None,
        }

    raise ValueError(
        f"no logP available for {drug_name!r}: it has no entry in "
        "EXPERIMENTAL_LOGP and no finite RDKit Crippen value was supplied. "
        "Refusing to substitute a default -- a silently-wrong logP drives "
        "binding, clearance, free fraction, permeability and BBB penetrance."
    )


def provenance_string(res: Dict) -> str:
    """One-line provenance for a :func:`resolve_logp` result (for profiles)."""
    if res["tier"] == "experimental":
        d = res.get("delta")
        moved = (f"; supersedes RDKit Crippen {res['computed_logp']:+.3f} "
                 f"(delta {d:+.3f})") if d is not None else ""
        return (f"logP={res['logp']:+.3f} MEASURED (tier: experimental); "
                f"{res['source']}{moved}; confidence: measured")
    return (f"logP={res['logp']:+.3f} RDKit Crippen MolLogP (tier: "
            f"rdkit_crippen); confidence: computed ESTIMATE -- {res['note']}")


if __name__ == "__main__":
    from .breast_library import BREAST_DRUGS
    from .library import GBM_DRUGS

    print("=" * 78)
    print("lipophilicity.resolve_logp self-test: measured logP > RDKit > raise")
    print("=" * 78)

    # 1. Both libraries must be FULLY accounted for: every drug is either in the
    #    measured table or in the explicit gap list. Silence is not allowed.
    all_drugs = set(GBM_DRUGS) | set(BREAST_DRUGS)
    accounted = set(EXPERIMENTAL_LOGP) | set(NO_EXPERIMENTAL_LOGP)
    missing = sorted(all_drugs - accounted)
    assert not missing, f"library drugs neither measured nor explained: {missing}"
    stray = sorted(accounted - all_drugs)
    assert not stray, f"table entries for drugs in neither library: {stray}"
    overlap = sorted(set(EXPERIMENTAL_LOGP) & set(NO_EXPERIMENTAL_LOGP))
    assert not overlap, f"drug both measured and 'no data': {overlap}"

    n_meas = len(all_drugs & set(EXPERIMENTAL_LOGP))
    print(f"\n{len(all_drugs)} library drugs: {n_meas} measured, "
          f"{len(all_drugs) - n_meas} left on RDKit (each with a stated reason).")

    # 2. Every measured entry carries a real citation and a CID.
    for name, rec in EXPERIMENTAL_LOGP.items():
        assert isinstance(rec["logp"], float), name
        assert math.isfinite(rec["logp"]), name
        assert rec.get("pubchem_cid"), f"{name}: no PubChem CID"
        assert len(rec.get("source", "")) > 40, f"{name}: citation too thin"
    print("every measured entry has a CID and a citation  OK")

    # 3. Tier behaviour.
    r = resolve_logp("methotrexate", computed_logp=0.268)
    assert r["tier"] == "experimental" and r["logp"] == -1.85, r
    assert abs(r["delta"] - (-1.85 - 0.268)) < 1e-12
    print("\nmethotrexate:", provenance_string(r))

    r = resolve_logp("paclitaxel", computed_logp=3.196)
    assert r["tier"] == "rdkit_crippen" and r["logp"] == 3.196, r
    print("\npaclitaxel  :", provenance_string(r)[:120], "...")

    # Unknown drug with a computed value still resolves (computed tier).
    r = resolve_logp("some_new_molecule", computed_logp=1.5)
    assert r["tier"] == "rdkit_crippen" and r["logp"] == 1.5

    # 4. FAIL LOUDLY: nothing measured, nothing computed.
    for bad in (None, float("nan"), float("inf")):
        try:
            resolve_logp("some_new_molecule", computed_logp=bad)
        except ValueError:
            pass
        else:
            raise SystemExit(f"expected ValueError for computed_logp={bad!r}")
    print("\nfail-loudly tier raises when neither tier is available  OK")

    # 5. The table is a lookup, not a mutable global.
    rec = experimental_logp("methotrexate")
    rec["logp"] = 999.0
    assert EXPERIMENTAL_LOGP["methotrexate"]["logp"] == -1.85, "table was mutated"

    print("\nAll lipophilicity self-tests passed.")
