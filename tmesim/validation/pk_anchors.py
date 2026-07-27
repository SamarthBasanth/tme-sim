"""REAL measured tissue-pharmacokinetics anchors for MULTIPLE drugs.

This module is **pure data + accessors**. It does NOT import or touch the
core/engine/molecular modules and it does NOT run any simulation. Its job is to
carry a rigorously-cited compendium of *real, published, quantitative*
tissue-PK measurements so the physics engine can be validated against reality
for more than the single doxorubicin anchor used in
:mod:`tmesim.validation.pk_reality`.

Every number here is a real measurement traced to a peer-reviewed paper with a
PubMed ID (PMID) that was verified against PubMed. No estimates are presented as
data; where a value is a range reported jointly for several drugs, or is an
in-vivo whole-organ ratio rather than a spatial gradient, that is flagged in the
entry's ``confidence`` and ``notes``.

Two DISTINCT physical quantities live here (do not mix them):

* ``perivascular_half_distance_um`` -- the distance from a tumor blood vessel at
  which drug (immunofluorescence) intensity falls to HALF of its perivascular
  value. This is a spatial tissue-penetration gradient. It is the class of
  measurement directly comparable to the sim's perivascular penetration outputs
  (``penetration_length_um`` / the single-vessel ``1/e`` distance) -- **after**
  the half-vs-1/e conversion documented below.

* ``tissue_to_plasma_ratio`` / ``csf_to_plasma_ratio`` -- a dimensionless
  concentration (or AUC) ratio between a tissue/CSF compartment and plasma. This
  is a blood-brain / blood-tumor *partition* quantity. It is comparable to the
  sim's structural ``DrugParams.bbb_penetrance`` field, **not** to a penetration
  length.

Half-distance vs 1/e (the conversion the validator must apply)
--------------------------------------------------------------
The perivascular anchors report the distance ``d_half`` at which intensity falls
to 50 %. For an exponential perivascular profile ``C(r) ~ exp(-r/L)`` the
e-folding length ``L`` (the sim's ``penetration_length_um``) relates to the
measured half-distance by::

    d_half = L * ln(2)  ~=  0.693 * L        =>        L = d_half / ln(2) ~= 1.443 * d_half
    d_(1/e) = L                                       (the sim's "1/e distance")

So a measured half-distance of 40-50 um corresponds to an L-equivalent of
~58-72 um, and a measured 25-75 um corresponds to L ~ 36-108 um. Compare the
sim's ``L`` / ``1/e`` numbers to these *converted* values, not to the raw
half-distances. (The existing ``pk_reality.py`` compares to the raw 40-50 um,
which is slightly conservative -- it treats the half-distance as if it were the
1/e distance, making the anchor if anything a touch too small relative to L;
the L-equivalent ~58-72 um leaves its order-of-magnitude verdict unchanged.)

Units: length um, concentration ratios dimensionless (fraction, not %).
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

# Quantity vocabulary (the ``quantity`` field of each anchor).
Q_PERIVASCULAR_HALF = "perivascular_half_distance_um"
Q_TISSUE_PLASMA = "tissue_to_plasma_ratio"
Q_CSF_PLASMA = "csf_to_plasma_ratio"

# Which quantities are a spatial tissue-penetration gradient (comparable to the
# sim's penetration length / 1/e distance) vs a partition ratio (comparable to
# DrugParams.bbb_penetrance).
PENETRATION_QUANTITIES = frozenset({Q_PERIVASCULAR_HALF})
PARTITION_QUANTITIES = frozenset({Q_TISSUE_PLASMA, Q_CSF_PLASMA})

# ln(2): half-distance -> e-folding length conversion factor is 1/ln(2).
LN2 = math.log(2.0)


# ---------------------------------------------------------------------------
# The compendium. Each dict is one real measurement.
#   drug        : lowercase common drug name
#   quantity    : one of the Q_* constants above
#   value       : the measured number (or a (lo, hi) tuple for a reported range)
#   units       : human-readable units string
#   model       : the biological system the number was measured in
#   method      : how it was measured
#   citation    : full human-readable citation
#   pmid        : PubMed ID (verified against PubMed)
#   doi         : DOI (for attribution)
#   confidence  : "high" / "medium" / "low" -- our read of how directly the
#                 number supports a single-drug quantitative comparison
#   notes       : caveats / comparability guidance
# ---------------------------------------------------------------------------
PK_ANCHORS: List[Dict] = [
    {
        "drug": "doxorubicin",
        "quantity": Q_PERIVASCULAR_HALF,
        "value": (40.0, 50.0),
        "units": "um (distance to 50% of perivascular intensity)",
        "model": "murine 16C and EMT6 tumors + human PC-3 prostate xenografts (mice)",
        "method": "quantitative immunofluorescence of doxorubicin autofluorescence "
                  "vs distance from CD31+ blood vessels",
        "citation": "Primeau AJ, Rendon A, Hedley D, Lilge L, Tannock IF. The "
                    "distribution of the anticancer drug doxorubicin in relation to "
                    "blood vessels in solid tumors. Clin Cancer Res 2005;11(24 Pt 1):8782-8.",
        "pmid": "16361566",
        "doi": "10.1158/1078-0432.CCR-05-1664",
        "confidence": "high",
        "notes": "The canonical perivascular-penetration anchor (already used in "
                 "pk_reality.py). Reported as the distance at which concentration "
                 "falls to HALF its perivascular value -> L-equivalent ~58-72 um "
                 "(divide by ln2). Directly comparable to the sim's penetration "
                 "length / 1/e distance.",
    },
    {
        "drug": "doxorubicin",
        "quantity": Q_PERIVASCULAR_HALF,
        "value": (25.0, 75.0),
        "units": "um (distance to 50% intensity, 10 min post-injection)",
        "model": "various human tumor xenografts in mice",
        "method": "quantitative immunofluorescence vs distance from tumor blood vessels",
        "citation": "Patel KJ, Tredan O, Tannock IF. Distribution of the anticancer "
                    "drugs doxorubicin, mitoxantrone and topotecan in tumors and "
                    "normal tissues. Cancer Chemother Pharmacol 2013;72(1):127-38.",
        "pmid": "23680920",
        "doi": "10.1007/s00280-013-2176-z",
        "confidence": "high",
        "notes": "Independent replication of the ~40-50 um doxorubicin gradient in a "
                 "multi-drug study; the ~25-75 um half-distance range is reported "
                 "jointly for the three drugs at 10 min. Same half-vs-1/e conversion "
                 "applies (L-equivalent ~36-108 um).",
    },
    {
        "drug": "mitoxantrone",
        "quantity": Q_PERIVASCULAR_HALF,
        "value": (25.0, 75.0),
        "units": "um (distance to 50% intensity, 10 min post-injection)",
        "model": "various human tumor xenografts in mice",
        "method": "quantitative immunofluorescence vs distance from tumor blood vessels",
        "citation": "Patel KJ, Tredan O, Tannock IF. Distribution of the anticancer "
                    "drugs doxorubicin, mitoxantrone and topotecan in tumors and "
                    "normal tissues. Cancer Chemother Pharmacol 2013;72(1):127-38.",
        "pmid": "23680920",
        "doi": "10.1007/s00280-013-2176-z",
        "confidence": "medium",
        "notes": "The ~25-75 um half-distance is reported as a shared range for "
                 "doxorubicin/mitoxantrone/topotecan at 10 min, not as a "
                 "mitoxantrone-specific single value; treat as order-of-magnitude. "
                 "Perivascular tissue gradient (comparable to sim penetration length).",
    },
    {
        "drug": "topotecan",
        "quantity": Q_PERIVASCULAR_HALF,
        "value": (25.0, 75.0),
        "units": "um (distance to 50% intensity, 10 min post-injection)",
        "model": "various human tumor xenografts in mice",
        "method": "quantitative immunofluorescence vs distance from tumor blood vessels",
        "citation": "Patel KJ, Tredan O, Tannock IF. Distribution of the anticancer "
                    "drugs doxorubicin, mitoxantrone and topotecan in tumors and "
                    "normal tissues. Cancer Chemother Pharmacol 2013;72(1):127-38.",
        "pmid": "23680920",
        "doi": "10.1007/s00280-013-2176-z",
        "confidence": "medium",
        "notes": "Shared ~25-75 um half-distance range for the three drugs at 10 min; "
                 "not a topotecan-specific point value. Perivascular tissue gradient.",
    },
    {
        "drug": "temozolomide",
        "quantity": Q_TISSUE_PLASMA,
        "value": 0.178,
        "units": "dimensionless (brain interstitium AUC / plasma AUC)",
        "model": "human peritumoral brain tissue (primary + metastatic brain tumor "
                 "patients), single 150 mg/m2 oral dose",
        "method": "intracerebral microdialysis; temozolomide by tandem mass "
                  "spectrometry; AUC ratio over 24 h",
        "citation": "Portnow J, Badie B, Chen M, Liu A, Blanchard S, Synold TW. The "
                    "neuropharmacokinetics of temozolomide in patients with resectable "
                    "brain tumors. Clin Cancer Res 2009;15(22):7092-8.",
        "pmid": "19861433",
        "doi": "10.1158/1078-0432.CCR-09-1349",
        "confidence": "high",
        "notes": "Brain extracellular-fluid : plasma AUC ratio = 17.8% (~0.2), the "
                 "widely-cited temozolomide brain-penetration figure; mean peak brain "
                 "ECF concentration 0.6 +/- 0.3 ug/mL. This is a BBB/partition ratio -- "
                 "compare to DrugParams.bbb_penetrance, NOT to penetration length.",
    },
    {
        "drug": "methotrexate",
        "quantity": Q_CSF_PLASMA,
        "value": 0.013,
        "units": "dimensionless (steady-state CSF / plasma concentration ratio)",
        "model": "children with acute lymphoblastic leukemia WITHOUT CNS involvement "
                 "(n=49), high-dose intravenous methotrexate",
        "method": "paired CSF and plasma antifolate concentrations during IV infusion "
                  "(257 infusions)",
        "citation": "Morse M, Savitch J, Balis F, Miser J, Feusner J, Reaman G, "
                    "Poplack D, Bleyer A. Altered central nervous system pharmacology "
                    "of methotrexate in childhood leukemia. J Clin Oncol 1985;3(1):19-24.",
        "pmid": "3855310",
        "doi": "10.1200/JCO.1985.3.1.19",
        "confidence": "high",
        "notes": "Mean steady-state CSF:plasma ratio 0.013 (SD 0.01) in patients "
                 "without meningeal disease (rose ~12x to 0.157 with overt CNS "
                 "leukemia). A blood-CSF partition ratio -- compare to bbb_penetrance. "
                 "Illustrates a poorly-CNS-penetrant small molecule.",
    },
    {
        "drug": "cisplatin",
        "quantity": Q_CSF_PLASMA,
        "value": (0.025, 0.029),
        "units": "dimensionless (peak CSF total-Pt / plasma total-Pt ratio)",
        "model": "patients with solid-type metastatic brain tumor, IV cisplatin "
                 "(60 and 100 mg/m2)",
        "method": "total platinum by atomic absorption spectroscopy; CSF via Ommaya "
                  "reservoir in lateral ventricle; peak CSF/plasma ratio",
        "citation": "Nakagawa H, Fujita T, Kubo S, et al. Difference in CDDP "
                    "penetration into CSF between selective intraarterial chemotherapy "
                    "in malignant glioma and intravenous or intracarotid administration "
                    "in metastatic brain tumor. Cancer Chemother Pharmacol 1996;37(4):317-26.",
        "pmid": "8548876",
        "doi": "10.1007/s002800050391",
        "confidence": "medium",
        "notes": "Intact blood-CSF barrier (IV route, solid metastases): peak "
                 "CSF:plasma total-Pt ratio ~2.5% (60 mg/m2) to ~2.9% (100 mg/m2). "
                 "Higher ratios in the same paper reflect disrupted barrier / "
                 "intra-arterial or post-operative-cavity sampling and are not used as "
                 "the systemic-partition anchor. Blood-CSF partition ratio (compare to "
                 "bbb_penetrance).",
    },
]


# ---------------------------------------------------------------------------
# Supporting references (context / directional evidence, not point anchors).
# ---------------------------------------------------------------------------
SUPPORTING_REFERENCES: List[Dict] = [
    {
        "citation": "Minchinton AI, Tannock IF. Drug penetration in solid tumours. "
                    "Nat Rev Cancer 2006;6(8):583-92.",
        "pmid": "16862189",
        "doi": "10.1038/nrc1893",
        "role": "Review establishing that limited tissue penetration (steep "
                "perivascular gradients) is a general cause of chemoresistance.",
    },
    {
        "citation": "Tannock IF, Lee CM, Tunggal JK, Cowan DSM, Egorin MJ. Limited "
                    "penetration of anticancer drugs through tumor tissue. Clin Cancer "
                    "Res 2002;8(3):878-84.",
        "pmid": "11895922",
        "doi": "",
        "role": "Multicellular-layer (~200 um) method: cisplatin, etoposide, "
                "gemcitabine, paclitaxel and vinblastine all penetrate slowly. "
                "Directional (slow penetration) support; not a clean distance value.",
    },
    {
        "citation": "Fujimori K, Covell DG, Fletcher JE, Weinstein JN. A "
                    "modeling analysis of monoclonal antibody percolation through "
                    "tumors: a binding-site barrier. J Nucl Med 1990;31(7):1191-8.",
        "pmid": "2362198",
        "doi": "",
        "role": "The 'binding-site barrier' -- antibodies penetrate poorly. "
                "Directional anchor for the small-molecule-vs-antibody ordering check "
                "(PMID not re-verified in this pass; cited from pk_reality.py).",
    },
]


# ---------------------------------------------------------------------------
# Accessors (pure functions over the data above).
# ---------------------------------------------------------------------------
def anchors_for(drug: str) -> List[Dict]:
    """Return all anchor dicts for ``drug`` (case-insensitive), in list order."""
    key = str(drug).strip().lower()
    return [a for a in PK_ANCHORS if a["drug"] == key]


def drugs() -> List[str]:
    """Return the sorted unique list of drugs that have at least one anchor."""
    return sorted({a["drug"] for a in PK_ANCHORS})


def anchors_by_quantity(quantity: str) -> List[Dict]:
    """Return all anchors whose ``quantity`` equals ``quantity``."""
    return [a for a in PK_ANCHORS if a["quantity"] == quantity]


def penetration_anchors() -> List[Dict]:
    """Anchors that are spatial tissue-penetration gradients -- the ones directly
    comparable to the sim's penetration length / single-vessel 1/e distance."""
    return [a for a in PK_ANCHORS if a["quantity"] in PENETRATION_QUANTITIES]


def partition_anchors() -> List[Dict]:
    """Anchors that are tissue/CSF:plasma partition ratios -- comparable to the
    sim's DrugParams.bbb_penetrance, NOT to a penetration length."""
    return [a for a in PK_ANCHORS if a["quantity"] in PARTITION_QUANTITIES]


def value_range(anchor: Dict) -> tuple:
    """Return ``(lo, hi)`` for an anchor whose ``value`` is a scalar or a range."""
    v = anchor["value"]
    if isinstance(v, (tuple, list)):
        return float(v[0]), float(v[1])
    return float(v), float(v)


def half_distance_to_efold_length(d_half_um: float) -> float:
    """Convert a measured perivascular HALF-distance (distance to 50% intensity)
    to the exponential e-folding length L used by the sim, assuming an
    exponential perivascular profile: ``L = d_half / ln(2)``."""
    return float(d_half_um) / LN2


def efold_length_equivalent(anchor: Dict) -> Optional[tuple]:
    """For a perivascular half-distance anchor, return the ``(lo, hi)`` L-equivalent
    (the value directly comparable to the sim's penetration length). Returns
    ``None`` for non-penetration (partition-ratio) anchors."""
    if anchor["quantity"] not in PENETRATION_QUANTITIES:
        return None
    lo, hi = value_range(anchor)
    return half_distance_to_efold_length(lo), half_distance_to_efold_length(hi)


if __name__ == "__main__":
    print(f"{len(PK_ANCHORS)} anchors across {len(drugs())} drugs: {', '.join(drugs())}")
    print("-" * 78)
    for a in PK_ANCHORS:
        lo, hi = value_range(a)
        val = f"{lo:g}" if lo == hi else f"{lo:g}-{hi:g}"
        extra = ""
        if a["quantity"] in PENETRATION_QUANTITIES:
            elo, ehi = efold_length_equivalent(a)
            extra = f"  (L-equiv {elo:.0f}-{ehi:.0f} um)"
        print(f"  {a['drug']:14s} {a['quantity']:28s} {val:>10s}  "
              f"PMID {a['pmid']:>9s}  [{a['confidence']}]{extra}")
    print("-" * 78)
    print(f"penetration-gradient anchors (compare to sim L / 1e): "
          f"{len(penetration_anchors())}")
    print(f"partition-ratio anchors (compare to bbb_penetrance): "
          f"{len(partition_anchors())}")
